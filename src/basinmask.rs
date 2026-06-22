//! WOA 0.25° basin mask → per-cell basin id on the mapping grid.
//!
//! Mirrors `WMO2024_interp_basins.m`: parse the text table (`Latitude, Longitude,
//! Basin_<depth>m`), take the surface column (`Basin_0m`), wrap longitudes `<20 → +360`,
//! and assign each 1° grid cell the basin of its nearest mask point by Euclidean distance
//! in (lat, lon) degrees (MATLAB `knnsearch`).
//!
//! The mask points form a regular 0.25° lattice (ocean only), so instead of a kd-tree
//! (which struggles with the ~1440 points sharing each latitude) we index into a dense
//! lattice and find the nearest *filled* node by an expanding ring search with an exact
//! Euclidean stop test. This reproduces nearest-neighbour exactly and removes a dependency.

use anyhow::{anyhow, bail, Result};
use ndarray::Array2;
use std::path::Path;

use crate::GridDef;

/// Column index of `Basin_0m` (surface) in the `.msk` table: Lat, Lon, Basin_0m, …
const SURFACE_COL: usize = 2;
/// Mask lattice spacing (degrees).
const STEP: f64 = 0.25;
/// Empty-cell sentinel in the dense lattice.
const EMPTY: i16 = i16::MIN;

struct Lattice {
    grid: Vec<i16>, // nlat_l * nlon_l, row-major; EMPTY where no mask point
    lat_min: f64,
    lon_min: f64,
    nlat_l: usize,
    nlon_l: usize,
}

impl Lattice {
    fn idx(&self, il: usize, jl: usize) -> i16 {
        self.grid[il * self.nlon_l + jl]
    }

    /// Nearest filled node to (qlat, qlon), exact Euclidean.
    fn nearest(&self, qlat: f64, qlon: f64) -> i16 {
        let ci = ((qlat - self.lat_min) / STEP).round() as i64;
        let cj = ((qlon - self.lon_min) / STEP).round() as i64;
        let mut best_d2 = f64::INFINITY;
        let mut best_b = 0i16;
        let mut r: i64 = 0;
        let rmax = (self.nlat_l + self.nlon_l) as i64 + 2;
        loop {
            for di in -r..=r {
                for dj in -r..=r {
                    if di.abs() != r && dj.abs() != r {
                        continue; // ring border only
                    }
                    let il = ci + di;
                    let jl = cj + dj;
                    if il < 0 || jl < 0 || il as usize >= self.nlat_l || jl as usize >= self.nlon_l {
                        continue;
                    }
                    let b = self.idx(il as usize, jl as usize);
                    if b != EMPTY {
                        let plat = self.lat_min + il as f64 * STEP;
                        let plon = self.lon_min + jl as f64 * STEP;
                        let d2 = (plat - qlat).powi(2) + (plon - qlon).powi(2);
                        if d2 < best_d2 {
                            best_d2 = d2;
                            best_b = b;
                        }
                    }
                }
            }
            // any node beyond radius r is ≥ (r+0.5)*STEP away; stop once that exceeds best.
            if best_d2.is_finite() && (r as f64 + 0.5) * STEP > best_d2.sqrt() {
                break;
            }
            r += 1;
            if r > rmax {
                break; // safety (won't trigger for a non-empty lattice)
            }
        }
        best_b
    }
}

fn build_lattice(path: impl AsRef<Path>) -> Result<Lattice> {
    let text = std::fs::read_to_string(path.as_ref())?;
    let mut lines = text.lines();
    let _comment = lines.next(); // "#..."
    let _header = lines.next();

    let mut pts: Vec<(f64, f64, i16)> = Vec::new();
    let (mut lat_min, mut lat_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut lon_min, mut lon_max) = (f64::INFINITY, f64::NEG_INFINITY);

    for line in lines {
        if line.is_empty() {
            continue;
        }
        let mut it = line.split(',');
        let lat: f64 = it.next().ok_or_else(|| anyhow!("missing lat"))?.trim().parse()?;
        let mut lon: f64 = it.next().ok_or_else(|| anyhow!("missing lon"))?.trim().parse()?;
        let surf = it.nth(SURFACE_COL - 2).unwrap_or("").trim();
        if surf.is_empty() {
            continue;
        }
        let basin: i16 = surf.parse()?;
        if lon < 20.0 {
            lon += 360.0; // match basinLong wrap
        }
        lat_min = lat_min.min(lat);
        lat_max = lat_max.max(lat);
        lon_min = lon_min.min(lon);
        lon_max = lon_max.max(lon);
        pts.push((lat, lon, basin));
    }
    if pts.is_empty() {
        bail!("no basin rows parsed");
    }

    let nlat_l = ((lat_max - lat_min) / STEP).round() as usize + 1;
    let nlon_l = ((lon_max - lon_min) / STEP).round() as usize + 1;
    let mut grid = vec![EMPTY; nlat_l * nlon_l];
    for (lat, lon, basin) in pts {
        let il = ((lat - lat_min) / STEP).round() as usize;
        let jl = ((lon - lon_min) / STEP).round() as usize;
        grid[il * nlon_l + jl] = basin;
    }
    Ok(Lattice { grid, lat_min, lon_min, nlat_l, nlon_l })
}

/// Parse the basin mask and assign a basin id to every `[lat, lon]` cell of `grid`.
pub fn read_basin_id(path: impl AsRef<Path>, grid: &GridDef) -> Result<Array2<i16>> {
    let lat = build_lattice(path)?;
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    let mut out = Array2::<i16>::zeros((nlat, nlon));
    for j in 0..nlat {
        let qlat = grid.lat[j];
        for i in 0..nlon {
            out[[j, i]] = lat.nearest(qlat, grid.lon[i]);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    // Set OHC_BASINMASK to basinmask_04.msk. Skipped if unset.
    fn path() -> Option<PathBuf> {
        std::env::var_os("OHC_BASINMASK").map(PathBuf::from)
    }

    #[test]
    fn basin_assignment_matches_reference() {
        let Some(p) = path() else { return };
        let grid = GridDef::mapping();
        let b = read_basin_id(p, &grid).unwrap();
        assert_eq!(b.dim(), (180, 360));
        // ground-truth (brute-force nearest in the Python reference):
        // index j=lat+89.5, i=lon-20.5
        assert_eq!(b[[100, 200]], 2); // central Pacific  (10.5, 220.5)
        assert_eq!(b[[120, 310]], 1); // N Atlantic       (30.5, 330.5)
        assert_eq!(b[[69, 60]], 3); //  Indian           (-20.5, 80.5)
        assert_eq!(b[[130, 130]], 2); // NW Pacific       (40.5, 150.5)
    }
}

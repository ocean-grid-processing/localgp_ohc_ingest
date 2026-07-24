//! Minimal classic-NetCDF reader for `etopo60.cdf`.
//!
//! `etopo60.cdf` is classic NetCDF (CDF v1, big-endian): vars `ETOPO60X` (360, degrees_east
//! 20.5…379.5), `ETOPO60Y` (180, degrees_north −89.5…89.5), `ROSE` (relief, meters, negative
//! below sea level) with dims `[ETOPO60Y, ETOPO60X]` → read row-major as `[lat, lon]`. The grid
//! is bit-identical to the mapping grid, so we assert that rather than regridding.
//!
//! Purpose-built (no general NetCDF support), pure Rust. Handles CDF v1 (u32 offsets) and
//! v2 (i64 offsets).

use anyhow::{anyhow, bail, Result};
use byteorder::{BigEndian, ReadBytesExt};
use ndarray::Array2;
use std::collections::HashMap;
use std::path::Path;

use crate::GridDef;

const NC_FLOAT: u32 = 5;
const NC_DOUBLE: u32 = 6;
const NC_DIMENSION: u32 = 0x0A;
const NC_VARIABLE: u32 = 0x0B;
const NC_ATTRIBUTE: u32 = 0x0C;

fn nctype_size(t: u32) -> usize {
    match t {
        1 | 2 => 1, // byte/char
        3 => 2,     // short
        4 | 5 => 4, // int / float
        6 => 8,     // double
        _ => 0,
    }
}

struct VarMeta {
    nctype: u32,
    shape: Vec<usize>,
    begin: u64,
}

struct Cdf {
    buf: Vec<u8>,
    vars: HashMap<String, VarMeta>,
}

struct Cursor<'a> {
    buf: &'a [u8],
    off: usize,
}
impl<'a> Cursor<'a> {
    fn u32(&mut self) -> Result<u32> {
        let v = (&self.buf[self.off..self.off + 4]).read_u32::<BigEndian>()?;
        self.off += 4;
        Ok(v)
    }
    fn u64(&mut self) -> Result<u64> {
        let v = (&self.buf[self.off..self.off + 8]).read_u64::<BigEndian>()?;
        self.off += 8;
        Ok(v)
    }
    fn name(&mut self) -> Result<String> {
        let n = self.u32()? as usize;
        let s = String::from_utf8_lossy(&self.buf[self.off..self.off + n]).to_string();
        self.off += n + ((4 - n % 4) % 4); // pad to 4
        Ok(s)
    }
    /// skip an att_list (tag, nelems, then name/type/nelems/values*)
    fn skip_att_list(&mut self) -> Result<()> {
        let tag = self.u32()?;
        let n = self.u32()?;
        if tag != NC_ATTRIBUTE && !(tag == 0 && n == 0) {
            bail!("expected attribute list, got tag {tag}");
        }
        for _ in 0..n {
            let _name = self.name()?;
            let t = self.u32()?;
            let ne = self.u32()? as usize;
            let nb = ne * nctype_size(t);
            self.off += nb + ((4 - nb % 4) % 4);
        }
        Ok(())
    }
}

impl Cdf {
    fn open(path: impl AsRef<Path>) -> Result<Self> {
        let buf = std::fs::read(path.as_ref())?;
        if &buf[0..3] != b"CDF" {
            bail!("not a classic NetCDF file");
        }
        let version = buf[3];
        let wide = version >= 2; // v2/v5 use 64-bit begin offsets
        let mut c = Cursor { buf: &buf, off: 4 };
        let _numrecs = c.u32()?;

        // dim_list
        let mut dim_sizes = Vec::new();
        let tag = c.u32()?;
        let nd = c.u32()?;
        if tag == NC_DIMENSION {
            for _ in 0..nd {
                let _nm = c.name()?;
                dim_sizes.push(c.u32()? as usize);
            }
        } else if !(tag == 0 && nd == 0) {
            bail!("expected dimension list, got tag {tag}");
        }

        // global att_list
        c.skip_att_list()?;

        // var_list
        let mut vars = HashMap::new();
        let tag = c.u32()?;
        let nv = c.u32()?;
        if tag == NC_VARIABLE {
            for _ in 0..nv {
                let name = c.name()?;
                let nda = c.u32()? as usize;
                let mut shape = Vec::with_capacity(nda);
                for _ in 0..nda {
                    let dimid = c.u32()? as usize;
                    shape.push(dim_sizes[dimid]);
                }
                c.skip_att_list()?;
                let nctype = c.u32()?;
                let _vsize = c.u32()?;
                let begin = if wide { c.u64()? } else { c.u32()? as u64 };
                vars.insert(name, VarMeta { nctype, shape, begin });
            }
        } else if !(tag == 0 && nv == 0) {
            bail!("expected variable list, got tag {tag}");
        }

        let off = c.off;
        let _ = off;
        Ok(Cdf { buf, vars })
    }

    fn read_var_f64(&self, name: &str) -> Result<(Vec<f64>, Vec<usize>)> {
        let v = self
            .vars
            .get(name)
            .ok_or_else(|| anyhow!("variable {name} not found"))?;
        let count: usize = v.shape.iter().product();
        let start = v.begin as usize;
        let sz = nctype_size(v.nctype);
        let mut slice = &self.buf[start..start + count * sz];
        let mut out = Vec::with_capacity(count);
        match v.nctype {
            NC_FLOAT => {
                for _ in 0..count {
                    out.push(slice.read_f32::<BigEndian>()? as f64);
                }
            }
            NC_DOUBLE => {
                for _ in 0..count {
                    out.push(slice.read_f64::<BigEndian>()?);
                }
            }
            other => bail!("unsupported nc type {other} for {name}"),
        }
        Ok((out, v.shape.clone()))
    }
}

/// Read `etopo60.cdf` and return bathymetry as `[lat, lon]` (meters, negative below sea
/// level), after asserting its grid equals the mapping grid.
pub fn read_etopo(path: impl AsRef<Path>, grid: &GridDef) -> Result<Array2<f64>> {
    let cdf = Cdf::open(path)?;

    let (x, _) = cdf.read_var_f64("ETOPO60X")?;
    let (y, _) = cdf.read_var_f64("ETOPO60Y")?;
    assert_grid(&x, &grid.lon, "ETOPO60X/lon")?;
    assert_grid(&y, &grid.lat, "ETOPO60Y/lat")?;

    let (rose, shape) = cdf.read_var_f64("ROSE")?;
    // ROSE dims are [ETOPO60Y, ETOPO60X] = [lat, lon]; row-major read is already [lat, lon].
    if shape != vec![grid.nlat(), grid.nlon()] {
        bail!("ROSE shape {:?} != [nlat, nlon] = [{}, {}]", shape, grid.nlat(), grid.nlon());
    }
    Array2::from_shape_vec((grid.nlat(), grid.nlon()), rose)
        .map_err(|e| anyhow!("reshape ROSE: {e}"))
}

fn assert_grid(got: &[f64], want: &[f64], label: &str) -> Result<()> {
    if got.len() != want.len() {
        bail!("{label} length {} != {}", got.len(), want.len());
    }
    for (a, b) in got.iter().zip(want) {
        if (a - b).abs() > 1e-6 {
            bail!("{label} mismatch: {a} vs {b}");
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    // etopo60.cdf ships in the crate under data/, so this test always runs (no env needed).
    // OHC_ETOPO overrides the path for an out-of-tree copy.
    fn etopo_path() -> PathBuf {
        std::env::var_os("OHC_ETOPO")
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/data/etopo60.cdf")))
    }

    #[test]
    fn etopo_matches_grid_and_reference() {
        let p = etopo_path();
        let grid = GridDef::mapping();
        let e = read_etopo(p, &grid).unwrap();
        assert_eq!(e.dim(), (180, 360)); // [lat, lon]
        // Reference cell from the Python parse: etopo[lat=100, lon=200] ≈ -4827.465 m.
        assert!((e[[100, 200]] - (-4827.46533203125)).abs() < 1e-3);
        let mn = e.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = e.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!((mn - (-7473.0)).abs() < 1.0);
        assert!((mx - 5731.0).abs() < 1.0);
    }
}

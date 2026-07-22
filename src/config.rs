//! Configuration.
//!
//! [`RunConfig`] holds the static constants + paths for a product/environment (one per
//! `config.toml`). The per-run **slice** — which single layer, over which months — is
//! supplied separately as a [`Slice`] and is mandatory on the command line. One run
//! processes exactly one layer; batching across layers is the job scheduler's job.

use anyhow::{bail, Result};
use serde::Deserialize;
use std::path::PathBuf;

/// One mapped pressure layer, bounds in dbar (shallow `top`, deep `bottom`).
#[derive(Debug, Clone, Copy, Deserialize)]
pub struct LayerSpec {
    pub top: i32,
    pub bottom: i32,
}

impl LayerSpec {
    pub fn tag(&self) -> String {
        format!("{}_{}", self.top, self.bottom)
    }
}

/// Static constants + paths for a run (everything except the per-run slice; see `Slice`).
#[derive(Debug, Clone, Deserialize)]
pub struct RunConfig {
    /// run identifier, set per-run via `--tag` (this default is a placeholder)
    #[serde(default = "default_tag")]
    pub run_tag: String,
    /// e.g. "potentialTemperature"
    pub var_name: String,
    /// e.g. "SpaceTimeTrend"
    pub model_name: String,
    /// keep cells with lat in [lo, hi]
    pub latitude_range_to_keep: [f64; 2],
    /// basin ids to drop (land + marginal seas)
    pub basins_to_remove: Vec<i16>,
    /// uniform bathymetry floor (m): flag cells shallower than this for every layer
    /// (`bed_above_floor`). `None` = no floor. WMO/GCOS product uses 300.
    #[serde(default)]
    pub bathy_floor_m: Option<f64>,
    /// dir holding the FullField mean `.mat` files
    pub dir_mean: PathBuf,
    /// dir holding the LocalCondSim ensemble `.mat` files
    pub dir_ensemble: PathBuf,
    /// where the zarr stores are written
    pub dir_out: PathBuf,
    /// path to etopo60.cdf
    pub etopo_path: PathBuf,
    /// path to basinmask_04.msk
    pub basinmask_path: PathBuf,
    #[serde(default = "default_cp0")]
    pub cp0: f64,
    #[serde(default = "default_rho0")]
    pub rho0: f64,
}

fn default_cp0() -> f64 { crate::consts::CP0 }
fn default_rho0() -> f64 { crate::consts::RHO0 }
fn default_tag() -> String { "UNSET".into() }

impl RunConfig {
    /// Constant defaults (current LocalGP conventions). The run tag and paths are
    /// placeholders — set them via `--tag` / config.toml / env vars.
    pub fn defaults() -> Self {
        RunConfig {
            run_tag: default_tag(), // required via --tag
            var_name: "potentialTemperature".into(),
            model_name: "SpaceTimeTrend".into(),
            latitude_range_to_keep: [-64.5, 64.5],
            basins_to_remove: vec![0, 5, 6, 7, 8, 9, 53],
            bathy_floor_m: None,
            dir_mean: PathBuf::from("."),
            dir_ensemble: PathBuf::from("."),
            dir_out: PathBuf::from("."),
            etopo_path: PathBuf::from("etopo60.cdf"),
            basinmask_path: PathBuf::from("basinmask_04.msk"),
            cp0: crate::consts::CP0,
            rho0: crate::consts::RHO0,
        }
    }

    fn fname_prefix(&self, cond_sim: bool) -> String {
        let cs = if cond_sim { "LocalCondSim" } else { "" };
        format!("{}FullField{}{}", self.var_name, cs, self.model_name)
    }

    /// Full path to one monthly `.mat` for a layer.
    /// `cond_sim=false` → FullField mean; `true` → LocalCondSim ensemble.
    pub fn mat_path(&self, layer: &LayerSpec, year: i32, month: u32, cond_sim: bool) -> PathBuf {
        let fname = format!(
            "{}_{}_{:02}_{}.mat",
            self.fname_prefix(cond_sim),
            layer.tag(),
            month,
            year
        );
        let dir = if cond_sim { &self.dir_ensemble } else { &self.dir_mean };
        dir.join(fname)
    }

    /// zarr store directory for a layer.
    pub fn store_path(&self, layer: &LayerSpec) -> PathBuf {
        self.dir_out
            .join(format!("ohc_{}_plev{}.zarr", self.run_tag, layer.tag()))
    }
}

/// One unit of work: a single layer over an (inclusive) year range and a month subset.
#[derive(Debug, Clone)]
pub struct Slice {
    pub layer: LayerSpec,
    pub years: [i32; 2],
    pub months: Vec<u32>,
}

impl Slice {
    /// Monthly (year, month) axis, day-15 implied.
    pub fn time_axis(&self) -> Vec<(i32, u32)> {
        let mut v = Vec::new();
        for y in self.years[0]..=self.years[1] {
            for &m in &self.months {
                v.push((y, m));
            }
        }
        v
    }

    /// CF time values: days since day-15 of the first month of the slice.
    pub fn time_days_since_start(&self) -> (Vec<f64>, (i32, u32)) {
        let axis = self.time_axis();
        let (y0, m0) = axis[0];
        let base = days_from_civil(y0, m0, 15);
        let days = axis
            .iter()
            .map(|&(y, m)| (days_from_civil(y, m, 15) - base) as f64)
            .collect();
        (days, (y0, m0))
    }
}

/// Days since 1970-01-01 (proleptic Gregorian). Hinnant's algorithm.
pub fn days_from_civil(y: i32, m: u32, d: u32) -> i64 {
    let y = if m <= 2 { y - 1 } else { y } as i64;
    let era = if y >= 0 { y } else { y - 399 } / 400;
    let yoe = (y - era * 400) as i64;
    let m = m as i64;
    let d = d as i64;
    let doy = (153 * (if m > 2 { m - 3 } else { m + 9 }) + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146097 + doe - 719468
}

/// Parse `--years` / `OHC_YEARS`: `"2016"` → `[2016,2016]`; `"2016:2018"`/`"2016-2018"` → `[2016,2018]`.
pub fn parse_years(s: &str) -> Result<[i32; 2]> {
    let parts: Vec<&str> = s.split(|c| c == ':' || c == '-').collect();
    match parts.len() {
        1 => {
            let y: i32 = parts[0].trim().parse()?;
            Ok([y, y])
        }
        2 => Ok([parts[0].trim().parse()?, parts[1].trim().parse()?]),
        _ => bail!("bad --years value: {s:?}"),
    }
}

/// Parse `--months` / `OHC_MONTHS`: `"8"` → `[8]`; `"1,2,3"`; `"1:3"`/`"1-3"` → `[1,2,3]`. Validates 1..=12.
pub fn parse_months(s: &str) -> Result<Vec<u32>> {
    let months: Vec<u32> = if s.contains(',') {
        s.split(',')
            .map(|p| p.trim().parse::<u32>())
            .collect::<std::result::Result<_, _>>()?
    } else {
        let parts: Vec<&str> = s.split(|c| c == ':' || c == '-').collect();
        match parts.len() {
            1 => vec![parts[0].trim().parse()?],
            2 => {
                let a: u32 = parts[0].trim().parse()?;
                let b: u32 = parts[1].trim().parse()?;
                (a..=b).collect()
            }
            _ => bail!("bad --months value: {s:?}"),
        }
    };
    if months.iter().any(|&m| !(1..=12).contains(&m)) {
        bail!("months out of range 1..=12: {months:?}");
    }
    Ok(months)
}

/// Parse `--layer` / `OHC_LAYER`: exactly one `top-bottom` (inner separator `-`, `_`, or `:`),
/// e.g. `"15-20"`, `"300_700"`, `"700:1850"`. Rejects multiple layers.
pub fn parse_layer(s: &str) -> Result<LayerSpec> {
    if s.contains(',') {
        bail!("--layer takes exactly one layer; got a list: {s:?}");
    }
    let parts: Vec<&str> = s.split(|c| c == '-' || c == '_' || c == ':').collect();
    if parts.len() != 2 {
        bail!("bad --layer value {s:?}, expected top-bottom (one layer)");
    }
    Ok(LayerSpec {
        top: parts[0].trim().parse()?,
        bottom: parts[1].trim().parse()?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_scope_overrides() {
        assert_eq!(parse_years("2016").unwrap(), [2016, 2016]);
        assert_eq!(parse_years("2004:2025").unwrap(), [2004, 2025]);
        assert_eq!(parse_months("8").unwrap(), vec![8]);
        assert_eq!(parse_months("1-3").unwrap(), vec![1, 2, 3]);
        assert!(parse_months("13").is_err());
        let l = parse_layer("15-20").unwrap();
        assert_eq!((l.top, l.bottom), (15, 20));
        assert_eq!(parse_layer("300_700").unwrap().bottom, 700);
        assert!(parse_layer("15-20,300-700").is_err()); // no lists
        assert!(parse_layer("15").is_err());
    }

    #[test]
    fn filenames_match_localgp_convention() {
        let c = RunConfig::defaults();
        let l = LayerSpec { top: 15, bottom: 20 };
        assert!(c
            .mat_path(&l, 2016, 8, false)
            .to_string_lossy()
            .ends_with("potentialTemperatureFullFieldSpaceTimeTrend_15_20_08_2016.mat"));
        assert!(c.mat_path(&l, 2016, 8, true).to_string_lossy().ends_with(
            "potentialTemperatureFullFieldLocalCondSimSpaceTimeTrend_15_20_08_2016.mat"
        ));
    }

    #[test]
    fn slice_time_axis_and_days() {
        let s = Slice { layer: LayerSpec { top: 15, bottom: 20 }, years: [2016, 2016], months: vec![8] };
        assert_eq!(s.time_axis(), vec![(2016, 8)]);
        let (days, base) = s.time_days_since_start();
        assert_eq!(days, vec![0.0]);
        assert_eq!(base, (2016, 8));

        let s2 = Slice { layer: LayerSpec { top: 15, bottom: 20 }, years: [2004, 2004], months: vec![1, 2] };
        let (days, _) = s2.time_days_since_start();
        assert_eq!(days, vec![0.0, 31.0]); // Jan15 -> Feb15 2004
    }
}

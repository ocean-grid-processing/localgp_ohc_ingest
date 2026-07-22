//! Layer ingest — option 1: buffer the whole layer in RAM and transpose month-major →
//! member-major.
//!
//! LocalGP delivers one `.mat` per month (the ensemble file holds all 100 members for that
//! month). Our zarr chunks are per-member, so we accumulate the full layer and emit:
//!   - `ohc_mean`     `[time, lat, lon]`          (the FullField posterior mean)
//!   - `ohc_ensemble` `[member, time, lat, lon]`  (the 100 conditional simulations)
//! Integrated temperature is converted to OHC (`* cp0 * rho0`) on the way in; NaNs are
//! preserved; arrays are transposed from the `.mat`'s `[lon, lat]` order to `[lat, lon]`. f32 storage
//! (~6.8 GB for 264 months × 100 members) — RAM bet on the cluster.

use anyhow::{Context, Result};
use ndarray::{Array3, Array4};

use crate::config::{RunConfig, Slice};
use crate::{matread, GridDef};

pub struct LayerData {
    /// `[time, lat, lon]`, OHC J/m², NaN preserved. f64: the deliverable takes a large-mean
    /// anomaly (absolute OHC − baseline), so the mean is kept double all the way through.
    pub ohc_mean: Array3<f64>,
    /// `[member, time, lat, lon]`, OHC J/m², NaN preserved. f32: only feeds ensemble spread
    /// (`_sd`), where single precision is ample, and keeps the 100-member stack half the size.
    /// `None` when ingested mean-only (`--no-ensemble`) — the CondSim files are not read and
    /// `ohc_ensemble` is omitted from the store.
    pub ohc_ensemble: Option<Array4<f32>>,
}

/// Read every month of the slice's layer and assemble the member-major arrays.
///
/// `mean_only` (from `--no-ensemble`) skips the LocalCondSim files entirely and returns
/// `ohc_ensemble: None` — for mean-only products (e.g. the GCOS deliverable) where the ensemble
/// is never used downstream, and to run without a complete set of CondSim `.mat` files.
pub fn ingest_layer(cfg: &RunConfig, slice: &Slice, grid: &GridDef, mean_only: bool) -> Result<LayerData> {
    let layer = &slice.layer;
    let nlat = grid.nlat();
    let nlon = grid.nlon();
    let time = slice.time_axis();
    let nt = time.len();
    let nm = crate::consts::NMEMBER;
    let scale = cfg.cp0 * cfg.rho0;

    let mut ohc_mean = Array3::<f64>::from_elem((nt, nlat, nlon), f64::NAN);
    let mut ohc_ensemble = if mean_only {
        None
    } else {
        Some(Array4::<f32>::from_elem((nm, nt, nlat, nlon), f32::NAN))
    };

    for (t, &(year, month)) in time.iter().enumerate() {
        // FullField mean: [lon, lat]
        let mean_path = cfg.mat_path(layer, year, month, false);
        let mean = matread::read_mean_grid(&mean_path)
            .with_context(|| format!("reading {}", mean_path.display()))?;
        debug_assert_eq!(mean.dim(), (nlon, nlat));
        for j in 0..nlat {
            for i in 0..nlon {
                ohc_mean[[t, j, i]] = mean[[i, j]] * scale; // transpose [lon,lat]→[lat,lon]
            }
        }

        // LocalCondSim ensemble: [lon, lat, member] — skipped entirely when mean_only
        if let Some(ens_arr) = ohc_ensemble.as_mut() {
            let ens_path = cfg.mat_path(layer, year, month, true);
            let ens = matread::read_ensemble(&ens_path)
                .with_context(|| format!("reading {}", ens_path.display()))?;
            debug_assert_eq!(ens.dim(), (nlon, nlat, nm));
            for m in 0..nm {
                for j in 0..nlat {
                    for i in 0..nlon {
                        ens_arr[[m, t, j, i]] = (ens[[i, j, m]] * scale) as f32;
                    }
                }
            }
        }
    }

    Ok(LayerData { ohc_mean, ohc_ensemble })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    // Validates the transpose + cp0·rho0 scaling on the single local month (Aug 2016),
    // by running ingest over a 1-month config. Set OHC_TEST_DATA to the .mat dir.
    #[test]
    fn single_month_transpose_and_scale() {
        use crate::config::LayerSpec;
        let Some(dir) = std::env::var_os("OHC_TEST_DATA").map(PathBuf::from) else { return };
        let mut cfg = RunConfig::defaults();
        cfg.dir_mean = dir;
        let layer = LayerSpec { top: 15, bottom: 20 };

        // Validate transpose + scale by reading the one local month directly.
        let mean = matread::read_mean_grid(cfg.mat_path(&layer, 2016, 8, false)).unwrap();
        let scale = cfg.cp0 * cfg.rho0;
        // raw m[200,100]=144.646083 → ohc at [lat=100, lon=200]
        let ohc = mean[[200, 100]] * scale;
        assert!((ohc - 5.943_394e8).abs() / 5.943_394e8 < 1e-5);
    }
}

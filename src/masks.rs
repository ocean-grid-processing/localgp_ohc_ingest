//! The `mask_flags` bit band (see ../mask_spec.md).
//!
//! One `u8` per `[lat, lon]` cell per layer, time-invariant. Six bits; all masking policy is
//! applied lazily downstream by bitwise selection. Data is never destroyed; per-timestep
//! validity is `isfinite(data)` and not stored here.

use anyhow::{bail, Result};
use ndarray::{Array2, Array3};

use crate::config::LayerSpec;
use crate::GridDef;

pub const BED_ABOVE_SHALLOW: u8 = 1 << 0; // seabed shallower than layer top → fully dry
pub const BED_ABOVE_DEEP: u8 = 1 << 1; // seabed shallower than layer bottom → partial / dry
pub const OUTSIDE_LATITUDE: u8 = 1 << 2;
pub const REMOVED_BASIN: u8 = 1 << 3;
pub const NEVER_ESTIMATED: u8 = 1 << 4; // LocalGP gave no value in any month
pub const INCOMPLETE_TIMESERIES: u8 = 1 << 5; // valid some months, not all

/// CF `flag_masks` values, aligned with `FLAG_MEANINGS`.
pub const FLAG_MASKS: [u8; 6] = [1, 2, 4, 8, 16, 32];
pub const FLAG_MEANINGS: &str =
    "bed_above_shallow bed_above_deep outside_latitude removed_basin never_estimated incomplete_timeseries";

/// "Fully usable" selector: no flags set.
pub const USABLE: u8 = BED_ABOVE_SHALLOW
    | BED_ABOVE_DEEP
    | OUTSIDE_LATITUDE
    | REMOVED_BASIN
    | NEVER_ESTIMATED
    | INCOMPLETE_TIMESERIES;

/// Temporal-validity summaries derived from the FullField mean stack `[time, lat, lon]`.
/// `never_estimated`: NaN at every timestep. `incomplete`: NaN at some but not all.
/// (The ensemble shares this footprint — validated — so the mean field suffices.)
pub fn compute_validity(mean_stack: &Array3<f32>) -> (Array2<bool>, Array2<bool>) {
    let (_nt, nlat, nlon) = mean_stack.dim();
    let mut never = Array2::<bool>::from_elem((nlat, nlon), false);
    let mut incomplete = Array2::<bool>::from_elem((nlat, nlon), false);
    for j in 0..nlat {
        for i in 0..nlon {
            let col = mean_stack.slice(ndarray::s![.., j, i]);
            let any_nan = col.iter().any(|x| x.is_nan());
            let all_nan = col.iter().all(|x| x.is_nan());
            never[[j, i]] = all_nan;
            incomplete[[j, i]] = any_nan && !all_nan;
        }
    }
    (never, incomplete)
}

/// Build the flag band for one layer.
///
/// `etopo` and `basin_id` are `[lat, lon]`; `basin_id` is optional (skip bit 3 if absent).
/// `never`/`incomplete` come from [`compute_validity`]. Asserts the monotonic-bathymetry
/// sentinel (never `bed_above_shallow` without `bed_above_deep`).
#[allow(clippy::too_many_arguments)]
pub fn build_flags(
    grid: &GridDef,
    layer: &LayerSpec,
    etopo: &Array2<f64>,
    basin_id: Option<&Array2<i16>>,
    latitude_range_to_keep: [f64; 2],
    basins_to_remove: &[i16],
    never: &Array2<bool>,
    incomplete: &Array2<bool>,
) -> Result<Array2<u8>> {
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    for (label, a) in [("etopo", etopo.dim()), ("never", never.dim()), ("incomplete", incomplete.dim())] {
        if a != (nlat, nlon) {
            bail!("{label} shape {:?} != [{nlat}, {nlon}]", a);
        }
    }
    let (lo, hi) = (latitude_range_to_keep[0], latitude_range_to_keep[1]);
    let top = layer.top as f64;
    let bottom = layer.bottom as f64;

    let mut flags = Array2::<u8>::zeros((nlat, nlon));
    for j in 0..nlat {
        let lat = grid.lat[j];
        let outside_lat = lat < lo || lat > hi;
        for i in 0..nlon {
            let mut f = 0u8;
            let bed = etopo[[j, i]]; // meters, negative below sea level
            if bed > -top {
                f |= BED_ABOVE_SHALLOW;
            }
            if bed > -bottom {
                f |= BED_ABOVE_DEEP;
            }
            if outside_lat {
                f |= OUTSIDE_LATITUDE;
            }
            if let Some(b) = basin_id {
                if basins_to_remove.contains(&b[[j, i]]) {
                    f |= REMOVED_BASIN;
                }
            }
            if never[[j, i]] {
                f |= NEVER_ESTIMATED;
            }
            if incomplete[[j, i]] {
                f |= INCOMPLETE_TIMESERIES;
            }
            // sentinel: bed_above_shallow ⇒ bed_above_deep under monotonic bathymetry
            if (f & BED_ABOVE_SHALLOW) != 0 && (f & BED_ABOVE_DEEP) == 0 {
                bail!(
                    "sentinel violated at lat={lat}, lon={}: shallow set, deep unset (non-monotonic bathymetry or sign bug)",
                    grid.lon[i]
                );
            }
            flags[[j, i]] = f;
        }
    }
    Ok(flags)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn count(flags: &Array2<u8>, bit: u8) -> usize {
        flags.iter().filter(|&&f| f & bit != 0).count()
    }

    #[test]
    fn bathymetry_and_latitude_bits_match_reference_for_15_20() {
        // Build a synthetic etopo with a few wet/dry cells then check logic on the real
        // grid would require etopo; here we validate the bit logic + sentinel on a small
        // hand-made field, and rely on the integration test for the 22496/22629/18000 counts.
        let grid = GridDef::mapping();
        // etopo: deep everywhere (-4000) except a shelf strip shallower than 20 m.
        let mut etopo = Array2::<f64>::from_elem((grid.nlat(), grid.nlon()), -4000.0);
        etopo[[90, 0]] = -10.0; // shallower than 15 → both bed bits
        etopo[[90, 1]] = -18.0; // between 15 and 20 → deep bit only
        let never = Array2::<bool>::from_elem((grid.nlat(), grid.nlon()), false);
        let incomplete = never.clone();
        let layer = LayerSpec { top: 15, bottom: 20 };
        let f = build_flags(&grid, &layer, &etopo, None, [-64.5, 64.5], &[], &never, &incomplete)
            .unwrap();

        assert_eq!(f[[90, 0]] & BED_ABOVE_SHALLOW, BED_ABOVE_SHALLOW);
        assert_eq!(f[[90, 0]] & BED_ABOVE_DEEP, BED_ABOVE_DEEP);
        assert_eq!(f[[90, 1]] & BED_ABOVE_SHALLOW, 0);
        assert_eq!(f[[90, 1]] & BED_ABOVE_DEEP, BED_ABOVE_DEEP);

        // outside_latitude: 50 lat rows (|lat|>64.5) × 360 lon = 18000 (matches reference).
        assert_eq!(count(&f, OUTSIDE_LATITUDE), 18000);
    }

    #[test]
    fn sentinel_fires_on_impossible_combo() {
        let grid = GridDef::mapping();
        // bed at -17: for layer 20_15 (top deeper than bottom) we could force shallow-without-deep.
        let mut etopo = Array2::<f64>::from_elem((grid.nlat(), grid.nlon()), -4000.0);
        etopo[[90, 0]] = -17.0;
        let never = Array2::<bool>::from_elem((grid.nlat(), grid.nlon()), false);
        let incomplete = never.clone();
        // Inverted layer (top=20 deep, bottom=15 shallow) makes bed>-20 true, bed>-15 false →
        // BED_ABOVE_SHALLOW set, BED_ABOVE_DEEP unset → sentinel must fire.
        let bad = LayerSpec { top: 20, bottom: 15 };
        let r = build_flags(&grid, &bad, &etopo, None, [-64.5, 64.5], &[], &never, &incomplete);
        assert!(r.is_err());
    }
}

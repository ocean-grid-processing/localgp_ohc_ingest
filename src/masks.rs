//! The `mask_flags` bit band (see ../mask_spec.md).
//!
//! One `u8` per `[lat, lon]` cell per layer, time-invariant. Eight bits; all masking policy is
//! applied lazily downstream by bitwise selection. Data is never destroyed; per-timestep
//! validity is `isfinite(data)` and not stored here.

use anyhow::{bail, Result};
use ndarray::{Array2, Array3, Array4};

use crate::config::LayerSpec;
use crate::GridDef;

pub const BED_ABOVE_SHALLOW: u8 = 1 << 0; // seabed shallower than layer top → fully dry
pub const BED_ABOVE_DEEP: u8 = 1 << 1; // seabed shallower than layer bottom → partial / dry
pub const OUTSIDE_LATITUDE: u8 = 1 << 2;
pub const REMOVED_BASIN: u8 = 1 << 3;
pub const NEVER_ESTIMATED: u8 = 1 << 4; // LocalGP gave no value in any month
pub const INCOMPLETE_TIMESERIES: u8 = 1 << 5; // valid some months, not all
pub const BED_ABOVE_FLOOR: u8 = 1 << 6; // seabed shallower than a fixed floor depth, uniform across layers (policy)
pub const ENSEMBLE_INCOMPLETE: u8 = 1 << 7; // some CondSim member is NaN-in-time here (mean may be fine)

/// CF `flag_masks` values, aligned with `FLAG_MEANINGS`.
pub const FLAG_MASKS: [u8; 8] = [1, 2, 4, 8, 16, 32, 64, 128];
pub const FLAG_MEANINGS: &str =
    "bed_above_shallow bed_above_deep outside_latitude removed_basin never_estimated incomplete_timeseries bed_above_floor ensemble_incomplete";

/// "Fully usable" selector: no flags set.
pub const USABLE: u8 = BED_ABOVE_SHALLOW
    | BED_ABOVE_DEEP
    | OUTSIDE_LATITUDE
    | REMOVED_BASIN
    | NEVER_ESTIMATED
    | INCOMPLETE_TIMESERIES
    | BED_ABOVE_FLOOR
    | ENSEMBLE_INCOMPLETE;

/// Temporal-validity summaries derived from the FullField mean stack `[time, lat, lon]`.
/// `never_estimated`: NaN at every timestep. `incomplete`: NaN at some but not all.
/// (The ensemble shares this footprint — validated — so the mean field suffices.)
pub fn compute_validity(mean_stack: &Array3<f64>) -> (Array2<bool>, Array2<bool>) {
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

/// Union of member NaN footprints for the `ensemble_incomplete` bit `[lat, lon]`.
///
/// True where **any** CondSim member is NaN at **any** timestep — i.e. the member contribution to
/// the original's mean∪members mask (`msk |= isnan(sum(member,3))` over all members). The mean's
/// own NaNs are already carried by `never`/`incomplete`; this adds the cells that only some
/// members drop (the flaky deep layers). `None` in `build_flags` when ingested `--no-ensemble`.
pub fn compute_ensemble_incomplete(ensemble: &Array4<f64>) -> Array2<bool> {
    let (nm, _nt, nlat, nlon) = ensemble.dim();
    let mut inc = Array2::<bool>::from_elem((nlat, nlon), false);
    for m in 0..nm {
        for (idx, &v) in ensemble.index_axis(ndarray::Axis(0), m).indexed_iter() {
            if v.is_nan() {
                inc[[idx.1, idx.2]] = true; // idx = (t, lat, lon)
            }
        }
    }
    inc
}

/// Build the flag band for one layer.
///
/// `etopo` and `basin_id` are `[lat, lon]`; `basin_id` is optional (skip bit 3 if absent).
/// `never`/`incomplete` come from [`compute_validity`]. Asserts the monotonic-bathymetry
/// sentinel (never `bed_above_shallow` without `bed_above_deep`).
///
/// `bathy_floor_m` (optional) sets `bed_above_floor` where the seabed is shallower than a fixed
/// floor depth, applied uniformly to every layer regardless of its own bounds. `None` = no floor
/// (the per-layer bed bits alone). Used by the WMO/GCOS product (floor = 300 m).
///
/// `ens_incomplete` (optional, from [`compute_ensemble_incomplete`]) sets `ensemble_incomplete`
/// where some CondSim member is NaN-in-time — the member half of the original's mean∪members
/// mask. `None` when ingested `--no-ensemble` (mean-only: the bit stays unset).
#[allow(clippy::too_many_arguments)]
pub fn build_flags(
    grid: &GridDef,
    layer: &LayerSpec,
    etopo: &Array2<f64>,
    basin_id: Option<&Array2<i16>>,
    latitude_range_to_keep: [f64; 2],
    basins_to_remove: &[i16],
    bathy_floor_m: Option<f64>,
    never: &Array2<bool>,
    incomplete: &Array2<bool>,
    ens_incomplete: Option<&Array2<bool>>,
) -> Result<Array2<u8>> {
    let (nlat, nlon) = (grid.nlat(), grid.nlon());
    for (label, a) in [("etopo", etopo.dim()), ("never", never.dim()), ("incomplete", incomplete.dim())] {
        if a != (nlat, nlon) {
            bail!("{label} shape {:?} != [{nlat}, {nlon}]", a);
        }
    }
    if let Some(ei) = ens_incomplete {
        if ei.dim() != (nlat, nlon) {
            bail!("ens_incomplete shape {:?} != [{nlat}, {nlon}]", ei.dim());
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
            if let Some(floor) = bathy_floor_m {
                if bed > -floor {
                    f |= BED_ABOVE_FLOOR;
                }
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
            if let Some(ei) = ens_incomplete {
                if ei[[j, i]] {
                    f |= ENSEMBLE_INCOMPLETE;
                }
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
        let f = build_flags(&grid, &layer, &etopo, None, [-64.5, 64.5], &[], None, &never, &incomplete, None)
            .unwrap();

        assert_eq!(f[[90, 0]] & BED_ABOVE_SHALLOW, BED_ABOVE_SHALLOW);
        assert_eq!(f[[90, 0]] & BED_ABOVE_DEEP, BED_ABOVE_DEEP);
        assert_eq!(f[[90, 1]] & BED_ABOVE_SHALLOW, 0);
        assert_eq!(f[[90, 1]] & BED_ABOVE_DEEP, BED_ABOVE_DEEP);

        // outside_latitude: 50 lat rows (|lat|>64.5) × 360 lon = 18000 (matches reference).
        assert_eq!(count(&f, OUTSIDE_LATITUDE), 18000);
    }

    #[test]
    fn bathy_floor_masks_shelves_uniformly() {
        // A 15_20 layer whose own bed bits keep everything (deep ocean), plus two shelf cells at
        // 150 m and 400 m. A 300 m floor must flag the 150 m cell but not the 400 m one, and it
        // must fire independently of the (unset) per-layer bed bits.
        let grid = GridDef::mapping();
        let mut etopo = Array2::<f64>::from_elem((grid.nlat(), grid.nlon()), -4000.0);
        etopo[[90, 0]] = -150.0; // shallower than 300 → floor bit
        etopo[[90, 1]] = -400.0; // deeper than 300 → no floor bit
        let never = Array2::<bool>::from_elem((grid.nlat(), grid.nlon()), false);
        let incomplete = never.clone();
        let layer = LayerSpec { top: 15, bottom: 20 };

        // No floor: neither cell flagged (both deeper than 20 m).
        let f0 = build_flags(&grid, &layer, &etopo, None, [-64.5, 64.5], &[], None, &never, &incomplete, None)
            .unwrap();
        assert_eq!(count(&f0, BED_ABOVE_FLOOR), 0);

        // Floor = 300 m: only the 150 m cell, and it does not touch the per-layer bed bits.
        let f = build_flags(&grid, &layer, &etopo, None, [-64.5, 64.5], &[], Some(300.0), &never, &incomplete, None)
            .unwrap();
        assert_eq!(f[[90, 0]] & BED_ABOVE_FLOOR, BED_ABOVE_FLOOR);
        assert_eq!(f[[90, 0]] & BED_ABOVE_DEEP, 0); // 150 m is deeper than the 20 m layer bottom
        assert_eq!(f[[90, 1]] & BED_ABOVE_FLOOR, 0);
        assert_eq!(count(&f, BED_ABOVE_FLOOR), 1);
    }

    #[test]
    fn ensemble_incomplete_unions_member_nans() {
        // A member NaN at (t, lat, lon) marks that cell; a fully-finite column stays clear.
        let grid = GridDef::mapping();
        let (nlat, nlon) = (grid.nlat(), grid.nlon());
        let mut ens = Array4::<f64>::from_elem((3, 4, nlat, nlon), 1.0); // 3 members, 4 timesteps
        ens[[2, 1, 90, 5]] = f64::NAN;                                   // member 2, t=1, one cell
        let ei = compute_ensemble_incomplete(&ens);
        assert!(ei[[90, 5]]);
        assert!(!ei[[90, 6]]);

        // build_flags sets the bit there and nowhere finite; None leaves it unset everywhere.
        let etopo = Array2::<f64>::from_elem((nlat, nlon), -4000.0);
        let never = Array2::<bool>::from_elem((nlat, nlon), false);
        let incomplete = never.clone();
        let layer = LayerSpec { top: 15, bottom: 20 };
        let with = build_flags(&grid, &layer, &etopo, None, [-64.5, 64.5], &[], None,
                               &never, &incomplete, Some(&ei)).unwrap();
        assert_eq!(with[[90, 5]] & ENSEMBLE_INCOMPLETE, ENSEMBLE_INCOMPLETE);
        assert_eq!(count(&with, ENSEMBLE_INCOMPLETE), 1);
        let without = build_flags(&grid, &layer, &etopo, None, [-64.5, 64.5], &[], None,
                                  &never, &incomplete, None).unwrap();
        assert_eq!(count(&without, ENSEMBLE_INCOMPLETE), 0);
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
        let r = build_flags(&grid, &bad, &etopo, None, [-64.5, 64.5], &[], None, &never, &incomplete, None);
        assert!(r.is_err());
    }
}

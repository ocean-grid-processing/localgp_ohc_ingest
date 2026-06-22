//! Grid geometry — analytic spherical cell areas.
//!
//! Spherical cell area on a sphere of radius `EARTH_RADIUS_M`. For a cell centred at latitude
//! φ spanning ±0.5° in lat and ±0.5° in lon, area = R² · Δλ · (sin φ_top − sin φ_bot). Area
//! depends only on latitude.

use ndarray::Array2;

use crate::consts::EARTH_RADIUS_M;
use crate::GridDef;

/// Per-cell area in m², indexed `[lat, lon]`.
pub fn cell_area(grid: &GridDef) -> Array2<f64> {
    let nlat = grid.nlat();
    let nlon = grid.nlon();
    // assume uniform 1° spacing (half-cell = 0.5°)
    let dlon_rad = 1.0_f64.to_radians();
    let half = 0.5_f64;
    let r2 = EARTH_RADIUS_M * EARTH_RADIUS_M;

    let mut a = Array2::<f64>::zeros((nlat, nlon));
    for (j, &lat) in grid.lat.iter().enumerate() {
        let s_hi = (lat + half).to_radians().sin();
        let s_lo = (lat - half).to_radians().sin();
        let area = r2 * dlon_rad * (s_hi - s_lo);
        for i in 0..nlon {
            a[[j, i]] = area;
        }
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_area_is_earth_surface() {
        let g = GridDef::mapping();
        let a = cell_area(&g);
        let total: f64 = a.sum();
        // Earth surface ≈ 5.1006e14 m² (validated against the Python reference: 5.100645e14).
        assert!((total - 5.100_645e14).abs() / 5.100_645e14 < 1e-4);
    }

    #[test]
    fn cells_shrink_toward_poles() {
        let g = GridDef::mapping();
        let a = cell_area(&g);
        let equator = a[[90, 0]]; // lat ≈ 0.5
        let polar = a[[0, 0]]; // lat ≈ -89.5
        assert!(polar < equator);
    }
}

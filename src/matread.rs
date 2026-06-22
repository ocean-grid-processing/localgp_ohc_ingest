//! Minimal MAT v7 reader for LocalGP `fullFieldGrid` output.
//!
//! LocalGP writes MATLAB v7 files: a 128-byte text header, then data elements where the
//! payload is a single zlib-compressed (`miCOMPRESSED`) `miMATRIX` holding one variable,
//! `fullFieldGrid` — `[lon,lat]` f64 for the FullField mean, `[lon,lat,member]` f64 for the
//! LocalCondSim ensemble. We only need that one double array, so this is a purpose-built
//! parser (no general MAT support), pure Rust: `flate2` to inflate + a small v5 element walk.
//!
//! Data is stored MATLAB column-major (Fortran order); we preserve that when building the
//! ndarray so logical indexing is `[lon, lat(, member)]`. NaNs are preserved verbatim.

use anyhow::{anyhow, bail, Result};
use byteorder::{LittleEndian, ReadBytesExt};
use ndarray::{Array2, Array3, ShapeBuilder};
use std::io::Read;
use std::path::Path;

// MAT data element types we care about.
const MI_INT32: u32 = 5;
const MI_SINGLE: u32 = 7;
const MI_DOUBLE: u32 = 9;
const MI_MATRIX: u32 = 14;
const MI_COMPRESSED: u32 = 15;

/// A parsed array: dims as stored (MATLAB order) + column-major data.
pub struct MatArray {
    pub name: String,
    pub dims: Vec<usize>,
    /// length == product(dims), MATLAB column-major (Fortran) order.
    pub data: Vec<f64>,
}

/// A single MAT data element (tag + payload slice), with the offset of the next element.
struct Element<'a> {
    dtype: u32,
    payload: &'a [u8],
    next: usize,
}

/// Read one element tag+payload at `off`. Handles the v5 "small element" format
/// (non-zero upper 16 bits of the type word ⇒ length there, 4 inline data bytes).
fn read_element(buf: &[u8], off: usize) -> Result<Element<'_>> {
    if off + 8 > buf.len() {
        bail!("element header past end of buffer at {off}");
    }
    let type_raw = (&buf[off..off + 4]).read_u32::<LittleEndian>()?;
    let small_len = type_raw >> 16;
    if small_len != 0 {
        // small element: type in low 16 bits, <=4 data bytes inline
        let dtype = type_raw & 0xFFFF;
        let nbytes = small_len as usize;
        let start = off + 4;
        Ok(Element {
            dtype,
            payload: &buf[start..start + nbytes],
            next: off + 8, // small element is always 8 bytes total
        })
    } else {
        let nbytes = (&buf[off + 4..off + 8]).read_u32::<LittleEndian>()? as usize;
        let start = off + 8;
        if start + nbytes > buf.len() {
            bail!("element payload past end ({start}+{nbytes} > {})", buf.len());
        }
        let pad = (8 - (nbytes % 8)) % 8; // elements align to 8 bytes
        Ok(Element {
            dtype: type_raw,
            payload: &buf[start..start + nbytes],
            next: start + nbytes + pad,
        })
    }
}

/// Decode a numeric payload (miDOUBLE / miSINGLE) into Vec<f64>.
fn payload_to_f64(dtype: u32, payload: &[u8]) -> Result<Vec<f64>> {
    match dtype {
        MI_DOUBLE => {
            let n = payload.len() / 8;
            let mut v = Vec::with_capacity(n);
            let mut rdr = payload;
            for _ in 0..n {
                v.push(rdr.read_f64::<LittleEndian>()?);
            }
            Ok(v)
        }
        MI_SINGLE => {
            let n = payload.len() / 4;
            let mut v = Vec::with_capacity(n);
            let mut rdr = payload;
            for _ in 0..n {
                v.push(rdr.read_f32::<LittleEndian>()? as f64);
            }
            Ok(v)
        }
        other => bail!("unsupported real-data type {other} (expected double/single)"),
    }
}

/// Parse the body of a `miMATRIX` element into a `MatArray`.
fn parse_matrix(body: &[u8]) -> Result<MatArray> {
    // subelements: array flags, dimensions (miINT32), name (miINT8), real part (pr)
    let flags = read_element(body, 0)?; // array flags — class lives here; not needed for doubles
    let dims_el = read_element(body, flags.next)?;
    if dims_el.dtype != MI_INT32 {
        bail!("dimensions element not miINT32 (got {})", dims_el.dtype);
    }
    let mut dims = Vec::new();
    let mut d = dims_el.payload;
    while !d.is_empty() {
        dims.push(d.read_i32::<LittleEndian>()? as usize);
    }
    let name_el = read_element(body, dims_el.next)?;
    let name = String::from_utf8_lossy(name_el.payload).to_string();
    let data_el = read_element(body, name_el.next)?;
    let data = payload_to_f64(data_el.dtype, data_el.payload)?;

    let expected: usize = dims.iter().product();
    if data.len() != expected {
        bail!(
            "data length {} != product(dims)={} for var {name:?}",
            data.len(),
            expected
        );
    }
    Ok(MatArray { name, dims, data })
}

/// Read the (single) `fullFieldGrid` variable from a LocalGP `.mat`.
pub fn read_fullfield(path: impl AsRef<Path>) -> Result<MatArray> {
    let raw = std::fs::read(path.as_ref())?;
    if raw.len() < 128 || &raw[0..6] != b"MATLAB" {
        bail!("not a MAT v5/v7 file: {}", path.as_ref().display());
    }
    let body = &raw[128..]; // after text header
    let el = read_element(body, 0)?;
    match el.dtype {
        MI_COMPRESSED => {
            let mut inflated = Vec::new();
            flate2::read::ZlibDecoder::new(el.payload).read_to_end(&mut inflated)?;
            let m = read_element(&inflated, 0)?;
            if m.dtype != MI_MATRIX {
                bail!("compressed element does not contain a matrix (got {})", m.dtype);
            }
            parse_matrix(m.payload)
        }
        MI_MATRIX => parse_matrix(el.payload),
        other => Err(anyhow!("first element is type {other}, expected matrix/compressed")),
    }
}

/// FullField mean → `Array2` indexed `[lon, lat]` (column-major preserved).
pub fn read_mean_grid(path: impl AsRef<Path>) -> Result<Array2<f64>> {
    let m = read_fullfield(&path)?;
    if m.dims.len() != 2 {
        bail!("expected 2D mean grid, got dims {:?}", m.dims);
    }
    let (nlon, nlat) = (m.dims[0], m.dims[1]);
    Array2::from_shape_vec((nlon, nlat).f(), m.data)
        .map_err(|e| anyhow!("reshape mean grid: {e}"))
}

/// LocalCondSim ensemble → `Array3` indexed `[lon, lat, member]` (column-major preserved).
pub fn read_ensemble(path: impl AsRef<Path>) -> Result<Array3<f64>> {
    let m = read_fullfield(&path)?;
    if m.dims.len() != 3 {
        bail!("expected 3D ensemble, got dims {:?}", m.dims);
    }
    let (nlon, nlat, nm) = (m.dims[0], m.dims[1], m.dims[2]);
    Array3::from_shape_vec((nlon, nlat, nm).f(), m.data)
        .map_err(|e| anyhow!("reshape ensemble: {e}"))
}

// ---------------------------------------------------------------------------
// Tests against the local 15_20 Aug-2016 sample. Set OHC_TEST_DATA to the dir
// holding the two .mat files (e.g. postprocesser/data). Skipped if unset.
// Ground-truth values locked via a Python reference parse of the same files.
// ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn data_dir() -> Option<PathBuf> {
        std::env::var_os("OHC_TEST_DATA").map(PathBuf::from)
    }

    const MEAN: &str = "potentialTemperatureFullFieldSpaceTimeTrend_15_20_08_2016.mat";
    const ENS: &str = "potentialTemperatureFullFieldLocalCondSimSpaceTimeTrend_15_20_08_2016.mat";

    #[test]
    fn mean_grid_matches_reference() {
        let Some(dir) = data_dir() else { return };
        let g = read_mean_grid(dir.join(MEAN)).unwrap();
        assert_eq!(g.dim(), (360, 180));
        let finite = g.iter().filter(|x| x.is_finite()).count();
        assert_eq!(finite, 31287);
        assert_eq!(g.iter().filter(|x| x.is_nan()).count(), 33513);
        assert!((g[[200, 100]] - 144.646083).abs() < 1e-5);
        assert!(g[[0, 0]].is_nan()); // Antarctica land corner
    }

    #[test]
    fn ensemble_shape_and_shared_footprint() {
        let Some(dir) = data_dir() else { return };
        let e = read_ensemble(dir.join(ENS)).unwrap();
        assert_eq!(e.dim(), (360, 180, 100));
        // All 100 members share an identical NaN footprint (the key design assumption).
        let m0 = e.index_axis(ndarray::Axis(2), 0).mapv(|x| x.is_nan());
        for k in 1..100 {
            let mk = e.index_axis(ndarray::Axis(2), k).mapv(|x| x.is_nan());
            assert_eq!(m0, mk, "member {k} footprint differs from member 0");
        }
        assert_eq!(
            e.index_axis(ndarray::Axis(2), 0)
                .iter()
                .filter(|x| x.is_finite())
                .count(),
            31287
        );
    }
}

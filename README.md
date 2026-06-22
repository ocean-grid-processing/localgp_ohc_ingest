# ohc_ingest

Rust ingest that turns LocalGP `.mat` output into clean, raw ocean-heat-content (OHC) grids
plus a bit-band mask, serialized as one zarr store per layer (one chunk file per ensemble
member). It is the heavy data-processing front half of a pipeline whose plotting and
consumer-specific analysis are deferred to a later Python stage that consumes these stores.

This is a from-scratch reimplementation of the data-processing core of the legacy MATLAB
`WMO2024_*` post-processor (in `../../postprocesser/code_to_Katie/WMO_2026/`).

## What it does (scope)

For each mapped pressure layer:

1. Read the LocalGP `.mat` files (FullField mean + 100-member LocalCondSim ensemble), month
   by month.
2. Convert integrated temperature to OHC (`× cp0 × rho0`), preserving NaNs.
3. Transpose month-major input → member-major arrays (buffer the whole layer in RAM).
4. Derive ancillary grids (`etopo`, `basin_id`, `cell_area`) and the `mask_flags` bit band.
5. Write a per-layer zarr v3 store. **No masks are applied to the data** — masking is carried
   in `mask_flags` and applied lazily downstream.

Out of scope (later, in Python): combined layers, anomalies, area integrals, trends,
plotting, and the per-collaborator NetCDF exports.

## Design docs

- [`../mask_spec.md`](../mask_spec.md) — the `mask_flags` bit band
- [`../zarr_schema.md`](../zarr_schema.md) — store layout
- [`../implementation_plan.md`](../implementation_plan.md) — build plan

## Inputs (all confirmed against the real files)

| input | notes |
|---|---|
| LocalGP `.mat` | MATLAB v7 (zlib-compressed); var `fullFieldGrid`, `[lon,lat]` mean / `[lon,lat,100]` ensemble |
| `etopo60.cdf` | Ferret 1° bathymetry; vars `ETOPO60X/ETOPO60Y/ROSE`; grid identical to mapping grid |
| `basinmask_04.msk` | WOA 0.25° basin table; nearest-neighbour to the 1° grid, surface column |

## Output

One store per layer: `ohc_<run>_plev<top>_<bottom>.zarr`, containing `ohc_mean`
`(time,lat,lon)`, `ohc_ensemble` `(member,time,lat,lon)` chunked one file per member,
`mask_flags` `(lat,lon)` (with CF `flag_masks`/`flag_meanings`), and `etopo` / `basin_id` /
`cell_area` ancillaries. zarr v3, `bytes`+`gzip` codecs (pure Rust), xarray-readable.

## Build

```bash
cargo build --release          # native
docker build -t ohc_ingest .   # container (pure Rust, no system deps)
```

## Cluster deployment (static binary)

For clusters with neither Docker nor Rust, build a fully static `x86_64-unknown-linux-musl`
binary locally and copy the single file up (no libc/runtime deps — the crate is pure Rust):

```bash
docker build -f Dockerfile.static --target bin --output type=local,dest=dist .
file dist/ohc_ingest          # → statically linked
scp dist/ohc_ingest cluster:~/bin/
```

Then run it directly on the cluster (one layer per job — a scheduler job array fans out the
full cube), with paths via env or a `config.toml`:

```bash
OHC_DIR_MEAN=... OHC_DIR_ENSEMBLE=... OHC_DIR_OUT=... \
OHC_ETOPO=.../etopo60.cdf OHC_BASINMASK=.../basinmask_04.msk \
  ./ohc_ingest --tag OP20260110 --layer 15-20 --years 2004:2025 --months 1:12
```

## Run

**One run processes exactly one layer.** The per-run slice (`--layer`, `--years`, `--months`)
is required; static constants + paths come from a config file or path env vars.

```bash
# env paths, Aug 2016 of the 15–20 dbar layer, tagged OP20260110:
OHC_DIR_MEAN=... OHC_DIR_ENSEMBLE=... OHC_DIR_OUT=... \
OHC_ETOPO=.../etopo60.cdf OHC_BASINMASK=.../basinmask_04.msk \
  ./target/release/ohc_ingest --tag OP20260110 --layer 15-20 --years 2016 --months 8

# from a config file (see config.example.toml), full record of one layer:
./target/release/ohc_ingest config.toml --tag OP20260110 --layer 300-700 --years 2004:2025 --months 1:12
```

### Required per-run slice (CLI > env)

| flag | env | examples |
|---|---|---|
| `--tag` | `OHC_TAG` | `OP20260110` (run identifier; labels store + metadata) |
| `--layer` | `OHC_LAYER` | `15-20`, `300_700`, `700:1850` (exactly one) |
| `--years` | `OHC_YEARS` | `2016`, `2004:2025` |
| `--months` | `OHC_MONTHS` | `8`, `1:3`, `1,6,12` |

The output store is named `ohc_<tag>_plev<top>_<bottom>.zarr`.

To build the whole cube, run one invocation per layer (e.g. a scheduler job array) — there is
deliberately no multi-layer mode, since each layer is an independent store.

## Test

Data-backed unit tests are gated on env vars (skipped if unset), with ground-truth locked
from the `15_20` Aug-2016 sample:

```bash
OHC_TEST_DATA=/path/to/postprocesser/data \
OHC_ETOPO=/path/to/data/etopo60.cdf \
OHC_BASINMASK=/path/to/data/basinmask_04.msk \
  cargo test
```

## Publish an ME4OH submission

The zarr store is the source of truth; `publish.py` projects it to a compliant `.nc` submission
— collapsing the selected mask bits to NaN, converting J/m² → TJ/m² and the time axis to days
since 1900-01-01, and writing `DATA(LONGITUDE, LATITUDE, TIME)` under the ME4OH filename. This
keeps the Rust binary pure (zarr only); NetCDF emission lives here in Python.

```bash
python scripts/publish.py /path/ohc_<tag>_plev15_20.zarr \
    --experiment B --product LocalGP --out submissions/
# -> submissions/OHC_<Y0>_<Y1>_lev15_20_expB_LocalGP.nc
```

Mask presets: `me4oh` (default) applies only physical/validity bits (so we submit the honest,
maximal valid field and let the assessment define the common domain); `wmo` applies all bits
(our latitude/basin-cropped product). Use `--levels LOW,HIGH` to set the filename's layer bounds
in meters (e.g. `--levels 0,286.6`) when the store's bounds aren't the submission bounds, and
`--anomaly` to subtract the per-cell time mean. Requires `netCDF4` in addition to the verify deps.

## Verify a written store

```bash
python scripts/verify_store.py /path/to/ohc_<tag>_plev15_20.zarr DIR_MEAN DIR_ENSEMBLE
```
Compares the **entire** store against the LocalGP `.mat` files (mean and ensemble from their
separate directories) — all timesteps of `ohc_mean` and all 100 members at all timesteps of
`ohc_ensemble` — checking every grid point matches exactly (after `cp0*rho0` and the lon/lat
transpose). It loads the full ensemble into memory once (~6.8 GB for the 264-month record), so
run it where there's RAM (inside the job allocation). Uses `scipy.io.loadmat` as an independent reader (not our
Rust parser), so it's a genuine oracle, and opens the store via xarray/zarr — the same path the
downstream consumer will use. Requires `xarray`, `zarr` (>=3 for the v3 store), `numpy`, `scipy`;
a recent xarray is needed for zarr-v3 support:

```bash
# fresh env:
conda create -n ohc -c conda-forge python=3.12 "xarray>=2025.1" "zarr>=3" scipy numpy
# or into an existing env (numpy/scipy usually already present):
pip install "zarr>=3" "xarray>=2024.10"
```

For a pinned environment, build the cross-check image (`Dockerfile.crosscheck`):

```bash
docker image build -f Dockerfile.crosscheck -t ohc_verify .
docker container run --rm -v /host/out:/out:ro \
  -v /host/FullField:/in_mean:ro -v /host/FullFieldLocalCondSim:/in_ens:ro \
  ohc_verify /out/ohc_<tag>_plev15_20.zarr /in_mean /in_ens
```

## Status & known check-points

All modules drafted and cross-validated in Python against the local sample. The crate was
**not** compiled in the authoring environment (no Rust toolchain there), so first `cargo build`
on the cluster is the real smoke test. Two spots most likely to need a small tweak:

- `zarrwrite.rs`: confirm your `zarr`/`xarray` version reads the emitted v3 metadata
  (`dimension_names`, `gzip` codec).

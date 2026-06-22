# ohc_ingest

`ohc_ingest` turns our LocalGP ocean-heat-content (OHC) mapping output into a clean, analysis-
ready store, and then into an ME4OH-protocol submission. It produces our group's contribution to
the MapEval4OceanHeat (ME4OH) mapping-method intercomparison; the cross-group assessment that
ingests every group's submission is a separate tool (`me4oh_assess`).

The pipeline has a Rust core and a thin Python edge:

```
LocalGP .mat  ──ohc_ingest (Rust)──▶  per-layer zarr store  ──publish.py──▶  ME4OH .nc submission
 (mean +                              (raw OHC + bit-band mask,                (DATA[, DATA_SD],
  ensemble)                            nothing masked out)                      mask applied as NaN)
```

The Rust binary is deliberately pure-Rust (no C deps) so it builds to a single static binary
for clusters with no Docker/Rust; all NetCDF work lives in the Python step, where libnetcdf is
already available.

## How this reflects the ME4OH protocol

The submission format is defined by the ME4OH protocol
(`MapEval4OceanHeat_protocol_Giglio_etal2023.pdf`). The numbers baked into this code come from
there:

- **Grid:** 1×1°, `lon 20.5…379.5`, `lat −89.5…89.5` — the protocol mandates this so no
  regridding is needed for the intercomparison. (Conveniently it's also LocalGP's native grid.)
- **Units:** OHC density in TJ/m², with `cp0 = 3989.244 J/kg/K`, `rho0 = 1030 kg/m³`.
- **Time:** monthly, days since 1900-01-01.
- **Submission file:** one NetCDF per layer, `DATA(LONGITUDE, LATITUDE, TIME)`, named
  `OHC_<Y0>_<Y1>_lev<low>_<high>_exp<X>_<product>.nc`.
- **Mask:** the protocol's `ocean_mask` is "points not-NaN at all timesteps," so a submission
  communicates "don't use this cell" only via NaN. `publish.py` is where our richer internal
  mask collapses to that NaN convention.

The store the Rust core writes is *not* a protocol artifact — it's our internal representation
(raw data + a bit-band mask, masking nothing). `publish.py` projects it down to the protocol.

## What the Rust core does (scope)

For one mapped layer:

1. Read the LocalGP `.mat` files (FullField mean + 100-member LocalCondSim ensemble), month by
   month.
2. Convert integrated temperature to OHC (`× cp0 × rho0`), preserving NaNs.
3. Transpose month-major input → member-major arrays (buffer the whole layer in RAM).
4. Derive ancillary grids (`etopo`, `basin_id`, `cell_area`) and the `mask_flags` bit band.
5. Write a per-layer zarr v3 store. **No masks are applied to the data** — masking is carried in
   `mask_flags` and applied lazily downstream.

Out of scope here: combined depth layers, anomalies, area integrals, trends, plotting, and the
cross-group assessment — all of which live in the Python consumer / `me4oh_assess`.

## Inputs

| input | notes |
|---|---|
| LocalGP `.mat` | MATLAB v7 (zlib-compressed); variable `fullFieldGrid`, `[lon,lat]` for the mean and `[lon,lat,100]` for the ensemble. Mean and ensemble live in **separate directories**. |
| `etopo60.cdf` | 1° bathymetry (classic NetCDF); vars `ETOPO60X`/`ETOPO60Y`/`ROSE`; its grid is identical to the mapping grid (asserted, not regridded). |
| `basinmask_04.msk` | WOA 0.25° basin table; nearest-neighbour to the 1° grid, surface column. |

## Output (the zarr store)

`ohc_<tag>_plev<top>_<bottom>.zarr`, containing:

- `ohc_mean` `(time, lat, lon)` — the posterior-mean OHC, J/m², NaN preserved.
- `ohc_ensemble` `(member, time, lat, lon)` — the 100 conditional simulations, chunked one
  file per member.
- `mask_flags` `(lat, lon)` — the bit band, with CF `flag_masks`/`flag_meanings` (see
  [`../mask_spec.md`](../mask_spec.md)).
- `etopo`, `basin_id`, `cell_area` `(lat, lon)` — ancillaries.

zarr v3, `bytes`+`gzip` codecs (pure Rust), xarray-readable. Layout details in
[`../zarr_schema.md`](../zarr_schema.md).

## The mask bit band

`mask_flags` is a `uint8` per cell carrying one bit per *reason* a cell might be excluded, so
the store never destroys data — masking is a downstream choice (and, for a submission, collapses
to NaN in `publish.py`). The bits:

| bit | value | name | meaning |
|---|---|---|---|
| 0 | 1 | `bed_above_shallow` | seafloor shallower than the layer's shallow edge (layer entirely in rock) |
| 1 | 2 | `bed_above_deep` | seafloor shallower than the layer's deep edge (seabed cuts through the layer) |
| 2 | 4 | `outside_latitude` | cell outside the kept latitude band |
| 3 | 8 | `removed_basin` | cell in a dropped basin (marginal / enclosed seas) |
| 4 | 16 | `never_estimated` | LocalGP produced no value in any month |
| 5 | 32 | `incomplete_timeseries` | valid in some months but not all |

Bits 0/1/4/5 are physical/validity reasons; bits 2/3 are policy reasons — `publish.py`'s presets
use exactly that split. Full definitions, the selector conventions, and the monotonic-bathymetry
sentinel are in [`../mask_spec.md`](../mask_spec.md).

## Build

```bash
cargo build --release          # native
docker build -t ohc_ingest .   # container (pure Rust, no system deps)
```

### Static binary for the cluster

For clusters with neither Docker nor Rust, build a fully static `x86_64-unknown-linux-musl`
binary and copy the single file up (no libc/runtime deps):

```bash
docker build -f Dockerfile.static --target bin --output type=local,dest=dist .
file dist/ohc_ingest          # → statically linked
scp dist/ohc_ingest cluster:~/bin/
```

## Run

**One run processes exactly one layer.** The per-run slice is required; static constants +
paths come from a `config.toml` (positional arg) or path env vars. To build the whole cube, run
one invocation per layer (e.g. a scheduler job array) — there is deliberately no multi-layer
mode, since each layer is an independent store.

```bash
# from a config file (see config.example.toml), full record of one layer:
./ohc_ingest config.toml --tag OP20260110 --layer 0-286.6 --years 2004:2025 --months 1:12

# or with paths from env instead of a config file:
OHC_DIR_MEAN=... OHC_DIR_ENSEMBLE=... OHC_DIR_OUT=... \
OHC_ETOPO=.../etopo60.cdf OHC_BASINMASK=.../basinmask_04.msk \
  ./ohc_ingest --tag OP20260110 --layer 0-286.6 --years 2016 --months 8
```

Note: when run from a scheduler, pass `config.toml` explicitly and use **absolute paths**
(inside it too) — a job's working directory is not guaranteed. The binary prints the resolved
config + paths on startup so a missing config is obvious.

### Required per-run slice (CLI > env)

| flag | env | examples |
|---|---|---|
| `--tag` | `OHC_TAG` | `OP20260110` (run identifier; labels the store + metadata) |
| `--layer` | `OHC_LAYER` | `0-286.6`, `300_700`, `700:1850` (exactly one) |
| `--years` | `OHC_YEARS` | `2016`, `2004:2025` |
| `--months` | `OHC_MONTHS` | `8`, `1:3`, `1,6,12` |

## Test

Data-backed unit tests are gated on env vars (skipped if unset), with ground-truth locked from
a single sample month/layer:

```bash
OHC_TEST_DATA=/dir/with/sample/mat/files \
OHC_ETOPO=/path/etopo60.cdf \
OHC_BASINMASK=/path/basinmask_04.msk \
  cargo test
```

## Python edge (publish + verify)

These scripts share one environment; build it once:

```bash
conda create -n ohc -c conda-forge python=3.12 "xarray>=2024.10" "zarr>=3" scipy numpy netCDF4
# or, into an existing env:  pip install -r scripts/requirements-crosscheck.txt
```

### publish.py — make the ME4OH submission

Projects a store to a compliant `.nc`: collapses the selected mask bits to NaN, converts
J/m² → TJ/m² and the time axis to days-since-1900, and writes `DATA(LONGITUDE, LATITUDE, TIME)`
under the ME4OH filename. It also adds `DATA_SD` (ensemble 1σ — the protocol's "associated
uncertainties, where available"); `--no-uncertainty` skips it (and the full-ensemble read).

```bash
python scripts/publish.py /path/ohc_<tag>_plev0_286.6.zarr \
    --experiment B --product LocalGP --out submissions/
# -> submissions/OHC_<Y0>_<Y1>_lev0_286.6_expB_LocalGP.nc
```

Mask presets: `me4oh` (default) applies only physical/validity bits, so we submit the honest,
maximal valid field and let the assessment define the common domain; `wmo` applies all bits
(our latitude/basin-cropped product). `--levels LOW,HIGH` sets the filename's layer bounds in
meters when they differ from the store's; `--anomaly` subtracts the per-cell time mean.

### Verification (two independent round-trips against the `.mat`)

```bash
# the store vs the upstream .mat (every grid point, all members):
python scripts/verify_store.py   STORE.zarr        DIR_MEAN DIR_ENSEMBLE

# the published .nc vs the upstream .mat (end-to-end pipeline):
python scripts/verify_publish.py SUBMISSION.nc     DIR_MEAN DIR_ENSEMBLE
```

Both use `scipy.io.loadmat` as an independent reader (not our Rust parser), so they're genuine
oracles. Comparisons hold to a small float32-scale tolerance, not bit-for-bit: the TJ/m²
conversion is done in float32, so a last-ULP (~1e-7 relative) difference is expected and
harmless; a real bug (transpose flip, unit error, wrong month) would still be caught. Both load
sizeable arrays, so run them inside the job allocation.

For a pinned, reproducible verify environment there is also `Dockerfile.crosscheck`:

```bash
docker image build -f Dockerfile.crosscheck -t ohc_verify .
docker container run --rm -v /host/out:/out:ro \
  -v /host/FullField:/in_mean:ro -v /host/FullFieldLocalCondSim:/in_ens:ro \
  ohc_verify /out/ohc_<tag>_plev0_286.6.zarr /in_mean /in_ens
```

## Layout

```
ohc_ingest/
├── src/                 Rust core: .mat + etopo + basinmask readers, masks, ingest, zarr writer
├── scripts/             Python edge: publish.py, verify_store.py, verify_publish.py
├── config.example.toml  constants + paths template
├── Dockerfile           runtime container
├── Dockerfile.static    static musl binary for the cluster
└── Dockerfile.crosscheck  pinned env for the verify scripts
```

Design notes: [`../mask_spec.md`](../mask_spec.md) (the mask bit band) and
[`../zarr_schema.md`](../zarr_schema.md) (store layout).

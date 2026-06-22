#!/usr/bin/env python3
"""Full round-trip check: does every grid point of the zarr match the upstream .mat?

Compares the entire store against the LocalGP .mat files — all timesteps of ohc_mean, and all
100 members at all timesteps of ohc_ensemble — for exact equality (after cp0*rho0 and the
lon/lat transpose). Uses scipy.io.loadmat as an independent reader (not our Rust parser), and
opens the store via xarray/zarr (the same path the downstream consumer uses).

    python verify_store.py STORE.zarr DIR_MEAN DIR_ENSEMBLE

Access pattern: the store is chunked one file per ensemble member, while the .mat are one file
per month, so we load the full ensemble into memory once (~6.8 GB for the 264-month record) and
stream the .mat. Run it where there's enough RAM (e.g. inside the job allocation).

Requires: xarray, zarr (>=3 for the v3 store), numpy, scipy. Recommended env:
    conda create -n ohc -c conda-forge python=3.12 "xarray>=2025.1" "zarr>=3" scipy numpy
"""
import datetime
import os
import sys

import numpy as np
import xarray as xr
from scipy.io import loadmat


def compare(expected, got):
    expected = expected.astype("float32")
    got = got.astype("float32")
    same_nan = np.array_equal(np.isnan(expected), np.isnan(got))
    finite = ~np.isnan(expected)
    max_diff = float(np.abs(expected[finite] - got[finite]).max()) if finite.any() else 0.0
    return same_nan, max_diff


def main(store, dir_mean, dir_ensemble):
    ds = xr.open_zarr(store, consolidated=False, decode_times=False)

    cp0 = ds.attrs["cp0"]
    rho0 = ds.attrs["rho0"]
    top = ds.attrs["layer_top"]
    bottom = ds.attrs["layer_bottom"]
    var = ds.attrs["var_name"]
    model = ds.attrs["model_name"]
    scale = cp0 * rho0

    units = ds["time"].attrs["units"]            # "days since YYYY-MM-15"
    y0, m0, d0 = (int(x) for x in units.split("since")[1].strip().split("-"))
    base = datetime.date(y0, m0, d0)
    times = np.asarray(ds["time"].values)
    nt = len(times)
    print("checking %s plev%d_%d — %d timesteps, %d members"
          % (ds.attrs.get("mapped_fields_tag"), top, bottom, nt, ds.sizes["member"]))

    # Read each side once: full store into memory, then stream the month-major .mat.
    print("loading full store into memory (~%.1f GB)…"
          % (ds["ohc_ensemble"].size * 4 / 1e9))
    zmean_all = ds["ohc_mean"].values            # [time, lat, lon]
    zens_all = ds["ohc_ensemble"].values         # [member, time, lat, lon]

    worst_mean = 0.0
    worst_ens = 0.0
    for t in range(nt):
        dt = base + datetime.timedelta(days=int(round(float(times[t]))))
        year, month = dt.year, dt.month
        stem = "%sFullField%%s%s_%d_%d_%02d_%d.mat" % (var, model, top, bottom, month, year)
        mean_path = os.path.join(dir_mean, stem % "")
        ens_path = os.path.join(dir_ensemble, stem % "LocalCondSim")

        # mean: loadmat gives [lon, lat] -> [lat, lon]
        mat_mean = loadmat(mean_path)["fullFieldGrid"].T * scale
        ok, md = compare(mat_mean, zmean_all[t])
        assert ok, "ohc_mean NaN footprint differs at %04d-%02d" % (year, month)
        assert md == 0.0, "ohc_mean differs at %04d-%02d (max %g)" % (year, month, md)
        worst_mean = max(worst_mean, md)

        # ensemble: [lon, lat, member] -> [member, lat, lon]
        mat_ens = np.transpose(loadmat(ens_path)["fullFieldGrid"], (2, 1, 0)) * scale
        ok, md = compare(mat_ens, zens_all[:, t])
        assert ok, "ohc_ensemble NaN footprint differs at %04d-%02d" % (year, month)
        assert md == 0.0, "ohc_ensemble differs at %04d-%02d (max %g)" % (year, month, md)
        worst_ens = max(worst_ens, md)

        if (t + 1) % 24 == 0 or t == nt - 1:
            print("  checked %d/%d timesteps (through %04d-%02d)" % (t + 1, nt, year, month))

    print("PASS — %d timesteps × %d members; ohc_mean max diff=%g, ohc_ensemble max diff=%g"
          % (nt, ds.sizes["member"], worst_mean, worst_ens))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("error: expected 3 arguments (got %d)\n\nusage: verify_store.py STORE.zarr DIR_MEAN DIR_ENSEMBLE"
                 % (len(sys.argv) - 1))
    main(sys.argv[1], sys.argv[2], sys.argv[3])

#!/usr/bin/env python3
"""Pipeline round-trip check: does the published .nc match the upstream .mat?

The companion to verify_store.py, but one stage later: it confirms the published submission
still equals the LocalGP .mat after the whole ingest+publish chain (mask applied,
J/m^2 -> TJ/m^2, time re-referenced to 1900). Uses scipy.io.loadmat as an independent reader.

    python verify_publish.py SUBMISSION.nc DIR_MEAN DIR_ENSEMBLE [--no-sd] [--ensemble]

`--ensemble` also checks the sibling `OHCENS_<...>.nc` (the full per-member ensemble written by
`publish.py --ensemble`, located by swapping the `OHC_` filename prefix) member-by-member
against the `.mat` ensemble.

The mean (DATA) is float64 end to end by default, so it matches a float64 recompute to ~1e-12;
published with `--dtype float32` it holds only to ~1e-6 (the float32 /1e12 rounding, and 1e12
isn't exactly representable in float32). The ensemble (DATA_SD, OHCENS) stays float32, so those
checks use the looser float32 tolerance. The DATA tolerance adapts to the stored dtype.
Requires: xarray, numpy, scipy, netCDF4.
"""
import argparse
import os

import numpy as np
import xarray as xr
from scipy.io import loadmat

TERA = 1e12
SD_RTOL = 1e-3
ENS_RTOL = 1e-6    # OHCENS members are float32; absorbs the float32 /1e12 rounding


def expected_data(mat_lonlat, cp0, rho0, out_dtype):
    """Recompute the published DATA in float64 (ingest stores the mean f64) and cast to the
    submission's stored dtype — float64 by default, float32 if published with --dtype float32."""
    v = (mat_lonlat.astype(np.float64) * cp0 * rho0) / TERA   # J/m^2 -> TJ/m^2, f64
    return v.astype(out_dtype)


def expected_sd(ens_lonlatmember, cp0, rho0):
    s = np.float32(ens_lonlatmember * cp0 * rho0)      # [lon, lat, member] f32
    sd = np.std(s, axis=2, ddof=1)                     # ensemble 1-sigma, f32
    return np.float32(sd.astype(np.float64) / TERA)


def expected_ens(ens_lonlatmember, cp0, rho0):
    """Per-member published value: same cast chain as expected_data, kept 3-D."""
    s = np.float32(ens_lonlatmember * cp0 * rho0)      # [lon, lat, member] f32
    return np.float32(s.astype(np.float64) / TERA)     # [lon, lat, member] TJ/m^2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("submission")
    ap.add_argument("dir_mean")
    ap.add_argument("dir_ensemble")
    ap.add_argument("--no-sd", action="store_true", help="skip the DATA_SD check")
    ap.add_argument("--ensemble", action="store_true",
                    help="also check the sibling OHCENS_<...>.nc ensemble file member-by-member")
    args = ap.parse_args()

    ds = xr.open_dataset(args.submission, decode_times=False)
    if ds["DATA"].ndim != 3:
        raise SystemExit(
            "error: %s has %d-D DATA, expected the 3-D OHC_ submission. Pass the OHC_ submission "
            "file, not the OHCENS_ ensemble file; --ensemble finds the OHCENS_ sibling itself."
            % (os.path.basename(args.submission), ds["DATA"].ndim))
    a = ds.attrs
    cp0, rho0 = a["cp0"], a["rho0"]
    var, model, layer = a["var_name"], a["model_name"], a["mapped_layer"]
    has_sd = ("DATA_SD" in ds) and not args.no_sd

    data = ds["DATA"].values                            # [lon, lat, time], TJ/m^2
    data_dtype = ds["DATA"].dtype
    data_rtol = 1e-12 if np.dtype(data_dtype) == np.float64 else 1e-6
    data_sd = ds["DATA_SD"].values if has_sd else None

    data_ens = None
    if args.ensemble:
        base = os.path.basename(args.submission)
        if not base.startswith("OHC_"):
            raise SystemExit("error: cannot derive the OHCENS_ name from %r (expected an OHC_ submission)" % base)
        ens_nc = os.path.join(os.path.dirname(args.submission), "OHCENS_" + base[len("OHC_"):])
        if not os.path.exists(ens_nc):
            raise SystemExit("error: --ensemble given but %s not found "
                             "(run publish.py --ensemble first)" % os.path.basename(ens_nc))
        eda = xr.open_dataset(ens_nc, decode_times=False)["DATA"]
        if eda.ndim != 4:
            raise SystemExit("error: %s DATA is %d-D, expected 4-D (MEMBER,LON,LAT,TIME)"
                             % (os.path.basename(ens_nc), eda.ndim))
        data_ens = eda.values   # [MEMBER, LON, LAT, TIME]

    days = np.round(ds["TIME"].values).astype("timedelta64[D]")
    dates = np.datetime64("1900-01-01") + days
    nt = len(dates)
    print("checking %s — %d timesteps, layer %s, uncertainty=%s, ensemble=%s"
          % (a.get("product"), nt, layer, has_sd, args.ensemble))

    worst_data = 0.0
    worst_sd = 0.0
    worst_ens = 0.0
    for t in range(nt):
        year = int(dates[t].astype("datetime64[Y]").astype(int) + 1970)
        month = int(dates[t].astype("datetime64[M]").astype(int) % 12 + 1)
        stem = "%sFullField%%s%s_%s_%02d_%d.mat" % (var, model, layer, month, year)

        mean_path = os.path.join(args.dir_mean, stem % "")
        exp = expected_data(loadmat(mean_path)["fullFieldGrid"], cp0, rho0, data_dtype)
        got = data[:, :, t]
        finite = np.isfinite(got)
        assert np.all(np.isfinite(exp[finite])), "DATA finite where .mat is NaN at %04d-%02d" % (year, month)
        if finite.any():
            md = float(np.abs(got[finite] - exp[finite]).max())
            scale = float(np.abs(exp[finite]).max())
            assert md <= data_rtol * scale, \
                "DATA differs at %04d-%02d (max %g, field scale %g)" % (year, month, md, scale)
            worst_data = max(worst_data, md)

        if has_sd or args.ensemble:
            ens_path = os.path.join(args.dir_ensemble, stem % "LocalCondSim")
            ens_mat = loadmat(ens_path)["fullFieldGrid"]    # [lon, lat, member]
            if has_sd:
                exp_s = expected_sd(ens_mat, cp0, rho0)
                got_s = data_sd[:, :, t]
                finite = np.isfinite(got_s)
                rel = np.abs(got_s[finite] - exp_s[finite]) / (np.abs(exp_s[finite]) + 1e-30)
                mr = float(rel.max()) if finite.any() else 0.0
                assert mr < SD_RTOL, "DATA_SD differs at %04d-%02d (max rel %g)" % (year, month, mr)
                worst_sd = max(worst_sd, mr)
            if args.ensemble:
                exp_e = np.transpose(expected_ens(ens_mat, cp0, rho0), (2, 0, 1))  # [member,lon,lat]
                got_e = data_ens[:, :, :, t]                                       # [MEMBER,LON,LAT]
                finite = np.isfinite(got_e)
                if finite.any():
                    me = float(np.abs(got_e[finite] - exp_e[finite]).max())
                    scale = float(np.abs(exp_e[finite]).max())
                    assert me <= ENS_RTOL * scale, \
                        "OHCENS differs at %04d-%02d (max %g, field scale %g)" % (year, month, me, scale)
                    worst_ens = max(worst_ens, me)

        if (t + 1) % 24 == 0 or t == nt - 1:
            print("  checked %d/%d (through %04d-%02d)" % (t + 1, nt, year, month))

    msg = "PASS — %d timesteps; DATA max diff=%g" % (nt, worst_data)
    if has_sd:
        msg += ", DATA_SD max rel diff=%g" % worst_sd
    if args.ensemble:
        msg += ", OHCENS max diff=%g" % worst_ens
    print(msg)


if __name__ == "__main__":
    main()

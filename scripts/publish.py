#!/usr/bin/env python3
"""Project an ohc_ingest zarr store to an ME4OH-compliant NetCDF submission.

The zarr store is our source of truth (raw OHC + a bit-band mask, nothing masked out). A
compliant submission can only say "don't use this point" via NaN, so this step collapses the
selected mask bits to NaN, converts J/m^2 -> TJ/m^2 and the time axis to days since 1900-01-01,
and writes DATA(LONGITUDE, LATITUDE, TIME) under the ME4OH filename.

    python publish.py STORE.zarr --experiment B --product LocalGP \
        [--preset me4oh|wmo] [--levels LOW,HIGH] [--anomaly] [--ensemble] [--out DIR]

Mask presets (see ../mask_spec.md):
  me4oh (default) = physical/validity bits only (never_estimated, incomplete_timeseries,
                    bed_above_shallow, bed_above_deep) — submit the honest, maximal valid
                    field and let the assessment define the common domain.
  wmo             = all bits (adds outside_latitude, removed_basin) — our cropped product.

--ensemble additionally writes the full conditional-simulation ensemble as a sibling file
  OHCENS_<...>.nc with DATA(MEMBER, LONGITUDE, LATITUDE, TIME) — same mask, units, and time
  axis — for downstream uses that derive per-member quantities before collapsing to a spread.
  It is NOT an ME4OH submission (different filename, extra dimension). Reads all members.

Requires: xarray, zarr>=3, numpy, netCDF4.
"""
import argparse
import datetime
import os

import numpy as np
import xarray as xr

# mask bit values (mask_spec.md)
BITS = {
    "bed_above_shallow": 1,
    "bed_above_deep": 2,
    "outside_latitude": 4,
    "removed_basin": 8,
    "never_estimated": 16,
    "incomplete_timeseries": 32,
}
PRESETS = {
    "me4oh": ["never_estimated", "incomplete_timeseries", "bed_above_shallow", "bed_above_deep"],
    "wmo": list(BITS),
}
TERA = 1e12


def preset_mask_value(preset):
    v = 0
    for name in PRESETS[preset]:
        v |= BITS[name]
    return v


def fmt_lev(x):
    xf = float(x)
    return str(int(xf)) if xf == int(xf) else ("%g" % xf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("store")
    ap.add_argument("--experiment", required=True)
    ap.add_argument("--product", default=None, help="product name for the filename (default: run tag)")
    ap.add_argument("--preset", default="me4oh", choices=list(PRESETS))
    ap.add_argument("--levels", default=None, help="LOW,HIGH meters for the filename (default: store layer bounds)")
    ap.add_argument("--anomaly", action="store_true", help="subtract the per-cell time mean before writing")
    ap.add_argument("--no-uncertainty", action="store_true",
                    help="skip the ensemble standard-deviation field DATA_SD (reads all members)")
    ap.add_argument("--ensemble", action="store_true",
                    help="also write the full ensemble as OHCENS_<...>.nc, DATA(MEMBER,LON,LAT,TIME)")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()

    ds = xr.open_zarr(args.store, consolidated=False)  # decodes time -> datetime64
    g = ds.attrs
    product = args.product or g["mapped_fields_tag"]

    # --- collapse the selected mask bits to NaN ---
    mval = preset_mask_value(args.preset)
    masked = xr.DataArray((ds["mask_flags"].values.astype("uint8") & mval) != 0,
                          dims=("lat", "lon"))
    data = ds["ohc_mean"].where(~masked) / TERA        # [time, lat, lon], TJ/m^2
    if args.anomaly:
        data = data - data.mean("time")

    # --- ensemble 1-sigma (the protocol's "associated uncertainties, where available") ---
    # ddof=1 (sample standard deviation); this reads all ensemble members.
    include_sd = not args.no_uncertainty
    sd = (ds["ohc_ensemble"].std("member", ddof=1).where(~masked) / TERA) if include_sd else None

    # --- time -> days since 1900-01-01 ---
    t = ds["time"].values                              # datetime64
    days1900 = (t - np.datetime64("1900-01-01T00:00:00")) / np.timedelta64(1, "D")
    years = t.astype("datetime64[Y]").astype(int) + 1970
    y0, y1 = int(years.min()), int(years.max())

    # --- layer bounds (meters) for the filename ---
    if args.levels:
        low, high = (s.strip() for s in args.levels.split(","))
    else:
        low, high = fmt_lev(g["layer_top"]), fmt_lev(g["layer_bottom"])

    # --- compliant dataset: DATA(LONGITUDE, LATITUDE, TIME) [+ optional DATA_SD] ---
    def to_lon_lat_time(da):
        return da.transpose("lon", "lat", "time").values.astype("float32")

    data_vars = {"DATA": (("LONGITUDE", "LATITUDE", "TIME"), to_lon_lat_time(data))}
    if include_sd:
        data_vars["DATA_SD"] = (("LONGITUDE", "LATITUDE", "TIME"), to_lon_lat_time(sd))
    out = xr.Dataset(
        data_vars,
        coords={
            "LONGITUDE": ("LONGITUDE", ds["lon"].values),
            "LATITUDE": ("LATITUDE", ds["lat"].values),
            "TIME": ("TIME", days1900.astype("float64")),
        },
    )
    out["LONGITUDE"].attrs = {"units": "degrees_east", "axis": "X"}
    out["LATITUDE"].attrs = {"units": "degrees_north", "axis": "Y"}
    out["TIME"].attrs = {"units": "days since 1900-01-01 00:00:00",
                         "calendar": "proleptic_gregorian", "axis": "T"}
    out["DATA"].attrs = {"units": "TJ/m^2", "long_name": "ocean heat content density"}
    if include_sd:
        out["DATA_SD"].attrs = {
            "units": "TJ/m^2",
            "long_name": "ocean heat content density, ensemble standard deviation (1-sigma)",
            "comment": "std across %d conditional-simulation members (ddof=1)" % ds.sizes["member"],
        }
    out.attrs = {
        "Conventions": "CF-1.8",
        "product": product,
        "experiment": args.experiment,
        "period": "%d_%d" % (y0, y1),
        "layer_m": "%s_%s" % (low, high),
        "source": g.get("source", ""),
        "var_name": g["var_name"],
        "model_name": g["model_name"],
        "mapped_layer": "%d_%d" % (int(g["layer_top"]), int(g["layer_bottom"])),
        "cp0": g["cp0"], "rho0": g["rho0"],
        "mask_preset": args.preset,
        "mask_applied": " ".join(PRESETS[args.preset]),
        "created": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    if include_sd:
        out.attrs["ensemble_size"] = int(ds.sizes["member"])

    fname = "OHC_%d_%d_lev%s_%s_exp%s_%s.nc" % (y0, y1, low, high, args.experiment, product)
    path = os.path.join(args.out, fname)
    chunk_enc = {"zlib": True, "complevel": 4, "_FillValue": np.float32(np.nan)}
    enc = {"DATA": dict(chunk_enc)}
    if include_sd:
        enc["DATA_SD"] = dict(chunk_enc)
    out.to_netcdf(path, engine="netcdf4", format="NETCDF4", encoding=enc)
    print("wrote", path, "(%d timesteps, preset=%s, uncertainty=%s)"
          % (len(days1900), args.preset, include_sd))

    # --- optional: the full ensemble as a member-dimensioned sibling file ---
    if args.ensemble:
        ens = (ds["ohc_ensemble"].where(~masked) / TERA).astype("float32")
        ens = ens.transpose("member", "lon", "lat", "time").rename(
            {"member": "MEMBER", "lon": "LONGITUDE", "lat": "LATITUDE", "time": "TIME"})
        ens = ens.assign_coords(MEMBER=ds["member"].values,
                                LONGITUDE=ds["lon"].values,
                                LATITUDE=ds["lat"].values,
                                TIME=days1900.astype("float64"))
        eds = ens.to_dataset(name="DATA")
        eds["LONGITUDE"].attrs = {"units": "degrees_east", "axis": "X"}
        eds["LATITUDE"].attrs = {"units": "degrees_north", "axis": "Y"}
        eds["TIME"].attrs = {"units": "days since 1900-01-01 00:00:00",
                             "calendar": "proleptic_gregorian", "axis": "T"}
        eds["MEMBER"].attrs = {"long_name": "conditional-simulation member"}
        eds["DATA"].attrs = {"units": "TJ/m^2",
                             "long_name": "ocean heat content density (per ensemble member)"}
        eds.attrs = dict(out.attrs)
        eds.attrs["ensemble_size"] = int(ds.sizes["member"])
        eds.attrs["note"] = ("full conditional-simulation ensemble for per-member downstream "
                             "analysis; NOT a single-field ME4OH submission")

        ename = "OHCENS_%d_%d_lev%s_%s_exp%s_%s.nc" % (y0, y1, low, high, args.experiment, product)
        epath = os.path.join(args.out, ename)
        nlon, nlat, ntime = len(ds["lon"]), len(ds["lat"]), len(days1900)
        eenc = {"DATA": {"zlib": True, "complevel": 4, "_FillValue": np.float32(np.nan),
                         "chunksizes": (1, nlon, nlat, ntime)}}
        eds.to_netcdf(epath, engine="netcdf4", format="NETCDF4", encoding=eenc)
        print("wrote", epath, "(ensemble: %d members, preset=%s)"
              % (ds.sizes["member"], args.preset))


if __name__ == "__main__":
    main()

# Mask bit band (`mask_flags`)

`ohc_ingest` never deletes data. Instead of NaN-ing cells out, it records — per grid cell — one
bit for each independent *reason* the cell might be excluded. Downstream code chooses which
reasons to honor with a bitwise test, and `publish.py` collapses a chosen set to NaN when it
writes a submission. This document defines the bits.

## The variable

| field | value |
|---|---|
| name | `mask_flags` |
| dtype | `uint8` |
| dims | `(lat, lon)` — time-invariant; one per layer store |
| fill | none; every cell is defined |

CF attributes on the variable:

```
flag_masks    = 1, 2, 4, 8, 16, 32, 64
flag_meanings = "bed_above_shallow bed_above_deep outside_latitude removed_basin never_estimated incomplete_timeseries bed_above_floor"
```

## Bits

| bit | value | name | set when | derived from |
|:--:|:--:|---|---|---|
| 0 | 1 | `bed_above_shallow` | seafloor shallower than the layer's shallow edge (layer entirely in rock) | `etopo > −shallow_edge` |
| 1 | 2 | `bed_above_deep` | seafloor shallower than the layer's deep edge (seabed cuts through the layer) | `etopo > −deep_edge` |
| 2 | 4 | `outside_latitude` | cell outside the kept latitude band | `lat ∉ [lo, hi]` |
| 3 | 8 | `removed_basin` | cell in a dropped basin (marginal / enclosed seas) | `basin_id ∈ remove-list` |
| 4 | 16 | `never_estimated` | LocalGP produced no value in any month | `all t: ¬isfinite` |
| 5 | 32 | `incomplete_timeseries` | valid in some months but not all | `(any t: ¬isfinite) ∧ ¬never_estimated` |
| 6 | 64 | `bed_above_floor` | seafloor shallower than a fixed floor depth, applied uniformly to every layer | `etopo > −floor` |
| 7 | 128 | (reserved) | | |

Bits 0/1/4/5 are **physical / validity** reasons; bits 2/3/6 are **policy** reasons. Bits 4 and 5
together encode the temporal-validity state: neither set = always valid, bit 5 = sometimes
valid, bit 4 = never valid. Per-timestep validity itself is not stored — it's `isfinite(data)`.

`bed_above_floor` differs from `bed_above_deep`: the deep bit is per-layer (against that layer's
own bottom edge), while the floor bit is a single depth applied to **all** layers alike. It is set
only when `bathy_floor_m` is configured (`None` = off). The WMO/GCOS product uses a 300 m floor so
every layer shares one open-ocean footprint (matches the original `max_pressure_to_keep = 300`).

## Selecting cells

"Fully usable" is "no bit set":

```
USABLE = 0b01111111   (= 127)
keep = (mask_flags & USABLE) == 0
```

Drop bits from the selector to relax a policy: a global (un-cropped) integral uses
`USABLE & ~outside_latitude`; tolerating partial-depth cells drops `bed_above_deep`; including
marginal seas drops `removed_basin`. `publish.py` ships two presets — `me4oh` honors only the
physical bits (0,1,4,5), `wmo` honors all seven (adds `outside_latitude`, `removed_basin`,
`bed_above_floor`).

## Sentinel (free correctness check)

Under monotonic bathymetry `bed_above_shallow ⟹ bed_above_deep`, so `bit 0 set with bit 1
unset` is impossible. If it ever occurs, suspect a non-monotonic cell, a grid misalignment, or
a sign error. The ingest asserts it never fires.

## Conventions

- **"above" means shallower** (smaller depth, higher in the water column): `bed_above_shallow`
  means the seafloor is shallower than the layer's shallow edge.
- **Bits may co-occur** — each is an independent reason (a dry shelf cell can be both
  `bed_above_shallow` and `never_estimated`; a Mediterranean cell both `removed_basin` and
  `outside_latitude`).
- **Depth bounds (dbar) are compared to `etopo` (m) as if equal** — ~1% error near the surface,
  a few % at depth, negligible at a 1° land/sea boundary.
- **`incomplete_timeseries` is relative to the stored time axis.** If a consumer subsets time,
  the authoritative per-month truth is `isfinite(data)`.

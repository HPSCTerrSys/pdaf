(gracetws)=
# GRACE/TWS Data Assimilation #

This page describes the assimilation of Terrestrial Water Storage (TWS)
anomaly observations derived from GRACE/GRACE-FO satellite data into the
eCLM/CLM5 land-surface model using the PDAF-OMI framework.

## Configuration ##

TWS DA is enabled by setting
[`CLM:update_tws = 1`](enkfpf:clm:update_tws) in `enkfpf.par`.  The
following additional parameters control its behaviour:

| Parameter | Section | Purpose |
|-----------|---------|---------|
| [`CLM:update_tws`](enkfpf:clm:update_tws) | `[CLM]` | Enable/disable TWS DA |
| [`DA:state_setup`](enkfpf:da:state_setup) | `[DA]` | State vector composition |
| [`DA:TWS_smoother`](enkfpf:da:tws_smoother) | `[DA]` | Instantaneous vs. monthly-mean water storage |
| [`DA:max_inc`](enkfpf:da:max_inc) | `[DA]` | Maximum increment as fraction of state value |
| [`DA:exclude_greenland`](enkfpf:da:exclude_greenland) | `[DA]` | Exclude Greenland from state vector |
| [`DA:set_zero_start`](enkfpf:da:set_zero_start) | `[DA]` | Reset running average before first GRACE obs |

## State Vector ##

Only hydrologically active gridcells (columns where
`col%hydrologically_active` is `.TRUE.`) contribute to the TWS state
vector.  Greenland can optionally be excluded via
[`DA:exclude_greenland`](enkfpf:da:exclude_greenland).

The composition of the state vector is controlled by
[`DA:state_setup`](enkfpf:da:state_setup):

| `state_setup` | Contents per gridcell |
|:---:|---|
| `0` | Liquid + ice summed per soil layer (all layers down to bedrock), snow water equivalent, surface water |
| `1` | Single pre-aggregated TWS value |
| `2` | Surface soil moisture (layer 1), root-zone moisture (layer 4), deep moisture (layer 13), snow water equivalent |

## Temporal Averaging ##

GRACE observations represent monthly averages.  The parameter
[`DA:TWS_smoother`](enkfpf:da:tws_smoother) selects which CLM fields
are used to populate the state vector:

- `TWS_smoother = 0`: Instantaneous CLM fields are used.  This is
  appropriate when the observation file already contains monthly-mean
  GRACE anomalies referenced to a model climatology.
- `TWS_smoother = 1`: Monthly running-mean CLM fields (computed
  internally by eCLM) are used.  The reset of the running average is
  controlled by a flag inside the observation file; the first reset can
  be advanced by [`DA:set_zero_start`](enkfpf:da:set_zero_start).

## Observation Operator ##

The observation operator (implemented in `obs_GRACE_pdafomi.F90`) maps
the state vector to an observed TWS value by summing all water-storage
components for each gridcell and subtracting the temporal mean read from
a reference file (`mean_filename`).  This converts absolute water
storage to a TWS anomaly comparable to GRACE observations.

The PDAF-OMI localization radius for GRACE observations is specified in
units of grid cells (not km, unlike the SM localization radius).

## Increment Capping ##

To prevent unphysically large updates,
[`DA:max_inc`](enkfpf:da:max_inc) caps the analysis increment for each
water-storage component:

```
if |Δx| > max_inc × x  →  Δx = sign(max_inc × x, Δx)
```

The default value of `1.0` limits the increment to 100 % of the
current state value (i.e. the updated value cannot become negative
through this mechanism).

## Safety Checks ##

- Columns that are not hydrologically active are never included in the
  state vector.
- The running temporal mean is only reset when the observation file
  requests it (or once at start-up via `DA:set_zero_start`).
- PDAF-OMI deallocates observation arrays after each analysis step to
  avoid memory leaks in multi-step runs.

## Configuration Examples ##

### Minimal GRACE DA setup ###

```ini
[CLM]
update_tws = 1

[DA]
state_setup    = 0
TWS_smoother   = 0
max_inc        = 1.0
set_zero_start = 0
```

### Monthly-mean TWS with Greenland excluded ###

```ini
[CLM]
update_tws = 1

[DA]
state_setup        = 0
TWS_smoother       = 1
max_inc            = 0.5
exclude_greenland  = 1
set_zero_start     = 720
```

`set_zero_start = 720` resets the running monthly mean 720 CLM time
steps before the first GRACE observation (e.g. one month at 30-minute
time steps).

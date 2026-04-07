(lstda)=
# Land Surface Temperature Data Assimilation #

Land Surface Temperature data assimilation (LST-DA) in TSMP-PDAF enables
the assimilation of remotely-sensed land surface temperature (LST)
observations into the eCLM land-surface model. The analysis step updates
one or more temperature variables in the CLM state, and auxiliary masking
and increment-clipping parameters allow the update to be restricted to
physically meaningful regimes.

## Configuration ##

LST-DA is controlled by the following parameters in the [`[CLM]` section of
`enkfpf.par`](enkfpf:clm):

- [`CLM:update_T`](enkfpf:clm:update_T): selects which CLM temperature
  variables are placed in the state vector and updated after the analysis
  step.
- [`CLM:statevec_max_layer`](enkfpf:clm:statevec_max_layer): limits the
  number of `t_soisno` layers included in the state vector for options 2
  and 3 (default 25, i.e. all layers).
- [`CLM:T_mask_snow`](enkfpf:clm:T_mask_snow): optionally masks out
  columns with snow cover from the temperature update.
- [`CLM:increment_type`](enkfpf:clm:increment_type): switches between
  multiplicative (default) and additive increments.
- [`CLM:T_max_increment`](enkfpf:clm:T_max_increment): clips additive
  increments to a maximum magnitude (only relevant when
  `CLM:increment_type=1`).
- [`CLM:T_mask_T`](enkfpf:clm:T_mask_T): masks out columns whose
  surface-layer soil temperature falls below a threshold close to
  freezing.

When LST-DA is active (`CLM:update_T != 0`), the temperature state
variables replace the soil water content (SWC) in the state vector.
Simultaneous assimilation of SWC and LST in the same PDAF update step is
not supported.

## State Vector ##

The state vector content depends on `CLM:update_T`. For options 2–5
the state vector is defined at the gridcell level (one value per grid
cell); for option 1 it is defined at the patch level.

| `update_T` | State vector variables                                                                        |
|:-----------|:----------------------------------------------------------------------------------------------|
| `1`        | `t_grnd` (per patch) + `t_veg` (per patch)                                                    |
| `2`        | `t_skin` + `t_soisno` (`min(nlevgrnd, statevec_max_layer)` layers) + `t_veg` (per gridcell)  |
| `3`        | Same as `2`, plus `t_grnd` (per gridcell)                                                     |
| `4`        | Same as `2`, plus `t_h2osfc` (per gridcell)                                                   |
| `5`        | Same as `3`, plus `t_h2osfc` (per gridcell)                                                   |

## Observation Operator ##

The observation operator maps the CLM state to a simulated LST:

- **Option 1**: The simulated LST is computed from `t_grnd` and
  `t_veg` using the radiometric mixing formula of [Kustas & Anderson
  (2009)](https://doi.org/10.1016/j.agrformet.2009.05.016) (Eq. 7),
  weighted by leaf area index (LAI).
- **Options 2–5**: The skin temperature `t_skin` is used directly as
  the simulated LST. Observations are indexed at the gridcell level
  (one observation per grid cell).

## Increment Application ##

After the PDAF analysis step, the updated state vector values are written
back to the CLM temperature arrays. Two increment modes are available,
selected via [`CLM:increment_type`](enkfpf:clm:increment_type):

**Multiplicative increment** (`CLM:increment_type=0`, default):
Each temperature variable is scaled by the ratio of the posterior to the
prior gridcell-mean skin temperature:

```
t_update = t_prior × (TSKIN_out / TSKIN_in)
```

**Additive increment** (`CLM:increment_type=1`):
The difference between posterior and prior gridcell-mean TSKIN is added
directly to each temperature variable:

```
t_update = t_prior + (TSKIN_out - TSKIN_in)
```

When `CLM:increment_type=1`, the increment is clipped to
`±CLM:T_max_increment` (default 5 K). Increments that exceed this
threshold are replaced by the clipped value and a warning counter is
incremented; the total number of clipped increments is printed after each
analysis step.

## Masking ##

Two independent masking conditions can suppress the update for individual
columns or grid cells:

**Snow masking** (`CLM:T_mask_snow`): When set to `1`, columns with a
snow depth exceeding `CLM:T_mask_snow_depth` (default 1 mm) are excluded
from the temperature update. This avoids applying a bare-soil LST
increment to snow-covered grid cells. The threshold can be adjusted via
[`CLM:T_mask_snow_depth`](enkfpf:clm:T_mask_snow_depth).

**Freeze masking** (`CLM:T_mask_T`): The update is suppressed whenever
the soil/snow temperature of the surface layer (`t_soisno(:,1)`) falls
below `T_freeze + CLM:T_mask_T`. With the default value of `0`, this
masks all grid cells at or below the freezing point (273.15 K).

Both masks are evaluated independently for each patch or column.

## Safety Checks ##

- NaN checks are applied to all updated temperature variables; a warning
  is printed to stdout if a NaN value is detected.
- Clipped additive increments are counted per analysis step and the total
  is printed for `t_skin`, `t_soisno`, `t_veg`, and `t_grnd` separately.
- The update is skipped entirely for a given patch/column if the
  gridcell-mean TSKIN change falls below `1e-7 K` (numerical tolerance).

## Configuration Examples ##

### Assimilate LST into ground and vegetation temperature (option 1) ###

```text
[CLM]
update_T       = 1
increment_type = 0
```

The simulated LST is computed from `t_grnd` and `t_veg` via the
[Kustas & Anderson (2009)](https://doi.org/10.1016/j.agrformet.2009.05.016)
radiometric mixing formula. Both variables are updated with a
multiplicative increment.

### Assimilate LST into skin and soil temperature (option 2) ###

```text
[CLM]
update_T       = 2
increment_type = 0
T_mask_snow    = 1
T_mask_T       = 0.0
```

TSKIN is used as the simulated LST. After analysis, `t_skin`, all
`t_soisno` layers, and `t_veg` are scaled by the TSKIN increment factor.
Snow-covered columns are excluded from the update.

### Assimilate LST with additive increment and clipping (option 2) ###

```text
[CLM]
update_T         = 2
increment_type   = 1
T_max_increment  = 3.0
T_mask_T         = 2.0
```

The TSKIN increment is applied additively and clipped to ±3 K. Grid cells
whose surface-layer temperature falls below 275.15 K (freezing + 2 K) are
excluded.

### Assimilate LST including ground temperature (option 3) ###

```text
[CLM]
update_T       = 3
increment_type = 0
T_mask_snow    = 1
```

Like option 2, but `t_grnd` is additionally included in the state vector
and updated. Snow-covered columns are masked out.

### Assimilate LST including surface water temperature (option 4) ###

```text
[CLM]
update_T       = 4
increment_type = 0
T_mask_snow    = 1
```

Like option 2, but `t_h2osfc` is additionally included in the state vector
and updated. Snow-covered columns are masked out.

### Assimilate LST including ground and surface water temperature (option 5) ###

```text
[CLM]
update_T       = 5
increment_type = 0
T_mask_snow    = 1
```

Like option 3, but `t_h2osfc` is additionally included in the state vector
and updated. Snow-covered columns are masked out.

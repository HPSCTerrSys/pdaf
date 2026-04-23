(snowda)=
# Snow Data Assimilation #

Snow data assimilation (Snow-DA) in TSMP-PDAF enables the assimilation
of snow-related observations—snow depth (SD) or snow water equivalent
(SWE)—into the eCLM land-surface model.

## Configuration ##

Snow-DA is controlled by two parameters in the [`[CLM]` section of
`enkfpf.par`](enkfpf:clm):

- [`CLM:update_snow`](enkfpf:clm:update_snow): selects which snow
  variable(s) are placed in the state vector and which CLM layer
  variables are updated after the analysis step.
- [`CLM:update_snow_repartitioning`](enkfpf:clm:update_snow_repartitioning):
  selects the method used to redistribute the updated bulk snow
  quantity across the individual CLM snow layers.

When Snow-DA is active (`CLM:update_snow != 0`), the snow state
variables replace the soil water content (SWC) in the state
vector. Simultaneous assimilation of SWC and snow in the same PDAF
update step is not supported.

## State Vector ##

The state vector size and content depend on `CLM:update_snow`. Snow
variables are always stored at the grid-cell level (one value per grid
cell); snow depth and SWE are taken from the first column of each grid
cell.

| `update_snow` | State vector size | Variables in state vector                  |
|:--------------|:------------------|:-------------------------------------------|
| `1`           | 1 × n_gridcells   | SD (`snow_depth_col`)                      |
| `2`           | 1 × n_gridcells   | SWE (`h2osno_col`)                         |
| `3`           | 2 × n_gridcells   | SD (`snow_depth_col`) + SWE (`h2osno_col`) |
| `4`           | 2 × n_gridcells   | SD (`snow_depth_col`) + SWE (`h2osno_col`) |
| `5`           | 2 × n_gridcells   | SD (`snow_depth_col`) + SWE (`h2osno_col`) |
| `6`           | 2 × n_gridcells   | SD (`snow_depth_col`) + SWE (`h2osno_col`) |
| `7`           | 2 × n_gridcells   | SD (`snow_depth_col`) + SWE (`h2osno_col`) |

For cases 3–7, SD occupies the first block of the state vector and SWE
the second block (offset by `clm_varsize`).

## Snow Layer Variable Update ##

After the PDAF analysis step updates the bulk snow quantities in the
state vector, the individual CLM snow layer variables must be adjusted
to remain physically consistent. The following layer-indexed arrays
(indexed over active snow layers `snlsno(j)+1 : 0`) may be updated:

| CLM variable                  | Description                                  |
|-------------------------------|----------------------------------------------|
| `h2osoi_ice_col`              | Ice content per snow layer (kg m⁻²)          |
| `h2osoi_liq_col`              | Liquid water content per snow layer (kg m⁻²) |
| `dz` (column layer thickness) | Thickness of each snow layer (m)             |

The redistribution method is selected via
`CLM:update_snow_repartitioning`.

### Repartitioning methods 1 and 2 ###

Methods 1 and 2 are only available for `CLM:update_snow=1` or `2`.
They operate by computing the posterior SWE and distributing the
SWE gain/loss across the snow layers, then deriving layer thickness
changes from the local snow density.

**Method 1** (DART-style): The entire SWE change is applied to the
bottom snow layer (`i=0`). All other layers are left unchanged. Layer
thickness is adjusted based on the local snow density. This approach
follows the DART/CLM implementation
(<https://github.com/NCAR/DART/blob/main/models/clm/dart_to_clm.f90>).

**Method 2**: The SWE change is distributed across all active snow
layers proportionally to their current water content. Each layer
receives a fraction `(h2osoi_liq + h2osoi_ice) / SWE_prior` of the
total SWE change.

For both methods, the posterior SWE (`h2osno_po`) is derived from:
- For `update_snow=1`: `h2osno_po = snow_depth_out × ρ_avg × frac_sno`,
  where `ρ_avg = min(800, h2osno / snowdp)` is the layer-averaged snow
  density.
- For `update_snow=2`: `h2osno_po = h2osno_out` directly from the state
  vector.

### Repartitioning method 3 (default) ###

Method 3 uses a multiplicative increment to scale the layer
variables. The scaling factor is the ratio of the posterior to the
prior bulk snow value. It is available for all values of
`CLM:update_snow`.

The following table summarises which layer variables are scaled and by
which increment for each `update_snow` option:

| `update_snow` | `h2osoi_ice` scaled by | `h2osoi_liq` scaled by | `dz` scaled by |
|:---:|:---:|:---:|:---:|
| `1` | `SD_out / SD_in` | — | — |
| `2` | `SWE_out / SWE_in` | — | — |
| `3` | `SWE_out / SWE_in` | — | — |
| `4` | `SWE_out / SWE_in` | `SWE_out / SWE_in` | `SD_out / SD_in` |
| `5` | `SWE_out / SWE_in` | — | `SD_out / SD_in` |
| `6` | `SWE_out / SWE_in` | `SWE_out / SWE_in` | — |
| `7` | `SWE_out / SWE_in` | — | — |

`SD_in` / `SD_out`: prior / posterior snow depth
`SWE_in` / `SWE_out`: prior / posterior snow water equivalent

For `update_snow=1`, the column-level SWE (`h2osno`) is additionally
updated after the layer loop to match the adjusted layer ice contents.

The increment update is only applied when both the prior and posterior
bulk snow values exceed `1e-6` to avoid numerical instabilities with
near-zero snow amounts.

## Observation Files ##

Snow-DA uses CLM-type observation files (see [Observation files](obs)).

- For `CLM:update_snow=1` and `3`–`7`, the observation variable should
  be snow depth (`SNOWDP` / `SNOW_DEPTH`).
- For `CLM:update_snow=2`, the observation variable should be SWE.

Observations are expected at the grid-cell level (one observation per
grid cell, no layer index). In the TSMP-PDAF observation operator
(`init_dim_obs_f_pdaf.F90`), when `CLM:update_snow != 0`, the
observation index `obs_index_p` is set to the grid-cell index rather
than a layer-specific state vector index.

## Safety Checks ##

The following protective measures are applied during the Snow-DA update:

- If the posterior SD (case 1) or SWE (case 2) is negative or zero,
  a warning is printed to stdout and no update is applied for that
  column.
- For cases 4–7, the increment update is skipped unless all four
  quantities (prior SD, posterior SD, prior SWE, posterior SWE)
  exceed `1e-6`.
- NaN checks are applied to all updated layer variables; a warning is
  printed if NaN values are detected.
- The SWC masking flag [`CLM:swc_mask_snow`](enkfpf:clm) is
  independent of Snow-DA and only affects SWC updates.

## Configuration Examples ##

### Assimilate snow depth (SD) only ###

```text
[CLM]
update_snow                = 1
update_snow_repartitioning = 3
```

The state vector contains one SD value per grid cell. After the PDAF
update, `snow_depth_col` is written back to CLM and `h2osoi_ice` in
each snow layer is scaled by `SD_out / SD_in`.

### Assimilate SWE only ###

```text
[CLM]
update_snow                = 2
update_snow_repartitioning = 3
```

The state vector contains one SWE value per grid cell. After the PDAF
update, `h2osno_col` is written back to CLM and `h2osoi_ice` in each
snow layer is scaled by `SWE_out / SWE_in`.

### Assimilate SD, update ice and layer thickness ###

```text
[CLM]
update_snow                = 5
update_snow_repartitioning = 3
```

Both SD and SWE are placed in the state vector. After the PDAF update,
`h2osoi_ice` is scaled by the SWE increment and layer thickness `dz`
is scaled by the SD increment.

### Assimilate SD, update ice, liquid water, and layer thickness ###

```text
[CLM]
update_snow                = 4
update_snow_repartitioning = 3
```

Both SD and SWE are placed in the state vector. After the PDAF update,
`h2osoi_ice` and `h2osoi_liq` are scaled by the SWE increment, and
layer thickness `dz` is scaled by the SD increment.

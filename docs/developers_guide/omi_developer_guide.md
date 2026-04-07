(omi:developer-guide)=
# Developer Guide: Adding a New Observation Type to the OMI Interface

This page describes how to add a new observation type — generically
called C below — to the TSMP-PDAF OMI (Observation Model Interface)
framework. It covers the architecture, the role of each source file,
and the changes required, using the existing GRACE-DA and SM-DA
implementations as reference.

## Architecture Overview

Observation handling sits between the PDAF library and the eCLM model
state in three layers:

```
PDAF library (src/)
    │  calls callback routines at each analysis step
    ▼
callback_obs_pdafomi.F90        ← single routing hub
    │  dispatches to one routine per observation type
    ▼
obs_GRACE_pdafomi.F90           ← one module per observation type
obs_SM_pdafomi.F90
obs_C_pdafomi.F90  (new)
    │  read observations, build H, apply H to state vector
    ▼
eCLM state vector (via enkf_clm_mod, mod_assimilation, …)
```

Each observation-type module is self-contained. The only coupling
between modules is through the central routing file and the shared
`mod_assimilation` module.

## Files Involved

To add a new observation type `C`, touch exactly **seven** things:

| File                                           | Action                                                                |
|------------------------------------------------|-----------------------------------------------------------------------|
| `interface/framework/obs_C_pdafomi.F90`        | **Create** (new file, modelled on the template)                       |
| `interface/framework/callback_obs_pdafomi.F90` | Add calls for C to every callback routine                             |
| `interface/framework/mod_assimilation.F90`     | Declare `cradius_C`, `sradius_C`, and any other C-specific parameters |
| `interface/framework/mod_read_obs.F90`         | Add `'C'` case to `update_obs_type`                                   |
| `interface/model/eclm/enkf_clm_mod_5.F90`      | Add `clmupdate_C` flag (parallel to `clmupdate_tws`, `clmupdate_swc`) |
| `interface/framework/init_pdaf.F90`            | Import `assim_C` and set it from `clmupdate_C`                        |
| `interface/framework/init_pdaf_parse.F90`      | Import `rms_obs_C` and add a `parse` call for it                      |

The template at `templates/omi/obs_OBSTYPE_pdafomi_TEMPLATE.F90` is a
good starting skeleton. The existing SM and GRACE modules show how the
template is applied in this codebase.

---

## Step 1 — Create `obs_C_pdafomi.F90`

### Module-level variables

At minimum, declare `assim_C` (a logical switch controlled by
`callback_obs_pdafomi`), `rms_obs_C` (a constant fallback error
standard deviation), and the PDAFomi types `thisobs` (of type `obs_f`,
thread-safe) and `thisobs_l` (of type `obs_l`, declared
`THREADPRIVATE`).

`assim_C` and `rms_obs_C` are wired up externally:
- **`init_pdaf.F90`** imports `assim_C` and sets it from `clmupdate_C`
  (i.e. `assim_C = (clmupdate_C /= 0)`), so the observation type is
  enabled automatically when the corresponding model-state update flag
  is active. Commented-out placeholders mark the insertion point.
- **`init_pdaf_parse.F90`** imports `rms_obs_C` and populates it via a
  `parse` call that reads `rms_obs_C` from `enkfpf.par`. Commented-out
  placeholders mark the insertion point here as well.

Additional module-level allocatable arrays — e.g. a cached
climatological mean needed for an anomaly operator — follow the same
pattern as in `obs_GRACE_pdafomi.F90`.

### The `obs_f` type

`obs_f` carries all observation metadata and is populated inside
`init_dim_obs_C`. Mandatory fields:

| Field                   | Type                   | Meaning                                               |
|-------------------------|------------------------|-------------------------------------------------------|
| `thisobs%doassim`       | `INTEGER`              | 1 = assimilate, 0 = skip                              |
| `thisobs%disttype`      | `INTEGER`              | distance metric for localization (see below)          |
| `thisobs%ncoord`        | `INTEGER`              | number of spatial coordinates (usually 2)             |
| `thisobs%id_obs_p(:,:)` | `INTEGER, ALLOCATABLE` | state-vector index for each process-local observation |

Frequently useful optional fields: 
- `icoeff_p`: bilinear interpolation weights
- `obs_err_type`: 0 = Gaussian, 1 = Laplace
- `inno_omit`: innovation outlier rejection threshold

**`disttype` choices:**
- `0` — Cartesian distance; coordinates in any consistent unit (GRACE
  uses integer grid-cell indices, so `cradius_GRACE` is in number of
  grid cells).
- `3` — geographic distance; coordinates in degrees, radii in km (used
  by SM).

Use `disttype=3` for point observations distributed geographically (C,
SM) and `disttype=0` for products already mapped to grid-cell indices
(e.g. GRACE).

---

### Routine 1 — `init_dim_obs_C`

Called once per analysis step before the filter loop. Responsibilities:

1. Set the mandatory `thisobs` fields.
2. Read the global observation vector from the NetCDF file using
   `read_obs_nc_type` from `mod_read_obs`.
3. Handle the `dim_obs == 0` case (no observations scheduled for this
   step): allocate dummy arrays of size 1 and call
   `PDAFomi_gather_obs` before returning.
4. For each global observation that falls within the PE-local domain,
   record the state-vector index in `thisobs%id_obs_p` and populate
   the local arrays `obs_p`, `ivar_obs_p`, and `ocoord_p`.
5. Call `PDAFomi_gather_obs` to finalize observation bookkeeping across
   MPI ranks.

**Reading observations:** Both SM and GRACE use `read_obs_nc_type`,
which checks the `type_clm` NetCDF variable and returns `dim_obs = 0`
if the file does not contain observations of the requested type.

**Who reads the file:** SM reads only on `mype_filter==0` and
broadcasts; GRACE reads on all filter processes because the spatial
averaging operator needs each PE to know its contribution. For point
observations like C, the SM pattern (read + broadcast) is simpler.

**Snapping observations to the state vector:** For point observations
(SM, C), loop over all global observations and all PE-local gridcells
and record the state-vector index for matching pairs. For GRACE the
logic is inverted: each observation aggregates multiple gridcells, so
the loop iterates over gridcells and maps each one to the observation
it contributes to.

---

### Routine 2 — `obs_op_C`

Applies the observation operator H, mapping the PE-local model state to
the observation space. For a simple point observation (SM pattern):
allocate a local output array, extract `state_p(obs_index_p(i))` for
each local observation, then call `PDAFomi_gather_obsstate` to assemble
the global result in `ostate`.

For aggregating operators (GRACE pattern), sum all gridcell
contributions mapping to the same observation, apply any transformation
(e.g. subtract temporal mean for TWS anomaly), and then call
`PDAFomi_gather_obsstate`. Follow the SM pattern for C unless your
observations represent spatial averages over a large footprint.

---

### Routines 3 & 4 — Localization (required for LESTKF/LETKF)

**`init_dim_obs_l_C`** counts observations within the localization
radius of the current local analysis domain. It is a thin wrapper
around `PDAFomi_init_dim_obs_l`, passing `thisobs_l`, `thisobs`, the
local domain coordinates `coords_l` (in the same units and `disttype`
as `ocoord_p`), and the localization parameters `cradius_C` and
`sradius_C`.

**`localize_covar_C`** is required only for the LEnKF. It is an
equally thin wrapper around `PDAFomi_localize_covar` with the same
localization parameters.

---

## Step 2 — Update `callback_obs_pdafomi.F90`

This file dispatches PDAF callbacks to each observation module. Add C
in the same pattern as GRACE and SM. The commented-out 
`!USE obs_C_pdafomi` placeholders show exactly where to insert new
entries.

For each of the eight callback subroutines, add a `USE obs_C_pdafomi`
statement and one call to the corresponding C routine:

- **`init_dim_obs_pdafomi`**: set `assim_C = .true.`, call
  `init_dim_obs_C`, and add `dim_obs_C` to the running total. The
  `assim_*` flags exist to disable an observation type for debugging;
  keep `assim_C = .true.` here because the module itself returns
  `dim_obs = 0` when no observations are present. (This logic is
  planned to be complemented in the future by setting `assim_C` from
  the EnKF input file)
- **`obs_op_pdafomi`**: call `obs_op_C`. Call order does not affect
  correctness because offsets in `ostate` are determined by the order
  of `init_dim_obs_*` calls above.
- **`init_dim_obs_l_pdafomi`**: call `init_dim_obs_l_C`. The variable
  `dim_obs_l` is incremented inside each call; PDAF initializes it to
  0 before invoking the callback.
- **`localize_covar_pdafomi`**: call `localize_covar_C`.
- **Remaining callbacks** (`add_obs_err`, `init_obscovar`, `prodRinvA`,
  `prodRinvA_l`, `deallocate_obs`): each needs a corresponding C
  call, all of which are thin wrappers around the matching PDAFomi
  generic routine, identical in structure to the existing GRACE and SM
  entries.

---

## Step 3 — Add Parameters to `mod_assimilation.F90`

Declare `cradius_C` and `sradius_C` next to the existing GRACE and
SM radius parameters. These are read from the `enkfpf.par` input file
via the existing namelist mechanism; add corresponding entries to the
`[DA]` section parser in the same module.

---

## Step 4 — Register the New Type in `mod_read_obs.F90`

The routine `update_obs_type` reads the `type_clm` string from the
observation file and sets the `clmupdate_*` flags for the next
analysis cycle. Add a `'C'` case that sets the appropriate flags. If
no existing `clmupdate_*` flag covers your new variable, declare a new
`clmupdate_C` in `enkf_clm_mod.F90` and wire it to the state-vector
fill and scatter logic in the eCLM interface routines (`eclm/`).

---

## Key Design Differences Between GRACE-DA and SM-DA

### Localization distance units

| Type  | `disttype`    | Coordinate arrays          | Radius unit |
|-------|---------------|----------------------------|-------------|
| GRACE | 0 (Cartesian) | integer grid-cell indices  | grid cells  |
| SM    | 3 (haversine) | degrees longitude/latitude | km          |

### Observation operator complexity

- **SM**: direct extraction from `state_p` using a precomputed index.
- **GRACE**: aggregating — all gridcell values within the GRACE
  footprint are summed and averaged, then the temporal mean is
  subtracted to produce a TWS anomaly. `thisobs%id_obs_p(1, g)` maps
  each gridcell `g` to its observation (inverse of the usual
  convention).

Follow the SM pattern for C point observations; follow the GRACE
pattern if observations represent spatial averages over a large
footprint.

### Coverage checks and error specification

GRACE discards observations where fewer than a minimum number of
gridcells lie within the support radius (important near coastlines). SM
and C point observations do not need this, but consider rejecting
observations over water bodies or outside the CLM domain.

Error specification is controlled by `multierr`: constant
(`multierr=0`), per-observation from `obserr_clm` (`multierr=1`), or
full covariance matrix from `obscov_clm` (`multierr=2`).

---

## Assimilating Multiple Observation Types at the Same Timestep

The framework generates a state vector for each type individually before the assimilation, some things would need to be adapted when mutliple observation types are assimilated at the same timestep. Currently, one observation file only consists of one observation type. As SM observations are usually assimilated at noon and GRACE observations are assimilated at the end of the month at midnight, this should not provide any problems.

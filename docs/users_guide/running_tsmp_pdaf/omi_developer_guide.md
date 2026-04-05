(omi:developer-guide)=
# Developer Guide: Adding a New Observation Type to the OMI Interface

This page is aimed at a software developer who wants to add a new
observation type — for example Land Surface Temperature (LST) — to the
TSMP-PDAF OMI (Observation Model Interface) framework.  It describes
the architecture of the existing interface, explains the roles of each
source file, and gives a step-by-step walkthrough of every change
required, using the existing GRACE-DA and SM-DA implementations as
concrete reference points.

## Architecture Overview

The observation handling sits between the PDAF library and the eCLM
model state.  Three layers are involved:

```
PDAF library (src/)
    │  calls call-back routines at each analysis step
    ▼
callback_obs_pdafomi.F90        ← single routing hub
    │  dispatches to one routine per observation type
    ▼
obs_GRACE_pdafomi.F90           ← one module per observation type
obs_SM_pdafomi.F90
obs_LST_pdafomi.F90  (new)
    │  read observations, build H, apply H to state vector
    ▼
eCLM state vector (via enkf_clm_mod, mod_assimilation, …)
```

Each observation-type module is independent and self-contained.  The
only coupling between modules is through the central routing file and
the shared `mod_assimilation` module where observation radii and
similar parameters live.

## Files Involved

To add a new observation type `LST` you need to touch exactly **five**
things:

| File                                           | Action                                                                      |
|------------------------------------------------|-----------------------------------------------------------------------------|
| `interface/framework/obs_LST_pdafomi.F90`      | **Create** (new file, modelled on the template)                             |
| `interface/framework/callback_obs_pdafomi.F90` | Add calls for LST to every callback routine                                 |
| `interface/framework/mod_assimilation.F90`     | Declare `cradius_LST`, `sradius_LST`, and any other LST-specific parameters |
| `interface/framework/mod_read_obs.F90`         | Add `'LST'` case to `update_obs_type`                                       |
| `interface/framework/enkf_clm_mod.F90`         | Add `clmupdate_LST` flag (parallel to `clmupdate_tws`, `clmupdate_swc`)     |

The template at `templates/omi/obs_OBSTYPE_pdafomi_TEMPLATE.F90` is a
good starting skeleton, but the existing SM and GRACE modules show
how the template is actually applied in this codebase.

---

## Step 1 — Create `obs_LST_pdafomi.F90`

### Module-level variables

Every observation module declares at minimum:

```fortran
MODULE obs_LST_pdafomi
  USE PDAFomi, ONLY: obs_f, obs_l
  IMPLICIT NONE
  SAVE
  PUBLIC

  LOGICAL :: assim_LST       ! controlled by callback_obs_pdafomi
  REAL    :: rms_obs_LST     ! constant observation error std dev (fallback)

  TYPE(obs_f), TARGET, PUBLIC :: thisobs    ! full obs data (set in init_dim_obs)
  TYPE(obs_l), TARGET, PUBLIC :: thisobs_l  ! local obs data (set by PDAFomi)
  !$OMP THREADPRIVATE(thisobs_l)
CONTAINS
  ...
END MODULE obs_LST_pdafomi
```

Additional module-level arrays — for example to cache the LST
climatology mean needed for an anomaly operator, or coordinate arrays
for the spatial matching — can be added as `ALLOCATABLE` module
variables, following the pattern in `obs_GRACE_pdafomi.F90`:

```fortran
! GRACE example — cache temporal mean loaded once per run
real, allocatable :: tws_temp_mean_d(:)
```

### The `obs_f` type

`obs_f` is the PDAF-OMI data type that carries all observation
metadata.  You set its fields inside `init_dim_obs_LST`.  The
mandatory fields are:

| Field                   | Type                   | Meaning                                               |
|-------------------------|------------------------|-------------------------------------------------------|
| `thisobs%doassim`       | `INTEGER`              | 1 = assimilate, 0 = skip                              |
| `thisobs%disttype`      | `INTEGER`              | distance metric for localization (see below)          |
| `thisobs%ncoord`        | `INTEGER`              | number of spatial coordinates (usually 2)             |
| `thisobs%id_obs_p(:,:)` | `INTEGER, ALLOCATABLE` | state-vector index for each process-local observation |

Optional fields that are frequently useful:

| Field                   | Default | Meaning                                                                      |
|-------------------------|---------|------------------------------------------------------------------------------|
| `thisobs%icoeff_p(:,:)` | —       | bilinear interpolation weights (SM uses this)                                |
| `thisobs%obs_err_type`  | 0       | 0 = Gaussian, 1 = Laplace                                                    |
| `thisobs%inno_omit`     | 0.0     | omit obs whose squared innovation exceeds this factor × obs variance         |
| `thisobs%infile`        | —       | custom flag used here to signal "no obs in file" (not a standard PDAF field) |

**`disttype` choices and units:**
- `0` — Cartesian distance, coordinates in any consistent unit (GRACE
  uses integer grid-cell indices, so `cradius_GRACE` is in grid cells)
- `3` — geographic distance via the haversine formula, coordinates in
  degrees (SM uses this, so `cradius_SM` is in km)

Choose `disttype=3` for point observations spread across the globe
(LST, SM) and `disttype=0` if your observations are already mapped to
integer grid-cell coordinates (coarser products like GRACE).

---

### Routine 1 — `init_dim_obs_LST`

This is the most complex routine.  It is called once per analysis
step, before the filter loop.  Its responsibilities are:

1. Set the mandatory `thisobs` fields.
2. Read the global observation vector from the NetCDF file.
3. Determine which observations fall in the PE-local domain and build
   the index array `thisobs%id_obs_p`.
4. Assemble the PE-local arrays `obs_p`, `ivar_obs_p`, `ocoord_p`.
5. Call `PDAFomi_gather_obs` to exchange information across MPI ranks.

**Reading observations from file**

Both SM and GRACE use the shared reader `read_obs_nc_type` from
`mod_read_obs`:

```fortran
USE mod_read_obs, ONLY: read_obs_nc_type
character(len=20) :: obs_type_name
obs_type_name = 'LST'
write(current_obs_filename, '(a, i5.5)') trim(obs_filename)//'.',  step
call read_obs_nc_type(current_obs_filename, obs_type_name, &
                      dim_obs, obs_g, lon_obs, lat_obs, layer_obs, &
                      dr_obs, obserr, obscov)
```

The reader checks the `type_clm` NetCDF variable in the file and
returns `dim_obs = 0` if the file does not contain observations of the
requested type.  Always handle the `dim_obs == 0` case — it is the
normal situation on time steps when no LST observations are scheduled:

```fortran
if (dim_obs == 0) then
  dim_obs_p = 0
  ! allocate dummy arrays of size 1 (required by PDAFomi_gather_obs)
  ALLOCATE(obs_p(1), ivar_obs_p(1), ocoord_p(2,1), thisobs%id_obs_p(1,1))
  thisobs%infile = 0
  CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
       thisobs%ncoord, cradius_LST, dim_obs)
  DEALLOCATE(obs_g, obs_p, ocoord_p, ivar_obs_p)
  return
end if
thisobs%infile = 1
```

Note the difference between GRACE and SM here: GRACE calls
`read_obs_nc_type` on **all** filter processes because the GRACE
spatial-averaging operator needs each PE to know which observations
it contributes to.  SM reads only on `mype_filter==0` and then
broadcasts `dim_obs` and the observation arrays.  For point
observations like LST, broadcasting (as SM does) is simpler.

**Snapping observations to the state vector**

For point observations (SM, LST), you loop over all global
observations `i` and all gridcells `g` on the PE-local domain, and
record the state-vector index for each observation that falls in the
local domain.  SM uses the helper arrays `longxy`, `latixy` (integer
grid indices) for the matching and stores the result in `obs_index_p`
from `mod_assimilation`:

```fortran
! SM pattern (simplified)
cnt = 0
do i = 1, dim_obs
  do g = begg, endg
    if (longxy_obs(i) == longxy(g-begg+1) .and. &
        latixy_obs(i) == latixy(g-begg+1)) then
      cnt = cnt + 1
      obs_index_p(cnt) = state_clm2pdaf_p(c, layer_obs(i))  ! CLMFIVE
      obs_p(cnt)       = obs_g(i)
      ivar_obs_p(cnt)  = 1.0 / (rms_obs_LST**2)
      ocoord_p(1,cnt)  = lon_obs(i)
      ocoord_p(2,cnt)  = lat_obs(i)
    end if
  end do
end do
dim_obs_p = cnt
```

For GRACE the matching is inverted: each observation aggregates
**multiple** gridcells, so the loop runs over gridcells and then
checks which observation each gridcell contributes to.
`thisobs%id_obs_p` maps each gridcell to the observation it belongs
to (1-to-many mapping), which is why the GRACE operator calls
`PDAFomi_gather_obsstate` with a gridcell-level sum rather than a
direct state element extraction.

**Calling `PDAFomi_gather_obs`**

This PDAFomi routine finalises the observation bookkeeping across MPI
ranks.  It must be called after all PE-local arrays are populated:

```fortran
CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
     thisobs%ncoord, cradius_LST, dim_obs)
```

After this call, `thisobs%dim_obs_p` holds the count of
process-local observations for use in `obs_op_LST`.

---

### Routine 2 — `obs_op_LST`

This routine applies the observation operator H: it maps the
PE-local model state `state_p` to the space of all observations
`ostate`.  For a simple point observation with no horizontal
interpolation (the SM case):

```fortran
SUBROUTINE obs_op_LST(dim_p, dim_obs, state_p, ostate)
  USE mod_assimilation, ONLY: obs_index_p
  USE PDAFomi_obs_f,    ONLY: PDAFomi_gather_obsstate
  IMPLICIT NONE
  INTEGER, INTENT(in) :: dim_p, dim_obs
  REAL,    INTENT(in)    :: state_p(dim_p)
  REAL,    INTENT(inout) :: ostate(dim_obs)
  REAL, ALLOCATABLE :: ostate_p(:)
  INTEGER :: i

  IF (thisobs%dim_obs_p > 0) THEN
    ALLOCATE(ostate_p(thisobs%dim_obs_p))
  ELSE
    ALLOCATE(ostate_p(1))
  END IF

  DO i = 1, thisobs%dim_obs_p
    ostate_p(i) = state_p(obs_index_p(i))   ! direct extraction
  END DO

  CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)
  DEALLOCATE(ostate_p)
END SUBROUTINE obs_op_LST
```

For GRACE, the operator is more complex: it must sum all gridcell
contributions that map to the same observation, subtract the temporal
mean, and then call `PDAFomi_gather_obsstate`.  If your new variable
requires a similar aggregation or a nonlinear transformation
(e.g. converting a model skin temperature to a brightness
temperature), implement that here.

The key invariant is that after `PDAFomi_gather_obsstate` returns,
`ostate` contains the model equivalent for **all** observations,
placed at the correct offsets determined by the order of
`init_dim_obs_*` calls in `callback_obs_pdafomi.F90`.

---

### Routines 3 & 4 — Localization (required for LESTKF/LETKF)

If you use a domain-localized filter (LESTKF or LETKF), you must
implement two additional routines:

**`init_dim_obs_l_LST`** — counts observations within the
localization radius of the current local analysis domain.  This is a
thin wrapper around `PDAFomi_init_dim_obs_l`:

```fortran
SUBROUTINE init_dim_obs_l_LST(domain_p, step, dim_obs, dim_obs_l)
  USE PDAFomi,        ONLY: PDAFomi_init_dim_obs_l
  USE mod_assimilation, ONLY: cradius_LST, locweight, sradius_LST
  ...
  ! coords_l: (ncoord) array with coordinates of local analysis domain
  CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
       locweight, cradius_LST, sradius_LST, dim_obs_l)
END SUBROUTINE
```

The local analysis domain coordinates `coords_l` must be in the same
units and with the same `disttype` as the observation coordinates set
in `init_dim_obs_LST`.

**`localize_covar_LST`** — required only for the LEnKF; a thin
wrapper around `PDAFomi_localize_covar`:

```fortran
SUBROUTINE localize_covar_LST(dim_p, dim_obs, HP_p, HPH, coords_p)
  USE PDAFomi, ONLY: PDAFomi_localize_covar
  USE mod_assimilation, ONLY: cradius_LST, locweight, sradius_LST
  ...
  CALL PDAFomi_localize_covar(thisobs, locweight, cradius_LST, &
       sradius_LST, coords_p, HP_p, HPH)
END SUBROUTINE
```

---

## Step 2 — Update `callback_obs_pdafomi.F90`

This file is the single hub that dispatches PDAF callbacks to each
observation module.  Add calls for LST in the same pattern as GRACE
and SM.  The commented-out `!USE obs_ST_pdafomi` lines are placeholder
reminders from the original code showing exactly where to add new
entries.

**`init_dim_obs_pdafomi`:**

```fortran
USE obs_LST_pdafomi, ONLY: assim_LST, init_dim_obs_LST
INTEGER :: dim_obs_LST = 0
assim_LST = .true.
IF (assim_LST) CALL init_dim_obs_LST(step, dim_obs_LST)
dim_obs = dim_obs_GRACE + dim_obs_SM + dim_obs_LST
```

Keep `assim_LST = .true.` here unconditionally, because the
observation module itself returns `dim_obs = 0` when the file
contains no LST observations.  The `assim_*` flags exist mainly to
let the user disable an observation type at compile/run time for
debugging.

**`obs_op_pdafomi`:**

```fortran
USE obs_LST_pdafomi, ONLY: obs_op_LST
CALL obs_op_LST(dim_p, dim_obs, state_p, ostate)
```

The order of `obs_op_*` calls does not matter because the offsets
within `ostate` are determined by the order of `init_dim_obs_*` calls
above.

**`init_dim_obs_l_pdafomi`:**

```fortran
USE obs_LST_pdafomi, ONLY: init_dim_obs_l_LST
CALL init_dim_obs_l_LST(domain_p, step, dim_obs, dim_obs_l)
```

Note that `dim_obs_l` is incremented inside each `init_dim_obs_l_*`
call, so the calls must be additive and the variable should be
initialized to 0 before the first call (this is handled by PDAF
itself before invoking this callback).

**`localize_covar_pdafomi`:**

```fortran
USE obs_LST_pdafomi, ONLY: localize_covar_LST
CALL localize_covar_LST(dim_p, dim_obs, HP_p, HPH, coords_p)
```

The `coords_p` array (PE-local state-vector coordinates) is allocated
in `localize_covar_pdafomi` and shared across all `localize_covar_*`
calls.  Make sure it is populated consistently with the coordinate
system expected by your `localize_covar_LST`.

**The remaining callbacks** (`add_obs_err_pdafomi`,
`init_obscovar_pdafomi`, `prodRinvA_pdafomi`, `prodRinvA_l_pdafomi`,
`deallocate_obs_pdafomi`) each need a corresponding LST call, which
are all thin wrappers around the matching PDAFomi generic routine,
identical in structure to the GRACE and SM entries.

---

## Step 3 — Add Parameters to `mod_assimilation.F90`

Declare the localization radii and any other LST-specific scalar
parameters next to the existing GRACE and SM entries (lines 276–279):

```fortran
REAL :: cradius_LST    ! cut-off radius for LST localization (km, for disttype=3)
REAL :: sradius_LST    ! support radius for 5th-order polynomial localization
```

These are read from the `enkfpf.par` input file via the existing
namelist mechanism.  Add corresponding entries to the `[DA]` section
parser in the same module.

---

## Step 4 — Register the New Type in `mod_read_obs.F90`

The routine `update_obs_type` (line 955) reads the `type_clm` string
from the observation file and sets the `clmupdate_*` flags accordingly.
The framework calls this between analysis steps to configure the model
state updates for the next cycle.  Add a `'LST'` case:

```fortran
case ('LST')
    clmupdate_tws     = 0
    clmupdate_swc     = 0
    clmupdate_T       = 1    ! or whichever flag controls LST updates
    clmupdate_texture = 0
```

If no existing `clmupdate_*` flag fits your new variable, declare a
new one (`clmupdate_LST`) in `enkf_clm_mod.F90` and add the
corresponding state-vector filling/scattering logic in the eCLM
interface routines (`eclm/`).

---

## Step 5 — Prepare the Observation NetCDF Files

The reader `read_obs_nc_type` (line 81 of `mod_read_obs.F90`) expects
the following variables in each observation file:

| NetCDF name | Fortran type | Dimension | Required? |
|---|---|---|---|
| `dim_obs` | dimension | — | yes |
| `obs_clm` | `REAL(dim_obs)` | observation values | yes |
| `type_clm` | `CHARACTER(20)(dim_obs)` | observation type string, e.g. `"LST"` | yes |
| `lon` | `REAL(dim_obs)` | longitude (degrees east) | yes |
| `lat` | `REAL(dim_obs)` | latitude (degrees north) | yes |
| `layer` | `INTEGER(dim_obs)` | soil/canopy layer index (1-based) | yes |
| `dr` | `REAL(1)` | search radius for matching (degrees) | yes |
| `obserr_clm` | `REAL(dim_obs)` | per-observation error std dev | if `multierr=1` |
| `obscov_clm` | `REAL(dim_obs,dim_obs)` | full error covariance | if `multierr=2` |

File naming follows the convention `<obs_filename>.<NNNNN>` where
`<NNNNN>` is the zero-padded PDAF time step (five digits), consistent
with the pattern used in `init_dim_obs_GRACE` and `init_dim_obs_SM`:

```fortran
write(current_obs_filename, '(a, i5.5)') trim(obs_filename)//'.', step
```

The `obs_filename` stem is a single shared prefix for all observation
types; the reader distinguishes them using `type_clm`.  This means a
file containing GRACE observations at one time step cannot
simultaneously contain SM or LST observations — the `type_clm` field
must be homogeneous within a file.

---

## Key Design Differences Between GRACE-DA and SM-DA

Understanding these differences helps avoid common pitfalls when
designing a new operator.

### Localization distance units

| Type  | `disttype`    | Coordinate arrays          | Radius unit |
|-------|---------------|----------------------------|-------------|
| GRACE | 0 (Cartesian) | integer grid-cell indices  | grid cells  |
| SM    | 3 (haversine) | degrees longitude/latitude | km          |

If you use `disttype=3`, pass true geographic coordinates
(degrees) in `ocoord_p` and set `cradius_*` in km.

### Observation operator complexity

- **SM** (`obs_op_SM`): trivial — a direct extraction from
  `state_p` using a precomputed index `obs_index_p(i)`.
- **GRACE** (`obs_op_GRACE`): aggregating — all gridcell values
  within the GRACE footprint are summed and averaged; the temporal
  mean is then subtracted to produce a TWS anomaly.  The mapping from
  gridcell to observation is stored in `thisobs%id_obs_p(1, g)` (where
  `g` is the gridcell index, not the observation index), which is the
  inverse of the usual convention.

For LST-DA with point observations, follow the SM pattern.  If LST
observations represent spatial averages over a footprint (e.g. thermal
infrared with a large IFOV), follow the GRACE pattern.

### Who reads the observation file

- **SM**: only `mype_filter==0` reads; all other ranks allocate
  arrays and receive data via `MPI_Bcast`.
- **GRACE**: all ranks call `read_obs_nc_type` independently (each
  gets `dim_obs` and the global observation arrays), then compute
  their PE-local contribution to the spatial average.  This avoids a
  broadcast of potentially large temporary arrays.

### Handling observations with insufficient model coverage

GRACE applies a minimum-coverage check: an observation is discarded
if fewer than `numPoints = ⌈π(dr/0.11)²/2⌉` model gridcells lie
within its support radius.  This is important for GRACE footprints
near coastlines.  SM and LST observations are point-like and do not
need this check, but you should decide whether your observation is
hydrologically "observable" in the model (e.g. reject observations
over water bodies or outside the CLM domain) and add the corresponding
guard in `init_dim_obs_LST`.

### Observation error specification

All three modes controlled by `multierr` (from `mod_assimilation`) are
available:

- `multierr=0`: constant error `rms_obs_LST` for all observations
- `multierr=1`: per-observation error read from `obserr_clm` in the
  NetCDF file
- `multierr=2`: full covariance matrix `obscov_clm`; PDAF then calls
  `prodRinvA` with the precomputed matrix inverse

GRACE supports all three modes because the TWS error covariance is
often spatially correlated.  SM currently uses only modes 0 and 1.

---

## Checklist Summary

When adding a new observation type `LST`:

- [ ] Create `interface/framework/obs_LST_pdafomi.F90`
  - [ ] Implement `init_dim_obs_LST` (read obs, snap to state vector, call `PDAFomi_gather_obs`)
  - [ ] Implement `obs_op_LST` (apply H, call `PDAFomi_gather_obsstate`)
  - [ ] Implement `init_dim_obs_l_LST` (needed for LESTKF/LETKF)
  - [ ] Implement `localize_covar_LST` (needed for LEnKF)
  - [ ] Implement `add_obs_err_LST`, `init_obscovar_LST`, `prodRinvA_LST`, `prodRinvA_l_LST`, `deallocate_obs_LST`
- [ ] In `callback_obs_pdafomi.F90`: add `USE obs_LST_pdafomi` and calls in all 8 callback subroutines
- [ ] In `mod_assimilation.F90`: declare `cradius_LST`, `sradius_LST`, read from namelist
- [ ] In `mod_read_obs.F90`: add `'LST'` case in `update_obs_type`
- [ ] In `enkf_clm_mod.F90`: add `clmupdate_LST` flag (if needed) and wire it to state-vector fill/scatter logic
- [ ] Prepare observation NetCDF files with `type_clm = "LST"` and the required variables
- [ ] Add `cradius_LST` and `sradius_LST` to `enkfpf.par`

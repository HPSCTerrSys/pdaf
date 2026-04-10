(omi:tsmp-pdaf-with-pdaf-omi)=
# TSMP-PDAF with PDAF-OMI

PDAF-OMI for multivariate data assimiliaton.
It is design to handle different observation types (currently soil moisture and TWS) automatically.
Additional to current observation files, the observation type has to be included (SM, GRACE).

See also create observation script for details: https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp-pdaf-observation-scripts/-/tree/main/omi_obs_data?ref_type=heads

Both global and local filters can be used. To enable multi-scale data assimilation, different localization radii for different observation types can be passed. Note that the localization radius for SM is currently in km and for GRACE in #gridcells.

If questions arise contact ewerdwalbesloh@geod.uni-bonn.de

## Observation File Format

The reader `read_obs_nc_type` expects the following variables in each
observation file:

| NetCDF name  | Type                           | Required?        |
|--------------|--------------------------------|------------------|
| `dim_obs`    | dimension                      | yes              |
| `obs_clm`    | `REAL(dim_obs)`                | yes              |
| `type_clm`   | `CHARACTER(20)(dim_obs)`       | yes              |
| `lon`        | `REAL(dim_obs)`                | yes              |
| `lat`        | `REAL(dim_obs)`                | yes              |
| `layer`      | `INTEGER(dim_obs)`             | yes              |
| `dr`         | `REAL(1)`                      | yes              |
| `obserr_clm` | `REAL(dim_obs)`                | if `multierr=1`  |
| `obscov_clm` | `REAL(dim_obs,dim_obs)`        | if `multierr=2`  |

Files follow the naming convention `<obs_filename>.<NNNNN>` (five-digit
zero-padded PDAF step). The `type_clm` field must be homogeneous within
a file — a file containing GRACE observations cannot simultaneously
contain SM or C observations.

**Note**: Each entry `type_clm` should have the length `20`. This
means the actual types should be padded with whitespaces (before or
after). So, currently possible entries would be
- `"GRACE               "`
- `"SM                  "`
- `"LST                 "`

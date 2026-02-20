(omi:tsmp-pdaf-with-pdaf-omi)=
# TSMP-PDAF with PDAF-OMI

PDAF-OMI for multivariate data assimiliaton.
It is design to handle different observation types (currently soil moisture and TWS) automatically.
Additional to current observation files, the observation type has to be included (SM, GRACE).

See also create observation script for details: https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp-pdaf-observation-scripts/-/tree/main/omi_obs_data?ref_type=heads

Both global and local filters can be used. To enable multi-scale data assimilation, different localization radii for different observation types can be passed. Note that the localization radius for SM is currently in km and for GRACE in #gridcells.

The framework generates a state vector for each type individually before the assimilation, some things would need to be adapted when mutliple observation types are assimilated at the same timestep. Currently, one observation file only consists of one observation type. As SM observations are usually assimilated at noon and GRACE observations are assimilated at the end of the month at midnight, this should not provide any problems.

If questions arise contact ewerdwalbesloh@geod.uni-bonn.de

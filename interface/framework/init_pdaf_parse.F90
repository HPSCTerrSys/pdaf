!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz and Guowei He (Forschungszentrum Juelich GmbH)
!
!This file is part of TSMP-PDAF
!
!TSMP-PDAF is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!TSMP-PDAF is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU LesserGeneral Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
!init_pdaf_parse.F90: TSMP-PDAF implementation of routine
!                     'init_pdaf_parse' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_pdaf_parse.F90 1442 2013-10-04 10:35:19Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf_parse - Parse command line options for PDAF
!
! !INTERFACE:
SUBROUTINE init_pdaf_parse()

! !DESCRIPTION:
! This routine calls the command line parser to initialize
! variables for the data assimilation with PDAF.
! Using the parser is optional and shows one possibility
! to modify the variables of the compiled program. An
! alternative to this might be Fortran namelist files.
!
! !REVISION HISTORY:
! 2011-05 - Lars Nerger - Initial code extracted from init_pdaf
! Later revisions - see svn log
!
! !USES:
  USE parser, &           ! Parser function
       ONLY: parse
  USE mod_parallel_pdaf, &     ! Parallelization variables
       ONLY: mype_world
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, toffset, &
       rms_obs, model_error, model_err_amp, incremental, type_forget, &
       forget, rank_analysis_enkf, locweight, cradius, &
       sradius, filename, type_trans, dim_obs, &
       type_sqrt, obs_filename, dim_lag
  USE mod_assimilation, ONLY: temp_mean_filename
  USE mod_assimilation, ONLY: use_omi

  use mod_assimilation,&
       only: cradius_GRACE, sradius_GRACE, &
       cradius_SM, sradius_SM
  use mod_assimilation, only: cradius LST
  use mod_assimilation, only: sradius_LST
#ifdef CLMFIVE
  use obs_GRACE_pdafomi, only: rms_obs_GRACE
  use obs_SM_pdafomi, only: rms_obs_SM
  use obs_LST_pdafomi, only: rms_obs_LST
#endif
#if defined CLMSA
#ifdef CLMFIVE
  use enkf_clm_mod, only: use_omi_model
#endif
#endif


  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: init_pdaf
! Calls: parse
!EOP

! Local variables
  CHARACTER(len=32) :: handle  ! handle for command line parser


! **********************************
! *** Parse command line options ***
! **********************************

  ! Settings for model and time stepping
  handle = 'model_error'             ! Control application of model error
  CALL parse(handle, model_error)
  handle = 'model_err_amp'           ! Amplitude of model error
  CALL parse(handle, model_err_amp)

  ! Observation settings
  handle = 'delt_obs'                ! Time step interval between filter analyses
  CALL parse(handle, delt_obs)
  handle = 'toffset'                ! Offset in time steps
  CALL parse(handle, toffset)
  handle = 'rms_obs'                 ! Assumed uniform RMS error of the observations
  CALL parse(handle, rms_obs)

#ifdef CLMFIVE
  rms_obs_GRACE = rms_obs  ! backward compatibility
  handle = 'rms_obs_GRACE'          ! RMS error for GRACE observations
  CALL parse(handle, rms_obs_GRACE)
  rms_obs_SM = rms_obs              ! backward compatibility
  handle = 'rms_obs_SM'             ! RMS error for SM observations
  CALL parse(handle, rms_obs_SM)
  rms_obs_LST = rms_obs              ! backward compatibility
  handle = 'rms_obs_LST'             ! RMS error for LST observations
  CALL parse(handle, rms_obs_LST)
  ! rms_obs_C = rms_obs              ! backward compatibility
  ! handle = 'rms_obs_C'             ! RMS error for C observations
  ! CALL parse(handle, rms_obs_C)
#endif

  handle = 'dim_obs'                 ! Number of observations
  CALL parse(handle, dim_obs)

  ! General settings for PDAF
  handle = 'screen'                  ! set verbosity of PDAF
  CALL parse(handle, screen)
  handle = 'dim_ens'                 ! set ensemble size/rank of covar matrix
  CALL parse(handle, dim_ens)
  handle = 'filtertype'              ! Choose filter algorithm
  CALL parse(handle, filtertype)
  handle = 'subtype'                 ! Set subtype of filter
  CALL parse(handle, subtype)
  handle = 'incremental'             ! Set whether to use incremental updating
  CALL parse(handle, incremental)
  handle = 'use_omi'                 ! Set whether to use OMI interface
  CALL parse(handle, use_omi)
#if defined CLMSA
#ifdef CLMFIVE
  use_omi_model = use_omi       ! Set variable for use in interface/model routines
#endif
#endif

  ! Filter-specific settings
  handle = 'type_trans'              ! Type of ensemble transformation in SEIK/ETKF/LSEIK/LETKF
  CALL parse(handle, type_trans)
  handle = 'rank_analysis_enkf'      ! Set rank for pseudo inverse in EnKF
  CALL parse(handle, rank_analysis_enkf)
  handle = 'type_forget'             ! Set type of forgetting factor
  CALL parse(handle, type_forget)
  handle = 'forget'                  ! Set forgetting factor
  CALL parse(handle,forget)
  handle = 'type_sqrt'               ! Set type of transformation square-root (SEIK-sub4, ESTKF)
  CALL parse(handle, type_sqrt)

  ! Settings for localization in LSEIK/LETKF
  handle = 'local_range'             ! For backward compatibility
  CALL parse(handle, cradius)
  handle = 'cradius'                 ! Set cut-off radius in grid points for observation domain
  CALL parse(handle, cradius)
  handle = 'locweight'               ! Set type of localizating weighting
  CALL parse(handle, locweight)
  sradius = cradius                  ! By default use cradius as support radius
  handle = 'srange'                  ! For backward compatibility
  CALL parse(handle, sradius)
  handle = 'sradius'                 ! Set support radius in grid points
             ! for 5th-order polynomial or radius for 1/e in exponential weighting
  CALL parse(handle, sradius)

  ! Settings for different observation types
  cradius_GRACE = cradius  ! For backward compatibility
  handle = 'cradius_GRACE'          ! Set cut-off radius for GRACE observations
  call parse(handle, cradius_GRACE)
  sradius_GRACE = sradius  ! For backward compatibility
  handle = 'sradius_GRACE'          ! Set support radius for GRACE observations
  call parse(handle, sradius_GRACE)
  cradius_SM = cradius              ! For backward compatibility
  handle = 'cradius_SM'             ! Set cut-off radius for SM observations
  call parse(handle, cradius_SM)
  sradius_SM = sradius              ! For backward compatibility
  handle = 'sradius_SM'             ! Set support radius for SM observations
  call parse(handle, sradius_SM)
  cradius_LST = cradius              ! For backward compatibility
  handle = 'cradius_LST'             ! Set cut-off radius for LST observations
  call parse(handle, cradius_LST)
  sradius_LST = sradius              ! For backward compatibility
  handle = 'sradius_LST'             ! Set support radius for LST observations
  call parse(handle, sradius_LST)

  ! Setting for file output
  handle = 'filename'                ! Set name of output file
  CALL parse(handle, filename)

  ! *** user defined observation filename *** !
  handle = 'obs_filename'
  call parse(handle, obs_filename)

  ! *** Yorck: user defined filename for temporal mean of TWS to be subtracted in observation operator *** !
  handle = 'temp_mean_filename'
  call parse(handle, temp_mean_filename)

  !kuw: add smoother support
  handle = 'smoother_lag'
  call parse(handle, dim_lag)
  !kuw end

END SUBROUTINE init_pdaf_parse

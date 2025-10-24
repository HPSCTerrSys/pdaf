!> callback_obs_pdafomi
!!
!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!!   When adding an observation type, one has to add one module
!!   obs_OBSTYPE_pdafomi (based on the template obs_OBSTYPE_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files.
!!
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

    use enkf_clm_mod, only: clmupdate_swc, clmupdate_tws

    ! Include functions for different observations
    USE obs_GRACE_pdafomi, ONLY: assim_GRACE, init_dim_obs_GRACE
    USE obs_SM_pdafomi, ONLY: assim_SM, init_dim_obs_SM
    !USE obs_ST_pdafomi, ONLY: assim_C, init_dim_obs_C

    use mod_assimilation, only: screen
    USE mod_parallel_pdaf, only: mype_world

    IMPLICIT NONE

  ! *** Arguments ***
    INTEGER, INTENT(in)  :: step     !< Current time step
    INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

  ! *** Local variables ***
    INTEGER :: dim_obs_GRACE ! Observation dimensions
    INTEGER :: dim_obs_SM ! Observation dimensions
    !INTEGER :: dim_obs_C ! Observation dimensions


  ! *********************************************
  ! *** Initialize full observation dimension ***
  ! *********************************************

    ! Initialize number of observations
    dim_obs_GRACE = 0
    dim_obs_SM = 0
    !dim_obs_C = 0


    assim_SM = .true.
    assim_GRACE = .true.



    ! Call observation-specific routines
    ! The routines are independent, so it is not relevant
    ! in which order they are called

    if (mype_world==0 .and. screen > 2) then
      write(*,*)'Call dimension initialization'
    end if

    IF (assim_GRACE) CALL init_dim_obs_GRACE(step, dim_obs_GRACE)
    IF (assim_SM) CALL init_dim_obs_SM(step, dim_obs_SM)
    !IF (assim_C) CALL init_dim_obs_C(step, dim_obs_C)

    dim_obs = dim_obs_GRACE + dim_obs_SM! + dim_obs_C

  END SUBROUTINE init_dim_obs_pdafomi



  !-------------------------------------------------------------------------------
  !> Call-back routine for obs_op
  !!
  !! This routine calls the observation-specific
  !! routines obs_op_TYPE.
  !!
  SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

    ! Include functions for different observations
    USE obs_GRACE_pdafomi, ONLY: obs_op_GRACE
    USE obs_SM_pdafomi, ONLY: obs_op_SM
    !USE obs_C_pdafomi, ONLY: obs_op_C

    use mod_assimilation, only: screen
    USE mod_parallel_pdaf, only: mype_world

    IMPLICIT NONE

  ! *** Arguments ***
    INTEGER, INTENT(in) :: step                 !< Current time step
    INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
    REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state


  ! ******************************************************
  ! *** Apply observation operator H on a state vector ***
  ! ******************************************************

    ! The order of these calls is not relevant as the setup
    ! of the overall observation vector is defined by the
    ! order of the calls in init_dim_obs_pdafomi

    if (mype_world==0 .and. screen > 2) then
      write(*,*)'Call observation operators'
    end if

    CALL obs_op_GRACE(dim_p, dim_obs, state_p, ostate)
    CALL obs_op_SM(dim_p, dim_obs, state_p, ostate)
    !CALL obs_op_C(dim_p, dim_obs, state_p, ostate)

  END SUBROUTINE obs_op_pdafomi



  !-------------------------------------------------------------------------------
  !> Call-back routine for init_dim_obs_l
  !!
  !! This routine calls the routine PDAFomi_init_dim_obs_l
  !! for each observation type
  !!
  SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

    ! Include functions for different observations
    USE obs_GRACE_pdafomi, ONLY: init_dim_obs_l_GRACE
    USE obs_SM_pdafomi, ONLY: init_dim_obs_l_SM
    !USE obs_C_pdafomi, ONLY: init_dim_obs_l_C

    IMPLICIT NONE

  ! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step       !< Current time step
    INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
    INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector


  ! **********************************************
  ! *** Initialize local observation dimension ***
  ! **********************************************

    ! Call init_dim_obs_l specific for each observation

    CALL init_dim_obs_l_GRACE(domain_p, step, dim_obs, dim_obs_l)
    CALL init_dim_obs_l_SM(domain_p, step, dim_obs, dim_obs_l)
    !CALL init_dim_obs_l_C(domain_p, step, dim_obs, dim_obs_l)

  END SUBROUTINE init_dim_obs_l_pdafomi



  !-------------------------------------------------------------------------------
  !> Call-back routine for localize_covar
  !!
  !! This routine calls the routine PDAFomi_localize_covar
  !! for each observation type to apply covariance
  !! localization in the LEnKF.
  !!
  SUBROUTINE localize_covar_pdafomi(dim_p, dim_obs, HP_p, HPH)

    ! Include functions for different observations
    USE obs_GRACE_pdafomi, ONLY: localize_covar_GRACE
    USE obs_SM_pdafomi, ONLY: localize_covar_SM
    !USE obs_C_pdafomi, ONLY: localize_covar_C

    ! Include information on model grid
    IMPLICIT NONE

  ! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< number of observations
    REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH

  ! *** local variables ***
    INTEGER :: i, j, cnt               ! Counters
    REAL, ALLOCATABLE :: coords_p(:,:) ! Coordinates of PE-local state vector entries


  ! **********************
  ! *** INITIALIZATION ***
  ! **********************

    ! Initialize coordinate array

    ALLOCATE(coords_p(2, dim_p))


  ! *************************************
  ! *** Apply covariance localization ***
  ! *************************************

    ! Call localize_covar specific for each observation

    CALL localize_covar_GRACE(dim_p, dim_obs, HP_p, HPH, coords_p)
    CALL localize_covar_SM(dim_p, dim_obs, HP_p, HPH, coords_p)
    !CALL localize_covar_C(dim_p, dim_obs, HP_p, HPH, coords_p)


  ! ****************
  ! *** Clean up ***
  ! ****************

    DEALLOCATE(coords_p)

  END SUBROUTINE localize_covar_pdafomi

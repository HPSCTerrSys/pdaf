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

! Author: Yorck Ewerdwalbesloh

#ifdef CLMFIVE
  SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

    use enkf_clm_mod, only: clmupdate_swc, clmupdate_tws
    use enkf_clm_mod, only: clmupdate_T

    ! Include functions for different observations
    USE obs_GRACE_pdafomi, ONLY: assim_GRACE
    USE obs_GRACE_pdafomi, ONLY: init_dim_obs_GRACE
    USE obs_SM_pdafomi, ONLY: assim_SM
    USE obs_SM_pdafomi, ONLY: init_dim_obs_SM
    USE obs_LST_pdafomi, ONLY: assim_LST
    USE obs_LST_pdafomi, ONLY: init_dim_obs_LST
    !USE obs_ST_pdafomi, ONLY: assim_C
    !USE obs_ST_pdafomi, ONLY: init_dim_obs_C

    use mod_assimilation, only: screen
    USE mod_parallel_pdaf, only: mype_world

    IMPLICIT NONE

  ! *** Arguments ***
    INTEGER, INTENT(in)  :: step     !< Current time step
    INTEGER, INTENT(out) :: dim_obs  !< Dimension of full observation vector

  ! *** Local variables ***
    INTEGER :: dim_obs_GRACE ! Observation dimensions
    INTEGER :: dim_obs_SM ! Observation dimensions
    INTEGER :: dim_obs_LST ! Observation dimensions
    !INTEGER :: dim_obs_C ! Observation dimensions


  ! *********************************************
  ! *** Initialize full observation dimension ***
  ! *********************************************

    ! Initialize number of observations
    dim_obs_GRACE = 0
    dim_obs_SM = 0
    dim_obs_LST = 0
    !dim_obs_C = 0


    ! Call observation-specific routines
    ! The routines are independent, so it is not relevant
    ! in which order they are called

    if (mype_world==0 .and. screen > 2) then
      write(*,*)'Call dimension initialization'
    end if

#ifdef PDAF_DEBUG
    if (mype_world==0) then
      write(*,*)'PDAF-OMI-DEBUG: assim_GRACE=', assim_GRACE
      write(*,*)'PDAF-OMI-DEBUG: assim_SM=', assim_SM
      write(*,*)'PDAF-OMI-DEBUG: assim_LST=', assim_LST
      ! write(*,*)'PDAF-OMI assim_C=', assim_C
    end if
#endif

    IF (assim_GRACE) CALL init_dim_obs_GRACE(step, dim_obs_GRACE)
    IF (assim_SM) CALL init_dim_obs_SM(step, dim_obs_SM)
    IF (assim_LST) CALL init_dim_obs_LST(step, dim_obs_LST)
    !IF (assim_C) CALL init_dim_obs_C(step, dim_obs_C)

    dim_obs = dim_obs_GRACE + dim_obs_SM + dim_obs_LST! + dim_obs_C

  END SUBROUTINE init_dim_obs_pdafomi



  !-------------------------------------------------------------------------------
  !> Call-back routine for obs_op
  !!
  !! This routine calls the observation-specific
  !! routines obs_op_TYPE.
  !!
  SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

    ! Include functions for different observations
    USE obs_GRACE_pdafomi, ONLY: assim_GRACE
    USE obs_GRACE_pdafomi, ONLY: obs_op_GRACE
    USE obs_SM_pdafomi, ONLY: assim_SM
    USE obs_SM_pdafomi, ONLY: obs_op_SM
    USE obs_LST_pdafomi, ONLY: assim_LST
    USE obs_LST_pdafomi, ONLY: obs_op_LST
    !USE obs_C_pdafomi, ONLY: assim_C
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

    IF (assim_GRACE) CALL obs_op_GRACE(dim_p, dim_obs, state_p, ostate)
    IF (assim_SM) CALL obs_op_SM(dim_p, dim_obs, state_p, ostate)
    IF (assim_LST) CALL obs_op_LST(dim_p, dim_obs, state_p, ostate)
    !IF (assim_C) CALL obs_op_C(dim_p, dim_obs, state_p, ostate)

  END SUBROUTINE obs_op_pdafomi



  !-------------------------------------------------------------------------------
  !> Call-back routine for init_dim_obs_l
  !!
  !! This routine calls the routine PDAFomi_init_dim_obs_l
  !! for each observation type
  !!
  SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

    ! Include functions for different observations
    USE obs_GRACE_pdafomi, ONLY: assim_GRACE
    USE obs_GRACE_pdafomi, ONLY: init_dim_obs_l_GRACE
    USE obs_SM_pdafomi, ONLY: assim_SM
    USE obs_SM_pdafomi, ONLY: init_dim_obs_l_SM
    USE obs_LST_pdafomi, ONLY: assim_LST
    USE obs_LST_pdafomi, ONLY: init_dim_obs_l_LST
    !USE obs_C_pdafomi, ONLY: assim_C
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

    IF (assim_GRACE) CALL init_dim_obs_l_GRACE(domain_p, step, dim_obs, dim_obs_l)
    IF (assim_SM) CALL init_dim_obs_l_SM(domain_p, step, dim_obs, dim_obs_l)
    IF (assim_LST) CALL init_dim_obs_l_LST(domain_p, step, dim_obs, dim_obs_l)
    !IF (assim_C) CALL init_dim_obs_l_C(domain_p, step, dim_obs, dim_obs_l)

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
    USE obs_GRACE_pdafomi, ONLY: assim_GRACE
    USE obs_GRACE_pdafomi, ONLY: localize_covar_GRACE
    USE obs_SM_pdafomi, ONLY: assim_SM
    USE obs_SM_pdafomi, ONLY: localize_covar_SM
    USE obs_LST_pdafomi, ONLY: assim_LST
    USE obs_LST_pdafomi, ONLY: localize_covar_LST
    !USE obs_C_pdafomi, ONLY: assim_C
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

    IF (assim_GRACE) CALL localize_covar_GRACE(dim_p, dim_obs, HP_p, HPH, coords_p)
    IF (assim_SM) CALL localize_covar_SM(dim_p, dim_obs, HP_p, HPH, coords_p)
    IF (assim_LST) CALL localize_covar_LST(dim_p, dim_obs, HP_p, HPH, coords_p)
    !IF (assim_C) CALL localize_covar_C(dim_p, dim_obs, HP_p, HPH, coords_p)


  ! ****************
  ! *** Clean up ***
  ! ****************

    DEALLOCATE(coords_p)

  END SUBROUTINE localize_covar_pdafomi


  SUBROUTINE add_obs_err_pdafomi(step, dim_obs, C)

    USE obs_GRACE_pdafomi, ONLY: assim_GRACE
    USE obs_GRACE_pdafomi, ONLY: add_obs_err_GRACE
    USE obs_SM_pdafomi, ONLY: assim_SM
    USE obs_SM_pdafomi, ONLY: add_obs_err_SM
    USE obs_LST_pdafomi, ONLY: assim_LST
    USE obs_LST_pdafomi, ONLY: add_obs_err_LST
    !USE obs_C_pdafomi, ONLY: assim_C
    !USE obs_C_pdafomi, ONLY: add_obs_err_C

    IMPLICIT NONE

    INTEGER, INTENT(in) :: step                ! Current time step
    INTEGER, INTENT(in) :: dim_obs             ! Dimension of obs. vector
    REAL, INTENT(inout) :: C(dim_obs, dim_obs) ! Matrix to that the observation
                                              !    error covariance matrix is added


    INTEGER :: i          ! index of observation component
    REAL :: variance_obs  ! variance of observations
    IF (assim_GRACE) CALL add_obs_err_GRACE(step, dim_obs, C)
    IF (assim_SM) CALL add_obs_err_SM(step, dim_obs, C)
    IF (assim_LST) CALL add_obs_err_LST(step, dim_obs, C)
    !IF (assim_C) CALL add_obs_err_C(step, dim_obs, C)

  END SUBROUTINE add_obs_err_pdafomi


  SUBROUTINE init_obscovar_pdafomi(step, dim_obs, dim_obs_p, covar, m_state_p, isdiag)

    USE obs_GRACE_pdafomi, ONLY: assim_GRACE
    USE obs_GRACE_pdafomi, ONLY: init_obscovar_GRACE
    USE obs_SM_pdafomi, ONLY: assim_SM
    USE obs_SM_pdafomi, ONLY: init_obscovar_SM
    USE obs_LST_pdafomi, ONLY: assim_LST
    USE obs_LST_pdafomi, ONLY: init_obscovar_LST
    !USE obs_C_pdafomi, ONLY: assim_C
    !USE obs_C_pdafomi, ONLY: init_obscovar_C
    IMPLICIT NONE

    INTEGER, INTENT(in) :: step                ! Current time step
    INTEGER, INTENT(in) :: dim_obs             ! Dimension of observation vector
    INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of observation vector
    REAL, INTENT(out) :: covar(dim_obs, dim_obs) ! Observation error covariance matrix
    REAL, INTENT(in)  :: m_state_p(dim_obs_p)  ! PE-local observation vector
    LOGICAL, INTENT(out) :: isdiag             ! Whether the observation error covar. matrix is diagonal

    integer :: i
    REAL :: variance_obs


    IF (assim_GRACE) CALL init_obscovar_GRACE(step, dim_obs, dim_obs_p, covar, m_state_p, isdiag)
    IF (assim_SM) CALL init_obscovar_SM(step, dim_obs, dim_obs_p, covar, m_state_p, isdiag)
    IF (assim_LST) CALL init_obscovar_LST(step, dim_obs, dim_obs_p, covar, m_state_p, isdiag)
    !IF (assim_C) CALL init_obscovar_C(step, dim_obs, dim_obs_p, covar, m_state_p, isdiag)

  END SUBROUTINE init_obscovar_pdafomi


  SUBROUTINE prodRinvA_pdafomi(step, dim_obs_p, rank, obs_p, A_p, C_p)

    use obs_GRACE_pdafomi, ONLY: assim_GRACE
    use obs_GRACE_pdafomi, ONLY: prodRinvA_GRACE
    use obs_SM_pdafomi, ONLY: assim_SM
    use obs_SM_pdafomi, ONLY: prodRinvA_SM
    use obs_LST_pdafomi, ONLY: assim_LST
    use obs_LST_pdafomi, ONLY: prodRinvA_LST
    !use obs_C_pdafomi, ONLY: assim_C
    !use obs_C_pdafomi, ONLY: prodRinvA_C

    implicit none
    INTEGER, INTENT(in) :: step                ! Current time step
    INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
    INTEGER, INTENT(in) :: rank                ! Rank of initial covariance matrix
    REAL, INTENT(in)    :: obs_p(dim_obs_p)    ! PE-local vector of observations
    REAL, INTENT(in)    :: A_p(dim_obs_p,rank) ! Input matrix from analysis routine
    REAL, INTENT(out)   :: C_p(dim_obs_p,rank) ! Output matrix

    IF (assim_GRACE) CALL prodRinvA_GRACE(step, dim_obs_p, rank, obs_p, A_p, C_p)
    IF (assim_SM) CALL prodRinvA_SM(step, dim_obs_p, rank, obs_p, A_p, C_p)
    IF (assim_LST) CALL prodRinvA_LST(step, dim_obs_p, rank, obs_p, A_p, C_p)
    !IF (assim_C) CALL prodRinvA_C(step, dim_obs_p, rank, obs_p, A_p, C_p)

  END SUBROUTINE prodRinvA_pdafomi


  SUBROUTINE prodRinvA_l_pdafomi(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

    use obs_GRACE_pdafomi, ONLY: assim_GRACE
    use obs_GRACE_pdafomi, ONLY: prodRinvA_l_GRACE
    use obs_SM_pdafomi, ONLY: assim_SM
    use obs_SM_pdafomi, ONLY: prodRinvA_l_SM
    use obs_LST_pdafomi, ONLY: assim_LST
    use obs_LST_pdafomi, ONLY: prodRinvA_l_LST
    !use obs_C_pdafomi, ONLY: assim_C
    !use obs_C_pdafomi, ONLY: prodRinvA_l_C

    implicit none

    INTEGER, INTENT(in) :: domain_p             ! Current local analysis domain
    INTEGER, INTENT(in) :: step                 ! Current time step
    INTEGER, INTENT(in) :: dim_obs_l            ! Dimension of local observation vector
    INTEGER, INTENT(in) :: rank                 ! Rank of initial covariance matrix
    REAL, INTENT(in)    :: obs_l(dim_obs_l)     ! Local vector of observations
    REAL, INTENT(inout) :: A_l(dim_obs_l, rank) ! Input matrix from analysis routine
    REAL, INTENT(out)   :: C_l(dim_obs_l, rank) ! Output matrix

    IF (assim_GRACE) CALL prodRinvA_l_GRACE(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
    IF (assim_SM) CALL prodRinvA_l_SM(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
    IF (assim_LST) CALL prodRinvA_l_LST(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)
    !IF (assim_C) CALL prodRinvA_l_C(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

  END SUBROUTINE prodRinvA_l_pdafomi

  SUBROUTINE deallocate_obs_pdafomi()

    use obs_GRACE_pdafomi, ONLY: assim_GRACE
    use obs_GRACE_pdafomi, ONLY: deallocate_obs_GRACE
    use obs_SM_pdafomi, ONLY: assim_SM
    use obs_SM_pdafomi, ONLY: deallocate_obs_SM
    use obs_LST_pdafomi, ONLY: assim_LST
    use obs_LST_pdafomi, ONLY: deallocate_obs_LST
    !use obs_C_pdafomi, ONLY: assim_C
    !use obs_C_pdafomi, ONLY: deallocate_obs_C

    implicit none

    IF (assim_GRACE) CALL deallocate_obs_GRACE()
    IF (assim_SM) CALL deallocate_obs_SM()
    IF (assim_LST) CALL deallocate_obs_LST()
    !IF (assim_C) CALL deallocate_obs_C()

  END SUBROUTINE deallocate_obs_pdafomi

#endif

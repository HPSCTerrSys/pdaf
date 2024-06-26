! Copyright (c) 2014-2024 Paul Kirchgessner
!
! This file is part of PDAF.
!
! PDAF is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! PDAF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with PDAF.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id$
!BOP
!
! !ROUTINE: PDAF_lnetf_update --- Control analysis update of the LNETF
!
! !INTERFACE:
SUBROUTINE  PDAF_lnetf_update(step, dim_p, dim_obs_f, dim_ens, &
     state_p, Uinv, ens_p, type_forget, noise_type, noise_amp, &
     U_obs_op, U_init_dim_obs, U_init_obs_l, U_likelihood_l, &
     U_init_n_domains_p, U_init_dim_l, U_init_dim_obs_l, U_g2l_state, U_l2g_state, &
     U_g2l_obs, U_prepoststep, screen, subtype, &
     dim_lag, sens_p, cnt_maxlag, flag)

! !DESCRIPTION:
! Routine to control the analysis update of the LNETF.
!
! The analysis is performed by first preparing several
! global quantities on the PE-local domain, like the
! observed part of the state ensemble for all local
! analysis domains on the PE-local state domain.
! Then the analysis (PDAF\_lnetf\_analysis) is performed within
! a loop over all local analysis domains in the PE-local 
! state domain. In this loop, the local state and 
! observation dimensions are initialized and the global 
! state ensemble is restricted to the local analysis domain.
! In addition, the routine U\_prepoststep is called prior
! to the analysis and after the resampling outside of
! the loop over the local domains to allow the user
! to access the ensemble information.
!
! Variant for domain decomposition.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2014-05 - Paul Kirchgessner - Initial code based on LETKF
! Later revisions - see svn log
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  USE mpi
  USE PDAF_timer, &
       ONLY: PDAF_timeit, PDAF_time_temp
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filter, &
       ONLY: obs_member, type_trans, type_winf, limit_winf, &
       forget, inloop, member_save, debug
  USE PDAF_mod_filtermpi, &
       ONLY: mype, dim_ens_l, npes_filter, COMM_filter, MPIerr
  USE PDAF_analysis_utils, &
       ONLY: PDAF_print_domain_stats, PDAF_init_local_obsstats, PDAF_incr_local_obsstats, &
       PDAF_print_local_obsstats

  IMPLICIT NONE

! !ARGUMENTS:
! ! Variable naming scheme:
! !   suffix _p: Denotes a full variable on the PE-local domain
! !   suffix _l: Denotes a local variable on the current analysis domain
! !   suffix _f: Denotes a full variable of all observations required for the
! !              analysis loop on the PE-local domain
  INTEGER, INTENT(in) :: step         ! Current time step
  INTEGER, INTENT(in) :: dim_p        ! PE-local dimension of model state
  INTEGER, INTENT(out) :: dim_obs_f   ! PE-local dimension of observation vector
  INTEGER, INTENT(in) :: dim_ens      ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)        ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens, dim_ens)      ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens) ! PE-local ensemble matrix
  INTEGER, INTENT(in) :: type_forget  ! Select the type of forgetting factor
                ! 0 : Inflate forecast ensemble on all analysis domains
                ! 1 : Inflate forecast ensemble on observed domains only
                ! 2 : Inflate analysis ensemble on all analysis domains
                ! 3 : Inflate analysis ensemble on observed domains only
                ! 4 : Inflate paricle weights according to N_eff/N>=forget
  INTEGER, INTENT(in) :: noise_type   ! Type of pertubing noise
  REAL, INTENT(in) :: noise_amp       ! Amplitude of noise
  INTEGER, INTENT(in) :: screen       ! Verbosity flag
  INTEGER, INTENT(in) :: subtype      ! Filter subtype
  INTEGER, INTENT(inout) :: dim_lag   ! Status flag
  REAL, INTENT(inout) :: sens_p(dim_p, dim_ens, dim_lag) ! PE-local smoother ensemble
  INTEGER, INTENT(inout) :: cnt_maxlag ! Count number of past time steps for smoothing
  INTEGER, INTENT(inout) :: flag      ! Status flag

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_obs_op, &    ! Observation operator
       U_init_n_domains_p, & ! Provide number of local analysis domains
       U_init_dim_l, &       ! Init state dimension for local ana. domain
       U_init_dim_obs, &     ! Initialize dimension of observation vector
       U_init_dim_obs_l, &   ! Initialize dim. of obs. vector for local ana. domain
       U_init_obs_l, &       ! Init. observation vector on local analysis domain
       U_g2l_state, &        ! Get state on local ana. domain from global state
       U_l2g_state, &        ! Init full state from state on local analysis domain
       U_g2l_obs, &          ! Restrict full obs. vector to local analysis domain
       U_likelihood_l, &     ! Compute observation likelihood for an ensemble member
       U_prepoststep         ! User supplied pre/poststep routine

! !CALLING SEQUENCE:
! Called by: PDAF_put_state_lnetf
! Calls: U_prepoststep
! Calls: U_init_n_domains_p
! Calls: U_init_dim_obs
! Calls: U_obs_op
! Calls: U_init_dim_l
! Calls: U_init_dim_obs_l
! Calls: U_g2l_state
! Calls: U_l2g_state
! Calls: PDAF_generate_rndmat
! Calls: PDAF_lnetf_analysis
! Calls: PDAF_lnetf_smootherT
! Calls: PDAF_smoother_lnetf
! Calls: PDAF_inflate_ens
! Calls: PDAF_timeit
! Calls: PDAF_memcount
! Calls: MPI_reduce
!EOP

! *** local variables ***
  INTEGER :: i, member, row        ! Counters
  INTEGER :: domain_p              ! Counter for local analysis domain
  INTEGER, SAVE :: allocflag = 0   ! Flag whether first time allocation is done
  INTEGER, SAVE :: allocflag_l = 0 ! Flag whether first time allocation is done
  REAL    :: invdimens             ! Inverse global ensemble size
  INTEGER :: minusStep             ! Time step counter
  INTEGER :: n_domains_p           ! number of PE-local analysis domains
  REAL, ALLOCATABLE :: HX_f(:,:)   ! HX for PE-local ensemble
  REAL, ALLOCATABLE :: HXbar_f(:)  ! PE-local observed mean state
  REAL, ALLOCATABLE :: TA_l(:,:)   ! Local ensemble transform matrix
  REAL, ALLOCATABLE :: HX_noinfl_f(:,:) ! HX for smoother (without inflation)
  REAL, ALLOCATABLE :: TA_noinfl_l(:,:) ! TA for smoother (without inflation)
  REAL, ALLOCATABLE :: rndmat(:,:) ! random rotation matrix for ensemble trans.
  ! Variables on local analysis domain
  INTEGER :: dim_l                 ! State dimension on local analysis domain
  INTEGER :: dim_obs_l             ! Observation dimension on local analysis domain
  REAL, ALLOCATABLE :: ens_l(:,:)  ! State ensemble on local analysis domain
  REAL, ALLOCATABLE :: state_l(:)  ! Mean state on local analysis domain
  REAL :: invforget                ! inverse forgetting factor
  REAL, ALLOCATABLE :: n_eff(:)    ! Effective sample size for each local domain
  LOGICAL, ALLOCATABLE :: MASK(:)  ! Mask for effective sample sizes > 0
  REAL :: max_n_eff_l, min_n_eff_l ! PE-local min/max. effective ensemble sizes
  REAL :: max_n_eff, min_n_eff     ! Global min/max. effective ensemble sizes
  INTEGER :: cnt_small_svals       ! Counter for small values
  INTEGER :: subtype_dummy         ! Dummy variable to avoid compiler warning
  REAL :: avg_n_eff_l, avg_n_eff   ! Average effective sample size


! *************************************
! *** Prestep for forecast ensemble ***
! *************************************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- START'

  ! Initialize variable to prevent compiler warning
  subtype_dummy = subtype

  CALL PDAF_timeit(5, 'new')
  minusStep = - step  ! Indicate forecast by negative time step number
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a, i7)') 'PDAF', 'Call pre-post routine after forecast; step ', step
  ENDIF
  CALL U_prepoststep(minusStep, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Uinv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF ', '--- duration of prestep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Analysis ', ('-', i = 1, 55)
  END IF


! **************************************
! *** Preparation for local analysis ***
! **************************************

#ifndef PDAF_NO_UPDATE
  CALL PDAF_timeit(3, 'new')  ! Time for assimilation 
  CALL PDAF_timeit(4, 'new')  ! global preparation

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(3) dim_lag     ', dim_lag
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(4) noise_type  ', noise_type
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(5) type_forget ', type_forget
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(6) type_trans  ', type_trans
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_int(7) type_winf   ', type_winf

     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_real(1) forget     ', forget
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_real(2) limit_winf ', limit_winf
     WRITE (*,*) '++ PDAF-debug PDAF_lnetf_update', debug, &
          'Configuration: param_real(3) noise amp. ', noise_amp
  END IF

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call init_n_domains'

  ! Query number of analysis domains for the local analysis
  ! in the PE-local domain
  CALL PDAF_timeit(42, 'new')
  CALL U_init_n_domains_p(step, n_domains_p)
  CALL PDAF_timeit(42, 'old')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  n_domains_p', n_domains_p

  ! Initialize effective sample size
  ALLOCATE(n_eff(n_domains_p))
  n_eff = 0.0  !initialize 
  ALLOCATE(MASK(n_domains_p))

  IF (screen > 0) THEN
     IF (mype == 0) THEN
        WRITE (*, '(a, i7, 3x, a)') &
             'PDAF ', step, 'Assimilating observations - LNETF analysis using T-matrix'
     END IF
     IF (screen<3) THEN
        CALL PDAF_print_domain_stats(n_domains_p)
     ELSE
        WRITE (*, '(a, 5x, a, i6, a, i10)') &
             'PDAF', '--- PE-domain:', mype, ' number of analysis domains:', n_domains_p
     END IF
  END IF


! *** Local analysis: initialize global quantities ***

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call init_dim_obs'

  ! Get observation dimension for all observations required 
  ! for the loop of local analyses on the PE-local domain.
  CALL PDAF_timeit(43, 'new')
  CALL U_init_dim_obs(step, dim_obs_f)
  CALL PDAF_timeit(43, 'old')

  IF (screen > 2) THEN
     WRITE (*, '(a, 5x, a, i6, a, i10)') &
          'PDAF', '--- PE-Domain:', mype, &
          ' dimension of PE-local full obs. vector', dim_obs_f
  END IF

  IF (dim_lag>0) THEN
     CALL PDAF_timeit(44, 'new')

     ! For the smoother get observed uninflated ensemble
     ALLOCATE(HX_noinfl_f(dim_obs_f, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_f * dim_ens)

     CALL PDAF_timeit(27, 'new') !Apply obs_op 

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_letkf_update -- call obs_op', dim_ens, 'times for smoother'

     ENSS: DO member = 1,dim_ens
        ! Store member index to make it accessible with PDAF_get_obsmemberid
        obs_member = member

        CALL U_obs_op(step, dim_p, dim_obs_f, ens_p(:, member), HX_noinfl_f(:, member))
     END DO ENSS

     CALL PDAF_timeit(27, 'old') 
     CALL PDAF_timeit(44, 'old')
  END IF

  CALL PDAF_timeit(51, 'new')

  ! Apply covariance inflation to global ens (only if prior infl is chosen)
  ! returns the full ensemble with unchanged mean, but inflated covariance
  IF (type_forget==0 .OR. type_forget==1) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: PDAF_letkf_update', debug, &
          'Inflate ensemble: type_forget, forget', type_forget, forget

     CALL PDAF_timeit(14, 'new') !Apply forgetting factor
     CALL PDAF_inflate_ens(dim_p, dim_ens, state_p, ens_p, forget)
     CALL PDAF_timeit(14, 'old')
  ENDIF

  CALL PDAF_timeit(51, 'old')

  ! HX = [Hx_1 ... Hx_(r+1)] for full DIM_OBS_F region on PE-local domain
  ALLOCATE(HX_f(dim_obs_f, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_f * dim_ens)

  CALL PDAF_timeit(44, 'new')

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call obs_op', dim_ens, 'times'

  ENS: DO member = 1,dim_ens
     ! Store member index to make it accessible with PDAF_get_obsmemberid
     obs_member = member

     CALL PDAF_timeit(12, 'new') !Apply obs_op 
     CALL U_obs_op(step, dim_p, dim_obs_f, ens_p(:, member), HX_f(:, member))
     CALL PDAF_timeit(12, 'old') 
  END DO ENS

  CALL PDAF_timeit(44, 'old')

  CALL PDAF_timeit(51, 'new')
  CALL PDAF_timeit(11, 'new')

  ! *** Compute mean forecast state
  state_p = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_p
        state_p(row) = state_p(row) + invdimens * ens_p(row, member)
     END DO
  END DO

  CALL PDAF_timeit(11, 'old')

  ! *** Compute mean state of ensemble on PE-local observation space 
  ALLOCATE(HXbar_f(dim_obs_f))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_obs_f)

  HXbar_f = 0.0
  invdimens = 1.0 / REAL(dim_ens)
  DO member = 1, dim_ens
     DO row = 1, dim_obs_f
        HXbar_f(row) = HXbar_f(row) + invdimens * HX_f(row, member)
     END DO
  END DO

  ! Generate orthogonal random matrix with eigenvector (1,...,1)^T
  CALL PDAF_timeit(13, 'new')

  ALLOCATE(rndmat(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens**2)

  IF (type_trans==0) THEN
     IF (screen > 0 .AND. mype == 0) &
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- Initialize random transformation'
     CALL PDAF_generate_rndmat(dim_ens, rndmat, 2)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  rndmat', rndmat
  ELSE
     IF (screen > 0 .AND. mype == 0) &
          WRITE (*, '(a, 5x, a)') 'PDAF', '--- Initialize deterministic transformation'
     rndmat = 0.0
     DO i = 1, dim_ens
        rndmat(i,i) = 1.0
     END DO
  END IF

  CALL PDAF_timeit(13, 'old')
  CALL PDAF_timeit(51, 'old')

  CALL PDAF_timeit(4, 'old')


! ************************
! *** Perform analysis ***
! ************************

  CALL PDAF_timeit(6, 'new')

  ! Initialize counters for statistics on local observations
  CALL PDAF_init_local_obsstats()

!$OMP PARALLEL default(shared) private(dim_l, dim_obs_l, ens_l, state_l, TA_l, TA_noinfl_l, flag)

  ! Allocate ensemble transform matrix
  ALLOCATE(TA_l(dim_ens, dim_ens))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)
  TA_l = 0.0

  ! For smoother: Allocate ensemble transform matrix without inflation
  IF (dim_lag>0) THEN
     ALLOCATE(TA_noinfl_l(dim_ens, dim_ens))
     IF (allocflag == 0) CALL PDAF_memcount(3, 'r', dim_ens*dim_ens)

     TA_noinfl_l = 0.0
  ELSE
     ALLOCATE(TA_noinfl_l(1, 1))
  END IF

  ! initialize number of small singular values
  cnt_small_svals = 0

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- Enter local analysis loop'

!$OMP BARRIER
!$OMP DO firstprivate(cnt_maxlag) lastprivate(cnt_maxlag) schedule(runtime)
  localanalysis: DO domain_p = 1, n_domains_p    

     ! Set flag that we are in the local analysis loop
     inloop = .true.

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug: ', debug, &
             'PDAF_letkf_update -- local analysis for domain_p', domain_p
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call init_dim_l'
     END IF

     ! local state dimension
     CALL PDAF_timeit(45, 'new')
     CALL U_init_dim_l(step, domain_p, dim_l)
     CALL PDAF_timeit(45, 'old')

     IF (debug>0) THEN
        WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  dim_l', dim_l
        WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- call init_dim_obs_l'
     END IF

     ! Get observation dimension for local domain
     CALL PDAF_timeit(9, 'new')
     dim_obs_l = 0
     CALL U_init_dim_obs_l(domain_p, step, dim_obs_f, dim_obs_l)
     CALL PDAF_timeit(9, 'old')

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  dim_obs_l', dim_obs_l

     CALL PDAF_timeit(51, 'new')

     ! Gather statistical information on local observations
     CALL PDAF_incr_local_obsstats(dim_obs_l)

     ! Allocate arrays for local analysis domain
     ALLOCATE(ens_l(dim_l, dim_ens))
     ALLOCATE(state_l(dim_l))
     IF (allocflag_l == 0) CALL PDAF_memcount(3, 'r', dim_l*dim_ens + dim_l)
     CALL PDAF_timeit(51, 'old')

     CALL PDAF_timeit(15, 'new')

     ! state ensemble and mean state on current analysis domain
     DO member = 1, dim_ens
        ! Store member index to make it accessible with PDAF_get_memberid
        member_save = member

        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_letkf_update -- call g2l_state for ensemble member', member
           if (member==1) &
                WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_letkf_update --    Note: if ens_l is incorrect check user-defined indices in g2l_state!'
        END IF

        CALL U_g2l_state(step, domain_p, dim_p, ens_p(:, member), dim_l, &
             ens_l(:, member))

        IF (debug>0) &
             WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  ens_l', ens_l(:,member)
     END DO

     ! Store member index to make it accessible with PDAF_get_memberid
     member_save = 0

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_letkf_update -- call g2l_state for ensemble mean'

     CALL U_g2l_state(step, domain_p, dim_p, state_p, dim_l, &
          state_l)

     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  meanens_l', state_l

     CALL PDAF_timeit(15, 'old')


     ! LNETF analysis
     IF (dim_obs_l == 0) THEN
        ! UNOBSERVED DOMAIN

        CALL PDAF_timeit(51, 'new')
        CALL PDAF_timeit(7, 'new')

        ! Depending on type_forget, inflation on unobserved domain has to be inverted or applied here
        IF (type_forget==1) THEN
           ! prior inflation NOT on unobserved domains - take it back!
           invforget = 1.0/forget
           CALL PDAF_inflate_ens(dim_l, dim_ens, state_l, ens_l, invforget)
        ELSEIF (type_forget==2) THEN 
           ! analysis inflation ALSO on unobserved domains - add it!
           invforget = forget
           CALL PDAF_inflate_ens(dim_l, dim_ens, state_l, ens_l, invforget)
        ENDIF

        CALL PDAF_timeit(7, 'old')
        CALL PDAF_timeit(51, 'old')

     ELSE     
        ! OBSERVED DOMAIN

        CALL PDAF_timeit(7, 'new')

        ! only necessary if there are observations
        CALL PDAF_lnetf_analysis(domain_p, step, dim_l, dim_obs_f, dim_obs_l, &
             dim_ens, ens_l, HX_f, HXbar_f, rndmat, U_g2l_obs, &
             U_init_obs_l, U_likelihood_l, screen, type_forget, forget, &
             type_winf, limit_winf, cnt_small_svals, n_eff(domain_p), TA_l, flag)

        CALL PDAF_timeit(7, 'old')
        CALL PDAF_timeit(17, 'new')

        IF (dim_lag>0) THEN
           ! Compute transform matrix for smoother
           CALL PDAF_lnetf_smootherT(domain_p, step, dim_obs_f, dim_obs_l, &
                dim_ens, HX_noinfl_f, rndmat, U_g2l_obs, U_init_obs_l, U_likelihood_l, &
                screen, TA_noinfl_l, flag)
        END IF

        CALL PDAF_timeit(17, 'old')

     ENDIF 

     CALL PDAF_timeit(16, 'new')

     ! re-initialize full state ensemble on PE and mean state from local domain
     DO member = 1, dim_ens
        member_save = member

        IF (debug>0) then
           WRITE (*,*) '++ PDAF-debug: ', debug, &
                'PDAF_letkf_update -- call l2g_state for ensemble member', member
           WRITE (*,*) '++ PDAF-debug PDAF_letkf_update:', debug, '  ens_l', ens_l(:,member)
        END IF

        CALL U_l2g_state(step, domain_p, dim_l, ens_l(:, member), dim_p, ens_p(:,member))
     END DO
    
     CALL PDAF_timeit(16, 'old')
     CALL PDAF_timeit(51, 'new')
     CALL PDAF_timeit(18, 'new')

     ! *** Perform smoothing of past ensembles ***
     IF (dim_lag>0) THEN
        CALL PDAF_smoother_lnetf(domain_p, step, dim_p, dim_l, dim_ens, &
             dim_lag, TA_noinfl_l, ens_l, sens_p, cnt_maxlag, &
             U_g2l_state, U_l2g_state, screen)
     END IF

     CALL PDAF_timeit(18, 'old')

     ! clean up
     DEALLOCATE(ens_l, state_l)
     CALL PDAF_timeit(51, 'old')

     ! Set allocflag
     allocflag_l = 1

  END DO localanalysis

  IF (debug>0 .and. n_domains_p>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- End of local analysis loop'

  ! Set flag that we are not in the local analysis loop
  inloop = .false.

  ! Clean up arrays allocated in parallel
  DEALLOCATE(TA_l)
  DEALLOCATE(TA_noinfl_l)

!$OMP END PARALLEL

  CALL PDAF_timeit(6, 'old')


  ! *****************************************
  ! *** Perturb particles by adding noise ***
  ! *****************************************

  CALL PDAF_timeit(23, 'new')

  IF (noise_type>0) THEN
     IF (debug>0) &
          WRITE (*,*) '++ PDAF-debug: ', debug, &
          'PDAF_letkf_update -- add noise to particles'

     CALL PDAF_pf_add_noise(dim_p, dim_ens, state_p, ens_p, noise_type, noise_amp, screen)
  END IF

  CALL PDAF_timeit(23, 'old')

  ! Initialize mask array for effective ensemble size
  MASK = (n_eff > 0.0)

  ! *** Print statistics for local analysis to the screen ***
  CALL PDAF_print_local_obsstats(screen)

  IF (npes_filter>1) THEN
     ! Min/max effective sample sizes
     max_n_eff_l = MAXVAL(n_eff)
     CALL MPI_Reduce(max_n_eff_l, max_n_eff, 1, MPI_REALTYPE, MPI_MAX, &
          0, COMM_filter, MPIerr)
     min_n_eff_l = MINVAL(n_eff, MASK)
     CALL MPI_Reduce(min_n_eff_l, min_n_eff, 1, MPI_REALTYPE, MPI_MIN, &
          0, COMM_filter, MPIerr)

     ! Average effective sample size
     avg_n_eff_l = SUM(n_eff)/n_domains_p
     CALL MPI_Reduce(avg_n_eff_l, avg_n_eff, 1, MPI_REALTYPE, MPI_SUM, &
          0, COMM_filter, MPIerr)
     avg_n_eff = avg_n_eff / REAL(npes_filter)
  ELSE
     ! Min/max effective ensemble sizes
     max_n_eff = MAXVAL(n_eff)
     min_n_eff = MINVAL(n_eff, MASK)

     ! Average effective sample sizes
     avg_n_eff = SUM(n_eff)/n_domains_p
  END IF
 
  CALL PDAF_timeit(3, 'old')
 
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 8x, a, 9x, f7.1)') &
         'PDAF', 'Minimal  effective ensemble size:', min_n_eff
     WRITE (*, '(a, 8x, a, 9x, f7.1)') &
          'PDAF', 'Maximal  effective ensemble size:', max_n_eff
     WRITE (*, '(a, 8x, a, 9x, f7.1)') &
          'PDAF', 'Average  effective ensemble size:', avg_n_eff
     WRITE (*, '(a, 8x, a, 11x, i6)') &
          'PDAF', 'Number of small singular values:', &
          cnt_small_svals 

     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- analysis/re-init duration:', PDAF_time_temp(3), 's'
     END IF
  END IF
 
! *** Clean up from local analysis update ***
  DEALLOCATE(HX_f,  HXbar_f, rndmat)
  IF (dim_lag>0) DEALLOCATE(HX_noinfl_f)
#else
  WRITE (*,'(/5x,a/)') &
       '!!! WARNING: ANALYSIS STEP IS DEACTIVATED BY PDAF_NO_UPDATE !!!' 
#endif
 

! *** Poststep for analysis ensemble ***
  CALL PDAF_timeit(5, 'new')
  IF (mype == 0 .AND. screen > 0) THEN
     WRITE (*, '(a, 5x, a)') 'PDAF', 'Call pre-post routine after analysis step'
  ENDIF
  CALL U_prepoststep(step, dim_p, dim_ens, dim_ens_l, dim_obs_f, &
       state_p, Uinv, ens_p, flag)
  CALL PDAF_timeit(5, 'old')

  IF (mype == 0 .AND. screen > 0) THEN
     IF (screen > 1) THEN
        WRITE (*, '(a, 5x, a, F10.3, 1x, a)') &
             'PDAF', '--- duration of poststep:', PDAF_time_temp(5), 's'
     END IF
     WRITE (*, '(a, 55a)') 'PDAF Forecast ', ('-', i = 1, 55)
  END IF


! ********************
! *** Finishing up ***
! ********************

  IF (allocflag == 0) allocflag = 1

#ifndef PDAF_NO_UPDATE
  DEALLOCATE(n_eff)
#endif

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_letkf_update -- END'

END SUBROUTINE PDAF_lnetf_update

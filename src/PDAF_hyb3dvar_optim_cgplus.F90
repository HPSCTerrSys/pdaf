! Copyright (c) 2004-2024 Lars Nerger
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
! !ROUTINE: PDAF_hyb3dvar_optim_cgplus --- Optimization with CG+ for Hyb3dVar
!
! !INTERFACE:
SUBROUTINE PDAF_hyb3dvar_optim_cgplus(step, dim_p, dim_ens, dim_cv_par_p, dim_cv_ens_p, &
     dim_obs_p, ens_p, obs_p, dy_p, v_par_p, v_ens_p, U_prodRinvA, &
     U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, U_obs_op_lin, U_obs_op_adj, &
     opt_parallel, beta_3dvar, screen)

! !DESCRIPTION:
! Optimization routine for ensemble 3D-Var using the CG+ solver
!
! Variant for domain decomposed states.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2021-03 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE PDAF_timer, &
       ONLY: PDAF_timeit
  USE PDAF_memcounting, &
       ONLY: PDAF_memcount
  USE PDAF_mod_filtermpi, &
       ONLY: mype, comm_filter, npes_filter
  USE PDAF_mod_filter, &
       ONLY: method_cgplus_var, irest_cgplus_var, eps_cgplus_var, debug

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                  ! Current time step
  INTEGER, INTENT(in) :: dim_p                 ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens               ! ensemble size
  INTEGER, INTENT(in) :: dim_cv_par_p          ! Size of control vector (parameterized)
  INTEGER, INTENT(in) :: dim_cv_ens_p          ! Size of control vector (ensemble)
  INTEGER, INTENT(in) :: dim_obs_p             ! PE-local dimension of observation vector
  REAL, INTENT(in) :: ens_p(dim_p, dim_ens)    ! PE-local state ensemble
  REAL, INTENT(in)  :: obs_p(dim_obs_p)        ! Vector of observations
  REAL, INTENT(in)  :: dy_p(dim_obs_p)         ! Background innovation
  REAL, INTENT(inout) :: v_par_p(dim_cv_par_p) ! Control vector (parameterized part)
  REAL, INTENT(inout) :: v_ens_p(dim_cv_ens_p) ! Control vector (ensemble part)
  INTEGER, INTENT(in) :: screen                ! Verbosity flag
  INTEGER, INTENT(in) :: opt_parallel          ! Whether to use a decomposed control vector
  REAL, INTENT(in) :: beta_3dvar               ! Hybrid weight

! ! External subroutines 
! ! (PDAF-internal names, real names are defined in the call to PDAF)
  EXTERNAL :: U_prodRinvA, &   ! Provide product R^-1 A
       U_cvt, &                ! Apply control vector transform matrix to control vector (parameterized)
       U_cvt_adj, &            ! Apply adjoint control vector transform matrix (parameterized)
       U_cvt_ens, &            ! Apply control vector transform matrix to control vector (ensemble)
       U_cvt_adj_ens, &        ! Apply adjoint control vector transform matrix (ensemble)
       U_obs_op_lin, &         ! Linearized observation operator
       U_obs_op_adj            ! Adjoint observation operator

! !CALLING SEQUENCE:
! Called by: PDAF_3dvar_analysis_cvt
! Calls: PDAF_timeit
! Calls: PDAF_memcount
!EOP

! *** local variables ***
  INTEGER, SAVE :: allocflag = 0       ! Flag whether first time allocation is done
  INTEGER :: optiter                   ! Additional iteration counter
  INTEGER :: dim_cv_p                  ! Full size of control vector
  REAL :: J_tot                        ! Cost function
  REAL, ALLOCATABLE :: gradJ_p(:)      ! PE-local part of gradient of J
  REAL, ALLOCATABLE :: v_p(:)          ! PE-local full control vector

  ! Variables for CG+
  INTEGER :: method=2
  INTEGER :: irest=5
  REAL :: eps=1.0e-5
  INTEGER :: iprint(2), iflag, icall, mp, lp, i
  REAL, ALLOCATABLE :: d(:), gradJ_old_p(:), w(:)
  REAL :: tlev
  LOGICAL :: finish, update_J
  INTEGER :: iter, nfun
  COMMON /cgdd/    mp,lp
  COMMON /runinf/  iter,nfun


! **********************
! *** INITIALIZATION ***
! **********************

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_optim_CGPLUS -- START'

  ! Initialize overall dimension of control vector
  dim_cv_p = dim_cv_par_p + dim_cv_ens_p

  ! Settings for CG+
  method = method_cgplus_var  ! (1) Fletcher-Reeves, (2) Polak-Ribiere, (3) positive Polak-Ribiere (default=2)
  irest = irest_cgplus_var    ! (0) no restarts; (1) restart every n steps (default=1)
  EPS = eps_cgplus_var        ! Convergence constant (default=1.0e-5)
  icall = 0
  iflag = 0
  FINISH = .FALSE.
  update_J = .TRUE.
  optiter = 1

  ! Set verbosity of solver
  IF (screen>0 .AND. screen<2) THEN
     iprint(1) = -1
  ELSEIF (screen==2) THEN
     iprint(1) = 0
     IF (mype>0) iprint(1) = -1
  ELSE
     iprint(1) = 0
  END IF
  iprint(2) = 0  
  if (debug>0) iprint(1) = 0

  ! Allocate arrays
  ALLOCATE(v_p(dim_cv_p))
  ALLOCATE(d(dim_cv_p), w(dim_cv_p))
  ALLOCATE(gradJ_p(dim_cv_p), gradJ_old_p(dim_cv_p))
  IF (allocflag == 0) CALL PDAF_memcount(3, 'r', 4*dim_cv_p)

  IF (debug>0) THEN
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_CGPLUS', debug, &
          'Solver config: method  ', method
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_CGPLUS', debug, &
          'Solver config: restarts', irest
     WRITE (*,*) '++ PDAF-debug PDAF_hyb3dvar_optim_CGPLUS', debug, &
          'Solver config: EPS     ', EPS
  END IF

  ! Initialize numbers
  v_p = 0.0


! ***************************
! ***   Iterative solving ***
! ***************************

  IF (mype==0 .AND. screen > 0) &
       WRITE (*, '(a, 5x, a)') 'PDAF', '--- OPTIMIZE' 

  minloop: DO

! ********************************
! ***   Evaluate cost function ***
! ********************************

     IF (update_J) THEN
        CALL PDAF_timeit(53, 'new')
        CALL PDAF_hyb3dvar_costf_cvt(step, optiter, dim_p, dim_ens, &
             dim_cv_p, dim_cv_par_p, dim_cv_ens_p, dim_obs_p, ens_p, obs_p, &
             dy_p, v_par_p, v_ens_p, v_p, J_tot, gradJ_p, &
             U_prodRinvA, U_cvt, U_cvt_adj, U_cvt_ens, U_cvt_adj_ens, &
             U_obs_op_lin, U_obs_op_adj, opt_parallel, beta_3dvar)
        CALL PDAF_timeit(53, 'old')
     END IF


! ***************************
! ***   Optimize with CG+ ***
! ***************************

     CALL PDAF_timeit(54, 'new')
     IF (opt_parallel==0) THEN
        CALL CGFAM(dim_cv_p, v_p, J_tot, gradJ_p, d, gradJ_old_p, IPRINT, EPS, W,  &
             iflag, IREST, METHOD, FINISH)
     ELSE
        CALL CGFAM_mpi(dim_cv_p, v_p, J_tot, gradJ_p, d, gradJ_old_p, IPRINT, EPS, W,  &
             iflag, IREST, METHOD, FINISH, comm_filter, npes_filter)
     END IF
     CALL PDAF_timeit(54, 'old')


     ! *** Check exit status ***

     ! iflag=
     !    0 : successful termination
     !    1 : return to evaluate F and G
     !    2 : return with a new iterate, try termination test
     !   -i : error

     update_J = .TRUE.
     IF (iflag <= 0 .OR. icall > 10000) EXIT minloop
     IF (iflag == 1 ) icall = icall + 1
     IF (iflag == 2) THEN

        ! Termination Test.
        tlev = eps*(1.0 + ABS(J_tot))
        i=0
        checktest: DO
           i = i + 1
           IF(i > dim_cv_p) THEN
              FINISH = .TRUE.
              update_J = .FALSE.
              EXIT checktest
           ENDIF
           IF(ABS(gradJ_p(i)) > tlev) THEN
              update_J = .FALSE.
              EXIT checktest
           ENDIF
        END DO checktest

     ENDIF

     ! Increment loop counter
     optiter = optiter+1

  END DO minloop


! ********************
! *** Finishing up ***
! ********************

  ! Initialize two partial control vectors
  v_par_p = v_p(1 : dim_cv_par_p)
  v_ens_p = v_p(dim_cv_par_p+1 : dim_cv_p)

  DEALLOCATE(gradJ_p, v_p)
  DEALLOCATE(d, gradJ_old_p, w)

  IF (allocflag == 0) allocflag = 1

  IF (debug>0) &
       WRITE (*,*) '++ PDAF-debug: ', debug, 'PDAF_hyb3dvar_optim_CGPLUS -- END'

END SUBROUTINE PDAF_hyb3dvar_optim_cgplus

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
! !ROUTINE: PDAF_lknetf_init --- PDAF-internal initialization of localized KNETF
!
! !INTERFACE:
SUBROUTINE PDAF_lknetf_init(subtype, param_int, dim_pint, param_real, dim_preal, &
     ensemblefilter, fixedbasis, verbose, outflag)

! !DESCRIPTION:
! Initialization of LKNETF within PDAF. Performed are:\\
!   - initialize filter-specific parameters\\
!   - print screen information on filter configuration.
!
! !  This is a core routine of PDAF and
!    should not be changed by the user   !
!
! !REVISION HISTORY:
! 2017-08 - Lars Nerger - Initial code based on code for LNETF
! Later revisions - see svn log
!
! !USES:
  USE PDAF_mod_filter, &
       ONLY: incremental, dim_ens, forget, type_forget, type_trans, &
       dim_lag, localfilter, type_hyb, hyb_g, hyb_k

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(inout) :: subtype             ! Sub-type of filter
  INTEGER, INTENT(in) :: dim_pint               ! Number of integer parameters
  INTEGER, INTENT(inout) :: param_int(dim_pint) ! Integer parameter array
  INTEGER, INTENT(in) :: dim_preal              ! Number of real parameters 
  REAL, INTENT(inout) :: param_real(dim_preal)  ! Real parameter array
  LOGICAL, INTENT(out) :: ensemblefilter        ! Is the chosen filter ensemble-based?
  LOGICAL, INTENT(out) :: fixedbasis            ! Does the filter run with fixed error-space basis?
  INTEGER, INTENT(in) :: verbose                ! Control screen output
  INTEGER, INTENT(inout):: outflag              ! Status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init_filters
!EOP

! *** local variables ***
  REAL :: param_real_dummy    ! Dummy variable to avoid compiler warning


! ****************************
! *** INITIALIZE VARIABLES ***
! ****************************

  ! Initialize variable to prevent compiler warning
  param_real_dummy = param_real(1)

  ! Size of lag considered for smoother
  IF (dim_pint>=3) THEN
     IF (param_int(3) > 0) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-WARNING: LKNETKF does not yet support smoothing!'
        outflag = 11
     END IF
     dim_lag = 0
  END IF

  ! Whether incremental updating is performed
  ! We do not have incremental updating for LETKF!
  if (dim_pint>=4) THEN
     incremental = param_int(4)
     IF (param_int(4) /= 0 .AND. param_int(4) /= 1) THEN
        WRITE (*,'(/5x, a/)') &
             'PDAF-ERROR(10): LETKF does not yet support incremental updating!'
        outflag = 10
     END IF
  END IF

  ! Store type for forgetting factor
  IF (dim_pint >= 5) THEN
     type_forget = param_int(5)
  END IF

  ! Type of ensemble transformation
  IF (dim_pint >= 6) THEN     
     type_trans = param_int(6)
  END IF

  ! Type of hybrid weight
  type_hyb = 1
  IF (dim_pint >= 7) THEN     
     type_hyb = param_int(7)
  END IF

  ! Store prescribed hybrid weight
  hyb_g = 1.0
  IF (dim_preal >= 2) THEN
     hyb_g = param_real(2)
  END IF

  ! Store scale factor for hybrid weight (type_hyb=3 and 4)
  hyb_k = dim_ens
  IF (dim_preal >= 3) THEN
     hyb_k = param_real(3)
  END IF

  ! Define whether filter is mode-based or ensemble-based
  ensemblefilter = .TRUE.

  ! Define whether filter is domain localized
  localfilter = 1

  ! Initialize flag for fixed-basis filters
  fixedbasis = .FALSE.


! *********************
! *** Screen output ***
! *********************

  filter_pe2: IF (verbose == 1) THEN
  
     WRITE(*, '(/a, 4x, a)') 'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++  Local Hybrid Kalman-Nonlinear Ensemble Transform Filter  +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++                                                           +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++                Domain-localized LKNETF by                 +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++ L. Nerger, QJRMS, 148 (2022) 620-640, doi:10.1002/qj.4221 +++'
     WRITE(*, '(a, 4x, a)')  'PDAF', '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

     ! *** General output ***
     WRITE (*, '(/a, 4x, a)') 'PDAF', 'LKNETF configuration'
     WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'filter sub-type = ', subtype
     IF (subtype == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HNK) 2-step LKNETF: NETF before LETKF'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HKN) 2-step LKNETF: LETKF before NETF'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> (HSync) LKNETF synchronous'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> offline mode - 2-step LKNETF: NETF before LETKF'

        ! Reset subtype
        subtype = 0
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(2): No valid sub type!'
        outflag = 3
     END IF
     IF (type_trans == 0) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Transform ensemble including product with random matrix'
     ELSE IF (type_trans == 1) THEN
        WRITE (*, '(a, 12x, a)') 'PDAF', '--> Deterministic symmetric ensemble transformation'
     ELSE
        WRITE (*,'(/1x, a/)') &
             'PDAF-ERROR(9): Invalid setting for ensemble transformation!'
        outflag = 9
     END IF
     IF (incremental == 1) &
          WRITE (*, '(a, 12x, a)') 'PDAF', '--> Perform incremental updating'
     IF (type_forget == 0) THEN
        WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> prior inflation, forgetting factor:', forget
     ELSEIF (type_forget == 1) THEN
        WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> prior inflation on observed domains, forgetting factor: ', forget
     ELSEIF (type_forget == 2) THEN
        WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> posterior inflation, forgetting factor:', forget
     ELSEIF (type_forget == 3) THEN
        WRITE (*, '(a, 12x, a, f5.2)') 'PDAF', '--> posterior inflation on observed domains, forgetting factor: ', forget
     ELSE
        WRITE (*, '(/5x, a/)') 'PDAF-ERROR(8): Invalid type of forgetting factor!'
        outflag = 8
     ENDIF
     WRITE (*, '(a, 10x, a, i1)') 'PDAF', 'hybridization type = ', type_hyb
     WRITE (*, '(a, 12x, a, es10.2)') 'PDAF','--> hybrid weight input (gamma):', hyb_g
     IF (type_hyb==3 .OR. type_hyb==4) &
          WRITE (*, '(a, 12x, a, es10.2)') 'PDAF','--> hybrid norm (kappa):', hyb_k
     IF (type_hyb == 0) THEN
        WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_fix: fixed hybrid weight', hyb_g
     ELSEIF (type_hyb == 1) THEN
        WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_lin: (1 - N_eff/N_e)*', hyb_g
     ELSEIF (type_hyb == 2) THEN
        WRITE(*, '(a, 12x, a, f8.3)') 'PDAF', '--> use gamma_alpha: hybrid weight from N_eff/N>=', hyb_g
     ELSEIF (type_hyb == 3) THEN
        WRITE(*, '(a, 12x, a, f8.3, a, f8.3)') &
             'PDAF', '--> use gamma_ska: 1 - min(s,k)/sqrt(', hyb_k, ') with N_eff/N>=', hyb_g
     ELSEIF (type_hyb == 4) THEN
        WRITE(*, '(a, 12x, a, f8.3, a, f8.3)') &
             'PDAF', '--> use gamma_sklin: 1 - min(s,k)/sqrt(', hyb_k, ') >= 1-N_eff/N>=', hyb_g
     END IF
     WRITE (*, '(a, 12x, a, i5)') 'PDAF', '--> ensemble size N:', dim_ens

  END IF filter_pe2

END SUBROUTINE PDAF_lknetf_init

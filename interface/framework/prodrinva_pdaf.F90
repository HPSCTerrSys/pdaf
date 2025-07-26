!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule(Forschungszentrum Juelich GmbH)
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
!prodrinva_pdaf.F90: TSMP-PDAF implementation of routine
!                    'prodrinva_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: prodrinva_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: prodRinvA_pdaf --- Compute product of inverse of R with some matrix
!
! !INTERFACE:
SUBROUTINE prodRinvA_pdaf(step, dim_obs_p, rank_dim_ens, obs_p, A_p, C_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to compute the product of the inverse of 
! the observation error covariance matrix with
! the observed ensemble perturbations (SEIK/ETKF/ESTKF).
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
   USE mod_assimilation, &
        ONLY: rms_obs
   USE mod_assimilation, &
        ONLY: dim_obs, obscov_inv
  
   use mod_read_obs, only: multierr,clm_obserr, pressure_obserr
   use mod_read_obs, only: clm_obscov, vec_useObs, vec_useObs_global

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                ! Current time step
  INTEGER, INTENT(in) :: dim_obs_p           ! PE-local dimension of obs. vector
  INTEGER, INTENT(in) :: rank_dim_ens        ! Ensemble size in case of ETKF filter else 
                                             ! rank of initial covariance matrix for ESTKF and SEIK filter
  REAL, INTENT(in)    :: obs_p(dim_obs_p)    ! PE-local vector of observations
  ! The second dimension of input and ouput matrix is dim_ens for ETKF, 
  ! while it is rank for the ESTKF and SEIK filter
  REAL, INTENT(in)    :: A_p(dim_obs_p,rank_dim_ens) ! Input matrix from SEIK_ANALYSIS
  REAL, INTENT(out)   :: C_p(dim_obs_p,rank_dim_ens) ! Output matrix

! !CALLING SEQUENCE:
! Called by: PDAF_seik_analysis        (as U_prodRinvA)
! Called by: PDAF_seik_analysis_newT   (as U_prodRinvA)
! Called by: PDAF_etkf_analysis        (as U_prodRinvA)
! Called by: PDAF_estkf_analysis       (as U_prodRinvA)
!EOP

! *** local variables ***
  INTEGER :: i, j       ! index of observation component
  REAL :: ivariance_obs ! inverse of variance of the observations
  REAL :: clm_obserr_model(dim_obs_p) ! errors of observations in the model domain
  REAL :: obscov_inv_l(dim_obs_p,dim_obs_p) ! local inverse covariance matrix
  integer :: countR, countC
  logical :: vec_useObs_p(dim_obs)

! **********************
! *** INITIALIZATION ***
! **********************

  WRITE (*,*) 'TEMPLATE prodrinva_pdaf.F90: Implement multiplication here!'

  SELECT CASE (multierr)
  CASE (0)
  ! *** initialize numbers
  ivariance_obs = 1.0 / rms_obs ** 2


! *************************************
! ***                -1             ***
! ***           C = R   A           ***
! ***                               ***
! *** The inverse observation error ***
! *** covariance matrix is not      ***
! *** computed explicitely.         ***
! *************************************

#ifdef PDAF_DEBUG
  print *, 'prodrinva_pdaf: covariance matrix: case 0'
#endif

  DO j = 1, rank_dim_ens !rank
     DO i = 1, dim_obs_p
        C_p(i, j) = ivariance_obs * A_p(i, j)
     END DO
  END DO

  CASE (1)   ! Diagonal covariance matrix, read from observation file
    print *, 'prodrinva_pdaf: Observation covariance matrix: case 1'    
    clm_obserr_model = pack(clm_obserr,vec_useObs)
    DO j = 1, rank_dim_ens
      DO i = 1, dim_obs_p
        C_p(i, j) = 1.0/clm_obserr_model(i) * A_p(i, j)
      END DO
    END DO

  CASE (2)	! Fully occupied covariance matrix, read from observation file
    print *, 'prodrinva_pdaf: Observation covariance matrix: case 2' 
    ! Find out observations that shoud be considered, obscov_inv contains all used observations
    j = 1
    DO i = 1, size(vec_useObs_global,1)
      if(vec_useObs_global(i).eqv..true.) then
        if(vec_useObs(i).eqv..true.) then 
          vec_useObs_p(j) = .true.
        else
          vec_useObs_p(j) = .false.
        end if
        j = j+1
      end if
    END DO

    ! Local covariance matrix
    countR = 1
    countC = 1
    do i = 1, size(obscov_inv,1)
      if(vec_useObs_p(i)) then
        do j = 1, size(obscov_inv,2)
          if(vec_useObs_p(j)) then
            obscov_inv_l(countR,countC) = obscov_inv(i,j)
            countC = countC+1;
          end if
        end do
        countC = 1
        countR = countR+1;
      end if
    end do

    C_p = matmul(obscov_inv_l,A_p)

    IF(allocated(obscov_inv))DEALLOCATE(obscov_inv)
  END SELECT
  if(allocated(clm_obserr))deallocate(clm_obserr)
  if(allocated(clm_obscov))deallocate(clm_obscov)

END SUBROUTINE prodRinvA_pdaf

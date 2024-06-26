! Copyright (c) 2012-2024 Lars Nerger, lars.nerger@awi.de
!
! This routine is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License
! as published by the Free Software Foundation, either version
! 3 of the License, or (at your option) any later version.
!
! This code is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this software.  If not, see <http://www.gnu.org/licenses/>.
!
!$Id$
!BOP
!
! !ROUTINE: PDAF_diag_EnsStats --- Compute ensemble statistics
!
! !INTERFACE:
SUBROUTINE PDAF_diag_ensstats(dim, dim_ens, element, &
     state, ens, skewness, kurtosis, status)

! !DESCRIPTION:
! This routine computes the higher-order ensemble statistics (skewness and 
! kurtosis). Inputs are the ensemble array and the state vector about which
! the statistics are computed (usually the ensemble mean). In addition, the
! index of the element has to be specified for which the statistics are
! computed. If this is 0, the mean statistics over all elements are computed.
!
! The definition used for kurtosis follows that used by Lawson and Hansen,
! Mon. Wea. Rev. 132 (2004) 1966.
!
! !REVISION HISTORY:
! 2012-09 - Lars Nerger - Initial code for SANGOMA based on PDAF
! 2013-11 - L. Nerger - Adaption to SANGOMA data model
! 2016-05 - Lars Nerger - Back-porting to PDAF
!
! !USES:
! Include definitions for real type of different precision
! (Defines BLAS/LAPACK routines and MPI_REALTYPE)
#include "typedefs.h"

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: dim               ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens           ! Ensemble size
  INTEGER, INTENT(in) :: element           ! ID of element to be used
       ! If element=0, mean values over all elements are computed
  REAL, INTENT(in)    :: state(dim)        ! State vector
  REAL, INTENT(in)    :: ens(dim, dim_ens) ! State ensemble
  REAL, INTENT(out)   :: skewness          ! Skewness of ensemble
  REAL, INTENT(out)   :: kurtosis          ! Kurtosis of ensemble
  INTEGER, INTENT(out) :: status           ! Status flag (0=success)
!EOP

! *** local variables ***
  INTEGER :: i, elem     ! Counters
  REAL :: m2, m3, m4     ! Statistical moments


! **************************
! *** Compute statistics ***
! **************************

  IF (element > 0 .AND. element <= dim) THEN

     m2 = 0.0
     m3 = 0.0
     m4 = 0.0
     do i=1,dim_ens
        m2 = m2 + (ens(element, i) - state(element))**2
        m3 = m3 + (ens(element, i) - state(element))**3
        m4 = m4 + (ens(element, i) - state(element))**4 
     end do
     m2 = m2 / real(dim_ens)
     m3 = m3 / real(dim_ens)
     m4 = m4 / real(dim_ens)

     skewness = m3 / sqrt(m2)**3
     kurtosis = m4 / m2**2 -3.0

     ! Set status flag for success
     status = 0

  ELSE IF (element == 0) THEN

     skewness = 0.0
     kurtosis = 0.0

     DO elem = 1, dim
        m2 = 0.0
        m3 = 0.0
        m4 = 0.0
        do i=1,dim_ens
           m2 = m2 + (ens(elem, i) - state(elem))**2
           m3 = m3 + (ens(elem, i) - state(elem))**3
           m4 = m4 + (ens(elem, i) - state(elem))**4 
        end do
        m2 = m2 / real(dim_ens)
        m3 = m3 / real(dim_ens)
        m4 = m4 / real(dim_ens)

        skewness = skewness + m3 / sqrt(m2)**3
        kurtosis = kurtosis + m4 / m2**2 -3.0

     END DO

     skewness = skewness / REAL(dim)
     kurtosis = kurtosis / REAL(dim)

     ! Set status flag for success
     status = 0

  ELSE

     ! Choice of 'element' not valid
     status = 1

  END IF


END SUBROUTINE PDAF_diag_ensstats

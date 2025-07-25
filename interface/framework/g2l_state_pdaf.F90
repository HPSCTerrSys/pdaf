!-------------------------------------------------------------------------------------------
!Copyright (c) 2013-2016 by Wolfgang Kurtz, Guowei He and Mukund Pondkule (Forschungszentrum Juelich GmbH)
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
!g2l_state_pdaf.F90: TSMP-PDAF implementation of routine
!                    'g2l_state_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: g2l_state_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: g2l_state_pdaf --- Restrict a model state to a local analysis domain
!
! !INTERFACE:
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during the loop over all
! local analysis domains in PDAF\_lseik\_update
! before the analysis on a single local analysis 
! domain.  It has to project the full PE-local 
! model state onto the current local analysis 
! domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_tsmp, ONLY: tag_model_parflow, &
       tag_model_clm, model
  USE mod_tsmp, &
       ONLY: nx_local, ny_local
#if defined CLMSA
#if defined CLMFIVE
  USE decompMod, ONLY: get_proc_bounds
#else  
  USE decompMod, ONLY: get_proc_bounds_atm
#endif
  use enkf_clm_mod, only: hactiveg_levels, num_layer, state_setup, num_hactiveg_patch, hactiveg_patch, clm_varsize_tws
  USE enkf_clm_mod, ONLY: clmupdate_tws
#endif

  USE iso_c_binding, ONLY: c_loc

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step
  INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          ! Local state dimension
  REAL, TARGET, INTENT(in)    :: state_p(dim_p) ! PE-local full state vector 
  REAL, TARGET, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis domain

  INTEGER :: i, n_domain, nshift_p
  INTEGER :: sub, j, g
  INTEGER :: begg, endg   ! per-proc gridcell ending gridcell indices
! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update    (as U_g2l_state)
! Called by: PDAF_letkf_update    (as U_g2l_state)
! Called by: PDAF_lestkf_update   (as U_g2l_state)
!EOP

! *************************************
! *** Initialize local state vector ***
! *************************************
#ifndef CLMSA
  if (model == tag_model_parflow) then
     n_domain = nx_local * ny_local
     DO i = 0, dim_l-1
        nshift_p = domain_p + i * n_domain
        state_l(i+1) = state_p(nshift_p)
     ENDDO
  else  if (model == tag_model_clm) then
     state_l(dim_l) = state_p(domain_p)
  end if
  !call g2l_state(domain_p, c_loc(state_p), dim_l, c_loc(state_l))
#else
  ! beg and end gridcell for atm
#if defined CLMFIVE
  call get_proc_bounds(begg, endg)
#else  
  call get_proc_bounds_atm(begg, endg)
#endif

  if (clmupdate_tws.eq.1) then

  
    ! first depth dependent variables --> liq and ice (together in statevector or not)
    ! dim_l is number of layers of gridcell + 3 (3 for other compartments that are added to the statevector)
    ! lets loop over known layers 
    
    ! first do canopy water to check how many variables are for this gridcell available next to snow and surface water as well as the layers 
    ! compartments to subtract from dim_l to get layers variables:
  
    if (clm_varsize_tws(5).ne.0) then
      sub=3
    else
      sub=2
    end if
  
  
    select case (state_setup)
    case(0) ! liq and ice seperated
      g = hactiveg_levels(domain_p,1)
      do i = 1, (dim_l-sub)/2 ! two entries for liq and ice seperated
        do j = 1, num_layer(i) ! i is the layer that we are in right now
          if (g==hactiveg_levels(j,i)) then ! if the counter is the gridcell of the local domain, we know the position in the statevector
  
            if (i == 1) then ! if first layer
              state_l(i) = state_p(j)     ! first liquid water as it is first in the statevector
              state_l(i+(dim_l-3)/2) = state_p(j+clm_varsize_tws(1))
            else
              state_l(i) = state_p(j + sum(num_layer(1:i-1)))
              state_l(i+(dim_l-3)/2) = state_p(j + sum(num_layer(1:i-1)) + clm_varsize_tws(1))
            end if
  
          end if
        end do
      end do
  
      do j = 1, num_layer(1)
        if (g==hactiveg_levels(j,1)) then
          if (sub==3) then
            state_l(dim_l-2) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2))
            state_l(dim_l-1) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3))
            state_l(dim_l) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3) + clm_varsize_tws(4))
          else
            state_l(dim_l-1) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2))
            state_l(dim_l) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3))
          end if
        end if
      end do
  
    case(1)
  
      g = hactiveg_levels(domain_p,1)
      do i = 1, dim_l-sub ! liq and ice added up
        do j = 1, num_layer(i) ! i is the layer that we are in right now
          if (g==hactiveg_levels(j,i)) then ! if the counter is the gridcell of the local domain, we know the position in the statevector
  
            if (i == 1) then ! if first layer
              state_l(i) = state_p(j)     ! first liquid water as it is first in the statevector
            else
              state_l(i) = state_p(j + sum(num_layer(1:i-1)))
            end if
  
          end if
        end do
      end do
  
      do j = 1, num_layer(1)
        if (g==hactiveg_levels(j,1)) then
          if (sub==3) then
            state_l(dim_l-2) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2))
            state_l(dim_l-1) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3))
            state_l(dim_l) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3) + clm_varsize_tws(4))
          else
            state_l(dim_l-1) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2))
            state_l(dim_l) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3))
          end if
        end if
      end do
    
    case(2) ! only tws in statevector
      g = hactiveg_levels(domain_p,1)
      do j = 1, num_layer(1)
        if (g==hactiveg_levels(j,1)) then
          state_l(1) = state_p(j)
        end if
      end do
  
    case(3)
  
      g = hactiveg_levels(domain_p,1)
  
      do j = 1, num_layer(1)
        if (g==hactiveg_levels(j,1)) then
          state_l(1) = state_p(j)
          state_l(2) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2))
        end if
      end do
  
    case(4)
  
      g = hactiveg_levels(domain_p,1)
  
      do j = 1, num_layer(1)
        if (g==hactiveg_levels(j,1)) then
          state_l(1) = state_p(j)
          state_l(dim_l) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2))
        end if
      end do
  
      if (dim_l==3) then
        do j = 1, num_layer(8)
          if (g==hactiveg_levels(j,8)) then
            state_l(2) = state_p(j + clm_varsize_tws(1))
          end if
        end do
      end if
  
    case(5)
  
      g = hactiveg_levels(domain_p,1)
  
      do j = 1, num_layer(1)
        if (g==hactiveg_levels(j,1)) then
          state_l(1) = state_p(j)
          state_l(dim_l) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3))
        end if
      end do
  
      if (dim_l>=3) then
        do j = 1, num_layer(4)
          if (g==hactiveg_levels(j,4)) then
            state_l(2) = state_p(j + clm_varsize_tws(1))
          end if
        end do
      end if
  
      if (dim_l>=4) then
        do j = 1, num_layer(13)
          if (g==hactiveg_levels(j,13)) then
            state_l(3) = state_p(j + clm_varsize_tws(1) + clm_varsize_tws(2))
          end if
        end do
      end if
  
    end select
  else  

  n_domain = endg - begg + 1
  DO i = 0, dim_l-1
     nshift_p = domain_p + i * n_domain
     state_l(i+1) = state_p(nshift_p)
  ENDDO

  end if
#endif
  
END SUBROUTINE g2l_state_pdaf

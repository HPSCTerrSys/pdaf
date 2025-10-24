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
!init_dim_l_pdaf.F90: TSMP-PDAF implementation of routine
!                     'init_dim_l_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_dim_l_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_dim_l_pdaf --- Set dimension of local model state
!
! !INTERFACE:
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: LSEIK/LETKF/LESTKF
!
! The routine is called during analysis step
! in the loop over all local analysis domain.
! It has to set the dimension of local model
! state on the current analysis domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_tsmp, ONLY: tag_model_parflow, &
       tag_model_clm, model
  USE mod_tsmp, &
       ONLY: init_dim_l_pfl
  use clm_varpar   , only : nlevsoi
#ifdef CLMSA
  USE enkf_clm_mod, &
       ONLY: hactiveg_levels, num_layer, state_setup, num_hactiveg_patch, hactiveg_patch, clm_varsize_tws
  USE enkf_clm_mod, ONLY: clmupdate_tws
  USE enkf_clm_mod, ONLY: init_dim_l_clm
#endif
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: step     ! Current time step
  INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    ! Local state dimension

  integer :: g, i, count

! !CALLING SEQUENCE:
! Called by: PDAF_lseik_update   (as U_init_dim_l)
! Called by: PDAF_lestkf_update  (as U_init_dim_l)
! Called by: PDAF_letkf_update   (as U_init_dim_l)
!EOP

! ****************************************
! *** Initialize local state dimension ***
! ****************************************

#if defined PARFLOW_STAND_ALONE
  ! Set the size of the local analysis domain
  call init_dim_l_pfl(dim_l)
#endif

#if defined COUP_OAS_PFL
  if (model == tag_model_parflow) then
     ! Set the size of the local analysis domain
     call init_dim_l_pfl(dim_l)
  end if
  if (model == tag_model_clm) then
     ! Set the size of the local analysis domain
     dim_l = 1
  end if
#endif

#if defined CLMSA
  ! Set the size of the local analysis domain
  ! for clm stand alone mode only

  if (clmupdate_tws/=1) then
  call init_dim_l_clm(domain_p, dim_l)
  end if

  ! Set the size of the local analysis domain
  ! for clm stand alone mode only
  if (clmupdate_tws==1) then

    ! initialize local state dimension --> different for each domain as a domain is a gridcell, number of layers in statevector differs per gridcell
    ! --> check with begg and endg as well as domain_p and hactiveg_levels what is going on
    ! hactiveg_levels in level 1 at position domain_p gives gridcell index of domain_p --> go through rest of hactiveg_levels and check for number of layers
    dim_l = 0
    g = hactiveg_levels(domain_p,1)

    select case(state_setup)
    case(0)
      do i = 1,nlevsoi
        do count = 1, num_layer(i)
          if (g==hactiveg_levels(count,i)) then
            dim_l = dim_l+2
          end if
        end do
      end do

      ! snow and surface water
      dim_l = dim_l+2

      ! canopy water
      if (clm_varsize_tws(5)/=0) then
        dim_l = dim_l+1
      end if

    case(1)
      do i = 1,nlevsoi
        do count = 1, num_layer(i)
          if (g==hactiveg_levels(count,i)) then
            dim_l = dim_l+1
          end if
        end do
      end do

      ! snow and surface water
      dim_l = dim_l+2

      if (clm_varsize_tws(5)/=0) then
        dim_l = dim_l+1
      end if

    case(2) ! only TWS in statevector

      dim_l=1

    case(3) ! sum over all soil layers and snow in statevector

      dim_l=2

    case(4)

      dim_l=2

      do count = 1, num_layer(8)
        if (g==hactiveg_levels(count,8)) then
          dim_l = dim_l+1
        end if
      end do

    case(5)

      dim_l=2

      do count = 1, num_layer(4)
        if (g==hactiveg_levels(count,4)) then
          dim_l = dim_l+1
        end if
      end do

      do count = 1, num_layer(13)
        if (g==hactiveg_levels(count,13)) then
          dim_l = dim_l+1
        end if
      end do

    end select

  end if
#endif

END SUBROUTINE init_dim_l_pdaf

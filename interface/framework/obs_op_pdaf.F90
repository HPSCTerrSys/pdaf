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
!obs_op_pdaf.F90: TSMP-PDAF implementation of routine
!                 'obs_op_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: obs_op_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: obs_op_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_pdaf(step, dim_p, dim_obs_p, state_p, m_state_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/EnKF/ETKF/ESTKF
!
! The routine is called during the analysis step.
! It has to perform the operation of the
! observation operator acting on a state vector.
! For domain decomposition, the action is on the
! PE-local sub-domain of the state and has to
! provide the observed sub-state for the PE-local
! domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
   USE mod_assimilation, &
        ONLY: obs_index_p, obs_p
   USE mod_parallel_pdaf, ONLY: mype_world
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
   USE mod_assimilation, ONLY: sc_p
#endif
#endif
   USE mod_assimilation, ONLY: obs_interp_indices_p
   USE mod_assimilation, ONLY: obs_interp_weights_p
   use mod_tsmp, &
       only: obs_interp_switch, &
       soilay, &
       soilay_fortran, &
       nz_glob, &
       da_crns_depth_tol, &
       crns_flag
!      tcycle

   USE, INTRINSIC :: iso_c_binding, ONLY: C_F_POINTER

#if defined CLMSA
   USE enkf_clm_mod, &
        ONLY : clm_varsize, clm_paramarr, clmupdate_swc, clmupdate_T, clmcrns_bd
   USE enkf_clm_mod, &
        ONLY : clmupdate_lai, clm_begp, clm_endp, clm_patch2gc, clm_patchwt, &
               clm_lai_patch_itype, clmt_printensemble, clmupdate_lai_params, &
               lai_mode3_npool
   USE pftconMod, only : pftcon
#ifdef CLMFIVE
   USE clm_instMod, &
     ONLY : soilstate_inst
#endif
#endif
  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step               ! Currrent time step
  INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
  REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state
  integer :: i, j, k, z, n
  integer :: icorner
  logical :: lpointobs       !If true: no special observation; use point observation
! !CALLING SEQUENCE:
! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT   (as U_obs_op)
! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
!EOP

! *** local variables ***

! hcp test with hardcoding variable declaration
real, dimension(:), allocatable :: soide !soil depth
!real(8), dimension(0:12), parameter :: &
! soide=(/0.d0,  0.02d0,  0.05d0,  0.1d0,  0.17d0, 0.3d0,  0.5d0, &
!                0.8d0,   1.3d0,   2.d0,  3.d0, 5.d0,  12.d0/) !soil depth

real :: tot, avesm, avelai, avesm_temp, Dp, slatop_val
integer :: nsc
character(len=48) :: fn_hx
! end of hcp

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
! Variables used in crns version 2
REAL :: weights_r1(920), weights_r2(920), weights_r3(920)
Real :: weights_layer(8)
Integer :: nweights(8)
Real  :: d86_r1, d86_r2, d86_r3
REAL, PARAMETER :: r1 = 1.0
REAL, PARAMETER :: r2 = 20.0
REAL, PARAMETER :: r3 = 85.0
REAL :: bd, y
REAL :: sum_r1, sum_r2, sum_r3, totw
#endif
#endif


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

! If no special observation operator is compiled, use point observations
lpointobs = .true.

! LAI assimilation assuming statevec contains leafc, slatop, dsladlai for each patch
! two component op to get to LAI obs
! 1 calculate tlai from statevec per patch
! state_p contains leafc of varsize then slatop of varsize and then dsadlai of varsize
! 2 average patches according to weight per gridcell
if (clmupdate_lai==2) then
  lpointobs = .false. ! disable general obs_op
  ! Note: This fails if there ever is an observation pointing to a gridcell with 0 patches
  ! this should not happen, but would cause div by zero.

  do i = 1, dim_obs_p ! loop over all observations
    avelai = 0.0 ! use avelai as grid cell average collecter
    do j = 1, clm_endp-clm_begp+1 ! loop over all patches
      write(*,*) 'DEBUG LAI : obs_index_p(i), clm_patch2gc(j)', obs_index_p(i), clm_patch2gc(j)
      if (obs_index_p(i)==clm_patch2gc(j)) then
        write(*,*) 'DEBUG LAI : state dsladlai', state_p(j+2*clm_varsize)
        if (state_p(j+2*clm_varsize)>0.0) then
          avelai = avelai + clm_patchwt(j) * ((state_p(j+1*clm_varsize)&
            *(exp(state_p(j)*state_p(j+2*clm_varsize)) - 1.0))/state_p(j+2*clm_varsize)) ! formula for tlai from leafc
        else
          avelai = avelai +  clm_patchwt(j) *(state_p(j+1*clm_varsize)*state_p(j)) ! 2nd formula
        endif ! dsladlai decider
      endif ! this patch is in gridcell of observation
      write(*,*) 'DEBUG LAI : average ', avelai
    enddo ! all patches for obs checked
    m_state_p(i) = avelai ! == lai of gridcell where the observation is
  enddo ! all observations
endif ! clmupdate_lai == 2 : m_state_p now contains lai for each gridcell with an observation.

! LAI assimilation assuming statevec contains leafc, livestemc, deadstemc
! for each vegetated patch (bare-ground PFT excluded).
! H uses slatop from state (update_lai_params=1) or pftcon; dsladlai from pftcon.
if (clmupdate_lai==3) then
  lpointobs = .false.

  do i = 1, dim_obs_p
    avelai = 0.0
    do j = 1, clm_varsize
      if (obs_index_p(i)==clm_patch2gc(j)) then
        if (clmupdate_lai_params == 1) then
          slatop_val = state_p(j + lai_mode3_npool*clm_varsize)
        else
          slatop_val = pftcon%slatop(clm_lai_patch_itype(j))
        end if
        if (pftcon%dsladlai(clm_lai_patch_itype(j)) > 0.0) then
          avelai = avelai + clm_patchwt(j) * ((slatop_val&
            *(exp(state_p(j)*pftcon%dsladlai(clm_lai_patch_itype(j))) - 1.0))&
            /pftcon%dsladlai(clm_lai_patch_itype(j)))
        else
          avelai = avelai + clm_patchwt(j) * (slatop_val*state_p(j))
        endif
      endif
    enddo
    m_state_p(i) = avelai
  enddo

#ifdef PDAF_DEBUG
  IF(clmt_printensemble == step .OR. clmt_printensemble == -1) THEN
    WRITE(fn_hx, "(a,i5.5,a,i5.5,a)") "lai_hx_", mype_world, "_", step, ".txt"
    OPEN(unit=72, file=fn_hx, action="write", status="replace")
    WRITE(72,"(a,i8)") "# mode3 H(x) at PDAF step ", step
    WRITE(72,"(a,i8)") "# dim_obs_p ", dim_obs_p
    DO i = 1, dim_obs_p
      WRITE(72,"(a,i8,a)") "# obs_index ", obs_index_p(i)
      WRITE(72,"(es22.15)") m_state_p(i)
    END DO
    CLOSE(72)
  END IF
#endif
endif ! clmupdate_lai == 3



#if defined CLMSA
if (clmupdate_T==1) then

  lpointobs = .false.

  DO i = 1, dim_obs_p
    ! Equation for LST computation from Kustas2009, Eq(7)
    ! http://dx.doi.org/10.1016/j.agrformet.2009.05.016
    !
    ! Comment: Fractional vegetation cover (Eq(8) from Kustas2009)
    ! currently implemented with simplified settings: Vegetation
    ! clumping parameter `Omega=1`; radiometer view angle `phi=0`
     m_state_p(i) &
    = (exp(-0.5*clm_paramarr(obs_index_p(i))) &
                     *state_p(obs_index_p(i))**4 &
       +(1.-exp(-0.5*clm_paramarr(obs_index_p(i)))) &
                     *state_p(clm_varsize+obs_index_p(i))**4)**0.25
  END DO
!  write(*,*) 'now is cycle ', tcycle
!  write(*,*) 'obs_p =', obs_p(:)
!  write(*,*) 'model LST', m_state_p(:)
!  write(*,*) 'TG', state_p(obs_index_p(:))
!  write(*,*) 'TV', state_p(clm_varsize+obs_index_p(:))
end if
#endif


#ifndef CLMSA
#ifndef OBS_ONLY_CLM
 if (crns_flag==1) then
    !Schroen et al HESS 2017 modelled CRNS averaging
    lpointobs = .false.
     call C_F_POINTER(soilay,soilay_fortran,[nz_glob])
     Allocate(soide(0:nz_glob))
     soide(0)=0.d0
     do i=1,nz_glob
       soide(i)=soide(i-1)+soilay_fortran(nz_glob-i+1)
     end do
     do i = 1, dim_obs_p

       !Initial average soil moisture for 1st iteration
       avesm=0.d0
       do j=1,nz_glob
            avesm=avesm+(soide(j)-soide(j-1))*state_p(sc_p(j,i))/soide(nz_glob)
       end do
       avesm_temp=0.d0

       !iteration
       do while (abs(avesm-avesm_temp)/avesm >= da_crns_depth_tol)
          !Averaging, conventional profile, Schroen et al HESS 2017 Eq. (3)
          avesm_temp=avesm
          Dp=0.058d0/(avesm+0.0829d0)

          !Sum weight*soil_moisture
          avesm=0.d0; nsc=nz_glob
          do j=1,nz_glob
             if ((soide(j-1)<Dp).AND.(Dp<=soide(j))) then
               nsc=j
             end if
          end do
          do j=1, nsc-1
              avesm=avesm+(1.d0-0.5d0*(soide(j)+soide(j-1))/Dp)*(soide(j)-soide(j-1)) &
                    *state_p(sc_p(j,i))/Dp
          end do
          avesm=avesm+(1.d0-0.5d0*(Dp+soide(nsc-1))/Dp)*(Dp-soide(nsc-1)) &
             *state_p(sc_p(nsc,i))/Dp

          !Sum weight
          tot=0.d0
          do j=1, nsc-1
              tot =   tot+(1.d0-0.5d0*(soide(j)+soide(j-1))/Dp)*(soide(j)-soide(j-1)) &
                                                    /Dp
          end do
          tot  =  tot+(1.d0-0.5d0*(Dp+soide(nsc-1))/Dp)*(Dp-soide(nsc-1)) &
                                               /Dp

          avesm=avesm/tot
       end do
       m_state_p(i)=avesm
     end do
     deallocate(soide)
 end if
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
#ifdef CLMFIVE
 if (crns_flag==2) then
   lpointobs = .false.
   ! CRNS implementation based on Schrön et al. 2017 using
   ! d86 for 3 different radius values
   DO i = 1, dim_obs_p
     ! Bulk density average for the 8 considered layers
     IF (clmcrns_bd > 0.0) THEN
       bd = clmcrns_bd
     ELSE
       bd = 0.0
       DO j = 1, 8
         bd = bd + soilstate_inst%bd_col(obs_index_p(i),j) ! bulk density
       END DO
         bd = bd / 8.0 * 0.001 ! average and convert from kg/m^3 to g/cm^3
     END IF
     ! CRNS observed value
     y = obs_p(i) ! CRNS observation
     ! Penetration depth calculations D86(bd, r, y)
     d86_r1 = (1/bd*(8.321+0.14249*(0.96655+exp(-0.01*r1))*(20+y)/(0.0429+y)))
     d86_r2 = (1/bd*(8.321+0.14249*(0.96655+exp(-0.01*r2))*(20+y)/(0.0429+y)))
     d86_r3 = (1/bd*(8.321+0.14249*(0.96655+exp(-0.01*r3))*(20+y)/(0.0429+y)))
     ! Then calculate the weights for thin (1mm) slices of the layers for 85cm
     sum_r1 = 0.0
     sum_r2 = 0.0
     sum_r3 = 0.0
     DO j = 1, 920 ! depth in mm but in calculation used in cm:
       weights_r1(j) = exp(-2*(j/10)/d86_r1)
       weights_r2(j) = exp(-2*(j/10)/d86_r2)
       weights_r3(j) = exp(-2*(j/10)/d86_r3)

       sum_r1 = sum_r1 + weights_r1(j)
       sum_r2 = sum_r2 + weights_r2(j)
       sum_r3 = sum_r3 + weights_r3(j)
     END DO
     ! Normalize the weights:
     weights_r1(:) = weights_r1(:) / sum_r1
     weights_r2(:) = weights_r2(:) / sum_r2
     weights_r3(:) = weights_r3(:) / sum_r3
     ! assign average weights to each layer
     nweights(:) = 0
     weights_layer(:) = 0.0
     ! z index for different layers, manually here, could be done better with model layer depth
     DO j = 1, 920
       IF (j > 680) then
         z = 8
       ELSE IF (j < 680 .and. j > 480) then
         z = 7
       ELSE IF (j < 480 .and. j > 320) then
         z = 6
       ELSE IF (j < 320 .and. j > 200) then
         z = 5
       ELSE IF (j < 200 .and. j > 120) then
         z = 4
       ELSE IF (j < 120 .and. j > 60) then
         z = 3
       ELSE IF (j < 60 .and. j > 20) then
         z = 2
       ELSE IF (j < 20) then
         z = 1
       END IF
       weights_layer(z) = weights_layer(z) + weights_r1(j) + weights_r2(j) + weights_r3(j)
       nweights(z) = nweights(z) + 1
     END DO
     ! Normalize the weights
     totw = 0.0
     DO j = 1, 8
        weights_layer(j) = weights_layer(j) / (nweights(j) * 3)
        totw = totw + weights_layer(j)
     END DO
     weights_layer(:) = weights_layer(:) / totw
     ! Finally use the weights to calculate the weighted average of the state variable
     avesm = 0.0
     DO j = 1, 8
       avesm = avesm + weights_layer(j) * state_p(obs_index_p(i) + (j-1))
       ! This assumes that obs_index_p(i) for obs i is the index of
       ! the first layer of the gridcell where obs i is
     END DO
     ! Assign new average as the state variable
     m_state_p(i) = avesm
     ! end loop over observations
   END DO
 end if
#endif
#endif
#endif


 if(obs_interp_switch == 1) then

      lpointobs = .false.

      do i = 1, dim_obs_p

          m_state_p(i) = 0
          do icorner = 1, 4
              m_state_p(i) = m_state_p(i) + state_p(obs_interp_indices_p(i,icorner)) * obs_interp_weights_p(i,icorner)
          end do

      end do

  end if

  if(lpointobs) then

  DO i = 1, dim_obs_p
     m_state_p(i) = state_p(obs_index_p(i))
  END DO

  end if

END SUBROUTINE obs_op_pdaf

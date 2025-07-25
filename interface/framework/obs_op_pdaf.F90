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
      ONLY: obs_index_p, &
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
      depth_obs_p, & !obs_p
      sc_p, &
#endif
#endif
      obs_interp_indices_p, &
      obs_interp_weights_p, &
      obs_id_p, &
      screen, &
      tws_temp_mean_d 
    USE mod_read_obs, ONLY: crns_flag, vec_numPoints_global !clm_obs
    use mod_tsmp, &
      only: obs_interp_switch, &
      soilay, &
      soilay_fortran, &
      nz_glob
  !      tcycle
  
    USE, INTRINSIC :: iso_c_binding
  
#if defined CLMSA
    USE enkf_clm_mod, & 
        ONLY : clm_varsize, clm_paramarr, clmupdate_swc, clmupdate_T, clmupdate_tws, clm_varsize_tws, state_setup, &
        num_layer, hactiveg_levels, num_hactiveg_patch, hactiveg_patch, remove_mean
    Use mod_read_obs, only: vec_useObs, vec_useObs_global
    use mod_parallel_pdaf, &
        only: mype_filter, comm_filter, &
        mpi_integer, mpi_double_precision, mpi_in_place, mpi_sum, &
        mype_world, mpi_2integer, mpi_maxloc
    use clm_varpar   , only : nlevsoi
    use decompMod , only : get_proc_bounds
    use clm_varcon, only: spval
    use clm_varctl    , only: inst_suffix
#endif
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_instMod, only : waterstate_inst
    IMPLICIT NONE
  
  ! !ARGUMENTS:
    INTEGER, INTENT(in) :: step               ! Currrent time step
    INTEGER, INTENT(in) :: dim_p              ! PE-local dimension of state
    INTEGER, INTENT(in) :: dim_obs_p          ! Dimension of observed state
    REAL, INTENT(in)    :: state_p(dim_p)     ! PE-local model state
    REAL, INTENT(out) :: m_state_p(dim_obs_p) ! PE-local observed state
    integer :: i, j, k, g
    integer :: icorner
    logical :: lpointobs       !If true: no special observation; use point observation
    ! !CALLING SEQUENCE:
    ! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT   (as U_obs_op)
    ! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
    !EOP
  
    ! *** local variables ***
  
    ! hcp test with hardcoding variable declaration
    real(8), dimension(:), allocatable :: soide !soil depth
    !real(8), dimension(0:12), parameter :: &
    ! soide=(/0.d0,  0.02d0,  0.05d0,  0.1d0,  0.17d0, 0.3d0,  0.5d0, &
    !                0.8d0,   1.3d0,   2.d0,  3.d0, 5.d0,  12.d0/) !soil depth
  
    real(8) :: tot, avesm, sum
    integer :: nsc
    ! end of hcp 
  
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: ierror
    integer :: obs_point    ! which observation is seen by which point?
  
    integer :: gridcell_index
  
    REAL:: m_state_sum(size(vec_useObs_global)) ! sum up all model grid cells and variables which correspond to an observation
    REAL:: m_state_sum_global(size(vec_useObs_global)) ! sum up all model grid cells and variables which correspond to an observation
  
    integer :: count
  
    REAL, allocatable :: tws_from_statevector(:)
  
  
  
    ! *********************************************
    ! *** Perform application of measurement    ***
    ! *** operator H on vector or matrix column ***
    ! *********************************************
  
    ! If no special observation operator is compiled, use point observations
    lpointobs = .true.
  
#if defined CLMSA
    if (clmupdate_T.EQ.1) then
  
      lpointobs = .false.
  
      DO i = 1, dim_obs_p
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
    endif
  
    if (clmupdate_tws.eq.1) then
      m_state_sum(:) = 0
      lpointobs = .false.
  
      call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  
  
      if (allocated(tws_from_statevector)) deallocate(tws_from_statevector)
      allocate(tws_from_statevector(begg:endg))
  
      tws_from_statevector(begg:endg) = spval
  
      select case(state_setup)
      case(0)
  
        do j = 1,nlevsoi
  
          do count = 1, num_layer(j)
  
            g = hactiveg_levels(count,j)
  
            if (j==1) then
              tws_from_statevector(g) = 0._r8
            end if
  
            if (j==1) then
  
              ! liq
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count)
  
              ! ice
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1))
  
            else
  
              ! liq
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count+sum(num_layer(1:j-1)))
  
              ! ice
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count+sum(num_layer(1:j-1)) + clm_varsize_tws(1))
  
            end if
  
            if (j == 1) then
  
              ! snow
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2))
  
              ! surface water
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2)+ clm_varsize_tws(3))
  
            end if 
  
          end do
  
        end do
  
        do count = 1, num_hactiveg_patch
  
          g = hactiveg_patch(count)
  
          ! canopy water
          tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2)+ clm_varsize_tws(3)+ clm_varsize_tws(4))
  
        end do
  
  
      case(1)
  
        do j = 1,nlevsoi
  
          do count = 1, num_layer(j)
  
            g = hactiveg_levels(count,j)
  
            if (j==1) then
              tws_from_statevector(g) = 0._r8
            end if
  
            if (j==1) then
  
              ! liq + ice
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count)
  
            else
  
              ! liq + ice
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count+sum(num_layer(1:j-1)))
  
            end if
  
            if (j == 1) then
  
              ! snow
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2))
  
              ! surface water
              tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2)+ clm_varsize_tws(3))
  
            end if 
  
          end do
  
        end do
  
        do count = 1, num_hactiveg_patch
  
          g = hactiveg_patch(count)
  
          ! canopy water
          tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2)+ clm_varsize_tws(3)+ clm_varsize_tws(4))
  
        end do
  
  
      case(2)
  
        do count = 1, num_layer(1)
  
          g = hactiveg_levels(count,1)
          
          tws_from_statevector(g) = state_p(count)
  
        end do
  
      case(3)
  
        do count = 1,num_layer(1)
  
          g = hactiveg_levels(count,1)
  
          tws_from_statevector(g) = state_p(count)
  
          tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2))
  
        end do
  
      case(4)
   
        do count = 1,num_layer(1)
  
          g = hactiveg_levels(count,1)
  
          tws_from_statevector(g) = state_p(count)
  
          tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2))
  
        end do
  
        do count = 1,num_layer(8)
  
          g = hactiveg_levels(count,8)
  
          tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1))
  
        end do
  
      case(5)
   
        do count = 1,num_layer(1)
  
          g = hactiveg_levels(count,1)
  
          tws_from_statevector(g) = state_p(count)
  
          tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3))
  
        end do
  
        do count = 1,num_layer(4)
  
          g = hactiveg_levels(count,4)
  
          tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1))
  
        end do
  
        do count = 1,num_layer(13)
  
          g = hactiveg_levels(count,13)
  
          tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2))
  
        end do
  
      end select
  
  
  
      do count = 1, num_layer(1)
  
        g = hactiveg_levels(count,1)
  
        obs_point = obs_id_p(g)
  
        if (obs_point /= 0) then
          ! now, the gridcell that was looked upon has been added to the sum for its corresponging observations to reproduce it. However, GRACE measures anomalies 
          ! (TWS changes). Due to this reason, a mean per gridcell has to be removed from this sum. The value from the mean corresponds to the mean per gridcell in an
          ! reference run with unperturbed forcings and surface data.
          !print*, 'difference TWS and reproduced (', g , ') = ', TWS(g)-tws_from_statevector(g)
          if (tws_temp_mean_d(g).ne.spval .and. tws_from_statevector(g).ne.spval) then
  
            if (remove_mean.eq.0) then
  
              m_state_sum(obs_point) = m_state_sum(obs_point) + tws_from_statevector(g)-tws_temp_mean_d(g)
  
            else
  
              m_state_sum(obs_point) = m_state_sum(obs_point) + tws_from_statevector(g)
  
            end if
          else if (tws_temp_mean_d(g).eq.spval .and. .not. tws_from_statevector(g).eq.spval) then
            print*, "error, tws temporal mean is spval and reproduced values is not spval for g = ", g
            print*, "reproduced = ", tws_from_statevector(g)
            stop
          else if (.not. tws_temp_mean_d(g).eq.spval .and. tws_from_statevector(g).eq.spval) then
            print*, "error, tws temporal mean is not spval and reproduced values is spvalfor g = ", g
            print*, "temp_mean = ", tws_temp_mean_d(g)
            stop
          end if
    
        end if
        
      end do
  
      call mpi_allreduce(m_state_sum, m_state_sum_global, size(vec_useObs_global), mpi_double_precision, mpi_sum, comm_filter, ierror)
  
      m_state_sum_global = m_state_sum_global/vec_numPoints_global
  
      m_state_p = pack(m_state_sum_global, vec_useObs)
  
      if (mype_filter==0) then
        print *, "m_state_global = ", m_state_sum_global
      end if
  
    end if
#endif
  
  
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
    if (crns_flag.EQ.1) then
  
      lpointobs = .false.
  
        call C_F_POINTER(soilay,soilay_fortran,[nz_glob])
        Allocate(soide(0:nz_glob))
        soide(0)=0.d0
        do i=1,nz_glob
          soide(i)=soide(i-1)+soilay_fortran(nz_glob-i+1) 
        enddo
        do i = 1, dim_obs_p
          nsc= size(sc_p(i)%scol_obs_in(:))
          avesm=0.d0
          do j=1, nsc-1
              avesm=avesm+(1.d0-0.5d0*(soide(j)+soide(j-1))/depth_obs_p(i))*(soide(j)-soide(j-1)) &
                    *state_p(sc_p(i)%scol_obs_in(j))/depth_obs_p(i)
          enddo
          avesm=avesm+(1.d0-0.5d0*(depth_obs_p(i)+soide(nsc-1))/depth_obs_p(i))*(depth_obs_p(i)-soide(nsc-1)) &
                *state_p(sc_p(i)%scol_obs_in(nsc))/depth_obs_p(i)
          tot=0.d0
          do j=1, nsc-1
              tot=tot+(1.d0-0.5d0*(soide(j)+soide(j-1))/depth_obs_p(i))*(soide(j)-soide(j-1)) &
                  /depth_obs_p(i)
          enddo
          tot=tot+(1.d0-0.5d0*(depth_obs_p(i)+soide(nsc-1))/depth_obs_p(i))*(depth_obs_p(i)-soide(nsc-1)) &
            /depth_obs_p(i)
  
          avesm=avesm/tot
          m_state_p(i)=avesm
        enddo
        deallocate(soide)
    end if
#endif
#endif
  
    if(obs_interp_switch == 1) then
  
      lpointobs = .false.
  
      do i = 1, dim_obs_p
  
          m_state_p(i) = 0
          do icorner = 1, 4
              m_state_p(i) = m_state_p(i) + state_p(obs_interp_indices_p(i,icorner)) * obs_interp_weights_p(i,icorner)
          enddo
  
      enddo
  
    end if
  
    if(lpointobs) then
  
      DO i = 1, dim_obs_p
          m_state_p(i) = state_p(obs_index_p(i))
      END DO
        
    end if
  
  END SUBROUTINE obs_op_pdaf
  
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
!obs_op_f_pdaf.F90: TSMP-PDAF implementation of routine
!                   'obs_op_f_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: obs_op_f_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: obs_op_f_pdaf --- Implementation of observation operator
!
! !INTERFACE:
SUBROUTINE obs_op_f_pdaf(step, dim_p, dim_obs_f, state_p, m_state_f)

  ! !DESCRIPTION:
  ! User-supplied routine for PDAF.
  ! Used in the filters: LSEIK/LETKF/LESTKF
  !
  ! The routine is called in PDAF\_X\_update
  ! before the loop over all local analysis domains
  ! is entered.  The routine has to perform the 
  ! operation of the observation operator acting on 
  ! a state vector.  The full vector of all 
  ! observations required for the localized analysis
  ! on the PE-local domain has to be initialized.
  ! This is usually data on the PE-local domain plus 
  ! some region surrounding the PE-local domain. 
  ! This data is gathered by MPI operations. The 
  ! gathering has to be done here, since in the loop 
  ! through all local analysis domains, no global
  ! MPI operations can be performed, because the 
  ! number of local analysis domains can vary from 
  ! PE to PE.
  !
  ! !REVISION HISTORY:
  ! 2013-09 - Lars Nerger - Initial code
  ! Later revisions - see svn log
  !
  ! !USES:
  USE mod_assimilation, &
       ONLY: obs_index_p, local_dims_obs, obs_id_p, obs_nc2pdaf_deprecated, &
       var_id_obs, dim_obs_p
  USE mod_assimilation, ONLY: tws_temp_mean_d 
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, comm_filter, MPI_DOUBLE, &
       MPI_DOUBLE_PRECISION, MPI_INT, MPI_SUM
  !USE mod_read_obs, & 
  !     ONLY: var_id_obs_nc 
#ifdef CLMSA
  use decompMod , only : get_proc_bounds
  use enkf_clm_mod, only: clmupdate_tws, clm_varsize_tws, state_setup, &
    num_layer, hactiveg_levels, num_hactiveg_patch, hactiveg_patch, remove_mean
  Use mod_read_obs, only: vec_useObs_global, vec_numPoints_global
  use clm_varcon, only: spval
  use clm_varpar   , only : nlevsoi
  use shr_kind_mod, only: r8 => shr_kind_r8
#endif
  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER, INTENT(in) :: step                 ! Current time step
  INTEGER, INTENT(in) :: dim_p                ! PE-local dimension of state
  INTEGER, INTENT(in) :: dim_obs_f            ! Dimension of observed state
  REAL, INTENT(in)    :: state_p(dim_p)       ! PE-local model state
  REAL, INTENT(inout) :: m_state_f(dim_obs_f) ! PE-local observed state

  ! !CALLING SEQUENCE:
  ! Called by: PDAF_lseik_update   (as U_obs_op)
  ! Called by: PDAF_lestkf_update  (as U_obs_op)
  ! Called by: PDAF_letkf_update   (as U_obs_op)
  !EOP

  ! local variables
  INTEGER :: ierror, max_var_id
  INTEGER :: i                         ! Counter
  INTEGER :: j                         ! Counter
  INTEGER :: g                         ! Counter
  REAL, ALLOCATABLE :: m_state_tmp(:)  ! Temporary process-local state vector
  INTEGER, ALLOCATABLE :: obs_nc2pdaf_deprecated_p_tmp(:)
  INTEGER, ALLOCATABLE :: obs_nc2pdaf_deprecated_tmp(:)

#ifdef CLMSA
  REAL:: m_state_sum(size(vec_useObs_global)) ! sum up all model grid cells and variables which correspond to an observation
  REAL:: m_state_sum_global(size(vec_useObs_global)) ! sum up all model grid cells and variables which correspond to an observation

  integer :: count

  REAL, allocatable :: tws_from_statevector(:)

  integer :: begp, endp   ! per-proc beginning and ending pft indices
  integer :: begc, endc   ! per-proc beginning and ending column indices
  integer :: begl, endl   ! per-proc beginning and ending landunit indices
  integer :: begg, endg   ! per-proc gridcell ending gridcell indices
  integer :: obs_point    ! which observation is seen by which point?
#endif

  ! *********************************************
  ! *** Perform application of measurement    ***
  ! *** operator H on vector or matrix column ***
  ! *********************************************


  if (clmupdate_tws.ne.1) then

  ! Check local observation dimension
  if (.not. local_dims_obs(mype_filter+1) == dim_obs_p) then
    print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR in local observation dimension"
    print *, "mype_filter=", mype_filter
    print *, "local_dims_obs(mype_filter+1)=", local_dims_obs(mype_filter+1)
    print *, "dim_obs_p=", dim_obs_p
    call abort_parallel()
  end if

  ! Initialize process-local observed state
  ALLOCATE(m_state_tmp(dim_obs_p))
  allocate(obs_nc2pdaf_deprecated_p_tmp (dim_obs_p))

  if(allocated(obs_nc2pdaf_deprecated)) deallocate(obs_nc2pdaf_deprecated)
  allocate(obs_nc2pdaf_deprecated(dim_obs_f))
  if(allocated(obs_nc2pdaf_deprecated_tmp)) deallocate(obs_nc2pdaf_deprecated_tmp)
  allocate(obs_nc2pdaf_deprecated_tmp(dim_obs_f))

  DO i = 1, dim_obs_p
     m_state_tmp(i) = state_p(obs_index_p(i))
     obs_nc2pdaf_deprecated_p_tmp(i)  = obs_id_p(obs_index_p(i))
  END DO
  
  !print *,'local_dims_obs(mype_filter+1) ', local_dims_obs(mype_filter+1)
  !print *,'dim_obs_p ', dim_obs_p

  ! Gather full observed state using local_dims_obs, local_disp_obs

  ! gather local observed states of different sizes in a vector
  CALL mpi_allgatherv(m_state_tmp, dim_obs_p, &
       MPI_DOUBLE_PRECISION, m_state_f, local_dims_obs, local_disp_obs, &
       MPI_DOUBLE_PRECISION, comm_filter, ierror)

  ! gather obs_nc2pdaf_deprecated_p
  CALL mpi_allgatherv(obs_nc2pdaf_deprecated_p_tmp, dim_obs_p, &
       MPI_INT, obs_nc2pdaf_deprecated, local_dims_obs, local_disp_obs, &  
       MPI_INT, comm_filter, ierror)

  ! At this point OBS_NC2PDAF_DEPRECATED should be the same as OBS_PDAF2NC from
  ! INIT_DIM_OBS_PDAF / INIT_DIM_OBS_F_PDAF
  do i = 1, dim_obs_f
    if(.not. obs_nc2pdaf_deprecated(i) .eq. obs_pdaf2nc(i)) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR in observation index arrays"
      print *, "i=", i
      print *, "obs_nc2pdaf_deprecated(i)=", obs_nc2pdaf_deprecated(i)
      print *, "obs_pdaf2nc(i)=", obs_pdaf2nc(i)
      call abort_parallel()
    end if
  end do

  ! Then OBS_NC2PDAF_DEPRECATED is inverted in the following lines

  ! resort obs_nc2pdaf_deprecated_p
  do i=1,dim_obs_f
     obs_nc2pdaf_deprecated_tmp(i) = obs_nc2pdaf_deprecated(i)
  enddo
  do i=1,dim_obs_f
    ! print *,'obs_nc2pdaf_deprecated_tmp(i) ', obs_nc2pdaf_deprecated_tmp(i)
     obs_nc2pdaf_deprecated(obs_nc2pdaf_deprecated_tmp(i)) = i
  enddo

  ! At this point OBS_NC2PDAF_DEPRECATED should be the same as OBS_NC2PDAF from
  ! INIT_DIM_OBS_PDAF / INIT_DIM_OBS_F_PDAF
  do i = 1, dim_obs_f
    if(.not. obs_nc2pdaf_deprecated(i) .eq. obs_nc2pdaf(i)) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR in observation index arrays"
      print *, "i=", i
      print *, "obs_nc2pdaf_deprecated(i)=", obs_nc2pdaf_deprecated(i)
      print *, "obs_nc2pdaf(i)=", obs_nc2pdaf(i)
      call abort_parallel()
    end if
  end do

  ! Clean up
  DEALLOCATE(m_state_tmp,obs_nc2pdaf_deprecated_p_tmp,obs_nc2pdaf_deprecated_tmp)

  else

#ifdef CLMSA
    m_state_sum(:) = 0

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

    m_state_f = pack(m_state_sum_global, vec_useObs_global)

    if (mype_filter==0) then
      print *, "m_state_global = ", m_state_sum_global
    end if

#endif
  end if

END SUBROUTINE obs_op_f_pdaf

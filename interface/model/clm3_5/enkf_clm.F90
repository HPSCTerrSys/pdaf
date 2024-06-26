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
!enkf_clm.F90: Wrapper functions for CLM
!-------------------------------------------------------------------------------------------

#include <misc.h>
#include <preproc.h>

!> @author Wolfgang Kurtz, Guowei He
!> @date 25.02.2022
!> @brief CLM3.5 Initialization routine for TSMP-PDAF
!> @param[in] finname CLM input file name
!> @details
!> Initialization routines for the land model analogous to
!> clm3_5/src/main/program_off.F90
subroutine clm_init(finname, pdaf_id, pdaf_max, mype) bind(C,name="clm_init")
  use iso_C_binding
  use enkf_clm_mod, only: &
#if (defined COUP_OAS_COS || defined COUP_OAS_PFL)
      kl_comm, &
#elif (defined CLMSA)
      COMM_model_clm, &
#else
      mpi_running, mpicom_glob, ier, &
#endif
#if defined CLMSA
      define_clm_statevec, &
#endif
      spmd_init, mpicom, comp_id, masterproc, &
      clmprefixlen, nlfilename, &
      ESMF_Initialize, &
      control_setNL, &
      log_print, &
      eccen, mvelpp, lambm0, obliqr, obliq, &
      iyear_AD, nmvelp, &
      shr_orb_params, SHR_ORB_UNDEF_REAL, &
      clm_init0, clm_init1, clm_init2, &
      atmdrv, atmdrv_init, &
      mct_world_init

  implicit none

!  character(c_char),target   :: finname
  character(kind=c_char,len=1),dimension(100),intent(in) :: finname
  integer(c_int), intent(in) :: pdaf_id !unused: only for clm5_0
  integer(c_int), intent(in) :: pdaf_max !unused: only for clm5_0
  integer(c_int), intent(in) :: mype
  integer(c_int) :: counter
  !character(100),pointer :: pchar

  !call c_f_pointer(c_loc(finname),pchar)
  !s_len = index(pchar,c_null_char)-1
  !nlfilename = trim(pchar(1:s_len))

  loop_string: do counter=1,clmprefixlen
        nlfilename(counter:counter) = finname(counter)
  end do loop_string

  ! -----------------------------------------------------------------
  ! Initialize MPI
  ! -----------------------------------------------------------------

#if (defined COUP_OAS_COS || defined COUP_OAS_PFL)
  call oas_clm_init
  call spmd_init(kl_comm)
  call mct_world_init(1,kl_comm,mpicom,comp_id)
#elif (defined CLMSA)
  call spmd_init(COMM_model_clm)
  call mct_world_init(1,COMM_model_clm,mpicom,comp_id)
#else
  call mpi_initialized (mpi_running, ier)
  if (.not. mpi_running) call mpi_init(ier)
  mpicom_glob = MPI_COMM_WORLD
  call spmd_init(mpicom_glob)
  call mct_world_init(1,mpicom_glob,mpicom,comp_id)
#endif

  ! call t_startf('init')

  ! -----------------------------------------------------------------
  ! Initialize ESMF (needed for time-manager)
  ! -----------------------------------------------------------------

  call ESMF_Initialize()

  ! -----------------------------------------------------------------
  ! Initialize timing library, and set full path to namelist
  ! -----------------------------------------------------------------

  call control_setNL( nlfilename )     ! Set namelist
  ! call t_initf(nlfilename, LogPrint=masterproc, Mpicom=mpicom, &
  !              MasterTask=masterproc)

  ! -----------------------------------------------------------------
  ! Initialize Orbital parameters
  ! -----------------------------------------------------------------

  ! obliq, eccen and nmvelp are determined based on value of iyear_AD

  if (masterproc) then
     log_print = .true.
  else
     log_print = .false.
  end if
  iyear_AD = 1950
  obliq    = SHR_ORB_UNDEF_REAL
  eccen    = SHR_ORB_UNDEF_REAL
  nmvelp   = SHR_ORB_UNDEF_REAL
  call shr_orb_params (iyear_AD, eccen, obliq, nmvelp, obliqr, &
                       lambm0, mvelpp, log_print)

  ! -----------------------------------------------------------------
  ! Initialize land model
  ! -----------------------------------------------------------------

  call clm_init0()
  call clm_init1()
  call clm_init2()

  ! -----------------------------------------------------------------
  ! Initialize "external" atmospheric forcing
  ! -----------------------------------------------------------------

  !if (masterproc) write (6,*) 'Attempting to set up atmospheric grid '
  call atmdrv_init()
  !if (masterproc) write (6,*) 'Successfully set up atmospheric grid '


#if defined CLMSA
  call define_clm_statevec(mype)
#endif

end subroutine clm_init


subroutine clm_advance(ntstep, tstartcycle, mype) bind(C,name="clm_advance")
  use enkf_clm_mod
  use iso_c_binding
  use mod_clm_statistics
  integer(c_int),intent(in) :: ntstep
  integer(c_int),intent(in) :: tstartcycle
  integer(c_int),intent(in) :: mype
  integer :: counter=0

   do counter=1,ntstep
     nstep = get_nstep()
     call atmdrv(nstep)

#if defined COUP_OAS_PFL
    !  received fields from ParFlow
     call receive_fld_2pfl(nstep-1)
    !
#endif

     call clm_run1()


#if defined COUP_OAS_COS
     ! send fields to Cosmo
     call send_fld_2cos
     !
#endif

#if defined COUP_OAS_PFL
    !  send flux fields to ParFlow
  call send_fld_2pfl
    !
#endif

     call clm_run2()

    ! Determine if time to stop

    if(.NOT.is_last_step()) call advance_timestep()
    !if(counter.eq.ntstep) call write_clm_statistics()

!  end do
  end do

#if defined CLMSA
  ! Calling PDAF Function to set state vector before assimiliation
  call set_clm_statevec(tstartcycle, mype)
#endif

end subroutine clm_advance


subroutine clm_finalize() bind(C,name="clm_finalize")
  use enkf_clm_mod
#if (defined BGL)
       call print_stack_size()
#endif

    if (masterproc) then
       write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
       call flush(6)      !CPS
    endif

#if (defined COUP_OAS_COS || defined COUP_OAS_PFL)
    ! Let OASIS finalize run
    call oas_clm_finalize
#else
    ! Finalize ESMF
!    call ESMF_Finalize()
#endif


!endif
end subroutine clm_finalize

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
!next_observation_pdaf.F90: TSMP-PDAF implementation of routine
!                           'next_observation_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: next_observation_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: next_observation_pdaf --- Initialize information on next observation
!
! !INTERFACE:
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

  ! !DESCRIPTION:
  ! User-supplied routine for PDAF.
  ! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
  !
  ! The subroutine is called before each forecast phase
  ! by PDAF\_get\_state. It has to initialize the number 
  ! of time steps until the next available observation 
  ! (nsteps) and the current model time (time). In 
  ! addition the exit flag (exit) has to be initialized.
  ! It indicates if the data assimilation process is 
  ! completed such that the ensemble loop in the model 
  ! routine can be exited.
  !
  ! The routine is called by all processes. 
  !
  ! !REVISION HISTORY:
  ! 2013-09 - Lars Nerger - Initial code
  ! Later revisions - see svn log
  !
  ! !USES:
    USE mod_assimilation, &
         ONLY: delt_obs, toffset, screen, da_interval_variable
    USE mod_parallel_pdaf, &
         ONLY: mype_world
    USE mod_tsmp, &
         ONLY: total_steps
    USE mod_assimilation, &
         ONLY: obs_filename
    use mod_read_obs, &
         only: check_n_observationfile, check_n_observationfile_da_interval, check_n_observationfile_set_zero
    use clm_time_manager, &
         only: get_nstep
    use enkf_clm_mod, &
         only: da_interval
    use clm_varcon, only: set_averaging_to_zero, ispval
    IMPLICIT NONE
  
  ! !ARGUMENTS:
    INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step
    INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
    INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
    REAL, INTENT(out)    :: time     ! Current model (physical) time
  
  ! !CALLING SEQUENCE:
  ! Called by: PDAF_get_state   (as U_next_obs)
  !EOP
  
    !kuw: local variables
    integer :: counter,no_obs=0, nstep
    character (len = 110) :: fn
    !kuw end
    
    time = 0.0    ! Not used in fully-parallel implementation variant
    doexit = 0
    nstep = get_nstep()
    nsteps = delt_obs

    if (mype_world==0 .and. screen > 2) then
      write(*,*) 'TSMP-PDAF (in next_observation_pdaf.F90) total_steps: ',total_steps
    end if
    ! Read steps until next observation from current observation file
    if (stepnow.eq.toffset) then
      set_averaging_to_zero = 0
      if (mype_world==0 .and. screen > 2) then
        write(*,*)'next_observation_pdaf: da_interval from enkfpf.par'
      end if
    else
      write(fn, '(a, i5.5)') trim(obs_filename)//'.', stepnow
      call check_n_observationfile_da_interval(fn,da_interval_variable)
      if (da_interval_variable.ne.ispval) then
        da_interval = da_interval_variable
      end if
      call check_n_observationfile_set_zero(fn, set_averaging_to_zero)
    end if
    if (mype_world==0 .and. screen > 2) then
      write(fn, '(a, i5.5)') trim(obs_filename)//'.', stepnow+delt_obs
      write(*,*)'next_observation_pdaf: fn = ', fn
      write(*,*)'da_interval (in next_observation_pdaf):',da_interval
    end if
    if (set_averaging_to_zero.ne.ispval) then
      set_averaging_to_zero = set_averaging_to_zero+nstep
    end if
    if (mype_world==0 .and. screen > 2) then
      write(*,*) 'set_averaging_to_zero (in next_observation_pdaf):',set_averaging_to_zero
    end if
  
  END SUBROUTINE next_observation_pdaf
  
  
  
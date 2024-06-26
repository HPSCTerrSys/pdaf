!> Main program
!!
!! The program only serves to be able to compile
!! the PDAF online template routines for testing
!! their consistency.
!!
!! The program shows the setup for the fully-parallel
!! implementation variant of PDAF.
!! 
!! In the online implementation with a real model
!! this driver program is replaced by the actual
!! model code.
!!
!! __Revision history:__
!! * 2021-11 - Lars Nerger - Initial code
!!
PROGRAM MAIN

  USE mpi                      ! MPI
  USE mod_parallel_pdaf, &     ! Parallelization
       ONLY: init_parallel, finalize_parallel, &
       n_modeltasks, npes_world, mype_world
  USE mod_model, &             ! Module provided by model code
       ONLY: step_final        

  IMPLICIT NONE

! local variables
  INTEGER :: istep       ! Counter


! *** Initialize MPI ***

  ! If the model itself is parallelized this step is done by the model

  CALL init_parallel() ! initializes MPI

  ! FOR TESTING: 
  n_modeltasks = npes_world

  IF (mype_world==0) THEN
     WRITE (*,*) '**********************************************************************'
     WRITE (*,*) '*   THIS IS A TEST PROGRAM TO CHECK THE TEMPLATE CODE CONSISTENCY    *'
     WRITE (*,*) '*                   Run this program with:                           *'
     WRITE (*,*) '*          mpirun -np NENS ./PDAF_online -dim_ens NENS               *'
     WRITE (*,*) '* with ensemble size NENS (=2 is good for testing, =1 does not work) *'
     WRITE (*,*) '**********************************************************************'
  END IF

  
! *** Initialize MPI communicators for PDAF (model, filter and coupling) ***

  ! This step is always inserted directly after the MPI initialization

  CALL init_parallel_pdaf(0, 1)

  
  ! MODEL: Here the model would perform its initialization


! *** Initialize PDAF ***

  ! This step is always inserted after the model initialization
  ! is complete and just before the time stepping starts

  CALL init_pdaf()


! *** Ensemble forecasting and analysis steps ***

  ! MODEL: In the real model this is the time stepping loop of the model
  timesteps: DO istep = 1, step_final

     ! MODEL: Here the model code would compute the time stepping


     ! *** Perform analysis ***

     ! This step is inserted in the time stepping loop
     ! usually just before the 'end do'

     CALL assimilate_pdaf()

  ENDDO timesteps


! *** Finalize PDAF - print memory and timing information ***

  ! This step can be inserted at the end of the model code
  ! before the MPI parallelization is finalized

  CALL finalize_pdaf(0)


! *** Terminate MPI

  ! If the model itself is parallelized this step is done by the model

  CALL finalize_parallel()

END PROGRAM MAIN

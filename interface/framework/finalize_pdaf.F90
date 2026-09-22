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
!finalize_pdaf.F90: TSMP-PDAF implementation of routine
!                   'finalize_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: finalize_pdaf.F90 1442 2013-10-04 10:35:19Z lnerger $
!BOP
!
! !ROUTINE: finalize_pdaf --- Finalize PDAF
!
! !INTERFACE:
SUBROUTINE finalize_pdaf()

! !DESCRIPTION:
! This routine call MPI_finalize
!
! !REVISION HISTORY:
! 2004-11 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, &      ! Variables for assimilation
       ONLY: dim_state_p_count, obs_p, &
             obs_index_p, &
             obs_index_l, global_to_local, &
             local_dims_obs, &
             local_disp_obs
  USE mod_parallel_pdaf, &
       ONLY: local_npes_model, mype_world
#ifdef CLMSA
#ifdef CLMFIVE
  USE obs_LST_pdafomi, &
       ONLY: t_LST_readobs, t_LST_snapping, t_LST_gather, &
             t_LST_total, t_LST_calls
#endif
#endif

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main program
!EOP

! *** Show allocated memory for PDAF ***
  IF (mype_world==0) CALL PDAF_print_info(10)

! *** Print PDAF timings onto screen ***
  IF (mype_world==0) CALL PDAF_print_info(3)
  IF (mype_world==0) CALL PDAF_print_info(5)
  IF (mype_world==0) CALL PDAF_print_info(1)

#ifdef CLMSA
#ifdef CLMFIVE
! *** Print LST-OMI sub-timers ***
  IF (mype_world==0 .AND. t_LST_calls > 0) THEN
    WRITE (*,'(a)') ''
    WRITE (*,'(a)') 'TSMP-PDAF    LST-OMI init_dim_obs sub-timers (rank 0, accumulated)'
    WRITE (*,'(a)') 'TSMP-PDAF    --------------------------------------------------------'
    WRITE (*,'(a,i6)')     'TSMP-PDAF      Number of calls:              ', t_LST_calls
    WRITE (*,'(a,f12.3,a)') 'TSMP-PDAF      Obs I/O + broadcast:     ', t_LST_readobs, ' s'
    WRITE (*,'(a,f12.3,a)') 'TSMP-PDAF      Snapping loops + MPI:    ', t_LST_snapping, ' s'
    WRITE (*,'(a,f12.3,a)') 'TSMP-PDAF      PDAFomi_gather_obs:      ', t_LST_gather, ' s'
    WRITE (*,'(a,f12.3,a)') 'TSMP-PDAF      Total init_dim_obs_LST:  ', t_LST_total, ' s'
    WRITE (*,'(a,f12.3,a)') 'TSMP-PDAF      Unaccounted (domain_lim):', &
      t_LST_total - t_LST_readobs - t_LST_snapping - t_LST_gather, ' s'
  END IF
#endif
#endif

! *** Deallocate PDAF arrays ***
  CALL PDAF_deallocate()

! *** Deallocate PDAF framework arrays ***
  if (allocated(local_npes_model)) deallocate(local_npes_model)
  IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  if (allocated(dim_state_p_count)) deallocate (dim_state_p_count)
  ! M.Pondkule: deallocating variables used in data assimilation
  ! with letkf filter
  IF (ALLOCATED(obs_index_l)) DEALLOCATE(obs_index_l)
  IF (ALLOCATED(global_to_local)) DEALLOCATE(global_to_local)
  IF (ALLOCATED(local_dims_obs)) DEALLOCATE(local_dims_obs)
  IF (ALLOCATED(local_disp_obs)) DEALLOCATE(local_disp_obs)

! *** Finalize parallel MPI region - if not done by model ***
!  CALL finalize_parallel()

END SUBROUTINE finalize_pdaf

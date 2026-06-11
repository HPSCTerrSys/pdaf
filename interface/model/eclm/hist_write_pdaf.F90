!-------------------------------------------------------------------------------------------
!Copyright (c) 2026 by HPSCTerrSys (Forschungszentrum Juelich GmbH)
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
!GNU Lesser General Public License for more details.
!
!You should have received a copy of the GNU Lesser General Public License
!along with TSMP-PDAF.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------
! hist_write_pdaf.F90: Write CLM history file snapshots around a PDAF DA update.
!
! Two C-callable entry points:
!   clm_hist_write_da_before  -- call before update_clm(); captures forecast state
!   clm_hist_write_da_after   -- call after  update_clm(); captures analysis state
!
! Both set the da_tape_phase label (used in the filename), populate the dedicated
! DA tape buffers with the current model state, and trigger hist_htapes_wrapup
! with da_call=.true. so that only the DA tape is written.
!-------------------------------------------------------------------------------------------

#if defined CLMSA

! ---------------------------------------------------------------------------
! Internal helper: set phase, update DA tape buffer, and write.
! ---------------------------------------------------------------------------
subroutine clm_hist_write_da(phase)

  use decompMod,   only : get_proc_bounds, bounds_type
  use clm_varpar,  only : nlevgrnd
  use clm_instMod, only : soilstate_inst
  use histFileMod, only : hist_set_da_tape_phase, hist_update_hbuf_da, &
                          hist_htapes_wrapup

  implicit none

  character(len=*), intent(in) :: phase

  type(bounds_type) :: bounds

  call get_proc_bounds(bounds)
  call hist_set_da_tape_phase(phase)
  call hist_update_hbuf_da(bounds)
  call hist_htapes_wrapup( &
       rstwr      = .false., &
       nlend      = .false., &
       bounds     = bounds, &
       watsat_col = soilstate_inst%watsat_col(bounds%begc:bounds%endc, 1:nlevgrnd), &
       sucsat_col = soilstate_inst%sucsat_col(bounds%begc:bounds%endc, 1:nlevgrnd), &
       bsw_col    = soilstate_inst%bsw_col(bounds%begc:bounds%endc, 1:nlevgrnd),    &
       hksat_col  = soilstate_inst%hksat_col(bounds%begc:bounds%endc, 1:nlevgrnd),  &
       da_call    = .true.)

end subroutine clm_hist_write_da

! ---------------------------------------------------------------------------
! C-callable: snapshot before the DA update (forecast / prior state)
! ---------------------------------------------------------------------------
subroutine clm_hist_write_da_before() bind(C, name="clm_hist_write_da_before")
  implicit none
  call clm_hist_write_da('bef')
end subroutine clm_hist_write_da_before

! ---------------------------------------------------------------------------
! C-callable: snapshot after the DA update (analysis / posterior state)
! ---------------------------------------------------------------------------
subroutine clm_hist_write_da_after() bind(C, name="clm_hist_write_da_after")
  implicit none
  call clm_hist_write_da('aft')
end subroutine clm_hist_write_da_after

#endif

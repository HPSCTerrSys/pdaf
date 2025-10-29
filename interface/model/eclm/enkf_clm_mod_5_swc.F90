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
!enkf_clm_mod_5_swc.F90: Module for CLM - SWC observations
!-------------------------------------------------------------------------------------------
module enkf_clm_mod_swc

  implicit none

  public

  contains

  subroutine define_clm_statevec_swc(mype)
    use decompMod , only : get_proc_bounds
    use clm_varpar   , only : nlevsoi
    use clm_varcon , only : ispval
    use ColumnType , only : col

    implicit none

    integer,intent(in) :: mype

    integer :: i
    integer :: c
    integer :: g
    integer :: cc

    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices


    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

#ifdef PDAF_DEBUG
    WRITE(*,"(a,i5,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a)") &
      "TSMP-PDAF mype(w)=", mype, " define_clm_statevec, CLM5-bounds (g,l,c,p)----",&
      begg,",",endg,",",begl,",",endl,",",begc,",",endc,",",begp,",",endp," -------"
#endif

    clm_begg     = begg
    clm_endg     = endg
    clm_begc     = begc
    clm_endc     = endc
    clm_begp     = begp
    clm_endp     = endp

    ! Soil Moisture DA: State vector index arrays

      ! 1) COL/GRC: CLM->PDAF
      IF (allocated(state_clm2pdaf_p)) deallocate(state_clm2pdaf_p)
      allocate(state_clm2pdaf_p(begc:endc,nlevsoi))
      do i=1,nlevsoi
        do c=clm_begc,clm_endc
          ! Default: inactive
          state_clm2pdaf_p = ispval
        end do
      end do

      ! All column variables in state vector
      if(clmstatevec_allcol==1) then

        ! Only hydrologically active columns
        if(clmstatevec_only_active == 1) then

          cc = 0

          do i=1,nlevsoi
            ! Only take into account layers above input maximum layer
            if(i<=clmstatevec_max_layer) then

              do c=clm_begc,clm_endc
                ! Only take into account hydrologically active columns
                ! and layers above bedrock
                if(col%hydrologically_active(c) .and. i<=col%nbedrock(c)) then
                  cc = cc + 1
                  state_clm2pdaf_p(c,i) = cc
                end if
              end do

            end if
          end do

          ! All column variables in state vector simplifying the indexing
        else

          do i=1,nlevsoi
            do c=clm_begc,clm_endc
              state_clm2pdaf_p(c,i) = (c - clm_begc + 1) + (i - 1)*(clm_endc - clm_begc + 1)
            end do
          end do

        end if

      ! Gridcell values or averages in state vector
      else

        ! Only hydrologically active columns
        if(clmstatevec_only_active==1) then

          cc = 0

          do i=1,nlevsoi
            ! Only layers above max_layer
            if(i<=clmstatevec_max_layer) then

              do g=clm_begg,clm_endg

                newgridcell = .true.

                do c=clm_begc,clm_endc
                  if(col%gridcell(c) == g) then
                    ! All (hydrologically active / above bedrock)
                    ! column-layer pairs that belong to a gridcell
                    ! point to the state vector index of the
                    if(col%hydrologically_active(c) .and. i<=col%nbedrock(c)) then
                      if(newgridcell) then
                        ! Update the index if first col found for grc,
                        ! otherwise reproduce previous index
                        cc = cc + 1
                        newgridcell = .false.
                      end if
                      state_clm2pdaf_p(c,i) = cc
                    end if
                  end if
                end do

              end do
            end if
          end do
        else
          do i=1,nlevsoi
            do c=clm_begc,clm_endc
              ! All columns in a gridcell are assigned the updated
              ! gridcell-SWC
              state_clm2pdaf_p(c,i) = (col%gridcell(c) - clm_begg + 1) + (i - 1)*(clm_endg - clm_begg + 1)
            end do
          end do
        end if

      end if

      ! 2) COL/GRC: STATEVECSIZE
      if(clmstatevec_only_active==1) then
        ! Use iterator cc for setting state vector size.
        !
        ! Set `clm_varsize`, even though it is currently not used
        ! for `clmupdate_swc.eq.1`
        clm_varsize      =  cc
        clm_statevecsize =  cc
      else
        if(clmstatevec_allcol==1) then
          ! #cols * #levels
          clm_varsize      =  (endc-begc+1) * nlevsoi
          clm_statevecsize =  (endc-begc+1) * nlevsoi
        else
          ! #grcs * #levels
          clm_varsize      =  (endg-begg+1) * nlevsoi
          clm_statevecsize =  (endg-begg+1) * nlevsoi
        end if
      end if

      ! 3) COL/GRC: PDAF->CLM
      IF (allocated(state_pdaf2clm_c_p)) deallocate(state_pdaf2clm_c_p)
      allocate(state_pdaf2clm_c_p(clm_statevecsize))
      IF (allocated(state_pdaf2clm_j_p)) deallocate(state_pdaf2clm_j_p)
      allocate(state_pdaf2clm_j_p(clm_statevecsize))

      ! Defaults
      do cc=1,clm_statevecsize
        state_pdaf2clm_c_p(cc) = ispval
        state_pdaf2clm_j_p(cc) = ispval
      end do

      do cc=1,clm_statevecsize

        lay: do i=1,nlevsoi
          do c=clm_begc,clm_endc
            if (state_clm2pdaf_p(c,i) == cc) then
              ! Set column index and then exit loop
              state_pdaf2clm_c_p(cc) = c
              state_pdaf2clm_j_p(cc) = i
              exit lay
            end if
          end do
        end do lay

#ifdef PDAF_DEBUG
        ! Check that all state vectors have been assigned c, i
        if(state_pdaf2clm_c_p(cc) == ispval) then
          write(*,*) 'cc: ', cc
          error stop "state_pdaf2clm_c_p not set at cc"
        end if
        if(state_pdaf2clm_j_p(cc) == ispval) then
          write(*,*) 'cc: ', cc
          error stop "state_pdaf2clm_j_p not set at cc"
        end if
#endif
      end do

  end subroutine define_clm_statevec_swc


  subroutine set_clm_statevec_swc()
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_varpar   , only : nlevsoi
    use ColumnType , only : col
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    real(r8), pointer :: swc(:,:)
    integer :: j,g,cc,c
    integer :: n_c

    cc = 0

    swc   => waterstate_inst%h2osoi_vol_col

      ! write swc values to state vector
      if (clmstatevec_colmean==1) then

        do cc = 1, clm_statevecsize

          clm_statevec(cc) = 0.0
          n_c = 0

          ! Get gridcell and layer
          g = col%gridcell(state_pdaf2clm_c_p(cc))
          j = state_pdaf2clm_j_p(cc)

          ! Loop over all columns
          do c=clm_begc,clm_endc
            ! Select columns in gridcell g
            if(col%gridcell(c)==g) then
              ! Select hydrologically active columns
              if(col%hydrologically_active(c)) then
                ! Add active column to swc-sum
                clm_statevec(cc) = clm_statevec(cc) + swc(c,j)
                n_c = n_c + 1
              end if
            end if
          end do

          if(n_c == 0) then
            write(*,*) "WARNING: Gridcell g=", g
            write(*,*) "WARNING: Layer    j=", j
            write(*,*) "Grid cell g at layer j without hydrologically active column! Setting SWC as in gridcell mode."
            clm_statevec(cc) = swc(state_pdaf2clm_c_p(cc), state_pdaf2clm_j_p(cc))
          else
            ! Normalize sum to average
            clm_statevec(cc) = clm_statevec(cc) / real(n_c, r8)
          end if

          ! Save prior column mean state vector for computing
          ! increment in updating the state vector
          clm_statevec_orig(cc) = clm_statevec(cc)

        end do

      else
        do cc = 1, clm_statevecsize
          clm_statevec(cc) = swc(state_pdaf2clm_c_p(cc), state_pdaf2clm_j_p(cc))
        end do
      end if

  end subroutine set_clm_statevec_swc


  subroutine update_clm_swc(tstartcycle, mype)
    use clm_varpar   , only : nlevsoi
    use shr_kind_mod , only : r8 => shr_kind_r8
    use ColumnType , only : col
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_varcon      , only : denh2o, denice, watmin
    use clm_varcon      , only : ispval
    use clm_varcon      , only : spval

    implicit none

    integer,intent(in) :: tstartcycle
    integer,intent(in) :: mype

    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: watsat(:,:)

    real(r8), pointer :: dz(:,:)          ! layer thickness depth (m)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)
    real(r8), pointer :: snow_depth(:)
    real(r8)  :: rliq,rice
    real(r8)  :: watmin_check      ! minimum soil moisture for checking clm_statevec (mm)
    real(r8)  :: watmin_set        ! minimum soil moisture for setting swc (mm)
    real(r8)  :: swc_update        ! updated SWC in loop

    integer :: i,j,cc
    character (len = 31) :: fn2    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn3    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn4    !TSMP-PDAF: function name for state vector outpu

    logical :: swc_zero_before_update

    cc = 0
    swc_zero_before_update = .false.

    swc   => waterstate_inst%h2osoi_vol_col
    watsat => soilstate_inst%watsat_col
    dz            => col%dz
    h2osoi_liq    => waterstate_inst%h2osoi_liq_col
    h2osoi_ice    => waterstate_inst%h2osoi_ice_col

    snow_depth => waterstate_inst%snow_depth_col ! snow height of snow covered area (m)

        ! Set minimum soil moisture for checking the state vector and
        ! for setting minimum swc for CLM
        if(clmwatmin_switch==3) then
          ! CLM3.5 type watmin
          watmin_check = 0.00
          watmin_set = 0.05
        else if(clmwatmin_switch==5) then
          ! CLM5.0 type watmin
          watmin_check = watmin
          watmin_set = watmin
        else
          ! Default
          watmin_check = 0.0
          watmin_set = 0.0
        end if

        ! cc = 0
        do i=1,nlevsoi
          ! CLM3.5: iterate over grid cells
          ! CLM5.0: iterate over columns
          ! do j=clm_begg,clm_endg
            do j=clm_begc,clm_endc

              ! If snow is masked, update only, when snow depth is less than 1mm
              if( (.not. clmswc_mask_snow) .or. snow_depth(j) < 0.001 ) then
              ! Update only those SWCs that are not excluded by ispval
              if(state_clm2pdaf_p(j,i) /= ispval) then

                if(swc(j,i)==0.0) then
                  swc_zero_before_update = .true.

                  ! Zero-SWC leads to zero denominator in computation of
                  ! rliq/rice, therefore setting rliq/rice to special
                  ! value
                  rliq = spval
                  rice = spval
                else
                  swc_zero_before_update = .false.

                  rliq = h2osoi_liq(j,i)/(dz(j,i)*denh2o*swc(j,i))
                  rice = h2osoi_ice(j,i)/(dz(j,i)*denice*swc(j,i))
                  !h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
                end if

                if (clmstatevec_colmean==1) then
                  ! If there is no significant increment, do not
                  ! implement any update / check.
                  !
                  ! Note: Computing the absolute difference here,
                  ! because the whole state vector should be soil
                  ! moistures. For variables with very small values in
                  ! the state vector, this would have to be adapted
                  ! (e.g. to relative difference).
                  if( abs(clm_statevec(state_clm2pdaf_p(j,i)) - clm_statevec_orig(state_clm2pdaf_p(j,i))) <= 1.0e-7) then
                    cycle
                  end if

                  ! Update SWC column value with the increment-factor
                  ! of the state vector update (state vector updates
                  ! are means of cols in grc)
                  swc_update = swc(j,i) * clm_statevec(state_clm2pdaf_p(j,i)) / clm_statevec_orig(state_clm2pdaf_p(j,i))
                else
                  ! Update SWC with updated state vector
                  swc_update = clm_statevec(state_clm2pdaf_p(j,i))
                end if

                if(swc_update<=watmin_check) then
                  swc(j,i) = watmin_set
                else if(swc_update>=watsat(j,i)) then
                  swc(j,i) = watsat(j,i)
                else
                  swc(j,i)   = swc_update
                endif

                if (isnan(swc(j,i))) then
                  swc(j,i) = watmin_set
                  print *, "WARNING: swc at j,i is nan: ", j, i
                endif

                if(swc_zero_before_update) then
                  ! This case should not appear for hydrologically
                  ! active columns/layers, where always: swc > watmin
                  !
                  ! If you want to make sure that no zero SWCs appear in
                  ! the code, comment out the error stop

#ifdef PDAF_DEBUG
                  ! error stop "ERROR: Update of zero-swc"
                  print *, "WARNING: Update of zero-swc"
                  print *, "WARNING: Any new H2O added to h2osoi_liq(j,i) with j,i = ", j, i
#endif
                  h2osoi_liq(j,i) = swc(j,i) * dz(j,i)*denh2o
                  h2osoi_ice(j,i) = 0.0
                else
                  ! update liquid water content
                  h2osoi_liq(j,i) = swc(j,i) * dz(j,i)*denh2o*rliq
                  ! update ice content
                  h2osoi_ice(j,i) = swc(j,i) * dz(j,i)*denice*rice
                end if

              end if
              end if
              ! cc = cc + 1
            end do
        end do

#ifdef PDAF_DEBUG
        IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN

            ! TSMP-PDAF: For debug runs, output the state vector in files
            WRITE(fn3, "(a,i5.5,a,i5.5,a)") "h2osoi_liq", mype, ".update.", tstartcycle, ".txt"
            OPEN(unit=71, file=fn3, action="write")
            WRITE (71,"(es22.15)") h2osoi_liq(:,:)
            CLOSE(71)

            ! TSMP-PDAF: For debug runs, output the state vector in files
            WRITE(fn4, "(a,i5.5,a,i5.5,a)") "h2osoi_ice", mype, ".update.", tstartcycle, ".txt"
            OPEN(unit=71, file=fn4, action="write")
            WRITE (71,"(es22.15)") h2osoi_ice(:,:)
            CLOSE(71)

            ! TSMP-PDAF: For debug runs, output the state vector in files
            WRITE(fn2, "(a,i5.5,a,i5.5,a)") "swcstate_", mype, ".update.", tstartcycle, ".txt"
            OPEN(unit=71, file=fn2, action="write")
            WRITE (71,"(es22.15)") swc(:,:)
            CLOSE(71)

        END IF
#endif

  end subroutine update_clm_swc


end module enkf_clm_mod_swc

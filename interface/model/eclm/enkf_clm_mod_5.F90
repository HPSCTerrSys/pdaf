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
!enkf_clm_mod.F90: Module for CLM
!-------------------------------------------------------------------------------------------

module enkf_clm_mod

  use iso_c_binding

! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8, SHR_KIND_CL

! !ARGUMENTS:
    implicit none

#if (defined CLMSA)
  integer :: COMM_model_clm
  integer :: clm_statevecsize
  integer :: clm_paramsize !hcp: Size of CLM parameter vector (f.e. LAI)
  integer :: clm_varsize
  integer :: clm_begg,clm_endg
  integer :: clm_begl,clm_endl
  integer :: clm_begc,clm_endc
  integer :: clm_begp,clm_endp
  real(r8),allocatable :: clm_statevec(:)
  real(r8),allocatable :: clm_statevec_orig(:)
  real(r8),allocatable :: clm_statevec_original_input(:) ! orginal values in statevector so that I can also access them in the update
  integer,allocatable :: state_pdaf2clm_c_p(:)
  integer,allocatable :: state_pdaf2clm_j_p(:)
  ! clm_paramarr: Contains LAI used in obs_op_pdaf for computing model
  ! LST in LST assimilation (clmupdate_T)
  real(r8),allocatable :: clm_paramarr(:)  !hcp CLM parameter vector (f.e. LAI)
  integer, allocatable :: state_clm2pdaf_p(:,:) !Index of column in hydraulic active state vector (nlevsoi,endc-begc+1)
  integer(c_int),bind(C,name="clmupdate_swc")     :: clmupdate_swc
  integer(c_int),bind(C,name="clmupdate_T")     :: clmupdate_T  ! by hcp
  integer(c_int),bind(C,name="clmupdate_texture") :: clmupdate_texture
  integer(c_int),bind(C,name="clmprint_swc")      :: clmprint_swc
#endif
  integer(c_int),bind(C,name="clmprint_et")       :: clmprint_et
  integer(c_int),bind(C,name="clmstatevec_allcol")       :: clmstatevec_allcol
  integer(c_int),bind(C,name="clmstatevec_colmean")       :: clmstatevec_colmean
  integer(c_int),bind(C,name="clmstatevec_only_active")  :: clmstatevec_only_active
  integer(c_int),bind(C,name="clmstatevec_max_layer")  :: clmstatevec_max_layer
  integer(c_int),bind(C,name="clmt_printensemble")       :: clmt_printensemble
  integer(c_int),bind(C,name="clmwatmin_switch")         :: clmwatmin_switch
  real(c_double),bind(C,name="clmcrns_bd")      :: clmcrns_bd

  ! Yorck
  integer(c_int),bind(C,name="clmupdate_tws") :: clmupdate_tws
  integer(c_int),bind(C,name="exclude_greenland") :: exclude_greenland
  real(r8),bind(C,name="da_interval") :: da_interval
  integer, dimension(1:5) :: clm_varsize_tws
  real(r8),bind(C,name="max_inc") :: max_inc
  integer(c_int),bind(C,name="TWS_smoother") :: TWS_smoother
  integer(c_int),bind(C,name="state_setup") :: state_setup !0: liq and ice seperated in statevector, 1: liq and ice together in statevector, 2: raw TWS values in statevector (just for testing)
  integer(c_int),bind(C,name="update_snow") :: update_snow !0: scripts from Lukas, 1: simple factor of old and new snow multiplied with old values
  integer(c_int),bind(C,name="remove_mean") :: remove_mean
  integer, allocatable :: num_layer(:)
  integer, allocatable :: num_layer_columns(:)

  real(r8), allocatable :: tws_temp_mean(:,:) ! temporal mean for TWS
  real(r8), allocatable :: lon_temp_mean(:,:) ! corresponding longitude
  real(r8), allocatable :: lat_temp_mean(:,:) ! corresponding latitude

  real(r8), allocatable :: tws_temp_mean_vector(:) ! temporal mean for TWS, in vector form, sorted just as sub-domain

  integer :: num_hactiveg, num_hactivec, num_hactiveg_patch, num_hactivep

  integer, allocatable :: hactiveg_levels(:,:)     ! hydrolocial active filter for all levels (gridcell) 
  integer, allocatable :: hactivec_levels(:,:)     ! hydrolocial active filter for all levels (column) 
  integer, allocatable :: hactivep(:)     ! hydrolocial active filter (patches)
  integer, allocatable :: hactiveg_patch(:)     ! hydrolocial active filter (patches)
  integer, allocatable :: gridcell_state(:)

  character(c_char),dimension(100),bind(C,name="mean_filename") :: mean_filename

  ! end Yorck

  integer  :: nstep     ! time step index
  real(r8) :: dtime     ! time step increment (sec)
  integer  :: ier       ! error code

  character(kind=c_char),dimension(100),bind(C,name="outdir"),target :: outdir

  logical  :: log_print    ! true=> print diagnostics
  real(r8) :: eccf         ! earth orbit eccentricity factor
  logical  :: mpi_running  ! true => MPI is initialized
  integer  :: mpicom_glob  ! MPI communicator

  character(len=SHR_KIND_CL) :: nlfilename = " "
  integer :: ierror, lengths_of_types, i
  logical :: flag
  integer(c_int),bind(C,name="clmprefixlen") :: clmprefixlen
  integer :: COMM_couple_clm    ! CLM-version of COMM_couple
                                ! (currently not used for clm5_0)
  logical :: newgridcell        !only clm5_0

  contains

#if defined CLMSA
  subroutine define_clm_statevec(mype)
    use shr_kind_mod, only: r8 => shr_kind_r8
    use decompMod , only : get_proc_bounds
    use clm_varpar   , only : nlevsoi
    use clm_varcon , only : ispval
    use clm_varcon, only: spval
    use GridcellType, only: grc
    use ColumnType , only : col
    use PatchType, only: patch

    implicit none

    integer,intent(in) :: mype

    integer :: i
    integer :: j
    integer :: jj
    integer :: c
    integer :: g
    integer :: p
    integer :: cg
    integer :: cc
    integer :: cccheck
    integer :: fa
    integer :: fg

    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices

    logical, allocatable :: found(:)

    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)

    lon   => grc%londeg
    lat   => grc%latdeg


    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

#ifdef PDAF_DEBUG
    WRITE(*,"(a,i5,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a,i10,a)") &
      "TSMP-PDAF mype(w)=", mype, " define_clm_statevec, CLM5-bounds (g,l,c,p)----",&
      begg,",",endg,",",begl,",",endl,",",begc,",",endc,",",begp,",",endp," -------"
#endif

    clm_begg     = begg
    clm_endg     = endg
    clm_begl     = begl
    clm_endl     = endl
    clm_begc     = begc
    clm_endc     = endc
    clm_begp     = begp
    clm_endp     = endp

    if (allocated(found)) deallocate(found)
    allocate(found(clm_begg:clm_endg))

    ! Soil Moisture DA: State vector index arrays
    if(clmupdate_swc.eq.1) then

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
      if(clmstatevec_allcol.eq.1) then

        ! Only hydrologically active columns
        if(clmstatevec_only_active .eq. 1) then

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
        if(clmstatevec_only_active.eq.1) then

          cc = 0

          do i=1,nlevsoi
            ! Only layers above max_layer
            if(i<=clmstatevec_max_layer) then

              do g=clm_begg,clm_endg

                newgridcell = .true.

                do c=clm_begc,clm_endc
                  if(col%gridcell(c) == g) then
                    ! All hydrologically active columns above
                    ! bedrock in a gridcell point to the state
                    ! vector index of the gridcell
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
      if(clmstatevec_only_active.eq.1) then
        ! Use iterator cc for setting state vector size.
        !
        ! Set `clm_varsize`, even though it is currently not used
        ! for `clmupdate_swc.eq.1`
        clm_varsize      =  cc
        clm_statevecsize =  cc
      else
        if(clmstatevec_allcol.eq.1) then
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

    endif

    if(clmupdate_swc.eq.2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    endif

    if(clmupdate_texture.eq.1) then
        clm_statevecsize = clm_statevecsize + 2*((endg-begg+1)*nlevsoi)
    endif

    if(clmupdate_texture.eq.2) then
        clm_statevecsize = clm_statevecsize + 3*((endg-begg+1)*nlevsoi)
    endif

    !hcp LST DA
    if(clmupdate_T.eq.1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    endif
    !end hcp

    if (clmupdate_tws.eq.1) then

      ! first we build a filter to determine which columns are active / are not active
      ! we also build a gridcell filter for gridcell averges
      num_hactiveg = 0
      num_hactivec = 0

      found(clm_begg:clm_endg) = .false.

      allocate(num_layer(1:nlevsoi))
      num_layer(1:nlevsoi) = 0

      allocate(num_layer_columns(1:nlevsoi))
      num_layer_columns(1:nlevsoi) = 0

      do c = clm_begc, clm_endc ! find out hydrological active cells

        g = col%gridcell(c) ! gridcell of column

        if ((exclude_greenland.eq.0) .or. (.not.(lon(g)<330 .and. lon(g)>180 .and. lat(g)>55))) then

          if (col%hydrologically_active(c)) then

            if (.not. found(g)) then ! if the gridcell is not found before

              found(g) = .true.

              do j = 1,nlevsoi
                ! get number in layers
    
                if (j<=col%nbedrock(c)) then
                  num_layer(j) = num_layer(j) + 1
                end if
    
              end do

              num_hactiveg = num_hactiveg + 1

            end if

            do j = 1,nlevsoi
              ! get number in layers

              if (j<=col%nbedrock(c)) then
                num_layer_columns(j) = num_layer_columns(j) + 1
              end if

            end do

            num_hactivec = num_hactivec + 1

          end if 
        end if

      end do


      found(clm_begg:clm_endg) = .false.
      num_hactiveg_patch = 0
      num_hactivep = 0
      do p = clm_begp, clm_endp
        c = patch%column(p)
        g = col%gridcell(c)

        if ((exclude_greenland.eq.0) .or. (.not.(lon(g)<330 .and. lon(g)>180 .and. lat(g)>55))) then

          if (col%hydrologically_active(c) .and. patch%active(p)) then
            if (.not. found(g)) then ! if the gridcell is not found before

              found(g) = .true.

              num_hactiveg_patch = num_hactiveg_patch+1

            end if

            num_hactivep = num_hactivep + 1

          end if

        end if

      end do

      allocate(hactiveg_levels(1:num_hactiveg,1:nlevsoi))
      allocate(hactivec_levels(1:num_hactivec,1:nlevsoi))
      allocate(hactiveg_patch(1:num_hactiveg_patch))
      allocate(hactivep(1:num_hactivep))

      ! now we fill these things with the columns and gridcells so that we can access all active things later on
      do j = 1,nlevsoi
        found(clm_begg:clm_endg) = .false. ! has to be inside the for lopp, else, the hactiveg_levels is only filled for the first level
        fa = 0
        fg = 0
        do c = clm_begc, clm_endc
          
          g = col%gridcell(c) ! gridcell of column

          if ((exclude_greenland.eq.0) .or. (.not.(lon(g)<330 .and. lon(g)>180 .and. lat(g)>55))) then

            if (col%hydrologically_active(c)) then

              if (.not. found(g)) then ! if the gridcell is not found before

                found(g) = .true.

                if (j<=col%nbedrock(c)) then
                  fg = fg+1
                  hactiveg_levels(fg,j) = g
                end if

              end if

              if (j<=col%nbedrock(c)) then
                fa = fa + 1
                hactivec_levels(fa,j) = c
              end if
              
            end if

          end if

        end do
      end do

      found(clm_begg:clm_endg) = .false.
      fa = 0
      fg = 0
      do p = clm_begp, clm_endp
        c = patch%column(p)
        g = col%gridcell(c) ! gridcell of column

        if ((exclude_greenland.eq.0) .or. (.not.(lon(g)<330 .and. lon(g)>180 .and. lat(g)>55))) then
          if (col%hydrologically_active(c) .and. patch%active(p)) then

            if (.not. found(g)) then ! if the gridcell is not found before
              found(g) = .true.
              fg = fg+1
              hactiveg_patch(fg) = g
            end if

            fa = fa+1
            hactivep(fa) = p

          end if
        end if

      end do

      if (allocated(found)) deallocate(found)

      ! now lets find out the dimension of the state vector

      ! first h2osoi_liq and h2osoi_ice
      clm_varsize_tws(:) = 0

      clm_statevecsize = 0

      select case (state_setup)
      case(0)
        do j = 1,nlevsoi
          clm_varsize_tws(1) = clm_varsize_tws(1) + num_layer(j)
          clm_statevecsize = clm_statevecsize + num_layer(j)

          clm_varsize_tws(2) = clm_varsize_tws(2) + num_layer(j)
          clm_statevecsize = clm_statevecsize + num_layer(j)
        end do

        ! snow
        clm_varsize_tws(3) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)

        ! surface water
        clm_varsize_tws(4) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)

        ! canopy water
        clm_varsize_tws(5) = num_hactiveg_patch
        clm_statevecsize = clm_statevecsize + num_hactiveg_patch

      case(1)

        do j = 1,nlevsoi
          clm_varsize_tws(1) = clm_varsize_tws(1) + num_layer(j)
          clm_statevecsize = clm_statevecsize + num_layer(j)

          clm_varsize_tws(2) = 0
          clm_statevecsize = clm_statevecsize + 0
        end do

        ! snow
        clm_varsize_tws(3) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)

        ! surface water
        clm_varsize_tws(4) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)

        ! canopy water
        clm_varsize_tws(5) = num_hactiveg_patch
        clm_statevecsize = clm_statevecsize + num_hactiveg_patch

      case(2)

        clm_varsize_tws(1) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)
        clm_varsize_tws(2) = 0
        clm_varsize_tws(3) = 0
        clm_varsize_tws(4) = 0
        clm_varsize_tws(5) = 0

      case(3) ! only sum over all soil layers and snow in state vector, maybe I will add other compartments too

        clm_varsize_tws(1) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)
        clm_varsize_tws(2) = 0
        
        ! snow
        clm_varsize_tws(3) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)
        clm_varsize_tws(4) = 0
        clm_varsize_tws(5) = 0

      case(4) ! sum over upper layer (1-7), sum over bottom layers (8-nlevsoi), snow

        clm_varsize_tws(1) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)
        clm_varsize_tws(2) = num_layer(8)
        clm_statevecsize = clm_statevecsize + num_layer(8)
        
        ! snow
        clm_varsize_tws(3) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)
        clm_varsize_tws(4) = 0
        clm_varsize_tws(5) = 0

      case(5) ! one variable for surface soil moisture (upper 10 cm), one for root zone soil moisture (until 200 cm), one for everything underneath and one for snow

        ! upper three layers for surface soil moisture (layer three is in a depth of 9 cm)

        clm_varsize_tws(1) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)

        ! layer 4 to 12 for root zone soil moisture (layer 12 is in a depth of 208 cm)

        clm_varsize_tws(2) = num_layer(4)
        clm_statevecsize = clm_statevecsize + num_layer(4)

        ! layer 13 to 20 for deep soil moisture

        clm_varsize_tws(3) = num_layer(13)
        clm_statevecsize = clm_statevecsize + num_layer(13)

        ! one variable for snow --> caution !!! normally snow is var 3, now it is var 4!!!

        clm_varsize_tws(4) = num_layer(1)
        clm_statevecsize = clm_statevecsize + num_layer(1)
        
      end select
      

    end if

#ifdef PDAF_DEBUG
    ! Debug output of clm_statevecsize
    WRITE(*, '(a,x,a,i5,x,a,i10)') "TSMP-PDAF-debug", "mype(w)=", mype, "define_clm_statevec: clm_statevecsize=", clm_statevecsize
#endif

    !write(*,*) 'clm_statevecsize is ',clm_statevecsize
    IF (allocated(clm_statevec)) deallocate(clm_statevec)
    if ((clmupdate_swc.ne.0) .or. (clmupdate_T.ne.0) .or. (clmupdate_texture.ne.0) .or. (clmupdate_tws.eq.1)) then
      !hcp added condition
      allocate(clm_statevec(clm_statevecsize))
    end if

    ! Allocate statevector-duplicate for saving original column mean
    ! values used in computing increments during updating the state
    ! vector in column-mean-mode.
    IF (allocated(clm_statevec_orig)) deallocate(clm_statevec_orig)
    if (clmupdate_swc.ne.0 .and. clmstatevec_colmean.ne.0) then
      allocate(clm_statevec_orig(clm_statevecsize))
    end if

    if (clmupdate_tws.eq.1) then
      IF (allocated(clm_statevec_original_input)) deallocate(clm_statevec_original_input)
      allocate(clm_statevec_original_input(1:clm_statevecsize))

      IF (allocated(gridcell_state)) deallocate(gridcell_state)
      allocate(gridcell_state(1:clm_statevecsize))
    end if

    !write(*,*) 'clm_paramsize is ',clm_paramsize
    if (allocated(clm_paramarr)) deallocate(clm_paramarr)         !hcp
    if ((clmupdate_T.ne.0)) then  !hcp
      error stop "Not implemented clmupdate_T.NE.0"
    end if

  end subroutine define_clm_statevec

  subroutine cleanup_clm_statevec()

    implicit none

    ! Deallocate arrays from `define_clm_statevec`
    IF (allocated(clm_statevec)) deallocate(clm_statevec)
    IF (allocated(state_pdaf2clm_c_p)) deallocate(state_pdaf2clm_c_p)
    IF (allocated(state_pdaf2clm_j_p)) deallocate(state_pdaf2clm_j_p)
    IF (allocated(state_clm2pdaf_p)) deallocate(state_clm2pdaf_p)
    IF (allocated(clm_statevec_original_input)) deallocate(clm_statevec_original_input)
    IF (allocated(gridcell_state)) deallocate(gridcell_state)

  end subroutine cleanup_clm_statevec

  subroutine set_clm_statevec(tstartcycle, mype)
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_varpar   , only : nlevsoi
    ! use clm_varcon, only: nameg, namec
    ! use GetGlobalValuesMod, only: GetGlobalWrite
    use ColumnType , only : col
    use PatchType, only: patch
    use GridcellType, only: grc
    use clm_varcon, only: spval
    use clm_varctl, only: inst_suffix
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    integer,intent(in) :: tstartcycle
    integer,intent(in) :: mype
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)
    real(r8), pointer :: h2osoi_liq(:,:) ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:) ! ice lens (kg/m2)
    real(r8), pointer :: h2osno(:) ! snow water (mm)
    real(r8), pointer :: h2osfc(:) ! surface water
    real(r8), pointer :: h2ocan(:) ! canopy water
    real(r8), pointer :: TWS(:)

    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)

    real(r8), pointer :: tws_state(:)
    real(r8), pointer :: h2osoi_liq_state(:,:)
    real(r8), pointer :: h2osoi_ice_state(:,:)
    real(r8), pointer :: h2osno_state(:)

    real(r8), pointer :: watsat(:,:)

    integer :: i,j,jj,g,c,cc=0,offset=0
    integer :: n_c
    character (len = 34) :: fn    !TSMP-PDAF: function name for state vector output
    character (len = 34) :: fn2    !TSMP-PDAF: function name for swc output

    integer :: count,count_columns, count_patch, p, l, k
    real(r8) :: avg_sum
    real(r8) :: avg_sum_ice
    real(r8) :: avg_sum_patch
    integer :: avg_divide
    integer :: avg_divide_patch

    character (len = 110) :: filename_temp
    
    swc   => waterstate_inst%h2osoi_vol_col
    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

    TWS => waterstate_inst%tws_hactive

    tws_state => waterstate_inst%tws_state_before
    h2osoi_liq_state => waterstate_inst%h2osoi_liq_state_before
    h2osoi_ice_state => waterstate_inst%h2osoi_ice_state_before
    h2osno_state => waterstate_inst%h2osno_state_before

    lon   => grc%londeg
    lat   => grc%latdeg

    watsat => soilstate_inst%watsat_col

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle + 1 .OR. clmt_printensemble < 0) THEN

      IF(clmupdate_swc.NE.0) THEN
        ! TSMP-PDAF: Debug output of CLM swc
        WRITE(fn2, "(a,i5.5,a,i5.5,a)") "swcstate_", mype, ".integrate.", tstartcycle + 1, ".txt"
        OPEN(unit=71, file=fn2, action="write")
        WRITE (71,"(es22.15)") swc(:,:)
        CLOSE(71)
      END IF

    END IF
#endif

    select case (TWS_smoother)
    case(0)

      !print*, 'instanteneous values in statevector'

      h2osoi_liq => waterstate_inst%h2osoi_liq_col 
      h2osoi_ice => waterstate_inst%h2osoi_ice_col
      h2osno => waterstate_inst%h2osno_col 
      h2osfc => waterstate_inst%h2osfc_col
      h2ocan => waterstate_inst%h2ocan_patch
      TWS => waterstate_inst%tws_hactive

    case default

      !print*, 'mean values over one month in statevector'

      h2osoi_liq => waterstate_inst%h2osoi_liq_col_mean 
      h2osoi_ice => waterstate_inst%h2osoi_ice_col_mean
      h2osno => waterstate_inst%h2osno_col_mean  
      h2osfc => waterstate_inst%h2osfc_col_mean
      h2ocan => waterstate_inst%h2ocan_patch_mean
      TWS => waterstate_inst%tws_hactive_mean

    end select

    ! calculate shift when CRP data are assimilated
    if(clmupdate_swc.eq.2) then
      error stop "Not implemented clmupdate_swc.eq.2"
    endif

    if(clmupdate_swc.ne.0) then
      ! write swc values to state vector
      if (clmstatevec_colmean.eq.1) then

        do cc = 1, clm_statevecsize

          clm_statevec(cc) = 0.0
          n_c = 0

          ! Get gridcell and layer
          g = col%gridcell(state_pdaf2clm_c_p(cc))
          j = state_pdaf2clm_j_p(cc)

          ! Loop over all columns
          do c=clm_begc,clm_endc
            ! Select columns in gridcell g
            if(col%gridcell(c).eq.g) then
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
    endif

    !hcp  LAI
    if(clmupdate_T.eq.1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    endif
    !end hcp  LAI

    ! write average swc to state vector (CRP assimilation)
    if(clmupdate_swc.eq.2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    endif

    ! write texture values to state vector (if desired)
    if(clmupdate_texture.ne.0) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
          clm_statevec(cc+1*clm_varsize+offset) = psand(j,i)
          clm_statevec(cc+2*clm_varsize+offset) = pclay(j,i)
          if(clmupdate_texture.eq.2) then
            !incl. organic matter values
            clm_statevec(cc+3*clm_varsize+offset) = porgm(j,i)
          end if
          cc = cc + 1
        end do
      end do
    endif


    if (clmupdate_tws.eq.1) then

      if (remove_mean.eq.1) then

        if (.not. allocated(tws_temp_mean_vector)) then

          do j = 1,100
            filename_temp(j:j) = mean_filename(j)
          end do
          filename_temp = trim(filename_temp)
          call read_temp_mean_model(filename_temp)


          if (allocated(tws_temp_mean_vector)) DEALLOCATE(tws_temp_mean_vector)
          ALLOCATE(tws_temp_mean_vector(clm_begg:clm_endg))
          tws_temp_mean_vector(:) = spval

          !this process only need the sub domain information
          do j = clm_begg,clm_endg
              ! find lon and lat in the file that corresponds to that of the grid point of the sub process
              outer: do l = 1,size(lon_temp_mean,1)
                do k=1,size(lon_temp_mean,2)
                    if (lon_temp_mean(l,k).eq.lon(j) .and. lat_temp_mean(l,k).eq.lat(j)) then
                      tws_temp_mean_vector(j) = tws_temp_mean(l,k)
                      exit outer
                    end if
                end do
              end do outer

              if (lon(j).ne.lon_temp_mean(l,k) .or. lat(j).ne.lat_temp_mean(l,k)) then
                print *, "Attention: distributing model mean to clumps does not work properly"
                print *, "idx_lon= ",l, "idx_lat= ",k
                print *, "lon(j)= ", lon(j),"lon_temp_mean(idx_lon)= ",lon_temp_mean(l,k)
                print *, "lat(j)= ", lat(j),"lat_temp_mean(idx_lat)= ",lat_temp_mean(l,k)
                stop
              end if
          end do

          deallocate(lon_temp_mean)
          deallocate(lat_temp_mean)
          deallocate(tws_temp_mean)

        end if

      end if

      select case (state_setup)

      case(0) ! all compartments, liq and ice water indidually

        if (inst_suffix=='_0000' .and. clm_begc==1) then
          print*, "Filling up state vector with all compartments, liq and ice water indidually"
        end if

        cc = 1

        do j = 1,nlevsoi 

          do count = 1, num_layer(j)

            g = hactiveg_levels(count,j)

            avg_sum = 0
            avg_sum_ice = 0
            avg_divide = 0

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osoi_liq(c,j)
                avg_sum_ice = avg_sum_ice + h2osoi_ice(c,j)

                avg_divide = avg_divide+1


              end if

            end do

            clm_statevec(cc) = avg_sum/avg_divide
            clm_statevec(cc+clm_varsize_tws(1)) = avg_sum_ice/avg_divide

            clm_statevec_original_input(cc) = avg_sum/avg_divide
            clm_statevec_original_input(cc+clm_varsize_tws(1)) = avg_sum_ice/avg_divide

            gridcell_state(cc) = g
            gridcell_state(cc+clm_varsize_tws(1)) = g

            h2osoi_liq_state(g,j) = clm_statevec(cc)
            h2osoi_ice_state(g,j) = clm_statevec(cc+clm_varsize_tws(1))

            avg_sum = 0
            avg_divide = 0
            if (j==1) then
              ! snow
              avg_sum = 0
              avg_divide = 0
              do count_columns = 1,num_layer_columns(j)
                c = hactivec_levels(count_columns,j)

                if (g==col%gridcell(c)) then

                  avg_sum = avg_sum + h2osno(c)

                  avg_divide = avg_divide+1


                end if

              end do

              clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = avg_sum/avg_divide
              clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = avg_sum/avg_divide

              gridcell_state(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = g

              h2osno_state(g) = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2))

              ! surface water
              avg_sum = 0
              avg_divide = 0
              do count_columns = 1,num_layer_columns(j)
                c = hactivec_levels(count_columns,j)

                if (g==col%gridcell(c)) then

                  avg_sum = avg_sum + h2osfc(c)

                  avg_divide = avg_divide+1


                end if

              end do

              clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) = avg_sum/avg_divide
              clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) = avg_sum/avg_divide

              gridcell_state(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) = g

            end if

            cc = cc+1

          end do
        end do



      case(1) ! all compartments, sum of ice and liq soil water to overcome balancing errors due to different partitioning of water caused by different temperature

        if (inst_suffix=='_0000' .and. clm_begc==1) then
          print*, "Filling up state vector with all compartments, sum of liq and ice"
        end if
        
        cc = 1

        do j = 1,nlevsoi 

          do count = 1, num_layer(j)

            g = hactiveg_levels(count,j)

            avg_sum = 0
            avg_sum_ice = 0
            avg_divide = 0

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

                avg_divide = avg_divide+1

              end if

            end do

            clm_statevec(cc) = avg_sum/avg_divide
            clm_statevec_original_input(cc) = avg_sum/avg_divide

            gridcell_state(cc) = g

            h2osoi_liq_state(g,j) = clm_statevec(cc)

            avg_sum = 0
            avg_divide = 0
            if (j==1) then
              ! snow
              avg_sum = 0
              avg_divide = 0
              do count_columns = 1,num_layer_columns(j)
                c = hactivec_levels(count_columns,j)

                if (g==col%gridcell(c)) then

                  avg_sum = avg_sum + h2osno(c)

                  avg_divide = avg_divide+1


                end if

              end do

              clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = avg_sum/avg_divide
              clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = avg_sum/avg_divide

              gridcell_state(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = g

              h2osno_state(g) = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2))

              ! surface water
              avg_sum = 0
              avg_divide = 0
              do count_columns = 1,num_layer_columns(j)
                c = hactivec_levels(count_columns,j)

                if (g==col%gridcell(c)) then

                  avg_sum = avg_sum + h2osfc(c)

                  avg_divide = avg_divide+1


                end if

              end do

              clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) = avg_sum/avg_divide
              clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) = avg_sum/avg_divide

              gridcell_state(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) = g

            end if

            cc = cc+1

          end do
        end do



      case(2) ! only TWS in statevector

        if (inst_suffix=='_0000' .and. clm_begc==1) then
          print*, "Filling up state vector with TWS"
        end if

        cc = 1
        do count = 1, num_layer(1)

          g = hactiveg_levels(count,1)

          if (remove_mean.eq.0) then

            clm_statevec(cc) = TWS(g)
            clm_statevec_original_input(cc) = TWS(g)

          else

            clm_statevec(cc) = TWS(g)-tws_temp_mean_vector(g)
            clm_statevec_original_input(cc) = TWS(g)-tws_temp_mean_vector(g)

          end if

          gridcell_state(cc) = g

          tws_state(g) = clm_statevec(cc)
          
          cc = cc+1

        end do

      case(3) ! sum over all soil layers and snow in statevector

        if (inst_suffix=='_0000' .and. clm_begc==1) then
          print*, "Filling up state vector with sum over all soil layers and snow"
        end if

        cc = 1

        do count = 1, num_layer(1)

          clm_statevec(cc) = 0

          g = hactiveg_levels(count,1)

          do j = 1, nlevsoi

            avg_sum = 0
            avg_divide = 0

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

                avg_divide = avg_divide+1

              end if

            end do

            if (avg_divide.ne.0) then

              clm_statevec(cc) = clm_statevec(cc) + avg_sum/avg_divide

            end if

          end do

          clm_statevec_original_input(cc) = clm_statevec(cc)
          h2osoi_liq_state(g,1) = clm_statevec(cc)

          cc = cc + 1

        end do


        ! snow

        cc = 1

        do count = 1, num_layer(1)

          g = hactiveg_levels(count,1)

          avg_sum = 0
          avg_divide = 0
          do count_columns = 1,num_layer_columns(1)
            c = hactivec_levels(count_columns,1)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osno(c)

              avg_divide = avg_divide+1


            end if

          end do

          clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = avg_sum/avg_divide
          clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = avg_sum/avg_divide

          h2osno_state(g) = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2))


          cc = cc+1


        end do


      case(4) ! sum of liq and ice water in the upper layers (1-7, until some gridcells have first bedrock layer in depth 8), sum underneath and snow in state vector

        if (inst_suffix=='_0000' .and. clm_begc==1) then
          print*, "Filling up state vector with sum in the upper layers, sum underneath and snow"
        end if

        cc = 1

        do count = 1, num_layer(1)

          clm_statevec(cc) = 0

          g = hactiveg_levels(count,1)

          do j = 1, 7

            avg_sum = 0
            avg_divide = 0

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

                avg_divide = avg_divide+1

              end if

            end do

            if (avg_divide.ne.0) then

              clm_statevec(cc) = clm_statevec(cc) + avg_sum/avg_divide

            end if

          end do

          clm_statevec_original_input(cc) = clm_statevec(cc)
          h2osoi_liq_state(g,1) = clm_statevec(cc)

          cc = cc + 1

        end do

        cc = 1

        do count = 1, num_layer(8)

          clm_statevec(cc+clm_varsize_tws(1)) = 0

          g = hactiveg_levels(count,8)

          do j = 8, nlevsoi

            avg_sum = 0
            avg_divide = 0

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

                avg_divide = avg_divide+1

              end if

            end do

            if (avg_divide.ne.0) then

              clm_statevec(cc+clm_varsize_tws(1)) = clm_statevec(cc+clm_varsize_tws(1)) + avg_sum/avg_divide

            end if

          end do

          clm_statevec_original_input(cc+clm_varsize_tws(1)) = clm_statevec(cc+clm_varsize_tws(1))
          h2osoi_liq_state(g,2) = clm_statevec(cc+clm_varsize_tws(1))

          cc = cc + 1

        end do

        ! snow

        cc = 1

        do count = 1, num_layer(1)

          g = hactiveg_levels(count,1)

          avg_sum = 0
          avg_divide = 0
          do count_columns = 1,num_layer_columns(1)
            c = hactivec_levels(count_columns,1)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osno(c)

              avg_divide = avg_divide+1


            end if

          end do

          clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = avg_sum/avg_divide
          clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = avg_sum/avg_divide

          h2osno_state(g) = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2))


          cc = cc+1


        end do

      case(5) ! sum of liq and ice water in surface soil mositure (1-3), root zone (4-12), sum underneath (13-20) and snow in state vector

        if (inst_suffix=='_0000' .and. clm_begc==1) then
          print*, "Filling up state vector with sum  over surface, root zone, 'groundwater', snow"
        end if

        cc = 1

        do count = 1, num_layer(1)

          clm_statevec(cc) = 0

          g = hactiveg_levels(count,1)

          do j = 1, 3

            avg_sum = 0
            avg_divide = 0

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

                avg_divide = avg_divide+1

              end if

            end do

            if (avg_divide.ne.0) then

              clm_statevec(cc) = clm_statevec(cc) + avg_sum/avg_divide

            end if

          end do

          clm_statevec_original_input(cc) = clm_statevec(cc)
          h2osoi_liq_state(g,1) = clm_statevec(cc)

          cc = cc + 1

        end do

        cc = 1

        do count = 1, num_layer(4)

          clm_statevec(cc+clm_varsize_tws(1)) = 0

          g = hactiveg_levels(count,4)

          do j = 4, 12

            avg_sum = 0
            avg_divide = 0

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

                avg_divide = avg_divide+1

              end if

            end do

            if (avg_divide.ne.0) then

              clm_statevec(cc+clm_varsize_tws(1)) = clm_statevec(cc+clm_varsize_tws(1)) + avg_sum/avg_divide

            end if

          end do

          clm_statevec_original_input(cc+clm_varsize_tws(1)) = clm_statevec(cc+clm_varsize_tws(1))
          h2osoi_liq_state(g,2) = clm_statevec(cc+clm_varsize_tws(1))

          cc = cc + 1

        end do

        cc = 1

        do count = 1, num_layer(13)

          clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = 0

          g = hactiveg_levels(count,13)

          do j = 13, nlevsoi

            avg_sum = 0
            avg_divide = 0

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

                avg_divide = avg_divide+1

              end if

            end do

            if (avg_divide.ne.0) then

              clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) + avg_sum/avg_divide

            end if

          end do

          clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2))
          h2osoi_liq_state(g,3) = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2))

          cc = cc + 1

        end do


        ! snow

        cc = 1

        do count = 1, num_layer(1)

          g = hactiveg_levels(count,1)

          avg_sum = 0
          avg_divide = 0
          do count_columns = 1,num_layer_columns(1)
            c = hactivec_levels(count_columns,1)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osno(c)

              avg_divide = avg_divide+1


            end if

          end do

          clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) = avg_sum/avg_divide
          clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) = avg_sum/avg_divide

          h2osno_state(g) = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3))


          cc = cc+1


        end do

      end select


    end if

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle + 1 .OR. clmt_printensemble < 0) THEN
      ! TSMP-PDAF: For debug runs, output the state vector in files
      WRITE(fn, "(a,i5.5,a,i5.5,a)") "clmstate_", mype, ".integrate.", tstartcycle + 1, ".txt"
      OPEN(unit=71, file=fn, action="write")
      DO i = 1, clm_statevecsize
        WRITE (71,"(es22.15)") clm_statevec(i)
      END DO
      CLOSE(71)
    END IF
#endif

  end subroutine set_clm_statevec

  subroutine update_clm(tstartcycle, mype) bind(C,name="update_clm")
    use clm_varpar   , only : nlevsoi
    use clm_time_manager  , only : update_DA_nstep
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
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    real(r8), pointer :: dz(:,:)          ! layer thickness depth (m)
    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)
    real(r8)  :: rliq,rice
    real(r8)  :: watmin_check      ! minimum soil moisture for checking clm_statevec (mm)
    real(r8)  :: watmin_set        ! minimum soil moisture for setting swc (mm)
    real(r8)  :: swc_update        ! updated SWC in loop

    integer :: i,j,jj,g,cc=0,offset=0
    character (len = 31) :: fn    !TSMP-PDAF: function name for state vector outpu
    character (len = 31) :: fn2    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn3    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn4    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn5    !TSMP-PDAF: function name for state vector outpu
    character (len = 32) :: fn6    !TSMP-PDAF: function name for state vector outpu

    logical :: swc_zero_before_update = .false.

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN
      ! TSMP-PDAF: For debug runs, output the state vector in files
      WRITE(fn, "(a,i5.5,a,i5.5,a)") "clmstate_", mype, ".update.", tstartcycle, ".txt"
      OPEN(unit=71, file=fn, action="write")
      DO i = 1, clm_statevecsize
        WRITE (71,"(es22.15)") clm_statevec(i)
      END DO
      CLOSE(71)
    END IF
#endif

    swc   => waterstate_inst%h2osoi_vol_col
    watsat => soilstate_inst%watsat_col
    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

    dz            => col%dz
    h2osoi_liq    => waterstate_inst%h2osoi_liq_col
    h2osoi_ice    => waterstate_inst%h2osoi_ice_col

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN

      IF(clmupdate_swc.NE.0) THEN
        ! TSMP-PDAF: For debug runs, output the state vector in files
        WRITE(fn5, "(a,i5.5,a,i5.5,a)") "h2osoi_liq", mype, ".bef_up.", tstartcycle, ".txt"
        OPEN(unit=71, file=fn5, action="write")
        WRITE (71,"(es22.15)") h2osoi_liq(:,:)
        CLOSE(71)

        ! TSMP-PDAF: For debug runs, output the state vector in files
        WRITE(fn6, "(a,i5.5,a,i5.5,a)") "h2osoi_ice", mype, ".bef_up.", tstartcycle, ".txt"
        OPEN(unit=71, file=fn6, action="write")
        WRITE (71,"(es22.15)") h2osoi_ice(:,:)
        CLOSE(71)
      END IF

    END IF
#endif

    ! calculate shift when CRP data are assimilated
    if(clmupdate_swc.eq.2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    endif

    ! CLM5: Update the Data Assimulation time-step to the current time
    ! step, since DA has been done. Used by CLM5 to skip BalanceChecks
    ! directly after the DA step.
    call update_DA_nstep()

    ! write updated swc back to CLM
    if(clmupdate_swc.ne.0) then

        ! Set minimum soil moisture for checking the state vector and
        ! for setting minimum swc for CLM
        if(clmwatmin_switch.eq.3) then
          ! CLM3.5 type watmin
          watmin_check = 0.00
          watmin_set = 0.05
        else if(clmwatmin_switch.eq.5) then
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

              ! Update only those SWCs that are not excluded by ispval
              if(state_clm2pdaf_p(j,i) .ne. ispval) then

                if(swc(j,i).eq.0.0) then
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

                if (clmstatevec_colmean.eq.1) then
                  ! Update SWC column value with the increment-factor
                  ! of the state vector update (state vector updates
                  ! are means of cols in grc)
                  swc_update = swc(j,i) * clm_statevec(state_clm2pdaf_p(j,i)) / clm_statevec_orig(state_clm2pdaf_p(j,i))
                else
                  ! Update SWC with updated state vector
                  swc_update = clm_statevec(state_clm2pdaf_p(j,i))
                end if

                if(swc_update.le.watmin_check) then
                  swc(j,i) = watmin_set
                else if(swc_update.ge.watsat(j,i)) then
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
              ! cc = cc + 1
            end do
        end do

#ifdef PDAF_DEBUG
        IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN

          IF(clmupdate_swc.NE.0) THEN
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

        END IF
#endif

    endif

    !hcp: TG, TV
    if(clmupdate_T.EQ.1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    endif
    ! end hcp TG, TV

    !! update liquid water content
    !do j=clm_begg,clm_endg
    !  do i=1,nlevsoi
    !    h2osoi_liq(j,i) = swc(j,i) * dz(j,i)*denh2o
    !  end do
    !end do

    ! write updated texture back to CLM
    if(clmupdate_texture.ne.0) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
          psand(j,i) = clm_statevec(cc+1*clm_varsize+offset)
          pclay(j,i) = clm_statevec(cc+2*clm_varsize+offset)
          if(clmupdate_texture.eq.2) then
            ! incl. organic matter
            porgm(j,i) = clm_statevec(cc+3*clm_varsize+offset)
          end if
          cc = cc + 1
        end do
      end do
      call clm_correct_texture
      call clm_texture_to_parameters
    endif

    if (clmupdate_tws.eq.1) then

      call clm_update_tws

    end if

  end subroutine update_clm

  subroutine clm_update_tws()
    use clm_instMod
    use clm_varpar   , only : nlevsoi, nlevsno
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varcon, only: spval, watmin, denh2o, denice, averaging_var
    use ColumnType         , only : col
    use LandunitType, only: lun
    use clm_varctl, only: inst_suffix
    implicit none

    integer, pointer :: snl(:) ! negative number of snow layers
    real(r8), pointer :: h2osno(:) ! snow water (mm)
    real(r8), pointer :: h2osoi_ice(:,:) ! ice lens (kg/m2)
    real(r8), pointer :: h2osoi_liq(:,:) ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_vol(:,:)

    real(r8), pointer :: TWS(:)

    real(r8), pointer :: h2osoi_liq_mean(:,:) ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice_mean(:,:) ! ice lens (kg/m2)
    real(r8), pointer :: h2osno_mean(:) ! snow water (mm)

    real(r8), pointer :: h2osoi_liq_inc(:,:) ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice_inc(:,:) ! ice lens (kg/m2)

    real(r8), pointer :: h2osno_inc(:) ! snow water (mm)
    real(r8), pointer :: snow_depth(:) ! snow water (mm)
    real(r8), pointer :: dz(:,:) ! snow water (mm)
    real(r8), pointer :: zi(:,:) ! snow water (mm)
    real(r8), pointer :: z(:,:) ! snow water (mm)

    real(r8), pointer :: watsat(:,:)

    real(r8), pointer :: forc_t(:)

    real(r8), pointer :: forc_wind(:)

    real(r8), pointer :: frac_iceold(:,:)

    real(r8), pointer :: tws_state(:)

    real(r8), pointer :: h2osoi_liq_state(:,:)
    real(r8), pointer :: h2osoi_ice_state(:,:)
    real(r8), pointer :: h2osno_state(:)

    

    ! Local variables:
    integer :: c, j, fc,cc, l,p,g, temp, count, count_columns                ! indices
    real(r8) :: inc, var_temp, inc_1, inc_2, inc_col, ratio, inc_ice !increment

    real(r8) :: scale

    real(r8) :: rsnow(clm_begc:clm_endc)
    real(r8) :: snowden, frac_swe, frac_liq, frac_ice
    real(r8) :: gain_h2osno, gain_h2oliq, gain_h2oice, gain_dzsno

    real(r8) :: t_for_bifall_degC  ! temperature to use in bifall equation (deg C)
    real(r8) :: bifall ! bulk density of newly fallen dry snow [kg/m3]


    real(r8) :: avg_sum
    real(r8) :: avg_sum_ice
    real(r8) :: avg_sum_patch
    integer :: avg_divide
    integer :: avg_divide_patch

    real(r8) :: mult_liq(clm_begc:clm_endc)
    real(r8) :: mult_ice(clm_begc:clm_endc)

    real(r8) :: lok_liq(clm_begc:clm_endc,1:nlevsoi)
    real(r8) :: lok_ice(clm_begc:clm_endc,1:nlevsoi)
    real(r8) :: lok_vol(clm_begc:clm_endc,1:nlevsoi)


    select case (TWS_smoother)
    case(0)
      h2osoi_liq_mean => waterstate_inst%h2osoi_liq_col 
      h2osoi_ice_mean => waterstate_inst%h2osoi_ice_col
      h2osno_mean => waterstate_inst%h2osno_col 
    case default
      h2osoi_liq_mean => waterstate_inst%h2osoi_liq_col_mean 
      h2osoi_ice_mean => waterstate_inst%h2osoi_ice_col_mean
      h2osno_mean => waterstate_inst%h2osno_col_mean     
    end select

    h2osoi_liq => waterstate_inst%h2osoi_liq_col 
    h2osoi_ice => waterstate_inst%h2osoi_ice_col 
    h2osno => waterstate_inst%h2osno_col
    snl => col%snl     
    dz         => col%dz
    zi         => col%zi
    z          => col%z
    snow_depth => waterstate_inst%snow_depth_col
    h2osoi_vol => waterstate_inst%h2osoi_vol_col
    watsat => soilstate_inst%watsat_col

    h2osoi_liq_inc => waterstate_inst%h2osoi_liq_col_inc
    h2osoi_ice_inc => waterstate_inst%h2osoi_ice_col_inc
    h2osno_inc => waterstate_inst%h2osno_col_inc

    forc_t      => atm2lnd_inst%forc_t_downscaled_col !  atmospheric temperature (Kelvin)
    forc_wind   => atm2lnd_inst%forc_wind_grc         !  atmospheric wind speed (m/s)
    frac_iceold =>  waterstate_inst%frac_iceold_col   !  fraction of ice relative to the tot water

    TWS => waterstate_inst%tws_hactive

    tws_state => waterstate_inst%tws_state_after

    h2osoi_liq_state => waterstate_inst%h2osoi_liq_state_after
    h2osoi_ice_state => waterstate_inst%h2osoi_ice_state_after
    h2osno_state => waterstate_inst%h2osno_state_after


     ! now all variables are updated. Restrictions have to be introduced to ensure that the model is still running correctly 

    ! set averaging factor to zero
    averaging_var = 0
    
    do j = 1,nlevsoi
      do count = 1,num_layer_columns(j)
        c = hactivec_levels(count,j)

        h2osoi_liq_inc(c,j) = h2osoi_liq(c,j)
        h2osoi_ice_inc(c,j) = h2osoi_ice(c,j)

        if (j==1) then
          h2osno_inc(c) = h2osno(c)
        end if
      end do
    end do



    cc = 1

    if (state_setup.eq.0 .or. state_setup.eq.1) then

      do j = 1,nlevsoi

        do count = 1, num_layer(j)

          g = hactiveg_levels(count,j)

          inc = clm_statevec(cc)-clm_statevec_original_input(cc)

          if (abs(inc)>1.e-10_r8) then

            if (state_setup.eq.0) then
              inc_ice = clm_statevec(cc+clm_varsize_tws(1))
            end if

            do count_columns = 1, num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (col%gridcell(c)==g) then

                select case(state_setup)
                case(0)

                  inc_col = inc

                  if (inc_col/=inc_col) then
                    inc_col = 0.0
                  end if

                  ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                  ! so that the direction of the increment is right
                  if (abs(inc_col).gt.max_inc*h2osoi_liq(c,j)) then
                    inc_col = sign(max_inc*h2osoi_liq(c,j),inc_col)
                  end if

                  h2osoi_liq(c,j) = h2osoi_liq(c,j) + inc_col

                  if (h2osoi_liq(c,j).lt.watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  inc_col = inc_ice-clm_statevec_original_input(cc+clm_varsize_tws(1))

                  if (inc_col/=inc_col) then
                    inc_col = 0.0
                  end if

                  if (abs(inc_col).gt.max_inc*h2osoi_ice(c,j)) then
                    inc_col = sign(max_inc*h2osoi_ice(c,j),inc_col)
                  end if

                  h2osoi_ice(c,j) = h2osoi_ice(c,j) + inc_col

                  if (h2osoi_ice(c,j).lt.0) then
                    h2osoi_ice(c,j) = 0._r8
                  end if

                case(1)

                  inc_col = inc

                  var_temp = h2osoi_liq(c,j)+h2osoi_ice(c,j)

                  if (abs(inc_col).gt.max_inc*var_temp) then
                    inc_col = sign(max_inc*var_temp,inc_col)
                  end if

                  inc_1 = inc_col*(h2osoi_liq(c,j)/var_temp)

                  h2osoi_liq(c,j) = h2osoi_liq(c,j)+inc_1

                  if (h2osoi_liq(c,j)<watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  inc_2 = inc_col*(h2osoi_ice(c,j)/var_temp)

                  h2osoi_ice(c,j) = h2osoi_ice(c,j)+inc_2

                  if (h2osoi_ice(c,j) < 0) then
                    h2osoi_ice(c,j) = 0._r8
                  end if

                end select

                h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)

                if (h2osoi_vol(c,j)-watsat(c,j)>0.00001 .and. j>1) then

                  if (h2osoi_ice(c,j) == 0) then

                    h2osoi_liq(c,j) = h2osoi_vol(c,j)*dz(c,j)*denh2o

                  else

                    var_temp = watsat(c,j) / h2osoi_vol(c,j)

                    h2osoi_liq(c,j) = h2osoi_liq(c,j)*var_temp

                    if (h2osoi_liq(c,j) < watmin) then
                      h2osoi_liq(c,j) = watmin
                    end if

                    h2osoi_ice(c,j) = (watsat(c,j)*dz(c,j)*denice)-h2osoi_liq(c,j)*denice/denh2o

                  end if

                  h2osoi_vol(c,j) = watsat(c,j)

                end if

              end if

            end do

          end if

          cc = cc + 1

        end do

      end do


      ! update snow

      cc = 1

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        inc = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) - clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) ! save increment for gridcell

        if (abs(inc)>1.e-10_r8) then

          do count_columns = 1, num_layer_columns(1)

            c = hactivec_levels(count_columns,1)

            if (col%gridcell(c)==g) then

              !ratio = h2osno_mean(c)/clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) ! ratio of column value to averaged gridcell value (input in state vector)
              !inc_col = inc*ratio - h2osno_mean(c) ! increment of column is then the gridcell increment scaled with this ratio
              inc_col = inc

              if (inc_col/=inc_col) then
                inc_col = 0.0
              end if

              if (abs(inc_col)>500) then
                inc_col = sign(500._r8,inc_col)
              end if

              select case (update_snow)
              ! Tests with snow DA, scripts adapted from Lukas Strebel
              case(0)

                if (inc_col.ne.0._r8) then

                  if (snl(c) < 0) then ! snow layers in the column

                    h2osno(c) = h2osno(c) + inc_col

                    if (h2osno(c)>10000) then
                      h2osno(c) = 10000
                    end if

                    do j=0,snl(c)+1,-1 ! iterate through the snow layers

                      ! snow density prior for each layer
                      if (dz(c,j)>0.0_r8) then
                        snowden = (h2osoi_liq(c,j) + h2osoi_ice(c,j)) / dz(c,j)
                      else
                        snowden = 0.0_r8
                      endif

                      ! fraction of SWE in each active layers
                      if(rsnow(c).gt.0.0_r8) then
                        frac_swe = (h2osoi_liq(c,j) + h2osoi_ice(c,j)) / rsnow(c)
                      else
                        frac_swe = 0.0_r8 ! no fraction SWE if no snow is present in column
                      end if ! end SWE fraction if

                      ! fraction of liquid and ice
                      if ((h2osoi_liq(c,j) + h2osoi_ice(c,j)).gt.0.0_r8) then
                        frac_liq = h2osoi_liq(c,j) / (h2osoi_liq(c,j) + h2osoi_ice(c,j))
                        frac_ice = 1.0_r8 - frac_liq
                      else
                        frac_liq = 0.0_r8
                        frac_ice = 0.0_r8
                      end if

                      ! SWE adjustment per layer 
                      ! assumes identical layer distribution of liq and ice than before DA (frac_*)
                      gain_h2osno = (h2osno(c) - rsnow(c)) * frac_swe
                      gain_h2oliq = gain_h2osno * frac_liq
                      gain_h2oice = gain_h2osno * frac_ice

                      ! layer level adjustments
                      if (snowden.gt.0.0_r8) then
                        gain_dzsno = gain_h2osno / snowden
                      else
                        gain_dzsno = 0.0_r8
                      end if
                      h2osoi_liq(c,j) = h2osoi_liq(c,j) + gain_h2oliq
                      h2osoi_ice(c,j) = h2osoi_ice(c,j) + gain_h2oice


                      ! Adjust snow layer dimensions so that CLM5 can calculate compaction / aggregation
                      ! in the DART code dzsno is adjusted directly but in CLM5 dzsno is local and diagnostic
                      ! i.e. calculated / assigned from frac_sno and dz(:, snow_layer) in SnowHydrologyMod
                      ! therefore we adjust dz(:, snow_layer) here

                      dz(c,j) = dz(c,j) + gain_dzsno
                      ! mid point and interface adjustments
                      ! i.e. zsno (col%z(:, snow_layers)) and zisno (col%zi(:, snow_layers))
                      ! DART version the sum goes from ilevel:nlevsno to fit with our indexing:
                      zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8
                      ! In DART the check is ilevel == nlevsno but here
                      
                      if (j.eq.0) then
                        z(c,j) = zi(c,j-1) / 2.0_r8
                      else
                        z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
                      end if


                    end do

                    ! Update the total snow depth to match updates to layers for active snow layers                
                    snow_depth(c) = sum(dz(c,snl(c)+1:0))
                    h2osno(c) = sum(h2osoi_ice(c,snl(c)+1:0)+h2osoi_liq(c,snl(c)+1:0))
                    
                  end if

                end if

              case (1) !update with factor of old and new snow

                if (snl(c) < 0) then ! snow layers in the column

                  ! if (inc_col>1000._r8) then
                  !   inc_col = 1000._r8
                  ! end if

                  inc_col = h2osno(c)+inc_col

                  if (inc_col>10000) then
                    inc_col = 10000
                  end if

                  scale = inc_col/h2osno(c)
                  h2osno(c) = inc_col
                  
                  do j=0,snl(c)+1,-1
                    h2osoi_liq(c,j) = h2osoi_liq(c,j)*scale
                    h2osoi_ice(c,j) = h2osoi_ice(c,j)*scale
                    dz(c,j) = dz(c,j)*scale
                    zi(c,j) = zi(c,j)*scale
                    z(c,j) = z(c,j)*scale
                  end do
                  zi(c,snl(c)) = zi(c,snl(c))*scale
                  snow_depth(c) = snow_depth(c)*scale

                end if

              end select

              ! snow negative
              if (h2osno(c) < 0._r8) then
                if (snl(c)<0) then

                  do j=0,snl(c)+1,-1
                      h2osoi_liq(c,j) = 0.0_r8
                      h2osoi_ice(c,j) = 0.00000001_r8
                      dz(c,j)  = 0.00000001_r8  
                      zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8                 
                      if (j.eq.0) then
                        z(c,j) = zi(c,j-1) / 2.0_r8
                      else
                        z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
                      end if                 
                  end do

                else

                  h2osoi_liq(c,0) = 0.0_r8
                  h2osoi_ice(c,0) = 0.00000001_r8
                  dz(c,0)  = 0.00000001_r8
                  zi(c,-1) = dz(c,0)*-1.0_r8 
                  z(c,0) = zi(c,-1) / 2.0_r8

                end if

                snow_depth(c) = sum(dz(c,-nlevsno+1:0))
                h2osno(c) = sum(h2osoi_ice(c,-nlevsno+1:0))
              end if

            end if

          end do

        end if

        cc = cc+1

      end do




    elseif (state_setup.eq.2) then

      cc = 1

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        inc = clm_statevec(cc)-clm_statevec_original_input(cc)

        if (abs(inc)>1.e-10_r8) then

          ! update soil water
          do j = 1,nlevsoi

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (col%gridcell(c)==g) then

                inc_col = inc*h2osoi_liq_mean(c,j)/clm_statevec_original_input(cc)

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_liq(c,j)) then
                  inc_col = sign(max_inc*h2osoi_liq(c,j),inc_col)
                end if

                h2osoi_liq(c,j) = h2osoi_liq(c,j) + inc_col

                if (h2osoi_liq(c,j).lt.watmin) then
                  h2osoi_liq(c,j) = watmin
                end if




                inc_col = inc*h2osoi_ice_mean(c,j)/clm_statevec_original_input(cc)

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_ice(c,j)) then
                  inc_col = sign(max_inc*h2osoi_ice(c,j),inc_col)
                end if

                h2osoi_ice(c,j) = h2osoi_ice(c,j) + inc_col

                if (h2osoi_ice(c,j).lt.0) then
                  h2osoi_ice(c,j) = 0._r8
                end if



                h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)

                if (j>1 .and. h2osoi_vol(c,j)-watsat(c,j)>0.000001) then

                  var_temp = watsat(c,j) / h2osoi_vol(c,j)

                  h2osoi_liq(c,j) = h2osoi_liq(c,j)*var_temp

                  if (h2osoi_liq(c,j) < watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  h2osoi_ice(c,j) = (watsat(c,j)*dz(c,j)*denice)-h2osoi_liq(c,j)*denice/denh2o
                  if (abs(h2osoi_ice(c,j))<1.e-10_r8) then !numerics, if h2osoiice was zero before, it is -something e-16 after the previous calculation, so something marginal negative
                    h2osoi_ice(c,j) = 0._r8
                  end if

                  h2osoi_vol(c,j) = watsat(c,j)

                end if

              end if

            end do

          end do

          ! update snow

          do count_columns = 1,num_layer_columns(1)

            c = hactivec_levels(count_columns,1)

            if (col%gridcell(c)==g) then

              inc_col = inc*h2osno_mean(c)/clm_statevec_original_input(cc)

              if (abs(inc_col)>500) then
                inc_col = sign(500._r8,inc_col)
              end if

              if (inc_col/=inc_col) then
                inc_col = 0.0
              end if

              if (snl(c) < 0) then ! snow layers in the column

                inc_col = h2osno(c)+inc_col

                if (inc_col>10000) then
                  inc_col = 10000
                end if

                scale = inc_col/h2osno(c)
                h2osno(c) = inc_col
                
                do j=0,snl(c)+1,-1
                  h2osoi_liq(c,j) = h2osoi_liq(c,j)*scale
                  h2osoi_ice(c,j) = h2osoi_ice(c,j)*scale
                  dz(c,j) = dz(c,j)*scale
                  zi(c,j) = zi(c,j)*scale
                  z(c,j) = z(c,j)*scale
                end do
                zi(c,snl(c)) = zi(c,snl(c))*scale
                snow_depth(c) = snow_depth(c)*scale

              end if

              ! snow negative
              if (h2osno(c) < 0._r8) then
                if (snl(c)<0) then

                  do j=0,snl(c)+1,-1
                      h2osoi_liq(c,j) = 0.0_r8
                      h2osoi_ice(c,j) = 0.00000001_r8
                      dz(c,j)  = 0.00000001_r8  
                      zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8                 
                      if (j.eq.0) then
                        z(c,j) = zi(c,j-1) / 2.0_r8
                      else
                        z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
                      end if                 
                  end do

                else

                  h2osoi_liq(c,0) = 0.0_r8
                  h2osoi_ice(c,0) = 0.00000001_r8
                  dz(c,0)  = 0.00000001_r8
                  zi(c,-1) = dz(c,0)*-1.0_r8 
                  z(c,0) = zi(c,-1) / 2.0_r8

                end if

                snow_depth(c) = sum(dz(c,-nlevsno+1:0))
                h2osno(c) = sum(h2osoi_ice(c,-nlevsno+1:0))
              end if


            end if

          end do

        end if
        cc = cc+1

      end do




    elseif (state_setup.eq.3) then

      ! ! soil water
      cc = 1

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        inc = clm_statevec(cc)-clm_statevec_original_input(cc)

        if (inc/=inc) then
          inc = 0.0
        end if

        if (abs(inc)>1.e-10_r8) then

          ! update soil water
          do j = 1,nlevsoi

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (col%gridcell(c)==g) then

                var_temp = h2osoi_liq(c,j)+h2osoi_ice(c,j)

                inc_col = inc*(h2osoi_liq_mean(c,j)+h2osoi_ice_mean(c,j))/clm_statevec_original_input(cc)

                if (abs(inc_col).gt.max_inc*var_temp) then
                  inc_col = sign(max_inc*var_temp,inc_col)
                end if

                h2osoi_liq(c,j) = h2osoi_liq(c,j) + inc_col*(h2osoi_liq(c,j)/var_temp)
                h2osoi_ice(c,j) = h2osoi_ice(c,j) + inc_col*(h2osoi_ice(c,j)/var_temp)

                if (h2osoi_liq(c,j).lt.watmin) then
                  h2osoi_liq(c,j) = watmin
                end if

                if (h2osoi_ice(c,j).lt.0) then
                  h2osoi_ice(c,j) = 0
                end if

                h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)

                if (j>1 .and. h2osoi_vol(c,j)>watsat(c,j)) then

                  var_temp = watsat(c,j) / h2osoi_vol(c,j)

                  h2osoi_liq(c,j) = h2osoi_liq(c,j)*var_temp

                  if (h2osoi_liq(c,j) < watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  h2osoi_ice(c,j) = (watsat(c,j)*dz(c,j)*denice)-h2osoi_liq(c,j)*denice/denh2o
                  if (abs(h2osoi_ice(c,j))<1.e-10_r8) then !numerics, if h2osoiice was zero before, it is -something e-16 after the previous calculation, so something marginal negative
                    h2osoi_ice(c,j) = 0._r8
                  end if

                  h2osoi_vol(c,j) = watsat(c,j)

                end if

              end if

            end do

          end do

        end if

        cc = cc+1

      end do

      ! ! snow

      ! cc = 1

      ! do count = 1, num_layer(1)

      !   g = hactiveg_levels(count,1)

      !   inc = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) - clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) ! save increment for gridcell

      !   if (abs(inc)>1.e-10_r8) then

      !     do count_columns = 1, num_layer_columns(1)

      !       c = hactivec_levels(count_columns,1)

      !       if (col%gridcell(c)==g) then

      !         inc_col = inc*h2osno(c)/clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2))

      !         if (inc_col/=inc_col) then
      !           inc_col = 0.0
      !         end if

      !         select case (update_snow)
      !         ! Tests with snow DA, scripts adapted from Lukas Strebel
      !         case(0)

      !           if (inc_col.ne.0._r8) then

      !             if (snl(c) < 0) then ! snow layers in the column

      !               h2osno(c) = h2osno(c) + inc_col

      !               do j=0,snl(c)+1,-1 ! iterate through the snow layers

      !                 ! snow density prior for each layer
      !                 if (dz(c,j)>0.0_r8) then
      !                   snowden = (h2osoi_liq(c,j) + h2osoi_ice(c,j)) / dz(c,j)
      !                 else
      !                   snowden = 0.0_r8
      !                 endif

      !                 ! fraction of SWE in each active layers
      !                 if(rsnow(c).gt.0.0_r8) then
      !                   frac_swe = (h2osoi_liq(c,j) + h2osoi_ice(c,j)) / rsnow(c)
      !                 else
      !                   frac_swe = 0.0_r8 ! no fraction SWE if no snow is present in column
      !                 end if ! end SWE fraction if

      !                 ! fraction of liquid and ice
      !                 if ((h2osoi_liq(c,j) + h2osoi_ice(c,j)).gt.0.0_r8) then
      !                   frac_liq = h2osoi_liq(c,j) / (h2osoi_liq(c,j) + h2osoi_ice(c,j))
      !                   frac_ice = 1.0_r8 - frac_liq
      !                 else
      !                   frac_liq = 0.0_r8
      !                   frac_ice = 0.0_r8
      !                 end if

      !                 ! SWE adjustment per layer 
      !                 ! assumes identical layer distribution of liq and ice than before DA (frac_*)
      !                 gain_h2osno = (h2osno(c) - rsnow(c)) * frac_swe
      !                 gain_h2oliq = gain_h2osno * frac_liq
      !                 gain_h2oice = gain_h2osno * frac_ice

      !                 ! layer level adjustments
      !                 if (snowden.gt.0.0_r8) then
      !                   gain_dzsno = gain_h2osno / snowden
      !                 else
      !                   gain_dzsno = 0.0_r8
      !                 end if
      !                 h2osoi_liq(c,j) = h2osoi_liq(c,j) + gain_h2oliq
      !                 h2osoi_ice(c,j) = h2osoi_ice(c,j) + gain_h2oice


      !                 ! Adjust snow layer dimensions so that CLM5 can calculate compaction / aggregation
      !                 ! in the DART code dzsno is adjusted directly but in CLM5 dzsno is local and diagnostic
      !                 ! i.e. calculated / assigned from frac_sno and dz(:, snow_layer) in SnowHydrologyMod
      !                 ! therefore we adjust dz(:, snow_layer) here

      !                 dz(c,j) = dz(c,j) + gain_dzsno
      !                 ! mid point and interface adjustments
      !                 ! i.e. zsno (col%z(:, snow_layers)) and zisno (col%zi(:, snow_layers))
      !                 ! DART version the sum goes from ilevel:nlevsno to fit with our indexing:
      !                 zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8
      !                 ! In DART the check is ilevel == nlevsno but here
                      
      !                 if (j.eq.0) then
      !                   z(c,j) = zi(c,j-1) / 2.0_r8
      !                 else
      !                   z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
      !                 end if


      !               end do

      !               ! Update the total snow depth to match updates to layers for active snow layers                
      !               snow_depth(c) = sum(dz(c,snl(c)+1:0))
      !               h2osno(c) = sum(h2osoi_ice(c,snl(c)+1:0)+h2osoi_liq(c,snl(c)+1:0))
                    
      !             end if

      !           end if

      !         case (1) !update with factor of old and new snow

      !           if (snl(c) < 0) then ! snow layers in the column

      !             ! if (inc_col>1000._r8) then
      !             !   inc_col = 1000._r8
      !             ! end if

      !             inc_col = h2osno(c)+inc_col
      !             scale = inc_col/h2osno(c)
      !             h2osno(c) = inc_col
                  
      !             do j=0,snl(c)+1,-1
      !               h2osoi_liq(c,j) = h2osoi_liq(c,j)*scale
      !               h2osoi_ice(c,j) = h2osoi_ice(c,j)*scale
      !               dz(c,j) = dz(c,j)*scale
      !               zi(c,j) = zi(c,j)*scale
      !               z(c,j) = z(c,j)*scale
      !             end do
      !             zi(c,snl(c)) = zi(c,snl(c))*scale
      !             snow_depth(c) = snow_depth(c)*scale

      !           end if

      !         end select

      !         ! snow negative
      !         if (h2osno(c) < 0._r8) then
      !           if (snl(c)<0) then

      !             do j=0,snl(c)+1,-1
      !                 h2osoi_liq(c,j) = 0.0_r8
      !                 h2osoi_ice(c,j) = 0.00000001_r8
      !                 dz(c,j)  = 0.00000001_r8  
      !                 zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8                 
      !                 if (j.eq.0) then
      !                   z(c,j) = zi(c,j-1) / 2.0_r8
      !                 else
      !                   z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
      !                 end if                 
      !             end do

      !           else

      !             h2osoi_liq(c,0) = 0.0_r8
      !             h2osoi_ice(c,0) = 0.00000001_r8
      !             dz(c,0)  = 0.00000001_r8
      !             zi(c,-1) = dz(c,0)*-1.0_r8 
      !             z(c,0) = zi(c,-1) / 2.0_r8

      !           end if

      !           snow_depth(c) = sum(dz(c,-nlevsno+1:0))
      !           h2osno(c) = sum(h2osoi_ice(c,-nlevsno+1:0))
      !         end if

      !       end if

      !     end do

      !   end if

      !   cc = cc+1

      ! end do

    elseif (state_setup.eq.4) then

      
      ! soil water

      cc = 1

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        inc = clm_statevec(cc)-clm_statevec_original_input(cc)

        if (abs(inc)>1.e-10_r8) then

          ! update soil water
          do j = 1,7

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (col%gridcell(c)==g) then

                inc_col = inc*h2osoi_liq_mean(c,j)/clm_statevec_original_input(cc)

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_liq(c,j)) then
                  inc_col = sign(max_inc*h2osoi_liq(c,j),inc_col)
                end if

                h2osoi_liq(c,j) = h2osoi_liq(c,j) + inc_col

                if (h2osoi_liq(c,j).lt.watmin) then
                  h2osoi_liq(c,j) = watmin
                end if




                inc_col = inc*h2osoi_ice_mean(c,j)/clm_statevec_original_input(cc)

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_ice(c,j)) then
                  inc_col = sign(max_inc*h2osoi_ice(c,j),inc_col)
                end if

                h2osoi_ice(c,j) = h2osoi_ice(c,j) + inc_col

                if (h2osoi_ice(c,j).lt.0) then
                  h2osoi_ice(c,j) = 0._r8
                end if



                h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)

                if (j>1 .and. h2osoi_vol(c,j)>watsat(c,j)) then

                  var_temp = watsat(c,j) / h2osoi_vol(c,j)

                  h2osoi_liq(c,j) = h2osoi_liq(c,j)*var_temp

                  if (h2osoi_liq(c,j) < watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  h2osoi_ice(c,j) = (watsat(c,j)*dz(c,j)*denice)-h2osoi_liq(c,j)*denice/denh2o
                  if (abs(h2osoi_ice(c,j))<1.e-10_r8) then !numerics, if h2osoiice was zero before, it is -something e-16 after the previous calculation, so something marginal negative
                    h2osoi_ice(c,j) = 0._r8
                  end if

                  h2osoi_vol(c,j) = watsat(c,j)

                end if

              end if

            end do
          end do

        end if

        cc = cc+1

      end do



      cc = 1

      do count = 1, num_layer(8)

        g = hactiveg_levels(count,8)

        inc = clm_statevec(cc+clm_varsize_tws(1))-clm_statevec_original_input(cc+clm_varsize_tws(1))

        if (abs(inc)>1.e-10_r8) then

          ! update soil water
          do j = 8,nlevsoi

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (col%gridcell(c)==g) then

                inc_col = inc*h2osoi_liq_mean(c,j)/clm_statevec_original_input(cc+clm_varsize_tws(1))

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_liq(c,j)) then
                  inc_col = sign(max_inc*h2osoi_liq(c,j),inc_col)
                end if

                h2osoi_liq(c,j) = h2osoi_liq(c,j) + inc_col

                if (h2osoi_liq(c,j).lt.watmin) then
                  h2osoi_liq(c,j) = watmin
                end if




                inc_col = inc*h2osoi_ice_mean(c,j)/clm_statevec_original_input(cc+clm_varsize_tws(1))

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_ice(c,j)) then
                  inc_col = sign(max_inc*h2osoi_ice(c,j),inc_col)
                end if

                h2osoi_ice(c,j) = h2osoi_ice(c,j) + inc_col

                if (h2osoi_ice(c,j).lt.0) then
                  h2osoi_ice(c,j) = 0._r8
                end if



                h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)

                if (j>1 .and. h2osoi_vol(c,j)>watsat(c,j)) then

                  var_temp = watsat(c,j) / h2osoi_vol(c,j)

                  h2osoi_liq(c,j) = h2osoi_liq(c,j)*var_temp

                  if (h2osoi_liq(c,j) < watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  h2osoi_ice(c,j) = (watsat(c,j)*dz(c,j)*denice)-h2osoi_liq(c,j)*denice/denh2o
                  if (abs(h2osoi_ice(c,j))<1.e-10_r8) then !numerics, if h2osoiice was zero before, it is -something e-16 after the previous calculation, so something marginal negative
                    h2osoi_ice(c,j) = 0._r8
                  end if

                  h2osoi_vol(c,j) = watsat(c,j)

                end if

              end if

            end do
          end do

        end if

        cc = cc+1

      end do



      ! snow

      cc = 1

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        inc = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) - clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)) ! save increment for gridcell
        
        if (abs(inc)>1.e-10_r8) then

          do count_columns = 1, num_layer_columns(1)

            c = hactivec_levels(count_columns,1)

            if (col%gridcell(c)==g) then

              inc_col = inc*h2osno_mean(c)/clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2))

              if (inc_col/=inc_col) then
                inc_col = 0.0
              end if

              if (abs(inc_col)>500) then
                inc_col = sign(500._r8,inc_col)
              end if

              select case (update_snow)
              ! Tests with snow DA, scripts adapted from Lukas Strebel
              case(0)

                if (inc_col.ne.0._r8) then

                  if (snl(c) < 0) then ! snow layers in the column

                    h2osno(c) = h2osno(c) + inc_col

                    do j=0,snl(c)+1,-1 ! iterate through the snow layers

                      ! snow density prior for each layer
                      if (dz(c,j)>0.0_r8) then
                        snowden = (h2osoi_liq(c,j) + h2osoi_ice(c,j)) / dz(c,j)
                      else
                        snowden = 0.0_r8
                      endif

                      ! fraction of SWE in each active layers
                      if(rsnow(c).gt.0.0_r8) then
                        frac_swe = (h2osoi_liq(c,j) + h2osoi_ice(c,j)) / rsnow(c)
                      else
                        frac_swe = 0.0_r8 ! no fraction SWE if no snow is present in column
                      end if ! end SWE fraction if

                      ! fraction of liquid and ice
                      if ((h2osoi_liq(c,j) + h2osoi_ice(c,j)).gt.0.0_r8) then
                        frac_liq = h2osoi_liq(c,j) / (h2osoi_liq(c,j) + h2osoi_ice(c,j))
                        frac_ice = 1.0_r8 - frac_liq
                      else
                        frac_liq = 0.0_r8
                        frac_ice = 0.0_r8
                      end if

                      ! SWE adjustment per layer 
                      ! assumes identical layer distribution of liq and ice than before DA (frac_*)
                      gain_h2osno = (h2osno(c) - rsnow(c)) * frac_swe
                      gain_h2oliq = gain_h2osno * frac_liq
                      gain_h2oice = gain_h2osno * frac_ice

                      ! layer level adjustments
                      if (snowden.gt.0.0_r8) then
                        gain_dzsno = gain_h2osno / snowden
                      else
                        gain_dzsno = 0.0_r8
                      end if
                      h2osoi_liq(c,j) = h2osoi_liq(c,j) + gain_h2oliq
                      h2osoi_ice(c,j) = h2osoi_ice(c,j) + gain_h2oice


                      ! Adjust snow layer dimensions so that CLM5 can calculate compaction / aggregation
                      ! in the DART code dzsno is adjusted directly but in CLM5 dzsno is local and diagnostic
                      ! i.e. calculated / assigned from frac_sno and dz(:, snow_layer) in SnowHydrologyMod
                      ! therefore we adjust dz(:, snow_layer) here

                      dz(c,j) = dz(c,j) + gain_dzsno
                      ! mid point and interface adjustments
                      ! i.e. zsno (col%z(:, snow_layers)) and zisno (col%zi(:, snow_layers))
                      ! DART version the sum goes from ilevel:nlevsno to fit with our indexing:
                      zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8
                      ! In DART the check is ilevel == nlevsno but here
                      
                      if (j.eq.0) then
                        z(c,j) = zi(c,j-1) / 2.0_r8
                      else
                        z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
                      end if


                    end do

                    ! Update the total snow depth to match updates to layers for active snow layers                
                    snow_depth(c) = sum(dz(c,snl(c)+1:0))
                    h2osno(c) = sum(h2osoi_ice(c,snl(c)+1:0)+h2osoi_liq(c,snl(c)+1:0))
                    
                  end if

                end if

              case (1) !update with factor of old and new snow

                if (snl(c) < 0) then ! snow layers in the column

                  ! if (inc_col>1000._r8) then
                  !   inc_col = 1000._r8
                  ! end if

                  inc_col = h2osno(c)+inc_col

                  if (inc_col>10000) then
                    inc_col = 10000
                  end if

                  scale = inc_col/h2osno(c)
                  h2osno(c) = inc_col
                  
                  do j=0,snl(c)+1,-1
                    h2osoi_liq(c,j) = h2osoi_liq(c,j)*scale
                    h2osoi_ice(c,j) = h2osoi_ice(c,j)*scale
                    dz(c,j) = dz(c,j)*scale
                    zi(c,j) = zi(c,j)*scale
                    z(c,j) = z(c,j)*scale
                  end do
                  zi(c,snl(c)) = zi(c,snl(c))*scale
                  snow_depth(c) = snow_depth(c)*scale

                end if

              end select

              ! snow negative
              if (h2osno(c) < 0._r8) then
                if (snl(c)<0) then

                  do j=0,snl(c)+1,-1
                      h2osoi_liq(c,j) = 0.0_r8
                      h2osoi_ice(c,j) = 0.00000001_r8
                      dz(c,j)  = 0.00000001_r8  
                      zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8                 
                      if (j.eq.0) then
                        z(c,j) = zi(c,j-1) / 2.0_r8
                      else
                        z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
                      end if                 
                  end do

                else

                  h2osoi_liq(c,0) = 0.0_r8
                  h2osoi_ice(c,0) = 0.00000001_r8
                  dz(c,0)  = 0.00000001_r8
                  zi(c,-1) = dz(c,0)*-1.0_r8 
                  z(c,0) = zi(c,-1) / 2.0_r8

                end if

                snow_depth(c) = sum(dz(c,-nlevsno+1:0))
                h2osno(c) = sum(h2osoi_ice(c,-nlevsno+1:0))
              end if

            end if

          end do

        end if

        cc = cc+1

      end do

    elseif (state_setup.eq.5) then   !--> I use this as default

      
      ! soil water

      cc = 1

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        inc = clm_statevec(cc)-clm_statevec_original_input(cc)

        if (abs(inc)>1.e-10_r8) then

          ! update soil water
          do j = 1,3

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (col%gridcell(c)==g) then

                inc_col = inc*h2osoi_liq_mean(c,j)/clm_statevec_original_input(cc)

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_liq(c,j)) then
                  inc_col = sign(max_inc*h2osoi_liq(c,j),inc_col)
                end if

                h2osoi_liq(c,j) = h2osoi_liq(c,j) + inc_col

                if (h2osoi_liq(c,j).lt.watmin) then
                  h2osoi_liq(c,j) = watmin
                end if




                inc_col = inc*h2osoi_ice_mean(c,j)/clm_statevec_original_input(cc)

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_ice(c,j)) then
                  inc_col = sign(max_inc*h2osoi_ice(c,j),inc_col)
                end if

                h2osoi_ice(c,j) = h2osoi_ice(c,j) + inc_col

                if (h2osoi_ice(c,j).lt.0) then
                  h2osoi_ice(c,j) = 0._r8
                end if



                h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)

                if (j>1 .and. h2osoi_vol(c,j)>watsat(c,j)) then
                !if (h2osoi_vol(c,j)>watsat(c,j)) then ! check if soil balancing error comes from this

                  var_temp = watsat(c,j) / h2osoi_vol(c,j)

                  h2osoi_liq(c,j) = h2osoi_liq(c,j)*var_temp

                  if (h2osoi_liq(c,j) < watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  h2osoi_ice(c,j) = (watsat(c,j)*dz(c,j)*denice)-h2osoi_liq(c,j)*denice/denh2o
                  if (abs(h2osoi_ice(c,j))<1.e-10_r8) then !numerics, if h2osoiice was zero before, it is -something e-16 after the previous calculation, so something marginal negative
                    h2osoi_ice(c,j) = 0._r8
                  end if

                  h2osoi_vol(c,j) = watsat(c,j)

                end if

              end if

            end do
          end do

        end if

        cc = cc+1

      end do



      cc = 1

      do count = 1, num_layer(4)

        g = hactiveg_levels(count,4)

        inc = clm_statevec(cc+clm_varsize_tws(1))-clm_statevec_original_input(cc+clm_varsize_tws(1))

        if (abs(inc)>1.e-10_r8) then

          ! update soil water
          do j = 4,12

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (col%gridcell(c)==g) then

                inc_col = inc*h2osoi_liq_mean(c,j)/clm_statevec_original_input(cc+clm_varsize_tws(1))

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_liq(c,j)) then
                  inc_col = sign(max_inc*h2osoi_liq(c,j),inc_col)
                end if

                h2osoi_liq(c,j) = h2osoi_liq(c,j) + inc_col

                if (h2osoi_liq(c,j).lt.watmin) then
                  h2osoi_liq(c,j) = watmin
                end if




                inc_col = inc*h2osoi_ice_mean(c,j)/clm_statevec_original_input(cc+clm_varsize_tws(1))

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_ice(c,j)) then
                  inc_col = sign(max_inc*h2osoi_ice(c,j),inc_col)
                end if

                h2osoi_ice(c,j) = h2osoi_ice(c,j) + inc_col

                if (h2osoi_ice(c,j).lt.0) then
                  h2osoi_ice(c,j) = 0._r8
                end if



                h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)

                if (h2osoi_vol(c,j)-watsat(c,j)>0.00001) then

                  var_temp = watsat(c,j) / h2osoi_vol(c,j)

                  h2osoi_liq(c,j) = h2osoi_liq(c,j)*var_temp

                  if (h2osoi_liq(c,j) < watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  h2osoi_ice(c,j) = (watsat(c,j)*dz(c,j)*denice)-h2osoi_liq(c,j)*denice/denh2o
                  if (abs(h2osoi_ice(c,j))<1.e-10_r8) then !numerics, if h2osoiice was zero before, it is -something e-16 after the previous calculation, so something marginal negative
                    h2osoi_ice(c,j) = 0._r8
                  end if

                  h2osoi_vol(c,j) = watsat(c,j)

                end if

              end if

            end do
          end do

        end if

        cc = cc+1

      end do

      cc = 1

      do count = 1, num_layer(13)

        g = hactiveg_levels(count,13)

        inc = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2))-clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2))

        if (abs(inc)>1.e-10_r8) then

          ! update soil water
          do j = 13,nlevsoi

            do count_columns = 1,num_layer_columns(j)

              c = hactivec_levels(count_columns,j)

              if (col%gridcell(c)==g) then

                inc_col = inc*h2osoi_liq_mean(c,j)/clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2))

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_liq(c,j)) then
                  inc_col = sign(max_inc*h2osoi_liq(c,j),inc_col)
                end if

                h2osoi_liq(c,j) = h2osoi_liq(c,j) + inc_col

                if (h2osoi_liq(c,j).lt.watmin) then
                  h2osoi_liq(c,j) = watmin
                end if




                inc_col = inc*h2osoi_ice_mean(c,j)/clm_statevec_original_input(cc+clm_varsize_tws(1))

                if (inc_col/=inc_col) then
                  inc_col = 0.0
                end if

                ! if increment larger than maximal increment, adjust it to maximal increment with the sign of old increment
                ! so that the direction of the increment is right
                if (abs(inc_col).gt.max_inc*h2osoi_ice(c,j)) then
                  inc_col = sign(max_inc*h2osoi_ice(c,j),inc_col)
                end if

                h2osoi_ice(c,j) = h2osoi_ice(c,j) + inc_col

                if (h2osoi_ice(c,j).lt.0) then
                  h2osoi_ice(c,j) = 0._r8
                end if



                h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)

                if (h2osoi_vol(c,j)-watsat(c,j)>0.00001) then

                  var_temp = watsat(c,j) / h2osoi_vol(c,j)

                  h2osoi_liq(c,j) = h2osoi_liq(c,j)*var_temp

                  if (h2osoi_liq(c,j) < watmin) then
                    h2osoi_liq(c,j) = watmin
                  end if

                  h2osoi_ice(c,j) = (watsat(c,j)*dz(c,j)*denice)-h2osoi_liq(c,j)*denice/denh2o
                  if (abs(h2osoi_ice(c,j))<1.e-10_r8) then !numerics, if h2osoiice was zero before, it is -something e-16 after the previous calculation, so something marginal negative
                    h2osoi_ice(c,j) = 0._r8
                  end if

                  h2osoi_vol(c,j) = watsat(c,j)

                end if

              end if

            end do
          end do

        end if

        cc = cc+1

      end do



      ! snow

      cc = 1

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        inc = clm_statevec(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) - clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3)) ! save increment for gridcell

        if (abs(inc)>1.e-10_r8) then

          do count_columns = 1, num_layer_columns(1)

            c = hactivec_levels(count_columns,1)

            if (col%gridcell(c)==g) then

              inc_col = inc*h2osno_mean(c)/clm_statevec_original_input(cc+clm_varsize_tws(1)+clm_varsize_tws(2)+clm_varsize_tws(3))

              if (inc_col/=inc_col) then
                inc_col = 0.0
              end if

              if (abs(inc_col).gt.max_inc*h2osno(c)) then
                inc_col = sign(max_inc*h2osno(c),inc_col)
              end if

              select case (update_snow)
              ! Tests with snow DA, scripts adapted from Lukas Strebel
              case(0)

                if (inc_col.ne.0._r8) then

                  if (snl(c) < 0) then ! snow layers in the column

                    h2osno(c) = h2osno(c) + inc_col

                    do j=0,snl(c)+1,-1 ! iterate through the snow layers

                      ! snow density prior for each layer
                      if (dz(c,j)>0.0_r8) then
                        snowden = (h2osoi_liq(c,j) + h2osoi_ice(c,j)) / dz(c,j)
                      else
                        snowden = 0.0_r8
                      endif

                      ! fraction of SWE in each active layers
                      if(rsnow(c).gt.0.0_r8) then
                        frac_swe = (h2osoi_liq(c,j) + h2osoi_ice(c,j)) / rsnow(c)
                      else
                        frac_swe = 0.0_r8 ! no fraction SWE if no snow is present in column
                      end if ! end SWE fraction if

                      ! fraction of liquid and ice
                      if ((h2osoi_liq(c,j) + h2osoi_ice(c,j)).gt.0.0_r8) then
                        frac_liq = h2osoi_liq(c,j) / (h2osoi_liq(c,j) + h2osoi_ice(c,j))
                        frac_ice = 1.0_r8 - frac_liq
                      else
                        frac_liq = 0.0_r8
                        frac_ice = 0.0_r8
                      end if

                      ! SWE adjustment per layer 
                      ! assumes identical layer distribution of liq and ice than before DA (frac_*)
                      gain_h2osno = (h2osno(c) - rsnow(c)) * frac_swe
                      gain_h2oliq = gain_h2osno * frac_liq
                      gain_h2oice = gain_h2osno * frac_ice

                      ! layer level adjustments
                      if (snowden.gt.0.0_r8) then
                        gain_dzsno = gain_h2osno / snowden
                      else
                        gain_dzsno = 0.0_r8
                      end if
                      h2osoi_liq(c,j) = h2osoi_liq(c,j) + gain_h2oliq
                      h2osoi_ice(c,j) = h2osoi_ice(c,j) + gain_h2oice


                      ! Adjust snow layer dimensions so that CLM5 can calculate compaction / aggregation
                      ! in the DART code dzsno is adjusted directly but in CLM5 dzsno is local and diagnostic
                      ! i.e. calculated / assigned from frac_sno and dz(:, snow_layer) in SnowHydrologyMod
                      ! therefore we adjust dz(:, snow_layer) here

                      dz(c,j) = dz(c,j) + gain_dzsno
                      ! mid point and interface adjustments
                      ! i.e. zsno (col%z(:, snow_layers)) and zisno (col%zi(:, snow_layers))
                      ! DART version the sum goes from ilevel:nlevsno to fit with our indexing:
                      zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8
                      ! In DART the check is ilevel == nlevsno but here
                      
                      if (j.eq.0) then
                        z(c,j) = zi(c,j-1) / 2.0_r8
                      else
                        z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
                      end if


                    end do

                    ! Update the total snow depth to match updates to layers for active snow layers                
                    snow_depth(c) = sum(dz(c,snl(c)+1:0))
                    h2osno(c) = sum(h2osoi_ice(c,snl(c)+1:0)+h2osoi_liq(c,snl(c)+1:0))
                    
                  end if

                end if

              case (1) !update with factor of old and new snow

                if (snl(c) < 0) then ! snow layers in the column

                  ! if (inc_col>1000._r8) then
                  !   inc_col = 1000._r8
                  ! end if

                  inc_col = h2osno(c)+inc_col
                  scale = inc_col/h2osno(c)
                  h2osno(c) = inc_col
                  
                  do j=0,snl(c)+1,-1
                    h2osoi_liq(c,j) = h2osoi_liq(c,j)*scale
                    h2osoi_ice(c,j) = h2osoi_ice(c,j)*scale
                    dz(c,j) = dz(c,j)*scale
                    zi(c,j) = zi(c,j)*scale
                    z(c,j) = z(c,j)*scale
                  end do
                  zi(c,snl(c)) = zi(c,snl(c))*scale
                  snow_depth(c) = snow_depth(c)*scale

                end if

              end select

              ! snow negative
              if (h2osno(c) < 0._r8) then
                if (snl(c)<0) then

                  do j=0,snl(c)+1,-1
                      h2osoi_liq(c,j) = 0.0_r8
                      h2osoi_ice(c,j) = 0.00000001_r8
                      dz(c,j)  = 0.00000001_r8  
                      zi(c,j-1) = sum(dz(c,j:0))*-1.0_r8                 
                      if (j.eq.0) then
                        z(c,j) = zi(c,j-1) / 2.0_r8
                      else
                        z(c,j) = sum(zi(c,j-1:j)) / 2.0_r8
                      end if                 
                  end do

                else

                  h2osoi_liq(c,0) = 0.0_r8
                  h2osoi_ice(c,0) = 0.00000001_r8
                  dz(c,0)  = 0.00000001_r8
                  zi(c,-1) = dz(c,0)*-1.0_r8 
                  z(c,0) = zi(c,-1) / 2.0_r8

                end if

                snow_depth(c) = sum(dz(c,-nlevsno+1:0))
                h2osno(c) = sum(h2osoi_ice(c,-nlevsno+1:0))
              end if

            end if

          end do

        end if

        cc = cc+1

      end do


    end if


    do j = 1,nlevsoi
      do count = 1,num_layer_columns(j)
        c = hactivec_levels(count,j)

        h2osoi_liq_inc(c,j) = h2osoi_liq(c,j)-h2osoi_liq_inc(c,j)
        h2osoi_ice_inc(c,j) = h2osoi_ice(c,j)-h2osoi_ice_inc(c,j)

        if (j==1) then
          h2osno_inc(c) = h2osno(c)-h2osno_inc(c)
        end if
      end do
    end do




    ! fill state after assimilation (the same script as filling the state vector before the assimilation, now only after the assimilation)

    select case (state_setup)

    case(0) ! all compartments, liq and ice water indidually

      do j = 1,nlevsoi 

        do count = 1, num_layer(j)

          g = hactiveg_levels(count,j)

          avg_sum = 0
          avg_sum_ice = 0
          avg_divide = 0

          do count_columns = 1,num_layer_columns(j)

            c = hactivec_levels(count_columns,j)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osoi_liq(c,j)
              avg_sum_ice = avg_sum_ice + h2osoi_ice(c,j)

              avg_divide = avg_divide+1


            end if

          end do

          h2osoi_liq_state(g,j) = avg_sum/avg_divide
          h2osoi_ice_state(g,j) = avg_sum_ice/avg_divide

          avg_sum = 0
          avg_divide = 0
          if (j==1) then
            ! snow
            avg_sum = 0
            avg_divide = 0
            do count_columns = 1,num_layer_columns(j)
              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osno(c)

                avg_divide = avg_divide+1


              end if

            end do

            h2osno_state(g) = avg_sum/avg_divide

          end if

        end do
      end do



    case(1) ! all compartments, sum of ice and liq soil water to overcome balancing errors due to different partitioning of water caused by different temperature

      do j = 1,nlevsoi 

        do count = 1, num_layer(j)

          g = hactiveg_levels(count,j)

          avg_sum = 0
          avg_sum_ice = 0
          avg_divide = 0

          do count_columns = 1,num_layer_columns(j)

            c = hactivec_levels(count_columns,j)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

              avg_divide = avg_divide+1

            end if

          end do

          h2osoi_liq_state(g,j) = avg_sum/avg_divide

          avg_sum = 0
          avg_divide = 0
          if (j==1) then
            ! snow
            avg_sum = 0
            avg_divide = 0
            do count_columns = 1,num_layer_columns(j)
              c = hactivec_levels(count_columns,j)

              if (g==col%gridcell(c)) then

                avg_sum = avg_sum + h2osno(c)

                avg_divide = avg_divide+1


              end if

            end do

            h2osno_state(g) = avg_sum/avg_divide

          end if

        end do
      end do



    case(2) ! only TWS in statevector

      do count = 1, num_layer(1)

        h2osoi_liq_state(g,1) = hactiveg_levels(count,1)

        h2osoi_liq_state(g,1) = TWS(g)

      end do

    case(3) ! sum over all soil layers and snow in statevector

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        h2osoi_liq_state(g,1) = 0

        do j = 1, nlevsoi

          avg_sum = 0
          avg_sum_ice = 0
          avg_divide = 0

          do count_columns = 1,num_layer_columns(j)

            c = hactivec_levels(count_columns,j)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

              avg_divide = avg_divide+1

            end if

          end do

          if (avg_divide.ne.0) then
            h2osoi_liq_state(g,1) = h2osoi_liq_state(g,1) + avg_sum/avg_divide
          end if

        end do

      end do


      ! snow

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        avg_sum = 0
        avg_divide = 0
        do count_columns = 1,num_layer_columns(1)
          c = hactivec_levels(count_columns,1)

          if (g==col%gridcell(c)) then

            avg_sum = avg_sum + h2osno(c)

            avg_divide = avg_divide+1


          end if

        end do

        h2osno_state(g) = avg_sum/avg_divide

      end do

    case(4)

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        h2osoi_liq_state(g,1) = 0

        do j = 1, 7

          avg_sum = 0
          avg_divide = 0

          do count_columns = 1,num_layer_columns(j)

            c = hactivec_levels(count_columns,j)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

              avg_divide = avg_divide+1

            end if

          end do

          if (avg_divide.ne.0) then

            h2osoi_liq_state(g,1) = h2osoi_liq_state(g,1) + avg_sum/avg_divide

          end if

        end do

      end do

      do count = 1, num_layer(8)

        g = hactiveg_levels(count,8)

        h2osoi_liq_state(g,2) = 0

        do j = 8, nlevsoi

          avg_sum = 0
          avg_divide = 0

          do count_columns = 1,num_layer_columns(j)

            c = hactivec_levels(count_columns,j)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

              avg_divide = avg_divide+1

            end if

          end do

          if (avg_divide.ne.0) then

            h2osoi_liq_state(g,2) = h2osoi_liq_state(g,2) + avg_sum/avg_divide

          end if

        end do

      end do

      ! snow

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        avg_sum = 0
        avg_divide = 0
        do count_columns = 1,num_layer_columns(1)
          c = hactivec_levels(count_columns,1)

          if (g==col%gridcell(c)) then

            avg_sum = avg_sum + h2osno(c)

            avg_divide = avg_divide+1


          end if

        end do

        h2osno_state(g) = avg_sum/avg_divide

      end do

    case(5)   --> I use this as default

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        h2osoi_liq_state(g,1) = 0

        do j = 1, 3

          avg_sum = 0
          avg_divide = 0

          do count_columns = 1,num_layer_columns(j)

            c = hactivec_levels(count_columns,j)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

              avg_divide = avg_divide+1

            end if

          end do

          if (avg_divide.ne.0) then

            h2osoi_liq_state(g,1) = h2osoi_liq_state(g,1) + avg_sum/avg_divide

          end if

        end do

      end do

      do count = 1, num_layer(4)

        g = hactiveg_levels(count,4)

        h2osoi_liq_state(g,2) = 0

        do j = 4, 12

          avg_sum = 0
          avg_divide = 0

          do count_columns = 1,num_layer_columns(j)

            c = hactivec_levels(count_columns,j)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

              avg_divide = avg_divide+1

            end if

          end do

          if (avg_divide.ne.0) then

            h2osoi_liq_state(g,2) = h2osoi_liq_state(g,2) + avg_sum/avg_divide

          end if

        end do

      end do

      do count = 1, num_layer(13)

        g = hactiveg_levels(count,13)

        h2osoi_liq_state(g,3) = 0

        do j = 13, nlevsoi

          avg_sum = 0
          avg_divide = 0

          do count_columns = 1,num_layer_columns(j)

            c = hactivec_levels(count_columns,j)

            if (g==col%gridcell(c)) then

              avg_sum = avg_sum + h2osoi_liq(c,j) + h2osoi_ice(c,j)

              avg_divide = avg_divide+1

            end if

          end do

          if (avg_divide.ne.0) then

            h2osoi_liq_state(g,3) = h2osoi_liq_state(g,3) + avg_sum/avg_divide

          end if

        end do

      end do


      ! snow

      do count = 1, num_layer(1)

        g = hactiveg_levels(count,1)

        avg_sum = 0
        avg_divide = 0
        do count_columns = 1,num_layer_columns(1)
          c = hactivec_levels(count_columns,1)

          if (g==col%gridcell(c)) then

            avg_sum = avg_sum + h2osno(c)

            avg_divide = avg_divide+1

          end if

        end do

        if (avg_divide.ne.0) then

          h2osno_state(g) = avg_sum/avg_divide

        end if

      end do

    end select


  end subroutine

  subroutine clm_correct_texture()

    use clm_varpar   , only : nlevsoi
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_instMod, only : soilstate_inst
    use CNSharedParamsMod, only : CNParamsShareInst

    implicit none

    integer :: c,lev
    real(r8) :: clay,sand,ttot
    real(r8) :: orgm
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col

    do c = clm_begg,clm_endg
      do lev = 1,nlevsoi
         clay = pclay(c,lev)
         sand = psand(c,lev)

         if(sand.le.0.0) sand = 1.0
         if(clay.le.0.0) clay = 1.0

         ttot = sand + clay
         if(ttot.gt.100) then
             sand = sand/ttot * 100.0
             clay = clay/ttot * 100.0
         end if

         pclay(c,lev) = clay
         psand(c,lev) = sand

         if(clmupdate_texture.eq.2) then
           orgm = (porgm(c,lev) / CNParamsShareInst%organic_max) * 100.0
           if(orgm.le.0.0) orgm = 0.0
           if(orgm.ge.100.0) orgm = 100.0
           porgm(c,lev) = (orgm / 100.0)* CNParamsShareInst%organic_max
         end if

       end do
    end do
  end subroutine clm_correct_texture

  subroutine clm_texture_to_parameters()
    use clm_varpar   , only : nlevsoi
    use clm_varcon   , only : zsoi, secspday
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_instMod, only : soilstate_inst
    use FuncPedotransferMod , only : pedotransf, get_ipedof
    use CNSharedParamsMod, only : CNParamsShareInst
    implicit none

    real(r8)           :: om_tkm         = 0.25_r8      
    ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8)           :: om_watsat_lake = 0.9_r8       
    ! porosity of organic soil
    real(r8)           :: om_hksat_lake  = 0.1_r8       
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat_lake = 10.3_r8      
    ! saturated suction for organic matter (Letts, 2000)
    real(r8)           :: om_b_lake      = 2.7_r8       
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000) (lake)
    real(r8)           :: om_watsat                     
    ! porosity of organic soil
    real(r8)           :: om_hksat                      
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat                     
    ! saturated suction for organic matter (mm)(Letts, 2000)
    real(r8)           :: om_csol        = 2.5_r8       
    ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8)           :: om_tkd         = 0.05_r8      
    ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8)           :: om_b                          
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8)           :: zsapric        = 0.5_r8
    ! depth (m) that organic matter takes on characteristics of sapric peat
    real(r8)           :: pcalpha        = 0.5_r8       ! percolation threshold
    real(r8)           :: pcbeta         = 0.139_r8     ! percolation exponent
    real(r8)           :: perc_frac    ! "percolating" fraction of organic soil
    real(r8)           :: perc_norm    ! normalize to 1 when 100% organic soil
    real(r8)           :: uncon_hksat  ! series conductivity of mineral/organic soil
    real(r8)           :: uncon_frac   ! fraction of "unconnected" soil
    real(r8)           :: bd           ! bulk density of dry soil material [kg/m^3]
    real(r8)           :: tkm          ! mineral conductivity 
    real(r8)           :: xksat        ! maximum hydraulic conductivity of soil [mm/s]
    
    integer  :: ipedof,c,lev
    real(r8) :: clay,sand,om_frac
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm   => soilstate_inst%cellorg_col

    do c = clm_begg,clm_endg
      do lev = 1,nlevsoi
         clay = pclay(c,lev)
         sand = psand(c,lev)
         om_frac = porgm(c,lev) / CNParamsShareInst%organic_max

         ipedof=get_ipedof(0)
         call pedotransf(ipedof, sand, clay, &
                         soilstate_inst%watsat_col(c,lev), &
                         soilstate_inst%bsw_col(c,lev), &
                         soilstate_inst%sucsat_col(c,lev), &
                         xksat)

         om_watsat         = max(0.93_r8 - 0.1_r8   *(zsoi(lev)/zsapric), 0.83_r8)
         om_b              = min(2.7_r8  + 9.3_r8   *(zsoi(lev)/zsapric), 12.0_r8)
         om_sucsat         = min(10.3_r8 - 0.2_r8   *(zsoi(lev)/zsapric), 10.1_r8)
         om_hksat          = max(0.28_r8 - 0.2799_r8*(zsoi(lev)/zsapric), xksat)

         soilstate_inst%bd_col(c,lev) = &
              (1._r8 - soilstate_inst%watsat_col(c,lev))*2.7e3_r8
         
         soilstate_inst%watsat_col(c,lev) = &
              (1._r8 - om_frac) * soilstate_inst%watsat_col(c,lev) + om_watsat*om_frac
         
         tkm = (1._r8 - om_frac) * (8.80_r8*sand+2.92_r8*clay)/(sand+clay)+om_tkm*om_frac 
         ! W/(m K)
         soilstate_inst%bsw_col(c,lev) = &
              (1._r8-om_frac) * (2.91_r8 + 0.159_r8*clay) + om_frac*om_b   
         
         soilstate_inst%sucsat_col(c,lev) = &
              (1._r8-om_frac) * soilstate_inst%sucsat_col(c,lev) + om_sucsat*om_frac  
         
         soilstate_inst%hksat_min_col(c,lev) = xksat

         ! perc_frac is zero unless perf_frac greater than percolation threshold
         if (om_frac > pcalpha) then
            perc_norm = (1._r8 - pcalpha)**(-pcbeta)
            perc_frac = perc_norm*(om_frac - pcalpha)**pcbeta
         else
            perc_frac = 0._r8
         endif

         ! uncon_frac is fraction of mineral soil plus fraction of
         ! "nonpercolating" organic soil
         uncon_frac = (1._r8-om_frac)+(1._r8-perc_frac)*om_frac

         ! uncon_hksat is series addition of mineral/organic
         ! conductivites
         if (om_frac < 1._r8) then
            uncon_hksat = uncon_frac/((1._r8-om_frac)/xksat &
                 +((1._r8-perc_frac)*om_frac)/om_hksat)
         else
            uncon_hksat = 0._r8
         end if

         soilstate_inst%hksat_col(c,lev)  = uncon_frac*uncon_hksat + &
             (perc_frac*om_frac)*om_hksat

         soilstate_inst%tkmg_col(c,lev)   = tkm ** (1._r8 - &
             soilstate_inst%watsat_col(c,lev))           

         soilstate_inst%tksatu_col(c,lev) = &
             soilstate_inst%tkmg_col(c,lev)*0.57_r8**soilstate_inst%watsat_col(c,lev)

         soilstate_inst%tkdry_col(c,lev)  = &
             ((0.135_r8*soilstate_inst%bd_col(c,lev) + 64.7_r8) / &
             (2.7e3_r8 - 0.947_r8*soilstate_inst%bd_col(c,lev)) &
             )*(1._r8-om_frac) + om_tkd*om_frac  

         soilstate_inst%csol_col(c,lev)   = &
             ((1._r8-om_frac)*(2.128_r8*sand+2.385_r8*clay) / (sand+clay) + &
             om_csol*om_frac)*1.e6_r8

         soilstate_inst%watdry_col(c,lev) = &
             soilstate_inst%watsat_col(c,lev) * &
             (316230._r8/soilstate_inst%sucsat_col(c,lev)  &
             ) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 

         soilstate_inst%watopt_col(c,lev) = & 
             soilstate_inst%watsat_col(c,lev) * &
             (158490._r8/soilstate_inst%sucsat_col(c,lev)  &
             ) ** (-1._r8/soilstate_inst%bsw_col(c,lev)) 

         ! secspday (day/sec)
         soilstate_inst%watfc_col(c,lev) = &
             soilstate_inst%watsat_col(c,lev) * &
             (0.1_r8 / (soilstate_inst%hksat_col(c,lev)*secspday) &
             )**(1._r8/(2._r8*soilstate_inst%bsw_col(c,lev)+3._r8))

      end do
    end do
  end subroutine clm_texture_to_parameters

  ! subroutine  average_swc_crp(profdat,profave)
  !   use clm_varcon  , only : zsoi

  !   implicit none

  !   real(r8),intent(in)  :: profdat(10)
  !   real(r8),intent(out) :: profave

  !   error stop "Not implemented average_swc_crp"
  ! end subroutine average_swc_crp
#endif

  subroutine domain_def_clm(lon_clmobs, lat_clmobs, dim_obs, &
                            longxy, latixy, longxy_obs, latixy_obs)
    use spmdMod,   only : npes, iam
    use domainMod, only : ldomain
    use decompMod, only : get_proc_total, get_proc_bounds, ldecomp

    implicit none
    real, intent(in) :: lon_clmobs(:)
    real, intent(in) :: lat_clmobs(:)
    integer, intent(in) :: dim_obs
    integer, allocatable, intent(inout) :: longxy(:)
    integer, allocatable, intent(inout) :: latixy(:)
    integer, allocatable, intent(inout) :: longxy_obs(:)
    integer, allocatable, intent(inout) :: latixy_obs(:)
    integer :: ni, nj, ii, jj, kk, cid, ier, ncells, nlunits
    integer :: ncols, counter
    integer :: npatches, ncohorts
    real :: minlon, minlat, maxlon, maxlat
    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices


    lon => ldomain%lonc
    lat => ldomain%latc
    ni = ldomain%ni
    nj = ldomain%nj

    ! get total number of gridcells, landunits,
    ! columns, patches and cohorts on processor

    call get_proc_total(iam, ncells, nlunits, ncols, npatches, ncohorts)

    ! beg and end gridcell
    call get_proc_bounds(begg=begg, endg=endg)
    
   !print *,'ni, nj ', ni, nj
   !print *,'cells per processor ', ncells
   !print *,'begg, endg ', begg, endg

    ! allocate vector with size of elements in x directions * size of elements in y directions
    allocate(longxy(ncells), stat=ier)
    allocate(latixy(ncells), stat=ier)

    ! initialize vector with zero values
    longxy(:) = 0
    latixy(:) = 0
  
    ! fill vector with index values
    counter = 1
    do ii = 1, nj
      do jj = 1, ni
        cid = (ii-1)*ni + jj
        do kk = begg, endg
          if(cid == ldecomp%gdc2glo(kk)) then
            latixy(counter) = ii
            longxy(counter) = jj
            counter = counter + 1
          end if
        end do
      end do
    end do

    ! set intial values for max/min of lon/lat
    minlon = 999
    minlat = 999
    maxlon = -999
    maxlat = -999

    ! looping over all cell centers to get min/max longitude and latitude
    minlon = MINVAL(lon(:) + 180)
    maxlon = MAXVAL(lon(:) + 180)
    minlat = MINVAL(lat(:) + 90)
    maxlat = MAXVAL(lat(:) + 90)

    allocate(longxy_obs(dim_obs), stat=ier)
    allocate(latixy_obs(dim_obs), stat=ier)

    do i = 1, dim_obs
       if(((lon_clmobs(i) + 180) - minlon) /= 0 .and. &
         ((lat_clmobs(i) + 90) - minlat) /= 0) then
          longxy_obs(i) = ceiling(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon)) !+ 1
          latixy_obs(i) = ceiling(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat)) !+ 1
          !print *,'longxy_obs(i) , latixy_obs(i) ', longxy_obs(i) , latixy_obs(i)
        else if(((lon_clmobs(i) + 180) - minlon) == 0 .and. &
                ((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs(i) = 1
          latixy_obs(i) = 1
       else if(((lon_clmobs(i) + 180) - minlon) == 0) then
          longxy_obs(i) = 1
          latixy_obs(i) = ceiling(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat))
       else if(((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs(i) = ceiling(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon))
          latixy_obs(i) = 1
       endif
    end do
    ! deallocate temporary arrays
    !deallocate(longxy)
    !deallocate(latixy)
    !deallocate(longxy_obs)
    !deallocate(latixy_obs)

  end subroutine domain_def_clm

  !> @author  Mukund Pondkule, Johannes Keller
  !> @date    27.03.2023
  !> @brief   Set indices of grid cells with lon/lat smaller than observation locations
  !> @details
  !>    This routine sets the indices of grid cells with lon/lat
  !>    smaller than observation locations.
  subroutine get_interp_idx(lon_clmobs, lat_clmobs, dim_obs, longxy_obs_floor, latixy_obs_floor)

    USE domainMod, ONLY: ldomain
    ! USE decompMod, ONLY: get_proc_total, get_proc_bounds_atm, adecomp
    ! USE spmdMod,   ONLY: npes, iam ! number of processors for clm and processor number
   ! USE mod_read_obs, ONLY: lon_clmobs, lat_clmobs
    ! USE clmtype,    ONLY : clm3
   ! USE mod_assimilation, ONLY: dim_obs
!#endif

    implicit none
    real, intent(in) :: lon_clmobs(:)
    real, intent(in) :: lat_clmobs(:)
    integer, intent(in) :: dim_obs
    integer, allocatable, intent(inout) :: longxy_obs_floor(:)
    integer, allocatable, intent(inout) :: latixy_obs_floor(:)
    integer :: i

    integer :: ni, nj
    ! integer :: ii, jj, kk, cid
    integer :: ier
    ! integer :: ncells, nlunits,ncols, npfts
    ! integer :: counter

    real :: minlon, minlat, maxlon, maxlat
    real(r8), pointer :: lon(:)
    real(r8), pointer :: lat(:)
    ! integer :: begg, endg   ! per-proc gridcell ending gridcell indices

    lon => ldomain%lonc
    lat => ldomain%latc
    ni = ldomain%ni
    nj = ldomain%nj

    ! ! get total number of gridcells, landunits,
    ! ! columns and pfts for any processor
    ! call get_proc_total(iam, ncells, nlunits, ncols, npfts)

    ! ! beg and end gridcell for atm
    ! call get_proc_bounds_atm(begg, endg)

    ! set intial values for max/min of lon/lat
    minlon = 999
    minlat = 999
    maxlon = -999
    maxlat = -999

    ! looping over all cell centers to get min/max longitude and latitude
    minlon = MINVAL(lon(:) + 180)
    maxlon = MAXVAL(lon(:) + 180)
    minlat = MINVAL(lat(:) + 90)
    maxlat = MAXVAL(lat(:) + 90)

    if(allocated(longxy_obs_floor)) deallocate(longxy_obs_floor)
    allocate(longxy_obs_floor(dim_obs), stat=ier)
    if(allocated(latixy_obs_floor)) deallocate(latixy_obs_floor)
    allocate(latixy_obs_floor(dim_obs), stat=ier)
    do i = 1, dim_obs
       if(((lon_clmobs(i) + 180) - minlon) /= 0 .and. ((lat_clmobs(i) + 90) - minlat) /= 0) then
          longxy_obs_floor(i) = floor(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon)) !+ 1
          latixy_obs_floor(i) = floor(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat)) !+ 1
          !print *,'longxy_obs(i) , latixy_obs(i) ', longxy_obs(i) , latixy_obs(i)
       else if(((lon_clmobs(i) + 180) - minlon) == 0 .and. ((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs_floor(i) = 1
          latixy_obs_floor(i) = 1
       else if(((lon_clmobs(i) + 180) - minlon) == 0) then
          longxy_obs_floor(i) = 1
          latixy_obs_floor(i) = floor(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat))
       else if(((lat_clmobs(i) + 90) - minlat) == 0) then
          longxy_obs_floor(i) = floor(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon))
          latixy_obs_floor(i) = 1
       endif
    end do

  end subroutine get_interp_idx

#if defined CLMSA
  subroutine init_clm_l_size(dim_l)
    use clm_varpar   , only : nlevsoi

    implicit none

    integer, intent(out) :: dim_l
    integer              :: nshift

    if(clmupdate_swc.eq.1) then
      dim_l = nlevsoi
      nshift = nlevsoi
    endif

    if(clmupdate_swc.eq.2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
      ! dim_l = nlevsoi + 1
      ! nshift = nlevsoi + 1
    endif

    if(clmupdate_texture.eq.1) then
      dim_l = 2*nlevsoi + nshift
    endif

    if(clmupdate_texture.eq.2) then
      dim_l = 3*nlevsoi + nshift
    endif

    if (clmupdate_tws.eq.1) then
      select case (state_setup)
      case(0)
        dim_l = 2*nlevsoi+3
      case(2)
        dim_l = 1
      case default 
        dim_l = 1*nlevsoi+3
      end select
    end if

  end subroutine init_clm_l_size
#endif

  !> @author Yorck Ewerdwalbesloh
  !> @date 05.09.2023
  !> @brief reading TWS temporal mean model file
  !> @param[in] temp_mean_filename Name of mean file
  !> @details
  !> This subroutine reads a provided temporal mean model file
  subroutine read_temp_mean_model(temp_mean_filename)
    
    use netcdf
    use mod_read_obs, only: check
    implicit none
    integer :: ncid, dim_lon, dim_lat, lon_varid, lat_varid, tws_varid
    character (len = *), parameter :: dim_lon_name = "lsmlon"
    character (len = *), parameter :: dim_lat_name = "lsmlat"
    character (len = *), parameter :: lon_name = "longitude"
    character (len = *), parameter :: lat_name = "latitude"
    character (len = *), parameter :: tws_name = "TWS"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid_lon, dimid_lat, status
    integer :: haserr
    character (len = *), intent(in) :: temp_mean_filename

    !print *, "Read temporal mean of CLM OL run"

    call check(nf90_open(temp_mean_filename, nf90_nowrite, ncid))
    call check(nf90_inq_dimid(ncid, dim_lon_name, dimid_lon))
    call check(nf90_inq_dimid(ncid, dim_lat_name, dimid_lat))
    call check(nf90_inquire_dimension(ncid, dimid_lon, recorddimname, dim_lon))
    call check(nf90_inquire_dimension(ncid, dimid_lat, recorddimname, dim_lat))
    
    if(allocated(lon_temp_mean))deallocate(lon_temp_mean)
    if(allocated(lat_temp_mean))deallocate(lat_temp_mean)
    if(allocated(tws_temp_mean))deallocate(tws_temp_mean)

    allocate(tws_temp_mean(dim_lon,dim_lat))
    allocate(lon_temp_mean(dim_lon,dim_lat))
    allocate(lat_temp_mean(dim_lon,dim_lat))

    call check( nf90_inq_varid(ncid, lon_name, lon_varid))
    call check(nf90_get_var(ncid, lon_varid, lon_temp_mean))

    call check( nf90_inq_varid(ncid, lat_name, lat_varid))
    call check(nf90_get_var(ncid, lat_varid, lat_temp_mean))

    call check( nf90_inq_varid(ncid, tws_name, tws_varid))
    call check(nf90_get_var(ncid, tws_varid, tws_temp_mean))

    call check( nf90_close(ncid) )

  end subroutine

  ! subroutine check(status)
  
  !   use netcdf
  !   integer, intent ( in) :: status

  !   if(status /= nf90_noerr) then
  !      print *, trim(nf90_strerror(status))
  !      stop "Stopped"
  !   end if
  ! end subroutine check

end module enkf_clm_mod


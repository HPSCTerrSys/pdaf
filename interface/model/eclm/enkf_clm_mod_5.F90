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

  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
  use, intrinsic :: IEEE_ARITHMETIC, only: ieee_is_nan

! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8, SHR_KIND_CL

! !ARGUMENTS:
    implicit none

    public

#if (defined CLMSA)
  integer :: COMM_model_clm
  integer :: clm_statevecsize
  integer :: clm_paramsize !hcp: Size of CLM parameter vector (f.e. LAI)
  integer :: clm_varsize
  integer :: clm_begg,clm_endg
  integer :: clm_begc,clm_endc
  integer :: clm_begp,clm_endp
  real(r8),allocatable :: clm_statevec(:)
  real(r8),allocatable :: clm_statevec_orig(:)
  integer,allocatable :: state_pdaf2clm_c_p(:)
  integer,allocatable :: state_pdaf2clm_j_p(:)
  integer,allocatable :: state_loc2clm_c_p(:)
  ! clm_paramarr: Contains LAI used in obs_op_pdaf for computing model
  ! LST in LST assimilation (clmupdate_T)
  real(r8),allocatable :: clm_paramarr(:)  !hcp CLM parameter vector (f.e. LAI)
  integer, allocatable :: state_clm2pdaf_p(:,:) !Index of column in hydraulic active state vector (nlevsoi,endc-begc+1)
  real(r8),allocatable :: clm_patch2gc(:) ! index array to simplify patch to gridcell averaging
  real(r8),allocatable :: clm_patchwt(:) ! weight array to simplify patch to gridcell averaging
  integer,allocatable :: clm_lai_patch_idx(:)   ! mode-3: CLM patch index per state slot
  integer,allocatable :: clm_lai_patch_itype(:) ! mode-3: PFT type per state slot
  integer(c_int),bind(C,name="clmupdate_swc")     :: clmupdate_swc
  integer(c_int),bind(C,name="clmupdate_T")     :: clmupdate_T  ! by hcp
  integer(c_int),bind(C,name="clmupdate_texture") :: clmupdate_texture
  integer(c_int),bind(C,name="clmprint_swc")      :: clmprint_swc
  integer(c_int),bind(C,name="clmupdate_lai")    :: clmupdate_lai
  integer(c_int),bind(C,name="clmupdate_lai_params") :: clmupdate_lai_params
  real(c_double),bind(C,name="clmupdate_lai_incr_w") :: clmupdate_lai_incr_w

  ! Yorck
  integer(c_int),bind(C,name="clmupdate_tws") :: clmupdate_tws
  integer(c_int),bind(C,name="exclude_greenland") :: exclude_greenland
  integer, dimension(1:5) :: clm_varsize_tws
  real(r8),bind(C,name="max_inc") :: max_inc
  integer(c_int),bind(C,name="TWS_smoother") :: TWS_smoother
  integer(c_int),bind(C,name="state_setup") :: state_setup
  integer(c_int),bind(C,name="set_zero_start") :: set_zero_start
  integer, allocatable :: num_layer(:)
  integer, allocatable :: num_layer_columns(:)
  integer :: num_hactiveg, num_hactivec

  integer, allocatable :: hactiveg_levels(:,:)     ! hydrolocial active filter for all levels (gridcell)
  integer, allocatable :: hactivec_levels(:,:)     ! hydrolocial active filter for all levels (column)
  integer, allocatable :: gridcell_state(:)

  logical :: first_cycle = .TRUE.
  logical :: use_omi_model = .FALSE.

  ! Mode-3 state layout: 3 carbon pools; +2 PFT params when update_lai_params=1
  ! (kmax is private/protected in PhotosynthesisMod and cannot be updated from PDAF)
  integer, parameter :: lai_mode3_npool = 3
  integer, parameter :: lai_mode3_nparam = 2   ! slatop, medlynslope


  ! OMI --> I want to update the observation type after each observation comes in.
  ! problem: the observation type is updated before the update of the assimilation
  ! this causes the wrong observation type to be used in the update

  ! idea: the state vector is newly initilized for each assimilation time step for the current variable
  ! use a new variable, e.g. obs_type_update_tws, and make it equal to clmupdate_tws at the beginning of the initialization
  ! this variable is then used in the update of the state vector

  integer :: obs_type_update_swc = 0
  integer :: obs_type_update_tws = 0
  integer :: obs_type_update_T = 0
  integer :: obs_type_update_texture = 0



  ! end Yorck

#endif
  integer(c_int),bind(C,name="clmprint_et")       :: clmprint_et
  integer(c_int),bind(C,name="clmprint_inc")       :: clmprint_inc
  integer(c_int),bind(C,name="clmstatevec_allcol")       :: clmstatevec_allcol
  integer(c_int),bind(C,name="clmstatevec_colmean")       :: clmstatevec_colmean
  integer(c_int),bind(C,name="clmstatevec_only_active")  :: clmstatevec_only_active
  integer(c_int),bind(C,name="clmstatevec_max_layer")  :: clmstatevec_max_layer
  integer(c_int),bind(C,name="clmt_printensemble")       :: clmt_printensemble
  integer(c_int),bind(C,name="clmwatmin_switch")         :: clmwatmin_switch
  integer(c_int),bind(C,name="clmswc_mask_snow")            :: clmswc_mask_snow
  real(c_double),bind(C,name="clmcrns_bd")      :: clmcrns_bd

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
                                ! (currently not used for eclm)
  logical :: newgridcell        !only eclm

  ! Arrays for LAI Assimilation
  real(r8),allocatable :: tlai(:)
  real(r8),allocatable :: tlai_orig(:)
  real(r8),allocatable :: lai_wtsum2_gc(:)

  contains

#if defined CLMSA
  subroutine define_clm_statevec(mype)
    use decompMod , only : get_proc_bounds
    use clm_varpar   , only : nlevsoi
    use PatchType  , only : patch
    use pftconMod  , only : noveg
    use clm_varcon, only: set_averaging_to_zero
    use PDAF_interfaces_module, only: PDAF_reset_dim_p

    implicit none

    integer,intent(in) :: mype

    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: i, cc


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

    clm_statevecsize = 0
    clm_varsize      = 0

    ! check which observation type should be assimilated and design the state vector accordingly
    ! I will introduce functions for each observation type so that other types can easily be included

    ! make update variable equal to the clmupdate variable
    obs_type_update_swc = clmupdate_swc
    obs_type_update_tws = clmupdate_tws
    obs_type_update_T   = clmupdate_T
    obs_type_update_texture = clmupdate_texture

    ! soil water content observations - case 1
    if(clmupdate_swc==1) then
      call define_clm_statevec_swc(mype)
    end if

    ! soil water content observations - case 2
    if(clmupdate_swc==2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    end if

    ! texture observations - case 1
    if(clmupdate_texture==1) then
      clm_statevecsize = clm_statevecsize + 2*((endg-begg+1)*nlevsoi)
    end if

    ! texture observations - case 2
    if(clmupdate_texture==2) then
      clm_statevecsize = clm_statevecsize + 3*((endg-begg+1)*nlevsoi)
    end if

    ! LAI assimilation
    if(clmupdate_lai==1) then
        ! Allocate array to store total leaf area index per patch
        if (allocated(tlai)) deallocate(tlai)
        allocate(tlai(clm_begp:clm_endp))
        if (allocated(tlai_orig)) deallocate(tlai_orig)
        allocate(tlai_orig(clm_begp:clm_endp))
        if (allocated(lai_wtsum2_gc)) deallocate(lai_wtsum2_gc)
        allocate(lai_wtsum2_gc(clm_begg:clm_endg))

        ! #patches values per grid-cell ! clm_varsize not accurate due to variable #patches/grid cell
        ! set sizes
        clm_varsize = 1 ! each grid cell gets 1 value of LAI
        ! Currently no combination of SWC and LAI DA

        clm_statevecsize = (clm_endg - clm_begg + 1) ! statevector is the size of the num of gridcells
        !(clm_endp-clm_begp+1) ! So like this if lai is set it takes priority

        if (clmupdate_lai_params==1) then ! add parameters to the state vector.
            ! clm_varsize = (clm_endg - clm_begg+1) ! Currently no combination of SWC and LAI DA
            ! clm_statevecsize = (clm_endg - clm_begg+1)*1 ! 1 parameters
            clm_varsize      =  (clm_endp-clm_begp+1)+1 ! Currently no combination of SWC and LAI DA
            clm_statevecsize =  (clm_endp-clm_begp+1)+1 ! 1 parameters
            ! So like this if lai is set it takes priority
        end if

        ! if (clmupdate_lai_params==3) then
        !     clm_varsize      =  (clm_endp-clm_begp+1)
        !     clm_statevecsize =  3*clm_varsize
        ! end if

    end if

    ! for obs_op version allocate additional helper index array
    if(clmupdate_lai==2) then
      clm_varsize      =  (clm_endp-clm_begp+1) ! equals number of patches
      clm_statevecsize =  3*clm_varsize ! 3 var/params per patch

      if (allocated(clm_patch2gc)) deallocate(clm_patch2gc)
      allocate(clm_patch2gc(clm_varsize))
      if (allocated(clm_patchwt)) deallocate(clm_patchwt)
      allocate(clm_patchwt(clm_varsize))

    endif

    ! LAI assimilation mode 3: leafc, livestemc, deadstemc per vegetated patch
    if(clmupdate_lai==3) then
      clm_varsize = 0
      do i=clm_begp,clm_endp
        if (patch%itype(i) /= noveg) clm_varsize = clm_varsize + 1
      end do
      clm_statevecsize = lai_mode3_npool*clm_varsize
      if (clmupdate_lai_params == 1) then
        clm_statevecsize = (lai_mode3_npool + lai_mode3_nparam)*clm_varsize
      end if

      if (allocated(clm_patch2gc)) deallocate(clm_patch2gc)
      allocate(clm_patch2gc(clm_varsize))
      if (allocated(clm_patchwt)) deallocate(clm_patchwt)
      allocate(clm_patchwt(clm_varsize))
      if (allocated(clm_lai_patch_idx)) deallocate(clm_lai_patch_idx)
      allocate(clm_lai_patch_idx(clm_varsize))
      if (allocated(clm_lai_patch_itype)) deallocate(clm_lai_patch_itype)
      allocate(clm_lai_patch_itype(clm_varsize))

      cc = 1
      do i=clm_begp,clm_endp
        if (patch%itype(i) /= noveg) then
          clm_lai_patch_idx(cc) = i
          clm_lai_patch_itype(cc) = patch%itype(i)
          clm_patch2gc(cc) = real(patch%gridcell(i), r8)
          clm_patchwt(cc) = patch%wtgcell(i)
          cc = cc + 1
        end if
      end do
    endif

    !hcp LST DA
    if(clmupdate_T==1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    end if
    !end hcp

    ! TWS observations
    if (clmupdate_tws==1) then
      call define_clm_statevec_tws(mype)
    end if

    !

      ! Include your own state vector definitions here for different variables (ET, LAI, etc.)

    !

    ! reset PDAF dimensions for multivariate assimilation, only not for first call as PDAF did not initalize yet
    ! non-first-cycle is only called in OMI simulations from clm_advance
    if (.not. first_cycle) then
        call PDAF_reset_dim_p(clm_statevecsize,ierror)
    end if

    if (first_cycle) then
      ! possibility to assimilate GRACE not in the first month -->
      ! enkfpf.par file has information set_zero_start where the
      ! running average should be resetted
      !
      ! This is usually one month prior to the first GRACE
      ! observation. If it is not included in the file, it is resetted
      ! when the first GRACE observation is assimilated.
      !
      ! Afterwards, the normal set_zero information inside the
      ! observation file is used (see next_observation_pdaf for
      ! details).
      if (set_zero_start/=0) then
        set_averaging_to_zero = set_zero_start
      end if
    end if

    first_cycle = .FALSE.


#ifdef PDAF_DEBUG
    ! Debug output of clm_statevecsize
    WRITE(*, "(a,x,a,i5,x,a,i10)") "TSMP-PDAF-debug", "mype(w)=", mype, "define_clm_statevec: clm_statevecsize=", clm_statevecsize
#endif

    IF (allocated(clm_statevec)) deallocate(clm_statevec)
    if ((clmupdate_swc/=0) .or. (clmupdate_T/=0) .or. (clmupdate_texture/=0) .or. (clmupdate_tws/=0) .or. (clmupdate_lai/=0)) then
      !hcp added condition
      allocate(clm_statevec(clm_statevecsize))
    end if


    ! Allocate statevector-duplicate for saving the forecast state
    ! vector. The forecast values may then be used during the
    ! update-phase, when mapping analysis increments back to CLM (used
    ! in SWC DA column-mean mode or LAI DA mode).
    IF (allocated(clm_statevec_orig)) deallocate(clm_statevec_orig)
    if ( (clmupdate_swc/=0 .and. clmstatevec_colmean/=0) .or. clmupdate_tws/=0  .or. (clmupdate_lai/=0)) then
      allocate(clm_statevec_orig(clm_statevecsize))
    end if

    !write(*,*) 'clm_paramsize is ',clm_paramsize
    if (allocated(clm_paramarr)) deallocate(clm_paramarr)         !hcp
    if ((clmupdate_T/=0)) then  !hcp
      error stop "Not implemented clmupdate_T.NE.0"
    end if

    if (allocated(gridcell_state)) deallocate(gridcell_state)
    allocate(gridcell_state(clm_statevecsize))


  end subroutine define_clm_statevec





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
          state_clm2pdaf_p(c,i) = ispval
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
          write(*,*) "cc: ", cc
          error stop "state_pdaf2clm_c_p not set at cc"
        end if
        if(state_pdaf2clm_j_p(cc) == ispval) then
          write(*,*) "cc: ", cc
          error stop "state_pdaf2clm_j_p not set at cc"
        end if
#endif
      end do

  end subroutine define_clm_statevec_swc

  !> @author Yorck Ewerdwalbesloh, Anne Springer
  !> @date 29.10.2025
  !> @brief define the state vector for TWS assimilation
      !> @param[in] mype MPI rank
  subroutine define_clm_statevec_tws(mype)
    use shr_kind_mod, only: r8 => shr_kind_r8
    use decompMod , only : get_proc_bounds
    use clm_varpar   , only : nlevsoi
    use ColumnType , only : col
    use PatchType, only: patch
    use GridcellType, only: grc

    implicit none

    integer,intent(in) :: mype

    integer :: i
    integer :: j
    integer :: c
    integer :: g
    integer :: cc

    integer :: fa, fg

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
    clm_begc     = begc
    clm_endc     = endc
    clm_begp     = begp
    clm_endp     = endp


    ! first we build a filter to determine which columns are active / are not active
    ! we also build a gridcell filter for gridcell averges

    num_hactiveg = 0
    num_hactivec = 0

    if (allocated(found)) deallocate(found)
    allocate(found(clm_begg:clm_endg))
    found(clm_begg:clm_endg) = .false.

    if (allocated(num_layer)) deallocate(num_layer)
    allocate(num_layer(1:nlevsoi))
    num_layer(1:nlevsoi) = 0

    if (allocated(num_layer_columns)) deallocate(num_layer_columns)
    allocate(num_layer_columns(1:nlevsoi))
    num_layer_columns(1:nlevsoi) = 0

    do c = clm_begc, clm_endc ! find hydrological active cells

      g = col%gridcell(c) ! gridcell of column

      ! greenland can be excluded from the statevector
      if ((exclude_greenland==0) .or. &
          (.not.(lon(g)<330 .and. lon(g)>180 .and. lat(g)>55))) then

        ! if the column is hydrologically active, add it, if the corresponding
        ! gridcell is not found before, add also the gridcell
        if (col%hydrologically_active(c)) then

          if (.not. found(g)) then ! if the gridcell is not found before

            found(g) = .true.

            do j = 1,nlevsoi ! add gridcell for each layer until bedrock
              ! get number in layers

              if (j<=col%nbedrock(c)) then
                num_layer(j) = num_layer(j) + 1
              end if

            end do

            num_hactiveg = num_hactiveg + 1

          end if

          do j = 1,nlevsoi ! add column for each layer until bedrock
            ! get number in layers

            if (j<=col%nbedrock(c)) then
              num_layer_columns(j) = num_layer_columns(j) + 1
            end if

          end do

          num_hactivec = num_hactivec + 1

        end if
      end if
    end do


    ! allocate matrices that will hold the index of each hydrologically active component (gridcell, column, patch) in each depth
    if (allocated(hactiveg_levels)) deallocate(hactiveg_levels)
    if (allocated(hactivec_levels)) deallocate(hactivec_levels)
    allocate(hactiveg_levels(1:num_hactiveg,1:nlevsoi))
    allocate(hactivec_levels(1:num_hactivec,1:nlevsoi))

    ! now we fill these things with the columns and gridcells so that we can access all active things later on

    do j = 1,nlevsoi

      ! has to be inside the for lopp, else, the hactiveg_levels is only filled for the first level
      found(clm_begg:clm_endg) = .false.
      fa = 0
      fg = 0

      do c = clm_begc, clm_endc
        g = col%gridcell(c) ! gridcell of column

        if ((exclude_greenland==0) .or. (.not.(lon(g)<330 .and. lon(g)>180 .and. lat(g)>55))) then

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


    if (allocated(found)) deallocate(found)

    ! now we have an array for the columns and gridcells of interest that we can
    ! use when we fill the statevector and distribute the update
    ! now lets find out the dimension of the state vector

    clm_varsize_tws(:) = 0

    clm_statevecsize = 0

    select case (state_setup)
    case(0)
      ! all compartments in state vector, liq and ice added up to not run into balancing errors due to different partioning
      ! in different ensemble members in these variables even if the total amount of water is similar

      do j = 1,nlevsoi
        clm_varsize_tws(1) = clm_varsize_tws(1) + num_layer(j)
        clm_statevecsize = clm_statevecsize + num_layer(j)

        clm_varsize_tws(2) = 0 ! as liq + ice is added up, we do not need another variable for ice
        clm_statevecsize = clm_statevecsize + 0
      end do

      ! snow
      clm_varsize_tws(3) = num_layer(1)
      clm_statevecsize = clm_statevecsize + num_layer(1)

      ! surface water
      clm_varsize_tws(4) = num_layer(1)
      clm_statevecsize = clm_statevecsize + num_layer(1)

      ! canopy water
      clm_varsize_tws(5) = 0
      clm_statevecsize = clm_statevecsize + 0

    case(1) ! TWS in statevector

      clm_varsize_tws(1) = num_layer(1)
      clm_statevecsize = clm_statevecsize + num_layer(1)
      clm_varsize_tws(2) = 0
      clm_varsize_tws(3) = 0
      clm_varsize_tws(4) = 0
      clm_varsize_tws(5) = 0

    case(2) ! snow and soil moisture aggregated over surface, root zone and deep soil moisture in state vector

      clm_varsize_tws(1) = num_layer(1)
      clm_statevecsize = clm_statevecsize + num_layer(1)

      clm_varsize_tws(2) = num_layer(4)
      clm_statevecsize = clm_statevecsize + num_layer(4)

      clm_varsize_tws(3) = num_layer(13)
      clm_statevecsize = clm_statevecsize + num_layer(13)

      ! one variable for snow
      clm_varsize_tws(4) = num_layer(1)
      clm_statevecsize = clm_statevecsize + num_layer(1)

    case default

      error stop "Unsupported state_setup"

    end select

  end subroutine define_clm_statevec_tws


  subroutine cleanup_clm_statevec()

    implicit none

    ! Deallocate arrays from `define_clm_statevec`
    IF (allocated(clm_statevec)) deallocate(clm_statevec)
    IF (allocated(state_pdaf2clm_c_p)) deallocate(state_pdaf2clm_c_p)
    IF (allocated(state_pdaf2clm_j_p)) deallocate(state_pdaf2clm_j_p)
    IF (allocated(state_clm2pdaf_p)) deallocate(state_clm2pdaf_p)
    IF (allocated(clm_statevec_orig)) deallocate(clm_statevec_orig)
    IF (allocated(clm_patch2gc)) deallocate(clm_patch2gc)
    IF (allocated(clm_patchwt)) deallocate(clm_patchwt)
    IF (allocated(clm_lai_patch_idx)) deallocate(clm_lai_patch_idx)
    IF (allocated(clm_lai_patch_itype)) deallocate(clm_lai_patch_itype)
    IF (allocated(lai_wtsum2_gc)) deallocate(lai_wtsum2_gc)
    IF (allocated(tlai)) deallocate(tlai)
    IF (allocated(tlai_orig)) deallocate(tlai_orig)

  end subroutine cleanup_clm_statevec

  subroutine set_clm_statevec(tstartcycle, mype)
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_instMod, only : bgc_vegetation_inst, canopystate_inst
    use clm_varpar   , only : nlevsoi
    ! use clm_varcon, only: nameg, namec
    ! use GetGlobalValuesMod, only: GetGlobalWrite
    use ColumnType , only : col
    use PatchType  , only : patch
    use pftconMod  , only : pftcon
    use shr_kind_mod, only: r8 => shr_kind_r8

    implicit none
    integer,intent(in) :: tstartcycle
    integer,intent(in) :: mype
    real(r8), pointer :: swc(:,:)
    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)
    real(r8), pointer :: leafc(:)
    real(r8), pointer :: livestemc(:)
    real(r8), pointer :: deadstemc(:)
    integer :: i, j, jj, g, c, cc, offset, igc, itype
    integer :: n_c
    character (len = 34) :: fn    !TSMP-PDAF: function name for state vector output
    character (len = 34) :: fn2    !TSMP-PDAF: function name for swc output
    character (len = 48) :: fn3    !TSMP-PDAF: function name for lai_leafc_* integrate output

    cc = 0
    offset = 0

    swc   => waterstate_inst%h2osoi_vol_col
    psand => soilstate_inst%cellsand_col
    pclay => soilstate_inst%cellclay_col
    porgm => soilstate_inst%cellorg_col
    leafc => bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_patch


#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle + 1 .OR. clmt_printensemble == -1) THEN

      IF(clmupdate_swc/=0) THEN
        ! TSMP-PDAF: Debug output of CLM swc
        WRITE(fn2, "(a,i5.5,a,i5.5,a)") "swcstate_", mype, ".integrate.", tstartcycle + 1, ".txt"
        OPEN(unit=71, file=fn2, action="write")
        WRITE (71,"(es22.15)") swc(:,:)
        CLOSE(71)
      END IF

    END IF
#endif

    if(clmupdate_swc==1) then
      call set_clm_statevec_swc
    end if

    ! calculate shift when CRP data are assimilated
    if(clmupdate_swc==2) then
      error stop "Not implemented clmupdate_swc.eq.2"
    end if

    !> @author  Lukas Strebel, Haojin Zhao
    !> @date    04.2026
    ! LAI assimilation
    ! Case 1: LAI all patches
    ! Case 1: Transformations in Set and Update functions
    if(clmupdate_lai==1) then
      clm_statevec(:) = 0._r8 ! reset statevector to 0 because we add to it for the average
      lai_wtsum2_gc(:) = 0._r8
      do i=clm_begp,clm_endp
        lai_wtsum2_gc(patch%gridcell(i)) = lai_wtsum2_gc(patch%gridcell(i)) &
          + patch%wtgcell(i)**2
      end do
       do i=clm_begp,clm_endp ! iterate through patches
           ! update the leaf area index based on leafC and SLA
           ! Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923.
           ! (slatop(ivt(p))*(exp(leafc(p)*dsladlai(ivt(p))) - 1._r8))/dsladlai(ivt(p))
           if (pftcon%dsladlai(patch%itype(i)) > 0._r8) then
             tlai(i) = (pftcon%slatop(patch%itype(i))*(exp(leafc(i)*pftcon%dsladlai(patch%itype(i))) - 1._r8)) &
               /pftcon%dsladlai(patch%itype(i))
           else
             tlai(i) = pftcon%slatop(patch%itype(i)) * leafc(i)
           end if
           tlai(i) = max(0._r8, tlai(i)) ! don't allow negative LAI
           tlai_orig(i) = tlai(i)
           ! Map CLM gridcell index (begg:endg) to PE-local state index (1:clm_statevecsize)
           igc = patch%gridcell(i) - clm_begg + 1
           if (igc >= 1 .and. igc <= clm_statevecsize) then
             clm_statevec(igc) = clm_statevec(igc) + patch%wtgcell(i)*tlai(i)
             write(*,'(A,I6,A,I6,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4)') &
                 "LAI1 DEBUG setvec ", i, ";", igc, ";", patch%wtgcell(i), ";", &
                 tlai(i), ";", patch%wtgcell(i)*tlai(i), ";", clm_statevec(igc), ";", &
                 canopystate_inst%tlai_patch(i), ";", leafc(i), ";", &
                 pftcon%dsladlai(patch%itype(i)), ";", pftcon%slatop(patch%itype(i))
           else
             write(*,*) 'DEBUG LAI : igc out of bounds in LAI mapping: ', igc, &
               ' valid range: 1-', clm_statevecsize, ' i=', i
           end if
       end do
       if (allocated(clm_statevec_orig)) then
         do igc = 1, min(clm_endg - clm_begg + 1, clm_statevecsize)
           clm_statevec_orig(igc) = clm_statevec(igc)
         end do
       end if

#ifdef PDAF_DEBUG
      IF(clmt_printensemble == tstartcycle + 1 .OR. clmt_printensemble == -1) THEN
        IF(clmupdate_lai/=0) THEN
          ! Same pattern as swcstate_*.integrate.*: one file, consecutive es22.15 blocks.
          ! Order: leafc(patch) tlai_calc(patch) tlai_clm(patch) L_gc(gridcell) lai_wtsum2_gc(gridcell)
          WRITE(fn3, "(a,i5.5,a,i5.5,a)") "lai_leafc_", mype, ".integrate.", tstartcycle + 1, ".txt"
          OPEN(unit=71, file=fn3, action="write", status="replace")
          WRITE(71,"(a)") "# blocks: leafc(patch) tlai_calc(patch) tlai_clm(patch) L_gc(gridcell) lai_wtsum2_sq(gridcell)"
          WRITE(71,"(es22.15)") leafc(clm_begp:clm_endp)
          WRITE(71,"(es22.15)") tlai(clm_begp:clm_endp)
          WRITE(71,"(es22.15)") canopystate_inst%tlai_patch(clm_begp:clm_endp)
          WRITE(71,"(es22.15)") clm_statevec(1:clm_statevecsize)
          WRITE(71,"(es22.15)") lai_wtsum2_gc(clm_begg:clm_endg)
          CLOSE(71)
        END IF
      END IF
#endif

    endif

    ! Case 2: statevec with leafc, slatop, dsladlai for transform in obs_op
    if (clmupdate_lai==2) then
      cc = 1
      do i=clm_begp,clm_endp
        clm_statevec(cc) = leafc(i)
        clm_statevec(cc + 1*clm_varsize) = pftcon%slatop(patch%itype(i))
        clm_statevec(cc + 2*clm_varsize) = pftcon%dsladlai(patch%itype(i))
        clm_patch2gc(cc) = patch%gridcell(i)
        clm_patchwt(cc) = patch%wtgcell(i)
        cc = cc + 1
      enddo
      write(*,*) 'DEBUG LAI : statevec ', clm_statevec(:)
      write(*,*) 'DEBUG LAI : patches ', clm_patch2gc(:), clm_patchwt(:)
    endif

    ! Case 3: leafc, livestemc, deadstemc per vegetated patch; H(x) in obs_op
    if (clmupdate_lai==3) then
      livestemc => bgc_vegetation_inst%cnveg_carbonstate_inst%livestemc_patch
      deadstemc => bgc_vegetation_inst%cnveg_carbonstate_inst%deadstemc_patch
      cc = 1
      do j = 1, clm_varsize
        i = clm_lai_patch_idx(j)
        clm_statevec(cc) = leafc(i)
        clm_statevec(cc + clm_varsize) = livestemc(i)
        clm_statevec(cc + 2*clm_varsize) = deadstemc(i)
        cc = cc + 1
      end do
      if (allocated(clm_statevec_orig)) then
        clm_statevec_orig(:) = clm_statevec(:)
      end if

      if (clmupdate_lai_params == 1) then
        cc = 1
        do j = 1, clm_varsize
          itype = clm_lai_patch_itype(j)
          clm_statevec(cc + lai_mode3_npool    *clm_varsize) = pftcon%slatop(itype)
          clm_statevec(cc + (lai_mode3_npool+1)*clm_varsize) = pftcon%medlynslope(itype)
          cc = cc + 1
        end do
      end if

#ifdef PDAF_DEBUG
      IF(clmt_printensemble == tstartcycle + 1 .OR. clmt_printensemble == -1) THEN
        call write_lai_da_mode3_debug('integrate', tstartcycle + 1, mype)
      END IF
#endif
    endif

    if (clmupdate_lai_params==1 .and. clmupdate_lai==1) then
      cc = 1
      do i=clm_begp,clm_endp
          clm_statevec(cc+1*clm_varsize+offset) = pftcon%slatop(patch%itype(i))
          cc = cc + 1
      end do
    endif

    ! if (clmupdate_lai_params==3) then
    !   cc = 1
    !   do i=clm_begp,clm_endp
    !     clm_statevec(cc+1*clm_varsize+offset) = pftcon%slatop(patch%itype(i))
    !     clm_statevec(cc+2*clm_varsize+offset) = params_inst%kmax(patch%itype(i),1)
    !     clm_statevec(cc+3*clm_varsize+offset) = params_inst%kmax(patch%itype(i),2)
    !     cc = cc + 1
    !   end do
    ! endif

    !hcp  LAI
    if(clmupdate_T==1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    end if
    !end hcp  LAI

    ! write average swc to state vector (CRP assimilation)
    if(clmupdate_swc==2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    end if

    ! write texture values to state vector (if desired)
    if(clmupdate_texture/=0) then
      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
          clm_statevec(cc+1*clm_varsize+offset) = psand(j,i)
          clm_statevec(cc+2*clm_varsize+offset) = pclay(j,i)
          if(clmupdate_texture==2) then
            !incl. organic matter values
            clm_statevec(cc+3*clm_varsize+offset) = porgm(j,i)
          end if
          cc = cc + 1
        end do
      end do
    end if

    if (clmupdate_tws==1) then
      call set_clm_statevec_tws
    end if

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle + 1 .OR. clmt_printensemble == -1) THEN
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

  !> @author Yorck Ewerdwalbesloh, Anne Springer
  !> @date 29.10.2025
  !> @brief set the state vector for TWS assimilation
  subroutine set_clm_statevec_tws()
    use clm_instMod, only : soilstate_inst, waterstate_inst
    use clm_varpar   , only : nlevsoi
    use ColumnType , only : col
    use shr_kind_mod, only: r8 => shr_kind_r8
    use PatchType, only: patch
    use GridcellType, only: grc
    implicit none

    integer :: j,g,cc,count,c,count_c
    integer :: n_c
    real(r8) :: avg_sum

    real(r8), pointer :: liq(:,:) ! liquid water (kg/m2)
    real(r8), pointer :: ice(:,:) ! ice lens (kg/m2)
    real(r8), pointer :: snow(:) ! snow water (mm)
    real(r8), pointer :: sfc(:) ! surface water
    real(r8), pointer :: can(:) ! canopy water
    real(r8), pointer :: TWS(:) ! TWS

    real(r8), pointer :: tws_state(:)
    real(r8), pointer :: liq_state(:,:)
    real(r8), pointer :: ice_state(:,:)
    real(r8), pointer :: snow_state(:)

    cc = 0

    tws_state => waterstate_inst%tws_state_before
    liq_state => waterstate_inst%h2osoi_liq_state_before
    ice_state => waterstate_inst%h2osoi_ice_state_before
    snow_state => waterstate_inst%h2osno_state_before

    select case (TWS_smoother)

    case(0) ! instantaneous values
      liq => waterstate_inst%h2osoi_liq_col
      ice => waterstate_inst%h2osoi_ice_col
      snow => waterstate_inst%h2osno_col
      sfc => waterstate_inst%h2osfc_col
      can => waterstate_inst%h2ocan_patch
      TWS => waterstate_inst%tws_hactive
    case(1) ! monthly means
      liq => waterstate_inst%h2osoi_liq_col_mean
      ice => waterstate_inst%h2osoi_ice_col_mean
      snow => waterstate_inst%h2osno_col_mean
      sfc => waterstate_inst%h2osfc_col_mean
      can => waterstate_inst%h2ocan_patch_mean
      TWS => waterstate_inst%tws_hactive_mean
    case default

      error stop "Unsupported TWS_smoother"

    end select

    select case (state_setup)
      case(0) ! all compartments in state vector
        cc = 1
        do j = 1,nlevsoi
          do count = 1, num_layer(j)
            g = hactiveg_levels(count,j)

            clm_statevec(cc) = 0.0
            n_c = 0

            do count_c = 1,num_layer_columns(j) ! get average over liq+ice for all columns in gridcell g
              c = hactivec_levels(count_c,j)
              if (g==col%gridcell(c)) then
                clm_statevec(cc) = clm_statevec(cc) + liq(c,j) + ice(c,j)
                n_c = n_c + 1
              end if
            end do

            clm_statevec(cc) = clm_statevec(cc) / real(n_c, r8)
            clm_statevec_orig(cc) = clm_statevec(cc)

            liq_state(g,j) = clm_statevec(cc)

            gridcell_state(cc) = g

            if (j==1) then ! for the first layer, we can fill the other compartments

              ! snow
              clm_statevec(cc+sum(clm_varsize_tws(1:2))) = 0.0
              n_c = 0

              do count_c = 1,num_layer_columns(j)
                c = hactivec_levels(count_c,j)
                if (g==col%gridcell(c)) then
                  clm_statevec(cc+sum(clm_varsize_tws(1:2))) = clm_statevec(cc+sum(clm_varsize_tws(1:2))) + snow(c)
                  n_c = n_c + 1
                end if
              end do

              clm_statevec(cc+sum(clm_varsize_tws(1:2))) = clm_statevec(cc+sum(clm_varsize_tws(1:2))) / real(n_c, r8)
              clm_statevec_orig(cc+sum(clm_varsize_tws(1:2))) = clm_statevec(cc+sum(clm_varsize_tws(1:2)))

              snow_state(g) = clm_statevec(cc+sum(clm_varsize_tws(1:2)))
              gridcell_state(cc+sum(clm_varsize_tws(1:2))) = g

              ! surface water
              clm_statevec(cc+sum(clm_varsize_tws(1:3))) = 0.0
              n_c = 0

              do count_c = 1,num_layer_columns(j)
                c = hactivec_levels(count_c,j)
                if (g==col%gridcell(c)) then
                  clm_statevec(cc+sum(clm_varsize_tws(1:3))) = clm_statevec(cc+sum(clm_varsize_tws(1:3))) + sfc(c)
                  n_c = n_c + 1
                end if

              end do

              clm_statevec(cc+sum(clm_varsize_tws(1:3))) = clm_statevec(cc+sum(clm_varsize_tws(1:3))) / real(n_c, r8)
              clm_statevec_orig(cc+sum(clm_varsize_tws(1:3))) = clm_statevec(cc+sum(clm_varsize_tws(1:3)))
              gridcell_state(cc+sum(clm_varsize_tws(1:3))) = g
            end if

            cc = cc+1
          end do
        end do

      case(1) ! TWS in state vector

        cc = 1

        do count = 1, num_layer(1)
          g = hactiveg_levels(count,1)
          clm_statevec(cc) = TWS(g)
          clm_statevec_orig(cc) = TWS(g)
          tws_state(g) = clm_statevec(cc)
          gridcell_state(cc) = g
          cc = cc+1
        end do

      case(2) ! snow and soil moisture aggregated over surface, root zone and deep soil moisture in state vector

        cc = 1

        ! surface soil moisture
        do count = 1, num_layer(1)
          clm_statevec(cc) = 0
          g = hactiveg_levels(count,1)

          do j = 1, 3 ! surface SM
            avg_sum = 0
            n_c = 0
            do count_c = 1,num_layer_columns(j)
              c = hactivec_levels(count_c,j)
              if (g==col%gridcell(c)) then
                avg_sum = avg_sum + liq(c,j) + ice(c,j)
                n_c = n_c+1
              end if

            end do

            if (n_c/=0) then

              clm_statevec(cc) = clm_statevec(cc) + avg_sum / real(n_c, r8)

            end if

          end do

          clm_statevec_orig(cc) = clm_statevec(cc)
          liq_state(g,1) = clm_statevec(cc)

          gridcell_state(cc) = g

          cc = cc + 1

        end do

        cc = 1

        ! root zone soil moisture
        do count = 1, num_layer(4)
          clm_statevec(cc+clm_varsize_tws(1)) = 0
          g = hactiveg_levels(count,4)
          do j = 4, 12
            avg_sum = 0
            n_c = 0
            do count_c = 1,num_layer_columns(j)
              c = hactivec_levels(count_c,j)
              if (g==col%gridcell(c)) then
                avg_sum = avg_sum + liq(c,j) + ice(c,j)
                n_c = n_c+1
              end if
            end do

            if (n_c/=0) then

              clm_statevec(cc+clm_varsize_tws(1)) = clm_statevec(cc+clm_varsize_tws(1)) + avg_sum / real(n_c, r8)

            end if

          end do

          clm_statevec_orig(cc+clm_varsize_tws(1)) = clm_statevec(cc+clm_varsize_tws(1))
          liq_state(g,2) = clm_statevec(cc+clm_varsize_tws(1))
          gridcell_state(cc+clm_varsize_tws(1)) = g

          cc = cc + 1

        end do

        cc = 1

        ! deep soil moisture
        do count = 1, num_layer(13)
          clm_statevec(cc+sum(clm_varsize_tws(1:2))) = 0
          g = hactiveg_levels(count,13)
          do j = 13, nlevsoi
            avg_sum = 0
            n_c = 0
            do count_c = 1,num_layer_columns(j)
              c = hactivec_levels(count_c,j)
              if (g==col%gridcell(c)) then
                avg_sum = avg_sum + liq(c,j) + ice(c,j)
                n_c = n_c+1
              end if
            end do

            if (n_c/=0) then

              clm_statevec(cc+sum(clm_varsize_tws(1:2))) = clm_statevec(cc+sum(clm_varsize_tws(1:2))) + avg_sum / real(n_c, r8)

            end if

          end do

          clm_statevec_orig(cc+sum(clm_varsize_tws(1:2))) = clm_statevec(cc+sum(clm_varsize_tws(1:2)))
          liq_state(g,3) = clm_statevec(cc+sum(clm_varsize_tws(1:2)))
          gridcell_state(cc+sum(clm_varsize_tws(1:2))) = g

          cc = cc + 1

        end do

        cc = 1

        ! snow
        do count = 1, num_layer(1)

          g = hactiveg_levels(count,1)

          clm_statevec(cc+sum(clm_varsize_tws(1:3))) = 0
          n_c = 0
          do count_c = 1,num_layer_columns(1)
            c = hactivec_levels(count_c,1)
            if (g==col%gridcell(c)) then
              clm_statevec(cc+sum(clm_varsize_tws(1:3))) = clm_statevec(cc+sum(clm_varsize_tws(1:3))) + snow(c)
              n_c = n_c+1
            end if
          end do

          clm_statevec(cc+sum(clm_varsize_tws(1:3))) = clm_statevec(cc+sum(clm_varsize_tws(1:3))) / real(n_c, r8)
          clm_statevec_orig(cc+sum(clm_varsize_tws(1:3))) = clm_statevec(cc+sum(clm_varsize_tws(1:3)))

          snow_state(g) = clm_statevec(cc+sum(clm_varsize_tws(1:3)))
          gridcell_state(cc+sum(clm_varsize_tws(1:3))) = g

          cc = cc+1

        end do

      case default

      error stop "Unsupported state_setup"

      end select


  end subroutine set_clm_statevec_tws




  subroutine lai_patch_tlai_blend(tlai_fc, lai_gc_fc, lai_gc_an, tlai_out)
    use shr_kind_mod , only : r8 => shr_kind_r8

    implicit none

    real(r8), intent(in) :: tlai_fc
    real(r8), intent(in) :: lai_gc_fc
    real(r8), intent(in) :: lai_gc_an
    real(r8), intent(out) :: tlai_out

    real(r8) :: incr_lai_gc
    real(r8) :: tlai_add
    real(r8) :: tlai_mult
    real(r8) :: w_add

    incr_lai_gc = lai_gc_an - lai_gc_fc
    tlai_add = tlai_fc + incr_lai_gc
    if (abs(lai_gc_fc) > 1.0e-12_r8) then
      tlai_mult = tlai_fc * (lai_gc_an / lai_gc_fc)
    else
      tlai_mult = tlai_add
    end if
    w_add = real(clmupdate_lai_incr_w, r8)
    tlai_out = w_add * tlai_add + (1._r8 - w_add) * tlai_mult
    tlai_out = max(0._r8, tlai_out)

  end subroutine lai_patch_tlai_blend


  subroutine write_lai_da_mode3_debug(phase, cycle, mype)
    use clm_instMod, only : bgc_vegetation_inst, canopystate_inst

    implicit none

    character(len=*), intent(in) :: phase
    integer, intent(in) :: cycle
    integer, intent(in) :: mype

    real(r8), pointer :: leafc(:)
    real(r8), pointer :: livestemc(:)
    real(r8), pointer :: deadstemc(:)
    real(r8), pointer :: leafn(:)
    real(r8), pointer :: livestemn(:)
    real(r8), pointer :: deadstemn(:)

    character(len=48) :: fn
    integer :: j, i

    IF(clmt_printensemble /= cycle .AND. clmt_printensemble /= -1) RETURN

    leafc => bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_patch
    livestemc => bgc_vegetation_inst%cnveg_carbonstate_inst%livestemc_patch
    deadstemc => bgc_vegetation_inst%cnveg_carbonstate_inst%deadstemc_patch
    leafn => bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_patch
    livestemn => bgc_vegetation_inst%cnveg_nitrogenstate_inst%livestemn_patch
    deadstemn => bgc_vegetation_inst%cnveg_nitrogenstate_inst%deadstemn_patch

    WRITE(fn, "(a,i5.5,a,a,a,i5.5,a)") "lai_da_", mype, ".", trim(phase), ".", cycle, ".txt"
    OPEN(unit=71, file=fn, action="write", status="replace")
    WRITE(71,"(a)") "# mode3 blocks: clm_statevec clm_lai_patch_idx clm_lai_patch_itype clm_patch2gc clm_patchwt"
    IF (clmupdate_lai_params == 1) THEN
      WRITE(71,"(a)") "# param blocks (if update_lai_params=1): slatop medlynslope per slot"
    END IF
    WRITE(71,"(a)") "#            leafc livestemc deadstemc leafn livestemn deadstemn tlai_clm"
    DO i = 1, clm_statevecsize
      WRITE(71,"(es22.15)") clm_statevec(i)
    END DO
    DO j = 1, clm_varsize
      WRITE(71,"(i12)") clm_lai_patch_idx(j)
    END DO
    DO j = 1, clm_varsize
      WRITE(71,"(i12)") clm_lai_patch_itype(j)
    END DO
    DO j = 1, clm_varsize
      WRITE(71,"(es22.15)") clm_patch2gc(j)
    END DO
    DO j = 1, clm_varsize
      WRITE(71,"(es22.15)") clm_patchwt(j)
    END DO
    DO j = 1, clm_varsize
      i = clm_lai_patch_idx(j)
      WRITE(71,"(es22.15)") leafc(i)
    END DO
    DO j = 1, clm_varsize
      i = clm_lai_patch_idx(j)
      WRITE(71,"(es22.15)") livestemc(i)
    END DO
    DO j = 1, clm_varsize
      i = clm_lai_patch_idx(j)
      WRITE(71,"(es22.15)") deadstemc(i)
    END DO
    DO j = 1, clm_varsize
      i = clm_lai_patch_idx(j)
      WRITE(71,"(es22.15)") leafn(i)
    END DO
    DO j = 1, clm_varsize
      i = clm_lai_patch_idx(j)
      WRITE(71,"(es22.15)") livestemn(i)
    END DO
    DO j = 1, clm_varsize
      i = clm_lai_patch_idx(j)
      WRITE(71,"(es22.15)") deadstemn(i)
    END DO
    DO j = 1, clm_varsize
      i = clm_lai_patch_idx(j)
      WRITE(71,"(es22.15)") canopystate_inst%tlai_patch(i)
    END DO
    CLOSE(71)

  end subroutine write_lai_da_mode3_debug


  subroutine update_clm(tstartcycle, mype) bind(C,name="update_clm")
    use clm_time_manager  , only : update_DA_nstep
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_instMod, only : waterstate_inst
    use clm_instMod, only : bgc_vegetation_inst, canopystate_inst
    use pftconMod       , only : pftcon
    use PatchType       , only : patch

    implicit none

    integer,intent(in) :: tstartcycle
    integer,intent(in) :: mype

    real(r8), pointer :: h2osoi_liq(:,:)  ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)
    real(r8), pointer :: leafc(:)         ! leaf carbon of patch
    real(r8), pointer :: leafn(:)         ! leaf nitrogen of patch
    real(r8), pointer :: livestemc(:)
    real(r8), pointer :: deadstemc(:)
    real(r8), pointer :: livestemn(:)
    real(r8), pointer :: deadstemn(:)

    real(r8) :: lai_gc_an
    real(r8) :: lai_gc_fc

    integer :: i
    integer :: j
    integer :: itype
    integer :: cc
    integer :: offset
    integer :: igc
    character (len = 31) :: fn    !TSMP-PDAF: function name for state vector output
    character (len = 31) :: fn2    !TSMP-PDAF: function name for state vector output
    character (len = 32) :: fn5    !TSMP-PDAF: function name for state vector output
    character (len = 32) :: fn6    !TSMP-PDAF: function name for state vector output
    character (len = 32) :: fn7    !TSMP-PDAF: function name for state vector output

    logical :: swc_zero_before_update

    swc_zero_before_update = .false.

    cc = 1
    offset = 0

#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble == -1) THEN
      ! TSMP-PDAF: For debug runs, output the state vector in files
      WRITE(fn, "(a,i5.5,a,i5.5,a)") "clmstate_", mype, ".update.", tstartcycle, ".txt"
      OPEN(unit=71, file=fn, action="write")
      DO i = 1, clm_statevecsize
        WRITE (71,"(es22.15)") clm_statevec(i)
      END DO
      CLOSE(71)
    END IF
#endif

    h2osoi_liq    => waterstate_inst%h2osoi_liq_col
    h2osoi_ice    => waterstate_inst%h2osoi_ice_col

    leafc => bgc_vegetation_inst%cnveg_carbonstate_inst%leafc_patch
    leafn => bgc_vegetation_inst%cnveg_nitrogenstate_inst%leafn_patch
    livestemc => bgc_vegetation_inst%cnveg_carbonstate_inst%livestemc_patch
    deadstemc => bgc_vegetation_inst%cnveg_carbonstate_inst%deadstemc_patch
    livestemn => bgc_vegetation_inst%cnveg_nitrogenstate_inst%livestemn_patch
    deadstemn => bgc_vegetation_inst%cnveg_nitrogenstate_inst%deadstemn_patch



#ifdef PDAF_DEBUG
    IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble == -1) THEN

      IF(obs_type_update_swc/=0) THEN
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

      IF(clmupdate_lai/=0) THEN
        ! TSMP-PDAF: For debug runs, output the state vector in files
        WRITE(fn5, "(a,i5.5,a,i5.5,a)") "leafc", mype, ".bef_up.", tstartcycle, ".txt"
        OPEN(unit=71, file=fn5, action="write")
        WRITE (71,"(es22.15)") leafc(:)
        CLOSE(71)
      END IF

    END IF
#endif

    ! calculate shift when CRP data are assimilated
    if(obs_type_update_swc==2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
    end if

    ! CLM5: Update the Data Assimulation time-step to the current time
    ! step, since DA has been done. Used by CLM5 to skip BalanceChecks
    ! directly after the DA step.
    call update_DA_nstep()

    ! write updated swc back to CLM
    if(obs_type_update_swc/=0) then
      call update_clm_swc(tstartcycle, mype)
    end if

    !hcp: TG, TV
    if(obs_type_update_T==1) then
      error stop "Not implemented: clmupdate_T.eq.1"
    end if
    ! end hcp TG, TV

    ! write updated texture back to CLM
    if(obs_type_update_texture/=0) then
      call update_clm_texture(tstartcycle, mype)
    end if

    if (obs_type_update_tws==1) then
      call clm_update_tws
    end if



    ! LAI assimilation:
    ! Case 1: gridcell-mean LAI state; map analysis back to patches via
    ! additive/multiplicative blend (CLM:update_lai_incr_w).
    if(clmupdate_lai==1) then
      do i=clm_begp,clm_endp
        igc = patch%gridcell(i) - clm_begg + 1
        if (igc >= 1 .and. igc <= clm_statevecsize .and. allocated(clm_statevec_orig)) then
          lai_gc_fc = clm_statevec_orig(igc)
          lai_gc_an = clm_statevec(igc)
          call lai_patch_tlai_blend(tlai_orig(i), lai_gc_fc, lai_gc_an, tlai(i))
        else
          tlai(i) = 0._r8
        endif
        ! update leafc based on Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923
        ! reformulate the equation to solve for leafc
        ! => leafc(p) = log(((tlai(p) * dsladlai(ivt(p)))/slatop(ivt(p))) + 1._r8) / dsladlai(ivt(p))
        if (tlai(i) > 0._r8) then ! invalid log protection
          if (pftcon%dsladlai(patch%itype(i)) > 0._r8) then
            if (pftcon%slatop(patch%itype(i)) > 0._r8) then
              leafc(i) = log(((tlai(i) * pftcon%dsladlai(patch%itype(i))) &
                / pftcon%slatop(patch%itype(i))) + 1.0_r8) &
                / pftcon%dsladlai(patch%itype(i))
            else
              leafc(i) = 0._r8
            endif
          else
              ! Forward: tlai = slatop * leafc  (dsladlai == 0) => leafc = tlai / slatop
              if (pftcon%slatop(patch%itype(i)) > 0._r8) then
                leafc(i) = tlai(i) / pftcon%slatop(patch%itype(i))
              else
                leafc(i) = 0._r8
              end if
          endif
        else
          leafc(i) = 0._r8
        endif
        leafc(i) = max(0._r8, leafc(i)) ! Don't allow negative leaf carbon
        if (pftcon%leafcn(patch%itype(i)) > 0._r8) then
          leafn(i) = leafc(i) / pftcon%leafcn(patch%itype(i))
        else
          leafn(i) = 0._r8
        endif
        ! DEBUG ONLY
        ! print *, "LAI1 DEBUG leafc,tlai(wt) AFTER: ", i, ";", leafc(i), ";", tlai(i)
      end do

#ifdef PDAF_DEBUG
        IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble == -1) THEN
          ! TSMP-PDAF: lai & leafc  (legacy column vector)
          WRITE(fn7, "(a,i5.5,a,i5.5,a)") "lai_leafc_", mype, ".update.", tstartcycle, ".txt"
          OPEN(unit=71, file=fn7, action="write", status="replace")
          WRITE(71,"(a)") "# patch  leafc_aft  tlai_patch_aft"
          DO i = clm_begp, clm_endp
            WRITE(71,"(i8,2es22.15)") i, leafc(i), tlai(i)
          END DO
          CLOSE(71)
        END IF
#endif
    endif

    ! Case 2: statevec contains leafc, slatop, dsladlai
    ! update of diagnostic tlai left to be done by clm5 itself.
    if(clmupdate_lai==2) then
      cc = 1
      do i=clm_begp,clm_endp
        leafc(i) = clm_statevec(cc)
        leafc(i) = max(0._r8, leafc(i))
        if (pftcon%leafcn(patch%itype(i)) > 0._r8) then
          leafn(i) = leafc(i) / pftcon%leafcn(patch%itype(i))
        else
          leafn(i) = 0._r8
        end if
        if(clmupdate_lai_params==2) then
          pftcon%slatop(patch%itype(i)) = clm_statevec(cc+1*clm_varsize)
          pftcon%dsladlai(patch%itype(i)) = clm_statevec(cc+2*clm_varsize)
        endif

        cc = cc + 1
      enddo

#ifdef PDAF_DEBUG
        IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble < 0) THEN
          ! TSMP-PDAF: For debug runs, output the state vector in files
          WRITE(fn7, "(a,i5.5,a,i5.5,a)") "leafc_", mype, ".update.", tstartcycle, ".txt"
          OPEN(unit=71, file=fn7, action="write")
          WRITE (71,"(es22.15)") leafc(:)
          CLOSE(71)
        END IF
#endif

    endif

    ! Case 3: leafc, livestemc, deadstemc per vegetated patch
    if(clmupdate_lai==3) then
      cc = 1
      do j = 1, clm_varsize
        i = clm_lai_patch_idx(j)
        itype = clm_lai_patch_itype(j)
        leafc(i) = max(0._r8, clm_statevec(cc))
        livestemc(i) = max(0._r8, clm_statevec(cc + clm_varsize))
        deadstemc(i) = max(0._r8, clm_statevec(cc + 2*clm_varsize))
        if (pftcon%leafcn(itype) > 0._r8) then
          leafn(i) = leafc(i) / pftcon%leafcn(itype)
        else
          leafn(i) = 0._r8
        end if
        if (pftcon%livewdcn(itype) > 0._r8) then
          livestemn(i) = livestemc(i) / pftcon%livewdcn(itype)
        else
          livestemn(i) = 0._r8
        end if
        if (pftcon%deadwdcn(itype) > 0._r8) then
          deadstemn(i) = deadstemc(i) / pftcon%deadwdcn(itype)
        else
          deadstemn(i) = 0._r8
        end if
        cc = cc + 1
      end do

      if (clmupdate_lai_params == 1) then
        cc = 1
        do j = 1, clm_varsize
          itype = clm_lai_patch_itype(j)
          pftcon%slatop(itype) = clm_statevec(cc + lai_mode3_npool    *clm_varsize)
          pftcon%medlynslope(itype) = clm_statevec(cc + (lai_mode3_npool+1)*clm_varsize)
          cc = cc + 1
        end do
      end if

#ifdef PDAF_DEBUG
      IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble == -1) THEN
        call write_lai_da_mode3_debug('update', tstartcycle, mype)
      END IF
#endif
    endif

    if (clmupdate_lai_params==1 .and. clmupdate_lai==1) then
      cc = 1
      do i=clm_begp,clm_endp
        pftcon%slatop(patch%itype(i)) = clm_statevec(cc+1*clm_varsize+offset)
        ! DEBUG ONLY
        print *, "LST DEBUG PARAM,i: ", pftcon%slatop(patch%itype(i)), ",", i
        cc = cc + 1
      end do
    endif

    ! if (clmupdate_lai_params==3) then
    !   cc = 1
    !   do i=clm_begp,clm_endp
    !     pftcon%slatop(patch%itype(i)) = clm_statevec(cc+1*clm_varsize+offset)
    !     params_inst%kmax(patch%itype(i),1) = clm_statevec(cc+2*clm_varsize+offset)
    !     params_inst%kmax(patch%itype(i),2) = clm_statevec(cc+3*clm_varsize+offset)

    !     ! DEBUG ONLY
    !     print *, "LST DEBUG PARAM slatop,i: ", pftcon%slatop(patch%itype(i)), ",", i
    !     print *, "LST DEBUG PARAM kmax sun,i: ", params_inst%kmax(patch%itype(i),1), ",", i
    !     print *, "LST DEBUG PARAM kmax shade,i: ", params_inst%kmax(patch%itype(i),2), ",", i

    !     cc = cc + 1
    !   end do

    ! end if

  end subroutine update_clm


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
    real(r8), pointer :: liq_inc(:,:), ice_inc(:,:), snow_inc(:)
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

    liq_inc => waterstate_inst%h2osoi_liq_col_inc
    ice_inc => waterstate_inst%h2osoi_ice_col_inc
    snow_inc => waterstate_inst%h2osno_col_inc


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

        do i = 1,nlevsoi
          do j = clm_begc,clm_endc

            liq_inc(j,i) = h2osoi_liq(j,i)
            ice_inc(j,i) = h2osoi_ice(j,i)

            if (i==1) then
              snow_inc(j) = waterstate_inst%h2osno_col(j)
            end if

          end do
        end do

        ! cc = 0
        do i=1,nlevsoi
          ! CLM3.5: iterate over grid cells
          ! CLM5.0: iterate over columns
          ! do j=clm_begg,clm_endg
            do j=clm_begc,clm_endc

              ! If snow is masked, update only, when snow depth is less than 1mm
              if( (clmswc_mask_snow == 0) .or. snow_depth(j) < 0.001 ) then
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
                end if

                if (ieee_is_nan(swc(j,i))) then
                  swc(j,i) = watmin_set
                  print *, "WARNING: swc at j,i is nan: ", j, i
                end if

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
        IF(clmt_printensemble == tstartcycle .OR. clmt_printensemble == -1) THEN

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

    do i = 1,nlevsoi
      do j=clm_begc,clm_endc

        liq_inc(j,i) = h2osoi_liq(j,i)-liq_inc(j,i)
        ice_inc(j,i) = h2osoi_ice(j,i)-ice_inc(j,i)

        if (i==1) then
          snow_inc(j) = waterstate_inst%h2osno_col(j)-snow_inc(j)
        end if

      end do
    end do

  end subroutine update_clm_swc


  subroutine update_clm_texture(tstartcycle, mype)
    use clm_varpar   , only : nlevsoi
    use shr_kind_mod , only : r8 => shr_kind_r8
    use clm_instMod, only : soilstate_inst

    implicit none

    integer,intent(in) :: tstartcycle
    integer,intent(in) :: mype

    integer :: i,j,cc,offset

    real(r8), pointer :: psand(:,:)
    real(r8), pointer :: pclay(:,:)
    real(r8), pointer :: porgm(:,:)

    cc = 0
    offset = 0

    psand   => soilstate_inst%cellsand_col
    pclay   => soilstate_inst%cellclay_col
    porgm   => soilstate_inst%cellorg_col

      cc = 1
      do i=1,nlevsoi
        do j=clm_begg,clm_endg
          psand(j,i) = clm_statevec(cc+1*clm_varsize+offset)
          pclay(j,i) = clm_statevec(cc+2*clm_varsize+offset)
          if(clmupdate_texture==2) then
            ! incl. organic matter
            porgm(j,i) = clm_statevec(cc+3*clm_varsize+offset)
          end if
          cc = cc + 1
        end do
      end do
      call clm_correct_texture
      call clm_texture_to_parameters

  end subroutine update_clm_texture


  !> @author Yorck Ewerdwalbesloh, Anne Springer
  !> @date 29.10.2025
  !> @brief update water storages from TWS assimilation
  subroutine clm_update_tws()

    use clm_varpar   , only : nlevsoi, nlevsno
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varcon, only: watmin, denh2o, denice, averaging_var
    use ColumnType, only : col
    use clm_instMod, only : soilstate_inst, waterstate_inst

    implicit none

    integer :: c, j

    real(r8), pointer :: liq(:,:), ice(:,:), swc(:,:)
    real(r8), pointer :: dz(:,:), zi(:,:), z(:,:)
    real(r8), pointer :: watsat(:,:), snow(:), snow_depth(:)
    integer, pointer  :: snl(:)

    real(r8), pointer :: liq_mean(:,:), ice_mean(:,:), snow_mean(:)
    real(r8), pointer :: liq_inc(:,:), ice_inc(:,:), snow_inc(:)
    real(r8), pointer :: TWS(:)

    ! Pointer declarations
    call assign_pointers()

    ! set averaging factor to zero
    averaging_var = 0

    ! Initialize increments
    call initialize_increments()

    select case (state_setup)
    case(0) ! all compartments in state vector
      call update_state_0()
    case(1) ! TWS in state vector
      call update_state_1()
    case(2) ! snow and soil moisture aggregated over surface, root zone and deep soil moisture in state vector
      call update_state_2()
    case default
      error stop "Unsupported state_setup"
    end select

    call finalize_increments()

    contains

    subroutine assign_pointers()

      liq         => waterstate_inst%h2osoi_liq_col
      ice         => waterstate_inst%h2osoi_ice_col
      swc         => waterstate_inst%h2osoi_vol_col
      dz          => col%dz
      zi          => col%zi
      z           => col%z
      watsat      => soilstate_inst%watsat_col
      snow        => waterstate_inst%h2osno_col
      snow_depth  => waterstate_inst%snow_depth_col
      snl         => col%snl

      select case (TWS_smoother)
      case(0)
        liq_mean => waterstate_inst%h2osoi_liq_col
        ice_mean => waterstate_inst%h2osoi_ice_col
        snow_mean => waterstate_inst%h2osno_col
      case default
        liq_mean => waterstate_inst%h2osoi_liq_col_mean
        ice_mean => waterstate_inst%h2osoi_ice_col_mean
        snow_mean => waterstate_inst%h2osno_col_mean
      end select

      liq_inc => waterstate_inst%h2osoi_liq_col_inc
      ice_inc => waterstate_inst%h2osoi_ice_col_inc
      snow_inc => waterstate_inst%h2osno_col_inc

      TWS => waterstate_inst%tws_hactive

    end subroutine assign_pointers

    subroutine initialize_increments()
      integer :: j, count, c
      ! set increment values for inc output file to old values, then they can be updated with new values
      liq_inc(:,:) = 0.0
      ice_inc(:,:) = 0.0
      snow_inc(:) = 0.0
      do j = 1,nlevsoi
        do count = 1,num_layer_columns(j)
          c = hactivec_levels(count,j)

          liq_inc(c,j) = liq(c,j)
          ice_inc(c,j) = ice(c,j)

          if (j==1) then
            snow_inc(c) = snow(c)
          end if
        end do
      end do
    end subroutine initialize_increments

    subroutine finalize_increments()
      integer :: j, count, c
      do j = 1,nlevsoi
        do count = 1,num_layer_columns(j)
          c = hactivec_levels(count,j)

          liq_inc(c,j) = liq(c,j)-liq_inc(c,j)
          ice_inc(c,j) = ice(c,j)-ice_inc(c,j)

          if (j==1) then
            snow_inc(c) = snow(c)-snow_inc(c)
          end if
        end do
      end do
    end subroutine finalize_increments




    subroutine update_state_0()
      real(r8) :: inc
      integer :: c, j, count, g, cc, count_c

      ! Update soil
      cc = 1
      do j = 1,nlevsoi
        do count = 1, num_layer(j)
          g = hactiveg_levels(count,j)
          inc = clm_statevec(cc)-clm_statevec_orig(cc)

          if (abs(inc)>1.e-10_r8) then
            do count_c = 1, num_layer_columns(j)
              c = hactivec_levels(count_c,j)
              if (col%gridcell(c)==g) then

                call update_soil_layer(c, j, inc, liq_mean(c,j), ice_mean(c,j), clm_statevec_orig(cc))

              end if
            end do
          end if
          cc = cc+1
        end do
      end do

      ! update snow
      cc = 1
      do count = 1, num_layer(1)
        g = hactiveg_levels(count,1)
        inc = clm_statevec(cc+sum(clm_varsize_tws(1:2))) - clm_statevec_orig(cc+sum(clm_varsize_tws(1:2)))
        if (abs(inc)>1.e-10_r8) then
          do count_c = 1, num_layer_columns(1)
            c = hactivec_levels(count_c,1)
            if (col%gridcell(c)==g) then

              if (snl(c) < 0) then
                call scale_snow(c, inc*snow_mean(c)/clm_statevec_orig(cc+sum(clm_varsize_tws(1:2))))
              end if

              if (snow(c) < 0._r8) then ! close snow layers if snow is negative
                call zero_snow_layers(c)
              end if

            end if
          end do
        end if
        cc = cc+1
      end do

    end subroutine update_state_0



    subroutine update_state_1()
      real(r8) :: inc
      integer :: c, j, count, g, cc, count_c

      cc = 1

      do count = 1, num_layer(1)
        g = hactiveg_levels(count,1)
        inc = clm_statevec(cc)-clm_statevec_orig(cc)

        if (abs(inc)>1.e-10_r8) then
          ! update soil water
          do j = 1,nlevsoi

            do count_c = 1,num_layer_columns(j)
              c = hactivec_levels(count_c,j)
              if (col%gridcell(c)==g) then
                call update_soil_layer(c, j, inc, liq_mean(c,j), ice_mean(c,j), clm_statevec_orig(cc))
              end if
            end do

          end do

          ! update snow

          do count_c = 1,num_layer_columns(1)
            c = hactivec_levels(count_c,1)
            if (col%gridcell(c)==g) then
              if (snl(c) < 0) then ! snow layers in the column
                call scale_snow(c, inc*snow_mean(c)/clm_statevec_orig(cc))
              end if
              if (snow(c) < 0._r8) then ! close snow layers if snow is negative
                call zero_snow_layers(c)
              end if
            end if
          end do
        end if

        cc = cc+1

      end do

    end subroutine update_state_1


    subroutine update_state_2()
      real(r8) :: inc
      integer :: c, j, count, g, cc, count_c

      ! surface soil moisture
      cc = 1
      do count = 1, num_layer(1)
        g = hactiveg_levels(count,1)
        inc = clm_statevec(cc)-clm_statevec_orig(cc)
        if (abs(inc)>1.e-10_r8) then
          do j = 1,3
            do count_c = 1,num_layer_columns(j)
              c = hactivec_levels(count_c,j)
              if (col%gridcell(c)==g) then

                call update_soil_layer(c, j, inc, liq_mean(c,j), ice_mean(c,j), clm_statevec_orig(cc))

              end if
            end do
          end do
        end if
        cc = cc+1
      end do

      ! root zone soil moisture
      cc = 1
      do count = 1, num_layer(4)
        g = hactiveg_levels(count,4)
        inc = clm_statevec(cc+clm_varsize_tws(1))-clm_statevec_orig(cc+clm_varsize_tws(1))
        if (abs(inc)>1.e-10_r8) then
          do j = 4,12
            do count_c = 1,num_layer_columns(j)
              c = hactivec_levels(count_c,j)
              if (col%gridcell(c)==g) then

                call update_soil_layer(c, j, inc, liq_mean(c,j), ice_mean(c,j), clm_statevec_orig(cc+clm_varsize_tws(1)))

              end if
            end do
          end do
        end if
        cc = cc+1
      end do

      ! deep soil moisture
      cc = 1
      do count = 1, num_layer(13)
        g = hactiveg_levels(count,13)
        inc = clm_statevec(cc+sum(clm_varsize_tws(1:2)))-clm_statevec_orig(cc+sum(clm_varsize_tws(1:2)))
        if (abs(inc)>1.e-10_r8) then
          do j = 13,nlevsoi
            do count_c = 1,num_layer_columns(j)
              c = hactivec_levels(count_c,j)
              if (col%gridcell(c)==g) then

                call update_soil_layer(c, j, inc, liq_mean(c,j), ice_mean(c,j), clm_statevec_orig(cc+sum(clm_varsize_tws(1:2))))

              end if
            end do
          end do
        end if
        cc = cc+1
      end do

      ! update snow
      cc = 1
      do count = 1, num_layer(1)
        g = hactiveg_levels(count,1)
        inc = clm_statevec(cc+sum(clm_varsize_tws(1:3))) - clm_statevec_orig(cc+sum(clm_varsize_tws(1:3)))
        if (abs(inc)>1.e-10_r8) then
          do count_c = 1, num_layer_columns(1)
            c = hactivec_levels(count_c,1)
            if (col%gridcell(c)==g) then

              if (snl(c) < 0) then
                call scale_snow(c, inc*snow_mean(c)/clm_statevec_orig(cc+sum(clm_varsize_tws(1:3))))
              end if

              if (snow(c) < 0._r8) then ! close snow layers if snow is negative
                call zero_snow_layers(c)
              end if

            end if
          end do
        end if
        cc = cc+1
      end do

    end subroutine update_state_2



    subroutine update_soil_layer(c, j, inc, liq_mean_cj, ice_mean_cj, vec_orig)
      real(r8), intent(in) :: inc, liq_mean_cj, ice_mean_cj, vec_orig
      integer, intent(in) :: c, j
      real(r8) :: inc_col, var_temp

      inc_col = inc * liq_mean_cj / vec_orig
      if (inc_col /= inc_col) inc_col = 0.0_r8
      if (abs(inc_col) > max_inc * liq(c,j)) inc_col = sign(max_inc * liq(c,j), inc_col)
      liq(c,j) = max(watmin, liq(c,j) + inc_col)

      inc_col = inc * ice_mean_cj / vec_orig
      if (inc_col /= inc_col) inc_col = 0.0_r8
      if (abs(inc_col) > max_inc * ice(c,j)) inc_col = sign(max_inc * ice(c,j), inc_col)
      ice(c,j) = max(0._r8, ice(c,j) + inc_col)

      swc(c,j) = liq(c,j)/(dz(c,j)*denh2o) + ice(c,j)/(dz(c,j)*denice)

      if (j > 1 .and. swc(c,j) - watsat(c,j) > 0.00001_r8) then
        var_temp = watsat(c,j) / swc(c,j)
        liq(c,j) = max(watmin, liq(c,j) * var_temp)
        ice(c,j) = max(0._r8, watsat(c,j)*dz(c,j)*denice - liq(c,j)*denice/denh2o)
        swc(c,j) = watsat(c,j)
        if (abs(ice(c,j)) < 1.e-10_r8) ice(c,j) = 0._r8
      end if
    end subroutine update_soil_layer


    subroutine scale_snow(c, inc)
      integer, intent(in) :: c
      real(r8), intent(in) :: inc
      real(r8) :: scale, inc_col
      integer :: j

      inc_col=inc
      if (inc_col /= inc_col) inc_col = 0.0

      if (abs(inc_col)>max_inc*snow(c)) inc_col = sign(max_inc*snow(c),inc_col)
      inc_col = snow(c) + inc_col
      scale = inc_col / snow(c)
      snow(c) = inc_col
      do j = 0, snl(c)+1, -1
        liq(c,j) = liq(c,j)*scale
        ice(c,j) = ice(c,j)*scale
        dz(c,j)  = dz(c,j)*scale
        zi(c,j)  = zi(c,j)*scale
        z(c,j)   = z(c,j)*scale
      end do

      zi(c,snl(c)) = zi(c,snl(c))*scale
      snow_depth(c) = snow_depth(c)*scale
    end subroutine scale_snow


    subroutine zero_snow_layers(c)
      integer, intent(in) :: c
      integer :: j

      if (snl(c) < 0) then
        do j = 0, snl(c)+1, -1
          liq(c,j) = 0.0_r8
          ice(c,j) = 1.e-8_r8
          dz(c,j) = 1.e-8_r8
          zi(c,j-1) = sum(dz(c,j:0)) * -1._r8
          if (j == 0) then
            z(c,j) = zi(c,j-1) / 2._r8
          else
            z(c,j) = sum(zi(c,j-1:j)) / 2._r8
          end if
        end do
      else
        liq(c,0) = 0.0_r8
        ice(c,0) = 1.e-8_r8
        dz(c,0) = 1.e-8_r8
        zi(c,-1) = dz(c,0) * -1._r8
        z(c,0) = zi(c,-1) / 2._r8
      end if

      snow_depth(c) = sum(dz(c,-nlevsno+1:0))
      snow(c) = sum(ice(c,-nlevsno+1:0))
    end subroutine zero_snow_layers

  end subroutine clm_update_tws

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

         if(sand<=0.0) sand = 1.0
         if(clay<=0.0) clay = 1.0

         ttot = sand + clay
         if(ttot>100) then
             sand = sand/ttot * 100.0
             clay = clay/ttot * 100.0
         end if

         pclay(c,lev) = clay
         psand(c,lev) = sand

         if(clmupdate_texture==2) then
           orgm = (porgm(c,lev) / CNParamsShareInst%organic_max) * 100.0
           if(orgm<=0.0) orgm = 0.0
           if(orgm>=100.0) orgm = 100.0
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

    real(r8), parameter           :: om_tkm         = 0.25_r8
    ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8), parameter           :: om_watsat_lake = 0.9_r8
    ! porosity of organic soil
    real(r8), parameter           :: om_hksat_lake  = 0.1_r8
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8), parameter           :: om_sucsat_lake = 10.3_r8
    ! saturated suction for organic matter (Letts, 2000)
    real(r8), parameter           :: om_b_lake      = 2.7_r8
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000) (lake)
    real(r8)           :: om_watsat
    ! porosity of organic soil
    real(r8)           :: om_hksat
    ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8)           :: om_sucsat
    ! saturated suction for organic matter (mm)(Letts, 2000)
    real(r8), parameter           :: om_csol        = 2.5_r8
    ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8), parameter           :: om_tkd         = 0.05_r8
    ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8)           :: om_b
    ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8), parameter           :: zsapric        = 0.5_r8
    ! depth (m) that organic matter takes on characteristics of sapric peat
    real(r8), parameter           :: pcalpha        = 0.5_r8       ! percolation threshold
    real(r8), parameter           :: pcbeta         = 0.139_r8     ! percolation exponent
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
         end if

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
    real(r8) :: dlon_span, dlat_span
    real(r8), parameter :: tol_degen = 1.0e-6_r8 
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
    if(allocated(longxy)) deallocate(longxy)
    allocate(longxy(ncells))
    if(allocated(latixy)) deallocate(latixy)
    allocate(latixy(ncells))

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

    ! Degenerate horizontal domain (e.g. 1x1 grid): min==max => zero span;
    ! Avoid dividing by (maxlon-minlon) or (maxlat-minlat).
    dlon_span = maxval(lon(:) + 180._r8) - minval(lon(:) + 180._r8)
    dlat_span = maxval(lat(:) + 90._r8) - minval(lat(:) + 90._r8)

    if(allocated(longxy_obs)) deallocate(longxy_obs)
    allocate(longxy_obs(dim_obs))
    if(allocated(latixy_obs)) deallocate(latixy_obs)
    allocate(latixy_obs(dim_obs))

    do i = 1, dim_obs
       if (dlon_span <= tol_degen .and. dlat_span <= tol_degen) then
          longxy_obs(i) = 1
          latixy_obs(i) = 1
       else if (dlon_span <= tol_degen) then
          longxy_obs(i) = 1
          latixy_obs(i) = ceiling(((lat_clmobs(i) + 90) - minlat) * nj / (maxlat - minlat))
       else if (dlat_span <= tol_degen) then
          longxy_obs(i) = ceiling(((lon_clmobs(i) + 180) - minlon) * ni / (maxlon - minlon))
          latixy_obs(i) = 1
       else if(((lon_clmobs(i) + 180) - minlon) /= 0 .and. &
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
       end if
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
    allocate(longxy_obs_floor(dim_obs))
    if(allocated(latixy_obs_floor)) deallocate(latixy_obs_floor)
    allocate(latixy_obs_floor(dim_obs))
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
       end if
    end do

  end subroutine get_interp_idx

#if defined CLMSA
  !> @author  Johannes Keller, adaptations by Yorck Ewerdwalbesloh
  !> @date    24.04.2025
  !> @brief   Set number of local analysis domains N_DOMAINS_P
  !> @details
  !>    This routine sets N_DOMAINS_P, the number of local analysis domains.
  subroutine init_n_domains_clm(n_domains_p)

    use decompMod, only : get_proc_bounds
    use clm_varcon      , only : ispval
    use ColumnType , only : col

    implicit none

    integer, intent(out) :: n_domains_p
    integer :: domain_p
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: begc, endc   ! per-proc beginning and ending column indices

    integer :: g
    integer :: c
    integer :: cc

    ! TODO: remove unnecessary calls of get_proc_bounds (use clm_begg,
    ! clm_endg, etc)
    call get_proc_bounds(begg=begg, endg=endg, begc=begc, endc=endc)

    if(clmupdate_swc==1) then
      if(clmstatevec_allcol==1) then
        ! Each column is a local domain
        ! -> DIM_L: number of layers in column
        n_domains_p = endc - begc + 1
      else
        ! Each gridcell is a local domain
        ! -> DIM_L: number of layers in gridcell
        n_domains_p = endg - begg + 1
      end if
    else
      ! Process-local number of gridcells Default, possibly not tested
      ! for other updates except SWC
      n_domains_p = endg - begg + 1
    end if

    NOGRACE: if (clmupdate_tws/=1) then

    ! If only_active: Use clm2pdaf to check which columsn/gridcells
    ! are inside. Possibly: number of columns/gridcells reduced by
    ! hydrologically inactive columns/gridcells.
    !
    ! Also: Set state_loc2clm_c_p: Returns the CLM-column c for the
    ! local domain domain_p (from the column, the gridcell can be
    ! derived)

    ! Allocate state_loc2clm_c_p with preliminary n_domains_p
    IF (allocated(state_loc2clm_c_p)) deallocate(state_loc2clm_c_p)
    allocate(state_loc2clm_c_p(n_domains_p))
    do domain_p=1,n_domains_p
      state_loc2clm_c_p(domain_p) = ispval
    end do

    if(clmstatevec_only_active == 1) then

      ! Reset n_domains_p
      n_domains_p = 0
      domain_p = 0

      if(clmstatevec_allcol == 1) then
        ! COLUMNS

        ! Each hydrologically active layer is a local domain
        ! -> DIM_L: number of layers in hydrologically active column
        do c=clm_begc,clm_endc
          ! Skip state vector loop and directly check if column is
          ! hydrologically active
          if(col%hydrologically_active(c)) then
            domain_p = domain_p + 1
            n_domains_p = n_domains_p + 1
            state_loc2clm_c_p(domain_p) = c
          end if
        end do

      else
        ! GRIDCELLS

        ! For gridcells
        do g = clm_begg,clm_endg

          ! Search the state vector for col in grc
          do cc = 1,clm_statevecsize

            if (col%gridcell(state_pdaf2clm_c_p(cc)) == g) then
              ! Set local domain index
              domain_p = domain_p + 1
              ! Set new number of local domains
              n_domains_p = n_domains_p + 1
              ! Set CLM-column-index corresponding to local domain
              state_loc2clm_c_p(domain_p) = state_pdaf2clm_c_p(cc)
              ! Exit state vector loop, when fitting column is found
              exit
            end if
          end do

        end do

      end if

    else

      ! Set state_loc2clm_c_p for non-excluding hydrologically
      ! inactive cols/grcs
      if(clmstatevec_allcol == 1) then
        ! COLUMNS
        do domain_p=1,n_domains_p
          state_loc2clm_c_p(domain_p) = clm_begc + domain_p - 1
        end do
      else
        ! GRIDCELLS
        do domain_p=1,n_domains_p
          state_loc2clm_c_p(domain_p) = clm_begg + domain_p - 1
        end do
      end if

    end if

    ! Possibly: Warning when final n_domains_p actually excludes
    ! hydrologically inactive gridcells

    else NOGRACE

      n_domains_p = num_hactiveg

    end if NOGRACE

  end subroutine init_n_domains_clm


  !> @author  Wolfgang Kurtz, Johannes Keller, Yorck Ewerdwalbesloh
  !> @date    20.11.2017
  !> @brief   Set local state vector dimension DIM_L local PDAF filters
  !> @details
  !>    This routine sets DIM_L, the local state vector dimension.
  subroutine init_dim_l_clm(domain_p, dim_l)
    use clm_varpar   , only : nlevsoi
    use ColumnType , only : col

    implicit none

    integer, intent(in)  :: domain_p
    integer, intent(out) :: dim_l
    integer              :: nshift

    integer :: g, i, count

    if(clmupdate_swc==1) then
      if(clmstatevec_only_active == 1) then
        ! Compare nlevsoi to clmstatevec_max_layer and bedrock if
        ! "hydrologically active" is turned on
        dim_l = min(nlevsoi, clmstatevec_max_layer, col%nbedrock(state_loc2clm_c_p(domain_p)))
        nshift = min(nlevsoi, clmstatevec_max_layer, col%nbedrock(state_loc2clm_c_p(domain_p)))
      else
        dim_l = nlevsoi
        nshift = nlevsoi
      end if
    end if

    if(clmupdate_swc==2) then
      error stop "Not implemented: clmupdate_swc.eq.2"
      ! dim_l = nlevsoi + 1
      ! nshift = nlevsoi + 1
    end if

    if(clmupdate_texture==1) then
      dim_l = 2*nlevsoi + nshift
    end if

    if(clmupdate_texture==2) then
      dim_l = 3*nlevsoi + nshift
    end if

    if (clmupdate_tws==1) then
      dim_l = 0
      g = hactiveg_levels(domain_p,1)

      select case(state_setup)
      case(0) ! subdivided setup
        do i = 1,nlevsoi
          do count = 1, num_layer(i)
            if (g==hactiveg_levels(count,i)) then
              ! I could also check with col%nbedrock but then I would need the column index and not the gridcell index
              dim_l = dim_l+1
            end if
          end do
        end do

        ! snow and surface water
        dim_l = dim_l+2

        if (clm_varsize_tws(5)/=0) then
          dim_l = dim_l+1
        end if

      case(1) ! only TWS in statevector

        dim_l=1

      case(2) ! aggregated setup

        dim_l=2

        do count = 1, num_layer(4)
          if (g==hactiveg_levels(count,4)) then
            dim_l = dim_l+1
          end if
        end do

        do count = 1, num_layer(13)
          if (g==hactiveg_levels(count,13)) then
            dim_l = dim_l+1
          end if
        end do
      case default

        error stop "Unsupported state_setup"

      end select
    end if

  end subroutine init_dim_l_clm

  !> @author  Wolfgang Kurtz, Johannes Keller, Yorck Ewerdwalbesloh
  !> @date    20.11.2017
  !> @brief   Set local state vector STATE_L from global state vector STATE_P
  !> @details
  !>    This routine sets STATE_L, the local state vector.
  !>
  !>    Source is STATE_P, the global (PE-local) state vector.
  subroutine g2l_state_clm(domain_p, dim_p, state_p, dim_l, state_l)

    use ColumnType , only : col

    implicit none

    INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
    INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
    INTEGER, INTENT(in) :: dim_l          ! Local state dimension
    REAL, TARGET, INTENT(in)    :: state_p(dim_p) ! PE-local full state vector
    REAL, TARGET, INTENT(out)   :: state_l(dim_l) ! State vector on local analysis d

    INTEGER :: i
    INTEGER :: n_domain
    INTEGER :: nshift_p

    integer :: sub, g, j

    ! call init_n_domains_clm(n_domain)

    ! DO i = 0, dim_l-1
    !   nshift_p = domain_p + i * n_domain
    !   state_l(i+1) = state_p(nshift_p)
    ! ENDDO
    NOGRACE: if (clmupdate_tws/=1) then
    ! Column index inside gridcell index domain_p
    DO i = 1, dim_l
      ! Column index from DOMAIN_P via STATE_LOC2CLM_C_P
      ! Layer index: i
      state_l(i) = state_p(state_clm2pdaf_p(state_loc2clm_c_p(domain_p),i))
    END DO
    else NOGRACE

      if (clm_varsize_tws(5)/=0) then
        sub=3
      else
        sub=2
      end if

      select case (state_setup)
      case(0) ! all compartmens in state vector
        g = hactiveg_levels(domain_p,1)
        do i = 1, dim_l-sub
          do j = 1, num_layer(i)
            if (g==hactiveg_levels(j,i)) then
              if (i == 1) then
                state_l(i) = state_p(j)
              else
                state_l(i) = state_p(j + sum(num_layer(1:i-1)))
              end if
            end if
          end do
        end do

        do j = 1, num_layer(1)
          if (g==hactiveg_levels(j,1)) then
            if (sub==3) then
              state_l(dim_l-2) = state_p(j + sum(clm_varsize_tws(1:2)))
              state_l(dim_l-1) = state_p(j + sum(clm_varsize_tws(1:3)))
              state_l(dim_l) = state_p(j + sum(clm_varsize_tws(1:4)))
            else
              state_l(dim_l-1) = state_p(j + sum(clm_varsize_tws(1:2)))
              state_l(dim_l) = state_p(j + sum(clm_varsize_tws(1:3)))
            end if
          end if
        end do

      case(1) ! only TWS in statevector

        g = hactiveg_levels(domain_p,1)
        do j = 1, num_layer(1)
          if (g==hactiveg_levels(j,1)) then
            state_l(1) = state_p(j)
          end if
        end do

      case(2) ! ! snow and soil moisture aggregated over surface, root zone and deep soil moisture in state vector

        g = hactiveg_levels(domain_p,1)

        do j = 1, num_layer(1)
          if (g==hactiveg_levels(j,1)) then
            state_l(1) = state_p(j) ! surface SM
            ! snow, same indexing as clm_varsize_tws(2:3) = 0 when only surface layers present
            state_l(dim_l) = state_p(j + sum(clm_varsize_tws(1:3)))
          end if
        end do

        if (dim_l>=3) then
          do j = 1, num_layer(4)
            if (g==hactiveg_levels(j,4)) then
              state_l(2) = state_p(j + clm_varsize_tws(1)) ! root zone SM
            end if
          end do
        end if

        if (dim_l>=4) then
          do j = 1, num_layer(13)
            if (g==hactiveg_levels(j,13)) then
              state_l(3) = state_p(j + sum(clm_varsize_tws(1:2))) ! deep SM
            end if
          end do
        end if

      case default

        error stop "Unsupported state_setup"

      end select

    end if NOGRACE

  end subroutine g2l_state_clm

  !> @author  Wolfgang Kurtz, Johannes Keller, Yorck Ewerdwalbesloh
  !> @date    20.11.2017
  !> @brief   Update global state vector STATE_P from local state vector STATE_L
  !> @details
  !>    This routine updates STATE_P, the global (PE-local) state vector.
  !>
  !>    Source is STATE_L, the local vector.
  subroutine l2g_state_clm(domain_p, dim_l, state_l, dim_p, state_p)

    use ColumnType , only : col

    implicit none

    INTEGER, INTENT(in) :: domain_p       ! Current local analysis domain
    INTEGER, INTENT(in) :: dim_l          ! Local state dimension
    INTEGER, INTENT(in) :: dim_p          ! PE-local full state dimension
    REAL, TARGET, INTENT(in)    :: state_l(dim_l) ! State vector on local analysis domain
    REAL, TARGET, INTENT(inout) :: state_p(dim_p) ! PE-local full state vector

    INTEGER :: i
    INTEGER :: n_domain
    INTEGER :: nshift_p

    integer :: sub, j, g

    ! ! beg and end gridcell for atm
    ! call init_n_domains_clm(n_domain)

    ! DO i = 0, dim_l-1
    !   nshift_p = domain_p + i * n_domain
    !   state_p(nshift_p) = state_l(i+1)
    ! ENDDO
    NOGRACE: if (clmupdate_tws==0) then
    ! Column index inside gridcell index domain_p
    DO i = 1, dim_l
      ! Column index from DOMAIN_P via STATE_LOC2CLM_C_P
      ! Layer index i
      state_p(state_clm2pdaf_p(state_loc2clm_c_p(domain_p),i)) = state_l(i)
    END DO
    else NOGRACE

      if (clm_varsize_tws(5)/=0) then
        sub=3
      else
        sub=2
      end if

      select case (state_setup)
      case(0) ! all compartments in state vector
        g = hactiveg_levels(domain_p,1)
        do i = 1, dim_l-sub
          do j = 1, num_layer(i) ! i is the layer that we are in right now
            ! if the counter is the gridcell of the local domain, we know the position in the statevector
            if (g==hactiveg_levels(j,i)) then
              if (i == 1) then ! if first layer
                state_p(j)  = state_l(i)    ! first liquid water as it is first in the statevector
              else
                state_p(j + sum(num_layer(1:i-1))) = state_l(i)
              end if
            end if
          end do
        end do

        do j = 1, num_layer(1)
          if (g==hactiveg_levels(j,1)) then

            if (sub==3) then
              state_p(j + sum(clm_varsize_tws(1:2))) = state_l(dim_l-2)
              state_p(j + sum(clm_varsize_tws(1:3))) = state_l(dim_l-1)
              state_p(j + sum(clm_varsize_tws(1:4))) = state_l(dim_l)
            else
              state_p(j + sum(clm_varsize_tws(1:2))) = state_l(dim_l-1)
              state_p(j + sum(clm_varsize_tws(1:3))) = state_l(dim_l)
            end if

          end if
        end do

      case(1) ! TWS in statevector

        g = hactiveg_levels(domain_p,1)
        do j = 1, num_layer(1)
          if (g==hactiveg_levels(j,1)) then
            state_p(j) = state_l(1)
          end if
        end do

      case(2) ! snow and soil moisture aggregated over surface, root zone and deep soil moisture in state vector

        g = hactiveg_levels(domain_p,1)

        do j = 1, num_layer(1)
          if (g==hactiveg_levels(j,1)) then

            state_p(j) = state_l(1)
            state_p(j + sum(clm_varsize_tws(1:3))) = state_l(dim_l)

          end if
        end do

        do j = 1, num_layer(4)
          if (g==hactiveg_levels(j,4)) then

            state_p(j + clm_varsize_tws(1)) = state_l(2)

          end if
        end do

        do j = 1, num_layer(13)
          if (g==hactiveg_levels(j,13)) then

            state_p(j + sum(clm_varsize_tws(1:2))) = state_l(3)

          end if
        end do

      case default

        error stop "Unsupported state_setup"

      end select

    end if NOGRACE

  end subroutine l2g_state_clm
#endif

end module enkf_clm_mod


!> PDAF-OMI observation module for type LST observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! OBSTYPE = LST
!!
!! __Observation type LST:__
!! The observation type LST are Land Surface Temperature observations.
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF
!! usually by callback_obs_pdafomi.F90
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data types (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!! Three observation operator modes are supported via clmupdate_T:
!! * clmupdate_T==1: Kustas (2009) formula combining ground (TG) and
!!                  vegetation (TV) temperature:
!!                  LST = (exp(-tau/2)*TG^4 + (1-exp(-tau/2))*TV^4)^0.25
!!                  where tau = clm_paramarr(obs_index) (vegetation optical depth)
!! * clmupdate_T==2,3: Direct TSKIN mapping (simulated LST = TSKIN)
!!
!! The state vector for LST is patch-based (not column-based as for SM).
!! For clmupdate_T==1 the state vector contains TG (first) and TV (second
!! half, offset by clm_varsize).
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_LST \n
!!           Count number of process-local and full observations;
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_LST \n
!!           observation operator to get full observation vector of this type.
!!
!! __Revision history:__
!! * 2024 - Initial code based on obs_SM_pdafomi.F90 (Yorck Ewerdwalbesloh)
!!          adapted for LST-DA (clmupdate_T)
!!

! Author: Based on obs_SM_pdafomi.F90 by Yorck Ewerdwalbesloh;
!         adapted for LST-DA (clmupdate_T)

#ifdef CLMFIVE
MODULE obs_LST_pdafomi

    USE mod_parallel_pdaf, &
         ONLY: mype_filter    ! Rank of filter process
    USE PDAFomi, &
         ONLY: obs_f, obs_l   ! Declaration of observation data types

    IMPLICIT NONE
    SAVE

    PUBLIC

    ! Variables which are inputs to the module (usually set in init_pdaf)
    LOGICAL :: assim_LST       !< Whether to assimilate this data type
    REAL    :: rms_obs_LST     !< Observation error standard deviation (for constant errors)

    ! longitude and latitude of grid cells and observation cells
    INTEGER, ALLOCATABLE :: longxy(:), latixy(:), longxy_obs(:), latixy_obs(:)

    ! Module-local obs-to-state index array (avoids conflict with global obs_index_p
    ! from mod_assimilation when SM and LST are both active simultaneously)
    INTEGER, ALLOCATABLE :: obs_index_p_LST(:)

  ! *********************************************************
  ! *** Data type obs_f defines the full observations by  ***
  ! *** internally shared variables of the module         ***
  ! *********************************************************

  ! Declare instances of observation data types used here
    TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
    TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

  !$OMP THREADPRIVATE(thisobs_l)


  !-------------------------------------------------------------------------------

  CONTAINS

  !> Initialize information on the module-type observation
  !!
  !! The routine is called by each filter process.
  !! at the beginning of the analysis step before
  !! the loop through all local analysis domains.
  !!
  !! It has to count the number of observations of the
  !! observation type handled in this module according
  !! to the current time step for all observations
  !! required for the analyses in the loop over all local
  !! analysis domains on the PE-local state domain.
  !!
SUBROUTINE init_dim_obs_LST(step, dim_obs)
  
  USE mpi, ONLY: MPI_INTEGER
  USE mpi, ONLY: MPI_DOUBLE_PRECISION
  USE mpi, ONLY: MPI_SUM
  USE mpi, ONLY: MPI_IN_PLACE
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, comm_filter, npes_filter, abort_parallel, &
       mype_world
  USE mod_assimilation, &
       ONLY: obs_filename, &
       obs_pdaf2nc, obs_nc2pdaf, &
       local_dims_obs, &
       local_disp_obs, &
       longxy_obs_floor, latixy_obs_floor, &
       screen, cradius_LST

  USE PDAFomi, &
    ONLY: PDAFomi_gather_obs, pi

  Use mod_read_obs, &
    only: multierr, read_obs_nc_type

  use enkf_clm_mod, only: state_clm2pdaf_p

  use enkf_clm_mod, only: domain_def_clm

  use mod_parallel_pdaf, &
    only: abort_parallel

  use shr_kind_mod, only: r8 => shr_kind_r8

  use GridcellType, only: grc

  use clm_varcon, only: ispval

  use decompMod , only : get_proc_bounds, get_proc_global

  use PatchType , only : patch

  use enkf_clm_mod, only: get_interp_idx

  use mod_tsmp, only: da_print_obs_index

  IMPLICIT NONE

  ! *** Arguments ***
  INTEGER, INTENT(in)    :: step       !< Current time step
  INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

  ! *** Local variables ***
  INTEGER :: i, p, g, pg            ! Counters
  INTEGER :: cnt                    ! Counter
  INTEGER :: dim_obs_p              ! Number of process-local observations
  REAL, ALLOCATABLE :: obs_p(:)     ! PE-local observation vector
  REAL, ALLOCATABLE :: obs_g(:)     ! Global observation vector
  REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
  REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates
  character (len = 110) :: current_observation_filename


  character(len=20) :: obs_type_name ! name of observation type (e.g. GRACE, SM, ST, ...)

  REAL, ALLOCATABLE :: lon_obs(:)
  REAL, ALLOCATABLE :: lat_obs(:)
  INTEGER, ALLOCATABLE :: layer_obs(:)
  REAL, ALLOCATABLE :: dr_obs(:)
  REAL, ALLOCATABLE :: obserr(:)
  REAL, ALLOCATABLE :: obscov(:,:)

  integer :: begp, endp   ! per-proc beginning and ending pft indices
  integer :: begc, endc   ! per-proc beginning and ending column indices
  integer :: begl, endl   ! per-proc beginning and ending landunit indices
  integer :: begg, endg   ! per-proc gridcell ending gridcell indices

  integer :: numg         ! total number of gridcells across all processors
  integer :: numl         ! total number of landunits across all processors
  integer :: numc         ! total number of columns across all processors
  integer :: nump         ! total number of pfts across all processors

  real(r8), pointer :: lon(:)
  real(r8), pointer :: lat(:)

  integer :: ierror

  real :: deltax, deltay

  logical :: is_use_dr
  logical :: obs_snapped     !Switch for checking multiple observation counts
  logical :: newgridcell

  INTEGER :: sum_dim_obs_p

  character (len = 27) :: fn    !TSMP-PDAF: function name for obs_index_p output



  ! *********************************************
  ! *** Initialize full observation dimension ***
  ! *********************************************

  IF (mype_filter==0) &
    WRITE (*,*) 'Assimilate observations - obs type LST'

  IF (assim_LST) thisobs%doassim = 1

  ! Geographic distance with haversine formula
  thisobs%disttype = 3
  thisobs%ncoord = 2


  ! **********************************
  ! *** Read PE-local observations ***
  ! **********************************


  obs_type_name = 'LST'

  ! now call function to get observations

  if(mype_filter==0 .and. screen > 2) then
    write(*,*)'load observations from type LST'
  end if
  write(current_observation_filename, '(a, i5.5)') trim(obs_filename)//'.', step


  if (mype_filter == 0) then
    call read_obs_nc_type(current_observation_filename, obs_type_name, &
      dim_obs, obs_g, lon_obs, lat_obs, layer_obs, &
      dr_obs, obserr, obscov)
  end if

  call mpi_bcast(dim_obs, 1, MPI_INTEGER, 0, comm_filter, ierror)

  ! check if file contains observations of type LST

  if (dim_obs == 0) then
    if (mype_filter==0 .and. screen > 2) then
      write(*,*)'TSMP-PDAF mype(w) =', mype_world, &
        ': No observations of type LST found in file ', &
        trim(current_observation_filename)
    end if
    dim_obs_p = 0
    if (allocated(obs_p)) deallocate(obs_p)
    if (allocated(ivar_obs_p)) deallocate(ivar_obs_p)
    if (allocated(ocoord_p)) deallocate(ocoord_p)
    if (allocated(thisobs%id_obs_p)) deallocate(thisobs%id_obs_p)
    ALLOCATE(obs_p(1))
    ALLOCATE(ivar_obs_p(1))
    ALLOCATE(ocoord_p(2, 1))
    ALLOCATE(thisobs%id_obs_p(1, 1))
    thisobs%id_obs_p(1, 1) = 0
    thisobs%infile = 0
    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
      thisobs%ncoord, cradius_LST, dim_obs)
    if(mype_filter==0) DEALLOCATE(obs_g)
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)
    return
  end if

  call mpi_bcast(multierr, 1, MPI_INTEGER, 0, comm_filter, ierror)

  ! Allocate observation arrays for non-root procs
  ! ----------------------------------------------
  if (mype_filter /= 0) then ! for all non-master proc
    if(allocated(obs_g)) deallocate(obs_g)
    allocate(obs_g(dim_obs))
    if(allocated(lon_obs)) deallocate(lon_obs)
    allocate(lon_obs(dim_obs))
    if(allocated(lat_obs)) deallocate(lat_obs)
    allocate(lat_obs(dim_obs))
    if(allocated(dr_obs)) deallocate(dr_obs)
    allocate(dr_obs(dim_obs))
    if(allocated(layer_obs)) deallocate(layer_obs)
    allocate(layer_obs(dim_obs))
    if(multierr==1) then
      if(allocated(obserr)) deallocate(obserr)
      allocate(obserr(dim_obs))
    end if
  end if

  call mpi_bcast(obs_g, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  if(multierr==1) call mpi_bcast(obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  call mpi_bcast(lon_obs, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  call mpi_bcast(lat_obs, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  call mpi_bcast(dr_obs,  dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
  call mpi_bcast(layer_obs, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)



  if (mype_filter==0 .and. screen > 2) then
    write(*,*)'Done: load observations from type LST'
  end if


  thisobs%infile = 1
  call domain_def_clm(lon_obs, lat_obs, dim_obs, longxy, latixy, longxy_obs, latixy_obs)

  ! Obtain CLM lon/lat information
  lon => grc%londeg
  lat => grc%latdeg

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)

  dim_obs_p = 0
  is_use_dr = .true.

  ! id_obs_p: placeholder to satisfy PDAFomi internal check;
  ! actual obs-to-state mapping is in obs_index_p_LST
  if(allocated(thisobs%id_obs_p)) deallocate(thisobs%id_obs_p)
  allocate(thisobs%id_obs_p(1, 1))
  thisobs%id_obs_p(1, 1) = 1

  ! *** Count PE-local observations ***
  ! LST uses patch loop: one observation per gridcell, assigned to
  ! the first patch found in that gridcell.
  do i = 1, dim_obs
    obs_snapped = .false.
    do g = begg,endg
      newgridcell = .true.
      do p = begp,endp
        pg = patch%gridcell(p)
        if(pg == g) then
          if(newgridcell) then

            if(is_use_dr) then
              if(lon(g)>180) then
                deltax = abs(lon(g)-lon_obs(i)-360)
              else
                deltax = abs(lon(g)-lon_obs(i))
              end if
              deltay = abs(lat(g)-lat_obs(i))
            end if

            if(((is_use_dr).and.(deltax<=dr_obs(1)).and.(deltay<=dr_obs(1))).or. &
              ((.not. is_use_dr).and.(longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1)))) then
              dim_obs_p = dim_obs_p + 1

              if(obs_snapped) then
                print *, "TSMP-PDAF mype(w)=", mype_world, &
                  ": ERROR Observation snapped at multiple grid cells."
                print *, "i=", i
                call abort_parallel()
              end if
              obs_snapped = .true.
            end if

            newgridcell = .false.

          end if
        end if
      end do
    end do
  end do

  if(screen > 2) then
    print *, "TSMP-PDAF mype(w)=", mype_world, &
      ": init_dim_obs_LST: dim_obs_p=", dim_obs_p
  end if

  ! Initialize OMI arrays
  IF (ALLOCATED(ivar_obs_p)) DEALLOCATE(ivar_obs_p)
  ALLOCATE(ivar_obs_p(dim_obs_p))
  IF (ALLOCATED(ocoord_p)) DEALLOCATE(ocoord_p)
  ALLOCATE(ocoord_p(2, dim_obs_p))

  ! Gather and check PE-local observation dimensions
  call mpi_allreduce(dim_obs_p, sum_dim_obs_p, 1, MPI_INTEGER, MPI_SUM, &
    comm_filter, ierror)

  if(.not. sum_dim_obs_p == dim_obs) then
    print *, "TSMP-PDAF mype(w)=", mype_world, &
      ": ERROR Sum of PE-local observation dimensions"
    print *, "sum_dim_obs_p=", sum_dim_obs_p
    print *, "dim_obs=", dim_obs
    call abort_parallel()
  end if

  IF (ALLOCATED(local_dims_obs)) DEALLOCATE(local_dims_obs)
  ALLOCATE(local_dims_obs(npes_filter))
  call mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dims_obs, 1, MPI_INTEGER, &
    comm_filter, ierror)

  IF (ALLOCATED(local_disp_obs)) DEALLOCATE(local_disp_obs)
  ALLOCATE(local_disp_obs(npes_filter))
  local_disp_obs(1) = 0
  do i = 2, npes_filter
    local_disp_obs(i) = local_disp_obs(i-1) + local_dims_obs(i-1)
  end do

  if(mype_filter==0 .and. screen > 2) then
    print *, "TSMP-PDAF mype(w)=", mype_world, &
      ": init_dim_obs_LST: local_disp_obs=", local_disp_obs
  end if

  ! Write index mapping obs_pdaf2nc / obs_nc2pdaf (used for debug output)
  if(allocated(obs_pdaf2nc)) deallocate(obs_pdaf2nc)
  allocate(obs_pdaf2nc(dim_obs))
  obs_pdaf2nc = 0
  if(allocated(obs_nc2pdaf)) deallocate(obs_nc2pdaf)
  allocate(obs_nc2pdaf(dim_obs))
  obs_nc2pdaf = 0

  cnt = 1
  do i = 1, dim_obs
    obs_snapped = .true.
    do g = begg,endg
      newgridcell = .true.
      do p = begp,endp
        pg = patch%gridcell(p)
        if(pg == g) then
          if(newgridcell) then

            if(is_use_dr) then
              if(lon(g)>180) then
                deltax = abs(lon(g)-lon_obs(i)-360)
              else
                deltax = abs(lon(g)-lon_obs(i))
              end if
              deltay = abs(lat(g)-lat_obs(i))
            end if

            if(((is_use_dr).and.(deltax<=dr_obs(1)).and.(deltay<=dr_obs(1))).or. &
              ((.not. is_use_dr).and.(longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1)))) then
              if(state_clm2pdaf_p(p,1)==ispval) then
                obs_snapped = .false.
                cycle
              end if
              obs_pdaf2nc(local_disp_obs(mype_filter+1)+cnt) = i
              obs_nc2pdaf(i) = local_disp_obs(mype_filter+1)+cnt
              cnt = cnt + 1
              obs_snapped = .true.
            end if

            newgridcell = .false.

          end if
        end if
      end do
    end do

    if(.not. obs_snapped) then
      print *, "TSMP-PDAF mype(w)=", mype_world, &
        ": ERROR observations exist at non-active gridcells."
      print *, "Observation-index in NetCDF-file: i=", i
      call abort_parallel()
    end if
  end do

  call mpi_allreduce(MPI_IN_PLACE,obs_pdaf2nc,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)
  call mpi_allreduce(MPI_IN_PLACE,obs_nc2pdaf,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)

  if(mype_filter==0 .and. screen > 2) then
    print *, "TSMP-PDAF mype(w)=", mype_world, &
      ": init_dim_obs_LST: obs_pdaf2nc=", obs_pdaf2nc
  end if

  ! Write process-local observation arrays
  ! Use module-local obs_index_p_LST to avoid conflict with global obs_index_p
  IF (ALLOCATED(obs_index_p_LST)) DEALLOCATE(obs_index_p_LST)
  ALLOCATE(obs_index_p_LST(dim_obs_p))
  IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
  ALLOCATE(obs_p(dim_obs_p))

  cnt = 1

  do i = 1, dim_obs

    do g = begg,endg
      newgridcell = .true.

      do p = begp,endp

        pg = patch%gridcell(p)

        if(pg == g) then

          if(newgridcell) then

            if(is_use_dr) then
              if(lon(g)>180) then
                deltax = abs(lon(g)-lon_obs(i)-360)
              else
                deltax = abs(lon(g)-lon_obs(i))
              end if
              deltay = abs(lat(g)-lat_obs(i))
            end if

            if(((is_use_dr).and.(deltax<=dr_obs(1)).and.(deltay<=dr_obs(1))).or. &
              ((.not. is_use_dr).and.(longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1)))) then

              ! Convert observation coordinates to radians for haversine distance
              if(thisobs%disttype==3) then
                ocoord_p(1,cnt) = lon_obs(i) * pi / 180.0
                ocoord_p(2,cnt) = lat_obs(i) * pi / 180.0
              else
                ocoord_p(1,cnt) = lon_obs(i)
                ocoord_p(2,cnt) = lat_obs(i)
              end if

              ! Set state vector index for this patch.
              ! LST uses patch (not column), layer 1 = surface temperature.
              obs_index_p_LST(cnt) = state_clm2pdaf_p(p,1)

              obs_p(cnt) = obs_g(i)
              if(multierr==1) ivar_obs_p(cnt) = 1.0/(obserr(i)*obserr(i))
              if(multierr==0) ivar_obs_p(cnt) = 1.0/(rms_obs_LST*rms_obs_LST)
              cnt = cnt + 1

            end if

            newgridcell = .false.

          end if

        end if

      end do
    end do

  end do

#ifdef PDAF_DEBUG
  IF (da_print_obs_index > 0) THEN
    WRITE(fn, "(a,i5.5,a,i5.5,a)") "obs_index_p_LST_", mype_world, ".", step, ".txt"
    OPEN(unit=72, file=fn, action="write")
    DO i = 1, dim_obs_p
      WRITE (72,"(i10)") obs_index_p_LST(i)
    END DO
    CLOSE(72)
  END IF
#endif

  ! ****************************************
  ! *** Gather global observation arrays ***
  ! ****************************************

  CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
    thisobs%ncoord, cradius_LST, dim_obs)

  ! ********************
  ! *** Finishing up ***
  ! ********************

  DEALLOCATE(obs_g)
  DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

END SUBROUTINE init_dim_obs_LST



!-------------------------------------------------------------------------------
!> Implementation of observation operator for LST
!!
!! Applies the observation operator H for LST. Three modes:
!!   clmupdate_T==1: Kustas (2009) Eq.(7)
!!                   LST = (exp(-tau/2)*TG^4 + (1-exp(-tau/2))*TV^4)^0.25
!!   clmupdate_T==2,3: direct TSKIN mapping, LST = state_p(obs_index)
!!
SUBROUTINE obs_op_LST(dim_p, dim_obs, state_p, ostate)

  use enkf_clm_mod, only: clm_paramarr, clm_varsize, clmupdate_T

  use PDAFomi_obs_f, only: PDAFomi_gather_obsstate

  IMPLICIT NONE

  INTEGER, INTENT(in) :: dim_p
  INTEGER, INTENT(in) :: dim_obs
  REAL, INTENT(in)    :: state_p(dim_p)
  REAL, INTENT(inout) :: ostate(dim_obs)

  real, allocatable :: ostate_p(:)
  integer :: i

  IF (thisobs%dim_obs_p > 0) THEN
    if(allocated(ostate_p)) deallocate(ostate_p)
    ALLOCATE(ostate_p(thisobs%dim_obs_p))
  ELSE
    if(allocated(ostate_p)) deallocate(ostate_p)
    ALLOCATE(ostate_p(1))
  END IF

  if(clmupdate_T == 1) then

    ! Kustas (2009) Eq.(7): combined ground/vegetation brightness temperature
    ! tau    = clm_paramarr(idx): vegetation optical depth (LAI-based)
    ! TG     = state_p(obs_index_p_LST(i)):              ground temperature
    ! TV     = state_p(clm_varsize + obs_index_p_LST(i)): vegetation temperature
    DO i = 1, thisobs%dim_obs_p
      ostate_p(i) = &
        ( exp(-0.5*clm_paramarr(obs_index_p_LST(i))) &
        * state_p(obs_index_p_LST(i))**4 &
        + (1.0 - exp(-0.5*clm_paramarr(obs_index_p_LST(i)))) &
        * state_p(clm_varsize + obs_index_p_LST(i))**4 )**0.25
    END DO

  else

    ! clmupdate_T==2 or 3: simulated LST equals TSKIN directly
    DO i = 1, thisobs%dim_obs_p
      ostate_p(i) = state_p(obs_index_p_LST(i))
    END DO

  end if

  CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)

  deallocate(ostate_p)

END SUBROUTINE obs_op_LST



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! Called during the loop over all local analysis domains.
!! Uses patch-based localization domain coordinates.
!!
SUBROUTINE init_dim_obs_l_LST(domain_p, step, dim_obs, dim_obs_l)

  USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l, pi

  USE mod_assimilation, &
    ONLY: cradius_LST, locweight, sradius_LST, screen

  USE enkf_clm_mod, ONLY: state_loc2clm_p_p

  use shr_kind_mod, only: r8 => shr_kind_r8

  use decompMod , only : get_proc_bounds

  USE GridcellType, ONLY: grc
  USE PatchType, ONLY: patch

  use clm_varcon, only: spval

  use mod_parallel_pdaf, ONLY: mype_world

  IMPLICIT NONE

  INTEGER, INTENT(in)    :: domain_p
  INTEGER, INTENT(in)    :: step
  INTEGER, INTENT(in)    :: dim_obs
  INTEGER, INTENT(inout) :: dim_obs_l

  REAL :: coords_l(2)

  real(r8), pointer :: lon(:)
  real(r8), pointer :: lat(:)

  integer :: pg_l

  lon => grc%londeg
  lat => grc%latdeg

  ! **********************************************
  ! *** Initialize local observation dimension ***
  ! **********************************************

  if(thisobs%infile==1) then

    ! Get gridcell of local analysis domain via patch index
    pg_l = patch%gridcell(state_loc2clm_p_p(domain_p))

    if(lon(pg_l) > 180) then
      coords_l(1) = lon(pg_l) - 360.0
    else
      coords_l(1) = lon(pg_l)
    end if
    coords_l(2) = lat(pg_l)

    if(thisobs%disttype==3) then
      coords_l(1) = coords_l(1) * pi / 180.0
      coords_l(2) = coords_l(2) * pi / 180.0
    end if

  else

    coords_l(1) = spval
    coords_l(2) = spval

  end if

  ! For disttype=3, cradius and sradius are in km; multiply by 1000 for meters
  if(thisobs%disttype==3) then
    CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
      locweight, cradius_LST*1000.0, sradius_LST*1000.0, dim_obs_l)
  else
    CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
      locweight, cradius_LST, sradius_LST, dim_obs_l)
  end if

END SUBROUTINE init_dim_obs_l_LST



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF on the module-type observation
!!
!! Uses patch-based state vector coordinates.
!!
SUBROUTINE localize_covar_LST(dim_p, dim_obs, HP_p, HPH, coords_p)

  USE PDAFomi, ONLY: PDAFomi_localize_covar

  USE mod_assimilation, &
    ONLY: cradius_LST, locweight, sradius_LST

  use enkf_clm_mod, only: state_pdaf2clm_p_p

  use shr_kind_mod, only: r8 => shr_kind_r8

  USE GridcellType, ONLY: grc
  USE PatchType, ONLY: patch

  IMPLICIT NONE

  INTEGER, INTENT(in) :: dim_p
  INTEGER, INTENT(in) :: dim_obs
  REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)
  REAL, INTENT(inout) :: HPH(dim_obs, dim_obs)
  REAL, INTENT(inout) :: coords_p(:,:)

  integer :: i, pg_i

  real(r8), pointer :: lon(:)
  real(r8), pointer :: lat(:)

  lon => grc%londeg
  lat => grc%latdeg

  ! Assign geographic coordinates from patch gridcells to state vector elements
  do i = 1,dim_p
    pg_i = patch%gridcell(state_pdaf2clm_p_p(i))
    if(lon(pg_i) > 180) then
      coords_p(1,i) = lon(pg_i) - 360.0
    else
      coords_p(1,i) = lon(pg_i)
    end if
    coords_p(2,i) = lat(pg_i)
  end do

  CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, cradius_LST, sradius_LST, &
    coords_p, HP_p, HPH)

END SUBROUTINE localize_covar_LST


subroutine add_obs_err_LST(step, dim_obs, C)

  USE mod_parallel_pdaf, ONLY: npes_filter

  use PDAFomi, only: obsdims

  implicit none
  INTEGER, INTENT(in) :: step
  INTEGER, INTENT(in) :: dim_obs
  REAL, INTENT(inout) :: C(dim_obs,dim_obs)

  integer :: i, pe, cnt
  INTEGER, ALLOCATABLE :: id_start(:)
  INTEGER, ALLOCATABLE :: id_end(:)

  ALLOCATE(id_start(npes_filter), id_end(npes_filter))

  pe = 1
  id_start(1) = 1
  IF (thisobs%obsid>1) id_start(1) = id_start(1) + sum(obsdims(1, 1:thisobs%obsid-1))
  id_end(1) = id_start(1) + obsdims(1,thisobs%obsid) - 1
  DO pe = 2, npes_filter
    id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1,thisobs%obsid:))
    IF (thisobs%obsid>1) id_start(pe) = id_start(pe) + sum(obsdims(pe,1:thisobs%obsid-1))
    id_end(pe) = id_start(pe) + obsdims(pe,thisobs%obsid) - 1
  END DO

  cnt = 1
  DO pe = 1, npes_filter
    DO i = id_start(pe), id_end(pe)
      C(i,i) = C(i,i) + 1.0/thisobs%ivar_obs_f(cnt)
      cnt = cnt + 1
    end do
  end do

  DEALLOCATE(id_start, id_end)

end subroutine add_obs_err_LST


subroutine init_obscovar_LST(step, dim_obs, dim_obs_p, covar, m_state_p, isdiag)

  USE mod_parallel_pdaf, ONLY: npes_filter

  use PDAFomi, only: obsdims, map_obs_id

  implicit none
  INTEGER, INTENT(in) :: step
  INTEGER, INTENT(in) :: dim_obs
  INTEGER, INTENT(in) :: dim_obs_p
  REAL, INTENT(inout) :: covar(dim_obs, dim_obs)
  REAL, INTENT(in)    :: m_state_p(dim_obs_p)
  LOGICAL, INTENT(inout) :: isdiag

  integer :: i, pe, cnt
  INTEGER, ALLOCATABLE :: id_start(:)
  INTEGER, ALLOCATABLE :: id_end(:)

  ALLOCATE(id_start(npes_filter), id_end(npes_filter))

  pe = 1
  id_start(1) = 1
  IF (thisobs%obsid>1) id_start(1) = id_start(1) + sum(obsdims(1, 1:thisobs%obsid-1))
  id_end(1) = id_start(1) + obsdims(1,thisobs%obsid) - 1
  DO pe = 2, npes_filter
    id_start(pe) = id_start(pe-1) + SUM(obsdims(pe-1,thisobs%obsid:))
    IF (thisobs%obsid>1) id_start(pe) = id_start(pe) + sum(obsdims(pe,1:thisobs%obsid-1))
    id_end(pe) = id_start(pe) + obsdims(pe,thisobs%obsid) - 1
  END DO

  cnt = 1
  IF (thisobs%obsid-1 > 0) cnt = cnt + SUM(obsdims(:,1:thisobs%obsid-1))
  DO pe = 1, npes_filter
    DO i = id_start(pe), id_end(pe)
      map_obs_id(i) = cnt
      cnt = cnt + 1
    END DO
  END DO

  cnt = 1
  DO pe = 1, npes_filter
    DO i = id_start(pe), id_end(pe)
      covar(i, i) = covar(i, i) + 1.0/thisobs%ivar_obs_f(cnt)
      cnt = cnt + 1
    ENDDO
  ENDDO

  isdiag = .TRUE.

  DEALLOCATE(id_start, id_end)

end subroutine init_obscovar_LST


subroutine prodRinvA_LST(step, dim_obs_p, rank, obs_p, A_p, C_p)

  INTEGER, INTENT(in) :: step
  INTEGER, INTENT(in) :: dim_obs_p
  INTEGER, INTENT(in) :: rank
  REAL, INTENT(in)    :: obs_p(dim_obs_p)
  REAL, INTENT(in)    :: A_p(dim_obs_p,rank)
  REAL, INTENT(inout) :: C_p(dim_obs_p,rank)

  INTEGER :: i, j, off

  off = thisobs%off_obs_f

  do j = 1, rank
    do i = 1, thisobs%dim_obs_f
      C_p(i+off, j) = thisobs%ivar_obs_f(i) * A_p(i+off, j)
    END DO
  end do

end subroutine prodRinvA_LST


subroutine prodRinvA_l_LST(domain_p, step, dim_obs_l, rank, obs_l, A_l, C_l)

  use shr_kind_mod, only: r8 => shr_kind_r8
  USE mod_assimilation, ONLY: cradius_LST, locweight, sradius_LST
  use pdafomi, only: PDAFomi_observation_localization_weights

  implicit none

  INTEGER, INTENT(in) :: domain_p
  INTEGER, INTENT(in) :: step
  INTEGER, INTENT(in) :: dim_obs_l
  INTEGER, INTENT(in) :: rank
  REAL, INTENT(in)    :: obs_l(dim_obs_l)
  REAL, INTENT(inout) :: A_l(dim_obs_l, rank)
  REAL, INTENT(out)   :: C_l(dim_obs_l, rank)

  INTEGER :: verbose
  INTEGER, SAVE :: domain_save = -1
  REAL, ALLOCATABLE :: weight(:)
  INTEGER :: i, j, off, idummy

  off = thisobs_l%off_obs_l
  idummy = dim_obs_l

  IF ((domain_p <= domain_save .OR. domain_save < 0) .AND. mype_filter==0) THEN
    verbose = 1
  ELSE
    verbose = 0
  END IF
  domain_save = domain_p

  IF (verbose == 1) THEN
    WRITE (*, '(8x, a, f12.3)') &
      '--- Use global rms for LST observations of ', rms_obs_LST
    WRITE (*, '(8x, a, 1x)') &
      '--- Domain localization'
    WRITE (*, '(12x, a, 1x, f12.2)') &
      '--- Local influence radius', cradius_LST
    IF (locweight > 0) THEN
      WRITE (*, '(12x, a)') &
        '--- Use distance-dependent weight for observation errors'
    END IF
  ENDIF

  ALLOCATE(weight(thisobs_l%dim_obs_l))
  call PDAFomi_observation_localization_weights(thisobs_l, thisobs, rank, A_l, &
    weight, verbose)

  do j = 1, rank
    do i = 1, thisobs_l%dim_obs_l
      C_l(i+off,j) = thisobs_l%ivar_obs_l(i) * weight(i) * A_l(i+off, j)
    end do
  end do

  deallocate(weight)

end subroutine prodRinvA_l_LST


subroutine deallocate_obs_LST()

  USE PDAFomi, ONLY: PDAFomi_deallocate_obs
  USE PDAFomi_obs_l, ONLY: obs_l_all, firstobs

  implicit none

  if(mype_filter==0) then
    WRITE (*,*) 'Deallocating observations type LST'
  end if
  call PDAFomi_deallocate_obs(thisobs)

  if(allocated(thisobs_l%id_obs_l))   deallocate(thisobs_l%id_obs_l)
  if(allocated(thisobs_l%ivar_obs_l)) deallocate(thisobs_l%ivar_obs_l)
  if(allocated(thisobs_l%distance_l)) deallocate(thisobs_l%distance_l)
  if(allocated(thisobs_l%cradius_l))  deallocate(thisobs_l%cradius_l)
  if(allocated(thisobs_l%sradius_l))  deallocate(thisobs_l%sradius_l)
  if(allocated(thisobs_l%dist_l_v))   deallocate(thisobs_l%dist_l_v)
  if(allocated(thisobs_l%cradius))    deallocate(thisobs_l%cradius)
  if(allocated(thisobs_l%sradius))    deallocate(thisobs_l%sradius)

  if(allocated(obs_l_all)) deallocate(obs_l_all)

  firstobs = 0

  if(allocated(obs_index_p_LST)) deallocate(obs_index_p_LST)

end subroutine deallocate_obs_LST


END MODULE obs_LST_pdafomi
#endif

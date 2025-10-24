!> PDAF-OMI observation module for type SM observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! OBSTYPE = SM
!!
!! __Observation type SM:__
!! The observation type SM are soil moisture observations
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
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_OBSTYPE \n
!!           Count number of process-local and full observations;
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_OBSTYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routines, which are required if filters
!! with localization are used:
!! * init_dim_obs_l_OBSTYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_OBSTYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_SM_pdafomi

    USE mod_parallel_pdaf, &
         ONLY: mype_filter    ! Rank of filter process
    USE PDAFomi, &
         ONLY: obs_f, obs_l   ! Declaration of observation data types

    IMPLICIT NONE
    SAVE

    PUBLIC

    ! Variables which are inputs to the module (usually set in init_pdaf)
    LOGICAL :: assim_SM        !< Whether to assimilate this data type
    REAL    :: rms_obs_SM      !< Observation error standard deviation (for constant errors)

    INTEGER, ALLOCATABLE :: longxy(:), latixy(:), longxy_obs(:), latixy_obs(:) ! longitude and latitude of grid cells and observation cells

    ! One can declare further variables, e.g. for file names which can
    ! be use-included in init_pdaf() and initialized there.


  ! *********************************************************
  ! *** Data type obs_f defines the full observations by  ***
  ! *** internally shared variables of the module         ***
  ! *********************************************************

  ! Relevant variables that can be modified by the user:
  !   TYPE obs_f
  !      ---- Mandatory variables to be set in INIT_DIM_OBS ----
  !      INTEGER :: doassim                    ! Whether to assimilate this observation type
  !      INTEGER :: disttype                   ! Type of distance computation to use for localization
  !                                            ! (0) Cartesian, (1) Cartesian periodic
  !                                            ! (2) simplified geographic, (3) geographic haversine function
  !      INTEGER :: ncoord                     ! Number of coordinates use for distance computation
  !      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
  !
  !      ---- Optional variables - they can be set in INIT_DIM_OBS ----
  !      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
  !      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
  !
  !      ---- Variables with predefined values - they can be changed in INIT_DIM_OBS  ----
  !      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
  !      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs.
  !                                           ! or (0) obs. restricted to those relevant for a process domain
  !      REAL :: inno_omit=0.0                ! Omit obs. if squared innovation larger this factor times
  !                                           !     observation variance
  !      REAL :: inno_omit_ivar=1.0e-12       ! Value of inverse variance to omit observation
  !   END TYPE obs_f

  ! Data type obs_l defines the local observations by internally shared variables of the module

  ! ***********************************************************************

  ! Declare instances of observation data types used here
  ! We use generic names here, but one could rename the variables
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
  !! The following four variables have to be initialized in this routine
  !! * thisobs\%doassim     - Whether to assimilate this type of observations
  !! * thisobs\%disttype    - type of distance computation for localization with this observaton
  !! * thisobs\%ncoord      - number of coordinates used for distance computation
  !! * thisobs\%id_obs_p    - array with indices of module-type observation in process-local state vector
  !!
  !! Optional is the use of
  !! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
  !! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
  !! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
  !! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
  !!                          (default: 1=use global full observations)
  !! * thisobs\%inno_omit   - Omit obs. if squared innovation larger this factor times observation variance
  !!                          (default: 0.0, omission is active if >0)
  !! * thisobs\%inno_omit_ivar - Value of inverse variance to omit observation
  !!                          (default: 1.0e-12, change this if this value is not small compared to actual obs. error)
  !!
  !! Further variables are set when the routine PDAFomi_gather_obs is called.
  !!
    SUBROUTINE init_dim_obs_SM(step, dim_obs)

      USE mod_parallel_pdaf, &
          ONLY: mype_filter, comm_filter, npes_filter, abort_parallel, &
          mpi_integer, mpi_double_precision, mpi_in_place, mpi_sum, &
          mype_world

      USE mod_assimilation, &
          ONLY: obs_index_p, obs_filename, &
          obs_interp_indices_p, &
          obs_interp_weights_p, &
          obs_pdaf2nc, obs_nc2pdaf, &
          local_dims_obs, &
          local_disp_obs, &
          ! dim_obs_p, &
          longxy_obs_floor, latixy_obs_floor, &
          screen, cradius_SM


      USE PDAFomi, &
           ONLY: PDAFomi_gather_obs, pi
      USE mod_assimilation, &
           ONLY: obs_filename, screen

      Use mod_read_obs, &
        only: multierr, read_obs_nc_type

      use enkf_clm_mod, only: clmstatevec_allcol, clmstatevec_only_active, clmstatevec_max_layer, state_clm2pdaf_p

      use enkf_clm_mod, only: domain_def_clm

      use mod_parallel_pdaf, &
        only: abort_parallel

      use shr_kind_mod, only: r8 => shr_kind_r8

      use GridcellType, only: grc

      use clm_varcon, only: spval

      USE mod_parallel_pdaf, &
       ONLY: mype_world

      use decompMod , only : get_proc_bounds, get_proc_global

      use ColumnType , only : col

      use mod_tsmp, only: obs_interp_switch

      use enkf_clm_mod, only: get_interp_idx

      use mod_tsmp, only: point_obs

      use mod_tsmp, only: da_print_obs_index

      IMPLICIT NONE

  ! *** Arguments ***
      INTEGER, INTENT(in)    :: step       !< Current time step
      INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

  ! *** Local variables ***
      INTEGER :: i, j, c, g, cg            ! Counters
      INTEGER :: cnt                       ! Counters
      INTEGER :: cnt_interp
      INTEGER :: dim_obs_p                 ! Number of process-local observations
      REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
      REAL, ALLOCATABLE :: obs_g(:)        ! Global observation vector
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
      integer, pointer :: mycgridcell(:)

      integer :: ierror, cnt_p

      real :: deltax, deltay

      logical :: is_use_dr
      logical :: obs_snapped     !Switch for checking multiple observation counts
      logical :: newgridcell

      INTEGER :: sum_dim_obs_p

      real    :: sum_interp_weights

      character (len = 27) :: fn    !TSMP-PDAF: function name for obs_index_p output



  ! *********************************************
  ! *** Initialize full observation dimension ***
  ! *********************************************

      IF (mype_filter==0) &
        WRITE (*,*) 'Assimilate observations - obs type soil moisture'

      ! Store whether to assimilate this observation type (used in routines below)

      IF (assim_SM) thisobs%doassim = 1
      ! Specify type of distance computation

      thisobs%disttype = 3 ! 0=Cartesian, 3=geographic, distance computed with haversine formula

      ! Number of coordinates used for distance computation
      ! The distance compution starts from the first row
      thisobs%ncoord = 2


  ! **********************************
  ! *** Read PE-local observations ***
  ! **********************************


      obs_type_name = 'SM'

      ! now call function to get observations

      if (mype_filter==0 .and. screen > 2) then
        write(*,*)'load observations from type SM'
      end if
      write(current_observation_filename, '(a, i5.5)') trim(obs_filename)//'.', step


      if (mype_filter == 0) then
        call read_obs_nc_type(current_observation_filename, obs_type_name, dim_obs, obs_g, lon_obs, lat_obs, layer_obs, dr_obs, obserr, obscov)
      end if

      call mpi_bcast(dim_obs, 1, MPI_INTEGER, 0, comm_filter, ierror)

      ! check if file contains observations from type SM

      if (dim_obs == 0) then
        if (mype_filter==0 .and. screen > 2) then
          write(*,*)'TSMP-PDAF mype(w) =', mype_world, ': No observations of type SM found in file ', trim(current_observation_filename)
        end if
        dim_obs_p = 0
        ALLOCATE(obs_p(1))
        ALLOCATE(ivar_obs_p(1))
        ALLOCATE(ocoord_p(2, 1))
        ALLOCATE(thisobs%id_obs_p(1, 1))
        thisobs%infile=0
        CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
           thisobs%ncoord, cradius_SM, dim_obs)
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
        write(*,*)'Done: load observations from type SM'
      end if


      thisobs%infile = 1
      call domain_def_clm(lon_obs, lat_obs, dim_obs, longxy, latixy, longxy_obs, latixy_obs)

      ! Interpolation of measured states: Save the indices of the
      ! nearest grid points
      if (obs_interp_switch == 1) then
         ! Get the floored values for latitudes and longitudes
         call get_interp_idx(lon_obs, lat_obs, dim_obs, longxy_obs_floor, latixy_obs_floor)
      end if

#ifdef CLMFIVE
      ! Obtain CLM lon/lat information
      lon   => grc%londeg
      lat   => grc%latdeg
      ! Obtain CLM column-gridcell information
      mycgridcell => col%gridcell
#else
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      mycgridcell => clm3%g%l%c%gridcell
#endif



    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    dim_obs_p = 0

    is_use_dr = .true.

    if(allocated(thisobs%id_obs_p)) deallocate(thisobs%id_obs_p)
    allocate(thisobs%id_obs_p(1,endg-begg+1))
    thisobs%id_obs_p(1, :) = 0

    cnt = 1

    do i = 1, dim_obs
      obs_snapped = .false.
      do g = begg, endg

        newgridcell = .true.

        do c = begc,endc

          cg = mycgridcell(c)

          if(cg == g) then

            if(newgridcell) then

              if(is_use_dr) then
                  if (lon(g)>180) then
                    deltax = abs(lon(g)-lon_obs(i)-360)
                  else
                    deltax = abs(lon(g)-lon_obs(i))
                  end if
                  deltay = abs(lat(g)-lat_obs(i))
              end if

              ! Assigning observations to grid cells according to
              ! snapping distance or index arrays
              if(((is_use_dr).and.(deltax<=dr_obs(1)).and.(deltay<=dr_obs(1))).or.((.not. is_use_dr).and.(longxy_obs(i) == longxy(cnt)) .and. (latixy_obs(i) == latixy(cnt)))) then

                  dim_obs_p = dim_obs_p + 1
                  ! Use index array for setting the correct state vector index in `obs_id_p`
                  thisobs%id_obs_p(1,state_clm2pdaf_p(c,layer_obs(i))) = i

                  if (obs_snapped) then
                    print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR Observation snapped at multiple grid cells."
                    print *, "i=", i
                    if (is_use_dr) then
                      print *, "lon_obs(i)=", lon_obs(i)
                      print *, "lat_obs(i)=", lat_obs(i)
                    end if
                    call abort_parallel()
                  end if

                  ! Set observation as counted
                  obs_snapped = .true.
                  newgridcell = .false.

                  cnt=cnt+1
              end if
            end if
          end if
        end do
      end do
    end do

    if (screen > 2) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dim_obs_p=", dim_obs_p
    end if


    ! initalize OMI arrays

    IF (ALLOCATED(ivar_obs_p)) DEALLOCATE(ivar_obs_p)
    ALLOCATE(ivar_obs_p(dim_obs_p))

    IF (ALLOCATED(ocoord_p)) DEALLOCATE(ocoord_p)
    ALLOCATE(ocoord_p(2, dim_obs_p))


    ! Dimension of full observation vector
    ! ------------------------------------

    ! add and broadcast size of PE-local observation dimensions using mpi_allreduce
    call mpi_allreduce(dim_obs_p, sum_dim_obs_p, 1, MPI_INTEGER, MPI_SUM, &
        comm_filter, ierror)

    ! Check sum of dimensions of PE-local observation vectors against
    ! dimension of full observation vector
    if (.not. sum_dim_obs_p == dim_obs) then
      print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR Sum of PE-local observation dimensions"
      print *, "sum_dim_obs_p=", sum_dim_obs_p
      print *, "dim_obs=", dim_obs
      call abort_parallel()
    end if

    !  Gather PE-local observation dimensions and displacements in arrays
    ! ----------------------------------------------------------------

    ! Allocate array of PE-local observation dimensions
    IF (ALLOCATED(local_dims_obs)) DEALLOCATE(local_dims_obs)
    ALLOCATE(local_dims_obs(npes_filter))

    ! Gather array of PE-local observation dimensions
    call mpi_allgather(dim_obs_p, 1, MPI_INTEGER, local_dims_obs, 1, MPI_INTEGER, &
        comm_filter, ierror)

    ! Allocate observation displacement array local_disp_obs
    IF (ALLOCATED(local_disp_obs)) DEALLOCATE(local_disp_obs)
    ALLOCATE(local_disp_obs(npes_filter))

    ! Set observation displacement array local_disp_obs
    local_disp_obs(1) = 0
    do i = 2, npes_filter
      local_disp_obs(i) = local_disp_obs(i-1) + local_dims_obs(i-1)
    end do

    if (mype_filter==0 .and. screen > 2) then
        print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: local_disp_obs=", local_disp_obs
    end if


    ! Write index mapping array NetCDF->PDAF
    ! --------------------------------------
    ! Set index mapping `obs_pdaf2nc` between observation order in
    ! NetCDF input and observation order in pdaf as determined by domain
    ! decomposition.

    ! Use-case: Correct index order in loops over NetCDF-observation
    ! file input arrays.

    ! Trivial example: The order in the NetCDF file corresponds exactly
    ! to the order in the domain decomposition in PDAF, e.g. for a
    ! single PE per component model run.

    ! Non-trivial example: The first observation in the NetCDF file is
    ! not located in the domain/subgrid of the first PE. Rather, the
    ! second observation in the NetCDF file (`i=2`) is the only
    ! observation (`cnt = 1`) in the subgrid of the first PE
    ! (`mype_filter = 0`). This leads to a non-trivial index mapping,
    ! e.g. `obs_pdaf2nc(1)==2`:
    !
    ! i = 2
    ! cnt = 1
    ! mype_filter = 0
    !
    ! obs_pdaf2nc(local_disp_obs(mype_filter+1)+cnt) = i
    !-> obs_pdaf2nc(local_disp_obs(1)+1) = 2
    !-> obs_pdaf2nc(1) = 2


    if (allocated(obs_pdaf2nc)) deallocate(obs_pdaf2nc)
    allocate(obs_pdaf2nc(dim_obs))
    obs_pdaf2nc = 0
    if (allocated(obs_nc2pdaf)) deallocate(obs_nc2pdaf)
    allocate(obs_nc2pdaf(dim_obs))
    obs_nc2pdaf = 0


    if (point_obs==1) then

      cnt = 1
      do i = 1, dim_obs
        ! Many processes may not contain the observation / do not need
        ! to snap it, so default true
        obs_snapped = .true.

        do g = begg,endg
          newgridcell = .true.
          do c = begc,endc
            cg =   mycgridcell(c)
            if(cg == g) then
              if(newgridcell) then

                if(is_use_dr) then
                  if (lon(g)>180) then
                    deltax = abs(lon(g)-lon_obs(i)-360)
                  else
                    deltax = abs(lon(g)-lon_obs(i))
                  end if
                  deltay = abs(lat(g)-lat_obs(i))
                end if

                if(((is_use_dr).and.(deltax<=dr_obs(1)).and.(deltay<=dr_obs(1))).or.((.not. is_use_dr).and.(longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1)))) then
#ifdef CLMFIVE
                  if(state_clm2pdaf_p(c,1)==ispval) then
                    ! `ispval`: column not in state vector, most likely
                    ! because it is hydrologically inactive

                    ! Observation not snapped, even though location is
                    ! right!
                    obs_snapped = .false.

                    ! Do not use this column for snapping an
                    ! observation, instead cycle to next column
                    cycle
                  end if
#endif
                  obs_pdaf2nc(local_disp_obs(mype_filter+1)+cnt) = i
                  obs_nc2pdaf(i) = local_disp_obs(mype_filter+1)+cnt
                  cnt = cnt + 1

                  ! Observation snapped at location (possibly
                  ! overwriting a false from inactive column before)
                  obs_snapped = .true.
                end if

                newgridcell = .false.

              end if
            end if
          end do
        end do

        ! Warning, when an observation has not been snapped to any
        ! active gridcell.
        if(.not. obs_snapped) then
          print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR observations exist at non-active gridcells."
          print *, "Consider removing the following observation from the observation files."
          print *, "Observation-index in NetCDF-file: i=", i
          call abort_parallel()
        end if
      end do

    end if

    ! collect values from all PEs, by adding all PE-local arrays (works
    ! since only the subsection belonging to a specific PE is non-zero)
    call mpi_allreduce(MPI_IN_PLACE,obs_pdaf2nc,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)
    call mpi_allreduce(MPI_IN_PLACE,obs_nc2pdaf,dim_obs,MPI_INTEGER,MPI_SUM,comm_filter,ierror)

    if (mype_filter==0 .and. screen > 2) then
        print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: obs_pdaf2nc=", obs_pdaf2nc
    end if

    ! Write process-local observation arrays
    ! --------------------------------------
    IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
    ALLOCATE(obs_p(dim_obs_p))
    IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
    ALLOCATE(obs_index_p(dim_obs_p))
    if(obs_interp_switch == 1) then
        ! Array for storing indices from states that are interpolated to observation locations
        IF (ALLOCATED(obs_interp_indices_p)) DEALLOCATE(obs_interp_indices_p)
        ALLOCATE(obs_interp_indices_p(dim_obs_p, 4)) ! Later 8 for 3D / ParFlow
        IF (ALLOCATED(obs_interp_weights_p)) DEALLOCATE(obs_interp_weights_p)
        ALLOCATE(obs_interp_weights_p(dim_obs_p, 4)) ! Later 8 for 3D / ParFlow
    end if


    cnt = 1

    do i = 1, dim_obs

      do g = begg,endg
        newgridcell = .true.

        do c = begc,endc

          cg =   mycgridcell(c)

          if(cg == g) then

            if(newgridcell) then

              if(is_use_dr) then
                if (lon(g)>180) then
                  deltax = abs(lon(g)-lon_obs(i)-360)
                else
                  deltax = abs(lon(g)-lon_obs(i))
                end if
                deltay = abs(lat(g)-lat_obs(i))
              end if

              if(((is_use_dr).and.(deltax<=dr_obs(1)).and.(deltay<=dr_obs(1))).or.((.not. is_use_dr).and.(longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1)))) then
                if (thisobs%disttype/=3) then ! if haversine formula in distance calculation, the coordinates have to be converted to radians
                  ocoord_p(1,cnt) = lon_obs(i)
                  ocoord_p(2,cnt) = lat_obs(i)
                else
                  ocoord_p(1,cnt) = lon_obs(i) * pi / 180.0
                  ocoord_p(2,cnt) = lat_obs(i) * pi / 180.0
                end if

                ! Different settings of observation-location-index in
                ! state vector depending on the method of state
                ! vector assembling.
                if(clmstatevec_allcol==1) then
#ifdef CLMFIVE
                  if(clmstatevec_only_active==1) then

                    ! Error if observation deeper than clmstatevec_max_layer
                    if(layer_obs(i) > min(clmstatevec_max_layer, col%nbedrock(c))) then
                      print *, "TSMP-PDAF mype(w)=", mype_world, ": ERROR observation layer deeper than clmstatevec_max_layer or bedrock."
                      print *, "i=", i
                      print *, "c=", c
                      print *, "layer_obs(i)=", layer_obs(i)
                      print *, "col%nbedrock(c)=", col%nbedrock(c)
                      print *, "clmstatevec_max_layer=", clmstatevec_max_layer
                      call abort_parallel()
                    end if
                    obs_index_p(cnt) = state_clm2pdaf_p(c,layer_obs(i))
                  else
                    obs_index_p(cnt) = c-begc+1 + ((endc-begc+1) * (layer_obs(i)-1))
                  end if
#else
                   obs_index_p(cnt) = c-begc+1 + ((endc-begc+1) * (clmobs_layer(i)-1))
#endif
                else
#ifdef CLMFIVE
                  if(clmstatevec_only_active==1) then
                    obs_index_p(cnt) = state_clm2pdaf_p(c,layer_obs(i))
                  else
                    obs_index_p(cnt) = g-begg+1 + ((endg-begg+1) * (layer_obs(i)-1))
                  end if
#else
                   obs_index_p(cnt) = g-begg+1 + ((endg-begg+1) * (clmobs_layer(i)-1))
#endif
                end if

                !write(*,*) 'obs_index_p(',cnt,') is',obs_index_p(cnt)
                obs_p(cnt) = obs_g(i)
                if(multierr==1) ivar_obs_p(cnt) = 1/(obserr(i)*obserr(i))
                if(multierr==0) ivar_obs_p(cnt) = 1/(rms_obs_SM*rms_obs_SM)
                cnt = cnt + 1
              end if

              newgridcell = .false.

            end if

          end if

        end do
      end do

    end do

    if(obs_interp_switch==1) then
      ! loop over all obs and save the indices of the nearest grid
      ! points to array obs_interp_indices_p and save the distance
      ! weights to array obs_interp_weights_p (later normalized)
      cnt = 1
      do i = 1, dim_obs
          cnt_interp = 0
          do g = begg,endg
              ! First: latitude and longitude smaller than observation location
              if((longxy_obs_floor(i) == longxy(g-begg+1)) .and. (latixy_obs_floor(i) == latixy(g-begg+1))) then

                  obs_interp_indices_p(cnt, 1) = g-begg+1 + ((endg-begg+1) * (layer_obs(i)-1))
                  obs_interp_weights_p(cnt, 1) = sqrt(abs(lon(g)-lon_obs(i)) * abs(lon(g)-lon_obs(i)) + abs(lat(g)-lat_obs(i)) * abs(lat(g)-lat_obs(i)))
                  cnt_interp = cnt_interp + 1
              end if
              ! Second: latitude larger than observation location, longitude smaller than observation location
              if((longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs_floor(i) == latixy(g-begg+1))) then
                  obs_interp_indices_p(cnt, 2) = g-begg+1 + ((endg-begg+1) * (layer_obs(i)-1))
                  obs_interp_weights_p(cnt, 2) =sqrt(abs(lon(g)-lon_obs(i)) * abs(lon(g)-lon_obs(i)) + abs(lat(g)-lat_obs(i)) * abs(lat(g)-lat_obs(i)))
                  cnt_interp = cnt_interp + 1
              end if
              ! Third: latitude smaller than observation location, longitude larger than observation location
              if((longxy_obs_floor(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1))) then
                  obs_interp_indices_p(cnt, 3) = g-begg+1 + ((endg-begg+1) * (layer_obs(i)-1))
                  obs_interp_weights_p(cnt, 3) = sqrt(abs(lon(g)-lon_obs(i)) * abs(lon(g)-lon_obs(i)) + abs(lat(g)-lat_obs(i)) * abs(lat(g)-lat_obs(i)))
                  cnt_interp = cnt_interp + 1
              end if
              ! Fourth: latitude and longitude larger than observation location
              if((longxy_obs(i) == longxy(g-begg+1)) .and. (latixy_obs(i) == latixy(g-begg+1))) then
                  obs_interp_indices_p(cnt, 4) = g-begg+1 + ((endg-begg+1) * (layer_obs(i)-1))
                  obs_interp_weights_p(cnt, 4) = sqrt(abs(lon(g)-lon_obs(i)) * abs(lon(g)-lon_obs(i)) + abs(lat(g)-lat_obs(i)) * abs(lat(g)-lat_obs(i)))
                  cnt_interp = cnt_interp + 1
              end if
              ! Check if all four corners are found
              if(cnt_interp == 4) then
                  cnt = cnt + 1
                  ! exit
              end if
          end do
      end do

      do i = 1, dim_obs

          ! Sum of distance weights
          sum_interp_weights = sum(obs_interp_weights_p(i, :))

          do j = 1, 4
              ! Normalize distance weights
                obs_interp_weights_p(i, j) = obs_interp_weights_p(i, j) / sum_interp_weights
            end do
      end do

    end if

#ifdef PDAF_DEBUG
  IF (da_print_obs_index > 0) THEN
    ! TSMP-PDAF: For debug runs, output the state vector in files
    WRITE(fn, "(a,i5.5,a,i5.5,a)") "obs_index_p_", mype_world, ".", step, ".txt"
    OPEN(unit=71, file=fn, action="write")
    DO i = 1, dim_obs_p
      WRITE (71,"(i10)") obs_index_p(i)
    END DO
    CLOSE(71)
  END IF
#endif


  ! ****************************************
  ! *** Gather global observation arrays ***
  ! ****************************************

      CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
           thisobs%ncoord, cradius_SM, dim_obs)

  ! ********************
  ! *** Finishing up ***
  ! ********************

      ! Deallocate all local arrays
      DEALLOCATE(obs_g)
      DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    END SUBROUTINE init_dim_obs_SM



  !-------------------------------------------------------------------------------
  !> Implementation of observation operator
  !!
  !! This routine applies the full observation operator
  !! for the type of observations handled in this module.
  !!
  !! One can choose a proper observation operator from
  !! PDAFOMI_OBS_OP or add one to that module or
  !! implement another observation operator here.
  !!
  !! The routine is called by all filter processes.
  !!
    SUBROUTINE obs_op_SM(dim_p, dim_obs, state_p, ostate)

        use mod_assimilation, only: obs_index_p

        use PDAFomi_obs_f, only: PDAFomi_gather_obsstate

        IMPLICIT NONE

        ! *** Arguments ***
        INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
        INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
        REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
        REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state

        real, allocatable :: ostate_p(:)

        integer :: i



        ! ******************************************************
        ! *** Apply observation operator H on a state vector ***
        ! ******************************************************

        IF (thisobs%dim_obs_p>0) THEN
            if (allocated(ostate_p)) deallocate(ostate_p)
            ALLOCATE(ostate_p(thisobs%dim_obs_p))
        ELSE
            if (allocated(ostate_p)) deallocate(ostate_p)
            ALLOCATE(ostate_p(1))
        END IF


        DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = state_p(obs_index_p(i))
        END DO

        ! *** Global: Gather full observed state vector
        CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)

        deallocate(ostate_p)


    END SUBROUTINE obs_op_SM



  !-------------------------------------------------------------------------------
  !> Initialize local information on the module-type observation
  !!
  !! The routine is called during the loop over all local
  !! analysis domains. It has to initialize the information
  !! about local observations of the module type. It returns
  !! number of local observations of the module type for the
  !! current local analysis domain in DIM_OBS_L and the full
  !! and local offsets of the observation in the overall
  !! observation vector.
  !!
  !! This routine calls the routine PDAFomi_init_dim_obs_l
  !! for each observation type. The call allows to specify a
  !! different localization radius and localization functions
  !! for each observation type and  local analysis domain.
  !!
    SUBROUTINE init_dim_obs_l_SM(domain_p, step, dim_obs, dim_obs_l)

      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l, pi

      ! Include localization radius and local coordinates
      USE mod_assimilation, &
           ONLY: cradius_SM, locweight, sradius_SM, screen

      USE enkf_clm_mod, ONLY: state_loc2clm_c_p, clmstatevec_allcol, clmstatevec_only_active

      use shr_kind_mod, only: r8 => shr_kind_r8

      use decompMod , only : get_proc_bounds


#ifdef CLMFIVE
      USE GridcellType, ONLY: grc
      USE ColumnType, ONLY : col
#else
      USE clmtype, ONLY : clm3
#endif
      use clm_varcon, only: spval

      use mod_parallel_pdaf, &
           ONLY: mype_world



      IMPLICIT NONE

  ! *** Arguments ***
      INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
      INTEGER, INTENT(in)  :: step         !< Current time step
      INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
      INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

      REAL :: coords_l(2)      ! Coordinates of local analysis domain

      real(r8), pointer :: lon(:)
      real(r8), pointer :: lat(:)
      integer, pointer :: mycgridcell(:) !Pointer for CLM3.5/CLM5.0 col->gridcell index arrays

#ifdef CLMFIVE
      ! Obtain CLM lon/lat information
      lon   => grc%londeg
      lat   => grc%latdeg
      ! Obtain CLM column-gridcell information
      mycgridcell => col%gridcell
#else
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      mycgridcell => clm3%g%l%c%gridcell
#endif


  ! **********************************************
  ! *** Initialize local observation dimension ***
  ! **********************************************
    ! count observations within a radius

    if (thisobs%infile==1) then


      if (clmstatevec_allcol==0 .and. clmstatevec_only_active==0) then
        if (lon(state_loc2clm_c_p(domain_p))>180) then
          coords_l(1) = lon(state_loc2clm_c_p(domain_p)) - 360.0
        else
          coords_l(1) = lon(state_loc2clm_c_p(domain_p))
        end if

      else

        ! get coords_l --> coordinates of local analysis domain
        if (lon(mycgridcell(state_loc2clm_c_p(domain_p)))>180) then
        ! if SM should be assimilated. Else, state_loc2clm_c_p is not allocated
          coords_l(1) = lon(mycgridcell(state_loc2clm_c_p(domain_p))) - 360.0
        else
          coords_l(1) = lon(mycgridcell(state_loc2clm_c_p(domain_p)))
        end if
        coords_l(2) = lat(mycgridcell(state_loc2clm_c_p(domain_p)))

      end if

      if (thisobs%disttype==3) then ! if haversine formula in distance calculation, the coordinates have to be converted to radians
        coords_l(1) = coords_l(1) * pi / 180.0
        coords_l(2) = coords_l(2) * pi / 180.0
      end if

    else

      coords_l(1) = spval
      coords_l(2) = spval

    end if

    ! for disttype=3, the cradius and sradius have to passed in meters, so I multiply by 1000 to be able to put it in km in the input file

    if (thisobs%disttype==3) then
      CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
          locweight, cradius_SM*1000.0, sradius_SM*1000.0, dim_obs_l)
    else
      CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
          locweight, cradius_SM, sradius_SM, dim_obs_l)
    end if


    END SUBROUTINE init_dim_obs_l_SM



  !-------------------------------------------------------------------------------
  !> Perform covariance localization for local EnKF on the module-type observation
  !!
  !! The routine is called in the analysis step of the localized
  !! EnKF. It has to apply localization to the two matrices
  !! HP and HPH of the analysis step for the module-type
  !! observation.
  !!
  !! This routine calls the routine PDAFomi_localize_covar
  !! for each observation type. The call allows to specify a
  !! different localization radius and localization functions
  !! for each observation type.
  !!
    SUBROUTINE localize_covar_SM(dim_p, dim_obs, HP_p, HPH, coords_p)

      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_localize_covar

      ! Include localization radius and local coordinates
      USE mod_assimilation, &
           ONLY: cradius_SM, locweight, sradius_SM

      use enkf_clm_mod, only: state_pdaf2clm_c_p

      use shr_kind_mod, only: r8 => shr_kind_r8

      USE GridcellType, ONLY: grc
      USE ColumnType, ONLY : col

      IMPLICIT NONE

  ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
      REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
      REAL, INTENT(inout) :: coords_p(:,:)         !< Coordinates of state vector elements

      integer :: i

      real(r8), pointer :: lon(:)
      real(r8), pointer :: lat(:)
      integer, pointer :: mycgridcell(:) !Pointer for CLM3.5/CLM5.0 col->gridcell index arrays

#ifdef CLMFIVE
      ! Obtain CLM lon/lat information
      lon   => grc%londeg
      lat   => grc%latdeg
      ! Obtain CLM column-gridcell information
      mycgridcell => col%gridcell
#else
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
      mycgridcell => clm3%g%l%c%gridcell
#endif


  ! *************************************
  ! *** Apply covariance localization ***
  ! *************************************

      do i = 1,dim_p
        if (lon(mycgridcell(state_pdaf2clm_c_p(i)))>180) then
          coords_p(1,i) = lon(mycgridcell(state_pdaf2clm_c_p(i))) - 360.0
        else
          coords_p(1,i) = lon(mycgridcell(state_pdaf2clm_c_p(i)))
        end if
        coords_p(2,i) = lat(mycgridcell(state_pdaf2clm_c_p(i)))
      end do



      CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, cradius_SM, sradius_SM, &
           coords_p, HP_p, HPH)

    END SUBROUTINE localize_covar_SM


  END MODULE obs_SM_pdafomi






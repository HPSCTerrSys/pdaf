!> PDAF-OMI observation module for type GRACE observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! OBSTYPE = GRACE
!!
!! __Observation type GRACE:__
!! The observation type GRACE are TWSA observations (gridded)
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
MODULE obs_GRACE_pdafomi

    USE mod_parallel_pdaf, &
         ONLY: mype_filter    ! Rank of filter process
    USE PDAFomi, &
         ONLY: obs_f, obs_l   ! Declaration of observation data types

    IMPLICIT NONE
    SAVE

    PUBLIC

    ! Variables which are inputs to the module (usually set in init_pdaf)
    LOGICAL :: assim_GRACE        !< Whether to assimilate this data type
    REAL    :: rms_obs_GRACE      !< Observation error standard deviation (for constant errors)
    logical, allocatable :: vec_useObs(:)
    integer, allocatable :: vec_numPoints_global(:) ! vector of number of points for each GRACE observation, same dimension as observation vector
    logical, allocatable :: vec_useObs_global(:) ! vector that tells if an observation of used (1) or not (0), same dimension as observation vector, global

    real, allocatable :: tws_temp_mean(:,:) ! temporal mean for TWS
    real, allocatable :: tws_temp_mean_d(:) ! temporal mean for TWS, vectorized with the same bounds as local process
    real, allocatable :: lon_temp_mean(:,:) ! corresponding longitude
    real, allocatable :: lat_temp_mean(:,:) ! corresponding latitude

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
    SUBROUTINE init_dim_obs_GRACE(step, dim_obs)

      USE PDAFomi, &
           ONLY: PDAFomi_gather_obs
      USE mod_assimilation, &
           ONLY: filtertype, cradius_GRACE, obs_filename, temp_mean_filename, screen

        use mod_read_obs, only: read_obs_nc_type, domain_def_clm, multierr

      use enkf_clm_mod, only: num_layer, hactiveg_levels

      use mod_parallel_pdaf, &
        only: mpi_integer, mpi_sum, mpi_2integer, mpi_maxloc, comm_filter

      use shr_kind_mod, only: r8 => shr_kind_r8

      use GridcellType, only: grc

      use clm_varcon, only: spval

      USE mod_parallel_pdaf, &
       ONLY: mype_world

      use decompMod , only : get_proc_bounds

      IMPLICIT NONE

  ! *** Arguments ***
      INTEGER, INTENT(in)    :: step       !< Current time step
      INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

  ! *** Local variables ***
      INTEGER :: i, j, count_points, c, l, k     ! Counters
      INTEGER :: cnt, cnt0                 ! Counters
      INTEGER :: dim_obs_p                 ! Number of process-local observations
      REAL, ALLOCATABLE :: obs_field(:,:)  ! Observation field read from file
      REAL, ALLOCATABLE :: obs_p(:)        ! PE-local observation vector
      REAL, ALLOCATABLE :: obs_g(:)        ! Global observation vector
      REAL, ALLOCATABLE :: ivar_obs_p(:)   ! PE-local inverse observation error variance
      REAL, ALLOCATABLE :: ocoord_p(:,:)   ! PE-local observation coordinates
      CHARACTER(len=2) :: stepstr          ! String for time step
      character (len = 110) :: current_observation_filename


      character(len=20) :: obs_type_name ! name of observation type (e.g. GRACE, SM, ST, ...)

      integer :: numPoints ! minimum number of points so that the GRACE observation is used
      integer, allocatable :: vec_numPoints(:) ! number of model grid cells that are in a radius of dr around the GRACE observation
      INTEGER, allocatable :: in_mpi(:,:), out_mpi(:,:)
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

      real(r8), pointer :: lon(:)
      real(r8), pointer :: lat(:)

      real(r8) :: pi

      integer :: ierror, cnt_p

      real :: deltax, deltay



  ! *********************************************
  ! *** Initialize full observation dimension ***
  ! *********************************************

      !IF (mype_filter==0) &
      IF (mype_filter==0) &
        WRITE (*,*) 'Assimilate observations - obs type GRACE'

      ! Store whether to assimilate this observation type (used in routines below)

      IF (assim_GRACE) thisobs%doassim = 1
      ! Specify type of distance computation
      thisobs%disttype = 0   ! 0=Cartesian

      ! Number of coordinates used for distance computation
      ! The distance compution starts from the first row
      thisobs%ncoord = 2


  ! **********************************
  ! *** Read PE-local observations ***
  ! **********************************

      ! read observations from nc file --> call function in mod_read_obs. Idea: when you have multiple observation types in one file,
      ! also pass the observation type, here 'GRACE' to the function. You have to give each observation in the file an information which type
      ! it is. This way, the output in this function is only the GRACE observation (or soil moisture,...; dependent on what you want to implement)


      obs_type_name = 'GRACE'

      ! now call function to get observations

      if (mype_filter==0 .and. screen > 2) then
        write(*,*)'load observations from type GRACE'
      end if
      write(current_observation_filename, '(a, i5.5)') trim(obs_filename)//'.', step
      call read_obs_nc_type(current_observation_filename, obs_type_name, dim_obs, obs_g, lon_obs, lat_obs, layer_obs, dr_obs, obserr, obscov)
      if (mype_filter==0 .and. screen > 2) then
        write(*,*)'Done: load observations from type GRACE'
      end if

      if (dim_obs == 0) then
        if (mype_filter==0 .and. screen > 2) then
          write(*,*)'TSMP-PDAF mype(w) =', mype_world, ': No observations of type GRACE found in file ', trim(current_observation_filename)
        end if
        dim_obs_p = 0
        ALLOCATE(obs_p(1))
        ALLOCATE(ivar_obs_p(1))
        ALLOCATE(ocoord_p(2, 1))
        ALLOCATE(thisobs%id_obs_p(1, 1))
        thisobs%infile=0
        CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
           thisobs%ncoord, cradius_GRACE, dim_obs)
        return
      end if
      thisobs%infile=1

  ! ***********************************************************
  ! *** Count available observations for the process domain ***
  ! *** and initialize index and coordinate arrays.         ***
  ! ***********************************************************

    call domain_def_clm(lon_obs, lat_obs, dim_obs, longxy, latixy, longxy_obs, latixy_obs)

    IF (ALLOCATED(vec_useObs)) DEALLOCATE(vec_useObs)
    ALLOCATE(vec_useObs(dim_obs))
    IF (ALLOCATED(vec_numPoints)) DEALLOCATE(vec_numPoints)
    ALLOCATE(vec_numPoints(dim_obs))
    vec_numPoints=0
    IF (ALLOCATED(vec_numPoints_global)) DEALLOCATE(vec_numPoints_global)
    ALLOCATE(vec_numPoints_global(dim_obs))
    IF (ALLOCATED(vec_useObs_global)) DEALLOCATE(vec_useObs_global)
    ALLOCATE(vec_useObs_global(dim_obs))
    vec_useObs_global = .true.
    IF (ALLOCATED(in_mpi)) DEALLOCATE(in_mpi)
    ALLOCATE(in_mpi(2,dim_obs))
    IF (ALLOCATED(out_mpi)) DEALLOCATE(out_mpi)
    ALLOCATE(out_mpi(2,dim_obs))


    ! additions for GRACE assimilation, it can be the case that not enough CLM gridpoints lie in the neighborhood of a GRACE observation
    ! if this is the case, the GRACE observations cannot be reproduced in a satisfactory manner and is not used in the assimilation
    ! count grdicells that are in a certain radius. This effect is especially present when the applied GRACE resolution is high or for
    ! observation lying directly at the coast

    lon   => grc%londeg
    lat   => grc%latdeg

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)


    do i = 1, dim_obs

        count_points = 0

        do c = 1, num_layer(1)

            j = hactiveg_levels(c,1)
            deltax = abs(longxy(c)-longxy_obs(i))
            deltay = abs(latixy(c)-latixy_obs(i))

            if((sqrt(real(deltax)**2 + real(deltay)**2)) < dr_obs(1)/0.11) then
                count_points = count_points+1
            end if

        end do

        vec_numPoints(i) = count_points

    end do

    ! get vec_numPoints from all processes and add them up together via mpi_allreduce
    ! then we have for each observation the number of all points that lie within their neighborhood, not only from one process

    call mpi_allreduce(vec_numPoints,vec_numPoints_global, dim_obs, mpi_integer, mpi_sum, comm_filter, ierror)

    ! only observations should be used that "see" enough gridcells
    pi = 3.14159265358979323846
    numPoints = int(ceiling((pi*(dr_obs(1)/0.11)**2)/2))
    if (screen > 2) then
        if (mype_filter==0) then
            print *, "Minimum number of points for using one observation is ", numPoints
        end if
    end if

    vec_useObs_global = merge(vec_useObs_global,.false.,vec_numPoints_global>=numPoints)
    vec_useObs = vec_useObs_global

    if (screen > 2) then
        if (mype_filter==0) then
           print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_f_pdaf: vec_useObs_global=", vec_useObs_global
           print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_f_pdaf: vec_numPoints_global=", vec_numPoints_global
        end if
    end if

    ! The observation should be reproduced by the process that has the most grid points inside the observation radius
    ! However, in the observation operator, all points will be used by using mpi

    in_mpi(1,:) = vec_numPoints
    in_mpi(2,:) = mype_filter

    call mpi_allreduce(in_mpi,out_mpi, dim_obs, mpi_2integer, mpi_maxloc, comm_filter, ierror)

    vec_useObs = merge(vec_useObs,.false.,out_mpi(2,:)==mype_filter)

    IF (ALLOCATED(in_mpi)) DEALLOCATE(in_mpi)
    IF (ALLOCATED(out_mpi)) DEALLOCATE(out_mpi)

    dim_obs_p = count(vec_useObs)

    if(allocated(thisobs%id_obs_p)) deallocate(thisobs%id_obs_p)
    ALLOCATE(thisobs%id_obs_p(1, begg:endg))
    thisobs%id_obs_p(1, :) = 0

    do i = 1, dim_obs

        if (vec_useObs_global(i)) then

            do c = 1, num_layer(1)

                j = hactiveg_levels(c,1)
                deltax = abs(longxy(c)-longxy_obs(i))
                deltay = abs(latixy(c)-latixy_obs(i))

                if((sqrt(real(deltax)**2 + real(deltay)**2))<=dr_obs(1)/0.11) then
                    thisobs%id_obs_p(1, j) = i
                end if

            end do

        end if

    end do



    IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
    ALLOCATE(obs_p(dim_obs_p))

    IF (ALLOCATED(ivar_obs_p)) DEALLOCATE(ivar_obs_p)
    ALLOCATE(ivar_obs_p(dim_obs_p))

    IF (ALLOCATED(ocoord_p)) DEALLOCATE(ocoord_p)
    ALLOCATE(ocoord_p(2, dim_obs_p))


    if (multierr==1) ivar_obs_p = pack(1/obserr, vec_useObs)

    cnt_p = 1
    do i = 1, dim_obs
        if (vec_useObs(i)) then
            obs_p(cnt_p) = obs_g(i)
            ocoord_p(1,cnt_p) = longxy_obs(i)
            ocoord_p(2,cnt_p) = latixy_obs(i)
            cnt_p = cnt_p+1
        end if
    end do

    dim_obs = count(vec_useObs_global)
  ! ****************************************
  ! *** Gather global observation arrays ***
  ! ****************************************

      CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
           thisobs%ncoord, cradius_GRACE, dim_obs)

  ! ********************
  ! *** Finishing up ***
  ! ********************

        ! load temporal mean of TWS and vectorize for process from begg to endg --> only has to be done once as the process does not change bounds

        if (.not. allocated(tws_temp_mean_d)) then
            call read_temp_mean_model(temp_mean_filename)

            if (allocated(tws_temp_mean_d)) DEALLOCATE(tws_temp_mean_d)
            ALLOCATE(tws_temp_mean_d(begg:endg))
            tws_temp_mean_d(:) = spval

            do j = begg,endg
                outer: do l = 1,size(lon_temp_mean,1)
                    do k=1,size(lon_temp_mean,2)
                        if (lon_temp_mean(l,k)==lon(j) .and. lat_temp_mean(l,k)==lat(j)) then
                            tws_temp_mean_d(j) = tws_temp_mean(l,k)
                            exit outer
                        end if
                    end do
                end do outer

                if (lon(j)/=lon_temp_mean(l,k) .or. lat(j)/=lat_temp_mean(l,k)) then
                    print *, "Attention: distributing model mean to clumps does not work properly"
                    print *, "idx_lon= ",l, "idx_lat= ",k
                    print *, "lon(j)= ", lon(j),"lon_temp_mean(idx_lon)= ",lon_temp_mean(l,k)
                    print *, "lat(j)= ", lat(j),"lat_temp_mean(idx_lat)= ",lat_temp_mean(l,k)
                    stop
                end if
            end do

            deallocate(tws_temp_mean)
            deallocate(lon_temp_mean)
            deallocate(lat_temp_mean)

        end if

      ! Deallocate all local arrays
      DEALLOCATE(obs_g)
      DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    END SUBROUTINE init_dim_obs_GRACE



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
    SUBROUTINE obs_op_GRACE(dim_p, dim_obs, state_p, ostate)

        use enkf_clm_mod, &
                only: clm_varsize_tws, state_setup, num_layer, hactiveg_levels

        use mod_parallel_pdaf, &
                only: comm_filter, mpi_double_precision, mpi_sum

        use clm_varpar   , only : nlevsoi

        use decompMod , only : get_proc_bounds

        use clm_varcon, only: spval

        use shr_kind_mod, only: r8 => shr_kind_r8

        use mod_assimilation, only: screen

        use PDAFomi_obs_f, only: PDAFomi_gather_obsstate

        IMPLICIT NONE

        ! *** Arguments ***
        INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
        INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
        REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
        REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state

        REAL, allocatable :: tws_from_statevector(:)
        real, allocatable :: ostate_p(:)

        integer :: begp, endp   ! per-proc beginning and ending pft indices
        integer :: begc, endc   ! per-proc beginning and ending column indices
        integer :: begl, endl   ! per-proc beginning and ending landunit indices
        integer :: begg, endg   ! per-proc gridcell ending gridcell indices
        integer :: ierror
        integer :: obs_point    ! which observation is seen by which point?

        REAL:: m_state_sum(size(vec_useObs_global))
        REAL:: m_state_sum_global(size(vec_useObs_global))

        integer :: count, j, g

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

        if (thisobs%infile == 1) then ! as I need also tasks with dim_obs_p==0 to reproduce GRACE observations
            ! I do not check whether dim_obs_p>0 but if I want to assimilate GRACE. A check still has to be included,
            ! else there will be segmentation faults

            m_state_sum(:) = 0

            call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

            if (allocated(tws_from_statevector)) deallocate(tws_from_statevector)
            allocate(tws_from_statevector(begg:endg))

            tws_from_statevector(begg:endg) = spval

            select case(state_setup)
            case(0) ! all compartments included in state vector

                do j = 1,nlevsoi
                    do count = 1, num_layer(j)
                        g = hactiveg_levels(count,j)
                        if (j==1) then
                            tws_from_statevector(g) = 0._r8
                        end if
                        if (j==1) then
                            ! liq + ice
                            tws_from_statevector(g) = tws_from_statevector(g) + state_p(count)
                            ! snow
                            tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2))
                            !surface water
                            tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2)+ clm_varsize_tws(3))
                        else
                            ! liq + ice
                            tws_from_statevector(g) = tws_from_statevector(g) + state_p(count+sum(num_layer(1:j-1)))
                        end if
                    end do
                end do

                ! do count = 1, num_hactiveg_patch
                !     g = hactiveg_patch(count)
                !     ! canopy water
                !     tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2)+ clm_varsize_tws(3)+ clm_varsize_tws(4))
                ! end do


            case(1) ! TWS in statevector

                do count = 1, num_layer(1)
                    g = hactiveg_levels(count,1)
                    tws_from_statevector(g) = state_p(count)
                end do

            case(2) ! snow and soil moisture aggregated over surface, root zone and deep soil moisture in state vector

                do count = 1,num_layer(1)
                    g = hactiveg_levels(count,1)
                    tws_from_statevector(g) = state_p(count)
                    ! snow added for first layer
                    tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2) + clm_varsize_tws(3))
                end do

                do count = 1,num_layer(4)
                    g = hactiveg_levels(count,4)
                    tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1))
                end do

                do count = 1,num_layer(13)
                    g = hactiveg_levels(count,13)
                    tws_from_statevector(g) = tws_from_statevector(g) + state_p(count + clm_varsize_tws(1) + clm_varsize_tws(2))
                end do

            end select



            ! subtract mean TWS to obtain TWSA model values
            do count = 1, num_layer(1)
                g = hactiveg_levels(count,1)
                obs_point = thisobs%id_obs_p(1, g)
                if (obs_point /= 0) then
                    if (tws_temp_mean_d(g)/=spval .and. tws_from_statevector(g)/=spval) then
                        m_state_sum(obs_point) = m_state_sum(obs_point) + tws_from_statevector(g)-tws_temp_mean_d(g)
                    else if (tws_temp_mean_d(g)==spval .and. .not. tws_from_statevector(g)==spval) then
                        print*, "error, tws temporal mean is spval and reproduced values is not spval for g = ", g
                        print*, "reproduced = ", tws_from_statevector(g)
                        stop
                    else if (.not. tws_temp_mean_d(g)==spval .and. tws_from_statevector(g)==spval) then
                        print*, "error, tws temporal mean is not spval and reproduced values is spval for g = ", g
                        print*, "temp_mean = ", tws_temp_mean_d(g)
                        stop
                    end if
                end if
            end do

            ! now get the sum of m_state_sum (sum over all TWSA values for the gridcells for one process) over all processes
            call mpi_allreduce(m_state_sum, m_state_sum_global, size(vec_useObs_global), mpi_double_precision, mpi_sum, comm_filter, ierror)
            m_state_sum_global = m_state_sum_global/vec_numPoints_global

            ostate_p = pack(m_state_sum_global, vec_useObs)

        end if

        ! *** Global: Gather full observed state vector
        CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)
        deallocate(ostate_p)
        if (screen>2 .and. mype_filter==0 .and. thisobs%dim_obs_f>0) then
            write(*,*)'m_state_sum_global = ', m_state_sum_global
        end if

    END SUBROUTINE obs_op_GRACE



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
    SUBROUTINE init_dim_obs_l_GRACE(domain_p, step, dim_obs, dim_obs_l)

      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l

      ! Include localization radius and local coordinates
      USE mod_assimilation, &
           ONLY: cradius_GRACE, locweight, sradius_GRACE

        use clm_varcon, only: spval

      IMPLICIT NONE

  ! *** Arguments ***
      INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
      INTEGER, INTENT(in)  :: step         !< Current time step
      INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
      INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

      REAL :: coords_l(2)      ! Coordinates of local analysis domain


  ! **********************************************
  ! *** Initialize local observation dimension ***
  ! **********************************************
    ! count observations within a radius

    if (thisobs%infile==1) then

      ! get coords_l --> coordinates of local analysis domain, I can just set this to the rotated CLM coordinates
      coords_l(1) = longxy(domain_p)
      coords_l(2) = latixy(domain_p)

    else
      coords_l(1) = spval
      coords_l(2) = spval
    end if

      CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
           locweight, cradius_GRACE, sradius_GRACE, dim_obs_l)

    END SUBROUTINE init_dim_obs_l_GRACE



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
    SUBROUTINE localize_covar_GRACE(dim_p, dim_obs, HP_p, HPH, coords_p)

      ! Include PDAFomi function
      USE PDAFomi, ONLY: PDAFomi_localize_covar

      ! Include localization radius and local coordinates
      USE mod_assimilation, &
           ONLY: cradius_GRACE, locweight, sradius_GRACE

      use enkf_clm_mod, only: gridcell_state

      IMPLICIT NONE

  ! *** Arguments ***
      INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
      INTEGER, INTENT(in) :: dim_obs               !< Dimension of observation vector
      REAL, INTENT(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
      REAL, INTENT(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
      REAL, INTENT(inout) :: coords_p(:,:)         !< Coordinates of state vector elements

      integer :: i


  ! *************************************
  ! *** Apply covariance localization ***
  ! *************************************

      do i = 1,dim_p
        coords_p(1,i) = longxy(gridcell_state(i))
        coords_p(2,i) = latixy(gridcell_state(i))
      end do



      CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, cradius_GRACE, sradius_GRACE, &
           coords_p, HP_p, HPH)

    END SUBROUTINE localize_covar_GRACE



    subroutine read_temp_mean_model(temp_mean_filename)

        use netcdf, only: nf90_open
        use netcdf, only: nf90_inq_dimid
        use netcdf, only: nf90_inquire_dimension
        use netcdf, only: nf90_inq_varid
        use netcdf, only: nf90_get_var
        use netcdf, only: nf90_close
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

    end subroutine read_temp_mean_model

  END MODULE obs_GRACE_pdafomi





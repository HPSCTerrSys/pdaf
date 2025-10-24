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
!mod_read_obs.F90: Module for reading of observation files
!-------------------------------------------------------------------------------------------

module mod_read_obs
  use iso_C_binding, only: c_int, c_ptr, c_loc

  implicit none

  public

  integer, allocatable :: idx_obs_nc(:)
  integer, allocatable :: x_idx_obs_nc(:)
  integer, allocatable :: y_idx_obs_nc(:)
  integer, allocatable :: z_idx_obs_nc(:)
  integer, allocatable :: var_id_obs_nc(:,:)
  real, allocatable :: x_idx_interp_d_obs_nc(:)
  real, allocatable :: y_idx_interp_d_obs_nc(:)

  integer(c_int), allocatable,target :: idx_obs_pf(:)
  integer(c_int), allocatable,target :: x_idx_obs_pf(:)
  integer(c_int), allocatable,target :: y_idx_obs_pf(:)
  integer(c_int), allocatable,target :: z_idx_obs_pf(:)
  integer(c_int), allocatable,target :: ind_obs_pf(:)
  type(c_ptr),bind(C,name="tidx_obs") :: ptr_tidx_obs
  type(c_ptr),bind(C,name="xidx_obs") :: ptr_xidx_obs
  type(c_ptr),bind(C,name="yidx_obs") :: ptr_yidx_obs
  type(c_ptr),bind(C,name="zidx_obs") :: ptr_zidx_obs
  type(c_ptr),bind(C,name="ind_obs")  :: ptr_ind_obs

  !kuw: obs variables for clm
  real, allocatable :: clmobs_lon(:)
  real, allocatable :: clmobs_lat(:)
  integer, allocatable :: clmobs_layer(:)
  real, allocatable :: clmobs_dr(:) ! snapping distance for clm obs
  real, allocatable :: clm_obs(:)
  real, allocatable :: clm_obserr(:)
  !kuw end

  real, allocatable :: pressure_obs(:)
  real, allocatable :: pressure_obserr(:)

  ! Flag: Use vector of observation errors in observation file
  integer :: multierr=0
  integer :: dim_nx, dim_ny
!  integer :: crns_flag=0   !hcp
  real, allocatable :: dampfac_state_time_dependent_in(:)
  real, allocatable :: dampfac_param_time_dependent_in(:)
contains


  !> @author Yorck Ewerdwalbesloh
  !> @date 17.03.2025
  !> @brief Read NetCDF observation file for different observation types to be able to use one full observation files with several types
  !> @param[in] Name of observation file, Name of observation type
  !> @param[inout] Full observation dimension, full observation vector, uncertainty information, coordinates (lon and lat)
  !> @details
  !> This subroutine reads the observation file and return the data

  subroutine read_obs_nc_type(current_observation_filename, current_observation_type, dim_obs_g, obs_g, lon_obs_g, lat_obs_g, layer_obs_g, dr_obs_g, obserr_g, obscov_g)
    use netcdf, only: nf90_max_name
    use netcdf, only: nf90_open
    use netcdf, only: nf90_nowrite
    use netcdf, only: nf90_inq_dimid
    use netcdf, only: nf90_inquire_dimension
    use netcdf, only: nf90_inq_varid
    use netcdf, only: nf90_noerr
    use netcdf, only: nf90_get_var
    use mod_assimilation, only: screen
    implicit none

    integer :: ncid
    integer :: dimid, status
    integer :: haserr

    character (len = *), intent(in) :: current_observation_filename
    character (len = *), intent(in) :: current_observation_type
    INTEGER, INTENT(inout) :: dim_obs_g    !< Dimension of full observation vector
    REAL, allocatable, INTENT(inout) :: obs_g(:)
    REAL, allocatable, INTENT(inout) :: lon_obs_g(:)
    REAL, allocatable, INTENT(inout) :: lat_obs_g(:)
    INTEGER, allocatable, INTENT(inout) :: layer_obs_g(:)
    REAL, allocatable, INTENT(inout) :: dr_obs_g(:)
    REAL, allocatable, INTENT(inout) :: obserr_g(:)
    REAL, allocatable, INTENT(inout) :: obscov_g(:,:)

    integer :: clmobs_varid, dr_varid,  clmobs_lon_varid,  clmobs_lat_varid,  &
          clmobs_layer_varid, clmobserr_varid, clmobscov_varid, obstype_varid

    character (len = *), parameter :: dim_name = "dim_obs"
    character (len = *), parameter :: obs_name   = "obs_clm"
    character (len = *), parameter :: dr_name    = "dr"
    character (len = *), parameter :: lon_name   = "lon"
    character (len = *), parameter :: lat_name   = "lat"
    character (len = *), parameter :: layer_name = "layer"
    character (len = *), parameter :: obserr_name   = "obserr_clm"
    character (len = *), parameter :: obscov_name   = "obscov_clm"
    character (len = *), parameter :: type_name   = "type_clm"
    character(len = nf90_max_name) :: RecordDimName

    character(len=20), allocatable :: obs_type(:)
    integer, allocatable :: indices(:)

    integer :: has_obs_clm, dim_obs



    call check( nf90_open(current_observation_filename, nf90_nowrite, ncid) )
    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))

    has_obs_clm = nf90_inq_varid(ncid, obs_name, clmobs_varid)

    if(has_obs_clm == nf90_noerr) then

        if(allocated(obs_type))   deallocate(obs_type)
        allocate(obs_type(dim_obs))

        call check(nf90_inq_varid(ncid, type_name, obstype_varid))
        call check(nf90_get_var(ncid, obstype_varid, obs_type))

        if (trim(obs_type(1)) /= trim(current_observation_type)) then

            dim_obs_g = 0

            if(allocated(obs_g))   deallocate(obs_g)
            allocate(obs_g(dim_obs_g))

            if(allocated(lon_obs_g))   deallocate(lon_obs_g)
            allocate(lon_obs_g(dim_obs_g))

            if(allocated(lat_obs_g))   deallocate(lat_obs_g)
            allocate(lat_obs_g(dim_obs_g))

            if(allocated(layer_obs_g))   deallocate(layer_obs_g)
            allocate(layer_obs_g(dim_obs_g))

            if(allocated(dr_obs_g))   deallocate(dr_obs_g)
            allocate(dr_obs_g(dim_obs_g))

        else

            if(allocated(obs_g))   deallocate(obs_g)
            allocate(obs_g(dim_obs))

            if(allocated(lon_obs_g))   deallocate(lon_obs_g)
            allocate(lon_obs_g(dim_obs))

            if(allocated(lat_obs_g))   deallocate(lat_obs_g)
            allocate(lat_obs_g(dim_obs))

            if(allocated(layer_obs_g))   deallocate(layer_obs_g)
            allocate(layer_obs_g(dim_obs))

            if(allocated(dr_obs_g))   deallocate(dr_obs_g)
            allocate(dr_obs_g(dim_obs))

            call check(nf90_get_var(ncid, clmobs_varid, obs_g))

            call check( nf90_inq_varid(ncid, lon_name, clmobs_lon_varid) )
            call check( nf90_get_var(ncid, clmobs_lon_varid, lon_obs_g) )

            call check( nf90_inq_varid(ncid, lat_name, clmobs_lat_varid) )
            call check( nf90_get_var(ncid, clmobs_lat_varid, lat_obs_g) )

            haserr = nf90_inq_varid(ncid, layer_name, clmobs_layer_varid)
            if (haserr == nf90_noerr) then
                call check( nf90_get_var(ncid, clmobs_layer_varid, layer_obs_g) )
            else
                layer_obs_g(:) = 1
            end if

            call check( nf90_inq_varid(ncid, dr_name, dr_varid) )
            call check( nf90_get_var(ncid, dr_varid, dr_obs_g) )

            haserr = nf90_inq_varid(ncid, obserr_name, clmobserr_varid)

            if(haserr == nf90_noerr) then

                multierr = 1

                if(allocated(obserr_g))   deallocate(obserr_g)
                allocate(obserr_g(dim_obs))

                call check(nf90_get_var(ncid, clmobserr_varid, obserr_g))

            end if

            dim_obs_g = dim_obs

        end if

    end if

  end subroutine read_obs_nc_type





  !> @author Wolfgang Kurtz, Guowei He, Mukund Pondkule
  !> @date 03.03.2023
  !> @brief Read NetCDF observation file
  !> @param[in] current_observation_filename Name of observation file
  !> @details
  !> This subroutine reads the observation file and stores the data in the
  !> corresponding arrays.
  subroutine read_obs_nc(current_observation_filename)
    USE mod_assimilation, &
        ! ONLY: obs_p, obs_index_p, dim_obs, obs_filename, screen
        ONLY: dim_obs, screen
    use mod_parallel_pdaf, &
         only: mype_world !, mpi_info_null
    ! use mod_parallel_pdaf, &
    !      only: comm_filter
    use mod_tsmp, &
        only: point_obs, obs_interp_switch, is_dampfac_state_time_dependent, &
        is_dampfac_param_time_dependent, crns_flag
    use netcdf, only: nf90_max_name
    use netcdf, only: nf90_open
    use netcdf, only: nf90_nowrite
    use netcdf, only: nf90_inq_dimid
    use netcdf, only: nf90_inquire_dimension
    use netcdf, only: nf90_inq_varid
    use netcdf, only: nf90_get_var
    use netcdf, only: nf90_noerr
    use netcdf, only: nf90_strerror
    use netcdf, only: nf90_close
    implicit none
    integer :: ncid
    character (len = *), parameter :: dim_name = "dim_obs"
    integer :: var_id_varid !, x, y
    integer :: damp_state_varid
    integer :: damp_param_varid
    ! integer :: comm, omode, info
    character (len = *), parameter :: dim_nx_name = "dim_nx"
    character (len = *), parameter :: dim_ny_name = "dim_ny"
    character (len = *), parameter :: var_id_name = "var_id"
    character (len = *), parameter :: damp_state_name = "dampfac_state"
    character (len = *), parameter :: damp_param_name = "dampfac_param"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status
    integer :: haserr
    integer :: has_damping_state
    integer :: has_damping_param
    ! This is the name of the data file we will read.
    character (len = *), intent(in) :: current_observation_filename

    ! ParFlow
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
    integer :: pres_varid,presserr_varid, &
        idx_varid,  x_idx_varid, y_idx_varid,  z_idx_varid, &
        depth_varid, &
        x_idx_interp_d_varid, y_idx_interp_d_varid
    character (len = *), parameter :: pres_name = "obs_pf"
    character (len = *), parameter :: presserr_name = "obserr_pf"
    character (len = *), parameter :: idx_name = "idx"
    character (len = *), parameter :: x_idx_name = "ix"
    character (len = *), parameter :: y_idx_name = "iy"
    character (len = *), parameter :: z_idx_name = "iz"
    character (len = *), parameter :: depth_name = "depth"
    character (len = *), parameter :: x_idx_interp_d_name = "ix_interp_d"
    character (len = *), parameter :: y_idx_interp_d_name = "iy_interp_d"
    integer :: has_obs_pf
    integer :: has_depth
#endif
#endif

    ! CLM
#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
    integer :: clmobs_varid, dr_varid,  clmobs_lon_varid,  clmobs_lat_varid,  &
        clmobs_layer_varid, clmobserr_varid
    character (len = *), parameter :: obs_name   = "obs_clm"
    character (len = *), parameter :: dr_name    = "dr"
    character (len = *), parameter :: lon_name   = "lon"
    character (len = *), parameter :: lat_name   = "lat"
    character (len = *), parameter :: layer_name = "layer"
    character (len = *), parameter :: obserr_name   = "obserr_clm"
    integer :: has_obs_clm
#endif
#endif

    if (screen > 2) then
        print *, "TSMP-PDAF mype(w)=", mype_world, ": read_obs_nc"
        print *, "TSMP-PDAF mype(w)=", mype_world, ": current_observation_filename=", current_observation_filename
    end if

    ! Observation file dimension
    ! --------------------------
    call check( nf90_open(current_observation_filename, nf90_nowrite, ncid) )
    !call check( nf90_open_par(current_observation_filename, nf90_nowrite, comm_filter, mpi_info_null, ncid) )

    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))
    if (screen > 2) then
        print *, "TSMP-PDAF mype(w)=", mype_world, ": dim_obs=", dim_obs
    end if

    ! Multi-scale data assimilation
    ! ----------------------------
    ! Not point observations, see TSMP-PDAF manual entry for input `point_obs`
    if(point_obs == 0) then
        call check(nf90_inq_dimid(ncid, dim_nx_name, dimid))
        call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_nx))
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": dim_nx=", dim_nx
        end if

        call check(nf90_inq_dimid(ncid, dim_ny_name, dimid))
        call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_ny))
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": dim_ny=", dim_ny
        end if

        if(allocated(var_id_obs_nc)) deallocate(var_id_obs_nc)
        allocate(var_id_obs_nc(dim_ny, dim_nx))
        !allocate(var_id_obs_nc(dim_obs))

        call check(nf90_inq_varid(ncid, var_id_name, var_id_varid))
        call check(nf90_get_var(ncid, var_id_varid, var_id_obs_nc))
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": var_id_obs_nc=", var_id_obs_nc
        end if
    end if

    ! Damping factors
    ! ---------------
    ! Input of flexible damping factors (could be different for each
    ! update step)
    has_damping_state = nf90_inq_varid(ncid, damp_state_name, damp_state_varid)

    if(has_damping_state == nf90_noerr) then

      is_dampfac_state_time_dependent = 1

      if(allocated(dampfac_state_time_dependent_in)) deallocate(dampfac_state_time_dependent_in)
      allocate(dampfac_state_time_dependent_in(1))

      call check(nf90_get_var(ncid, damp_state_varid, dampfac_state_time_dependent_in))
      if (screen > 2) then
        print *, "TSMP-PDAF mype(w)=", mype_world, ": dampfac_state_time_dependent_in=", dampfac_state_time_dependent_in(1)
      end if

    end if

    has_damping_param = nf90_inq_varid(ncid, damp_param_name, damp_param_varid)

    if(has_damping_param == nf90_noerr) then

      is_dampfac_param_time_dependent = 1

      if(allocated(dampfac_param_time_dependent_in)) deallocate(dampfac_param_time_dependent_in)
      allocate(dampfac_param_time_dependent_in(1))

      call check(nf90_get_var(ncid, damp_param_varid, dampfac_param_time_dependent_in))
      if (screen > 2) then
        print *, "TSMP-PDAF mype(w)=", mype_world, ": dampfac_param_time_dependent_in=", dampfac_param_time_dependent_in(1)
      end if

    end if

#ifndef CLMSA
#ifndef OBS_ONLY_CLM
    ! ParFlow observations
    ! --------------------
    has_obs_pf = nf90_inq_varid(ncid, pres_name, pres_varid)

    if(has_obs_pf == nf90_noerr) then

        if(allocated(pressure_obs)) deallocate(pressure_obs)
        allocate(pressure_obs(dim_obs))

        call check(nf90_get_var(ncid, pres_varid, pressure_obs))
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": pressure_obs=", pressure_obs
        end if

        if(allocated(idx_obs_nc))   deallocate(idx_obs_nc)
        allocate(idx_obs_nc(dim_obs))

        call check( nf90_inq_varid(ncid, idx_name, idx_varid) )
        status =  nf90_get_var(ncid, idx_varid, idx_obs_nc)
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": idx_obs_nc=", idx_obs_nc
            print *, "TSMP-PDAF mype(w)=", mype_world, ": status=", status
            print *, "TSMP-PDAF mype(w)=", mype_world, ": nf90_strerror(status)=", nf90_strerror(status)
        end if

        !check, if observation errors are present in observation file
        haserr = nf90_inq_varid(ncid, presserr_name, presserr_varid)
        if(haserr == nf90_noerr) then
            multierr = 1
            !hcp pressure_obserr must be reallocated because dim_obs is not necessary
            !the same for every obs file.
            if(allocated(pressure_obserr)) deallocate(pressure_obserr)
            allocate(pressure_obserr(dim_obs))
            !hcp fin
            call check(nf90_get_var(ncid, presserr_varid, pressure_obserr))
            if (screen > 2) then
                print *, "TSMP-PDAF mype(w)=", mype_world, ": pressure_obserr=", pressure_obserr
            end if
        endif

        !has_depth = nf90_inq_varid(ncid, depth_name, depth_varid)
        !if(has_depth == nf90_noerr) then
        !    crns_flag = 1
        !end if

        ! Read the surface pressure and idxerature data from the file.
        ! Since we know the contents of the file we know that the data
        ! arrays in this program are the correct size to hold all the data.

        if(allocated(x_idx_obs_nc)) deallocate(x_idx_obs_nc)
        allocate(x_idx_obs_nc(dim_obs))

        call check( nf90_inq_varid(ncid, X_IDX_NAME, x_idx_varid) )
        call check( nf90_get_var(ncid, x_idx_varid, x_idx_obs_nc) )
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": x_idx_obs_nc=", x_idx_obs_nc
        end if

        if(allocated(y_idx_obs_nc)) deallocate(y_idx_obs_nc)
        allocate(y_idx_obs_nc(dim_obs))

        call check( nf90_inq_varid(ncid, Y_IDX_NAME, y_idx_varid) )
        call check( nf90_get_var(ncid, y_idx_varid, y_idx_obs_nc) )
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": y_idx_obs_nc=", y_idx_obs_nc
        end if

        if(allocated(z_idx_obs_nc)) deallocate(z_idx_obs_nc)
        allocate(z_idx_obs_nc(dim_obs))

        call check( nf90_inq_varid(ncid, Z_IDX_NAME, z_idx_varid) )
        call check( nf90_get_var(ncid, z_idx_varid, z_idx_obs_nc) )
        !hcp
        if  (crns_flag == 1) then
            z_idx_obs_nc(:)=1
            !if ((maxval(z_idx_obs_nc).NE.1) .OR. (minval(z_idx_obs_nc).NE.1)) then
            !   write(*,*) 'For crns average mode parflow obs layer iz must be 1'
            !   stop
            !endif
        endif
        !end hcp
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": z_idx_obs_nc=", z_idx_obs_nc
        end if

        ! Read observation distances to input observation grid point
        if (obs_interp_switch == 1) then
            if(allocated(x_idx_interp_d_obs_nc)) deallocate(x_idx_interp_d_obs_nc)
            allocate(x_idx_interp_d_obs_nc(dim_obs))

            call check( nf90_inq_varid(ncid, X_IDX_INTERP_D_NAME, x_idx_interp_d_varid) )
            call check( nf90_get_var(ncid, x_idx_interp_d_varid, x_idx_interp_d_obs_nc) )
            if (screen > 2) then
                print *, "TSMP-PDAF mype(w)=", mype_world, ": x_idx_interp_d_obs_nc=", x_idx_interp_d_obs_nc
            end if

            if(allocated(y_idx_interp_d_obs_nc)) deallocate(y_idx_interp_d_obs_nc)
            allocate(y_idx_interp_d_obs_nc(dim_obs))

            call check( nf90_inq_varid(ncid, Y_IDX_INTERP_D_NAME, y_idx_interp_d_varid) )
            call check( nf90_get_var(ncid, y_idx_interp_d_varid, y_idx_interp_d_obs_nc) )
            if (screen > 2) then
                print *, "TSMP-PDAF mype(w)=", mype_world, ": y_idx_interp_d_obs_nc=", y_idx_interp_d_obs_nc
            end if
        end if

    end if
#endif
#endif

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
    ! CLM observations
    ! ----------------
    has_obs_clm = nf90_inq_varid(ncid, obs_name, clmobs_varid)
    if(has_obs_clm == nf90_noerr) then

        if(allocated(clm_obs))   deallocate(clm_obs)
        allocate(clm_obs(dim_obs))

        call check(nf90_get_var(ncid, clmobs_varid, clm_obs))
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": clm_obs=", clm_obs
        end if

        !check, if observation errors are present in observation file
        haserr = nf90_inq_varid(ncid, obserr_name, clmobserr_varid)
        if(haserr == nf90_noerr) then
            multierr = 1
            if(allocated(clm_obserr)) deallocate(clm_obserr)
            allocate(clm_obserr(dim_obs))
            call check(nf90_get_var(ncid, clmobserr_varid, clm_obserr))
            if (screen > 2) then
                print *, "TSMP-PDAF mype(w)=", mype_world, ": clm_obserr=", clm_obserr
            end if
        endif

        ! Read the longitude latidute data from the file.

        if(allocated(clmobs_lon))   deallocate(clmobs_lon)
        allocate(clmobs_lon(dim_obs))

        call check( nf90_inq_varid(ncid, lon_name, clmobs_lon_varid) )
        call check( nf90_get_var(ncid, clmobs_lon_varid, clmobs_lon) )
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": clmobs_lon=", clmobs_lon
        end if

        if(allocated(clmobs_lat))   deallocate(clmobs_lat)
        allocate(clmobs_lat(dim_obs))

        call check( nf90_inq_varid(ncid, lat_name, clmobs_lat_varid) )
        call check( nf90_get_var(ncid, clmobs_lat_varid, clmobs_lat) )
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": clmobs_lat=", clmobs_lat
        end if

        if(allocated(clmobs_layer))   deallocate(clmobs_layer)
        allocate(clmobs_layer(dim_obs))

        haserr = nf90_inq_varid(ncid, layer_name, clmobs_layer_varid)
        if(haserr == nf90_noerr) then
            call check( nf90_get_var(ncid, clmobs_layer_varid, clmobs_layer) )
        else
            ! Default layer is 1
            clmobs_layer(:)=1   !hcp for LST DA
        end if
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": clmobs_layer=", clmobs_layer
        end if

        if(allocated(clmobs_dr))   deallocate(clmobs_dr)
        allocate(clmobs_dr(2))

        call check( nf90_inq_varid(ncid, dr_name, dr_varid) )
        call check( nf90_get_var(ncid, dr_varid, clmobs_dr) )
        if (screen > 2) then
            print *, "TSMP-PDAF mype(w)=", mype_world, ": clmobs_dr=", clmobs_dr
        end if

    end if
#endif
#endif

    call check( nf90_close(ncid) )
    if (screen > 2) then
        print *, "TSMP-PDAF mype(w)=", mype_world, "*** SUCCESS reading observation file ", current_observation_filename, "! "
    end if

  end subroutine read_obs_nc
  !mp end

  !> @author Wolfgang Kurtz, Guowei He, Mukund Pondkule
  !> @date 03.03.2023
  !> @brief Read observation index arrays for C-code
  !> @param[out] no_obs Number of observations
  !> @details
  !> This subroutine reads the observation index arrays for usage in
  !> the enkf_parflow.c for groundwater masking.
  !>
  !> Only used, when ParFlow is one of the component models.
  !>
  !> Index is for ParFlow-type observations
  !>
  !> Only used in `enkf_parflow.c` with `pf_gwmasking=2`.
  !>
  !> Outputs:
  !> --------
  !> Number of observations in `no_obs`.
  !>
  !> Index arrays that are set from NetCDF observation file:
  !> - `tidx_obs`
  !> - `xidx_obs`
  !> - `yidx_obs`
  !> - `zidx_obs`
  !> - `ind_obs`
  subroutine get_obsindex_currentobsfile(no_obs) bind(c,name='get_obsindex_currentobsfile')
    USE mod_tsmp, ONLY: tcycle
    USE mod_assimilation, only: obs_filename
    use netcdf, only: nf90_max_name
    use netcdf, only: nf90_open
    use netcdf, only: nf90_nowrite
    use netcdf, only: nf90_inq_dimid
    use netcdf, only: nf90_inquire_dimension
    use netcdf, only: nf90_inq_varid
    use netcdf, only: nf90_get_var
    use netcdf, only: nf90_close

    implicit none
    integer, intent(out) :: no_obs
    character (len = 110) :: filename
    integer :: ncid, varid
    character (len = *), parameter :: dim_name   = "dim_obs"
    character (len = *), parameter :: idx_name   = "idx"
    character (len = *), parameter :: x_idx_name = "ix"
    character (len = *), parameter :: y_idx_name = "iy"
    character (len = *), parameter :: z_idx_name = "iz"
    character (len = *), parameter :: gwind_name = "gw_indicator"
    character(len = nf90_max_name) :: RecordDimName
    integer :: dimid, status
    integer :: haserr

    write(filename, '(a, i5.5)') trim(obs_filename)//'.', tcycle

    if(allocated(idx_obs_pf))   deallocate(idx_obs_pf)
    if(allocated(x_idx_obs_pf)) deallocate(x_idx_obs_pf)
    if(allocated(y_idx_obs_pf)) deallocate(y_idx_obs_pf)
    if(allocated(z_idx_obs_pf)) deallocate(z_idx_obs_pf)
    if(allocated(ind_obs_pf))   deallocate(ind_obs_pf)

    call check( nf90_open(filename, nf90_nowrite, ncid) )

    call check(nf90_inq_dimid(ncid, dim_name, dimid))
    call check(nf90_inquire_dimension(ncid, dimid, recorddimname, no_obs))

    allocate(idx_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, idx_name, varid) )
    status =  nf90_get_var(ncid, varid, idx_obs_pf)

    allocate(x_idx_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, X_IDX_NAME, varid) )
    call check( nf90_get_var(ncid, varid, x_idx_obs_pf) )

    allocate(y_idx_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, Y_IDX_NAME, varid) )
    call check( nf90_get_var(ncid, varid, y_idx_obs_pf) )

    allocate(z_idx_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, Z_IDX_NAME, varid) )
    call check( nf90_get_var(ncid, varid, z_idx_obs_pf) )

    allocate(ind_obs_pf(no_obs))
    call check( nf90_inq_varid(ncid, gwind_name, varid) )
    call check( nf90_get_var(ncid, varid, ind_obs_pf) )

    ptr_tidx_obs = c_loc(idx_obs_pf)
    ptr_xidx_obs = c_loc(x_idx_obs_pf)
    ptr_yidx_obs = c_loc(y_idx_obs_pf)
    ptr_zidx_obs = c_loc(z_idx_obs_pf)
    ptr_ind_obs  = c_loc(ind_obs_pf)

    call check( nf90_close(ncid) )

  end subroutine get_obsindex_currentobsfile

  !> @author Wolfgang Kurtz, Guowei He, Mukund Pondkule
  !> @date 03.03.2023
  !> @brief Deallocation of observation arrays
  !> @details
  !> This subroutine deallocates the observation arrays used in
  !> subroutine `read_obs_nc`.
  subroutine clean_obs_nc()

    USE mod_assimilation, ONLY: filtertype

    implicit none
   ! if(allocated(idx_obs_nc))deallocate(idx_obs_nc)
    if(allocated(pressure_obs))deallocate(pressure_obs)
    !if(allocated(pressure_obserr))deallocate(pressure_obserr)
    !if(allocated(x_idx_obs_nc))deallocate(x_idx_obs_nc)
    !if(allocated(y_idx_obs_nc))deallocate(y_idx_obs_nc)
    !if(allocated(z_idx_obs_nc))deallocate(z_idx_obs_nc)
    !kuw: clean clm observations
    IF (.NOT. filtertype == 5 .AND. .NOT. filtertype == 7 .AND. .NOT. filtertype == 8) THEN
      ! For LETKF, LESTKF, LEnKF lat/lon are used
      if(allocated(clmobs_lon))deallocate(clmobs_lon)
      if(allocated(clmobs_lat))deallocate(clmobs_lat)
    END IF
    if(allocated(clm_obs))deallocate(clm_obs)
    if(allocated(clmobs_layer))deallocate(clmobs_layer)
    if(allocated(clmobs_dr))deallocate(clmobs_dr)
    !if(allocated(clm_obserr))deallocate(clm_obserr)
    !kuw end
  end subroutine clean_obs_nc

  !> @author Wolfgang Kurtz, Guowei He, Mukund Pondkule
  !> @date 03.03.2023
  !> @brief Deallocation of observation index arrays
  !> @details
  !> This subroutine deallocates the observation index arrays used in
  !> subroutine `get_obsindex_currentobsfile`.
  !>
  !> Only used in `enkf_parflow.c` with `pf_gwmasking=2`.
  subroutine clean_obs_pf() bind(c,name='clean_obs_pf')
    implicit none
    if(allocated(idx_obs_pf))deallocate(idx_obs_pf)
    if(allocated(x_idx_obs_pf))deallocate(x_idx_obs_pf)
    if(allocated(y_idx_obs_pf))deallocate(y_idx_obs_pf)
    if(allocated(z_idx_obs_pf))deallocate(z_idx_obs_pf)
    if(allocated(ind_obs_pf))   deallocate(ind_obs_pf)
  end subroutine clean_obs_pf

  !> @author Wolfgang Kurtz, Guowei He
  !> @date 21.03.2022
  !> @brief Return number of observations from file
  !> @param[in] fn Filename of the observation file
  !> @param[out] nn number of observations in `fn`
  !> @details
  !>     Reads the content of the variable (!) named `no_obs` from
  !>     NetCDF file `fn`.
  !>
  !>     Uses  subroutines from the NetCDF module.
  !>
  !>     The result is returned in `nn`.
  !>
  !>     The result is used to decide if the next observation file is
  !>     used or not.
  subroutine check_n_observationfile(fn,nn)
      use netcdf, only: nf90_max_name, nf90_open, nf90_nowrite, &
          nf90_inq_varid, nf90_get_var, nf90_close

      implicit none

      character(len=*),intent(in) :: fn
      integer, intent(out)        :: nn

      integer :: ncid, varid, status !,dimid
      character (len = *), parameter :: varname = "no_obs"

      !character (len = *), parameter :: dim_name = "dim_obs"
      !character(len = nf90_max_name) :: recorddimname

      call check(nf90_open(fn, nf90_nowrite, ncid))
      !call check(nf90_inq_dimid(ncid, dim_name, dimid))
      !call check(nf90_inquire_dimension(ncid, dimid, recorddimname, nn))
      call check( nf90_inq_varid(ncid, varname, varid) )
      call check( nf90_get_var(ncid, varid, nn) )
      call check(nf90_close(ncid))

  end subroutine check_n_observationfile

  !> @author Yorck Ewerdwalbesloh, Johannes Keller
  !> @date 11.09.2023
  !> @brief Return data assimilation interval from file
  !> @param[in] fn Filename of the observation file
  !> @param[out] aa new da_interval (number of time steps until next assimilation time step)
  !> @details
  !>     Reads the content of the variable name `da_interval` from NetCDF
  !>     file `fn` using subroutines from the NetCDF module.
  !>     The result is returned in `aa`.
  !>
  !>     The result is used to adapt the da_interval until the next observation file.
  !>
  !>     Adapted for TSMP-PDAF by Johannes Keller, 28.05.2025.
  subroutine check_n_observationfile_da_interval(fn,aa)
    use shr_kind_mod, only: r8 => shr_kind_r8
    use netcdf, only: nf90_max_name, nf90_open, nf90_nowrite, &
      nf90_inq_varid, nf90_get_var, nf90_close, nf90_noerr

    implicit none

    character(len=*),intent(in) :: fn
    real, intent(out)        :: aa


    integer :: ncid, varid, status !,dimid
    character (len = *), parameter :: varname = "da_interval"
    real(r8) :: dtime ! land model time step (sec)

    !character (len = *), parameter :: dim_name = "dim_obs"
    !character(len = nf90_max_name) :: recorddimname

    call check(nf90_open(fn, nf90_nowrite, ncid))
    !call check(nf90_inq_dimid(ncid, dim_name, dimid))
    !call check(nf90_inquire_dimension(ncid, dimid, recorddimname, nn))

    call check( nf90_inq_varid(ncid, varname, varid))
    call check( nf90_get_var(ncid, varid, aa) )
    call check(nf90_close(ncid))

  end subroutine check_n_observationfile_da_interval

  !> @author Wolfgang Kurtz, Guowei He, Mukund Pondkule
  !> @date 03.03.2023
  !> @brief Error handling for netCDF commands
  !> @param[in] status netCDF command status
  !> @details
  !> This subroutine checks the status of a netCDF command and prints
  !> an error message if necessary.
  subroutine check(status)

    use netcdf, only: nf90_noerr
    use netcdf, only: nf90_strerror
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped from NetCDF error handling subroutine check"
    end if
  end subroutine check


  ! !> @author Yorck Ewerdwalbesloh
  !     !> @date 11.09.2023
  !     !> @brief Return data assimilation interval from file
  !     !> @param[in] fn Filename of the observation file
  !     !> @param[out] nn number of hours until next assimilation time step
  !     !> @details
  !     !>     Reads the content of the variable name `da_interval_variable` from NetCDF
  !     !>     file `fn` using subroutines from the NetCDF module.
  !     !>     The result is returned in `nn`.
  !     !>
  !     !>     The result is used to decide if the next observation file is
  !     !>     used or not.
  !   subroutine check_n_observationfile_da_interval(fn,nn)
  !       use shr_kind_mod, only: r8 => shr_kind_r8
  !       use netcdf, only: nf90_max_name, nf90_open, nf90_nowrite, &
  !           nf90_inq_varid, nf90_get_var, nf90_close, nf90_noerr
  !       use clm_varcon, only: ispval
  !       use clm_time_manager, only : get_step_size

  !       implicit none

  !       character(len=*),intent(in) :: fn
  !       real, intent(out)        :: nn


  !       integer :: ncid, varid, status !,dimid
  !       character (len = *), parameter :: varname = "da_interval_variable"
  !       real(r8) :: dtime ! land model time step (sec)

  !       !character (len = *), parameter :: dim_name = "dim_obs"
  !       !character(len = nf90_max_name) :: recorddimname

  !       dtime = get_step_size()

  !       call check(nf90_open(fn, nf90_nowrite, ncid))
  !       !call check(nf90_inq_dimid(ncid, dim_name, dimid))
  !       !call check(nf90_inquire_dimension(ncid, dimid, recorddimname, nn))
  !       status = nf90_inq_varid(ncid, varname, varid)
  !       if (status == nf90_noerr) then
  !           call check(nf90_inq_varid(ncid, varname, varid))
  !           call check( nf90_get_var(ncid, varid, nn) )
  !           call check(nf90_close(ncid))
  !           ! at this point: half hourly time steps, this is adjusted here. In the GRACE files, da_interval is set up as hours
  !           ! --> is adjusted using information from inside CLM
  !           nn = nn*INT(3600/dtime)
  !       else
  !           nn = ispval
  !       end if

  !   end subroutine check_n_observationfile_da_interval



    !> @author Yorck Ewerdwalbesloh
    !> @date 04.12.2023
    !> @brief Return set zero interval for running mean of model variables from file
    !> @param[in] fn Filename of the observation file
    !> @param[out] nn number of hours until setting zero
    !> @details
    !>     Reads the content of the variable name `set_zero` from NetCDF
    !>     file `fn` using subroutines from the NetCDF module.
    !>     The result is returned in `nn`.
    !>
    !>     The result is used to reset the running average of state variables.
    subroutine check_n_observationfile_set_zero(fn,nn)
        use shr_kind_mod, only: r8 => shr_kind_r8
        use netcdf, only: nf90_max_name, nf90_open, nf90_nowrite, &
            nf90_inq_varid, nf90_get_var, nf90_close, nf90_noerr
        use clm_varcon, only: ispval
        use clm_time_manager, only : get_step_size

        implicit none

        character(len=*),intent(in) :: fn
        integer, intent(out)        :: nn

        integer :: ncid, varid, status !,dimid
        character (len = *), parameter :: varname = "set_zero"
        real(r8) :: dtime ! land model time step (sec)

        !character (len = *), parameter :: dim_name = "dim_obs"
        !character(len = nf90_max_name) :: recorddimname

        dtime = get_step_size()

        call check(nf90_open(fn, nf90_nowrite, ncid))
        !call check(nf90_inq_dimid(ncid, dim_name, dimid))
        !call check(nf90_inquire_dimension(ncid, dimid, recorddimname, nn))
        status = nf90_inq_varid(ncid, varname, varid)
        if (status == nf90_noerr) then
            call check(nf90_inq_varid(ncid, varname, varid))
            call check( nf90_get_var(ncid, varid, nn) )
            call check(nf90_close(ncid))
            ! at this point: half hourly time steps, this is adjusted here. In the GRACE files, set_zero is set up as hours
            ! --> is adjusted using information from inside CLM
            if (nn/=ispval) then
                nn = nn*INT(3600/dtime)
            end if
        else
            nn = ispval
        end if

    end subroutine check_n_observationfile_set_zero


    subroutine check_n_observationfile_next_type(fn, obs_type_str)
        use netcdf, only: nf90_max_name
        use netcdf, only: nf90_open
        use netcdf, only: nf90_nowrite
        use netcdf, only: nf90_inq_dimid
        use netcdf, only: nf90_inquire_dimension
        use netcdf, only: nf90_inq_varid
        use netcdf, only: nf90_noerr
        use netcdf, only: nf90_get_var
        use netcdf, only: nf90_close
        use mod_assimilation, only: screen
        implicit none

        character(len=*), intent(in)  :: fn
        character(len=*), intent(out) :: obs_type_str

        integer :: ncid, status, obstype_varid, dimid
        integer :: dim_obs
        character (len = *), parameter :: dim_name = "dim_obs"
        character (len = *), parameter :: type_name   = "type_clm"
        character(len=20), allocatable :: obs_type_lok(:)
        character(len = nf90_max_name) :: RecordDimName


        call check( nf90_open(fn, nf90_nowrite, ncid) )
        call check(nf90_inq_dimid(ncid, dim_name, dimid))
        call check(nf90_inquire_dimension(ncid, dimid, recorddimname, dim_obs))

        if(allocated(obs_type_lok))   deallocate(obs_type_lok)
        allocate(obs_type_lok(dim_obs))

        obs_type_str = ''

        status = nf90_inq_varid(ncid, "type_clm", obstype_varid)
        if (status == nf90_noerr) then
            call check(nf90_inq_varid(ncid, "type_clm", obstype_varid))
            call check(nf90_get_var(ncid, obstype_varid, obs_type_lok))
            obs_type_str = trim(obs_type_lok(1))
        end if
        call check(nf90_close(ncid))

        if(allocated(obs_type_lok))   deallocate(obs_type_lok)

    end subroutine check_n_observationfile_next_type


    subroutine update_obs_type(obs_type_str)
        use enkf_clm_mod, only: clmupdate_tws, clmupdate_swc, clmupdate_T, clmupdate_texture
        use mod_parallel_pdaf, only: abort_parallel
        implicit none

        character(len=*), intent(in) :: obs_type_str

        select case (trim(adjustl(obs_type_str)))
        case ('GRACE')
            clmupdate_tws     = 1
            clmupdate_swc     = 0
            clmupdate_T       = 0
            clmupdate_texture = 0

        case ('SM')
            clmupdate_tws     = 0
            clmupdate_swc     = 1
            clmupdate_T       = 0
            clmupdate_texture = 0

        case default
            write(*,*) 'ERROR: Unknown obs_type_str in update_obs_type:', trim(obs_type_str)
            call abort_parallel()
        end select
    end subroutine update_obs_type




    subroutine domain_def_clm(lon_clmobs, lat_clmobs, dim_obs, &
        longxy, latixy, longxy_obs, latixy_obs)

        use spmdMod,   only : npes, iam
        use domainMod, only : ldomain, lon1d, lat1d
        use decompMod, only : get_proc_total, get_proc_bounds, ldecomp
        use GridcellType, only: grc
        use shr_kind_mod, only: r8 => shr_kind_r8
        use enkf_clm_mod, only: hactiveg_levels, num_hactiveg
        !USE mod_parallel_pdaf, &
        !   ONLY: mpi_2integer, mpi_minloc
        USE mod_parallel_pdaf, &
            ONLY: comm_filter, npes_filter, abort_parallel, &
            mpi_integer, mpi_double_precision, mpi_in_place, mpi_sum, &
            mype_world, mpi_2integer, mpi_minloc, mype_filter
        real, intent(in) :: lon_clmobs(:)
        real, intent(in) :: lat_clmobs(:)
        integer, intent(in) :: dim_obs
        integer, allocatable, intent(inout) :: longxy(:)
        integer, allocatable, intent(inout) :: latixy(:)
        integer, allocatable, intent(inout) :: longxy_obs(:)
        integer, allocatable, intent(inout) :: latixy_obs(:)
        integer :: ni, nj, ii, jj, kk, cid, ier, ncells, nlunits, &
        ncols, npatches, ncohorts, counter, i, g, ll
        real :: minlon, minlat, maxlon, maxlat
        real(r8), pointer :: lon(:)
        real(r8), pointer :: lat(:)

        real(r8), allocatable :: longxy_obs_lokal(:), latixy_obs_lokal(:)

        INTEGER, allocatable :: in_mpi_(:,:), out_mpi_(:,:)

        integer :: begg, endg   ! per-proc gridcell ending gridcell indices

        real(r8) :: lat1, lon1, lat2, lon2, a, c, R, pi

        real(r8) :: dist
        real(r8), allocatable :: min_dist(:)
        integer, allocatable :: min_g(:)

        integer :: ierror

        integer :: lok_lon, lok_lat



        lon   => grc%londeg
        lat   => grc%latdeg


        ni = ldomain%ni
        nj = ldomain%nj

        ! get total number of gridcells, landunits,
        ! columns, patches and cohorts on processor

        call get_proc_total(iam, ncells, nlunits, ncols, npatches, ncohorts)

        ! beg and end gridcell
        call get_proc_bounds(begg=begg, endg=endg)

        if (allocated(longxy)) deallocate(longxy)
        if (allocated(latixy)) deallocate(latixy)
        allocate(longxy(num_hactiveg), stat=ier)
        allocate(latixy(num_hactiveg), stat=ier)


        longxy(:) = 0
        latixy(:) = 0


        counter = 1
        do ii = 1, nj
            do jj = 1, ni
                cid = (ii-1)*ni + jj
                do ll = 1, num_hactiveg
                    kk = hactiveg_levels(ll,1)
                    if(cid == ldecomp%gdc2glo(kk)) then
                        latixy(counter) = ii
                        longxy(counter) = jj
                        counter = counter + 1
                    end if
                end do
            end do
        end do

        if (allocated(min_dist)) deallocate(min_dist)
        allocate(min_dist(dim_obs))
        min_dist(:) = huge(1.0d0)

        if (allocated(min_g)) deallocate(min_g)
        allocate(min_g(dim_obs))

        R = 6371.0
        pi = 3.14159265358979323846
        do i = 1, dim_obs
            do g = begg, endg

                ! check distance from each grid point to observation location --> take the coordinate in local system that equals
                ! the one of the closest coordinate
                lat1 = lat(g) * pi / 180.0
                lon1 = lon(g) * pi / 180.0
                lat2 = lat_clmobs(i) * pi / 180.0
                lon2 = lon_clmobs(i) * pi / 180.0

                a = sin((lat2 - lat1) / 2)**2 + cos(lat1) * cos(lat2) * sin((lon2 - lon1) / 2)**2
                c = 2 * atan2(sqrt(a), sqrt(1 - a))
                dist = R * c

                if (dist < min_dist(i)) then
                    min_dist(i) = dist
                    min_g(i) = g
                end if
            end do
        end do


        IF (ALLOCATED(in_mpi_)) DEALLOCATE(in_mpi_)
        ALLOCATE(in_mpi_(2,dim_obs))
        IF (ALLOCATED(out_mpi_)) DEALLOCATE(out_mpi_)
        ALLOCATE(out_mpi_(2,dim_obs))

        in_mpi_(1,:) = int(ceiling(min_dist))
        in_mpi_(2,:) = min_g

        if (allocated(longxy_obs_lokal)) deallocate(longxy_obs_lokal)
        if (allocated(latixy_obs_lokal)) deallocate(latixy_obs_lokal)

        allocate(longxy_obs_lokal(dim_obs))
        allocate(latixy_obs_lokal(dim_obs))

        do i =1, dim_obs
            outer: do ii = 1, nj
                do jj = 1, ni
                    cid = (ii-1)*ni + jj
                    do kk = begg, endg
                        if (kk == in_mpi_(2,i)) then
                            if(cid == ldecomp%gdc2glo(kk)) then
                                if (min_dist(i)<30) then
                                    latixy_obs_lokal(i) = ii
                                    longxy_obs_lokal(i) = jj
                                else
                                    longxy_obs_lokal(i) = -9999
                                    latixy_obs_lokal(i) = -9999
                                end if
                            exit outer
                            end if
                        end if
                    end do
                end do
            end do outer
        end do


        if (allocated(longxy_obs)) deallocate(longxy_obs)
        if (allocated(latixy_obs)) deallocate(latixy_obs)
        allocate(longxy_obs(dim_obs), stat=ier)
        allocate(latixy_obs(dim_obs), stat=ier)

        in_mpi_(2,:) = longxy_obs_lokal
        call mpi_allreduce(in_mpi_,out_mpi_, dim_obs, mpi_2integer, mpi_minloc, comm_filter, ierror)
        longxy_obs(:) = out_mpi_(2,:)

        in_mpi_(2,:) = latixy_obs_lokal
        call mpi_allreduce(in_mpi_,out_mpi_, dim_obs, mpi_2integer, mpi_minloc, comm_filter, ierror)
        latixy_obs(:) = out_mpi_(2,:)

        deallocate(longxy_obs_lokal)
        deallocate(latixy_obs_lokal)
        deallocate(in_mpi_)
        deallocate(out_mpi_)
        deallocate(min_dist)
        deallocate(min_g)


        if (mype_filter == 0) then
            print*, "longxy_obs = ", longxy_obs
            print*, "latixy_obs = ", latixy_obs
        end if

    end subroutine domain_def_clm


end module mod_read_obs

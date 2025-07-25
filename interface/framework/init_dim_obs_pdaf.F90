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
!init_dim_obs_pdaf.F90: TSMP-PDAF implementation of routine
!                       'init_dim_obs_pdaf' (PDAF online coupling)
!-------------------------------------------------------------------------------------------

!$Id: init_dim_obs_pdaf.F90 1441 2013-10-04 10:33:42Z lnerger $
!BOP
!
! !ROUTINE: init_dim_obs_pdaf --- Compute number of observations
!
! !INTERFACE:
SUBROUTINE init_dim_obs_pdaf(step, dim_obs_p)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/EnKF/ETKF/ESTKF
!
! The routine is called at the beginning of each
! analysis step.  It has to initialize the size of
! the observation vector according to the current
! time step for the PE-local domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
!   USE mod_assimilation, &
!        ONLY : nx, ny, local_dims, obs_p, obs_index_p, &
!        coords_obs, local_dims_obs
!   USE mod_parallel_pdaf, &
!        ONLY: mype_filter, npes_filter, COMM_filter, MPI_INTEGER, &
!        MPIerr, MPIstatus
   USE mod_parallel_pdaf, &
      ONLY: mype_filter, comm_filter, npes_filter, abort_parallel, &
      mpi_integer, mpi_double_precision, mpi_in_place, mpi_sum, &
      mype_world, mpi_2integer, mpi_maxloc
   USE mod_assimilation, &
      ONLY: obs_p, obs_index_p, dim_obs, obs_filename, &
      obs, &
      obs_interp_indices_p, &
      obs_interp_weights_p, &
      pressure_obserr_p, clm_obserr_p, &
      local_dims_obs, &
      ! dim_obs_p, &
      obs_id_p, &
      obscov, obscov_inv, filtertype, &
      tws_temp_mean_d, &
      temp_mean_filename, &
#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
!hcp 
!CLMSA needs the physical  coordinates of the elements of state vector 
!and observation array.        
      longxy, latixy, longxy_obs, latixy_obs, &
      longxy_obs_floor, latixy_obs_floor, &
!hcp end
#endif
#endif
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
      depth_obs_p, &
      sc_p, idx_obs_nc_p, &
#endif
#endif
      var_id_obs, maxlon, minlon, maxlat, &
      minlat, maxix, minix, maxiy, miniy, lon_var_id, ix_var_id, lat_var_id, iy_var_id, &
      screen
   Use mod_read_obs, &
      only: idx_obs_nc, pressure_obs, pressure_obserr, multierr, &
      read_obs_nc, clean_obs_nc, x_idx_obs_nc, y_idx_obs_nc, &
      z_idx_obs_nc, &
      x_idx_interp_d_obs_nc, y_idx_interp_d_obs_nc, &
      clm_obs, &
      var_id_obs_nc, dim_nx, dim_ny, &
      clmobs_lon, clmobs_lat, clmobs_layer, clmobs_dr, clm_obserr, &
      crns_flag, depth_obs, dampfac_state_time_dependent_in, dampfac_param_time_dependent_in, &
      clm_obscov, vec_numPoints_global, vec_useObs_global, &
      vec_useObs, lon_temp_mean, lat_temp_mean, tws_temp_mean, read_temp_mean_model
   use mod_tsmp, &
      only: idx_map_subvec2state_fortran, tag_model_parflow, enkf_subvecsize, &
      nx_glob, ny_glob, nz_glob, &
#ifndef CLMSA
#ifndef OBS_ONLY_CLM
      xcoord, ycoord, zcoord, xcoord_fortran, ycoord_fortran, &
      zcoord_fortran, &
#endif
#endif
      tag_model_clm, point_obs, obs_interp_switch, is_dampfac_state_time_dependent, &
      dampfac_state_time_dependent, is_dampfac_param_time_dependent, dampfac_param_time_dependent, model

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
!kuw
   use shr_kind_mod, only: r8 => shr_kind_r8
#ifdef CLMFIVE
   use GridcellType, only: grc
   use ColumnType         , only : col
   use clm_varcon, only: spval
#else  
   USE clmtype,                  ONLY : clm3
#endif  
   use decompMod , only : get_proc_bounds, get_proc_global
   
!kuw end
!hcp
!use the subroutine written by Mukund "domain_def_clm" to evaluate longxy,
!latixy, longxy_obs, latixy_obs
   USE enkf_clm_mod, only: get_interp_idx, clmupdate_tws, num_layer, hactiveg_levels
!hcp end
#endif
#endif

   USE, INTRINSIC :: iso_c_binding

   IMPLICIT NONE
   ! !ARGUMENTS:
   INTEGER, INTENT(in)  :: step       ! Current time step
   INTEGER, INTENT(out) :: dim_obs_p  ! Dimension of observation vector
   ! !CALLING SEQUENCE:
   ! Called by: PDAF_seik_analysis, PDAF_seik_analysis_newT    (as U_init_dim_obs)
   ! Called by: PDAF_enkf_analysis_rlm, PDAF_enkf_analysis_rsm
   ! Called by: PDAF_etkf_analysis, PDAF_etkf_analysis_T
   ! Called by: PDAF_estkf_analysis, PDAF_estkf_analysis_fixed
   !EOP

   ! *** Local variables
   integer :: ierror
   INTEGER :: max_var_id
   INTEGER :: sum_dim_obs_p
   INTEGER :: i,j,k,count_points, countC, countR, countObs  ! Counters
   INTEGER :: count_interp ! Counter for interpolation grid cells
   INTEGER :: m,l          ! Counters
   INTEGER :: idx         ! Computed Index
   logical :: is_multi_observation_files
   character (len = 110) :: current_observation_filename
   integer,allocatable :: local_dis(:),local_dim(:)
   integer :: k_count,nsc !hcp
   real    :: sum_interp_weights
   integer :: c

#ifndef PARFLOW_STAND_ALONE
#ifndef OBS_ONLY_PARFLOW
   real(r8), pointer :: lon(:)
   real(r8), pointer :: lat(:)
   ! pft: "plant functional type"
   integer :: begp, endp   ! per-proc beginning and ending pft indices
   integer :: begc, endc   ! per-proc beginning and ending column indices
   integer :: begl, endl   ! per-proc beginning and ending landunit indices
   integer :: begg, endg   ! per-proc gridcell ending gridcell indices
   integer :: numg         ! total number of gridcells across all processors
   integer :: numl         ! total number of landunits across all processors
   integer :: numc         ! total number of columns across all processors
   integer :: nump         ! total number of pfts across all processors
   real    :: deltax, deltay
   !real    :: deltaxy, y1 , x1, z1, x2, y2, z2, R, dist, deltaxy_max
   integer :: numPoints ! minimum number of points so that the GRACE observation is used
   integer, allocatable :: vec_numPoints(:) ! number of model grid cells that are in a radius of dr around the GRACE observation
   INTEGER, allocatable :: in_mpi(:,:), out_mpi(:,:)
   INTEGER, ALLOCATABLE :: ipiv(:)
   real, ALLOCATABLE :: work(:)
   logical :: is_use_dr
#endif
#endif

   ! ****************************************
   ! *** Initialize observation dimension ***
   ! ****************************************

   ! Read observation file
   ! ---------------------

   ! Default: no local damping factors
   is_dampfac_state_time_dependent = 0
   is_dampfac_param_time_dependent = 0

   !  if I'm root in filter, read the nc file
   is_multi_observation_files = .true.
   if (is_multi_observation_files) then
         ! Set name of current NetCDF observation file
         write(current_observation_filename, '(a, i5.5)') trim(obs_filename)//'.', step
   else
         ! Single NetCDF observation file (currently NOT used)
         write(current_observation_filename, '(a, i5.5)') trim(obs_filename)
   end if

   if (mype_filter .eq. 0) then
         ! Read current NetCDF observation file
         call read_obs_nc(current_observation_filename)
   end if

   ! Broadcast first variables
   ! -------------------------
   ! Dimension of observation vector
   call mpi_bcast(dim_obs, 1, MPI_INTEGER, 0, comm_filter, ierror)
   ! Switch for vector of observation errors
   call mpi_bcast(multierr, 1, MPI_INTEGER, 0, comm_filter, ierror)
   ! broadcast crns_flag
   call mpi_bcast(crns_flag, 1, MPI_INTEGER, 0, comm_filter, ierror)
   ! broadcast dim_ny and dim_nx
   if(point_obs.eq.0) then
      call mpi_bcast(dim_nx, 1, MPI_INTEGER, 0, comm_filter, ierror)
      call mpi_bcast(dim_ny, 1, MPI_INTEGER, 0, comm_filter, ierror)
   endif
   ! broadcast damping factor flags
   call mpi_bcast(is_dampfac_state_time_dependent, 1, MPI_INTEGER, 0, comm_filter, ierror)
   call mpi_bcast(is_dampfac_param_time_dependent, 1, MPI_INTEGER, 0, comm_filter, ierror)

   ! broadcast dampfac_state_time_dependent_in
   if(is_dampfac_state_time_dependent.eq.1) then

      if (mype_filter .ne. 0) then ! for all non-master proc
         if(allocated(dampfac_state_time_dependent_in)) deallocate(dampfac_state_time_dependent_in)
         allocate(dampfac_state_time_dependent_in(1))
      end if

      if (screen > 2) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": Before setting dampfac_state_time_dependent"
      end if

      call mpi_bcast(dampfac_state_time_dependent_in, 1, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      if (screen > 2) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dampfac_state_time_dependent_in=", dampfac_state_time_dependent_in
      end if

      ! Set C-version of dampfac_state_time_dependent with value read from obsfile
      dampfac_state_time_dependent = dampfac_state_time_dependent_in(1)

      if (screen > 2) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dampfac_state_time_dependent=", dampfac_state_time_dependent
      end if

   end if

   ! broadcast dampfac_param_time_dependent_in
   if(is_dampfac_param_time_dependent.eq.1) then

      if (mype_filter .ne. 0) then ! for all non-master proc
         if(allocated(dampfac_param_time_dependent_in)) deallocate(dampfac_param_time_dependent_in)
         allocate(dampfac_param_time_dependent_in(1))
      end if

      if (screen > 2) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": Before setting dampfac_param_time_dependent"
      end if

      call mpi_bcast(dampfac_param_time_dependent_in, 1, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
      if (screen > 2) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dampfac_param_time_dependent_in=", dampfac_param_time_dependent_in
      end if

      ! Set C-version of dampfac_param_time_dependent with value read from obsfile
      dampfac_param_time_dependent = dampfac_param_time_dependent_in(1)

      if (screen > 2) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dampfac_param_time_dependent=", dampfac_param_time_dependent
      end if

   end if

! Allocate observation arrays for non-root procs
! ----------------------------------------------
   if (mype_filter .ne. 0) then ! for all non-master proc

      if(allocated(clm_obs)) deallocate(clm_obs)
      allocate(clm_obs(dim_obs))
      if(allocated(clmobs_lon)) deallocate(clmobs_lon)
      allocate(clmobs_lon(dim_obs))
      if(allocated(clmobs_lat)) deallocate(clmobs_lat)
      allocate(clmobs_lat(dim_obs))
      if(allocated(clmobs_dr)) deallocate(clmobs_dr)
      allocate(clmobs_dr(2))
      if(allocated(clmobs_layer)) deallocate(clmobs_layer)
      allocate(clmobs_layer(dim_obs))
      if(point_obs.eq.0) then
            if(allocated(var_id_obs_nc)) deallocate(var_id_obs_nc)
            allocate(var_id_obs_nc(dim_ny, dim_nx))
      endif
      if(multierr.eq.1) then 
            if(allocated(clm_obserr)) deallocate(clm_obserr)
            allocate(clm_obserr(dim_obs))
      end if
      if(multierr.eq.2) then
            if(allocated(clm_obscov)) deallocate(clm_obscov)
            allocate(clm_obscov(dim_obs, dim_obs))
      end if


   end if


   call mpi_bcast(clm_obs, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
   if(multierr.eq.1) call mpi_bcast(clm_obserr, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
   if(multierr.eq.2) call mpi_bcast(clm_obscov, dim_obs*dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
   call mpi_bcast(clmobs_lon, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
   call mpi_bcast(clmobs_lat, dim_obs, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
   call mpi_bcast(clmobs_dr,  2, MPI_DOUBLE_PRECISION, 0, comm_filter, ierror)
   call mpi_bcast(clmobs_layer, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)
   if(point_obs.eq.0) call mpi_bcast(var_id_obs_nc, dim_obs, MPI_INTEGER, 0, comm_filter, ierror)


! CLM grid information
! --------------------
! Results used only in `localize_covar_pdaf` for LEnKF
! Calling could be restricted to LEnKF
!hcp
!use the subroutine written by Mukund "domain_def_clm" to evaluate longxy,
!latixy, longxy_obs, latixy_obs
! Index arrays of longitudes and latitudes

   ! if exist CLM-type obs
   if(model .eq. tag_model_clm) then 

      ! Obtain general CLM index information
#ifdef CLMFIVE
      lon   => grc%londeg
      lat   => grc%latdeg
#else      
      lon   => clm3%g%londeg
      lat   => clm3%g%latdeg
#endif
      call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
      call get_proc_global(numg, numl, numc, nump)
   end if
!hcp end

   ! Number of observations in process-local domain
   ! ----------------------------------------------
   ! Additionally `obs_id_p` is set (the NetCDF index of the
   ! observation corresponding to the state index in the local domain)
   dim_obs_p = 0



   ! Switch for how to check index of CLM observations
   ! True: Use snapping distance between long/lat on CLM grid
   ! False: Use index arrays from `domain_def_clm`
   is_use_dr = .true.

   if(model .eq. tag_model_clm) then


      if (ALLOCATED(vec_useObs_global).eqv..false.) then
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
         ! count grdicells that are in a certain radius
         do i = 1, dim_obs
            count_points = 0
            ! only take gridcells into account that have at least one hydrological active column
            do c = 1, num_layer(1)
               j = hactiveg_levels(c,1)
               if (lon(j)<180)then
                  deltax = abs(lon(j)-clmobs_lon(i))
               else
                  deltax = abs((lon(j)-360)-clmobs_lon(i))
               end if
               deltay = abs(lat(j)-clmobs_lat(i))
               if((deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2))) then
                  count_points = count_points+1
               end if
            end do           
            vec_numPoints(i) = count_points
         end do
         ! get vec_numPoints from all processes and add them up together via mpi_allreduce
         call mpi_allreduce(vec_numPoints,vec_numPoints_global, dim_obs, mpi_integer, mpi_sum, comm_filter, ierror)
         ! only observations should be used that "see" enough gridcells
         ! numPoints equation from Anne
         numPoints = int(ceiling((clmobs_dr(1)*2/0.11 * clmobs_dr(2)*2/0.11)/2.0)) 
         if (screen > 2) then
            if (mype_filter==0) then
            print *, "Minimum number of points for using one observation is ", numPoints
            end if
         end if
         vec_useObs_global = merge(vec_useObs_global,.false.,vec_numPoints_global.ge.numPoints)
         vec_useObs = vec_useObs_global
   

         ! The observation should be reproduced by the process that has the most grid points inside the observation radius
         ! However, in the observation operator, all points will be used by using mpi
         ! obs_id_p could be of good use for this as it assigns each grid point in the local process a observation (or none if
         ! there isn't any in the radius)

         ! The corresponding process for each observation is found out using the following code

         in_mpi(1,:) = vec_numPoints
         in_mpi(2,:) = mype_filter
         call mpi_allreduce(in_mpi,out_mpi, dim_obs, mpi_2integer, mpi_maxloc, comm_filter, ierror)

         vec_useObs = merge(vec_useObs,.false.,out_mpi(2,:).eq.mype_filter)

         IF (ALLOCATED(in_mpi)) DEALLOCATE(in_mpi)
         IF (ALLOCATED(out_mpi)) DEALLOCATE(out_mpi)

         ! vec_useObs_global(1) = .false.
         ! vec_useObs_global(2) = .true.

         dim_obs_p = count(vec_useObs)
         ! if (mype_filter == 0) then
         !    dim_obs_p = COUNT(vec_useObs_global)
         !    vec_useObs = vec_useObs_global
         ! else
         !    dim_obs_p = 0
         !    vec_useObs(:) = .false.
         ! end if

         if(allocated(obs_id_p)) deallocate(obs_id_p)
         allocate(obs_id_p(begg:endg))
         obs_id_p(:) = 0

         !IF (ALLOCATED(obs_index_p)) DEALLOCATE(obs_index_p)
         !ALLOCATE(obs_index_p(sum(vec_numPoints)))
         !count_points = 1
         do i = 1, dim_obs
            if (vec_useObs_global(i)) then
               do c = 1, num_layer(1)
                  j = hactiveg_levels(c,1)
                  if (lon(j)<180)then
                     deltax = abs(lon(j)-clmobs_lon(i))
                  else
                     deltax = abs((lon(j)-360)-clmobs_lon(i))
                  end if
                  deltay = abs(lat(j)-clmobs_lat(i))
                  if((deltax.le.clmobs_dr(1)).and.(deltay.le.clmobs_dr(2))) then
                     obs_id_p(j) = i
                     !obs_index_p(count_points) = j
                     !count_points = count_points+1
                  end if
               end do
            end if
         end do   

      else

         dim_obs_p = count(vec_useObs)
            ! if (mype_filter == 0) then
            !    dim_obs_p = COUNT(vec_useObs_global)
            ! else
            !    dim_obs_p = 0
            ! end if

      end if

      dim_obs = count(vec_useObs_global)
   end if

   if (screen > 2) then
      if (mype_filter==0) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_f_pdaf: vec_useObs_global=", vec_useObs_global
         print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_f_pdaf: vec_numPoints_global=", vec_numPoints_global
      end if
   end if

   if (screen > 2) then
         print *, "TSMP-PDAF mype(w)=", mype_world, ": init_dim_obs_pdaf: dim_obs_p=", dim_obs_p
   end if

   ! add and broadcast size of local observation dimensions using mpi_allreduce 
   call mpi_allreduce(dim_obs_p, sum_dim_obs_p, 1, MPI_INTEGER, MPI_SUM, &
         comm_filter, ierror) 




   IF (ALLOCATED(obs)) DEALLOCATE(obs)
   ALLOCATE(obs(dim_obs))
   IF (ALLOCATED(obs_p)) DEALLOCATE(obs_p)
   ALLOCATE(obs_p(dim_obs_p))
   obs = pack(clm_obs,vec_useObs_global)
   obs_p = pack(clm_obs,vec_useObs)
   


   if(multierr.eq.1) then
      if (allocated(clm_obserr_p)) deallocate(clm_obserr_p)
      allocate(clm_obserr_p(dim_obs_p))
   endif
   if (multierr.eq.2) then
      IF (ALLOCATED(obscov)) DEALLOCATE(obscov)
      ALLOCATE(obscov(COUNT(vec_useObs_global),COUNT(vec_useObs_global)))
      countR = 1
      countC = 1
      do i = 1, size(clm_obscov,1)
         if(vec_useObs_global(i)) then
            do j = 1, size(clm_obscov,2)
               if (vec_useObs_global(j)) then
                  obscov(countR,countC) = clm_obscov(i,j)
                  countC = countC + 1
               end if
            end do
            countC = 1
            countR = countR + 1
         end if
      end do 

      if (filtertype.ne.2) then
         print *, 'Compute inverse of observation covariance matrix'
         IF (ALLOCATED(obscov_inv)) DEALLOCATE(obscov_inv)
         ALLOCATE(obscov_inv(dim_obs,dim_obs))
         ALLOCATE(ipiv(dim_obs))
         ALLOCATE(work(dim_obs))
         obscov_inv = obscov
         !LU factorization
         call dgetrf(dim_obs,dim_obs,obscov_inv,dim_obs,ipiv,ierror)
         ! Inverse using LU factorization
         call dgetri(dim_obs, obscov_inv, dim_obs, ipiv, work, dim_obs, ierror)
         if (ierror /= 0) then
               stop 'init_dim_obs_pdaf:  inversion failed!'
         end if
         IF (ALLOCATED(ipiv)) DEALLOCATE(ipiv)
         IF (ALLOCATED(work)) DEALLOCATE(work)
      end if

   end if    

   count_points = 1
   countObs = 1 


   !  clean up the temp data from nc file
   ! ------------------------------------
   call clean_obs_nc()
   
   ! Read temporal mean TWS from model for observation operator, only for GRACE data assimilation
   if (clmupdate_tws.eq.1) then
      ! do it only in the first call of this routine

      if (.not. allocated(tws_temp_mean_d)) then

         ! fill tws_temp mean, lat_temp_mean and lon_temp_mean
         call read_temp_mean_model(temp_mean_filename)

         if (allocated(tws_temp_mean_d)) DEALLOCATE(tws_temp_mean_d)
         ALLOCATE(tws_temp_mean_d(begg:endg))
         tws_temp_mean_d(:) = spval

         !this process only need the sub domain information
         do j = begg,endg
            ! find lon and lat in the file that corresponds to that of the grid point of the sub process
            outer3: do l = 1,size(lon_temp_mean,1)
               do k=1,size(lon_temp_mean,2)
                  if (lon_temp_mean(l,k).eq.lon(j) .and. lat_temp_mean(l,k).eq.lat(j)) then
                     tws_temp_mean_d(j) = tws_temp_mean(l,k)
                     exit outer3
                  end if
               end do
            end do outer3

            if (lon(j).ne.lon_temp_mean(l,k) .or. lat(j).ne.lat_temp_mean(l,k)) then
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
   end if

END SUBROUTINE init_dim_obs_pdaf
   
   
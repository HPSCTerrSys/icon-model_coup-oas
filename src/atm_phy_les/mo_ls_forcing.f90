! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! This module initializes and applies the large-scale forcing for idealized simulations
! This module assumes that the model grid is FLAT
! 2013-JUNE-04: AT THIS STAGE LS FORCING WILL WORK IN RESTART MODE ONLY IF ITS CALLED EVERY
! DYNAMIC TIMESTEP. TO MAKE IT WORK "SMOOTHLY" ADD_VAR HAS TO WORK ON 1D VARS

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_ls_forcing

  USE mo_kind,                ONLY: wp, vp
  USE mo_io_units,            ONLY: find_next_free_unit
  USE mo_exception,           ONLY: message, finish, message_text
  USE mo_impl_constants,      ONLY: success, max_char_length
  USE mo_nonhydro_types,      ONLY: t_nh_prog, t_nh_diag, t_nh_metrics
  USE mo_model_domain,        ONLY: t_patch
  USE mo_parallel_config,     ONLY: nproma
  USE mo_loopindices,         ONLY: get_indices_c
  USE mo_nh_vert_interp_les,  ONLY: vert_intp_linear_1d, &
    &                               vertical_derivative
  USE mo_ls_forcing_nml,      ONLY: is_subsidence_moment, is_subsidence_heat, is_ls_forcing, &
    &                               is_advection, is_advection_uv,is_advection_tq,           &
    &                               is_nudging, is_nudging_uv, is_nudging_t, is_nudging_q,   &
    &                               is_geowind, is_rad_forcing, dt_relax_uv, dt_relax_t,     &
    &                               dt_relax_q, is_ls_coriolis, theta_nudging, &
    &                               nudge_start_height_uv, nudge_full_height_uv,  &
    &                               nudge_start_height_t, nudge_full_height_t,  &
    &                               nudge_start_height_q, nudge_full_height_q, is_sim_rad
  USE mo_fortran_tools,       ONLY: init
  USE mo_scm_nml,             ONLY: i_scm_netcdf, lscm_ls_forcing_ini, lscm_icon_ini, &
    &                               lat_scm, lon_scm
  USE mo_read_interface,      ONLY: nf
  USE mo_netcdf
  USE mo_mpi,                 ONLY: get_my_global_mpi_id
  USE mo_statistics,          ONLY: levels_horizontal_mean
  USE mo_dynamics_config,     ONLY: lcoriolis
  USE mo_math_constants,      ONLY: deg2rad
  USE mo_grid_config,         ONLY: grid_angular_velocity
  USE mo_fortran_tools,       ONLY: assert_acc_host_only
  USE mo_physical_constants,  ONLY: grav,rd_o_cpd,vtmpc1
  USE mo_run_config,          ONLY: iqv, iqc, iqi

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_ls_forcing, apply_ls_forcing, apply_ls_forc_nudge_qv, apply_ls_forc_nudge_uvt, apply_ls_sfc_forcing

  !Anurag Dipankar, MPIM (2013-May): variables for large-scale (LS)
  !forcing to be used in periodic domain or other relveant cases.

  !Dependent on z and time only
  REAL(wp), SAVE, ALLOCATABLE ::  &
    w_ls(:,:),            & !subsidence          [m/s]
    u_geo  (:,:),         & !u-geostrophic wind  [m/s]
    v_geo  (:,:),         & !v-geostrophic wind  [m/s]
    ddt_temp_hadv_ls(:,:),& !LS horizontal advective tendency for temp [K/s]
    ddt_temp_rad(:,:),    & !radiative tendency for temp [k/s]
    ddt_qv_hadv_ls(:,:),  & !LS horizontal advective tendency for qv [1/s]
    ddt_u_hadv_ls(:,:),   & !LS horizontal advective tendency for u [m/s^2]
    ddt_v_hadv_ls(:,:),   & !LS horizontal advective tendency for m [m/s^2]
    u_nudg(:,:),          & !LS nudging u    [m/s]
    v_nudg(:,:),          & !LS nudging v    [m/s]
    theta_nudg(:,:),      & !LS nudging temp [K]
    qv_nudg(:,:),         & !LS nudging qv   [kg/kg]
    bnd_sfc_lat_flx(:),   & !surface Latent heat flux[W/m2]
    bnd_sfc_sens_flx(:),  & !surface sensible heat flux[W/m2]
    bnd_ts(:),            & !Surface temperature[K]
    bnd_tg(:),            & !Skin temperature[K]
    bnd_qvs(:),           & !Surface specific vapor humidity[kg/kg]
    bnd_Ch(:),            & !Drag coefficient for heat
    bnd_Cm(:),            & !Drag coefficient for momentum
    bnd_Cq(:),            & !Drag coefficient for moisture
    bnd_ustar(:),         & !Friction velocity[m/s]
    bnd_lat(:),           & !Latitude (degrees)
    bnd_lon(:),           & !Longitude (degrees)
    tempf(:,:),           &
    tempf_f(:,:),         &
    tempf_sf(:)

  ! To make sure that LS forcing is not applied before it is read in
  REAL(wp), SAVE :: dt_forcing = 0._wp
  REAL(wp), SAVE :: dt_nudging = 0._wp

  INTEGER,  SAVE :: nt      ! number of time steps in init_SCM.nc


  CONTAINS


  !>
  !! init_ls_forcing
  !!------------------------------------------------------------------------
  !! Initialize large-scale forcings from input file
  !!------------------------------------------------------------------------
  !! Initial release by Anurag Dipankar, MPI-M (2013-May-30)
  !! Martin Koehler, DWD (2022-May-16): switch fluxes after input to downward-positive
  !!                                    convention of ICON (attention LES uses opposite)

  SUBROUTINE init_ls_forcing(p_metrics)

    TYPE(t_nh_metrics),INTENT(in), TARGET :: p_metrics

    CHARACTER(len=*), PARAMETER :: routine = 'mo_ls_forcing:init_ls_forcing'

    REAL(wp), ALLOCATABLE, DIMENSION(:) :: zz, zw, zu, zv, z_dt_temp_adv, &
                                         & z_dt_temp_rad, z_dt_qv_adv,z_dt_u_adv,z_dt_v_adv
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw_nc, zu_nc, zv_nc, z_dt_temp_adv_nc, &
                                         & z_dt_temp_rad_nc, z_dt_qv_adv_nc,ztheta_nc,z_qv_nc, &
                                         & z_dt_u_adv_nc,z_dt_v_adv_nc, zz_nc
    REAL(wp), ALLOCATABLE, DIMENSION(:) :: time_nc

    REAL(wp), ALLOCATABLE, DIMENSION(:) :: ztheta, z_qv
    REAL(wp)  :: end_time
    INTEGER   :: iunit, ist, nk, jk, nlev, i, nskip, n, t0
    INTEGER   :: varid                  !< id of variable in netcdf file
    INTEGER   :: fileid                 !< id number of netcdf file
    INTEGER   :: dimid                  !< id number of dimension
    INTEGER   :: nf_status, nf_status2  !< return status of netcdf function

    CHARACTER(len=max_char_length) :: time_unit, time_unit_short
    REAL(wp)                       :: time_unit_in_sec
    REAL(wp)                       :: temp_nf(1)

    nlev   = SIZE(p_metrics%z_mc,2)

!------------------------------------------------------------------------------
! READ LS FORCING FILE
!------------------------------------------------------------------------------

    IF (is_ls_forcing) THEN

      IF (i_scm_netcdf > 0) THEN

        IF (i_scm_netcdf == 2) THEN

!------------------------------------------------------------------------------
! Forcing with DEPHY unified SCM NETCDF file

          CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
            & TRIM(routine)//'   File init_SCM.nc cannot be opened')

          CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), TRIM(routine)//' lev')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nk), routine)

          CALL nf (nf90_inq_dimid(fileid, 't0', dimid), routine)
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = t0), routine)

          CALL nf (nf90_inq_dimid(fileid, 'time', dimid), TRIM(routine)//' nt')    ! dimension 'time'
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt), routine)

          WRITE(message_text,'(a,i6,a,i6)') 'SCM LS forcing input file: nk, ', nk, ', nt ', nt
          CALL message (routine,message_text)

          ALLOCATE( zz_nc(nk,nt), time_nc(nt), zw_nc(nk,nt), zu_nc(nk,nt), zv_nc(nk,nt), z_dt_temp_adv_nc(nk,nt), &
                  & z_dt_temp_rad_nc(nk,nt), z_dt_qv_adv_nc(nk,nt), z_dt_u_adv_nc(nk,nt), z_dt_v_adv_nc(nk,nt) )
          zz_nc=0.0_wp ; time_nc=0.0_wp ; zw_nc=0.0_wp ; zu_nc=0.0_wp ; zv_nc=0.0_wp ; z_dt_temp_adv_nc=0.0_wp
          z_dt_temp_rad_nc=0.0_wp ; z_dt_qv_adv_nc=0.0_wp ; z_dt_u_adv_nc=0.0_wp ; z_dt_v_adv_nc=0.0_wp

          ALLOCATE( w_ls(nlev,nt), u_geo(nlev,nt), v_geo(nlev,nt), ddt_temp_hadv_ls(nlev,nt), &
                  & ddt_temp_rad(nlev,nt), ddt_qv_hadv_ls(nlev,nt),                           &
                  & ddt_u_hadv_ls(nlev,nt), ddt_v_hadv_ls(nlev,nt),                           &
                  & tempf(nk,t0), tempf_f(nk,nt), tempf_sf(nt) )
          w_ls=0.0_wp ; u_geo=0.0_wp ; v_geo=0.0_wp ; ddt_temp_hadv_ls=0.0_wp
          ddt_temp_rad=0.0_wp ; ddt_qv_hadv_ls=0.0_wp ; ddt_u_hadv_ls=0.0_wp ; ddt_v_hadv_ls=0.0_wp

          ALLOCATE( bnd_sfc_lat_flx(nt), bnd_sfc_sens_flx(nt), bnd_ts(nt), bnd_qvs(nt), &
                  & bnd_Ch(nt), bnd_Cm(nt), bnd_Cq(nt), bnd_ustar(nt), bnd_tg(nt),  &
                  & bnd_lat(nt), bnd_lon(nt))
          bnd_sfc_lat_flx=0.0_wp ; bnd_sfc_sens_flx=0.0_wp ; bnd_ts=0.0_wp ; bnd_qvs=0.0_wp
          bnd_Ch=0.0_wp ; bnd_Cm=0.0_wp ; bnd_Cq=0.0_wp ; bnd_ustar=0.0_wp ; bnd_tg=0.0_wp
          bnd_lat=0.0_wp; bnd_lon=0.0_wp

          CALL nf (nf90_inq_varid (fileid, 'zh_forc', varid ), routine)
          CALL nf (nf90_get_var   (fileid, varid, tempf_f)   , routine)
          zz_nc=tempf_f(nk:1:-1,:)

          !Check if the file is written in descending order
          IF(zz_nc(1,1) < zz_nc(nk,1)) CALL finish ( routine, 'Write LS forcing data in descending order!')

          CALL nf(nf90_inq_varid(fileid, 'time', varid), TRIM(routine)//' time')
          CALL nf(nf90_get_var(fileid, varid, time_nc), routine)
          time_unit = ''               !necessary for formatting
          CALL nf(nf90_get_att(fileid, varid, 'units', time_unit), TRIM(routine)//' units')

          time_unit_short = time_unit(1:INDEX(time_unit,"since")-2)  ! e.g. "minutes since 2020-1-1 00:00:00"
          IF ( get_my_global_mpi_id() == 0 ) THEN
            WRITE(*,*) 'time unit in init_SCM.nc, long: ', TRIM(time_unit), '   short: ', TRIM(time_unit_short), &
              & '     times: ', time_nc
          END IF

          SELECT CASE (time_unit_short)
          CASE ('s', 'sec', 'second', 'seconds')
            time_unit_in_sec = 1
          CASE ('min', 'minute', 'minutes')
            time_unit_in_sec = 60
          CASE ('h', 'hour', 'hours')
            time_unit_in_sec = 3600
          CASE ('d', 'day', 'days')
            time_unit_in_sec = 84600
          CASE DEFAULT
            CALL finish ( routine, 'time unit in init_SCM.nc not s/min/h!')
          END SELECT

          time_nc = time_nc * time_unit_in_sec  ! conversion to [s] from min/hours/days

          IF (nt<=1) THEN
            dt_forcing = -999._wp
          ELSE
            dt_forcing = time_nc(2) - time_nc(1)
            DO i=2,size(time_nc)-1
              IF ((time_nc(i+1)-time_nc(i)) /= dt_forcing) THEN
                CALL finish(routine,'input timesteps of equal length needed')
              END IF
            END DO
          ENDIF

          nf_status = nf90_inq_varid (fileid, 'wa', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'w not available in init_SCM.nc.  It will be set to 0.')
            zw_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid,tempf_f)
            zw_nc=tempf_f(nk:1:-1,:)
          END IF
!         write(*,*) 'wLS:',zw_nc(:,1)

          nf_status = nf90_inq_varid (fileid, 'ug', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'ug not available in init_SCM.nc.  It will be set to 0.')
            zu_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_f)
            zu_nc=tempf_f(nk:1:-1,:)
          END IF
!         write(*,*) 'ug:',zu_nc(:,1)

          nf_status = nf90_inq_varid (fileid, 'vg', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'vg not available in init_SCM.nc.  It will be set to 0.')
            zv_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_f)
            zv_nc=tempf_f(nk:1:-1,:)
          END IF
!         write(*,*) 'vg:',zv_nc(:,1)

          nf_status = nf90_inq_varid (fileid, 'tnta_adv', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_temp_adv_nc not available in init_SCM.nc.  It will be set to 0.')
            z_dt_temp_adv_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_f)
            z_dt_temp_adv_nc=tempf_f(nk:1:-1,:)
          END IF
!         write(*,*) 'dTadv:',z_dt_temp_adv_nc(:,1)

          nf_status  = nf90_inq_varid (fileid, 'temp_rad', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_temp_rad_nc not available in init_SCM.nc.  It will be set to 0.')
            z_dt_temp_rad_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_f)
            z_dt_temp_rad_nc=tempf_f(nk:1:-1,:)
          END IF
!         write(*,*) 'dTrad:',z_dt_temp_rad_nc(:,1)

          nf_status  = nf90_inq_varid     (fileid, 'tnqv_adv',varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_qv_adv_nc not available in init_SCM.nc It will be set to 0.')
            z_dt_qv_adv_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid,varid, tempf_f)
            z_dt_qv_adv_nc=tempf_f(nk:1:-1,:)
          END IF
!         write(*,*) 'dQVadv:',z_dt_qv_adv_nc(:,1)

          nf_status = nf90_inq_varid (fileid, 'tnua_adv', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_u_adv_nc not available in init_SCM.nc.  It will be set to 0.')
            z_dt_u_adv_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_f)
            z_dt_u_adv_nc=tempf_f(nk:1:-1,:)
          END IF
!         write(*,*) 'dUadv:',z_dt_u_adv_nc(:,1)

          nf_status = nf90_inq_varid (fileid, 'tnva_adv', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'z_dt_v_adv_nc not available in init_SCM.nc.  It will be set to 0.')
            z_dt_v_adv_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_f)
            z_dt_v_adv_nc=tempf_f(nk:1:-1,:)
          END IF
!         write(*,*) 'dVadv:',z_dt_v_adv_nc(:,1)

          nf_status = nf90_inq_varid (fileid, 'hfls', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'bnd_sfc_lat_flx_nc not available in init_SCM.nc.  It will be set to 0.')
            bnd_sfc_lat_flx=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_sfc_lat_flx=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_sfc_lat_flx:',bnd_sfc_lat_flx(:)

          nf_status = nf90_inq_varid (fileid, 'hfss', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'bnd_sfc_sens_flx not available in init_SCM.nc.  It will be set to 0.')
            bnd_sfc_sens_flx=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_sfc_sens_flx=tempf_sf(:)
          END IF
!          write(*,*) 'bnd_sfc_sens_flx:',bnd_sfc_sens_flx(:)

          nf_status = nf90_inq_varid (fileid, 'ts_forc', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'bnd_ts not available in init_SCM.nc.  It will be set to 400.')
            bnd_ts=400.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_ts=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_ts:',bnd_ts(:)

          nf_status = nf90_inq_varid (fileid, 'tskin', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'bnd_tg not available in init_SCM.nc.  It will be set to 400.')
            bnd_tg=400.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_tg=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_tg:',bnd_tg(:)

          ! Use either ts_forc or tskin as inputs for the variables bnd_tg and bnd_ts
          IF ( bnd_ts(1) == 400.0_wp )  bnd_ts = bnd_tg
          IF ( bnd_tg(1) == 400.0_wp )  bnd_tg = bnd_ts
!         write(*,*) 'new bnd_ts and bnd_tg:',bnd_ts(:),bnd_tg(:)

          nf_status  = nf90_inq_varid     (fileid, 'qvs', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'qvs not available in init_SCM.nc.  It will be set to 0.')
            bnd_qvs=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_qvs=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_qvs:',bnd_qvs(:)

          nf_status  = nf90_inq_varid (fileid, 'Ch', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'Ch not available in init_SCM.nc.  It will be set to 0.')
            bnd_Ch=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_Ch=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_Ch:',bnd_Ch(:)

          nf_status  = nf90_inq_varid (fileid, 'Cm', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'Cm not available in init_SCM.nc.  It will be set to 0.')
            bnd_Cm=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_Cm=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_Cm:',bnd_Cm(:)

          nf_status = nf90_inq_varid (fileid, 'Cq', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'Cq not available in init_SCM.nc.  It will be set to 0.')
            bnd_Cq=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_Cq=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_Cq:',bnd_Cq(:)

          nf_status = nf90_inq_varid (fileid, 'ustar', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'ustar not available in init_SCM.nc.  It will be set to 0.')
            bnd_ustar=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_ustar=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_ustar:',bnd_ustar(:)

          nf_status = nf90_inq_varid (fileid, 'lat', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'bnd_lat not available in init_SCM.nc.  It will be set to 0.')
            bnd_lat=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_lat=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_lat:',bnd_lat(:)

          nf_status = nf90_inq_varid (fileid, 'lon', varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'bnd_lon not available in init_SCM.nc.  It will be set to 0.')
            bnd_lon=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf_sf)
            bnd_lon=tempf_sf(:)
          END IF
!         write(*,*) 'bnd_lon:',bnd_lon(:)

          CALL nf (nf90_close(fileid), routine)

          ! switch fluxes from DEPHY format to ICON convention of downward positive
          bnd_sfc_sens_flx = - bnd_sfc_sens_flx
          bnd_sfc_lat_flx  = - bnd_sfc_lat_flx

          DO n = 1 , nt
            !Now perform interpolation to grid levels assuming:
            !a) linear interpolation
            !b) Beyond the last Z level the values are linearly extrapolated
            !c) Assuming model grid is flat-NOT on sphere

            CALL vert_intp_linear_1d(zz_nc(:,n),zw_nc           (:,n),p_metrics%z_mc(1,:,1),w_ls            (:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),zu_nc           (:,n),p_metrics%z_mc(1,:,1),u_geo           (:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),zv_nc           (:,n),p_metrics%z_mc(1,:,1),v_geo           (:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_temp_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_temp_hadv_ls(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_temp_rad_nc(:,n),p_metrics%z_mc(1,:,1),ddt_temp_rad    (:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_qv_adv_nc  (:,n),p_metrics%z_mc(1,:,1),ddt_qv_hadv_ls  (:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_u_adv_nc   (:,n),p_metrics%z_mc(1,:,1),ddt_u_hadv_ls   (:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_v_adv_nc   (:,n),p_metrics%z_mc(1,:,1),ddt_v_hadv_ls   (:,n))
          END DO

          DEALLOCATE( zz_nc, time_nc, zw_nc, zu_nc, zv_nc, z_dt_temp_adv_nc, z_dt_temp_rad_nc, &
                   &  z_dt_qv_adv_nc, z_dt_u_adv_nc, z_dt_v_adv_nc, tempf, tempf_f, tempf_sf )


!------------------------------------------------------------------------------
! Forcing with normal NETCDF file (including real ICON NWP forcing)

        ELSE

          CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
            & TRIM(routine)//'   File init_SCM.nc cannot be opened')

          CALL nf (nf90_inq_dimid(fileid, 'lev' , dimid), TRIM(routine)//' lev')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nk), routine)

          CALL nf (nf90_inq_dimid(fileid, 'nt', dimid), TRIM(routine)//' nt')    ! dimension 'nt'
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt), routine)

          WRITE(message_text,'(a,i6,a,i6)') 'SCM LS forcing input file: nk, ', nk, ', nt ', nt
          CALL message (routine,message_text)

          ALLOCATE( zz_nc(nk,nt), time_nc(nt), zw_nc(nk,nt), zu_nc(nk,nt), zv_nc(nk,nt), z_dt_temp_adv_nc(nk,nt), &
                    z_dt_temp_rad_nc(nk,nt), z_dt_qv_adv_nc(nk,nt), z_dt_u_adv_nc(nk,nt), z_dt_v_adv_nc(nk,nt) )
          zz_nc=0.0_wp;time_nc=0.0_wp;zw_nc=0.0_wp;zu_nc=0.0_wp;zv_nc=0.0_wp;z_dt_temp_adv_nc=0.0_wp
          z_dt_temp_rad_nc=0.0_wp;z_dt_qv_adv_nc=0.0_wp;z_dt_u_adv_nc=0.0_wp;z_dt_v_adv_nc=0.0_wp

          ALLOCATE( w_ls(nlev,nt), u_geo(nlev,nt), v_geo(nlev,nt), ddt_temp_hadv_ls(nlev,nt), &
                    ddt_temp_rad(nlev,nt), ddt_qv_hadv_ls(nlev,nt),&
                    ddt_u_hadv_ls(nlev,nt),ddt_v_hadv_ls(nlev,nt) )
          w_ls=0.0_wp;u_geo=0.0_wp;v_geo=0.0_wp;ddt_temp_hadv_ls=0.0_wp
          ddt_temp_rad=0.0_wp;ddt_qv_hadv_ls=0.0_wp;ddt_u_hadv_ls=0.0_wp;ddt_v_hadv_ls=0.0_wp

          ALLOCATE( bnd_sfc_lat_flx(nt), bnd_sfc_sens_flx(nt), bnd_ts(nt), bnd_qvs(nt),&
                    bnd_Ch(nt), bnd_Cm(nt),bnd_Cq(nt), bnd_ustar(nt),bnd_tg(nt), &
                    bnd_lat(nt),bnd_lon(nt))
          bnd_sfc_lat_flx=0.0_wp; bnd_sfc_sens_flx=0.0_wp; bnd_ts=0.0_wp; bnd_qvs=0.0_wp
          bnd_Ch=0.0_wp; bnd_Cm=0.0_wp;bnd_Cq=0.0_wp; bnd_ustar=0.0_wp;bnd_tg=0.0_wp
          bnd_lat=0.0_wp;bnd_lon=0.0_wp

          CALL nf (nf90_inq_varid(fileid, 'height', varid), TRIM(routine)//' height')
          CALL nf (nf90_get_var(fileid, varid, zz_nc), routine)
         !write(*,*) 'zz_nc',zz_nc
          !Check if the file is written in descending order
          IF(zz_nc(1,1) < zz_nc(nk,1)) CALL finish ( routine, 'Write LS forcing data in descending order!')


          CALL nf(nf90_inq_varid(fileid, 'time', varid), TRIM(routine)//' time')
          CALL nf(nf90_get_var(fileid, varid, time_nc), routine)
          time_unit = ''               !necessary for formatting
          CALL nf(nf90_get_att(fileid, varid, 'units', time_unit), TRIM(routine)//' units')

          time_unit_short = time_unit(1:INDEX(time_unit,"since")-2)  ! e.g. "minutes since 2020-1-1 00:00:00"
          IF ( get_my_global_mpi_id() == 0 ) THEN
            WRITE(*,*) 'time unit in init_SCM.nc, long: ', TRIM(time_unit), '   short: ', TRIM(time_unit_short), &
              & '     times: ', time_nc
          END IF

          SELECT CASE (time_unit_short)
          CASE ('s', 'sec', 'second', 'seconds')
            time_unit_in_sec = 1
          CASE ('min', 'minute', 'minutes')
            time_unit_in_sec = 60
          CASE ('h', 'hour', 'hours')
            time_unit_in_sec = 3600
          CASE ('d', 'day', 'days')
            time_unit_in_sec = 84600
          CASE DEFAULT
            CALL finish ( routine, 'time unit in init_SCM.nc not s/min/h!')
          END SELECT

          time_nc = time_nc * time_unit_in_sec     ! conversion to [s] from min/hours/days

          IF (nt<=1) THEN
            dt_forcing = -999._wp
          ELSE
            dt_forcing = time_nc(2) - time_nc(1)
            DO i=2,size(time_nc)-1
              IF ((time_nc(i+1)-time_nc(i)) /= dt_forcing) THEN
                CALL finish(routine,'input timesteps of equal length needed')
              END IF
            END DO
          ENDIF

          CALL nf (nf90_inq_varid(fileid, 'wLS', varid), TRIM(routine)//' wLS')
          CALL nf (nf90_get_var(fileid, varid,zw_nc), routine)
         !write(*,*) 'wLS',zw_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'uGEO', varid), TRIM(routine)//' uGEO')
          CALL nf (nf90_get_var(fileid, varid,zu_nc), routine)
         !write(*,*) 'uGEO',zu_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'vGEO', varid), TRIM(routine)//' vGEO')
          CALL nf (nf90_get_var(fileid, varid,zv_nc), routine)
         !write(*,*) 'vGEO',zv_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'dTadv', varid), TRIM(routine)//' dTadv')
          CALL nf (nf90_get_var(fileid, varid,z_dt_temp_adv_nc), routine)
         !write(*,*) 'dTadv',z_dt_temp_adv_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'dTrad', varid), TRIM(routine)//' dTrad')
          CALL nf (nf90_get_var(fileid, varid,z_dt_temp_rad_nc), routine)
         !write(*,*) 'dTrad',z_dt_temp_rad_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'dQVadv', varid), TRIM(routine)//' dQVadv')
          CALL nf (nf90_get_var(fileid, varid,z_dt_qv_adv_nc), routine)
         !write(*,*) 'dQVadv',z_dt_qv_adv_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'dUadv', varid), TRIM(routine)//' dUadv')
          CALL nf (nf90_get_var(fileid, varid,z_dt_u_adv_nc), routine)
         !write(*,*) 'dUadv',z_dt_u_adv_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'dVadv', varid), TRIM(routine)//' dVadv')
          CALL nf (nf90_get_var(fileid, varid,z_dt_v_adv_nc), routine)
         !write(*,*) 'dVadv',z_dt_v_adv_nc(:,1)

          CALL nf (nf90_inq_varid(fileid, 'sfc_lat_flx', varid), TRIM(routine)//' sfc_lat_flx')
          CALL nf (nf90_get_var(fileid, varid, bnd_sfc_lat_flx), routine)
         !write(*,*) 'bnd_sfc_lat_flx',bnd_sfc_lat_flx(:)

          CALL nf (nf90_inq_varid(fileid, 'sfc_sens_flx', varid), TRIM(routine)//' sfc_sens_flx')
          CALL nf (nf90_get_var(fileid, varid, bnd_sfc_sens_flx), routine)
         !write(*,*) 'bnd_sfc_sens_flx',bnd_sfc_sens_flx(:)

          CALL nf (nf90_inq_varid(fileid, 'ts', varid), TRIM(routine)//' ts')
          CALL nf (nf90_get_var(fileid, varid, bnd_ts), routine)
         !write(*,*) 'bnd_ts',bnd_ts(:)

          CALL nf (nf90_inq_varid(fileid, 'tg', varid), TRIM(routine)//' tg')
          CALL nf (nf90_get_var(fileid, varid, bnd_tg), routine)
         !write(*,*) 'bnd_tg',bnd_tg(:)

          CALL nf (nf90_inq_varid(fileid, 'qvs', varid), TRIM(routine)//' qvs')
          CALL nf (nf90_get_var(fileid, varid, bnd_qvs), routine)
         !write(*,*) 'bnd_qvs',bnd_qvs(:)

          CALL nf (nf90_inq_varid(fileid, 'Ch', varid), TRIM(routine)//' Ch')
          CALL nf (nf90_get_var(fileid, varid, bnd_Ch), routine)
         !write(*,*) 'bnd_Ch',bnd_Ch(:)

          CALL nf (nf90_inq_varid(fileid, 'Cm', varid), TRIM(routine)//' Cm')
          CALL nf (nf90_get_var(fileid, varid, bnd_Cm), routine)
         !write(*,*) 'bnd_Cm',bnd_Cm(:)

          CALL nf (nf90_inq_varid(fileid, 'Cq', varid), TRIM(routine)//' Cq')
          CALL nf (nf90_get_var(fileid, varid, bnd_Cq), routine)
         !write(*,*) 'bnd_Cq',bnd_Cq(:)

          CALL nf (nf90_inq_varid(fileid, 'ustar', varid), TRIM(routine)//' ustar')
          CALL nf (nf90_get_var(fileid, varid, bnd_ustar), routine)
         !write(*,*) 'bnd_ustar',bnd_ustar(:)

          CALL nf (nf90_inq_varid(fileid, 'latitude', varid), TRIM(routine)//' lat')
          CALL nf (nf90_get_var(fileid, varid, temp_nf), routine)
          bnd_lat(:)=temp_nf(1)
         !write(*,*) 'bnd_lat',bnd_lat(:)

          CALL nf (nf90_inq_varid(fileid, 'longitude', varid), TRIM(routine)//' lon')
          CALL nf (nf90_get_var(fileid, varid, temp_nf), routine)
          bnd_lon(:)=temp_nf(1)
         !write(*,*) 'bnd_lon',bnd_lon(:)

          CALL nf (nf90_close(fileid), routine)

          ! switch fluxes from old input format to ICON convention of downward positive
          IF ( .NOT. lscm_icon_ini ) THEN
            bnd_sfc_sens_flx = - bnd_sfc_sens_flx
            bnd_sfc_lat_flx  = - bnd_sfc_lat_flx
          ENDIF

          DO n = 1 , nt
              !Now perform interpolation to grid levels assuming:
              !a) linear interpolation
              !b) Beyond the last Z level the values are linearly extrapolated
              !c) Assuming model grid is flat-NOT on sphere

              CALL vert_intp_linear_1d(zz_nc(:,n),zw_nc(:,n),p_metrics%z_mc(1,:,1),w_ls(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),zu_nc(:,n),p_metrics%z_mc(1,:,1),u_geo(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),zv_nc(:,n),p_metrics%z_mc(1,:,1),v_geo(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_temp_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_temp_hadv_ls(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_temp_rad_nc(:,n),p_metrics%z_mc(1,:,1),ddt_temp_rad(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_qv_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_qv_hadv_ls(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_u_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_u_hadv_ls(:,n))
              CALL vert_intp_linear_1d(zz_nc(:,n),z_dt_v_adv_nc(:,n),p_metrics%z_mc(1,:,1),ddt_v_hadv_ls(:,n))
          END DO

          DEALLOCATE( zz_nc, time_nc, zw_nc, zu_nc, zv_nc, z_dt_temp_adv_nc, z_dt_temp_rad_nc, &
                   &  z_dt_qv_adv_nc, z_dt_u_adv_nc, z_dt_v_adv_nc)

        ENDIF

!------------------------------------------------------------------------------
! Forcing with ASCII file

      ELSE

        !Open formatted file to read ls forcing data
        iunit = find_next_free_unit(10,20)
        OPEN (unit=iunit,file='ls_forcing.dat',access='SEQUENTIAL', &
              form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)

        IF (ist/=success) CALL finish(routine, 'open ls_forcing.dat failed')

        !Read the input file till end. The order of file assumed is:
        !Z(m) - u_geo(m/s) - v_geo(m/s) - w_ls(m/s) - dt_temp_rad(K/s) - ddt_qv_hadv_ls(1/s) - ddt_temp_hadv_ls(K/s)

        !Skip the first line and read next 2 lines with information about vertical and time levels
        READ(iunit,*,IOSTAT=ist)                       !skip
        IF (ist == success) READ(iunit,*,IOSTAT=ist)nk,dt_forcing,end_time !vertical levels, forcing interval, end of forcing data
        IF(ist/=success) &
          CALL finish(routine, &
            'Must provide vertical level, forcing interval, end of forcing time info in forcing file')

        nt = INT(end_time/dt_forcing)+1

        IF (nt>1) THEN
          READ(iunit,*,IOSTAT=ist)nskip !lines to skip between successive time levels
          IF (ist/=success) &
            CALL finish(routine, &
              'Time levels > 1 so must provide number lines to skip between successive time levels in forcing file')
        END IF

        IF (nt<=1) dt_forcing=-999._wp

        ALLOCATE( zz(nk), zw(nk), zu(nk), zv(nk), z_dt_temp_adv(nk),                        &
                & z_dt_temp_rad(nk), z_dt_qv_adv(nk), z_dt_u_adv(nk), z_dt_v_adv(nk) )

        ALLOCATE( w_ls(nlev,nt), u_geo(nlev,nt), v_geo(nlev,nt), ddt_temp_hadv_ls(nlev,nt), &
                & ddt_temp_rad(nlev,nt), ddt_qv_hadv_ls(nlev,nt),                           &
                & ddt_u_hadv_ls(nlev,nt),ddt_v_hadv_ls(nlev,nt) )

        DO n = 1,nt

          DO jk = nk , 1, -1
            READ(iunit,*,IOSTAT=ist)zz(jk),zu(jk),zv(jk),zw(jk),z_dt_temp_rad(jk), &
                & z_dt_qv_adv(jk),z_dt_temp_adv(jk),z_dt_u_adv(jk),z_dt_v_adv(jk)
            IF(ist/=success) CALL finish(routine, 'something wrong in forcing.dat')
          END DO

          !Skip lines
          IF(nt>1)THEN
            DO i = 1 , nskip
              READ(iunit,*,IOSTAT=ist)
            END DO
          END IF

          !Check if the file is written in descending order
          IF (zz(1) < zz(nk)) &
            CALL finish(routine, 'Write LS forcing data in descending order!')

          !Now perform interpolation to grid levels assuming:
          !a) linear interpolation
          !b) Beyond the last Z level the values are linearly extrapolated
          !c) Assuming model grid is flat-NOT on sphere

          CALL vert_intp_linear_1d(zz,zw,p_metrics%z_mc(1,:,1),w_ls(:,n))
          CALL vert_intp_linear_1d(zz,zu,p_metrics%z_mc(1,:,1),u_geo(:,n))
          CALL vert_intp_linear_1d(zz,zv,p_metrics%z_mc(1,:,1),v_geo(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_temp_adv,p_metrics%z_mc(1,:,1),ddt_temp_hadv_ls(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_temp_rad,p_metrics%z_mc(1,:,1),ddt_temp_rad(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_qv_adv,p_metrics%z_mc(1,:,1),ddt_qv_hadv_ls(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_u_adv,p_metrics%z_mc(1,:,1),ddt_u_hadv_ls(:,n))
          CALL vert_intp_linear_1d(zz,z_dt_v_adv,p_metrics%z_mc(1,:,1),ddt_v_hadv_ls(:,n))

        END DO !n

        CLOSE(iunit)

        DEALLOCATE( zz, zw, zu, zv, z_dt_temp_adv, z_dt_temp_rad, z_dt_qv_adv,z_dt_u_adv,z_dt_v_adv )

      ENDIF

      WRITE(message_text,*) dt_forcing
      CALL message('Time varying LS forcing read in:',message_text)

    END IF !is_ls_forcing


!------------------------------------------------------------------------------
! READ NUDGING FILE
!------------------------------------------------------------------------------

    IF(is_nudging) THEN

      IF(i_scm_netcdf > 0) THEN

        IF(i_scm_netcdf == 2) THEN

!------------------------------------------------------------------------------
! Nudging with DEPHY unified SCM NETCDF file
! ... Attention: the DEPHY format doesn't yet support time dependent prognostic
!                variables (T, qv, ...).  They are dimension with t0 not time.

          CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
            & TRIM(routine)//'   File init_SCM.nc cannot be opened')

          CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), TRIM(routine)//' lev')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nk), routine)

          CALL nf (nf90_inq_dimid        (fileid, 't0' , dimid), routine)
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = t0) , routine)

          CALL nf (nf90_inq_dimid        (fileid, 'time', dimid), TRIM(routine)//' nt')    ! dimension 'time'
          CALL nf (nf90_inquire_dimension(fileid, dimid , len = nt), routine)

          WRITE(message_text,'(a,i5,a,i5)') 'SCM nudging input: nk', nk, '  nt', nt
          CALL message (routine,message_text)

          ALLOCATE( zz_nc(nk,nt), time_nc(nt), zu_nc(nk,nt), zv_nc(nk,nt), ztheta_nc(nk,nt), z_qv_nc(nk,nt), &
                  & tempf(nk,nt) )
          zz_nc=0.0_wp ; time_nc=0.0_wp ; zu_nc=0.0_wp ; zv_nc=0.0_wp ; ztheta_nc=0.0_wp ; z_qv_nc=0.0_wp

          ALLOCATE( u_nudg(nlev,nt), v_nudg(nlev,nt), theta_nudg(nlev,nt), qv_nudg(nlev,nt) )
          u_nudg=0.0_wp ; v_nudg=0.0_wp ; theta_nudg=0.0_wp ; qv_nudg=0.0_wp

          CALL nf (nf90_inq_varid (fileid, 'zh_forc', varid), TRIM(routine)//' height')
          CALL nf (nf90_get_var   (fileid, varid,tempf), routine)
          zz_nc=tempf(nk:1:-1,:)

          !Check if the file is written in descending order
          IF(zz_nc(1,1) < zz_nc(nk,1)) CALL finish ( routine, 'Write LS forcing data in descending order!')

          CALL nf(nf90_inq_varid(fileid, 'time', varid), TRIM(routine)//' time')
          CALL nf(nf90_get_var(fileid, varid, time_nc), routine)
          time_unit = ''               !necessary for formatting
          CALL nf(nf90_get_att(fileid, varid, 'units', time_unit), TRIM(routine)//' units')

          time_unit_short = time_unit(1:INDEX(time_unit,"since")-2)  ! e.g. "minutes since 2020-1-1 00:00:00"
          IF ( get_my_global_mpi_id() == 0 ) THEN
            WRITE(*,*) 'time unit in init_SCM.nc, long: ', TRIM(time_unit), '   short: ', TRIM(time_unit_short), &
              & '     times: ', time_nc
          END IF

          SELECT CASE (time_unit_short)
          CASE ('s', 'sec', 'second', 'seconds')
            time_unit_in_sec = 1
          CASE ('min', 'minute', 'minutes')
            time_unit_in_sec = 60
          CASE ('h', 'hour', 'hours')
            time_unit_in_sec = 3600
          CASE ('d', 'day', 'days')
            time_unit_in_sec = 84600
          CASE DEFAULT
            CALL finish ( routine, 'time unit in init_SCM.nc not s/min/h!')
          END SELECT

          time_nc = time_nc * time_unit_in_sec  ! conversion to [s] from min/hours/days

          IF (nt<=1) THEN
            dt_nudging = -999._wp
          ELSE
            dt_nudging = time_nc(2) - time_nc(1)
            DO i=2,size(time_nc)-1
              IF ((time_nc(i+1)-time_nc(i)) /= dt_nudging) THEN
                CALL finish(routine,'input timesteps of equal length needed')
              END IF
            END DO
          ENDIF

          IF ( dt_relax_uv <= 0.0_wp ) THEN   ! dt_relax_uv not defined by namelist input
            CALL nf (nf90_inq_varid (fileid, 'dt_relax_uv', varid), TRIM(routine)//' dt_relax_uv')
            CALL nf (nf90_get_var   (fileid, varid, temp_nf), routine)
            dt_relax_uv = temp_nf(1)
            write(*,*) 'dt_relax_uv [s]: ', dt_relax_uv
          END IF
          IF ( dt_relax_t <= 0.0_wp ) THEN   ! dt_relax_t not defined by namelist input
            CALL nf (nf90_inq_varid (fileid, 'dt_relax_t', varid), TRIM(routine)//' dt_relax_t')
            CALL nf (nf90_get_var   (fileid, varid, temp_nf), routine)
            dt_relax_t = temp_nf(1)
            write(*,*) 'dt_relax_t [s]: ', dt_relax_t
          END IF
          IF ( dt_relax_q <= 0.0_wp ) THEN   ! dt_relax_q not defined by namelist input
            CALL nf (nf90_inq_varid (fileid, 'dt_relax_q', varid), TRIM(routine)//' dt_relax_q')
            CALL nf (nf90_get_var   (fileid, varid, temp_nf), routine)
            dt_relax_q = temp_nf(1)
            write(*,*) 'dt_relax_q [s]: ', dt_relax_q
          END IF

          nf_status = nf90_inq_varid (fileid, 'ua_nud',  varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'ua_nud not available in init_SCM.nc.  It will be set to 0.')
            zu_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf)
            zu_nc=tempf(nk:1:-1,:)
          END IF

          nf_status = nf90_inq_varid (fileid, 'va_nud',  varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'va_nud not available in init_SCM.nc.  It will be set to 0.')
            zv_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf)
            zv_nc=tempf(nk:1:-1,:)
          END IF

          IF (theta_nudging) THEN
            nf_status = nf90_inq_varid (fileid, 'theta_nud',  varid)
            IF (nf_status /= nf90_noerr) THEN
              CALL message(routine,'theta_nud not available in init_SCM.nc.  It will be set to 0.')
              ztheta_nc=0.0_wp
            ELSE
              CALL message(routine,'Using potential temperature nudging with theta_nud.')
              nf_status2 = nf90_get_var(fileid, varid, tempf)
              ztheta_nc=tempf(nk:1:-1,:)
            END IF
          ELSE
            ! T nudging
            nf_status  = nf90_inq_varid (fileid, 'ta_nud',  varid)
            IF (nf_status /= nf90_noerr) THEN
              CALL message(routine,'ta_nud not available in init_SCM.nc.  It will be set to 0.')
              ztheta_nc=0.0_wp
            ELSE
              CALL message(routine,'Using temperature nudging with ta_nud.')
              nf_status2 = nf90_get_var(fileid, varid, tempf)
              ztheta_nc=tempf(nk:1:-1,:)
            END IF
          ENDIF

          nf_status = nf90_inq_varid (fileid, 'qv_nud',  varid)
          IF (nf_status /= nf90_noerr) THEN
            CALL message(routine,'qv_nud not available in init_SCM.nc.  It will be set to 0.')
            z_qv_nc=0.0_wp
          ELSE
            nf_status2 = nf90_get_var(fileid, varid, tempf)
            z_qv_nc=tempf(nk:1:-1,:)
          END IF

         !write(*,*) 'uLS:' , zu_nc(:,1)
         !write(*,*) 'vLS:' , zv_nc(:,1)
         !write(*,*) 'qvLS:', z_qv_nc(:,1)
         !write(*,*) 'thLS:', ztheta_nc(:,1)

          CALL nf (nf90_close(fileid), routine)

          DO n = 1 , nt
            !Now perform interpolation to grid levels assuming:
            !a) linear interpolation
            !b) Beyond the last Z level the values are linearly extrapolated
            !c) Assuming model grid is flat-NOT on sphere

            CALL vert_intp_linear_1d(zz_nc(:,n),zu_nc(:,n)    ,p_metrics%z_mc(1,:,1),u_nudg(:,n)    )
            CALL vert_intp_linear_1d(zz_nc(:,n),zv_nc(:,n)    ,p_metrics%z_mc(1,:,1),v_nudg(:,n)    )
            CALL vert_intp_linear_1d(zz_nc(:,n),ztheta_nc(:,n),p_metrics%z_mc(1,:,1),theta_nudg(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_qv_nc(:,n)  ,p_metrics%z_mc(1,:,1),qv_nudg(:,n)   )
          END DO

          DEALLOCATE( zz_nc, time_nc, zu_nc, zv_nc, ztheta_nc, z_qv_nc, tempf )

!------------------------------------------------------------------------------
! Nudging with normal NETCDF file (including real ICON NWP forcing)

        ELSE

          CALL nf (nf90_open('init_SCM.nc', NF90_NOWRITE, fileid), &
            & TRIM(routine)//'   File init_SCM.nc cannot be opened')

          CALL nf (nf90_inq_dimid(fileid, 'lev', dimid), TRIM(routine)//' lev')
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nk), routine)

          CALL nf (nf90_inq_dimid(fileid, 'nt', dimid), TRIM(routine)//' nt')   ! dimension 'nt'
          CALL nf (nf90_inquire_dimension(fileid, dimid, len = nt), routine)

          WRITE(message_text,'(a,i5,a,i5)') 'SCM nudging input: nk', nk, '  nt', nt
          CALL message(routine,message_text)

          ALLOCATE( zz_nc(nk,nt), time_nc(nt), zu_nc(nk,nt), zv_nc(nk,nt), ztheta_nc(nk,nt), z_qv_nc(nk,nt) )
          zz_nc=0.0_wp ; time_nc=0.0_wp ; zu_nc=0.0_wp ; zv_nc=0.0_wp ; ztheta_nc=0.0_wp ; z_qv_nc=0.0_wp

          ALLOCATE( u_nudg(nlev,nt), v_nudg(nlev,nt), theta_nudg(nlev,nt), qv_nudg(nlev,nt) )
          u_nudg=0.0_wp ; v_nudg=0.0_wp ; theta_nudg=0.0_wp ; qv_nudg=0.0_wp

          CALL nf (nf90_inq_varid(fileid, 'height', varid), TRIM(routine)//' height')
          CALL nf (nf90_get_var  (fileid, varid, zz_nc), routine)
!         write(*,*) 'zz_nc',zz_nc
          !Check if the file is written in descending order
          IF(zz_nc(1,1) < zz_nc(nk,1)) CALL finish ( routine, 'Write LS forcing data in descending order!')

          CALL nf(nf90_inq_varid(fileid, 'time', varid), TRIM(routine)//' time')
          CALL nf(nf90_get_var  (fileid, varid, time_nc), routine)
          time_unit = ''               !necessary for formatting
          CALL nf(nf90_get_att(fileid, varid, 'units', time_unit), TRIM(routine)//' units')

          time_unit_short = time_unit(1:INDEX(time_unit,"since")-2)  ! e.g. "minutes since 2020-1-1 00:00:00"
          IF ( get_my_global_mpi_id() == 0 ) THEN
            WRITE(*,*) 'time unit in init_SCM.nc, long: ', TRIM(time_unit), '   short: ', TRIM(time_unit_short), &
              & '     times: ', time_nc
          END IF

          SELECT CASE (time_unit_short)
          CASE ('s', 'sec', 'second', 'seconds')
            time_unit_in_sec = 1
          CASE ('min', 'minute', 'minutes')
            time_unit_in_sec = 60
          CASE ('h', 'hour', 'hours')
            time_unit_in_sec = 3600
          CASE ('d', 'day', 'days')
            time_unit_in_sec = 84600
          CASE DEFAULT
            CALL finish ( routine, 'time unit in init_SCM.nc not s/min/h!')
          END SELECT

          time_nc = time_nc * time_unit_in_sec  ! conversion to [s] from min/hours/days

          IF (nt<=1) THEN
            dt_nudging = -999._wp
          ELSE
            dt_nudging = time_nc(2) - time_nc(1)
            DO i=2,size(time_nc)-1
              IF ((time_nc(i+1)-time_nc(i)) /= dt_nudging) THEN
                CALL finish(routine,'input timesteps of equal length needed')
              END IF
            END DO
          ENDIF

          IF ( dt_relax_uv <= 0.0_wp ) THEN   ! dt_relax_uv not defined by namelist input
            CALL nf (nf90_inq_varid (fileid, 'dt_relax_uv', varid), TRIM(routine)//' dt_relax_uv')
            CALL nf (nf90_get_var   (fileid, varid, temp_nf), routine)
            dt_relax_uv = temp_nf(1)
            write(*,*) 'dt_relax_uv [s]: ', dt_relax_uv
          END IF
          IF ( dt_relax_t <= 0.0_wp ) THEN   ! dt_relax_t not defined by namelist input
            CALL nf (nf90_inq_varid (fileid, 'dt_relax_t', varid), TRIM(routine)//' dt_relax_t')
            CALL nf (nf90_get_var   (fileid, varid, temp_nf), routine)
            dt_relax_t = temp_nf(1)
            write(*,*) 'dt_relax_t [s]: ', dt_relax_t
          END IF
          IF ( dt_relax_q <= 0.0_wp ) THEN   ! dt_relax_q not defined by namelist input
            CALL nf (nf90_inq_varid (fileid, 'dt_relax_q', varid), TRIM(routine)//' dt_relax_q')
            CALL nf (nf90_get_var   (fileid, varid, temp_nf), routine)
            dt_relax_q = temp_nf(1)
            write(*,*) 'dt_relax_q [s]: ', dt_relax_q
          END IF

          CALL nf (nf90_inq_varid (fileid, 'uIN',  varid), TRIM(routine)//' uIN')
          CALL nf (nf90_get_var   (fileid, varid,zu_nc) , routine)

          CALL nf (nf90_inq_varid (fileid, 'vIN',  varid), TRIM(routine)//' vIN')
          CALL nf (nf90_get_var   (fileid, varid, zv_nc), routine)

          CALL nf (nf90_inq_varid (fileid, 'qvIN', varid), TRIM(routine)//' qvIN')
          CALL nf (nf90_get_var   (fileid, varid, z_qv_nc), routine)

          CALL nf (nf90_inq_varid (fileid, 'thIN', varid), TRIM(routine)//' thIN')
          CALL nf (nf90_get_var   (fileid, varid, ztheta_nc), routine)

!         write(*,*) 'uLS',zu_nc(:,1)
!         write(*,*) 'vLS',zv_nc(:,1)
!         write(*,*) 'qvLS',z_qv_nc(:,1)
!         write(*,*) 'thLS',ztheta_nc(:,1)

          CALL nf (nf90_close(fileid), routine)

          DO n = 1 , nt
            !Now perform interpolation to grid levels assuming:
            !a) linear interpolation
            !b) Beyond the last Z level the values are linearly extrapolated
            !c) Assuming model grid is flat-NOT on sphere

            CALL vert_intp_linear_1d(zz_nc(:,n),zu_nc(:,n)    ,p_metrics%z_mc(1,:,1),u_nudg(:,n)    )
            CALL vert_intp_linear_1d(zz_nc(:,n),zv_nc(:,n)    ,p_metrics%z_mc(1,:,1),v_nudg(:,n)    )
            CALL vert_intp_linear_1d(zz_nc(:,n),ztheta_nc(:,n),p_metrics%z_mc(1,:,1),theta_nudg(:,n))
            CALL vert_intp_linear_1d(zz_nc(:,n),z_qv_nc(:,n)  ,p_metrics%z_mc(1,:,1),qv_nudg(:,n)   )
          END DO

          DEALLOCATE( zz_nc, time_nc, zu_nc, zv_nc, ztheta_nc, z_qv_nc )

        ENDIF

!------------------------------------------------------------------------------
! Nudging with formatted ASCII file

      ELSE

        !Open formatted file to read nudging data
        !The order of file assumed is: Z(m) - tau(s) - u(m/s) - v(m/s) - pot. temp(K) - qv(kg/kg)
        iunit = find_next_free_unit(10,20)
        OPEN (unit=iunit,file='nudging.dat',access='SEQUENTIAL', &
              form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)

        IF(ist/=success)THEN
          CALL finish (routine, 'open nudging.dat failed')
        ENDIF


        !Skip the first line and read next 2 lines with information about vertical and time levels
        READ(iunit,*,IOSTAT=ist)                         !skip
        READ(iunit,*,IOSTAT=ist)nk,dt_nudging,end_time   !vertical levels,time levels and interval
        IF (ist/=success) CALL finish(routine, &
          'Must provide vertical level, time level and time interval info in nudging file')

        IF (nt>1) THEN
          READ(iunit,*,IOSTAT=ist)nskip              !lines to skip between successive time levels
          IF (ist/=success) CALL finish(routine, &
            'Time levels > 1 so must provide number lines to skip between successive time levels in nudging file')
        END IF
        IF (nt<=1) dt_nudging=-999._wp

        ALLOCATE( zz(nk), zu(nk), zv(nk), ztheta(nk), z_qv(nk) )

        ALLOCATE( u_nudg(nlev,nt), v_nudg(nlev,nt), theta_nudg(nlev,nt), qv_nudg(nlev,nt))

        DO n = 1 , nt

          DO jk = nk , 1, -1
            IF ( dt_relax_t <= 0.0_wp ) THEN   ! dt_relax not defined by namelist input
              READ(iunit,*,IOSTAT=ist)zz(jk),dt_relax_t,zu(jk),zv(jk),ztheta(jk),z_qv(jk)
              IF (ist/=success) CALL finish(routine, 'something wrong in nudging.dat')
            END IF
          END DO

          !Skip lines
          IF(nt>1)THEN
            DO i = 1 , nskip
              READ(iunit,*,IOSTAT=ist)
            END DO
          END IF

          !Check if the file is written in descending order
          IF (zz(1) < zz(nk)) &
            CALL finish(routine, 'Write nuding data in descending order!')

          !Now perform interpolation to grid levels assuming:
          !a) linear interpolation
          !b) Beyond the last Z level the values are linearly extrapolated
          !c) Assuming model grid is flat-NOT on sphere

          CALL vert_intp_linear_1d(zz,zu    ,p_metrics%z_mc(1,:,1),u_nudg(:,n))
          CALL vert_intp_linear_1d(zz,zv    ,p_metrics%z_mc(1,:,1),v_nudg(:,n))
          CALL vert_intp_linear_1d(zz,ztheta,p_metrics%z_mc(1,:,1),theta_nudg(:,n))
          CALL vert_intp_linear_1d(zz,z_qv  ,p_metrics%z_mc(1,:,1),qv_nudg(:,n))

        END DO

        CLOSE(iunit)

        DEALLOCATE( zz, zu, zv, ztheta, z_qv )

      END IF

      WRITE(message_text,*) dt_nudging
      CALL message('Time varying nudging read in with time step:', message_text)

    END IF !is_nudging


  END SUBROUTINE init_ls_forcing

  !>
  !! apply_ls_forcing
  !!------------------------------------------------------------------------
  !! Apply large-scale forcing: called in the end from mo_nh_interface_nwp
  !! It uses the most updated u,v to calculate the large-scale subsidence induced
  !! advective tendencies. All other tendencies don't need any computation.
  !! All tendencies are then accumulated in the slow physics tendency terms
  !!
  !!------------------------------------------------------------------------

  SUBROUTINE apply_ls_forcing(p_patch, p_metrics, curr_sim_time,       & !in
                              p_prog, p_diag, qv, rl_start, rl_end,    & !in
                              ddt_u_ls, ddt_v_ls,                      & !out
                              ddt_temp_ls, ddt_qv_ls,                  & !out
                              ddt_temp_subs_ls, ddt_qv_subs_ls,        & !output
                              ddt_temp_adv_ls, ddt_qv_adv_ls,          & !output
                              ddt_u_adv_ls, ddt_v_adv_ls, wsub,        & !output
                              temp_nudge, u_nudge, v_nudge, q_nudge)     !output

    TYPE(t_patch),INTENT(inout),     TARGET :: p_patch
    TYPE(t_nh_metrics),INTENT(in),   TARGET :: p_metrics
    TYPE(t_nh_prog),INTENT(in),      TARGET :: p_prog
    TYPE(t_nh_diag),INTENT(in),      TARGET :: p_diag
    INTEGER,        INTENT(in)              :: rl_start
    INTEGER,        INTENT(in)              :: rl_end
    REAL(wp),       INTENT(in)              :: curr_sim_time
    REAL(wp),       INTENT(in)              :: qv(:,:,:)
    !Individual LS tendencies for profile output
    REAL(wp),       INTENT(out)             :: ddt_u_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_v_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_subs_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_subs_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_temp_adv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_qv_adv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_u_adv_ls(:)
    REAL(wp),       INTENT(out)             :: ddt_v_adv_ls(:)
    REAL(wp),       INTENT(out)             :: wsub(:)
    REAL(wp),       INTENT(out)              :: temp_nudge(:)
    REAL(wp),       INTENT(out)              :: u_nudge(:)
    REAL(wp),       INTENT(out)              :: v_nudge(:)
    REAL(wp),       INTENT(out)              :: q_nudge(:,:)


    CHARACTER(len=*), PARAMETER :: routine = 'mo_ls_forcing:apply_ls_forcing'
    REAL(wp) :: thetain(nproma,p_patch%nlev,p_patch%nblks_c)
    REAL(wp) :: inv_no_gb_cells
    REAL(wp), DIMENSION(p_patch%nlev)  :: u_gb, v_gb, theta_gb, qv_gb
    REAL(wp), DIMENSION(p_patch%nlev+1):: u_gb_hl, v_gb_hl, theta_gb_hl, qv_gb_hl
    REAL(wp), DIMENSION(p_patch%nlev)  :: inv_dz, exner_gb
    REAL(wp), DIMENSION(p_patch%nlev)  :: z_ddt_t_rad, z_ugeo, z_vgeo

    REAL(wp) :: int_weight
    INTEGER  :: i_startblk, i_endblk, jk, nlev, jb, jc
    INTEGER  :: n_curr, n_next, i_startidx, i_endidx
    REAL(wp) :: corio_lat_scm

    nlev      = p_patch%nlev
    inv_no_gb_cells = 1._wp / REAL(p_patch%n_patch_cells_g,wp)

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !Find where in the array of LS forcing and nudging current time stands
    !and do linear interpolation in time
    IF (dt_forcing > 0._wp) THEN
      n_curr     = FLOOR(curr_sim_time/dt_forcing)+1
      n_next     = n_curr+1
      int_weight = curr_sim_time/dt_forcing-n_curr+1
      !ibd, exception for last time step
      IF (int_weight == 0._wp) THEN
        n_next = n_curr
      END IF

      ! Christopher Moseley: temporary catch
      IF (int_weight < 0 .OR. int_weight > 1) THEN
        WRITE(message_text,'(a,2(i0,", "),g0)') &
          'INTERPOLATION ERROR in LS forcing: ', n_curr,n_next,int_weight
        CALL message(routine, message_text)
      END IF

      wsub            = w_ls            (:,n_curr)*(1.-int_weight) + w_ls            (:,n_next)*int_weight
      ddt_qv_adv_ls   = ddt_qv_hadv_ls  (:,n_curr)*(1.-int_weight) + ddt_qv_hadv_ls  (:,n_next)*int_weight
      ddt_u_adv_ls    = ddt_u_hadv_ls   (:,n_curr)*(1.-int_weight) + ddt_u_hadv_ls   (:,n_next)*int_weight
      ddt_v_adv_ls    = ddt_v_hadv_ls   (:,n_curr)*(1.-int_weight) + ddt_v_hadv_ls   (:,n_next)*int_weight
      ddt_temp_adv_ls = ddt_temp_hadv_ls(:,n_curr)*(1.-int_weight) + ddt_temp_hadv_ls(:,n_next)*int_weight
      z_ddt_t_rad     = ddt_temp_rad    (:,n_curr)*(1.-int_weight) + ddt_temp_rad    (:,n_next)*int_weight
      z_ugeo          = u_geo           (:,n_curr)*(1.-int_weight) + u_geo           (:,n_next)*int_weight
      z_vgeo          = v_geo           (:,n_curr)*(1.-int_weight) + v_geo           (:,n_next)*int_weight
      lat_scm         = bnd_lat         (n_curr)  *(1.-int_weight) + bnd_lat         (n_next)  *int_weight
      lon_scm         = bnd_lon         (n_curr)  *(1.-int_weight) + bnd_lon         (n_next)  *int_weight
    ELSE
      wsub            = w_ls(:,1)
      ddt_qv_adv_ls   = ddt_qv_hadv_ls(:,1)
      ddt_temp_adv_ls = ddt_temp_hadv_ls(:,1)
      ddt_u_adv_ls    = ddt_u_hadv_ls(:,1)
      ddt_v_adv_ls    = ddt_v_hadv_ls(:,1)
      z_ddt_t_rad     = ddt_temp_rad(:,1)
      z_ugeo          = u_geo(:,1)
      z_vgeo          = v_geo(:,1)
      lat_scm         = bnd_lat(1)
      lon_scm         = bnd_lon(1)
    END IF


! 3D forcing variable calculations

    !0) Initialize all passed ddt's to 0
!$OMP PARALLEL
    CALL init(ddt_u_ls   , lacc=.FALSE.)
    CALL init(ddt_v_ls   , lacc=.FALSE.)
    CALL init(ddt_temp_ls, lacc=.FALSE.)
    call init(ddt_qv_ls  , lacc=.FALSE.)

    !use theta instead of temperature for subsidence (Anurag, Christopher)
!$OMP BARRIER
!$OMP DO PRIVATE(jc,jb,jk,i_startidx,i_endidx)
    DO jb = i_startblk,i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
                            i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1 , nlev
        DO jc = i_startidx, i_endidx
          thetain(jc,jk,jb)  = p_diag%temp(jc,jk,jb)/p_prog%exner(jc,jk,jb)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

    !Global means use > 50% of CPU time!!!
    IF (is_subsidence_moment .or. is_ls_coriolis) THEN
      CALL levels_horizontal_mean(p_diag%u  , p_patch%cells%area, p_patch%cells%owned, u_gb    (1:nlev))
      CALL levels_horizontal_mean(p_diag%v  , p_patch%cells%area, p_patch%cells%owned, v_gb    (1:nlev))
    ENDIF
    IF (is_subsidence_heat) THEN
      CALL levels_horizontal_mean(thetain   , p_patch%cells%area, p_patch%cells%owned, theta_gb(1:nlev))
      CALL levels_horizontal_mean(qv        , p_patch%cells%area, p_patch%cells%owned, qv_gb   (1:nlev))
    ENDIF
    CALL levels_horizontal_mean(p_prog%exner, p_patch%cells%area, p_patch%cells%owned, exner_gb(1:nlev))

    IF (is_subsidence_moment.OR.is_subsidence_heat) THEN
      inv_dz(:) = 1._wp / p_metrics%ddqz_z_full(1,:,1)
    END IF


    !1a) Horizontal interpolation and their advective tendency - momentum

    IF (is_subsidence_moment) THEN
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),u_gb,p_metrics%z_ifc(1,:,1),u_gb_hl)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),v_gb,p_metrics%z_ifc(1,:,1),v_gb_hl)

      ddt_u_ls =  ddt_u_ls - wsub*vertical_derivative(u_gb_hl,inv_dz)
      ddt_v_ls =  ddt_v_ls - wsub*vertical_derivative(v_gb_hl,inv_dz)
    END IF

    !1b) Horizontal interpolation and their advective tendency - temperature and water vapor

    IF (is_subsidence_heat) THEN
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),theta_gb,p_metrics%z_ifc(1,:,1),theta_gb_hl)
      CALL vert_intp_linear_1d(p_metrics%z_mc(1,:,1),qv_gb   ,p_metrics%z_ifc(1,:,1),qv_gb_hl)

      ddt_temp_subs_ls = -wsub*vertical_derivative(theta_gb_hl,inv_dz)  ! vertical advective tendency for theta
      ddt_qv_subs_ls   = -wsub*vertical_derivative(qv_gb_hl   ,inv_dz)  ! vertical advective tendenciy for qv

      ddt_temp_ls = ddt_temp_ls + ddt_temp_subs_ls
      ddt_qv_ls   = ddt_qv_ls   + ddt_qv_subs_ls
    END IF

    !2) Horizontal advective forcing: at present only for tracers and temperature
    IF (is_advection) THEN
      IF(is_advection_tq)THEN
        ddt_qv_ls   = ddt_qv_ls   + ddt_qv_adv_ls
        ddt_temp_ls = ddt_temp_ls + ddt_temp_adv_ls
      ENDIF
      IF (is_advection_uv) THEN
        ddt_u_ls = ddt_u_ls + ddt_u_adv_ls
        ddt_v_ls = ddt_v_ls + ddt_v_adv_ls
      ENDIF
    END IF

    !3) Coriolis and geostrophic wind
    ! Coriolis parameter based on time-varying latitude read from large-scale forcing
    corio_lat_scm = 2._wp*grid_angular_velocity*SIN(lat_scm * deg2rad)

    ! If dynamical core is calculating Coriolis term on beta plain, update latitude from LS forcing
    IF (lcoriolis) THEN
      p_patch%cells%f_c(:,:) = corio_lat_scm
      p_patch%edges%f_e(:,:) = corio_lat_scm
      p_patch%verts%f_v(:,:) = corio_lat_scm
    ENDIF

    ! If pressure gradient term is calculated from geostrophic wind, add here
    IF (is_geowind) THEN
      ddt_u_ls = ddt_u_ls - corio_lat_scm * z_vgeo
      ddt_v_ls = ddt_v_ls + corio_lat_scm * z_ugeo
    END IF

    ! If dynamics does not calculate Coriolis term, calculate here based on domain-average winds
    IF (is_ls_coriolis) THEN
      ddt_u_ls = ddt_u_ls + corio_lat_scm * v_gb
      ddt_v_ls = ddt_v_ls - corio_lat_scm * u_gb
    ENDIF

    !4) Radiative forcing
    IF (is_rad_forcing) THEN
      ddt_temp_ls = ddt_temp_ls + z_ddt_t_rad
    END IF

    !5) Nudging
    IF (is_nudging) THEN

      WRITE(message_text,*) dt_nudging
      CALL message('dt_nudging (time-step input data) =',message_text)

      ! multiple time steps available separated by dt_nudging
      IF (dt_nudging > 0._wp) THEN
        n_curr = FLOOR(curr_sim_time/dt_nudging)+1
        n_next = n_curr+1
        int_weight = curr_sim_time/dt_nudging-n_curr+1

        !ibd, exception for last time step
        IF (int_weight.eq.0._wp) THEN
          n_next = n_curr
        END IF

        ! stop if end of data reached
        IF (n_next .GT. nt) THEN
          CALL finish ( routine, 'end of SCM input file for nudging' )
        END IF

        ! Interpolation error catch
        IF (int_weight.LT.0 .OR. int_weight.GT.1) THEN
          WRITE (message_text,'(a,2(i0,", "),g0)') &
               'INTERPOLATION ERROR in nudging:', n_curr,n_next,int_weight
          CALL message(routine, message_text)
        END IF

        IF (is_nudging_uv) THEN
          u_nudge         = u_nudg(:,n_curr)*(1.-int_weight) &
                        & + u_nudg(:,n_next)*    int_weight
          v_nudge         = v_nudg(:,n_curr)*(1.-int_weight) &
                        & + v_nudg(:,n_next)*    int_weight
        ENDIF
        IF(is_nudging_t) THEN
          temp_nudge      = theta_nudg(:,n_curr)*(1.-int_weight) &
                        & + theta_nudg(:,n_next)*    int_weight
        ENDIF
        IF(is_nudging_q) THEN
          q_nudge(:,1)    = qv_nudg   (:,n_curr)*(1.-int_weight) &
                        & + qv_nudg   (:,n_next)*    int_weight
        ENDIF
      ! single time step available
      ELSE
        IF (is_nudging_uv) THEN
          u_nudge         = u_nudg(:,1)
          v_nudge         = v_nudg(:,1)
        ENDIF
        IF (is_nudging_t) THEN
          temp_nudge      = theta_nudg(:,1)
        ENDIF
        IF (is_nudging_q) THEN
          q_nudge(:,1)    = qv_nudg   (:,1)
        ENDIF
      END IF

      q_nudge(:,2) = 0.0_wp         ! qc undefined
      q_nudge(:,3) = 0.0_wp         ! qi undefined

    END IF

    !Convert theta tendency to temp at once
    ddt_temp_ls = exner_gb * ddt_temp_ls
    IF (theta_nudging) THEN
      ! only necessary if not nudging against temperature profile
      temp_nudge  = exner_gb * temp_nudge
    ENDIF

  END SUBROUTINE apply_ls_forcing

  !-------------------------------------------------------------------------------
  ! Apply large-scale forcing with/without nudging to moisture in case SCM/LES is
  ! run with large-scale forcing, and nudging is enabled.
  !-------------------------------------------------------------------------------


  SUBROUTINE apply_ls_forc_nudge_qv(i_startidx, i_endidx, &
       & dt_loc, nqtendphy, pdtime, kstart_moist, kend,   &
       & geopot_agl, ddt_tracer_ls, q_nudge, tracer, lacc)


    INTEGER,  INTENT(IN)    :: i_startidx, i_endidx ! block start/end indices
    INTEGER,  INTENT(IN)    :: nqtendphy    !< number of water species for which physical tendencies are stored
    REAL(wp), INTENT(IN)    :: pdtime       !< time step
    REAL(wp), INTENT(IN)    :: dt_loc       !< (advective) time step applicable to local grid level
    INTEGER,  INTENT(IN)    :: kstart_moist !< vertical start index
    INTEGER,  INTENT(IN)    :: kend         !< vertical end index
    REAL(wp), INTENT(IN)    :: geopot_agl(:,:)    !< geopotential, height above ground
    REAL(wp), INTENT(IN)    :: ddt_tracer_ls(:,:) !< tendency of tracers due to LS forcing
    REAL(wp), INTENT(IN)    :: q_nudge(:,:)       !< humidity tendency due to qv nudging
    REAL(wp), INTENT(INOUT) :: tracer(:,:,:)      !< tracer array (qv, qi, qc ...)
    LOGICAL,  INTENT(IN), OPTIONAL :: lacc  !< If true, use openacc

    !Local variables
    REAL(wp) :: nudgecoeff
    REAL(wp) :: z_ddt_q_nudge
    INTEGER  :: jt, jk, jc

    CALL assert_acc_host_only("apply_ls_forc_nudge_qv", lacc)
    DO jt=1, nqtendphy  ! qv,qc,qi
      DO jk = kstart_moist, kend
!DIR$ IVDEP
        DO jc = i_startidx, i_endidx

          ! add q nudging (T, U, V nudging is in mo_nh_interface_nwp)
          IF ( is_nudging_q ) THEN

            ! linear nudging profile between "start" and "full" heights - prevent sfc layer instability
            IF ( nudge_full_height_q == nudge_start_height_q ) THEN
              IF ( geopot_agl(jc,jk)/grav >= nudge_start_height_q ) THEN
                nudgecoeff = 1.0_wp
              ELSE
                nudgecoeff = 0.0_wp
              ENDIF
            ELSE
              nudgecoeff = ( geopot_agl(jc,jk)/grav - nudge_start_height_q ) / &
                   & ( nudge_full_height_q          - nudge_start_height_q )
              nudgecoeff = MAX( MIN( nudgecoeff, 1.0_wp ), 0.0_wp )
            END IF
            ! analytic implicit: (q,n+1 - q,n) / dt = (q,nudge - q,n) / dt_relax * exp(-dt/dt_relax)
            z_ddt_q_nudge =                                                          &
                 &  - ( tracer(jc,jk,jt) - q_nudge(jk,jt) ) &
                 &  / dt_relax_q * exp(-dt_loc/dt_relax_q) * nudgecoeff
          ELSE
            z_ddt_q_nudge = 0.0_wp
          END IF

          tracer(jc,jk,jt) = MAX(0._wp, tracer(jc,jk,jt)   &
               &             + pdtime*ddt_tracer_ls(jk,jt) &
               &             + pdtime*z_ddt_q_nudge)

        ENDDO
      ENDDO
    END DO

  END SUBROUTINE apply_ls_forc_nudge_qv

  !--------------------------------------------------------------------------------
  ! Apply large-scale tendency with/without nudging to temperature and wind in case
  ! SCM/LES is run with large-scale forcing, and nudging is enabled.
  !--------------------------------------------------------------------------------

  SUBROUTINE  apply_ls_forc_nudge_uvt(i_startidx, i_endidx, nlev,&
       &  ddt_u_ls,ddt_v_ls,ddt_temp_ls,temp, &
       &  theta_v, u, v,qv, z_qsum, &
       &  ddt_tracer_ls, temp_nudge, u_nudge, v_nudge, &
       &  dt_loc,&
       &  geopot_agl,ddt_temp_sim_rad, &
       &  ddt_exner_phy, z_ddt_u_tot,z_ddt_v_tot, z_ddt_temp, z_ddt_alpha, lacc)


    INTEGER,  INTENT(IN)    :: i_startidx, i_endidx ! block start/end indices
    INTEGER,  INTENT(IN)    :: nlev                 ! number of vertical levels
    REAL(wp), INTENT(IN)    :: ddt_u_ls(:)          ! large-scale forcing u wind
    REAL(wp), INTENT(IN)    :: ddt_v_ls(:)          ! large-scale forcing v wind
    REAL(wp), INTENT(IN)    :: ddt_temp_ls(:)       ! large-scale forcing temperature
    REAL(wp), INTENT(IN)    :: temp(:,:)            ! temperature
    REAL(wp), INTENT(IN)    :: theta_v(:,:)         ! virtual potential temperature
    REAL(wp), INTENT(IN)    :: u(:,:)               ! zonal wind
    REAL(wp), INTENT(IN)    :: v(:,:)               ! meridional wind
    REAL(wp), INTENT(IN)    :: qv(:,:)              ! mixing ratio vapor
    REAL(wp), INTENT(IN)    :: z_qsum(:,:)          ! summand of virtual increment
    REAL(wp), INTENT(IN)    :: ddt_tracer_ls(:,:)   ! tracer tendency due to large-scale forcing
    REAL(wp), INTENT(IN)    :: temp_nudge(:)        ! nudging tendency temperature
    REAL(wp), INTENT(IN)    :: u_nudge(:)           ! nudging tendency u wind
    REAL(wp), INTENT(IN)    :: v_nudge(:)           ! nudging tendency v wind
    REAL(wp), INTENT(IN)    :: geopot_agl(:,:)      ! geopotential height above ground
    REAL(wp), INTENT(IN)    :: ddt_temp_sim_rad(:,:)! temp tendency due to radiation
    REAL(wp), INTENT(IN)    :: dt_loc               !< (advective) time step applicable to local grid level
    REAL(vp), INTENT(INOUT) :: ddt_exner_phy(:,:)   ! tendency of exner due to physics processes
    REAL(wp), INTENT(INOUT) :: z_ddt_u_tot(:,:)     ! horizontal wind tendencies
    REAL(wp), INTENT(INOUT) :: z_ddt_v_tot(:,:)     ! horizontal wind tendencies
    REAL(wp), INTENT(INOUT) :: z_ddt_temp(:,:)      ! temperature tendency
    REAL(wp), INTENT(INOUT) :: z_ddt_alpha(:,:)     ! tendency of virtual increment
    LOGICAL,  INTENT(IN), OPTIONAL :: lacc  !< If true, use openacc

    !Local variables
    REAL(wp) :: nudgecoeff_uv, nudgecoeff_t
    REAL(wp) :: z_ddt_q_nudge
    INTEGER :: jt, jk, jc

    CALL assert_acc_host_only("apply_ls_forc_nudge_uvt", lacc)
    DO jk = 1, nlev
      DO jc = i_startidx, i_endidx

        ! add u/v/T forcing tendency
        z_ddt_u_tot(jc,jk)   = z_ddt_u_tot(jc,jk)       &
             &                     + ddt_u_ls(jk)

        z_ddt_v_tot(jc,jk)   = z_ddt_v_tot(jc,jk)       &
             &                     + ddt_v_ls(jk)

        z_ddt_temp(jc,jk)       = z_ddt_temp(jc,jk)           &
             + ddt_temp_ls(jk)

        ! simplified radiation scheme
        IF(is_sim_rad) THEN
          z_ddt_temp(jc,jk)     = z_ddt_temp(jc,jk)           &
               + ddt_temp_sim_rad(jc,jk)
        ENDIF

        ! linear nudging profile between "start" and "full" heights - prevent sfc layer instability
        IF ( nudge_full_height_uv == nudge_start_height_uv ) THEN
          IF ( geopot_agl(jc,jk)/grav >= nudge_start_height_uv ) THEN
            nudgecoeff_uv = 1.0_wp
          ELSE
            nudgecoeff_uv = 0.0_wp
          ENDIF
        ELSE
          nudgecoeff_uv = ( geopot_agl(jc,jk)/grav - nudge_start_height_uv ) / &
               & ( nudge_full_height_uv                   - nudge_start_height_uv )
          nudgecoeff_uv = MAX( MIN( nudgecoeff_uv, 1.0_wp ), 0.0_wp )
        END IF

        ! linear nudging profile between "start" and "full" heights - prevent sfc layer instability
        IF ( nudge_full_height_t == nudge_start_height_t ) THEN
          IF (geopot_agl(jc,jk)/grav >= nudge_start_height_t ) THEN
            nudgecoeff_t = 1.0_wp
          ELSE
            nudgecoeff_t = 0.0_wp
          ENDIF
        ELSE
          nudgecoeff_t = ( geopot_agl(jc,jk)/grav - nudge_start_height_t ) / &
               & ( nudge_full_height_t                   - nudge_start_height_t )
          nudgecoeff_t = MAX( MIN( nudgecoeff_t, 1.0_wp ), 0.0_wp )
        END IF

        ! add u/v/T nudging
        IF ( is_nudging_uv ) THEN
          ! explicit:          (u,n+1 - u,n) / dt = (u,nudge - u,n)   / dt_relax
          ! implicit:          (u,n+1 - u,n) / dt = (u,nudge - u,n+1) / dt_relax
          ! analytic implicit: (u,n+1 - u,n) / dt = (u,nudge - u,n)   / dt_relax * exp(-dt/dt_relax)
          z_ddt_u_tot(jc,jk) = z_ddt_u_tot(jc,jk)       &
               &  - ( u(jc,jk) - u_nudge(jk) ) / dt_relax_uv * exp(-dt_loc/dt_relax_uv) &
               &  * nudgecoeff_uv

          z_ddt_v_tot(jc,jk) = z_ddt_v_tot(jc,jk)       &
               &  - ( v(jc,jk) - v_nudge(jk) ) / dt_relax_uv * exp(-dt_loc/dt_relax_uv) &
               &  * nudgecoeff_uv

        END IF

        ! attention: T nudging results in dt-step oscillation/instability!
        ! q nudging done in tracer_add_phytend in mo_util_phys

        IF ( is_nudging_t ) THEN
          ! explicit:          (T,n+1 - T,n) / dt = (T,nudge - T,n)   / dt_relax
          ! implicit:          (T,n+1 - T,n) / dt = (T,nudge - T,n+1) / dt_relax
          ! analytic implicit: (T,n+1 - T,n) / dt = (T,nudge - T,n)   / dt_relax * exp(-dt/dt_relax)
          z_ddt_temp(jc,jk)     = z_ddt_temp(jc,jk)           &
               &  - ( temp(jc,jk) - temp_nudge(jk) ) / dt_relax_t &
               &  * exp(-dt_loc/dt_relax_t) * nudgecoeff_t
        END IF

        ! Convert temperature tendency into Exner function tendency
        z_ddt_alpha(jc,jk) = z_ddt_alpha(jc,jk)                          &
             &                + vtmpc1 * ddt_tracer_ls(jk,iqv) &
             &                - ddt_tracer_ls(jk,iqc)          &
             &                - ddt_tracer_ls(jk,iqi)

        ddt_exner_phy(jc,jk) = rd_o_cpd / theta_v(jc,jk)           &
             &                             * (z_ddt_temp(jc,jk)                             &
             &                             *(1._wp + vtmpc1*qv(jc,jk)&
             &                             - z_qsum(jc,jk))                                 &
             &                             + temp(jc,jk) * z_ddt_alpha(jc,jk) )

      END DO  ! jc
    END DO  ! jk


  END SUBROUTINE apply_ls_forc_nudge_uvt
!-----------------------------------------------------------------------------------

  !>
  !! apply_ls_sfc_forcing
  !!------------------------------------------------------------------------
  !! Apply large-scale surface forcing: surface variables from forcing are
  !! read and interpolated here. This has to occur before turbulence
  !! calculations, hence two calls to apply_ls_*_forcing.
  !!------------------------------------------------------------------------

  SUBROUTINE apply_ls_sfc_forcing(p_patch, p_metrics, curr_sim_time,   & !in
                              rl_start, rl_end,                        & !in
                              fc_sfc_lat_flx,fc_sfc_sens_flx,fc_ts,    & !output
                              fc_tg, fc_qvs,fc_Ch,fc_Cq,fc_Cm,fc_ustar)  !output

    TYPE(t_patch),INTENT(inout),     TARGET :: p_patch
    TYPE(t_nh_metrics),INTENT(in),   TARGET :: p_metrics
    INTEGER,        INTENT(in)              :: rl_start
    INTEGER,        INTENT(in)              :: rl_end
    REAL(wp),       INTENT(in)              :: curr_sim_time
    !Individual LS tendencies for profile output
    REAL(wp),       INTENT(out)             :: fc_sfc_lat_flx
    REAL(wp),       INTENT(out)             :: fc_sfc_sens_flx
    REAL(wp),       INTENT(out)             :: fc_ts              !surface temperature
    REAL(wp),       INTENT(out)             :: fc_tg              !skin temperature
    REAL(wp),       INTENT(out)             :: fc_qvs
    REAL(wp),       INTENT(out)             :: fc_Ch
    REAL(wp),       INTENT(out)             :: fc_Cq
    REAL(wp),       INTENT(out)             :: fc_Cm
    REAL(wp),       INTENT(out)             :: fc_ustar

    CHARACTER(len=*), PARAMETER :: routine = 'mo_ls_forcing:apply_ls_sfc_forcing'

    REAL(wp) :: int_weight
    INTEGER  :: i_startblk, i_endblk, jk, nlev, jb, jc
    INTEGER  :: n_curr, n_next, i_startidx, i_endidx

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !Find where in the array of LS forcing and nudging current time stands
    !and do linear interpolation in time
    IF (dt_forcing > 0._wp) THEN
      n_curr     = FLOOR(curr_sim_time/dt_forcing)+1
      n_next     = n_curr+1
      int_weight = curr_sim_time/dt_forcing-n_curr+1
      !ibd, exception for last time step
      IF (int_weight == 0._wp) THEN
        n_next = n_curr
      END IF

      ! Christopher Moseley: temporary catch
      IF (int_weight < 0 .OR. int_weight > 1) THEN
        WRITE(message_text,'(a,2(i0,", "),g0)') &
          'INTERPOLATION ERROR in LS sfc forcing: ', n_curr,n_next,int_weight
        CALL message(routine, message_text)
      END IF

      fc_sfc_lat_flx  = bnd_sfc_lat_flx (n_curr)  *(1.-int_weight) + bnd_sfc_lat_flx (n_next)  *int_weight
      fc_sfc_sens_flx = bnd_sfc_sens_flx(n_curr)  *(1.-int_weight) + bnd_sfc_sens_flx(n_next)  *int_weight
      fc_ts           = bnd_ts          (n_curr)  *(1.-int_weight) + bnd_ts          (n_next)  *int_weight
      fc_tg           = bnd_tg          (n_curr)  *(1.-int_weight) + bnd_tg          (n_next)  *int_weight
      fc_qvs          = bnd_qvs         (n_curr)  *(1.-int_weight) + bnd_qvs         (n_next)  *int_weight
      fc_Ch           = bnd_Ch          (n_curr)  *(1.-int_weight) + bnd_Ch          (n_next)  *int_weight
      fc_Cm           = bnd_Cm          (n_curr)  *(1.-int_weight) + bnd_Cm          (n_next)  *int_weight
      fc_Cq           = bnd_Cq          (n_curr)  *(1.-int_weight) + bnd_Cq          (n_next)  *int_weight
      fc_ustar        = bnd_ustar       (n_curr)  *(1.-int_weight) + bnd_ustar       (n_next)  *int_weight
    ELSE
      fc_sfc_lat_flx  = bnd_sfc_lat_flx(1)
      fc_sfc_sens_flx = bnd_sfc_sens_flx(1)
      fc_ts           = bnd_ts(1)
      fc_tg           = bnd_tg(1)
      fc_qvs          = bnd_qvs(1)
      fc_Ch           = bnd_Ch(1)
      fc_Cm           = bnd_Cm(1)
      fc_Cq           = bnd_Cq(1)
      fc_ustar        = bnd_ustar(1)
    END IF

    ! apply_ls_forcing has been called, tell SCM setup
    lscm_ls_forcing_ini = .TRUE.

  END SUBROUTINE apply_ls_sfc_forcing
!-----------------------------------------------------------------------------------
END MODULE mo_ls_forcing

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

! This module prepares aerosol for the use in radiation

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_nwp_aerosol

! ICON infrastructure
  USE mo_kind,                    ONLY: wp, rp
  USE mo_exception,               ONLY: finish, message, message_text
  USE mo_model_domain,            ONLY: t_patch
  USE mo_grid_config,             ONLY: n_dom, nroot, l_scm_mode
  USE mo_ext_data_types,          ONLY: t_external_data
  USE mo_nonhydro_types,          ONLY: t_nh_diag
  USE mo_nwp_phy_types,           ONLY: t_nwp_phy_diag
  USE mo_parallel_config,         ONLY: nproma
  USE mo_loopindices,             ONLY: get_indices_c
  USE mo_impl_constants,          ONLY: min_rlcell_int, SUCCESS, &
                                    &   iss, iorg, ibc, iso4, idu, n_camsaermr
  USE mo_impl_constants_grf,      ONLY: grf_bdywidth_c
  USE mo_physical_constants,      ONLY: rd, grav, cpd, rdv, o_m_rdv
  USE mo_math_constants,          ONLY: pi
  USE mo_reader_cams,             ONLY: t_cams_reader
  USE mo_interpolate_time,        ONLY: t_time_intp, t_time_intp_monthlyclim, t_time_intp_transient
  USE mo_io_units,                ONLY: filename_max
  USE mo_fortran_tools,           ONLY: init, set_acc_host_or_device, assert_acc_device_only
  USE mo_util_string,             ONLY: int2string, associate_keyword, t_keyword_list, with_keywords
! ICON configuration
  USE mo_atm_phy_nwp_config,      ONLY: atm_phy_nwp_config, i2daero_dust, i2daero_seas, i2daero_anthro,     &
    &                                   icpl_aero_conv
  USE mo_run_config,              ONLY: iqv
  USE mo_thdyn_functions,         ONLY: sat_pres_water
  USE mo_radiation_config,        ONLY: irad_aero, iRadAeroConstKinne, iRadAeroKinne, iRadAeroCAMSclim,     &
                                    &   iRadAeroCAMStd, iRadAeroVolc, iRadAeroKinneVolc, iRadAeroART,       &
                                    &   iRadAeroKinneVolcSP, iRadAeroKinneSP, iRadAeroTegen,                &
                                    &   iRadAeroExternal, cams_aero_filename
  USE mo_cpl_aerosol_microphys,   ONLY: calc_cdnc_from_cloud_num
  USE mo_nwp_tuning_config,       ONLY: tune_sc_eis
! External infrastruture
  USE mtime,                      ONLY: datetime, timedelta, newDatetime, newTimedelta,       &
                                    &   operator(+), deallocateTimedelta, deallocateDatetime
! Aerosol-specific
  USE mo_aerosol_util,            ONLY: aerdis
  USE mo_bc_aeropt_kinne,         ONLY: read_bc_aeropt_kinne, set_bc_aeropt_kinne
  USE mo_bc_aeropt_volc,          ONLY: read_bc_aeropt_volc, add_bc_aeropt_volc
  USE mo_bc_aeropt_splumes,       ONLY: add_bc_aeropt_splumes, cloud_num_scaling_factor
  USE mo_coupling_config,         ONLY: is_coupled_to_aero
  USE mo_bcs_time_interpolation,  ONLY: t_time_interpolation_weights,         &
    &                                   calculate_time_interpolation_weights
  USE mo_io_config,               ONLY: var_in_output
#ifdef __ICON_ART
  USE mo_aerosol_util,            ONLY: tegen_scal_factors
  USE mo_art_radiation_interface, ONLY: art_rad_aero_interface
#endif
  USE mo_scm_nml,                 ONLY: lat_scm

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_aerosol'

  PUBLIC :: nwp_aerosol_interface
  PUBLIC :: nwp_aerosol_cleanup
  PUBLIC :: nwp_aerosol_init
  PUBLIC :: nwp_aerosol_alloc
  PUBLIC :: nwp_aerosol_dealloc
  PUBLIC :: cams_reader
  PUBLIC :: cams_intp

  TYPE(t_cams_reader),      ALLOCATABLE, TARGET :: cams_reader(:)
  CLASS(t_time_intp),       ALLOCATABLE         :: cams_intp(:)

  ! Local memory for aerosol fields.
  ! These are only allocated for certain aerosol options in the scope of the radiation.
  ! If available, global memory fields from prm_diag are used instead
  REAL(wp), TARGET, ALLOCATABLE :: &
    &  locmem_od_lw(:,:,:,:), & !< LW optical thickness of aerosols
    &  locmem_od_sw(:,:,:,:), & !< SW aerosol optical thickness
    &  locmem_g_sw (:,:,:,:), & !< SW aerosol asymmetry factor
    &  locmem_ssa_sw(:,:,:,:)   !< SW aerosol single scattering albedo

CONTAINS

  SUBROUTINE nwp_aerosol_alloc

    ALLOCATE(cams_reader(n_dom))

    SELECT CASE (irad_aero)
    CASE (iRadAeroCAMSclim)
      ALLOCATE(t_time_intp_monthlyclim :: cams_intp(n_dom))
    CASE (iRadAeroCAMStd)
      ALLOCATE(t_time_intp_transient :: cams_intp(n_dom))
    CASE DEFAULT
    END SELECT

  END SUBROUTINE nwp_aerosol_alloc

  SUBROUTINE nwp_aerosol_dealloc
    INTEGER :: jg

    IF (ALLOCATED(cams_intp)) THEN
      DEALLOCATE(cams_intp)
    END IF

    IF (ALLOCATED(cams_reader)) THEN
      DO jg = 1, n_dom
        CALL cams_reader(jg)%deinit
      END DO

      DEALLOCATE(cams_reader)
    END IF

  END SUBROUTINE nwp_aerosol_dealloc

  !---------------------------------------------------------------------------------------
  !! This subroutine uploads CAMS aerosols and updates them once a day
  SUBROUTINE nwp_aerosol_init(mtime_datetime, p_patch)

    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                            !< Current datetime
    TYPE(t_patch), INTENT(in)           :: &
      &  p_patch
    ! Local variables
    INTEGER                             :: &
      &  jg                                        !< Domain index
    CHARACTER(LEN=filename_max)         :: &
      &  cams_aero_td_file                         !< CAMS file names

    CHARACTER(len=*), PARAMETER :: routine = modname // '::nwp_aerosol_init'



    jg     = p_patch%id

    cams_aero_td_file = generate_cams_filename(TRIM(cams_aero_filename), nroot, p_patch%level, p_patch%id)

    IF (ALLOCATED(cams_intp)) THEN
#if defined(_CRAYFTN)
      ! Cray Fortran cannot parse CALL cams_intp(jg)%init(...)
      SELECT TYPE (cams_intp)
      TYPE IS (t_time_intp_monthlyclim)
        CALL message  ('nwp_aerosol_init opening CAMS 3D climatology file: ', TRIM(cams_aero_td_file))
        CALL cams_reader(jg)%init(p_patch, TRIM(cams_aero_td_file))
      TYPE IS (t_time_intp_transient)
        CALL message  ('nwp_aerosol_init opening CAMS forecast file: ', TRIM(cams_aero_td_file))
        CALL cams_reader(jg)%init(p_patch, TRIM(cams_aero_td_file))
      CLASS DEFAULT
        CALL finish(routine, 'Internal error: unknown interpolator for CAMS aerosol.')
      END SELECT

      CALL init_interp(cams_intp(jg))
#else
      SELECT TYPE (cams_intp)
      TYPE IS (t_time_intp_monthlyclim)
        CALL message  ('nwp_aerosol_init opening CAMS 3D climatology file: ', TRIM(cams_aero_td_file))
        CALL cams_reader(jg)%init(p_patch, TRIM(cams_aero_td_file))
        CALL cams_intp(jg)%init(cams_reader(jg), mtime_datetime, '')

      TYPE IS (t_time_intp_transient)
        CALL message  ('nwp_aerosol_init opening CAMS forecast file: ', TRIM(cams_aero_td_file))
        CALL cams_reader(jg)%init(p_patch, TRIM(cams_aero_td_file))
        CALL cams_intp(jg)%init(cams_reader(jg), mtime_datetime, '')

      CLASS DEFAULT
        CALL finish(routine, 'Internal error: unknown interpolator for CAMS aerosol.')
      END SELECT
#endif
    END IF

    CONTAINS

      FUNCTION generate_cams_filename(filename_in, nroot, jlev, idom) RESULT(result_str)
        CHARACTER(filename_max)                 :: result_str
        CHARACTER(LEN=*), INTENT(in)            :: filename_in
        INTEGER,                     INTENT(in) :: nroot, jlev, idom
        ! Local variables
        TYPE (t_keyword_list), POINTER      :: keywords => NULL()

        CALL associate_keyword("<nroot>",    TRIM(int2string(nroot,"(i0)")),   keywords)
        CALL associate_keyword("<nroot0>",   TRIM(int2string(nroot,"(i2.2)")), keywords)
        CALL associate_keyword("<jlev>",     TRIM(int2string(jlev, "(i2.2)")), keywords)
        CALL associate_keyword("<idom>",     TRIM(int2string(idom, "(i2.2)")), keywords)

        result_str = TRIM(with_keywords(keywords, TRIM(filename_in)))
      END FUNCTION generate_cams_filename

#if defined(_CRAYFTN)
      SUBROUTINE init_interp (interp)
        CLASS(t_time_intp), INTENT(OUT) :: interp

        SELECT TYPE (interp)
        TYPE IS (t_time_intp_monthlyclim)
          CALL interp%init(cams_reader(jg), mtime_datetime, '')
        TYPE IS (t_time_intp_transient)
          CALL interp%init(cams_reader(jg), mtime_datetime, '')
        END SELECT
      END SUBROUTINE init_interp
#endif

  END SUBROUTINE nwp_aerosol_init

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_interface(mtime_datetime, pt_patch, ext_data, pt_diag, prm_diag,          &
    &                              zf, zh, dz, dt_rad,                                             &
    &                              inwp_radiation, nbands_lw, nbands_sw, wavenum1_sw, wavenum2_sw, &
    &                              zaeq1, zaeq2, zaeq3, zaeq4, zaeq5,                              &
    &                              od_lw, od_sw, ssa_sw, g_sw, lacc)
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//':nwp_aerosol_interface'

    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime          !< Current datetime
    TYPE(t_patch), TARGET, INTENT(in) :: &
      &  pt_patch                !< Grid/patch info
    TYPE(t_external_data), INTENT(inout) :: &
      &  ext_data                !< External data
    TYPE(t_nh_diag), TARGET, INTENT(inout) :: &
      &  pt_diag                 !< the diagnostic variables
    TYPE(t_nwp_phy_diag), INTENT(inout) :: &
      &  prm_diag                !< Physics diagnostics

    REAL(wp), INTENT(in) ::    &
      &  zf(:,:,:), zh(:,:,:), & !< model full/half layer height
      &  dz(:,:,:),            & !< Layer thickness
      &  dt_rad                  !< Radiation time step
    REAL(rp), POINTER, INTENT(in) :: &
      &  wavenum1_sw(:),       & !< Shortwave wavenumber lower band bounds
      &  wavenum2_sw(:)          !< Shortwave wavenumber upper band bounds
    REAL(wp), ALLOCATABLE, TARGET, INTENT(inout) :: &
      &  zaeq1(:,:,:),         & !< Tegen optical thicknesses       1: continental
      &  zaeq2(:,:,:),         & !< relative to 550 nm, including   2: maritime
      &  zaeq3(:,:,:),         & !< a vertical profile              3: desert
      &  zaeq4(:,:,:),         & !< for 5 different                 4: urban
      &  zaeq5(:,:,:)            !< aerosol species.                5: stratospheric background
    INTEGER, INTENT(in) ::     &
      &  inwp_radiation,       & !< Radiation scheme (1=rrtmg, 4=ecrad)
      &  nbands_lw, nbands_sw    !< Number of short and long wave bands
    REAL(wp), POINTER, INTENT(out) :: &
      &  od_lw(:,:,:,:),       & !< Longwave optical thickness
      &  od_sw(:,:,:,:),       & !< Shortwave optical thickness
      &  ssa_sw(:,:,:,:),      & !< Shortwave asymmetry factor
      &  g_sw(:,:,:,:)           !< Shortwave single scattering albedo

    LOGICAL, OPTIONAL, INTENT(in) :: lacc ! If true, use openacc
! Local variables
#ifdef __ICON_ART
    REAL(wp), ALLOCATABLE ::   &
      &  od_lw_art_vr(:,:,:),  & !< AOD LW (vertically reversed)
      &  od_sw_art_vr(:,:,:),  & !< AOD SW (vertically reversed)
      &  ssa_sw_art_vr(:,:,:), & !< SSA SW (vertically reversed)
      &  g_sw_art_vr(:,:,:)      !< Assymetry parameter SW (vertically reversed)
    INTEGER ::                 &
      &  jk_vr, jband            !< Loop indices
#endif
    REAL(wp) ::                &
      &  cloud_num_fac(nproma),& !< Scaling factor (simple plumes) for CDNC
      &  latitude(nproma),     & !< Geographical latitude
      &  time_weight             !< Weihting for temporal interpolation
    REAL(wp),  ALLOCATABLE ::  &
      &  cams(:,:,:,:)           !< CAMS climatology fields taken from external file
    INTEGER ::                 &
      &  jk, jc, jb, jt,       &
      &  jg,                   & !< Domain index
      &  rl_start, rl_end,     &
      &  i_startblk, i_endblk, &
      &  i_startidx, i_endidx, &
      &  istat,                & !< Error code
      &  imo1 , imo2             !< Month index (current and next month)
    LOGICAL :: lzacc

    jg     = pt_patch%id

    CALL set_acc_host_or_device(lzacc, lacc)

    SELECT CASE(irad_aero)
!---------------------------------------------------------------------------------------
! Tegen aerosol (+ART if chosen)
!---------------------------------------------------------------------------------------
      CASE(iRadAeroTegen, iRadAeroART)

        !$ACC DATA CREATE(latitude) IF(lzacc)

        ALLOCATE( zaeq1( nproma, pt_patch%nlev, pt_patch%nblks_c), &
          &       zaeq2( nproma, pt_patch%nlev, pt_patch%nblks_c), &
          &       zaeq3( nproma, pt_patch%nlev, pt_patch%nblks_c), &
          &       zaeq4( nproma, pt_patch%nlev, pt_patch%nblks_c), &
          &       zaeq5( nproma, pt_patch%nlev, pt_patch%nblks_c)  )
        !$ACC ENTER DATA CREATE(zaeq1, zaeq2, zaeq3, zaeq4, zaeq5) IF(lzacc)

        ! Outer two rows need dummy values as RRTM always starts at 1
        rl_start   = 1
        rl_end     = 2
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,jk,i_startidx,i_endidx)  ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR COLLAPSE(2)
          DO jk = 1, pt_patch%nlev
            DO jc = i_startidx,i_endidx
              zaeq1(jc,jk,jb) = 0._wp
              zaeq2(jc,jk,jb) = 0._wp
              zaeq3(jc,jk,jb) = 0._wp
              zaeq4(jc,jk,jb) = 0._wp
              zaeq5(jc,jk,jb) = 0._wp
            ENDDO
          ENDDO
          !$ACC END PARALLEL
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        ! Start at third row instead of fifth as two rows are needed by the reduced grid aggregation
        rl_start   = grf_bdywidth_c-1
        rl_end     = min_rlcell_int
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

        ! Calculate the weighting factor and month indices for temporal interpolation
        CALL get_time_intp_weights(mtime_datetime, imo1 , imo2, time_weight)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,latitude)  ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)

          !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
          !$ACC LOOP GANG VECTOR
          DO jc = i_startidx,i_endidx
            ! SCM: fixed latitude/longitude
            IF ( l_scm_mode ) THEN
              latitude(jc) = lat_scm * pi/180._wp
            ELSE
              latitude(jc) = pt_patch%cells%center(jc,jb)%lat
            ENDIF
          ENDDO
          !$ACC END PARALLEL

          CALL nwp_aerosol_tegen(i_startidx, i_endidx, pt_patch%nlev, pt_patch%nlevp1, prm_diag%k850(:,jb), &
            &                    pt_diag%temp(:,:,jb), pt_diag%pres(:,:,jb),                              &
            &                    pt_diag%pres_ifc(:,:,jb),  prm_diag%tot_cld(:,:,jb,iqv),                 &
            &                    ext_data%atm_td%aer_ss(:,jb,imo1),   ext_data%atm_td%aer_org(:,jb,imo1), &
            &                    ext_data%atm_td%aer_bc(:,jb,imo1),   ext_data%atm_td%aer_so4(:,jb,imo1), &
            &                    ext_data%atm_td%aer_dust(:,jb,imo1), ext_data%atm_td%aer_ss(:,jb,imo2),  &
            &                    ext_data%atm_td%aer_org(:,jb,imo2),  ext_data%atm_td%aer_bc(:,jb,imo2),  &
            &                    ext_data%atm_td%aer_so4(:,jb,imo2),  ext_data%atm_td%aer_dust(:,jb,imo2),&
            &                    prm_diag%pref_aerdis(:,jb),latitude(:), pt_diag%dpres_mc(:,:,jb), time_weight, &
            &                    prm_diag%aerosol(:,:,jb), prm_diag%aercl_ss(:,jb), prm_diag%aercl_or(:,jb), &
            &                    prm_diag%aercl_bc(:,jb), prm_diag%aercl_su(:,jb), prm_diag%aercl_du(:,jb), &
            &                    zaeq1(:,:,jb),zaeq2(:,:,jb),zaeq3(:,:,jb),zaeq4(:,:,jb),zaeq5(:,:,jb),lacc )

        ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

#if defined(__ECRAD) && defined(__ICON_ART)
        ! Replace Tegen selectively with ART aerosol
        IF ( irad_aero ==iRadAeroART ) THEN
#ifdef _OPENACC
          IF (lzacc) CALL finish(routine, "irad_aero==iRadAeroART is not ported to openACC.")
#endif
          IF (inwp_radiation == 4) THEN
            ! Allocations
            ALLOCATE(locmem_od_lw (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_lw), &
              &      locmem_od_sw (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw), &
              &      locmem_ssa_sw(nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw), &
              &      locmem_g_sw  (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw), &
              &      od_lw_art_vr (nproma,pt_patch%nlev,                 nbands_lw), &
              &      od_sw_art_vr (nproma,pt_patch%nlev,                 nbands_lw), &
              &      ssa_sw_art_vr(nproma,pt_patch%nlev,                 nbands_lw), &
              &      g_sw_art_vr  (nproma,pt_patch%nlev,                 nbands_lw), &
              &      STAT=istat)
            IF(istat /= SUCCESS) &
              &  CALL finish(routine, 'Allocation of od_lw, od_sw, ssa_sw, g_sw plus ART variants failed')
            od_lw  => locmem_od_lw
            od_sw  => locmem_od_sw
            ssa_sw => locmem_ssa_sw
            g_sw   => locmem_g_sw

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx,                              &
!$OMP            od_lw_art_vr,od_sw_art_vr,ssa_sw_art_vr,g_sw_art_vr, &
!$OMP            jc,jk,jk_vr,jband) ICON_OMP_DEFAULT_SCHEDULE
            DO jb = i_startblk,i_endblk
              CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
              IF (i_startidx>i_endidx) CYCLE

              CALL art_rad_aero_interface(zaeq1(:,:,jb),zaeq2(:,:,jb),       & !
                &                         zaeq3(:,:,jb),zaeq4(:,:,jb),       & !< Tegen aerosol
                &                         zaeq5(:,:,jb),                     & !
                &                         tegen_scal_factors%absorption,     & !
                &                         tegen_scal_factors%scattering,     & !< Tegen coefficients
                &                         tegen_scal_factors%asymmetry,      & !
                &                         pt_patch%id, jb, 1, pt_patch%nlev, & !< Indices domain, block, level
                &                         i_startidx, i_endidx,              & !< Indices nproma loop
                &                         nbands_lw,                         & !< Number of SW bands
                &                         nbands_sw,                         & !< Number of LW bands
                &                         od_lw_art_vr(:,:,:),               & !< OUT: Optical depth LW
                &                         od_sw_art_vr(:,:,:),               & !< OUT: Optical depth SW
                &                         ssa_sw_art_vr(:,:,:),              & !< OUT: SSA SW
                &                         g_sw_art_vr(:,:,:))                  !< OUT: Assymetry parameter SW


              DO jk = 1, pt_patch%nlev
                jk_vr = pt_patch%nlev+1-jk
! LONGWAVE
                DO jband = 1, nbands_lw
                  DO jc = i_startidx, i_endidx
                    od_lw(jc,jk,jb,jband) = od_lw_art_vr(jc,jk_vr,jband)
                  ENDDO !jc
                ENDDO !jband
! SHORTWAVE
                DO jband = 1, nbands_sw
                  DO jc = i_startidx, i_endidx
                    od_sw(jc,jk,jb,jband) = od_sw_art_vr(jc,jk_vr,jband)
                    ssa_sw(jc,jk,jb,jband) = ssa_sw_art_vr(jc,jk_vr,jband)
                    g_sw(jc,jk,jb,jband) = g_sw_art_vr(jc,jk_vr,jband)
                  ENDDO !jc
                ENDDO !jband
              ENDDO !jk
            ENDDO !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

            ! Deallocations
            DEALLOCATE(od_lw_art_vr, od_sw_art_vr, ssa_sw_art_vr, g_sw_art_vr, &
              &        STAT=istat)
            IF(istat /= SUCCESS) &
              &  CALL finish(routine, 'Deallocation of od_lw_art_vr, od_sw_art_vr, ssa_sw_art_vr, g_sw_art_vr failed')
          ENDIF
        ENDIF ! iRadAeroART
#endif

        !$ACC WAIT
        !$ACC END DATA

!---------------------------------------------------------------------------------------
! Kinne aerosol
!---------------------------------------------------------------------------------------
      CASE(iRadAeroConstKinne, iRadAeroKinne, iRadAeroVolc, iRadAeroKinneVolc, iRadAeroKinneVolcSP, iRadAeroKinneSP)
        !$ACC DATA CREATE(cloud_num_fac) IF(lzacc)

        rl_start   = grf_bdywidth_c-1
        rl_end     = min_rlcell_int
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

        ! Compatibility checks
#ifdef __ECRAD
        IF (inwp_radiation /= 4) THEN
          WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' only implemented for ecrad (inwp_radiation=4).'
          CALL finish(routine, message_text)
        ENDIF
#else
        WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' requires to compile with --enable-ecrad.'
        CALL finish(routine, message_text)
#endif

        ! Update Kinne aerosol from files once per day
        CALL nwp_aerosol_daily_update_kinne(mtime_datetime, pt_patch, dt_rad, inwp_radiation, &
          &                                 nbands_lw, nbands_sw)

        ! Allocations
        ALLOCATE(locmem_od_lw (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_lw)  , &
          &      locmem_od_sw (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      locmem_ssa_sw(nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      locmem_g_sw  (nproma,pt_patch%nlev,pt_patch%nblks_c,nbands_sw)  , &
          &      STAT=istat)
        IF(istat /= SUCCESS) &
          &  CALL finish(routine, 'Allocation of od_lw, od_sw, ssa_sw, g_sw failed')
        !$ACC ENTER DATA CREATE(locmem_od_lw, locmem_od_sw, locmem_ssa_sw, locmem_g_sw) IF(lzacc)
        od_lw  => locmem_od_lw
        od_sw  => locmem_od_sw
        ssa_sw => locmem_ssa_sw
        g_sw   => locmem_g_sw
        !$ACC ENTER DATA ATTACH(od_lw, od_sw, ssa_sw, g_sw) IF(lzacc)

        IF ( .NOT. ASSOCIATED(wavenum1_sw) .OR. .NOT. ASSOCIATED(wavenum2_sw) ) &
          &  CALL finish(routine, 'wavenum1 or wavenum2 not associated')
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jc,i_startidx,i_endidx,cloud_num_fac) ICON_OMP_DEFAULT_SCHEDULE
        DO jb = i_startblk,i_endblk
          CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
          IF (i_startidx>i_endidx) CYCLE

          CALL nwp_aerosol_kinne(mtime_datetime, zf(:,:,jb), zh(:,:,jb), dz(:,:,jb),   &
            &                    pt_patch%id, jb, i_startidx, i_endidx, pt_patch%nlev, &
            &                    nbands_lw, nbands_sw, wavenum1_sw(:), wavenum2_sw(:), &
            &                    od_lw(:,:,jb,:), od_sw(:,:,jb,:),                     &
            &                    ssa_sw(:,:,jb,:), g_sw(:,:,jb,:), cloud_num_fac(:),   &
            &                    lacc=lzacc)

          IF ( atm_phy_nwp_config(pt_patch%id)%scale_cdnc_mode /= 0 ) THEN
            ! scale the cdnc with the scaling factor:
            !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
            !$ACC LOOP GANG VECTOR
            DO jc = i_startidx,i_endidx
              prm_diag%cloud_num(jc,jb) = cloud_num_fac(jc) * ext_data%atm%cdnc(jc,jb)
            END DO
            !$ACC END PARALLEL
          ENDIF

          IF ( var_in_output(jg)%aod_550nm ) THEN
            CALL calc_aod550_kinne(i_startidx, i_endidx, pt_patch%nlev, od_sw(:,:,jb,10), &
              &                    prm_diag%aod_550nm(:,jb), lacc=lzacc)
          END IF

          IF (atm_phy_nwp_config(pt_patch%id)%icpl_aero_gscp == 3 .OR. icpl_aero_conv == 1) THEN
            CALL calc_cdnc_from_cloud_num(i_startidx, i_endidx, pt_patch%nlev, pt_diag%pres_sfc(:,jb),  &
                           &   pt_diag%pres(:,:,jb), prm_diag%acdnc(:,:,jb), prm_diag%cloud_num(:,jb),  &
                               prm_diag%k_inversion(:,jb), prm_diag%conv_eis(:,jb), tune_sc_eis, lacc=lzacc)
          ENDIF

        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        !$ACC WAIT
        !$ACC END DATA
      ! CAMS climatology/forecasted aerosols
      CASE(iRadAeroCAMSclim,iRadAeroCAMStd)

#ifdef _OPENACC
        IF (lzacc) THEN
          WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' not ported to OpenACC.'
          CALL finish(routine, message_text)
        ENDIF
#endif

        rl_start   = grf_bdywidth_c+1
        rl_end     = min_rlcell_int
        i_startblk = pt_patch%cells%start_block(rl_start)
        i_endblk   = pt_patch%cells%end_block(rl_end)

        ! Compatibility checks
#ifdef __ECRAD
        IF (inwp_radiation /= 4) THEN
          WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' only implemented for ecrad (inwp_radiation=4).'
          CALL finish(routine, message_text)
        ENDIF
#else
        WRITE(message_text,'(a,i2,a)') 'irad_aero = ', irad_aero,' requires to compile with --enable-ecrad.'
        CALL finish(routine, message_text)
#endif

        CALL nwp_aerosol_update_cams(mtime_datetime, pt_patch%id, cams)

        DO jt=1, n_camsaermr
!$OMP PARALLEL
!$OMP DO PRIVATE(jb,i_startidx,i_endidx) ICON_OMP_DEFAULT_SCHEDULE
          DO jb = i_startblk,i_endblk
            CALL get_indices_c(pt_patch,jb,i_startblk,i_endblk,i_startidx,i_endidx,rl_start,rl_end)
            IF (i_startidx>i_endidx) CYCLE
            pt_diag%camsaermr(:,:,jb,jt) = 0._wp

            CALL vinterp_cams(i_startidx, i_endidx, pt_diag%pres(:,:,jb),pt_diag%pres_ifc(:,pt_patch%nlev+1,jb),&
              &               cams_pres_in=cams(:,:,jb,n_camsaermr+1), cams=cams(:,:,jb,jt), &
              &               nlev=pt_patch%nlev, camsaermr=pt_diag%camsaermr(:,:,jb,jt))
          END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
        END DO ! jt aerosol loop

      CASE(iRadAeroExternal)

        IF ( .NOT. ASSOCIATED(prm_diag%od_lw) .OR. .NOT. ASSOCIATED(prm_diag%od_sw) .OR. &
          &  .NOT. ASSOCIATED(prm_diag%ssa_sw) .OR. .NOT. ASSOCIATED(prm_diag%g_sw) ) &
          &  CALL finish(routine, 'od_lw, od_sw, ssa_sw or g_sw not associated!')

        ! Externally specified aerosol is stored in prm_diag
        od_lw  => prm_diag%od_lw
        od_sw  => prm_diag%od_sw
        ssa_sw => prm_diag%ssa_sw
        g_sw   => prm_diag%g_sw
        !$ACC ENTER DATA ATTACH(od_lw, od_sw, ssa_sw, g_sw) IF(lzacc)

      CASE DEFAULT
        ! Currently continue as not all cases are ported to nwp_aerosol_interface yet
    END SELECT

  END SUBROUTINE nwp_aerosol_interface


  SUBROUTINE calc_aod550_kinne(i_startidx, i_endidx, nlev, od_sw_band10, aod_550nm, lacc)

    INTEGER,  INTENT(in)                :: &
      &  i_startidx, i_endidx, nlev             !< loop start and end indices (nproma, vertical)
    REAL(wp), INTENT(in)                :: &
      &  od_sw_band10(:,:)                      !< Shortwave optical thickness 10th band range (442 - 625nm)
    REAL(wp), INTENT(inout)             :: &
      &  aod_550nm(:)                           !< cloud droplet number concentration
    LOGICAL,  INTENT(in), OPTIONAL      :: &
      &  lacc                                   !< If true, use openacc

    INTEGER :: jk, jc

    ! Output AOD in SW band 550nm
    ! 10th band range (442 - 625nm) is used to output aod_550 nm
    ! due to a lack of spectrally resolved information for a particular wavelength

    CALL assert_acc_device_only("calc_aod550_kinne", lacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR
    DO jc = i_startidx, i_endidx
      aod_550nm(jc) = 0.0_wp
    ENDDO

    !$ACC LOOP SEQ
    DO jk = 1, nlev
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO jc = i_startidx, i_endidx
        aod_550nm(jc) = aod_550nm(jc) + od_sw_band10(jc,jk)
      ENDDO !jc
    ENDDO !jk
    !$ACC END PARALLEL

  END SUBROUTINE calc_aod550_kinne

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_daily_update_kinne(mtime_datetime, pt_patch, dt_rad, inwp_radiation, nbands_lw, nbands_sw)
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                    !< Current datetime
    TYPE(t_patch), TARGET, INTENT(in) :: &
      &  pt_patch                          !< Grid/patch info
    REAL(wp), INTENT(in) ::              &
      &  dt_rad                            !< Radiation time step
    INTEGER, INTENT(in) ::               &
      &  inwp_radiation,                 & !< Radiation scheme (1=rrtmg, 4=ecrad)
      &  nbands_lw, nbands_sw              !< Number of short and long wave bands
    ! Local variables
    TYPE(datetime), POINTER ::           &
      &  prev_radtime                      !< Datetime of previous radiation time step
    TYPE(timedelta), POINTER ::          &
      &  td_dt_rad                         !< Radiation time step

    td_dt_rad => newTimedelta('-',0,0,0,0,0, second=NINT(dt_rad), ms=0)
    prev_radtime => newDatetime(mtime_datetime + td_dt_rad)

    IF (prev_radtime%date%day /= mtime_datetime%date%day) THEN
      IF (inwp_radiation == 4) THEN
        IF (ANY(irad_aero == [iRadAeroKinne, iRadAeroKinneVolc])) &
            & CALL read_bc_aeropt_kinne(mtime_datetime, pt_patch, .TRUE., nbands_lw, nbands_sw, &
                                        opt_from_coupler = is_coupled_to_aero())
        IF (ANY(irad_aero == [iRadAeroVolc, iRadAeroKinneVolc, iRadAeroKinneVolcSP])) &
            & CALL read_bc_aeropt_volc(mtime_datetime, nbands_lw, nbands_sw)
      ENDIF
    ENDIF

    CALL deallocateTimedelta(td_dt_rad)
    CALL deallocateDatetime(prev_radtime)

  END SUBROUTINE nwp_aerosol_daily_update_kinne

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_kinne(mtime_datetime, zf, zh, dz, jg, jb, i_startidx, i_endidx, nlev, &
    &                          nbands_lw, nbands_sw, wavenum1_sw, wavenum2_sw,     &
    &                          od_lw, od_sw, ssa_sw, g_sw, cloud_num_fac, lacc)
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                      !< Current datetime
    REAL(wp), INTENT(in) ::                &
      &  zf(:,:), zh(:,:), dz(:,:)           !< model full/half layer height, layer thickness
    REAL(rp), INTENT(in) ::                &
      &  wavenum1_sw(:),                   & !< Shortwave wavenumber lower band bounds
      &  wavenum2_sw(:)                      !< Shortwave wavenumber upper band bounds
    INTEGER, INTENT(in) ::                 &
      &  jg, jb,                           & !< Domain and block index
      &  i_startidx,                       & !< Loop bound
      &  i_endidx,                         & !< Loop bound
      &  nlev,                             & !< Number of vertical levels
      &  nbands_lw, nbands_sw                !< Number of short and long wave bands
    REAL(wp), INTENT(out) ::               &
      &  od_lw(:,:,:), od_sw(:,:,:),       & !< LW/SW optical thickness
      &  ssa_sw(:,:,:), g_sw(:,:,:),       & !< SW asymmetry factor, SW single scattering albedo
      &  cloud_num_fac(:)                    !< Scaling factor for Cloud Droplet Number Concentration;
                                             !< cdnc is scaled if scale_cdnc_mode /= 0
    ! Local variables
    REAL(wp) ::                            &
      &  od_lw_vr (nproma,nlev,nbands_lw), & !< LW optical thickness of aerosols    (vertically reversed)
      &  od_sw_vr (nproma,nlev,nbands_sw), & !< SW aerosol optical thickness        (vertically reversed)
      &  g_sw_vr  (nproma,nlev,nbands_sw), & !< SW aerosol asymmetry factor         (vertically reversed)
      &  ssa_sw_vr(nproma,nlev,nbands_sw)    !< SW aerosol single scattering albedo (vertically reversed)
    TYPE(datetime), POINTER :: &
      &  mtime_local,                      & !< local time variable to get x_cdnc
      &  mtime_ref                           !< reference time in 2005 to get x_cdnc_ref

    INTEGER ::                             &
      &  jk, jc, jwl                         !< Loop indices
    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//':nwp_aerosol_kinne'
    LOGICAL, INTENT(in), OPTIONAL :: lacc   !< If true, use openacc
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(od_lw_vr, od_sw_vr, g_sw_vr, ssa_sw_vr) IF(lzacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO jwl = 1, nbands_lw
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = 1, nproma
          od_lw_vr(jc,jk,jwl)  = 0.0_wp
        END DO
      END DO
    END DO

    !$ACC LOOP SEQ
    DO jwl = 1, nbands_sw
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO jk = 1, nlev
        DO jc = 1, nproma
          od_sw_vr(jc,jk,jwl)  = 0.0_wp
          ssa_sw_vr(jc,jk,jwl) = 1.0_wp
          g_sw_vr(jc,jk,jwl)   = 0.0_wp
        END DO
      END DO
    END DO
    !$ACC END PARALLEL


    ! Tropospheric Kinne aerosol
    IF (ANY( irad_aero == (/iRadAeroConstKinne,iRadAeroKinne,iRadAeroKinneVolc, &
      &                     iRadAeroKinneVolcSP,iRadAeroKinneSP/) )) THEN
      CALL set_bc_aeropt_kinne(mtime_datetime, jg, i_startidx, i_endidx, nproma, nlev, jb, &
        &                      nbands_sw, nbands_lw, zf(:,:), dz(:,:),            &
        &                      od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),                 &
        &                      g_sw_vr (:,:,:), od_lw_vr(:,:,:),                  &
        &                      lacc=lzacc, opt_from_coupler = is_coupled_to_aero())

    ENDIF

    ! Volcanic stratospheric aerosols
    IF (ANY( irad_aero == (/iRadAeroVolc,iRadAeroKinneVolc,iRadAeroKinneVolcSP/) )) THEN
      CALL add_bc_aeropt_volc(mtime_datetime, jg, i_startidx, i_endidx, nproma, nlev, jb, &
        &                           nbands_sw, nbands_lw, zf(:,:), dz(:,:),            &
        &                           od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),                 &
        &                           g_sw_vr (:,:,:), od_lw_vr(:,:,:), lacc=lzacc       )
    END IF

    ! Simple plumes
    IF ( ANY( irad_aero == (/iRadAeroKinneVolcSP,iRadAeroKinneSP/) ) ) THEN
      ! Add the anthropogenic aerosol optical properties on top of the background aerosol
      CALL add_bc_aeropt_splumes(jg, i_startidx, i_endidx, nproma, nlev, jb,  &
        &                        nbands_sw, mtime_datetime,          &
        &                        zf(:,:), dz(:,:), zh(:,nlev+1),     & ! in
        &                        wavenum1_sw(:), wavenum2_sw(:),     & ! in
        &                        od_sw_vr(:,:,:), ssa_sw_vr(:,:,:),  & ! inout
        &                        g_sw_vr (:,:,:),                    & ! inout
        &                        lacc=lzacc                          )
    END IF

    IF (atm_phy_nwp_config(jg)%scale_cdnc_mode /= 0) THEN
      ! set the reference year as 2005:
      mtime_ref => newDatetime(mtime_datetime)
      mtime_ref%date%year = 2005

      ! Set the time variable to calculate x_cdnc to year 1850
      ! when constant cdnc scaling to year 1850 is required (e.g. in picontrol experiment type)
      mtime_local => newDatetime(mtime_datetime)
      IF (atm_phy_nwp_config(jg)%scale_cdnc_mode == 2) mtime_local%date%year = 1850

      ! Get the cloud droplet number scaling factor:
      CALL cloud_num_scaling_factor(jg, i_startidx, i_endidx, nproma, nlev, jb,  &
        &                        mtime_ref, mtime_local,             & ! in
        &                        zf(:,:), dz(:,:), zh(:,nlev+1),     & ! in
        &                        cloud_num_fac(:),                   & ! out
        &                        lacc=lzacc                          )

      CALL deallocateDatetime(mtime_ref)
      CALL deallocateDatetime(mtime_local)

    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    ! Vertically reverse the fields:
    !$ACC LOOP SEQ
    DO jk = 1, nlev
      !$ACC LOOP SEQ
      DO jwl = 1, nbands_lw
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nproma
          od_lw (jc,jk,jwl) = od_lw_vr (jc,nlev-jk+1,jwl)
        END DO
      END DO

      !$ACC LOOP SEQ
      DO jwl = 1, nbands_sw
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO jc = 1, nproma
          od_sw (jc,jk,jwl) = od_sw_vr (jc,nlev-jk+1,jwl)
          ssa_sw(jc,jk,jwl) = ssa_sw_vr(jc,nlev-jk+1,jwl)
          g_sw  (jc,jk,jwl) = g_sw_vr  (jc,nlev-jk+1,jwl)
        END DO
      END DO
    END DO
    !$ACC END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE nwp_aerosol_kinne


  !---------------------------------------------------------------------------------------
  !! This subroutine uploads CAMS aerosols mixing ratios 3D climatology and updates them
  SUBROUTINE nwp_aerosol_update_cams(mtime_datetime, jg, cams)

    TYPE(datetime), POINTER, INTENT(in)  :: &
      &  mtime_datetime                    !< Current datetime
    REAL(wp), ALLOCATABLE, INTENT(inout) :: &
      &  cams(:,:,:,:)                     !< CAMS fields taken from external file
    INTEGER, INTENT(in)                  :: &
      &  jg                                !< Domain index
    ! Local variables
    REAL(wp), ALLOCATABLE                :: &
      &  cams_dat(:,:,:,:)

    ALLOCATE(cams( nproma, cams_reader(jg)%nlev_cams, cams_reader(jg)%p_patch%nblks_c, n_camsaermr+1 ))
    cams(:,:,:,:) = 0.0_wp

    CALL cams_intp(jg)%intp(mtime_datetime, cams_dat, lacc=.FALSE.)

    cams(:,:,:,:) = cams_dat(:,:,:,:)

  END SUBROUTINE nwp_aerosol_update_cams

  !---------------------------------------------------------------------------------------------------------
  ! Vertically interpolate CAMS climatology onto ICON model levels. CAMS aerosol is provided
  ! as mixing ratios, and output is mixing ratio on ICON levels
  SUBROUTINE vinterp_cams(i_startidx, i_endidx, pres,sfcpres, cams_pres_in, cams, nlev, camsaermr )

    REAL(wp),      INTENT(in)      :: &
      &  cams(:,:),                      & !< CAMS fields taken from external file; mixing ratio [kg/kg]
      &  pres(:,:),                      & !< ICON diagnosed pressure, model layer centers [Pa]
      &  sfcpres(:),                     & !< ICON surface pressure, lowest model layer interface [Pa]
      &  cams_pres_in(:,:)                 !< CAMS climatology pressure taken from external file,
                                           !< on model layer centers [Pa]
    REAL(wp),      INTENT(inout)   :: &
      &  camsaermr(:,:)                    !< CAMS aerosols mixing ratios [kg/kg]
    INTEGER,       INTENT(in)      :: &
      &  nlev,                           & !< ICON number of vertical levels
      &  i_startidx,                     & !< Loop indices
      &  i_endidx                          !< Loop indices

    ! local variables
    INTEGER                        :: &
      &  jc, jk, jk1,                 &    !< Loop indices
      &  nk1                               !< number of vertical levels in original CAMS climatology data
    REAL(wp)                       :: rescale, sigmasfc

    CHARACTER(len=*), PARAMETER    :: &
      &  routine = modname//':vinterp_cams'

    nk1 = size(cams_pres_in,2)
    ! To get CAMS surface pressure, use knowledge of appropriate 'B' hybrid coefficient of second
    ! lowest model layer interface to reconstruct surface pressure from pressure on lowest model
    ! layer centre pres_mc(nlev), which is known.
    ! Coefficients can be found here: https://confluence.ecmwf.int/display/UDOC/Model+level+definitions
    !
    !nlev-1: layer interface ~~~~~~~~~~~~~~~~~ A=0, B=0.99XXXX (dependent on number of model layers)
    !
    !nlev: layer centre      - - - - - - - - - pres_mc(nlev)=0.5*(pres_ifc(nlev)+pres_ifc(nlev-1))
    !
    !nlev: layer interface   ~~~~~~~~~~~~~~~~~ A=0, B=1 (SURFACE)
    !
    IF (nk1==60) THEN
       sigmasfc=1._wp/(0.5_wp*1.997630_wp)
    ELSEIF (nk1==21) THEN
       sigmasfc=1._wp/(0.5_wp*1.992281_wp)
    ELSEIF (nk1==137) THEN
       sigmasfc=1._wp/(0.5_wp*1.997630_wp)
    ELSE
       WRITE(message_text,'(a,i2,a)') 'A CAMS climatology with ',nk1,' number of levels is currently not supported'
       CALL finish(routine, message_text)
    ENDIF

    ! loop over grid points
    DO jc = i_startidx, i_endidx
      rescale=sfcpres(jc)/(cams_pres_in(jc,nk1)*sigmasfc)
      ! loop over target ICON levels
      DO jk = 1, nlev
        IF (pres(jc,jk) .gt. cams_pres_in(jc,nk1)*rescale) THEN
           ! Current ICON layer is below CAMS surface: extrapolate
           camsaermr(jc,jk) = cams(jc,nk1)/rescale
        ENDIF
        DO jk1 = 1, nk1-1 ! loop over original CAMS levels
           IF (pres(jc,jk) .gt. cams_pres_in(jc,jk1)*rescale .and. pres(jc,jk).le. cams_pres_in(jc,jk1+1)*rescale) THEN
              ! Current ICON layer between two CAMS layers: linear interpolation
              camsaermr(jc,jk) = (cams(jc,jk1)+(pres(jc,jk)-cams_pres_in(jc,jk1)*rescale)* &
                   & (cams(jc,jk1+1)-cams(jc,jk1))/(cams_pres_in(jc,jk1+1)*rescale-cams_pres_in(jc,jk1)*rescale))/rescale
           ENDIF
        ENDDO
        IF (pres(jc,jk) .le. cams_pres_in(jc,1)*rescale) THEN
           ! Current ICON layer is above CAMS top: extrapolate
           camsaermr(jc,jk) = cams(jc,1)/rescale
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE vinterp_cams


  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_tegen ( istart, iend, nlev, nlevp1, k850, temp, pres, pres_ifc, qv,             &
    &                            aer_ss_mo1, aer_org_mo1, aer_bc_mo1, aer_so4_mo1, aer_dust_mo1,         &
    &                            aer_ss_mo2, aer_org_mo2, aer_bc_mo2, aer_so4_mo2, aer_dust_mo2,         &
    &                            pref_aerdis,latitude, dpres_mc, time_weight,                            &
    &                            aerosol, aercl_ss, aercl_or, aercl_bc, aercl_su, aercl_du, &
    &                            zaeq1,zaeq2,zaeq3,zaeq4,zaeq5,lacc )

    INTEGER,  INTENT(in)                :: &
      &  nlev, nlevp1,                     & !< Number of vertical full/half levels
      &  istart, iend,                     & !< Start and end index of jc loop
      &  k850(:)                             !< Index of 850 hPa layer
    REAL(wp), INTENT(in)                :: &
      &  temp(:,:), pres(:,:),             & !< temperature and pressure at full level
      &  pres_ifc(:,:), qv(:,:),           & !< pressure at half level, specific humidity
      &  aer_ss_mo1(:), aer_org_mo1(:),    & !< Month 1 climatology from extpar file (sea salt, organic)
      &  aer_bc_mo1(:), aer_so4_mo1(:),    & !< Month 1 climatology from extpar file (blck carbon, sulphate)
      &  aer_dust_mo1(:),                  & !< Month 1 climatology from extpar file (dust)
      &  aer_ss_mo2(:), aer_org_mo2(:),    & !< Month 2 climatology from extpar file (sea salt, organic)
      &  aer_bc_mo2(:), aer_so4_mo2(:),    & !< Month 2 climatology from extpar file (blck carbon, sulphate)
      &  aer_dust_mo2(:),                  & !< Month 2 climatology from extpar file (dust)
      &  pref_aerdis(:),                   & !< Reference pressure for vertical distribution of aerosol
      &  latitude(:),                      & !< geographical latitude
      &  dpres_mc(:,:),                    & !< pressure thickness
      &  time_weight                         !< Temporal weighting factor
    REAL(wp), TARGET, INTENT(inout)     :: &
      &  aerosol(:,:),                     & !< Aerosol field incl. temporal interpolation
      &  aercl_ss(:), aercl_or(:),         & !< Climatological fields for relaxation (i2daero_seas/i2daero_anthro > 0)
      &  aercl_bc(:), aercl_su(:),         & !< Climatological fields for relaxation (i2daero_anthro > 0)
      &  aercl_du(:)                         !< Climatological fields for relaxation (i2daero_dust > 0)
    REAL(wp), INTENT(inout)             :: &
      &  zaeq1(:,:), zaeq2(:,:),           & !< organics, sea salt
      &  zaeq3(:,:), zaeq4(:,:),           & !< dust, black carbon
      &  zaeq5(:,:)                          !< sulphate (incl. stratospheric background zstbga)
    LOGICAL, INTENT(in), OPTIONAL       :: &
      &  lacc                                !< If true, use openacc
! Local variables
    CHARACTER(len=*), PARAMETER         :: &
      &  routine = modname//':nwp_aerosol_tegen'
    REAL(wp)                            :: &
      &  zsign (nproma,nlevp1),            &
      &  zvdaes(nproma,nlevp1),            &
      &  zvdael(nproma,nlevp1),            &
      &  zvdaeu(nproma,nlevp1),            &
      &  zvdaed(nproma,nlevp1),            &
      &  zaeqdo (nproma), zaeqdn,          &
      &  zaequo (nproma), zaequn,          &
      &  zaeqlo (nproma), zaeqln,          &
      &  zaeqsuo(nproma), zaeqsun,         &
      &  zaeqso (nproma), zaeqsn,          &
      &  zptrop (nproma), zdtdz(nproma),   &
      &  zlatfac(nproma), zstrfac,         &
      &  zpblfac, zslatq, tunefac_pbl, rh, humidity_fac
    REAL(wp), POINTER                   :: &
      &  climaero_ss(:),                   &
      &  climaero_org(:),                  &
      &  climaero_bc(:),                   &
      &  climaero_so4(:),                  &
      &  climaero_du(:)
    REAL(wp), PARAMETER                 :: &
      & ztrbga = 0.03_wp  / (101325.0_wp - 19330.0_wp), &
      ! original value for zstbga of 0.045 is much higher than recently published climatologies
      & zstbga = 0.015_wp  / 19330.0_wp
    INTEGER                             :: &
      &  jc,jk                               !< Loop indices
    LOGICAL                             :: &
      &  lzacc                               !< non-optional version of lacc

    ! increase aerosol enhancement in stable PBLs with prognostic aersosol
    tunefac_pbl = MERGE(1._wp, 2._wp, i2daero_anthro == 0)

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(climaero_ss, climaero_org, climaero_bc, climaero_so4, climaero_du) &
    !$ACC   CREATE(zsign, zvdaes, zvdael, zvdaeu, zvdaed) &
    !$ACC   CREATE(zaeqdo, zaequo, zaeqlo, zaeqsuo, zaeqso, zptrop) &
    !$ACC   CREATE(zdtdz, zlatfac) IF(lzacc)

    IF (i2daero_dust == 1) THEN
      climaero_du => aercl_du(:)
    ELSE
      climaero_du => aerosol(:,idu)
    ENDIF

    IF (i2daero_seas == 1) THEN
      climaero_ss => aercl_ss(:)
    ELSE
      climaero_ss => aerosol(:,iss)
    ENDIF

    IF (i2daero_anthro == 1) THEN
      climaero_org => aercl_or(:)
      climaero_bc  => aercl_bc(:)
      climaero_so4 => aercl_su(:)
    ELSE
      climaero_org => aerosol(:,iorg)
      climaero_bc  => aerosol(:,ibc)
      climaero_so4 => aerosol(:,iso4)
    ENDIF

    !$ACC UPDATE DEVICE(climaero_ss, climaero_org, climaero_bc, climaero_so4, climaero_du)

!DIR$ IVDEP
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    DO jc = istart, iend
      climaero_ss (jc) = aer_ss_mo1  (jc) + ( aer_ss_mo2  (jc) - aer_ss_mo1  (jc) ) * time_weight
      climaero_org(jc) = aer_org_mo1 (jc) + ( aer_org_mo2 (jc) - aer_org_mo1 (jc) ) * time_weight
      climaero_bc (jc) = aer_bc_mo1  (jc) + ( aer_bc_mo2  (jc) - aer_bc_mo1  (jc) ) * time_weight
      climaero_so4(jc) = aer_so4_mo1 (jc) + ( aer_so4_mo2 (jc) - aer_so4_mo1 (jc) ) * time_weight
      climaero_du (jc) = aer_dust_mo1(jc) + ( aer_dust_mo2(jc) - aer_dust_mo1(jc) ) * time_weight
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR COLLAPSE(2)
    DO jk = 2, nlevp1
      DO jc = istart, iend
        zsign(jc,jk) = pres_ifc(jc,jk) / MAX(pref_aerdis(jc),0.95_wp*pres_ifc(jc,nlevp1))
      ENDDO
    ENDDO
    !$ACC END PARALLEL

    CALL aerdis ( &
      & kbdim  = nproma,      & !in
      & jcs    = istart,      & !in
      & jce    = iend,        & !in
      & klevp1 = nlevp1,      & !in
      & petah  = zsign(1,1),  & !in
      & pvdaes = zvdaes(1,1), & !out
      & pvdael = zvdael(1,1), & !out
      & pvdaeu = zvdaeu(1,1), & !out
      & pvdaed = zvdaed(1,1), & !out
      & lacc = lzacc)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR PRIVATE(jk, zslatq)
    DO jc = istart, iend
      ! top level
      zaeqso (jc) = zvdaes(jc,1) * aerosol(jc,iss)
      zaeqlo (jc) = zvdael(jc,1) * aerosol(jc,iorg)
      zaeqsuo(jc) = zvdael(jc,1) * aerosol(jc,iso4)
      zaequo (jc) = zvdaeu(jc,1) * aerosol(jc,ibc)
      zaeqdo (jc) = zvdaed(jc,1) * aerosol(jc,idu)

      ! tropopause pressure and PBL stability
      jk          = k850(jc)
      zslatq      = SIN(latitude(jc))**2
      zptrop(jc)  = 1.e4_wp + 2.e4_wp*zslatq ! 100 hPa at the equator, 300 hPa at the poles
      zdtdz(jc)   = (temp(jc,jk)-temp(jc,nlev-1))/(-rd/grav*                     &
        &           (temp(jc,jk)+temp(jc,nlev-1))*(pres(jc,jk)-pres(jc,nlev-1))/ &
        &           (pres(jc,jk)+pres(jc,nlev-1)))
      ! latitude-dependence of tropospheric background
      zlatfac(jc) = MAX(0.1_wp, 1._wp-MERGE(zslatq**3, zslatq, latitude(jc) > 0._wp))
    ENDDO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP SEQ
    DO jk = 1,nlev
      !$ACC LOOP GANG VECTOR PRIVATE(zaeqsn, zaeqln, zaeqsun, zaequn, zaeqdn, zstrfac, zpblfac, rh, humidity_fac)
      DO jc = istart, iend
        zaeqsn  = zvdaes(jc,jk+1) * aerosol(jc,iss)
        zaeqln  = zvdael(jc,jk+1) * aerosol(jc,iorg)
        zaeqsun = zvdael(jc,jk+1) * aerosol(jc,iso4)
        zaequn  = zvdaeu(jc,jk+1) * aerosol(jc,ibc)
        zaeqdn  = zvdaed(jc,jk+1) * aerosol(jc,idu)

        ! stratosphere factor: 1 in stratosphere, 0 in troposphere, width of transition zone 0.1*p_TP
        zstrfac = MIN(1._wp,MAX(0._wp,10._wp*(zptrop(jc)-pres(jc,jk))/zptrop(jc)))
        ! PBL stability factor; enhance organic, sulfate and black carbon aerosol for stable stratification;
        ! account for particle growth in nearly saturated air
        rh = qv(jc,jk)*pres(jc,jk)/((rdv+o_m_rdv*qv(jc,jk))*sat_pres_water(temp(jc,jk)))
        humidity_fac = MERGE(0._wp, 20._wp*MAX(0._wp,rh-0.8_wp), i2daero_anthro == 0)
        zpblfac = 1._wp + tunefac_pbl*MIN(1.5_wp,1.e2_wp*MAX(0._wp, zdtdz(jc) + grav/cpd)) + humidity_fac

        zaeq1(jc,jk) = (1._wp-zstrfac)*MAX(zpblfac*(zaeqln-zaeqlo(jc)), ztrbga*zlatfac(jc)*dpres_mc(jc,jk))
        zaeq2(jc,jk) = (1._wp-zstrfac)*(zaeqsn-zaeqso(jc))
        zaeq3(jc,jk) = (1._wp-zstrfac)*(zaeqdn-zaeqdo(jc))
        zaeq4(jc,jk) = (1._wp-zstrfac)*zpblfac*(zaequn-zaequo(jc))
        zaeq5(jc,jk) = (1._wp-zstrfac)*zpblfac*(zaeqsun-zaeqsuo(jc)) + zstrfac*zstbga*dpres_mc(jc,jk)

        zaeqso(jc)  = zaeqsn
        zaeqlo(jc)  = zaeqln
        zaeqsuo(jc) = zaeqsun
        zaequo(jc)  = zaequn
        zaeqdo(jc)  = zaeqdn

      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC END DATA

  END SUBROUTINE nwp_aerosol_tegen

  !---------------------------------------------------------------------------------------

  SUBROUTINE get_time_intp_weights(mtime_datetime, imo1 , imo2, time_weight)
    TYPE(datetime), POINTER, INTENT(in) :: &
      &  mtime_datetime                      !< Current datetime
    INTEGER, INTENT(out)                :: &
      &  imo1 , imo2                         !< Month indices for temporal interpolation
    REAL(wp), INTENT(out)               :: &
      &  time_weight
    ! Local variables
    TYPE(datetime), POINTER             :: &
      &  current_time_hours
    TYPE(t_time_interpolation_weights)  :: &
      &  current_time_interpolation_weights

    current_time_hours => newDatetime(mtime_datetime)
    current_time_hours%time%minute = 0
    current_time_hours%time%second = 0
    current_time_hours%time%ms = 0

    current_time_interpolation_weights = calculate_time_interpolation_weights(current_time_hours)

    imo1        = current_time_interpolation_weights%month1
    imo2        = current_time_interpolation_weights%month2
    time_weight = current_time_interpolation_weights%weight2

    CALL deallocateDatetime(current_time_hours)

  END SUBROUTINE get_time_intp_weights

  !---------------------------------------------------------------------------------------
  SUBROUTINE nwp_aerosol_cleanup(zaeq1, zaeq2, zaeq3, zaeq4, zaeq5, od_lw, od_sw, ssa_sw, g_sw, lacc)

    CHARACTER(len=*), PARAMETER :: &
      &  routine = modname//':nwp_aerosol_cleanup'

    REAL(wp), ALLOCATABLE, INTENT(inout) :: &
      &  zaeq1(:,:,:),         & !< Tegen optical thicknesses       1: continental
      &  zaeq2(:,:,:),         & !< relative to 550 nm, including   2: maritime
      &  zaeq3(:,:,:),         & !< a vertical profile              3: desert
      &  zaeq4(:,:,:),         & !< for 5 different                 4: urban
      &  zaeq5(:,:,:)            !< aerosol species.                5: stratospheric background
    REAL(wp), POINTER, INTENT(inout) :: &
      &  od_lw(:,:,:,:),       & !< Longwave optical thickness
      &  od_sw(:,:,:,:),       & !< Shortwave optical thickness
      &  ssa_sw(:,:,:,:),      & !< Shortwave asymmetry factor
      &  g_sw(:,:,:,:)           !< Shortwave single scattering albedo
    ! Local variables
    INTEGER :: istat
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    CALL assert_acc_device_only("nwp_aerosol_cleanup", lacc)

    !$ACC EXIT DATA DETACH(od_lw) IF(ASSOCIATED(od_lw))
    !$ACC EXIT DATA DETACH(od_sw) IF(ASSOCIATED(od_sw))
    !$ACC EXIT DATA DETACH(ssa_sw) IF(ASSOCIATED(ssa_sw))
    !$ACC EXIT DATA DETACH(g_sw) IF(ASSOCIATED(g_sw))
    NULLIFY(od_lw)
    NULLIFY(od_sw)
    NULLIFY(ssa_sw)
    NULLIFY(g_sw)

    !$ACC WAIT
    IF( ALLOCATED(zaeq1) ) THEN
      !$ACC EXIT DATA DELETE(zaeq1)
      DEALLOCATE(zaeq1, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of zaeq1 failed.')
    ENDIF
    IF( ALLOCATED(zaeq2) ) THEN
      !$ACC EXIT DATA DELETE(zaeq2)
      DEALLOCATE(zaeq2, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of zaeq2 failed.')
    ENDIF
    IF( ALLOCATED(zaeq3) ) THEN
      !$ACC EXIT DATA DELETE(zaeq3)
      DEALLOCATE(zaeq3, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of zaeq3 failed.')
    ENDIF
    IF( ALLOCATED(zaeq4) ) THEN
      !$ACC EXIT DATA DELETE(zaeq4)
      DEALLOCATE(zaeq4, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of zaeq4 failed.')
    ENDIF
    IF( ALLOCATED(zaeq5) ) THEN
      !$ACC EXIT DATA DELETE(zaeq5)
      DEALLOCATE(zaeq5, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of zaeq5 failed.')
    ENDIF

    IF( ALLOCATED(locmem_od_lw) ) THEN
      !$ACC EXIT DATA DELETE(locmem_od_lw)
      DEALLOCATE(locmem_od_lw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of locmem_od_lw failed.')
    ENDIF
    IF( ALLOCATED(locmem_od_sw) ) THEN
      !$ACC EXIT DATA DELETE(locmem_od_sw)
      DEALLOCATE(locmem_od_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of locmem_od_sw failed.')
    ENDIF
    IF( ALLOCATED(locmem_ssa_sw) ) THEN
      !$ACC EXIT DATA DELETE(locmem_ssa_sw)
      DEALLOCATE(locmem_ssa_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of locmem_ssa_sw failed.')
    ENDIF
    IF( ALLOCATED(locmem_g_sw) ) THEN
      !$ACC EXIT DATA DELETE(locmem_g_sw)
      DEALLOCATE(locmem_g_sw, STAT=istat)
      IF(istat /= SUCCESS) CALL finish(routine, 'Deallocation of locmem_g_sw failed.')
    ENDIF

  END SUBROUTINE nwp_aerosol_cleanup

END MODULE mo_nwp_aerosol

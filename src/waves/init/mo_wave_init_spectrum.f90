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

MODULE mo_wave_init_spectrum

  USE mo_kind,                 ONLY: wp
  USE mo_exception,            ONLY: message
  USE mo_impl_constants,       ONLY: MAX_CHAR_LENGTH, VNAME_LEN, min_rlcell
  USE mo_math_constants,       ONLY: rpi_2
  USE mo_model_domain,         ONLY: t_patch
  USE mo_loopindices,          ONLY: get_indices_c
  USE mo_parallel_config,      ONLY: nproma
  USE mo_wave_phy_util,        ONLY: jonswap_nonblk, fetch_law_nonblk
  USE fortran_support,         ONLY: t_ptr_3d_wp, t_ptr_2d3d
  USE mo_wave_types,           ONLY: t_wave_state, t_wesd
  USE mo_wave_config,          ONLY: t_wave_config, generate_filename
  USE mo_wave_constants,       ONLY: EMIN
  USE mo_io_units,             ONLY: filename_max
  USE mo_read_interface,       ONLY: openInputFile, closeFile, t_stream_id, on_cells, &
    &                                read_2D_1time
  USE mo_io_config,            ONLY: default_read_method
  USE mo_sync,                 ONLY: SYNC_C, sync_patch_array_mult
  USE mo_grid_config,          ONLY: nroot
  USE mo_dynamics_config,      ONLY: nnow
  USE mo_initwave_config,      ONLY: initwave_config
  USE mo_time_config,          ONLY: time_config
  USE mo_master_config,        ONLY: getModelBaseDir
  USE mo_post_op,              ONLY: inverse_post_op
  USE mo_timer,                ONLY: timer_start, timer_stop, timers_level
  USE mo_wave_timer,           ONLY: timer_wave_exch

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_wave_spectrum_analytic
  PUBLIC :: init_wave_spectrum_from_file

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_init_spectrum'

CONTAINS


  !>
  !! Initialisation of the wave spectrum by the analytic 1D JONSWAP spectrum
  !!
  !! Calculation of wind dependent initial spectrum from
  !! the fetch law and from the 1D JONSWAP spectrum. The minimum
  !! of wave energy is limited to FLMIN.
  !!
  SUBROUTINE init_wave_spectrum_analytic(p_patch, wave_config, sp10m, dir10m, wesd)
    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':init_wave_spectrum_analytic'

    TYPE(t_patch),               INTENT(IN)    :: p_patch
    TYPE(t_wave_config), TARGET, INTENT(IN)    :: wave_config
    REAL(wp),                    INTENT(IN)    :: sp10m(:,:)      !< 10m wind speed (m/s)
    REAL(wp),                    INTENT(IN)    :: dir10m(:,:)     !< wind direction (rad)
    TYPE(t_wesd),                INTENT(INOUT) :: wesd(:)         !< wave energy spectrum (m**2 s)

    INTEGER :: i_rlstart, i_rlend, i_startblk, i_endblk
    INTEGER :: i_startidx, i_endidx
    INTEGER :: jc,jb,jd,jf
    REAL(wp):: alphaj(nproma)                 ! Phillips constant (scaling parameter)
    REAL(wp):: fp(nproma)                     ! JONSWAP peak frequency (where spectrum has maximum energy)
    REAL(wp):: et(nproma, wave_config%nfreqs) ! JONSWAP spectrum
    REAL(wp):: st

    TYPE(t_wave_config), POINTER :: wc => NULL()

    ! save some paperwork
    wc => wave_config

    ! halo points must be included
    i_rlstart  = 1
    i_rlend    = min_rlcell
    i_startblk = p_patch%cells%start_block(i_rlstart)
    i_endblk   = p_patch%cells%end_block(i_rlend)

!$OMP PARALLEL
!$OMP DO PRIVATE(jb,jf,jd,jc,i_startidx,i_endidx,st,fp,alphaj,et)
    DO jb = i_startblk, i_endblk
      CALL get_indices_c( p_patch, jb, i_startblk, i_endblk,           &
        &                 i_startidx, i_endidx, i_rlstart, i_rlend)

      ! Calculate the peak frequency fp from a fetch law
      ! and set the Phillips constant alphaj
      CALL fetch_law_nonblk(              &
        &     i_startidx  = i_startidx,   & !in
        &     i_endidx    = i_endidx,     & !in
        &     fetch       = wc%fetch,     & !in
        &     fpmax       = wc%fm,        & !in
        &     sp10m       = sp10m(:,jb),  & !in
        &     fp          = fp(:),        & !out
        &     alphaj      = alphaj(:))      !out


      ! Set the 1D JONSWAP spectrum
      !
      CALL jonswap_nonblk(               &
        &   i_startidx = i_startidx,     & !in
        &   i_endidx   = i_endidx,       & !in
        &   freqs      = wc%freqs(:),    & !in
        &   gamma      = wc%GAMMA_wave,  & !in
        &   sa         = wc%SIGMA_A,     & !in
        &   sb         = wc%SIGMA_B,     & !in
        &   flmin      = wc%flmin,       & !in
        &   alphaj     = alphaj(:),      & !in
        &   fp         = fp(:),          & !in
        &   et         = et(:,:)         ) !out


      DO jf = 1,wc%nfreqs
        DO jd = 1,wc%ndirs
          DO jc = i_startidx, i_endidx
            st = rpi_2*MAX(0._wp, COS(wc%dirs(jd)-dir10m(jc,jb)) )**2
            IF (st < 0.1E-08_wp) st = 0._wp

            ! Wave model initialisation
            wesd(jf)%ptr(jc,jd,jb) = et(jc,jf) * st
            wesd(jf)%ptr(jc,jd,jb) = MAX(wesd(jf)%ptr(jc,jd,jb),EMIN)
          END DO  !jc
        END DO  !jd
      END DO  !jf
    END DO  !jb
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    CALL message(routine,'finished')

  END SUBROUTINE init_wave_spectrum_analytic


  !>
  ! -----------------------------------------------------------------------
  ! SUBROUTINE: init_wave_spectrum_from_file
  ! PURPOSE:
  !   This subroutine initializes the wave spectrum by reading data from
  !   an external state file specified in the wave configuration.
  !   The data is read for each component of the energy spectrum and
  !   assigned to the wave state structure (p_wave_state) for use in the
  !   simulation.
  ! -----------------------------------------------------------------------
  SUBROUTINE init_wave_spectrum_from_file(p_patch, wave_config, p_wave_state)
    TYPE(t_patch),       INTENT(IN   ) :: p_patch
    TYPE(t_wave_config), INTENT(IN   ) :: wave_config
    TYPE(t_wave_state),  INTENT(INOUT) :: p_wave_state

    TYPE(t_stream_id) :: stream_id

    CHARACTER(len=MAX_CHAR_LENGTH), PARAMETER ::  &
         &  routine = modname//':init_spectrum_from_file'

    CHARACTER(LEN=filename_max) :: filename
    CHARACTER(LEN=VNAME_LEN) :: freq_ind_str, dir_ind_str
    CHARACTER(LEN=VNAME_LEN) :: wesd_name
    INTEGER :: jg, jf, jd
    INTEGER :: ndirs
    TYPE(t_ptr_3d_wp):: wesd_ptr(wave_config%nfreqs)
    TYPE(t_ptr_2d3d) :: input_data

    ! 1. Setup parameters
    ndirs = wave_config%ndirs
    jg = p_patch%id


    ! 2. Build filename dynamically using generator
    filename = TRIM(generate_filename( &
      & initwave_config(jg)%initial_wave_spectrum_filename, &
      & getModelBaseDir(), nroot, 1, jg, time_config%tc_exp_startdate))

    ! Log which file we are reading
    CALL message(routine, 'Reading initial wave spectrum from file: '//TRIM(filename))

    ! 3. Open file and read data
    CALL openInputFile(stream_id, TRIM(filename), p_patch, default_read_method)


    DO jf = 1,wave_config%nfreqs
      write(freq_ind_str,'(I0.3)') jf

      DO jd = 1,wave_config%ndirs
        ! wesd
        write(dir_ind_str,'(I0.3)') jd
        wesd_name = 'wesd_f'//TRIM(freq_ind_str)//'_d'//TRIM(dir_ind_str)

        input_data%p_2d => p_wave_state%prog(nnow(jg))%wesd(jf)%ptr(:,jd,:)
        CALL read_2D_1time(stream_id, on_cells, wesd_name, input_data%p_2d)
        !
        ! field conversion to internally used SI units
        IF (.NOT. initwave_config(jg)%lskip_inv_post_op) THEN
          CALL inverse_post_op(TRIM(wesd_name), input_data%p_2d)
        ENDIF
      ENDDO
      !
      ! structure of type t_ptr_3d_wp required for halo synchronization below
      wesd_ptr(jf)%p => p_wave_state%prog(nnow(jg))%wesd(jf)%ptr(:,:,:)
    END DO !frequencies

    CALL closeFile(stream_id)

    IF (timers_level >= 5) CALL timer_start(timer_wave_exch)

    CALL sync_patch_array_mult(typ         = SYNC_C,          &
      &                        p_patch     = p_patch,         &
      &                        nfields     = SIZE(wesd_ptr),  &
      &                        f3din_arr   = wesd_ptr,        &
      &                        opt_varname = 'wesd_now',      &
      &                        lacc        = .FALSE.)

    IF (timers_level >= 5) CALL timer_stop(timer_wave_exch)

    CALL message(routine, 'finished')
    !
  END SUBROUTINE init_wave_spectrum_from_file

END MODULE mo_wave_init_spectrum

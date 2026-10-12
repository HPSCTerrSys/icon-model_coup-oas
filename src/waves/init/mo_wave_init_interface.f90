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

! Initialize wave energy spectrum
!
! The following options exist:
! - read energy spectrum from file
! - initialize energy spectrum analytically by wind-speed based parameterization
!   (such as JONSWAP)
!
!
MODULE mo_wave_init_interface

  USE mo_exception,            ONLY: message, message_text, finish
  USE mo_timer,                ONLY: timers_level, timer_start, timer_stop, timer_read_restart
  USE mo_wave_timer,           ONLY: timer_wave_reader
  USE mo_wave_constants,       ONLY: MODE_ANA, MODE_COLD
  USE mo_master_config,        ONLY: isInitFromRestart
  USE mo_grid_config,          ONLY: n_dom
  USE mo_wave_config,          ONLY: wave_config
  USE mo_initwave_config,      ONLY: initwave_config
  USE mo_run_config,           ONLY: ltestcase
  USE mo_dynamics_config,      ONLY: nnow
  USE mo_time_config,          ONLY: time_config
  USE mo_model_domain,         ONLY: t_patch
  USE mo_wave_types,           ONLY: t_wave_state
  USE mo_wave_ext_data_types,  ONLY: t_external_wave
  USE mo_wave_forcing_types,   ONLY: t_wave_forcing
  USE mo_wave_forcing,         ONLY: reader_wave_forcing
  USE mo_intp_data_strc,       ONLY: t_int_state
  USE mo_wave_adv_exp,         ONLY: init_analytic_forcing
  USE mo_wave_td_update,       ONLY: update_water_depth_and_grad
  USE mo_wave_init_spectrum,   ONLY: init_wave_spectrum_analytic, init_wave_spectrum_from_file
  USE mo_load_restart,         ONLY: read_restart_files

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_init_interface'

  PUBLIC :: init_wave

CONTAINS

  ! Interface for wave model initialization, in particular the
  ! wave energy spectrum.
  !
  ! The following options are available:
  ! init_mode = MODE_COLD
  !   initialize wave spectrum by analytic JONSWAP spectrum
  !   ltestcase=.FALSE.: velocity vector for initialization is read from file
  !   ltestcase=.TRUE. : velocity vector is set internally by analytic function
  !
  ! init_mode = MODE_ANA
  !   read wave spectrum from analysis or first guess file
  !   isInitFromRestart=.TRUE.: the ICON restart input machinery is (mis)used
  !                             for reading the spectrum
  !   isInitFromRestart=.FALSE.: TO BE IMPLEMENTED
  !                             the ICON input module is used for reading the spectrum
  !
  SUBROUTINE init_wave (p_patch, p_int_state, &
    &                   wave_state, wave_forcing_state, wave_ext_data)

    TYPE(t_patch),         INTENT   (IN) :: p_patch(:)
    TYPE(t_int_state),     INTENT   (IN) :: p_int_state(:)
    TYPE(t_wave_state),    INTENT(INOUT) :: wave_state(:)
    TYPE(t_wave_forcing),  INTENT(INOUT) :: wave_forcing_state(:)
    TYPE(t_external_wave), INTENT(INOUT) :: wave_ext_data(:)

    CHARACTER(LEN = *), PARAMETER :: routine = modname//':init_wave'

    ! local
    INTEGER :: jg

    DO jg = 1,n_dom

      IF (.NOT. p_patch(jg)%ldom_active) CYCLE

      SELECT CASE(initwave_config(jg)%init_mode)
      CASE (MODE_ANA)
        !
        ! Initialize wave spectrum from analysis/first guess
        !
        CALL message(routine,'Warmstart: Initialize wave spectrum from analysis/first guess')

        IF (isInitFromRestart()) THEN
          IF (timers_level > 4) CALL timer_start(timer_read_restart)
          !
          IF (p_patch(jg)%ldom_active) THEN
            CALL message(routine,'Warmstart: Initialize from restart file')
            CALL read_restart_files( p_patch(jg), n_dom)
          END IF
          !
          CALL message(routine,'normal exit from read_restart_files')
          !
          IF (timers_level > 4) CALL timer_stop(timer_read_restart)
        ELSE
          ! Initialize from First Guess or analysis file
          !
          CALL init_wave_spectrum_from_file(p_patch(jg), wave_config(jg), wave_state(jg))
          !
          !
          !
        ENDIF

      CASE (MODE_COLD)
        !
        ! Coldstart: Initialize by analytic JONSWAP spectrum
        !
        CALL message(routine,'Coldstart: Initialize wave spectrum by analytic JONSWAP spectrum')
        IF (ltestcase) THEN
          !-----------------------------------------------------------------------
          ! advection experiment
          CALL message(routine,'Coldstart: use analytic wind field for spectrum initialization')

          ! Analytic time-constant initialisation of the following forcing fields:
          ! u10m, v10m, sp10m, dir10m, sea_ice_c, ice_free_mask
          !
          CALL init_analytic_forcing(p_patch(jg), wave_config(jg), wave_forcing_state(jg))
        ELSE
          ! ltestcase=.FALSE.
          !
          IF (wave_config(jg)%lread_forcing) THEN
            IF (timers_level >= 5) CALL timer_start(timer_wave_reader)

            CALL message(routine,'Coldstart: read wind field for spectrum initialization from file')
            !
            ! get initial forcing data set (read from file and copy to forcing state vector)
            CALL reader_wave_forcing(jg)%update_forcing(                     &
              &     destination_time = time_config%tc_current_date,          & !in
              &     u10m             = wave_forcing_state(jg)%u10m,          & !out
              &     v10m             = wave_forcing_state(jg)%v10m,          & !out
              &     sp10m            = wave_forcing_state(jg)%sp10m,         & !out
              &     dir10m           = wave_forcing_state(jg)%dir10m,        & !out
              &     sic              = wave_forcing_state(jg)%sea_ice_c,     & !out
              &     slh              = wave_forcing_state(jg)%sea_level_c,   & !out
              &     uosc             = wave_forcing_state(jg)%usoce_c,       & !out
              &     vosc             = wave_forcing_state(jg)%vsoce_c,       & !out
              &     sp_osc           = wave_forcing_state(jg)%sp_soce_c,     & !out
              &     dir_osc          = wave_forcing_state(jg)%dir_soce_c,    & !out
              &     ice_free_mask_c  = wave_forcing_state(jg)%ice_free_mask_c) !out

            ! update depth and gradient
            CALL update_water_depth_and_grad(                              &
              &     p_patch          = p_patch(jg),                        & !in
              &     p_int_state      = p_int_state(jg),                    & !in
              &     bathymetry_c     = wave_ext_data(jg)%bathymetry_c,     & !in
              &     sea_level_c      = wave_forcing_state(jg)%sea_level_c, & !in
              &     depth_c          = wave_ext_data(jg)%depth_c,          & !out
              &     depth_e          = wave_ext_data(jg)%depth_e,          & !out
              &     geo_depth_grad_c = wave_ext_data(jg)%geo_depth_grad_c)   !out

            IF (timers_level >= 5) CALL timer_stop(timer_wave_reader)
          ELSE
            WRITE(message_text,'(a,a,a)') 'No forcing files specified, testcase run is assumed.'
            CALL message(routine, message_text)
          ENDIF ! lread_forcing
        ENDIF

        ! Initialisation of the wave spectrum by the analytic 1D JONSWAP spectrum
        CALL init_wave_spectrum_analytic(                                    &
          &     p_patch     = p_patch(jg),                                   & !in
          &     wave_config = wave_config(jg),                               & !in
          &     sp10m       = wave_forcing_state(jg)%sp10m(:,:),             & !in
          &     dir10m      = wave_forcing_state(jg)%dir10m(:,:),            & !in
          &     wesd        = wave_state(jg)%prog(nnow(jg))%wesd(:))           !out

      CASE DEFAULT
        CALL finish(routine, "Invalid operation mode!")
      END SELECT
    ENDDO

  END SUBROUTINE init_wave

END MODULE mo_wave_init_interface

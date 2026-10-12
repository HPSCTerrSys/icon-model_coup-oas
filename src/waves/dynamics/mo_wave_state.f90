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

! Contains the variables to set up the wave model.

MODULE mo_wave_state

  USE mo_kind,                      ONLY: wp
  USE mo_master_control,            ONLY: get_my_process_name
  USE mo_exception,                 ONLY: message, finish
  USE mo_parallel_config,           ONLY: nproma
  USE mo_model_domain,              ONLY: t_patch
  USE mo_grid_config,               ONLY: n_dom, l_limited_area, ifeedback_type
  USE mo_coupling_config,           ONLY: is_coupled_to_ocean
  USE mo_impl_constants,            ONLY: success, max_char_length, VNAME_LEN, TLEV_NNOW, &
    &                                     HINTP_TYPE_LONLAT_NNB, HINTP_TYPE_LONLAT_BCTR
  USE mo_math_constants,            ONLY: rad2deg
  USE mo_var_list,                  ONLY: add_var, add_ref, t_var_list_ptr
  USE mo_math_constants,            ONLY: pi2
  USE mo_var_list_register,         ONLY: vlr_add, vlr_del
  USE mo_var_groups,                ONLY: groups
  USE mo_cdi_constants,             ONLY: GRID_UNSTRUCTURED_CELL, GRID_CELL, &
    &                                     GRID_UNSTRUCTURED_EDGE, GRID_EDGE
  USE mo_cdi,                       ONLY: DATATYPE_FLT32, DATATYPE_FLT64, GRID_UNSTRUCTURED, &
    &                                     DATATYPE_PACK16, DATATYPE_INT
  USE mo_zaxis_type,                ONLY: ZA_SURFACE, ZA_FREQ_GENERIC, ZA_DIR_GENERIC, &
    &                                     ZA_DEPTH_BELOW_SEA
  USE mo_cf_convention,             ONLY: t_cf_var
  USE mo_grib2,                     ONLY: t_grib2_var, grib2_var, t_grib2_int_key, OPERATOR(+), &
    &                                     t_grib2_intarr_key
  USE mo_io_config,                 ONLY: lnetcdf_flt64_output
  USE mo_wave_io_config,            ONLY: t_wave_var_in_output
  USE mo_var_metadata,              ONLY: get_timelevel_string, create_hor_interp_metadata, post_op
  USE mo_var_metadata_types,        ONLY: CLASS_WAVE_SPECTRUM, POST_OP_SCALE
  USE mo_tracer_metadata,           ONLY: create_tracer_metadata
  USE mo_wave_types,                ONLY: t_wave_prog, t_wave_source, t_wave_const, &
    &                                     t_wave_diag_dyn, t_wave_diag_out, t_wave_diag_cpl, &
    &                                     t_wave_state, t_wave_state_lists
  USE mo_wave_config,               ONLY: t_wave_config, wave_config


  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_state'

  PUBLIC :: construct_wave_state    ! Constructor for the wave state
  PUBLIC :: destruct_wave_state     ! Destructor

  PUBLIC :: p_wave_state            ! state vector of wave variables
  PUBLIC :: p_wave_state_lists      ! lists for state vector of wave variables

  TYPE(t_wave_state),       TARGET, ALLOCATABLE :: p_wave_state(:)
  TYPE(t_wave_state_lists), TARGET, ALLOCATABLE :: p_wave_state_lists(:)

CONTAINS

  SUBROUTINE construct_wave_state(p_patch, n_timelevels, var_in_output)

    TYPE(t_patch),              INTENT(IN) :: p_patch(:)
    INTEGER,                    INTENT(IN) :: n_timelevels
    TYPE(t_wave_var_in_output), INTENT(IN) :: var_in_output(:) !< switches for optional diagnostics


    CHARACTER(len=max_char_length) :: listname
    CHARACTER(len=*), PARAMETER :: routine = modname//'::construct_wave_state'

    INTEGER :: ntl,  &! local number of timelevels
         ist,        &! status
         jg,         &! grid level counter
         jt           ! time level counter

    CALL message (routine, 'Construction of wave state started')

    ALLOCATE (p_wave_state(n_dom),p_wave_state_lists(n_dom), stat=ist)
    IF (ist /= success) THEN
       CALL finish(routine,'allocation for wave state failed')
    END IF

    DO jg = 1, n_dom

       ntl = n_timelevels

       ! As grid nesting is not called at every dynamics time step, an extra time
       ! level is needed for full-field interpolation and boundary-tendency calculation
       IF (n_dom > 1) THEN
          ntl = ntl + 1
       END IF

       IF (ifeedback_type == 1 .AND. jg > 1 .OR. l_limited_area .AND. jg == 1) ntl = ntl + 1

       ALLOCATE(p_wave_state(jg)%prog(1:ntl), STAT=ist)
       IF (ist/=SUCCESS) CALL finish(routine,                                   &
            'allocation of prognostic state array failed')

       ! create state list
       ALLOCATE(p_wave_state_lists(jg)%prog_list(1:ntl), STAT=ist)
       IF (ist/=SUCCESS) CALL finish(routine,                                   &
            'allocation of prognostic state list array failed')

       !
       ! Build lists for every timelevel
       !
       DO jt = 1, ntl
         WRITE(listname,'(a,i2.2,a,i2.2)') 'wave_state_prog_of_domain_',jg, &
             &                            '_and_timelev_',jt
         ! Build prog state list
         ! includes memory allocation
         CALL new_wave_state_prog_list(p_patch(jg), p_wave_state(jg)%prog(jt), &
              & p_wave_state_lists(jg)%prog_list(jt), &
              & listname, jt)
       END DO

       ! Build source state list
       ! includes memory allocation
       WRITE(listname,'(a,i2.2)') 'wave_state_source_of_domain_',jg
       CALL new_wave_state_source_list(&
            p_patch(jg), &
            p_wave_state(jg)%source, &
            p_wave_state_lists(jg)%source_list, &
            listname)

       ! Build const state list
       ! includes memory allocation
       WRITE(listname,'(a,i2.2)') 'wave_state_const_of_domain_',jg
       CALL new_wave_state_const_list(&
            p_patch(jg), &
            p_wave_state(jg)%const, &
            p_wave_state_lists(jg)%const_list, &
            listname)

       ! Build diag_dyn state list
       ! includes memory allocation
       WRITE(listname,'(a,i2.2)') 'wave_state_diag_dyn_of_domain_',jg
       CALL new_wave_state_diag_dyn_list(&
            p_patch(jg), &
            p_wave_state(jg)%diag_dyn, &
            p_wave_state_lists(jg)%diag_dyn_list, &
            listname)

       ! Build diag_out state list
       ! includes memory allocation
       WRITE(listname,'(a,i2.2)') 'wave_state_diag_out_of_domain_',jg
       CALL new_wave_state_diag_out_list(&
            p_patch(jg), &
            p_wave_state(jg)%diag_out, &
            p_wave_state_lists(jg)%diag_out_list, &
            listname)

       ! Build diag_cpl state list
       ! includes memory allocation
       WRITE(listname,'(a,i2.2)') 'wave_state_diag_cpl_of_domain_',jg
       CALL new_wave_state_diag_cpl_list(&
            p_patch(jg), &
            p_wave_state(jg)%diag_cpl, &
            p_wave_state_lists(jg)%diag_cpl_list, &
            listname, &
            var_in_output(jg))
    END DO

    CALL message (routine, 'wave state construction completed')

  END SUBROUTINE construct_wave_state


  SUBROUTINE new_wave_state_prog_list(p_patch, p_prog, p_prog_list, listname, timelev)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_prog),     INTENT(INOUT) :: p_prog
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_prog_list !< current prognostic state list
    CHARACTER(len=*),      INTENT(IN)    :: listname
    INTEGER,               INTENT(IN)    :: timelev

    CHARACTER(len=*), PARAMETER :: routine = modname//'::new_wave_state_prog_list'
    TYPE(t_cf_var)    :: cf_desc, new_cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: nblks_c       !< number of cell blocks to allocate
    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: ist

    CHARACTER(len=4)         :: suffix
    CHARACTER(len=VNAME_LEN) :: dir_ind_str, freq_ind_str
    CHARACTER(LEN=VNAME_LEN) :: wesd_container_name
    CHARACTER(LEN=VNAME_LEN) :: wesd_name

    TYPE(t_wave_config),               POINTER :: wc

    INTEGER :: shape3d_c(3), shape2d_c(2)
    INTEGER :: jf, jd

    INTEGER :: scaleFactor
    INTEGER :: scaled_wdsp1, scaled_wdsp2  ! scaled Grib2 waveDirectionSequenceParameters
    INTEGER :: scaled_wfsp1, scaled_wfsp2  ! scaled Grib2 waveFrequencySequenceParameters
    INTEGER :: scaleFactorArr(2), wdspArr(2), wfspArr(2)

    !determine size of arrays
    nblks_c = p_patch%nblks_c

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    shape3d_c = (/nproma, wc%ndirs, nblks_c/)
    shape2d_c = (/nproma, nblks_c/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF (lnetcdf_flt64_output) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    END IF

    ! Suffix (mandatory for time level dependent variables)
    suffix = get_timelevel_string(timelev)

    ! allocate wesd state
    ALLOCATE(p_prog%wesd(wc%nfreqs), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(routine, 'allocation of prognostic wesd state failed')

    ! scaling for waveDirectionSequenceParameters and waveFrequencySequenceParameters
    ! according to
    ! value = scaledValue * 10^{-scaleFactor}
    ! strictly speaking the scaleFactor refers to a power-of-10 exponent!
    !
    scaleFactor    = 6
    scaleFactorArr = (/scaleFactor, scaleFactor/)
    !
    ! scaled waveDirectionSequenceParameters
    ! hint: arithmetic sequence for direction calculation (see mo_wave_config)
    ! dirs(jd) =  0.5_wp*wc%delth + REAL(jd-1,wp) * wc%delth
    scaled_wdsp1 = NINT(rad2deg * 0.5_wp*wc%delth * 10**(scaleFactor))
    scaled_wdsp2 = NINT(rad2deg * wc%delth * 10**(scaleFactor))
    wdspArr      = (/scaled_wdsp1, scaled_wdsp2/)
    !
    ! scaled waveFrequencySequenceParameter
    ! hint: geometric sequence for frequency calculation (see mo_wave_config)
    ! freqs(jf) = wc%fr1 * wc%co**(jf-1)
    scaled_wfsp1 = NINT(wc%fr1 * 10**(scaleFactor))
    scaled_wfsp2 = NINT(wc%co * 10**(scaleFactor))
    wfspArr      = (/scaled_wfsp1, scaled_wfsp2/)

    !
    ! Register a field list and apply default settings
    !
    CALL vlr_add(p_prog_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &          model_type=get_my_process_name())

    DO jf = 1,wc%nfreqs
      write(freq_ind_str,'(I0.3)') jf

      ! wesd        wesd(jf)%ptr(nproma,ndirs,nblks_c)
      wesd_container_name = 'wesd_f'//TRIM(freq_ind_str)//suffix
      cf_desc    = t_cf_var(TRIM(wesd_container_name), 'm^2 s', &
        &                   'wave energy spectral density', datatype_flt)
      grib2_desc = grib2_var(10, 0, 42, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var( p_prog_list, wesd_container_name, p_prog%wesd(jf)%ptr, &
           & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC, cf_desc, grib2_desc,  &
           & ldims=shape3d_c,                                              &
           & var_class=CLASS_WAVE_SPECTRUM,                                &
           & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE. )

      ALLOCATE(p_prog%wesd(jf)%dir(wc%ndirs), STAT=ist)
      IF (ist/=SUCCESS) CALL finish(routine, &
        &                    'allocation of p_prog%wesd(jf)%dir failed')

      DO jd = 1, wc%ndirs
        write(dir_ind_str,'(I0.3)') jd
        wesd_name = 'wesd_f'//TRIM(freq_ind_str)//'_d'//TRIM(dir_ind_str)//suffix

        cf_desc    = t_cf_var(TRIM(wesd_name), 'm^2 s', 'wave energy spectral density', datatype_flt)
        new_cf_desc = t_cf_var(TRIM(wesd_name), 'm^2 s rad^-1', 'wave energy spectral density', &
          &                    datatype_flt)
        grib2_desc = grib2_var(10, 0, 42, ibits, GRID_UNSTRUCTURED, GRID_CELL)                         &
          &        + t_grib2_int_key   ("numberOfWaveDirections", wc%ndirs)                            &
          &        + t_grib2_int_key   ("typeOfWaveDirectionSequence", 2)                              &
          &        + t_grib2_int_key   ("waveDirectionNumber", jd)                                     &
          &        + t_grib2_int_key   ("numberOfWaveDirectionSequenceParameters", 2)                  &
          &        + t_grib2_intarr_key("scaleFactorOfWaveDirectionSequenceParameter", scaleFactorArr) &
          &        + t_grib2_intarr_key("scaledValueOfWaveDirectionSequenceParameter", wdspArr)        &
          &        + t_grib2_int_key   ("numberOfWaveFrequencies", wc%nfreqs)                          &
          &        + t_grib2_int_key   ("typeOfWaveFrequencySequence", 1)                              &
          &        + t_grib2_int_key   ("waveFrequencyNumber", jf)                                     &
          &        + t_grib2_int_key   ("numberOfWaveFrequencySequenceParameters", 2)                  &
          &        + t_grib2_intarr_key("scaleFactorOfWaveFrequencySequenceParameter", scaleFactorArr) &
          &        + t_grib2_intarr_key("scaledValueOfWaveFrequencySequenceParameter", wfspArr)

        CALL add_ref( p_prog_list, wesd_container_name,                          &
          &  TRIM(wesd_name), p_prog%wesd(jf)%dir(jd)%p_2d,                      &
          &  GRID_UNSTRUCTURED_CELL, ZA_SURFACE,                                 &
          &  cf_desc, grib2_desc,                                                &
          &  ldims=shape2d_c, opt_var_ref_pos=2, ref_idx=jd,                     &
          &  loutput=.TRUE., lrestart=.TRUE.,                                    &
          &  tlev_source=TLEV_NNOW,                                              &
          &  tracer_info=create_tracer_metadata(lis_tracer=.TRUE.,               &
          &                        name        = TRIM(wesd_name),                &
          &                        lfeedback   = .TRUE.,                         &
          &                        ihadv_tracer= 2,                              &
          &                        ivadv_tracer= 0 ),                            &
          &  var_class=CLASS_WAVE_SPECTRUM,                                      &
          &  post_op=post_op(POST_OP_SCALE, arg1=1._wp/pi2, new_cf=new_cf_desc), &
          &  in_group=groups("wave_spectrum","DWD_FG_WAVE_VARS"))
      ENDDO
    ENDDO

  END SUBROUTINE new_wave_state_prog_list



  SUBROUTINE new_wave_state_source_list(p_patch, p_source, p_source_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_source),   INTENT(INOUT) :: p_source
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: p_source_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    CHARACTER(len=*), PARAMETER :: routine = modname//'::new_wave_state_source_list'

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: nblks_c
    INTEGER :: jf
    INTEGER :: ist
    INTEGER :: shape4d_c(4), shape3d_c(3)
    CHARACTER(len=VNAME_LEN) :: sl_name, fl_name, llws_name
    CHARACTER(len=3) :: freq_ind_str

    TYPE(t_wave_config), POINTER :: wc

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(p_source%fl,  &
      &     p_source%sl,  &
      &     p_source%llws)

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    nblks_c = p_patch%nblks_c

    shape4d_c = (/nproma, wc%ndirs, wc%nfreqs, nblks_c/)
    shape3d_c = (/nproma, wc%ndirs, nblks_c/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(p_source_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &           model_type=get_my_process_name())


    ! fl          p_source%fl(nproma,ndirs,nfreqs,nblks_c)
    cf_desc    = t_cf_var('fl', '-', 'DIAG. MTRX OF FUNC. DERIVATIVE', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_source_list, 'fl', p_source%fl,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC, cf_desc, grib2_desc, &
         & ldims=shape4d_c, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! sl          p_source%sl(nproma,ndirs,nfreqs,nblks_c)
    cf_desc    = t_cf_var('sl', '-', 'TOTAL SOURCE FUNCTION', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_source_list, 'sl', p_source%sl,                      &
         & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC, cf_desc, grib2_desc, &
         & ldims=shape4d_c, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)

    ! llws        p_source%llws(nproma,ndirs,nfreqs,nblks_c)
    cf_desc    = t_cf_var('llws', '-', '1 where sinput is positive', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(p_source_list, 'llws', p_source%llws,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC, cf_desc, grib2_desc, &
         & ldims=shape4d_c,                                             &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.)


    ALLOCATE(p_source%sl_ptr(wc%nfreqs),           &
      &      p_source%fl_ptr(wc%nfreqs),           &
      &      p_source%llws_ptr(wc%nfreqs), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(routine, &
      &                    'allocation of sl_ptr, fl_ptr and llws_ptr failed')

    DO jf = 1, wc%nfreqs
      write(freq_ind_str,'(I0.3)') jf

      sl_name = 'sl_f'//TRIM(freq_ind_str)
      CALL add_ref(p_source_list, 'sl',                                     &
           & sl_name, p_source%sl_ptr(jf)%p_3d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC,                        &
           & t_cf_var(sl_name, '-',sl_name, datatype_flt),                  &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ldims=shape3d_c, opt_var_ref_pos=3, ref_idx=jf,                &
           & lrestart=.FALSE., loutput=.TRUE.)

      fl_name = 'fl_f'//TRIM(freq_ind_str)
      CALL add_ref(p_source_list, 'fl',                                     &
           & fl_name, p_source%fl_ptr(jf)%p_3d,                             &
           & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC,                        &
           & t_cf_var(fl_name, '-',fl_name, datatype_flt),                  &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ldims=shape3d_c, opt_var_ref_pos=3, ref_idx=jf,                &
           & lrestart=.FALSE., loutput=.TRUE.)

      llws_name = 'llws_f'//TRIM(freq_ind_str)
      CALL add_ref(p_source_list, 'llws',                                   &
           & llws_name, p_source%llws_ptr(jf)%p,                            &
           & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC,                        &
           & t_cf_var(llws_name, '-',llws_name, datatype_int),              &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & ldims=shape3d_c, opt_var_ref_pos=3, ref_idx=jf,                &
           & lrestart=.FALSE., loutput=.TRUE.,                              &
           & in_group=groups("wave_debug","DWD_FG_WAVE_VARS"))
    END DO

  END SUBROUTINE new_wave_state_source_list


  !>
  !! Variable list collecting constant coefficient fields for model integration.
  !! Fields are computed once during model initialization.
  !!
  SUBROUTINE new_wave_state_const_list(p_patch, const, const_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_const),    INTENT(INOUT) :: const
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: const_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    CHARACTER(len=*), PARAMETER :: routine = modname//'::new_wave_state_const_list'

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: shape1d_freq_p4(1), shape1d_dir_2(2)

    TYPE(t_wave_config),      POINTER :: wc

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(const%flminfr_tab, &
    &       const%IKP,         &
    &       const%IKP1,        &
    &       const%IKM,         &
    &       const%IKM1,        &
    &       const%K1W,         &
    &       const%K2W,         &
    &       const%K11W,        &
    &       const%K21W,        &
    &       const%JA1,         &
    &       const%JA2,         &
    &       const%AF11,        &
    &       const%FKLAP,       &
    &       const%FKLAP1,      &
    &       const%FKLAM,       &
    &       const%FKLAM1)

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    shape1d_freq_p4   = (/wc%nfreqs + 4/)
    shape1d_dir_2     = (/wc%ndirs, 2/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(const_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &         model_type=get_my_process_name())


    ! flminfr_tab      const%flminfr_tab(jmax,nfreqs)
    cf_desc    = t_cf_var('flminfr_tab', 'm^2 Hz-1', 'minimum allowed energy level', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'flminfr_tab', const%flminfr_tab,            &
         & GRID_UNSTRUCTURED_CELL, ZA_FREQ_GENERIC, cf_desc, grib2_desc,  &
         & ldims=(/wc%jmax, wc%nfreqs/),                                  &
         & lrestart=.FALSE., loutput=.FALSE.)

    ! for discrete approximation of nonlinear transfer
    !
    ! ikp              const%ikp(nfreqs + 4)
    cf_desc    = t_cf_var('IKP', '-', 'IKP', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'IKP', const%IKP,                      &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

    ! ikp1             const%ikp1(nfreqs + 4)
    cf_desc    = t_cf_var('IKP1', '-', 'IKP1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'IKP1', const%IKP1,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

    ! ikm              const%ikm(nfreqs + 4)
    cf_desc    = t_cf_var('IKM', '-', 'IKM', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'IKM', const%IKM,                      &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

    ! ikm1             const%ikm1(nfreqs + 4)
    cf_desc    = t_cf_var('IKM1', '-', 'IKM1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'IKM1', const%IKM1,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

    ! k1w              const%k1w(ndirs,2)
    cf_desc    = t_cf_var('K1W', '-', 'K1W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'K1W', const%K1W,                      &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_dir_2)

    ! k2w              const%k2w(ndirs,2)
    cf_desc    = t_cf_var('K2W', '-', 'K2W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'K2W', const%K2W,                      &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_dir_2)

    ! k11w             const%k11w(ndirs,2)
    cf_desc    = t_cf_var('K11W', '-', 'K11W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'K11W', const%K11W,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_dir_2)

    ! k21w             const%k21w(ndirs,2)
    cf_desc    = t_cf_var('K21W', '-', 'K21W', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'K21W', const%K21W,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_dir_2)

    ! ja1              const%ja1(ndirs,2)
    cf_desc    = t_cf_var('JA1', '-', 'JA1', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'JA1', const%JA1,                      &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_dir_2)

    ! ja2              const%ja2(ndirs,2)
    cf_desc    = t_cf_var('JA2', '-', 'JA2', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'JA2', const%JA2,                      &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_dir_2)

    ! af11             const%af11(nfreqs + 4)
    cf_desc    = t_cf_var('AF11', '-', 'AF11', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'AF11', const%AF11,                    &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

    ! fklap            const%fklap(nfreqs + 4)
    cf_desc    = t_cf_var('FKLAP', '-', 'FKLAP', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'FKLAP', const%FKLAP,                  &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

    ! fklap1           const%fklap1(nfreqs + 4)
    cf_desc    = t_cf_var('FKLAP1', '-', 'FKLAP1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'FKLAP1', const%FKLAP1,                &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

    ! fklam            const%fklam(nfreqs + 4)
    cf_desc    = t_cf_var('FKLAM', '-', 'FKLAM', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'FKLAM', const%FKLAM,                  &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

    ! fklam1           const%fklam1(nfreqs + 4)
    cf_desc    = t_cf_var('FKLAM1', '-', 'FKLAM1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(const_list, 'FKLAM1', const%FKLAM1,                &
         & GRID_UNSTRUCTURED_EDGE, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.FALSE.,                       &
         & ldims=shape1d_freq_p4)

  END SUBROUTINE new_wave_state_const_list


  !>
  !! Variable list collecting diagnostic fields which are required
  !! for time integration (affect solution).
  !!
  SUBROUTINE new_wave_state_diag_dyn_list(p_patch, diag_dyn, diag_dyn_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_diag_dyn), INTENT(INOUT) :: diag_dyn
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: diag_dyn_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    CHARACTER(len=*), PARAMETER :: routine = modname//'::new_wave_state_diag_dyn_list'

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: nblks_c, nblks_e
    INTEGER :: jg
    INTEGER :: shape2d_c(2)
    INTEGER :: shape3d_freq_c(3), shape3d_freq_e(3)

    TYPE(t_wave_config),      POINTER :: wc

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(diag_dyn%gv_c, &
    &       diag_dyn%gv_e, &
    &       diag_dyn%emean, &
    &       diag_dyn%emeanws, &
    &       diag_dyn%femean, &
    &       diag_dyn%hrms_frac, &
    &       diag_dyn%wbr_frac, &
    &       diag_dyn%wave_num_c, &
    &       diag_dyn%wave_num_e, &
    &       diag_dyn%f1mean, &
    &       diag_dyn%femeanws, &
    &       diag_dyn%akmean, &
    &       diag_dyn%xkmean, &
    &       diag_dyn%last_prog_freq_ind, &
    &       diag_dyn%ustar, &
    &       diag_dyn%z0, &
    &       diag_dyn%tauw, &
    &       diag_dyn%phiaw, &
    &       diag_dyn%tauhf, &
    &       diag_dyn%phihf, &
    &       diag_dyn%tm1, &
    &       diag_dyn%tm2)

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    jg      = p_patch%id
    nblks_c = p_patch%nblks_c
    nblks_e = p_patch%nblks_e

    shape2d_c         = (/nproma, nblks_c/)
    shape3d_freq_c    = (/nproma, wc%nfreqs, nblks_c/)
    shape3d_freq_e    = (/nproma, wc%nfreqs, nblks_e/)


    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(diag_dyn_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &         model_type=get_my_process_name())


    ! Wave group velocity
    cf_desc    = t_cf_var('gv_c', 'm s-1', 'group velocity at cells', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'gv_c', diag_dyn%gv_c,                   &
         & GRID_UNSTRUCTURED_CELL, ZA_FREQ_GENERIC, cf_desc, grib2_desc, &
         & ldims=shape3d_freq_c, in_group=groups("wave_phy_ext"),        &
         & lrestart=.FALSE., loutput=.TRUE.)

    cf_desc    = t_cf_var('gv_e', 'm s-1', 'group velocity at edges', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var(diag_dyn_list, 'gv_e', diag_dyn%gv_e,                   &
         & GRID_UNSTRUCTURED_EDGE, ZA_FREQ_GENERIC, cf_desc, grib2_desc, &
         & ldims=shape3d_freq_e, in_group=groups("wave_phy_ext"),        &
         & lrestart=.FALSE., loutput=.TRUE.)

    ! Wave physics group
    cf_desc    = t_cf_var('emean', 'm^2', 'total wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'emean', diag_dyn%emean,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('emeanws', 'm^2', 'wind sea wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'emeanws', diag_dyn%emeanws,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('femean', 's-1', 'mean frequency (m0/m-1)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'femean', diag_dyn%femean,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('hrms_frac', '-', 'square ratio (Hrms / Hmax)^2', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'hrms_frac', diag_dyn%hrms_frac,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c)

    cf_desc    = t_cf_var('wbr_frac', '-', 'fraction of breaking waves', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'wbr_frac', diag_dyn%wbr_frac,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('wave_num_c', '1/m', 'wave number at cell center', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'wave_num_c', diag_dyn%wave_num_c,       &
         & GRID_UNSTRUCTURED_CELL, ZA_FREQ_GENERIC, cf_desc, grib2_desc, &
         & ldims=shape3d_freq_c,                                         &
         & lrestart=.FALSE., loutput=.TRUE.)

    cf_desc    = t_cf_var('wave_num_e', 'm-1', 'wave number at edge midpoint', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_EDGE)
    CALL add_var(diag_dyn_list, 'wave_num_e', diag_dyn%wave_num_e,       &
         & GRID_UNSTRUCTURED_EDGE, ZA_FREQ_GENERIC, cf_desc, grib2_desc, &
         & ldims=shape3d_freq_e,                                         &
         & lrestart=.FALSE., loutput=.TRUE.)

    cf_desc    = t_cf_var('f1mean', 's-1', 'mean frequency based on f-moment (m1/m0)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'f1mean', diag_dyn%f1mean,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('femeanws', 's-1', 'wind sea mean frequency (m0/m-1)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'femeanws', diag_dyn%femeanws,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('akmean', '', 'Mean wavenumber based on SQRT(1/K)-moment, wm1', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'akmean', diag_dyn%akmean,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('xkmean', '', 'Mean wavenumber based on SQRT(K)-moment, wm2', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'xkmean', diag_dyn%xkmean,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))


    cf_desc    = t_cf_var('last_prog_freq_ind', '-', 'last frequency index of the prognostic range', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'last_prog_freq_ind', diag_dyn%last_prog_freq_ind, &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,            &
         & lrestart=.FALSE., loutput=.TRUE.,                                   &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('friction_velocity', 'm s-1', 'friction velocity', datatype_flt)
    grib2_desc = grib2_var(10, 0, 17, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'ustar', diag_dyn%ustar,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c , in_group=groups("wave_phy", "DWD_FG_WAVE_VARS"))

    cf_desc    = t_cf_var('roughness_length', 'm', 'Surface roughness length', datatype_flt)
    grib2_desc = grib2_var(2, 0, 1, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'z0', diag_dyn%z0,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c , in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('wave_stress', 'm2 s-2', 'wave stress divided by air density', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'tauw', diag_dyn%tauw,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.TRUE., loutput=.TRUE.,                         &
         & ldims=shape2d_c , in_group=groups("wave_phy","wave_fluxes"))

    cf_desc    = t_cf_var('integrated_energy_flux', 'W m-2', &
      &                   'integrated energy flux from wind into waves', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'phiaw', diag_dyn%phiaw,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c , in_group=groups("wave_phy","wave_fluxes"))

    cf_desc    = t_cf_var('tauhf', 'm2 s-2', 'high-fequency stress divided by air density', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'tauhf', diag_dyn%tauhf,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('phihf', 'W m-2', 'high-frequency energy flux into waves', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'phihf', diag_dyn%phihf,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    cf_desc    = t_cf_var('tm1', 's', 'Total m1 wave period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 25, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'tm1', diag_dyn%tm1,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    cf_desc    = t_cf_var('tm2', 's', 'Total m2 wave period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 28, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_dyn_list, 'tm2', diag_dyn%tm2,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

  END SUBROUTINE new_wave_state_diag_dyn_list


  !>
  !! Variable list of output diagnostics for scientific analysis (do not impact solution).
  !!
  SUBROUTINE new_wave_state_diag_out_list(p_patch, diag_out, diag_out_list, listname)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_diag_out), INTENT(INOUT) :: diag_out
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: diag_out_list
    CHARACTER(len=*),      INTENT(IN)    :: listname

    CHARACTER(len=*), PARAMETER :: routine = modname//'::new_wave_state_diag_out_list'

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: jf
    INTEGER :: ist
    INTEGER :: shape2d_c(2)
    INTEGER :: shape4d_c(4)
    INTEGER :: shape3d_dir_c(3)

    CHARACTER(len=3) :: freq_ind_str
    CHARACTER(len=VNAME_LEN) :: out_name

    TYPE(t_wave_config),      POINTER :: wc

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(diag_out%swell_mask,   &
      &     diag_out%beta,         &
      &     diag_out%drag,         &
      &     diag_out%ds,           &
      &     diag_out%hs_dir,       &
      &     diag_out%kp,           &
      &     diag_out%tauwn,        &
      &     diag_out%tmp,          &
      &     diag_out%tpp,          &
      &     diag_out%ds_sea,       &
      &     diag_out%emean_sea,    &
      &     diag_out%femean_sea,   &
      &     diag_out%f1mean_sea,   &
      &     diag_out%hs_sea,       &
      &     diag_out%hs_sea_dir,   &
      &     diag_out%kp_sea,       &
      &     diag_out%mp_sea,       &
      &     diag_out%m1_sea,       &
      &     diag_out%m2_sea,       &
      &     diag_out%pp_sea,       &
      &     diag_out%ds_swell,     &
      &     diag_out%emean_swell,  &
      &     diag_out%femean_swell, &
      &     diag_out%f1mean_swell, &
      &     diag_out%hs_swell,     &
      &     diag_out%hs_swell_dir, &
      &     diag_out%kp_swell,     &
      &     diag_out%mp_swell,     &
      &     diag_out%m1_swell,     &
      &     diag_out%m2_swell,     &
      &     diag_out%pp_swell,     &
      &     diag_out%steepness,    &
      &     diag_out%thp_adj,      &
      &     diag_out%sigma_f,      &
      &     diag_out%sigma_th,     &
      &     diag_out%qp_goda,      &
      &     diag_out%nu_f_LH,      &
      &     diag_out%bfis,         &
      &     diag_out%relw_fth,     &
      &     diag_out%kurtosis,     &
      &     diag_out%hmaxn,        &
      &     diag_out%hmax,         &
      &     diag_out%Tmax          &
      &    )

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    shape2d_c      = (/nproma, p_patch%nblks_c/)
    shape4d_c      = (/nproma, wc%ndirs, wc%nfreqs, p_patch%nblks_c/)
    shape3d_dir_c  = (/nproma, wc%ndirs, p_patch%nblks_c/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(diag_out_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &         model_type=get_my_process_name())

    ! beta             diag_out%beta(nproma,nblks_c)
    cf_desc    = t_cf_var('beta', '-', 'Charnock parameter', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'beta', diag_out%beta,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_debug"))

    ! drag             diag_out%drag(nproma,nblks_c)
    cf_desc    = t_cf_var('drag', '-', 'Drag coefficient', datatype_flt)
    grib2_desc = grib2_var(10, 0, 16, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'drag', diag_out%drag,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_debug"))

    ! ds               diag_out%ds(nproma,nblks_c)
    cf_desc    = t_cf_var('ds', 'deg', 'Total directional wave spread', datatype_flt)
    grib2_desc = grib2_var(10, 0, 31, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'ds', diag_out%ds,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! hs_dir           diag_out%hs_dir(nproma,nblks_c)
    cf_desc    = t_cf_var('hs_dir', 'deg', 'Total mean wave direction', datatype_flt)
    grib2_desc = grib2_var(10, 0, 14, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'hs_dir', diag_out%hs_dir,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! kp               diag_out%kp(nproma,nblks_c)
    cf_desc    = t_cf_var('kp', 'm-1', 'Total wave peak wavenumber', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'kp', diag_out%kp,                  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! tauwn            diag_out%tauwn(nproma,nblks_c)
    cf_desc    = t_cf_var('tauwn', '-', 'Normalized wave stress', datatype_flt)
    grib2_desc = grib2_var(10, 0, 19, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'tauwn', diag_out%tauwn,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_debug","wave_fluxes"))

    ! tmp              diag_out%tmp(nproma,nblks_c)
    cf_desc    = t_cf_var('tmp', 's', 'Total wave mean period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 15, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'tmp', diag_out%tmp,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! tpp              diag_out%tpp(nproma,nblks_c)
    cf_desc    = t_cf_var('tpp', 's', 'Total wave peak period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 34, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'tpp', diag_out%tpp,                &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    !
    ! wind sea
    !

    ! ds_sea           diag_out%ds_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('ds_sea', 'deg', 'Sea directional wave spread', datatype_flt)
    grib2_desc = grib2_var(10, 0, 32, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'ds_sea', diag_out%ds_sea,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c)

    ! emean_sea        diag_out%emean_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('emean_sea', 'm^2', 'Wind sea wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'emean_sea', diag_out%emean_sea,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    ! femean_sea       diag_out%femean_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('femean_sea', 's-1', 'Wind sea mean frequency (m0/m-1)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'femean_sea', diag_out%femean_sea,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    ! f1mean_sea       diag_out%f1mean_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('f1mean_sea', 's-1', 'Wind sea mean frequency (m1/m0)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'f1mean_sea', diag_out%f1mean_sea,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    ! hs_sea           diag_out%hs_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('hs_sea', 'm', 'Sea significant wave height', datatype_flt)
    grib2_desc = grib2_var(10, 0, 5, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'hs_sea', diag_out%hs_sea,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! hs_sea_dir       diag_out%hs_sea_dir(nproma,nblks_c)
    cf_desc    = t_cf_var('hs_sea_dir', 'deg', 'Sea mean wave direction', datatype_flt)
    grib2_desc = grib2_var(10, 0, 4, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'hs_sea_dir', diag_out%hs_sea_dir,  &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! kp_sea          diag_out%kp_sea(nproma,nblks_c)
    cf_desc   = t_cf_var('kp_sea', 'm-1', 'Sea wave peak wavenumber', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'kp_sea', diag_out%kp_sea,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! mp_sea           diag_out%mp_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('mp_sea', 's', 'Sea wave mean period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 6, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'mp_sea', diag_out%mp_sea,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! m1_sea           diag_out%m1_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('m1_sea', 's', 'Sea m1 wave period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'm1_sea', diag_out%m1_sea,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! m2_sea           diag_out%m2_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('m2_sea', 's', 'Sea m2 wave period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'm2_sea', diag_out%m2_sea,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! pp_sea           diag_out%pp_sea(nproma,nblks_c)
    cf_desc    = t_cf_var('pp_sea', 's', 'Sea wave peak period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 35, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'pp_sea', diag_out%pp_sea,          &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))


    !
    ! swell
    !

    ! ds_swell         diag_out%ds_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('ds_swell', 'deg', 'Swell directional wave spread', datatype_flt)
    grib2_desc = grib2_var(10, 0, 33, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'ds_swell', diag_out%ds_swell,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! emean_swell      diag_out%emean_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('emean_swell', 'm^2', 'Swell wave energy', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'emean_swell', diag_out%emean_swell, &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,  &
         & lrestart=.FALSE., loutput=.TRUE.,                         &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    ! femean_swell     diag_out%femean_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('femean_swell', 's-1', 'Swell mean frequency (m0/m-1)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'femean_swell', diag_out%femean_swell, &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
         & lrestart=.FALSE., loutput=.TRUE.,                           &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    ! f1mean_swell     diag_out%f1mean_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('f1mean_swell', 's-1', 'Swell mean frequency (m1/m0)', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'f1mean_swell', diag_out%f1mean_swell, &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
         & lrestart=.FALSE., loutput=.TRUE.,                           &
         & ldims=shape2d_c, in_group=groups("wave_phy"))

    ! hs_swell         diag_out%hs_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('hs_swell', 'm', 'Swell significant wave height', datatype_flt)
    grib2_desc = grib2_var(10, 0, 8, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'hs_swell', diag_out%hs_swell,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! hs_swell_dir     diag_out%hs_swell_dir(nproma,nblks_c)
    cf_desc    = t_cf_var('hs_swell_dir', 'deg', 'Swell mean wave direction', datatype_flt)
    grib2_desc = grib2_var(10, 0, 7, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'hs_swell_dir', diag_out%hs_swell_dir, &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
         & lrestart=.FALSE., loutput=.TRUE.,                           &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! kp_swell        diag_out%kp_swell(nproma,nblks_c)
    cf_desc   = t_cf_var('kp_swell', 'm-1', 'Swell wave peak wavenumber', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'kp_swell', diag_out%kp_swell,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! mp_swell         diag_out%mp_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('mp_swell', 's', 'Swell wave mean period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 9, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'mp_swell', diag_out%mp_swell,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! m1_swell         diag_out%m1_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('m1_swell', 's', 'Swell m1 wave period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'm1_swell', diag_out%m1_swell,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! m2_swell         diag_out%m2_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('m2_swell', 's', 'Swell m2 wave period', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'm2_swell', diag_out%m2_swell,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! pp_swell         diag_out%pp_swell(nproma,nblks_c)
    cf_desc    = t_cf_var('pp_swell', 's', 'Swell wave peak period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 36, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'pp_swell', diag_out%pp_swell,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_short"))

    ! swell_mask       diag_out%swell_mask(nproma,ndirs,nfreqs,nblks_c)
    cf_desc    = t_cf_var('swell_mask', '-', 'swell mask for tracers', datatype_int)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'swell_mask', diag_out%swell_mask,&
         & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC, cf_desc, grib2_desc, &
         & lcontainer=.TRUE., lrestart=.FALSE., loutput=.FALSE.,    &
         & ldims=shape4d_c)

    ALLOCATE(diag_out%swmask_ptr(wc%nfreqs), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(routine,               &
          'allocation of swmask_ptr failed')

    DO jf = 1, wc%nfreqs
      write(freq_ind_str,'(I0.3)') jf

      out_name = 'swmask_f'//TRIM(freq_ind_str)
      CALL add_ref(diag_out_list, 'swell_mask',                             &
           & out_name, diag_out%swmask_ptr(jf)%p,                           &
           & GRID_UNSTRUCTURED_CELL, ZA_DIR_GENERIC,                        &
           & t_cf_var(out_name, '-',out_name, datatype_int),                &
           & grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL), &
           & opt_var_ref_pos=3, ref_idx=jf, ldims=shape3d_dir_c,            &
           & lrestart=.TRUE., loutput=.TRUE.,                               &
           & in_group=groups("wave_debug","DWD_FG_WAVE_VARS"))
    END DO


    !
    ! extreme wave parameters
    !

    ! steepness        diag_out%steepness(nproma,nblks_c)
    cf_desc    = t_cf_var('steepness', '-', 'Wave steepness', datatype_flt)
    grib2_desc = grib2_var(10, 0, 192, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'steepness', diag_out%steepness,    &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! thp_adj          diag_out%thp_adj(nproma,nblks_c)
    cf_desc    = t_cf_var('thp_adj', '-', 'Adjusted peak direction', datatype_flt)
    grib2_desc = grib2_var(10, 0, 46, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'thp_adj', diag_out%thp_adj,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! sigma_f          diag_out%sigma_f(nproma,nblks_c)
    cf_desc    = t_cf_var('sigma_f', '-', 'Frequency bandwith', datatype_flt)
    grib2_desc = grib2_var(255,255,255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'sigma_f', diag_out%sigma_f,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! sigma_th         diag_out%sigma_th(nproma,nblks_c)
    cf_desc    = t_cf_var('sigma_th', '-', 'Directional bandwith', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'sigma_th', diag_out%sigma_th,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! qp_goda         diag_out%qp_goda(nproma,nblks_c)
    cf_desc    = t_cf_var('qp_goda', '-', 'Goda peakedness parameter', datatype_flt)
    grib2_desc = grib2_var(10, 0, 98, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'qp_goda', diag_out%qp_goda,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! nu_f_lh          diag_out%nu_f_lh(nproma,nblks_c)
    cf_desc    = t_cf_var('nu_f_LH', '-', 'Longuet-Higgins broadbandness parameter', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'nu_f_LH', diag_out%nu_f_LH,        &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! bfis             diag_out%bfis(nproma,nblks_c)
    cf_desc    = t_cf_var('bfis', '-', 'finite-depth Benjamin-Feir index', datatype_flt)
    grib2_desc = grib2_var(10, 0, 44, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'bfis', diag_out%bfis,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! relw_fth         diag_out%relw_fth(nproma,nblks_c)
    cf_desc    = t_cf_var('relw_fth', '-', 'Relative direction-frequency width', datatype_flt)
    grib2_desc = grib2_var(10, 0, 80, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'relw_fth', diag_out%relw_fth,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! kurtosis         diag_out%kurtosis(nproma,nblks_c)
    cf_desc    = t_cf_var('kurtosis', '-', 'Spectral kurtosis', datatype_flt)
    grib2_desc = grib2_var(10, 0, 43, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'kurtosis', diag_out%kurtosis,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! hmaxn            diag_out%hmaxn(nproma,nblks_c)
    cf_desc    = t_cf_var('hmaxn', '-', 'Normalised max. significant wave height', datatype_flt)
    grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'hmaxn', diag_out%hmaxn,            &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! The grib2 reference is for the envelope-max individual wave height
    ! hmax             diag_out%hmax(nproma,nblks_c)
    cf_desc    = t_cf_var('hmax', 'm', 'Max. significant wave height', datatype_flt)
    grib2_desc = grib2_var(10, 0, 93, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'hmax', diag_out%hmax,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

    ! Tmax             diag_out%Tmax(nproma,nblks_c)
    cf_desc    = t_cf_var('Tmax', 's', 'Max. wave period', datatype_flt)
    grib2_desc = grib2_var(10, 0, 23, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_out_list, 'Tmax', diag_out%Tmax,              &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_extreme"))

  END SUBROUTINE new_wave_state_diag_out_list


  !>
  !! Variable list of diagnostics for coupling and output, if requested.
  !!
  SUBROUTINE new_wave_state_diag_cpl_list(p_patch, diag_cpl, diag_cpl_list, listname, var_in_output)

    TYPE(t_patch),         INTENT(IN)    :: p_patch
    TYPE(t_wave_diag_cpl), INTENT(INOUT) :: diag_cpl
    TYPE(t_var_list_ptr),  INTENT(INOUT) :: diag_cpl_list
    CHARACTER(len=*),      INTENT(IN)    :: listname
    TYPE(t_wave_var_in_output), INTENT(IN)    :: &  !< optional diagnostic switches
      &  var_in_output

    CHARACTER(len=*), PARAMETER :: routine = modname//'::new_wave_state_diag_cpl_list'

    TYPE(t_cf_var)    :: cf_desc
    TYPE(t_grib2_var) :: grib2_desc

    INTEGER :: ibits         !< "entropy" of horizontal slice
    INTEGER :: datatype_flt  !< floating point accuracy in NetCDF output
    INTEGER :: shape2d_c(2)
    INTEGER :: shape3d_depth_c(3)

    TYPE(t_wave_config),      POINTER :: wc

    !------------------------------
    ! Ensure that all pointers have a defined association status
    !------------------------------
    NULLIFY(diag_cpl%hs,             &
      &     diag_cpl%u_stokes,       &
      &     diag_cpl%v_stokes,       &
      &     diag_cpl%last_idx_depth, &
      &     diag_cpl%u3d_stokes,     &
      &     diag_cpl%v3d_stokes,     &
      &     diag_cpl%tauoc_x,        &
      &     diag_cpl%tauoc_y,        &
      &     diag_cpl%tauoc,          &
      &     diag_cpl%phioc           &
      &    )

    ! pointer to wave_config(jg) to save some paperwork
    wc => wave_config(p_patch%id)

    shape2d_c       = (/nproma, p_patch%nblks_c/)
    shape3d_depth_c = (/nproma, wc%oce_stokes_nlev, p_patch%nblks_c/)

    ibits = DATATYPE_PACK16   ! "entropy" of horizontal slice

    IF ( lnetcdf_flt64_output ) THEN
      datatype_flt = DATATYPE_FLT64
    ELSE
      datatype_flt = DATATYPE_FLT32
    ENDIF

    CALL vlr_add(diag_cpl_list, TRIM(listname), patch_id=p_patch%id, lrestart=.TRUE., &
      &         model_type=get_my_process_name())


    IF (is_coupled_to_ocean() .OR. var_in_output%hs) THEN
      ! hs               diag_cpl%hs(nproma,nblks_c)
      cf_desc    = t_cf_var('hs', 'm', 'Total significant wave height', datatype_flt)
      grib2_desc = grib2_var(10, 0, 3, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(diag_cpl_list, 'hs', diag_cpl%hs,                  &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
           & hor_interp=create_hor_interp_metadata(                   &
           &    hor_intp_type=HINTP_TYPE_LONLAT_BCTR,                 &
           &    fallback_type=HINTP_TYPE_LONLAT_NNB),                 &
           & lrestart=.FALSE., loutput=.TRUE.,                        &
           & ldims=shape2d_c, in_group=groups("wave_short"))
    ENDIF

    cf_desc    = t_cf_var('u_stokes', 'ms-1', 'U-component surface Stokes drift', datatype_flt)
    grib2_desc = grib2_var(10, 0, 21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_cpl_list, 'u_stokes', diag_cpl%u_stokes,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_stokes"))

    cf_desc    = t_cf_var('v_stokes', 'ms-1', 'V-component surface Stokes drift', datatype_flt)
    grib2_desc = grib2_var(10, 0, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var(diag_cpl_list, 'v_stokes', diag_cpl%v_stokes,      &
         & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
         & lrestart=.FALSE., loutput=.TRUE.,                        &
         & ldims=shape2d_c, in_group=groups("wave_stokes"))


    IF (var_in_output%last_idx_depth .OR. &
      & var_in_output%u3d_stokes     .OR. &
      & var_in_output%v3d_stokes)   THEN

      IF (TRIM(wc%oce_vct_filename) == "") THEN
        CALL finish(routine, "file name for ocean vertical interfaces table (namelist oce_vct_filename) is not defined")
      END IF

      cf_desc    = t_cf_var('last_idx_depth', '-', 'last index of depth layer', datatype_int)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(diag_cpl_list, 'last_idx_depth', diag_cpl%last_idx_depth, &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
           & lrestart=.FALSE., loutput=.TRUE.,                           &
           & ldims=shape2d_c)

      cf_desc    = t_cf_var('u3d_stokes', 'ms-1', 'U-component of 3d Stokes drift', datatype_flt)
      grib2_desc = grib2_var(10, 0, 21, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(diag_cpl_list, 'u3d_stokes', diag_cpl%u3d_stokes,          &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA, cf_desc, grib2_desc, &
           & lrestart=.FALSE., loutput=.TRUE.,                                &
           & ldims=shape3d_depth_c)

      cf_desc    = t_cf_var('v3d_stokes', 'ms-1', 'V-component of 3d Stokes drift', datatype_flt)
      grib2_desc = grib2_var(10, 0, 22, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(diag_cpl_list, 'v3d_stokes', diag_cpl%v3d_stokes,          &
           & GRID_UNSTRUCTURED_CELL, ZA_DEPTH_BELOW_SEA, cf_desc, grib2_desc, &
           & lrestart=.FALSE., loutput=.TRUE.,                                &
           & ldims=shape3d_depth_c)
    END IF

    IF (is_coupled_to_ocean() .OR. var_in_output%tauoc_x    &
      &                       .OR. var_in_output%tauoc_y    &
      &                       .OR. var_in_output%tauoc      &
      &                       .OR. var_in_output%phioc)  THEN

      ! wave stress normalised by roair: bare x wave stress has code (10,0,90) and units N/m^2
      cf_desc    = t_cf_var('tauoc_x', '(m/s)^2', 'Wave-to-ocean stress x-component', datatype_int)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(diag_cpl_list, 'tauoc_x', diag_cpl%tauoc_x,           &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
           & lrestart=.FALSE., loutput=.TRUE.,                           &
           & ldims=shape2d_c, in_group=groups("wave_fluxes"))

      ! wave stress normalised by roair: bare y stress has code (10,0,91) and units N/m^2
      cf_desc    = t_cf_var('tauoc_y', '(m/s)^2', 'Wave-to-ocean stress y-component', datatype_int)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(diag_cpl_list, 'tauoc_y', diag_cpl%tauoc_y,           &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc,    &
           & lrestart=.FALSE., loutput=.TRUE.,                           &
           & ldims=shape2d_c, in_group=groups("wave_fluxes"))

      ! wave stress normalised by roair:  normalised stress has code (10,0,84)
      cf_desc    = t_cf_var('tauoc', '(m/s)^2', 'Wave-to-ocean stress', datatype_int)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(diag_cpl_list, 'tauoc', diag_cpl%tauoc,            &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
           & lrestart=.FALSE., loutput=.TRUE.,                        &
           & ldims=shape2d_c, in_group=groups("wave_fluxes"))

      ! wave-to-ocean energy flux divided by roair: normalised flux has code (10,0,85)
      cf_desc    = t_cf_var('phioc', 'kg/s^3', 'Wave-to-ocean energy flux', datatype_int)
      grib2_desc = grib2_var(255, 255, 255, ibits, GRID_UNSTRUCTURED, GRID_CELL)
      CALL add_var(diag_cpl_list, 'phioc', diag_cpl%phioc,            &
           & GRID_UNSTRUCTURED_CELL, ZA_SURFACE, cf_desc, grib2_desc, &
           & lrestart=.FALSE., loutput=.TRUE.,                        &
           & ldims=shape2d_c, in_group=groups("wave_fluxes"))
    END IF

  END SUBROUTINE new_wave_state_diag_cpl_list


  !>
  !! Destruction of wave-specific variable lists and memory deallocation
  !!
  SUBROUTINE destruct_wave_state ()

    INTEGER :: ist
    INTEGER :: jg, jt
    CHARACTER(len=*), PARAMETER :: routine = modname//'::destruct_wave_state'

    DO jg = 1, n_dom
      ! delete prognostic state list elements
      DO jt = 1, SIZE(p_wave_state_lists(jg)%prog_list(:))
        CALL vlr_del(p_wave_state_lists(jg)%prog_list(jt))
      ENDDO

      ! delete source state list elements
      CALL vlr_del(p_wave_state_lists(jg)%source_list)

      ! delete const state list elements
      CALL vlr_del(p_wave_state_lists(jg)%const_list)

      ! delete diagnostics state list elements
      CALL vlr_del(p_wave_state_lists(jg)%diag_dyn_list)
      CALL vlr_del(p_wave_state_lists(jg)%diag_out_list)
      CALL vlr_del(p_wave_state_lists(jg)%diag_cpl_list)

      ! deallocate state lists and arrays
      DEALLOCATE(p_wave_state_lists(jg)%prog_list, stat=ist)
      IF (ist /= success) THEN
        CALL finish(routine,'deallocation for prog_list array failed')
      END IF

      DEALLOCATE(p_wave_state(jg)%prog, stat=ist)
      IF (ist /= success) THEN
        CALL finish(routine,'deallocation of prognostic state array failed')
      END IF
    ENDDO

    ! deallocate states
    DEALLOCATE(p_wave_state, p_wave_state_lists, stat=ist)
    IF (ist /= success) THEN
      CALL finish(routine,'deallocation for p_wave_state failed')
    END IF

  END SUBROUTINE destruct_wave_state


END MODULE mo_wave_state

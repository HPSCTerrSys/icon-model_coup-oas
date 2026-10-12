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

! Waves config.
! The content is mostly adopted from the WAM 4.5.4.

MODULE mo_wave_config

  USE mo_kind,                 ONLY: wp, vp
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_impl_constants,       ONLY: max_dom, SUCCESS, MAX_CHAR_LENGTH
  USE mo_math_constants,       ONLY: pi2, rad2deg, dbl_eps
  USE mo_physical_constants,   ONLY: grav, rhoh2o
  USE mo_wave_constants,       ONLY: EX_TAIL
  USE mo_fortran_tools,        ONLY: DO_DEALLOCATE
  USE mo_io_units,             ONLY: filename_max, find_next_free_unit
  USE mo_util_string,          ONLY: t_keyword_list, associate_keyword, with_keywords, &
    &                                int2string
  USE mo_mpi,                  ONLY: p_io, p_comm_work, my_process_is_stdio, &
    &                                p_bcast, p_comm_work_test
  USE mo_parallel_config,      ONLY: p_test_run
  USE mo_run_config,           ONLY: msg_level
  USE mtime,                   ONLY: datetime

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_config'

  ! types
  PUBLIC :: t_wave_config

  ! objects
  PUBLIC :: wave_config

  ! subroutines
  PUBLIC :: configure_wave
  PUBLIC :: generate_filename

  TYPE t_wave_config
    INTEGER  :: ndirs    ! number of directions.
    INTEGER  :: nfreqs   ! number of frequencies.

    INTEGER  :: Tlength  ! duration of sample time series

    REAL(wp) :: fr1      ! first frequency [hz].
    REAL(wp) :: co       ! frequency ratio
    INTEGER  :: iref     ! frequency bin number of reference frequency

    REAL(wp) :: umax     ! maximum wind speed in flminfr table
    INTEGER  :: jmax     ! number of wind speed bins in flminfr table
    REAL(wp) :: flmin    ! minimum of energy in flminfr table
    REAL(wp) :: delu     ! wind speed increment in flminfr table

    REAL(wp) :: alpha      ! phillips' parameter  (not used if iopti = 1)
    REAL(wp) :: fm         ! peak frequency (hz) and/or maximum frequency
    REAL(wp) :: gamma_wave ! overshoot factor
    REAL(wp) :: sigma_a    ! left peak width
    REAL(wp) :: sigma_b    ! right peak width
    REAL(wp) :: fetch             ! fetch in metres used for initialisation of spectrum
    REAL(wp) :: fetch_min_energy  ! fetch in meters used for calculation of minimum allowed energy level

    REAL(wp) :: roair   ! air density
    REAL(wp) :: rnuair  ! kinematic air viscosity
    REAL(wp) :: rnuairm ! kinematic air viscosity for momentum transfer
    REAL(wp) :: rowater ! water density
    REAL(wp) :: xeps
    REAL(wp) :: xinveps

    REAL(wp) :: xkappa   ! von karman constant.
    REAL(wp) :: xnlev    ! windspeed ref. level.
    REAL(wp) :: betamax  ! parameter for wind input (ecmwf cy45r1).
    REAL(wp) :: zalp     ! shifts growth curve (ecmwf cy45r1).
    REAL(wp) :: alpha_ch ! minimum charnock constant (ecmwf cy45r1)

    CHARACTER(LEN=filename_max) :: oce_vct_filename ! name of ocean vertical coordinates file

    REAL(wp) :: depth        ! ocean depth (m) if not 0, then constant depth
    REAL(wp) :: depth_min    ! allowed minimum of model depth (m)
    REAL(wp) :: depth_max    ! allowed maximum of model depth (m)
    INTEGER  :: stokes_method ! 1 - calculation of Stokes profile from the full spectrum
                              ! 2 - calculation based on Breivik (2016)
    REAL(wp) :: stokes_depth ! maximum of Stokes layer depth (m)
    REAL(wp) :: stokes_th    ! individual Stokes layer thickness (m)

    INTEGER  :: niter_smooth ! number of smoothing iterations for wave bathymetry
                             ! if 0 then no smoothing
    INTEGER  :: nsubs_refrac ! number of susteps in wave refraction

    INTEGER  :: jtot_tauhf ! dimension of wtauhf, must be odd
    REAL(wp) :: x0tauhf    ! lowest limit for integration in tau_phi_hf: x0 *(g/ustar)

    CHARACTER(LEN=filename_max) :: forc_file_prefix ! prefix of forcing file name
                                           ! the real file name will be constructed as:
                                           ! forc_file_prefix+'_wind.nc' for U and V 10 meter wind (m/s)
                                           ! forc_file_prefix+'_ice.nc'  for sea ice concentration (fraction of 1)
                                           ! forc_file_prefix+'_slh.nc'  for sea level height (m)
                                           ! forc_file_prefix+'_osc.nc'  for U and V ocean surface currents (m/s)

    LOGICAL :: linput_sf1      ! if .TRUE., calculate wind input source function term, first call
    LOGICAL :: linput_sf2      ! if .TRUE., calculate wind input source function term, second call
    LOGICAL :: ldissip_sf      ! if .TRUE., calculate dissipation source function term
    LOGICAL :: lwave_brk_sf    ! if .TRUE., calculate wave breaking dissipation source function term
    LOGICAL :: lnon_linear_sf  ! if .TRUE., calculate non linear source function term
    LOGICAL :: lbottom_fric_sf ! if .TRUE., calculate bottom_friction source function term
    LOGICAL :: lwave_stress1   ! if .TRUE., calculate wave stress, first call
    LOGICAL :: lwave_stress2   ! if .TRUE., calculate wave stress, second call

    REAL(wp) :: peak_lat, peak_lon ! geographical location (deg) of wind peak value
    REAL(wp) :: peak_wsp10         ! peak value of 10m wind spead at (peak_lat, peak_lon)
    REAL(wp) :: dir_wsp10          ! flow direction measured clockwise from true north (deg)
    REAL(wp) :: r_wsp10            ! physical radius of wind field (m)

    REAL(wp) :: impl_fac       !! implicitness factor for total source function time integration
                               !! impl_fac=0.5 : second order Crank-Nicholson/trapezoidal scheme
                               !! impl_fac=1   : first order Euler backward scheme
                               !! valid range: 0.5 <= impl_fac <= 1

    ! derived variables and fields
    !

    REAL(wp) ::            &
      &  delth,            & ! angular increment of spectrum [rad].
      &  mo_tail,          & ! mo  tail factor.
      &  mm1_tail,         & ! m-1 tail factor.
      &  mp1_tail,         & ! m+1 tail factor.
      &  mp2_tail            ! m+2 tail factor.

    REAL(vp) ::            &
      &  acl1,             & ! weight in angular grid for interpolation,
                             ! wave no. 3 ("1+lambda" term).
      &  acl2,             & ! weight in angular grid for interpolation,
                             ! wave no. 4 ("1-lambda" term).
      &  cl11,             & ! 1.-acl1.
      &  cl21,             & ! 1.-acl2.
      &  dal1,             & ! 1./acl1.
      &  dal2,             & ! 1./acl2.
      &  frh(30)             ! tail frequency ratio **5

    REAL(wp), ALLOCATABLE :: &
      &  freqs(:),         & ! frequencies (1:nfreqs) of wave spectrum [hz]
      &  dfreqs(:),        & ! frequency interval (1:nfreqs)
      &  dfreqs_freqs(:),  & ! dfreqs * freqs
      &  dfreqs_freqs2(:), & ! dfreqs * freqs * freqs
      &  dirs(:),          & ! directions (1:ndirs) of wave spectrum [rad]
      &  dfim(:),          & ! mo  integration weights.
      &  dfimofr(:),       & ! m-1 integration weights.
      &  dfim_fr(:),       & ! m+1 integration weights.
      &  dfim_fr2(:),      & ! m+2 integration weights.
      &  sin_dir(:),       & ! sine of direction
      &  cos_dir(:),       & ! cosine of direction
      &  rhowg_dfim(:),    & ! momentum and energy flux weights.
      &  wtauhf(:)           ! integration weight for tau_phi_hf

     REAL(wp), ALLOCATABLE :: &
      &  oce_ifc(:),        & ! depth of ocean vertical cell interfaces (m)
      &  oce_mc(:),         & ! depth of ocean vertical cell midpoints (m)
      &  oce_lay_th(:),     & ! ocean cell thickness (m)
      &  oce_stokes_ifc(:), & ! depth of Stokes interfaces (m)
      &  oce_stokes_mc(:)     ! depth of Stokes midpoints (m)

    INTEGER, ALLOCATABLE :: &
      &  dir_neig_ind(:,:)   ! index of direction neighbor (2,1:ndirs)

    LOGICAL :: lread_forcing ! set to .TRUE. if a forcing file prefix has been specified (forc_file_prefix)

    ! For precomputed index list of coastal edges (see mo_wave_ext_data_init:init_coastedge_list, including allocation)
    INTEGER :: n_coastedges
    INTEGER, ALLOCATABLE :: idx_coastedges(:), blk_coastedges(:)

    REAL(wp), ALLOCATABLE :: orient_coastedges(:) ! corresponding edge orientation

    INTEGER :: oce_nifc ! number of ocean vertical interfaces
    INTEGER :: oce_nlev ! number of ocean vertical layers
    INTEGER :: oce_stokes_nifc ! number of Stokes interfaces
    INTEGER :: oce_stokes_nlev ! number of Stokes vertical layers

  CONTAINS
    !
    ! calculate directions for wave spectrum
    PROCEDURE :: set_directions      => wave_config__set_directions
    !
    ! calculate frequencies and intervals for wave spectrum
    PROCEDURE :: set_frequencies     => wave_config__set_frequencies
    !
    ! destruct wave_config object
    PROCEDURE :: destruct            => wave_config_destruct
  END type t_wave_config

  TYPE(t_wave_config), TARGET:: wave_config(max_dom)

CONTAINS

  !>
  !! deallocate memory used by object of type t_wave_config
  !!
  SUBROUTINE wave_config_destruct(me)
    CLASS(t_wave_config) :: me
    CHARACTER(*), PARAMETER :: routine = modname//'::wave_config_destruct'

    CALL DO_DEALLOCATE(me%freqs)
    CALL DO_DEALLOCATE(me%dfreqs)
    CALL DO_DEALLOCATE(me%dfreqs_freqs)
    CALL DO_DEALLOCATE(me%dfreqs_freqs2)
    CALL DO_DEALLOCATE(me%dirs)
    CALL DO_DEALLOCATE(me%sin_dir)
    CALL DO_DEALLOCATE(me%cos_dir)
    CALL DO_DEALLOCATE(me%DFIM)
    CALL DO_DEALLOCATE(me%DFIMOFR)
    CALL DO_DEALLOCATE(me%DFIM_FR)
    CALL DO_DEALLOCATE(me%DFIM_FR2)
    CALL DO_DEALLOCATE(me%RHOWG_DFIM)
    CALL DO_DEALLOCATE(me%dir_neig_ind)
    CALL DO_DEALLOCATE(me%wtauhf)
    CALL DO_DEALLOCATE(me%idx_coastedges)
    CALL DO_DEALLOCATE(me%blk_coastedges)
    CALL DO_DEALLOCATE(me%orient_coastedges)
    CALL DO_DEALLOCATE(me%oce_ifc)
    CALL DO_DEALLOCATE(me%oce_mc)
    CALL DO_DEALLOCATE(me%oce_lay_th)
    CALL DO_DEALLOCATE(me%oce_stokes_ifc)
    CALL DO_DEALLOCATE(me%oce_stokes_mc)


  END SUBROUTINE wave_config_destruct


  !>
  !! Read the number of ocean vertical cell interfaces from oce_vct_file
  !!
  SUBROUTINE read_oce_vct(oce_vct_file, oce_nifc, oce_ifc)

    CHARACTER(LEN=*),      INTENT(IN)    :: oce_vct_file
    INTEGER,               INTENT(INOUT) :: oce_nifc
    REAL(wp), ALLOCATABLE, INTENT(INOUT) :: oce_ifc(:)

    CHARACTER(*),PARAMETER :: routine = modname//'::read_oce_vct'

    INTEGER :: ist, iunit, jk, ik
    INTEGER :: mpi_comm

    IF (my_process_is_stdio()) THEN
      iunit = find_next_free_unit(10,20)
      OPEN (unit=iunit,file=TRIM(oce_vct_file),access='SEQUENTIAL', &
        &  form='FORMATTED', action='READ', status='OLD', IOSTAT=ist)
      IF (ist/=success) THEN
        CALL finish (routine, 'open vertical coordinate table file failed')
      END IF

      READ (iunit,*,IOSTAT=ist) oce_nifc

      IF (ist/=success) THEN

        CALL finish (routine, 'reading number of ocean interface levels failed')

      ELSE

        ALLOCATE(oce_ifc(oce_nifc), stat=ist)
        IF (ist/=SUCCESS) CALL finish(routine, "allocation for oce_ifc of type REAL failed")

        DO jk=1, oce_nifc

          READ (iunit,*,IOSTAT=ist) ik, oce_ifc(jk)
          IF (ist/=success) THEN
            CALL finish (routine, 'reading of ocean interface levels failed')
          END IF

        END DO

        CALL message(routine, 'ocean vertical coordinate table file successfully read')

      END IF

      CLOSE(iunit)
    END IF

    IF (p_test_run) THEN
      mpi_comm = p_comm_work_test
    ELSE
      mpi_comm = p_comm_work
    END IF

    CALL p_bcast(oce_nifc,  p_io, mpi_comm)
    IF (.NOT.my_process_is_stdio()) THEN
      ALLOCATE(oce_ifc(oce_nifc), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for oce_ifc of type REAL failed")
    END IF
    CALL p_bcast(oce_ifc,  p_io, mpi_comm)

  END SUBROUTINE read_oce_vct

  !>
  !! Calculate interfaces and midpoint depths of ocean levels and
  !! within the Stokes layer according to the ocean vertical coordinates table.
  !! The Stokes layer is limited by the namelist parameter stokes_depth,
  !! if given as > 0, or by the depth of the last ocean level.
  !!
  SUBROUTINE oce_stokes_levels(wave_config)

    TYPE(t_wave_config), TARGET, INTENT(INOUT) :: wave_config

    TYPE(t_wave_config), POINTER :: wc => NULL()

    INTEGER :: jk, ist

    CHARACTER(len=*), PARAMETER ::  &
      &  routine = modname//':oce_stokes_levels'

    wc => wave_config

    ! read ocean vertical coordinates table
    IF (TRIM(wc%oce_vct_filename) /= "") THEN

      CALL read_oce_vct(oce_vct_file = wc%oce_vct_filename, &
        &                   oce_nifc = wc%oce_nifc, &
        &                    oce_ifc = wc%oce_ifc)

      wc%oce_nlev = wc%oce_nifc-1 ! number of ocean levels

      ALLOCATE(wc%oce_mc(wc%oce_nlev),wc%oce_lay_th(wc%oce_nlev), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for ocean vct of type REAL failed")

      CALL message ('','')
      CALL message (':-----------------------------------------------------------','')
      WRITE(message_text,'(i5)') wc%oce_nlev
      CALL message (' Run with ocean layers, number of layers: ',message_text)
      CALL message (':-----------------------------------------------------------','')
      CALL message (' index, depth of center, depth of interface, layer thickness','')
      CALL message (':-----------------------------------------------------------','')

      DO jk=1, wc%oce_nlev

        wc%oce_mc(jk) = 0.5_wp * (wc%oce_ifc(jk) + wc%oce_ifc(jk+1))
        wc%oce_lay_th(jk) = wc%oce_ifc(jk+1) - wc%oce_ifc(jk)

        WRITE(message_text,'(i6,f17.1,f19.1,f15.1)') jk, wc%oce_mc(jk), wc%oce_ifc(jk), wc%oce_lay_th(jk)
        CALL message ('',message_text)

      END DO

      WRITE(message_text,'(i6,a17,f19.1,a16)') jk, '-',  wc%oce_ifc(jk), '-'
      CALL message ('',message_text)
      CALL message (':-----------------------------------------------------------','')

      ! calculate depth of each Stokes level (interface and midpoint) within wc%stokes_depth
      IF (wc%stokes_depth > 0._wp) THEN
        ! find the last interface index according to wc%stokes_depth
        wc%oce_stokes_nifc = MINLOC(ABS(wc%oce_ifc-wc%stokes_depth), dim=1)
      ELSE
        ! use full depth of ocean table
        wc%oce_stokes_nifc = wc%oce_nifc
      END IF

      wc%oce_stokes_nlev = wc%oce_stokes_nifc - 1

      ALLOCATE(wc%oce_stokes_ifc(wc%oce_stokes_nifc),wc%oce_stokes_mc(wc%oce_stokes_nlev), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for stokes_ifc field of type REAL failed")

      CALL message ('','')
      CALL message (':--- Run with Stokes layer ---------------------------------','')
      WRITE(message_text,'(f10.5)') wc%stokes_depth
      CALL message (' depth of Stokes layer (m)',message_text)
      CALL message (':-----------------------------------------------------------','')
      CALL message (' index, depth of center, depth of interface                 ','')
      CALL message (':-----------------------------------------------------------','')

      DO jk = 1, wc%oce_stokes_nifc
        wc%oce_stokes_ifc(jk) = wc%oce_ifc(jk)
      END DO

      DO jk = 1, wc%oce_stokes_nlev
        wc%oce_stokes_mc(jk) = wc%oce_mc(jk)

        WRITE(message_text,'(i6,f17.1,f19.1)') jk, wc%oce_stokes_mc(jk), wc%oce_stokes_ifc(jk)
        CALL message ('',message_text)
      END DO

      WRITE(message_text,'(i6,a17,f19.1)') jk, '-',  wc%oce_stokes_ifc(jk)
      CALL message ('',message_text)
      CALL message (':-----------------------------------------------------------','')

    END IF

  END SUBROUTINE oce_stokes_levels

  !>
  !! setup the waves model
  !!
  !! Setup of additional waves control variables and constant fields
  !! which depend on the waves-NAMELIST and potentially other namelists.
  !! This routine is called, after all namelists have been read and a
  !! synoptic consistency check has been done.
  !!
  SUBROUTINE configure_wave(n_dom)

    INTEGER, INTENT(IN) :: n_dom    !< number of domains

    INTEGER :: j, jd     ! loop index
    INTEGER :: jg        ! patch ID
    INTEGER :: ist       ! error status

    TYPE(t_wave_config), POINTER :: wc =>NULL()     ! convenience pointer

    REAL(wp) :: X0, FF, F, DF, CONST1

    CHARACTER(*), PARAMETER :: routine = modname//'::configure_waves'

    CALL message (' ','')
    CALL message ('!','')
    CALL message ('!!','')
    CALL message ('!!!!','')
    CALL message ('!!!!!!','')
    CALL message ('!!!!!!!!','')
    CALL message ('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!','')
    CALL message ('!     ICON-WAVES IS STILL UNDER DEVELOPMENT                              !','')
    CALL message ('!     THIS CODE IS FOR TECHNICAL TESTING ONLY                            !','')
    CALL message ('!     THE WAVE PHYSICS IS UNDER EVALUATION                               !','')
    CALL message ('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!','')
    CALL message ('!!!!!!!!','')
    CALL message ('!!!!!!','')
    CALL message ('!!!!','')
    CALL message ('!!','')
    CALL message ('!','')
    CALL message (' ','')


    DO jg=1,n_dom

      ! convenience pointer
      wc => wave_config(jg)


      ! reading of external forcing data yes/no
      ! set to .TRUE. if a forcing file prefix has been specified
      IF (TRIM(wave_config(jg)%forc_file_prefix) /= '') THEN
        wc%lread_forcing = .TRUE.
      ELSE
        wc%lread_forcing = .FALSE.
      ENDIF

      ALLOCATE(                             &
        &      wc%dfim         (wc%nfreqs), &
        &      wc%dfimofr      (wc%nfreqs), &
        &      wc%dfim_fr      (wc%nfreqs), &
        &      wc%dfim_fr2     (wc%nfreqs), &
        &      wc%rhowg_dfim   (wc%nfreqs), &
        &      wc%wtauhf   (wc%jtot_tauhf), &
        &      stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for fields of type REAL failed")

      ALLOCATE(wc%dir_neig_ind (2,wc%ndirs), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for fields of type INTEGER failed")

      ! calculate interface and midpoint depths of ocean levels and
      ! within the Stokes layer
      CALL oce_stokes_levels(wc)

      ! calculate wind speed interval in the flminfr table
      wc%delu = wc%umax/REAL(wc%jmax,wp)

      !
      ! configuration of spectral setup
      !
      CALL message ('  ','')
      CALL message (':----------------------------------------------------------','')
      WRITE(message_text,'(a,i4)') ' domain ', jg
      CALL message ('Frequencies and directions of wave spectrum',message_text)
      WRITE(message_text,'(a,i4)') 'Number of directions  = ', wc%ndirs
      CALL message ('  ',message_text)
      WRITE(message_text,'(a,i4)') 'Number of frequencies = ', wc%nfreqs
      CALL message ('  ',message_text)
      CALL message (':----------------------------------------------------------','')

      ! calculate directions for wave spectrum
      !
      CALL wc%set_directions()

      ! calculate frequencies and intervals for wave spectrum
      !
      CALL wc%set_frequencies()


      ! calculate direction neighbor index
      DO jd = 1,wc%ndirs
        IF (jd == 1) THEN
          wc%dir_neig_ind(1,jd) = wc%ndirs
          wc%dir_neig_ind(2,jd) = jd+1
        ELSEIF (jd == wc%ndirs) THEN
          wc%dir_neig_ind(1,jd) = jd-1
          wc%dir_neig_ind(2,jd) = 1
        ELSE
          wc%dir_neig_ind(1,jd) = jd-1
          wc%dir_neig_ind(2,jd) = jd+1
        END IF
      END DO


      ! MO  TAIL FACTOR.
      wc%MO_TAIL  = - wc%DELTH / (EX_TAIL + 1.0_wp) * wc%freqs(wc%nfreqs)

      ! M-1 TAIL FACTOR.
      wc%MM1_TAIL = - wc%DELTH / EX_TAIL

      ! M+1 TAIL FACTOR.
      wc%MP1_TAIL = - wc%DELTH / (EX_TAIL + 2.0_wp) * wc%freqs(wc%nfreqs)**2

      ! M+2 TAIL FACTOR.
      wc%MP2_TAIL = - wc%DELTH / (EX_TAIL + 3.0_wp) * wc%freqs(wc%nfreqs)**3


      ! calculate freqs_dirs parameters
      ! MO INTEGRATION WEIGHTS
      wc%DFIM = wc%dfreqs * wc%DELTH

      ! M-1 INTEGRATION WEIGHTS.
      wc%DFIMOFR = wc%DFIM / wc%freqs

      ! M+1 INTEGRATION WEIGHTS.
      wc%DFIM_FR = wc%dfreqs_freqs * wc%DELTH

      ! M+2 INTEGRATION WEIGHTS.
      wc%DFIM_FR2 = wc%dfreqs_freqs2 * wc%DELTH

      ! MOMENTUM AND ENERGY FLUX WEIGHTS.
      wc%RHOWG_DFIM(:) = rhoh2o * grav * wc%DELTH * LOG(wc%CO) * wc%freqs(:)
      !
      wc%RHOWG_DFIM(1)         = 0.5_wp * wc%RHOWG_DFIM(1)
      wc%RHOWG_DFIM(wc%nfreqs) = 0.5_wp * wc%RHOWG_DFIM(wc%nfreqs)

      ! initialisation of wtauhf
      wc%wtauhf(:) = 0._wp

      X0 = 0.005_wp
      FF = EXP(wc%XKAPPA / (X0 + wc%ZALP))
      F = wc%ALPHA_CH * X0**2 * FF -1.0_wp

      j = 1
      DO WHILE(ABS(F)>dbl_eps .AND. j<30)
        FF = EXP(wc%XKAPPA / (X0 + wc%ZALP))
        F = wc%ALPHA_CH * X0**2 * FF -1.0_wp
        DF = wc%ALPHA_CH * FF *(2.0_wp * X0 - wc%XKAPPA * (X0 / (X0 + wc%ZALP))**2)
        X0 = X0 - F / DF
        j = j + 1
      END DO
      wc%X0TAUHF = X0

      CONST1 = (wc%BETAMAX/wc%XKAPPA**2)/3.0_wp

      ! Simpson integration weights
      wc%WTAUHF(1) = CONST1
      DO J = 2, wc%jtot_tauhf-1,2
        wc%WTAUHF(J) = 4.0_wp * CONST1
        wc%WTAUHF(J+1) = 2.0_wp * CONST1
      END DO
      wc%WTAUHF(wc%jtot_tauhf) = CONST1

    END DO

  END SUBROUTINE configure_wave


  FUNCTION generate_filename(input_filename, model_base_dir, nroot, jlev, idom, &
    &                        file_datetime, filenumber)  RESULT(result_str)
    CHARACTER(len=*),         INTENT(IN)   :: input_filename, model_base_dir
    INTEGER,                  INTENT(IN)  :: nroot, jlev, idom
    TYPE(datetime), OPTIONAL, INTENT(IN)  :: file_datetime       !< optional argument for <y><m><d><h><min><sec> for backward compatibility
    INTEGER,        OPTIONAL, INTENT(IN)  :: filenumber          !< optional, for <num> wave current_setup

    CHARACTER(LEN=MAX_CHAR_LENGTH) :: str
    CHARACTER(len=MAX_CHAR_LENGTH) :: result_str
    TYPE(t_keyword_list), POINTER  :: keywords => NULL()

    ! Existing ICON placeholders
    CALL associate_keyword("<path>",   TRIM(model_base_dir),             keywords)
    CALL associate_keyword("<nroot>",  TRIM(int2string(nroot,"(i0)")),   keywords)
    CALL associate_keyword("<nroot0>", TRIM(int2string(nroot,"(i2.2)")), keywords)
    CALL associate_keyword("<jlev>",   TRIM(int2string(jlev, "(i2.2)")), keywords)
    CALL associate_keyword("<idom>",   TRIM(int2string(idom, "(i2.2)")), keywords)
    CALL associate_keyword("<dom>",    TRIM(int2string(idom, "(i2.2)")), keywords)

    ! Date/time placeholders from file_datetime
    IF (PRESENT(file_datetime)) THEN
      WRITE(str,'(i4)')   file_datetime%date%year
      CALL associate_keyword("<y>",   TRIM(str), keywords)
      WRITE(str,'(i2.2)') file_datetime%date%month
      CALL associate_keyword("<m>",   TRIM(str), keywords)
      WRITE(str,'(i2.2)') file_datetime%date%day
      CALL associate_keyword("<d>",   TRIM(str), keywords)
      WRITE(str,'(i2.2)') file_datetime%time%hour
      CALL associate_keyword("<h>",   TRIM(str), keywords)
      WRITE(str,'(i2.2)') file_datetime%time%minute
      CALL associate_keyword("<min>", TRIM(str), keywords)
      WRITE(str,'(i2.2)') file_datetime%time%second
      CALL associate_keyword("<sec>", TRIM(str), keywords)
    ENDIF

    ! Optional file number placeholder <num>
    IF (PRESENT(filenumber)) THEN
     CALL associate_keyword("<num>", TRIM(int2string(filenumber, "(i4.4)")), keywords)
    ELSE
     CALL associate_keyword("<num>", "0001", keywords)
    END IF
    ! replace keywords in "input_filename"
    result_str = TRIM(with_keywords(keywords, TRIM(input_filename)))

  END FUNCTION generate_filename


  !< Computes equally spaced directions for wave spectrum discretization.
  !
  SUBROUTINE wave_config__set_directions(me)
    CLASS(t_wave_config) :: me
    ! local
    INTEGER :: jd, ist
    REAL(wp):: startdir
    CHARACTER(*), PARAMETER :: routine = modname//'::set_directions'

    IF (.NOT.ALLOCATED(me%dirs)) THEN
      ALLOCATE(me%dirs(me%ndirs), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for 'dirs' of type REAL failed")
    ENDIF
    !
    IF (.NOT.ALLOCATED(me%sin_dir)) THEN
      ALLOCATE(me%sin_dir(me%ndirs), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for 'sin_dir' of type REAL failed")
    END IF
    !
    IF (.NOT.ALLOCATED(me%cos_dir)) THEN
      ALLOCATE(me%cos_dir(me%ndirs), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for 'cos_dir' of type REAL failed")
    ENDIF

    me%delth = pi2 / REAL(me%ndirs,wp)
    !
    startdir = 0.5_wp * me%delth
    DO jd = 1, me%ndirs
      me%dirs(jd)    = startdir + REAL(jd-1,wp) * me%delth
      me%sin_dir(jd) = SIN(me%dirs(jd))
      me%cos_dir(jd) = COS(me%dirs(jd))
    ENDDO

    IF (msg_level >= 10) THEN
      CALL message ('  ','')
      CALL message ('  ','Directions [Degree]: ')
      DO jd = 1, me%ndirs
        WRITE(message_text,'(i3,f10.5)') jd, me%dirs(jd)*rad2deg
        CALL message ('  ',message_text)
      ENDDO
    ENDIF

  END SUBROUTINE wave_config__set_directions


  !< Computes frequency bins for wave spectrum discretization
  !
  SUBROUTINE wave_config__set_frequencies(me)
    CLASS(t_wave_config) :: me
    ! local
    INTEGER :: jf, ist
    REAL(wp):: co1
    CHARACTER(*), PARAMETER :: routine = modname//'::set_frequencies'

    IF (.NOT.ALLOCATED(me%freqs)) THEN
      ALLOCATE(me%freqs(me%nfreqs), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for 'freqs' of type REAL failed")
    ENDIF
    !
    IF (.NOT.ALLOCATED(me%dfreqs)) THEN
      ALLOCATE(me%dfreqs(me%nfreqs), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for 'dfreqs' of type REAL failed")
    ENDIF
    !
    IF (.NOT.ALLOCATED(me%dfreqs_freqs)) THEN
      ALLOCATE(me%dfreqs_freqs(me%nfreqs), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for 'dfreqs_freqs' of type REAL failed")
    ENDIF
    !
    IF (.NOT.ALLOCATED(me%dfreqs_freqs2)) THEN
      ALLOCATE(me%dfreqs_freqs2(me%nfreqs), stat=ist)
      IF (ist/=SUCCESS) CALL finish(routine, "allocation for 'dfreqs_freqs2' of type REAL failed")
    ENDIF


    ! frequencies
    DO jf = 1,me%nfreqs
      me%freqs(jf) = me%co**(-me%iref + jf)* me%fr1
    ENDDO

    ! intervals
    co1 = 0.5_wp * (me%co - 1.0_wp)
    me%dfreqs(1)         = co1 * me%freqs(1)
    DO jf = 2,me%nfreqs-1
      me%dfreqs(jf) = co1 * (me%freqs(jf) +  me%freqs(jf-1))
    ENDDO
    me%dfreqs(me%nfreqs) = co1 * me%freqs(me%nfreqs-1)
    !
    ! f**2, f**3
    DO jf = 1,me%nfreqs
      me%dfreqs_freqs(jf)  = me%dfreqs(jf) * me%freqs(jf)
      me%dfreqs_freqs2(jf) = me%dfreqs_freqs(jf) * me%freqs(jf)
    ENDDO

    IF (msg_level >= 10) THEN
      CALL message ('  ','')
      CALL message ('  ','Frequencies [Hz], intervals [Hz]: ')
      DO jf = 1, me%nfreqs
        WRITE(message_text,'(i3,2f12.7)') jf, me%freqs(jf), me%dfreqs(jf)
        CALL message ('  ',message_text)
      ENDDO
    ENDIF

  END SUBROUTINE wave_config__set_frequencies

END MODULE mo_wave_config

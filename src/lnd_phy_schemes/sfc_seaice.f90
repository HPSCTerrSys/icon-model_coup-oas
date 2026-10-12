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

! The main program unit of the sea-ice parameterization scheme for NWP.
!
! It contains a number of procedures.
! SUBROUTINE seaice_init_nwp
! initializes the sea-ice scheme and performs some consistency checks.
! SUBROUTINE seaice_timestep_nwp
! advances prognostic variables of the sea-ice scheme one time step.
! SUBROUTINE seaice_coldinit_nwp
! performs a cold start of the sea-ice scheme.
! SUBROUTINE seaice_coldinit_albsi_nwp
! performs a cold start initialization of prognostic sea-ice albedo.
!
! The present sea-ice parameterization scheme is a bulk thermodynamic (no rheology) scheme
! intended for use in NWP and similar applications.
! The scheme is based on a self-similar parametric representation (assumed shape) of the
! evolving temperature profile within the ice and on the integral heat budget of the ice slab.
! The scheme carries ordinary differential equations (in time)
! for the ice surface temperature and the ice thickness.
! An explicit Euler scheme is used for time advance.
! In the current configuration of the scheme, snow over sea ice is not treated explicitly.
! The effect of snow above the ice is accounted for implicitly (parametrically).
! To this end, a rate equation for the sea-ice surface albedo
! with respect to (diffuse) solar radiation is used.
! The rate equation contains
! the relaxation terms that drive sea-ice albedo towards its equilibrium value,
! and "albedo source term" due to precipitation that accounts for the increase
! of albedo after snowfalls.
! The equilibrium albedo is a function of the sea-ice surface temperature.
! Optionally, the sea-ice albedo may be treated diagnostically using
! a temperature-dependent equilibrium albedo.
! For the "sea water" type ICON grid boxes, the snow thickness is set to zero and
! the snow surface temperature is set equal to the ice surface temperature
! (both temperatures are set equal to the fresh-water freezing point if the ice is absent).
!
! Prognostic equations for the ice thickness and the ice surface temperature are solved
! for the ICON grid boxes with the ice fraction
! (area fraction of a given model grid box of the type "sea water" that is covered by ice)
! that exceeds a threshold value of 0.0015. Otherwise, the grid box is treated as ice-free.
! The ice fraction is determined on the basis of observational data
! by the data assimilation scheme and is kept constant over the entire model forecast period.
! However, if the ice melts away during the forecast, the ice fraction is reset to zero.
! This is done within SUBROUTINE update_idx_lists_sea.
! If the ICON grid box is set ice-free during the initialization,
! no ice is created over the forecast period.
! If observational data indicate open water conditions for a given ICON grid box,
! residual ice from the previous model run is removed,
! i.e., the ice thickness is set to zero and
! the ice surface temperature is set to the fresh-water freezing point.
! The newly formed ice has the surface temperature equal to the salt-water freezing point
! and the thickness from 0.1 m to 0.5 m depending on the ice fraction.
! The new ice is formed instantaneously if the data assimilation scheme
! indicates the presence of ice in a given ICON grid box
! but there was no ice in that grid box at the end of the previous model run.
! Prognostic ice thickness is limited by a maximum value of 3 m and a minimum value of 0.05 m.
! Constant values of the ice density, ice molecular heat conductivity, specific heat of ice,
! the latent heat of fusion, and the salt-water freezing point are used.
!
! In uncoupled runs (no interaction of the ice slab with the ocean boundary layer beneath the ice),
! a simple ad hoc formulation of the heat flux from water to ice is used.
! Optionally, adaptive tuning of the parameter(s) of the temperature profile within the ice
! and of the heat flux from water to ice can be used that is based on the assimilation increments
! of the atmopsheric surface layer quantities.
!
! A detailed description of the sea-ice scheme is given in
! Mironov, D., B. Ritter, J.-P. Schulz, M. Buchhold, M. Lange, and E. Machulskaya, 2012:
! Parameterization of sea and lake ice in numerical weather prediction models
! of the German Weather Service.
! Tellus A, 64, 17330. doi:10.3402/tellusa.v64i0.17330
!
! The present sea-ice scheme (with minor modifications)
! is also implemented into the NWP models GME and COSMO
! (see Mironov et al. 2012, for details).
!
! Adaptive parameter tuning in NWP models is described in
! Zaengl, G., 2023:
! Adaptive tuning of uncertain parameters in a numerical weather prediction model
! based upon data assimilation.
! Q. J. R. Meteorol. Soc., 1-20. doi:https://doi.org/10.1002/qj.4535
!
!
!
! Lines embraced with "!_tmp>" and "!_tmp<" contain temporary parts of the code.
! Lines embraced/marked with "!_dev>" and "!_dev<" may be replaced
! as improved formulations are developed and tested.
! Lines embraced/marked with "!_cdm>" and "!_cdm<" are DM's comments that may be helpful to a user.
! Lines embraced/marked with "!_dbg>" and "!_dbg<" are used for debugging purposes only.
! Lines starting with "!_nu" are not used.

MODULE sfc_seaice

  USE mo_kind, ONLY:      &
                   &  wp  !< KIND-type parameter for real variables

  USE mo_exception, ONLY:           &
                        &  finish , &  !< external procedure, finishes model run and reports the reason
                        &  message     !< external procedure, sends a message (error, warning, etc.)

!_cdm>
! Note that ki is equal to 2.1656 in ICON, but is 2.29 in COSMO and GME.
!_cdm<
  USE mo_physical_constants, ONLY:                      &
                                 & tf_fresh => tmelt  , &  !< fresh-water freezing point [K]
                                 &             alf    , &  !< latent heat of fusion [J/kg]
                                 &             rhoi   , &  !< density of ice [kg/m^3]
                                 & rhos_def => rhos   , &  !< default snow density [kg/m^3]
                                 &             ci     , &  !< specific heat of ice [J/(kg K)]
                                 &             cs     , &  !< specific heat of snow [J/(kg K)]
                                 &             ki          !< molecular heat conductivity of ice [J/(m s K)]

  USE mo_lnd_nwp_config,     ONLY:                      &
                                 & lprog_albsi        , &  !< sea-ice albedo is computed prognostically
                                 & lbottom_hflux      , &  !< use parameterization for heat flux through sea ice bottom
                                 & frsi_min           , &  !< minimum sea-ice fraction  [-]
                                 & hice_min           , &  !< minimum sea-ice thickness [m]
                                 & hice_max           , &  !< maximum sea-ice thickness [m]
                                 & albsi_snow_max     , &  !< maximum albedo of snow over sea ice [-]
                                 & albsi_snow_min     , &  !< minimum albedo of snow over sea ice [-]
                                 & albsi_max          , &  !< maximum albedo of sea ice [-]
                                 & albsi_min          , &  !< minimum albedo of sea ice [-]
                                 & lsnow_on_seaice    , &  !< consider snow on seaice
                                 & tf_salt                 !< salt-water freezing point [K]

  USE mo_coupling_config,    ONLY: is_coupled_to_ocean     !< TRUE for coupled ocean-atmosphere runs

  USE mo_lnd_nwp_config,     ONLY: lcuda_graph_lnd

  USE mo_fortran_tools,      ONLY: set_acc_host_or_device  !< Routine can be run on CPU and on GPU

  IMPLICIT NONE

  PRIVATE

  REAL (wp), PARAMETER ::                             &
                       &  hice_ini_min = 0.1_wp     , &  !< minimum thickness of the newly formed sea ice [m]
                       &  hice_ini_max = 0.5_wp     , &  !< maximum thickness of the newly formed sea ice [m]
                       &  csi_lin      = 0.5_wp     , &  !< shape factor for linear temperature profile [-]
                       &  phiipr0_lin  = 1.0_wp     , &  !< derivative (at zeta=0) of the linear
                                                         !< temperature profile shape function [-]
                       &  csidp_nlin   = 2.0_wp     , &  !< disposable parameter for non-linear
                                                         !< temperature profile shape factor [-]
                       &  csidp_nlin_d =              &  !< derived disposable parameter for non-linear
                                                         !< temperature profile shape factor [-]
                       &  (1._wp+csidp_nlin)/12._wp , &
                       &  cmaxearg     = 1.0E+02_wp , &  !< maximum value of the EXP function argument (security constant) [-]
                       &  csmall       = 1.0E-05_wp      !< small number (security constant) [-]

  REAL (wp), PARAMETER ::                             &
    &  taualbsi_min   = 7.0_wp*86400._wp            , &  !< minimum relaxation time scale for sea-ice albedo [s]
    &  taualbsi_max   = 21.0_wp*86400._wp           , &  !< maximum relaxation time scale for sea-ice albedo [s]
    &  t_taualbsi_min = 268.15_wp                   , &  !< lower bound of the temperature range in the interpolation
                                                         !< formula for the sea-ice albedo relaxation time scale [K]
    &  rdelt_taualbsi =                               &  !< reciprocal of the temperature range in the interpolation
    &    1._wp/(tf_fresh-t_taualbsi_min)            , &  !< formula for the sea-ice albedo relaxation time scale [K^{-1}]
    &  c2_albsi_snow  = 136.6_wp/tf_fresh           , &  !< constant in the expression for albedo
                                                         !< of snow over sea ice [K^{-1}]
    &  c_tausi_snow   = 1._wp/5._wp                 , &  !< constant used to define the relaxation time scale towards
                                                         !< the equilibrium albedo of snow over sea ice [(kg/m^2)^{-1}]
                                                         !< (corresponds to 5 mm snow water equaivalent precipitated over
                                                         !< an e-folding time scale)
    &  t_albsi_snow_max  = 272.95_wp                     !< upper bound of the temperature range over which relaxation
                                                         !< towards snow-overice albedo is applied [K]

  REAL (wp), PARAMETER :: csidp_nlin_sonsi = 2.2_wp !< Phi_*i adjusted to MOSAiC data (hice < 1.4m) for hice_max = 6m [-].
  REAL (wp), PARAMETER :: csidp_nlin_sonsi_d = (1._wp+csidp_nlin_sonsi)/12._wp !< Height-dependent part of the ice shape factor [-].
  REAL (wp), PARAMETER :: cssdp_nlin = 1.3_wp !< Phi_*s from MOSAiC data [-].
  REAL (wp), PARAMETER :: hsnow_max = 0.6_wp !< Limit for snow temperature profile from MOSAiC data [m].
  REAL (wp), PARAMETER :: cssdp_nlin_d = (1._wp+cssdp_nlin)/12._wp !< Height-dependent part of the snow shape factor [-].
  REAL (wp), PARAMETER :: snow_frac_scale = 0.05_wp !< Scaling height for full snow cover (for albedo) [m].
  REAL (wp), PARAMETER :: r_alb_snow_si_scale = 0.5_wp !< Reciprocal scale temperature for snow albedo [K**-1].
  REAL (wp), PARAMETER :: r_alb_seaice_scale = 0.35_wp !< Reciprocal scale temperature for ice albedo [K**-1].

  !>
  !! Optical characteristics of sea ice.
  !!
  !! A storage for an n-band approximation of the exponential decay law
  !! for the flux of solar radiation is allocated.
  !! A maximum value of the wave-length bands is currently set equal to two.
  !!

  INTEGER, PARAMETER ::                      &
                     &  nband_optic_max = 2  !< maximum number of wave-length bands in the decay law
                                             !< for the solar radiation flux [-]

  TYPE opticpar_seaice
    INTEGER ::                                       &
            &  nband_optic                              !< number of wave-length bands [-]
    REAL (wp) ::                                     &
              &  frac_optic(nband_optic_max)       , &  !< fractions of total solar radiation flux for different bands [-]
              &  extincoef_optic(nband_optic_max)       !< extinction coefficients for different bands [1/m]
  END TYPE opticpar_seaice

  ! One-band approximation for opaque sea ice.
  ! The use of large extinction coefficient prevents the penetration of solar radiation
  ! into the ice interior, i.e. the volumetric character of the solar radiation heating is ignored.
  TYPE (opticpar_seaice), PARAMETER ::                             &
                                    &  opticpar_seaice_opaque =    &
                                    &  opticpar_seaice(1,          &
                                    &  (/1._wp, 0._wp/),           &
                                    &  (/1.0E+07_wp, 1.E+10_wp/))

  ! Minimum values of the sea-ice fraction and of the sea-ice thickness are used outside
  ! "sfc_seaice" to compose an index list of grid boxes where sea ice is present.
  PUBLIC ::                       &
         &  hice_ini_min             , & ! parameter
         &  hice_ini_max             , & ! parameter
         &  seaice_init_nwp          , & ! procedure
         &  seaice_coldinit_nwp      , & ! procedure
         &  seaice_coldinit_albsi_nwp, & ! procedure
         &  seaice_timestep_nwp      , & ! procedure
         &  alb_seaice_equil

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

CONTAINS

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Prognostic variables of the sea-ice parameterization scheme are initialized
  !! and some consistency checks are performed.
  !!
  !! The procedure arguments are arrays (vectors) of the sea-ice fraction, and
  !! of the sea-ice scheme prognostic variables.
  !! The vector length is equal to the number of grid boxes (within a given block)
  !! where sea ice is present, i.e. where the sea-ice fraction exceeds its minimum value.
  !! First, the sea-ice fraction is checked.
  !! If a value less than a minimum threshold value is found
  !! (indicating that a grid box is declared as partially ice-covered but no sea ice
  !! should be present), an error message is sent and the model abort is called.
  !! Next, "new" ice is formed in the grid boxes where "old" value of the sea-ice thickness
  !! is less than a minimum threshold value (i.e. there was no ice at the end
  !! of the the previous model run). The newly formed ice has the surface temperature equal
  !! to the salt-water freezing point and the thickness of 0.1 m to 0.5 m
  !! depending on the ice fraction.
  !! For the grid boxes where new ice is created, prognostic sea-ice albedo is set equal
  !! to its equilibrium value (function of sea-ice surface temperature).
  !! Then, the ice thickness is limited from above and below
  !! by the maximum and minimum threshold values,
  !! and the ice surface temperature is limited from above by the fresh-water freezing point.
  !! These security measures (taken for each grid box where sea ice is present)
  !! are required to avoid non-allowable values of prognostic variables
  !! that may occure due to a loss of accuracy during the model IO
  !! (e.g. due to GRIB encoding and decoding).
  !! Finally, the snow thickness is set to zero, and
  !! the snow surface temperature is set equal to the ice surface temperature
  !! (recall that snow over sea ice is not treated explicitly).
  !!

  SUBROUTINE seaice_init_nwp (                                  &
                          &  linit_hice,                        &
                          &  nsigb,                             &
                          &  frsi,                              &
                          &  tice_p, hice_p, tsnow_p, hsnow_p,  &
                          &  albsi_p,                           &
                          &  tice_n, hice_n, tsnow_n, hsnow_n,  &
                          &  albsi_n,                           &
                          &  lacc)

    IMPLICIT NONE

    ! Procedure arguments

    LOGICAL, INTENT(IN) ::        &
                        &  linit_hice     !< TRUE if initialization of hice is required
                                          !< for new seaice points.
                                          !< FALSE if updated information on hice is provided
                                          !< (e.g. in coupled atmosphere-ocean runs)

    LOGICAL, OPTIONAL, INTENT(IN) :: lacc !< if true, use OpenACC

    INTEGER, INTENT(IN) ::        &
                        &  nsigb  !< Array (vector) dimension
                                  !< (equal to the number of grid boxes within a block
                                  !< where the sea ice is present)

    REAL(wp), DIMENSION(:), INTENT(IN)    ::         &
                                          &  frsi       !< sea-ice fraction [-]

    REAL(wp), DIMENSION(:), INTENT(INOUT) ::           &
                                          &  tice_p  , &  !< temperature of ice upper surface at previous time level [K]
                                          &  hice_p  , &  !< ice thickness at previous time level [m]
                                          &  tsnow_p , &  !< temperature of snow upper surface at previous time level [K]
                                          &  hsnow_p , &  !< snow thickness at previous time level [m]
                                          &  albsi_p , &  !< sea-ice albedo at previous time level [-]
                                          &  tice_n  , &  !< temperature of ice upper surface at new time level [K]
                                          &  hice_n  , &  !< ice thickness at new time level [m]
                                          &  tsnow_n , &  !< temperature of snow upper surface at new time level [K]
                                          &  hsnow_n , &  !< snow thickness at new time level [m]
                                          &  albsi_n      !< sea-ice albedo at new time level [-]

    ! Local variables

    LOGICAL :: lzacc !< non-optional version of lacc

    CHARACTER(len=256) ::           &
                       &  nameerr , &  !< name of procedure where an error occurs
                       &  texterr      !< error/warning message text

    INTEGER ::      &
            &  isi  !< DO loop index

    REAL(wp) :: frsi_err !< used for check if lowest sea-ice fraction is lower than allowed value frsi_min


    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    CALL set_acc_host_or_device(lzacc, lacc)

    frsi_err = frsi_min
    !$ACC PARALLEL ASYNC(1) DEFAULT(PRESENT) REDUCTION(MIN: frsi_err) IF(lzacc)
    !$ACC LOOP GANG VECTOR REDUCTION(MIN: frsi_err)
    DO isi=1, nsigb
      ! Find minimum sea-ice fraction
      IF (frsi(isi) < frsi_err) THEN
        frsi_err = frsi(isi)
      END IF
    END DO
    !$ACC END PARALLEL
    !$ACC WAIT

    ! Call model abort if errors are encountered
    IF (frsi_err < frsi_min) THEN
      ! Send an error message
      WRITE(nameerr,*) "MODULE sfc_seaice, SUBROUTINE seaice_init_nwp"
      WRITE(texterr,*) "Sea-ice fraction ", frsi_err,                        &
                    &  " is less than a minimum threshold value ", frsi_min, &
                    &  " Call model abort."
      CALL message(TRIM(nameerr), TRIM(texterr), all_print=.TRUE.)

      ! Call model abort
      WRITE(nameerr,*) "sfc_seaice:seaice_init_nwp"
      WRITE(texterr,*) "error in sea-ice fraction"
      CALL finish(TRIM(nameerr), TRIM(texterr))
    END IF

    ! Loop over grid boxes where sea ice is present
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    !$ACC LOOP GANG VECTOR
    GridBoxesWithSeaIce: DO isi=1, nsigb

      IF ( linit_hice ) THEN
        ! Here we assume that new seaice points are characterized by (fr_seaice>0, hice_p<hice_min)
        ! Create new ice
        IF( hice_p(isi) < (hice_min-csmall) ) THEN
          hice_p(isi) = hice_ini_min + frsi(isi) * (hice_ini_max-hice_ini_min)
          hsnow_p(isi) = 0._wp
          tice_p(isi) = tf_salt
          tsnow_p(isi) = tice_p(isi)
          ! Set sea-ice albedo to its equilibrium value
          ! (only required if sea-ice albedo is treated prognostically)
          IF ( lprog_albsi ) THEN
            albsi_p(isi) = alb_seaice_equil( tsnow_p(isi), 0._wp )
          ENDIF
        ENDIF

      ELSE ! linit_hice=.FALSE.
        !
        ! Here we assume that an educated guess for ice thickness is provided,
        ! with hice>=hice_min (e.g. in coupled mode) in addition to the ice fraction.
        ! Hence, we skip the initialization of ice thickness but keep
        ! the initialization of ice temperature and albedo for new ice points.
        !
        ! In this case the identification of new ice points by the above condition
        ! (fr_seaice>0, hice_p<hice_min) fails and so does any identification via tice.
        ! Hence, new ice points must be identified before entering
        ! this routine by comparing the old and new seaice fraction fields and
        ! by creating a "NEW seaice points only" index list.
        !
        ! For this IF condition to work correctly, seaice_init_nwp must be called for
        ! NEW seaice points only rather than for ALL seaice points
        ! (see example in mo_nwp_sfc_utils:process_sst_and_seaice).
        ! Otherwise ice temperature and albedo will erroneously be reset
        ! to default values for all seaice points.
        !
        tice_p(isi) = tf_salt
        tsnow_p(isi) = tice_p(isi)
        ! Set sea-ice albedo to its equilibrium value
        ! (only required if sea-ice albedo is treated prognostically)
        IF ( lprog_albsi ) THEN
          albsi_p(isi) = alb_seaice_equil( tsnow_p(isi), hsnow_p(isi) )
        ENDIF
      ENDIF


      ! In general we assume that new seaice points are characterized by
      ! ( fr_seaice>0, h_ice_p=0 ). However, it may happen that h_ice_p
      ! has already been adjusted consistently by the data assimilation process.
      ! In that case, we fail to identify new sea-ice points by the above condition,
      ! and we miss the initialization of the prognostic seaice albedo. Thus,
      ! the following statement is added.
      IF ( lprog_albsi .AND. albsi_p(isi) <= 0._wp) THEN
        albsi_p(isi) = alb_seaice_equil( tsnow_p(isi), hsnow_p(isi) )
      ENDIF

      ! Take security measures
      hice_p(isi) = MAX(MIN(hice_p(isi), hice_max), hice_min)
      tice_p(isi) = MIN(tice_p(isi), tf_fresh)
      tsnow_p(isi) = MIN(tsnow_p(isi), tf_fresh)

      ! Set variables at new time level
      tice_n(isi)  = tice_p(isi)
      hice_n(isi)  = hice_p(isi)
      tsnow_n(isi) = tsnow_p(isi)
      hsnow_n(isi) = hsnow_p(isi)
      IF ( lprog_albsi ) THEN
        albsi_n(isi) = albsi_p(isi)
      ENDIF

    END DO GridBoxesWithSeaIce
    !$ACC END PARALLEL

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_init_nwp


!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Prognostic variables of the sea-ice scheme are advanced one time step.
  !!
  !! Ordinary differential equations (in time) for the ice surface temperature and
  !! the ice thickness are solved using an explicit Euler scheme for time advance.
  !! The sea-ice surface albedo with respect to solar radiation is determined
  !! by solving a rate equation (if the sea-ice albedo is treated diagnostically,
  !! no albedo calculations are performed in the present routine).
  !! The shape factor for the temperature profile within the ice and the derivative of the
  !! temperature profile shape function at the underside of the ice are functions of the ice
  !! thickness (see Mironov et al. 2012, for details).
  !! Optionally, constant values of the shape factor and of the shape-function derivative
  !! corresponding to the linear temperature profile within the ice
  !! (cf. the sea-ice schemes of GME and COSMO) can be used
  !! (as there is no logical switch to activate this option, changes in the code should be made).
  !! In the regime of ice growth or melting from below, the solution may become spurious
  !! when the ice thickness is small and/or the model time step is large.
  !! In such a case, a quasi-steady heat transfer through the ice is assumed
  !! and the differential equation for the ice surface temperature
  !! is reduced to an algebraic relation.
  !! In the current configuration of the sea-ice scheme, snow over ice is not treated explicitly.
  !! The effect of snow is accounted for implicitely through changes
  !! of the sea-ice albedo with respect to solar radiation.
  !! For the "sea water" type grid boxes, the snow thickness is set to zero and
  !! the snow surface temperature is set equal to the ice surface temperature.
  !! Prognostic ice thickness is limited by a maximum value of 3 m and a minimum value of 0.05 m.
  !! In "uncoupled" runs, i.e., where the interaction between the sea ice and the sea water
  !! beneath is not considered, an ad hoc formulation for the heat flux from water to ice is used
  !! that serves to quench the ice growth rate as the ice thickness approaches its maximum value.
  !! No ice is created during the forecast period.
  !! If the ice melts away during the forecast (i.e., the ice becomes thinner than 0.05 m),
  !! the ice tickness is set to zero and
  !! the ice surface temperatures is set to the fresh-water freezing point.
  !! The procedure arguments are arrays (vectors)
  !! of the sea-ice scheme prognostic variables,
  !! of the components of the heat balance at the ice upper surface
  !! (i.e., the fluxes of sensible and latent heat, the net flux of long-wave radiation,
  !! and the net flux of solar radiation with due regard for the ice surface albedo),
  !  of the precipitation rates of snow and rain,
  !! and of the sea-ice surface albedo with respect to solar radiation.
  !! Fluxes are positive when directed downward.
  !! The vector length is equal to the number of model grid boxes
  !! (within a given block) where sea ice is present
  !! (i.e., where the sea-ice fraction exceeds its minimum threshold value).
  !! The time tendecies of the ice thickness, the ice surface temperature,
  !! the snow thickness and the snow surface temperature are also computed.
  !! These are optional arguments of the procedure.
  !!

  SUBROUTINE seaice_timestep_nwp (                                      &
                              &  dtime,                                 &
                              &  nsigb,                                 &
                              &  qsen, qlat, qlwrnet, qsolnet,          &
                              &  snow_rate, rain_rate, fac_bottom_hflx, &
                              &  tice_p, hice_p, tsnow_p, hsnow_p,      &
                              &  albsi_p,                               &
                              &  tice_n, hice_n, tsnow_n, hsnow_n,      &
                              &  condhf, meltpot, albsi_n,              &
                              &  opt_dticedt, opt_dhicedt, opt_dtsnowdt,&
                              &  opt_dhsnowdt                           )

    IMPLICIT NONE

    ! Procedure arguments

    REAL(wp), INTENT(IN) ::        &
                         &  dtime  !< model time step [s]

    INTEGER, INTENT(IN) ::        &
                        &  nsigb  !< number of grid boxes within a block
                                  !< where the sea ice is present (<=nproma)

    REAL(wp), DIMENSION(:), INTENT(IN) ::           &
                                       &  qsen    , &  !< sensible heat flux at the surface [W/m^2]
                                       &  qlat    , &  !< latent heat flux at the surface [W/m^2]
                                       &  qlwrnet , &  !< net long-wave radiation flux at the surface [W/m^2]
                                       &  qsolnet      !< net solar radiation flux at the surface [W/m^2]

    REAL(wp), DIMENSION(:), INTENT(IN) ::             &
                                       &  snow_rate , &  !< snow rate (convecive + grid-scale) [kg/(m^2 s)]
                                       &  rain_rate      !< rain rate (convecive + grid-scale) [kg/(m^2 s)]

    REAL(wp), DIMENSION(:), OPTIONAL, INTENT(IN) ::   &
                                       &  fac_bottom_hflx  !< tuning factor for heat flux from water to ice [-]

    REAL(wp), DIMENSION(:), INTENT(IN) ::           &
                                       &  tice_p  , &  !< temperature of ice upper surface at previous time level [K]
                                       &  hice_p  , &  !< ice thickness at previous time level [m]
                                       &  tsnow_p , &  !< temperature of snow upper surface at previous time level [K]
                                       &  hsnow_p , &  !< snow thickness at previous time level [m]
                                       &  albsi_p      !< sea-ice albedo at previous time level [-]

    REAL(wp), DIMENSION(:), INTENT(OUT) ::          &
                                       &  tice_n  , &  !< temperature of ice upper surface at new time level [K]
                                       &  hice_n  , &  !< ice thickness at new time level [m]
                                       &  tsnow_n , &  !< temperature of snow upper surface at new time level [K]
                                       &  hsnow_n , &  !< snow thickness at new time level [m]
                                       &  condhf  , &  !< conductive heat flux within the sea ice
                                                       !< just above the ice lower boundary [W/m^2]
                                       &  meltpot , &  !< melt potential at top [W/m^2]
                                       &  albsi_n      !< sea-ice albedo at new time level [-]

    REAL(wp), DIMENSION(:), INTENT(OUT), OPTIONAL ::    &
                                       &  opt_dticedt , &  !< time tendency of ice surface temperature [K/s]
                                       &  opt_dhicedt , &  !< time tendency of ice thickness [m/s]
                                       &  opt_dtsnowdt, &  !< time tendency of snow surface temperature [K/s]
                                       &  opt_dhsnowdt     !< time tendency of snow thickness [m/s]

    INTEGER :: isi
    LOGICAL :: lhave_dtsnowdt
    LOGICAL :: lhave_dhsnowdt

    IF (lsnow_on_seaice) THEN
      CALL seaice_timestep_snow_nwp ( &
          & dtime=dtime, &
          & nsigb=nsigb, &
          & qsen=qsen, &
          & qlat=qlat, &
          & qlwrnet=qlwrnet, &
          & qsolnet=qsolnet, &
          & snow_rate=snow_rate, &
          & rain_rate=rain_rate, &
          & fac_bottom_hflx=fac_bottom_hflx, &
          & tice_p=tice_p, &
          & hice_p=hice_p, &
          & tsnow_p=tsnow_p, &
          & hsnow_p=hsnow_p, &
          & albsi_p=albsi_p, &
          & tice_n=tice_n, &
          & hice_n=hice_n, &
          & tsnow_n=tsnow_n, &
          & hsnow_n=hsnow_n, &
          & albsi_n=albsi_n, &
          & condhf=condhf, &
          & meltpot=meltpot, &
          & opt_dticedt=opt_dticedt, &
          & opt_dhicedt=opt_dhicedt, &
          & opt_dtsnowdt=opt_dtsnowdt, &
          & opt_dhsnowdt=opt_dhsnowdt &
        )
    ELSE
      CALL seaice_timestep_nosnow_nwp ( &
          & dtime=dtime, &
          & nsigb=nsigb, &
          & qsen=qsen, &
          & qlat=qlat, &
          & qlwrnet=qlwrnet, &
          & qsolnet=qsolnet, &
          & snow_rate=snow_rate, &
          & rain_rate=rain_rate, &
          & fac_bottom_hflx=fac_bottom_hflx, &
          & tice_p=tice_p, &
          & hice_p=hice_p, &
          & albsi_p=albsi_p, &
          & tice_n=tice_n, &
          & hice_n=hice_n, &
          & albsi_n=albsi_n, &
          & condhf=condhf, &
          & meltpot=meltpot, &
          & opt_dticedt=opt_dticedt, &
          & opt_dhicedt=opt_dhicedt &
        )

      lhave_dhsnowdt = PRESENT(opt_dhsnowdt)
      lhave_dtsnowdt = PRESENT(opt_dtsnowdt)

      !$ACC DATA ASYNC(1) NO_CREATE(opt_dhsnowdt, opt_dtsnowdt) PRESENT(nsigb)
      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
      DO isi = 1, nsigb
        tsnow_n(isi) = tice_n(isi)
        hsnow_n(isi) = 0._wp
        IF (lhave_dhsnowdt) opt_dhsnowdt(isi) = 0._wp
        IF (lhave_dtsnowdt) opt_dtsnowdt(isi) = (tice_n(isi) - tice_p(isi)) / dtime
      END DO
      !$ACC END DATA
    END IF

  END SUBROUTINE seaice_timestep_nwp

  SUBROUTINE seaice_timestep_nosnow_nwp (                               &
                              &  dtime,                                 &
                              &  nsigb,                                 &
                              &  qsen, qlat, qlwrnet, qsolnet,          &
                              &  snow_rate, rain_rate, fac_bottom_hflx, &
                              &  tice_p, hice_p, albsi_p,               &
                              &  tice_n, hice_n, condhf, meltpot,       &
                              &  albsi_n,                               &
                              &  opt_dticedt, opt_dhicedt               )

    IMPLICIT NONE

    ! Procedure arguments

    REAL(wp), INTENT(IN) ::        &
                         &  dtime  !< model time step [s]

    INTEGER, INTENT(IN) ::        &
                        &  nsigb  !< number of grid boxes within a block
                                  !< where the sea ice is present (<=nproma)

    REAL(wp), DIMENSION(:), INTENT(IN) ::           &
                                       &  qsen    , &  !< sensible heat flux at the surface [W/m^2]
                                       &  qlat    , &  !< latent heat flux at the surface [W/m^2]
                                       &  qlwrnet , &  !< net long-wave radiation flux at the surface [W/m^2]
                                       &  qsolnet      !< net solar radiation flux at the surface [W/m^2]

    REAL(wp), DIMENSION(:), INTENT(IN) ::             &
                                       &  snow_rate , &  !< snow rate (convecive + grid-scale) [kg/(m^2 s)]
                                       &  rain_rate      !< rain rate (convecive + grid-scale) [kg/(m^2 s)]

    REAL(wp), DIMENSION(:), OPTIONAL, INTENT(IN) ::   &
                                       &  fac_bottom_hflx  !< tuning factor for heat flux from water to ice [-]

    REAL(wp), DIMENSION(:), INTENT(IN) ::           &
                                       &  tice_p  , &  !< temperature of ice upper surface at previous time level [K]
                                       &  hice_p  , &  !< ice thickness at previous time level [m]
                                       &  albsi_p      !< sea-ice albedo at previous time level [-]

    REAL(wp), DIMENSION(:), INTENT(OUT) ::          &
                                       &  tice_n  , &  !< temperature of ice upper surface at new time level [K]
                                       &  hice_n  , &  !< ice thickness at new time level [m]
                                       &  condhf  , &  !< conductive heat flux within the sea ice
                                                       !< just above the ice lower boundary [W/m^2]
                                       &  meltpot , &  !< melt potential at top [W/m^2]
                                       &  albsi_n      !< sea-ice albedo at new time level [-]

    REAL(wp), DIMENSION(:), INTENT(OUT), OPTIONAL ::    &
                                       &  opt_dticedt , &  !< time tendency of ice surface temperature [K/s]
                                       &  opt_dhicedt      !< time tendency of ice thickness [m/s]

    ! Derived parameters
    ! (combinations of physical constants encountered several times in the code)

    REAL (wp), PARAMETER ::                                   &
                         &  r_rhoici     = 1._wp/(rhoi*ci)  , &
                         &  ki_o_rhoici  = ki*r_rhoici      , &
                         &  r_rhoialf    = 1._wp/(rhoi*alf) , &
                         &  ki_o_rhoialf = ki*r_rhoialf     , &
                         &  ci_o_alf     = ci/alf

    ! Local variables
    ! Need to allocate these variables to nproma size for the CUDA graphs
    !  because nsigb may change from step to step and we need to make sure
    !  these arrays are large enough for any nsigb at the captured time step
    REAL(wp), DIMENSION(SIZE(qsen)) ::       &
                                &  dticedt , &  !< time tendency of ice surface temperature [K/s]
                                &  dhicedt , &  !< time tendency of ice thickness [m/s]
                                &  thetaipr0    !< Deriavtive (at z=0) of ice temperature [K/m]

    INTEGER ::      &
            &  isi  !< DO loop index

    REAL (wp) ::                &
              &  qsoliw       , &  !< solar radiation flux at the ice-water interface (positive downward) [W/m^2]
              &  qatm         , &  !< "total atmospheric heat flux" for the ice slab  (positive downward) [W/m^2]
                                   !< (sum of the sensible heat flux, latent heat flux and the net flux
                                   !< of long-waver radiation at the ice upper surface,
                                   !< and the difference of solar radiation fluxes
                                   !< at the upper and lower surfaces of the ice slab)
              &  qwat         , &  !< heat flux from water to ice (positive downward) [W/m^2]
              &  csice        , &  !< shape factor for the temperature profile within the ice [-]
              &  phiipr0      , &  !< derivative (at zeta=0) of the temperature profile shape function [-]
              &  hice_thrshld , &  !< threshold value of the ice tickness to switch between a quasi-equilibrium
                                   !< model of heat transfer through the ice and a complete model [m]
              &  rti          , &  !< dimensionless parameter [-]
              &  r_dtime      , &  !< reciprocal of the time step [1/s]
              &  strg_1       , &  !< temporary storage variable
              &  strg_2            !< temporary storage variable

    REAL (wp) ::                &
              &  albsi_e      , &  !< equilibrium sea-ice albedo [-]
              &  albsi_snow_e , &  !< equilibrium albedo of snow over sea ice [-]
              &  c1_albsi_snow, &  !< for use in the expression for albedo of snow over sea ice
              &  taualbsi     , &  !< relaxation time scale for sea-ice albedo [s]
              &  rtaualbsisn  , &  !< reciprocal of relaxation time scale for snow-over-ice albedo [s^{-1}]
              &  albsi_e_wghtd     !< weighted equilibrium albedo (storage variable) [-]

    LOGICAL ::   lis_coupled_to_ocean !< TRUE for coupled ocean-atmosphere runs (copy for vectorisation)
    LOGICAL ::   lpres_fac_hflux   !< TRUE if fac_bottom_hflx is provided as input (from NWP interface only)

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Reciprocal of the time step
    r_dtime = 1._wp/dtime

    ! For use in the expression for albedo of snow over sea ice
    c1_albsi_snow  = 1._wp-albsi_snow_min/albsi_snow_max

    ! If fac_bottom_hflx is provided, adaptive tuning of the parameter(s)
    ! of the temperature profile within the ice and of the heat flux from water to ice is used
    lpres_fac_hflux = PRESENT(fac_bottom_hflx)

    ! for vectorisation
    lis_coupled_to_ocean = is_coupled_to_ocean()

    !$ACC DATA CREATE(dticedt, dhicedt, thetaipr0) &
    !$ACC   PRESENT(nsigb) &
    !$ACC   NO_CREATE(fac_bottom_hflx) ASYNC(1)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR &
    !$ACC   PRIVATE(qsoliw, qatm, qwat, csice, phiipr0, hice_thrshld, rti, strg_1, strg_2) &
    !$ACC   PRIVATE(albsi_e, albsi_snow_e, taualbsi, rtaualbsisn, albsi_e_wghtd)
    ! Loop over grid boxes where sea ice is present
    GridBoxesWithSeaIce: DO isi=1, nsigb

      ! Compute solar radiation flux at the ice-water interface (positive downward)
      qsoliw = qsolnet(isi)*(                                                                 &
             & opticpar_seaice_opaque%frac_optic(1)                                           &
             & *EXP(-MIN(opticpar_seaice_opaque%extincoef_optic(1)*hice_p(isi), cmaxearg)) +  &
             & opticpar_seaice_opaque%frac_optic(2)                                           &
             & *EXP(-MIN(opticpar_seaice_opaque%extincoef_optic(2)*hice_p(isi), cmaxearg)) )

      ! Compute total atmospheric heat flux for the ice slab  (positive downward)
      qatm = qsen(isi) + qlat(isi) + qlwrnet(isi) + qsolnet(isi) - qsoliw

      ! Provision is made to account for the heat flux from water to ice (upward flux is negative).
      ! This is the ocean heat flux just below the ice,
      ! not including the latent flux from melting and freezing.
      ! In the case of coupled icon atmosphere-ocean runs,
      ! this flux has to be passed from the ocean if available.
      ! For uncoupled runs, an ad hoc formulation for the heat flux from water to ice is used
      ! that basically serves to quench the ice growth rate
      ! as the ice thickness approaches its maximum value.

      IF (lbottom_hflux) THEN
        ! Derivative (at zeta=0) of the temperature profile shape function.
        IF (lpres_fac_hflux) THEN
          ! Adaptive tuning of phiipr0 and of bottom heat flux
          phiipr0 = phiipr0_lin - fac_bottom_hflx(isi)*hice_p(isi)/hice_max
        ELSE
          ! A constant value is used (cf. linear temperature profile).
          phiipr0 = phiipr0_lin
        END IF
        ! Heat flux from water to ice is limited from above (no flux from ice to water is allowed).
        qwat = (1._wp-hice_p(isi)/hice_max-phiipr0)*ki*(tf_salt-tice_p(isi))/MAX(hice_p(isi),hice_min)
        qwat = MIN(qwat, 0._wp)
      ELSE
        ! Derivative (at zeta=0) of the temperature profile shape function
        phiipr0 = 1._wp - hice_p(isi)/hice_max
        ! Heat flux from water to ice is set to zero.
        qwat = 0._wp
      END IF

      ! Compute temperature profile shape factor and temporary storage variables
      ! (to recover linear temperature profile, set csice=csi_lin)
      csice = csi_lin - csidp_nlin_d*hice_p(isi)/hice_max
      rti = ci_o_alf*(tice_p(isi)-tf_salt)
      strg_1 = rti*(1.5_wp-2.0_wp*csice)
      strg_2 = 1._wp + strg_1

      FreezingMeltingRegime: IF( tice_p(isi)>=(tf_fresh-csmall) .AND. qatm>0._wp ) THEN

        ! Melting from above

        ! Set the ice surface temperature equal to the fresh-water freezing point
        tice_n(isi) = tf_fresh

        IF ( lis_coupled_to_ocean ) THEN
          ! Coupling to icon-o:
          ! Heat flux within the ice just above the ice-water intrface, phiipr0 is computed above
          ! Note that for coupled runs, lbottom_hflux=.FALSE. and the heat flux
          ! from water to ice qwat is set to zero (qwat should be provided by the ocean model).
          condhf(isi) = ki*phiipr0*(tice_n(isi)-tf_salt)/hice_p(isi)
          meltpot(isi) = (qatm-qwat)/strg_2
          ! No change of ice thickness
          hice_n(isi) = hice_p(isi)
        ELSE
          ! Compute the rate of ice melting (note the sign of heat fluxes)
          dhicedt(isi) = -(qatm-qwat)*r_rhoialf/strg_2
          ! Update the ice thickness
          hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)
        END IF

      ELSE FreezingMeltingRegime

        ! Freezing or melting from below

        ! Compute threshold value of the ice thickness
        ! Note that an expression in parentheses should be multiplied with
        ! MAX(1._wp, ABS(2._wp*csice*rti)) (cf. the code of the lake parameterization scheme FLake).
        ! However, |2*csice*(ci/alf)*(tice-tf_salt)| < 1
        ! at all conceivable values of the ice surface temperature.
        hice_thrshld = SQRT(phiipr0*ki_o_rhoici*dtime/csice)

        IF( hice_p(isi)<hice_thrshld) THEN

          ! Use a quasi-equilibrium model of heat transfer through the ice

          IF ( lis_coupled_to_ocean ) THEN
            ! Coupling to icon-o:
            ! Heat flux (phiipr0 is computed above)
            condhf(isi) = ki*phiipr0*(tice_p(isi)-tf_salt)/hice_p(isi)
            meltpot(isi) = 0._wp
            ! No change of ice thickness
            hice_n(isi) = hice_p(isi)
          ELSE

            ! Compute the time-rate-of-change of the ice thickness (note the sign of heat fluxes)
            dhicedt(isi) = -(qatm-qwat)*r_rhoialf/strg_2
            ! Update the ice thickness
            hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)

          END IF

          ! Compute the (updated) ice surface temperature (note the sign of heat fluxes)
          tice_n(isi) = tf_salt + (qatm+qwat*strg_1)*hice_n(isi)/(phiipr0*ki*strg_2)

        ELSE

          ! Use a complete model of heat transfer through the ice

          thetaipr0(isi) = phiipr0*(tice_p(isi)-tf_salt)/hice_p(isi)

          ! Compute the time-rate-of-change of the ice surface temperature
          ! (note the sign of heat fluxes)
          dticedt(isi) = ( (qatm+strg_1*qwat)*r_rhoici - ki_o_rhoici*thetaipr0(isi)*strg_2 )  &
                     & /(csice*hice_p(isi))
          ! Update the ice surface temperature
          tice_n(isi) = tice_p(isi) + dtime*dticedt(isi)

          IF ( lis_coupled_to_ocean ) THEN
            ! Coupling to icon-o:
            ! Heat flux (phiipr0 is computed above)
            condhf(isi) = ki*phiipr0*(tice_p(isi)-tf_salt)/hice_p(isi)
            meltpot(isi) = 0._wp
            ! No change of ice thickness
            hice_n(isi) = hice_p(isi)
          ELSE
            ! Compute the time-rate-of-change of the ice thickness (note the sign of heat fluxes)
            dhicedt(isi) = -ki_o_rhoialf*thetaipr0(isi) + qwat*r_rhoialf
            ! Update the ice thickness
            hice_n(isi) = hice_p(isi) + dtime*dhicedt(isi)
          END IF

        END IF

      END IF FreezingMeltingRegime

      ! Remove too thin ice or impose security constraints
      IF( hice_n(isi)<hice_min ) THEN
        ! Remove too thin ice
        hice_n(isi) = 0._wp
        ! Set the ice surface temperature equal to the fresh water freezing point
        tice_n(isi) = tf_fresh
      ELSE
        ! Limit the ice thickness from above (security)
        hice_n(isi) = MIN(hice_n(isi), hice_max)
        ! Limit the ice surface temperature from above (security)
        tice_n(isi) = MIN(tice_n(isi), tf_fresh)
      END IF

      ! Compute tendencies (for eventual use outside the sea-ice scheme program units)
      dticedt(isi)  = (tice_n(isi)-tice_p(isi))*r_dtime
      dhicedt(isi)  = (hice_n(isi)-hice_p(isi))*r_dtime

    END DO GridBoxesWithSeaIce

    ! Compute sea-ice albedo through a rate equation
    PrognosticSeaIceAlbedo: IF ( lprog_albsi ) THEN

      ! Loop over grid boxes where sea ice is present
      !$ACC LOOP GANG(STATIC: 1) VECTOR PRIVATE(albsi_e, albsi_snow_e, taualbsi, rtaualbsisn, albsi_e_wghtd)
      DO isi=1, nsigb

        ! Equilibrium sea-ice albedo (function of sea-ice surface temperature)
        albsi_e = alb_seaice_equil( tice_n(isi), 0._wp )

        ! Equilibrium albedo of snow over sea ice (function of sea-ice surface temperature)
        albsi_snow_e = albsi_snow_max * ( 1.0_wp - c1_albsi_snow                   &
          &                           * EXP(-c2_albsi_snow*(tf_fresh-tice_n(isi))) )

        ! Relaxation time scale for sea-ice albedo
        ! Interpolate linearly between maximum and minimum relaxation time scales
        ! over a given temperaure range
        taualbsi = taualbsi_max + (taualbsi_min-taualbsi_max)  &
                 * ((tice_n(isi)-t_taualbsi_min)*rdelt_taualbsi)
        ! Limit relaxation time scale from below and from above
        taualbsi = MIN(taualbsi_max,MAX(taualbsi,taualbsi_min))
        ! Use temperature-dependent relaxation time scale
        ! if sea-ice albedo tends to decrease,
        ! and a maximum time scale otherwise
        taualbsi = MERGE( taualbsi, taualbsi_max, (albsi_p(isi)>albsi_e) )

        ! Reciprocal of the relaxation time scale for snow-over-ice albedo
        ! Relaxation towards snow-over-ice albedo is only applied if
        ! albedo tends to increase, and
        ! the sea-ice surface temperature is not too close to the freezing point
        rtaualbsisn = snow_rate(isi)*c_tausi_snow*MERGE( 1._wp, 0._wp,  &
          &           ((albsi_p(isi)<albsi_snow_e).AND.(tice_n(isi)<t_albsi_snow_max)) )

        ! Weighted equilibrium albedo
        albsi_e_wghtd = (albsi_e+taualbsi*rtaualbsisn*albsi_snow_e)  &
          &           / (1._wp+taualbsi*rtaualbsisn)

        ! Relax sea-ice albedo towards equilibrium value
        albsi_n(isi) = albsi_e_wghtd+(albsi_p(isi)-albsi_e_wghtd)  &
          &          * EXP(-dtime*(1._wp/taualbsi+rtaualbsisn))

      END DO

    ENDIF PrognosticSeaIceAlbedo
    !$ACC END PARALLEL

    ! Store time tendencies (optional)
    IF (PRESENT(opt_dticedt)) THEN
#ifdef _OPENACC
        CALL finish ('seaice_timestep_nwp', 'OpenACC version currently does not support the optional argument opt_dticedt')
#endif
      opt_dticedt(1:nsigb)  = dticedt(1:nsigb)
      IF (nsigb < SIZE(opt_dticedt)) THEN
        opt_dticedt(nsigb+1:) = 0._wp
      ENDIF
    ENDIF
    IF (PRESENT(opt_dhicedt)) THEN
#ifdef _OPENACC
        CALL finish ('seaice_timestep_nwp', 'OpenACC version currently does not support the optional argument opt_dhicedt')
#endif
      opt_dhicedt(1:nsigb)  = dhicedt(1:nsigb)
      IF (nsigb < SIZE(opt_dhicedt)) THEN
        opt_dhicedt(nsigb+1:) = 0._wp
      ENDIF
    ENDIF

    IF (.NOT. lcuda_graph_lnd) THEN
      !$ACC WAIT(1)
    END IF
    !$ACC END DATA
    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_timestep_nosnow_nwp

  SUBROUTINE seaice_timestep_snow_nwp( &
        & dtime, nsigb, qsen, qlat, qlwrnet, qsolnet, snow_rate, rain_rate, fac_bottom_hflx, &
        & tice_p, hice_p, tsnow_p, hsnow_p, albsi_p, tice_n, hice_n, tsnow_n, hsnow_n, albsi_n, &
        & condhf, meltpot, opt_dticedt, opt_dhicedt, opt_dtsnowdt, opt_dhsnowdt &
      )

    IMPLICIT NONE

    !> model time step [s]
    REAL(wp), INTENT(IN) :: dtime
    !> number of grid boxes within a block where the sea ice is present (<=nproma)
    INTEGER, INTENT(IN) :: nsigb

    REAL(wp), INTENT(IN) :: qsen(:) !< sensible heat flux at the surface [W/m^2]
    REAL(wp), INTENT(IN) :: qlat(:) !< latent heat flux at the surface [W/m^2]
    REAL(wp), INTENT(IN) :: qlwrnet(:) !< net long-wave radiation flux at the surface [W/m^2]
    REAL(wp), INTENT(IN) :: qsolnet(:) !< net solar radiation flux at the surface [W/m^2]

    REAL(wp), INTENT(IN) :: snow_rate(:) !< snow rate (convective + grid-scale) [kg/(m^2 s)]
    REAL(wp), INTENT(IN) :: rain_rate(:) !< rain rate (convective + grid-scale) [kg/(m^2 s)]

    !> tuning factor for heat flux from water to ice [-]
    REAL(wp), OPTIONAL, INTENT(IN) :: fac_bottom_hflx(:)

    REAL(wp), INTENT(IN) :: tice_p(:) !< temperature of ice upper surface at previous time level [K]
    REAL(wp), INTENT(IN) :: hice_p(:) !< ice thickness at previous time level [m]
    REAL(wp), INTENT(IN) :: tsnow_p(:) !< temperature of snow upper surface at previous time level [K]
    REAL(wp), INTENT(IN) :: hsnow_p(:) !< snow thickness at previous time level [m]
    REAL(wp), INTENT(IN) :: albsi_p(:) !< sea-ice albedo at previous time level [-]

    REAL(wp), INTENT(OUT) :: tice_n(:) !< temperature of ice upper surface at new time level [K]
    REAL(wp), INTENT(OUT) :: hice_n(:) !< ice thickness at new time level [m]
    REAL(wp), INTENT(OUT) :: tsnow_n(:) !< temperature of snow upper surface at new time level [K]
    REAL(wp), INTENT(OUT) :: hsnow_n(:) !< snow thickness at new time level [m]
    REAL(wp), INTENT(OUT) :: albsi_n(:) !< sea-ice albedo at new time level [-]

    !> conductive heat flux within the sea ice just above the ice lower boundary [W/m^2]
    REAL(wp), INTENT(OUT) :: condhf(:)
    !> melt potential at top [W/m^2]
    REAL(wp), INTENT(OUT) :: meltpot(:)

    REAL(wp), INTENT(OUT), OPTIONAL :: opt_dticedt(:) !< time tendency of ice surface temperature [K/s]
    REAL(wp), INTENT(OUT), OPTIONAL :: opt_dhicedt(:) !< time tendency of ice thickness [m/s]
    REAL(wp), INTENT(OUT), OPTIONAL :: opt_dtsnowdt(:) !< time tendency of snow surface temperature [K/s]
    REAL(wp), INTENT(OUT), OPTIONAL :: opt_dhsnowdt(:) !< time tendency of snow thickness [m/s]


    LOGICAL :: lis_coupled_to_ocean !< Simulation is coupled to ocean
    LOGICAL :: lpres_fac_hflux !< Heat-flux tuning factor is present
    LOGICAL :: lpres_dticedt !< opt_dticedt is present
    LOGICAL :: lpres_dhicedt !< opt_dhicedt is present
    LOGICAL :: lpres_dtsnowdt !< opt_dtsnowdt is present
    LOGICAL :: lpres_dhsnowdt !< opt_dhsnowdt is present

    REAL(wp) :: r_dtime !< Reciprocal of dtime [1/s]

    REAL(wp) :: qsoliw !< Solar flux at ice-water interface [W/m^2]
    REAL(wp) :: qatm !< Total atmospheric forcing [W/m^2]
    REAL(wp) :: qwat !< Heat flux from water to the melting/freezing zone [W/m^2]

    REAL(wp) :: phiipr0 !< Phi_i'(0), derivative of ice temperature profile [-]
    REAL(wp) :: phiipr1 !< Phi_i'(1), derivative of ice temperature profile [-]
    REAL(wp) :: phispr0 !< Phi_s'(0), derivative of snow temperature profile [-]
    REAL(wp) :: csice !< Ice shape factor [-]
    REAL(wp) :: cssnow !< Snow shape factor [-]
    REAL(wp) :: dcsicedh !< Height derivative of ice shape factor [1/m]
    REAL(wp) :: dcssnowdh !< Height derivative of snow shape factor [1/m]
    REAL(wp) :: dphiipr1dh !< d Phi_i'(1,h) / dh [1/m]
    REAL(wp) :: dphispr0dh !< d Phi_s'(0,h) / dh [1/m]

    REAL(wp) :: rhos !< Snow density [kg/m^3]
    REAL(wp) :: drhosdt !< Time derivative of snow density [kg/(m^3 s)]
    REAL(wp) :: ks !< Snow conductivity [W/(K m)]
    REAL(wp) :: dksdt !< Time derivative of snow conductivity [W/(K m s)]

    REAL(wp) :: theta_s !< Snow surface temperature [K]
    REAL(wp) :: theta_i !< Ice surface temperature [K]
    REAL(wp) :: a_i !< Ice temperature coefficient (theta_s) [-]
    REAL(wp) :: b_i !< Ice temperature offset [K]
    REAL(wp) :: denom !< Denominator of a_i [W/K]
    REAL(wp) :: e_coeff !< Time derivative coefficient of a_i theta_s + b_i (dhsnowdt) [1/m]
    REAL(wp) :: f_coeff !< Time derivative offset of a_i theta_s + b_i [1/s]

    REAL(wp) :: c_eff !< Effective heat capacity of snow-ice slab [J/(m^2 K)]
    REAL(wp) :: c_bnd !< Effective heat capacity of snow-ice boundary [J/(m^2 K)]
    REAL(wp) :: heat_net !< Total heat absorbed by the snow-ice sheet in timestep [J/m^2]

    REAL(wp) :: meltpot_top !< Top melt potential [W/m^2]
    REAL(wp) :: meltpot_bot !< Bottom melt potential [W/m^2]
    REAL(wp) :: flux_bot !< Mass flux at bottom surface [kg/(m^2 s)]

    REAL(wp) :: dhicedt !< Time derivative of ice height [m/s]
    REAL(wp) :: dhsnowdt !< Time derivative of snow height [m/s]
    REAL(wp) :: dhicedt_pred !< Time derivative of ice height w/o corrections [m/s]

    REAL(wp) :: snow_frac !< Diagnosed snow fraction for albedo [-]
    REAL(wp) :: albsi_snow_e !< Equilibrium snow albedo [-]
    REAL(wp) :: albsi_ice_e !< Equilibrium ice albedo [-]
    REAL(wp) :: albsi_e !< Equilibrium albedo [-]
    REAL(wp) :: taualbsi !< Albedo relaxation timescale [s]
    REAL(wp) :: rtaualbsisn !< Reciprocal of snow-fall albedo timescale [1/s]
    REAL(wp) :: albsi_e_wghtd !< Weighted equilibrium albedo after snow fall [-]

    INTEGER :: isi

    r_dtime = 1._wp/dtime
    lis_coupled_to_ocean = is_coupled_to_ocean()

    ! If fac_bottom_hflx is provided, adaptive tuning of the parameter(s)
    ! of the temperature profile within the ice and of the heat flux from water to ice is used
    lpres_fac_hflux = PRESENT(fac_bottom_hflx)

    lpres_dticedt = PRESENT(opt_dticedt)
    lpres_dhicedt = PRESENT(opt_dhicedt)
    lpres_dtsnowdt = PRESENT(opt_dtsnowdt)
    lpres_dhsnowdt = PRESENT(opt_dhsnowdt)

    !$ACC DATA ASYNC(1) &
    !$ACC   PRESENT(qsen, qlat, qlwrnet, qsolnet, snow_rate) &
    !$ACC   PRESENT(tice_p, hice_p, tsnow_p, hsnow_p, albsi_p) &
    !$ACC   PRESENT(tice_n, hice_n, tsnow_n, hsnow_n, albsi_n) &
    !$ACC   PRESENT(condhf, meltpot) &
    !$ACC   NO_CREATE(fac_bottom_hflx, opt_dticedt, opt_dhicedt, opt_dtsnowdt, opt_dhsnowdt)

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG(STATIC: 1) VECTOR &
    !$ACC   PRIVATE(qsoliw, qatm, qwat, phiipr0, phiipr1, phispr0, csice, cssnow) &
    !$ACC   PRIVATE(rhos, drhosdt, ks, dksdt, theta_s, theta_i, a_i, b_i, c_eff, c_bnd) &
    !$ACC   PRIVATE(denom, heat_net, meltpot_top, meltpot_bot, flux_bot) &
    !$ACC   PRIVATE(dhicedt_pred, dhicedt, dhsnowdt) &
    !$ACC   PRIVATE(dcsicedh, dcssnowdh, dphiipr1dh, dphispr0dh, e_coeff, f_coeff)
    ! Loop over grid boxes where sea ice is present
    GridBoxesWithSeaIce: DO isi=1, nsigb

      ! Compute solar radiation flux at the ice-water interface (positive downward)
      qsoliw = qsolnet(isi)*(                                                                 &
             & opticpar_seaice_opaque%frac_optic(1)                                           &
             & *EXP(-MIN(opticpar_seaice_opaque%extincoef_optic(1)*hice_p(isi), cmaxearg)) +  &
             & opticpar_seaice_opaque%frac_optic(2)                                           &
             & *EXP(-MIN(opticpar_seaice_opaque%extincoef_optic(2)*hice_p(isi), cmaxearg)) )

      ! Compute total atmospheric heat flux for the ice slab  (positive downward)
      qatm = qsen(isi) + qlat(isi) + qlwrnet(isi) + qsolnet(isi) - qsoliw

      ! Provision is made to account for the heat flux from water to ice (upward flux is negative).
      ! This is the ocean heat flux just below the ice,
      ! not including the latent flux from melting and freezing.
      ! In the case of coupled icon atmosphere-ocean runs,
      ! this flux has to be passed from the ocean if available.
      ! For uncoupled runs, an ad hoc formulation for the heat flux from water to ice is used
      ! that basically serves to quench the ice growth rate
      ! as the ice thickness approaches its maximum value.
      IF (lbottom_hflux) THEN
        ! Derivative (at zeta=0) of the temperature profile shape function.
        IF (lpres_fac_hflux) THEN
          ! Adaptive tuning of phiipr0 and of bottom heat flux
          phiipr0 = phiipr0_lin - fac_bottom_hflx(isi)*hice_p(isi)/hice_max
        ELSE
          ! A constant value is used (cf. linear temperature profile).
          phiipr0 = phiipr0_lin
        END IF
        ! Heat flux from water to ice is limited from above (no flux from ice to water is allowed).
        qwat = (1._wp-hice_p(isi)/hice_max-phiipr0)*ki*(tf_salt-tice_p(isi))/MAX(hice_p(isi),hice_min)
        qwat = MIN(qwat, 0._wp)
      ELSE
        ! Derivative (at zeta=0) of the temperature profile shape function
        phiipr0 = 1._wp - MIN(hice_p(isi) / hice_max, 1._wp)
        ! Heat flux from water to ice is set to zero.
        qwat = 0._wp
      END IF

      csice = csi_lin - csidp_nlin_sonsi_d * MIN(hice_p(isi) / hice_max, 1._wp)
      cssnow = csi_lin - cssdp_nlin_d * MIN(hsnow_p(isi) / hsnow_max, 0.9_wp)
      phiipr1 = 1._wp + csidp_nlin_sonsi * MIN(hice_p(isi) / hice_max, 1._wp)
      phispr0 = 1._wp - MIN(hsnow_p(isi) / hsnow_max, 0.9_wp)

      ! Fixed snow density.
      rhos = rhos_def
      drhosdt = 0._wp

      ! Snow conductivity parametrization from TERRA.
      ks = 2.22_wp*EXP(1.88_wp*LOG(rhos/rhoi))
      dksdt = 0._wp

      ! Coefficients for theta_i = a_i theta_s + b_i.
      denom = ki * hsnow_p(isi) * phiipr1 + ks * hice_p(isi) * phispr0
      a_i = ks * hice_p(isi) * phispr0 / denom
      b_i = (1._wp - a_i) * tf_salt
      theta_s = MIN(tf_fresh, tsnow_p(isi))
      theta_i = a_i * theta_s + b_i

      ! Ensure that theta_i is between theta_s and tf_salt.
      theta_i = MIN(MAX(MIN(theta_s, tf_salt), theta_i), MAX(theta_s, tf_salt))

      ! effective heat capacity of the snow-ice slab.
      c_eff = rhoi * ci * hice_p(isi) * a_i * csice + rhos * cs * hsnow_p(isi) * (cssnow + a_i * (1._wp - cssnow))

      ! boundary heat capacity, associated with ice top and snow bottom
      c_bnd = rhoi * ci * hice_p(isi) * csice + rhos * cs * hsnow_p(isi) * (1._wp - cssnow)

      ! Bottom melt potential (equal to heat flux through the ice bottom)
      meltpot_bot = ki * (theta_i - tf_salt) / hice_p(isi) * phiipr0
      flux_bot = (meltpot_bot - qwat) / alf

      dphispr0dh = -1._wp / hsnow_max
      dcssnowdh = -cssdp_nlin_d / hsnow_max
      dphiipr1dh = csidp_nlin_sonsi / hice_max
      dcsicedh = -csidp_nlin_sonsi_d / hice_max

      IF (hsnow_p(isi) > 0._wp) THEN
        dhicedt_pred = -flux_bot / rhoi

        ! da_idt theta_s + db_idt =: E dhsnowdt + F
        e_coeff = ((1._wp - a_i) * hice_p(isi) * ks * dphispr0dh - a_i * ki * phiipr1) / denom
        f_coeff = ((1._wp - a_i) * phispr0 * (dhicedt_pred * ks + hice_p(isi) * dksdt) &
            & - a_i * hsnow_p(isi) * ki * dphiipr1dh * dhicedt_pred) / denom

        ! heat absorbed during time step.
        heat_net = dtime * (qatm - qwat &
            & - c_bnd * (theta_s - tf_salt) * (e_coeff / rhos * (snow_rate(isi) - drhosdt * hsnow_p(isi)) + f_coeff) &
            & + (ci * (csice + hice_p(isi) * dcsicedh) * (theta_i - tf_salt) - alf) * flux_bot &
            & + cs * (1._wp - cssnow - hsnow_p(isi) * dcssnowdh) * (theta_s - theta_i) * snow_rate(isi) &
            & + cs * hsnow_p(isi)**2 * (theta_s - theta_i) * dcssnowdh * drhosdt)

        tsnow_n(isi) = MIN(theta_s + heat_net / c_eff, tf_fresh)
        tice_n(isi) = a_i * tsnow_n(isi) + b_i

        heat_net = heat_net - c_eff * (tsnow_n(isi) - theta_s)

        meltpot_top = MAX(0._wp, heat_net) / ( &
            & 1._wp + cs / alf * (1._wp - cssnow - hsnow_p(isi) * dcssnowdh) * ((1._wp - a_i) * tf_fresh - b_i) &
            & - c_bnd * (tf_fresh - tf_salt) * e_coeff / (rhos * alf)) * r_dtime
      ELSE ! no snow
        heat_net = dtime * (qatm - qwat &
            & + (ci * (csice + hice_p(isi) * dcsicedh) * (theta_i - tf_salt) - alf) * flux_bot &
            & + ci * rhoi * ki / (rhos * ks) * (theta_s - tf_salt) * phiipr1 * snow_rate(isi))

        tice_n(isi) = MIN(theta_s + heat_net / c_eff, tf_fresh)
        tsnow_n(isi) = tice_n(isi)

        heat_net = heat_net - c_eff * (tice_n(isi) - theta_s)

        meltpot_top = MAX(0._wp, heat_net) / ( &
            & 1._wp + ci / alf * (1._wp - csice - hice_p(isi) * dcsicedh) * (tf_fresh - tf_salt)) * r_dtime
      END IF

      IF (lis_coupled_to_ocean) THEN
        dhicedt = 0._wp
        dhsnowdt = 0._wp
        hsnow_n(isi) = hsnow_p(isi)
        hice_n(isi) = hice_p(isi)
        meltpot(isi) = meltpot_top
        condhf(isi) = meltpot_bot
      ELSE
        dhsnowdt = MAX(-hsnow_p(isi) * r_dtime, snow_rate(isi) / rhos - meltpot_top / (alf * rhos) - drhosdt / rhos * hsnow_p(isi))
        meltpot_top = MAX(0._wp, meltpot_top - snow_rate(isi) * alf + dhsnowdt * alf * rhos + alf * drhosdt * hsnow_p(isi))
        dhicedt = MAX(-hice_p(isi) * r_dtime, -meltpot_top / (alf * rhoi) - flux_bot / rhoi)

        hsnow_n(isi) = MIN(MAX(0._wp, hsnow_p(isi) + dhsnowdt * dtime), hsnow_max)
        hice_n(isi) = MIN(MAX(0._wp, hice_p(isi) + dhicedt * dtime), hice_max)
      END IF

      IF (lpres_dticedt) opt_dticedt(isi) = (tice_n(isi) - tice_p(isi)) * r_dtime
      IF (lpres_dhicedt) opt_dhicedt(isi) = dhicedt
      IF (lpres_dtsnowdt) opt_dtsnowdt(isi) = (tsnow_n(isi) - tsnow_p(isi)) * r_dtime
      IF (lpres_dhsnowdt) opt_dhsnowdt(isi) = dhsnowdt

      IF (lprog_albsi) THEN
        snow_frac = MAX(0._wp, MIN(hsnow_n(isi)/snow_frac_scale, 1._wp))

        albsi_snow_e = albsi_snow_max &
            & - (albsi_snow_max - albsi_snow_min) * EXP(-r_alb_snow_si_scale * (tf_fresh-tsnow_n(isi)))
        albsi_ice_e = albsi_max &
            & - (albsi_max - albsi_min) * EXP(-r_alb_seaice_scale * (tf_fresh-tsnow_n(isi)))
        albsi_e = (1._wp - snow_frac) * albsi_ice_e + snow_frac * albsi_snow_e

        taualbsi = taualbsi_max + (taualbsi_min-taualbsi_max) &
            & * ((tsnow_n(isi)-t_taualbsi_min) * rdelt_taualbsi)
        taualbsi = MIN(taualbsi_max, MAX(taualbsi, taualbsi_min))

        ! Use temperature-dependent relaxation time scale
        ! if sea-ice albedo tends to decrease,
        ! and a maximum time scale otherwise
        taualbsi = MERGE( taualbsi, taualbsi_max, (albsi_p(isi)>albsi_e) )

        rtaualbsisn = MERGE( snow_rate(isi)*c_tausi_snow, 0._wp,  &
            &           ((albsi_p(isi)<albsi_snow_max).AND.(tsnow_n(isi)<t_albsi_snow_max)) )

        albsi_e_wghtd = (albsi_e + taualbsi * rtaualbsisn * albsi_snow_max) / (1._wp + taualbsi * rtaualbsisn)

        albsi_n(isi) =  albsi_e_wghtd &
            & + (albsi_p(isi) - albsi_e_wghtd) &
            &   * EXP( - dtime * (1._wp / taualbsi + rtaualbsisn))
      END IF
    END DO GridBoxesWithSeaIce
    !$ACC END PARALLEL

    !$ACC END DATA

  END SUBROUTINE seaice_timestep_snow_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Coldstart for sea-ice parameterization scheme.
  !!
  !! Coldstart for sea-ice parameterization scheme. Sea-ice surface temperature and sea-ice
  !! thickness are initialized with meaningful values.
  !! Note that an estimate of the sea-ice temperature is required for the cold start and is
  !! assumed to be available. The only option at the time being is to use the IFS skin
  !! temperature for the cold start initialization of t_ice. Since an estimate of the
  !! ice thickness h_ice is generally not available, h_ice is initialized with a
  !! meaningful constant value (1m).
  !! Note that only "sea" grid boxes are initialized;
  !! the "lake" grid boxes are left intact.
  !!

  SUBROUTINE seaice_coldinit_nwp (                              &
                          &  nswgb,                             &
                          &  frice_thrhld,                      &
                          &  frsi,                              &
                          &  temp_in,                           &
                          &  tice_p, hice_p, tsnow_p, hsnow_p,  &
                          &  albsi_p,                           &
                          &  tice_n, hice_n, tsnow_n, hsnow_n,  &
                          &  albsi_n                            &
                          &  )

    IMPLICIT NONE

    ! Procedure arguments

    INTEGER, INTENT(IN)  ::                &
                        &  nswgb      !< number of "sea" grid boxes within a block (<=nproma)

    REAL(wp), INTENT(IN) :: frice_thrhld     !< fraction threshold for creating a sea grid point

    REAL(wp), DIMENSION(:), INTENT(IN)    ::           &
                                          &  frsi    , &  !< sea-ice fraction [-]
                                          &  temp_in      !< meaningfull guess of ice surface temperature [K]
                                                          !  e.g. tskin from IFS


    REAL(wp), DIMENSION(:), INTENT(INOUT) ::           &
                                          &  tice_p  , &  !< temperature of ice upper surface at previous time level [K]
                                          &  hice_p  , &  !< ice thickness at previous time level [m]
                                          &  tsnow_p , &  !< temperature of snow upper surface at previous time level [K]
                                          &  hsnow_p , &  !< snow thickness at previous time level [m]
                                          &  albsi_p , &  !< sea-ice albedo at previous time level [-]
                                          &  tice_n  , &  !< temperature of ice upper surface at new time level [K]
                                          &  hice_n  , &  !< ice thickness at new time level [m]
                                          &  tsnow_n , &  !< temperature of snow upper surface at new time level [K]
                                          &  hsnow_n , &  !< snow thickness at new time level [m]
                                          &  albsi_n      !< sea-ice albedo at new time level [-]

    ! Local variables

    INTEGER ::      &
            &  isi  !< DO loop index


    REAL(wp), PARAMETER :: h_ice_coldstart = 1.0_wp   ! sea-ice thickness for cold start [m]

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Loop over all grid boxes
    DO isi=1, nswgb

      ! Note that we make use of >= instead of > in order to be consistent
      ! with the seaice index list generation routine
      IF ( frsi(isi) >= frice_thrhld ) THEN  ! ice point

        hice_p(isi)  = h_ice_coldstart             ! constant ice thickness of 1m
        tice_p(isi)  = temp_in(isi)                ! some proper estimate (here: tskin from IFS)
        tice_p(isi)  = MIN(tice_p(isi), tf_fresh)  ! security
        tsnow_p(isi) = tice_p(isi)                ! snow temperature is equal to ice temperature
        hsnow_p(isi) = 0._wp                      ! snow over ice is not treated explicitly
        IF ( lprog_albsi ) THEN                   ! set sea-ice albedo to its equilibrium value
          albsi_p(isi) = alb_seaice_equil( tice_p(isi), hsnow_p(isi) )
        ENDIF

        ! Set variables at new time level
        tice_n(isi)  = tice_p(isi)
        hice_n(isi)  = hice_p(isi)
        tsnow_n(isi) = tsnow_p(isi)
        hsnow_n(isi) = hsnow_p(isi)
        IF ( lprog_albsi ) THEN                   ! set sea-ice albedo to its equilibrium value
          albsi_n(isi) = albsi_p(isi)
        ENDIF
      ENDIF

    END DO ! isi

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_coldinit_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Cold start initialization of prognostic sea-ice albedo.
  !!
  !! This routine is used when
  !! cold start initialization of prognostic sea-ice albedo is necessary,
  !! whereas the sea-ice surface temperature and the sea-ice thickness
  !! should not be (re-)initilaized.
  !! This occurs, for example, if the sea-ice scheme has already been used,
  !! but the sea-ice albedo from previous runs is not available.
  !! The sea-ice albedo is initialized with its equilibrium value
  !! that is a function of sea-ice surface temperature.
  !! Note that only "sea" grid boxes are initialized;
  !! the "lake" grid boxes are left intact.
  !!

  SUBROUTINE seaice_coldinit_albsi_nwp (              &
                                    &  nswgb,         &
                                    &  frice_thrhld,  &
                                    &  frsi,          &
                                    &  tsnow_p,       &
                                    &  hsnow_p,       &
                                    &  albsi_p,       &
                                    &  albsi_n        &
                                    &  )

    IMPLICIT NONE

    ! Procedure arguments

    INTEGER, INTENT(IN) ::          &
                        &  nswgb      !< number of "see" grid boxes within a block (<=nproma)

    REAL(wp), INTENT(IN) :: frice_thrhld     !< fraction threshold for creating a sea grid point

    REAL(wp), DIMENSION(:), INTENT(IN)    ::           &
                                          &  frsi    , &  !< sea-ice fraction [-]
                                          &  tsnow_p , &  !< temperature of upper surface at previous time level [K]
                                          &  hsnow_p      !< height of snow at previous time level [m]


    REAL(wp), DIMENSION(:), INTENT(INOUT) ::           &
                                          &  albsi_p , &  !< sea-ice albedo at previous time level [-]
                                          &  albsi_n      !< sea-ice albedo at new time level [-]
    ! Local variables

    INTEGER ::      &
            &  isi  !< DO loop index

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------

    ! Loop over sea-water grid boxes
    DO isi=1, nswgb

      ! Note that we make use of >= instead of > in order to be consistent
      ! with the seaice index list generation routine
      IF ( frsi(isi) >= frice_thrhld ) THEN  ! ice point

        ! set sea-ice albedo to its equilibrium value
        albsi_p(isi) = alb_seaice_equil( tsnow_p(isi), hsnow_p(isi) )

        ! set albedo at new time level
        albsi_n(isi) = albsi_p(isi)

      ENDIF

    END DO ! isi

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END SUBROUTINE seaice_coldinit_albsi_nwp

!234567890023456789002345678900234567890023456789002345678900234567890023456789002345678900234567890

  !>
  !! Equilibrium sea-ice albedo is computed
  !! as function of the sea-ice surface temperature.
  !!

  REAL (wp) FUNCTION alb_seaice_equil ( t_top, h_snow )
    !$ACC ROUTINE SEQ

    IMPLICIT NONE

    ! Procedure arguments

    REAL(wp), INTENT(IN) :: t_top       !< temperature of ice/snow upper surface [K]
    REAL(wp), INTENT(IN) :: h_snow      !< snow height [m]

    REAL(wp) :: snow_frac
    REAL(wp) :: albsi_snow_e
    REAL(wp) :: albsi_ice_e

    !===============================================================================================
    !  Start calculations
    !-----------------------------------------------------------------------------------------------


    ! albsi_max and albsi_min are defined in the namelist (mo_lnd_nwp_nml).
    ! A derived constant 0.35 is equal to 95.6/tf_fresh, where tf_fresh=273.15 K,
    ! and has a dimensions of K^{-1}.
    IF (lsnow_on_seaice) THEN
      snow_frac = MAX(0._wp, MIN(h_snow/snow_frac_scale, 1._wp))
      albsi_snow_e = albsi_snow_max &
          & - (albsi_snow_max - albsi_snow_min) * EXP(-r_alb_snow_si_scale * (tf_fresh-t_top))
      albsi_ice_e = albsi_max &
          & - (albsi_max - albsi_min) * EXP(-r_alb_seaice_scale * (tf_fresh-t_top))
      alb_seaice_equil = (1._wp - snow_frac) * albsi_ice_e + snow_frac * albsi_snow_e
    ELSE
      alb_seaice_equil = albsi_max-(albsi_max-albsi_min)*EXP(-0.35_wp*(tf_fresh-t_top))
    END IF

    !-----------------------------------------------------------------------------------------------
    !  End calculations
    !===============================================================================================

  END FUNCTION alb_seaice_equil

END MODULE sfc_seaice

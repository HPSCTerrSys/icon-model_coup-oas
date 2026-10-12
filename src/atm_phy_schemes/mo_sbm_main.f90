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

! Description:
!
! this is the spectral-bin microphysics scheme based on the hebrew university
! cloud model (hucm), originally formulated and coded by alexander khain
!
! the wrf bin microphysics scheme (fast sbm or fsbm) solves equations for four
! size distribution FUNCTIONs: aerosols, drop (including rain drops), snow and
! graupel/hail (from which mass mixing ratio qna, qc, qr, qs, qg/qh and
! their number concentrations are calculated).

! the scheme is generally written in cgs units. in the updated scheme (fsbm-2)
! the users can choose either graupel or hail to describe dense particles
! (see the 'hail_opt' switch). by default, the 'hail_opt = 1' is used.
! hail particles have larger terminal velocity than graupel per mass bin.
! 'hail_opt' is recommended to be used in simulations of continental clouds
! systems. the graupel option may lead to better results in simulations of
! maritime convection.

! the aerosol spectrum in fsbm-2 is approximated by 3-lognormal size distribution
! representing smallest aerosols (nucleation mode), intermediate-size
! (accumulation mode) and largest aerosols (coarse mode). the bc/ic for aerosols
! ,as well as aerosols vertical distribution profile -- are set from within the
! fsbm-2 scheme (see the 'domain_id' parameter). the domain_id forces bc to be applied
! for the parent domain only.

! the user can set the liquid water content threshold (lwc) in which rimed snow
! is being transferred to hail/graupel (see 'alcr' parameter).
! the default value is alcr = 0.5 [g/m3]. increasing this value will result
! in an increase of snow mass content, and a decrease in hail/graupel mass
! contents.

! we thank and acknowledge contribution from jiwen fan (pnnl), alexander rhyzkov
! (cimms/nssl), jeffery snyder (cimms/nssl), jimy dudhia (ncar) and dave gill! (ncar).

! useful references:
! -------------------
!   khain, a. p., and i. sednev, 1996: simulation of precipitation formation in
! the eastern mediterranean coastal zone using a spectral microphysics cloud
! ensemble model. atmospheric research, 43: 77-110;
!   khain, a. p., a. pokrovsky and m. pinsky, a. seifert, and v. phillips, 2004:
! effects of atmospheric aerosols on deep convective clouds as seen from
! simulations using a spectral microphysics mixed-phase cumulus cloud model
! part 1: model description. j. atmos. sci 61, 2963-2982);
!   khain a. p. and m. pinsky, 2018: physical processes in clouds and cloud
! modeling. cambridge university press. 642 pp
!   shpund, j., a. khain, and d. rosenfeld, 2019: effects of sea spray on the
! dynamics and microphysics of an idealized tropical cyclone. j. atmos. sci., 0,
! https://doi.org/10.1175/jas-d-18-0270.1 (a preliminary description of the
! updated fsbm-2 scheme)

! when using the fsbm-2 version please cite:
! -------------------------------------------
! shpund, j., khain, a., lynn, b., fan, j., han, b., ryzhkov, a., snyder, j.,
! dudhia, j. and gill, d., 2019. simulating a mesoscale convective system using wrf
! with a new spectral bin microphysics: 1: hail vs graupel.
! journal of geophysical research: atmospheres.
! note by jiwen fan
! (1) the main SUBROUTINE is fast_sbm where all the microphysics processes are
!     called
! (2) for aerosol setup, seach "aerosol setup", where one can set up aerosol sd,
! composition information (molecular weight, ions, and density). for sd, there
! is a choice for a lognormal distribution, or read from an observed sd.
! (3) my postdoc yuwei zhang has added cloud related diagnostics (mainly process
! rates) and added an option to read in the observed sd. observed sd data should be
! processed following! a format of the file "ccn_size_33bin.dat" which is in
! size (cm), dn (# cm-3), and dndlogd for 33bins

MODULE mo_sbm_main

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, txt => message_text
  USE mo_run_config,      ONLY: iqb_water_start, iqb_water_end,    &
       &                      iqb_snow_start, iqb_snow_end,      &
       &                      iqb_graupel_start, iqb_graupel_end,&
       &                      iqb_ccn_start, iqb_ccn_end
  USE mo_thdyn_functions, ONLY: sat_pres_water, sat_pres_ice
  USE mo_math_constants,  ONLY: pi
  USE mo_sbm_util,        ONLY: &
    use_cloud_base_nuc, t_nucl_drop_min, isign_3point,             &
    coeff_remaping, ratio_icew_min, rw_pw_min, rw_pw_ri_pi_min,    &
    ventpl_max, epsil, chucm, col, nkr, icemax, ibreakup, ima,     &
    jbreak, krmin_breakup, pkij, qkj, convert_micro2advect,        &
    mwaero, ions, krdrop, ncondcoll, rccn, rlec, ro_solute,     &
    iceprocs, stick_param1,stick_param2, vr1, vr2, vr3, vr4, vr5,  &
    riec, rsec, rgec, rhec, ro1bl, ro2bl, ro3bl, ro4bl, ro5bl,     &
    hail_opt, melt_on,breakup_rain_spont_on,                       &
    snow_breakup_on, kr_snow_min, ikr_spon_break,xl,xi,xs,xh,xg,   &
    xl_mg,xs_mg,xg_mg,xh_mg,prob,gain_var_new,nnd,krice,           &
    alcr,alcr_g, cwll_all, cwls_all,cwsl_all,cwss_all,             &
    cwlg_all,cwgl_all,cwlh_all,cwhl_all,cwsg_all,cwgs_all,         &
    cwsh_all, cwhs_all, cwgg_all, cwhh_all, p_z_up, p_z_del,       &
    t_nucl_ice_min,isign_tq_icenucl,delsupice_max,                 &
    i_break_method, icempl, kr_icempl,                             &
    & rhoh2o,  & ! [kg/m3]  density of liquid water
    & alv,     & ! [J/kg]   latent heat of vaporization
    & als,     & ! [J/kg]   latent heat of sublimation
    & alf,     & ! [J/kg]   latent heat for fusion (freezing)
    & cpd,     & ! [J/K/kg] specific heat at constant pressure
    & cvd,     & ! [J/K/kg] specific heat at constant volume
    & rhoice,  & ! [kg/m3]  density of pure ice
    & tmelt,   & ! [K]      melting temperature of ice/snow
    & rv,      & ! [J/K/kg] gas constant for water vapor
    & grav,    & ! [m/s2] av. gravitational acceleration
    & rd,      & ! [J/K/kg] gas constant (dry air)
    & ccnconstarr, use_ccn_const, &
    & latheatfac1, latheatfac2, latheatfac3, latheatfac4, &
    & tune_fall, tune_long_relax, &
    & usenkrp1a,usenkrp1b,usenkrp2,usenkrp3, & ! parameters for 34th pseodo bin
    & tune_melt_factor, & !artificial factor for increasing melt rate of snow, graupel and hail
    & snha2ha, positive_t_coll ! allow collisions with ice also at positive temperatures

  IMPLICIT NONE
  PRIVATE
  PUBLIC fast_sbm
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sbm_main'

  REAL(KIND=wp), PARAMETER :: zdv=2.22e-5_wp           ! molecular diffusion coefficient for water vapour
  REAL(KIND=wp), PARAMETER :: rv_cgs=rv*10000.0_wp     ! cgs: cm^2/s^2/K
  REAL(KIND=wp), PARAMETER :: alv_cgs=alv*10000.0_wp   ! cgs: cm^2/s^2
  REAL(KIND=wp), PARAMETER :: als_cgs=als*10000.0_wp   ! cgs: cm^2/s^2
  REAL(KIND=wp), PARAMETER :: al1=alv/cpd              ! K
  REAL(KIND=wp), PARAMETER :: al2=als/cpd              ! K
  REAL(KIND=wp), PARAMETER :: zdv_cgs=zdv*10000.0_wp   ! cgs: cm^2/s
  REAL(KIND=wp), PARAMETER :: eps_diff=1.0e-16_wp      ! minimal sum over all bins, above which we do diffusional growth
  REAL(KIND=wp), PARAMETER :: eps_diff_bin=1.0e-10_wp  ! minimal maximum PSD bin value above which we do diffusional growth
  REAL(KIND=wp), PARAMETER :: temp_no_diff=213.15_wp   ! minimum temperature for nucleation or diffisional growth
  REAL(KIND=wp), PARAMETER :: latheat_freez=alf/cvd    ! latent heat release of freezing [J/g] (was before alf/cpd=334000/1004=334)
  REAL(KIND=wp), PARAMETER :: temp_stick_lim=233.15_wp ! below that limit, ice-ice do not stick in collission
  REAL(KIND=wp), PARAMETER :: epsdel=0.1e-03_wp
  REAL(KIND=wp), PARAMETER :: aa1_my=2.53e12_wp        ! 2.53e11 Pa = 2.53e12 dynes/cm/cm (1Pa=10dynes/cm/cm)
  REAL(KIND=wp), PARAMETER :: bb1_my=5.42e3_wp
  REAL(KIND=wp), PARAMETER :: aa2_my=3.41e13_wp
  REAL(KIND=wp), PARAMETER :: bb2_my=6.13e3_wp
  REAL(KIND=wp), PARAMETER :: pzero=1.01325e6_wp       ! reference pressure (cgs), i.e. dynes/cm/cm
  REAL(KIND=wp), PARAMETER :: cf_my=2.4e3_wp           ! coefficient of thermal conductivity of air (mks), i.e. g*cm/sec/sec/sec/kelvin
  REAL(KIND=wp), PARAMETER :: const_dmdt=12.566372_wp
  REAL(KIND=wp), PARAMETER :: bfreezmax=0.66_wp
  REAL(KIND=wp), PARAMETER :: ttcoal=233.15_wp
  REAL(KIND=wp), PARAMETER :: afreezmy=0.3333e-04_wp
  REAL(KIND=wp), PARAMETER :: bfreezmy=0.66_wp
  INTEGER, PARAMETER       :: krfreez=21
  REAL(KIND=wp) :: lh_ce_1, ce_bf, ce_af, cldnucl_af, cldnucl_bf, del_cldnucl_sum, del_ce_sum !processes rates (golbal)

  CONTAINS

  SUBROUTINE fast_sbm (dt                &!in:    dt
                      ,dz8w              &!in:    vertical layer thickness
                      ,rho_phy           &!in:    density
                      ,p_phy             &!in:    pressure
                      ,pi_phy            &!in:    exner
                      ,w                 &!in:    updraft velocity
                      ,qv_old            &!in     cloud water before advection
                      ,th_phy            &!inout: theta. Check how to update prognostic theta_v ???
                      ,qv                &!inout: specific humidity
                      ,chem_new          &!inout: 99 mass bins
                      ,prec_r_sbm        &!inout: 1 time step precip. (mm/sec): input: 0, output can go further to the model
                      ,qc                &!inout: cloud water:                   input: 0
                      ,qr                &!inout: rain water:                    input: 0
                      ,qnc               &!inout: cloud water concentration:     input: 0
                      ,qnr               &!inout: rain water concentration:      input: 0
                      ,qna               &!inout: ccn concentration:             input: 0
                      ,lh_rate           &!inout: rate 1:                        input: 0, output can    go further to the model
                      ,ce_rate           &!inout: rate 2:                        input: 0, output can    go further to the model
                      ,cldnucl_rate      &!inout: rate 3:                        input: 0, output can    go further to the model
                      ,its,ite, kts,kte  &!in:    subdomain indeces
                      ,diag_satur_ba     & !inout: diagnostic supersaturation before advection
                      ,diag_satur_aa     & !inout: diagnostic supersaturation after advection
                      ,diag_satur_am     & !inout: diagnostic supersaturation after microphysics
                      ,diag_supsat_out   & !inout: diagnostic supersaturation after cond_evap subroutines
                      ,temp_old &
                      ,temp_new,reff,reffc,reffr &
                      ,qi &
                      ,qs &
                      ,qg &
                      ,qni &
                      ,qns &
                      ,qng &
                      ,prec_s_sbm &
                      ,prec_g_sbm)
    IMPLICIT NONE
    INTEGER,INTENT(IN)                            :: its,ite,kts,kte
    REAL(KIND=wp), INTENT(IN)                     :: dt
    REAL(KIND=wp), DIMENSION(:,:), INTENT(IN)     :: w, qv_old, &
                 temp_old, temp_new, dz8w, p_phy,pi_phy,rho_phy
    REAL(KIND=wp), DIMENSION(:,:,:),INTENT(INOUT) :: chem_new
    REAL(KIND=wp), DIMENSION(:,:), INTENT(INOUT)  :: qc, qnc,  &
                 qr,qnr,qna,reff,reffc,reffr,qi, qs, qg, qni,  &
                 qns, qng, lh_rate,ce_rate,cldnucl_rate,       &
                 th_phy, qv,diag_satur_ba, diag_satur_aa,      &
                 diag_satur_am, diag_supsat_out
    REAL(KIND=wp), INTENT(INOUT), DIMENSION(:) :: prec_r_sbm,  &
                 prec_s_sbm,prec_g_sbm
    INTEGER :: kr,ikl,ice,i,k, flag_condevap,krr,              &
                 isym1,isym2(icemax),isym3,isym4,isym5,        &
                 is_this_cloudbase, kz_cloud_base(its:ite),    &
                 k_found, microphys_flag, ncondcoll2
    REAL(KIND=wp), DIMENSION(its:ite, kts:kte)    :: qv_oldtmp,&
                 t_new,t_old,zcgs,rhocgs,pcgs
    REAL(KIND=wp) :: mom_2c,mom_3c,mom_2r,mom_3r,mom_2,mom_3,  &
                 supsat_out, sup2_old, dtcond, dtnew, dtcoll,  &
                 tt,qq,tta,qqa,pp,div1,div2,         &
                 del1in,del2in,z_full,           &
                 fmax1,fmax2(icemax),fmax3,fmax4,fmax5
    REAL(KIND=wp), DIMENSION(nkr) :: fccn, fccn_nucl, ff1in,   &
                 ff1r,ff1r_d,ff3in,ff4in,ff5in,ff3r,ff4r,ff5r
    REAL(KIND=wp),DIMENSION (nkr,icemax) :: ff2in,ff2r
    REAL(KIND=wp) :: zcgs_z(kts:kte),pcgs_z(kts:kte),          &
                 rhocgs_z(kts:kte),ffx_z(kts:kte,nkr),         &
                 vrx(kts:kte,nkr),vr1_z(nkr), factor_p,        &
                 vr1_z3d(nkr,its:ite,kts:kte),vr2_z(nkr,icemax), &
                 vr3_z(nkr),vr4_z(nkr), vr5_z(nkr),            &
                 vr3_z3d(nkr,its:ite,kts:kte),                 &
                 vr4_z3d(nkr,its:ite,kts:kte),vr5_z3d(nkr,its:ite,kts:kte), &
                 w_stag_my, satur_diag, pcgs_loc, rhocgs_loc,  &
                 rho_phy_loc, delta_tt, delta_qq, divfac

    IF ( iceprocs==0 ) THEN
      qi=0.0_wp
      qs=0.0_wp
      qg=0.0_wp
      qni=0.0_wp
      qns=0.0_wp
      qng=0.0_wp
      prec_s_sbm=0.0_wp
      prec_g_sbm=0.0_wp
    END IF

    dtcond = MAX(tune_long_relax*dt/REAL(ncondcoll),1.0_wp)
    dtcoll = dtcond
    ncondcoll2 = INT(dt/dtcond)

!   dt ncondcoll      dtcond ncondcoll2
!   6  6         -->  1      1
!   3  6         -->  1      3
!   1  6         -->  1      1
!   6  3         -->  2      3
!   1  3         -->  1      1

    DO k = kts,kte
      DO i = its,ite
        qv_oldtmp(i,k) = qv_old(i,k)
        t_old(i,k) = temp_old(i,k)
        t_new(i,k) = temp_new(i,k)
        rhocgs(i,k) = rho_phy(i,k)*0.001_wp
      END DO
    END DO

    ! ... drops
    krr=0
    DO kr = iqb_water_start, iqb_water_end
      krr=krr+1
      divfac = 1.0_wp / (col*xl(krr)*xl(krr)*3.0_wp)
      DO k = kts,kte
        DO i = its,ite
          chem_new(i,k,kr)=chem_new(i,k,kr)*rhocgs(i,k)*divfac ! input from model: kg/kg. inside: #/(g*cm^3)
        END DO
      END DO
    END DO
    ! ... snow
    IF ( iceprocs==1 ) THEN
      krr=0
      DO kr = iqb_snow_start, iqb_snow_end
        krr=krr+1
        divfac = 1.0_wp / (col*xs(krr)*xs(krr)*3.0_wp)
        DO k = kts,kte
          DO i = its,ite
            chem_new(i,k,kr)=chem_new(i,k,kr)*rhocgs(i,k)*divfac
          END DO
        END DO
      END DO
    END IF
    ! ... aerosols
    krr=0
    DO kr = iqb_ccn_start, iqb_ccn_end
      krr=krr+1
      ! In the original mixed code there was .../col here, but no *col in the output of chem_new. This was wrong. Email Koby
      IF ( use_ccn_const == 1 ) THEN
        DO k = kts,kte
          DO i = its,ite
            chem_new(i,k,kr) = ccnconstarr(k,krr)*rhocgs(i,k)/1000.0_wp !use time-constant ccn profile
          END DO                                                                !i.e. initialize it again every microphysical time step
        END DO
      ELSE
        DO k = kts,kte
          DO i = its,ite
            chem_new(i,k,kr) = chem_new(i,k,kr)*rhocgs(i,k)/1000.0_wp   !input from model: #/kg. inside: #/cm^3
          END DO                                                                !i.e. initialize it again every microphysical time step
        END DO
      END IF
    END DO
    ! ... hail or graupel [same registry adresses]
    IF ( iceprocs==1 ) THEN
      IF (hail_opt == 1) THEN
        krr=0
        DO kr=iqb_graupel_start, iqb_graupel_end
          krr=krr+1
          divfac = 1.0_wp / (col*xh(krr)*xh(krr)*3.0_wp)
          DO k = kts,kte
            DO i = its,ite
              chem_new(i,k,kr)=chem_new(i,k,kr)*rhocgs(i,k)*divfac
            END DO
          END DO
        END DO
      ELSE
        krr=0
        DO kr=iqb_graupel_start, iqb_graupel_end
          krr=krr+1
          divfac = 1.0_wp / (col*xg(krr)*xg(krr)*3.0_wp)
          DO k = kts,kte
            DO i = its,ite
              chem_new(i,k,kr)=chem_new(i,k,kr)*rhocgs(i,k)*divfac
            END DO
          END DO
        END DO
      END IF
    END IF

    DO i = its,ite
      z_full=0.0_wp
      DO k = kte,kts,-1 ! was DO k = kts,kte before reversing z axis
        pcgs(i,k)=p_phy(i,k)*10.0_wp !Pa --> dyn/cm^2
        zcgs(i,k)=z_full+0.5_wp*dz8w(i,k)*100.0_wp
        z_full=z_full+dz8w(i,k)*100.0_wp
      END DO
    END DO

    ! find cloud base using sup.sat over water:
    DO i = its,ite
      kz_cloud_base(i) = 0
      k_found = 0
      DO k = kte,kts,-1 !loop from 65 up to 2
        satur_diag=0.0_wp
        CALL satcalc_water(t_new(i,k),qv(i,k),satur_diag,rho_phy(i,k)) !supersaturation over water
        IF (k < kte) THEN
          w_stag_my = 50.0_wp*(w(i,k)+w(i,k+1)) !updraft in cm/s
        ELSE
          w_stag_my = 100.0_wp*w(i,k) !updraft in cm/s
        END IF
        ! look for cloud base at the lowest 3km:
        IF (satur_diag > 0.0_wp .AND. w_stag_my > 0.1_wp*1.0e2_wp .AND. k_found == 0 .AND. k < kte .AND. zcgs(i,k) < 3.0_wp*1.0e5_wp) THEN
          kz_cloud_base(i) = k ! k-level index of cloud base
          k_found = 1
        END IF
      END DO
    END DO

    DO k = kts,kte
      DO i = its,ite
        pcgs_loc=pcgs(i,k)
        rhocgs_loc=rhocgs(i,k)
        rho_phy_loc=rho_phy(i,k)

        ! ... correcting look-up-table terminal velocities (stronger with height):
        factor_p = SQRT(1.0e6_wp/pcgs_loc)        !10^6 dyn/cm2: The fall vel. depends on 1/SQRT(pressure)
        vr1_z(1:nkr) =  vr1(1:nkr)*factor_p*tune_fall    ! water
        vr1_z3d(1:nkr,i,k) = vr1(1:nkr)*factor_p*tune_fall ! water
        IF ( iceprocs==1 ) THEN
          vr2_z(1:nkr,1) = vr2(1:nkr,1)*factor_p  ! crystal 1 plates (check?)
          vr2_z(1:nkr,2) = vr2(1:nkr,2)*factor_p  ! crystal 2 dendrites (check?)
          vr2_z(1:nkr,3) = vr2(1:nkr,3)*factor_p  ! crystal 3 columns (check?)
          vr3_z(1:nkr) = vr3(1:nkr)*factor_p     ! snow
          vr4_z(1:nkr) = vr4(1:nkr)*factor_p*tune_fall     ! graupel
          vr5_z(1:nkr) = vr5(1:nkr)*factor_p*tune_fall     ! hail
          vr3_z3d(1:nkr,i,k) = vr3(1:nkr)*factor_p ! snow
          vr4_z3d(1:nkr,i,k) = vr4(1:nkr)*factor_p*tune_fall ! graupel
          vr5_z3d(1:nkr,i,k) = vr5(1:nkr)*factor_p*tune_fall ! hail
        END IF
        ! ... droplet / drops
        krr = 0
        DO kr = iqb_water_start, iqb_water_end
          krr = krr + 1
          ff1r(krr) = chem_new(i,k,kr) ! #/(g*cm^3)
          IF (ff1r(krr) < 0.0_wp) ff1r(krr) = 0.0_wp
        END DO
        ! ... ccn
        krr = 0
        DO kr=iqb_ccn_start, iqb_ccn_end
          krr = krr + 1
          fccn(krr) = chem_new(i,k,kr)
          IF (fccn(krr) < 0.0_wp) fccn(krr) = 0.0_wp
        END DO
        ! ... nucleated ccn - not implemented here

        DO kr=1,nkr
          fccn_nucl(kr)=0.0_wp
        END DO

        lh_ce_1 = 0.0_wp
        ce_bf = 0.0_wp
        ce_af = 0.0_wp
        cldnucl_af = 0.0_wp
        cldnucl_bf = 0.0_wp
        del_cldnucl_sum = 0.0_wp; del_ce_sum = 0.0_wp

        IF ( iceprocs==1 ) THEN
          ! no explicit ice crystals in fsbm (plates, dendrites, columns)
          ! temporarly used for ice nucleation, but THEN transformed to snow bins (and not advected)
          ff2r(:,1) = 0.0_wp
          ff2r(:,2) = 0.0_wp
          ff2r(:,3) = 0.0_wp

          ! ... snow
          krr=0
          DO kr=iqb_snow_start, iqb_snow_end
            krr=krr+1
            ff3r(krr)=chem_new(i,k,kr) ! #/(g*cm^3)
            IF (ff3r(krr) < 0.0_wp) ff3r(krr) = 0.0_wp
          END DO

          ! ... hail or graupel: the chem_new bins are inserted here to graupel ff4r or hail ff5r PSD:
          IF (hail_opt == 1) THEN
            krr=0
            DO kr=iqb_graupel_start, iqb_graupel_end
              krr=krr+1
              ff5r(krr) = chem_new(i,k,kr) ! #/(g*cm^3)
              IF (ff5r(krr) < 0.0_wp) ff5r(krr) = 0.0_wp
              ff4r(krr) = 0.0_wp
            END DO
          ELSE
            krr=0
            DO kr=iqb_graupel_start, iqb_graupel_end
              krr=krr+1
              ff4r(krr) = chem_new(i,k,kr) ! #/(g*cm^3)
              IF (ff4r(krr) < 0.0_wp) ff4r(krr) = 0.0_wp
                ff5r(krr) = 0.0_wp
            END DO
          END IF
        END IF

!       check mass conservation:
!       tot_water(i,k,1)=0.0_wp
!       CALL total_water(qv(i,k),ff1r,ff3r,ff4r,ff5r,tot_water(i,k,1),rhocgs_loc)

        ! +---------------------------------------------+
        ! nucleation, condensation, collisions
        ! +---------------------------------------------+
        satur_diag=0.0_wp
        CALL satcalc_water(t_old(i,k),qv_oldtmp(i,k),satur_diag,rho_phy_loc)
        diag_satur_ba(i,k)=satur_diag
        CALL satcalc_water(t_new(i,k),qv(i,k),satur_diag,rho_phy_loc)
        diag_satur_aa(i,k)=satur_diag

        tt=t_old(i,k)  !temperature before advection
        qq=qv_oldtmp(i,k) !specific humidity before advection kg/kg
        IF (qq <= 0.0_wp) qq = 1.e-10_wp
        pp=pcgs_loc
        tta=t_new(i,k) !temperature after advection
        qqa=qv(i,k)    !specific humidity after advection
        qqa=MAX(qqa,1.0e-10_wp)

        delta_tt=(tta-tt)/ncondcoll2
        delta_qq=(qqa-qq)/ncondcoll2
        dtnew = 0.0_wp

        DO ikl=1,ncondcoll2 !main loop over condensation/collisions substep, for given i,k
          dtcond = min(tune_long_relax*dt-dtnew,dtcond)
          dtnew = dtnew + dtcond
          sup2_old=0.0_wp
          CALL satcalc_ice(tt,qq,sup2_old,rho_phy_loc)
          tt=tt+delta_tt
          qq=qq+delta_qq
          del1in=0.0_wp
          CALL satcalc_water(tt,qq,del1in,rho_phy_loc)
          del2in=0.0_wp
          CALL satcalc_ice(tt,qq,del2in,rho_phy_loc)
          div1=del1in+1.0_wp
          div2=del2in+1.0_wp

          microphys_flag=0

          IF (((del1in > 0.0_wp) .OR. (del2in > 0.0_wp) .OR. &
            (sum(ff1r)+sum(ff3r)+sum(ff4r)+sum(ff5r) > eps_diff)) .AND. &
            .NOT. ((div1 >= div2 .AND. tt <= 272.0_wp) .OR. (tt < temp_no_diff))) THEN !nucleation+diffusional growth

            microphys_flag=1

            ! +------------------------------------------+
            ! droplet or ice nucleation :
            ! +------------------------------------------+
            ff1in(:) = ff1r(:)

            IF ( iceprocs==1 ) THEN
              DO kr=1,nkr
                DO ice=1,icemax !3 crystal types
                  ff2in(kr,ice)=ff2r(kr,ice)
                END DO
              END DO
            END IF
            ! This is for nucleated mass rate diagnostics - was in warm code but not here:
            cldnucl_bf = 3.0_wp*col*( sum(ff1in*(xl**2.0_wp)) )/rhocgs_loc

            !max supsat calculation method:
            is_this_cloudbase = 0 !simple use of local supsat as in warm code
            w_stag_my = 0.0_wp
            IF (use_cloud_base_nuc > 0) THEN !find max supsat (~20m) above c.base
                                             !for better estimation of nucleated drops
                                             !concentration. Currently off
              IF (kz_cloud_base(i) == k .AND. col*sum(ff1in*xl) < 5.0_wp) THEN ! condition of low LWC (no rain)
                is_this_cloudbase = 1
              END IF
              IF (k < kte) THEN
                w_stag_my = 50.0_wp*(w(i,k)+w(i,k+1)) !updraft in cm/s
              ELSE
                w_stag_my = 100.0_wp*w(i,k) !updraft in cm/s
              END IF
            END IF

            ! water and ice nucleation:
            CALL nucleation_main(ff1in,ff2in,fccn,fccn_nucl,tt,qq,pcgs_loc &
                         ,del1in,del2in,sup2_old,w_stag_my,is_this_cloudbase,rho_phy_loc)

            !This is meaningless nucleated mass rate diagnostics that was in warm code but not here:
            !It would be interesting diagnostics IF number concentration (set xl insteadof xl**2) nucleation rate change is calculated.
            cldnucl_af = 3.0_wp*col*( sum(ff1in*(xl**2.0_wp)) )/rhocgs_loc
            del_cldnucl_sum = del_cldnucl_sum + (cldnucl_af - cldnucl_bf)

            DO kr=1,nkr
              ff1r(kr)=ff1in(kr)
              IF ( iceprocs==1 ) THEN
                DO ice=1,icemax
                  ff3r(kr) = ff3r(kr) + ff2in(kr,ice) ! adding 3 crystal types to the snow PSD
                  ff2in(kr,ice) = 0.0_wp
                  ff2r(kr,ice) = 0.0_wp
                END DO
              END IF
            END DO
            ! Here, after nucleation there is no ice ff2r PSD. It is moved to snow PSD.
            ! graupel/hail PSD can exist as input to FAST_SBM, and is not changed in nucleation above

            fmax1=0.0_wp
            IF ( iceprocs==1 ) THEN
              fmax2=0.0_wp
              fmax3=0.0_wp
              fmax4=0.0_wp
              fmax5=0.0_wp
            END IF
            DO kr=1,nkr !find max PSDs values
              ff1in(kr)=ff1r(kr)
              fmax1=MAX(ff1r(kr),fmax1)
              IF ( iceprocs==1 ) THEN
                ff3in(kr)=ff3r(kr)
                fmax3=MAX(ff3r(kr),fmax3)
                ff4in(kr)=ff4r(kr)
                fmax4=MAX(ff4r(kr),fmax4)
                ff5in(kr)=ff5r(kr)
                fmax5=MAX(ff5r(kr),fmax5)
                DO ice=1,icemax
                  ff2in(kr,ice)=ff2r(kr,ice)
                  fmax2(ice)=MAX(ff2r(kr,ice),fmax2(ice))
                END DO
              END IF
            END DO

            isym1=0
            isym2=0
            isym3=0
            isym4=0
            isym5=0
            !check IF there is "something" of given type at the grid point:
            IF (fmax1 > eps_diff_bin) isym1 = 1
            IF ( iceprocs==1 ) THEN
              IF (fmax2(1) > eps_diff_bin) isym2(1) = 1
              IF (fmax2(2) > eps_diff_bin) isym2(2) = 1
              IF (fmax2(3) > eps_diff_bin) isym2(3) = 1
              IF (fmax3 > eps_diff_bin)    isym3 = 1
              IF (fmax4 > eps_diff_bin)    isym4 = 1
              IF (fmax5 > eps_diff_bin)    isym5 = 1
            END IF

            !This is water mass rate due to diffusional growth diagnostics:
            ce_bf = 3.0_wp*col*( SUM(ff1r*(xl**2.0_wp)) )/rhocgs_loc    ! total water mass before

            flag_condevap=0
            IF (isym1==1 .AND. ((tt-tmelt)>-0.187_wp .OR. (sum(isym2)==0 .AND. isym3==0 .AND. isym4==0 .AND. isym5==0) )) THEN
              !when there is no ice, condevap_water is called
              !when there is ice, but T>0, condevap_water is called
              !since here (with the current simple melting scheme)
              !ice diffusional growth is not calculated for T>0:
              flag_condevap=1
            ELSE IF ( iceprocs==1 ) THEN
              IF ( tt > temp_no_diff) THEN !supersaturation formulas are valid (no diff growth below -60C - no humidity there)
                !when there is only ice, and T<0, condevap_ice is called:
                IF (isym1==0 .AND. (tt-tmelt)<-0.187_wp .AND. (sum(isym2)>1 .OR. isym3==1 .OR. isym4==1 .OR. isym5==1)) THEN
                  flag_condevap=2
                END IF
                !when there are both water+ice (mixed), and T<0, condevap_mixed is called:
                !when there is ice, but T>0, condevap_water is called since here ice diffusional growth is not calculated for T>0:
                IF (isym1==1 .AND. (tt-tmelt)<-0.187_wp .AND. (sum(isym2)>1 .OR. isym3==1 .OR. isym4==1 .OR. isym5==1)) THEN
                    flag_condevap=3
                END IF
                ! we do not condensate/evaporate cirrus particles below -60C, they behave as passive tracers
              END IF
            END IF
            ! Note: the formula of Rogers&Yau is not valid for -0.187<t<0, but we ignore this problem here

            supsat_out=0.0_wp
            IF ( flag_condevap == 1 ) THEN ! ... only warm phase - diffusional growth

              CALL condevap_water(tt,qq,pp,vr1_z(:) &
                              ,del1in,del2in,div1,div2 &
                              ,ff1r,ff1in &
                              ,dtcond,isym1 &
                              ,isym2,isym3,isym4,isym5,i,k &
                              ,supsat_out,rho_phy_loc)

            ELSE IF ( flag_condevap == 2 ) THEN ! ... only ice phase - nucleation + diffusional growth

              CALL condevap_ice(tt,qq,pp,vr2_z(:,:),vr3_z(:),vr4_z(:),vr5_z(:) &
                              ,del1in,del2in,div1,div2 &
                              ,ff2r,ff2in &
                              ,ff3r,ff3in &
                              ,ff4r,ff4in &
                              ,ff5r,ff5in &
                              ,dtcond,isym1,isym2,isym3,isym4,isym5,i,k &
                              ,supsat_out,rho_phy_loc)

            ELSE IF ( flag_condevap == 3 ) THEN ! ... mixed phase - nucleation + diffusional growth

              CALL condevap_mixed(tt,qq,pp,vr1_z(:),vr2_z(:,:),vr3_z(:),vr4_z(:),vr5_z(:) &
                              ,del1in,del2in,div1,div2 &
                              ,ff1r,ff1in &
                              ,ff2r,ff2in &
                              ,ff3r,ff3in &
                              ,ff4r,ff4in &
                              ,ff5r,ff5in,dtcond &
                              ,isym1,isym2,isym3,isym4,isym5,i,k &
                              ,supsat_out,rho_phy_loc)
            END IF

            !This is water mass rate due to diffusional growth diagnostics:
            ce_af = 3.0_wp*col*( sum(ff1r*(xl**2.0_wp)) )/rhocgs_loc
            del_ce_sum = del_ce_sum + (ce_af - ce_bf)

          END IF ! microphys_flag

          ! After nucleation and condevap subroutines there is graupel in ff4r or hail in ff5r, depending on hail_opt=0/1

          ! +----------------------------------+
          ! collision-coallescnce
          ! +----------------------------------+
          IF (SUM(ff1r)+SUM(ff3r)+SUM(ff4r)+SUM(ff5r) > eps_diff) THEN
            IF ( tt >= temp_stick_lim ) THEN !ice-ice collisions do not cause coagulation (do not stick)
              CALL coallescence (ff1r,ff2r,ff3r,ff4r,ff5r,tt,qq,pp,rhocgs_loc,dtcoll,i,k)
            END IF
          END IF

          t_new(i,k) = tt
          qv(i,k) = qq !kg/kg !updated only in case microphysics was active

        END DO ! end of ncondcoll loop

        ! +-------------------------------- +
        ! immediate freezing
        ! +---------------------------------+
        IF ( t_new(i,k) < tmelt .AND. iceprocs==1 ) THEN
          CALL freezing (ff1r,ff2r,ff3r,ff4r,ff5r,t_new(i,k),dt,rhocgs_loc)
          ! after nucleation and diffusional growth that were separate for each
          ! of the ice crystal types, to save computational time, the 3 ice types are
          ! added to the snow PSD. Similarly there is an option to add hail to graupel psd below
          ! (Note that freezing subroutine vreates hail but not graupel)
          DO kr = 1,nkr
            DO ice = 1,icemax
              ff3r(kr) = ff3r(kr) + ff2r(kr,ice)
              ff2r(kr,ice) = 0.0_wp
            END DO
            IF ( hail_opt == 0 ) THEN
              ff4r(kr) = ff4r(kr) + ff5r(kr)
              ff5r(kr) = 0.0_wp
            END IF
          END DO
        END IF

        ! --------------------------------------------------------------+
        ! jiwen fan melting (simplified melting with melting rates depending on hydometeor type and size)
        ! --------------------------------------------------------------+
        IF (melt_on == 1 .AND. t_new(i,k) > tmelt .AND. iceprocs == 1 ) THEN
          CALL melting(ff1r,ff2r,ff3r,ff4r,ff5r,t_new(i,k),dt,rhocgs_loc)
        END IF

        ! +---------------------------+
        ! spontanaous rain breakup
        ! +---------------------------+
        IF ( breakup_rain_spont_on == 1 .AND. ( sum(ff1r) > nkr*1.0e-15_wp ) ) THEN
          ff1r_d(:) = ff1r(:)
          CALL breakup_rain_spont (dt ,ff1r_d)
          ff1r(:) = ff1r_d(:)
        END IF

        ! ----------------------------+
        ! ... snow breakup
        ! ----------------------------+
        IF ( snow_breakup_on == 1 .AND. iceprocs == 1 .AND. sum(ff3r(kr_snow_min:nkr)) > epsil ) THEN
          CALL breakup_snow_spont(ff3r)
        END IF

        !calculate supersaturation after microphysics:
        CALL satcalc_water(t_new(i,k),qv(i,k),satur_diag,rho_phy_loc)
        diag_satur_am(i,k)=satur_diag
        diag_supsat_out(i,k)=supsat_out

        ! ... process rate (integrated)
        lh_rate(i,k) = lh_rate(i,k) +  lh_ce_1/dt
        ce_rate(i,k) = ce_rate(i,k) +  del_ce_sum/dt
        cldnucl_rate(i,k) = cldnucl_rate(i,k) + del_cldnucl_sum/dt

        ! update temperature at the end of mp
        th_phy(i,k) = t_new(i,k)/pi_phy(i,k)

        ! ... drops
        krr = 0
        DO kr = iqb_water_start, iqb_water_end
          krr = krr+1
          chem_new(i,k,kr) = ff1r(krr)
        END DO
        ! ... ccn
        krr = 0
        DO kr=iqb_ccn_start, iqb_ccn_end
          krr=krr+1
          chem_new(i,k,kr)=fccn(krr)
        END DO

        IF ( iceprocs==1 ) THEN
          ! ... snow
          krr = 0
          DO kr=iqb_snow_start, iqb_snow_end
            krr=krr+1
            chem_new(i,k,kr)=ff3r(krr)
          END DO
          ! ... hail/ graupel
          IF (hail_opt == 1)THEN
            krr = 0
            DO kr=iqb_graupel_start, iqb_graupel_end
              krr=krr+1
              chem_new(i,k,kr) = ff5r(krr)
            END DO
          ELSE
            krr = 0
            DO kr=iqb_graupel_start, iqb_graupel_end
              krr=krr+1
              chem_new(i,k,kr) = ff4r(krr)
            END DO
          END IF
        END IF
      END DO
    END DO

    ! +-----------------------------+
    ! hydrometeor sedimentation
    ! +-----------------------------+
    DO i = its,ite

    ! ... drops ...
      DO k = kts,kte
        rhocgs_z(k)=rhocgs(i,k)
        pcgs_z(k)=pcgs(i,k)
        zcgs_z(k)=zcgs(i,k)
        vrx(k,:)=vr1_z3d(:,i,k)
        krr = 0
        DO kr=iqb_water_start, iqb_water_end
          krr=krr+1
          ffx_z(k,krr)=chem_new(i,k,kr)/rhocgs(i,k)
        END DO
      END DO
      CALL sedimentation(ffx_z,vrx,rhocgs_z,zcgs_z,dt,kts,kte)
      DO k = kts,kte
        krr = 0
        DO kr=iqb_water_start, iqb_water_end
          krr=krr+1
          chem_new(i,k,kr)=ffx_z(k,krr)*rhocgs(i,k)
        END DO
      END DO

      IF ( iceprocs==1 ) THEN
        ! ... snow ...
        DO k = kts,kte
          rhocgs_z(k)=rhocgs(i,k)
          pcgs_z(k)=pcgs(i,k)
          zcgs_z(k)=zcgs(i,k)
          vrx(k,:)=vr3_z3d(:,i,k)
          krr=0
          DO kr=iqb_snow_start, iqb_snow_end
            krr=krr+1
            ffx_z(k,krr)=chem_new(i,k,kr)/rhocgs(i,k)
          END DO
        END DO

        CALL sedimentation(ffx_z,vrx,rhocgs_z,zcgs_z,dt,kts,kte)
        DO k = kts,kte
          krr=0
          DO kr=iqb_snow_start, iqb_snow_end
            krr=krr+1
            chem_new(i,k,kr)=ffx_z(k,krr)*rhocgs(i,k)
          END DO
        END DO

        ! ... hail or graupel ...
        DO k = kts,kte
          rhocgs_z(k)=rhocgs(i,k)
          pcgs_z(k)=pcgs(i,k)
          zcgs_z(k)=zcgs(i,k)
          IF (hail_opt == 1) THEN
            vrx(k,:) = vr5_z3d(:,i,k)
          ELSE
            vrx(k,:) = vr4_z3d(:,i,k)
          END IF
          krr=0
          DO kr=iqb_graupel_start, iqb_graupel_end
            krr=krr+1
            ffx_z(k,krr)=chem_new(i,k,kr)/rhocgs(i,k)
          END DO
        END DO

        CALL sedimentation(ffx_z,vrx,rhocgs_z,zcgs_z,dt,kts,kte)
        DO k = kts,kte
          krr=0
          DO kr=iqb_graupel_start, iqb_graupel_end
            krr=krr+1
            chem_new(i,k,kr)=ffx_z(k,krr)*rhocgs(i,k)
          END DO
        END DO
      END IF
    END DO

    ! ... output block
    DO k = kts,kte
      DO i = its,ite
        qc(i,k) = 0.0_wp
        qr(i,k) = 0.0_wp
        qnc(i,k) = 0.0_wp
        qnr(i,k) = 0.0_wp
        IF ( iceprocs==1 ) THEN
          qi(i,k) = 0.0_wp
          qs(i,k) = 0.0_wp
          qg(i,k) = 0.0_wp
          qni(i,k) = 0.0_wp
          qns(i,k) = 0.0_wp
          qng(i,k) = 0.0_wp
        END IF
        qna(i,k) = 0.0_wp
        reff(i,k) = 0.0_wp
        reffc(i,k) = 0.0_wp
        reffr(i,k) = 0.0_wp
        mom_2c = 0.0_wp
        mom_3c = 0.0_wp
        mom_2r = 0.0_wp
        mom_3r = 0.0_wp
        mom_2 = 0.0_wp
        mom_3 = 0.0_wp

        ! ... drop output
        krr = 0
        DO kr = iqb_water_start, iqb_water_end
          krr=krr+1
          IF (krr < krdrop)THEN
            qc(i,k) = qc(i,k) &
            + (1.0_wp/rhocgs(i,k))*col*chem_new(i,k,kr)*xl(krr)*xl(krr)*3.0_wp  ! [qc]=kg/kg, [chem_new1-33]=#/(gr*cm^3)
            qnc(i,k) = qnc(i,k) &
            + col*chem_new(i,k,kr)*xl(krr)*3.0/rhocgs(i,k)*1000.0_wp        ! [qnc]=#/kg, [chem_new1-33]=#/(gr*cm^3)

            mom_2c=mom_2c+3.0_wp*xl(krr)*chem_new(i,k,kr)*(3.0_wp*xl(krr)/(4.0_wp*pi))**(2.0_wp/3.0_wp)
            mom_3c=mom_3c+3.0_wp*xl(krr)*chem_new(i,k,kr)*(3.0_wp*xl(krr)/(4.0_wp*pi))
          ELSE
            qr(i,k) = qr(i,k) &
            + (1.0_wp/rhocgs(i,k))*col*chem_new(i,k,kr)*xl(krr)*xl(krr)*3.0_wp  ! [qr]=kg/kg, [chem_new1-33]=#/(gr*cm^3)
            qnr(i,k) = qnr(i,k) &
            + col*chem_new(i,k,kr)*xl(krr)*3.0_wp/rhocgs(i,k)*1000.0_wp        ! [qnr]=#/kg, [chem_new1-33]=#/(gr*cm^3)

            mom_2r=mom_2r+3.0_wp*xl(krr)*chem_new(i,k,kr)*(3.0_wp*xl(krr)/(4.0_wp*pi))**(2.0_wp/3.0_wp)
            mom_3r=mom_3r+3.0_wp*xl(krr)*chem_new(i,k,kr)*(3.0_wp*xl(krr)/(4.0_wp*pi))
          END IF
          mom_2=mom_2+3.0_wp*xl(krr)*chem_new(i,k,kr)*(3.0_wp*xl(krr)/(4.0_wp*pi))**(2.0_wp/3.0_wp)
          mom_3=mom_3+3.0_wp*xl(krr)*chem_new(i,k,kr)*(3.0_wp*xl(krr)/(4.0_wp*pi))
        END DO

        IF (qc(i,k) > 1.0e-6_wp .AND. mom_2c > 0.0_wp) THEN
          reffc(i,k)=(mom_3c/mom_2c)*1.0e4_wp
        END IF
        IF (qr(i,k) > 1.0e-6_wp .AND. mom_2r > 0.0_wp) THEN
          reffr(i,k)=(mom_3r/mom_2r)*1.0e4_wp
        END IF
        IF (qc(i,k) + qr(i,k) > 1.0e-6_wp .AND. mom_2 > 0.0_wp) THEN
          reff(i,k)=(mom_3/mom_2)*1.0e4_wp
        END IF

        krr=0
        IF ( iceprocs==1 ) THEN
          ! ... snow output
          krr=0
          DO kr=iqb_snow_start, iqb_snow_end
            krr=krr+1
            IF (krr <= krice)THEN
              qi(i,k) = qi(i,k) &
                      + (1.0_wp/rhocgs(i,k))*col*chem_new(i,k,kr)*xs(krr)*xs(krr)*3.0_wp
              qni(i,k)= qni(i,k) &
                      + col*chem_new(i,k,kr)*xs(krr)*3.0/rhocgs(i,k)*1000.0_wp ! #/kg
            ELSE
              qs(i,k) = qs(i,k) &
                      + (1.0_wp/rhocgs(i,k))*col*chem_new(i,k,kr)*xs(krr)*xs(krr)*3.0_wp
              qns(i,k)= qns(i,k) &
                      + col*chem_new(i,k,kr)*xs(krr)*3.0/rhocgs(i,k)*1000.0_wp ! #/kg
            END IF
          END DO
          ! ... hail / graupel output
          krr=0
          DO kr=iqb_graupel_start, iqb_graupel_end
            krr=krr+1
            ! ... hail or graupel: the output used in ICON is named qg, even IF we put hail inside
            IF (hail_opt == 1)THEN
              qg(i,k)=qg(i,k) &
              +(1.0_wp/rhocgs(i,k))*col*chem_new(i,k,kr)*xh(krr)*xh(krr)*3.0_wp
              qng(i,k)=qng(i,k) &
              +col*chem_new(i,k,kr)*xh(krr)*3.0_wp/rhocgs(i,k)*1000.0_wp ! #/kg
            ELSE
              qg(i,k)=qg(i,k) &
              +(1.0_wp/rhocgs(i,k))*col*chem_new(i,k,kr)*xg(krr)*xg(krr)*3.0_wp
              qng(i,k)=qng(i,k) &
              +col*chem_new(i,k,kr)*xg(krr)*3.0_wp/rhocgs(i,k)*1000.0_wp ! #/kg
            END IF
          END DO
        END IF

        ! ... aerosols output
        krr = 0
        DO kr = iqb_ccn_start, iqb_ccn_end
          krr = krr + 1
          qna(i,k) = qna(i,k) + col*chem_new(i,k,kr)/rhocgs(i,k)*1000.0_wp              ! [qna]=#/kg, [chem_new34-66]=#/(cm^3)
        END DO
      END DO
    END DO

    DO i = its,ite
      prec_r_sbm(i) = 0.0_wp
      IF ( iceprocs==1 ) THEN
        prec_s_sbm(i) = 0.0_wp
        prec_g_sbm(i) = 0.0_wp
      END IF
      krr=0
      DO kr=iqb_water_start, iqb_water_end
        krr=krr+1
        prec_r_sbm(i) = prec_r_sbm(i)+10.0_wp*(3.0_wp/ro1bl(krr))*col*vr1_z3d(krr,i,kte)* &
                        chem_new(i,kte,kr)*xl(krr)*xl(krr) !lowest level
      END DO

      IF ( iceprocs==1 ) THEN
        krr=0
        DO kr=iqb_snow_start, iqb_snow_end
          krr=krr+1
          prec_s_sbm(i) = prec_s_sbm(i)+10.0_wp*(3.0_wp/ro1bl(krr))*col*vr3_z3d(krr,i,kte)* & !ro1bl since we look at mm of water eqvivalent
                          chem_new(i,kte,kr)*xs(krr)*xs(krr)
        END DO
        IF (hail_opt == 1) THEN
          krr=0
          DO kr=iqb_graupel_start, iqb_graupel_end
            krr=krr+1
            prec_g_sbm(i) = prec_g_sbm(i)+10.0_wp*(3.0_wp/ro1bl(krr))*col*vr5_z3d(krr,i,kte)* &
                            chem_new(i,kte,kr)*xh(krr)*xh(krr)
          END DO
        ELSE
          krr=0
          DO kr=iqb_graupel_start, iqb_graupel_end
            krr=krr+1
            prec_g_sbm(i) = prec_g_sbm(i)+10.0_wp*(3.0_wp/ro1bl(krr))*col*vr4_z3d(krr,i,kte)* &
                            chem_new(i,kte,kr)*xg(krr)*xg(krr)
          END DO
        END IF
      END IF
    END DO

    IF (convert_micro2advect) THEN
      ! ... drops
      krr=0
      DO kr=iqb_water_start, iqb_water_end
        krr=krr+1
        DO k = kts,kte
          DO i = its,ite
            chem_new(i,k,kr)=chem_new(i,k,kr)/rhocgs(i,k)*col*xl(krr)*xl(krr)*3.0_wp
            IF (qc(i,k)+qr(i,k) < epsil) chem_new(i,k,kr)=0.0_wp
          END DO
        END DO
      END DO
      ! ... snow
      IF ( iceprocs==1 ) THEN
        krr=0
        DO kr=iqb_snow_start, iqb_snow_end
          krr=krr+1
          DO k = kts,kte
            DO i = its,ite
              chem_new(i,k,kr)=chem_new(i,k,kr)/rhocgs(i,k)*col*xs(krr)*xs(krr)*3.0_wp
              IF (qs(i,k) <  epsil) chem_new(i,k,kr)=0.0_wp
            END DO
          END DO
        END DO
      END IF
      ! ... ccn
      krr=0
      DO kr=iqb_ccn_start, iqb_ccn_end
        krr=krr+1
        DO k = kts,kte
          DO i = its,ite
            chem_new(i,k,kr)=chem_new(i,k,kr)/rhocgs(i,k)*1000.0_wp
          END DO
        END DO
      END DO
      ! ... hail / graupel
      IF ( iceprocs==1 ) THEN
        IF (hail_opt == 1) THEN
          krr=0
          DO kr=iqb_graupel_start, iqb_graupel_end
            krr=krr+1
            DO k = kts,kte
              DO i = its,ite
                chem_new(i,k,kr)=chem_new(i,k,kr)/rhocgs(i,k)*col*xh(krr)*xh(krr)*3.0_wp
                IF (qg(i,k) < epsil) chem_new(i,k,kr) = 0.0_wp
              END DO
            END DO
          END DO
        ELSE
          krr=0
          DO kr=iqb_graupel_start, iqb_graupel_end
            krr=krr+1
            DO k = kts,kte
              DO i = its,ite
                chem_new(i,k,kr)=chem_new(i,k,kr)/rhocgs(i,k)*col*xg(krr)*xg(krr)*3.0_wp
                IF (qg(i,k) < epsil) chem_new(i,k,kr) = 0.0_wp
              END DO
            END DO
          END DO
        END IF
      END IF
    END IF

    RETURN
  END SUBROUTINE fast_sbm

  SUBROUTINE total_water(qv,ff1r,ff3r,ff4r,ff5r,tot_water,rhocgs)
    IMPLICIT NONE
    REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: ff1r,ff3r,ff4r,ff5r
    REAL(KIND=wp), INTENT(INOUT)            :: tot_water
    REAL(KIND=wp), INTENT(IN)               :: qv,rhocgs
    INTEGER :: kr, krr

    tot_water=qv

    ! ... drop output
    krr = 0
    DO kr = iqb_water_start, iqb_water_end
      krr=krr+1
      tot_water = tot_water+(1.0_wp/rhocgs)*col*ff1r(krr)*xl(krr)*xl(krr)*3.0_wp  ! [qc]=kg/kg, [chem_new1-33]=#/(gr*cm^3)
    END DO
    IF ( iceprocs==1 ) THEN
      ! ... snow output
      krr=0
      DO kr=iqb_snow_start, iqb_snow_end
        krr=krr+1
        tot_water = tot_water+(1.0_wp/rhocgs)*col*ff3r(krr)*xs(krr)*xs(krr)*3.0_wp
      END DO
      ! ... hail / graupel output
      krr=0
      DO kr=iqb_graupel_start, iqb_graupel_end
        krr=krr+1
        ! ... hail or graupel: the output used in ICON is named qg, even IF we put hail inside
        IF (hail_opt == 1) THEN
          tot_water = tot_water+(1.0_wp/rhocgs)*col*ff5r(krr)*xh(krr)*xh(krr)*3.0_wp
        ELSE
          tot_water = tot_water+(1.0_wp/rhocgs)*col*ff4r(krr)*xg(krr)*xg(krr)*3.0_wp
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE total_water

  SUBROUTINE satcalc_water(t,q,satur_diag,rho) ! saturation calculation
    IMPLICIT NONE
    REAL(KIND=wp), INTENT(IN)   :: t, rho, q   !rho [kg/m^3], q [kg/kg]
    REAL(KIND=wp), INTENT(INOUT):: satur_diag
    REAL(KIND=wp) :: es1n, ew1n, qv

    qv=q
    IF (qv <= 0.0_wp) qv = 1.e-10_wp
    es1n=10.0_wp*sat_pres_water(t) !water vapor pressure at saturation over water [dynes/cm^2] (1Pa=10dynes/cm^2)
    ew1n=10.0_wp*qv*rho*rv*t       !water vapor pressure over water [dynes/cm^2] (1Pa=10dynes/cm^2) !rv [J/K/kg]=[m^2/K/s^2]
    satur_diag=ew1n/es1n-1.0_wp

    RETURN
  END SUBROUTINE satcalc_water

  SUBROUTINE satcalc_ice(t,q,satur_diag,rho) ! saturation calculation
    IMPLICIT NONE
    REAL(KIND=wp), INTENT(IN)   :: t, rho, q !rho [kg/m^3], q [kg/kg]
    REAL(KIND=wp), INTENT(INOUT):: satur_diag
    REAL(KIND=wp) :: es2n, ew1n, qv

    qv=q
    IF (qv <= 0.0_wp) qv = 1.e-10_wp
    es2n=10.0_wp*sat_pres_ice(t)
    ew1n=10.0_wp*qv*rho*rv*t
    satur_diag=ew1n/es2n-1.0_wp

    RETURN
  END SUBROUTINE satcalc_ice

  SUBROUTINE sedimentation(chem_new,vr1,rhocgs,zcgs,dt,kts,kte)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(INOUT) :: chem_new(:,:)
    REAL(KIND=wp),INTENT(IN) :: vr1(:,:),rhocgs(:),zcgs(:),dt
    INTEGER,INTENT(IN) :: kts,kte
    INTEGER :: k,kr
    REAL(KIND=wp) :: tfall,dtfall,vfall(kte),dwflux(kte)
    INTEGER :: ifall,n,nsub
    ! falling fluxes for each kind of cloud particles: c.g.s. unit
    ! adapted from gsfc code for hucm
    ! the flux at k=1 is assumed to be the ground so flux(1) is the
    ! flux into the ground. dwflux(1) is at the lowest half level where
    ! q(1) etc are defined. the formula for flux(1) uses q(1) etc which
    ! is actually half a grid level above it. this is what is meant by
    ! an upstream method. upstream in this case is above because the
    ! velocity is downwards.
    ! use upstream method (vfall is positive)

    DO kr=1,nkr
      ifall=0
      DO k = kts,kte
      IF (chem_new(k,kr).GE.epsil) ifall=1   ! IF there is some mass to fall
      END DO
      IF (ifall.EQ.1) THEN
        tfall=1.e10_wp
        DO k=kts,kte
          ! [ks] vfall(k) = vr1(k,kr)*SQRT(1.e6/pcgs(k))
          vfall(k) = vr1(k,kr) ! ... [ks] : the pressure effect is taken into account at the beggining of the calculations
          tfall=MIN(tfall,zcgs(k)/(vfall(k)+epsil))   ! set tfall=z/v
        END DO
        IF (tfall.GE.1.e10_wp) stop
        nsub=(INT(2.0_wp*dt/tfall)+1)
        dtfall=dt/nsub         ! dt used in sedimintation, which is smaller than model dt to obey CFL criterion

        DO n=1,nsub            ! loop over sedimintation sub steps
          DO k=kte,kts+1,-1    !loop from 65 up to 2, each level -(lower-upper)/upper. Used to be DO k=kts,kte-1
            dwflux(k)=-(rhocgs(k)*vfall(k)*chem_new(k,kr)-&         ! flux difference (divergence) at level k
                        rhocgs(k-1)*vfall(k-1)*chem_new(k-1,kr))/ &
                        (rhocgs(k)*(zcgs(k-1)-zcgs(k)))
          END DO
          ! no z above top, so use the same deltaz
          dwflux(kts)=-(rhocgs(kts)*vfall(kts)*chem_new(kts,kr))/ & ! flux difference (divergence) at top
                       (rhocgs(kts)*(zcgs(kts)-zcgs(kts+1)))
          DO k=kts,kte
            chem_new(k,kr)=chem_new(k,kr)+dwflux(k)*dtfall !if the droplets fall into some level, add them to the psd there
          END DO
        END DO
      END IF
    END DO

    RETURN
  END SUBROUTINE sedimentation

  SUBROUTINE condevap_water (tt,qq,pp,vr1,del1n,del2n,div1,div2,ff1,psi1,dtcond, &
                             isym1,isym2,isym3,isym4,isym5,iin,kin, supsat_out,rho)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(INOUT) :: tt, qq, del1n, del2n, div1, div2, ff1(nkr), supsat_out, psi1(nkr)
    REAL(KIND=wp),INTENT(IN) :: pp, vr1(nkr), rho, dtcond
    INTEGER, INTENT(INOUT) :: isym1,isym2(icemax),isym3,isym4,isym5
    INTEGER, INTENT(IN) :: iin, kin
    INTEGER k,kr,itime,kcond,idrop
    REAL(KIND=wp) :: ror,gam,dt,dtt,b6,del1,del2,del1s,del2s, &
     & timenew,timerev,sfn11,sfn12,sfnl,sfni,b5l,b5i,b7l,b7i,dopl,dopi,rw,ri,qw,pw, &
     & pi,dtnewl,d1n,ff1_old(nkr),supintw(nkr),dsupintw(nkr),tpn,tps,qpn,qps,told,qold, &
     & rmasslaa,rmasslbb,es1n,es2n,ew1n,oper2,ar1,delmassl1, &
     & fi1(nkr),b11_my(nkr),fl1(nkr),sfndummy(3), &
     & del1_d,del2_d,rw_d,pw_d,ri_d,pi_d,d1n_d,d2n_d
    DATA gam /1.e-4_wp/

    oper2(ar1)=0.622_wp/(0.622_wp+0.378_wp*ar1)/ar1

    ror=rho*0.001_wp !total density in cgs

    sfndummy = 0.0_wp
    b11_my = 0.0_wp
    kcond=0

    dt=dtcond
    dtt=dtcond
    b6=0.0_wp

    DO kr=1,nkr
      ff1_old(kr)=ff1(kr)
      supintw(kr)=0.0_wp
      dsupintw(kr)=0.0_wp
    END DO

    tpn=tt
    qpn=qq
    DO kr=1,nkr
      fi1(kr)=ff1(kr)
    END DO

    ! warm mp (condensation or evaporation) (begin)
    timenew=0.0_wp
    itime=0
    told = tpn
    qold = qpn
    sfnl = 0.0_wp
    sfn11 = 0.0_wp

56  itime = itime+1
    timerev = dt-timenew
    del1 = del1n
    del2 = del2n
    del1s = del1n
    del2s = del2n
    tps = tpn
    qps = qpn

    IF (isym1 == 1) THEN !water droplets exist
      fl1 = 0.0_wp
      ! calculation of f_1 in (3.10) khain&sednev,1996:
      CALL condevap_dmdt_coef(xl,tps,pp,vr1,rlec,ro1bl,b11_my,1,1,fl1)
      sfndummy(1)=sfn11
      ! calculation of coefficients in equations for supersaturation
      CALL condevap_dsdt_coef(fi1,xl,sfndummy,b11_my,1./ror,1)
      sfn11 = sfndummy(1)
    END IF

    sfn12 = 0.0_wp
    sfnl = sfn11 + sfn12
    sfni = 0.0_wp
    b5l=bb1_my/tps/tps
    b5i=bb2_my/tps/tps
    b7l=b5l*b6
    b7i=b5i*b6
    dopl=1.0_wp+del1s
    dopi=1.0_wp+del2s

    ! after each substep, the new sup sat over water and ice are calculated. These values are used to
    ! calculate the new temperature and humidity at each substep. Therefore sup sat over ice is calculated even
    ! when there is no ice. The procedure can in future be simplified by assuming a linear changes of T and Q

    IF (.not. latheatfac1) THEN
      rw=(oper2(qps)+b5l*al1)*dopl*sfnl !R1 coeff in (3.11) in Khain&Sednev,1996, al1=Lw/Cp
      ri=(oper2(qps)+b5l*al2)*dopl*sfni !R2 coeff in (3.11) in Khain&Sednev,1996, al2=Li/Cp
    ELSE
      rw=(oper2(qps)+b5l*al1*cpd/cvd)*dopl*sfnl !R1 coeff in (3.11) in Khain&Sednev,1996, al1=Lw/Cp
      ri=(oper2(qps)+b5l*al2*cpd/cvd)*dopl*sfni !R2 coeff in (3.11) in Khain&Sednev,1996, al2=Li/Cp
    END IF
    IF (.not. latheatfac1) THEN
      pw=(oper2(qps)+b5i*al1)*dopi*sfnl !P1 coeff in (3.11) in Khain&Sednev,1996
      pi=(oper2(qps)+b5i*al2)*dopi*sfni !P2 coeff in (3.11) in Khain&Sednev,1996
    ELSE
      pw=(oper2(qps)+b5i*al1*cpd/cvd)*dopi*sfnl !P1 coeff in (3.11) in Khain&Sednev,1996
      pi=(oper2(qps)+b5i*al2*cpd/cvd)*dopi*sfni !P2 coeff in (3.11) in Khain&Sednev,1996
    END IF
    qw=b7l*dopl
    IF ( rw .NE. rw .OR. pw .NE. pw ) THEN
      PRINT*, 'nan in condevap_water',psi1
      CALL finish(TRIM(modname),"fatal error in condevap_water (rw or pw are nan), model stop")
    END IF

    ! every ncond substep can be still too large for droplets growth. we want that during diffusional growth,
    ! the droplet radius change will be less THEN few bins. therefore ncond substep is further divided by up
    ! to kcond=10. in case that after 10 steps, time will be still lower than ncond sub step, the last step
    ! will have this delta time. note that this is done without updates of t,q,s
    kcond=10 ! kcond is a flag
    IF (del1n >= 0.0_wp) kcond=11 ! is it diffusional growth or evaporation --> kcond
    IF (kcond == 11) THEN
      dtnewl = dt
      dtnewl = MIN(dtnewl,timerev) ! timerev are small substeps within ncond sub step
      timenew = timenew + dtnewl
      dtt = dtnewl

      IF (dtt < 0.0_wp) CALL finish(TRIM(modname),"fatal error in condevap_water-del1n>0:(dtt<0), model stop")
      del1_d = del1
      del2_d = del2
      rw_d = rw
      pw_d = pw
      ri_d = ri
      pi_d = pi

      ! solving the equation for 2 supersaturations del1n,del2n
      CALL condevap_supsat_eqn(del1_d,del2_d,del1n,del2n, &
                  rw_d,pw_d,ri_d,pi_d, &
                  dtt,d1n_d,d2n_d,0.0_wp,0.0_wp, &
                  isym1,isym2,isym3,isym4,isym5)
      del1 = del1_d
      del2 = del2_d
      rw = rw_d
      pw = pw_d
      ri = ri_d
      pi = pi_d
      d1n = d1n_d

      IF (isym1 == 1) THEN ! water exists
        idrop = isym1
        !bin mass change and remapping - (3.14) in khain&sednev,1996:
        CALL condevap_mass_eqn_and_remap1(xl, b11_my, fi1, psi1, d1n, &
                    isym1, 1, 1, idrop, iin, kin)
      END IF

      IF ((del1 .GT. 0 .AND. del1n .LT. 0) .AND. ABS(del1n) .GT. epsdel) THEN
        CALL finish(TRIM(modname),"fatal error in condevap_water-1 (del1.GT.0.AND.del1n.LT.0), model stop")
      END IF

    ELSE
      ! evaporation - only water
      ! in case : kcond.NE.11
      dtnewl = dt
      dtnewl = MIN(dtnewl,timerev)
      timenew = timenew + dtnewl
      dtt = dtnewl
      IF (dtt < 0.0_wp) CALL finish(TRIM(modname),"fatal error in condevap_water-del1n<0:(dtt<0), model stop")

      del1_d = del1
      del2_d = del2
      rw_d = rw
      pw_d = pw
      ri_d = ri
      pi_d = pi

      ! solving the equation for 2 supersaturations del1n,del2n. Check, why division of the code for diff
      ! growth and evap is needed (???ok):
      CALL condevap_supsat_eqn(del1_d,del2_d,del1n,del2n, &
                 rw_d,pw_d,ri_d,pi_d, dtt,d1n_d,d2n_d,0.0_wp,0.0_wp, &
                 isym1,isym2,isym3,isym4,isym5)
      del1 = del1_d
      del2 = del2_d
      rw = rw_d
      pw = pw_d
      ri = ri_d
      pi = pi_d
      d1n = d1n_d

      IF (isym1 == 1) THEN
        idrop = isym1
        !bin mass change and remapping - (3.14) in khain&sednev,1996:
        CALL condevap_mass_eqn_and_remap1(xl, b11_my,fi1, psi1, d1n, &
                 isym1, 1, 1, idrop, iin, kin)
      END IF

      IF ((del1 .LT. 0 .AND. del1n .GT. 0) .AND. ABS(del1n) .GT. epsdel) THEN
        CALL finish(TRIM(modname),"fatal error in condevap_water-2 (del1.LT.0.AND.del1n.GT.0), model stop")
      END IF
    END IF

    ! masses:
    rmasslbb=0.
    rmasslaa=0.

    ! ... before jernewf (only water)
    DO k=1,nkr
      rmasslbb = rmasslbb+fi1(k)*xl(k)*xl(k) !calculate cloud water content before diffusional
                                             !growth ncond substep, xl(k) !mass of bin k, named xl outside
    END DO
    rmasslbb = rmasslbb*col*3.0_wp/ror !result: cloud water content before diffusional growth ncond substep
    IF (rmasslbb .LE. 0.0_wp) rmasslbb=0.0_wp
    ! ... after jernewf (only water)

    DO k=1,nkr
      rmasslaa=rmasslaa+psi1(k)*xl(k)*xl(k)
    END DO
    rmasslaa=rmasslaa*col*3.0_wp/ror !result: cloud water content after diffusional growth ncond substep

    IF (rmasslaa.LE.0.0_wp) rmasslaa=0.0_wp

    delmassl1 = rmasslaa - rmasslbb !cloud mass change during substep (new-old)
    qpn = qps - delmassl1 !new specific humidity

    tpn = tps + al1*delmassl1*cpd/cvd ! new temperature after substep (latent heat release or evaporation)

    IF (ABS(al1*delmassl1) > 10.0_wp )THEN !temperature change during one substep (latent heat release or evaporation)
      WRITE(txt,'(a,i3,a,i3,a,i3,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5, &
&                      a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5,a,e12.5)') &
&                      ' i=',iin,' k=',kin,' delmassl1=',delmassl1,' del1n=',del1n,' del2n=',del2n, &
&                      ' del1=',del1,' del2=',del2,' d1n=',d1n,' rw=',rw,' pw=',pw,' ri=',ri,' pi=',pi, &
&                      ' dt=',dt,' tps=',tps,' qps=',qps,' before min(fi1)=',MINVAL(fi1), &
&                      ' before max(fi1)=',MAXVAL(fi1),' before min(psi1)=',MINVAL(psi1), &
&                      ' before max(psi1)=',MAXVAL(psi1)
      CALL message('', TRIM(txt))
      CALL finish(TRIM(txt),"fatal error in condevap_water-in (ABS(al1*delmassl1) > 10.0_wp), model stop")
    END IF

    ! ... supersaturation (only water)
    es1n=10.0_wp*sat_pres_water(tpn) !10.0_wp because Pa-->dynes/cm^2
    es2n=10.0_wp*sat_pres_ice(tpn)   !dynes/cm^2
    ew1n=10.0_wp*qpn*rho*rv*tpn      !dynes/cm^2

    IF (es1n == 0.0_wp) THEN
      del1n=0.5_wp
      div1=1.5_wp
    ELSE
      div1 = ew1n/es1n
      del1n = ew1n/es1n-1.0_wp
    END IF
    IF (es2n == 0.0_wp)THEN
      del2n=0.5_wp
      div2=1.5_wp
    ELSE
      del2n = ew1n/es2n-1.0_wp
      div2 = ew1n/es2n
    END IF

    ! calculation of the full integral over the sup sat (dsupintw), without ncond substeps, and
    ! THEN perform a single remapping, instead of doing remapping after each substep
    IF (isym1 == 1) THEN
      DO kr=1,nkr
        supintw(kr)=supintw(kr)+b11_my(kr)*d1n
        dsupintw(kr)=dsupintw(kr)+b11_my(kr)*d1n
      END DO
    END IF

    ! ... repeate time step (only water: condensation or evaporation)
    IF (timenew.LT.dt) GOTO 56

    IF (isym1 == 1) THEN
      ! a single diffusional growth step using the integral of sup sat. We do it for water only,
      ! since exact collision initiation time is important, to prevent water PSD broadening
      ! and early rain initiation. Bin mass change and remapping - (3.14) in khain&sednev,1996:
      ! The difference between condevap_mass_eqn_and_remap1 and condevap_mass_eqn_and_remap2 is:
      ! In condevap_mass_eqn_and_remap1 (3.14) is summed up within the ncond loop, and the
      ! remapping is done ncond times, but only for t update and not for the update of the final PSD.
      ! In condevap_mass_eqn_and_remap2 (only for water) (3.14) is summed up at the end, after ncond
      ! loop, so that the remapping which updates the final PSD is performed only once:
      CALL condevap_mass_eqn_and_remap2 (xl,supintw, ff1_old,psi1, idrop, iin,kin)
    END IF

    rmasslaa=0.0_wp
    rmasslbb=0.0_wp

    DO k=1,nkr
      rmasslbb=rmasslbb+ff1_old(k)*xl(k)*xl(k)
    END DO
    rmasslbb=rmasslbb*col*3.0_wp/ror
    IF (rmasslbb.LT.0.0_wp) rmasslbb=0.0_wp

    DO k=1,nkr
      rmasslaa=rmasslaa+psi1(k)*xl(k)*xl(k)
    END DO
    rmasslaa=rmasslaa*col*3.0_wp/ror
    IF (rmasslaa.LT.0.0_wp) rmasslaa=0.0_wp
    delmassl1 = rmasslaa-rmasslbb

    !latent heat release:
    qpn = qold - delmassl1

    tpn = told + al1*delmassl1*cpd/cvd
    lh_ce_1 = lh_ce_1 + al1*delmassl1*cpd/cvd

    IF (ABS(al1*delmassl1) > 10.0_wp )THEN
      PRINT*,"condevap_water-out"," i=",iin,"kin",kin
      PRINT*,"del1n,del2n,rw,pw,ri,pi,dt"
      PRINT*,del1n,del2n,rw,pw,ri,pi,dtt
      PRINT*,"i=",iin,"kin",kin,"tps=",tps,"qps=",qps,"delmassl1",delmassl1
      PRINT*,"al1=",al1,rmasslbb,rmasslaa,"fi1",fi1,"psi1",psi1
      CALL finish(TRIM(modname),"fatal error 2 in condevap_water-out (ABS(al1*delmassl1) > 10.0_wp), model stop")
    END IF

    ! ... supersaturation
    es1n=10.0_wp*sat_pres_water(tpn) !dynes/cm^2
    es2n=10.0_wp*sat_pres_ice(tpn)   !dynes/cm^2
    ew1n=10.0_wp*qpn*rho*rv*tpn      !dynes/cm^2

    IF (es1n == 0.0_wp) THEN
      del1n=0.5_wp
      div1=1.5_wp
      CALL finish(TRIM(modname),"fatal error in condevap_water (es1n.EQ.0_wp), model stop")
    ELSE
      div1=ew1n/es1n
      del1n=ew1n/es1n-1.0_wp
    END IF
    IF (es2n .EQ. 0.0_wp)THEN
      del2n=0.5_wp
      div2=1.5_wp
      CALL finish(TRIM(modname),"fatal error in condevap_water (es2n.EQ.0_wp), model stop")
    ELSE
      del2n=ew1n/es2n-1.
      div2=ew1n/es2n
    END IF

    tt=tpn
    qq=qpn
    DO kr=1,nkr
      ff1(kr)=psi1(kr)
    END DO
    supsat_out=del1n

    RETURN
  END SUBROUTINE condevap_water

  SUBROUTINE condevap_ice &          ! deposition (like condensation) and sublimation (like evaporation)
                    & (tt,qq,pp,vr2,vr3,vr4,vr5,del1n,del2n,div1,div2 &
                    & ,ff2,psi2 &
                    & ,ff3,psi3 &
                    & ,ff4,psi4 &
                    & ,ff5,psi5 &
                    & ,dtcond,isym1,isym2,isym3,isym4,isym5 &
                    & ,iin,kin,supsat_out,rho)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(INOUT) :: tt, qq, del1n, del2n, div1, div2, &
                              & ff2(nkr,icemax), ff3(nkr), ff4(nkr), ff5(nkr), &
                              & psi2(nkr,icemax), psi3(nkr), psi4(nkr), psi5(nkr), supsat_out
    REAL(KIND=wp),INTENT(IN) :: pp, vr2(nkr,icemax), vr3(nkr), vr4(nkr), vr5(nkr), &
                              & rho, dtcond
    INTEGER, INTENT(INOUT) :: isym1,isym2(icemax),isym3,isym4,isym5
    INTEGER, INTENT(IN) :: iin, kin
    INTEGER k,kr,ice,itime,kcond,idrop
    REAL(kind=wp) :: ror, fi2(nkr,icemax), fi3(nkr), fi4(nkr), fi5(nkr), &
                   & dt,dtt,b6,del1,del2,del1s,del2s,timenew,timerev, &
                   & sfnl,sfni,sfnii1,sfn21,sfn22,sfn31,sfn41,sfn51, &
                   & b5l,b5i,b7l,b7i,dopl,dopi,operq,rw,ri,qw,pw, &
                   & pi,qi,dtnewl,d2n,tpn,tps,qpn,qps,rmassiaa,rmassibb, &
                   & es1n,es2n,ew1n,oper2,ar1,delmassi1, &
                   & b21_my(nkr,icemax),b31_my(nkr),b41_my(nkr),b51_my(nkr), &
                   & sfni1(icemax), fl1(nkr), sfndummy(3), fl3(nkr), fl4(nkr), fl5(nkr), &
                   & del1_d, del2_d, rw_d, pw_d, ri_d, pi_d, d1n_d, d2n_d

    oper2(ar1)=0.622_wp/(0.622_wp+0.378_wp*ar1)/ar1

    ror=rho*0.001_wp !total density in cgs

    b21_my = 0.0_wp
    b31_my = 0.0_wp
    b41_my = 0.0_wp
    b51_my = 0.0_wp

    sfndummy = 0.0_wp

    sfni1 = 0.0_wp
    sfn31 = 0.0_wp
    sfn41 = 0.0_wp
    sfn51 = 0.0_wp

    itime=0
    kcond=0

    dt=dtcond
    dtt=dtcond
    b6=0.0_wp

    tpn=tt
    qpn=qq

! only ice (condensation or evaporation) (begin)
    timenew = 0.0_wp
    itime = 0

46  itime = itime + 1

    timerev=dt-timenew

    del1=del1n
    del2=del2n
    del1s=del1n
    del2s=del2n
    tps=tpn
    qps=qpn
    DO kr=1,nkr
      fi3(kr)=psi3(kr)
      fi4(kr)=psi4(kr)
      fi5(kr)=psi5(kr)
      DO ice=1,icemax
        fi2(kr,ice)=psi2(kr,ice)
      END DO
    END DO

    IF (sum(isym2) > 0) THEN !ice exists
      fl1 = 0.0_wp
      ! calculation of coefficients in the eqn for ice diffusional growth (deposition/sublimation):
      CALL condevap_dmdt_coef (xi,tps,pp,vr2,riec,ro2bl,b21_my,3,2,fl1)
      ! calculation of coefficients in the eqn for sup sat over ice:
      CALL condevap_dsdt_coef (fi2,xi,sfni1,b21_my,1.0_wp/ror,icemax)
    END IF
    IF (isym3 == 1) THEN !snow
      fl3 = 0.0_wp
      CALL condevap_dmdt_coef (xs,tps,pp,vr3,rsec,ro3bl,b31_my,1,3,fl3)
      sfndummy(1) = sfn31
      CALL condevap_dsdt_coef(fi3,xs,sfndummy,b31_my,1.0_wp/ror,1)
      sfn31 = sfndummy(1)
    END IF
    IF (isym4 == 1) THEN !graupel
      fl4 = 0.0_wp
      CALL condevap_dmdt_coef(xg,tps,pp,vr4,rgec,ro4bl,b41_my,1,2,fl4)
      sfndummy(1) = sfn41
      CALL condevap_dsdt_coef(fi4,xg,sfndummy,b41_my,1./ror,1)
      sfn41 = sfndummy(1)
    END IF
    IF (isym5 == 1) THEN !hail
      fl5 = 0.0_wp
      CALL condevap_dmdt_coef(xh,tps,pp,vr5,rhec,ro5bl,b51_my,1,2,fl5)
      sfndummy(1) = sfn51
      CALL condevap_dsdt_coef(fi5,xh,sfndummy,b51_my,1.0_wp/ror,1)
      sfn51 = sfndummy(1)
    END IF

    ! note that the coefficients in condevap_dmdt_coef,condevap_dsdt_coef defer f
    ! ice/snow/graupel/hail because different shapes (capacities)

    sfnii1 = sfni1(1) + sfni1(2) + sfni1(3)
    sfn21 = sfnii1 + sfn31 + sfn41 + sfn51
    sfnl = 0.0_wp
    sfn22 = 0.0_wp
    sfni = sfn21 + sfn22

    b5l=bb1_my/tps/tps
    b5i=bb2_my/tps/tps
    b7l=b5l*b6
    b7i=b5i*b6
    dopl=1.0_wp+del1s
    dopi=1.0_wp+del2s
    operq=oper2(qps)

    ! after each substep, the new sup sat over water and ice are calculated. These values are used to
    ! calculate the new temperature and humidity at each substep. Therefore sup sat over ice is calculated even
    ! when there is no ice. The procedure can in future be simplified by assuming a linear changes of T and Q

    IF (.not. latheatfac1) THEN
      rw=(operq+b5l*al1)*dopl*sfnl !R1 coeff in (3.11) in Khain&Sednev,1996, al1=Lw/Cp
      ri=(operq+b5l*al2)*dopl*sfni !R2 coeff in (3.11) in Khain&Sednev,1996, al2=Li/Cp
    ELSE
      rw=(operq+b5l*al1*cpd/cvd)*dopl*sfnl !R1 coeff in (3.11) in Khain&Sednev,1996, al1=Lw/Cp
      ri=(operq+b5l*al2*cpd/cvd)*dopl*sfni !R2 coeff in (3.11) in Khain&Sednev,1996, al2=Li/Cp
    END IF
    IF (.not. latheatfac1) THEN
      pw=(oper2(qps)+b5i*al1)*dopi*sfnl !P1 coeff in (3.11) in Khain&Sednev,1996
      pi=(oper2(qps)+b5i*al2)*dopi*sfni !P2 coeff in (3.11) in Khain&Sednev,1996
    ELSE
      pw=(oper2(qps)+b5i*al1*cpd/cvd)*dopi*sfnl !P1 coeff in (3.11) in Khain&Sednev,1996
      pi=(oper2(qps)+b5i*al2*cpd/cvd)*dopi*sfni !P2 coeff in (3.11) in Khain&Sednev,1996
    END IF
    qw=b7l*dopl
    qi=b7i*dopi

    IF (rw.ne.rw .OR. pw.ne.pw) THEN
      PRINT*, 'nan in condevap_ice'
      CALL finish(TRIM(modname),"fatal error in condevap_ice (rw or pw are nan), model stop")
    END IF

    ! every ncond substep can be still too large for ice growth. We want that during
    ! diffusional growth, the ice size change will be less THEN few bins. Therefore
    ! ncond substep is further divided for steps not less than 0.4 sec. In case that
    ! after several steps, time will be still lower than ncond sub step, the last
    ! step will have this delta time. Note that this is done without updates of t,q,s

    kcond=20 ! kcond is a flag
    IF (del2n > 0.0_wp) kcond=21 ! is it deposition/sublimation? --> kcond

    ! ... (only ice)
    IF (kcond == 21) THEN
      ! ... only_ice: condensation
      dtnewl = dt
      dtnewl = MIN(dtnewl,timerev) ! timerev are small substeps within ncond sub step
      timenew = timenew + dtnewl
      dtt = dtnewl

      IF (dtt < 0.0_wp) CALL finish(TRIM(modname),"fatal error in condevap_ice-del2n>0:(dtt<0), model stop")

      del1_d = del1
      del2_d = del2
      rw_d = rw
      pw_d = pw
      ri_d = ri
      pi_d = pi

      ! solving the equation for 2 supersaturations del1n,del2n:
      CALL condevap_supsat_eqn(del1_d,del2_d,del1n,del2n, &
                        rw_d,pw_d,ri_d,pi_d, &
                        dtt,d1n_d,d2n_d,0.0_wp,0.0_wp, &
                        isym1,isym2,isym3,isym4,isym5)
      del1 = del1_d
      del2 = del2_d
      rw = rw_d
      pw = pw_d
      ri = ri_d
      pi = pi_d
      d2n = d2n_d

      IF (sum(isym2) > 0)THEN ! ice exists
        idrop = 0
        fl1 = 0.0_wp
        IF (isym2(1) == 1) THEN
          ! diffusional growth+remapping for crystals type 1:
          CALL condevap_mass_eqn_and_remap1(xi(:,1), b21_my(:,1), &
                     fi2(:,1), psi2(:,1), d2n, &
                     isym2(1), icemax, 1, idrop, iin, kin)
        END IF
        ! diffusional growth+remapping for crystals type 2:
        IF (isym2(2) == 1) THEN
          CALL condevap_mass_eqn_and_remap1(xi(:,2), b21_my(:,2), &
                     fi2(:,2), psi2(:,2), d2n, &
                     isym2(2), icemax, 2, idrop, iin ,kin)
        END IF
        ! diffusional growth+remapping for crystals type 3:
        IF (isym2(3) == 1) THEN
          CALL condevap_mass_eqn_and_remap1(xi(:,3), b21_my(:,3), &
                     fi2(:,3), psi2(:,3), d2n, &
                     isym2(3), icemax, 3, idrop, iin ,kin)
        END IF
      END IF

      IF (isym3 == 1) THEN
        idrop = 0
        fl3 = 0.0_wp
        ! diffusional growth+remapping for snow (why different?):
        CALL condevap_mass_eqn_and_remap1(xs, b31_my, &
                     fi3, psi3, d2n, &
                     isym3, 1, 3, idrop, iin ,kin)
      END IF

      IF (isym4 == 1) THEN
        idrop = 0
        fl4 = 0.0_wp
        ! diffusional growth+remapping for graupel:
        CALL condevap_mass_eqn_and_remap1(xg, b41_my, &
                     fi4, psi4, d2n, &
                     isym4, 1, 4, idrop, iin ,kin)
      END IF

      IF (isym5 == 1) THEN
        idrop = 0
        fl5 = 0.0_wp
        ! diffusional growth+remapping for hail:
        CALL condevap_mass_eqn_and_remap1(xh, b51_my, &
                     fi5, psi5, d2n, &
                     isym5, 1, 5, idrop, iin, kin)
      END IF

      IF ((del2.GT.0.AND.del2n.LT.0) .AND. ABS(del2n).GT.epsdel) THEN
        CALL finish(TRIM(modname),"fatal error in module_mp_fast_sbm (del2.GT.0.AND.del2n.LT.0), model stop")
      END IF

    ELSE ! it used to be not separated
      ! ... in case kcond.ne.21
      ! only ice: evaporation
      dtnewl = dt
      dtnewl = MIN(dtnewl,timerev)
      timenew = timenew + dtnewl
      dtt = dtnewl

      IF (dtt < 0.0_wp) CALL finish(TRIM(modname),"fatal error in condevap_ice-del2n<0:(dtt<0), model stop")

      del1_d = del1
      del2_d = del2
      rw_d = rw
      pw_d = pw
      ri_d = ri
      pi_d = pi

      ! solving the equation for 2 supersaturations del1n,del2n. Check, why division of the code for diff growth and evap is needed:
      CALL condevap_supsat_eqn(del1_d,del2_d,del1n,del2n, &
                                 rw_d,pw_d,ri_d,pi_d, &
                                 dtt,d1n_d,d2n_d,0.0_wp,0.0_wp, &
                                 isym1,isym2,isym3,isym4,isym5)
      del1 = del1_d
      del2 = del2_d
      rw = rw_d
      pw = pw_d
      ri = ri_d
      pi = pi_d
      d2n = d2n_d

      IF (sum(isym2) > 0) THEN
        idrop = 0
        fl1 = 0.0_wp
        ! ?? bin mass change and remapping - (3.14) in Khain&Sednev,1996:
        IF (isym2(1)==1)THEN
          CALL condevap_mass_eqn_and_remap1(xi(:,1), b21_my(:,1), &
                       fi2(:,1), psi2(:,1), d2n, &
                       isym2(1), icemax, 1, idrop, iin, kin)
        END IF
        IF (isym2(2)==1)THEN
          CALL condevap_mass_eqn_and_remap1(xi(:,2), b21_my(:,2), &
                       fi2(:,2), psi2(:,2), d2n, &
                       isym2(2), icemax, 2, idrop, iin, kin)
        END IF
        IF (isym2(3)==1)THEN
          CALL condevap_mass_eqn_and_remap1(xi(:,3), b21_my(:,3), &
                       fi2(:,3), psi2(:,3), d2n, &
                       isym2(3), icemax, 3, idrop, iin ,kin)
        END IF
      END IF

      IF (isym3 == 1) THEN
        ! ... snow
        idrop = 0
        fl3 = 0.0_wp
        CALL condevap_mass_eqn_and_remap1(xs, b31_my, &
                       fi3, psi3, d2n, &
                       isym3, 1, 3, idrop, iin ,kin)
      END IF

      IF (isym4 == 1) THEN
        ! ... graupels (only_ice: evaporation)
        ! ... new jerdfun
        idrop = 0
        fl4 = 0.0_wp
        CALL condevap_mass_eqn_and_remap1(xg, b41_my, &
                       fi4, psi4, d2n, &
                       isym4, 1, 4, idrop, iin ,kin)
      END IF

      IF (isym5 == 1) THEN
        ! ... hail (only_ice: evaporation)
        ! ... new jerdfun
        idrop = 0
        fl5 = 0.0_wp
        CALL condevap_mass_eqn_and_remap1(xh, b51_my, &
                       fi5, psi5, d2n, &
                       isym5, 1, 5, idrop, iin ,kin)
      END IF

      IF ((del2.LT.0.AND.del2n.GT.0) .AND. ABS(del2n).GT.epsdel) THEN
        CALL finish(TRIM(modname),"fatal error in condevap_ice (del2.LT.0.AND.del2n.GT.0), model stop")
      END IF

      ! in case : kcond.ne.21
    END IF

    ! masses
    rmassibb=0.0_wp
    rmassiaa=0.0_wp

    DO k=1,nkr
      DO ice = 1,icemax ! sum over each of crystal types
        rmassibb = rmassibb + fi2(k,ice)*xi(k,ice)*xi(k,ice)
      END DO
      rmassibb=rmassibb+fi3(k)*xs(k)*xs(k)
      rmassibb=rmassibb+fi4(k)*xg(k)*xg(k)
      rmassibb=rmassibb+fi5(k)*xh(k)*xh(k)
    END DO
    rmassibb=rmassibb*col*3.0_wp/ror        ! total mass of solid particles before diffusional growth
    IF (rmassibb.LT.0.0_wp) rmassibb=0.0_wp

    DO k=1,nkr
      DO ice =1,icemax
        rmassiaa=rmassiaa+psi2(k,ice)*xi(k,ice)*xi(k,ice)
      END DO
      rmassiaa=rmassiaa+psi3(k)*xs(k)*xs(k)
      rmassiaa=rmassiaa+psi4(k)*xg(k)*xg(k)
      rmassiaa=rmassiaa+psi5(k)*xh(k)*xh(k)
    END DO
    rmassiaa = rmassiaa*col*3.0_wp/ror ! total mass of solid particles after diffusional growth

    IF (rmassiaa.LT.0.0_wp) rmassiaa=0.0_wp

    delmassi1 = rmassiaa-rmassibb
    qpn = qps-delmassi1

    tpn = tps + al2*delmassi1*cpd/cvd ! new temperature after substep (latent heat release or evaporation)

    IF (ABS(al2*delmassi1) > 5.0_wp )THEN
      PRINT*,"condevap_ice-out"
      PRINT*,"i=",iin,"kin",kin,"del1n,del2n,d2n,rw,pw,ri,pi,dt"
      PRINT*,del1n,del2n,d2n,rw,pw,ri,pi,dtt
      PRINT*,"tps=",tps,"qps=",qps,"delmassi1",delmassi1
      PRINT*,"al2=",al2,rmassibb,rmassiaa,"fi2_1",fi2(:,1),"fi2_2",fi2(:,2),"fi2_3",fi2(:,3)
      PRINT*,"fi3",fi3,"fi4",fi4,"fi5",fi5,"psi2_1",psi2(:,1),"psi2_2",psi2(:,2)
      PRINT*,"psi2_3",psi2(:,3),"psi3",psi3,"psi4",psi4,"psi5",psi5
      IF (ABS(al2*delmassi1) > 5.0_wp )THEN
        CALL finish(TRIM(modname),"fatal error in condevap_ice-out (ABS(al2*delmassi1) > 5.0_wp), model stop")
      END IF
    END IF

    ! ... supersaturation
    es1n=10.0_wp*sat_pres_water(tpn) !dynes/cm^2
    es2n=10.0_wp*sat_pres_ice(tpn)   !dynes/cm^2
    ew1n=10.0_wp*qpn*rho*rv*tpn      !dynes/cm^2

    IF (es1n == 0.0_wp)THEN
      del1n=0.5_wp
      div1=1.5_wp
      PRINT*,'es1n condevap_ice = 0'
      CALL finish(TRIM(modname),"fatal error in condevap_ice (es1n.eq.0), model stop")
    ELSE
      div1=ew1n/es1n
      del1n=ew1n/es1n-1.0_wp
    END IF
    IF (es2n == 0.0_wp)THEN
      del2n=0.5_wp
      div2=1.5_wp
      PRINT*,'es2n condevap_ice = 0'
      CALL finish(TRIM(modname),"fatal error in condevap_ice (es2n.eq.0), model stop")
    ELSE
      del2n=ew1n/es2n-1.0_wp
      div2=ew1n/es2n
    END IF

    !  end of time splitting
    ! (only ice: condensation or evaporation)
    IF (timenew.LT.dt) GOTO 46

    tt=tpn
    qq=qpn
    ! final PSD of 3 ice crystal types, as well as snow, graupel, hail:
    DO kr=1,nkr
      DO ice=1,icemax
        ff2(kr,ice)=psi2(kr,ice)
      END DO
      ff3(kr)=psi3(kr)
      ff4(kr)=psi4(kr)
      ff5(kr)=psi5(kr)
    END DO

    supsat_out=del1n

    RETURN
  END SUBROUTINE condevap_ice

  SUBROUTINE condevap_mixed &
                    & (tt,qq,pp,vr1,vr2,vr3,vr4,vr5,del1n,del2n,div1,div2 &
                    & ,ff1,psi1 &
                    & ,ff2,psi2 &
                    & ,ff3,psi3 &
                    & ,ff4,psi4 &
                    & ,ff5,psi5 &
                    & ,dtcond,isym1,isym2,isym3,isym4,isym5 &
                    & ,iin,kin,supsat_out,rho)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(INOUT) :: tt, qq, del1n, del2n, div1, div2, &
                              & ff1(nkr), ff2(nkr,icemax), ff3(nkr), ff4(nkr), ff5(nkr), &
                              & psi1(nkr), psi2(nkr,icemax), psi3(nkr), psi4(nkr), psi5(nkr), supsat_out
    REAL(KIND=wp),INTENT(IN) :: pp, vr1(nkr), vr2(nkr,icemax), vr3(nkr), vr4(nkr), vr5(nkr), &
                              & rho, dtcond
    INTEGER, INTENT(INOUT) :: isym1,isym2(icemax),isym3,isym4,isym5
    INTEGER, INTENT(IN) :: iin, kin
    INTEGER k,kr,ice,itime,kcond,idrop
    REAL(kind=wp) :: ror, fi1(nkr), fi2(nkr,icemax), fi3(nkr), fi4(nkr), fi5(nkr), &
                   & dt,dtt,b6,del1,del2,del1s,del2s,timenew,timerev, &
                   & sfnl,sfni,sfnii1,sfn11,sfn12,sfn21,sfn22,sfn31,sfn41,sfn51, &
                   & b5l,b5i,b7l,b7i,dopl,dopi,rw,ri,qw,pw, &
                   & pi,qi,dtnewl,d1n,d2n,tpn,tps,qpn,qps,rmasslaa,rmasslbb,rmassiaa,rmassibb, &
                   & es1n,es2n,ew1n,oper2,ar1,deltaq1,delmassi1,delmassl1, &
                   & b11_my(nkr),b21_my(nkr,icemax),b31_my(nkr),b41_my(nkr),b51_my(nkr), &
                   & sfni1(icemax), fl1(nkr), sfndummy(3), fl3(nkr), fl4(nkr), fl5(nkr), &
                   & del1_d, del2_d, rw_d, pw_d, ri_d, pi_d, d1n_d, d2n_d

    oper2(ar1)=0.622_wp/(0.622_wp+0.378_wp*ar1)/ar1

    ror=rho*0.001_wp !total density in cgs

    b11_my = 0.0_wp
    b21_my = 0.0_wp
    b31_my = 0.0_wp
    b41_my = 0.0_wp
    b51_my = 0.0_wp

    sfn11 = 0.0_wp
    sfni1 = 0.0_wp
    sfn31 = 0.0_wp
    sfn41 = 0.0_wp
    sfn51 = 0.0_wp

    itime = 0
    timenew = 0.0_wp
    dt = dtcond
    dtt = dtcond
    b6=0.0_wp

    tpn=tt
    qpn=qq

16  itime = itime + 1

    IF ((tpn-tmelt).GE.-0.187_wp) GOTO 17
    timerev = dt - timenew
    del1 = del1n
    del2 = del2n
    del1s = del1n
    del2s = del2n

    tps = tpn
    qps = qpn
    DO kr = 1,nkr
      fi1(kr) = psi1(kr)
      fi3(kr) = psi3(kr)
      fi4(kr) = psi4(kr)
      fi5(kr) = psi5(kr)
      DO ice = 1,icemax
        fi2(kr,ice) = psi2(kr,ice)
      END DO
    END DO

    IF (isym1 == 1) THEN !water droplets exist
      fl1 = 0.0_wp
      ! calculation of F_1 in (3.10) Khain&Sednev,1996:
      CALL condevap_dmdt_coef (xl,tps,pp,vr1,rlec,ro1bl,b11_my,1,1,fl1)
      sfndummy(1) = sfn11
      ! calculation of coefficients in equations for supersaturation
      CALL condevap_dsdt_coef(fi1,xl,sfndummy,b11_my,1.0_wp/ror,1)
      sfn11 = sfndummy(1)
    END IF

    IF (sum(isym2) > 0) THEN !ice crystals
      fl1 = 0.0_wp
      CALL condevap_dmdt_coef (xi,tps,pp,vr2,riec,ro2bl,b21_my,3,2,fl1)
      CALL condevap_dsdt_coef (fi2,xi,sfni1,b21_my,1.0_wp/ror,icemax)
    END IF
    IF (isym3 == 1) THEN !snow
      fl3 = 0.0_wp
      CALL condevap_dmdt_coef (xs,tps,pp,vr3,rsec,ro3bl,b31_my,1,3,fl3)
      sfndummy(1) = sfn31
      CALL condevap_dsdt_coef(fi3,xs,sfndummy,b31_my,1.0_wp/ror,1)
      sfn31 = sfndummy(1)
    END IF
    IF (isym4 == 1) THEN !graupel
      fl4 = 0.0_wp
      CALL condevap_dmdt_coef(xg,tps,pp,vr4,rgec,ro4bl,b41_my,1,2,fl4)
      sfndummy(1) = sfn41
      CALL condevap_dsdt_coef(fi4,xg,sfndummy,b41_my,1.0_wp/ror,1)
      sfn41 = sfndummy(1)
    END IF
    IF (isym5 == 1) THEN !hail
      fl5 = 0.0_wp
      CALL condevap_dmdt_coef(xh,tps,pp,vr5,rhec,ro5bl,b51_my,1,2,fl5)
      sfndummy(1) = sfn51
      CALL condevap_dsdt_coef(fi5,xh,sfndummy,b51_my,1.0_wp/ror,1)
      sfn51 = sfndummy(1)
    END IF

    sfnii1 = sfni1(1) + sfni1(2) + sfni1(3)
    sfn21 = sfnii1 + sfn31 + sfn41 + sfn51
    sfn12 = 0.0_wp
    sfnl = sfn11 + sfn12
    sfn22 = 0.0_wp
    sfni = sfn21 + sfn22

    b5l=bb1_my/tps/tps
    b5i=bb2_my/tps/tps
    b7l=b5l*b6
    b7i=b5i*b6
    dopl=1.+del1s
    dopi=1.+del2s

    ! after each substep, the new sup sat over water and ice are calculated. These values are used to
    ! calculate the new temperature and humidity at each substep. Therefore sup sat over ice is calculated even
    ! when there is no ice. The procedure can in future be simplified by assuming a linear changes of T and Q

    IF (.not. latheatfac1) THEN
      rw=(oper2(qps)+b5l*al1)*dopl*sfnl !R1 coeff in (3.11) in Khain&Sednev,1996, al1=Lw/Cp
      ri=(oper2(qps)+b5l*al2)*dopl*sfni !R2 coeff in (3.11) in Khain&Sednev,1996, al2=Li/Cp
    ELSE
      rw=(oper2(qps)+b5l*al1*cpd/cvd)*dopl*sfnl !R1 coeff in (3.11) in Khain&Sednev,1996, al1=Lw/Cp
      ri=(oper2(qps)+b5l*al2*cpd/cvd)*dopl*sfni !R2 coeff in (3.11) in Khain&Sednev,1996, al2=Li/Cp
    END IF
    IF (.not. latheatfac1) THEN
      pw=(oper2(qps)+b5i*al1)*dopi*sfnl !P1 coeff in (3.11) in Khain&Sednev,1996
      pi=(oper2(qps)+b5i*al2)*dopi*sfni !P2 coeff in (3.11) in Khain&Sednev,1996
    ELSE
      pw=(oper2(qps)+b5i*al1*cpd/cvd)*dopi*sfnl !P1 coeff in (3.11) in Khain&Sednev,1996
      pi=(oper2(qps)+b5i*al2*cpd/cvd)*dopi*sfni !P2 coeff in (3.11) in Khain&Sednev,1996
    END IF
    qw=b7l*dopl
    qi=b7i*dopi

    IF (rw.ne.rw .OR. pw.ne.pw)THEN
      PRINT*, 'nan in condevap_mixed'
      CALL finish(TRIM(modname),"fatal error in condevap_mixed (rw or pw are nan), model stop")
    END IF

    ! del1 > 0, del2 < 0    (antibergeron mixed phase - kcond=50)
    ! del1 < 0 and del2 < 0 (evaporation mixed_phase - kcond=30)
    ! del1 > 0 and del2 > 0 (condensation mixed phase - kcond=31)
    ! del1 < 0, del2 > 0    (bergeron mixed phase - kcond=32)

    ! every ncond substep can be still too large for droplets growth. We want that during diffusional
    ! growth, the droplet radius change will be less THEN few bins. Therefore ncond substep is further
    ! divided for steps not less than 0.4 sec. In case that after several steps, time will be still
    ! lower than ncond sub step, the last step will have this delta time. Note that this is done without updates of t,q,s
    kcond=50 ! kcond is a flag
    IF (del1n .LT. 0.0_wp .AND. del2n .LT. 0.0_wp) kcond=30
    IF (del1n .GT. 0.0_wp .AND. del2n .GT. 0.0_wp) kcond=31 !AK 20220929: is it diffusional growth or evaporation? --> kcond
    IF (del1n .LT. 0.0_wp .AND. del2n .GT. 0.0_wp) kcond=32

    IF (kcond == 50) THEN
      dtnewl = dt
      dtnewl = MIN(dtnewl,timerev) !AK 20220929: timerev are small substeps within ncond sub step
      timenew = timenew + dtnewl
      dtt = dtnewl

      ! ... in case the anti-bregeron regime we do not call diffusional-growth
      !     (for some strange reason sup sat over water > sup sat over ice)
      PRINT*, "anti-bregeron regime, no diffu"
      PRINT*,  del1, del2, tt, qq, kin
      GOTO 17                           ! leave this SUBROUTINE
      ! in case : kcond = 50
    END IF
    IF (kcond == 31) THEN
      ! ... del1 > 0 and del2 > 0 (condensation mixed phase - kcond=31)
      ! ... condensation mixed phase:
      dtnewl = dt
      dtnewl = MIN(dtnewl,timerev)
      timenew = timenew + dtnewl
      dtt = dtnewl
    END IF
    IF (kcond == 30) THEN
      ! ... del1 < 0 and del2 < 0 (evaporation mixed_phase - kcond=30)
      ! ... evaporation mixed phase:
      dtnewl = dt
      dtnewl = MIN(dtnewl,timerev)
      timenew = timenew + dtnewl
      dtt = dtnewl
    END IF
    IF (kcond == 32) THEN
      ! ... IF (del1n < 0.0_wp .AND. del2n > 0.0_wp) kcond=32
      ! ... bergeron mixed phase:
      dtnewl = dt
      dtnewl = MIN(dtnewl,timerev)
      timenew = timenew + dtnewl
      dtt = dtnewl
    END IF

    IF (dtt < 0.0_wp) CALL finish(TRIM(modname),"fatal error in condevap_mixed:(dtt<0), model stop")

    del1_d = del1
    del2_d = del2
    rw_d = rw
    pw_d = pw
    ri_d = ri
    pi_d = pi

    ! solving the equation for 2 supersaturations del1n,del2n:
    CALL condevap_supsat_eqn(del1_d,del2_d,del1n,del2n, &
                      rw_d,pw_d,ri_d,pi_d, &
                      dtt,d1n_d,d2n_d,0.0_wp,0.0_wp, &
                      isym1,isym2,isym3,isym4,isym5)
    del1 = del1_d
    del2 = del2_d
    rw = rw_d
    pw = pw_d
    ri = ri_d
    pi = pi_d
    d1n = d1n_d
    d2n = d2n_d

    IF (isym1 == 1) THEN ! water exists
      idrop = isym1
      fl1 = 0.0_wp
      CALL condevap_mass_eqn_and_remap1(xl, b11_my, &
                   fi1, psi1, d1n, &
                   isym1, 1, 1, idrop, iin, kin)
    END IF

    IF (sum(isym2) > 0) THEN !crystals
      idrop = 0
      fl1 = 0.0_wp
      IF (isym2(1)==1) THEN ! columns
        CALL condevap_mass_eqn_and_remap1(xi(:,1), b21_my(:,1), &
                  fi2(:,1), psi2(:,1), d2n, &
                  isym2(1), icemax, 1, idrop, iin ,kin)
      END IF

      IF (isym2(2)==1)THEN ! plates
        CALL condevap_mass_eqn_and_remap1(xi(:,2), b21_my(:,2), &
                  fi2(:,2), psi2(:,2), d2n, &
                  isym2(2), icemax, 2, idrop, iin ,kin)
      END IF

      IF (isym2(3)==1)THEN ! dendrites
        CALL condevap_mass_eqn_and_remap1(xi(:,3), b21_my(:,3), &
                  fi2(:,3), psi2(:,3), d2n, &
                  isym2(3), icemax, 3, idrop, iin ,kin)
      END IF
    END IF

    IF (isym3 == 1) THEN !snow
      idrop = 0
      fl3 = 0.0_wp
      CALL condevap_mass_eqn_and_remap1(xs, b31_my, &
                  fi3, psi3, d2n, &
                  isym3, 1, 3, idrop, iin ,kin)
    END IF

    IF (isym4 == 1) THEN !graupels
      idrop = 0
      fl4 = 0.0_wp
      CALL condevap_mass_eqn_and_remap1(xg, b41_my, &
                  fi4, psi4, d2n, &
                  isym4, 1, 4, idrop, iin ,kin)
    END IF

    IF (isym5 == 1) THEN !hail
      idrop = 0
      fl5 = 0.0_wp
      CALL condevap_mass_eqn_and_remap1(xh, b51_my, &
                  fi5, psi5, d2n, &
                  isym5, 1, 5, idrop, iin ,kin)
    END IF

    rmasslbb=0.0_wp
    rmassibb=0.0_wp
    rmasslaa=0.0_wp
    rmassiaa=0.0_wp

    DO k=1,nkr
      rmasslbb=rmasslbb+fi1(k)*xl(k)*xl(k) !calculate cloud water content before diffusional
                                           !growth ncond substep, xl(k) !mass of bin k , named xl outside
      DO ice =1,icemax
        rmassibb=rmassibb+fi2(k,ice)*xi(k,ice)*xi(k,ice) !calculate cloud ice content before diffusional growth
                                                         !(deposition/sublimation) ncond substep
      END DO
      rmassibb=rmassibb+fi3(k)*xs(k)*xs(k)
      rmassibb=rmassibb+fi4(k)*xg(k)*xg(k)
      rmassibb=rmassibb+fi5(k)*xh(k)*xh(k)
    END DO
    rmassibb=rmassibb*col*3.0_wp/ror
    IF (rmassibb.LT.0.0_wp) rmassibb=0.0_wp
    rmasslbb=rmasslbb*col*3.0_wp/ror
    IF (rmasslbb.LT.0.0_wp) rmasslbb=0.0_wp

    DO k=1,nkr
      rmasslaa=rmasslaa+psi1(k)*xl(k)*xl(k)
      DO ice =1,icemax
        rmassiaa=rmassiaa+psi2(k,ice)*xi(k,ice)*xi(k,ice)
      END DO
      rmassiaa=rmassiaa+psi3(k)*xs(k)*xs(k)
      rmassiaa=rmassiaa+psi4(k)*xg(k)*xg(k)
      rmassiaa=rmassiaa+psi5(k)*xh(k)*xh(k)
    END DO
    rmassiaa=rmassiaa*col*3.0_wp/ror !result: cloud ice content after diffusional growth (deposition/sublimation) ncond substep
    IF (rmassiaa.LE.0.0_wp) rmassiaa=0.0_wp
    rmasslaa=rmasslaa*col*3.0_wp/ror !result: cloud water content after diffusional growth ncond substep

    IF (rmasslaa.LT.0.0_wp) rmasslaa=0.0_wp

    delmassl1=rmasslaa-rmasslbb !cloud water mass change during substep (new-old)
    delmassi1=rmassiaa-rmassibb !cloud ice mass change during substep (new-old)
    deltaq1=delmassl1+delmassi1
    qpn=qps-deltaq1 !new specific humidity

    tpn = tps + (al1*delmassl1+al2*delmassi1)*cpd/cvd ! new temperature after substep (latent heat release or evaporation)

    IF (ABS(al1*delmassl1+al2*delmassi1) > 5.0_wp )THEN
      PRINT*,"condevap_mixed-input"
      PRINT*,"pp",pp,"ror",ror,'vr1',vr1,'div1',div1,'div2',div2
      PRINT*,'rlec',rlec,'ro1bl',ro1bl,'const',col,iin,kin
      PRINT*,'isym1',isym1,'isym2',isym2,'isym3',isym3,'isym4',isym4,'isym5',isym5
      PRINT*,"condevap_mixed-out"
      PRINT*,"i=",iin,"kin",kin,"del1n,del2n,d1n,d2n,rw,pw,ri,pi,dt"
      PRINT*,del1n,del2n,d1n,d2n,rw,pw,ri,pi,dtt
      PRINT*,"tps=",tps,"tpn=",tpn,"qps=",qps,"delmassl1",delmassl1,"delmassi1",delmassi1
      PRINT*,"al2=",al2,"al1=",al1,rmasslaa,rmasslbb,rmassiaa,rmassibb
      PRINT*,"fi1",fi1,"fi3",fi3,"fi4",fi4,"fi5",fi5,"psi1",psi1,"psi3",psi3,"psi4",psi4,"psi5",psi5
      IF (ABS(al1*delmassl1+al2*delmassi1) > 5.0_wp )THEN
        CALL finish(TRIM(modname),"fatal error in condevap_mixed-out (ABS(al1*delmassl1+al2*delmassi1) > 5.0_wp), model stop")
      END IF
    END IF

    ! supersaturation over water and ice
    es1n=10.0_wp*sat_pres_water(tpn) !dynes/cm^2
    es2n=10.0_wp*sat_pres_ice(tpn)   !dynes/cm^2
    ew1n=10.0_wp*qpn*rho*rv*tpn      !dynes/cm^2

    IF (es1n == 0.0_wp) THEN
      del1n=0.5_wp
      div1=1.5_wp
      PRINT*,'es1n condevap_mixed = 0'
      CALL finish(TRIM(modname),"fatal error in condevap_mixed (es1n.eq.0), model stop")
    ELSE
      div1=ew1n/es1n
      del1n=ew1n/es1n-1.0_wp
    END IF
    IF (es2n == 0.0_wp) THEN
      del2n=0.5_wp
      div2=1.5_wp
      PRINT*,'es2n condevap_mixed = 0'
      CALL finish(TRIM(modname),"fatal error in condevap_mixed (es2n.eq.0), model stop")
    ELSE
      del2n=ew1n/es2n-1.0_wp
      div2=ew1n/es2n
    END IF

    ! At this point in condevap_water, there is a calculation
    ! of the full integral over the sup sat
    ! (dsupintw), without ncond substeps, and THEN
    ! perform a single remapping, instead of doing
    ! remapping after each substep.
    ! Here in condevap_mixed we don't need it

    ! end of time splitting

    IF (timenew < dt) GOTO 16
17  CONTINUE

    tt=tpn
    qq=qpn
    DO kr=1,nkr
      ff1(kr)=psi1(kr)
      DO ice=1,icemax
        ff2(kr,ice)=psi2(kr,ice)
      END DO
      ff3(kr)=psi3(kr)
      ff4(kr)=psi4(kr)
      ff5(kr)=psi5(kr)
    END DO

    supsat_out=del1n

    RETURN
  END SUBROUTINE condevap_mixed

  SUBROUTINE coallescence(ff1r, ff2r, ff3r, ff4r, ff5r, tt, qq, pp, rho, dtcoll, iin, kin) !cwgg,cwhh are missing

    ! lwf=liquid water fraction - this option comes from HUCM and is not implemented here
    ! use module_mp_sbm_collision,only:coll_xyy_lwf,coll_xyx_lwf,coll_xxx_lwf, &
    ! coll_xyz_lwf, coll_ice_stick_eff, coll_breakup, coll_xxy_lwf, coll_xxx_bott

    ! uses SUBROUTINEs:
    ! coll_xxx_bott: collision of large particle X with small particle X, resulting particle X (e.g. water+water->water)
    ! coll_xyy_bott: collision of large particle X with small particle Y, resulting particle Y (e.g. water+ice-> ice)
    ! coll_xyx_bott: collision of large particle X with small particle Y, resulting particle X (e.g. ice+water-> ice)
    ! coll_xyz_bott: collision of large particle X with small particle Y, resulting particle Z (e.g. water+graupel->hail)
    ! coll_xyz_lwf: is able to account water rimed fraction (liquid water fraction), but in this code is used without this option

    !coll   subroutine     indc
    !----   ----------     ----
    !cwll   coll_xxx_bott
    !
    !cwls   coll_xyz_lwf   0
    !cwsl   coll_xyx_lwf   1
    !cwsl   coll_xyz_lwf   1
    !
    !cwhs   coll_xyz_lwf   0
    !cwgs   coll_xyz_lwf   0
    !
    !cwlg   coll_xyy_lwf   0
    !cwgl   coll_xyx_lwf   1
    !cwlg   coll_xyz_lwf   0  effectively turned off
    !cwgl   coll_xyz_lwf   1  effectively turned off
    !
    !cwlh   coll_xyy_lwf   0
    !cwhl   coll_xyx_lwf   1
    !cwss   coll_xxx_lwf

    ! ff1r - PSD of water
    ! ff2r - PSD of crystals (plates, dendrites, columns)
    ! ff3r - PSD of snow
    ! ff4r - PSD of graupel
    ! ff5r - PSD of hail

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: iin,kin
    REAL(KIND=wp),INTENT(IN) :: pp,rho,dtcoll
    REAL(KIND=wp),INTENT(INOUT) :: ff1r(:),ff2r(:,:),ff3r(:),ff4r(:),ff5r(:),tt,qq
    INTEGER :: kr,ice,icol_drop,icol_snow,icol_graupel,icol_hail,icol_drop_brk, it, ndiv, kr1, nkf, do_more_collisions
    REAL(KIND=wp) :: g1(nkr),g2(nkr,icemax),g3(nkr),g4(nkr),g5(nkr), &
                  gdumb(nkr), cont_fin_drop,conc_icempl,deldrop,t_new, &
                  cont_fin_ice,conc_old,conc_new,cont_init_drop,alwc,pp_r, &
                  break_drop_bef,break_drop_aft,dtbreakup,break_drop_per, prdkrn, cwxx(nkr,nkr), &
                  g1p(nkr+1),x1p(nkr+1),x2p(nkr+1),g3p(nkr+1),g4p(nkr+1),g5p(nkr+1)
    REAL(KIND=wp), PARAMETER :: prdkrn1 = 1.0_wp, epsilbreak=0.001_wp

    icol_drop=0
    icol_snow=0
    icol_graupel=0
    icol_hail=0

    t_new = tt
    pp_r = pp

    IF ( iceprocs==1 ) THEN
      ! in case of ice, sticking efficiency depends on temperature (ice doesnt stick well at low T):
      CALL coll_ice_stick_eff(tt,qq,pp,prdkrn)
    END IF

    ! the calculation of collisions will be done for distribution funcions with respect to masses (mg)
    ! and therefore we do here translation from PSD (integral=concentration) to PSD (integral==mass content):

    DO kr=1,nkr !loop over all bins
      g1(kr) = ff1r(kr)*3.0_wp*xl(kr)*xl(kr)*1.0e3_wp        ! *1.e3 for g->mg
      IF ( iceprocs==1 ) THEN
        !in fast_sbm there are no crystals, so ff2r=0. Nucleated crystals (in condevap_ice) are transformed to snow ff3r
        !here is a general version which includes crystals (3 PSD) and graupel/hail PSD, instead of simply 2: for snow and graupel/hail:
        g2(kr,1)= ff2r(kr,1)*3.0_wp*xi(kr,1)*xi(kr,1)*1.0e3_wp ! 0 in this version
        g2(kr,2)= ff2r(kr,2)*3.0_wp*xi(kr,2)*xi(kr,2)*1.0e3_wp ! 0 in this version
        g2(kr,3)= ff2r(kr,3)*3.0_wp*xi(kr,3)*xi(kr,3)*1.0e3_wp ! 0 in this version
        g3(kr)  = ff3r(kr)*3.0_wp*xs(kr)*xs(kr)*1.e3_wp
        g4(kr)  = ff4r(kr)*3.0_wp*xg(kr)*xg(kr)*1.e3_wp ! 1 if hail_opt=0, 0 if hail_opt=1
        g5(kr)  = ff5r(kr)*3.0_wp*xh(kr)*xh(kr)*1.e3_wp ! 1 if hail_opt=1, 0 if hail_opt=0
      END IF

      ! check whether collision/breakup is possible for each type, i.e.
      ! if at least one bin is non-zero, we will call collision for this type:
      IF (g1(kr) > epsil) icol_drop=1
      IF ( iceprocs==1 ) THEN
        IF (g3(kr).GT.epsil) icol_snow = 1
        IF (g4(kr).GT.epsil) icol_graupel = 1
        IF (g5(kr).GT.epsil) icol_hail = 1
      END IF
    END DO
    ! here icol_snow=1 if there is some snow,
    ! here icol_graupel=1 if there is some graupel and hail_opt=0,
    ! here icol_hail=1 if there is some graupel and hail_opt=1

    IF ( iceprocs==1 ) THEN
      ! calculation of initial hydromteors content in g/cm**3 :
      cont_init_drop = 0.0_wp
      cont_init_drop = sum(g1(1:nkr))
      cont_init_drop=col*cont_init_drop*1.e-3_wp
      ! calculation of critical liquid water contnent which determines the
      ! type of resulting particle (graupel or snow): alwc in g/m**3
      alwc=cont_init_drop*1.e6_wp
    END IF

    ! ----------------------------
    ! ... drop-drop collisions
    ! ----------------------------

    IF (icol_drop == 1) THEN
      ! at this point, in case of rates_diag==1 there are dummy diagnostic calls to coll_xxx in operational warm version
      ! here we call the standard 'coll_xxx()' routine with the full spectrum (dsd+rsd)      ! in [g/cm3]
      ! input:  g1 - mass PSD, cwll - collision kernel liquid-liquid, xl_mg - mass bins:
      ! output: new g1:
      CALL kernel2d(cwxx,pp_r,"cwll")
      IF ( usenkrp1a .EQ. 0 ) THEN
        CALL coll_xxx_bott(g1,cwxx,xl_mg,1.0_wp,nkr) !cwll is reduced by coallescence efficiency in mo_sbm_util, even IF break-up is off
      ELSE
        CALL coll_add_pseodo_bin(x1p,g1p,nkf,xl_mg,g1)
        CALL coll_xxx_bott(g1p,cwxx,x1p,1.0_wp,nkf) !cwll is reduced by coallescence efficiency in mo_sbm_util, even IF break-up is off
        DO kr1=1,nkr
          g1(kr1)=g1p(kr1)
        END DO
        g1(nkr)=g1(nkr)+g1p(nkr+1)
      END IF

      ! --------------------------------------------------------
      ! ... probability of drop breakup after collision:
      ! --------------------------------------------------------

      ! There is an impirical distribution of fragments sizes after breakup
      ! Here we do iterative process to use this distribution with the restriction of mass conservation:

      icol_drop_brk=0
      DO kr=1,nkr !loop over all bins
        IF (kr > krmin_breakup .AND. g1(kr) > epsil) icol_drop_brk = 1
      END DO
      IF (ibreakup == 0) icol_drop_brk = 0

      IF (icol_drop_brk == 1) THEN
        ndiv = 1 !initial number of iterations, which will increase in case of mass conservation violation
        dtbreakup = dtcoll/ndiv ! timestep of breakup
10      CONTINUE
        DO it = 1,ndiv            ! ndiv=number of iterations
          IF (it == 1) THEN
            DO kr=1,nkr
              gdumb(kr)= g1(kr)*1.e-3_wp       !mg->g
            END DO
            break_drop_bef=0.0_wp
            DO kr=1,nkr
              break_drop_bef = break_drop_bef+g1(kr)*1.e-3_wp !total mass before in g
            END DO
          END IF
          CALL coll_breakup(gdumb, dtbreakup)
        END DO

        DO kr=1,nkr
          IF (gdumb(kr) < 0.0_wp) THEN
            IF (ndiv < 8) THEN
              ndiv = 2*ndiv   ! in this case calculate breakup again with twice more iterations
              GOTO 10
            ELSE              ! if ndiv is already large and still the final PSD is <0, give up:
              GOTO 11         ! in this case g1 will not be updated by gdumb
            END IF
          END IF
          IF (gdumb(kr) .NE. gdumb(kr)) THEN ! check for nan
            PRINT*,kr,gdumb(kr),xl(kr),it,ndiv, dtbreakup,gdumb
            CALL finish(TRIM(modname),"in coal_bott after coll_breakup - ff1r nan, model stop") !ff1r=g1/(3*xl*xl*1.e3) where g1 in mg
          END IF
          IF (gdumb(kr) .LT. 0.0_wp ) THEN ! check for < 0
            PRINT*,kr,gdumb(kr),xl(kr),it,ndiv, dtbreakup,gdumb
            CALL finish(TRIM(modname),"in coal_bott after coll_breakup - ff1r<0 , model stop")
          END IF
        END DO

        ! mass conservation:
        break_drop_aft=0.0_wp
        DO kr=1,nkr
          break_drop_aft=break_drop_aft+gdumb(kr)
        END DO
        break_drop_per=ABS(break_drop_aft/break_drop_bef-1.0_wp)

        IF ( break_drop_per > epsilbreak) THEN !check that the mass is not growing
          print*,'drop breakup details: ',ndiv,break_drop_bef,break_drop_aft,break_drop_per
          !CALL finish(TRIM(modname),"no mass conservation in breakup, model stop")
          ndiv=ndiv*2  ! in this case calculate breakup again with twice more iterations
          GOTO 10
        ELSE
          DO kr=1,nkr
            g1(kr) = gdumb(kr)*1.e3_wp !g->mg
          END DO
        END IF

        !clipping maximum g:
        DO kr=1,nkr
          IF ( g1(kr) .LT. epsil ) THEN
            g1(kr)=0.0_wp
          END IF
        END DO

      END IF ! IF icol_drop_brk.EQ.1
    END IF ! IF icol_drop.EQ.1

11  CONTINUE

    do_more_collisions=0
    IF (( positive_t_coll == 1 ) .AND. (iceprocs == 1)) THEN
      do_more_collisions=1 ! allow collisions with ice particles at positive temperatures
    ELSE IF (( positive_t_coll == 0 ) .AND. (iceprocs == 1) .AND. (tt <= tmelt)) THEN
      do_more_collisions=1 ! allow collisions with ice particles at negative temperatures only
    END IF

    IF ( do_more_collisions == 1 ) THEN

      IF (icol_drop == 1) THEN
        !        drop - snow = graupel/hail
        !        snow - drop = snow
        !        or
        !        snow - drop = graupel/hail
        IF (icol_snow == 1) THEN
          ! ------------------------------------
          ! Collisions between drops and snow
          ! ------------------------------------
          CALL kernel2d(cwxx,pp_r,"cwls")
          IF (hail_opt == 1) THEN       ! g1 (water) + g3 (snow) --> g5 (hail):
            IF ( usenkrp3 .EQ. 0 ) THEN
              ! when water is involved, sticking efficiency is prdkrn1=1
              CALL coll_xyz_lwf(g1,g3,g5,cwxx,xl_mg,xs_mg,prdkrn1,0,nkr) ! 0 means allow collision of same bins !cwls
            ELSE
              CALL coll_add_pseodo_3bin(x1p,x2p,g1p,g3p,g5p,nkf,xl_mg,xs_mg,g1,g3,g5)
              CALL coll_xyz_lwf(g1p,g3p,g5p,cwxx,x1p,x2p,prdkrn1,0,nkf) ! 0 means allow collision of same bins !cwls
              DO kr1=1,nkr
                g1(kr1)=g1p(kr1)
                g3(kr1)=g3p(kr1)
                g5(kr1)=g5p(kr1)
              END DO
              g1(nkr)=g1(nkr)+g1p(nkr+1)
              g3(nkr)=g3(nkr)+g3p(nkr+1)
              g5(nkr)=g5(nkr)+g5p(nkr+1)
            END IF
          ELSE                         ! g1 (water) + g3 (snow) --> g4 (graupel):
            IF ( usenkrp3 .EQ. 0 ) THEN
              ! when water is involved, sticking efficiency is prdkrn1=1
              CALL coll_xyz_lwf(g1,g3,g4,cwxx,xl_mg,xs_mg,prdkrn1,0,nkr) ! 0 means allow collision of same bins
            ELSE
              CALL coll_add_pseodo_3bin(x1p,x2p,g1p,g3p,g4p,nkf,xl_mg,xs_mg,g1,g3,g4)
              CALL coll_xyz_lwf(g1p,g3p,g4p,cwxx,x1p,x2p,prdkrn1,0,nkf) ! 0 means allow collision of same bins !cwls
              DO kr1=1,nkr
                g1(kr1)=g1p(kr1)
                g3(kr1)=g3p(kr1)
                g4(kr1)=g4p(kr1)
              END DO
              g1(nkr)=g1(nkr)+g1p(nkr+1)
              g3(nkr)=g3(nkr)+g3p(nkr+1)
              g4(nkr)=g4(nkr)+g4p(nkr+1)
            END IF
          END IF

          CALL kernel2d(cwxx,pp_r,"cwsl")
          IF (alwc < alcr) THEN         ! In case of low liquid water content, then: g3 (snow) + g1 (water) --> g3 (snow):
            IF ( usenkrp2 .EQ. 0 ) THEN
              ! when water is involved, sticking efficiency is prdkrn1=1
              CALL coll_xyx_lwf(g3,g1,cwxx,xs_mg,xl_mg,prdkrn1,1,nkr)  ! 1 means do not allow collision of same
                                                                       ! bins since water+snow collided above !cwsl
            ELSE
              CALL coll_add_pseodo_2bin(x1p,x2p,g3p,g1p,nkf,xs_mg,xl_mg,g3,g1)
              CALL coll_xyx_lwf(g3p,g1p,cwxx,x1p,x2p,prdkrn1,1,nkf)
              DO kr1=1,nkr
                g1(kr1)=g1p(kr1)
                g3(kr1)=g3p(kr1)
              END DO
              g1(nkr)=g1(nkr)+g1p(nkr+1)
              g3(nkr)=g3(nkr)+g3p(nkr+1)
            END IF
          ELSE                          ! In case of high liquid water content, then: g3 (snow) + g1 (water) --> hail or graupel:
            IF (hail_opt == 1) THEN
              IF ( usenkrp3 .EQ. 0 ) THEN
                ! when water is involved, sticking efficiency is prdkrn1=1
                CALL coll_xyz_lwf(g3,g1,g5,cwxx,xs_mg,xl_mg,prdkrn1,1,nkr) ! 1 means do not allow collision of same
                                                                           ! bins since water+snow collided above !cwsl
              ELSE
                CALL coll_add_pseodo_3bin(x1p,x2p,g3p,g1p,g5p,nkf,xs_mg,xl_mg,g3,g1,g5)
                CALL coll_xyz_lwf(g3p,g1p,g5p,cwxx,x1p,x2p,prdkrn1,1,nkf)
                DO kr1=1,nkr
                  g1(kr1)=g1p(kr1)
                  g3(kr1)=g3p(kr1)
                  g5(kr1)=g5p(kr1)
                END DO
                g1(nkr)=g1(nkr)+g1p(nkr+1)
                g3(nkr)=g3(nkr)+g3p(nkr+1)
                g5(nkr)=g5(nkr)+g5p(nkr+1)
              END IF
            ELSE
              IF ( usenkrp3 .EQ. 0 ) THEN
                ! when water is involved, sticking efficiency is prdkrn1=1
                CALL coll_xyz_lwf(g3,g1,g4,cwxx,xs_mg,xl_mg,prdkrn1,1,nkr) ! 1 means do not allow collision of same
                                                                           ! bins since water+snow collided above !cwsl
              ELSE
                CALL coll_add_pseodo_3bin(x1p,x2p,g3p,g1p,g4p,nkf,xs_mg,xl_mg,g3,g1,g4)
                CALL coll_xyz_lwf(g3p,g1p,g4p,cwxx,x1p,x2p,prdkrn1,1,nkf)
                DO kr1=1,nkr
                  g1(kr1)=g1p(kr1)
                  g3(kr1)=g3p(kr1)
                  g4(kr1)=g4p(kr1)
                END DO
                g1(nkr)=g1(nkr)+g1p(nkr+1)
                g3(nkr)=g3(nkr)+g3p(nkr+1)
                g4(nkr)=g4(nkr)+g4p(nkr+1)
              END IF
            END IF
          END IF
        END IF
      END IF

      IF (icol_graupel == 1) THEN
        ! ---------------------------------------
        ! Collisions between drops and graupel
        ! ---------------------------------------
        ! In this version of Fast-SBM, we have either Graupel or Hail as our dense particle.
        ! That is, in addition to aerosols and drops, we use additional 2 size distribution
        ! (SD; advected & mixed): continuous ice-crystals/snow and graupel or hail
        ! (total of 4 advected SD in Fast-SBM).
        ! This version is intended to describe two generic case: continental deep convection
        ! where hail is more abundant, and maritime deep convection where particles more like
        ! graupel are abundant. This is done simply for getting a compact/representative setup
        ! which is more affordable from computational perspective, that is- there is no problem
        ! to add graupel (or hail) as standalone additional SD.
        ! Therefore, in case we have either graupel or hail, once some of the collisions process
        ! are designed to transfer mass between these two distinct particles categories, we prevent
        ! it from happening, so we do not have mass/number leakage to a non-existed category.
        ! That means, when the user chooses to use with graupel, when graupel grows while colliding
        ! with drops, it remain in it's category as increased size of graupel. In reality, when
        ! the density of graupel increases it becomes hail.
        ! We have plans to update also the Full-SBM where all SD are explicitly represented and
        ! interacting (CCN, drops, 3xice-crystals, snow, graupel, hail - total of 8 SD).

        !     drops - graupel = graupel
        !     graupel - drops = graupel
        !     drops - graupel = hail (no transition in fsbm)
        !     graupel - drop = hail (no transition in fsbm)

        IF (alwc < alcr_g) THEN
          ! In case of low liquid water content, THEN: g1 (water) + g4 (graupel) --> g4 (graupel):
          CALL kernel2d(cwxx,pp_r,"cwlg")
          IF ( usenkrp2 .EQ. 0 ) THEN
            CALL coll_xyy_lwf(g1,g4,cwxx,xl_mg,xg_mg,prdkrn1,0,nkr) !cwlg
          ELSE
            CALL coll_add_pseodo_2bin(x1p,x2p,g1p,g4p,nkf,xl_mg,xg_mg,g1,g4)
            CALL coll_xyy_lwf(g1p,g4p,cwxx,x1p,x2p,prdkrn1,0,nkf)
            DO kr1=1,nkr
              g1(kr1)=g1p(kr1)
              g4(kr1)=g4p(kr1)
            END DO
            g1(nkr)=g1(nkr)+g1p(nkr+1)
            g4(nkr)=g4(nkr)+g4p(nkr+1)
          END IF
          ! ... for ice multiplication
          conc_old = 0.0_wp
          DO kr = kr_icempl,nkr
            conc_old = conc_old+col*g1(kr)/xl_mg(kr)
          END DO

          ! In case of low liquid water content, THEN: g4 (graupel) + g1 (water) --> g4 (graupel):
          CALL kernel2d(cwxx,pp_r,"cwgl")
          IF ( usenkrp2 .EQ. 0 ) THEN
            CALL coll_xyx_lwf(g4,g1,cwxx,xg_mg,xl_mg,prdkrn1,1,nkr) ! 1 means do not allow collision of same
                                                                    ! bins since water+snow collided above !cwgl
          ELSE
            CALL coll_add_pseodo_2bin(x1p,x2p,g4p,g1p,nkf,xg_mg,xl_mg,g4,g1)
            CALL coll_xyx_lwf(g4p,g1p,cwxx,x1p,x2p,prdkrn1,1,nkf)
            DO kr1=1,nkr
              g1(kr1)=g1p(kr1)
              g4(kr1)=g4p(kr1)
            END DO
            g1(nkr)=g1(nkr)+g1p(nkr+1)
            g4(nkr)=g4(nkr)+g4p(nkr+1)
          END IF

        ELSE ! effectively turned off in this version

          ! In case of high liquid water content, THEN: g1 (water) + g4 (graupel) --> g5 (hail):
          CALL kernel2d(cwxx,pp_r,"cwlg")
          IF ( usenkrp3 .EQ. 0 ) THEN
            CALL coll_xyz_lwf(g1,g4,g5,cwxx,xl_mg,xg_mg,prdkrn1,0,nkr) ! 1 means do not allow collision of same
                                                                       ! bins since water+snow collided above !cwlg
          ELSE
            CALL coll_add_pseodo_3bin(x1p,x2p,g1p,g4p,g5p,nkf,xl_mg,xg_mg,g1,g4,g5)
            CALL coll_xyz_lwf(g1p,g4p,g5p,cwxx,x1p,x2p,prdkrn1,0,nkf)
            DO kr1=1,nkr
              g1(kr1)=g1p(kr1)
              g4(kr1)=g4p(kr1)
              g5(kr1)=g5p(kr1)
            END DO
            g1(nkr)=g1(nkr)+g1p(nkr+1)
            g4(nkr)=g4(nkr)+g4p(nkr+1)
            g5(nkr)=g5(nkr)+g5p(nkr+1)
          END IF

          ! ... for ice multiplication
          conc_old = 0.0_wp
          DO kr = kr_icempl,nkr
            conc_old = conc_old+col*g1(kr)/xl_mg(kr)
          END DO

          ! In case of high liquid water content, THEN: g4 (graupel) + g1 (water) --> g5 (hail):
          CALL kernel2d(cwxx,pp_r,"cwgl")
          IF ( usenkrp3 .EQ. 0 ) THEN
            CALL coll_xyz_lwf(g4,g1,g5,cwxx,xg_mg,xl_mg,prdkrn1,1,nkr)  ! 1 means do not allow collision of same
                                                                        ! bins since water+snow collided above !cwgl
          ELSE
            CALL coll_add_pseodo_3bin(x1p,x2p,g4p,g1p,g5p,nkf,xg_mg,xl_mg,g4,g1,g5)
            CALL coll_xyz_lwf(g4p,g1p,g5p,cwxx,x1p,x2p,prdkrn1,1,nkf)
            DO kr1=1,nkr
              g1(kr1)=g1p(kr1)
              g4(kr1)=g4p(kr1)
              g5(kr1)=g5p(kr1)
            END DO
            g1(nkr)=g1(nkr)+g1p(nkr+1)
            g4(nkr)=g4(nkr)+g4p(nkr+1)
            g5(nkr)=g5(nkr)+g5p(nkr+1)
          END IF
        END IF
      END IF

      IF (icol_hail == 1) THEN
        ! ---------------------------------------
        ! Collisions between drops and hail
        ! ---------------------------------------

        ! drops - hail = hail
        ! hail - water = hail
        CALL kernel2d(cwxx,pp_r,"cwlh")
        IF ( usenkrp2 .EQ. 0 ) THEN
          CALL coll_xyy_lwf(g1,g5,cwxx,xl_mg,xh_mg,prdkrn1,0,nkr) !cwlh
        ELSE
          CALL coll_add_pseodo_2bin(x1p,x2p,g1p,g5p,nkf,xl_mg,xh_mg,g1,g5)
          CALL coll_xyy_lwf(g1p,g5p,cwxx,x1p,x2p,prdkrn1,0,nkf)
          DO kr1=1,nkr
            g1(kr1)=g1p(kr1)
            g5(kr1)=g5p(kr1)
          END DO
          g1(nkr)=g1(nkr)+g1p(nkr+1)
          g5(nkr)=g5(nkr)+g5p(nkr+1)
        END IF

        ! ... for ice multiplication
        conc_old = 0.0_wp
        DO kr = kr_icempl,nkr
          conc_old = conc_old+col*g1(kr)/xl_mg(kr)
        END DO

        CALL kernel2d(cwxx,pp_r,"cwhl")
        IF ( usenkrp2 .EQ. 0 ) THEN
          CALL coll_xyx_lwf(g5,g1,cwxx,xh_mg,xl_mg,prdkrn1,1,nkr) ! 1 means do not allow collision of same
                                                                  ! bins since water+snow collided above !cwhl
        ELSE
          CALL coll_add_pseodo_2bin(x1p,x2p,g5p,g1p,nkf,xh_mg,xl_mg,g5,g1)
          CALL coll_xyx_lwf(g5p,g1p,cwxx,x1p,x2p,prdkrn1,1,nkf)
          DO kr1=1,nkr
            g1(kr1)=g1p(kr1)
            g5(kr1)=g5p(kr1)
          END DO
          g1(nkr)=g1(nkr)+g1p(nkr+1)
          g5(nkr)=g5(nkr)+g5p(nkr+1)
        END IF
      END IF

      ! Hallet Mossop ice-multiplication:
      ! small mass added to snow first bin due to water-graupel or water-hail collisions
      ! this mass is negligible and therefore is not subtracted from other PSD to conserve mass
      IF ((icol_graupel == 1 .OR. icol_hail == 1) .AND. icempl == 1) THEN
        IF (tt .GE. 265.15_wp .AND. tt .LE. 270.15_wp) THEN
          ! ... ice-multiplication (h-m) (linear interpolation using "triangle" with maximum at -5C):
          conc_new = 0.0
          DO kr = kr_icempl,nkr
            conc_new=conc_new+col*g1(kr)/xl_mg(kr)
          END DO
          IF (tt .LE. 268.15_wp) THEN
            conc_icempl=(conc_old-conc_new)*4.e-3_wp*(265.15_wp-tt)/(265.15_wp-268.15_wp)
          END IF
          IF (tt .GT. 268.15_wp) THEN
            conc_icempl=(conc_old-conc_new)*4.e-3_wp*(270.15_wp-tt)/(270.15_wp-268.15_wp)
          END IF
          !g2_2(1)=g2_2(1)+conc_icempl*xi2_mg(1)/col ! no crystals in fast-sbm
          g3(1)=g3(1)+conc_icempl*xs_mg(1)/col       ! and therefore they are added as small snow
        END IF
      END IF

      ! -------------------------------------
      ! Collisions which do not include drops
      ! -------------------------------------

      ! new addition: graupel+snow-->graupel or hail+snow-->hail:
      IF (icol_snow == 1) THEN
        IF ((hail_opt == 1) .AND. (icol_hail == 1)) THEN          ! g5 (hail) + g3 (snow) --> g5 (hail):
          CALL kernel2d(cwxx,pp_r,"cwhs")
          IF ( usenkrp2 .EQ. 0 ) THEN
            ! water is not involved, sticking efficiency is prdkrn=1 only for tc>0
            CALL coll_xyx_lwf(g5,g3,cwxx,xh_mg,xs_mg,prdkrn,0,nkr) !cwhs
          ELSE
            CALL coll_add_pseodo_2bin(x1p,x2p,g5p,g3p,nkf,xh_mg,xs_mg,g5,g3)
            CALL coll_xyx_lwf(g5p,g3p,cwxx,x1p,x2p,prdkrn,0,nkf)
            DO kr1=1,nkr
              g3(kr1)=g3p(kr1)
              g5(kr1)=g5p(kr1)
            END DO
            g3(nkr)=g3(nkr)+g3p(nkr+1)
            g5(nkr)=g5(nkr)+g5p(nkr+1)
          END IF
        ELSE IF ((hail_opt == 0) .AND. (icol_graupel == 1)) THEN  ! g4 (graupel) + g3 (snow) --> g4 (graupel):
          CALL kernel2d(cwxx,pp_r,"cwgs")
          IF ( usenkrp2 .EQ. 0 ) THEN
            ! water is not involved, sticking efficiency is prdkrn=1 only for tc>0
            CALL coll_xyx_lwf(g4,g3,cwxx,xg_mg,xs_mg,prdkrn,0,nkr) !cwgs
          ELSE
            CALL coll_add_pseodo_2bin(x1p,x2p,g4p,g3p,nkf,xg_mg,xs_mg,g4,g3)
            CALL coll_xyx_lwf(g4p,g3p,cwxx,x1p,x2p,prdkrn,0,nkf)
            DO kr1=1,nkr
              g3(kr1)=g3p(kr1)
              g4(kr1)=g4p(kr1)
            END DO
            g3(nkr)=g3(nkr)+g3p(nkr+1)
            g4(nkr)=g4(nkr)+g4p(nkr+1)
          END IF
        END IF
      END IF

      ! new addition: snow+graupel-->snow or snow+hail-->snow:
      IF (icol_snow == 1) THEN
        IF ((hail_opt == 1) .AND. (icol_hail == 1)) THEN          ! g3 (snow) + g5 (hail) --> g3 (snow) or g5 (hail):
          CALL kernel2d(cwxx,pp_r,"cwsh")
          IF ( usenkrp2 .EQ. 0 ) THEN
            ! water is not involved, sticking efficiency is prdkrn=1 only for tc>0
            ! 1 means do not allow collision of same bins since hail+snow collided above
            IF ( snha2ha == 1 ) THEN
              CALL coll_xyy_lwf(g3,g5,cwxx,xs_mg,xh_mg,prdkrn,1,nkr) !cwsh
            ELSE
              CALL coll_xyx_lwf(g3,g5,cwxx,xs_mg,xh_mg,prdkrn,1,nkr) !cwsh
            END IF
          ELSE
            CALL coll_add_pseodo_2bin(x1p,x2p,g3p,g5p,nkf,xs_mg,xh_mg,g3,g5)
            IF ( snha2ha == 1 ) THEN
              CALL coll_xyy_lwf(g3p,g5p,cwxx,x1p,x2p,prdkrn,1,nkf)
            ELSE
              CALL coll_xyx_lwf(g3p,g5p,cwxx,x1p,x2p,prdkrn,1,nkf)
            END IF
            DO kr1=1,nkr
              g3(kr1)=g3p(kr1)
              g5(kr1)=g5p(kr1)
            END DO
            g3(nkr)=g3(nkr)+g3p(nkr+1)
            g5(nkr)=g5(nkr)+g5p(nkr+1)
          END IF
        ELSE IF ((hail_opt == 0) .AND. (icol_graupel == 1)) THEN  ! g3 (snow) + g4 (graupel) --> g3 (snow) or g4 (graupel):
          CALL kernel2d(cwxx,pp_r,"cwsg")
          IF ( usenkrp2 .EQ. 0 ) THEN
            ! water is not involved, sticking efficiency is prdkrn=1 only for tc>0
            ! 1 means do not allow collision of same bins since hail+snow collided above
            IF ( snha2ha == 1 ) THEN
              CALL coll_xyy_lwf(g3,g4,cwxx,xs_mg,xg_mg,prdkrn,1,nkr) !cwsg
            ELSE
              CALL coll_xyx_lwf(g3,g4,cwxx,xs_mg,xg_mg,prdkrn,1,nkr) !cwsg
            END IF
          ELSE
            CALL coll_add_pseodo_2bin(x1p,x2p,g3p,g4p,nkf,xs_mg,xg_mg,g3,g4)
            IF ( snha2ha == 1 ) THEN
              CALL coll_xyy_lwf(g3p,g4p,cwxx,x1p,x2p,prdkrn,1,nkf)
            ELSE
              CALL coll_xyx_lwf(g3p,g4p,cwxx,x1p,x2p,prdkrn,1,nkf)
            END IF
            DO kr1=1,nkr
              g3(kr1)=g3p(kr1)
              g4(kr1)=g4p(kr1)
            END DO
            g3(nkr)=g3(nkr)+g3p(nkr+1)
            g4(nkr)=g4(nkr)+g4p(nkr+1)
          END IF
        END IF
      END IF

      ! --------------------------------
      ! Collisions between snowflakes
      ! --------------------------------

      IF (icol_snow == 1) THEN
        CALL kernel2d(cwxx,pp_r,"cwss")
        IF ( usenkrp1b .EQ. 0 ) THEN
          CALL coll_xxx_lwf(g3,cwxx,xs_mg,prdkrn,nkr) !cwss
        ELSE
          CALL coll_add_pseodo_bin(x1p,g3p,nkf,xs_mg,g3)
          CALL coll_xxx_lwf(g3p,cwxx,x1p,prdkrn,nkf)
          DO kr1=1,nkr
            g3(kr1)=g3p(kr1)
          END DO
          g3(nkr)=g3(nkr)+g3p(nkr+1)
        END IF
      END IF

      ! -------------------------------------------
      ! new addition: graupel+graupel or hail+hail:
      ! -------------------------------------------
      IF (icol_graupel == 1) THEN
        CALL kernel2d(cwxx,pp_r,"cwgg")
        IF ( usenkrp1b .EQ. 0 ) THEN
          CALL coll_xxx_lwf(g4,cwxx,xg_mg,prdkrn,nkr) !cwgg
        ELSE
          CALL coll_add_pseodo_bin(x1p,g4p,nkf,xg_mg,g4)
          CALL coll_xxx_lwf(g4p,cwxx,x1p,prdkrn,nkf)
          DO kr1=1,nkr
            g4(kr1)=g4p(kr1)
          END DO
          g4(nkr)=g4(nkr)+g4p(nkr+1)
        END IF
      END IF
      IF (icol_hail == 1) THEN
        CALL kernel2d(cwxx,pp_r,"cwhh")
        IF ( usenkrp1b .EQ. 0 ) THEN
          CALL coll_xxx_lwf(g5,cwxx,xh_mg,prdkrn,nkr) !cwhh
        ELSE
          CALL coll_add_pseodo_bin(x1p,g5p,nkf,xh_mg,g5)
          CALL coll_xxx_lwf(g5p,cwxx,x1p,prdkrn,nkf)
          DO kr1=1,nkr
            g5(kr1)=g5p(kr1)
          END DO
          g5(nkr)=g5(nkr)+g5p(nkr+1)
        END IF
      END IF
    END IF

    IF ( iceprocs==1 ) THEN
      ! latent heat release during freezing because of collisions:
      cont_fin_drop=0.
      cont_fin_ice=0.
      DO kr=1,nkr
        cont_fin_drop=cont_fin_drop+g1(kr)
        cont_fin_ice=cont_fin_ice+g3(kr)+g4(kr)+g5(kr)
        DO ice=1,icemax
          cont_fin_ice=cont_fin_ice+g2(kr,ice)
        END DO
      END DO
      cont_fin_drop=col*cont_fin_drop*1.e-3_wp
      cont_fin_ice=col*cont_fin_ice*1.e-3_wp
      deldrop=cont_init_drop-cont_fin_drop ! [g/cm**3] ! reduction of water mass due to collision freezing
      ! riming temperature correction (rho in g/cm**3) :
      IF (t_new <= tmelt) THEN
        IF (deldrop > 0.0_wp) THEN
          t_new = t_new + latheat_freez*deldrop/rho    ! latheat_freez takes into account cp/cv bug fix (was: t_new + 320.*deldrop/rho)
        ELSE ! IF deldrop < 0 ! problem: somehow the mass of water increased after collision freezing
          IF (ABS(deldrop).GT.cont_init_drop*0.05_wp) THEN
            CALL finish(TRIM(modname),"fatal error in module_mp_fast_sbm (ABS(deldrop).GT.cont_init_drop), model stop")
          END IF
        END IF
      END IF
    END IF

    ! recalculation of density FUNCTION f1,f3,f4,f5 in  units [1/(g*cm**3)] :
    DO kr=1,nkr
      ff1r(kr)=g1(kr)/(3.0_wp*xl(kr)*xl(kr)*1.e3_wp)
      IF ((ff1r(kr) .NE. ff1r(kr)) .OR. ff1r(kr) < 0.0_wp) THEN
        PRINT*,"g1",g1
        CALL finish(TRIM(modname),"stop at end coal_bott - ff1r nan or ff1r < 0.0_wp, model stop")
      END IF
      IF ( iceprocs==1 ) THEN
        ff3r(kr)=g3(kr)/(3.0_wp*xs(kr)*xs(kr)*1.e3)
        IF((ff3r(kr) .ne. ff3r(kr)) .OR. ff3r(kr) < 0.0_wp) THEN
          CALL finish(TRIM(modname),"stop at end coal_bott - ff3r nan or ff3r < 0.0_wp, model stop")
        END IF
        IF(hail_opt == 0) THEN
          ff4r(kr)=g4(kr)/(3.0_wp*xg(kr)*xg(kr)*1.e3_wp)+g5(kr)/(3.0_wp*xh(kr)*xh(kr)*1.e3_wp)
          ff5r(kr)=0.0_wp
          IF((ff4r(kr) .ne. ff4r(kr)) .OR. ff4r(kr) < 0.0_wp) THEN
            CALL finish(TRIM(modname),"stop at end coal_bott - ff4r nan or ff4r < 0.0_wp, model stop")
          END IF
        ELSE
          ff5r(kr)=g4(kr)/(3.0_wp*xg(kr)*xg(kr)*1.e3_wp)+g5(kr)/(3.0_wp*xh(kr)*xh(kr)*1.e3_wp)
          ff4r(kr)=0.0_wp
          IF((ff5r(kr) .ne. ff5r(kr)) .OR. ff5r(kr) < 0.0_wp) THEN
            CALL finish(TRIM(modname),"stop at end coal_bott - ff5r nan or ff5r < 0.0_wp, model stop")
          END IF
        END IF
      END IF
    END DO

    IF (ABS(tt-t_new).GT.5.0_wp) THEN
      CALL finish(TRIM(modname),"fatal error in module_mp_warm_sbm del_t 5 k, model stop")
    END IF

    tt = t_new

    RETURN
  END SUBROUTINE coallescence

  SUBROUTINE coll_xxx_bott(g,ckxx,x,prdkrn,nkf)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nkf
    REAL(KIND=wp),INTENT(INOUT) :: g(:)
    REAL(KIND=wp),INTENT(IN) :: ckxx(:,:),x(:), prdkrn
    REAL(KIND=wp):: gmin,x01,x02,x03,gsi,gsj,gsk,gk,flux,x1
    INTEGER :: i,ix0,ix1,j,k,kp,flag_i

    gmin = epsil

    ! ix0 - lower limit of integration by i
    DO i=1,nkf-1
      ix0=i
      IF (g(i).GT.gmin) EXIT
    END DO
    ! here ix0 is the smallest bin which is above gmin

    IF (ix0.EQ.nkf-1) RETURN

    ! ix1 - upper limit of integration by i
    DO i=nkf-1,1,-1
      ix1=i
      IF (g(i).GT.gmin) EXIT
    END DO
    ! here ix1 is the largest bin which is above gmin

    ! ... collisions
    DO i=ix0,ix1 !loop over bins with non negligible psd
      IF (g(i).LE.gmin) CYCLE !skip collisions of this negligible bin
      flag_i=0
      DO j=i,ix1 !loop over bins of j>=i
        IF (g(j).LE.gmin) CYCLE !skip collisions of this negligible bin
        !here we perform collisions of 2 non-negligible bins i,j where j>=i
        k=ima(i,j) !k and k+1 are the resulting bins from collision of i,j
        kp=k+1
        x01=ckxx(i,j)*g(i)*g(j)*prdkrn !prdkrn is sticking efficiency (input of this subroutine, 1 for water-water)
        x02=min(x01,g(i)*x(j))
        IF (j.NE.k) x03=min(x02,g(j)*x(i)) ! case when j=i and THEN k=j+1, kp=k+1
        IF (j.EQ.k) x03=x02                  ! case when j>i and THEN k=j, kp=k+1
        gsi=x03/x(j) ! mass subtracted from bin i
        gsj=x03/x(i) ! mass subtracted from bin j

        !make sure that we do not cause bin<0 by colliding it with others:
        IF ( gsi > g(i) ) THEN
          gsi=g(i)
          g(i)=0.0_wp
          flag_i=1
        ELSE
          g(i)=g(i)-gsi  ! resulting mass of bin i
        END IF

        ! In case j==i the k bin is >j and therefore negative g(j) is not permitted
        ! However, in case j>i, the k bin =j and THEN it will be refilled just afterwads in the "flux" remapping, so g(j)<0 is
        ! permitted
        IF ( ( j .EQ. i ) .AND. ( gsj > g(j) ) ) THEN
          gsj=g(j)
          g(j)=0.0_wp
        ELSE
          g(j)=g(j)-gsj ! resulting mass of bin j. Might be temporary < 0 before it is refilled after remapping
        END IF

        gsk=gsi+gsj  ! mass which is to be added to bin k
        gk=g(k)+gsk    ! when j=/k needs to be limited (only for different hydro)
        flux=0.0_wp

        ! g(i), g(j) - PSD f of bins i,j.
        ! their parts sum up (=gsk) have to be added to bin gk.
        ! since it falls between k and kp=k+1, it is added to bins g(k) and g(kp=k+1) via type of remapping (see flux below)

        IF (gk.GT.gmin) THEN
          x1=LOG(g(kp)/gk+gmin) ! avoid LOG(1) --> x1=0
          flux=gsk/x1*(EXP(0.5_wp*x1)-EXP(x1*(0.5_wp-chucm(i,j))))
          flux=min(flux,gk,gsk)
          g(k)=gk-flux
          g(kp)=g(kp)+flux
        END IF

        IF (g(i) < 0.0_wp .OR. g(j) < 0.0_wp .OR. g(k) < 0.0_wp .OR. g(kp) < 0.0_wp) THEN
          PRINT*, 'i,j,k,kp',i,j,k,kp,'ix0,ix1',ix0,ix1
          PRINT*, 'g(i),g(j),g(k),g(kp)'
          WRITE(*,'(A,4D13.5)') g(i),g(j),g(k),g(kp)
          stop 'stop in collisions coll_xxx_bott'
        END IF

        IF ( flag_i .EQ. 1 ) THEN !g(i) is now zero, no need to collide it with  other bins j>=i
          EXIT
        END IF
      END DO
    END DO

    RETURN
  END SUBROUTINE coll_xxx_bott

  SUBROUTINE coll_xxx_lwf(g,ckxx,x,prdkrn,nkf)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nkf
    REAL(KIND=wp),INTENT(INOUT) :: g(:)
    REAL(KIND=wp),INTENT(IN) :: ckxx(:,:),x(:), prdkrn
    REAL(KIND=wp):: gmin,x01,x02,x03,gsi,gsj,gsk,gk,flux,x1
    INTEGER :: i,ix0,ix1,j,k,kp,flag_i

    gmin = epsil

    ! ix0 - lower limit of integration by i
    DO i=1,nkf-1
      ix0=i
      IF (g(i).GT.gmin) EXIT
    END DO

    IF (ix0.EQ.nkf-1) RETURN

    ! ix1 - upper limit of integration by i
    DO i=nkf-1,1,-1
      ix1=i
      IF (g(i).GT.gmin) EXIT
    END DO

    ! ... collisions
    DO i=ix0,ix1
      IF (g(i).LE.gmin) CYCLE
      flag_i=0
      DO j=i,ix1
        IF (g(j).LE.gmin) CYCLE
        k=ima(i,j) ! target bin after collision, maximum value: nkf-1
        kp=k+1     ! maximum value: nkf
        x01=ckxx(i,j)*g(i)*g(j)*prdkrn
        x02=min(x01,g(i)*x(j))
        IF (j.NE.k) x03=min(x02,g(j)*x(i))
        IF (j.EQ.k) x03=x02
        gsi=x03/x(j)
        gsj=x03/x(i)

        !make sure that we do not cause bin<0 by colliding it with others:
        IF ( gsi > g(i) ) THEN
          gsi=g(i)
          g(i)=0.0_wp
          flag_i=1
        ELSE
          g(i)=g(i)-gsi  ! resulting mass of bin i
        END IF

        ! In case j==i the k bin is >j and therefore negative g(j) is not permitted
        ! However, in case j>i, the k bin =j and THEN it will be refilled just afterwads in the "flux" remapping, so g(j)<0 is
        ! permitted
        IF ( ( j .EQ. i ) .AND. ( gsj > g(j) ) ) THEN
          gsj=g(j)
          g(j)=0.0_wp
        ELSE
          g(j)=g(j)-gsj ! resulting mass of bin j. Might be temporary < 0 before it is refilled after remapping
        END IF

        gsk=gsi+gsj
        gk=g(k)+gsk
        IF ((gsk.LE.gmin) .OR. (gk.LE.gmin)) CYCLE

        flux=0.0_wp
        x1=LOG(g(kp)/gk+gmin)
        flux=gsk/x1*(EXP(0.5_wp*x1)-EXP(x1*(0.5_wp-chucm(i,j))))
        flux=min(flux,gsk)
        flux=min(flux,gk)
        IF (kp.GE.nkf) flux=0.5_wp*flux ! reduce transfer rate of mass from bin 32 to 33, ask Koby why?
        g(k)=gk-flux
        g(k)=max(g(k),gmin)
        g(kp)=g(kp)+flux
        g(kp)=max(g(kp),gmin)

        IF ( flag_i .EQ. 1 ) THEN !g(i) is now zero, no need to collide it with  other bins j>=i
          EXIT
        END IF
      END DO
    END DO

    RETURN
  END SUBROUTINE coll_xxx_lwf

  SUBROUTINE coll_xyy_lwf (gx,gy,ckxy,x,y,prdkrn,indc,nkf) !drop+hail->hail or drop+graupel->graupel
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nkf, indc
    REAL(KIND=wp),INTENT(INOUT) :: gy(:),gx(:)
    REAL(KIND=wp),INTENT(IN) :: ckxy(:,:),x(:),y(:),prdkrn
    REAL(KIND=wp) :: gmin,x01,x02,x03,gsi,gsj,gsk,gk,flux,x1
    INTEGER :: j,jx0,jx1,i,iy0,iy1,jmin,k,kp,flag_i

    gmin = epsil

    ! jx0 - lower limit of integration by j
    DO j=1,nkf-1
      jx0=j
      IF (gx(j).GT.gmin) EXIT
    END DO

    IF (jx0.eq.nkf-1) RETURN

    ! jx1 - upper limit of integration by j
    DO j=nkf-1,jx0,-1
      jx1=j
      IF (gx(j).GT.gmin) EXIT
    END DO

    ! iy0 - lower limit of integration by i
    DO i=1,nkf-1
      iy0=i
      IF (gy(i).GT.gmin) EXIT
    END DO

    IF (iy0.eq.nkf-1) RETURN

    ! iy1 - upper limit of integration by i
    DO i=nkf-1,iy0,-1
      iy1=i
      IF (gy(i).GT.gmin) EXIT
    END DO

    ! xyy --> running over small y(i) e.g. graupel
    ! and colliding with big x(j) e.g. water --> new y(k) e.g. graupel

    ! collisions :
    DO i = iy0,iy1
      IF (gy(i).LE.gmin) CYCLE
      flag_i=0
      jmin = i
      IF (jmin.eq.nkf-1) RETURN
      IF (i.LT.jx0) jmin=jx0-indc
      !indc=1: same particles with the same index cant collide, indc=0: valid for different particles
      !bottom line: the loop over j starts from i+1 or IF this bin is too small, from jx0
      DO j=jmin+indc,jx1 !indc=0 in this subroutine
        IF (gx(j).LE.gmin) CYCLE
        k=ima(i,j)
        kp=k+1
        x01=ckxy(j,i)*gy(i)*gx(j)*prdkrn
        x02=min(x01,gy(i)*x(j))
        x03=min(x02,gx(j)*y(i))
        gsi=x03/x(j)
        gsj=x03/y(i)

        !make sure that we do not cause bin<0 by colliding it with others:
        IF ( gsi > gy(i) ) THEN
          gsi=gy(i)
          gy(i)=0.0_wp
          flag_i=1
        ELSE
          gy(i)=gy(i)-gsi  ! resulting mass of bin i
        END IF

        IF ( gsj > gx(j) ) THEN
          gsj=gx(j)
          gx(j)=0.0_wp
        ELSE
          gx(j)=gx(j)-gsj  ! resulting mass of bin j
        END IF

        gsk=gsi+gsj
        gk=gy(k)+gsk !e.g. big graupel at bin k, which later goes into gy
        IF ((gsk.LE.gmin) .OR. (gk.LE.gmin)) CYCLE

        flux=0.0_wp

        x1=LOG(gy(kp)/gk+gmin)
        flux=gsk/x1*(EXP(0.5_wp*x1)-EXP(x1*(0.5_wp-chucm(i,j))))
        flux=min(flux,gsk)
        flux=min(flux,gk)

        IF (kp.GT.nkf) flux=0.5_wp*flux

        gy(k)=gk-flux
        gy(k)=max(gy(k),gmin)
        gy(kp)=gy(kp)+flux
        gy(kp)=max(gy(kp),gmin)

        IF ( flag_i .EQ. 1 ) THEN !gy(i) is now zero, no need to collide it with  other bins j>=i
          EXIT
        END IF
      END DO ! cycle by j
    END DO ! cycle by i

    RETURN
  END SUBROUTINE coll_xyy_lwf

  SUBROUTINE coll_xyx_lwf (gx,gy,ckxy,x,y,prdkrn,indc,nkf) !e.g. graupel+water-->graupel
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nkf, indc
    REAL(KIND=wp),INTENT(INOUT) :: gy(:),gx(:)
    REAL(KIND=wp),INTENT(IN) :: ckxy(:,:),x(:),y(:),prdkrn
    REAL(KIND=wp) :: gmin,x01,x02,x03,gsi,gsj,gsk,gk,flux,x1
    INTEGER :: j, jx0, jx1, i, iy0, iy1, jmin, k, kp, flag_i

    gmin = epsil

    ! jx0 - lower limit of integration by j
    DO j=1,nkf-1
      jx0=j
      IF (gx(j).GT.gmin) EXIT
    END DO

    IF (jx0.eq.nkf-1) RETURN

    ! jx1 - upper limit of integration by j
    DO j=nkf-1,jx0,-1
      jx1=j
      IF (gx(j).GT.gmin) EXIT
    END DO

    ! iy0 - lower limit of integration by i
    DO i=1,nkf-1
      iy0=i
      IF (gy(i).GT.gmin) EXIT
    END DO

    IF (iy0.eq.nkf-1) RETURN
    ! iy1 - upper limit of integration by i
    DO i=nkf-1,iy0,-1
      iy1=i
      IF (gy(i).GT.gmin) EXIT
    END DO

    ! xyx --> running over small y(i) e.g. graupel
    ! and colliding with big x(j) e.g. water --> new x(k) e.g. graupel

    ! ... collisions :
    DO i=iy0,iy1
      IF (gy(i).LE.gmin) CYCLE
      flag_i=0
      jmin=i
      IF (jmin.eq.nkf-1) RETURN
      IF (i.LT.jx0) jmin=jx0-indc
      !indc=1: same particles with the same index cant collide, indc=0: valid for different particles
      !bottom line: the loop over j starts from i+1 or IF this bin is too small, from jx0
      DO j=jmin+indc,jx1 !indc=1 in this suboutine
        IF (gx(j).LE.gmin) CYCLE
        k=ima(i,j) !k and k+1 are the resulting bins from collision of i,j
        kp=k+1
        x01=ckxy(j,i)*gy(i)*gx(j)*prdkrn
        x02=min(x01,gy(i)*x(j))
        IF (j.ne.k) x03=min(x02,gx(j)*y(i))
        IF (j.eq.k) x03=x02
        gsi=x03/x(j) ! mass subtracted from bin i
        gsj=x03/y(i) ! mass subtracted from bin j

        !make sure that we do not cause bin<0 by colliding it with others:
        IF ( gsi > gy(i) ) THEN
          gsi=gy(i)
          gy(i)=0.0_wp
          flag_i=1
        ELSE
          gy(i)=gy(i)-gsi  ! resulting mass of bin i
        END IF

        ! In case j==i the k bin is >j and therefore negative g(j) is not permitted
        ! However, in case j>i, the k bin =j and THEN it will be refilled just afterwads
        ! in the "flux" remapping, so g(j)<0 is permitted
        IF ( ( j .EQ. i ) .AND. ( gsj > gx(j) ) ) THEN
          gsj=gx(j)
          gx(j)=0.0_wp
        ELSE
          gx(j)=gx(j)-gsj  ! resulting mass of bin j
        END IF

        gsk=gsi+gsj
        gk=gx(k)+gsk
        IF ((gsk.LE.gmin) .OR. (gk.LE.gmin)) CYCLE

        flux=0.0_wp
        x1=LOG(gx(kp)/gk+gmin)
        flux=gsk/x1*(EXP(0.5_wp*x1)-EXP(x1*(0.5_wp-chucm(i,j))))
        flux=min(flux,gsk)
        flux=min(flux,gk)

        IF (kp.GT.nkf) flux=0.5_wp*flux
        gx(k)=gk-flux
        gx(k)=max(gx(k),gmin)

        gx(kp)=gx(kp)+flux
        gx(kp)=max(gx(kp),gmin)

        IF ( flag_i .EQ. 1 ) THEN !gy(i) is now zero, no need to collide it with  other bins j>=i
          EXIT
        END IF
      END DO ! cycle by j
    END DO ! cycle by i

    RETURN
  END SUBROUTINE coll_xyx_lwf

  SUBROUTINE coll_xyz_lwf(gx,gy,gz,ckxy,x,y,prdkrn,indc,nkf)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nkf, indc
    REAL(KIND=wp),INTENT(INOUT) :: gx(:),gy(:),gz(:)
    REAL(KIND=wp),INTENT(IN) :: ckxy(:,:),x(:),y(:),prdkrn
    REAL(KIND=wp) :: gmin,x01,x02,x03,gsi,gsj,gsk,gk,flux,x1
    INTEGER :: j,jx0,jx1,i,iy0,iy1,jmin,k,kp,flag_i

    gmin = epsil

    ! jx0 - lower limit of integration by j
    DO j=1,nkf-1
      jx0=j
      IF (gx(j) .GT. gmin) EXIT
    END DO
    ! here jx0 is the smallest bin which is above gmin

    IF (jx0 .EQ. nkf-1) RETURN

    ! jx1 - upper limit of integration by j
    DO j=nkf-1,jx0,-1
      jx1=j
      IF (gx(j) .GT. gmin) EXIT
    END DO
    ! here jx1 is the largest bin which is above gmin

    ! iy0 - lower limit of integration by i
    DO i=1,nkf-1
      iy0=i
      IF (gy(i) .GT. gmin) EXIT
    END DO

    IF (iy0.eq.nkf-1) RETURN

    ! iy1 - upper limit of integration by i
    DO i=nkf-1,iy0,-1
      iy1=i
      IF (gy(i) .GT. gmin) EXIT
    END DO

    ! ... collisions
    DO i=iy0,iy1
      IF (gy(i) .LE. gmin) CYCLE
      flag_i=0
      jmin=i
      IF (jmin .EQ. nkf-1) RETURN
      IF (i .LT. jx0) jmin=jx0-indc
      !indc=1: same particles with the same index cant collide, indc=0: valid for different particles
      !bottom line: the loop over j starts from i+1 or IF this bin is too small, from jx0
      !indc=0/1 here, 1 IF collision of these 2 same bin particles occured already in another subroutine
      DO j=jmin+indc,jx1
        IF (gx(j) .LE. gmin) CYCLE
        k=ima(i,j) !k and k+1 are the resulting bins from collision of i,j
        kp=k+1
        x01=ckxy(j,i)*gy(i)*gx(j)*prdkrn
        x02=min(x01,gy(i)*x(j))
        x03=min(x02,gx(j)*y(i))
        gsi=x03/x(j) ! mass subtracted from bin i
        gsj=x03/y(i) ! mass subtracted from bin j

        !make sure that we do not cause bin<0 by colliding it with others:
        IF ( gsi > gy(i) ) THEN
          gsi=gy(i)
          gy(i)=0.0_wp
          flag_i=1
        ELSE
          gy(i)=gy(i)-gsi  ! resulting mass of bin i
        END IF

        IF ( gsj > gx(j) ) THEN
          gsj=gx(j)
          gx(j)=0.0_wp
        ELSE
          gx(j)=gx(j)-gsj  ! resulting mass of bin j
        END IF

        gsk=gsi+gsj
        gk=gz(k)+gsk
        IF ((gsk.LE.gmin) .OR. (gk.LE.gmin)) CYCLE

        flux=0.0_wp

        x1=LOG(gz(kp)/gk+gmin)

        flux=gsk/x1*(EXP(0.5_wp*x1)-EXP(x1*(0.5_wp-chucm(i,j))))
        flux=min(flux,gsk)
        flux=min(flux,gk)

        IF (kp.GT.nkf) flux=0.5_wp*flux

        gz(k)=gk-flux
        gz(k)=max(gz(k),gmin)

        gz(kp)=gz(kp)+flux
        gz(kp)=max(gz(kp),gmin)

        IF ( flag_i .EQ. 1 ) THEN !gy(i) is now zero, no need to collide it with  other bins j>=i
          EXIT
        END IF
      END DO ! cycle by j
    END DO ! cycle by i

    RETURN
  END SUBROUTINE coll_xyz_lwf

  SUBROUTINE coll_breakup (gt, dt)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: dt
    REAL(KIND=wp),INTENT(INOUT) :: gt(:)
    INTEGER :: jdiff, k, i, j
    REAL(KIND=wp) :: xt(nkr+1), ft(nkr), fa(nkr), &
                         dg(nkr), df(nkr), dbreak(jbreak), &
                         amweight(jbreak), gain, aloss
    ! jbreak=18, maximum target bin number for splinter after collision
    ! gt : mass distribution function [g]
    ! xt : mass of bin [g]
    ! dt : timestep in s
    ! in cgs
    DO j=1,nkr
      xt(j)=xl(j)
      ft(j)=gt(j)/xt(j)/xt(j)
    END DO
    !shift between coagulation and breakup grid
    jdiff=nkr-jbreak

    !initialization
    !shift to breakup grid
    fa = 0.0_wp
    DO k=1,jbreak
      fa(k)=ft(k+jdiff)
    END DO

    !breakup: bleck's first order method
    !pkij: gain coefficients [g^3*cm^3/s]
    !qkj : loss coefficients

    xt(nkr+1)=xt(nkr)*2.0_wp

    amweight = 0.0_wp
    dbreak = 0.0_wp
    DO k=1,jbreak
      gain=0.0_wp
      DO i=1,jbreak
        DO j=1,i
          gain=gain+fa(i)*fa(j)*pkij(k,i,j)
        END DO
      END DO
      aloss=0.0_wp
      DO j=1,jbreak
        aloss=aloss+fa(j)*qkj(k,j)
      END DO
      j=nkr-jbreak+k
      amweight(k)=2.0_wp/(xt(j+1)**2.0_wp-xt(j)**2.0_wp)
      dbreak(k)=amweight(k)*(gain-fa(k)*aloss)

      IF (dbreak(k) .NE. dbreak(k)) THEN
        PRINT*,dbreak(k),amweight(k),gain,fa(k),aloss
        PRINT*,dbreak,amweight,j,nkr,jbreak,k,fa,xt,gt
        CALL finish(TRIM(modname)," inside coll_breakup, nan, model stop")
      END IF
    END DO

    !shift rate to coagulation grid
    df = 0.0_wp
    DO j=1,jdiff
      df(j)=0.0_wp
    END DO

    DO j=1,jbreak
      df(j+jdiff)=dbreak(j)
    END DO

    !transformation to mass distribution function g(ln x)
    DO j=1,nkr
      dg(j)=df(j)*xt(j)*xt(j)
    END DO

    !time integration
    DO j=1,nkr
      gt(j)=gt(j)+dg(j)*dt
    END DO

    RETURN
  END SUBROUTINE coll_breakup

  SUBROUTINE kernel2d(cwxx, pp_r, kernel_typ)
    IMPLICIT NONE
    REAL(KIND=wp), INTENT(INOUT) :: pp_r, cwxx(nkr,nkr)
    CHARACTER (len=4), INTENT(IN) :: kernel_typ
    INTEGER :: i,j,p_z_ind

    p_z_ind=(NINT(pp_r/p_z_del)*p_z_del-p_z_up)/p_z_del+1

    IF ( kernel_typ == "cwll" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwll_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwls" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwls_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwsl" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwsl_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwss" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwss_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwlg" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwlg_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwgl" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwgl_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwlh" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwlh_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwhl" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwhl_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwsg" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwsg_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwgs" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwgs_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwsh" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwsh_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwhs" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwhs_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwgg" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwgg_all(i,j,p_z_ind)
        END DO
      END DO
    ELSE IF ( kernel_typ == "cwhh" ) THEN
      DO i=1,nkr
        DO j=1,nkr
          cwxx(i,j)=cwhh_all(i,j,p_z_ind)
        END DO
      END DO
    END IF

    RETURN
  END SUBROUTINE kernel2d

  SUBROUTINE condevap_dmdt_coef (xmass,tp,pp,vxl,rxec,roxbl,bxx_my,icedim,in_,fl1)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: icedim, in_
    REAL(KIND=wp),INTENT(IN) :: xmass(nkr,icedim), tp, pp, vxl(nkr,icedim), &
                                rxec(nkr,icedim), roxbl(nkr,icedim), fl1(nkr)
    REAL(KIND=wp),INTENT(INOUT) :: bxx_my(nkr,icedim)
    INTEGER :: kr, nskin(nkr), ice
    REAL(KIND=wp) :: ventplm(nkr), fd1(nkr,icemax),fk1(nkr,icemax), &
              al1_my(2),esat1(2), d_my, coeff_viscous, &
              shmidt_number, a, b, rvt, reinolds_number, reshm, ventpl, constl

    DO kr=1,nkr
      IF (in_==2 .AND. fl1(kr)==0.0_wp .OR. in_==6 .OR. in_==3 .AND. tp<tmelt) THEN
        nskin(kr) = 2
      ELSE !in_==1 or in_==6 or lef/=0
        nskin(kr) = 1
      END IF
    END DO

    ! constants for clausius-clapeyron equation :
    al1_my(1)=alv_cgs ! latent heat of vaporization
    al1_my(2)=als_cgs ! latent heat of sublimation
    d_my=zdv_cgs*(pzero/pp)*(tp/tmelt)**1.94_wp !coefficient of diffusion
    coeff_viscous=0.13_wp !coefficient of viscousity (cm*cm/sec)
    shmidt_number=coeff_viscous/d_my !shmidt number

    ! constants used for calculation of reinolds number
    a=2.0_wp*(3.0_wp/4.0_wp/pi)**(1.0_wp/3.0_wp)
    b=a/coeff_viscous

    rvt=rv_cgs*tp
    ! update the saturation vapor pressure
    esat1(1) = 10.0_wp*sat_pres_water(tp) !dynes/cm^2
    esat1(2) = 10.0_wp*sat_pres_ice(tp)   !dynes/cm^2
    DO kr=1,nkr
      ventplm(kr)=0.0_wp
    END DO
    DO ice=1,icedim
      DO kr=1,nkr
        ! reynolds numbers
        reinolds_number=b*vxl(kr,ice)*(xmass(kr,ice)/roxbl(kr,ice))**(1.0_wp/3.0_wp)
        reshm=SQRT(reinolds_number)*(shmidt_number**(1.0_wp/3.0_wp)) !Pruppacher and Klett 1997, p.212
        IF (reinolds_number<2.5_wp) THEN
          ventpl=1.0_wp+0.108_wp*reshm*reshm
          ventplm(kr)=ventpl
        ELSE
          ventpl=0.78_wp+0.308_wp*reshm
          ventplm(kr)=ventpl
        END IF
      END DO

      ! ventpl_max is given in micro.prm include file
      DO kr=1,nkr
        ventpl=ventplm(kr)
        IF (ventpl>ventpl_max) THEN
          ventpl=ventpl_max
          ventplm(kr)=ventpl
        END IF
        constl=const_dmdt*rxec(kr,ice)
        fd1(kr,ice)=rvt/d_my/esat1(nskin(kr))
        IF (.not. latheatfac2) THEN !no signficant effect, related to cp/cv bug
          fk1(kr,ice)=(al1_my(nskin(kr))/rvt-1.0_wp)*al1_my(nskin(kr))/cf_my/tp
        ELSE
          fk1(kr,ice)=(al1_my(nskin(kr))*(cpd/cvd)/rvt-1.0_wp)*al1_my(nskin(kr))*(cpd/cvd)/cf_my/tp
        END IF

        ! growth rate
        bxx_my(kr,ice)=ventpl*constl/(fk1(kr,ice)+fd1(kr,ice))
      END DO
    END DO

    RETURN
  END SUBROUTINE condevap_dmdt_coef

  SUBROUTINE condevap_dsdt_coef (fi1,x1,sfn11,bxx_my,cf,icedim)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: icedim
    REAL(KIND=wp),INTENT(IN) :: bxx_my(nkr,icedim), fi1(nkr,icedim), cf, x1(nkr,icedim)
    REAL(KIND=wp),INTENT(OUT) :: sfn11(icedim)
    INTEGER :: ice, kr
    REAL(KIND=wp) :: sfn11s, delm, fun, b11

    DO ice=1,icedim
      sfn11s=0.0_wp
      sfn11(ice)=cf*sfn11s
      DO kr=1,nkr
        ! delta-m
        delm=x1(kr,ice)*3.0_wp*col
        ! integral's expression
        fun=fi1(kr,ice)*delm
        ! values of integrals
        b11=bxx_my(kr,ice)
        sfn11s=sfn11s+fun*b11
      END DO
      ! cycle by kr
      ! correction
      sfn11(ice)=cf*sfn11s
    END DO
    ! cycle by ice
    RETURN
  END SUBROUTINE condevap_dsdt_coef

  SUBROUTINE condevap_supsat_eqn (del1,del2,del1n,del2n,rw,pw,ri,pi_, &
                            dt,del1int,del2int,dyn1,dyn2,  &
                            isym1,isym2,isym3,isym4,isym5)
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: isym1, isym2(:), isym3, isym4, isym5
    REAL(KIND=wp),INTENT(IN) :: dt, dyn1, dyn2, del1, del2
    REAL(KIND=wp),INTENT(INOUT) :: del1n,del2n,del1int,del2int,rw, pw, ri, pi_
    INTEGER ::  isymice, irw, iri
    REAL(KIND=wp) :: x, expm1, deter, expr, expp, a, alfa, beta, gama, g31, g32, g2, expb, expg, &
              c11, c21, c12, c22, a1del1n, a2del1n, a3del1n, a4del1n, a1del1int, a2del1int, &
              a3del1int, a4del1int, a1del2n, a2del2n, a3del2n , a4del2n, a1del2int, a2del2int, &
              a3del2int, a4del2int, a5del2int

    expm1(x)=x+x*x/2.0_wp+x*x*x/6.0_wp+x*x*x*x/24.0_wp+x*x*x*x*x/120.0_wp

    isymice = sum(isym2) + isym3 + isym4 + isym5
    irw = 1
    iri = 1
    IF (MAX(rw,pw,ri,pi_)<=rw_pw_ri_pi_min) THEN
      rw = 0.0_wp
      irw = 0
      pw = 0.0_wp
      ri = 0.0_wp
      iri = 0
      pi_ = 0.0_wp
      isym1 = 0
      isymice = 0
    ELSE
      IF (max(rw,pw)>rw_pw_min) THEN
        ! a zero can pass through, assign a minimum value
        IF (rw < rw_pw_min*rw_pw_min) THEN
          rw = 1.0e-20_wp
          irw = 0
        END IF
        IF (pw < rw_pw_min*rw_pw_min)THEN
          pw = 1.0e-20_wp
        END IF
        IF (max(pi_/pw,ri/rw)<=ratio_icew_min) THEN
          ! ... only water
          ri = 0.0_wp
          iri = 0
          pi_ = 0.0_wp
          isymice = 0
        END IF
        IF (min(pi_/pw,ri/rw)>1.0_wp/ratio_icew_min) THEN
          ! ... only ice
          rw = 0.0_wp
          irw = 0
          pw = 0.0_wp
          isym1 = 0
        END IF
      ELSE
        ! only ice
        rw = 0.0_wp
        irw = 0
        pw = 0.0_wp
        isym1 = 0
      END IF
    END IF
    IF (isymice == 0) THEN
      isym2 = 0
      isym3 = 0
      isym4 = 0
      isym5 = 0
    END IF
    deter=rw*pi_-pw*ri
    IF (irw == 0 .AND. iri == 0) THEN
      del1n=del1+dyn1*dt
      del2n=del2+dyn2*dt
      del1int=del1*dt+dyn1*dt*dt/2.0_wp
      del2int=del2*dt+dyn2*dt*dt/2.0_wp
      RETURN
    END IF
    ! solution of equation for supersaturation with different deter values
    IF (iri == 0) THEN
      ! ... only water                                       (start)
      expr=EXP(-rw*dt)
      IF (ABS(rw*dt)>1.0e-6_wp) THEN
        del1n=del1*expr+(dyn1/rw)*(1.0_wp-expr)
        del2n=pw*del1*expr/rw-pw*dyn1*dt/rw-pw*dyn1*expr/(rw*rw)+dyn2*dt+ &
              del2-pw*del1/rw+pw*dyn1/(rw*rw)
        del1int=-del1*expr/rw+dyn1*dt/rw+dyn1*expr/(rw*rw)+del1/rw-dyn1/(rw*rw)
        del2int=pw*del1*expr/(-rw*rw)-pw*dyn1*dt*dt/(2.0_wp*rw)+ &
                pw*dyn1*expr/(rw*rw*rw)+dyn2*dt*dt/2.0_wp+ &
                del2*dt-pw*del1*dt/rw+pw*dyn1*dt/(rw*rw)+ &
                pw*del1/(rw*rw)-pw*dyn1/(rw*rw*rw)
        RETURN
        ! in case abs(rw*dt)>1.0e-6
      ELSE
        ! in case abs(rw*dt)<=1.0e-6
        expr=expm1(-rw*dt)
        del1n=del1+del1*expr+(dyn1/rw)*(0.0_wp-expr)
        del2n=pw*del1*expr/rw-pw*dyn1*dt/rw-pw*dyn1*expr/(rw*rw)+dyn2*dt+del2
        del1int=-del1*expr/rw+dyn1*dt/rw+dyn1*expr/(rw*rw)
        del2int=pw*del1*expr/(-rw*rw)-pw*dyn1*dt*dt/(2.0_wp*rw)+ &
              pw*dyn1*expr/(rw*rw*rw)+dyn2*dt*dt/2.0_wp+ &
              del2*dt-pw*del1*dt/rw+pw*dyn1*dt/(rw*rw)
        RETURN
      END IF
       ! ... only water                                                    (end)
       ! in case ri==0.0_wp
    END IF
    IF (irw  ==  0) THEN
      ! ... only ice                                                    (start)
      expp=EXP(-pi_*dt)
      IF (ABS(pi_*dt)>1.0e-6_wp) THEN
        del2n = del2*expp+(dyn2/pi_)*(1.0_wp-expp)
        del2int = -del2*expp/pi_+dyn2*dt/pi_+dyn2*expp/(pi_*pi_)+del2/pi_-dyn2/(pi_*pi_)
        del1n = +ri*del2*expp/pi_-ri*dyn2*dt/pi_-ri*dyn2*expp/(pi_*pi_)+dyn1*dt+ &
                 del1-ri*del2/pi_+ri*dyn2/(pi_*pi_)
        del1int = -ri*del2*expp/(pi_*pi_)-ri*dyn2*dt*dt/(2.0_wp*pi_)+ &
                 ri*dyn2*expp/(pi_*pi_*pi_)+dyn1*dt*dt/2.0_wp+ &
                 del1*dt-ri*del2*dt/pi_+ri*dyn2*dt/(pi_*pi_)+ &
                 ri*del2/(pi_*pi_)-ri*dyn2/(pi_*pi_*pi_)
        RETURN
        ! in case abs(pi_*dt)>1.0e-6
      ELSE
        ! in case abs(pi_*dt)<=1.0e-6
        expp=expm1(-pi_*dt)
        del2n=del2+del2*expp-expp*dyn2/pi_
        del2int=-del2*expp/pi_+dyn2*dt/pi_+dyn2*expp/(pi_*pi_)
        del1n=+ri*del2*expp/pi_-ri*dyn2*dt/pi_-ri*dyn2*expp/(pi_*pi_)+dyn1*dt+del1
        del1int=-ri*del2*expp/(pi_*pi_)-ri*dyn2*dt*dt/(2.0_wp*pi_)+ &
                     ri*dyn2*expp/(pi_*pi_*pi_)+dyn1*dt*dt/2.0_wp+ &
                     del1*dt-ri*del2*dt/pi_+ri*dyn2*dt/(pi_*pi_)
        RETURN
      END IF
      ! ... only ice                                                      (end)
      ! in case rw==0.0_wp
    END IF
    IF (irw == 1 .AND. iri == 1) THEN
      a=(rw-pi_)*(rw-pi_)+4.0_wp*pw*ri
      IF (a < 0.0_wp) THEN
        PRINT*,   'in SUBROUTINE jersupsat: a < 0'
        WRITE (*,'(A,D13.5)')  'deter',        deter
        WRITE (*,'(A,4D13.5)') 'rw,pw,ri,pi_',  rw,pw,ri,pi_
        WRITE (*,'(A,3D13.5)') 'dt,dyn1,dyn2', dt,dyn1,dyn2
        WRITE (*,'(A,2D13.5)') 'del1,del2',    del1,del2
        PRINT*,   'stop 1905:a < 0'
        CALL finish(TRIM(modname),"fatal error: stop 1905:a < 0, model stop")
      END IF
      ! ... water and ice                                               (start)
      alfa=SQRT((rw-pi_)*(rw-pi_)+4.0_wp*pw*ri)
      ! beta is negative to the simple solution so it will decay
      beta=0.5_wp*(alfa+rw+pi_)
      gama=0.5_wp*(alfa-rw-pi_)
      g31=pi_*dyn1-ri*dyn2
      g32=-pw*dyn1+rw*dyn2
      g2=rw*pi_-ri*pw
      IF (g2 < 1.0e-20_wp) g2 = 1.0004e-11_wp*1.0003e-11_wp-1.0002e-11_wp*1.0001e-11_wp
      expb=EXP(-beta*dt)
      expg=EXP(gama*dt)

      IF (ABS(gama*dt)>1.0e-6_wp) THEN
        c11=(beta*del1-rw*del1-ri*del2-beta*g31/g2+dyn1)/alfa
        c21=(gama*del1+rw*del1+ri*del2-gama*g31/g2-dyn1)/alfa
        c12=(beta*del2-pw*del1-pi_*del2-beta*g32/g2+dyn2)/alfa
        c22=(gama*del2+pw*del1+pi_*del2-gama*g32/g2-dyn2)/alfa
        del1n=c11*expg+c21*expb+g31/g2
        del1int=c11*expg/gama-c21*expb/beta+(c21/beta-c11/gama)+g31*dt/g2
        del2n=c12*expg+c22*expb+g32/g2
        del2int=c12*expg/gama-c22*expb/beta+(c22/beta-c12/gama)+g32*dt/g2
        RETURN
        ! in case abs(gama*dt)>1.0e-6
      ELSE
        ! in case abs(gama*dt)<=1.0e-6
        IF (ABS(ri/rw)>epsil) THEN
          IF (ABS(rw/ri)>epsil) THEN
            alfa=SQRT((rw-pi_)*(rw-pi_)+4.0_wp*pw*ri)
            beta=0.5_wp*(alfa+rw+pi_)
            gama=0.5_wp*(alfa-rw-pi_)
            IF (gama < 0.5_wp*2.0e-10_wp) gama=0.5_wp*(2.002e-10_wp-2.001e-10_wp)
            expg=expm1(gama*dt)
            expb=EXP(-beta*dt)
            ! beta/alfa could be very close to 1 that why i transform it
            ! remember alfa-beta=gama
            c11=(beta*del1-rw*del1-ri*del2+dyn1)/alfa
            c21=(gama*del1+rw*del1+ri*del2-gama*g31/g2-dyn1)/alfa
            c12=(beta*del2-pw*del1-pi_*del2+dyn2)/alfa
            c22=(gama*del2+pw*del1+pi_*del2-gama*g32/g2-dyn2)/alfa

            a1del1n=c11
            a2del1n=c11*expg
            a3del1n=c21*expb
            a4del1n=g31/g2*(gama/alfa+(gama/alfa-1.0_wp)*expg)

            del1n=a1del1n+a2del1n+a3del1n+a4del1n

            a1del1int=c11*expg/gama
            a2del1int=-c21*expb/beta
            a3del1int=c21/beta
            a4del1int=g31/g2*dt*(gama/alfa)

            del1int=a1del1int+a2del1int+a3del1int+a4del1int

            a1del2n=c12
            a2del2n=c12*expg
            a3del2n=c22*expb
            a4del2n=g32/g2*(gama/alfa+ &
                       (gama/alfa-1.0_wp)* &
                       (gama*dt+gama*gama*dt*dt/2.0_wp))

            del2n=a1del2n+a2del2n+a3del2n+a4del2n

            a1del2int=c12*expg/gama
            a2del2int=-c22*expb/beta
            a3del2int=c22/beta
            a4del2int=g32/g2*dt*(gama/alfa)
            a5del2int=g32/g2*(gama/alfa-1.0_wp)*(gama*dt*dt/2.0_wp)

            del2int=a1del2int+a2del2int+a3del2int+a4del2int+a5del2int
            ! in case abs(rw/ri)>1e-12
          ELSE
            ! in case abs(rw/ri)<=1e-12
            x=-2.0_wp*rw*pi_+rw*rw+4.0_wp*pw*ri

            alfa=pi_*(1+(x/pi_)/2.0_wp-(x/pi_)*(x/pi_)/8.0_wp)
            beta=pi_+(x/pi_)/4.0_wp-(x/pi_)*(x/pi_)/16.0_wp+rw/2.0_wp
            gama=(x/pi_)/4.0_wp-(x/pi_)*(x/pi_)/16.0_wp-rw/2.0_wp

            expg=expm1(gama*dt)
            expb=EXP(-beta*dt)

            c11=(beta*del1-rw*del1-ri*del2+dyn1)/alfa
            c21=(gama*del1+rw*del1+ri*del2-gama*g31/g2-dyn1)/alfa
            c12=(beta*del2-pw*del1-pi_*del2+dyn2)/alfa
            c22=(gama*del2+pw*del1+pi_*del2-gama*g32/g2-dyn2)/alfa

            del1n=c11+c11*expg+c21*expb+g31/g2*(gama/alfa+(gama/alfa-1)*expg)
            del1int=c11*expg/gama-c21*expb/beta+(c21/beta)+g31/g2*dt*(gama/alfa)
            del2n=c12+c12*expg+c22*expb+g32/g2*(gama/alfa+(gama/alfa-1.0_wp)* &
                    (gama*dt+gama*gama*dt*dt/2.0_wp))
            del2int=c12*expg/gama-c22*expb/beta+(c22/beta)+g32/g2*dt*(gama/alfa)+ &
                    g32/g2*(gama/alfa-1.0_wp)*(gama*dt*dt/2.0_wp)
            ! in case abs(rw/ri)<=1e-12
          END IF
          ! alfa/beta 2
          ! in case abs(ri/rw)>1e-12
        ELSE
          ! in case abs(ri/rw)<=1e-12
          x=-2.0_wp*rw*pi_+pi_*pi_+4.0_wp*pw*ri

          alfa=rw*(1.0_wp+(x/rw)/2.0_wp-(x/rw)*(x/rw)/8.0_wp)
          beta=rw+(x/rw)/4.0_wp-(x/rw)*(x/rw)/16.0_wp+pi_/2.0_wp
          gama=(x/rw)/4.0_wp-(x/rw)*(x/rw)/16.0_wp-pi_/2.0_wp

          expg=expm1(gama*dt)
          expb=EXP(-beta*dt)

          c11=(beta*del1-rw*del1-ri*del2+dyn1)/alfa
          c21=(gama*del1+rw*del1+ri*del2-gama*g31/g2-dyn1)/alfa
          c12=(beta*del2-pw*del1-pi_*del2+dyn2)/alfa
          c22=(gama*del2+pw*del1+pi_*del2-gama*g32/g2-dyn2)/alfa

          del1n=c11+c11*expg+c21*expb+g31/g2*(gama/alfa+(gama/alfa-1.0_wp)*expg)
          del1int=c11*expg/gama-c21*expb/beta+(c21/beta)+g31/g2*dt*(gama/alfa)
          del2n=c12+c12*expg+c22*expb+g32/g2*(gama/alfa+ &
                   (gama/alfa-1.0_wp)*(gama*dt+gama*gama*dt*dt/2.0_wp))
          del2int=c12*expg/gama-c22*expb/beta+c22/beta+g32/g2*dt*(gama/alfa)+ &
                   g32/g2*(gama/alfa-1.0_wp)*(gama*dt*dt/2.0_wp)
          ! alfa/beta
          ! in case abs(ri/rw)<=1e-12
        END IF
        ! in case abs(gama*dt)<=1e-6
      END IF
      ! water and ice                                                 (end)
      ! in case isym1/=0.AND.isym2/=0
    END IF

    RETURN
  END SUBROUTINE condevap_supsat_eqn

  SUBROUTINE condevap_mass_eqn_and_remap1 (xmass,b21_my,fi2,psi2,del2n, &
                                  isym2,ind,itype,idrop,iin,kin)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: isym2, ind, itype, iin, kin
    INTEGER,INTENT(INOUT) :: idrop
    REAL(KIND=wp),INTENT(IN) :: b21_my(:), fi2(:), del2n, xmass(:)
    REAL(KIND=wp),INTENT(INOUT) :: psi2(:)
    INTEGER :: kr, nr
    REAL(KIND=wp) :: xin(nkr),fi2r(nkr),psi2r(nkr),d,ratexi,b,a,xir(nkr),xinr(nkr)

    DO kr=1,nkr
      psi2r(kr) = fi2(kr)
      fi2r(kr) = fi2(kr)
    END DO
    nr=nkr
    ! new size distribution functions                             (start)
    IF (isym2 == 1) THEN
      IF (ind==1 .AND. itype==1) THEN
        ! drop diffusional growth
        DO kr=1,nkr
          d=xmass(kr)**(1.0_wp/3.0_wp)
          ratexi=(2.0_wp/3.0_wp)*del2n*b21_my(kr)/d
          b=xmass(kr)**(2.0_wp/3.0_wp)
          a=b+ratexi                !mass change due to diffusional growth/evaporation of water
          IF (a<0.0_wp) THEN
            xin(kr)=1.0e-50_wp
          ELSE
            xin(kr)=a**(3.0_wp/2.0_wp)
          END IF
        END DO
        ! in case ind==1.AND.itype==1
      ELSE
        ! in case ind/=1.OR.itype/=1
        DO kr=1,nkr
          ratexi = del2n*b21_my(kr)
          xin(kr) = xmass(kr) + ratexi !mass change due to diffusional growth/evaporation of ice
        END DO
      END IF

      ! recalculation of size distribution FUNCTIONs                (start)
      DO kr=1,nkr
        xir(kr) = xmass(kr)
        xinr(kr) = xin(kr)
        fi2r(kr) = fi2(kr)
      END DO

      CALL remapping(nr,xir,fi2r,psi2r,xinr,idrop,iin,kin)

      DO kr=1,nkr
        IF (psi2r(kr)<0.0_wp) THEN
          PRINT*,    'stop 1506 : psi2r(kr)<0.0_wp, in condevap_mass_eqn_and_remap1'
          CALL finish(TRIM(modname),"fatal error in psi2r(kr)<0.0_wp, in condevap_mass_eqn_and_remap1, model stop")
        END IF
        psi2(kr) = psi2r(kr)
      END DO
      ! cycle by ice
      ! recalculation of size distribution FUNCTIONs                  (end)
      ! in case isym2/=0
    END IF
    ! new size distribution FUNCTIONs                               (end)

    RETURN
  END SUBROUTINE condevap_mass_eqn_and_remap1

  SUBROUTINE condevap_mass_eqn_and_remap2(xmass,b21_my,fi2,psi2,idrop,iin,kin)
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: idrop
    INTEGER,INTENT(IN) :: iin,kin
    REAL(KIND=wp),INTENT(IN) :: xmass(:), fi2(:), b21_my(:)
    REAL(KIND=wp),INTENT(INOUT) :: psi2(:)
    INTEGER :: nr, kr
    REAL(KIND=wp) :: d, ratexi, b, a, xir(nkr),fi2r(nkr),psi2r(nkr),xinr(nkr)

    nr=nkr
    xir = xmass
    fi2r = fi2
    psi2r = psi2

    ! new drop size distribution functions:

    ! drop diffusional growth
    DO kr=1,nkr
      d = xir(kr)**(1.0_wp/3.0_wp)
      ratexi = (2.0_wp/3.0_wp)*b21_my(kr)/d
      b = xir(kr)**(2.0_wp/3.0_wp)
      a = b+ratexi
      IF (a<0.0_wp) THEN
        xinr(kr) = 1.0e-50_wp
      ELSE
        xinr(kr) = a**(3.0_wp/2.0_wp)
      END IF
    END DO

    ! recalculation of size distribution functions:
    !calculate the new PSD (remapping):
    CALL remapping(nr,xir,fi2r,psi2r,xinr,idrop,iin,kin)

    psi2 = psi2r

    RETURN
  END SUBROUTINE condevap_mass_eqn_and_remap2

  SUBROUTINE remapping(nrx,rr,fi,psi,rn,idrop,iin,kin)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nrx, iin, kin
    INTEGER,INTENT(INOUT) :: idrop
    REAL(KIND=wp),INTENT(INOUT) :: psi(:), rn(:), fi(:), rr(:)
    INTEGER :: kmax, kr, i, k , nrxp, isign_diffusional_growth, nrx1, i3point_condevap
    REAL(KIND=wp) :: rntmp,rrtmp,rrp,rrm,rntmp2,rrtmp2,rrp2,rrm2, gn1,gn2, &
             gn3,gn1p,gmat,gmat2,cdrop(nrx),delta_cdrop(nrx),rrs(nrx+1),psinew(nrx+1), &
             psi_im,psi_i,psi_ip
    INTEGER,PARAMETER :: krdrop_remaping_min = 6, krdrop_remaping_max = 12

    nrxp = nrx + 1
    nrx1 = nkr

    DO i=1,nrx
      ! rn(i), g - new masses after condensation or evaporation
      IF (rn(i) < 0.0_wp) THEN
        rn(i) = 1.0e-50_wp
        fi(i) = 0.0_wp
      END IF
    END DO

    DO k=1,nrx
      rrs(k)=rr(k)
    END DO

    i3point_condevap = isign_3point

    IF (rn(1) < rrs(1)) THEN
      ! evaporation
      i3point_condevap = 0
      idrop = 0
      nrx1 = nrx
    END IF

    IF (idrop == 0) i3point_condevap = 0

    DO k=1,nrx
      psi(k)=0.0_wp
      cdrop(k)=0.0_wp
      delta_cdrop(k)=0.0_wp
      psinew(k)=0.0_wp
    END DO
    rrs(nrxp)=rrs(nrx)*1024.0_wp
    psinew(nrxp) = 0.0_wp

    isign_diffusional_growth = 0
    DO k=1,nrx
      IF (rn(k).NE.rr(k)) THEN
        isign_diffusional_growth = 1
        EXIT
      END IF
    END DO

    IF (isign_diffusional_growth == 1) THEN
      ! kovetz-olund method                                         (start)
      DO k=1,nrx1 ! nrx1-1
        IF (fi(k) > 0.0_wp) THEN
          IF (ABS(rn(k)-rr(k)) < 1.0e-16_wp) THEN
            psinew(k) = fi(k)*rr(k)
            cycle
          END IF

          i = 1
          DO WHILE (.not.(rrs(i) <= rn(k) .AND. rrs(i+1) >= rn(k)) .AND. i .LT. nrx1) ! was nrx1-1
            i = i + 1
          END DO

          IF (rn(k).LT.rrs(1)) THEN
            rntmp=rn(k)
            rrtmp=0.0_wp
            rrp=rrs(1)
            gmat2=(rntmp-rrtmp)/(rrp-rrtmp)
            psinew(1)=psinew(1)+fi(k)*rr(k)*gmat2
          ELSE
            rntmp=rn(k)
            rrtmp=rrs(i)
            rrp=rrs(i+1)
            gmat2=(rntmp-rrtmp)/(rrp-rrtmp)
            gmat=(rrp-rntmp)/(rrp-rrtmp)
            psinew(i)=psinew(i)+fi(k)*rr(k)*gmat
            psinew(i+1)=psinew(i+1)+fi(k)*rr(k)*gmat2
         END IF
        END IF
      END DO

      DO kr=1,nrx1
        psi(kr)=psinew(kr)
      END DO

      DO kr=nrx1+1,nrx
        psi(kr)=fi(kr)
      END DO
      ! kovetz-olund method                                 (end)

      ! calculation both new total drop concentrations(after ko) and new total drop masses (after ko)

      ! 3point method	                                         (start)
      IF (i3point_condevap == 1) THEN
        DO k=1,nrx1-1
          IF (fi(k) > 0.0_wp) THEN
            IF (ABS(rn(k)-rr(k)).LT.1.0e-16_wp) THEN
              psi(k) = fi(k)*rr(k)
              CYCLE
            END IF

            IF (rrs(2).LT.rn(k)) THEN
              i = 2
              DO WHILE (.not.(rrs(i) <= rn(k) .AND. rrs(i+1) >= rn(k)) .AND. i.LT.nrx1-1)
                i=i+1
              END DO
              rntmp=rn(k)
              rrtmp=rrs(i)
              rrp=rrs(i+1)
              rrm=rrs(i-1)
              rntmp2=rn(k+1)
              rrtmp2=rrs(i+1)
              rrp2=rrs(i+2)
              rrm2=rrs(i)
              gn1=(rrp-rntmp)*(rrtmp-rntmp)/(rrp-rrm)/(rrtmp-rrm)
              gn1p=(rrp2-rntmp2)*(rrtmp2-rntmp2)/(rrp2-rrm2)/(rrtmp2-rrm2)
              gn2=(rrp-rntmp)*(rntmp-rrm)/(rrp-rrtmp)/(rrtmp-rrm)
              gmat=(rrp-rntmp)/(rrp-rrtmp)
              gn3=(rrtmp-rntmp)*(rrm-rntmp)/(rrp-rrm)/(rrp-rrtmp)
              gmat2=(rntmp-rrtmp)/(rrp-rrtmp)
              psi_im = psi(i-1)+gn1*fi(k)*rr(k)
              psi_i = psi(i)+gn1p*fi(k+1)*rr(k+1)+(gn2-gmat)*fi(k)*rr(k)
              psi_ip = psi(i+1)+(gn3-gmat2)*fi(k)*rr(k)
              IF (psi_im > 0.0_wp) THEN
                IF (psi_ip > 0.0_wp) THEN
                  IF (i > 2) THEN
                    ! smoothing criteria
                    IF (psi_im > psi(i-2) .AND. psi_im < psi_i &
                      .AND. psi(i-2) < psi(i) .OR. psi(i-2) >= psi(i)) THEN
                      psi(i-1) = psi_im
                      psi(i) = psi(i) + fi(k)*rr(k)*(gn2-gmat)
                      psi(i+1) = psi_ip
                    END IF
                  END IF
                ELSE
                  EXIT
                END IF
              ELSE
                EXIT
              END IF
            END IF
          END IF
        END DO
      END IF
      ! 3 point method                                    (end)

      ! psi(k) - new hydrometeor size distribution FUNCTION
      DO k=1,nrx1
        psi(k)=psi(k)/rr(k)
      END DO
      DO k=nrx1+1,nrx
        psi(k)=fi(k)
      END DO
      IF (idrop == 1) THEN
        DO k=krdrop_remaping_min,krdrop_remaping_max
          cdrop(k)=3.0_wp*col*psi(k)*rr(k)
        END DO
        ! kmax - right boundary spectrum of drop sdf
        !(krdrop_remap_min =< kmax =< krdrop_remap_max)
        DO k=krdrop_remaping_max,krdrop_remaping_min,-1
          kmax=k
          IF (psi(k).GT.0.0_wp) EXIT
        END DO

        DO k=kmax-1,krdrop_remaping_min,-1
          IF (cdrop(k).GT.0.0_wp) THEN
            delta_cdrop(k)=cdrop(k+1)/cdrop(k)
            IF (delta_cdrop(k).LT.coeff_remaping) THEN
              cdrop(k)=cdrop(k)+cdrop(k+1)
              cdrop(k+1)=0.0_wp
            END IF
          END IF
        END DO

        DO k=krdrop_remaping_min,kmax
          psi(k)=cdrop(k)/(3.0_wp*col*rr(k))
        END DO
      END IF
      ! in case isign_diffusional_growth.NE.0
    ELSE
      ! in case isign_diffusional_growth.EQ.0
      DO k=1,nrx
        psi(k)=fi(k)
      END DO
    END IF

    DO kr=1,nrx
      IF (psi(kr) < 0.0_wp) THEN
        PRINT*, 'psi(kr)<0',' before EXIT'
        PRINT*, 'isign_diffusional_growt', isign_diffusional_growth
        PRINT*, 'i3point_condevap', i3point_condevap
        PRINT*, 'k,rr(k),rn(k),k=1,nrx', (k,rr(k),rn(k),k=1,nrx)
        PRINT*, 'k,rr(k),rn(k),fi(k),psi(k),k=1,nrx'
        WRITE (*,'(A,1X,I2,2X,4D13.5)') (k,rr(k),rn(k),fi(k),psi(k),k=1,nrx)
        PRINT*, idrop,iin,kin
        CALL finish(TRIM(modname),"fatal error in SUBROUTINE jernewf psi(kr)<0, < min, model stop")
      END IF
    END DO

    RETURN
  END SUBROUTINE remapping

  SUBROUTINE nucleation_main(psi1_r,psi2_r,fccnr_r,fccnr_nucl_r,tt,qq,pp,sup1,sup2,sup2_old, &
                             win,is_this_cloudbase,rho)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: is_this_cloudbase
    REAL(KIND=wp),INTENT(IN) :: pp,win,rho
    REAL(KIND=wp),INTENT(INOUT) :: psi1_r(:),psi2_r(:,:),fccnr_r(:),fccnr_nucl_r(:),tt,qq,sup1,sup2,sup2_old
    INTEGER :: kr, ice, k
    REAL(KIND=wp) :: dropconcn(nkr), tpn, qpn, tpc, ror, &
                   sum_ice,del2n,fi2(nkr,icemax),rmassiaa_nucl,rmassibb_nucl, &
                   delmassice_nucl,es1n,es2n,ew1n,fccnr_nucl(nkr),psi1(nkr),psi2(nkr,icemax),fccnr(nkr)

    ror=rho*0.001_wp !total density in cgs

    ! ... adjust the input
    DO kr=1,nkr
      psi1(kr) = psi1_r(kr)               !drop size distribution
      IF ( iceprocs==1 ) THEN
        DO ice=1,icemax
          psi2(kr,ice) = psi2_r(kr,ice)   !ice size distribution
        END DO
      END IF
    END DO
    fccnr = fccnr_r             !ccn size distribution
    fccnr_nucl = fccnr_nucl_r   !nucleated ccn size distribution. It will be used when we turn on regeneration

    IF ( iceprocs==1 ) THEN
      tpn = tt
      qpn = qq
    END IF

    ! -----------------------------
    ! ... drop nucleation
    ! -----------------------------
    tpc = tt - tmelt
    IF (sup1>0.0_wp .AND. tpc > t_nucl_drop_min) THEN !do nucleation IF S over water >0 and t>t_nucl_drop_min=-80 celcius
      IF (sum(fccnr) > 0.0_wp)THEN                       !if ccn concentration >0
        dropconcn = 0.0_wp
        CALL nucleation_water(psi1,fccnr,fccnr_nucl,tt,ror,sup1,dropconcn,pp,is_this_cloudbase,win)

        DO kr=1,nkr
          qq=qq-(psi1(kr)-psi1_r(kr))*col*3.0_wp*xl(kr)*xl(kr)/ror !reduction of qq due to nucleation
        END DO
      END IF

      ! ... when t<-38 transfer nucleated drops to ice-crystals via direct homogenous nucleation
      IF ((tpc <= -38.0_wp) .AND. (iceprocs == 1)) THEN
        sum_ice = 0.0_wp
        DO kr=1,nkr      ! in future it is better to separate: small ones (till krfreez) to ice, and larger ones to graupel
          psi2(kr,2) = psi2(kr,2) + psi1(kr) ! psi1 is the frozen water, which is added to ice (plates) psi2(kr,2)
          sum_ice = sum_ice + col*3.0_wp*xl(kr)*xl(kr)*psi1(kr)
          psi1(kr) = 0.0_wp
        END DO
        tt = tt + latheat_freez*sum_ice/ror   ! latent heat of freezing. cp/cv bug was taken into account
      END IF
    END IF

    ! -------------------------------
    ! ... crystals nucleation
    ! -------------------------------
    IF ( iceprocs==1 ) THEN
      del2n = 100.0_wp*sup2     !supsat over ice in %
      tpc = tt-tmelt

      ! inhomogeneous nucleation in case of t<0 and t>t_nucl_ice_min=-38 celcius  and S over ice >0:
      IF (tpc < 0.0_wp .AND. tpc >= t_nucl_ice_min .AND. del2n > 0.0_wp) THEN
        DO kr=1,nkr
          DO ice=1,icemax
            fi2(kr,ice)=psi2(kr,ice)
          END DO
        END DO

        CALL nucleation_ice (psi2,sup2,tt,sup2_old) !Meyers method

        IF (isign_tq_icenucl == 1) THEN
          rmassiaa_nucl=0.0_wp     ! ice mass after nucleation
          rmassibb_nucl=0.0_wp     ! ice mass before nucleation

          ! before ice crystal nucleation
          DO k=1,nkr
            DO ice=1,icemax
              rmassibb_nucl=rmassibb_nucl+fi2(k,ice)*xi(k,ice)*xi(k,ice)
            END DO
          END DO
          rmassibb_nucl = rmassibb_nucl*col*3.0_wp/ror          !ice mass before nucleation
          IF (rmassibb_nucl < 0.0_wp) rmassibb_nucl = 0.0_wp

          ! after ice crystal nucleation
          DO k=1,nkr
            DO ice=1,icemax
              rmassiaa_nucl=rmassiaa_nucl+psi2(k,ice)*xi(k,ice)*xi(k,ice)
            END DO
          END DO
          rmassiaa_nucl = rmassiaa_nucl*col*3.0_wp/ror           !ice mass after nucleation
          IF (rmassiaa_nucl < 0.0_wp) rmassiaa_nucl=0.0_wp

          delmassice_nucl = rmassiaa_nucl-rmassibb_nucl     !ice mass difference due to nucleation

          qpn = qq-delmassice_nucl       !water vapor change due to ice nucleation
          qq = qpn

          IF (.not. latheatfac3) THEN !no signficant effect, related to cp/cv bug
            tpn = tt + al2*delmassice_nucl !temperature change due to ice nucleation
          ELSE
            tpn = tt + al2*delmassice_nucl*cpd/cvd
          END IF
          tt = tpn

          es1n = 10.0_wp*sat_pres_water(tpn) !dynes/cm^2
          es2n = 10.0_wp*sat_pres_ice(tpn)   !dynes/cm^2
          ew1n=10.0_wp*qpn*rho*rv*tpn        !!dynes/cm^2

          sup1 = ew1n/es1n-1.0_wp         !sup sat over water due to ice nucleation
          sup2 = ew1n/es2n-1.0_wp         !sup sat over ice due to ice nucleation

        END IF ! in case isign_tq_icenucl==1
      END IF ! in case tpc < 0.0_wp .AND. tpc >= t_nucl_ice_min .AND. del2n > 0.0_wp
    END IF

    ! -----------------------------
    ! ... output
    ! -----------------------------
    DO kr=1,nkr
      psi1_r(kr) = psi1(kr)               !new psd of water
      IF ( iceprocs==1 ) THEN
        DO ice=1,icemax
          psi2_r(kr,ice) = psi2(kr,ice)   !new psd of ice crystals
        END DO
      END IF
    END DO
    fccnr_r = fccnr            !new psd of ccn
    fccnr_nucl_r = fccnr_nucl  !new psd of nucleated ccn

    RETURN
  END SUBROUTINE nucleation_main

  SUBROUTINE nucleation_water (psi1, fccnr, fccnr_nucl, tt, ror, sup1, & !called from nucleation_main
                            dropconcn, pp, is_this_cloudbase, win)
    ! psi1(kr), 1/g/cm3 - non conservative drop size distribution function
    ! fccnr(kr), 1/cm^3 - aerosol(ccn) non conservative, size distribution function
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: is_this_cloudbase
    REAL(KIND=wp),INTENT(IN) :: ror, pp, win
    REAL(KIND=wp),INTENT(INOUT) :: fccnr(:), fccnr_nucl(:), psi1(:), dropconcn(:), tt, sup1
    INTEGER :: imax, i, ncriti, kr
    REAL(KIND=wp) :: dx,ar2,rcriti,ccnconc(nkr),akoe,bkoe,rccn_minimum,dln1,dln2,rmassl_nucl

    dropconcn(:) = 0.0_wp
    imax = nkr !right ccn spectrum boundary
    DO i=imax,1,-1
      IF (fccnr(i) > 0.0_wp) THEN
        imax = i
        EXIT
      END IF
    END DO

    ncriti=0
    ! every iteration we will nucleate one bin, then we will check the new supersaturation and new rcriti
    DO WHILE (imax>=ncriti)
      ccnconc = 0.0_wp
      ! akoe & bkoe - constants in koehler equation
      akoe=3.3e-05_wp/tt
      bkoe = ions*4.3_wp/mwaero
      bkoe=bkoe*(4.0_wp/3.0_wp)*pi*ro_solute

      IF (use_cloud_base_nuc == 1) THEN ! currently 0 (off) in mo_sbm_utils
        IF (is_this_cloudbase == 1) THEN
          CALL nucleation_water_cbase_supsat (fccnr, tt, pp, win, rcriti)
        ELSE
          rcriti = (akoe/3.0_wp)*(4.0_wp/bkoe/sup1/sup1)**(1.0_wp/3.0_wp) !critical radius of "dry" aerosol (cm)
        END IF
      ELSE ! ismax_cloud_base==0
        rcriti=(akoe/3.0_wp)*(4.0_wp/bkoe/sup1/sup1)**(1.0_wp/3.0_wp) !critical radius of "dry" aerosol (cm)
      END IF

      IF (rcriti >= rccn(imax)) EXIT ! nothing to nucleate
      ! find the minimum bin to nucleate
      ncriti = imax
      DO WHILE (rcriti<=rccn(ncriti) .AND. ncriti>1)
        ncriti=ncriti-1
      END DO
      rccn_minimum = rccn(1)/10000.0_wp !minimum aerosol(ccn) radius
      ! calculation of ccnconc(ii)=fccnr(ii)*col - aerosol(ccn) bin concentrations, ii=imin,...,imax
      ! determination of ncriti   - number bin in which is located rcriti
      ! calculation of ccnconc(ncriti)=fccnr(ncriti)*dln1/(dln1+dln2),
      ! where,
      ! dln1=ln(rcriti)-ln(rccn_minimum)
      ! dln2=ln(rccn(1)-ln(rcriti)
      ! calculation of new value of fccnr(ncriti)

      ! update:
      ! the problem is that ncriti=1 and imax=2 may occur 2 times:
      ! a. rccn(1)<rcriti<rccn(2) --> we should clean part of the "1-2" bin: fccnr(imax)=fccnr(imax)*dln1/col
      ! b. rcriti<rccn(1)         --> we should clean the entire bin:  fccnr(imax)=0
      ! and there is an option when criti=1 and imax=1, THEN:
      ! c. rcriti<rccn(1)         --> we should clean part of the "0-1" bin: fccnr(imax)=fccnr(imax)*dln1/(dln1+dln2)
      !                               this bin is special and has a width of dln1+dln2 and not just col
      !---------------------------------------------------
      IF ((imax-1>ncriti) .OR. ((ncriti==1) .AND. (imax==2) .AND. (rcriti<rccn(1)))) THEN ! imax bin should be cleaned to zero
        ccnconc(imax) = col*fccnr(imax)
        fccnr_nucl(imax) = fccnr_nucl(imax) + fccnr(imax)
        fccnr(imax) = 0.0_wp
      ELSE IF ((ncriti==1) .AND. (imax==1) .AND. (rcriti<rccn(1))) THEN ! rcriti<rccn(1) we should clean part
                                                                        ! of the "0-1" bin rccn_minimum<-->rccn(1)
        dln1=LOG(rcriti/rccn_minimum)
        dln2=LOG(rccn(1)/rcriti)
        ccnconc(imax)=dln2*fccnr(imax)
        fccnr_nucl(imax) = fccnr_nucl(imax) + fccnr(imax)*(1.0_wp - (dln1/LOG(rccn(1)/rccn_minimum)))
        fccnr(imax)=fccnr(imax)*dln1/LOG(rccn(1)/rccn_minimum)
      ELSE IF ((ncriti==imax-1) .AND. (rcriti > rccn(imax-1))) THEN ! rccn(ncriti)<rcriti<rccn(imax) imax bin should
                                                                    ! be cleaned partially
        dln1=LOG(rcriti/rccn(imax-1))
        dln2=col-dln1
        ccnconc(imax)=dln2*fccnr(imax)
        fccnr_nucl(imax) = fccnr_nucl(imax) + fccnr(imax)*(1.0_wp - dln1/col)
        fccnr(imax)=fccnr(imax)*dln1/col
      ELSE
        WRITE (txt,*) 'rcriti,rccn1,ncriti,imax=',rcriti,rccn(1),ncriti,imax
        CALL finish(TRIM(modname),'ccn bins problem, model stop' )
      END IF

      ! calculate the mass change due to nucleation
      rmassl_nucl=0.0_wp
      IF (imax <= nkr-7) THEN ! we pass it to drops mass grid
        dropconcn(1) = dropconcn(1)+ccnconc(imax)
        rmassl_nucl = rmassl_nucl+ccnconc(imax)*xl(1)*xl(1)
      ELSE
        dropconcn(8-(nkr-imax)) = dropconcn(8-(nkr-imax))+ccnconc(imax)
        rmassl_nucl = rmassl_nucl + ccnconc(imax)*xl(8-(nkr-imax))*xl(8-(nkr-imax))
      END IF
      rmassl_nucl = rmassl_nucl*col*3.0_wp/ror
      ! prepering to check IF we need to nucleate the next bin
      imax = imax-1
      ! cycle imax>=ncriti
    END DO

    ! intergarting for including the previous nucleated drops
    IF (sum(dropconcn) > 0.0_wp)THEN
      DO kr = 1,8
        dx = 3.0_wp*col*xl(kr)
        psi1(kr) = psi1(kr)+dropconcn(kr)/dx
      END DO
    END IF

    RETURN
  END SUBROUTINE nucleation_water

  SUBROUTINE nucleation_water_cbase_supsat (fccnr,tt,pp,wbase,rcriti) ! currently not called since use_cloud_base_nuc=0
    ! fccnr(kr), 1/cm^3 - aerosol(ccn) non conservative, size
    !                     distribution function in point with x,z
    !                     coordinates, kr=1,...,nkr
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) ::  tt, pp, wbase
    REAL(KIND=wp),INTENT(INOUT) :: fccnr(:), rcriti
    INTEGER :: nr, nn, kr
    REAL(KIND=wp) :: pl(nkr), supmax(nkr), akoe, bkoe, c3, pr, ccnconact, dl1, dl2

    CALL cbase_supsat_max_coef(akoe,bkoe,c3,pp,tt)
    ! supmax calculation
    ! 'analytical estimation of droplet concentration at cloud base', eq.21, 2012
    ! calculation of right side hand of equation for s_max
    ! WHILE wbase>0, calculation pr

    pr = c3*wbase**(0.75_wp)

    ! calculation supersaturation in cloud base
    supmax = 999.0_wp
    pl = 0.0_wp
    nn = -1
    DO nr=2,nkr
      supmax(nr)=SQRT((4.0_wp*akoe**3.0_wp)/(27.0_wp*rccn(nr)**3.0_wp*bkoe))
      ! calculation ccnconact- the concentration of ccn that were activated
      ! following nucleation
      ! ccnconact=n from the paper
      ! 'analytical estimation of droplet concentration at cloud base', eq.19, 2012
      ! ccnconact, 1/cm^3- concentration of activated ccn = new droplet concentration
      ! ccnconact=fccnr(kr)*col
      ! col=ln2/3

      ccnconact=0.0_wp

      ! nr represents the number of bin in which rcriti is located
      ! from nr bin to nkr bin goes to droplets
      DO kr=nr,nkr
        ccnconact = ccnconact + col*fccnr(kr)
      END DO
      ! calculate lhs of equation for s_max
      ! when pl<pr ccn will activate
      pl(nr)=supmax(nr)*(SQRT(ccnconact))
      IF (pl(nr).LE.pr) THEN
        nn = nr
        EXIT
      END IF
    END DO ! nr

    IF (nn == -1) THEN
      PRINT*,"pr, wbase [cm/s], c3",pr,wbase,c3,"pl",pl
      CALL finish(TRIM(modname),'nn is not defined in cloud base routine, model stop' )
    END IF

    ! linear interpolation- finding radius criti of aerosol between
    ! bin number (nn-1) to (nn)
    ! 1) finding the difference between pl and pr in the left and right over the
    ! final bin.

    dl1 = ABS(pl(nn-1)-pr) ! left side in the final bin
    dl2 = ABS(pl(nn)-pr)   ! right side in the final bin

    ! 2) fining the left part of bin that will not activate
    !	  dln1=col*dl1/(dl2+dl1)
    ! 3) finding the right part of bin that activate
    !	  dln2=col-dln1
    ! 4) finding radius criti of aerosol- rcriti

    rcriti = rccn(nn-1)*EXP(col*dl1/(dl1+dl2))
    ! end linear interpolation

    RETURN
  END SUBROUTINE nucleation_water_cbase_supsat

  SUBROUTINE cbase_supsat_max_coef (akoe,bkoe,c3,pp,tt) ! currently not called since cbase_supsat_max_coef = 0
    ! akoe, cm- constant in koehler equation
    ! bkoe    - constant in koehler equation
    ! f, cm^-2*s-  from koehler equation
    ! c3 - coefficient depends on thermodynamical parameters
    ! pp, (dynes/cm/cm)- pressure
    ! tt, (k)- temperature
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: pp, tt
    REAL(KIND=wp),INTENT(OUT) :: akoe, bkoe, c3
    REAL(KIND=wp) :: alw1,sw,ro_w,hc,ew,ro_v,dv,ro_a,fl,fr,f,tpc,qv,a1,a2,c1,c2,g

    g=grav*100.0_wp
    tpc = tt-tmelt
    ! cgs :
    ! cp=1005.0d4 cm*cm/sec/sec/kelvin- specific heat capacity of moist air at constant pressure
    ! g=9.8d2 cm/sec/sec- acceleration of gravity
    ! rd=287.0d4 cm*cm/sec/sec/kelvin - individual gas constant for dry air
    ! al2_my=2.834d10 cm*cm/sec/sec - latent heat of sublimation
    ! al1_my,         cm*cm/sec/sec - latent heat of vaporization
    ! alw1=al1_my - alw1 depends on temperature
    ! alw1, [m^2/sec^2] -latent heat of vaporization-
    alw1 = -6.143419998e-2_wp*tpc**(3.0_wp)+1.58927_wp*tpc**(2.0_wp)-2.36418e3_wp*tpc+2.50079e6_wp
    ! alw1, [cm^2/sec^2]
    alw1 = alw1*10.0e3_wp
    ! sw, [n*m^-1] - surface tension of water-air interface
    IF (tpc.LT.-5.5_wp) THEN
      sw = 5.285e-11_wp*tpc**(6.0_wp)+6.283e-9_wp*tpc**(5.0_wp)+ &
         2.933e-7_wp*tpc**(4.0_wp)+6.511e-6_wp*tpc**(3.0_wp)+ &
         6.818e-5_wp*tpc**(2.0_wp)+1.15e-4_wp*tpc+7.593e-2_wp
    ELSE
      sw = -1.55e-4_wp*tpc+7.566165e-2_wp
    END IF
    ! sw, [g/sec^2]
    sw = sw*10.0e2_wp
    ! ro_w, [kg/m^3] - density of liquid water
    IF (tpc.LT.0.0_wp) THEN
      ro_w= -7.497e-9_wp*tpc**(6.0_wp)-3.6449e-7_wp*tpc**(5.0_wp) &
          -6.9987e-6_wp*tpc**(4.0_wp)+1.518e-4_wp*tpc**(3.0_wp) &
          -8.486e-3_wp*tpc**(2.0_wp)+6.69e-2_wp*tpc+9.9986e2_wp
    ELSE
      ro_w=(-3.932952e-10_wp*tpc**(5.0_wp)+1.497562e-7_wp*tpc**(4.0_wp) &
          -5.544846e-5_wp*tpc**(3.0_wp)-7.92221e-3_wp*tpc**(2.0_wp)+ &
          1.8224944e1_wp*tpc+9.998396e2_wp)/(1.0_wp+1.8159725e-2_wp*tpc)
    END IF
    ! ro_w, [g/cm^3]
    ro_w=ro_w*1.0e-3_wp
    ! hc, [kg*m/kelvin*sec^3] - coefficient of air heat conductivity
    hc=7.1128e-5_wp*tpc+2.380696e-2_wp
    ! hc, [g*cm/kelvin*sec^3]
    hc=hc*10.0e4_wp

    ! ew, water vapor pressure ! ... ks (kg/m2/sec)
    ew = 6.38780966e-9_wp*tpc**(6.0_wp)+2.03886313e-6_wp*tpc**(5.0_wp)+ &
       3.02246994e-4_wp*tpc**(4.0_wp)+2.65027242e-2_wp*tpc**(3.0_wp)+ &
       1.43053301_wp*tpc**(2.0_wp)+4.43986062e1_wp*tpc+6.1117675e2_wp

    ! ew, [g/cm*sec^2]
    ew=ew*10.0_wp
    ! akoe & bkoe - constants in koehler equation
    !ro_solute=2.16_wp
    akoe=2.0_wp*sw/(rv_cgs*ro_w*(tpc+tmelt))
    bkoe = ions*4.3_wp/mwaero
    bkoe=bkoe*(4.0_wp/3.0_wp)*pi*ro_solute

    ! ro_v, g/cm^3 - density of water vapor, calculate from equation of state for water vapor
    ro_v = ew/(rv_cgs*(tpc+tmelt))
    ! dv,  [cm^2/sec] - coefficient of diffusion
    ! 'pruppacher, h.r., klett, j.d., 1997. microphysics of clouds and precipitation page 503, eq. 13-3'
    dv = 0.211_wp*(pzero/pp)*((tpc+tmelt)/tmelt)**(1.94_wp)

    ! qv,  g/g- water vapor mixing ratio
    ! ro_a, g/cm^3 - density of air, from equation of state
    ro_a=pzero/((tpc+tmelt)*rd)

    ! f, s/m^2 - coefficient depending on thermodynamics parameters
    !            such as temperature, thermal conductivity of air, etc
    ! left side of f equation
    fl=(ro_w*alw1**(2.0_wp))/(hc*rv_cgs*(tpc+tmelt)**(2.0_wp))

    ! right side of f equation
    fr = ro_w*rv_cgs*(tpc+tmelt)/(ew*dv)
    f = fl + fr

    ! qv, g/g - water vapor mixing ratio
    qv=ro_v/ro_a

    ! a1,a2 -  terms from equation describing changes of supersaturation in an adiabatic cloud air parcel
    ! a1, [cm^-1] - constant
    ! a2, [-]     - constant
    ! Note: this code is note called (see above), but cp/cv bug should be considered below:
    IF (.not. latheatfac4) THEN
      a1=(g*alw1/(cpd*rv_cgs*(tpc+tmelt)**(2.0_wp)))-(g/(rd*(tpc+tmelt)))
      a2=(1.0_wp/qv)+(alw1**(2.0_wp))/(cpd*rv_cgs*(tpc+tmelt)**(2.0_wp))
    ELSE
      a1=(g*alw1/(cvd*rv_cgs*(tpc+tmelt)**(2.0_wp)))-(g/(rd*(tpc+tmelt)))
      a2=(1.0_wp/qv)+(alw1**(2.0_wp))/(cvd*rv_cgs*(tpc+tmelt)**(2.0_wp))
    END IF

    ! c1,c2,c3,c4- constant parameters
    c1=1.058_wp
    c2=1.904_wp
    c3=c1*(f*a1/3.0_wp)**(0.75_wp)*SQRT(3.0_wp*ro_a/(4.0_wp*pi*ro_w*a2))

    RETURN
  END SUBROUTINE cbase_supsat_max_coef

  SUBROUTINE nucleation_ice (psi2,sup2,tt,sup2_old) !called from nucleation_main
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(INOUT) :: psi2(:,:),tt,sup2,sup2_old
    REAL(KIND=wp) :: fi2(nkr,icemax), conci_bfnucl(icemax), conci_afnucl(icemax)
    REAL(KIND=wp),parameter :: a1 = -0.639_wp, b1 = 0.1296_wp, a2 = -2.8_wp, b2 = 0.262_wp, &
                                       temp1 = -5.0_wp, temp2 = -2.0_wp, temp3 = -20.0_wp
    REAL(KIND=wp),parameter::c1_mey = 1.0e-3_wp ! deposition-condensation. concentration of 1/litre in Meyers --> 1/cm^3
    REAL(KIND=wp),parameter::c2_mey = 0.0_wp ! contact-freezing - neglected
    INTEGER,parameter :: nrgi = 2
    INTEGER :: kr,ice,itype
    REAL(KIND=wp) :: tpc,del2n,del2nn,helek1,helek2,ff1bn,fact,dsup2n,deltacd,deltaf, &
                             addf,delconci_afnucl,tpcc,dx

    ! size distribution functions of crystals before ice nucleation
    DO kr=1,nkr
      DO ice=1,icemax
        fi2(kr,ice)=psi2(kr,ice)
      END DO
    END DO

    ! calculation concentration of crystals before ice nucleation
    DO ice=1,icemax
      conci_bfnucl(ice)=0.0_wp
      DO kr=1,nkr
        conci_bfnucl(ice)=conci_bfnucl(ice)+3.0_wp*col*psi2(kr,ice)*xi(kr,ice)
      END DO
    END DO

    ! type of ice with nucleation:
    ! depending on T range, different types of crystals are nucleated:
    tpc = tt-tmelt
    itype=0

    IF ((tpc>-4.0_wp).OR.(tpc<=-8.1_wp.AND.tpc>-12.7_wp).OR.(tpc<=-17.8_wp.AND.tpc>-22.4_wp)) THEN
      itype=2
    ELSE
      IF ((tpc<=-4.0_wp.AND.tpc>-8.1_wp).OR.(tpc<=-22.4_wp)) THEN
        itype=1
      ELSE
       itype=3
      END IF
    END IF

    ! 2 algorithms exist for ice nucleation - deposition condensation freezing (helek1) and contact freezing nucleation (helek2)
    ! crystal size distribution function:
    ice=itype
    IF (tpc < temp1) THEN
      del2n = 100.0_wp*sup2
      del2nn = del2n
      IF ( del2n > delsupice_max) del2nn = delsupice_max ! 59%
      helek1 = c1_mey*EXP(a1+b1*del2nn) !Meyers formula: activated ice particles vs S
    ELSE
      helek1 = 0.0_wp
    END IF

    IF (tpc < temp2) THEN
      tpcc = tpc
      IF (tpcc < temp3) tpcc = temp3
      helek2 = c2_mey*EXP(a2-b2*tpcc)  ! c2 is set to 0 instead of 1 (helek 2 is switched off)
    ELSE
      helek2 = 0.0_wp
    END IF

    ff1bn = helek1+helek2
    fact=1.0_wp
    !dsup2n = (sup2-sup2_old+dsupice_xyz)*100.0_wp !it was: Ds/Dt=ds/dt+uds/dx+vds/dy+wds/dz, but we keep only ds/dt
    dsup2n = (sup2-sup2_old)*100.0_wp ! change in sup sat over ice during ncond sub step
    sup2_old = sup2 ! we calculate sup2_old outside of the subroutine

    IF (dsup2n > delsupice_max) dsup2n = delsupice_max ! limit the sup sat change to 59%, otherwise too many crystals
                                                       ! will nucleate (corresponds to the allowed range of Meyers formula)

    deltacd = ff1bn*b1*dsup2n ! the number of nucleated crystals is proportional to the change in sup sat

    IF (deltacd>=ff1bn) deltacd=ff1bn ! change in ice particles is limitted to the max according Meyers

    IF (deltacd>0.0_wp) THEN
      deltaf=deltacd*fact
      ! concentration of ice crystals
      IF (conci_bfnucl(ice)<=helek1) THEN
        DO kr=1,nrgi-1
          dx=3.0_wp*xi(kr,ice)*col
          addf=deltaf/dx
          psi2(kr,ice)=psi2(kr,ice)+addf ! updated PSD of ice crystals
        END DO
      END IF
    END IF

    ! calculation of crystal concentration after ice nucleation
    DO ice=1,icemax
      conci_afnucl(ice)=0.0_wp
      DO kr=1,nkr
        conci_afnucl(ice)=conci_afnucl(ice)+ &
                3.0_wp*col*psi2(kr,ice)*xi(kr,ice) ! updated concentration of ice crystals
      END DO
      delconci_afnucl=ABS(conci_afnucl(ice)-conci_bfnucl(ice))
      IF (delconci_afnucl>10.0_wp) THEN
        PRINT*, 'in SUBROUTINE nucleation_ice, after nucleation'
        PRINT*, 'because delconci_afnucl > 10/cm^3'
        PRINT*, 'conci_bfnucl(ice),conci_afnucl(ice)'
        PRINT*, conci_bfnucl(ice),conci_afnucl(ice)
        PRINT*, 'deltacd,dsup2n,ff1bn,b1,sup2'
        PRINT*, deltacd,dsup2n,ff1bn,b1,sup2
        PRINT*, 'kr,   fi2(kr,ice),   psi2(kr,ice), kr=1,nkr'
        PRINT*, (kr,   fi2(kr,ice),   psi2(kr,ice), kr=1,nkr)
        PRINT*, 'stop 099 : delconci_afnucl(ice) > 10/cm^3'
        stop 099
      END IF
    END DO

    RETURN
  END SUBROUTINE nucleation_ice

  SUBROUTINE stick_eff_interpol (nh, h_tab, x_tab, h, x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nh
    REAL(KIND=wp), INTENT(IN) :: h_tab(nh), x_tab(nh), h
    REAL(KIND=wp), INTENT(INOUT) :: x
    INTEGER :: i, j

    IF (h > h_tab(1)) THEN
      x = x_tab(1)
      RETURN
    END IF

    IF (h < h_tab(nh)) THEN
      x = x_tab(nh)
      RETURN
    END IF

    DO i = 2,nh
      IF (h > h_tab(i)) THEN
        j = i-1
        x = x_tab(j)+(x_tab(i)-x_tab(j))/(h_tab(i)-h_tab(j))*(h-h_tab(j))
        RETURN
      END IF
    END DO

    RETURN
  END SUBROUTINE stick_eff_interpol

  SUBROUTINE coll_ice_stick_eff (tt,qq,pp,factor_t)
    ! Sticking efficiency is very important but uncertain
    ! Generally, as the temperature is lower the sticking efficiency is lower,
    ! but on the other hand dendrites collide efficiently at -15C
    ! Moreover, the more humid is the air, the higher is the sticking efficiency
    ! The complexity comes also from the fact that we neglect detailed sepration between the 3 ice types here
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: tt, pp
    REAL(KIND=wp),INTENT(INOUT) :: qq
    REAL(KIND=wp),INTENT(OUT) :: factor_t ! sticking efficiency. Collision kernel which includes ice is multiplied by factor_t
    REAL(KIND=wp) :: satq2, temp, epsf, tc, qs2, qq1, tc_min, tc_max, factor_max, &
                    factor_min, t, a, b, c, p, d, at, bt, ct, dt, t_tab(7), se_tab(7)

    satq2(t,p) = 3.80e3_wp*(10.0_wp**(9.76421_wp-2667.1_wp/t))/p
    temp(a,b,c,d,t) = d*t*t*t+c*t*t+b*t+a

    tc = tt - tmelt
    IF (tc > 0.0_wp) RETURN

    SELECT CASE (stick_param1)
    CASE(1)
      DATA at, bt, ct, dt /0.88333_wp,  0.0931878_wp,  0.0034793_wp,  4.5185186e-05_wp/
      IF (qq.LE.0.0_wp) qq = epsil
      epsf = 0.5_wp
      tc = tt - tmelt
      qs2    =satq2(tt,pp)
      qq1    =qq*(0.622_wp+0.378_wp*qs2)/(0.622_wp+0.378_wp*qq)/qs2
      IF (tc.GE.-6.0_wp) THEN
        factor_t = temp(at,bt,ct,dt,tc)*qq1
        IF (factor_t.LT.epsf) factor_t = epsf
        IF (factor_t.GT.1.0_wp) factor_t = 1.0_wp
      END IF
      !The following values are uncertain and influence the ammount of resulting large snow flakes.
      !in the range of -12C - -17C the crystals are dendrites which stick efficiently, therefore
      !large values are chosen below
      IF (stick_param2 == 0) THEN
        IF (tc.GE.-12.5_wp .AND. tc.LT.-6.0_wp) factor_t = 0.5_wp
        IF (tc.GE.-17.0_wp .AND. tc.LT.-12.5_wp) factor_t = 1.0_wp
        IF (tc.GE.-20.0_wp .AND. tc.LT.-17.0_wp) factor_t = 0.4_wp
      ELSE
        IF (tc.GE.-12.5_wp .AND. tc.LT.-6.0_wp) factor_t = 0.6_wp
        IF (tc.GE.-17.0_wp .AND. tc.LT.-12.5_wp) factor_t = 0.25_wp
        IF (tc.GE.-20.0_wp .AND. tc.LT.-17.0_wp) factor_t = 0.05_wp
      END IF
      IF (tc.LT.-20.0_wp) THEN
        tc_min = ttcoal-tmelt
        tc_max = -20.0_wp
        IF (stick_param2 == 0)THEN
          factor_max = 0.4_wp
          factor_min = 0.0_wp
        ELSE
          factor_max = 0.05_wp
          factor_min = 0.0_wp
        END IF
        factor_t = factor_min + (tc-tc_min)*(factor_max-factor_min)/(tc_max-tc_min)
      END IF
      IF (tc.LT.-40.0_wp) THEN
        factor_t = 0.0_wp
      END IF
      IF (factor_t > 1.0_wp) factor_t = 1.0_wp
      IF (tc.GE.0.0_wp) THEN
        factor_t = 1.0_wp
      END IF

    CASE(2)
      ! ... linear
      t_tab =  [0.0_wp, -0.813_wp, -5.26_wp, -10.13_wp, -14.63_wp, -20.02_wp, -40.0_wp ]
      se_tab = [10.0_wp**(-0.693_wp), 10.0_wp**(-0.72_wp), 10.0_wp**(-0.877_wp), &
                10.0_wp**(-1.050_wp),  10.0_wp**(-1.212_wp),  10.0_wp**(-1.401_wp),  10.0_wp**(-2.082_wp) ]
      CALL stick_eff_interpol (size(se_tab), t_tab, se_tab, tc, factor_t)
      IF (tc < -40.0_wp) factor_t = 0.0_wp
      IF ((factor_t > 1.0_wp) .OR. (tc > 0.0_wp)) factor_t = 1.0_wp
    END SELECT

    RETURN
  END SUBROUTINE coll_ice_stick_eff

  SUBROUTINE breakup_rain_spont (dtime_spon_break, ff1r)
    ! +-----------------------------------------------------------------------------+
    ! i_break_method=1: Spontaneous breakup according to Srivastava1971_JAS -
    ! Size distribution of raindrops generated by their breakup and coalescence
    ! i_break_method=2: Spontaneous breakup according to Kamra et al 1991 JGR -
    ! SPONTANEOUS BREAKUP OF CHARGED AND UNCHARGED WATER DROPS FREELY SUSPENDED IN A WIND TUNNEL
    ! Description of variables:
    ! FF1R(KR), 1/g/cm3 - non conservative drop size distribution
    ! XL(kr), g - Mass of liquid drops
    ! prob, dimensionless - probability for breakup
    ! dropconc_bf(kr), cm^-3 - drops concentration before breakup
    ! dropconc_af(kr), cm^-3 - drops concentration after breakup
    ! drops_break(kr), cm^-3 - concentration of breaking drops
    ! +-----------------------------------------------------------------------------+
    IMPLICIT NONE
    REAL(KIND=wp), INTENT(INOUT) :: ff1r(:)
    REAL(KIND=wp), INTENT(IN) :: dtime_spon_break
    INTEGER :: kr,i,imax,j
    REAL(KIND=wp) :: dm, tmp_1, tmp_2, tmp_3, &
                     dropconc_bf(nkr), dropconc_af(nkr), drops_break(nkr), psi1(nkr)

    IF (sum(ff1r) <= nkr*1.e-30_wp) RETURN

    imax=nkr
    DO i=nkr,1,-1
      imax=i
      IF (ff1r(i) > 0.0_wp) EXIT !find the largest existing drop size
    END DO

    !drops of size < ikr_spon_break do not break. In mo_sbm_util.f90 it is the largest bin before 0.3 cm
    !in other words, only drops >=0.3cm break
    IF (imax<ikr_spon_break) RETURN

    ! initialization:
    psi1(:)=ff1r(:)
    drops_break(:)=0.0_wp
    dropconc_bf(:)=0.0_wp

    ! b) calculation of concentration of raindrops in all bins
    DO kr=1,imax
      dm=3.0_wp*col*xl(kr)
      dropconc_bf(kr)=dropconc_bf(kr)+dm*psi1(kr)   !drop concentration before breakup
    END DO
    dropconc_af(:)=dropconc_bf(:)   !initialization

    ! c+d) calculation of number of breaking drops  and the concentration of drops remaining in particular bin
    DO kr=imax,ikr_spon_break,-1
      !dropconc_af(kr)=dropconc_bf(kr)/(1+prob(kr)*dtime_spon_break)
      tmp_1 = prob(kr)*dtime_spon_break
      tmp_2 = EXP(-tmp_1)  !exponential decrease with time of concentration in the breaking bin
      tmp_3 = dropconc_bf(kr)
      dropconc_af(kr) = tmp_2*tmp_3
      !dropconc_af(kr) = EXP(-dtime_spon_break*prob(kr))*dropconc_bf(kr)
      drops_break(kr) = dropconc_bf(kr)-dropconc_af(kr)   !the concentration of the broken drops in this bin
      !IF (dropconc_af(kr)<0.0_wp) stop 'spontaneous breakup'
    END DO

    ! e) recalculation of dsd in bin j using new concentration
    !        DO kr=ikr_spon_break,imax
    !           dm=3.0_wp*col*xl(kr)
    !           psi1(kr)=psi1(kr)-drops_break(kr)/dm
    !        END DO

    ! f+g) redistributing and calculations drops concentration over smaller (i<j) bins
    !
    SELECT CASE (i_break_method)
    CASE(1)
      DO j=ikr_spon_break,imax
        DO i=1,j-1
          !adding water to small bins, using gain_var_new from calculations in spontanous_init
          dropconc_af(i)=dropconc_af(i)+drops_break(j)*gain_var_new(j,i)*xl(i)
        END DO
      END DO
    CASE(2)
      DO j=ikr_spon_break,imax
        DO i=1,j-1
          dropconc_af(i)=dropconc_af(i)+drops_break(j)*nnd(j,i)
        END DO
      END DO
    END SELECT

    ! h) recalculation of dsd in bins kr using new concentrations
    DO kr=1,imax
      dm=3.0_wp*col*xl(kr)
      psi1(kr)=dropconc_af(kr)/dm
    END DO

    ff1r(:)=psi1(:)

    RETURN
  END SUBROUTINE breakup_rain_spont

  SUBROUTINE breakup_snow_spont(f)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(INOUT) :: f(:) !snow PSD
    REAL(KIND=wp) :: g(nkr), break_snow(nkr)
    INTEGER :: kr

    ! break_snow are probabilities of break-up. For 33 bins it exists only for bins 31,32,33
    ! bin 30 is ~1cm, kr_snow_min is bin 31
    DO kr=1,nkr-3
      break_snow(kr)=0.0_wp
    END DO
    break_snow(nkr-2)=0.004_wp
    break_snow(nkr-1)=0.012_wp
    break_snow(nkr)=0.02_wp

    DO kr=1,nkr
      g(kr)=f(kr)*3.0_wp*xs(kr)*xs(kr) !PSD --> masses
    END DO

    DO kr=nkr,kr_snow_min,-1
      g(kr-1) = g(kr-1)+g(kr)*break_snow(kr)
      f(kr-1) = g(kr-1)/3.0_wp/xs(kr-1)/xs(kr-1)

      g(kr) = g(kr)*(1.0_wp-break_snow(kr))
      f(kr) = g(kr)/3.0_wp/xs(kr)/xs(kr)
    END DO

    RETURN
  END SUBROUTINE breakup_snow_spont

  SUBROUTINE freezing(ff1,ff2,ff3,ff4,ff5,tin,dt,ro)
    IMPLICIT NONE
    REAL(KIND=wp), INTENT(INOUT) :: ff1(:),ff2(:,:),ff3(:),ff4(:),ff5(:)
    INTEGER kr,ice,ice_type
    REAL(KIND=wp):: dt,ro,pf,del_t,yk2,cfreez,sum_ice,tin,ttin, &
                    f1_max,f2_max,f3_max,f4_max,f5_max

    ttin=tin
    del_t=ttin-tmelt
    ice_type=2
    f1_max=0.0_wp
    f2_max=0.0_wp
    f3_max=0.0_wp
    f4_max=0.0_wp
    f5_max=0.0_wp
    DO kr=1,nkr
      f1_max=MAX(f1_max,ff1(kr))
      f3_max=MAX(f3_max,ff3(kr))
      f4_max=MAX(f4_max,ff4(kr))
      f5_max=MAX(f5_max,ff5(kr))
      DO ice=1,icemax
        f2_max=MAX(f2_max,ff2(kr,ice))
      END DO
    END DO
    !******************************* freezing ****************************
    IF (del_t.LT.0.AND.f1_max.ne.0) THEN
      sum_ice=0.0_wp
      cfreez = (bfreezmax-bfreezmy)/xl(nkr)
      !***************************** mass loop ***************************
      DO kr=1,nkr
        pf = xl(kr)*afreezmy*exp(-(bfreezmy+cfreez*xl(kr))*del_t)
        yk2 = ff1(kr)*(1.0_wp-exp(-pf*dt))
        ff1(kr) = ff1(kr)*exp(-pf*dt)
        IF (kr.LE.krfreez) THEN !bin 21
          ff2(kr,ice_type)=ff2(kr,ice_type)+yk2
        ELSE
          ff5(kr) = ff5(kr)+yk2
        END IF
        sum_ice=sum_ice+yk2*3.0_wp*xl(kr)*xl(kr)*col
      END DO
      !************************** new temperature *************************
      tin = ttin+latheat_freez*sum_ice/ro ! cp/cv bug was taken into accound (was 333.*sum_ice/ro)
    END IF

    RETURN
  END SUBROUTINE freezing

  SUBROUTINE melting(ff1,ff2,ff3,ff4,ff5,tin,dt,ro)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN)    :: dt,ro
    REAL(KIND=wp),INTENT(INOUT) :: ff1(:),ff2(:,:),ff3(:),ff4(:),ff5(:),tin
    INTEGER :: kr,ice
    REAL(KIND=wp) :: arg_m,sum_ice,ff_max,f2_max,f3_max,f4_max,f5_max, del_t,meltrate

    del_t = tin-tmelt
    f2_max=0.0_wp
    f3_max=0.0_wp
    f4_max=0.0_wp
    f5_max=0.0_wp
    DO kr=1,nkr
      f3_max=MAX(f3_max,ff3(kr))
      f4_max=MAX(f4_max,ff4(kr))
      f5_max=MAX(f5_max,ff5(kr))
      DO ice=1,icemax
        f2_max=MAX(f2_max,ff2(kr,ice))
      END DO
      ff_max=MAX(f2_max,f3_max,f4_max,f5_max)
    END DO
    sum_ice=0.
    IF (del_t.GE.0.AND.ff_max.ne.0) THEN
      DO kr = 1,nkr
        arg_m = 0.0_wp
        DO ice = 1,icemax
          IF (ice ==1) THEN
            IF (kr .LE. 10) THEN
              arg_m = arg_m+ff2(kr,ice)
              ff2(kr,ice) = 0.0_wp
            ELSE IF (kr .GT. 10 .AND. kr .LT. 18) THEN
              meltrate = 0.5_wp/50.0_wp
              arg_m=arg_m+ff2(kr,ice)*(meltrate*dt)
              ff2(kr,ice)=ff2(kr,ice)-ff2(kr,ice)*(meltrate*dt)
            ELSE
              meltrate = 0.683_wp/120.0_wp
              arg_m=arg_m+ff2(kr,ice)*(meltrate*dt)
              ff2(kr,ice)=ff2(kr,ice)-ff2(kr,ice)*(meltrate*dt)
            END IF
          END IF
          IF (ice ==2 .OR. ice ==3) THEN
            IF (kr .LE. 12) THEN
              arg_m = arg_m+ff2(kr,ice)
              ff2(kr,ice)=0.0_wp
            ELSE IF (kr .GT. 12 .AND. kr .LT. 20) THEN
              meltrate = 0.5_wp/50.0_wp
              arg_m=arg_m+ff2(kr,ice)*(meltrate*dt)
              ff2(kr,ice)=ff2(kr,ice)-ff2(kr,ice)*(meltrate*dt)
            ELSE
              meltrate = 0.683_wp/120.0_wp
              arg_m=arg_m+ff2(kr,ice)*(meltrate*dt)
              ff2(kr,ice)=ff2(kr,ice)-ff2(kr,ice)*(meltrate*dt)
            END IF
          END IF
        END DO  ! DO ice
        !... snow
        IF (kr .LE. 14) THEN
          arg_m = arg_m + ff3(kr)
          ff3(kr) = 0.0_wp
        ELSE IF (kr .GT. 14 .AND. kr .LT. 22) THEN
          meltrate = (0.5_wp/50.0_wp)*tune_melt_factor
          arg_m=arg_m+ff3(kr)*(meltrate*dt)
          ff3(kr)=ff3(kr)-ff3(kr)*(meltrate*dt)
        ELSE
          meltrate = (0.683_wp/120.0_wp)*tune_melt_factor
          arg_m=arg_m+ff3(kr)*(meltrate*dt)
          ff3(kr)=ff3(kr)-ff3(kr)*(meltrate*dt)
        END IF
        ! ... graupel/hail
        IF (kr .LE. 13) THEN
          arg_m = arg_m+ff4(kr)+ff5(kr)
          ff4(kr)=0.0_wp
          ff5(kr)=0.0_wp
        ELSE IF (kr .GT. 13 .AND. kr .LT. 23) THEN
          meltrate = (0.5_wp/50.0_wp)*tune_melt_factor
          arg_m=arg_m+(ff4(kr)+ff5(kr))*(meltrate*dt)
          ff4(kr)=ff4(kr)-ff4(kr)*(meltrate*dt)
          ff5(kr)=ff5(kr)-ff5(kr)*(meltrate*dt)
        ELSE
          meltrate = (0.683_wp/120.0_wp)*tune_melt_factor
          arg_m=arg_m+(ff4(kr)+ff5(kr))*(meltrate*dt)
          ff4(kr)=ff4(kr)-ff4(kr)*(meltrate*dt)
          ff5(kr)=ff5(kr)-ff5(kr)*(meltrate*dt)
        END IF

        ff1(kr) = ff1(kr) + arg_m
        sum_ice=sum_ice+arg_m*3.0_wp*xl(kr)*xl(kr)*col
      END DO

      tin = tin - latheat_freez*sum_ice/ro ! cp/cv bug was taken into account (was -333.*sum_ice/ro)
    END IF

    RETURN
  END SUBROUTINE melting

  SUBROUTINE coll_add_pseodo_bin(x1p,g1p,nkf,x,g1)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: x(:), g1(:)
    REAL(KIND=wp),INTENT(INOUT) :: x1p(:),g1p(:)
    INTEGER,INTENT(INOUT) :: nkf
    INTEGER :: kr1

    DO kr1=1,nkr
      x1p(kr1)=x(kr1)
      g1p(kr1)=g1(kr1)
    END DO
    x1p(nkr+1)=2*x(nkr)
    g1p(nkr+1)=0.0_wp
    nkf=nkr+1

    RETURN
  END SUBROUTINE coll_add_pseodo_bin

  SUBROUTINE coll_add_pseodo_2bin(x1p,x2p,g1p,g2p,nkf,xa,xb,ga,gb)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: xa(:), xb(:), ga(:), gb(:)
    REAL(KIND=wp),INTENT(INOUT) :: x1p(:), x2p(:), g1p(:), g2p(:)
    INTEGER,INTENT(INOUT) :: nkf
    INTEGER :: kr1

    DO kr1=1,nkr
      x1p(kr1)=xa(kr1)
      x2p(kr1)=xb(kr1)
      g1p(kr1)=ga(kr1)
      g2p(kr1)=gb(kr1)
    END DO
    x1p(nkr+1)=2*xa(nkr)
    x2p(nkr+1)=2*xb(nkr)
    g1p(nkr+1)=0.0_wp
    g2p(nkr+1)=0.0_wp
    nkf=nkr+1

    RETURN
  END SUBROUTINE coll_add_pseodo_2bin

  SUBROUTINE coll_add_pseodo_3bin(x1p,x2p,g1p,g2p,g3p,nkf,xa,xb,ga,gb,gc)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: xa(:), xb(:), ga(:), gb(:), gc(:)
    REAL(KIND=wp),INTENT(INOUT) :: x1p(:), x2p(:), g1p(:), g2p(:), g3p(:)
    INTEGER,INTENT(INOUT) :: nkf
    INTEGER :: kr1

    DO kr1=1,nkr
      x1p(kr1)=xa(kr1)
      x2p(kr1)=xb(kr1)
      g1p(kr1)=ga(kr1)
      g2p(kr1)=gb(kr1)
      g3p(kr1)=gc(kr1)
    END DO
    x1p(nkr+1)=2*xa(nkr)
    x2p(nkr+1)=2*xb(nkr)
    g1p(nkr+1)=0.0_wp
    g2p(nkr+1)=0.0_wp
    g3p(nkr+1)=0.0_wp
    nkf=nkr+1

    RETURN
  END SUBROUTINE coll_add_pseodo_3bin

END MODULE mo_sbm_main

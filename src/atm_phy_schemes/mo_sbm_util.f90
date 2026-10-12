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

MODULE mo_sbm_util

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, txt => message_text
  USE mo_physical_constants, ONLY: &
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
       & rd         ! [J/K/kg] gas constant (dry air)
  USE mo_mpi,                ONLY: my_process_is_stdio, p_bcast, p_comm_work, p_io
  USE mo_io_units,           ONLY: find_next_free_unit
  USE mo_nwp_tuning_config,  ONLY: tune_sbmccn
  USE mo_model_domain,       ONLY: t_patch
  USE mo_nonhydro_types,     ONLY: t_nh_prog
  USE mo_impl_constants,     ONLY: min_rlcell
  USE mo_loopindices,        ONLY: get_indices_c
  USE mo_atm_phy_nwp_config, ONLY: atm_phy_nwp_config
  USE mo_2mom_mcrph_driver,  ONLY: two_moment_mcrph_init
  USE mo_2mom_mcrph_setup,   ONLY: cfg_params
  USE mo_run_config,         ONLY: iqv, iqc, iqi, iqr, iqs, iqg, iqh, iqbin, nkr => iqb_length, msg_level, &
       &                           iqb_water_start, iqb_water_end,      &
       &                           iqb_snow_start, iqb_snow_end,        &
       &                           iqb_graupel_start, iqb_graupel_end,  &
       &                           iqb_ccn_start, iqb_ccn_end

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: &
       & qx_from_bins_diag,        & !Update qc,qr,qi,qs,qg using the mass-bins of SBM microphysics
       & sbm_data, coll_kernel, ccn_init_sbm, sbm_init, & !local subroutines
       & ibreakup,                 & ! 1, 0/1 probability of drop breakup after collision
       & krmin_breakup,            & ! 31, minimum bin number for spontaneous breakup
       & jbreak,                   & ! 18, maximum target bin number for splinter after collision
       & snow_breakup_on,          & ! 1, 0/1 snow breakup
       & kr_snow_min,              & ! start of snow break up at R=~1cm: 30 out of 33
       & breakup_rain_spont_on,    & ! 0, 0/1 spontanaous rain breakup (effective for 33-36 bins, i.e. not relevant when nkr=33)
       & i_break_method,           & ! 1, 1/2: Srivastava 1971 / Kamra et al 1991
       & convert_micro2advect,     & ! T, T/F convert the units of PSD (chem_new) back to the model
       & melt_on,                  & ! 1, 0/1 jiwen fan simplified melting with melting rates depending on hydometeor type and size
       & hail_opt,                 & ! 1, 0/1 use the last PSD for graupel or hail
       & iceprocs,                 & ! 0/1 turn on ice processes
       & stick_param1,stick_param2,& ! 1/2 parameters for ice and snow sticking efficiency
       & icempl,                   & ! 1, 0/1 Hallet Mossop ice-multiplication
       & kr_icempl,                & ! 9, Hallet Mossop ice-multiplication: graupel+drops ->ice splinters only for drops>20mic
       & icemax,                   & ! 3, num of ice categories: here 3 ice categ. are nucleated, but then transfered to snow PSD
       & krice, krdrop,            & ! 18/15, ice/drop bin which accounts as snow/rain for output
       & nkr,                      & ! 33, number of bins at each PSD
       & col,                      & ! ln(2)/3, resolution of logaitmic double-mass grid
       & alcr, alcr_g,             & ! 0.5/inf, LWC thresholds (in gr/m^3) for resulting hydrometeor after collision
       & ncondcoll,                & ! 3/1, substepping of condensation/collisions
       & epsil,                    & ! epsilon. Represents small values in the code
       & isign_3point,             & ! 1, diff. growth remapping for 3 bins, conserving mass, concentration and radar reflectivity
       & coeff_remaping,           & ! 0.0066667_wp, proximity to current bin, which cancels rempapping to distant bins
       & ventpl_max,               & ! 5, maximum ventilation coefficient (due to drop fall) in condevap_dmdt_coef subroutine
       & rw_pw_min,                & ! epsilon in condevap_supsat_eqn subroutine for water
       & rw_pw_ri_pi_min,          & ! epsilon in condevap_supsat_eqn subroutine for ice
       & ratio_icew_min,           & ! epsilon in condevap_supsat_eqn subroutine
       & use_cloud_base_nuc,       & ! 0, max supsat calculation method
       & t_nucl_drop_min,          & ! minimum temperature for water nucleation
       & t_nucl_ice_min,           & ! minimum temperature for inhomogeneous nucleation
       & isign_tq_icenucl,         & ! 1, 0/1 flag to update temperaute and mixing ratio after ice nucleation
       & delsupice_max,            & ! 59%, max limit of supersaturation during ice nucleation (Meyers formula)
       & xl,xi,xs,xg,xh,           & ! mass bin grids (center, not boundaries)
       & xl_mg,xs_mg,xg_mg,xh_mg,  & ! mass bin grids (center, not boundaries)
       & rlec, riec, rsec, rgec, rhec,        & ! electrostatic capacitence needed for diffusional growth
       & vr1, vr2, vr3, vr4, vr5,             & ! terminal velocities
       & ywll_1000mb, ywll_750mb, ywll_500mb, & ! collision kernals at 3 heights. Used only in utils
       & ro1bl, ro2bl, ro3bl, ro4bl, ro5bl,   & ! bulk density of hydrometeors as function of their mass
       & ima,                                 & ! target bin after collision (remaping done between ima and ima+1)
       & chucm,                               & ! "courant number" in coll_bott_remap subroutine
       & ecoalmassm,                          & ! coallessence efficiency of drops
       & prob, gain_var_new, nnd,             & ! probabilities in breakup_rain_spont_prob subroutine
       & dropradii,                           & ! radiuses of drop-bins
       & pkij, qkj,                           & ! breakup (bleck's 1st order): gain and loss coefficients
       & ikr_spon_break,                      & ! for drops size < ikr_spon_break - no spontaneous breakup
       & coll_turb_fact,                      & ! artificial increase of collision kernel due to turbulence
       & cwll_all,                            & ! collision kernel liquid-liquid
       & cwlg_all, cwlh_all, cwls_all,        & ! collision kernels water with g,h,s
       & cwgl_all, cwhl_all, cwsl_all,        & ! collision kernels g,h,s with water
       & cwsg_all, cwss_all,                  & ! collision kernels s-g, s-s
       & cwgs_all, cwsh_all, cwhs_all, cwgg_all, cwhh_all, & ! artificially added on 20240830, in future: read tables for these kerenels
       & p_z_down, p_z_up, p_z_del,           & ! 1000 and 300mb for kernels veritcal interpolation
       & fccnr_mar, fccnr_con,                & ! 2 theoretical options for initial aerosol distribution
       & rccn,                                & ! ccn bins radii (cm)
       & fccnr_obs,                           & ! empirical option for initial aerosol distribution: dndlogd for 33bins
       & ilognormal_modes_aerosol,            & ! 1, 0/1 for empirical/lognormal initial aerosol distributions
       & mwaero,                              & ! molecular mass of aerosol (we use NaCl here)
       & ions,                                & ! number of ions after solution (2 for NaCl)
       & ro_solute,                           & ! density of solid salt aerosol
       & z0in, zmin,                          & ! parameters for exp decrease of initial aerosol concentration with height (cm)
       & rhoh2o,                              & ! [kg/m3]  density of liquid water
       & alv,                                 & ! [J/kg]   latent heat of vaporization
       & als,                                 & ! [J/kg]   latent heat of sublimation
       & alf,                                 & ! [J/kg]   latent heat for fusion (freezing)
       & cpd,                                 & ! [J/K/kg] specific heat at constant pressure
       & cvd,                                 & ! [J/K/kg] specific heat at constant volume
       & rhoice,                              & ! [kg/m3]  density of pure ice
       & tmelt,                               & ! [K]      melting temperature of ice/snow
       & rv,                                  & ! [J/K/kg] gas constant for water vapor
       & grav,                                & ! [m/s2] av. gravitational acceleration
       & rd,                                  & ! [J/K/kg] gas constant (dry air)
       & ccnconstarr, use_ccn_const,          & ! option to use time-constant ccn profile i.e. initialize it  every micro. time step
       & latheatfac1,latheatfac2,latheatfac3,latheatfac4, & !uncertain but less important switches for cp/cv bug
       & tune_fall, tune_long_relax,          & ! tuning factors
       & usenkrp1a,usenkrp1b,usenkrp2,usenkrp3, & ! optional parameters for 34th pseodo bin
       & tune_melt_factor,                    & !artificial factor for increasing melt rate of snow, graupel and hail
       & snha2ha,                             & ! snow+hail/gra. collisions. 0: snow+hail/gra. --> snow. 1: snow+hail/gra. --> hail/gra.
       & positive_t_coll                        ! 1/0: allow / dont allow collisions with solid particles also at positive temperatures

 CHARACTER(LEN=*), PARAMETER :: modname = 'mo_sbm_util'

 !--------- IMPORTANT PARAMETERS ---------------! start
 INTEGER,PARAMETER :: iceprocs = 1           ! 1/0 allow/dont allow ice processes
 INTEGER,PARAMETER :: icempl=1,icemax=3,kr_icempl=9,krice=18,krdrop=15,jbreak=18,krmin_breakup=31   ! krdrop=bin 15 --> 50um
 INTEGER,PARAMETER :: ibreakup=0             ! probability of drop breakup after collision. 1 is not working and probably not needed for 33 bins
 INTEGER,PARAMETER :: snow_breakup_on=1      ! 1, 0/1 snow breakup
 INTEGER,PARAMETER :: kr_snow_min=31         ! start of snow break up at R=~1cm: 30 out of 33
 INTEGER,PARAMETER :: breakup_rain_spont_on=0 ! 0, 0/1 spontanaous rain breakup (effective for 33-36 bins, i.e. not relevant when nkr=33)
 INTEGER,PARAMETER :: i_break_method=1       ! 1/2: Srivastava 1971 / Kamra et al 1991
 LOGICAL, PARAMETER :: latheatfac1=.TRUE., latheatfac2=.TRUE., latheatfac3=.FALSE., latheatfac4=.TRUE. !latent heat release
                                                                                     ! factors related to cp/cv bug in ICON
 INTEGER,PARAMETER :: hail_opt=1             ! 0: use graupel PSD. 1: use hail PSD instead. ICON output is written to graupel anyway
 REAL(KIND=wp), PARAMETER :: alcr=0.5        ! gr/m^3. LWC<alcr: snow+water-->snow, LWC>alcr: snow+water-->graupel
 REAL(KIND=wp), PARAMETER :: alcr_g=100.0    ! forcing no transition from graupel to hail in this version
 INTEGER,PARAMETER :: ncondcoll=3            ! number of substeps for diffusional growth / collisions
 INTEGER,PARAMETER :: coll_turb_fact=1       ! 0/1: artificial increase of collision kernel due to turbulence
                                             ! In future real turbulence in ICON should be taken into account
 INTEGER,PARAMETER :: use_ccn_const=0        ! 1: use time-constant ccn profile (initialize it again every microphysical time step)
 REAL(KIND=wp), PARAMETER :: tune_fall=1.0   ! >1.0 artificial increase factor for graupel and hail fall speed
 REAL(KIND=wp), PARAMETER :: tune_long_relax=1.0 ! >1.0 artificial increase factor for relaxation time
 INTEGER,PARAMETER :: usenkrp1a=0,usenkrp1b=0,usenkrp2=0,usenkrp3=0 ! parameters for 34th pseodo bin. 0: 33 bins, 1: use pseodo 34 bin
 REAL(KIND=wp), PARAMETER :: tune_melt_factor=4.0 !default: 1.0, artificial factor for increasing melt rate of snow, graupel and hail
 INTEGER,PARAMETER :: stick_param1=1         ! 1/2 parameters for collision efficiency which includes ice, was 2
 INTEGER,PARAMETER :: stick_param2=1         ! 1/2 parameters for collision efficiency which includes ice, was 1
 INTEGER,PARAMETER :: snha2ha=1              ! snow+hail/graupel collisions. 0: snow+hail/graupel --> snow. 1: snow+hail/graupel --> hail/graupel
 INTEGER,PARAMETER :: positive_t_coll=1      ! 1/0: allow / dont allow collisions with solid particles also at positive temperatures
 !--------- IMPORTANT PARAMETERS ---------------! end

 LOGICAL, PARAMETER :: convert_micro2advect=.true.    ! must be true. false means ignore PSD changes
 INTEGER, PARAMETER :: melt_on=1             ! turn on(1) or off(0) melting
 REAL(KIND=wp), PARAMETER :: z0in=2.0e5, zmin=2.0e5   ! parameters for exp decrease of initial aerosol concentration with height (cm)
 REAL(KIND=wp), PARAMETER :: col=0.23105_wp  ! ln(2)/3: needed to convert 2^n mass grid to linear vs ln(raius)
 INTEGER,PARAMETER :: p_z_down=1050000, p_z_up=20000, p_z_del=10000 ! 1050 and 20mb for kernels veritcal interpolation, 10mb step
 REAL(KIND=wp), PARAMETER :: epsil=1.0e-16   ! negligible value
 INTEGER,PARAMETER :: isign_3point=1         ! diff. growth remapping for 3 bins, conserving mass, concentration and radar reflectivity
 REAL(KIND=wp), PARAMETER::coeff_remaping=0.0066667_wp !remapping limiter to decrease PSD numerical broadening
 REAL(KIND=wp), PARAMETER::ventpl_max=5.0_wp       ! maximum ventilation coefficient (due to drop fall) in condevap_dmdt_coef subroutine
 REAL(KIND=wp), PARAMETER::rw_pw_min=1.0e-10       ! epsilon in condevap_supsat_eqn subroutine for water
 REAL(KIND=wp), PARAMETER::rw_pw_ri_pi_min=1.0e-10 ! epsilon in condevap_supsat_eqn subroutine for ice
 REAL(KIND=wp), PARAMETER::ratio_icew_min=1.0e-4   ! epsilon in condevap_supsat_eqn subroutine
 INTEGER,PARAMETER :: use_cloud_base_nuc=0   ! max supsat calculation method:
                                             ! 0: simple use of local supsat as in warm code
                                             ! 1: find max supsat (~20m) above c.base for better
                                             ! estimation of nucleated drops concentration. It is currently off
                                             ! in the code since rain existence at cloud base at later stages
                                             ! can consume supersaturation, and the typical supsat maximum
                                             ! above cloud base does not exist
 REAL(KIND=wp), PARAMETER::t_nucl_drop_min=-80.0_wp   ! minimum temperature for water nucleation
 REAL(KIND=wp), PARAMETER::t_nucl_ice_min=-38.0_wp    ! minimum temperature for inhomogeneous nucleation
 INTEGER,PARAMETER :: isign_tq_icenucl=1     ! 1, 0/1 flag to update temperaute and mixing ratio after ice nucleation
 REAL(KIND=wp), PARAMETER::delsupice_max=59.0_wp   !clipping of supsat over ice for nucleation to delsupice_max<59% (Meyers formula)
 REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: xl_mg,xs_mg,xg_mg,xh_mg, & ! mass bins in mg
                                             xl,xs,xg,xh,             & ! mass bins
                                             rlec,rsec,rgec,rhec,     & ! capacities of particles during diffusional growth
                                             vr1,vr3,vr4,vr5            ! fall velocoties
 REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: riec,xi,vr2
 REAL(KIND=wp), ALLOCATABLE :: &             ! collision kernals at 3 heights
                         ywll_1000mb(:,:),ywll_750mb(:,:),ywll_500mb(:,:), &
                         ywlg_300mb(:,:),ywlg_500mb(:,:),ywlg_750mb(:,:),  &
                         ywlh_300mb(:,:),ywlh_500mb(:,:),ywlh_750mb(:,:),  &
                         ywls_300mb(:,:),ywls_500mb(:,:),ywls_750mb(:,:),  &
                         ywsg_300mb(:,:),ywsg_500mb(:,:),ywsg_750mb(:,:),  &
                         ywss_300mb(:,:),ywss_500mb(:,:),ywss_750mb(:,:),  &
                         cwlg_all(:,:,:),cwlh_all(:,:,:),cwls_all(:,:,:),  &
                         cwgl_all(:,:,:),cwhl_all(:,:,:),cwsl_all(:,:,:),  &
                         cwll_all(:,:,:),cwsg_all(:,:,:),cwss_all(:,:,:),  &
                         !artificially added (tables missing):
                         cwsh_all(:,:,:),cwhs_all(:,:,:),cwgs_all(:,:,:),cwgg_all(:,:,:),cwhh_all(:,:,:)
 REAL(KIND=wp), ALLOCATABLE :: ro1bl(:), ro2bl(:,:), ro3bl(:), ro4bl(:), ro5bl(:), ccnconstarr(:,:), &
                         chucm(:,:), ecoalmassm(:,:), prob(:),gain_var_new(:,:),nnd(:,:), &
                         dropradii(:),pkij(:,:,:),qkj(:,:), fccnr_mar(:),fccnr_con(:), fccnr_obs(:), &
                         rccn(:)
 INTEGER ::              ikr_spon_break
 INTEGER,ALLOCATABLE ::  ima(:,:)
 INTEGER,PARAMETER :: ilognormal_modes_aerosol = 1 !follow lognormal distribution for aerosol size distribution
                                                   ! in case of ilognormal_modes_aerosol = 0 ! read in a sd file from observation.
                                                   ! Currently the file name for the observed
                                                   ! sd is "ccn_size_33bin.dat", which is from the july 18 2017 "ENA" case.
 INTEGER,PARAMETER :: ions = 2                ! sea salt. Set ions = 3 for ammonium-sulfate
 REAL(KIND=wp), PARAMETER :: mwaero = 22.9 + 35.5 ! sea salt. ! mwaero = 115.0
 REAL(KIND=wp), PARAMETER :: ro_solute = 2.16 ! sea salt. Set ro_solute = 1.79 for ammonium-sulfate
 !-----------------------------------------------------------
 ! following are two type declarations to hold the values of equidistand lookup tables. Similar lookup tables are
 ! used in the segal-khain parameterization of ccn-activation and for determining the wet growth diameter of graupel.

 CONTAINS

  SUBROUTINE sbm_data(fccnr_con,fccnr_mar,fccnr_obs)
    IMPLICIT NONE
    INTEGER :: unitnr,i,j,kr,error
    REAL(KIND=wp) :: dlnr, ccnr(nkr)
    CHARACTER(LEN=256), PARAMETER :: dir_43 = "SBM_input_43", dir_33 = "SBM_input_33"
    CHARACTER(LEN=256) :: input_dir
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::sbm_data'
    REAL(KIND=wp) ,INTENT(INOUT) :: fccnr_con(:), fccnr_mar(:), fccnr_obs(:)

    IF (nkr == 33) input_dir = TRIM(dir_33)
    IF (nkr == 43) input_dir = TRIM(dir_43)
    CALL message(routine," fast sbm: initializing hujisbm ")
    CALL message(routine," fast sbm: ****** hujisbm ******* ")

    dlnr=LOG(2._wp)/(3._wp)
    ! +----------------------------------------------+
    ! lookuptable #2: electrostatic capacitence needed for
    ! diffusional growth. For spherical particle it =R, for others f(R)
    ! +----------------------------------------------+
    IF (.NOT. ALLOCATED(rlec)) ALLOCATE(rlec(nkr))
    IF (.NOT. ALLOCATED(riec)) ALLOCATE(riec(nkr,icemax))
    IF (.NOT. ALLOCATED(rsec)) ALLOCATE(rsec(nkr))
    IF (.NOT. ALLOCATED(rgec)) ALLOCATE(rgec(nkr))
    IF (.NOT. ALLOCATED(rhec)) ALLOCATE(rhec(nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'sbm_data : table-2 -- opening capacity33.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/capacity33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'sbm_data : table-2 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(modname,'error in sbm_data')
      END IF
      READ (unitnr,*) rlec,riec,rsec,rgec,rhec
    END IF

    CALL p_bcast(rlec, p_io, p_comm_work)
    CALL p_bcast(riec, p_io, p_comm_work)
    CALL p_bcast(rsec, p_io, p_comm_work)
    CALL p_bcast(rgec, p_io, p_comm_work)
    CALL p_bcast(rhec, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-2'
    CALL message(routine,TRIM(txt))
    ! +-----------------------------------------------+
    ! lookuptable #3: bin masses array for each hydrometeor
    ! +-----------------------------------------------+
    IF (.NOT. ALLOCATED(xl)) ALLOCATE(xl(nkr))
    IF (.NOT. ALLOCATED(xi)) ALLOCATE(xi(nkr,icemax))
    IF (.NOT. ALLOCATED(xs)) ALLOCATE(xs(nkr))
    IF (.NOT. ALLOCATED(xg)) ALLOCATE(xg(nkr))
    IF (.NOT. ALLOCATED(xh)) ALLOCATE(xh(nkr))
    IF (.NOT. ALLOCATED(xl_mg)) ALLOCATE(xl_mg(nkr))
    IF (.NOT. ALLOCATED(xs_mg)) ALLOCATE(xs_mg(nkr))
    IF (.NOT. ALLOCATED(xg_mg)) ALLOCATE(xg_mg(nkr))
    IF (.NOT. ALLOCATED(xh_mg)) ALLOCATE(xh_mg(nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'sbm_data : table-3 -- opening masses33.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/masses33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'sbm_data : table-3 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(modname,'error in sbm_data')
      END IF
      READ (unitnr,*) xl,xi,xs,xg,xh
    END IF

    CALL p_bcast(xl, p_io, p_comm_work)
    CALL p_bcast(xi, p_io, p_comm_work)
    CALL p_bcast(xs, p_io, p_comm_work)
    CALL p_bcast(xg, p_io, p_comm_work)
    CALL p_bcast(xh, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-3'
    CALL message(routine,TRIM(txt))
    ! +---------------------------------------------------+
    ! lookuptable #4: fall velocity as function mass for different hydrometeors, assuming 1000mb
    ! +---------------------------------------------------+
    IF (.NOT. ALLOCATED(vr1)) ALLOCATE(vr1(nkr))
    IF (.NOT. ALLOCATED(vr2)) ALLOCATE(vr2(nkr,icemax))
    IF (.NOT. ALLOCATED(vr3)) ALLOCATE(vr3(nkr))
    IF (.NOT. ALLOCATED(vr4)) ALLOCATE(vr4(nkr))
    IF (.NOT. ALLOCATED(vr5)) ALLOCATE(vr5(nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'sbm_data : table-4 -- opening termvels.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/termvels33_corrected.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'sbm_data : table-4 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(modname,'error in sbm_data')
      END IF
      READ (unitnr,*) vr1,vr2,vr3,vr4,vr5
    END IF

    CALL p_bcast(vr1, p_io, p_comm_work)
    CALL p_bcast(vr2, p_io, p_comm_work)
    CALL p_bcast(vr3, p_io, p_comm_work)
    CALL p_bcast(vr4, p_io, p_comm_work)
    CALL p_bcast(vr5, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-4'
    CALL message(routine,TRIM(txt))

    ! +-----------------------------------------------------------------+
    ! lookuptable #6: collision kernels depending on pressure: collision
    ! probability of 2 particles (not including sticking efficiency): drops-drops
    ! +-----------------------------------------------------------------+
    IF (.NOT. ALLOCATED(ywll_1000mb)) ALLOCATE(ywll_1000mb(nkr,nkr))
    IF (.NOT. ALLOCATED(ywll_750mb)) ALLOCATE(ywll_750mb(nkr,nkr))
    IF (.NOT. ALLOCATED(ywll_500mb)) ALLOCATE(ywll_500mb(nkr,nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'sbm_data : table-6 -- opening kernels_z.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/kernLL_z33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'sbm_data : table-6 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(modname,'error in sbm_data')
      END IF
      READ (unitnr,*) ywll_1000mb,ywll_750mb,ywll_500mb
    END IF
    DO i=1,nkr
       DO j=1,nkr
          IF (i > nkr .or. j > nkr) THEN
             ywll_1000mb(i,j) = 0.0
             ywll_750mb(i,j) =  0.0
             ywll_500mb(i,j) =  0.0
          END IF
       END DO
    END DO

    CALL p_bcast(ywll_1000mb, p_io, p_comm_work)
    CALL p_bcast(ywll_750mb, p_io, p_comm_work)
    CALL p_bcast(ywll_500mb, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-6'
    CALL message(routine,TRIM(txt))

    IF ( iceprocs == 1 ) THEN
      ! +-----------------------------------------------------------------------+
      ! lookuptable #7
      ! collisions kernels: collision probability of 2 particles (not including sticking efficiency):
      ! +-----------------------------------------------------------------------+
      ! ... drops - graupel
      IF (.NOT. ALLOCATED(ywlg_300mb)) ALLOCATE(ywlg_300mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywlg_500mb)) ALLOCATE(ywlg_500mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywlg_750mb)) ALLOCATE(ywlg_750mb(nkr,nkr))
      ! ... drops - hail
      IF (.NOT. ALLOCATED(ywlh_300mb)) ALLOCATE(ywlh_300mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywlh_500mb)) ALLOCATE(ywlh_500mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywlh_750mb)) ALLOCATE(ywlh_750mb(nkr,nkr))
      ! ... drops - snow
      IF (.NOT. ALLOCATED(ywls_300mb)) ALLOCATE(ywls_300mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywls_500mb)) ALLOCATE(ywls_500mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywls_750mb)) ALLOCATE(ywls_750mb(nkr,nkr))
      ! ... snow - graupel
      IF (.NOT. ALLOCATED(ywsg_300mb)) ALLOCATE(ywsg_300mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywsg_500mb)) ALLOCATE(ywsg_500mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywsg_750mb)) ALLOCATE(ywsg_750mb(nkr,nkr))
      ! ... snow - snow
      IF (.NOT. ALLOCATED(ywss_300mb)) ALLOCATE(ywss_300mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywss_500mb)) ALLOCATE(ywss_500mb(nkr,nkr))
      IF (.NOT. ALLOCATED(ywss_750mb)) ALLOCATE(ywss_750mb(nkr,nkr))

      ! ... kernels depending on pressure :

      ! ... drops - graupel
      IF (my_process_is_stdio()) THEN
        unitnr=find_next_free_unit(10,999)
        WRITE(txt, '(a,i2)') 'sbm_data : table-7 -- opening cklg 300mb 500mb 750mb files'
        CALL message(modname,trim(txt))
        OPEN(unitnr, file=TRIM(input_dir)//"/cklg_33_300mb_500mb_750mb.asc", status='old', form='formatted', iostat=error)
        IF (error /= 0) THEN
          WRITE (txt,*) 'sbm_data : table-7 cklg_33_300mb_500mb_750mb.asc not found'
          CALL message(modname,TRIM(txt))
          CALL finish(modname,'error in sbm_data')
        END IF
        READ (unitnr,*) ywlg_300mb,ywlg_500mb,ywlg_750mb
      END IF

      ! ... drops - hail
      IF (my_process_is_stdio()) THEN
        unitnr=find_next_free_unit(10,999)
        WRITE(txt, '(a,i2)') 'sbm_data : table-7 -- opening cklh 300mb 500mb 750mb files'
        CALL message(modname,trim(txt))
        OPEN(unitnr, file=TRIM(input_dir)//"/cklh_33_300mb_500mb_750mb.asc", status='old', form='formatted', iostat=error)
        IF (error /= 0) THEN
          WRITE (txt,*) 'sbm_data : table-7 cklh_33_300mb_500mb_750mb.asc not found'
          CALL message(modname,TRIM(txt))
          CALL finish(modname,'error in sbm_data')
        END IF
        READ (unitnr,*) ywlh_300mb,ywlh_500mb,ywlh_750mb
      END IF

      ! ... drops - snow
      IF (my_process_is_stdio()) THEN
        unitnr=find_next_free_unit(10,999)
        WRITE(txt, '(a,i2)') 'sbm_data : table-7 -- opening ckls 300mb 500mb 750mb files'
        CALL message(modname,trim(txt))
        OPEN(unitnr, file=TRIM(input_dir)//"/ckls_33_300mb_500mb_750mb.asc", status='old', form='formatted', iostat=error)
        IF (error /= 0) THEN
          WRITE (txt,*) 'sbm_data : table-7 ckls_33_300mb_500mb_750mb.asc not found'
          CALL message(modname,TRIM(txt))
          CALL finish(modname,'error in sbm_data')
        END IF
        READ (unitnr,*) ywls_300mb,ywls_500mb,ywls_750mb
      END IF

      ! ... snow - graupel
      IF (my_process_is_stdio()) THEN
        unitnr=find_next_free_unit(10,999)
        WRITE(txt, '(a,i2)') 'sbm_data : table-7 -- opening cksg 300mb 500mb 750mb files'
        CALL message(modname,trim(txt))
        OPEN(unitnr, file=TRIM(input_dir)//"/cksg_33_300mb_500mb_750mb.asc", status='old', form='formatted', iostat=error)
        IF (error /= 0) THEN
          WRITE (txt,*) 'sbm_data : table-7 cksg_33_300mb_500mb_750mb.asc not found'
          CALL message(modname,TRIM(txt))
          CALL finish(modname,'error in sbm_data')
        END IF
        READ (unitnr,*) ywsg_300mb,ywsg_500mb,ywsg_750mb
      END IF

      ! ... snow - snow
      IF (my_process_is_stdio()) THEN
        unitnr=find_next_free_unit(10,999)
        WRITE(txt, '(a,i2)') 'sbm_data : table-7 -- opening ckss 300mb 500mb 750mb files'
        CALL message(modname,trim(txt))
        OPEN(unitnr, file=TRIM(input_dir)//"/ckss_33_300mb_500mb_750mb.asc", status='old', form='formatted', iostat=error)
        IF (error /= 0) THEN
          WRITE (txt,*) 'sbm_data : table-7 ckss_33_300mb_500mb_750mb.asc not found'
          CALL message(modname,TRIM(txt))
          CALL finish(modname,'error in sbm_data')
        END IF
        READ (unitnr,*) ywss_300mb,ywss_500mb,ywss_750mb
      END IF

      CALL p_bcast(ywlg_300mb, p_io, p_comm_work)
      CALL p_bcast(ywlg_500mb, p_io, p_comm_work)
      CALL p_bcast(ywlg_750mb, p_io, p_comm_work)

      CALL p_bcast(ywlh_300mb, p_io, p_comm_work)
      CALL p_bcast(ywlh_500mb, p_io, p_comm_work)
      CALL p_bcast(ywlh_750mb, p_io, p_comm_work)

      CALL p_bcast(ywls_300mb, p_io, p_comm_work)
      CALL p_bcast(ywls_500mb, p_io, p_comm_work)
      CALL p_bcast(ywls_750mb, p_io, p_comm_work)

      CALL p_bcast(ywsg_300mb, p_io, p_comm_work)
      CALL p_bcast(ywsg_500mb, p_io, p_comm_work)
      CALL p_bcast(ywsg_750mb, p_io, p_comm_work)

      CALL p_bcast(ywss_300mb, p_io, p_comm_work)
      CALL p_bcast(ywss_500mb, p_io, p_comm_work)
      CALL p_bcast(ywss_750mb, p_io, p_comm_work)

      WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-7'
      CALL message(routine,trim(txt))
      ! +-----------------------------------------------------------------------+
    END IF

    ! +--------------------------------------------------------------+
    ! lookuptable #8: bulk density of hydrometeors as function of their mass
    ! +--------------------------------------------------------------+
    IF (.NOT. ALLOCATED(ro1bl)) ALLOCATE(ro1bl(nkr))
    IF (.NOT. ALLOCATED(ro2bl)) ALLOCATE(ro2bl(nkr,icemax))
    IF (.NOT. ALLOCATED(ro3bl)) ALLOCATE(ro3bl(nkr))
    IF (.NOT. ALLOCATED(ro4bl)) ALLOCATE(ro4bl(nkr))
    IF (.NOT. ALLOCATED(ro5bl)) ALLOCATE(ro5bl(nkr))

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'sbm_data : table-8 -- opening bulkdens.asc'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(input_dir)//"/bulkdens33.asc", status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'sbm_data : table-8 not found'
        CALL message(modname,TRIM(txt))
        CALL finish(modname,'error in sbm_data')
      END IF
      READ (unitnr,*) ro1bl,ro2bl,ro3bl,ro4bl,ro5bl
    END IF

    CALL p_bcast(ro1bl, p_io, p_comm_work)
    CALL p_bcast(ro2bl, p_io, p_comm_work)
    CALL p_bcast(ro3bl, p_io, p_comm_work)
    CALL p_bcast(ro4bl, p_io, p_comm_work)
    CALL p_bcast(ro5bl, p_io, p_comm_work)

    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading table-8'
    CALL message(routine,TRIM(txt))

    DO i=1,nkr
      xl_mg(i) = xl(i)*1.e3
      xs_mg(i) = xs(i)*1.e3
      xg_mg(i) = xg(i)*1.e3
      xh_mg(i) = xh(i)*1.e3
    END DO

    IF (.NOT. ALLOCATED(ima)) ALLOCATE(ima(nkr,nkr))
    IF (.NOT. ALLOCATED(chucm)) ALLOCATE(chucm(nkr,nkr))
    chucm  = 0.0_wp
    ima = 0
    CALL coll_bott_remap(xl, chucm, ima)
    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading "coll_bott_remap" '
    CALL message(routine,TRIM(txt))

    IF (.NOT. ALLOCATED(dropradii)) ALLOCATE(dropradii(nkr))
    DO kr=1,nkr
       dropradii(kr)=(3.0*xl(kr)/4.0/3.141593/1.0)**(1.0_wp/3.0_wp)   ! [cm]
    END DO

    ! +-------------------------+
    ! allocating aerosols array
    ! +-------------------------+
    IF (.NOT. ALLOCATED(rccn)) ALLOCATE(rccn(nkr))

    ! ... initializing the fccnr_mar and fccnr_con
    fccnr_con = 0.0
    fccnr_mar = 0.0
    fccnr_obs = 0.0
    rccn = 0.0

    IF (ilognormal_modes_aerosol == 1) THEN
       CALL lognormal_modes_aerosol(fccnr_con,fccnr_mar,col,xl,rccn,1) ! Defining fccnr_mar.
       CALL lognormal_modes_aerosol(fccnr_con,fccnr_mar,col,xl,rccn,2) ! Defining fccnr_con
       WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull definition of "lognormal_modes_aerosol" '
       CALL message(routine,TRIM(txt))
    ELSE ! read an observed sd with a format of aerosol size (cm), dn (#cm-3) and dndlogd for 33bins (jinwe fan)
       IF (my_process_is_stdio()) THEN
         unitnr=find_next_free_unit(10,999)
         WRITE(txt, '(a,i2)') 'sbm_data : opening CCN_size_33bin.dat'
         CALL message(modname,TRIM(txt))
         OPEN(unitnr, file=TRIM(input_dir)//"/CCN_size_33bin.dat", status='old', form='formatted', iostat=error)
         IF (error /= 0) THEN
           WRITE (txt,*) 'sbm_data : CCN_size_33bin.dat not found'
           CALL message(modname,TRIM(txt))
           CALL finish(modname,'error in sbm_data')
         END IF
         DO kr=1,nkr
           READ (unitnr,*) rccn(kr),ccnr(kr),fccnr_obs(kr) !---aerosol size (cm), dn (# cm-3) and dndlogd for 33bins
         END DO
         CALL message(routine,"fast_sbm_init: succesfull reading aerosol sd from observation")
       END IF
    END IF

    IF (.NOT. ALLOCATED(pkij)) ALLOCATE(pkij(jbreak,jbreak,jbreak))
    IF (.NOT. ALLOCATED(qkj)) ALLOCATE(qkj(jbreak,jbreak))
    IF (.NOT. ALLOCATED(ecoalmassm)) ALLOCATE(ecoalmassm(nkr,nkr))
    pkij = 0.0_wp
    qkj = 0.0_wp
    ecoalmassm = 0.0_wp
    ! rain collisional breakup probability, performed once in a run. Output: coallessence efficiency: ecoalmassm:
    CALL breakup_rain_prob(pkij,qkj,ecoalmassm,xl,jbreak,vr1)

    CALL p_bcast(pkij, p_io, p_comm_work)
    CALL p_bcast(qkj, p_io, p_comm_work)
    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading breakup_rain_prob" '
    CALL message(routine,TRIM(txt))

    ! The following collision kerenls are defined here and will be interpolated in coll_kerenl given the current pressure p_z:
    IF (.NOT. ALLOCATED(cwll_all)) ALLOCATE(cwll_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwlg_all)) ALLOCATE(cwlg_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwlh_all)) ALLOCATE(cwlh_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwls_all)) ALLOCATE(cwls_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwgl_all)) ALLOCATE(cwgl_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwhl_all)) ALLOCATE(cwhl_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwsl_all)) ALLOCATE(cwsl_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwsg_all)) ALLOCATE(cwsg_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwss_all)) ALLOCATE(cwss_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    ! artificially added:
    IF (.NOT. ALLOCATED(cwsh_all)) ALLOCATE(cwsh_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwhs_all)) ALLOCATE(cwhs_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwgs_all)) ALLOCATE(cwgs_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwgg_all)) ALLOCATE(cwgg_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))
    IF (.NOT. ALLOCATED(cwhh_all)) ALLOCATE(cwhh_all(nkr,nkr,(p_z_down-p_z_up)/p_z_del+1))

    cwll_all(:,:,:) = 0.0_wp
    cwlg_all(:,:,:) = 0.0_wp ; cwlh_all(:,:,:) = 0.0_wp ; cwls_all(:,:,:) = 0.0_wp
    cwgl_all(:,:,:) = 0.0_wp ; cwhl_all(:,:,:) = 0.0_wp ; cwsl_all(:,:,:) = 0.0_wp
    cwsg_all(:,:,:) = 0.0_wp ; cwss_all(:,:,:) = 0.0_wp
    ! artificially added on 20240830 (tables missing):
    cwsh_all(:,:,:) = 0.0_wp
    cwhs_all(:,:,:) = 0.0_wp
    cwgs_all(:,:,:) = 0.0_wp
    cwgg_all(:,:,:) = 0.0_wp
    cwhh_all(:,:,:) = 0.0_wp

    !+---+-----------------------------------------+
    IF (.NOT. ALLOCATED(prob)) ALLOCATE(prob(nkr))
    IF (.NOT. ALLOCATED(gain_var_new)) ALLOCATE(gain_var_new(nkr,nkr))
    IF (.NOT. ALLOCATED(nnd)) ALLOCATE(nnd(nkr,nkr))
    prob = 0.0
    gain_var_new = 0.0
    nnd = 0.0
    ! calculation of breakup probabilities, performed once in a run, and used later in breakup_rain_spont (missing in warm code):
    CALL breakup_rain_spont_prob(xl, prob, gain_var_new, nnd, ikr_spon_break) ! missing in warm_sbm, inform Koby
    WRITE(txt, '(a,i2)') 'fast_sbm_init : succesfull reading "breakup_rain_spont_prob" '

    RETURN
  END SUBROUTINE sbm_data

  SUBROUTINE coll_kernel(dtime_coal) ! The intepolation of kernels is performed now only at initialization (inform Koby)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: dtime_coal
    INTEGER :: i,j,cwturb,p_zi,p_z
    REAL(KIND=wp), PARAMETER :: p1=1.0e6,p2=0.75e6,p3=0.50e6,p4=0.3e6    !1000,750,500,300mb
    REAL(KIND=wp) :: dtimelnr

    dtimelnr = dtime_coal*LOG(2.0_wp)/(3.0_wp)

    ! vertical interpolation of collision kernal of water-water using kernals at 3 pressure levels (1000,750,500)
    ! performed within coll_kernel_interp subroutine
    DO p_zi=1,(p_z_down-p_z_up)/p_z_del+1 !1050 and 10mb for kernels veritcal interpolation, 1mb step
      p_z=p_z_up+(p_zi-1)*p_z_del
      DO i=1,nkr     ! bin number of hydrometeor 1
        DO j=1,nkr   ! bin number of hydrometeor 2
          ! water - water
          cwll_all(i,j,p_zi) = coll_kernel_interp(p_z,p1,p2,p3,ywll_1000mb(i,j),ywll_750mb(i,j),ywll_500mb(i,j))*dtimelnr
        END DO
      END DO

      ! ... ecoalmassm is from "breakiniit_ks"
      DO i=1,nkr
        DO j=1,nkr
          cwll_all(i,j,p_zi) = ecoalmassm(i,j)*cwll_all(i,j,p_zi) !sticking efficiency for water, reduced by ecoalmassm factor
        END DO
      END DO

      IF ( iceprocs == 1 ) THEN
        ! vertical interpolation of collision kernal of ice-ice&ice-water using kernals at 3 pressure levels (750,500,300mb)
        ! p_1=750mb is first pressure level band, p_3=300mb is last pressure level band, p_z is running pressure level

        DO i=1,nkr     ! bin number of hydrometeor 1
          DO j=1,nkr   ! bin number of hydrometeor 2
            cwlg_all(i,j,p_zi) = coll_kernel_interp(p_z,p2,p3,p4,ywlg_750mb(i,j),ywlg_500mb(i,j),ywlg_300mb(i,j))*dtimelnr
            cwlh_all(i,j,p_zi) = coll_kernel_interp(p_z,p2,p3,p4,ywlh_750mb(i,j),ywlh_500mb(i,j),ywlh_300mb(i,j))*dtimelnr
            cwls_all(i,j,p_zi) = coll_kernel_interp(p_z,p2,p3,p4,ywls_750mb(i,j),ywls_500mb(i,j),ywls_300mb(i,j))*dtimelnr
            cwsg_all(i,j,p_zi) = coll_kernel_interp(p_z,p2,p3,p4,ywsg_750mb(i,j),ywsg_500mb(i,j),ywsg_300mb(i,j))*dtimelnr
            cwss_all(i,j,p_zi) = coll_kernel_interp(p_z,p2,p3,p4,ywss_750mb(i,j),ywss_500mb(i,j),ywss_300mb(i,j))*dtimelnr
          END DO
        END DO
      END IF

      ! artificial increase of collision kernel due to turbulence:
      IF (coll_turb_fact == 1) THEN
        DO i=1,nkr
          DO j=1,nkr
            IF ((i < krdrop) .AND. (j < krdrop)) THEN
              cwturb=5
            ELSEIF ((i > krdrop) .AND. (j > krdrop)) THEN
              cwturb=2
            ELSE
              cwturb=3
            END IF
            cwll_all(i,j,p_zi)=cwturb*cwll_all(i,j,p_zi)
            cwlg_all(i,j,p_zi)=cwturb*cwlg_all(i,j,p_zi)
            cwlh_all(i,j,p_zi)=cwturb*cwlh_all(i,j,p_zi)
            cwls_all(i,j,p_zi)=cwturb*cwls_all(i,j,p_zi)
            cwsg_all(i,j,p_zi)=cwturb*cwsg_all(i,j,p_zi)
            cwss_all(i,j,p_zi)=cwturb*cwss_all(i,j,p_zi)
          END DO
        END DO
      END IF

      ! In future: import missing tables for cwgl, cwhl, cwsl, cwgs, cwhs, cwhs, cwgg, cwhh
      DO i=1,nkr   !larger particle
        DO j=1,nkr !smaller particle
          ! graupel - water
          cwgl_all(i,j,p_zi)=cwlg_all(j,i,p_zi) !ok
          ! hail - water
          cwhl_all(i,j,p_zi)=cwlh_all(j,i,p_zi) !ok
          ! snow - water
          cwsl_all(i,j,p_zi)=cwls_all(j,i,p_zi) !ok
          ! graupel - snow
          cwgs_all(i,j,p_zi)=cwsg_all(j,i,p_zi)
          ! snow - hail
          cwsh_all(i,j,p_zi)=cwsg_all(j,i,p_zi) !cwsh is missing
          ! hail - snow
          cwhs_all(i,j,p_zi)=cwsg_all(j,i,p_zi) !cwsh is missing
          ! grauel - graupel
          cwgg_all(i,j,p_zi)=cwsg_all(j,i,p_zi) !cwgg is missing
          ! hail - hail
          cwhh_all(i,j,p_zi)=cwsg_all(j,i,p_zi) !cwhh is missing
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE coll_kernel

  REAL(KIND=wp) FUNCTION coll_kernel_interp (p_z,p_1,p_2,p_3,ckern_1,ckern_2,ckern_3) ! Linear interpolation of Kernels vs height
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: p_1,p_2,p_3,ckern_1,ckern_2,ckern_3
    INTEGER,INTENT(IN) :: p_z
    IF (p_z>=p_1) coll_kernel_interp = ckern_1
    IF (p_z<=p_3) coll_kernel_interp = ckern_3
    IF (p_z<p_1 .AND. p_z>=p_2) coll_kernel_interp = ckern_2 + (ckern_1-ckern_2)*(p_z-p_2)/(p_1-p_2)
    IF (p_z<p_2 .AND. p_z>p_3) coll_kernel_interp = ckern_3 + (ckern_2-ckern_3)*(p_z-p_3)/(p_2-p_3)

    RETURN
  END FUNCTION coll_kernel_interp

  SUBROUTINE lognormal_modes_aerosol(fccnr_con,fccnr_mar,col,xl,rccn,itype)
    ! initial aerosol size distribution. Later vertical dependence is added, assuming the same distribution vs height
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: itype
    REAL(KIND=wp) ,INTENT(IN) :: xl(:), col
    REAL(KIND=wp) ,INTENT(OUT) :: fccnr_con(:), fccnr_mar(:), rccn(:)
    INTEGER :: kr
    REAL(KIND=wp)  :: ccncon1, ccncon2, ccncon3, radius_mean1, radius_mean2, radius_mean3,  &
                      sig1, sig2, sig3, fccnr_tmp(nkr), x0, r0, x0ccn, roccn(nkr), &
                      arg11,arg12,arg13,arg21,arg22,arg23,arg31,arg32,arg33, &
                      dnbydlogr_norm1,dnbydlogr_norm2,dnbydlogr_norm3
    REAL(KIND=wp) , PARAMETER :: rccn_max = 0.4e-4_wp         ! [cm]
    ! ... minimal radii for dry aerosol for the 3 log normal distribution
    REAL(KIND=wp) , PARAMETER :: rccn_min_3ln = 0.00048e-4_wp ! [cm]
    REAL(KIND=wp) , PARAMETER :: pi = 3.14159265_wp
    REAL(KIND=wp) , PARAMETER :: roccn0 = 1.0_wp

    ! note: rccn(1)  = 1.2 nm
    !       rccn(33) = 2.1 um
    x0ccn = xl(1)/(2.0**nkr)
    DO kr = nkr,1,-1
      roccn(kr) = roccn0 ! [g/cm3]
      x0 = x0ccn*2.0_wp**(kr) !previously was xccn(kr) = x0 (mass grid of ccn), but is not used here
      r0 = (3.0_wp*x0/4.0_wp/3.141593_wp/roccn(kr))**(1.0_wp/3.0_wp)
      rccn(kr) = r0
    END DO

    IF (itype == 1) THEN ! maritime regime

      ccncon1 = 340.000
      radius_mean1 = 0.00500e-04
      sig1 = 1.60000

      ccncon2 = 60.0000
      radius_mean2 = 0.03500e-04
      sig2 = 2.00000

      ccncon3 = 3.10000
      radius_mean3 = 0.31000e-04
      sig3 = 2.70000

    ELSE IF(itype == 2) THEN ! continental regime
      ccncon1 = 1000.000
      radius_mean1 = 0.00800e-04
      sig1 = 1.60000

      ccncon2 = 800.0000
      radius_mean2 = 0.03400e-04
      sig2 = 2.10000

      ccncon3 = 0.72000
      radius_mean3 = 0.46000e-04
      sig3 = 2.20000
    END IF

    fccnr_tmp = 0.0

    arg11 = ccncon1/(sqrt(2.0_wp*pi)*LOG(sig1))
    arg21 = ccncon2/(sqrt(2.0_wp*pi)*LOG(sig2))
    arg31 = ccncon3/(sqrt(2.0_wp*pi)*LOG(sig3))

    dnbydlogr_norm1 = 0.0
    dnbydlogr_norm2 = 0.0
    dnbydlogr_norm3 = 0.0
    DO kr = nkr,1,-1
      IF (rccn(kr) > rccn_min_3ln .AND. rccn(kr) < rccn_max)THEN
        arg12 = (LOG(rccn(kr)/radius_mean1))**2.0
        arg13 = 2.0_wp*((LOG(sig1))**2.0);
        dnbydlogr_norm1 = arg11*exp(-arg12/arg13)*(LOG(2.0)/3.0)
        arg22 = (LOG(rccn(kr)/radius_mean2))**2.0
        arg23 = 2.0_wp*((LOG(sig2))**2.0)
        dnbydlogr_norm2 = dnbydlogr_norm1 + arg21*exp(-arg22/arg23)*(LOG(2.0)/3.0)
        arg32 = (LOG(rccn(kr)/radius_mean3))**2.0
        arg33 = 2.0_wp*((LOG(sig3))**2.0)
        dnbydlogr_norm3 = dnbydlogr_norm2 + arg31*exp(-arg32/arg33)*(LOG(2.0)/3.0); !col=LOG(2.0)/3.0
        fccnr_tmp(kr) = dnbydlogr_norm3
      END IF
    END DO
    IF (itype == 1) fccnr_mar = fccnr_tmp/col !finally, fccnr_mar --> chem_new --> qna=col*sum(chem_new)
    IF (itype == 2) fccnr_con = fccnr_tmp/col !finally, fccnr_con --> chem_new --> qna=col*sum(chem_new)
    !Inform Koby: it was fccnr_con = scale_fa*fccnr_tmp here

    RETURN
  END SUBROUTINE lognormal_modes_aerosol

     ! +----------------------------------------------------+
  SUBROUTINE coll_bott_remap(xl, chucm, ima)
    ! rate of bin-shift fluxes during collisions in Bott, 2000 scheme
    ! kind of remapping in collisions scheme, which in coontrast to the
    ! remapping in the diff growth scheme which conserves mass and
    ! concentration, this one conserves only mass, but uses numerical scheme which
    ! has low numerical diffusion (kind of Courant nubmer idea)
    ! ima(i,j) - k-category number
    ! chucm(i,j) - courant number
    ! logarithmic grid distance(dlnr)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: xl(:)
    REAL(KIND=wp),INTENT(INOUT) :: chucm(:,:)
    INTEGER,INTENT(INOUT) :: ima(:,:)
    INTEGER :: k, kk, j, i
    REAL(KIND=wp) :: x0, xl_mg(nkr), dlnr

    xl_mg(1:nkr) = xl(1:nkr)*1.0e3_wp
    dlnr=LOG(2.0_wp)/(3.0_wp)
    DO i = 1,nkr
      DO j = i,nkr
        x0 = xl_mg(i) + xl_mg(j)
        IF (( i == nkr ) .OR. ( j == nkr )) THEN
          ima(i,j) = nkr
          chucm(i,j) = LOG(x0/xl_mg(nkr))/(3._wp*dlnr)
          goto 2000
        ELSE
          DO k = j,nkr
            kk = k
            IF (k == 1) goto 1000
            IF (xl_mg(k) >= x0 .AND. xl_mg(k-1) < x0) THEN
              chucm(i,j) = LOG(x0/xl_mg(k-1))/(3._wp*dlnr)
              IF (chucm(i,j) > 1.0_wp-1.e-8_wp) THEN
                chucm(i,j) = 0.0_wp
                kk = kk + 1
              END IF
              ima(i,j) = min(nkr-1,kk-1)
              goto 2000
            END IF
            1000 continue
          END DO
        END IF
        2000  continue
        chucm(j,i) = chucm(i,j)
        ima(j,i) = ima(i,j)
      END DO
    END DO

    RETURN
  END SUBROUTINE coll_bott_remap

  SUBROUTINE breakup_rain_prob(pkij,qkj,ecoalmassm,xl_r,jbreak,vr1)
    !...input variables
    !   gt    : mass distribution FUNCTION
    !   xt_mg : mass of bin in mg
    !...local variables
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: jbreak
    REAL(KIND=wp),INTENT(INOUT) :: ecoalmassm(:,:), pkij(:,:,:),qkj(:,:)
    REAL(KIND=wp),INTENT(IN) :: xl_r(:), vr1(:)
    CHARACTER(LEN=*), PARAMETER :: routine = TRIM(modname)//'::breakup_rain_prob'
    REAL(KIND=wp) :: vr1_d(nkr)
    INTEGER :: unitnr, error, ie,je,ke,i,j,k,ip,kp,jp,kq,jq
    CHARACTER*256 file_p, file_q

    ie = jbreak
    je = jbreak
    ke = jbreak

    !probability of break up after collision as function of drop mass
    IF (nkr == 43) file_p = 'SBM_input_43/'//'coeff_p43.dat'
    IF (nkr == 43) file_q = 'SBM_input_43/'//'coeff_q43.dat'
    IF (nkr == 33) file_p = 'SBM_input_33/'//'coeff_p_new_33.dat' ! new version 33 (taken from 43bins)
    IF (nkr == 33) file_q = 'SBM_input_33/'//'coeff_q_new_33.dat' ! new version 33 (taken from 43 bins)

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'breakup_rain_prob : opening coeff_p dat file'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(file_p), status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'breakup_rain_prob : coeff_p dat file not found'
        CALL message(modname,TRIM(txt))
        CALL finish(modname,'error in breakup_rain_prob')
      END IF
      DO k=1,ke
        DO i=1,ie
          DO j=1,i
            READ (unitnr,*) kp,ip,jp,pkij(kp,ip,jp) ! pkij=[g^3*cm^3/s]
          END DO
        END DO
      END DO
      CALL message(routine,"breakup_rain_prob: succesfull reading coeff_p dat file")
    END IF

    IF (my_process_is_stdio()) THEN
      unitnr=find_next_free_unit(10,999)
      WRITE(txt, '(a,i2)') 'breakup_rain_prob : opening coeff_q dat file'
      CALL message(modname,TRIM(txt))
      OPEN(unitnr, file=TRIM(file_q), status='old', form='formatted', iostat=error)
      IF (error /= 0) THEN
        WRITE (txt,*) 'breakup_rain_prob : coeff_q dat file not found'
        CALL message(modname,TRIM(txt))
        CALL finish(modname,'error in breakup_rain_prob')
      END IF
      DO k=1,ke
        DO j=1,je
          READ (unitnr,*) kq,jq,qkj(kq,jq)
        END DO
      END DO
      CALL message(routine,"breakup_rain_prob: succesfull reading coeff_q dat file")
    END IF

    vr1_d = vr1
    DO j=1,nkr
      DO i=1,nkr
        ecoalmassm(i,j)=ecoalmass(xl_r(i), xl_r(j), vr1_d)
      END DO
    END DO
    ! ... correction of coalescence efficiencies for drop collision kernels
    ! write Koby!? DO j=25,31
    ! write Koby!?   ecoalmassm(nkr,j)=0.1e-29_wp
    ! write Koby!? END DO

    RETURN
  END SUBROUTINE breakup_rain_prob

  REAL(wp) FUNCTION ecoalmass(x1, x2, vr1_breakup) !coalescence efficiency as FUNCTION of masses
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: vr1_breakup(nkr), x1, x2
    REAL(KIND=wp) :: rho, pi, akpi, deta, dksi

    rho=1.0_wp             ! [rho]=g/cm^3
    pi=3.1415927_wp
    akpi=6.0_wp/pi

    deta = (akpi*x1/rho)**(1.0_wp/3.0_wp)
    dksi = (akpi*x2/rho)**(1.0_wp/3.0_wp)

    ecoalmass = ecoaldiam(deta, dksi, vr1_breakup)

    RETURN
  END FUNCTION ecoalmass

  REAL(wp) FUNCTION ecoaldiam(deta,dksi,vr1_breakup) !coalescence efficiency as function of diameters
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: vr1_breakup(nkr),deta,dksi
    REAL(KIND=wp) :: dgr, dkl, q, qmin, qmax, e, x, e1, e2, sin1, cos1
    REAL(KIND=wp), PARAMETER :: eps=1.0e-30_wp,pi=3.1415927_wp

    dgr=max(deta,dksi)
    dkl=min(deta,dksi)
    q=0.5_wp*(0.5_wp*dkl+0.5_wp*dgr)
    qmin=250.0e-4_wp
    qmax=500.0e-4_wp

    IF (dkl<100.0e-4_wp) THEN
      e=1.0_wp
    ELSE IF (q<qmin) THEN
      e = ecoalochs(dgr,dkl, vr1_breakup)
    ELSE IF(q>=qmin.AND.q<qmax) THEN
      x=(q-qmin)/(qmax-qmin)
      sin1=sin(pi/2.0_wp*x)
      cos1=cos(pi/2.0_wp*x)
      e1=ecoalochs(dgr, dkl, vr1_breakup)
      e2=ecoallowlist(dgr, dkl, vr1_breakup)
      e=cos1**2*e1+sin1**2*e2
    ELSE IF(q>=qmax) THEN
      e=ecoallowlist(dgr, dkl, vr1_breakup)
    ELSE
      e=0.999_wp
    END IF
    ecoaldiam=max(min(1.0_wp,e),eps)

    RETURN
  END FUNCTION ecoaldiam

  REAL(wp) FUNCTION ecoallowlist(dgr,dkl,vr1_breakup) !coalescence efficiency for large drops (Low and List, 1982):
                                                              !experimental data + parametrization vs collision energy
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: vr1_breakup(nkr)
    REAL(KIND=wp),INTENT(INOUT) :: dgr, dkl
    REAL(KIND=wp) :: sigma, aka, akb, dstsc, st, sc, et, cke, qq0, qq1, qq2, ecl, w1, w2, dc
    REAL(KIND=wp), PARAMETER :: epsi=1.e-20_wp

    ! 1 j = 10^7 g cm^2/s^2
    sigma=72.8_wp  ! surface tension,[sigma]=g/s^2 (7.28e-2 n/m)
    aka=0.778_wp   ! empirical constant
    akb=2.61e-4_wp   ! empirical constant,[b]=2.61e6 m^2/j^2

    CALL collenergy(dgr,dkl,cke,st,sc,w1,w2,dc,vr1_breakup)

    dstsc=st-sc         ! diff. of surf. energies   [dstsc] = g*cm^2/s^2
    et=cke+dstsc        ! coal. energy,             [et]    =     "

    IF (et<50.0_wp) THEN    ! et < 5 uj (= 50 g*cm^2/s^2)
      qq0=1.0_wp+(dkl/dgr)
      qq1=aka/qq0**2
      qq2=akb*sigma*(et**2)/(sc+epsi)
      ecl=qq1*exp(-qq2)
    ELSE
      ecl=0.0_wp
    END IF
    ecoallowlist=ecl

    RETURN
  END FUNCTION ecoallowlist

  REAL(wp) FUNCTION ecoalochs(d_l,d_s,vr1_breakup) ! coalescence efficiency for small drops
                                                           ! (Beard and Ochs, 1984): experimental data
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: vr1_breakup(nkr), d_l, d_s
    REAL(KIND=wp) :: pi, sigma, r_s, r_l, p, vtl, vts, dv, weber_number, pa1, pa2, pa3, g
    REAL(KIND=wp), PARAMETER :: fpmin=1.e-30_wp

    pi=3.1415927_wp
    sigma=72.8_wp       ! surface tension [sigma] = g/s^2 (7.28e-2 n/m)
                       ! all in cgs (1 j = 10^7 g cm^2/s^2)
    r_s=0.5_wp*d_s
    r_l=0.5_wp*d_l
    p=r_s/r_l
    vtl=fall_veloc_drop_beard(d_l,vr1_breakup)
    vts=fall_veloc_drop_beard(d_s,vr1_breakup)
    dv=abs(vtl-vts)
    IF (dv<fpmin) dv=fpmin
    weber_number=r_s*dv**2/sigma
    pa1=1.0_wp+p
    pa2=1.0_wp+p**2
    pa3=1.0_wp+p**3
    g=2**(3.0_wp/2.0_wp)/(6.0_wp*pi)*p**4*pa1/(pa2*pa3)
    ecoalochs=0.767_wp-10.14_wp*weber_number**(0.5_wp)*g

    RETURN
  END FUNCTION ecoalochs

  SUBROUTINE collenergy(dgr,dkl,cke,st,sc,w1,w2,dc,vr1_breakup) !calculating the collision energy
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: vr1_breakup(nkr)
    REAL(KIND=wp),INTENT(INOUT) :: dgr, dkl, cke, st, sc, w1, w2, dc
    REAL(KIND=wp) :: pi, rho, sigma, ak10, dgka2, dgka3, v1, v2, dv, dgkb3
    REAL(KIND=wp), PARAMETER :: epsf = 1.e-30, fpmin = 1.e-30

    pi=3.1415927_wp
    rho=1.0_wp            ! water density,[rho]=g/cm^3
    sigma=72.8_wp         ! surf. tension,(h2o,20C)=7.28e-2 n/m
                         ! [sigma]=g/s^2
    ak10=rho*pi/12.0_wp
    dgr=max(dgr,epsf)
    dkl=max(dkl,epsf)
    dgka2=(dgr**2)+(dkl**2)
    dgka3=(dgr**3)+(dkl**3)
    IF (dgr/=dkl) THEN
      v1=fall_veloc_drop_beard(dgr,vr1_breakup)
      v2=fall_veloc_drop_beard(dkl,vr1_breakup)
      dv=(v1-v2)
      IF (dv<fpmin) dv=fpmin
      dv=dv**2
      IF (dv<fpmin) dv=fpmin
      dgkb3=(dgr**3)*(dkl**3)
      cke=ak10*dv*dgkb3/dgka3         ! collision energy [cke]=g*cm^2/s^2
    ELSE
      cke = 0.0_wp
    END IF
    st=pi*sigma*dgka2                 ! surf.energy (parent drop)
    sc=pi*sigma*dgka3**(2.0_wp/3.0_wp)  ! surf.energy (coal.system)
    w1=cke/(sc+epsf)                  ! weber number 1
    w2=cke/(st+epsf)                  ! weber number 2
    dc=dgka3**(1.0_wp/3.0_wp)           ! diam. of coal. system

    RETURN
  END SUBROUTINE collenergy

  REAL(wp) FUNCTION fall_veloc_drop_beard(diam,vr1_breakup) !calculating terminal velocity (Beard-formula)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN) :: vr1_breakup(nkr), diam
    INTEGER :: kr
    REAL(KIND=wp) :: aa

    aa   = diam/2.0_wp           ! radius in cm
    IF (aa <= dropradii(1)) fall_veloc_drop_beard=vr1_breakup(1)
    IF (aa > dropradii(nkr)) fall_veloc_drop_beard=vr1_breakup(nkr)
    DO kr=1,nkr-1
      IF (aa>dropradii(kr).AND.aa<=dropradii(kr+1)) THEN
        fall_veloc_drop_beard=vr1_breakup(kr+1)
      END IF
    END DO

    RETURN
  END FUNCTION fall_veloc_drop_beard

  SUBROUTINE ccn_init_sbm(chem_new,fccnr_con,fccnr_mar,fccnr_obs,xland,rhocgs,zcgs)
    IMPLICIT NONE
    REAL(KIND=wp), DIMENSION(:), INTENT(INOUT) :: chem_new
    REAL(KIND=wp), DIMENSION(:), INTENT(IN) :: fccnr_con, fccnr_mar, fccnr_obs
    REAL(KIND=wp), INTENT(IN) :: xland, rhocgs, zcgs
    REAL(KIND=wp) :: factz
    INTEGER :: kr,krr

    IF (zcgs .LE. zmin)THEN
      factz = 1.0
    ELSE
      factz=exp(-(zcgs-zmin)/z0in) !check that the height dependence is good like in mo_nwp_phy_init.f90 !
    END IF

    IF (ilognormal_modes_aerosol == 1)THEN
      ! ... generic ccn
      krr = 0
      DO kr = 1,nkr
        krr = krr + 1
        IF (xland > 0.9_wp) THEN !seems that my wk82 case is over sea
        !IF (xland < 0.9_wp) THEN ! pt to get continental ccn for my wk82 experiment
          chem_new(kr)=fccnr_con(krr)*factz*tune_sbmccn ![#/cm^3]
        ELSE
          chem_new(kr)=fccnr_mar(krr)*factz*tune_sbmccn ![#/cm^3]
        END IF
        !continental anyway:
        chem_new(kr)=fccnr_con(krr)*factz*tune_sbmccn ![#/cm^3]
        !maritime anyway:
        !chem_new(kr)=fccnr_mar(krr)*factz*tune_sbmccn ![#/cm^3]
      END DO
    ELSE
      ! ... ccn input from observation
      krr = 0
      DO kr = 1,nkr
        krr = krr + 1
        chem_new(kr) = fccnr_obs(krr)*factz ![#/cm^3]
      END DO
    END IF
    DO kr = 1,nkr
      chem_new(kr)=chem_new(kr)/(rhocgs/1000.0) ! local chem_new for ccn (used in sbm):
                                                ! [#/cm^3]. chem_new for ccn advected by the model is [#/kg]
    END DO
  END SUBROUTINE ccn_init_sbm

  SUBROUTINE sbm_init(p_patch, p_prog_now, fr_land, z_mc, dt_fast)
    TYPE(t_patch),    INTENT(IN)    :: p_patch
    TYPE(t_nh_prog),  INTENT(INOUT) :: p_prog_now         ! the prognostic variables
    REAL(wp),         INTENT(IN)    :: fr_land(:,:)       ! land fraction
    REAL(wp),         INTENT(IN)    :: z_mc(:,:,:)        ! height of model levels asl [m]
    REAL(wp),         INTENT(IN)    :: dt_fast            ! model time step

    INTEGER :: nlev
    INTEGER :: rl_start, rl_end
    INTEGER :: i_startblk, i_endblk    !> blocks
    INTEGER :: i_startidx, i_endidx    !! slices
    INTEGER :: jg, jb, jc, jk1, jk, krr
    INTEGER :: nshift                  !< shift with respect to global grid
    REAL(KIND=wp) :: fccnr_con(nkr), fccnr_mar(nkr), fccnr_obs(nkr)
    REAL(KIND=wp) :: zcgs, dtime_coal

    dtime_coal=MAX(dt_fast/REAL(ncondcoll),1.0_wp)
    rl_start = 1 ! Initialization should be done for all points
    rl_end   = min_rlcell
    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)
    jg = p_patch%id
    nlev   = p_patch%nlev !number of vertical levels
    nshift = p_patch%nshift_total

    IF (.NOT. atm_phy_nwp_config(jg)%lsbm_coupled) THEN !FALSE: use 2M for feedback and run uncoupled SBM, TRUE: use SBM feedback
      IF (tune_sbmccn < 1.0_wp) THEN
        atm_phy_nwp_config(jg)%cfg_2mom%ccn_type=6 ! maritime 2mom conditions
      ELSE IF (tune_sbmccn == 1.0_wp) THEN
        atm_phy_nwp_config(jg)%cfg_2mom%ccn_type=8 ! continental 2mom conditions
      END IF
      IF (jg == 1) CALL two_moment_mcrph_init(igscp=4, msg_level=msg_level, cfg_2mom=atm_phy_nwp_config(jg)%cfg_2mom)
    ELSE
      cfg_params = atm_phy_nwp_config(jg)%cfg_2mom
    END IF

    CALL sbm_data(fccnr_con,fccnr_mar,fccnr_obs)
    CALL coll_kernel(dtime_coal)

    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk=1,nlev
        jk1 = jk + nshift !nshift is 0 in WK case
        DO jc=i_startidx,i_endidx
          !in 2M height is defined similarly:  zf = 0.5_wp*(atmo%zh(i,k)+atmo%zh(i,k+1)), where zh=hhl. In SBM we directly use full levels
          !there is an option to subtract hag = p_metrics%z_mc(jc,jk,jb)-ext_data%atm%topography_c(jc,jb) but I am not sure we need it for
          !the ccn profile

          zcgs=z_mc(jc,jk1,jb)*100.0_wp

          CALL ccn_init_sbm( &
                         & chem_new=p_prog_now%tracer(jc,jk,jb,7+3*nkr+1:7+4*nkr), &
                         & fccnr_con=fccnr_con, &
                         & fccnr_mar=fccnr_mar, &
                         & fccnr_obs=fccnr_obs, &
                         & xland=fr_land(jc,jb),&
                         & rhocgs=0.001_wp*p_prog_now%rho(jc,jk1,jb),   &
                         & zcgs=zcgs)
        END DO
      END DO
    END DO

    ! Saving initial vertical ccn profile. The exact horizontal dependence is not important, so specific i_startidx, i_startblk are used
    IF ( use_ccn_const == 1 ) THEN
      IF (.NOT. ALLOCATED(ccnconstarr)) ALLOCATE(ccnconstarr(nlev,nkr))
      jb = i_startblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO krr=1,nkr
        DO jk=1,nlev
          ccnconstarr(jk,krr)=p_prog_now%tracer(i_startidx,jk,i_startblk,7+3*nkr+krr)
        END DO
      END DO
    END IF

    ! Special treatment of initial conditions for SBM microphysics < -------------------------------
    ! Moved here from atm_dyn_iconam/mo_initicon_utils.f90
    ! transfer all qx to qv in case of SBM.
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk=1,nlev
        DO jc=i_startidx,i_endidx
          p_prog_now%tracer(jc,jk,jb,iqv)=p_prog_now%tracer(jc,jk,jb,iqv) + &
                                          p_prog_now%tracer(jc,jk,jb,iqc) + &
                                          p_prog_now%tracer(jc,jk,jb,iqi) + &
                                          p_prog_now%tracer(jc,jk,jb,iqr) + &
                                          p_prog_now%tracer(jc,jk,jb,iqs)
          p_prog_now%tracer(jc,jk,jb,iqc)=0.0_wp
          p_prog_now%tracer(jc,jk,jb,iqi)=0.0_wp
          p_prog_now%tracer(jc,jk,jb,iqr)=0.0_wp
          p_prog_now%tracer(jc,jk,jb,iqs)=0.0_wp
        ENDDO
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE sbm_init

  SUBROUTINE breakup_rain_spont_prob(xl, prob, gain_var_new, nnd, ikr_spon_break)
    IMPLICIT NONE
    REAL(KIND=wp),INTENT(IN)  :: xl(:)
    REAL(KIND=wp),INTENT(OUT) :: prob(:), gain_var_new(:,:), nnd(:,:)
    REAL(KIND=wp) :: diameter(nkr), ratio_new, q_m, gain_var(nkr,nkr), xl_dp(nkr)
    INTEGER :: kr,i,j
    INTEGER, INTENT(INOUT) :: ikr_spon_break
    REAL(KIND=wp),PARAMETER :: gamma = 0.453_wp

    xl_dp = xl
    diameter(:) = dropradii(:)*2.0_wp*10.0_wp ! diameter in mm
    DO kr=1,nkr
      ikr_spon_break=kr
      IF (dropradii(kr)>=0.3) exit
    END DO

    WRITE (txt,*) 'ikr_spon_break=',ikr_spon_break
    CALL message(modname,TRIM(txt))

    IF (i_break_method==1) THEN
      DO kr=1,nkr
        prob(kr)=2.94e-7*exp(34.0_wp*dropradii(kr))
      END DO
    ELSE IF  (i_break_method==2) THEN
      DO kr=1,nkr
        prob(kr)=0.155e-3*exp(1.466_wp*10.0_wp*dropradii(kr))
      END DO
    END IF

    DO j=ikr_spon_break,nkr
      DO i=1,j-1
        !2 methods of breakup calculation (Sriavstava 1971 and Kamra et al, 1991):
        gain_var(j,i)=(145.37_wp/xl_dp(i))*(dropradii(i)/dropradii(j))*exp(-7.0_wp*dropradii(i)/dropradii(j))
        nnd(j,i)=gamma*exp(-gamma*diameter(i))/(1-exp(-gamma*diameter(j)))
      END DO
    END DO

    ! calculation the ratio that leads to mass conservation
    q_m = 0.0
    DO i=1,ikr_spon_break-1
      q_m = q_m + gain_var(ikr_spon_break,i)*xl_dp(i)**2;
    END DO
    ratio_new = q_m/xl_dp(ikr_spon_break)

    DO i=1,j-1
      DO j=ikr_spon_break,nkr
        gain_var_new(j,i) = gain_var(j,i)/ratio_new
      END DO
    END DO

    RETURN
  END SUBROUTINE breakup_rain_spont_prob

  SUBROUTINE qx_from_bins_diag(tracer, i_startidx, i_endidx, kstart, nlev)
    IMPLICIT NONE
    REAL(wp), DIMENSION(:,:,:), INTENT(INOUT), TARGET :: tracer
    INTEGER, INTENT(IN) :: i_startidx, i_endidx, kstart, nlev
    REAL(wp) :: qx
    INTEGER  :: i,k,iqb

    DO k = kstart, nlev
      DO i = i_startidx, i_endidx
        qx = 0.0_wp
        DO iqb = iqb_water_start,iqb_water_start-1+krdrop
            qx = qx + tracer(i,k,iqbin(iqb))
        END DO
        tracer(i,k,iqc)=qx

        qx = 0.0_wp
        DO iqb = krdrop+1, iqb_water_end
            qx = qx + tracer(i,k,iqbin(iqb))
        END DO
        tracer(i,k,iqr)=qx

        qx = 0.0_wp
        DO iqb = iqb_snow_start, iqb_snow_start-1+krice
            qx = qx + tracer(i,k,iqbin(iqb))
        END DO
        tracer(i,k,iqi)=qx

        qx = 0.0_wp
        DO iqb = iqb_snow_start+krice, iqb_snow_end
            qx = qx + tracer(i,k,iqbin(iqb))
        END DO
        tracer(i,k,iqs)=qx

        qx = 0.0_wp
        DO iqb = iqb_graupel_start, iqb_graupel_end
            qx = qx + tracer(i,k,iqbin(iqb))
        END DO
        tracer(i,k,iqg)=qx
      END DO
    END DO

    RETURN
  END SUBROUTINE

END MODULE mo_sbm_util

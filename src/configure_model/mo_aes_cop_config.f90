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

! Configuration of the parameterization for cloud optical properties,
! that is used in the AES physics package.

MODULE mo_aes_cop_config

  USE mo_exception            ,ONLY: message, print_value, finish
  USE mo_kind                 ,ONLY: wp
  USE mo_impl_constants       ,ONLY: max_dom
  USE mo_physical_constants   ,ONLY: tmelt
  USE mtime,                   ONLY: OPERATOR(>)
  USE mo_vertical_coord_table ,ONLY: vct_a

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::                     name   !< name for this unit

  ! configuration
  PUBLIC ::         aes_cop_config   !< user specified configuration parameters
  PUBLIC ::    init_aes_cop_config   !< allocate and initialize aes_cop_config
  PUBLIC ::    eval_aes_cop_config   !< evaluate aes_cop_config
  PUBLIC ::   print_aes_cop_config   !< print out

  !>
  !! Name of this unit
  !!
  CHARACTER(LEN=*), PARAMETER :: name = 'aes_cop'

  !>
  !! Configuration type containing parameters for the configuration of the cloud optical properties
  !!  and parameters of aes cloud microphysics, still used in mo_aes_convect_tables
  !!
  TYPE t_aes_cop_config
     !
     ! configuration parameters
     ! ------------------------
     !
     ! cloud droplet number concentration
     REAL(wp) :: cn1lnd   ! [1e6/m3] over land, p <= 100 hPa
     REAL(wp) :: cn2lnd   ! [1e6/m3] over land, p >= 800 hPa
     REAL(wp) :: cn1sea   ! [1e6/m3] over sea , p <= 100 hPa
     REAL(wp) :: cn2sea   ! [1e6/m3] over sea , p >= 800 hPa
     !
     ! cloud inhomogeneity factors
     REAL(wp) :: cinhomi            ! ice clouds
     REAL(wp) :: cinhoml_cfl        ! liquid water cumuliform clouds over land
     REAL(wp) :: cinhoml_cfo        ! liquid water cumuliform clouds over ocean
     REAL(wp) :: cinhoml_sf         ! liquid water stratiform clouds
     REAL(wp) :: cinhoml_lts_height ! height (m) used for computing the lower tropospheric stability (lts)
     REAL(wp) :: cinhoml_del1       !       del1 = ordinate value of atan2 function for transition from cfl/cfo to cs value
     REAL(wp) :: cinhoml_del2       ! lts - del2 = abcissae value of ..., lts = del2 -> (cfl/clo + cs)/2
     INTEGER  :: cinhoml_jk         ! level index for blending function of cinhoml cf and sf
     !
     ! freezing/deposition/sublimation for mo_aes_convect_tables:
     REAL(wp) :: cthomi   ! [K]      maximum temperature for homogeneous freezing
     REAL(wp) :: csecfrl  ! [kg/kg]  minimum in-cloud water mass mixing ratio in mixed phase clouds
     !
  END TYPE t_aes_cop_config

  !>
  !! Configuration state vectors, for multiple domains/grids.
  !!
  TYPE(t_aes_cop_config), TARGET :: aes_cop_config(max_dom)

CONTAINS

  !----

  !>
  !! Initialize the configuration state vector
  !!
  SUBROUTINE init_aes_cop_config

    !
    ! aes cloud microphysics active
    !
    ! cloud optical properties configuration
    ! --------------------------------------
    !
    ! cloud droplet number concentration
    aes_cop_config(:)% cn1lnd   =  20._wp
    aes_cop_config(:)% cn2lnd   = 180._wp
    aes_cop_config(:)% cn1sea   =  20._wp
    aes_cop_config(:)% cn2sea   =  80._wp
    !
    ! cloud inhomogeneity factors
    aes_cop_config(:)% cinhomi     = 0.80_wp
    aes_cop_config(:)% cinhoml_cfl = 0.40_wp
    aes_cop_config(:)% cinhoml_cfo = 0.40_wp
    aes_cop_config(:)% cinhoml_sf  = 0.80_wp
    aes_cop_config(:)% cinhoml_jk  = 0
    aes_cop_config(:)% cinhoml_lts_height = 3200._wp ! m
    aes_cop_config(:)% cinhoml_del1 =  2._wp
    aes_cop_config(:)% cinhoml_del2 = 20._wp
    !
    ! freezing/deposition/sublimation
    aes_cop_config(:)% cthomi   = tmelt-35.0_wp
    aes_cop_config(:)% csecfrl  = 1.5e-5_wp
    !
  END SUBROUTINE init_aes_cop_config

  !----

  !>
  !! Evaluate additional derived parameters
  !!
  SUBROUTINE eval_aes_cop_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg, jk
    CHARACTER(LEN=2)    :: cg
    !
    DO jg = 1,ng
       !
       IF (aes_cop_config(jg)% cinhoml_jk /= 0) THEN
          WRITE(cg,'(i0)') jg
          CALL finish('eval_aes_cop_config', &
               &      'aes_cop_config('//TRIM(cg)//')% cinhoml_jk must not be set in the namelist!')
       END IF
       !
       IF (aes_cop_config(jg)% cinhoml_sf /= aes_cop_config(jg)% cinhoml_cfo) THEN
          !
          aes_cop_config(jg)% cinhoml_jk = SIZE(vct_a)-1
          DO jk = 1,SIZE(vct_a)-1
             IF ((vct_a(jk)+vct_a(jk+1))*0.5_wp <= aes_cop_config(jg)% cinhoml_lts_height) THEN
                aes_cop_config(jg)% cinhoml_jk = jk
                EXIT
             END IF
          END DO
          !
       END IF
       !
    END DO
    !
  END SUBROUTINE eval_aes_cop_config

  !----

  !>
  !! Print out the user controlled configuration state
  !!
  SUBROUTINE print_aes_cop_config(ng)
    !
    INTEGER, INTENT(in) :: ng
    !
    INTEGER             :: jg
    CHARACTER(LEN=2)    :: cg
    !
    CALL message    ('','')
    CALL message    ('','========================================================================')
    CALL message    ('','')
    CALL message    ('','cloud optical properties configuration')
    CALL message    ('','======================================')
    CALL message    ('','')
    !
    DO jg = 1,ng
       !
       WRITE(cg,'(i0)') jg
       !
       CALL message    ('','For domain '//cg)
       CALL message    ('','------------')
       CALL message    ('','')
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cn1lnd   ',aes_cop_config(jg)% cn1lnd  )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cn2lnd   ',aes_cop_config(jg)% cn2lnd  )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cn1sea   ',aes_cop_config(jg)% cn1sea  )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cn2sea   ',aes_cop_config(jg)% cn2sea  )
       CALL message    ('','')
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhomi    ',aes_cop_config(jg)% cinhomi )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoml_cfl',aes_cop_config(jg)% cinhoml_cfl )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoml_cfo',aes_cop_config(jg)% cinhoml_cfo )
       CALL message    ('','')
       IF (aes_cop_config(jg)% cinhoml_sf /= aes_cop_config(jg)% cinhoml_cfo) THEN
          CALL message    ('','with modification to differentiate stratiform and cumuliform clouds:')
          CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoml_sf        ',aes_cop_config(jg)% cinhoml_sf )
          CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoml_del1      ',aes_cop_config(jg)% cinhoml_del1 )
          CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoml_del2      ',aes_cop_config(jg)% cinhoml_del2 )
          CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoml_lts_height',aes_cop_config(jg)% cinhoml_lts_height )
          CALL print_value('    aes_cop_config('//TRIM(cg)//')% cinhoml_jk        ',aes_cop_config(jg)% cinhoml_jk )
          CALL message    ('','')
       END IF
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% cthomi   ',aes_cop_config(jg)% cthomi  )
       CALL print_value('    aes_cop_config('//TRIM(cg)//')% csecfrl  ',aes_cop_config(jg)% csecfrl )
       CALL message    ('','')

       !
    END DO
    !
  END SUBROUTINE print_aes_cop_config

  !----

END MODULE mo_aes_cop_config

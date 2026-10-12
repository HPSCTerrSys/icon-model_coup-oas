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

! ---------------------------------------------------------------!
! Memory for splumes: Geometric factors depending on the grid
! only can be pre-calculated.
! ---------------------------------------------------------------

MODULE mo_bc_aeropt_splumes_memory

  USE mo_kind,                ONLY: dp, wp
  USE mo_model_domain,        ONLY: t_patch
  USE mo_impl_constants,      ONLY: SUCCESS
  USE mo_exception,           ONLY: message, finish
  USE mo_parallel_config,     ONLY: nproma
  USE mtime,                  ONLY: timedelta
  USE mo_time_config,         ONLY: time_config
  USE mo_var_list_register,   ONLY: vlr_add, vlr_del
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var, grib2_var
  USE mo_cdi,                 ONLY: DATATYPE_PACK16, DATATYPE_PACK24,  &
    &                               DATATYPE_FLT32,  DATATYPE_FLT64,   &
    &                               GRID_UNSTRUCTURED, cdiInqMissval
  USE mo_io_config,           ONLY: lnetcdf_flt64_output
  USE mo_gribout_config,      ONLY: gribout_config
  USE mo_cdi_constants,       ONLY: GRID_UNSTRUCTURED_CELL,    &
    &                               GRID_CELL
  USE mo_master_control,      ONLY: get_my_process_name
  USE mo_var_list,            ONLY: add_var, t_var_list_ptr

#include "add_var_acc_macro.inc"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: bc_spl_field, bc_spl_field_list       !< memory variables
  PUBLIC :: construct_bc_aeropt_splumes_memory    !< subroutine
  PUBLIC :: destruct_bc_aeropt_splumes_memory     !< subroutine
  PUBLIC :: t_bc_spl_field                        !< derived type
  PUBLIC :: nplumes, nfeatures, ntimes, nyears

  CHARACTER(len=*), PARAMETER :: thismodule = 'mo_bc_aeropt_splumes_memory'

  TYPE t_bc_spl_field
         ! fields depending on geometric features alone for simple plumes
    REAL(wp),POINTER ::       &
      & aer_sp_f1 (:,:,:)=>NULL(),    &!< contribution of feature 1
      & aer_sp_f2 (:,:,:)=>NULL(),    &!< contribution of feature 2
      & aer_sp_f3 (:,:,:)=>NULL(),    &!< contribution from feature 1 in natural background of Twomey effect
      & aer_sp_f4 (:,:,:)=>NULL()      !< contribution from feature 2 in natural background of Twomey effect
  END TYPE t_bc_spl_field

  TYPE(t_bc_spl_field), ALLOCATABLE, TARGET :: bc_spl_field(:) !< shape: (n_dom)
  TYPE(t_var_list_ptr), ALLOCATABLE  :: bc_spl_field_list(:)   !< shape: (n_dom)

  INTEGER, PARAMETER      ::     &
       nplumes   = 9            ,& !< Number of plumes
       nfeatures = 2            ,& !< Number of features per plume
       ntimes    = 52           ,& !< Number of times resolved per year (52 => weekly resolution)
       nyears    = 251             !< Number of years of available forcing

  REAL(dp), SAVE :: cdimissval

CONTAINS

  SUBROUTINE construct_bc_aeropt_splumes_memory(patch_array)
    TYPE(t_patch), INTENT(IN) :: patch_array(:)

    CHARACTER(len=16) :: listname_f
    INTEGER :: ndomain, ist, jg, nblks

    !---

    CALL message(thismodule,'Construction of memory stream for simple plume aerosols started.')

    cdimissval = cdiInqMissval()

    ! get actual number of domains
    ndomain = SIZE(patch_array)

    ALLOCATE(bc_spl_field(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(thismodule, &
       &'allocation of bc_spl_field for simple plumes failed')
    ALLOCATE(bc_spl_field_list(ndomain), STAT=ist)
    IF (ist/=SUCCESS) CALL finish(thismodule, &
       &'allocation of bc_spl_field_list for simple plumes failed')


    DO jg = 1,ndomain

      nblks = patch_array(jg)%nblks_c

      WRITE(listname_f,'(a,i2.2)') 'bc_spl_field_D',jg
      CALL new_bc_spl_field_list( jg,                   nproma,          nblks, listname_f, &
                                & bc_spl_field_list(jg),bc_spl_field(jg)                    )
    END DO

    CALL message(thismodule,'Construction of memory stream for simple plume aerosols finished.')

  END SUBROUTINE construct_bc_aeropt_splumes_memory

  SUBROUTINE destruct_bc_aeropt_splumes_memory

    INTEGER :: ndomain
    INTEGER :: jg
    INTEGER :: ist

    IF (ALLOCATED(bc_spl_field)) THEN
      CALL message(thismodule,'Destruction of memory stream for simple plume aerosols started.')
      ndomain = SIZE(bc_spl_field)

      DO jg = 1,ndomain
        CALL vlr_del(bc_spl_field_list(jg))
      END DO

      DEALLOCATE(bc_spl_field, STAT=ist)
      IF (ist/=SUCCESS) CALL finish(thismodule, &
         & 'deallocation of bc_spl_field array failed')

      CALL message(thismodule,'Destruction of memory stream for simple plume aerosols finished.')
    END IF

  END SUBROUTINE destruct_bc_aeropt_splumes_memory
  SUBROUTINE new_bc_spl_field_list( jg,                kproma,        kblks,     listname, &
                                  & bc_spl_field_list, bc_spl_field                        )

    INTEGER, INTENT(IN) :: jg !< domain index
    INTEGER, INTENT(IN) :: kproma, kblks
    CHARACTER(LEN=*), INTENT(IN) :: listname
    TYPE(t_var_list_ptr), INTENT(INOUT) :: bc_spl_field_list
    TYPE(t_bc_spl_field), INTENT(INOUT) :: bc_spl_field

    INTEGER           :: shape3d_aer_sp(3)
    REAL(wp)          :: initial_value
    TYPE(t_cf_var)    ::    cf_desc
    TYPE(t_grib2_var) :: grib2_desc
    INTEGER           :: datatype_flt, ibits, iextbits, ivarbits


    datatype_flt = MERGE(DATATYPE_FLT64, DATATYPE_FLT32, lnetcdf_flt64_output)
    ibits        = DATATYPE_PACK16
    iextbits     = DATATYPE_PACK24
    ivarbits     = MERGE(DATATYPE_PACK24, DATATYPE_PACK16, gribout_config(jg)%lgribout_24bit)

    CALL vlr_add(bc_spl_field_list, listname, patch_id=jg, lrestart=.FALSE., &
      &          model_type=get_my_process_name())

    ! weights for simple plume aerosols
    shape3d_aer_sp  = (/kproma, nplumes, kblks/)
    initial_value = -99999._wp
    cf_desc    = t_cf_var('aer_sp_f1','-', &
        & 'weight for simple plumes - contribution from feature 1', &
        & datatype_flt)
    grib2_desc = grib2_var(0,20,201, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( bc_spl_field_list, 'aer_sp_f1', bc_spl_field%aer_sp_f1,      &
                & GRID_UNSTRUCTURED_CELL, nplumes, cf_desc, grib2_desc,        &
                & ldims=shape3d_aer_sp,lrestart = .FALSE., lopenacc=.TRUE.,    &
                & initval=initial_value                                        )
    __acc_attach(bc_spl_field%aer_sp_f1)
    cf_desc    = t_cf_var('aer_sp_f2','-', &
        & 'weight for simple plumes - contribution from feature 2', &
        & datatype_flt)
    grib2_desc = grib2_var(0,20,202, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( bc_spl_field_list, 'aer_sp_f2', bc_spl_field%aer_sp_f2,      &
                & GRID_UNSTRUCTURED_CELL, nplumes, cf_desc, grib2_desc,        &
                & ldims=shape3d_aer_sp,lrestart = .FALSE., lopenacc=.TRUE.,    &
                & initval=initial_value                                        )
    __acc_attach(bc_spl_field%aer_sp_f2)
    cf_desc    = t_cf_var('aer_sp_f3','-', &
        & 'weight for simple plumes - contribution from feature 3', &
        & datatype_flt)
    grib2_desc = grib2_var(0,20,203, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( bc_spl_field_list, 'aer_sp_f3', bc_spl_field%aer_sp_f3,      &
                & GRID_UNSTRUCTURED_CELL, nplumes, cf_desc, grib2_desc,        &
                & ldims=shape3d_aer_sp,lrestart = .FALSE., lopenacc=.TRUE.,    &
                & initval=initial_value                                        )
    __acc_attach(bc_spl_field%aer_sp_f3)
    cf_desc    = t_cf_var('aer_sp_f4','-', &
        & 'weight for simple plumes - contribution from feature 4', &
        & datatype_flt)
    grib2_desc = grib2_var(0,20,204, ibits, GRID_UNSTRUCTURED, GRID_CELL)
    CALL add_var( bc_spl_field_list, 'aer_sp_f4', bc_spl_field%aer_sp_f4,      &
                & GRID_UNSTRUCTURED_CELL, nplumes, cf_desc, grib2_desc,        &
                & ldims=shape3d_aer_sp,lrestart = .FALSE., lopenacc=.TRUE.,    &
                & initval=initial_value                                        )
    __acc_attach(bc_spl_field%aer_sp_f4)

  END SUBROUTINE new_bc_spl_field_list

END MODULE mo_bc_aeropt_splumes_memory

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

MODULE mo_wave_grib2

  USE mo_exception,          ONLY: finish
  USE mo_cdi,                ONLY: cdiInqKeyInt, cdiDefKeyInt, CDI_NOERR, &
    &                              CDI_KEY_TYPEOFGENERATINGPROCESS,       &
    &                              CDI_KEY_PRODUCTDEFINITIONTEMPLATE
  USE mo_var_metadata_types, ONLY: t_var_metadata, CLASS_WAVE_SPECTRUM

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_grib2_pdt_wave_spectra

  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_wave_grib2'
CONTAINS

  !------------------------------------------------------------------------------------------------
  !> Set product definition template (PDT) for wave 2D spectra
  !
  !
  SUBROUTINE set_grib2_pdt_wave_spectra(vlistID, varID, info)
    INTEGER,               INTENT(IN) :: vlistID, varID
    TYPE(t_var_metadata),  INTENT(IN) :: info

    ! local
    INTEGER :: res
    INTEGER :: typeOfGeneratingProcess
    INTEGER :: productDefinitionTemplate

    CHARACTER(len=*), PARAMETER :: routine = modname//'::set_grib2_pdt_wave_spectra'
  !----------------------------------------------------------------

    ! Skip inapplicable fields
    IF ( info%var_class /= CLASS_WAVE_SPECTRUM ) RETURN

    res = cdiInqKeyInt(vlistID, varID, CDI_KEY_TYPEOFGENERATINGPROCESS, typeOfGeneratingProcess)
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when inquiring typeOfGeneratingProcess")
    ENDIF

    IF (typeOfGeneratingProcess == 4) THEN
      ! ensemble
      ! 'Individual ensemble forecast, control and perturbed, at a horizontal level
      ! or in a horizontal layer at a point in time for wave 2D spectra with
      ! frequencies and directions defined by formulae'
      productDefinitionTemplate = 102
    ELSE
      ! deterministic
      ! 'Analysis or forecast at a horizontal level or in a horizontal layer
      ! at a point in time for wave 2D spectra with frequencies and directions defined by formulae'
      productDefinitionTemplate = 101
    ENDIF

    res = cdiDefKeyInt(vlistID, varID, CDI_KEY_PRODUCTDEFINITIONTEMPLATE, productDefinitionTemplate)
    IF (res/=CDI_NOERR) THEN
      CALL finish(routine, "error when defining productDefinitionTemplate")
    ENDIF

  END SUBROUTINE set_grib2_pdt_wave_spectra

END MODULE mo_wave_grib2

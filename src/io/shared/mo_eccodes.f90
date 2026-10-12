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

!
! Wrappers for ecCodes interfaces.
!
! Notes:
!
! - Main objectives:
!   * Encapsulate error handling
!   * Bundle ecCodes GRIB keys in groups
!     (mainly, but not exclusively, in terms of the sections, in which GRIB metadata are subdivided)
!   * Confine "HAVE_ECCODES-ifdefing" to this one module
!
! - Certain sequences of subroutine calls repeat themselves again and again
!   inside the wrappers below. Nevertheless, they are not wrapped into auxiliary subroutines,
!   in order to adhere to a shallow nesting.
!   Subroutine calls (especially repeating ones) are unwelcome in ICON-NWP, as they cost some runtime.
!
! - The wrappers are low-level and hand over the values
!   of basically "pure" ecCodes GRIB and concept keys to the caller.
!   So this is not the place to apply higher-level logics or diagnostics (that is up to the callers).
!
! - Please try to avoid optional arguments.
!   99% of the wrapper arguments below are rather "lightweight" in terms of memory consumption,
!   so there is no compelling reason to take the trouble of optional arguments.
!
! - The procedures of this module are likely called by Work PEs, but not the Workroot PE.
!   As a consequence, we have to call "warning" instead of "message".
!   If a Work PE would call "message", nothing would be printed.
!

MODULE mo_eccodes

  USE mo_kind,           ONLY: i4, i8, sp, dp, wp
  USE mo_impl_constants, ONLY: SUCCESS
  USE mo_exception,      ONLY: finish, warning, message_text
#ifdef HAVE_ECCODES
  USE eccodes,           ONLY: kindOfInt, kindOfLong, kindOfFloat,                               &
    &                          CODES_SUCCESS, CODES_NULL, CODES_NULL_HANDLE, CODES_END_OF_FILE,  &
    &                          codes_get_error_string, codes_open_file, codes_close_file,        &
    &                          codes_release, codes_get_string, codes_get_int, codes_is_defined, &
    &                          codes_is_missing, codes_get_byte_array, codes_get_long,           &
    &                          codes_get_size_long, codes_set, codes_get_real4_array,            &
    &                          codes_read_from_file_int4, codes_new_from_message_int4,           &
    &                          codes_get_real4
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ECC_NULL_HANDLE
  PUBLIC :: ECC_MAX_SHORTNAME_LENGTH
  PUBLIC :: ECC_GRID_ELEMENT_CELL
  PUBLIC :: ECC_GRID_ELEMENT_VERTEX
  PUBLIC :: ECC_GRID_ELEMENT_EDGE
  PUBLIC :: ECC_NULL
  PUBLIC :: ecc_open_file
  PUBLIC :: ecc_close_file
  PUBLIC :: ecc_read_record
  PUBLIC :: ecc_new_from_record
  PUBLIC :: ecc_release
  PUBLIC :: ecc_get_info_on_product
  PUBLIC :: ecc_get_info_on_datetime
  PUBLIC :: ecc_get_info_on_horizontal_grid
  PUBLIC :: ecc_get_info_on_vertical_grid
  PUBLIC :: ecc_get_info_on_generating_centre
  PUBLIC :: ecc_get_info_on_generating_process
  PUBLIC :: ecc_get_info_on_tiles
  PUBLIC :: ecc_get_local_info
  PUBLIC :: ecc_get_values

  !------------
  ! Parameters
  !------------

  !> Module name
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_eccodes'

  !> Error message: ecCodes is not enabled
  CHARACTER(LEN=*), PARAMETER :: ECC_ERROR_MESSAGE_NO_API = "The Fortran interface of ecCodes is not enabled"

  !> Maximum length of value of ecCodes concept key: shortName
  INTEGER, PARAMETER :: ECC_MAX_SHORTNAME_LENGTH = 50

  !> Grid element identifiers:
  INTEGER, PARAMETER :: ECC_GRID_ELEMENT_CELL   = 1
  INTEGER, PARAMETER :: ECC_GRID_ELEMENT_VERTEX = 2
  INTEGER, PARAMETER :: ECC_GRID_ELEMENT_EDGE   = 3

#ifdef HAVE_ECCODES
  !> Eccodes: kind of standard integer
  INTEGER, PARAMETER :: ECC_kindOfInt = kindOfInt
  !> Eccodes: kind of long integer
  INTEGER, PARAMETER :: ECC_kindOfLong = kindOfLong
  !> Eccodes: kind of single precision
  INTEGER, PARAMETER :: ECC_kindOfFloat = kindOfFloat
  !> Eccodes: status success
  INTEGER, PARAMETER :: ECC_SUCCESS = CODES_SUCCESS
  !> Eccodes: status end of file
  INTEGER, PARAMETER :: ECC_END_OF_FILE = CODES_END_OF_FILE
  !> Eccodes: null
  INTEGER, PARAMETER :: ECC_NULL = INT(CODES_NULL)
  !> Eccodes: null handle
  INTEGER, PARAMETER :: ECC_NULL_HANDLE = INT(CODES_NULL_HANDLE)
#else
  INTEGER, PARAMETER :: ECC_kindOfInt   = i4
  INTEGER, PARAMETER :: ECC_kindOfLong  = i8
  INTEGER, PARAMETER :: ECC_kindOfFloat = sp
  INTEGER, PARAMETER :: ECC_SUCCESS     = SUCCESS
  INTEGER, PARAMETER :: ECC_END_OF_FILE = -1
  INTEGER, PARAMETER :: ECC_NULL        = -1
  INTEGER, PARAMETER :: ECC_NULL_HANDLE = -20
#endif

  !> GRIB "missing" for 1-byte-keys
  INTEGER, PARAMETER :: ECC_MISSING = 255

  !> Rescale factors
  REAL(wp), PARAMETER :: ecc_downscale_factor(0:10) = &
    & [1.0_wp, 1.0E-1_wp, 1.0E-2_wp, 1.0E-3_wp, 1.0E-4_wp, 1.0E-5_wp, 1.0E-6_wp, 1.0E-7_wp, 1.0E-8_wp, 1.0E-9_wp, 1.0E-10_wp]

  REAL(wp), PARAMETER :: ecc_upscale_factor(1:10) = &
    & [1.0E+1_wp, 1.0E+2_wp, 1.0E+3_wp, 1.0E+4_wp, 1.0E+5_wp, 1.0E+6_wp, 1.0E+7_wp, 1.0E+8_wp, 1.0E+9_wp, 1.0E+10_wp]

CONTAINS

  !>
  !! @brief Error handling
  !!
  SUBROUTINE ecc_error_handling(routine_of_occurrence, error)

    !-----------
    ! Arguments
    !-----------

    !> Routine where error occurred
    CHARACTER(LEN=*), INTENT(IN) :: routine_of_occurrence

    !> ecCodes error identifier
    INTEGER(KIND=ECC_kindOfInt), INTENT(IN) :: error

    !-----------------
    ! Local variables
    !-----------------

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_error_handling'

    !----------------------------

#ifdef HAVE_ECCODES

    ! Get error message from ecCodes that corresponds to error identifier
    CALL codes_get_error_string(error=error, error_message=message_text, status=ecc_status)
    IF (ecc_status /= ECC_SUCCESS) message_text = "Request of error string from ecCodes failed"

    CALL finish(modname//routine_of_occurrence, message_text)

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_error_handling

  !--------------------------------------------------------------

  !>
  !! @brief Open a GRIB file and return its handle
  !!
  SUBROUTINE ecc_open_file(grib_filename, ecc_mode, ecc_ifile)

    !-----------
    ! Arguments
    !-----------

    !> Name of GRIB file to be opened
    CHARACTER(LEN=*), INTENT(IN) :: grib_filename

    !> File mode: 'r' -> read, 'w' -> write
    CHARACTER(LEN=*), INTENT(IN) :: ecc_mode

    !> File identifier/handle from ecCodes
    INTEGER, INTENT(OUT) :: ecc_ifile

    !-----------------
    ! Local variables
    !-----------------

    !> Length of filename
    INTEGER :: grib_filename_length

    !> Local ecc_ifile
    INTEGER(KIND=ECC_kindOfInt) :: ecc_ifile_local

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Flag to indicate if file exists
    LOGICAL :: does_grib_file_exist

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_open_file'

    !----------------------------

    ! Initialize intent-out argument
    ecc_ifile = ECC_NULL_HANDLE

    ! Length of file name
    grib_filename_length = LEN_TRIM(grib_filename)

    ! Check arguments
    IF (grib_filename_length < 1) THEN
      CALL finish(modname//routine, "No file name provided")
    ELSEIF (.NOT. ((ecc_mode(1:1) == "r") .OR. (ecc_mode(1:1) == "w"))) THEN
      CALL finish(modname//routine, "Invalid ecCodes file mode")
    ENDIF

    ! Check if file exists
    INQUIRE(file=grib_filename(1:grib_filename_length), exist=does_grib_file_exist)

    IF (.NOT. does_grib_file_exist) THEN
      message_text = "The file: "//grib_filename(1:grib_filename_length)//" does not exist"
      CALL finish(modname//routine, message_text)
    END IF

#ifdef HAVE_ECCODES

    CALL codes_open_file(ifile=ecc_ifile_local, filename=grib_filename(1:grib_filename_length), &
      &                  mode=ecc_mode(1:1), status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_ifile = INT(ecc_ifile_local)

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_open_file

  !--------------------------------------------------------------

  !>
  !! @brief Close a GRIB file
  !!
  SUBROUTINE ecc_close_file(ecc_ifile)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of file to be closed
    INTEGER, INTENT(IN) :: ecc_ifile

    !-----------------
    ! Local variables
    !-----------------

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_close_file'

    !----------------------------

#ifdef HAVE_ECCODES

    CALL codes_close_file(ifile=INT(ecc_ifile, KIND=ECC_kindOfInt), status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_close_file

  !--------------------------------------------------------------

  !>
  !! @brief Read next message (record) from GRIB file
  !!
  SUBROUTINE ecc_read_record(ecc_ifile, ecc_record, ecc_nbytes, ecc_record_length, ecc_eof)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of file to be closed
    INTEGER, INTENT(IN) :: ecc_ifile

    !> Buffer to hold the GRIB message (GRIB record)
    INTEGER(KIND=i4), INTENT(INOUT) :: ecc_record(:)

    !> Length of GRIB message (GRIB record) in bytes
    INTEGER, INTENT(INOUT) :: ecc_nbytes

    !> Length of GRIB message (GRIB record) in multiples of 4 bytes
    INTEGER, INTENT(OUT) :: ecc_record_length

    !> Flag to indicate that end of GRIB file is reached
    LOGICAL, INTENT(OUT) :: ecc_eof

    !-----------------
    ! Local variables
    !-----------------

    !> Local length of GRIB message in bytes
    INTEGER(KIND=ECC_kindOfInt) :: ecc_nbytes_local

    !> Remainder for computation of GRIB message length
    INTEGER :: remainder

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_read_record'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_record_length = 0
    ecc_eof           = .FALSE.

#ifdef HAVE_ECCODES

    ! Initially assumed size of record in bytes
    ecc_nbytes_local = INT(ecc_nbytes, KIND=ECC_kindOfInt)

    CALL codes_read_from_file_int4(ifile=INT(ecc_ifile, KIND=ECC_kindOfInt), buffer=ecc_record, &
      &                            nbytes=ecc_nbytes_local, status=ecc_status)

    IF ((ecc_status /= ECC_SUCCESS) .AND. (ecc_status /= ECC_END_OF_FILE)) THEN
      CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)
    ELSEIF ((ecc_status /= ECC_SUCCESS) .AND. (ecc_status == ECC_END_OF_FILE)) THEN
      ecc_eof = .TRUE.
      RETURN
    ENDIF

    ! Return actual size of record in bytes
    ecc_nbytes = INT(ecc_nbytes_local)

    ! Compute "standard" length of GRIB message (multiples of 4 bytes)
    ecc_record_length = ecc_nbytes / 4

    ! The length of a GRIB message in bytes is not necessary a multiple of 4, therefore:
    remainder = MAX(0, ecc_nbytes - ecc_record_length * 4)

    IF (remainder > 0) ecc_record_length = ecc_record_length + 1

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_read_record

  !--------------------------------------------------------------

  !>
  !! @brief Get handle of GRIB message (GRIB record) from "raw" GRIB message
  !!
  SUBROUTINE ecc_new_from_record(ecc_record, ecc_msgid)

    !-----------
    ! Arguments
    !-----------

    !> The "raw" GRIB message (GRIB record)
    INTEGER(KIND=i4), INTENT(IN) :: ecc_record(:)

    !> GRIB message identifier/handle from ecCodes
    INTEGER, INTENT(OUT) :: ecc_msgid

    !-----------------
    ! Local variables
    !-----------------

    !> Local GRIB message identifier/handle
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_new_from_record'

    !----------------------------

    ! Initialize intent-out argument
    ecc_msgid = ECC_NULL_HANDLE

#ifdef HAVE_ECCODES

    CALL codes_new_from_message_int4(msgid=ecc_msgid_local, message=ecc_record, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ! Return the handle
    ecc_msgid = INT(ecc_msgid_local)

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_new_from_record

  !--------------------------------------------------------------

  !>
  !! @brief Release a GRIB message (GRIB record)
  !!
  SUBROUTINE ecc_release(ecc_msgid)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message to be released
    INTEGER, INTENT(IN) :: ecc_msgid

    !-----------------
    ! Local variables
    !-----------------

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_release'

    !----------------------------

#ifdef HAVE_ECCODES

    CALL codes_release(msgid=INT(ecc_msgid, KIND=ECC_kindOfInt), status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_release

  !--------------------------------------------------------------

  !>
  !! @brief Get values of product-related ecCodes GRIB and concept keys
  !!
  SUBROUTINE ecc_get_info_on_product(ecc_msgid, ecc_productDefinitionTemplateNumber,             &
       &                             ecc_discipline, ecc_parameterCategory, ecc_parameterNumber, &
       &                             ecc_shortName, ecc_shortName_length, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Value of GRIB key: productDefinitionTemplateNumber
    INTEGER, INTENT(OUT) :: ecc_productDefinitionTemplateNumber

    !> Value of product discipline
    INTEGER, INTENT(OUT) :: ecc_discipline

    !> Value of product category
    INTEGER, INTENT(OUT) :: ecc_parameterCategory

    !> Value of product number
    INTEGER, INTENT(OUT) :: ecc_parameterNumber

    !> Value of ecCodes concept key: shortName
    CHARACTER(LEN=*), INTENT(OUT) :: ecc_shortName

    !> Length of value of ecCodes concept key: shortName
    INTEGER, INTENT(OUT) :: ecc_shortName_length

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Flag to indicate if GRIB key is defined
    INTEGER(KIND=ECC_kindOfInt) :: ecc_is_defined

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_info_on_product'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_productDefinitionTemplateNumber = ECC_NULL
    ecc_discipline        = ECC_NULL
    ecc_parameterCategory = ECC_NULL
    ecc_parameterNumber   = ECC_NULL
    ecc_shortName         = " "
    ecc_shortName_length  = 0
    ecc_successful        = .FALSE.

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! Get identifier of the product definition template:
    ! (The key 'productDefinitionTemplateNumber' is an integral part of product definition section 4.
    ! Therefore, we do not check for its existence.)
    CALL codes_get_int(msgid=ecc_msgid_local, key='productDefinitionTemplateNumber', &
      &                value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_productDefinitionTemplateNumber = INT(ecc_integer_value)

    ! Parameter discipline:
    ! (The key 'discipline' is an integral part of indicator section 0.
    ! Therefore, we do not check for its existence.)
    CALL codes_get_int(msgid=ecc_msgid_local, key='discipline', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_discipline = INT(ecc_integer_value)

    ! Parameter category:
    ! (Parameter category and number are in principle template-specific keys.
    ! Therefore, we first check for their existence.)
    CALL codes_is_defined(msgid=ecc_msgid_local, key='parameterCategory', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ! If this GRIB key is not part of the metadata of the GRIB message (GRIB record),
    ! there is no need to continue and we can return
    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'parameterCategory' is not defined")
      RETURN
    ENDIF

    CALL codes_get_int(msgid=ecc_msgid_local, key='parameterCategory', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_parameterCategory = INT(ecc_integer_value)

    ! Parameter number:
    CALL codes_is_defined(msgid=ecc_msgid_local, key='parameterNumber', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'parameterNumber' is not defined")
      RETURN
    ENDIF

    CALL codes_get_int(msgid=ecc_msgid_local, key='parameterNumber', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_parameterNumber = INT(ecc_integer_value)

    ! Short name:
    ! (This ecCodes concept key as such should always be defined.)
    CALL codes_get_string(msgid=ecc_msgid_local, key='shortName', value=ecc_shortName, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_shortName_length = LEN_TRIM(ecc_shortName)

    ! A valid short name is mandatory,
    ! otherwise we have no idea which field we got.
    ! (If ecCodes cannot find a shortName value that matches the GRIB metadata,
    ! it should return "unknow".)
    ecc_successful = (ecc_shortName_length > 0)
    IF (ecc_successful) ecc_successful = ALL(["unknown", "UNKNOWN"] /= ecc_shortName(1:ecc_shortName_length))

    IF (.NOT. ecc_successful) CALL warning(modname//routine, "No valid 'shortName' found")

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_info_on_product

  !--------------------------------------------------------------

  !>
  !! @brief Get dataDate + dataTime (reference date and time) and validityDate + validityTime
  !!
  SUBROUTINE ecc_get_info_on_datetime(ecc_msgid, ecc_significanceOfReferenceTime, &
    &                                 ecc_dataDateTime, ecc_validityDateTime,     &
    &                                 ecc_dataDateTime_length, ecc_validityDateTime_length, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Meaning of reference date and time
    INTEGER, INTENT(OUT) :: ecc_significanceOfReferenceTime

    !> Reference date and time in format: 'YYYY-MM-DDThh:mm:ss'
    CHARACTER(LEN=*), INTENT(OUT) :: ecc_dataDateTime

    !> Validity date and time in format: 'YYYY-MM-DDThh:mm:ss'
    CHARACTER(LEN=*), INTENT(OUT) :: ecc_validityDateTime

    !> Length of reference date and time string
    INTEGER, INTENT(OUT) :: ecc_dataDateTime_length

    !> Length of validity date and time string
    INTEGER, INTENT(OUT) :: ecc_validityDateTime_length

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> For value of key 'dataDate' with format 'YYYYMMDD'
    CHARACTER(LEN=10) :: ecc_dataDate

    !> For value of key 'dataTime' with format 'hhmm'
    CHARACTER(LEN=6) :: ecc_dataTime

    !> For value of key 'validityDate' with format 'YYYYMMDD'
    CHARACTER(LEN=10) :: ecc_validityDate

    !> For value of key 'validityTime' with format 'hhmm'
    CHARACTER(LEN=6) :: ecc_validityTime

    !> Actual length of 'dataTime' and 'validityTime'
    INTEGER :: ecc_time_length

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_info_on_datetime'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_significanceOfReferenceTime = ECC_MISSING
    ecc_successful = .FALSE.

    ! LEN_TRIM('YYYY-MM-DDThh:mm:ss') -> 19
    IF (LEN(ecc_dataDateTime) < 19) THEN
      CALL finish(modname//routine, "Length of argument 'ecc_dataDateTime' is too short (LEN>=19 required)")
    ELSEIF (LEN(ecc_validityDateTime) < 19) THEN
      CALL finish(modname//routine, "Length of argument 'ecc_validityDateTime' is too short (LEN>=19 required)")
    ENDIF

    ecc_dataDateTime            = ' '
    ecc_validityDateTime        = ' '
    ecc_dataDateTime_length     = 0
    ecc_validityDateTime_length = 0

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! Get meaning of reference date and time:
    ! (The key 'significanceOfReferenceTime' is an integral part of identification section 1.
    ! Therefore, we do not check for its existence.)
    CALL codes_get_int(msgid=ecc_msgid_local, key='significanceOfReferenceTime', &
      &                value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_significanceOfReferenceTime = INT(ecc_integer_value)

    ! Get reference date and time:
    ! (These ecCodes concept keys as such should always be defined.)
    CALL codes_get_string(msgid=ecc_msgid_local, key='dataDate', value=ecc_dataDate, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    CALL codes_get_string(msgid=ecc_msgid_local, key='dataTime', value=ecc_dataTime, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ! 'dataTime' has format 'hhmm'. However, if its values is '0000', for instance,
    ! ecCodes may return just '0'. In such a case, we pad with '0's.
    ecc_time_length = LEN_TRIM(ecc_dataTime)
    ! ('ecc_time_length > 0' is just for safety reasons)
    IF ((ecc_time_length > 0) .AND. (ecc_time_length < 4)) THEN
      ecc_dataTime(ecc_time_length+1:4) = REPEAT('0', 4 - ecc_time_length)
    ENDIF

    ! Get validity date and time:
    ! (These ecCodes concept keys as such should always be defined.)
    CALL codes_get_string(msgid=ecc_msgid_local, key='validityDate', value=ecc_validityDate, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    CALL codes_get_string(msgid=ecc_msgid_local, key='validityTime', value=ecc_validityTime, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ! 'validityTime' has format 'hhmm'. However, if its values is '0000', for instance,
    ! ecCodes may return just '0'. In such a case, we pad with '0's.
    ecc_time_length = LEN_TRIM(ecc_validityTime)
    ! ('ecc_time_length > 0' is just for safety reasons)
    IF ((ecc_time_length > 0) .AND. (ecc_time_length < 4)) THEN
      ecc_validityTime(ecc_time_length+1:4) = REPEAT('0', 4 - ecc_time_length)
    ENDIF

    ! Transfer 'YYYYMMDD' and 'hhmm' into desired format: 'YYYY-MM-DDThh:mm:ss':

    ! Both reference and validity date are mandatory
    IF (ALL([LEN_TRIM(ecc_dataDate), LEN_TRIM(ecc_validityDate)] >= 8)) THEN

      ecc_dataDateTime = ecc_dataDate(1:4)//"-"//ecc_dataDate(5:6)//"-"//ecc_dataDate(7:8) &
        &              //"T"//ecc_dataTime(1:2)//":"//ecc_dataTime(3:4)//":00"

      ecc_validityDateTime = ecc_validityDate(1:4)//"-"//ecc_validityDate(5:6)//"-"//ecc_validityDate(7:8) &
        &                 //"T"//ecc_validityTime(1:2)//":"//ecc_validityTime(3:4)//":00"

      ! Lengths of the above strings (should be 19)
      ecc_dataDateTime_length     = LEN_TRIM(ecc_dataDateTime)
      ecc_validityDateTime_length = LEN_TRIM(ecc_validityDateTime)

      ! Indicate that inquiry was (probably) successful
      ecc_successful = .TRUE.

    ENDIF

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_info_on_datetime

  !--------------------------------------------------------------

  !>
  !! @brief Get values of horizontal-grid-related keys
  !!
  SUBROUTINE ecc_get_info_on_horizontal_grid(ecc_msgid, ecc_gridDefinitionTemplateNumber,  &
    &                                        ecc_numberOfDataPoints, ecc_numberOfGridUsed, &
    &                                        ecc_numberOfGridInReference, ecc_uuidOfHGrid, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Value of GRIB key: gridDefinitionTemplateNumber
    INTEGER, INTENT(OUT) :: ecc_gridDefinitionTemplateNumber

    !> Value of GRIB key: numberOfDataPoints
    INTEGER(KIND=i8), INTENT(OUT) :: ecc_numberOfDataPoints

    !> Value of GRIB key: numberOfGridUsed
    INTEGER, INTENT(OUT) :: ecc_numberOfGridUsed

    !> Value of GRIB key: numberOfGridInReference
    INTEGER, INTENT(OUT) :: ecc_numberOfGridInReference

    !> Value of GRIB key: uuidOfHGrid
    ! This value is a 128-bit UUID.
    ! ecCodes provides this value in different types/formats.
    ! Here, we store the 16 bytes of the UUID
    ! in a character array of size 16.
    CHARACTER(LEN=1), INTENT(OUT) :: ecc_uuidOfHGrid(16)

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> Local ecCodes long integer value
    INTEGER(KIND=ECC_kindOfLong) :: ecc_long_value

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_info_on_horizontal_grid'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_gridDefinitionTemplateNumber = ECC_NULL
    ecc_numberOfDataPoints           = INT(ECC_NULL, KIND=i8)
    ecc_numberOfGridUsed             = ECC_NULL
    ecc_numberOfGridInReference      = ECC_NULL
    ecc_uuidOfHGrid(:)               = ' '
    ecc_successful                   = .FALSE.

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! Get identifier of the grid definition template:
    ! (The key 'gridDefinitionTemplateNumber' is an integral part of grid definition section 3.
    ! Therefore, we do not check for its existence first.)
    CALL codes_get_int(msgid=ecc_msgid_local, key='gridDefinitionTemplateNumber', &
      &                value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_gridDefinitionTemplateNumber = INT(ecc_integer_value)

    ! The number of data points (i.e. grid points) is an integral part of grid definition section 3:
    CALL codes_get_long(msgid=ecc_msgid_local, key='numberOfDataPoints', value=ecc_long_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_numberOfDataPoints = INT(ecc_long_value, KIND=i8)

    ! The GRIB keys used in the following are unique to the general-unstructured-horizontal-grid template (3.101).
    ! Only if the horizontal grid is of this type, this subroutine may return "success"!
    IF (ecc_gridDefinitionTemplateNumber == 101) THEN

       ! Get the consecutive number of the specific grid (if it is an official grid):
       CALL codes_get_int(msgid=ecc_msgid_local, key='numberOfGridUsed', value=ecc_integer_value, status=ecc_status)

       IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

       ecc_numberOfGridUsed = INT(ecc_integer_value)

       ! Get the relevant grid element (cells, edges or vertices):
       CALL codes_get_int(msgid=ecc_msgid_local, key='numberOfGridInReference', value=ecc_integer_value, status=ecc_status)

       IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

       ecc_numberOfGridInReference = INT(ecc_integer_value)

       ! It has to be one of 1 (cells), 2 (vertices) or 3 (edges)
       IF (ALL([ECC_GRID_ELEMENT_CELL, ECC_GRID_ELEMENT_VERTEX, ECC_GRID_ELEMENT_EDGE] /= ecc_numberOfGridInReference)) THEN
         CALL warning(modname//routine, "The value of key 'numberOfGridInReference' is invalid (valid: 1, 2, or 3)")
         RETURN
       ENDIF

       ! Get the 128-bit fingerprint of the horizontal grid:
       CALL codes_get_byte_array(msgid=ecc_msgid_local, key='uuidOfHGrid', value=ecc_uuidOfHGrid, status=ecc_status)

       IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

       ! Finally, indicate that inquiry was (probably) successful
       ecc_successful = .TRUE.

    ENDIF ! IF (ecc_gridDefinitionTemplateNumber == 101)

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_info_on_horizontal_grid

  !--------------------------------------------------------------

  !>
  !! @brief Get values of keys associated with the vertical grid
  !!
  SUBROUTINE ecc_get_info_on_vertical_grid(ecc_msgid, ecc_genVertHeightCoords, ecc_NV, ecc_nlev,           &
    &                                      ecc_uuidOfVGrid, ecc_typeOfFirstFixedSurface,                   &
    &                                      ecc_valueOfFirstFixedSurface, ecc_firstFixedSurface_is_missing, &
    &                                      ecc_typeOfSecondFixedSurface, ecc_valueOfSecondFixedSurface,    &
    &                                      ecc_secondFixedSurface_is_missing, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> General vertical height coordinates?
    LOGICAL, INTENT(OUT) :: ecc_genVertHeightCoords

    !> Number of vertical coordinates for ecc_genVertHeightCoords = .FALSE.
    INTEGER, INTENT(OUT) :: ecc_NV

    !> Number of vertical levels for ecc_genVertHeightCoords = .TRUE.
    INTEGER, INTENT(OUT) :: ecc_nlev

    !> UUID of vertical grid for ecc_genVertHeightCoords = .TRUE.
    ! This value is a 128-bit UUID.
    ! ecCodes provides this value in different types/formats.
    ! Here, we store the 16 bytes of the UUID
    ! in a character array of size 16.
    CHARACTER(LEN=1), INTENT(OUT) :: ecc_uuidOfVGrid(16)

    !> Type of first fixed surface
    INTEGER, INTENT(OUT) :: ecc_typeOfFirstFixedSurface

    !> Value of first fixed surface
    REAL(wp), INTENT(OUT) :: ecc_valueOfFirstFixedSurface

    !> Flag to indicate if first fixed surface is missing
    LOGICAL, INTENT(OUT) :: ecc_firstFixedSurface_is_missing

    !> Type of second fixed surface
    INTEGER, INTENT(OUT) :: ecc_typeOfSecondFixedSurface

    !> Value of second fixed surface
    REAL(wp), INTENT(OUT) :: ecc_valueOfSecondFixedSurface

    !> Flag to indicate if second fixed surface is missing
    LOGICAL, INTENT(OUT) :: ecc_secondFixedSurface_is_missing

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> Local ecCodes float value
    REAL(KIND=ECC_kindOfFloat) :: ecc_float_value

    !> Flag to indicate if GRIB key is defined
    INTEGER(KIND=ECC_kindOfInt) :: ecc_is_defined

    !> Flag to indicate if specific GRIB keys are defined
    INTEGER(KIND=ECC_kindOfInt) :: ecc_is_nlev_defined, ecc_is_uuidOfVGrid_defined

    !> Local scale factor
    INTEGER(KIND=ECC_kindOfInt) :: ecc_scaleFactor

    !> Local scaled value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_scaledValue

    !> Flag to indicate if GRIB key is set to missing
    INTEGER(KIND=ECC_kindOfInt) :: ecc_is_missing

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_info_on_vertical_grid'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_genVertHeightCoords           = .FALSE.
    ecc_NV                            = ECC_NULL
    ecc_nlev                          = ECC_NULL
    ecc_uuidOfVGrid(:)                = ' '
    ecc_typeOfFirstFixedSurface       = ECC_MISSING
    ! Use default value of CDI library:
    ecc_valueOfFirstFixedSurface      = 0.0_wp
    ecc_firstFixedSurface_is_missing  = .TRUE.
    ecc_typeOfSecondFixedSurface      = ECC_MISSING
    ! Use default value of CDI library:
    ecc_valueOfSecondFixedSurface     = 0.0_wp
    ecc_secondFixedSurface_is_missing = .TRUE.
    ecc_successful                    = .FALSE.

    ecc_is_nlev_defined        = ECC_NULL
    ecc_is_uuidOfVGrid_defined = ECC_NULL

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! The key 'NV' is an integral part of product definition section 4.
    ! Therefore, we do not check if it exists.
    CALL codes_get_int(msgid=ecc_msgid_local, key='NV', value=ecc_integer_value, status=ecc_status)
    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_NV = INT(ecc_integer_value)

    ! The key 'nlev' is defined for general vertical height coordinates only.
    ! ecCodes assigns a float value to this key.
    CALL codes_is_defined(msgid=ecc_msgid_local, key='nlev', is_defined=ecc_is_nlev_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_nlev_defined == 1_ECC_kindOfInt) THEN

      CALL codes_get_real4(msgid=ecc_msgid_local, key='nlev', value=ecc_float_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_nlev = INT(ecc_float_value)

    ENDIF

    ! The key 'uuidOfVGrid' is defined for general vertical height coordinates only.
    CALL codes_is_defined(msgid=ecc_msgid_local, key='uuidOfVGrid', is_defined=ecc_is_uuidOfVGrid_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_uuidOfVGrid_defined == 1_ECC_kindOfInt) THEN

      CALL codes_get_byte_array(msgid=ecc_msgid_local, key='uuidOfVGrid', value=ecc_uuidOfVGrid, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ENDIF

    ! If both keys, 'nlev' and 'uuidOfVGrid', are defined,
    ! we should deal with general vertical height coordinates
    ecc_genVertHeightCoords = (ecc_is_nlev_defined == 1_ECC_kindOfInt) .AND. (ecc_is_uuidOfVGrid_defined == 1_ECC_kindOfInt)

    ! Check if fixed surface keys are defined (if not, return):

    CALL codes_is_defined(msgid=ecc_msgid_local, key='typeOfFirstFixedSurface', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'typeOfFirstFixedSurface' is not defined")
      RETURN
    ENDIF

    CALL codes_is_defined(msgid=ecc_msgid_local, key='typeOfSecondFixedSurface', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'typeOfSecondFixedSurface' is not defined")
      RETURN
    ENDIF

    ! Check if fixed surface types are not set to MISSING = 255

    CALL codes_is_missing(msgid=ecc_msgid_local, key='typeOfFirstFixedSurface', is_missing=ecc_is_missing, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_firstFixedSurface_is_missing = (ecc_is_missing == 1_ECC_kindOfInt)

    CALL codes_is_missing(msgid=ecc_msgid_local, key='typeOfSecondFixedSurface', is_missing=ecc_is_missing, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_secondFixedSurface_is_missing = (ecc_is_missing == 1_ECC_kindOfInt)

    ! Get types and values of fixed surfaces:

    IF (.NOT. ecc_firstFixedSurface_is_missing) THEN

      CALL codes_get_int(msgid=ecc_msgid_local, key='typeOfFirstFixedSurface', value=ecc_integer_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_typeOfFirstFixedSurface = INT(ecc_integer_value)

      CALL codes_is_missing(msgid=ecc_msgid_local, key='scaleFactorOfFirstFixedSurface', &
        &                   is_missing=ecc_is_missing, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      IF (ecc_is_missing == 0_ECC_kindOfInt) THEN

        CALL codes_get_int(msgid=ecc_msgid_local, key='scaleFactorOfFirstFixedSurface', &
          &                value=ecc_scaleFactor, status=ecc_status)

        IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

        ! Note that the scaled value is an unsigned integer of 4 bytes.
        ! In principle, its range of (absolute) values exceeds the range of absolute values
        ! of a signed 4-byte-Fortran-integer.
        ! However, we assume that a case of exeedance is rather unlikely.
        CALL codes_get_int(msgid=ecc_msgid_local, key='scaledValueOfFirstFixedSurface', &
          &                value=ecc_scaledValue, status=ecc_status)

        IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

        ! Original value = (scaled value) * 10**[-(scale factor)]
        !
        ! Note that the scale factor can be signed according to its ecCodes implementation.

        IF (ecc_scaleFactor >= 0_ECC_kindOfInt) THEN

          ecc_valueOfFirstFixedSurface = REAL(ecc_scaledValue, KIND=wp) * ecc_downscale_factor(INT(ecc_scaleFactor))

        ELSE

          ecc_valueOfFirstFixedSurface = REAL(ecc_scaledValue, KIND=wp) * ecc_upscale_factor(ABS(INT(ecc_scaleFactor)))

        ENDIF

      ENDIF ! IF (ecc_is_missing == 0_ECC_kindOfInt)

    ENDIF ! IF (.NOT. ecc_firstFixedSurface_is_missing)

    IF (.NOT. ecc_secondFixedSurface_is_missing) THEN

      CALL codes_get_int(msgid=ecc_msgid_local, key='typeOfSecondFixedSurface', value=ecc_integer_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_typeOfSecondFixedSurface = INT(ecc_integer_value)

      CALL codes_is_missing(msgid=ecc_msgid_local, key='scaleFactorOfSecondFixedSurface', &
        &                   is_missing=ecc_is_missing, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      IF (ecc_is_missing == 0_ECC_kindOfInt) THEN

        CALL codes_get_int(msgid=ecc_msgid_local, key='scaleFactorOfSecondFixedSurface', &
          &                value=ecc_scaleFactor, status=ecc_status)

        IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

        ! Note that the scaled value is an unsigned integer of 4 bytes.
        ! In principle, its range of (absolute) values exceeds the range of absolute values
        ! of a signed 4-byte-Fortran-integer.
        ! However, we assume that a case of exeedance is rather unlikely.
        CALL codes_get_int(msgid=ecc_msgid_local, key='scaledValueOfSecondFixedSurface', &
          &                value=ecc_scaledValue, status=ecc_status)

        IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

        ! Original value = (scaled value) * 10**[-(scale factor)]
        !
        ! Note that the scale factor can be signed according to its ecCodes implementation.

        IF (ecc_scaleFactor >= 0_ECC_kindOfInt) THEN

          ecc_valueOfSecondFixedSurface = REAL(ecc_scaledValue, KIND=wp) * ecc_downscale_factor(INT(ecc_scaleFactor))

        ELSE

          ecc_valueOfSecondFixedSurface = REAL(ecc_scaledValue, KIND=wp) * ecc_upscale_factor(ABS(INT(ecc_scaleFactor)))

        ENDIF

      ENDIF ! IF (ecc_is_missing == 0_ECC_kindOfInt)

    ENDIF ! IF (.NOT. ecc_secondFixedSurface_is_missing)

    ! Indicate that inquiry was (probably) successful
    ecc_successful = .TRUE.

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_info_on_vertical_grid

!--------------------------------------------------------------

!>
  !! @brief Get values of GRIB keys: centre, subCentre
  !!
  SUBROUTINE ecc_get_info_on_generating_centre(ecc_msgid, ecc_centre, ecc_subCentre, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Value of GRIB key: centre
    INTEGER, INTENT(OUT) :: ecc_centre

    !> Value of GRIB key: subCentre
    INTEGER, INTENT(OUT) :: ecc_subCentre

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_info_on_generating_centre'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_centre     = ECC_NULL
    ecc_subCentre  = ECC_NULL
    ecc_successful = .FALSE.

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! The centre keys are integral part of GRIB section 1.
    ! Therefore, we do not ensure that they are defined.

    CALL codes_get_int(msgid=ecc_msgid_local, key='centre', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_centre = INT(ecc_integer_value)

    CALL codes_get_int(msgid=ecc_msgid_local, key='subCentre', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_subCentre = INT(ecc_integer_value)

    ! Indicate that inquiry was (probably) successful
    ecc_successful = .TRUE.

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_info_on_generating_centre

  !--------------------------------------------------------------

  !>
  !! @brief Get values of GRIB keys associated with the generating process or production in general
  !!
  SUBROUTINE ecc_get_info_on_generating_process(ecc_msgid, ecc_typeOfProcessedData, ecc_typeOfGeneratingProcess, &
    &                                           ecc_backgroundProcess, ecc_generatingProcessIdentifier, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Value of GRIB key: typeOfProcessedData
    INTEGER, INTENT(OUT) :: ecc_typeOfProcessedData

    !> Value of GRIB key: typeOfGeneratingProcess
    INTEGER, INTENT(OUT) :: ecc_typeOfGeneratingProcess

    !> Value of GRIB key: backgroundProcess
    INTEGER, INTENT(OUT) :: ecc_backgroundProcess

    !> Value of GRIB key: generatingProcessIdentifier
    INTEGER, INTENT(OUT) :: ecc_generatingProcessIdentifier

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Flag to indicate if GRIB key is defined
    INTEGER(KIND=ECC_kindOfInt) :: ecc_is_defined

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_info_on_generating_process'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_typeOfProcessedData         = ECC_NULL
    ecc_typeOfGeneratingProcess     = ECC_NULL
    ecc_backgroundProcess           = ECC_NULL
    ecc_generatingProcessIdentifier = ECC_NULL
    ecc_successful                  = .FALSE.

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! The key 'typeOfProcessedData' is an integral part of GRIB section 1.
    ! Therefore, we do not ensure that it is defined.
    CALL codes_get_int(msgid=ecc_msgid_local, key='typeOfProcessedData', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_typeOfProcessedData = INT(ecc_integer_value)


    CALL codes_is_defined(msgid=ecc_msgid_local, key='typeOfGeneratingProcess', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'typeOfGeneratingProcess' is not defined")
      RETURN
    ENDIF

    CALL codes_get_int(msgid=ecc_msgid_local, key='typeOfGeneratingProcess', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_typeOfGeneratingProcess = INT(ecc_integer_value)


    CALL codes_is_defined(msgid=ecc_msgid_local, key='backgroundProcess', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'backgroundProcess' is not defined")
      RETURN
    ENDIF

    CALL codes_get_int(msgid=ecc_msgid_local, key='backgroundProcess', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_backgroundProcess = INT(ecc_integer_value)


    CALL codes_is_defined(msgid=ecc_msgid_local, key='generatingProcessIdentifier', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'generatingProcessIdentifier' is not defined")
      RETURN
    ENDIF

    CALL codes_get_int(msgid=ecc_msgid_local, key='generatingProcessIdentifier', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_generatingProcessIdentifier = INT(ecc_integer_value)

    ! Indicate that inquiry was (probably) successful
    ecc_successful = .TRUE.

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_info_on_generating_process

  !--------------------------------------------------------------

  !>
  !! @brief Get values of keys associated with data representation and bitmap
  !!
  SUBROUTINE ecc_get_info_on_data_representation(ecc_msgid, ecc_dataRepresentationTemplateNumber, &
    &                                            ecc_bitsPerValue, ecc_referenceValue,            &
    &                                            ecc_bitmapPresent, ecc_numberOfValues, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Value of GRIB key: dataRepresentationTemplateNumber
    INTEGER, INTENT(OUT) :: ecc_dataRepresentationTemplateNumber

    !> Value of GRIB key: bitsPerValue
    INTEGER, INTENT(OUT) :: ecc_bitsPerValue

    !> Value of GRIB key: referenceValue
    REAL(dp), INTENT(OUT) :: ecc_referenceValue

    !> Value of GRIB key: bitmapPresent
    INTEGER, INTENT(OUT) :: ecc_bitmapPresent

    !> Value of GRIB key: numberOfValues
    INTEGER(KIND=i8), INTENT(OUT) :: ecc_numberOfValues

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> Local ecCodes long integer value
    INTEGER(KIND=ECC_kindOfLong) :: ecc_long_value

    !> Local scale factor
    INTEGER(KIND=ECC_kindOfInt) :: ecc_scaleFactor

    !> Local scaled value
    REAL(KIND=ECC_kindOfFloat) :: ecc_scaledValue

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Flag to indicate if GRIB key is defined
    INTEGER(KIND=ECC_kindOfInt) :: ecc_is_defined

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_info_on_data_representation'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_dataRepresentationTemplateNumber = ECC_NULL
    ecc_bitsPerValue   = ECC_NULL
    ecc_bitmapPresent  = ECC_NULL
    ecc_numberOfValues = INT(ECC_NULL, KIND=i8)
    ecc_referenceValue = -999.0_dp
    ecc_successful     = .FALSE.

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! Get identifier of the data representation template:
    ! (The key 'dataRepresentationTemplateNumber' is an integral part of data representation section 5.
    ! Therefore, we do not check for its existence first.)

    CALL codes_get_int(msgid=ecc_msgid_local, key='dataRepresentationTemplateNumber', &
      &                value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_dataRepresentationTemplateNumber = INT(ecc_integer_value)

    ! Not every template for data representation section 5
    ! may contain the key 'bitsPerValue'.
    ! Therefore, we first check if it exists.

    CALL codes_is_defined(msgid=ecc_msgid_local, key='bitsPerValue', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'bitsPerValue' is not defined")
      RETURN
    ENDIF

    CALL codes_get_int(msgid=ecc_msgid_local, key='bitsPerValue', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_bitsPerValue = INT(ecc_integer_value)

    ! Get flag which indicates whether a bitmap is active (=> i.e. missing values are present)

    CALL codes_is_defined(msgid=ecc_msgid_local, key='bitmapPresent', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'bitmapPresent' is not defined")
      RETURN
    ENDIF

    CALL codes_get_int(msgid=ecc_msgid_local, key='bitmapPresent', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_bitmapPresent = INT(ecc_integer_value)

    ! Get reference value of the field

    CALL codes_is_defined(msgid=ecc_msgid_local, key='referenceValue', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'referenceValue' is not defined")
      RETURN
    ENDIF

    CALL codes_get_real4(msgid=ecc_msgid_local, key='referenceValue', value=ecc_scaledValue, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ! If the reference value is defined, the decimal scale factor is defined, too
    CALL codes_get_int(msgid=ecc_msgid_local, key='decimalScaleFactor', value=ecc_scaleFactor, status=ecc_status)

    ! Original value = (scaled value) * 10**[-(scale factor)]
    !
    ! Note that the scale factor can be signed according to its ecCodes implementation.

    IF (ecc_scaleFactor == 0_ECC_kindOfInt) THEN
      ecc_referenceValue = REAL(ecc_scaledValue, KIND=dp)
    ELSE
      ecc_referenceValue = REAL(ecc_scaledValue, KIND=dp) * 10.0_dp**(-INT(ecc_scaleFactor))
    ENDIF

    ! The number of coded values is an integral part of data representation section 5:
    CALL codes_get_long(msgid=ecc_msgid_local, key='numberOfValues', value=ecc_long_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_numberOfValues = INT(ecc_long_value, KIND=i8)

    ! Indicate that inquiry was (probably) successful
    ecc_successful = .TRUE.

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_info_on_data_representation

  !--------------------------------------------------------------

  !>
  !! @brief Get values of keys associated with tile-base products
  !!
  SUBROUTINE ecc_get_info_on_tiles(ecc_msgid, ecc_productDefinitionTemplateNumber, ecc_tileClassification, &
    &                              ecc_totalNumberOfTileAttributePairs, ecc_numberOfUsedSpatialTiles,      &
    &                              ecc_tileIndex, ecc_numberOfUsedTileAttributes, ecc_attributeOfTile,     &
    &                              ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Value of GRIB key: productDefinitionTemplateNumber
    INTEGER, INTENT(IN) :: ecc_productDefinitionTemplateNumber

    !> Value of GRIB key:
    INTEGER, INTENT(OUT) :: ecc_tileClassification

    !> Value of GRIB key:
    INTEGER, INTENT(OUT) :: ecc_totalNumberOfTileAttributePairs

    !> Value of GRIB key:
    INTEGER, INTENT(OUT) :: ecc_numberOfUsedSpatialTiles

    !> Value of GRIB key:
    INTEGER, INTENT(OUT) :: ecc_tileIndex

    !> Value of GRIB key:
    INTEGER, INTENT(OUT) :: ecc_numberOfUsedTileAttributes

    !> Value of GRIB key:
    INTEGER, INTENT(OUT) :: ecc_attributeOfTile

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Flag to indicate if GRIB key is defined
    INTEGER(KIND=ECC_kindOfInt) :: ecc_is_defined

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_info_on_tiles'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_tileClassification              = ECC_MISSING
    ecc_totalNumberOfTileAttributePairs = ECC_NULL
    ecc_numberOfUsedSpatialTiles        = ECC_NULL
    ecc_tileIndex                       = ECC_NULL
    ecc_numberOfUsedTileAttributes      = ECC_NULL
    ecc_attributeOfTile                 = ECC_MISSING
    ecc_successful                      = .FALSE.

    ! Check product-definition-template number
    IF (.NOT. ANY([55, 59, 40455, 40456] == ecc_productDefinitionTemplateNumber)) THEN
      CALL warning(modname//routine, &
        & "The value of key 'productDefinitionTemplateNumber' is invalid (valid: 55, 59, 40455 or 40456)")
      RETURN
    ENDIF

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! For safety reasons, we inquire if one of the GRIB keys,
    ! which is (so far) unique to tile-base products, exists in the metadate.
    CALL codes_is_defined(msgid=ecc_msgid_local, key='tileClassification', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'tileClassification' is not defined")
      RETURN
    ENDIF

    ! Unfortunately, the names of some GRIB keys differ between the official WMO templates
    ! and the local DWD templates.
    ! We start with the common keys:

    CALL codes_get_int(msgid=ecc_msgid_local, key='tileClassification', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_tileClassification = INT(ecc_integer_value)

    CALL codes_get_int(msgid=ecc_msgid_local, key='totalNumberOfTileAttributePairs', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_totalNumberOfTileAttributePairs = INT(ecc_integer_value)

    CALL codes_get_int(msgid=ecc_msgid_local, key='tileIndex', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_tileIndex = INT(ecc_integer_value)

    ! Keys with different names:

    SELECT CASE(ecc_productDefinitionTemplateNumber)
    CASE(55, 59)

      ! Official WMO template

      ! For safety reasons, we inquire if one of the three following keys really exists in the metadatas
      CALL codes_is_defined(msgid=ecc_msgid_local, key='numberOfUsedSpatialTiles', is_defined=ecc_is_defined, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
        CALL warning(modname//routine, "The key 'numberOfUsedSpatialTiles' is not defined")
        RETURN
      ENDIF

      CALL codes_get_int(msgid=ecc_msgid_local, key='numberOfUsedSpatialTiles', value=ecc_integer_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_numberOfUsedSpatialTiles = INT(ecc_integer_value)

      CALL codes_get_int(msgid=ecc_msgid_local, key='numberOfUsedTileAttributes', value=ecc_integer_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_numberOfUsedTileAttributes = INT(ecc_integer_value)

      CALL codes_get_int(msgid=ecc_msgid_local, key='attributeOfTile', value=ecc_integer_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_attributeOfTile = INT(ecc_integer_value)

    CASE(40455, 40456)

      ! Local DWD template

      ! For safety reasons, we inquire if one of the three following keys really exists in the metadatas
      CALL codes_is_defined(msgid=ecc_msgid_local, key='numberOfTiles', is_defined=ecc_is_defined, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
        CALL warning(modname//routine, "The key 'numberOfTiles' is not defined")
        RETURN
      ENDIF

      CALL codes_get_int(msgid=ecc_msgid_local, key='numberOfTiles', value=ecc_integer_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_numberOfUsedSpatialTiles = INT(ecc_integer_value)

      CALL codes_get_int(msgid=ecc_msgid_local, key='numberOfTileAttributes', value=ecc_integer_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_numberOfUsedTileAttributes = INT(ecc_integer_value)

      CALL codes_get_int(msgid=ecc_msgid_local, key='tileAttribute', value=ecc_integer_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_attributeOfTile = INT(ecc_integer_value)

    END SELECT

    ! Indicate that inquiry was (probably) successful
    ecc_successful = .TRUE.

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_info_on_tiles

  !--------------------------------------------------------------

  !>
  !! @brief Get values of keys associated with local-use GRIB section 2
  !!
  SUBROUTINE ecc_get_local_info(ecc_msgid, ecc_localNumberOfExperiment, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Value of GRIB key: localNumberOfExperiment
    INTEGER, INTENT(OUT) :: ecc_localNumberOfExperiment

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Local ecCodes integer value
    INTEGER(KIND=ECC_kindOfInt) :: ecc_integer_value

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Flag to indicate if GRIB key is defined
    INTEGER(KIND=ECC_kindOfInt) :: ecc_is_defined

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_local_info'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_localNumberOfExperiment = ECC_NULL
    ecc_successful              = .FALSE.

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! The local section is "volatile",
    ! so every GRIB key should be checked for existence
    ! before inquiring its value.

    CALL codes_is_defined(msgid=ecc_msgid_local, key='localNumberOfExperiment', is_defined=ecc_is_defined, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    IF (ecc_is_defined == 0_ECC_kindOfInt) THEN
      CALL warning(modname//routine, "The key 'localNumberOfExperiment' is not defined")
      RETURN
    ENDIF

    CALL codes_get_int(msgid=ecc_msgid_local, key='localNumberOfExperiment', value=ecc_integer_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_localNumberOfExperiment = INT(ecc_integer_value)

    ! Indicate that inquiry was (probably) successful
    ecc_successful = .TRUE.

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_local_info

  !--------------------------------------------------------------

  !>
  !! @brief Get data values, associated metadata and statistics
  !!
  SUBROUTINE ecc_get_values(ecc_msgid, ecc_missingValue, ecc_values, ecc_sizeOfValues, &
    &                       ecc_bitsPerValue, ecc_numberOfMissing, ecc_isUniform,      &
    &                       ecc_uniformValue, ecc_min, ecc_max, ecc_avg, ecc_successful)

    !-----------
    ! Arguments
    !-----------

    !> Identifier/handle of GRIB message (GRIB record)
    INTEGER, INTENT(IN) :: ecc_msgid

    !> Missing value (single precision)
    REAL(sp), INTENT(IN) :: ecc_missingValue

    !> Array of data values (single precision)
    REAL(sp), ALLOCATABLE, INTENT(INOUT) :: ecc_values(:)

    !> Size of array of data values
    INTEGER(KIND=i8), INTENT(OUT) :: ecc_sizeOfValues

    !> Precision of data values
    INTEGER, INTENT(OUT) :: ecc_bitsPerValue

    !> Number of missing values within level/layer
    INTEGER(KIND=i8), INTENT(OUT) :: ecc_numberOfMissing

    !> Flag to indicate whether field is uniform within level/layer
    LOGICAL, INTENT(OUT) :: ecc_isUniform

    !> Uniform field value
    REAL(dp), INTENT(OUT) :: ecc_uniformValue

    !> Field statistics: min, max and average
    REAL(wp), INTENT(OUT) :: ecc_min, ecc_max, ecc_avg

    !> Flag to indicate that inquiry was successful
    LOGICAL, INTENT(OUT) :: ecc_successful

    !-----------------
    ! Local variables
    !-----------------

    !> Local status flag
    INTEGER :: status

    !> Local message handle/identifier
    INTEGER(KIND=ECC_kindOfInt) :: ecc_msgid_local

    !> Local ecCodes long integer value
    INTEGER(KIND=ECC_kindOfLong) :: ecc_long_value

    !> Local ecCodes float value
    REAL(KIND=ECC_kindOfFloat) :: ecc_float_value

    !> Status identifier for ecCodes interface
    INTEGER(KIND=ECC_kindOfInt) :: ecc_status

    !> Definition number/identifier of data representation (packing)
    INTEGER :: ecc_dataRepresentationTemplateNumber

    !> Reference value of packed data
    REAL(dp) :: ecc_referenceValue

    !> Flag which indicates presence of bitmap (missing values)
    INTEGER :: ecc_bitmapPresent

    !> GRIB key holding the number of coded field values in level/layer
    INTEGER(KIND=i8) :: ecc_numberOfValues

    !> GRIB key holding the number of horizontal grid points
    INTEGER(KIND=i8) :: ecc_numberOfDataPoints

    !> Local flag to indicate that inquiry was successful
    LOGICAL :: ecc_successful_local

    !> Procedure name
    CHARACTER(LEN=*), PARAMETER :: routine = 'ecc_get_values'

    !----------------------------

    ! Initialize intent-out arguments
    ecc_sizeOfValues    = INT(ECC_NULL, KIND=i8)
    ecc_bitsPerValue    = ECC_NULL
    ecc_numberOfMissing = INT(ECC_NULL, KIND=i8)
    ecc_isUniform       = .FALSE.
    ecc_uniformValue    = -9.0E+33_wp
    ecc_min             = -9.0E+33_wp
    ecc_max             = -9.0E+33_wp
    ecc_avg             = -9.0E+33_wp
    ecc_successful      = .FALSE.

    ! Important note: If the field values turn out to be uniform within the level or layer
    ! (ecc_isUniform = .TRUE.), ecc_values will be returned unchanged!
    ! This means that ecc_values will not contain the uniform value (ecc_uniformValue)!
    ! This is for reasons of efficiency.

#ifdef HAVE_ECCODES

    ecc_msgid_local = INT(ecc_msgid, KIND=ECC_kindOfInt)

    ! First of all, we have to find out if there are missing values within the level or layer:

    CALL ecc_get_info_on_data_representation( &
      & ecc_msgid                            = ecc_msgid,                            & ! in
      & ecc_dataRepresentationTemplateNumber = ecc_dataRepresentationTemplateNumber, & ! out
      & ecc_bitsPerValue                     = ecc_bitsPerValue,                     & ! out
      & ecc_referenceValue                   = ecc_referenceValue,                   & ! out
      & ecc_bitmapPresent                    = ecc_bitmapPresent,                    & ! out
      & ecc_numberOfValues                   = ecc_numberOfValues,                   & ! out
      & ecc_successful                       = ecc_successful_local                  ) ! out

    IF (.NOT. ecc_successful_local) THEN
      CALL warning(modname//routine, "Unable to get info on data representation")
      RETURN
    ENDIF

    ! In order to know how many values are missing, ecCodes provides the derived concept key:
    ! numberOfMissing = numberOfDataPoints - numberOfValues
    ! Unfortunately, this key is only defined for a small subset of data-representation templates.
    ! Therefore, we try to compute it ourselves.
    ! For this, we need the value of the GRIB key numberOfDataPoints (which is an integral part of grid-definition section 3):
    CALL codes_get_long(msgid=ecc_msgid_local, key='numberOfDataPoints', value=ecc_long_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_numberOfDataPoints = INT(ecc_long_value, KIND=i8)

    IF (ecc_numberOfDataPoints < ecc_numberOfValues) THEN
      CALL warning(modname//routine, "The number of horizontal grid points is less than the number of data values")
      RETURN
    ELSEIF((ecc_bitmapPresent == 1) .AND. (ecc_numberOfValues /= ecc_numberOfDataPoints)) THEN
      ecc_numberOfMissing = MAX(0_i8, ecc_numberOfDataPoints - ecc_numberOfValues)
    ELSE
      ecc_numberOfMissing = 0_i8
    ENDIF

    ! Next, we have to find out if the field values are uniform within the level or layer:

    ! A field should have uniform values whithin a level/layer if:
    ! - there is only a reference value, but no data vector (bitsPerValue = 0)
    ! - there are no missing values
    ecc_isUniform = (ecc_bitsPerValue == 0) .AND. (ecc_numberOfMissing < 1_i8)
    IF (ecc_isUniform) ecc_uniformValue = ecc_referenceValue

    ! Get the size, which is necessary to hold the data vector
    CALL codes_get_size_long(msgid=ecc_msgid_local, key='values', size=ecc_long_value, status=ecc_status)

    IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

    ecc_sizeOfValues = INT(ecc_long_value, KIND=i8)

    IF (ecc_sizeOfValues < 1_i8) THEN
      ! If the size is zero, there is no reason to continue
      CALL warning(modname//routine, "The data vector is of size zero")
      RETURN
    ELSEIF (sp /= ECC_kindOfFloat) THEN
      ! The following two kinds have to be equal:
      ! - Kind of 'values' here in ICON: sp
      ! - Kind of 'values' in ecCodes:   kindOfFloat
      CALL finish(modname//routine, "sp /= ECC_kindOfFloat")
    ENDIF

    IF (.NOT. ecc_isUniform) THEN

      ! Getting the field values is necessary only if
      ! they are not uniform within the level or layer:

      ! Set the missing value (in case that a bitmap applies)
      CALL codes_set(msgid=ecc_msgid_local, key='missingValue', value=ecc_missingValue, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      IF (ALLOCATED(ecc_values)) THEN

        ! Return in case of size mismatch
        IF (SIZE(ecc_values) < ecc_sizeOfValues) THEN
          CALL warning(modname//routine, "Size of argument ecc_values is too small to hold the data vector")
          RETURN
        ENDIF

      ELSE

        ! Allocate argument 'ecc_values' for inquired size
        ALLOCATE(ecc_values(ecc_sizeOfValues), STAT=status)
        IF (status /= SUCCESS) CALL finish(modname//routine, "Allocation of ecc_values failed")

      ENDIF ! IF (ALLOCATED(ecc_values))

      ! Get data values
      CALL codes_get_real4_array(msgid=ecc_msgid_local, key='values', value=ecc_values, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ! Finally, get some field statistics: min, max, average
      CALL codes_get_real4(msgid=ecc_msgid_local, key='min', value=ecc_float_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_min = REAL(ecc_float_value, KIND=wp)

      CALL codes_get_real4(msgid=ecc_msgid_local, key='max', value=ecc_float_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_max = REAL(ecc_float_value, KIND=wp)

      CALL codes_get_real4(msgid=ecc_msgid_local, key='average', value=ecc_float_value, status=ecc_status)

      IF (ecc_status /= ECC_SUCCESS) CALL ecc_error_handling(routine_of_occurrence=routine, error=ecc_status)

      ecc_avg = REAL(ecc_float_value, KIND=wp)

    ELSE

      ! If the field is uniform, we just have to set the statistics:
      ecc_min = REAL(ecc_uniformValue, KIND=wp)
      ecc_max = ecc_min
      ecc_avg = ecc_min

    ENDIF ! IF (.NOT. ecc_isUniform)

    ! Indicate that inquiry was (probably) successful
    ecc_successful = .TRUE.

#else

    CALL finish(modname//routine, ECC_ERROR_MESSAGE_NO_API)

#endif

  END SUBROUTINE ecc_get_values

END MODULE mo_eccodes

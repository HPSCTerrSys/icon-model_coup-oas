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

! A list of requests for input.

! The use case for which this class has been developed is this:
!  1. A list is created with `myList = t_InputRequestList()`, and the names (only the names) of the requested fields are added with `request()`.
!  2. Files are read with readFile(). This reads all DATA ASSOCIATED with the requested variable names into memory (already distributing it to the worker PEs to keep memory footprint down).
!  3. User code inspects AND retrieves the fetched DATA.

MODULE mo_input_request_list
    USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_DOUBLE, C_NULL_PTR, C_NULL_CHAR, C_ASSOCIATED

    USE mo_cdi, ONLY: t_CdiIterator, cdiIterator_new, cdiIterator_nextField, cdiIterator_delete, cdiIterator_inqVTime, &
                    & cdiIterator_inqLevelType, cdiIterator_inqLevel, cdiIterator_inqGridId, cdiIterator_inqVariableName, &
                    & gridInqType, gridInqUuid, CDI_UNDEFID, ZAXIS_SURFACE, ZAXIS_GENERIC, ZAXIS_HYBRID, &
                    & ZAXIS_HYBRID_HALF, &
                    & ZAXIS_PRESSURE, ZAXIS_HEIGHT, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_ISENTROPIC, &
                    & ZAXIS_TRAJECTORY, ZAXIS_ALTITUDE, ZAXIS_SIGMA, ZAXIS_MEANSEA, ZAXIS_TOA, ZAXIS_SEA_BOTTOM, &
                    & ZAXIS_ATMOSPHERE, ZAXIS_CLOUD_BASE, ZAXIS_CLOUD_TOP, ZAXIS_ISOTHERM_ZERO, ZAXIS_SNOW, ZAXIS_LAKE_BOTTOM, &
                    & ZAXIS_SEDIMENT_BOTTOM, ZAXIS_SEDIMENT_BOTTOM_TA, ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_MIX_LAYER, &
                    & ZAXIS_REFERENCE, cdiIterator_inqTile, CDI_NOERR, CDI_EINVAL, GRID_UNSTRUCTURED, t_CdiParam, &
                    & cdiIterator_inqParamParts, gridInqNumber, gridInqPosition, cdiGribIterator_inqLongValue, t_CdiGribIterator, &
                    & cdiGribIterator_clone, cdiGribIterator_delete, cdiIterator_inqRTime, &
                    & cdiIterator_inqFiletype, FILETYPE_GRB, FILETYPE_GRB2, institutInq, institutInqNamePtr
    USE mo_dictionary, ONLY: t_dictionary
    USE mo_exception, ONLY: message, warning, finish, message_text
    USE mo_grid_config, ONLY: n_dom
    USE mo_impl_constants, ONLY: SUCCESS, vname_len
    USE mo_initicon_config, ONLY: lconsistency_checks
    USE mo_parallel_config, ONLY: use_omp_input, process_stride_pgrib
    USE mo_post_op, ONLY: inverse_post_op
    USE mo_input_container, ONLY: t_InputContainer, inputContainer_make
    USE mo_kind, ONLY: wp, dp, sp, i8, i4
    USE mo_lnd_nwp_config, ONLY: tile_list
    USE mo_nwp_sfc_tiles, ONLY: t_tile_att, t_tileinfo_icon, t_tileinfo_grb2, trivial_tile_att
    USE mo_math_types, ONLY: t_Statistics
    USE mo_model_domain, ONLY: t_patch
    USE mo_mpi, ONLY: my_process_is_mpi_workroot, p_bcast, process_mpi_root_id, p_comm_work,  &
      &               p_pe, p_isEqual, p_mpi_wtime, p_isend, p_irecv, p_wait, p_gather,       &
      &               my_process_is_work, num_work_procs, get_my_mpi_work_id, MPI_REQUEST_NULL
    USE mo_run_config, ONLY: msg_level
    USE mo_time_config, ONLY: time_config
    USE mo_util_string, ONLY: real2string, int2string, toCharacter, c2f_char, tolower, one_of
    USE mo_util_table, ONLY: t_table, initialize_table, add_table_column, set_table_entry, print_table, finalize_table
    USE mo_util_uuid_types, ONLY: t_uuid, uuid_string_length
    USE mo_util_uuid, ONLY: uuid_unparse, uuid_parse, char2uuid, OPERATOR(==)
    USE mtime, ONLY: datetime, timedelta, newDatetime, datetimeToString, newTimedelta, timedeltaToString, deallocateDatetime, &
                   & deallocateTimedelta, max_timedelta_str_len, max_datetime_str_len, OPERATOR(-), OPERATOR(+), OPERATOR(==)
    USE mo_eccodes, ONLY: ECC_NULL_HANDLE, ECC_MAX_SHORTNAME_LENGTH, ECC_GRID_ELEMENT_CELL, ECC_GRID_ELEMENT_EDGE, ECC_NULL, &
      &                   ecc_open_file, ecc_close_file, ecc_read_record, ecc_new_from_record, ecc_release,                  &
      &                   ecc_get_info_on_product, ecc_get_info_on_datetime, ecc_get_info_on_horizontal_grid,                &
      &                   ecc_get_info_on_vertical_grid, ecc_get_info_on_generating_centre,                                  &
      &                   ecc_get_info_on_generating_process, ecc_get_info_on_tiles, ecc_get_local_info, ecc_get_values
    USE mo_timer, ONLY: timer_start, timer_stop, timer_file_reading, timer_raw_data_distribution, timer_metadata_decoding, &
      &                 timer_raw_data_decompression, timer_data_distribution, timer_file_inventory

    IMPLICIT NONE

PUBLIC :: t_InputRequestList, InputRequestList_create

    TYPE :: t_InputRequestList
        PRIVATE
        TYPE(t_ListEntry), POINTER :: list(:)
        INTEGER :: variableCount

    CONTAINS
        PROCEDURE :: request => InputRequestList_request    !< Require that a variable be read.
        PROCEDURE :: requestMultiple => InputRequestList_requestMultiple    !< Require that a list of variables. Unlike request() this will request the
                                                                            ! trimmed strings (because it's impossible to pass an array of strings of different LEN).
        PROCEDURE :: readFile => InputRequestList_readFile  !< Scan a file for input data to satisfy the requests.

        PROCEDURE :: getLevels => InputRequestList_getLevels    !< Get the count AND height values (elevation/presure/whatever) of all the levels encountered IN the file.

        !> The `fetchXXX()` methods simply RETURN FALSE IF the DATA could NOT be fetched entirely.
        !> This IS an atomic operation: Either the entire output array IS overwritten OR it IS NOT touched at all.
        PROCEDURE :: fetch2d => InputRequestList_fetch2d
        PROCEDURE :: fetch3d => InputRequestList_fetch3d
        PROCEDURE :: fetchSurface => InputRequestList_fetchSurface  !No level given, fail IF there are several levels.
        PROCEDURE :: fetchTiled2d => InputRequestList_fetchTiled2d
        PROCEDURE :: fetchTiled3d => InputRequestList_fetchTiled3d
        PROCEDURE :: fetchTiledSurface => InputRequestList_fetchTiledSurface  !No level given, fail IF there are several levels.

        !> The `fetchRequiredXXX()` methods will CALL `finish()` IF there are holes IN the DATA that was READ.
        PROCEDURE :: fetchRequired2d => InputRequestList_fetchRequired2d
        PROCEDURE :: fetchRequired3d => InputRequestList_fetchRequired3d
        PROCEDURE :: fetchRequiredSurface => InputRequestList_fetchRequiredSurface
        PROCEDURE :: fetchRequiredTiled2d => InputRequestList_fetchRequiredTiled2d
        PROCEDURE :: fetchRequiredTiled3d => InputRequestList_fetchRequiredTiled3d
        PROCEDURE :: fetchRequiredTiledSurface => InputRequestList_fetchRequiredTiledSurface

        PROCEDURE :: printInventory => InputRequestList_printInventory
        PROCEDURE :: checkRuntypeAndUuids => InputRequestList_checkRuntypeAndUuids

        PROCEDURE :: findIconName => InputRequestList_findIconName    !< Retrieve a t_ListEntry for the given ICON variable name if it exists already.

        PROCEDURE :: destruct => InputRequestList_destruct  !< Destructor.

        PROCEDURE, PRIVATE :: checkRequests => InputRequestList_checkRequests   !< Check that all processes IN the communicator have the same view on which variables are needed.
        PROCEDURE, PRIVATE :: translateNames => InputRequestList_translateNames !< Recalculates the translatedNames of all list entries using the given dictionary.
        PROCEDURE, PRIVATE :: findTranslatedName => InputRequestList_findTranslatedName    !< As findIconName, but uses the translatedVarName.
        PROCEDURE, PRIVATE :: sendStopMessage => InputRequestList_sendStopMessage
        PROCEDURE, PRIVATE :: sendFieldMetadata => InputRequestList_sendFieldMetadata
        PROCEDURE, PRIVATE :: receiveFieldMetadata => InputRequestList_receiveFieldMetadata
        PROCEDURE, PRIVATE :: isRecordValid => InputRequestList_isRecordValid
        PROCEDURE, PRIVATE :: nextField => InputRequestList_nextField
        !> Conceptual copy of InputRequestList_readFile for parallel decoding
        !> via direct use of GRIB library 'ecCodes'
        PROCEDURE :: readFile_grib => InputRequestList_readFile_grib
        !> Conceptual copy of InputRequestList_isRecordValid for use in InputRequestList_readFile_grib
        PROCEDURE, PRIVATE :: isRecordValid_grib => InputRequestList_isRecordValid_grib
    END TYPE

PRIVATE

    !> CDI default missing value
    !  (9999 would be the default missing value of ecCodes)
    REAL(wp), PARAMETER :: MISSING_VALUE    = -9.0E+33_wp
    REAL(sp), PARAMETER :: MISSING_VALUE_SP = -9.0E+33_sp

    ! These objects are created via findDomainData(, , opt_lcreate = .TRUE.), which will already instanciate an empty container.
    ! On the I/O PE, it IS the job of InputRequestList_isRecordValid() to immediately add a MetadataCache, so that ANY DomainData
    ! object returned by findDomainData() CONTAINS both a valid InputContainer AND a valid MetadataCache.
    TYPE :: t_DomainData
        INTEGER :: jg
        TYPE(t_DomainData), POINTER :: next

        CLASS(t_InputContainer), POINTER :: container
        TYPE(t_MetadataCache), POINTER :: metadata  !< Some metadata connected with the variable, which IS only used for consistency checking AND printing of the inventory table.
        TYPE(t_Statistics) :: statistics
    END TYPE

    TYPE :: t_ListEntry
      !> The name as it has been requested.
      CHARACTER(len=:), ALLOCATABLE :: iconVarName
      !> The name as it is matched against the stored name. This has
      !! dictionary translation, trimming, and case canonicalization
      !! applied.
      CHARACTER(len=:), ALLOCATABLE :: translatedVarName
        TYPE(t_DomainData), POINTER :: domainData   !< A linked list of an InputContainer AND a MetadataCache for each domain. Only accessed via findDomainData().
    END TYPE

    TYPE :: t_MetadataCache
        CHARACTER(KIND = C_CHAR), POINTER :: rtime(:), vtime(:)
        INTEGER :: levelType, gridNumber, gridPosition, runClass, experimentId, generatingProcessType
        TYPE(t_CdiParam) :: param
        TYPE(t_uuid) :: gridUuid

    CONTAINS
        PROCEDURE :: equalTo => MetadataCache_equalTo
        PROCEDURE :: destruct => MetadataCache_destruct
    END TYPE

    !> Linked-list node containing one file-inventory element
    !  For use in InputRequestList_readFile_grib and InputRequestList_isRecordValid_grib
    TYPE :: t_FileInventoryElement
      !> No direct access to components outside this module
      PRIVATE
      !
      ! Metadata:
      !
      !> ICON-internal variable name
      CHARACTER(LEN=vname_len) :: iconVarName = " "
      !> Reference date and time in format:'YYYY-MM-DDThh:mm:ss'
      CHARACTER(LEN=max_datetime_str_len) :: dataDateTime = " "
      !> Validity date and time in format:'YYYY-MM-DDThh:mm:ss'
      CHARACTER(LEN=max_datetime_str_len) :: validityDateTime = " "
      !> Parameter triple
      INTEGER :: discipline        = ECC_NULL
      INTEGER :: parameterCategory = ECC_NULL
      INTEGER :: parameterNumber   = ECC_NULL
      !> CDI level type
      INTEGER :: levelType = ECC_NULL
      !> Consecutive number of horizontal grid (GRIB key: numberOfGridUsed)
      INTEGER :: gridNumber = ECC_NULL
      !> Identifier of horizontal-grid element (cell, edge or vertex, GRIB key: numberOfGridInReference)
      INTEGER :: gridPosition = ECC_NULL
      !> Identifier of background generating process (GRIB key: backgroundProcess)
      INTEGER :: runClass = ECC_NULL
      !> DWD-specific experiment identifier (GRIB key: localNumberOfExperiment)
      INTEGER :: experimentId = ECC_NULL
      !> Identifier of process used to create the product (GRIB key: typeOfGeneratingProcess)
      INTEGER :: generatingProcessType = ECC_NULL
      !> Fingerprint of horizontal grid (GRIB key: uuidOfHGrid)
      TYPE(t_uuid) :: gridUuid
      !
      !> Actual length of variable name
      INTEGER :: iconVarNameLength = 0
      !> Actual length of reference date and time
      INTEGER :: dataDateTimeLength = 0
      !> Actual length of validity date and time
      INTEGER :: validityDateTimeLength = 0
      !
      ! Statistics:
      !
      !> Minimum of field values
      REAL(wp) :: min = MISSING_VALUE
      !> Maximum of field values
      REAL(wp) :: max = MISSING_VALUE
      !> Average of field values
      REAL(wp) :: mean_sum = MISSING_VALUE
      !> Update counter for computation of field average
      INTEGER :: statistics_counter = 0
      !> Level/layer counter
      INTEGER :: level_counter = 0
      !> Counter for levels/layers with uniform field values
      INTEGER :: uniform_level_counter = 0
      !> Maximum tile index
      INTEGER :: maxTileId = -1
      !> Does the field exhibit missing values?
      LOGICAL :: missingValuesPresent = .FALSE.
      !> Did a Work PE that holds an instance of this type
      !> get a GRIB record for decoding?
      LOGICAL :: gotRecord = .FALSE.
      !> Is this GRIB record valid for further processing?
      LOGICAL :: isValid = .FALSE.
      !
      !> Next link in list
      TYPE(t_FileInventoryElement), POINTER :: next => NULL()
    CONTAINS
      PRIVATE
      !> Initialize metadata of an instance of this type
      PROCEDURE :: FileInventoryElement_init_metadata
      !> Initialize statistics of an instance of this type
      PROCEDURE :: FileInventoryElement_init_statistics
      GENERIC   :: init => FileInventoryElement_init_metadata, &
        &                  FileInventoryElement_init_statistics
      !> Initialization in case of a GRIB record
      !> that is not required and ignored
      PROCEDURE :: reject => FileInventoryElement_reject
      !> Update an instance of this type
      PROCEDURE :: update => FileInventoryElement_update
      !> Get metadata of an instance of this type
      PROCEDURE :: FileInventoryElement_get_metadata
      !> Get statistics of an instance of this type
      PROCEDURE :: FileInventoryElement_get_statistics
      GENERIC   :: get => FileInventoryElement_get_metadata, &
        &                 FileInventoryElement_get_statistics
      !> Reset content to default values
      PROCEDURE :: reset => FileInventoryElement_reset
      !> Overload the assignment operator
      PROCEDURE :: FileInventoryElement_assign
      GENERIC   :: ASSIGNMENT(=) => FileInventoryElement_assign
      !> Overload the equality operator: ==
      PROCEDURE :: FileInventoryElement_equal
      GENERIC   :: OPERATOR(==) => FileInventoryElement_equal
    END TYPE t_FileInventoryElement

    !> Linked list for file inventory
    !  For use in InputRequestList_readFile_grib and InputRequestList_isRecordValid_grib
    TYPE :: t_FileInventory
      !> No direct access to components outside this module
      PRIVATE
      !> Counter of inventory-element nodes
      INTEGER :: element_counter = 0
      !> Counter of records in GRIB file
      INTEGER :: record_counter = 0
      !> Counter of GRIB records, which are not required and ignored
      INTEGER :: rejection_counter = 0
      !> Linked list of inventory elements
      TYPE(t_FileInventoryElement), POINTER :: first_element => NULL()
      TYPE(t_FileInventoryElement), POINTER :: last_element  => NULL()
    CONTAINS
      PRIVATE
      !> Add element to inventory
      PROCEDURE :: add => FileInventory_add
      !> Print file inventory
      PROCEDURE :: print => FileInventory_print
      !> Clear the linked list
      PROCEDURE :: destruct => FileInventory_destruct
    END TYPE t_FileInventory

    CHARACTER(*), PARAMETER :: modname = "mo_input_request_list"
    LOGICAL, PARAMETER :: debugModule = .FALSE.

CONTAINS

    !Can't use a type constructor interface t_InputRequestList() since the cray compiler looses the list pointer while returning + assigning the function result.
    FUNCTION InputRequestList_create() RESULT(resultVar)
        TYPE(t_InputRequestList), POINTER :: resultVar

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_create"
        INTEGER :: error, i

        ALLOCATE(resultVar, STAT = error)
        if(error /= SUCCESS) CALL finish(routine, "error allocating memory")

        resultVar%variableCount = 0
        ALLOCATE(resultVar%list(8), STAT = error)
        if(error /= SUCCESS) CALL finish(routine, "error allocating memory")
        DO i = 1, SIZE(resultVar%list, 1)
            resultVar%list(i)%domainData => NULL()
        END DO
    END FUNCTION InputRequestList_create

    FUNCTION findDomainData(listEntry, jg, opt_lcreate) RESULT(resultVar)
        TYPE(t_ListEntry), POINTER, INTENT(INOUT) :: listEntry
        INTEGER, INTENT(IN) :: jg
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lcreate
        TYPE(t_DomainData), POINTER :: resultVar

        CHARACTER(*), PARAMETER :: routine = modname//":findDomainData"
        INTEGER :: error

        IF(.NOT.ASSOCIATED(listEntry)) CALL finish(routine, "assertion failed, listEntry IS NOT ASSOCIATED")

        ! Try to find a preexisting DomainData object.
        resultVar => listEntry%domainData
        DO WHILE (ASSOCIATED(resultVar))
            IF(resultVar%jg == jg) RETURN
            resultVar => resultVar%next
        END DO

        ! Nothing preexisting found, should we create a new one?
        if(PRESENT(opt_lcreate)) THEN
            IF(opt_lcreate) THEN
                ALLOCATE(resultVar, STAT = error)
                IF(error /= SUCCESS) CALL finish(routine, "error allocating memory")
                resultVar%jg = jg
                resultVar%next => listEntry%domainData
                resultVar%container => InputContainer_make()
                resultVar%metadata => NULL()
                CALL resultVar%statistics%reset()
                listEntry%domainData => resultVar
            END IF
        END IF
    END FUNCTION findDomainData

    SUBROUTINE InputRequestList_translateNames(me, opt_dict)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_dictionary), OPTIONAL, INTENT(IN) :: opt_dict

        INTEGER :: i

        DO i = 1, me%variableCount
          IF(PRESENT(opt_dict)) THEN
            me%list(i)%translatedVarName &
                 = tolower(TRIM(opt_dict%get(me%list(i)%iconVarName, &
                 &                           me%list(i)%iconVarName)))
          ELSE
            me%list(i)%translatedVarName = tolower(me%list(i)%iconVarName)
          END IF
        END DO
    END SUBROUTINE InputRequestList_translateNames

    FUNCTION InputRequestList_findIconName(me, fieldName, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: fieldName
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug
        TYPE(t_ListEntry), POINTER :: resultVar

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_findIconName"
        INTEGER :: i
        LOGICAL :: debugInfo
        CHARACTER(:), POINTER :: tempName

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        resultVar => NULL()
        DO i = 1, me%variableCount
            tempName => me%list(i)%iconVarName
            IF(fieldName == tempName) THEN
                IF(debugInfo) CALL message(routine, fieldName//" == "//tempName)
                resultVar => me%list(i)
                RETURN
            ELSE
                IF(debugInfo) CALL message(routine, fieldName//" /= "//tempName)
            END IF
        END DO
    END FUNCTION InputRequestList_findIconName

    FUNCTION InputRequestList_findTranslatedName(me, fieldName) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(INOUT) :: fieldName
        TYPE(t_ListEntry), POINTER :: resultVar

        INTEGER :: i

        fieldname = toLower(fieldName)
        resultVar => NULL()
        DO i = 1, me%variableCount
          IF(fieldName == me%list(i)%translatedVarName) THEN
            resultVar => me%list(i)
            RETURN
          END IF
        END DO
    END FUNCTION InputRequestList_findTranslatedName

    !XXX: This also ensures that the requests have been given in the same order on all processes. Not technically necessary, but easier to
    !     implement and I guess, if the order is not the same, that's a hint that there is a bug somewhere else.
    SUBROUTINE InputRequestList_checkRequests(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_checkRequests"
        INTEGER :: i, j, error, concatenatedSize, curSize, accumulatedSize
        CHARACTER(KIND = C_CHAR), ALLOCATABLE :: concatenatedNames(:)

        !compute the concatenation of all requested variables
        concatenatedSize = 0
        DO i = 1, me%variableCount
          concatenatedSize = concatenatedSize + LEN(me%list(i)%iconVarName) + 1
        END DO
        ALLOCATE(concatenatedNames(concatenatedSize), STAT = error)
        IF(error /= SUCCESS) CALL finish(routine, "memory allocation error")
        accumulatedSize = 0
        DO i = 1, me%variableCount
          curSize = LEN(me%list(i)%iconVarName)
            DO j = 1, curSize
                concatenatedNames(accumulatedSize + j) = me%list(i)%iconVarName(j:j)
            END DO
            concatenatedNames(accumulatedSize + curSize + 1) = C_NULL_CHAR
            accumulatedSize = accumulatedSize + curSize + 1
        END DO

        !check that all processes have the same concatenatedNames string
        IF(.NOT. p_isEqual(concatenatedNames, p_comm_work)) THEN
            print*, "process ", p_pe, " has the variable list: ", concatenatedNames
            CALL finish(routine, "not all processes have the same requests in their t_InputRequestList")
        END IF
    END SUBROUTINE InputRequestList_checkRequests

    SUBROUTINE InputRequestList_request(me, fieldName)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: fieldName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_request"
        INTEGER :: i, listSize, error
        TYPE(t_ListEntry), POINTER :: tempList(:), newEntry

        ! don't add a name twice
        IF(ASSOCIATED(me%findIconName(fieldName))) RETURN

        IF(debugModule .AND. my_process_is_mpi_workroot()) print*, 'Adding request for variable "', fieldName, '"'

        ! ensure space for the new container
        listSize = SIZE(me%list, 1)
        IF(me%variableCount == listSize) THEN
            ALLOCATE(tempList(2*listSize), STAT = error)
            if(error /= SUCCESS) CALL finish(routine, "error allocating memory")
            DO i = 1, listSize
                tempList(i) = me%list(i)
            END DO
            DO i = listSize + 1, 2*listSize
                tempList(i)%domainData => NULL()
            END DO
            DEALLOCATE(me%list)
            me%list => tempList
        END IF

        ! add the entry to our list
        me%variableCount = me%variableCount + 1
        newEntry => me%list(me%variableCount)

        newEntry%iconVarName = fieldName
    END SUBROUTINE InputRequestList_request

    SUBROUTINE InputRequestList_requestMultiple(me, fieldNames)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        CHARACTER(LEN = *), INTENT(IN) :: fieldNames(:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_requestMultiple"
        INTEGER :: i

        DO i = 1, SIZE(fieldNames, 1)
            CALL me%request(TRIM(fieldNames(i)))
        END DO
    END SUBROUTINE InputRequestList_requestMultiple

    SUBROUTINE fail(message, variableName, resultVar)
        CHARACTER(LEN = *), INTENT(IN) :: message
        CHARACTER(LEN = *), INTENT(in) :: variableName
        LOGICAL, INTENT(inout) :: resultVar

        IF(msg_level >= 1) print*, 'invalid record for variable "', variableName, '" encountered: '//message
        resultVar = .FALSE.
    END SUBROUTINE fail

    LOGICAL FUNCTION InputRequestList_isRecordValid(me, iterator, p_patch, level, tileId, variableName, lIsFg) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_CdiIterator) :: iterator
        TYPE(t_patch), INTENT(IN) :: p_patch
        REAL(dp), INTENT(OUT) :: level
        INTEGER, INTENT(OUT) :: tileId
        CHARACTER(:), ALLOCATABLE, INTENT(out) :: variableName
        LOGICAL, INTENT(IN) :: lIsFg

        INTEGER(KIND = C_INT) :: error, gridId, gridType, tileIndex, tileAttribute
        REAL(KIND = C_DOUBLE) :: levelValue
        TYPE(t_MetadataCache), POINTER :: metadata
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        TYPE(t_CdiGribIterator) :: gribIterator
        TYPE(datetime), POINTER :: tempTime, iniTime, startTime
        CHARACTER(:), POINTER :: vtimeString
        CHARACTER(max_datetime_str_len) :: debugDatetimeString
        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_isRecordValid"
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute
        TYPE(t_tileinfo_grb2) :: tileinfo_grb2
        TYPE(t_tileinfo_icon) :: tileinfo_icon
        CHARACTER(KIND = C_CHAR), POINTER :: instName(:)  ! institute name
        INTEGER :: instId                                 ! institute ID
        INTEGER :: generatingCenter, generatingSubCenter
        CHARACTER(KIND=C_CHAR), POINTER :: variableName_(:)

        resultVar = .TRUE.
        variableName_ => cdiIterator_inqVariableName(iterator)
        CALL c2f_char(variableName, variableName_)
        DEALLOCATE(variableName_)
        metadata => MetadataCache_create()
        CALL cdiIterator_inqParamParts(iterator, metadata%param%discipline, metadata%param%category, metadata%param%number)

        !Check the time.
        metadata%vtime => cdiIterator_inqVTime(iterator)
        metadata%rtime => cdiIterator_inqRTime(iterator)
        IF (.NOT. ASSOCIATED(metadata%vtime) .OR. .NOT. ASSOCIATED(metadata%rtime)) THEN
          CALL finish(routine, "Internal error!")
        END IF

        IF(lconsistency_checks) THEN
            vtimeString => toCharacter(metadata%vtime)
            tempTime => newDatetime(vtimeString)

            ALLOCATE(iniTime, STAT = error)
            IF(error /= SUCCESS) CALL fail("memory allocation failure", variableName, resultVar)
            iniTime = time_config%tc_startdate
            IF(lIsFg) THEN
                ! add timeshift to INI-datetime to get true starting time
                ALLOCATE(startTime, STAT = error)
                IF(error /= SUCCESS) CALL fail("memory allocation failure", variableName, resultVar)
                startTime = iniTime + time_config%timeshift%mtime_shift
                IF(.NOT.(tempTime == startTime)) THEN
                    CALL datetimeToString(startTime, debugDatetimeString)
                    CALL fail("vtime of first-guess field ("//vtimeString//") does not match model start time (" &
                             &//TRIM(debugDatetimeString)//")", variableName, resultVar)
                END IF
                DEALLOCATE(startTime)
            ELSE
                IF(.NOT.(tempTime == iniTime)) THEN
                    CALL datetimeToString(iniTime, debugDatetimeString)
                    CALL fail("vtime of analysis field ("//vtimeString//") does not match model initialization time (" &
                             &//TRIM(debugDatetimeString)//")", variableName, resultVar)
                END IF
            END IF
            DEALLOCATE(iniTime)
            DEALLOCATE(vtimeString)
            CALL deallocateDatetime(tempTime)
        END IF

        !We only check the primary (top) level (selector = 0). Usually that's the only one, but GRIB does allow a secondary lower boundary level.
        metadata%levelType = cdiIterator_inqLevelType(iterator, 0)
        SELECT CASE(metadata%levelType)
            !the level types that translate to a single height VALUE
            CASE(ZAXIS_SURFACE, ZAXIS_PRESSURE, ZAXIS_HEIGHT, ZAXIS_DEPTH_BELOW_SEA, ZAXIS_DEPTH_BELOW_LAND, ZAXIS_ALTITUDE, &
                &ZAXIS_REFERENCE, ZAXIS_SNOW)

                levelValue = 0.0
                error = cdiIterator_inqLevel(iterator, 1, outValue1 = levelValue)
                level = REAL(levelValue, dp)
                IF(error /= 0) CALL fail("cdiIterator_inqLevel() failed", variableName, resultVar)
                !TODO: check the zaxis UUID

            !the level types for special levels
            CASE(ZAXIS_TOA, ZAXIS_ATMOSPHERE, ZAXIS_CLOUD_BASE, ZAXIS_CLOUD_TOP, ZAXIS_ISOTHERM_ZERO, &
                &ZAXIS_MEANSEA, ZAXIS_SEA_BOTTOM, ZAXIS_LAKE_BOTTOM, ZAXIS_SEDIMENT_BOTTOM, ZAXIS_SEDIMENT_BOTTOM_TA, &
                &ZAXIS_SEDIMENT_BOTTOM_TW, ZAXIS_MIX_LAYER)

                level = REAL(-metadata%levelType, dp)

            !the known z-axis types that are NOT handled by this code
            CASE(ZAXIS_GENERIC)
                CALL fail("z-axis type ZAXIS_GENERIC is not implemented", variableName, resultVar)
            CASE(ZAXIS_HYBRID)
                CALL fail("z-axis type ZAXIS_HYBRID is not implemented", variableName, resultVar)
            CASE(ZAXIS_HYBRID_HALF)
                CALL fail("z-axis type ZAXIS_HYBRID_HALF is not implemented", variableName, resultVar)
            CASE(ZAXIS_ISENTROPIC)
                CALL fail("z-axis type ZAXIS_ISENTROPIC is not implemented", variableName, resultVar)
            CASE(ZAXIS_TRAJECTORY)
                CALL fail("z-axis type ZAXIS_TRAJECTORY is not implemented", variableName, resultVar)
            CASE(ZAXIS_SIGMA)
                CALL fail("z-axis type ZAXIS_SIGMA is not implemented", variableName, resultVar)

            !fallback to catch future expansions of the list of available z-axis types
            CASE DEFAULT
                CALL fail("unknown z-axis TYPE ("//int2string(metadata%levelType)//")", variableName, resultVar)
        END SELECT

        !Check the grid.
        gridId = cdiIterator_inqGridId(iterator)
        IF(gridId /= CDI_UNDEFID) THEN
            !XXX: I believe, it's enough sanity checking if we check the type and the size of the grid,
            !     I don't want to go into checking the lon/lat for all its vertices here...
            !     A test for the correct grid SIZE IS IMPLICIT IN the t_InputContainer when it selects the scatter pattern to USE.
            gridType = gridInqType(gridId)
            IF(gridType /= GRID_UNSTRUCTURED) THEN
                CALL fail("support for this gridtype is not implemented (CDI grid type = "//TRIM(int2string(gridType))//")", variableName, resultVar)
            ELSE
                CALL gridInqUuid(gridId, metadata%gridUuid%DATA)
                metadata%gridNumber = gridInqNumber(gridID)
                metadata%gridPosition = gridInqPosition(gridID)
            END IF
        ELSE
            CALL fail("couldn't inquire grid ID", variableName, resultVar)
        END IF

        error = cdiIterator_inqTile(iterator, tileIndex, tileAttribute);
        SELECT CASE(error)
            CASE(CDI_NOERR)
              tileinfo_grb2%idx = tileIndex
              tileinfo_grb2%att = tileAttribute
              this_att => tile_list%getTileAtt(tileinfo_grb2)
              tileinfo_icon = this_att%getTileinfo_icon()
              tileId = tileinfo_icon%idx

            CASE(CDI_EINVAL)
              !There IS no tile information connected to this field, so we USE the trivial tileId.
              tileinfo_icon = trivial_tile_att%getTileinfo_icon()
              tileId = tileinfo_icon%idx

            CASE DEFAULT
                CALL finish(routine, "unexpected error while reading tile information")
        END SELECT

        !Fetch some additional metadata.
        metadata%runClass = -1
        metadata%experimentId = -1
        metadata%generatingProcessType = -1
        IF(resultVar) THEN
            gribIterator = cdiGribIterator_clone(iterator)
            IF(C_ASSOCIATED(gribIterator%ptr)) THEN
                metadata%runClass = INT(cdiGribIterator_inqLongValue(gribIterator, "backgroundProcess"))
                metadata%generatingProcessType = INT(cdiGribIterator_inqLongValue(gribIterator, "typeOfGeneratingProcess"))
                !
                ! fetching GRIB2 key "localNumberOfExperiment" is restricted to input data generated by DWD
                generatingCenter = INT(cdiGribIterator_inqLongValue(gribIterator, "centre"))
                generatingSubCenter = INT(cdiGribIterator_inqLongValue(gribIterator, "subCentre"))
                instId = institutInq(generatingCenter, generatingSubcenter, '', '')
                instName => institutInqNamePtr(instId)
                !
                ! Check if instName is associated as CDI returns a NULL pointer if the center is not found within
                ! a CDI internal list (instituteDefaultEntries)
                IF (ASSOCIATED(instName)) THEN
                  IF (TRIM(toCharacter(instName))=="DWD") THEN
                    metadata%experimentId = INT(cdiGribIterator_inqLongValue(gribIterator, "localNumberOfExperiment"))
                  ENDIF
                ENDIF
                CALL cdiGribIterator_delete(gribIterator)
            END IF
        END IF

        !Check whether the metadata of this record IS consistent with the metadata we've already seen for this variable.
        listEntry => me%findTranslatedName(variableName)
        IF(.NOT.ASSOCIATED(listEntry)) THEN
            resultVar = .FALSE.    !We are NOT interested IN this variable.
            CALL metadata%destruct()
            DEALLOCATE(metadata)
            RETURN
        END IF
        IF(resultVar) domainData => findDomainData(listEntry, p_patch%id)
        IF(resultVar .AND. ASSOCIATED(domainData)) resultVar = metadata%equalTo(domainData%metadata)

        !Commit AND cleanup.
        IF(resultVar) THEN
            IF(ASSOCIATED(domainData)) THEN
                CALL metadata%destruct()
                DEALLOCATE(metadata)
            ELSE
                domainData => findDomainData(listEntry, p_patch%id, opt_lcreate = .TRUE.)
                domainData%metadata => metadata  !We don't have a metadata cache yet, so we just remember this one.
            END IF
        ELSE
            !The record was not valid.
            DEALLOCATE(variableName)
            CALL metadata%destruct()
            DEALLOCATE(metadata)
        END IF

    END FUNCTION InputRequestList_isRecordValid

    ! Format of the messages which are broadcasted by the following three FUNCTION:
    !
    ! REAL(dp) :: message(3)
    ! message(1) = length of variable NAME, zero IF this is a stop message
    ! message(2) = level
    ! message(3) = tileId
    !
    ! In the case that the variable NAME length is nonzero, this is followed by another message containing the NAME itself.
    ! Note: The broadcasts IN receiveFieldMetadata() are matched with the broadcasts within sendFieldMetadata() and sendStopMessage().
    SUBROUTINE InputRequestList_sendStopMessage(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_sendStopMessage"
        REAL(dp) :: message(3)

        message(1) = 0.0_dp
        message(2) = 0.0_dp
        message(3) = 0.0_dp
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)
    END SUBROUTINE InputRequestList_sendStopMessage

    SUBROUTINE InputRequestList_sendFieldMetadata(me, level, tileId, variableName)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        REAL(dp), INTENT(in) :: level
        INTEGER, INTENT(in) :: tileId
        CHARACTER(*), INTENT(IN) :: variableName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_sendFieldMetadata"
        REAL(dp) :: message(3)
        CHARACTER(:), POINTER :: tempName
        INTEGER :: error

        message(1) = REAL(LEN(variableName), dp)
        message(2) = level
        message(3) = REAL(tileId, dp)
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)

        IF(debugModule) print*, 'Reading field for variable "'//variableName//'"'
        CALL p_bcast(variableName, process_mpi_root_id, p_comm_work)
    END SUBROUTINE InputRequestList_sendFieldMetadata

    LOGICAL FUNCTION InputRequestList_receiveFieldMetadata(me, level, tileId, variableName) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        REAL(dp), INTENT(INOUT) :: level
        INTEGER, INTENT(INOUT) :: tileId
        CHARACTER(:), ALLOCATABLE, INTENT(out) :: variableName

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_receiveFieldMetadata"
        REAL(dp) :: message(3)
        CHARACTER(:), ALLOCATABLE :: tempName
        INTEGER :: error

        message(1) = 0.0_dp
        message(2) = 0.0_dp
        message(3) = 0.0_dp
        CALL p_bcast(message, process_mpi_root_id, p_comm_work)
        level = message(2)
        tileId = INT(message(3))
        resultVar = message(1) /= 0.0_dp
        IF (resultVar) THEN
          ALLOCATE(CHARACTER(LEN = INT(message(1))) :: variableName, STAT = error)
          IF(error /= SUCCESS) CALL finish(routine, "error allocating memory")
          CALL p_bcast(variableName, process_mpi_root_id, p_comm_work)
        END IF
    END FUNCTION InputRequestList_receiveFieldMetadata

    ! Find the next field that we are interested IN.
    ! This FUNCTION is collective: either all processes RETURN .TRUE. or all RETURN .FALSE. .
    ! When this FUNCTION returns FALSE, THEN there is no further field IN the file.
    !
    ! ignoredRecords IS NOT reset by this FUNCTION, it IS ONLY incremented
    LOGICAL FUNCTION InputRequestList_nextField(me, iterator, p_patch, level, tileId, variableName, ignoredRecords, lIsFg) &
    &RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_CdiIterator), INTENT(IN) :: iterator
        TYPE(t_patch), INTENT(IN) :: p_patch
        REAL(dp), INTENT(OUT) :: level
        INTEGER, INTENT(OUT) :: tileId
        CHARACTER(:), ALLOCATABLE, INTENT(OUT) :: variableName
        INTEGER, INTENT(INOUT) :: ignoredRecords
        LOGICAL, INTENT(IN) :: lIsFg

        resultVar = .FALSE.
        IF(my_process_is_mpi_workroot()) THEN
            ! Scan the file until we find a field that we are interested in.
            DO
                IF(cdiIterator_nextField(iterator) /= 0) THEN
                  IF (.NOT. use_omp_input) THEN
                    CALL me%sendStopMessage()
                  ENDIF
                    RETURN
                ELSE
                    IF(me%isRecordValid(iterator, p_patch, level, tileId, variableName, lIsFg)) THEN
                        IF(ASSOCIATED(me%findTranslatedName(variableName))) THEN
                          IF (.NOT. use_omp_input) THEN
                            !NEC: skip communcation here in VH_OMP case, but do in read routine
                            CALL me%sendFieldMetadata(level, tileId, variableName)
                          ENDIF
                          resultVar = .TRUE.
                          RETURN
                        END IF
                    ELSE
                        ignoredRecords = ignoredRecords + 1
                    END IF
                END IF
            END DO
        ELSE
            resultVar = me%receiveFieldMetadata(level, tileId, variableName)
        END IF
    END FUNCTION InputRequestList_nextField

    SUBROUTINE InputRequestList_readFile(me, p_patch, path, lIsFg, opt_dict)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        TYPE(t_patch), INTENT(IN) :: p_patch
        CHARACTER(LEN = *, KIND = C_CHAR), INTENT(IN) :: path
        TYPE(t_dictionary), OPTIONAL, INTENT(IN) :: opt_dict
        LOGICAL, INTENT(IN) :: lIsFg

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_readFile"
        TYPE(t_CdiIterator) :: iterator
        REAL(dp) :: level
        CHARACTER(KIND = C_CHAR), DIMENSION(:), POINTER :: vtime
        CHARACTER(len = :), ALLOCATABLE :: variableName, variableName_prev
        INTEGER :: i, tileId, recordsRead, recordsIgnored
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        REAL(dp) :: timer(5), savetime
        LOGICAL  :: ret, l_exist

        INTEGER :: iread
        iread = 0

        recordsRead = 0
        recordsIgnored = 0
        timer(1) = p_mpi_wtime()
        timer(2:5) = 0._dp

        CALL me%checkRequests() !sanity checks
        CALL me%translateNames(opt_dict)

        iterator%ptr = C_NULL_PTR
        IF(my_process_is_mpi_workroot()) THEN

            INQUIRE (FILE=path, EXIST=l_exist)
            IF (.NOT.l_exist) THEN
              CALL finish(TRIM(routine),'File is not found: '//TRIM(path))
            ENDIF

            iterator = cdiIterator_new(path)
            IF(.NOT. C_ASSOCIATED(iterator%ptr)) THEN
              CALL finish(routine, "can't open file "//'"'//path//'" for reading')
            END IF
            ! Check whether a Map file for translating fileInputName<=>internalName is required
            ! - a Map File is mandatory, if input is read in GRIB2-Format
            IF (((cdiIterator_inqFiletype(iterator) == FILETYPE_GRB)  .OR.   &
              &  (cdiIterator_inqFiletype(iterator) == FILETYPE_GRB2)) .AND. &
              &  .NOT. PRESENT(opt_dict)) THEN
              CALL finish( routine,                         &
                &  'dictionary missing. It is required when trying to read data in GRIB format.')
            END IF
        END IF
        IF (use_omp_input .AND. my_process_is_mpi_workroot()) THEN
           !NEC_RP: if masterprocess: use readField_omp routine that OMP parallelizes read, statistics and distribution
           DO
               savetime = p_mpi_wtime()
               ret = me%nextField(iterator, p_patch, level, tileId, variableName, recordsIgnored, lIsFg)
               timer(2) = timer(2) + p_mpi_wtime() - savetime
               IF (.NOT. ret) EXIT
               IF (ALLOCATED(variableName_prev)) THEN
               IF (variableName /= variableName_prev) THEN
                  CALL domainData%container%readField_omp(variableName_prev, level, tileId, timer, &
                       p_patch%id, iterator, domainData%statistics, -1)
                  iread = 0
               END IF
               END IF
               recordsRead = recordsRead + 1
               ! We have now found the next field that we are interested in.
               listEntry => me%findTranslatedName(variableName)
               IF(.NOT.ASSOCIATED(listEntry)) CALL finish(routine, &
                 "Assertion failed: Processes have different input request lists!")
               domainData => findDomainData(listEntry, p_patch%id, opt_lcreate = .TRUE.)
               iread = iread + 1
               CALL domainData%container%readField_omp(variableName, level, tileId, timer, &
                  p_patch%id, iterator, domainData%statistics, iread)
               CALL MOVE_ALLOC(variableName, variableName_prev)
           END DO
           IF (ALLOCATED(variableName_prev)) THEN
              CALL domainData%container%readField_omp(variableName_prev, level, tileId, timer, &
                p_patch%id, iterator, domainData%statistics, -2)
              DEALLOCATE(variableName_prev)
           END IF
           CALL me%sendStopMessage()
        ELSE
          !NEC_RP: all other processes use original code
          DO
            savetime = p_mpi_wtime()
            ret = me%nextField(iterator, p_patch, level, tileId, variableName, recordsIgnored, lIsFg)
            timer(2) = timer(2) + p_mpi_wtime() - savetime
            IF (.NOT. ret) EXIT
            recordsRead = recordsRead + 1
            ! We have now found the next field that we are interested IN.
            listEntry => me%findTranslatedName(variableName)
            IF(.NOT.ASSOCIATED(listEntry)) CALL finish(routine, "Assertion failed: Processes have different input request lists!")
            domainData => findDomainData(listEntry, p_patch%id, opt_lcreate = .TRUE.)
            CALL domainData%container%readField(variableName, level, tileId, timer, p_patch%id, iterator, domainData%statistics)
          END DO
        END IF

        timer(1) = p_mpi_wtime() - timer(1)
        IF(my_process_is_mpi_workroot()) THEN
            IF(msg_level > 4) THEN
              WRITE(0, *) routine//": READ "//TRIM(int2string(recordsRead))//" records from file '"//path//"', &
                         &ignoring "//TRIM(int2string(recordsIgnored))//" records"
              WRITE(0, '(3(a,f10.5),a)') ' Timer report: Total ', timer(1), ' s, Read metadata ', timer(2), &
                                       & ' s, Read data ', timer(3), ' s'
              WRITE(0, '(2(a,f10.5),a)') '               Compute statistics ', timer(4), ' s, Distribute data ', timer(5), 's'
            ENDIF
            CALL cdiIterator_delete(iterator)
        END IF
    END SUBROUTINE InputRequestList_readFile

    !>
    !! @brief Reading a GRIB file via direct use of ecCodes
    !!
    SUBROUTINE InputRequestList_readFile_grib(me, grib_file_path, lIsFg, dict, jg, ncells_global, nedges_global, &
      &                                       hgrid_uuid, vgrid_uuid, ana_incr_list, verify_hgrid_uuid,          &
      &                                       verify_vgrid_uuid, verify_ana_incr_list, verbose, timing, inventory)

      !-----------
      ! Arguments
      !-----------

      CLASS(t_InputRequestList), INTENT(INOUT) :: me

      !> Path of GRIB file
      CHARACTER(LEN=*),          INTENT(IN)    :: grib_file_path

      !> Type of products:
      !> - .TRUE.:  First guess
      !> - .FALSE.: Analysis (increments)
      LOGICAL,                   INTENT(IN)    :: lIsFg

      !> Dictionary: ICON variable names <=> ecCodes shortNames
      TYPE(t_dictionary),        INTENT(IN)    :: dict

      !> Patch index
      INTEGER,                   INTENT(IN)    :: jg

      !> Number of grid cells of the global patch
      !> (Long integer to be somewhat prepared for relatively large grids)
      INTEGER(KIND=i8),          INTENT(IN)    :: ncells_global

      !> Number of grid edges of the global patch
      !> (Long integer to be somewhat prepared for relatively large grids)
      INTEGER(KIND=i8),          INTENT(IN)    :: nedges_global

      !> UUID of horizontal grid of ICON
      TYPE(t_uuid),              INTENT(IN)    :: hgrid_uuid

      !> UUID of vertical grid of ICON
      TYPE(t_uuid),              INTENT(IN)    :: vgrid_uuid

      !> List of input fields (with internal variable names)
      !> which may be provided as analysis increments
      CHARACTER(LEN=*),          INTENT(IN)    :: ana_incr_list(:)

      !> Flag to indicate if UUID of horizontal grid has to be verified
      LOGICAL,                   INTENT(IN)    :: verify_hgrid_uuid

      !> Flag to indicate if UUID of vertical grid has to be verified
      LOGICAL,                   INTENT(IN)    :: verify_vgrid_uuid

      !> Flag to indicate whether to verify if metadata of ana_incr_list-fields
      !> conform to analysis increments
      LOGICAL,                   INTENT(IN)    :: verify_ana_incr_list

      !> Flag to indicate verbose messaging
      LOGICAL,                   INTENT(IN)    :: verbose

      !> Flag to indicate measurement of required times
      LOGICAL,                   INTENT(IN)    :: timing

      !> Take a file inventory?
      LOGICAL,                   INTENT(IN)    :: inventory

      !-----------------
      ! Local variables
      !-----------------

      !> (???)
      TYPE(t_ListEntry), POINTER :: listEntry

      !> (???)
      TYPE(t_DomainData), POINTER :: domainData

      !> Status identifier
      INTEGER :: status

      !> Identifier of this Work PE in MPI communicator: p_comm_work
      INTEGER :: my_mpi_work_id

      !> File identifier/handle from ecCodes
      INTEGER :: ecc_ifile

      !> Number/index of a GRIB message (aka GRIB record)
      INTEGER :: ecc_count

      !> GRIB message identifier/handle for interaction with ecCodes
      INTEGER :: ecc_msgid

      !> Length of GRIB message (GRIB record) in bytes
      INTEGER :: ecc_nbytes

      !> Assumed max. length of a GRIB message (GRIB record) in units: 1-byte and 4-bytes
      INTEGER :: ecc_max_record_length_in_byte, ecc_max_record_length

      !> Actual length of a GRIB message (GRIB record) in unit: 4-bytes
      INTEGER :: ecc_record_length

      !> Max. number of status requests for MPI_ISEND and MPI_IRECV
      INTEGER :: max_status_requests

      !> Status request for MPI_ISEND and MPI_IRECV
      INTEGER, ALLOCATABLE :: status_request(:)

      !> Loop index for status requests
      INTEGER :: jsr

      !> Index of status request of Workroot PE currently in use
      INTEGER :: idx_status_requests_curr

      !> For arguments of MPI_ISEND and MPI_IRECV
      INTEGER :: send_to_work_proc, recv_from_work_proc, send_tag, recv_tag, send_count, recv_count

      !> Array to hold GRIB messages
      INTEGER(KIND=i4), ALLOCATABLE :: ecc_record(:,:)

      !> Loop index for data scattering
      INTEGER :: jproc

      !> Value of level for parametric level types
      REAL(dp) :: level

      !> ICON-specific tile identifier
      INTEGER :: tileId

      !> ecCodes shortName
      ! (In contrast to other procedures of this module,
      ! we do not use allocatable strings here,
      ! as this would make things only more complicated and error-prone.)
      CHARACTER(LEN=ECC_MAX_SHORTNAME_LENGTH) :: variableName

      !> Length of variableName
      INTEGER :: variableNameLength

      !> Identifier of sub-grid type: cells or edges
      INTEGER :: subGridId

      !> Number of (global) grid points
      INTEGER(KIND=i8) :: gridSize

      !> Flag to indicate whether field is uniform within level/layer
      LOGICAL :: isUniform

      !> Uniform field value
      REAL(dp) :: uniformValue

      !> Data values defined on grid cells or edges (single precision)
      REAL(sp), ALLOCATABLE :: field(:)

      !> Number of grid elements (either value of ncells_global or nedges_global)
      INTEGER(KIND=i8) :: nelems_global

      !> Flag to indicate a valid GRIB message (GRIB record)
      LOGICAL :: found_match

      !> Flag to indicate that PE is Work Root PE
      LOGICAL :: process_is_workroot

      !> Flag to indicate whether a Work PE participates in decoding
      LOGICAL :: process_does_decoding

      !> Flag to indicate whether it is a Work PE's turn to do something
      LOGICAL :: process_is_next_to_distribute

      !> Size of array of data values
      INTEGER(KIND=i8) :: ecc_sizeOfValues

      !> Number of missing values within level/layer
      INTEGER(KIND=i8) :: ecc_numberOfMissing

      !> Min., max. and average of field values
      REAL(wp) :: ecc_min, ecc_max, ecc_avg

      !> Precision of data values
      INTEGER :: ecc_bitsPerValue

      !> Flag to indicate successful ecCodes inquiries
      LOGICAL :: ecc_successful

      !> Flag to indicate that the end of a GRIB file is reached
      LOGICAL :: ecc_eof

      !> Length of message prefix
      INTEGER :: message_prefix_length, message_prefix_base_length

      !> Prefix for messages
      CHARACTER(LEN=LEN_TRIM(grib_file_path)+50) :: message_prefix

      !> Buffer for distribution of the subsequent metadata among Work PEs
      REAL(dp) :: buffer_real_dp(9)

      !> Metadata currently in progress during distribution among Work PEs
      REAL(dp)         :: level_curr
      INTEGER          :: subGridId_curr
      INTEGER(KIND=i8) :: gridSize_curr
      INTEGER          :: variableNameLength_curr
      INTEGER          :: tileId_curr
      LOGICAL          :: isUniform_curr
      REAL(dp)         :: uniformValue_curr
      LOGICAL          :: found_match_curr
      CHARACTER(LEN=ECC_MAX_SHORTNAME_LENGTH) :: variableName_curr

      !> Record flag:status
      INTEGER(KIND=i4) :: record_flag_status, record_flag_status_curr

      !> File inventory element
      TYPE(t_FileInventoryElement) :: inventory_element

      !> File inventory
      TYPE(t_FileInventory) :: file_inventory

      !> Status flag
      LOGICAL :: successful_local

      !> Shift lower boundary of GRIB message array to store additional info flags
      INTEGER, PARAMETER :: lbound_for_record_flags = -2, &
        &                   idx_record_flag_status  = -2, &
        &                   idx_record_flag_count   = -1, &
        &                   idx_record_flag_length  = 0

      !> Record flag status identifiers
      INTEGER(KIND=i4), PARAMETER :: RECORD_FLAG_STATUS_GOTDATA = 0_i4, &
        &                            RECORD_FLAG_STATUS_EOF     = 1_i4

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":InputRequestList_readFile_grib"

      !----------------------------

      !-----------------------------------------------------------------------------------
      !  - Workroot PE reads data and distributes them to all other Work PEs for decoding
      !  - Other Work PEs decode data
      !  - Used library: ecCodes
      !  - Applicable to file format: GRIB2
      !-----------------------------------------------------------------------------------

      !
      ! Notes:
      !
      ! - The method used in this subroutine was originally developed within the COSMO model
      !   (see https://www.cosmo-model.org/content/default.htm).
      !
      ! - With "ecc_..." we mark variables that are closely related to ecCodes functionalities.
      !
      ! - The overarching (infra)structure of this subroutine is inherited from InputRequestList_readFile.
      !   With "(???)" we mark data structures and procedures whose purpose/functionality/etc.
      !   we could not decrypt from the uncommented source code alone.
      !
      ! - Unfortunately, we cannot use the elements 'metadata' and 'statistics' of derived type 't_DomainData', here.
      !   They are required for printing the file inventory, and have been designed for the Workroot PE
      !   being the sole process that decodes a GRIB record.
      !   As a consequence, we cannot use 'InputRequestList_printInventory' either,
      !   as it does access 'metadata' and 'statistics' explicitly.
      !   In order to nevertheless allow for a file inventory,
      !   new derived types 't_FileInventoryElement' and 't_FileInventory' were implemented.
      !   (Any attempt to adapt 'metadata' and 'statistics' to our needs - if possible at all -
      !   would have resulted in something considerably worse than the current solution.)
      !
      ! - In case of GRIB records, the metadata contain the information
      !   whether a field is uniform within a level/layer:
      !   - There is only a reference value, but no data vector (bitsPerValue = 0)
      !   - There is no bitmap, i.e. no missing values (bitmapPresent = 0)
      !   For instance, the specific cloud water content qc is set to zero above the height
      !   where moist physics are switched off (nonhydrostatic_nml: htop_moist_proc).
      !   In such cases, it is sufficient to distribute only the uniform value among the Work PEs,
      !   instead of distributing an entire horizontal field cross-section of uniform values.
      !   The more model layers between htop_moist_proc and the model top (sleve_nml: top_height),
      !   the larger the savings, in general. (So an ICON-Global setup may potentially benefit more
      !   from this than an ICON-D2 setup.)
      !
      ! - Regarding the precision of the data values packed in a GRIB message,
      !   a bit depth up to bitsPerValue = 32 is allowed here.
      !   However, the data are buffered in single precision for reasons of efficiency.
      !   This may potentially result in a loss of precision for 32 bit data.
      !
      ! - At many places, we have to use "warning" instead of "message" for logging,
      !   as a "message" triggered by other PEs than Workroot are not printed.
      !

      ! The task of this subroutine is meant for Work PEs only
      IF (.NOT. my_process_is_work()) RETURN

      !-------------
      ! Some checks
      !-------------

      IF (jg < 1) THEN
        CALL finish(routine, "Invalid patch index (jg)")
      ELSEIF (ncells_global < 1_i8) THEN
        CALL finish(routine, "Invalid number of grid cells (ncells_global)")
      ELSEIF (nedges_global < 1_i8) THEN
        CALL finish(routine, "Invalid number of grid edges (nedges_global)")
      ELSEIF (process_stride_pgrib > num_work_procs - 1) THEN
        ! process_stride_pgrib has to be element of range [1, num_work_procs - 1].
        ! The lower bound has already been checked in mo_parallel_config: check_parallel_configuration.
        ! As the total number of Work PEs, num_work_procs, is not yet known there, the upper bound is checked here.
        CALL finish(routine, "process_stride_pgrib > num_work_procs-1 is invalid")
      ENDIF

      !--------------
      ! Preparations
      !--------------

      ! Initialize length of message prefix
      message_prefix_length = 0

      ! Unchanging part of message prefix
      message_prefix = "[File: "//TRIM(grib_file_path)//"][GRIB message:"

      ! Length of unchanging part of message prefix
      message_prefix_base_length = LEN_TRIM(message_prefix)

      ! We assume that 'get_my_mpi_work_id()' identifies this Work PE
      ! within MPI communicator: p_comm_work
      my_mpi_work_id = get_my_mpi_work_id()

      ! Is this Work PE the Workroot PE?
      process_is_workroot = my_process_is_mpi_workroot()

      ! Does this Work PE participate in decoding GRIB messages?
      process_does_decoding = (.NOT. process_is_workroot) .AND. (MOD(my_mpi_work_id, process_stride_pgrib) == 0)

      ! Init counter of records in GRIB file
      ecc_count = 0

      ! Initialize ecCodes handles (for safety reasons)
      ecc_msgid = ECC_NULL_HANDLE

      ! Initialize metadata to be distributed among Work PEs
      record_flag_status = -999_i4
      level              = REAL(MISSING_VALUE, KIND=dp)
      subGridId          = -999
      gridSize           = 0_i8
      variableNameLength = 0
      tileId             = -999
      isUniform          = .FALSE.
      uniformValue       = REAL(MISSING_VALUE, KIND=dp)
      found_match        = .FALSE.
      variableName       = " "

      record_flag_status_curr = record_flag_status
      level_curr              = level
      subGridId_curr          = subGridId
      gridSize_curr           = gridSize
      variableNameLength_curr = variableNameLength
      tileId_curr             = tileId
      isUniform_curr          = isUniform
      uniformValue_curr       = uniformValue
      found_match_curr        = found_match
      variableName_curr       = variableName

      ! Check consistency of requested input variables between the involved PEs(???)
      CALL me%checkRequests()

      ! Translate ICON variable names into ecCodes shortNames or vice versa(???)
      CALL me%translateNames(dict)

      ! Initialize content of file-inventory element with default values, just to make sure
      IF (inventory) CALL inventory_element%reset()

      ! We need enough memory to hold one GRIB message (aka GRIB record):
      !  - bitsPerValue = 32 (4 bytes) is the max. allowed precision here for the time being
      !  - The data vector of a GRIB message has either ncells_global or nedges_global entries
      !    (or less, if a bitmap applies)
      !  - Safety margin for the metadata header: 1000 bytes
      !  => MAX(ncells_global, nedges_global) * 4 + 1000
      ! (Note: Although GRIB records may become larger and larger in the future, we cannot use type 'INTEGER(KIND=i8)'.
      ! This is because there is no corresponding MPI interface.)
      ecc_max_record_length_in_byte = INT(MAX(ncells_global, nedges_global)) * 4 + 1000
      ecc_max_record_length         = ecc_max_record_length_in_byte / 4

      ! Max. number of status requests for MPI_ISEND and MPI_IRECV:
      !  - 5 for Workroot PE
      !  - 1 for the other Work PEs
      max_status_requests = MERGE(5, 1, process_is_workroot)

      ! Allocate fields
      IF (process_is_workroot) THEN
        ALLOCATE(ecc_record(lbound_for_record_flags:ecc_max_record_length,max_status_requests), &
          &      status_request(max_status_requests), field(1), STAT=status)
      ELSEIF (process_does_decoding) THEN
        ALLOCATE(ecc_record(lbound_for_record_flags:ecc_max_record_length,max_status_requests), &
          &      status_request(max_status_requests), field(MAX(ncells_global, nedges_global)), STAT=status)
      ELSE
        ! Dummy allocation
        ALLOCATE(ecc_record(1,1), status_request(1), field(1), STAT=status)
      ENDIF
      IF (status /= SUCCESS) &
        & CALL finish(routine, "Allocation of ecc_record, status_request and field failed")

      ecc_record(:,:)   = 0_i4
      status_request(:) = MPI_REQUEST_NULL !???
      field(:)          = 0.0_sp

      ! Initialize flag that indicates whether the end of a GRIB file is reached
      ecc_eof = .FALSE.

      IF (process_is_workroot) THEN

        IF (timing) CALL timer_start(timer_file_reading)

        ! Workroot PE opens GRIB file for reading.
        ! (Note that the following subroutine would initiate a model abort itself,
        ! if something would go wrong with opening the file.)
        CALL ecc_open_file(grib_filename = grib_file_path, & ! in
          &                ecc_mode      = 'r',            & ! in
          &                ecc_ifile     = ecc_ifile       ) ! out

        ! Initialize arguments for MPI_ISEND
        idx_status_requests_curr = 1
        send_to_work_proc        = 0
        send_tag                 = 0
        send_count               = ecc_max_record_length + ABS(lbound_for_record_flags) + 1

        IF (timing) CALL timer_stop(timer_file_reading)

      ELSEIF (process_does_decoding) THEN

        IF (timing) CALL timer_start(timer_raw_data_distribution)

        ! All other Work PEs launch first call of MPI_IRECV (which is non-blocking => use MPI_WAIT)
        idx_status_requests_curr = 1
        recv_from_work_proc      = process_mpi_root_id
        recv_tag                 = my_mpi_work_id
        recv_count               = ecc_max_record_length + ABS(lbound_for_record_flags) + 1

        CALL p_irecv(t_buffer = ecc_record(:,idx_status_requests_curr), & ! inout
          &          p_source = recv_from_work_proc,                    & ! in
          &          p_tag    = recv_tag,                               & ! in
          &          p_count  = recv_count,                             & ! in
          &          comm     = p_comm_work,                            & ! in
          &          request  = status_request(idx_status_requests_curr)) ! optout

        IF (timing) CALL timer_stop(timer_raw_data_distribution)

      ENDIF ! IF (process_is_workroot)

      !-------------------------
      ! Precessing of GRIB file
      !-------------------------

      ! Loop over the messages (records) contained in the GRIB file:
      ! (It might be helpful for its understanding to keep in mind that its effective appearance
      ! is quite different for the Workroot PE, on the one hand, and the Worker PEs, on the other hand.)
      FILE_PROCESSING_LOOP: DO

        IF (process_is_workroot) THEN

          !--------------
          ! Workroot PE:
          !--------------

          ! Increment counter of read in GRIB messages (GRIB records)
          ecc_count = ecc_count + 1

          IF (timing) CALL timer_start(timer_raw_data_distribution)

          ! Wait for send buffer to be free
          CALL p_wait(status_request(idx_status_requests_curr))

          IF (timing) CALL timer_stop(timer_raw_data_distribution)

          IF (timing) CALL timer_start(timer_file_reading)

          IF (.NOT. ecc_eof) THEN

            ! Initially assumed size of record in bytes = max. record length in bytes defined above
            ecc_nbytes = ecc_max_record_length_in_byte

            ! Read next GRIB message (GRIB record) from file
            CALL ecc_read_record(ecc_ifile         = ecc_ifile,         & ! in
              &                  ecc_record        = ecc_record(1:ecc_max_record_length,idx_status_requests_curr), & ! inout
              &                  ecc_nbytes        = ecc_nbytes,        & ! inout
              &                  ecc_record_length = ecc_record_length, & ! out
              &                  ecc_eof           = ecc_eof            ) ! out

            ! Either received a record, or reached the end of the GRIB file?
            record_flag_status = MERGE(RECORD_FLAG_STATUS_EOF, RECORD_FLAG_STATUS_GOTDATA, ecc_eof)

          ELSE

            ! Reached end of GRIB file in some previous loop cycle
            record_flag_status = RECORD_FLAG_STATUS_EOF

          ENDIF ! IF (.NOT. ecc_eof)

          IF (timing) CALL timer_stop(timer_file_reading)

          IF (timing) CALL timer_start(timer_raw_data_distribution)

          ! Set record status flag (index: -2)
          ecc_record(idx_record_flag_status,idx_status_requests_curr) = record_flag_status

          ! Set record count flag (index: -1)
          ecc_record(idx_record_flag_count,idx_status_requests_curr) = INT(ecc_count, KIND=i4)

          ! Set record length "flag" (index: 0)
          ecc_record(idx_record_flag_length,idx_status_requests_curr) = INT(ecc_record_length, KIND=i4)

          ! Lenght of send buffer: length of GRIB record + number of record flags
          send_count = ecc_record_length + ABS(lbound_for_record_flags) + 1

          ! ID of the Work PE, which shall receive
          ! the GRIB message for decoding
          send_to_work_proc = send_to_work_proc + process_stride_pgrib

          ! For the tag, we use the Work PE ID, too
          send_tag = send_to_work_proc

          ! Send GRIB record to the Work PE
          CALL p_isend( &
            & t_buffer      = ecc_record(lbound_for_record_flags:ecc_record_length,idx_status_requests_curr), & ! inout
            & p_destination = send_to_work_proc,                      & ! in
            & p_tag         = send_tag,                               & ! in
            & p_count       = send_count,                             & ! optin
            & comm          = p_comm_work,                            & ! optin
            & request       = status_request(idx_status_requests_curr)) ! optinout

          ! Update index of status request for next send
          idx_status_requests_curr = idx_status_requests_curr + 1

          ! Reset to one if we exceed the max. number of status requests
          IF (idx_status_requests_curr > max_status_requests) idx_status_requests_curr = 1

          ! If all Work PEs got a GRIB record for decoding,
          ! we have to reset the destination of sending
          IF (.NOT. (send_to_work_proc < num_work_procs - process_stride_pgrib)) THEN

            send_to_work_proc = 0
            send_tag          = send_to_work_proc

            ! Free the rest of the status request buffers
            ! if we reached the end of the GRIB file
            IF (ecc_eof) THEN
              DO jsr = 1, max_status_requests
                CALL p_wait(status_request(jsr))
              END DO
            ENDIF

          ENDIF ! IF (.NOT. (send_to_work_proc < num_work_procs - 1))

          IF (timing) CALL timer_stop(timer_raw_data_distribution)

          ! Condition for cycling the processing loop:
          ! The Workroot PE will cycle the processing loop up to this point
          ! about '(num_work_procs - 1)/process_stride_pgrib' times more often
          ! than all the other Work PEs do, in order to send one GRIB message
          ! to the latter for decoding.
          ! Only if the other Work PEs got a GRIB message,
          ! the Workroot PE will not enter the following loop-cycle-condition,
          ! but will advance to the "all-to-all" distribution
          ! of the decoded GRIB messages further below.
          IF (MOD(ecc_count, (num_work_procs - 1)/process_stride_pgrib) /= 0) CYCLE FILE_PROCESSING_LOOP

        ELSEIF (process_does_decoding) THEN

          !---------------------------------------------------------------------------
          ! All the other Work PEs (or a subset of them if process_stride_pgrib > 1):
          !---------------------------------------------------------------------------

          IF (timing) CALL timer_start(timer_raw_data_distribution)

          ! A Work PE has to wait until it got
          ! a GRIB message from the Workroot PE
          CALL p_wait(status_request(idx_status_requests_curr))

          ! First, get flags from record buffer
          record_flag_status = ecc_record(idx_record_flag_status,idx_status_requests_curr)
          ecc_count          = INT(ecc_record(idx_record_flag_count,idx_status_requests_curr))
          ecc_record_length  = INT(ecc_record(idx_record_flag_length,idx_status_requests_curr))

          nelems_global = 0_i8
          found_match   = .FALSE.

          IF (timing) CALL timer_stop(timer_raw_data_distribution)

          SELECT CASE(record_flag_status)
          CASE(RECORD_FLAG_STATUS_GOTDATA)

            IF (timing) CALL timer_start(timer_metadata_decoding)

            ! Check count (index) and length of GRIB record
            IF (ecc_count < 1) THEN
              CALL finish(routine, "Invalid GRIB record count")
            ELSEIF (ecc_record_length < 1) THEN
              CALL finish(routine, "Invalid GRIB record length")
            ENDIF

            ! Start decoding the GRIB message:
            CALL ecc_new_from_record(ecc_record = ecc_record(1:ecc_record_length,idx_status_requests_curr), & ! in
              &                      ecc_msgid  = ecc_msgid                                                 ) ! out

            ! For messages
            IF (verbose) THEN
              ! Remove prefix extension from previous processing loop cycle
              message_prefix(message_prefix_base_length+1:) = " "
              ! Append extension for current loop cycle
              message_prefix        = message_prefix(1:message_prefix_base_length)//" "//TRIM(int2string(ecc_count))//"]"
              message_prefix_length = message_prefix_base_length + LEN_TRIM(message_prefix(message_prefix_base_length+1:))
            ENDIF

            ! Check metadata of GRIB message
            CALL me%isRecordValid_grib(jg                    = jg,                    & ! in
              &                        ecc_msgid             = ecc_msgid,             & ! in
              &                        lIsFg                 = lIsFg,                 & ! in
              &                        verbose               = verbose,               & ! in
              &                        inventory             = inventory,             & ! in
              &                        hgrid_uuid            = hgrid_uuid,            & ! in
              &                        vgrid_uuid            = vgrid_uuid,            & ! in
              &                        ana_incr_list         = ana_incr_list,         & ! in
              &                        verify_hgrid_uuid     = verify_hgrid_uuid,     & ! in
              &                        verify_vgrid_uuid     = verify_vgrid_uuid,     & ! in
              &                        verify_ana_incr_list  = verify_ana_incr_list,  & ! in
              &                        message_prefix        = message_prefix,        & ! inout
              &                        message_prefix_length = message_prefix_length, & ! inout
              &                        level                 = level,                 & ! out
              &                        tileId                = tileId,                & ! out
              &                        variableName          = variableName,          & ! out
              &                        variableNameLength    = variableNameLength,    & ! out
              &                        subGridId             = subGridId,             & ! out
              &                        gridSize              = gridSize,              & ! out
              &                        inventory_element     = inventory_element,     & ! out
              &                        found_match           = found_match            ) ! out

            IF (timing) CALL timer_stop(timer_metadata_decoding)

            IF (timing) CALL timer_start(timer_raw_data_decompression)

            IF (found_match) THEN

              ! Check if grid size from GRIB file fits grid size of ICON model
              IF (gridSize == ncells_global) THEN

                ! Data are defined on grid cells:
                nelems_global = ncells_global

                IF (subGridId /= ECC_GRID_ELEMENT_CELL) THEN

                  ! 'subGridId' = 'numberOfGridInReference' should be 1 for grid cells
                  IF (verbose) THEN
                    WRITE(message_text,'(A,I0)') &
                      & ": Data are defined on grid cells, so numberOfGridInReference should be 1, but is ", subGridId
                    CALL warning(routine, message_prefix(1:message_prefix_length)//message_text)
                  ENDIF

                  ! We reset its value, as it is required later on
                  ! for data distribution among Work PEs
                  subGridId = ECC_GRID_ELEMENT_CELL

                ENDIF ! IF (subGridId /= ECC_GRID_ELEMENT_CELL)

              ELSEIF (gridSize == nedges_global) THEN

                ! Data are defined on grid edges:
                nelems_global = nedges_global

                IF (subGridId /= ECC_GRID_ELEMENT_EDGE) THEN

                  ! 'subGridId' = 'numberOfGridInReference' should be 3 for grid edges
                  IF (verbose) THEN
                    WRITE(message_text,'(A,I0)') &
                      & ": Data are defined on grid edges, so numberOfGridInReference should be 3, but is ", subGridId
                    CALL warning(routine, message_prefix(1:message_prefix_length)//message_text)
                  ENDIF

                  ! We reset its value, as it is required later on
                  ! for data distribution among Work PEs
                  subGridId = ECC_GRID_ELEMENT_EDGE

                ENDIF ! IF (subGridId /= ECC_GRID_ELEMENT_EDGE)

              ELSE

                ! Invalid grid size
                WRITE(message_text,'(3(A,I0),A)') ": Grid size (", gridSize, ") neither fits number of cells (", &
                  & ncells_global, ") nor number of edges (", nedges_global, ") of global patch"
                CALL finish(routine, message_prefix(1:message_prefix_length)//message_text)

              ENDIF ! If valid grid size

              ! Get data values:
              ! (Important note: If the field values turn out to be uniform within the level or layer
              ! (ecc_isUniform = .TRUE.), the argument field will be returned unchanged by the following subroutine!
              ! This means that field will not contain the uniform values (ecc_uniformValue)!
              ! This is for reasons of efficiency.)
              CALL ecc_get_values(ecc_msgid           = ecc_msgid,           & ! in
                &                 ecc_missingValue    = MISSING_VALUE_SP,    & ! in
                &                 ecc_values          = field,               & ! inout
                &                 ecc_sizeOfValues    = ecc_sizeOfValues,    & ! out
                &                 ecc_bitsPerValue    = ecc_bitsPerValue,    & ! out
                &                 ecc_numberOfMissing = ecc_numberOfMissing, & ! out
                &                 ecc_isUniform       = isUniform,           & ! out
                &                 ecc_uniformValue    = uniformValue,        & ! out
                &                 ecc_min             = ecc_min,             & ! out
                &                 ecc_max             = ecc_max,             & ! out
                &                 ecc_avg             = ecc_avg,             & ! out
                &                 ecc_successful      = ecc_successful       ) ! out

              IF (.NOT. ecc_successful) THEN
                CALL finish(routine, message_prefix(1:message_prefix_length) &
                  & //": Unable to get data values")
              ELSEIF (ecc_sizeOfValues /= nelems_global) THEN
                CALL finish(routine, message_prefix(1:message_prefix_length) &
                  & //": Mismatch between grid size and size of data vector")
              ELSEIF (ecc_bitsPerValue > 32) THEN
                CALL finish(routine, message_prefix(1:message_prefix_length) &
                  & //": bitsPerValue > 32 are not supported")
              ENDIF

            ENDIF ! IF (found_match)

            IF (timing) CALL timer_stop(timer_raw_data_decompression)

            IF (timing) CALL timer_start(timer_metadata_decoding)

            ! Finally, release GRIB message
            CALL ecc_release(ecc_msgid=ecc_msgid)

            IF (timing) CALL timer_stop(timer_metadata_decoding)

          CASE(RECORD_FLAG_STATUS_EOF)

            ! There is nothing to do here in this case.

          CASE DEFAULT

            CALL finish(routine, "Invalid record status flag")

          END SELECT ! SELECT CASE(record_flag_status)

        ENDIF ! IF (process_is_workroot)

        ! At this stage every decoding Work PE has to have a record or an EOF,
        ! and should have decoded and checked the metadata and payload.
        ! Now, we go into the loop over all Work PEs to distribute the decoded data:

        IF (timing) CALL timer_start(timer_data_distribution)

        ! Note that the Workroot PE would be 'jproc = 0', so it does not distribute
        ! something in the following loop itself. Nevertheless, it has to participate
        ! in the collective "all-to-all" distribution as such, also in order to obtain necessary information
        ! about each decoded GRIB message.
        DISTRIBUTION_LOOP: DO jproc = 1, num_work_procs - 1

          ! Work PE that is next in line for distributing its data
          process_is_next_to_distribute = (jproc == my_mpi_work_id) .AND. process_does_decoding

          ! First, we have to distribute a number of metadata.
          ! (It seems that the type-bound procedures 'InputRequestList_sendFieldMetadata' and
          ! 'InputRequestList_receiveFieldMetadata' are intended for this purpose.
          ! However, we cannot use them, as they assume the Workroot PE as the sole sender.)
          IF (process_is_next_to_distribute) THEN

            buffer_real_dp(1) = REAL(record_flag_status, KIND=dp)
            buffer_real_dp(2) = MERGE(10.0_dp, -10.0_dp, found_match)
            buffer_real_dp(3) = REAL(subGridId, KIND=dp)
            buffer_real_dp(4) = REAL(gridSize, KIND=dp)
            buffer_real_dp(5) = REAL(tileId, KIND=dp)
            buffer_real_dp(6) = level
            buffer_real_dp(7) = REAL(variableNameLength, KIND=dp)
            buffer_real_dp(8) = MERGE(10.0_dp, -10.0_dp, isUniform)
            buffer_real_dp(9) = uniformValue

          ELSE

            buffer_real_dp(:) = -999.0_dp

          ENDIF ! IF (process_is_next_to_distribute)

          CALL p_bcast(buffer_real_dp, jproc, p_comm_work)

          IF (process_is_next_to_distribute) THEN

            ! The current distributer can just take its own values
            record_flag_status_curr = record_flag_status
            found_match_curr        = found_match
            subGridId_curr          = subGridId
            gridSize_curr           = gridSize
            tileId_curr             = tileId
            level_curr              = level
            variableNameLength_curr = variableNameLength
            isUniform_curr          = isUniform
            uniformValue_curr       = uniformValue

          ELSE

            ! All the others take what was broadcast by the current distributor
            record_flag_status_curr = NINT(buffer_real_dp(1), KIND=i4)
            found_match_curr        = (buffer_real_dp(2) > 0.0_dp)
            subGridId_curr          = NINT(buffer_real_dp(3))
            gridSize_curr           = NINT(buffer_real_dp(4), KIND=i8)
            tileId_curr             = NINT(buffer_real_dp(5))
            level_curr              = buffer_real_dp(6)
            variableNameLength_curr = NINT(buffer_real_dp(7))
            isUniform_curr          = (buffer_real_dp(8) > 0.0_dp)
            uniformValue_curr       = buffer_real_dp(9)

          ENDIF ! IF (process_is_next_to_distribute)

          ! The following depends on there are two status flags only:
          ! - RECORD_FLAG_STATUS_GOTDATA
          ! - RECORD_FLAG_STATUS_EOF
          IF (record_flag_status_curr == RECORD_FLAG_STATUS_EOF) THEN

            ! End-of-file occurred.
            ! (This should overwrite the initialization of 'ecc_eof' with .FALSE.
            ! on all Work PEs at the beginning. It does also "overwrite"
            ! the 'ecc_eof=.TRUE.'-state, the Workroot PE is already in,
            ! but that should do no harm.)
            ecc_eof = .TRUE.

            ! In this case, there is no more to be done,
            ! so we can just skip to the next Work PE,
            ! as there might still be some more valid GRIB records
            ! to be distributed
            CYCLE DISTRIBUTION_LOOP

          ELSEIF (.NOT. found_match_curr) THEN

            ! record_flag_status_curr == RECORD_FLAG_STATUS_GOTDATA & found_match_curr == .FALSE.

            ! This GRIB message is not required and ignored,
            ! so no more to be done here.
            ! We can skip to the next GRIB message distributed by the next Work PE
            CYCLE DISTRIBUTION_LOOP

          ENDIF ! IF (record_flag_status_curr == RECORD_FLAG_STATUS_EOF ...)

          ! If we reached this point,
          ! the current GRIB message is a valid one,
          ! which we have to process further now:

          ! We got the length of the variable name string,
          ! so now the distributor can broadcast the variable name string itself
          IF (process_is_next_to_distribute) THEN

            variableName_curr(1:variableNameLength_curr) = variableName(1:variableNameLength)

          ELSE

            variableName_curr(1:variableNameLength_curr) = " "

          ENDIF

          CALL p_bcast(variableName_curr(1:variableNameLength_curr), jproc, p_comm_work)

          ! (???):

          listEntry => me%findTranslatedName(variableName_curr(1:variableNameLength_curr))

          ! (Note: Within the distribution loop, we do not use 'message_prefix' for the message text,
          ! as the GRIB record counter, it refers to,
          ! does not match the count of the GRIB record that is currently distributed, in general.)
          IF(.NOT. ASSOCIATED(listEntry)) &
            & CALL finish(routine, "Assertion failed: Processes have different input request lists!")

          domainData => findDomainData(listEntry, jg, opt_lcreate=.TRUE.)

          ! Finally, the current Work PE 'jproc' tries to distribute its data to all other Work PEs
          IF (.NOT. process_does_decoding) THEN
            ! See allocation of field above
            nelems_global = 1_i8
          ELSEIF (subGridId_curr == ECC_GRID_ELEMENT_CELL) THEN
            ! Grid cells:
            nelems_global = ncells_global
          ELSEIF (subGridId_curr == ECC_GRID_ELEMENT_EDGE) THEN
            ! Grid edges:
            nelems_global = nedges_global
          ELSE
            ! (Note: Within the distribution loop, we do not use 'message_prefix'
            ! for the message text, as the GRIB record counter, it refers to,
            ! does not match the count of the GRIB record that is currently distributed, in general.)
            CALL finish(routine, "Invalid numberOfGridInReference within distribution!")
          ENDIF ! IF (subGridId_curr == ...)

          CALL domainData%container%distributeField_grib( &
            & jg                 = jg,                                           & ! in
            & gridSize           = gridSize_curr,                                & ! in
            & mpi_work_id_sender = jproc,                                        & ! in
            & variableName       = variableName_curr(1:variableNameLength_curr), & ! in
            & level              = level_curr,                                   & ! in
            & tileId             = tileId_curr,                                  & ! in
            & isUniform          = isUniform_curr,                               & ! in
            & uniformValue       = uniformValue_curr,                            & ! in
            & field              = field(1:nelems_global)                        ) ! in

          !-------------------------------
          ! File inventory and statistics
          !-------------------------------

          ! We assume that the following holds true at this point(!):
          ! - record_flag_status == RECORD_FLAG_STATUS_GOTDATA
          ! - found_match        == .TRUE.
          IF (inventory .AND. process_is_next_to_distribute) THEN

            IF (timing) CALL timer_start(timer_file_inventory)

            ! The following initialization of the inventory element with field statistics
            ! is the second and last part of its initialization got started in InputRequestList_isRecordValid_grib below.
            CALL inventory_element%init(min                  = ecc_min,                      & ! in
              &                         max                  = ecc_max,                      & ! in
              &                         mean                 = ecc_avg,                      & ! in
              &                         tileId               = tileId,                       & ! in
              &                         isUniform            = isUniform,                    & ! in
              &                         missingValuesPresent = (ecc_numberOfMissing > 0_i8), & ! in
              &                         successful           = successful_local              ) ! out
            IF (.NOT. successful_local) CALL finish(routine, "Initialization of statistics of inventory element failed!")

            IF (timing) CALL timer_stop(timer_file_inventory)

          ENDIF ! IF (inventory .AND. process_is_next_to_distribute)

        END DO DISTRIBUTION_LOOP

        ! All PEs have to run through the file-inventory update
        IF (inventory) THEN

          IF (timing) CALL timer_start(timer_file_inventory)

          CALL file_inventory%add(inventory_element = inventory_element,   & ! inout
            &                     my_rank           = my_mpi_work_id,      & ! in
            &                     root_rank         = process_mpi_root_id, & ! in
            &                     comm              = p_comm_work,         & ! in
            &                     comm_size         = num_work_procs,      & ! in
            &                     reset             = .TRUE.,              & ! in
            &                     successful        = successful_local     ) ! out
          IF (.NOT. successful_local) CALL finish(routine, "Integration of inventory elements into file inventory failed!")

          IF (timing) CALL timer_stop(timer_file_inventory)

        ENDIF ! IF (inventory)

        IF (timing) CALL timer_stop(timer_data_distribution)

        IF (ecc_eof) THEN

          ! All-to-all GRIB message distribution is over.
          ! If we received "end-of-file" during this distribution cycle,
          ! this is the signal to exit the file processing loop
          EXIT FILE_PROCESSING_LOOP

        ELSEIF (.NOT. process_is_workroot) THEN

          IF (timing) CALL timer_start(timer_raw_data_distribution)

          ! The Non-Workroot PEs have to launch
          ! the next MPI_IRECV-call for the next loop cycle.
          ! (Note: Actually, no value of the following 4 variables
          ! does change for a Work PE. Nevertheless, we "re-set" them here,
          ! in order to reduce the error potential.)
          idx_status_requests_curr = 1
          recv_from_work_proc      = process_mpi_root_id
          recv_tag                 = my_mpi_work_id
          recv_count               = ecc_max_record_length + ABS(lbound_for_record_flags) + 1

          CALL p_irecv(t_buffer = ecc_record(:,idx_status_requests_curr), & ! inout
            &          p_source = recv_from_work_proc,                    & ! in
            &          p_tag    = recv_tag,                               & ! in
            &          p_count  = recv_count,                             & ! in
            &          comm     = p_comm_work,                            & ! in
            &          request  = status_request(idx_status_requests_curr)) ! optout

          IF (timing) CALL timer_stop(timer_raw_data_distribution)

        ENDIF ! IF (ecc_eof)

      END DO FILE_PROCESSING_LOOP

      IF (process_is_workroot) THEN

        IF (timing) CALL timer_start(timer_file_reading)

        ! Workroot PE has to close the GRIB file
        CALL ecc_close_file(ecc_ifile=ecc_ifile)

        IF (timing) CALL timer_stop(timer_file_reading)

        ! Print file inventory and statistics:
        ! (Note that we do not time this under timer_file_inventory,
        ! as it may perhaps result in problems with the timer hierarchy.)
        IF (inventory) THEN
          CALL file_inventory%print(grib_file_path         = grib_file_path,      & ! in
            &                       lIsFg                  = lIsFg,               & ! in
            &                       jg                     = jg,                  & ! in
            &                       legend                 = .TRUE.,              & ! in
            &                       clear                  = .TRUE.,              & ! in
            &                       successful             = successful_local     ) ! out
          IF (.NOT. successful_local) CALL finish(routine, "Printing the file inventory failed!")
        ENDIF

      ENDIF ! IF (process_is_workroot)

      ! Clean-up:

      IF (ALLOCATED(ecc_record)) THEN
        DEALLOCATE(ecc_record, STAT=status)
        IF (status /= SUCCESS) CALL finish(routine, "Deallocation of ecc_record failed")
      ENDIF

      IF (ALLOCATED(status_request)) THEN
        DEALLOCATE(status_request, STAT=status)
        IF (status /= SUCCESS) CALL finish(routine, "Deallocation of status_request failed")
      ENDIF

      IF (ALLOCATED(field)) THEN
        DEALLOCATE(field, STAT=status)
        IF (status /= SUCCESS) CALL finish(routine, "Deallocation of field failed")
      ENDIF

    END SUBROUTINE InputRequestList_readFile_grib

    !>
    !! @brief Check if metadate of GRIB message are consistent with requests
    !!
    SUBROUTINE InputRequestList_isRecordValid_grib(me, jg, ecc_msgid, lIsFg, verbose, inventory, hgrid_uuid, vgrid_uuid,      &
      &                                            ana_incr_list, verify_hgrid_uuid, verify_vgrid_uuid, verify_ana_incr_list, &
      &                                            message_prefix, message_prefix_length, level, tileId, variableName,        &
      &                                            variableNameLength, subGridId, gridSize, inventory_element, found_match)

      !-----------
      ! Arguments
      !-----------

      CLASS(t_InputRequestList), INTENT(INOUT) :: me

      !> Value of p_patch%id
      INTEGER,                   INTENT(IN)    :: jg

      !> GRIB message handle for interaction with ecCodes
      INTEGER,                   INTENT(IN)    :: ecc_msgid

      !> Flag to indicate first guess
      LOGICAL,                   INTENT(IN)    :: lIsFg

      !> Flag to indicate verbose messaging
      LOGICAL,                   INTENT(IN)    :: verbose

      !> Take a file inventory?
      LOGICAL,                   INTENT(IN)    :: inventory

      !> UUID of horizontal grid of ICON
      TYPE(t_uuid),              INTENT(IN)    :: hgrid_uuid

      !> UUID of vertical grid of ICON
      TYPE(t_uuid),              INTENT(IN)    :: vgrid_uuid

      !> List of input fields (with internal variable names)
      !> which may be provided as analysis increments
      CHARACTER(LEN=*),          INTENT(IN)    :: ana_incr_list(:)

      !> Flag to indicate if UUID of horizontal grid has to be verified
      LOGICAL,                   INTENT(IN)    :: verify_hgrid_uuid

      !> Flag to indicate if UUID of vertical grid has to be verified
      LOGICAL,                   INTENT(IN)    :: verify_vgrid_uuid

      !> Flag to indicate whether to verify if metadata of ana_incr_list-fields
      !> conform to analysis increments
      LOGICAL,                   INTENT(IN)    :: verify_ana_incr_list

      !> Prefix for messages
      CHARACTER(LEN=*),          INTENT(INOUT) :: message_prefix

      !> Length of message prefix
      INTEGER,                   INTENT(INOUT) :: message_prefix_length

      !> Level value
      REAL(dp),                  INTENT(OUT)   :: level

      !> ICON-internal tile identifier
      INTEGER,                   INTENT(OUT)   :: tileId

      !> Value of ecCodes concept key 'shortName'
      CHARACTER(LEN=*),          INTENT(OUT)   :: variableName

      !> Length of variableName
      INTEGER,                   INTENT(OUT)   :: variableNameLength

      !> Identifier of grid-element type (cells, edges(, vertices))
      INTEGER,                   INTENT(OUT)   :: subGridId

      !> Number of (global) grid points
      INTEGER(KIND=i8),          INTENT(OUT)   :: gridSize

      !> File inventory element
      TYPE(t_FileInventoryElement), INTENT(OUT) :: inventory_element

      !> Flag to indicate if check result is positive
      LOGICAL,                   INTENT(OUT)   :: found_match

      !-----------------
      ! Local variables
      !-----------------

      !> Metadata cache pointer(???)
      TYPE(t_MetadataCache), POINTER :: metadata

      !> (???)
      TYPE(t_ListEntry), POINTER :: listEntry

      !> (???)
      TYPE(t_DomainData), POINTER :: domainData

      !> Tile handling
      TYPE(t_tile_att), POINTER :: tile_att
      TYPE(t_tileinfo_icon)     :: tileinfo_icon
      TYPE(t_tileinfo_grb2)     :: tileinfo_grb2

      !> Status identifier
      INTEGER :: status

      !> GRIB triple of variable
      INTEGER :: ecc_discipline, ecc_parameterCategory, ecc_parameterNumber

      !> Meaning of reference date and time
      INTEGER :: ecc_significanceOfReferenceTime

      !> Lengths of datetime strings
      INTEGER :: ecc_dataDateTime_length, ecc_validityDateTime_length

      !> First fixed surface variables
      INTEGER  :: ecc_typeOfFirstFixedSurface
      REAL(wp) :: ecc_valueOfFirstFixedSurface
      LOGICAL  :: ecc_firstFixedSurface_is_missing

      !> Second fixed surface variables
      INTEGER  :: ecc_typeOfSecondFixedSurface
      REAL(wp) :: ecc_valueOfSecondFixedSurface
      LOGICAL  :: ecc_secondFixedSurface_is_missing

      !> GRIB key holding the number of horizontal grid points
      INTEGER(KIND=i8) :: ecc_numberOfDataPoints

      !> Definition number/identifier of type of horizontal grid
      INTEGER :: ecc_gridDefinitionTemplateNumber

      !> Value of GRIB key: numberOfGridUsed
      INTEGER :: ecc_numberOfGridUsed

      !> Value of GRIB key: numberOfGridInReference
      INTEGER :: ecc_numberOfGridInReference

      !> Value of GRIB key: uuidOfHGrid
      ! This value is a 128-bit UUID.
      ! ecCodes provides this value in different types/formats.
      ! Here, we store the 16 bytes of the UUID
      ! in a character array of size 16.
      CHARACTER(LEN=1) :: ecc_uuidOfHGrid(16)

      !> Switch to indicate general vertical height coordinates
      LOGICAL :: ecc_genVertHeightCoords

      !> Value of GRIB key: NV
      INTEGER :: ecc_NV

      !> Value of GRIG key: nlev
      INTEGER :: ecc_nlev

      !> Value of GRIB key: uuidOfVGrid
      ! This value is a 128-bit UUID.
      ! ecCodes provides this value in different types/formats.
      ! Here, we store the 16 bytes of the UUID
      ! in a character array of size 16.
      CHARACTER(LEN=1) :: ecc_uuidOfVGrid(16)

      !> uuidOfVGrid in internal UUID format
      TYPE(t_uuid) :: uuidOfVGrid

      !> uuidOfHGrid in internal UUID format
      TYPE(t_uuid) :: uuidOfHGrid

      !> Definition number/identifier of type of product
      INTEGER :: ecc_productDefinitionTemplateNumber

      !> Values of tile GRIB keys:
      INTEGER :: ecc_tileClassification
      INTEGER :: ecc_totalNumberOfTileAttributePairs
      INTEGER :: ecc_numberOfUsedSpatialTiles
      INTEGER :: ecc_tileIndex
      INTEGER :: ecc_numberOfUsedTileAttributes
      INTEGER :: ecc_attributeOfTile

      !> Values of centre GRIB keys:
      INTEGER :: ecc_centre
      INTEGER :: ecc_subCentre

      !> Values of generating-process-related GRIB keys:
      INTEGER :: ecc_typeOfProcessedData
      INTEGER :: ecc_typeOfGeneratingProcess
      INTEGER :: ecc_backgroundProcess
      INTEGER :: ecc_generatingProcessIdentifier

      !> Value of local GRIB key: localNumberOfExperiment
      INTEGER :: ecc_localNumberOfExperiment

      !> Reference and validity date and time in format:'YYYY-MM-DDThh:mm:ss'
      CHARACTER(LEN=max_datetime_str_len) :: ecc_dataDateTime, ecc_validityDateTime

      !> Flag to indicate successful ecCodes inquiries
      LOGICAL :: ecc_successful

      !> Loop index
      INTEGER :: j

      !> Instances of mtime types
      TYPE(datetime), POINTER :: tempTime, iniTime, startTime

      !> String corresponding to mtime instance
      CHARACTER(LEN=max_datetime_str_len) :: debugDatetimeString

      !> Strings for UUIDs
      CHARACTER(LEN=uuid_string_length) :: expected_uuid, found_uuid

      !> Status flag
      LOGICAL :: successful_local

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":InputRequestList_isRecordValid_grib"

      !----------------------------

      !
      ! Notes:
      !
      ! - The notes at the beginning of subroutine InputRequestList_readFile_grib apply here, too.
      !
      ! - The output arguments of the ecCodes wrapper routines
      !   will be initialized with some value in any case
      !   (no matter if they return "success" or not).
      !
      ! - ICON treats missing dictionary entries in a special way (see, e.g., InputRequestList_translateNames):
      !   If one of the variables among those, which are necessary for the current model setup,
      !   has no entry in the dictionary, the internal dictionary representation is expanded by a new entry,
      !   where the shortName value is set equal to the ICON-internal variable name.
      !

      ! Initialize intent-out arguments
      level              = REAL(MISSING_VALUE, KIND=dp)
      tileId             = -999
      variableName       = ' '
      variableNameLength = 0
      subGridId          = -999
      gridSize           = 0_i8
      found_match        = .TRUE.

      ! Create/get pointer of/to metadata cache(???)
      metadata => MetadataCache_create()

      !-----------------------------------------
      ! Decoding of GRIB metadata using ecCodes
      !-----------------------------------------

      ! Get values of product-related GRIB keys and ecCodes concept keys
      CALL ecc_get_info_on_product( &
        & ecc_msgid                           = ecc_msgid,                           & ! in
        & ecc_productDefinitionTemplateNumber = ecc_productDefinitionTemplateNumber, & ! out
        & ecc_discipline                      = ecc_discipline,                      & ! out
        & ecc_parameterCategory               = ecc_parameterCategory,               & ! out
        & ecc_parameterNumber                 = ecc_parameterNumber,                 & ! out
        & ecc_shortName                       = variableName,                        & ! out
        & ecc_shortName_length                = variableNameLength,                  & ! out
        & ecc_successful                      = ecc_successful                       ) ! out

      ! If ecc_successful is .false., this means that ecCodes
      ! was unable to find a shortName value, whose definition matches
      ! the metadata of the GRIB message (GRIB record).
      IF (.NOT. ecc_successful) THEN
        IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length)//": Unable to get values of product keys")
        found_match = .FALSE.
      ENDIF

      ! Store triple in metadata cache
      metadata%param%discipline = INT(ecc_discipline, KIND=C_INT)
      metadata%param%category   = INT(ecc_parameterCategory, KIND=C_INT)
      metadata%param%number     = INT(ecc_parameterNumber, KIND=C_INT)

      ! Expand message prefix by variableName
      IF (verbose) THEN
        message_prefix        = message_prefix(1:message_prefix_length)//"[shortName: " &
          &                   //variableName(1:variableNameLength)//"]"
        message_prefix_length = message_prefix_length + variableNameLength + 13
      ENDIF

      ! Get date and time in format 'YYYY-MM-DDThh:mm:ss':
      ! - Reference date and time => ecc_dataDateTime
      ! - Validity date and time  => ecc_validityDateTime (in general equal to reference date and time + forecast time)
      CALL ecc_get_info_on_datetime( &
        & ecc_msgid                       = ecc_msgid,                       & ! in
        & ecc_significanceOfReferenceTime = ecc_significanceOfReferenceTime, & ! out
        & ecc_dataDateTime                = ecc_dataDateTime,                & ! out
        & ecc_validityDateTime            = ecc_validityDateTime,            & ! out
        & ecc_dataDateTime_length         = ecc_dataDateTime_length,         & ! out
        & ecc_validityDateTime_length     = ecc_validityDateTime_length,     & ! out
        & ecc_successful                  = ecc_successful                   ) ! out

      IF (.NOT. ecc_successful) THEN
        IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length) &
          & //": Unable to get reference and validity date and time")
        found_match = .FALSE.
      ENDIF

      ! Store dates and times in metadata cache
      ! ('rtime' and 'vtime' seem to be some kind of pointer. Therefore, we assume
      ! that memory should be allocated before we can store something in them???)
      ALLOCATE(metadata%rtime(ecc_dataDateTime_length), metadata%vtime(ecc_validityDateTime_length), STAT=status)
      IF (status /= SUCCESS) CALL finish(routine, message_prefix(1:message_prefix_length) &
        & //": Allocation of rtime and vtime failed")

      ! Nothing is said about the format of 'rtime' and 'vtime'.
      ! We just assume it is: 'YYYY-MM-DDThh:mm:ss'(???)
      DO j = 1, ecc_dataDateTime_length
        metadata%rtime(j) = ecc_dataDateTime(j:j)
      END DO
      DO j = 1, ecc_validityDateTime_length
        metadata%vtime(j) = ecc_validityDateTime(j:j)
      END DO

      ! Check if date and time in GRIB metadate fit model date and time
      ! (lconsistency_checks is a switch of namelist initicon_nml)
      IF (lconsistency_checks) THEN

        tempTime => newDatetime(ecc_validityDateTime(1:ecc_validityDateTime_length))

        ALLOCATE(iniTime, STAT=status)
        IF (status /= SUCCESS) CALL finish(routine, message_prefix(1:message_prefix_length) &
          & //": Allocation of iniTime failed")

        ! Get model start data and time
        iniTime = time_config%tc_startdate

        IF (lIsFg) THEN

          ! First guess:

          ALLOCATE(startTime, STAT=status)
          IF (status /= SUCCESS) CALL finish(routine, message_prefix(1:message_prefix_length) &
            & //": Allocation of startTime failed")

          ! IAU time shift(???)
          startTime = iniTime + time_config%timeshift%mtime_shift

          IF (.NOT. (tempTime == startTime)) THEN
            IF (verbose) THEN
              CALL datetimeToString(startTime, debugDatetimeString)
              message_text = message_prefix(1:message_prefix_length)//": validity date and time of first-guess field (" &
                &          //ecc_validityDateTime(1:ecc_validityDateTime_length)//") does not match model start time (" &
                &          //TRIM(debugDatetimeString)//")"
              CALL warning(routine, message_text)
            ENDIF
            found_match = .FALSE.
          ENDIF

          DEALLOCATE(startTime, STAT=status)
          IF (status /= SUCCESS) CALL finish(routine, message_prefix(1:message_prefix_length) &
            & //": Deallocation of startTime failed")

        ELSE

          ! Analysis (increments):

          IF (.NOT. (tempTime == iniTime)) THEN
            IF (verbose) THEN
              CALL datetimeToString(iniTime, debugDatetimeString)
              message_text = message_prefix(1:message_prefix_length)//": validity date and time of analysis field ("             &
                &          //ecc_validityDateTime(1:ecc_validityDateTime_length)//") does not match model initialization time (" &
                &          //TRIM(debugDatetimeString)//")"
              CALL warning(routine, message_text)
            ENDIF
            found_match = .FALSE.
          ENDIF

        ENDIF ! IF (lIsFg)

        DEALLOCATE(iniTime, STAT=status)
        IF (status /= SUCCESS) CALL finish(routine, message_prefix(1:message_prefix_length) &
          & //": Deallocation of iniTime failed")

        CALL deallocateDatetime(tempTime)

      ENDIF ! IF (lconsistency_checks)

      ! Evaluate vertical grid and level/layer information in GRIB metadate:
      !
      !  |- typeOfFirstFixedSurface in [0, 254]        -|
      ! -|                                              | -> Level
      !  |- typeOfSecondFixedSurface = 255 ('missing') -|
      !
      !  |- typeOfFirstFixedSurface in [0, 254]        -|
      ! -|                                              | -> Layer
      !  |- typeOfSecondFixedSurface in [0, 254]       -|
      !
      ! Notes:
      !
      ! - Currently, only the type of the first fixed surface is effectively evaluated.
      !
      ! - We do not use the ecCodes concept key 'typeOfLevel' for the following reasons:
      !   - The value of this key is of type string.
      !     Unfortunately, any kind of string operation are not welcome
      !     in operational NWP production at DWD.
      !   - As an ecCodes concept key, the situation of 'typeOfLevel'
      !     is the same as, for instance, the sitation of the concept key 'shortName'.
      !     The definition of its range of values may change from center to center.
      !     But in contrast to the 'shortName', there is no dictionary for the 'typeOfLevel'
      !     (mapping its values to CDI zaxis, for instance).
      !
      CALL ecc_get_info_on_vertical_grid( &
        & ecc_msgid                         = ecc_msgid,                         & ! in
        & ecc_genVertHeightCoords           = ecc_genVertHeightCoords,           & ! out
        & ecc_NV                            = ecc_NV,                            & ! out
        & ecc_nlev                          = ecc_nlev,                          & ! out
        & ecc_uuidOfVGrid                   = ecc_uuidOfVGrid(:),                & ! out
        & ecc_typeOfFirstFixedSurface       = ecc_typeOfFirstFixedSurface,       & ! out
        & ecc_valueOfFirstFixedSurface      = ecc_valueOfFirstFixedSurface,      & ! out
        & ecc_firstFixedSurface_is_missing  = ecc_firstFixedSurface_is_missing,  & ! out
        & ecc_typeOfSecondFixedSurface      = ecc_typeOfSecondFixedSurface,      & ! out
        & ecc_valueOfSecondFixedSurface     = ecc_valueOfSecondFixedSurface,     & ! out
        & ecc_secondFixedSurface_is_missing = ecc_secondFixedSurface_is_missing, & ! out
        & ecc_successful                    = ecc_successful                     ) ! out

      IF (.NOT. ecc_successful) THEN
        IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length) &
          & //": Unable to get values of vertical-grid-related keys")
        found_match = .FALSE.
      ELSEIF (ecc_firstFixedSurface_is_missing) THEN
        IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length) &
          & //": Type of first fixed surface is missing")
        found_match = .FALSE.
      ELSE

        ! Special treatment for the level:
        ! (This is an ecCodes transfer of the corresponding code in InputRequestList_isRecordValid above,
        ! by trying some kind of "reverse-engineering" of CDI ZAXIS, cdiIterator_inqLevelType and cdiIterator_inqLevel.)
        SELECT CASE(ecc_typeOfFirstFixedSurface)
        CASE(1, 100, 102, 103, 106, 114, 150, 160)
          ! WMO GRIB2 code table 4.5 - Fixed surface types and units:
          !-----------------------------------------------------------
          ! 1   "Ground or water surface (-)"
          ! 100 "Isobaric surface (Pa)"
          ! 102 "Specific altitude above mean sea level (m)"
          ! 103 "Specified height level above ground (m)"
          ! 106 "Depth below land surface (m)"
          ! 114 "Snow level (Numeric)"
          ! 150 "Generalized vertical height coordinate (-)"
          ! 160 "Depth below sea level (m)"
          !-----------------------------------------------------------
          ! Parametric level type, e.g., 2m or 10m height above ground level
          ! (Actually, "ground or water surface" is a non-parametric level type,
          ! but who knows why it was originally added here...)
          level = REAL(ecc_valueOfFirstFixedSurface, KIND=dp)
          ! While values are stored in SI units in GRIB,
          ! ICON seems to measure "depth below land surface" in mm (???)
          ! (See icon/src/io/shared/mo_name_list_output_zaxes.f90: setup_ml_axes_atmo and
          ! icon/externals/cdi/src/cdilib.c: logicalLevelValue2).
          ! In addition, CDI seems to multiply isobaric levels by the factor 1000, too(???)
          IF (ANY([100, 106] == ecc_typeOfFirstFixedSurface)) level = level * 1000.0_dp
        CASE(2, 3, 4, 8, 9, 10, 101, 162, 163, 164, 165, 166, 173, 175)
          ! WMO GRIB2 code table 4.5 - Fixed surface types and units:
          !-----------------------------------------------------------
          ! 2   "Cloud base level (-)"
          ! 3   "Level of cloud tops (-)"
          ! 4   "Level of 0 degC isotherm (-)"
          ! 8   "Nominal top of the atmosphere (-)"
          ! 9   "Sea bottom (-)"
          ! 10  "Entire atmosphere (-)"
          ! 101 "Mean sea level (-)"
          ! 162 "Lake or river bottom (-)"
          ! 163 "Bottom of sediment layer (-)"
          ! 164 "Bottom of thermally active sediment layer (-)"
          ! 165 "Bottom of sediment layer penetrated by thermal wave (-)"
          ! 166 "Mixing layer (-)"
          ! 173 "Top of snow over sea ice on sea, lake or river"
          ! 175 "Top surface of ice, under snow cover, on sea, lake or river"
          !-----------------------------------------------------------
          ! Non-parametric level types, e.g. mean sea level
          level = -REAL(ecc_typeOfFirstFixedSurface, KIND=dp)
        CASE DEFAULT
          IF (verbose) THEN
            WRITE(message_text,'(A,I0,A)')                                   &
              & ": typeOfFirstFixedSurface = ", ecc_typeOfFirstFixedSurface, &
              & " is not supported (see WMO GRIB2 code tabel 4.5 'Fixed surface types and units' for its meaning)"
            CALL warning(routine, message_prefix(1:message_prefix_length)//message_text)
          ENDIF
          found_match = .FALSE.
        END SELECT ! SELECT CASE(ecc_typeOfFirstFixedSurface)

      ENDIF ! IF (.NOT. ecc_successful) ...

      ! Important difference to the implementation in InputRequestList_isRecordValid:
      ! metadata%levelType does no longer contain the CDI level type,
      ! but the ecCodes/GRIB level type!
      metadata%levelType = ecc_typeOfFirstFixedSurface

      ! Verify UUIDs of vertical grid
      IF (verify_vgrid_uuid .AND. ecc_genVertHeightCoords) THEN
        ! Transform character array with grid uuid into internal representation
        ! of derived type 't_uuid'
        CALL char2uuid(string=ecc_uuidOfVGrid(:), uuid=uuidOfVGrid)
        IF (.NOT. (uuidOfVGrid == vgrid_uuid)) THEN
          ! Get UUIDs as human-readable strings
          CALL uuid_unparse(vgrid_uuid, expected_uuid)
          CALL uuid_unparse(uuidOfVGrid, found_uuid)
          message_text = message_prefix(1:message_prefix_length)//": vgrid_uuid = "//TRIM(expected_uuid) &
            &          //" /= uuidOfVGrid = "//TRIM(found_uuid)
          IF (verbose) CALL warning(routine, message_text)
          found_match = .FALSE.
        ENDIF ! IF (.NOT. (uuidOfVGrid == vgrid_uuid))
      ENDIF ! IF (verify_vgrid_uuid .AND. ecc_genVertHeightCoords)

      ! Evaluate the horizontal grid:
      !
      ! Notes:
      !
      ! - Exclusively supported horizontal grid: 101 "General unstructured grid"
      !

      CALL ecc_get_info_on_horizontal_grid( &
        & ecc_msgid                        = ecc_msgid,                        & ! in
        & ecc_gridDefinitionTemplateNumber = ecc_gridDefinitionTemplateNumber, & ! out
        & ecc_numberOfDataPoints           = ecc_numberOfDataPoints,           & ! out
        & ecc_numberOfGridUsed             = ecc_numberOfGridUsed,             & ! out
        & ecc_numberOfGridInReference      = ecc_numberOfGridInReference,      & ! out
        & ecc_uuidOfHGrid                  = ecc_uuidOfHGrid(:),               & ! out
        & ecc_successful                   = ecc_successful                    ) ! out

      IF (.NOT. ecc_successful) THEN
        IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length) &
          & //": Unable to get info on horizontal grid")
        found_match = .FALSE.
      ELSEIF (ecc_gridDefinitionTemplateNumber /= 101) THEN
        IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length) &
          & //": Unsupported type of horizontal grid")
        found_match = .FALSE.
      ENDIF

      ! Identifier of sub-grid type: cells or edges
      subGridId = ecc_numberOfGridInReference

      ! Number of data points = number of grid (element) points
      gridSize = ecc_numberOfDataPoints

      ! Transform character array with grid uuid into internal representation
      ! of derived type 't_uuid'
      CALL char2uuid(string=ecc_uuidOfHGrid(:), uuid=metadata%gridUuid)

      metadata%gridNumber   = ecc_numberOfGridUsed
      metadata%gridPosition = ecc_numberOfGridInReference

      ! Verify UUIDs of horizontal grid
      IF (verify_hgrid_uuid) THEN
        IF (.NOT. (metadata%gridUuid == hgrid_uuid)) THEN
          ! Get UUIDs as human-readable strings
          CALL uuid_unparse(hgrid_uuid, expected_uuid)
          CALL uuid_unparse(metadata%gridUuid, found_uuid)
          message_text = message_prefix(1:message_prefix_length)//": hgrid_uuid = "//TRIM(expected_uuid) &
            &          //" /= uuidOfHGrid = "//TRIM(found_uuid)
          IF (verbose) CALL warning(routine, message_text)
          found_match = .FALSE.
        ENDIF ! IF (.NOT. (metadata%gridUuid == hgrid_uuid))
      ENDIF ! IF (verify_hgrid_uuid)

      ! Evaluate tile information:

      ! Determine argument tileId
      SELECT CASE(ecc_productDefinitionTemplateNumber)
      CASE(55, 59, 40455, 40456)

        ! Tile-based field:
        ! - 55, 59:...... official WMO tile templates
        ! - 40455, 40456: local DWD tile templates
        CALL ecc_get_info_on_tiles( &
          & ecc_msgid                           = ecc_msgid,                           & ! in
          & ecc_productDefinitionTemplateNumber = ecc_productDefinitionTemplateNumber, & ! in
          & ecc_tileClassification              = ecc_tileClassification,              & ! out
          & ecc_totalNumberOfTileAttributePairs = ecc_totalNumberOfTileAttributePairs, & ! out
          & ecc_numberOfUsedSpatialTiles        = ecc_numberOfUsedSpatialTiles,        & ! out
          & ecc_tileIndex                       = ecc_tileIndex,                       & ! out
          & ecc_numberOfUsedTileAttributes      = ecc_numberOfUsedTileAttributes,      & ! out
          & ecc_attributeOfTile                 = ecc_attributeOfTile,                 & ! out
          & ecc_successful                      = ecc_successful                       ) ! out

        IF (.NOT. ecc_successful) THEN
          IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length) &
            & //": Unable to get info on tiles")
          found_match = .FALSE.
        ENDIF

        tileinfo_grb2%idx = ecc_tileIndex
        tileinfo_grb2%att = ecc_attributeOfTile

        tile_att => tile_list%getTileAtt(tileinfo_grb2)

        tileinfo_icon = tile_att%getTileinfo_icon()

        ! Finally, get ICON-internal tileId
        tileId = tileinfo_icon%idx

      CASE DEFAULT

        ! No tile-based field
        tileinfo_icon = trivial_tile_att%getTileinfo_icon()
        tileId        = tileinfo_icon%idx

      END SELECT

      ! Evaluate generating centre:

      CALL ecc_get_info_on_generating_centre( &
        & ecc_msgid      = ecc_msgid,     & ! in
        & ecc_centre     = ecc_centre,    & ! out
        & ecc_subCentre  = ecc_subCentre, & ! out
        & ecc_successful = ecc_successful ) ! out

      IF (.NOT. ecc_successful) THEN
        IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length) &
          & //": Unable to get info on generating centre")
        found_match = .FALSE.
      ENDIF

      ! Evaluate generating process:

      ! Initialize metadata about the generating environment and process
      metadata%runClass              = -1
      metadata%experimentId          = -1
      metadata%generatingProcessType = -1

      CALL ecc_get_info_on_generating_process( &
        & ecc_msgid                       = ecc_msgid,                       & ! in
        & ecc_typeOfProcessedData         = ecc_typeOfProcessedData,         & ! out
        & ecc_typeOfGeneratingProcess     = ecc_typeOfGeneratingProcess,     & ! out
        & ecc_backgroundProcess           = ecc_backgroundProcess,           & ! out
        & ecc_generatingProcessIdentifier = ecc_generatingProcessIdentifier, & ! out
        & ecc_successful                  = ecc_successful                   ) ! out

      IF (.NOT. ecc_successful) THEN
        IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length) &
          & //": Unable to get info on generating process")
        found_match = .FALSE.
      ENDIF

      metadata%runClass              = ecc_backgroundProcess
      metadata%generatingProcessType = ecc_typeOfGeneratingProcess

      ! Use of local GRIB key 'localNumberOfExperiment' is restricted
      ! to input data generated by DWD (centre = 78 / edzw)
      IF (ecc_centre == 78) THEN

        CALL ecc_get_local_info( &
          & ecc_msgid                   = ecc_msgid,                   & ! in
          & ecc_localNumberOfExperiment = ecc_localNumberOfExperiment, & ! out
          & ecc_successful              = ecc_successful               ) ! out

        IF (.NOT. ecc_successful) THEN
          IF (verbose) CALL warning(routine, message_prefix(1:message_prefix_length)//": Unable to get local info")
          found_match = .FALSE.
        ENDIF

        metadata%experimentId = ecc_localNumberOfExperiment

      ENDIF ! IF (ecc_centre == 78)

      !------------------------------------
      ! Internal matching of metadata(???)
      !------------------------------------

      ! Check whether metadata of this GRIB record is consistent
      ! with metadata we've already seen for this variable:

      ! Note that we actually do not know whether this is really necessary
      ! or makes any sense in this subroutine.
      ! It is inherited from InputRequestList_isRecordValid,
      ! which in turn was designed for the Workroot PE being the sole process
      ! that decodes GRIB messages.
      ! Here, each Work PE does it on its own.

      listEntry => me%findTranslatedName(variableName(1:variableNameLength))

      IF (.NOT. ASSOCIATED(listEntry)) THEN

        ! We are not interested in this variable:

        found_match = .FALSE.
        ! Destruct pointer of/to metadata cache(???)
        CALL metadata%destruct()
        DEALLOCATE(metadata, STAT=status)
        IF (status /= SUCCESS) CALL finish(routine, message_prefix(1:message_prefix_length) &
          & //": Deallocation of metadata failed")

      ELSEIF (found_match .AND. verify_ana_incr_list .AND. .NOT. lIsFg) THEN

        ! Check if metadata conform with analysis increments if they should (applies to init_mode = 5 (MODE_IAU) only):

        ! Important note:
        ! In case of sequential input-file processing, a list of analysis-increment fields
        ! is explicitly specified in atm_dyn_iconam/mo_initicon: read_dwdana (=> incrementsList).
        ! Here, in the parallel input-file processing, however, we try to do without such explicit lists,
        ! as they are a nightmare for maintenance.
        ! Instead, the variable group "ana_increment" is introduced, in order to mark analysis-increment fields
        ! as such at their add_var registration (following the guiding principle that all field-specific information
        ! is centralized in this one place).
        ! Unfortunately, this comes along with new problems for which ICON peculiarities bear the blame,
        ! as outlined by the following example:
        ! The number concentration of cloud droplets "qnc" is potentially among those fields,
        ! which may be provided as analysis increments in the ANA input file.
        ! qnc is required for the two-moment microphysics scheme (nwp_phy_nml: inwp_gscp = 4).
        ! The add_var registration of qnc is encapsulated by an IF-block in such a way,
        ! that it takes place only if the two-moment scheme is switched on.
        ! If another scheme is used, the info about qnc made in add_var is inaccessible,
        ! which includes its "ana_increment"-group membership.
        ! This in turn means, that qnc will not appear in "ana_incr_list".
        ! However, analysis increments of qnc may be contained in the ANA input file for whatever reason.
        ! For the following check this would mean that qnc could not be found inside ana_incr_list,
        ! so that its metadata would be checked for "typeOfGeneratingProcess = 0",
        ! which inevitably leads to a program abort.
        ! Now, to make a long story short:
        ! For the following check to work, we rely on the assumption that the above assignment:
        !
        !   listEntry => me%findTranslatedName(variableName(1:variableNameLength))
        !
        ! would leave "listEntry" unassociated in case of qnc, whereby the following check would be skipped.

        ! WMO GRIB2 code table 4.3 - Type of generating process:
        !--------------------------------------------------------
        ! 0   "Analysis"
        ! 20  "Analysis increment"
        ! 201 "Diff. analysis - first guess" (local DWD definition, which is equivalent to 20)
        !--------------------------------------------------------

        IF (one_of(TRIM(listEntry%iconVarName), ana_incr_list) > 0) THEN
          IF (ALL([20, 201] /= ecc_typeOfGeneratingProcess)) CALL finish(routine, message_prefix(1:message_prefix_length) &
            & //": Metadata do not conform with analysis increments")
        ELSE
          IF (ecc_typeOfGeneratingProcess /= 0) CALL finish(routine, message_prefix(1:message_prefix_length) &
            & //": Metadata do not conform with a full analysis field")
        ENDIF

      ENDIF ! IF (.NOT. ASSOCIATED(listEntry))

      !-----------------------------------
      ! Store metadata for file inventory
      !-----------------------------------

      ! After the final potential change of found_match, we can initialize the inventory
      ! (in terms of metadata that are independent of the level value(s) and the field values)
      IF (inventory .AND. found_match) THEN

        CALL char2uuid(string=ecc_uuidOfHGrid(:), uuid=uuidOfHGrid)

        CALL inventory_element%init(iconVarName           = TRIM(listEntry%iconVarName), & ! in
          &                         dataDateTime          = ecc_dataDateTime,            & ! in
          &                         validityDateTime      = ecc_validityDateTime,        & ! in
          &                         discipline            = ecc_discipline,              & ! in
          &                         parameterCategory     = ecc_parameterCategory,       & ! in
          &                         parameterNumber       = ecc_parameterNumber,         & ! in
          &                         levelType             = ecc_typeOfFirstFixedSurface, & ! in
          &                         gridNumber            = ecc_numberOfGridUsed,        & ! in
          &                         gridPosition          = ecc_numberOfGridInReference, & ! in
          &                         runClass              = ecc_backgroundProcess,       & ! in
          &                         experimentId          = ecc_localNumberOfExperiment, & ! in
          &                         generatingProcessType = ecc_typeOfGeneratingProcess, & ! in
          &                         gridUuid              = uuidOfHGrid,                 & ! in
          &                         successful            = successful_local             ) ! out

        IF (.NOT. successful_local) CALL finish(routine, message_prefix(1:message_prefix_length) &
          & //": Initialization of inventory element failed")

      ELSEIF (inventory) THEN

        ! In this case the GRIB record is not required and ignored
        CALL inventory_element%reject()

      ENDIF ! IF (found_match .AND. inventory) ...

      ! In order to allow for a proper inventory also in case of found_match == .FALSE.,
      ! we had to move the RETURN from the last but one if-block to here
      IF (.NOT. ASSOCIATED(listEntry)) RETURN

      ! Store metadata on variable in linked list(???):

      IF (found_match) THEN

        domainData => findDomainData(listEntry, jg)

        ! Does the examined GRIB message match one of the requests(???)
        ! Uncommented because it does not work for multiple reading and decoding PEs for unknown reasons.
!!!        IF (ASSOCIATED(domainData)) found_match = metadata%equalTo(domainData%metadata)

      ENDIF ! IF (found_match)

      ! Commit and clean-up:

      IF (found_match) THEN

        IF (ASSOCIATED(domainData)) THEN

          CALL metadata%destruct()
          DEALLOCATE(metadata, STAT=status)
          IF (status /= SUCCESS) CALL finish(routine, message_prefix(1:message_prefix_length) &
            & //": Deallocation of metadata failed")

        ELSE

          domainData => findDomainData(listEntry, jg, opt_lcreate=.TRUE.)

          ! There is no metadata cache, so this one is remembered(???)
          domainData%metadata => metadata

        ENDIF ! IF (ASSOCIATED(domainData))

      ELSE

        ! The record is not valid:

        variableName       = " "
        variableNameLength = 0

        CALL metadata%destruct()
        DEALLOCATE(metadata, STAT=status)
        IF (status /= SUCCESS) CALL finish(routine, message_prefix(1:message_prefix_length) &
          & //": Deallocation of metadata failed")

      ENDIF ! IF (found_match)

    END SUBROUTINE InputRequestList_isRecordValid_grib

    FUNCTION InputRequestList_getLevels(me, varName, jg, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug
        REAL(dp), POINTER :: resultVar(:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_getLevels"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch level data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar => NULL()
        IF(ASSOCIATED(domainData)) resultVar => domainData%container%getLevels()
    END FUNCTION InputRequestList_getLevels

    LOGICAL FUNCTION InputRequestList_fetch2d(me, varName, level, tile, jg, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), INTENT(IN) :: level
        INTEGER, INTENT(IN) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetch2d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER   :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) resultVar = domainData%container%fetch2d(level, tile, outData, opt_lDebug)
        IF(resultVar) THEN
          this_att => tile_list%getTileAtt(t_tileinfo_icon(tile))
          CALL inverse_post_op( &
            &   TRIM(varName//TRIM(this_att%getTileSuffix())), &
            &   outData)
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetch2d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetch2d

    LOGICAL FUNCTION InputRequestList_fetch3d(me, varName, tile, jg, outData, optLevelDimension, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetch3d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        TYPE(t_tile_att), POINTER   :: this_att  ! pointer to attribute
        LOGICAL :: debugInfo

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) resultVar = domainData%container%fetch3d(tile, outData, optLevelDimension, opt_lDebug)
        IF(resultVar .AND. varName /= 'smi' .AND. varName /= 'SMI') THEN   !SMI IS NOT IN the ICON variable lists, so we need to skip inverse postprocessing for it manually.
          this_att => tile_list%getTileAtt(t_tileinfo_icon(tile))
          CALL inverse_post_op( &
            &   TRIM(varName//TRIM(this_att%getTileSuffix())), &
            &   outData)
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetch3d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetch3d

    LOGICAL FUNCTION InputRequestList_fetchSurface(me, varName, tile, jg, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchSurface"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        REAL(dp), POINTER :: levels(:)
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) THEN
            levels => domainData%container%getLevels()
            SELECT CASE(SIZE(levels, 1))
                CASE(0)
                    resultVar = .FALSE.
                    IF(debugInfo) CALL message(routine, "no levels found")
                CASE(1)
                    resultVar = domainData%container%fetch2d(levels(1), tile, outData, opt_lDebug)
                    IF(debugInfo .AND. .NOT. resultVar) CALL message(routine, "InputContainer_fetch2d() returned an error")
                CASE DEFAULT
                    CALL finish(routine, "trying to read '"//varName//"' as a surface variable, but the file contains several &
                                         &levels of this variable")
            END SELECT
        END IF
        IF(resultVar) THEN
            this_att => tile_list%getTileAtt(t_tileinfo_icon(tile))
            CALL inverse_post_op( &
            &   TRIM(varName//TRIM(this_att%getTileSuffix())), &
            &   outData)
        END IF
    END FUNCTION InputRequestList_fetchSurface

    LOGICAL FUNCTION InputRequestList_fetchTiled2d(me, varName, level, jg, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), INTENT(IN) :: level
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiled2d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        INTEGER :: i
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) resultVar = domainData%container%fetchTiled2d(level, outData, opt_lDebug)
        IF(resultVar) THEN
            DO i = 1, SIZE(outData, 3)
                this_att => tile_list%getTileAtt(t_tileinfo_icon(i))
                CALL inverse_post_op(TRIM(varName//TRIM(this_att%getTileSuffix())), &
                  &                  outData(:,:,i))
            END DO
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetchTiled2d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetchTiled2d

    LOGICAL FUNCTION InputRequestList_fetchTiled3d(me, varName, jg, outData, optLevelDimension, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiled3d"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        INTEGER :: i
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) resultVar = domainData%container%fetchTiled3d(outData, optLevelDimension, opt_lDebug)
        IF(resultVar .AND. varName /= 'smi' .AND. varName /= 'SMI') THEN   !SMI IS NOT IN the ICON variable lists, so we need to skip inverse postprocessing for it manually.
            DO i = 1, SIZE(outData, 4)
                this_att => tile_list%getTileAtt(t_tileinfo_icon(i))
                CALL inverse_post_op(TRIM(varName//TRIM(this_att%getTileSuffix())), &
                  &                  outData(:,:,:,i))
            END DO
        ELSE IF(debugInfo) THEN
            CALL message(routine, "InputContainer_fetchTiled3d() returned an error")
        END IF
    END FUNCTION InputRequestList_fetchTiled3d

    LOGICAL FUNCTION InputRequestList_fetchTiledSurface(me, varName, jg, outData, opt_lDebug) RESULT(resultVar)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        LOGICAL, OPTIONAL, INTENT(IN) :: opt_lDebug

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchTiledSurface"
        TYPE(t_ListEntry), POINTER :: listEntry
        TYPE(t_DomainData), POINTER :: domainData
        REAL(dp), POINTER :: levels(:)
        INTEGER :: i
        LOGICAL :: debugInfo
        TYPE(t_tile_att), POINTER  :: this_att  ! pointer to attribute

        debugInfo = .FALSE.
        IF(PRESENT(opt_lDebug)) debugInfo = opt_lDebug

        listEntry => me%findIconName(varName, opt_lDebug)
        IF(.NOT. ASSOCIATED(listEntry)) THEN
            CALL finish(routine, 'attempt to fetch data for an input variable "'//varName//'" that has not been requested')
        END IF
        domainData => findDomainData(listEntry, jg)
        resultVar = ASSOCIATED(domainData)
        IF(resultVar) THEN
            levels => domainData%container%getLevels()
            SELECT CASE(SIZE(levels, 1))
                CASE(0)
                    resultVar = .FALSE.
                    IF(debugInfo) CALL message(routine, "no levels found")
                CASE(1)
                    resultVar = domainData%container%fetchTiled2d(levels(1), outData, opt_lDebug)
                    IF(debugInfo .AND. .NOT. resultVar) CALL message(routine, "InputContainer_fetch2d() returned an error")
                CASE DEFAULT
                    CALL finish(routine, "trying to read '"//varName//"' as a surface variable, but the file contains several &
                                         &levels of this variable")
            END SELECT
        END IF
        IF(resultVar) THEN
            DO i = 1, SIZE(outData, 3)
                this_att => tile_list%getTileAtt(t_tileinfo_icon(i))
                CALL inverse_post_op(TRIM(varName//TRIM(this_att%getTileSuffix())), &
                  &                  outData(:,:,i))
            END DO
        END IF
    END FUNCTION InputRequestList_fetchTiledSurface

    SUBROUTINE InputRequestList_fetchRequired2d(me, varName, level, tile, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), INTENT(in) :: level
        INTEGER, INTENT(in) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequired2d"

        IF(.NOT. me%fetch2d(varName, level, tile, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequired2d

    SUBROUTINE InputRequestList_fetchRequired3d(me, varName, tile, jg, outData, optLevelDimension)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(in) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequired3d"

        IF(.NOT. me%fetch3d(varName, tile, jg, outData, optLevelDimension)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequired3d

    SUBROUTINE InputRequestList_fetchRequiredSurface(me, varName, tile, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(in) :: tile, jg
        REAL(wp), INTENT(INOUT) :: outData(:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredSurface"

        IF(.NOT. me%fetchSurface(varName, tile, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredSurface

    SUBROUTINE InputRequestList_fetchRequiredTiled2d(me, varName, level, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        REAL(dp), INTENT(in) :: level
        INTEGER, INTENT(in) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiled2d"

        IF(.NOT. me%fetchTiled2d(varName, level, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiled2d

    SUBROUTINE InputRequestList_fetchRequiredTiled3d(me, varName, jg, outData, optLevelDimension)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:,:)
        INTEGER, OPTIONAL, INTENT(IN) :: optLevelDimension

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiled3d"

        IF(.NOT. me%fetchTiled3d(varName, jg, outData, optLevelDimension)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiled3d

    SUBROUTINE InputRequestList_fetchRequiredTiledSurface(me, varName, jg, outData)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: varName
        INTEGER, INTENT(IN) :: jg
        REAL(wp), INTENT(INOUT) :: outData(:,:,:)

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_fetchRequiredTiledSurface"

        IF(.NOT. me%fetchTiledSurface(varName, jg, outData)) THEN
            CALL finish(routine, 'data read for variable "'//varName//'" is incomplete')
        END IF
    END SUBROUTINE InputRequestList_fetchRequiredTiledSurface

    SUBROUTINE InputRequestList_checkRuntypeAndUuids(me, incrementVariables, gridUuids, lIsFg, lHardCheckUuids)
        CLASS(t_InputRequestList), INTENT(IN) :: me
        CHARACTER(*), INTENT(IN) :: incrementVariables(:)
        TYPE(t_uuid), INTENT(IN) :: gridUuids(:)    !< gridUuids(n_dom)
        LOGICAL, INTENT(IN) :: lIsFg, lHardCheckUuids

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_checkRuntypeAndUuids"
        INTEGER :: i, jg, expectedRuntype
        TYPE(t_ListEntry), POINTER :: curVar
        TYPE(t_DomainData), POINTER :: curDomain
        CHARACTER(:), POINTER :: varnameString
        CHARACTER(LEN = uuid_string_length) :: expectedUuid, foundUuid

        IF(.NOT.my_process_is_mpi_workroot()) CALL finish(routine, "assertion failed")
        IF(SIZE(gridUuids, 1) /= n_dom) CALL finish(routine, "assertion failed")
        DO jg = 1, n_dom
            DO i = 1, me%variableCount
                curVar => me%list(i)
                curDomain => findDomainData(curVar, jg)
                IF(.NOT.ASSOCIATED(curDomain)) CYCLE
                varnameString => curVar%iconVarName

                ! first check the TYPE of the generating process of the DATA
                IF(.NOT.lIsFg) THEN
                    IF(one_of(varnameString, incrementVariables) > 0) THEN
                        expectedRuntype = 201   ! analysis increment variables
                    ELSE
                        expectedRuntype = 0 ! analysis full variables
                    END IF
                    IF(expectedRuntype /= curDomain%metadata%generatingProcessType) THEN
                        CALL finish(routine, "detected wrong type of generating process on variable '"//varnameString//"', &
                                             &expected "//TRIM(int2string(expectedRuntype))//", &
                                             &found "//TRIM(int2string(curDomain%metadata%generatingProcessType)))
                    END IF
                END IF

                ! second check the UUID of the grid
                IF(.NOT.(gridUuids(jg) == curDomain%metadata%gridUuid)) THEN
                    CALL uuid_unparse(gridUuids(jg), expectedUuid)
                    CALL uuid_unparse(curDomain%metadata%gridUuid, foundUuid)
                    IF(lHardCheckUuids) THEN
                        CALL finish(routine, "detected wrong UUID of grid for variable '"//varnameString//"', &
                                             &expected "//expectedUuid//", &
                                             &found "//foundUuid)
                    ELSE
                        CALL message(routine, "warning: unexpected UUID of grid for variable '"//varnameString//"', &
                                             &expected "//expectedUuid//", &
                                             &found "//foundUuid)
                    END IF
                END IF
            END DO
        END DO
    END SUBROUTINE InputRequestList_checkRuntypeAndUuids

    SUBROUTINE InputRequestList_printInventory(me)
        CLASS(t_InputRequestList), INTENT(IN) :: me

        CHARACTER(*), PARAMETER :: routine = modname//":InputRequestList_printInventory"
        INTEGER :: i, jg, curRow, levelCount, tileCount
        LOGICAL :: lUntiledData
        TYPE(t_ListEntry), POINTER :: curVar
        TYPE(t_DomainData), POINTER :: curDomain
        CHARACTER(:), POINTER :: rtimeString, vtimeString
        TYPE(t_table) :: table
        CHARACTER(*), PARAMETER :: domainCol = "jg", &
                                 & variableCol = "variable", &
                                 & tripleCol = "triple", &
                                 & vtimeCol = "validity time", &
                                 & levelTypeCol = "levTyp", &
                                 & levelCountCol = "nlev", &
                                 & tileCountCol = "tileCnt", &
                                 & untiledCol = "untiled", &
                                 & runtypeCol = "runtype", &
                                 & vvmmCol = "vvmm", &
                                 & clasCol = "clas", &
                                 & expidCol = "expid", &
                                 & gridCol = "grid", &
                                 & rgridCol = "rgrid", &
                                 & minCol = "min", &
                                 & meanCol = "mean", &
                                 & maxCol = "max"
        CHARACTER(LEN = 3*3+2) :: parameterString
        TYPE(datetime), POINTER :: rtime, vtime
        TYPE(timedelta), POINTER :: forecastTime
        CHARACTER(len=max_timedelta_str_len) :: forecastTimeString

        CALL initialize_table(table)

        CALL add_table_column(table, domainCol)
        CALL add_table_column(table, variableCol)
        CALL add_table_column(table, tripleCol)
        CALL add_table_column(table, vtimeCol)
        CALL add_table_column(table, vvmmCol)
        CALL add_table_column(table, levelTypeCol)
        CALL add_table_column(table, levelCountCol)
        CALL add_table_column(table, tileCountCol)
        CALL add_table_column(table, untiledCol)
        CALL add_table_column(table, runtypeCol)
        CALL add_table_column(table, clasCol)
        CALL add_table_column(table, expidCol)
        CALL add_table_column(table, gridCol)
        CALL add_table_column(table, rgridCol)
        CALL add_table_column(table, minCol)
        CALL add_table_column(table, meanCol)
        CALL add_table_column(table, maxCol)

        IF(.NOT.my_process_is_mpi_workroot()) CALL finish(routine, "assertion failed")
        curRow = 1  !we can have zero to n_dom rows for each variable, so we can't USE the loop counter for the rows
        DO jg = 1, n_dom
            DO i = 1, me%variableCount
                curVar => me%list(i)
                curDomain => findDomainData(curVar, jg)
                IF(.NOT.ASSOCIATED(curDomain)) CYCLE
                CALL curDomain%container%getCounts(levelCount, tileCount, lUntiledData)

                !domain, NAME, AND triple columns
                CALL set_table_entry(table, curRow, domainCol, TRIM(int2string(curDomain%jg)))
                CALL set_table_entry(table, curRow, variableCol, curVar%iconVarName)
                WRITE(parameterString, '(3(I3,:,"."))') curDomain%metadata%param%discipline, curDomain%metadata%param%category, &
                &                                    curDomain%metadata%param%number
                CALL set_table_entry(table, curRow, tripleCol, parameterString)


                !date AND forecast time columns
                rtimeString => toCharacter(curDomain%metadata%rtime)
                vtimeString => toCharacter(curDomain%metadata%vtime)
                CALL set_table_entry(table, curRow, vtimeCol, vtimeString)

                rtime => newDatetime(rtimeString)
                vtime => newDatetime(vtimeString)

                forecastTime => newTimedelta("PT00H")  ! this 'initialization' IS necessary, IN order to correctly deal with timedelta=0.
                forecastTime = vtime - rtime
                CALL timedeltaToString(forecastTime, forecastTimeString)
                CALL set_table_entry(table, curRow, vvmmCol, TRIM(forecastTimeString))

                CALL deallocateDatetime(rtime)
                CALL deallocateDatetime(vtime)
                CALL deallocateTimedelta(forecastTime)
                DEALLOCATE(rtimeString)
                DEALLOCATE(vtimeString)


                !the simpler columns
                CALL set_table_entry(table, curRow, levelTypeCol, TRIM(int2string(curDomain%metadata%levelType)))
                CALL set_table_entry(table, curRow, levelCountCol, TRIM(int2string(levelCount)))
                IF(tileCount /= 0) CALL set_table_entry(table, curRow, tileCountCol, TRIM(int2string(tileCount)))
                IF(lUntiledData) THEN
                    CALL set_table_entry(table, curRow, untiledCol, "yes")
                ELSE
                    CALL set_table_entry(table, curRow, untiledCol, "no")
                END IF
                CALL set_table_entry(table, curRow, clasCol, TRIM(int2string(curDomain%metadata%runClass)))
                CALL set_table_entry(table, curRow, expidCol, TRIM(int2string(curDomain%metadata%experimentId)))
                IF(curDomain%metadata%generatingProcessType /= -1) THEN
                    CALL set_table_entry(table, curRow, runtypeCol, TRIM(int2string(curDomain%metadata%generatingProcessType)))
                END IF
                IF(curDomain%metadata%gridNumber /= -1) THEN
                    CALL set_table_entry(table, curRow, gridCol, TRIM(int2string(curDomain%metadata%gridNumber)))
                END IF
                IF(curDomain%metadata%gridPosition /= -1) THEN
                    CALL set_table_entry(table, curRow, rgridCol, TRIM(int2string(curDomain%metadata%gridPosition)))
                END IF
                CALL set_table_entry(table, curRow, minCol, TRIM(real2string(curDomain%statistics%MIN)))
                CALL set_table_entry(table, curRow, meanCol, TRIM(real2string(curDomain%statistics%mean)))
                CALL set_table_entry(table, curRow, MAXCol, TRIM(real2string(curDomain%statistics%MAX)))


                !next row
                curRow = curRow + 1
            END DO
        END DO

        CALL print_table(table, opt_delimiter = " | ")
        CALL finalize_table(table)

    END SUBROUTINE InputRequestList_printInventory


    SUBROUTINE InputRequestList_destruct(me)
        CLASS(t_InputRequestList), INTENT(INOUT) :: me
        INTEGER :: i
        TYPE(t_ListEntry), POINTER :: currentEntry
        TYPE(t_DomainData), POINTER :: domainData, domainDataTemp

        DO i = 1, me%variableCount
            currentEntry => me%list(i)
            domainData => currentEntry%domainData
            DO WHILE (ASSOCIATED(domainData))
                IF(ASSOCIATED(domainData%container)) THEN
                    CALL domainData%container%destruct()
                    DEALLOCATE(domainData%container)
                END IF
                IF(ASSOCIATED(domainData%metadata)) THEN
                    CALL domainData%metadata%destruct()
                    DEALLOCATE(domainData%metadata)
                END IF
                domainDataTemp => domainData%next
                DEALLOCATE(domainData)
                domainData => domainDataTemp
            END DO
        END DO
        DEALLOCATE(me%list)
    END SUBROUTINE InputRequestList_destruct

    FUNCTION MetadataCache_create() RESULT(resultVar)
        TYPE(t_MetadataCache), POINTER :: resultVar

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":MetadataCache_create"
        INTEGER :: error

        ALLOCATE(resultVar, STAT = error)
        IF(error /= success) CALL finish(routine, "memory allocation error")
        resultVar%vtime => NULL()
        resultVar%rtime => NULL()
    END FUNCTION MetadataCache_create

    LOGICAL FUNCTION MetadataCache_equalTo(me, other) RESULT(resultVar)
        CLASS(t_MetadataCache), INTENT(IN) :: me, other

        CHARACTER(LEN = *), PARAMETER :: routine = modname//":MetadataCache_create"

        INTEGER :: gridNumber, gridPosition, runClass, experimentId, generatingProcessType

        resultVar = .FALSE.

        !compare the time strings
        IF(.NOT.ASSOCIATED(me%rtime).OR..NOT.ASSOCIATED(other%rtime)) THEN
            CALL finish(routine, "internal error, please report this bug")
        END IF
        IF(SIZE(me%rtime) /= SIZE(other%rtime) .OR. ANY(me%rtime /= other%rtime)) THEN
            CALL message(routine, "inconsistent rtime detected")
            RETURN
        END IF

        IF(.NOT.ASSOCIATED(me%vtime).OR..NOT.ASSOCIATED(other%vtime)) THEN
            CALL finish(routine, "internal error, please report this bug")
        END IF
        IF(SIZE(me%vtime) /= SIZE(other%vtime) .OR. ANY(me%vtime /= other%vtime)) THEN
            CALL message(routine, "inconsistent vtime detected")
            RETURN
        END IF

        !compare the parameters
        IF(me%param%discipline /= other%param%discipline) THEN
            CALL message(routine, "inconsistent discipline detected")
            RETURN
        END IF
        IF(me%param%category /= other%param%category) THEN
            CALL message(routine, "inconsistent category detected")
            RETURN
        END IF
        IF(me%param%number /= other%param%number) THEN
            CALL message(routine, "inconsistent number detected")
            RETURN
        END IF

        !compare the other fields
        IF(me%levelType /= other%levelType) THEN
            CALL message(routine, "inconsistent level type detected")
            RETURN
        END IF

        IF(me%gridNumber /= other%gridNumber) THEN
            CALL message(routine, "inconsistent number of grids detected")
            RETURN
        END IF

        IF(me%gridPosition /= other%gridPosition) THEN
            CALL message(routine, "inconsistent grid index detected")
            RETURN
        END IF

        IF(me%runClass /= other%runClass) THEN
            CALL message(routine, "inconsistent run CLASS detected")
            RETURN
        END IF

        IF(me%experimentId /= other%experimentId) THEN
            CALL message(routine, "inconsistent experiment ID detected")
            RETURN
        END IF

        IF(me%generatingProcessType /= other%generatingProcessType) THEN
            CALL message(routine, "inconsistent type of generating process detected")
            RETURN
        END IF

        IF(.NOT.(me%gridUuid == other%gridUuid)) THEN
            CALL message(routine, "inconsistent UUID of grid detected")
            RETURN
        END IF

        resultVar = .TRUE.
    END FUNCTION MetadataCache_equalTo

    SUBROUTINE MetadataCache_destruct(me)
        CLASS(t_MetadataCache), INTENT(INOUT) :: me

        IF(ASSOCIATED(me%vtime)) DEALLOCATE(me%vtime)
        IF(ASSOCIATED(me%rtime)) DEALLOCATE(me%rtime)
    END SUBROUTINE MetadataCache_destruct

    !>
    !! @brief Initialize metadata of an instance of t_FileInventoryElement
    !!
    SUBROUTINE FileInventoryElement_init_metadata(me, iconVarName, dataDateTime, validityDateTime, discipline, &
      &                                           parameterCategory, parameterNumber, levelType, gridNumber,   &
      &                                           gridPosition, runClass, experimentId, generatingProcessType, &
      &                                           gridUuid, successful)

      !-----------
      ! Arguments
      !-----------

      !> Passed-object dummy argument
      CLASS(t_FileInventoryElement), INTENT(OUT) :: me

      !> ICON-internal variable name
      CHARACTER(LEN=*), INTENT(IN)  :: iconVarName

      !> Reference and validity date and time in format:'YYYY-MM-DDThh:mm:ss'
      CHARACTER(LEN=*), INTENT(IN)  :: dataDateTime, validityDateTime

      !> Parameter triple
      INTEGER,          INTENT(IN)  :: discipline, parameterCategory, parameterNumber

       !> GRIB level type
      INTEGER,          INTENT(IN)  :: levelType

      !> Consecutive number of horizontal grid (GRIB key: numberOfGridUsed)
      INTEGER,          INTENT(IN)  :: gridNumber

      !> Identifier of horizontal-grid element (cell, edge or vertex, GRIB key: numberOfGridInReference)
      INTEGER,          INTENT(IN)  :: gridPosition

      !> Identifier of background generating process (GRIB key: backgroundProcess)
      INTEGER,          INTENT(IN)  :: runClass

      !> DWD-specific experiment identifier (GRIB key: localNumberOfExperiment)
      INTEGER,          INTENT(IN)  :: experimentId

      !> Identifier of process used to create the product (GRIB key: typeOfGeneratingProcess)
      INTEGER,          INTENT(IN)  :: generatingProcessType

      !> Fingerprint of horizontal grid (GRIB key: uuidOfHGrid)
      TYPE(t_uuid),     INTENT(IN)  :: gridUuid

      !> Status flag
      LOGICAL,          INTENT(OUT) :: successful

      !-----------------
      ! Local variables
      !-----------------

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":FileInventoryElement_init_metadata"

      !----------------------------

      successful = .FALSE.

      me%iconVarNameLength      = LEN_TRIM(iconVarName)
      me%dataDateTimeLength     = LEN_TRIM(dataDateTime)
      me%validityDateTimeLength = LEN_TRIM(validityDateTime)

      IF (ANY([me%iconVarNameLength, me%dataDateTimeLength, me%validityDateTimeLength] < 1)) THEN
        CALL warning(routine, "One of iconVarName, dataDateTime or validityDateTime is empty")
        RETURN
      ENDIF

      me%iconVarName           = iconVarName(1:me%iconVarNameLength)
      me%dataDateTime          = dataDateTime(1:me%dataDateTimeLength)
      me%validityDateTime      = validityDateTime(1:me%validityDateTimeLength)

      me%discipline            = discipline
      me%parameterCategory     = parameterCategory
      me%parameterNumber       = parameterNumber
      me%levelType             = levelType
      me%gridNumber            = gridNumber
      me%gridPosition          = gridPosition
      me%runClass              = runClass
      me%experimentId          = experimentId
      me%generatingProcessType = generatingProcessType
      me%gridUuid              = gridUuid
      ! Default initialization for statistics
      me%min                   = MISSING_VALUE
      me%max                   = MISSING_VALUE
      me%mean_sum              = MISSING_VALUE
      me%statistics_counter    = 0
      me%level_counter         = 0
      me%uniform_level_counter = 0
      me%maxTileId             = -1
      me%missingValuesPresent  = .FALSE.
      me%isValid               = .FALSE.

      ! If this subroutine is called, we take it as evidence
      ! that the calling Work PE got a GRIB record for decoding
      me%gotRecord             = .TRUE.

      successful = .TRUE.

    END SUBROUTINE FileInventoryElement_init_metadata

    !>
    !! @brief Initialize statistics of an instance of t_FileInventoryElement
    !!
    SUBROUTINE FileInventoryElement_init_statistics(me, min, max, mean, tileId, isUniform, missingValuesPresent, successful)

      !-----------
      ! Arguments
      !-----------

      !> Passed-object dummy argument
      CLASS(t_FileInventoryElement), INTENT(INOUT) :: me

      !> Minimum, maximum and average of field values
      REAL(wp), INTENT(IN)  :: min, max, mean

      !> Tile index (if applicable)
      INTEGER,  INTENT(IN)  :: tileId

      !> Is field uniform within level/layer?
      LOGICAL,  INTENT(IN)  :: isUniform

      !> Are missing values present within level/layer?
      LOGICAL,  INTENT(IN)  :: missingValuesPresent

      !> Status flag
      LOGICAL,  INTENT(OUT) :: successful

      !-----------------
      ! Local variables
      !-----------------

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":FileInventoryElement_init_statistics"

      !----------------------------

      successful = .FALSE.

      ! Initialization of statistics takes place after initialization of metadata
      IF (.NOT. me%gotRecord) THEN
        CALL warning(routine, "Initialization of statistics can take place only after init. of metadata")
        RETURN
      ELSEIF (ANY([me%iconVarNameLength, me%dataDateTimeLength, me%validityDateTimeLength] < 1)) THEN
        CALL warning(routine, "One of iconVarName, dataDateTime or validityDateTime is empty")
        RETURN
      ELSEIF (ANY([me%statistics_counter, me%level_counter, me%uniform_level_counter] > 1)) THEN
        CALL warning(routine, "One of statistics_counter, level_counter or uniform_level_counter is > 1")
        RETURN
      ENDIF

      me%min                   = min
      me%max                   = max
      me%mean_sum              = mean
      me%statistics_counter    = 1
      me%level_counter         = 1
      me%uniform_level_counter = MERGE(1, 0, isUniform)
      me%missingValuesPresent  = (me%missingValuesPresent .OR. missingValuesPresent)
      me%maxTileId             = MERGE(tileId, -1, (tileId > 0))

      ! If this subroutine is called, we take it as evidence
      ! that the decoded GRIB record is valid for further processing
      me%isValid               = .TRUE.

      successful = .TRUE.

    END SUBROUTINE FileInventoryElement_init_statistics

    !>
    !! @brief Got GRIB record that is not required and ignored
    !!
    SUBROUTINE FileInventoryElement_reject(me)

      !> Passed-object dummy argument
      CLASS(t_FileInventoryElement), INTENT(INOUT) :: me

      !----------------------------

      ! Just to make sure
      CALL me%reset()

      me%gotRecord = .TRUE.
      me%isValid   = .FALSE.

    END SUBROUTINE FileInventoryElement_reject

    !>
    !! @brief Update field-value statistics of an instance of t_FileInventoryElement
    !!
    SUBROUTINE FileInventoryElement_update(me, inventory_element, successful)

      !-----------
      ! Arguments
      !-----------

      !> Passed-object dummy argument
      CLASS(t_FileInventoryElement), INTENT(INOUT) :: me

      !> New inventory element
      TYPE(t_FileInventoryElement),  INTENT(IN)    :: inventory_element

      !> Status flag
      LOGICAL,  INTENT(OUT) :: successful

      !-----------------
      ! Local variables
      !-----------------

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":FileInventoryElement_update"

      !----------------------------

      successful = .FALSE.

      ! Note that we will not check the content of 'me',
      ! as this subroutine should be called only after checking 'me' and 'inventory_element' for equality.
      ! Nevertheless, at least the content of one of the two should be checked.
      IF (.NOT. inventory_element%isValid) THEN
        CALL warning(routine, "New element is not qualified for update")
        RETURN
      ELSEIF (ANY([inventory_element%iconVarNameLength, inventory_element%dataDateTimeLength, &
        &      inventory_element%validityDateTimeLength] < 1)) THEN
        CALL warning(routine, "One of new iconVarName, dataDateTime or validityDateTime is empty")
        RETURN
      ELSEIF (ANY([inventory_element%statistics_counter, inventory_element%level_counter] < 1)) THEN
        CALL warning(routine, "One of new statistics_counter, or level_counter is < 1")
        RETURN
      ENDIF

      ! Update statistics
      me%min                   = MIN(me%min, inventory_element%min)
      me%max                   = MAX(me%max, inventory_element%max)
      me%mean_sum              = me%mean_sum + inventory_element%mean_sum
      me%statistics_counter    = me%statistics_counter + inventory_element%statistics_counter
      me%level_counter         = me%level_counter + inventory_element%level_counter
      me%uniform_level_counter = me%uniform_level_counter + inventory_element%uniform_level_counter
      me%maxTileId             = MAX(me%maxTileId, inventory_element%maxTileId)
      me%missingValuesPresent  = (me%missingValuesPresent .OR. inventory_element%missingValuesPresent)
      me%gotRecord             = inventory_element%gotRecord
      me%isValid               = inventory_element%isValid

      successful = .TRUE.

    END SUBROUTINE FileInventoryElement_update

    !>
    !! @brief Get metadata of an instance of t_FileInventoryElement
    !!
    SUBROUTINE FileInventoryElement_get_metadata(me, iconVarName, iconVarNameLength, dataDateTime, dataDateTimeLength,    &
      &                                          validityDateTime, validityDateTimeLength, discipline, parameterCategory, &
      &                                          parameterNumber, levelType, gridNumber, gridPosition, runClass,          &
      &                                          experimentId, generatingProcessType, gridUuid, successful)

      !-----------
      ! Arguments
      !-----------

      !> Passed-object dummy argument
      CLASS(t_FileInventoryElement), INTENT(IN) :: me

      !> ICON-internal variable name
      CHARACTER(LEN=vname_len), INTENT(OUT) :: iconVarName

      !> Actual length of variable name
      INTEGER,      INTENT(OUT) :: iconVarNameLength

      !> Reference and validity date and time in format:'YYYY-MM-DDThh:mm:ss'
      CHARACTER(LEN=max_datetime_str_len), INTENT(OUT) :: dataDateTime, validityDateTime

      !> Actual lengths of reference and validity date and time
      INTEGER,      INTENT(OUT) :: dataDateTimeLength, validityDateTimeLength

      !> Parameter triple
      INTEGER,      INTENT(OUT) :: discipline, parameterCategory, parameterNumber

       !> GRIB level type
      INTEGER,      INTENT(OUT) :: levelType

      !> Consecutive number of horizontal grid (GRIB key: numberOfGridUsed)
      INTEGER,      INTENT(OUT) :: gridNumber

      !> Identifier of horizontal-grid element (cell, edge or vertex, GRIB key: numberOfGridInReference)
      INTEGER,      INTENT(OUT) :: gridPosition

      !> Identifier of background generating process (GRIB key: backgroundProcess)
      INTEGER,      INTENT(OUT) :: runClass

      !> DWD-specific experiment identifier (GRIB key: localNumberOfExperiment)
      INTEGER,      INTENT(OUT) :: experimentId

      !> Identifier of process used to create the product (GRIB key: typeOfGeneratingProcess)
      INTEGER,      INTENT(OUT) :: generatingProcessType

      !> Fingerprint of horizontal grid (GRIB key: uuidOfHGrid)
      TYPE(t_uuid), INTENT(OUT) :: gridUuid

      !> Status flag
      LOGICAL,      INTENT(OUT) :: successful

      !-----------------
      ! Local variables
      !-----------------

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":FileInventoryElement_get_metadata"

      !----------------------------

      successful = .FALSE.

      IF (.NOT. me%isValid) THEN
        CALL warning(routine, "Inventory element has no valid content")
        RETURN
      ENDIF

      discipline            = me%discipline
      parameterCategory     = me%parameterCategory
      parameterNumber       = me%parameterNumber
      levelType             = me%levelType
      gridNumber            = me%gridNumber
      gridPosition          = me%gridPosition
      runClass              = me%runClass
      experimentId          = me%experimentId
      generatingProcessType = me%generatingProcessType
      gridUuid              = me%gridUuid

      IF (me%iconVarNameLength > 0) THEN
        iconVarName       = me%iconVarName(1:me%iconVarNameLength)
        iconVarNameLength = me%iconVarNameLength
      ELSE
        iconVarName       = " "
        iconVarNameLength = 0
        CALL warning(routine, "iconVarName is empty")
        RETURN
      ENDIF

      IF (me%dataDateTimeLength > 0) THEN
        dataDateTime       = me%dataDateTime(1:me%dataDateTimeLength)
        dataDateTimeLength = me%dataDateTimeLength
      ELSE
        dataDateTime       = " "
        dataDateTimeLength = 0
        CALL warning(routine, "dataDateTime is empty")
        RETURN
      ENDIF

      IF (me%validityDateTimeLength > 0) THEN
        validityDateTime       = me%validityDateTime(1:me%validityDateTimeLength)
        validityDateTimeLength = me%validityDateTimeLength
      ELSE
        CALL warning(routine, "validityDateTime is empty")
        validityDateTime       = " "
        validityDateTimeLength = 0
        RETURN
      ENDIF

      successful = .TRUE.

    END SUBROUTINE FileInventoryElement_get_metadata

    !>
    !! @brief Get statistics of an instance of t_FileInventoryElement
    !!
    SUBROUTINE FileInventoryElement_get_statistics(me, min, max, mean, nlev, ntiles, nunival, missval, successful)

      !-----------
      ! Arguments
      !-----------

      !> Passed-object dummy argument
      CLASS(t_FileInventoryElement), INTENT(IN) :: me

      !> Minimum, maximum and average of field values
      REAL(wp), INTENT(OUT) :: min, max, mean

      !> Number of layers or levels
      INTEGER,  INTENT(OUT) :: nlev

      !> Number of tiles
      INTEGER,  INTENT(OUT) :: ntiles

      !> Number of levels/layers with uniform field values
      INTEGER,  INTENT(OUT) :: nunival

      !> Missing-value flag
      LOGICAL,  INTENT(OUT) :: missval

      !> Status flag
      LOGICAL,  INTENT(OUT) :: successful

      !-----------------
      ! Local variables
      !-----------------

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":FileInventoryElement_get_statistics"

      !----------------------------

      successful = .FALSE.

      ! Initialize intent-out arguments
      min     = MISSING_VALUE
      max     = MISSING_VALUE
      mean    = MISSING_VALUE
      nlev    = -1
      ntiles  = -1
      nunival = -1
      missval = .FALSE.

      IF (.NOT. me%isValid) THEN
        CALL warning(routine, "Inventory element has no valid content")
        RETURN
      ELSEIF (ANY([me%iconVarNameLength, me%dataDateTimeLength, me%validityDateTimeLength] < 1)) THEN
        CALL warning(routine, "One of iconVarName, dataDateTime or validityDateTime is empty")
        RETURN
      ELSEIF (ANY([me%statistics_counter, me%level_counter] < 1)) THEN
        CALL warning(routine, "One of statistics_counter or level_counter is < 1")
        RETURN
      ENDIF

      min  = me%min
      max  = me%max
      mean = me%mean_sum / REAL(me%statistics_counter, KIND=wp)

      missval = me%missingValuesPresent

      IF (me%maxTileId < 1) THEN
        ! A non-tiled field
        ntiles  = 0
        nlev    = me%level_counter
        nunival = me%uniform_level_counter
      ELSE
        ! A tile field
        ntiles = me%maxTileId
        ! Each level should have been encountered ntiles times
        nlev   = me%level_counter / me%maxTileId
        ! Note that we assume that if the field of one tile is uniform within the level/layer,
        ! the fields of all other tiles are uniform, too!
        ! Counting per tile would be too cumbersome.
        nunival = me%uniform_level_counter / me%maxTileId
      ENDIF

      successful = .TRUE.

    END SUBROUTINE FileInventoryElement_get_statistics

    !>
    !! @brief Reset content of an instance of t_FileInventoryElement to default values
    !!
    SUBROUTINE FileInventoryElement_reset(me)

      !> Passed-object dummy argument
      CLASS(t_FileInventoryElement), INTENT(INOUT) :: me

      !----------------------------

      me%iconVarName            = " "
      me%dataDateTime           = " "
      me%validityDateTime       = " "
      me%discipline             = ECC_NULL
      me%parameterCategory      = ECC_NULL
      me%parameterNumber        = ECC_NULL
      me%levelType              = ECC_NULL
      me%gridNumber             = ECC_NULL
      me%gridPosition           = ECC_NULL
      me%runClass               = ECC_NULL
      me%experimentId           = ECC_NULL
      me%generatingProcessType  = ECC_NULL
      ! No default value for me%gridUuid
      me%iconVarNameLength      = 0
      me%dataDateTimeLength     = 0
      me%validityDateTimeLength = 0
      me%min                    = MISSING_VALUE
      me%max                    = MISSING_VALUE
      me%mean_sum               = MISSING_VALUE
      me%statistics_counter     = 0
      me%level_counter          = 0
      me%uniform_level_counter  = 0
      me%maxTileId              = -1
      me%missingValuesPresent   = .FALSE.
      me%gotRecord              = .FALSE.
      me%isValid                = .FALSE.

    END SUBROUTINE FileInventoryElement_reset

    !>
    !! @brief Assignment operator for instances of t_FileInventoryElement: lhs = rhs
    !!
    PURE SUBROUTINE FileInventoryElement_assign(lhs, rhs)

      !> Left-hand side of assignment
      CLASS(t_FileInventoryElement), INTENT(INOUT) :: lhs

      !> Right-hand side of assignment
      TYPE(t_FileInventoryElement), INTENT(IN) :: rhs

      !----------------------------

      ! No error handling in assignment operator

      lhs%discipline            = rhs%discipline
      lhs%parameterCategory     = rhs%parameterCategory
      lhs%parameterNumber       = rhs%parameterNumber
      lhs%levelType             = rhs%levelType
      lhs%gridNumber            = rhs%gridNumber
      lhs%gridPosition          = rhs%gridPosition
      lhs%runClass              = rhs%runClass
      lhs%experimentId          = rhs%experimentId
      lhs%generatingProcessType = rhs%generatingProcessType
      lhs%gridUuid              = rhs%gridUuid

      IF (rhs%statistics_counter > 0) THEN
        lhs%min                = rhs%min
        lhs%max                = rhs%max
        lhs%mean_sum           = rhs%mean_sum
        lhs%statistics_counter = rhs%statistics_counter
      ELSE
        lhs%min                = MISSING_VALUE
        lhs%max                = MISSING_VALUE
        lhs%mean_sum           = MISSING_VALUE
        lhs%statistics_counter = 0
      ENDIF

      IF (rhs%iconVarNameLength > 0) THEN
        lhs%iconVarName       = rhs%iconVarName(1:rhs%iconVarNameLength)
        lhs%iconVarNameLength = rhs%iconVarNameLength
      ELSE
        lhs%iconVarName       = " "
        lhs%iconVarNameLength = 0
      ENDIF

      IF (rhs%dataDateTimeLength > 0) THEN
        lhs%dataDateTime       = rhs%dataDateTime(1:rhs%dataDateTimeLength)
        lhs%dataDateTimeLength = rhs%dataDateTimeLength
      ELSE
        lhs%dataDateTime       = " "
        lhs%dataDateTimeLength = 0
      ENDIF

      IF (rhs%validityDateTimeLength > 0) THEN
        lhs%validityDateTime       = rhs%validityDateTime(1:rhs%validityDateTimeLength)
        lhs%validityDateTimeLength = rhs%validityDateTimeLength
      ELSE
        lhs%validityDateTime       = " "
        lhs%validityDateTimeLength = 0
      ENDIF

      lhs%level_counter         = MAX(0, rhs%level_counter)
      lhs%uniform_level_counter = MAX(0, rhs%uniform_level_counter)
      lhs%maxTileId             = rhs%maxTileId
      lhs%missingValuesPresent  = rhs%missingValuesPresent
      lhs%gotRecord             = rhs%gotRecord
      lhs%isValid               = rhs%isValid

    END SUBROUTINE FileInventoryElement_assign

    !>
    !! @brief Check if two instances of t_FileInventoryElement are equal: equal = (arg1 == arg2)
    !!
    PURE FUNCTION FileInventoryElement_equal(arg1, arg2) RESULT(equal)

      !> First inventory element
      CLASS(t_FileInventoryElement), INTENT(IN) :: arg1

      !> Second inventory element
      CLASS(t_FileInventoryElement), INTENT(IN) :: arg2

      !> Result
      LOGICAL :: equal

      !----------------------------

      ! No error handling in binary operator

      equal = (arg1%discipline             == arg2%discipline)             .AND. &
        &     (arg1%parameterCategory      == arg2%parameterCategory)      .AND. &
        &     (arg1%parameterNumber        == arg2%parameterNumber)        .AND. &
        &     (arg1%levelType              == arg2%levelType)              .AND. &
        &     (arg1%gridNumber             == arg2%gridNumber)             .AND. &
        &     (arg1%gridPosition           == arg2%gridPosition)           .AND. &
        &     (arg1%runClass               == arg2%runClass)               .AND. &
        &     (arg1%experimentId           == arg2%experimentId)           .AND. &
        &     (arg1%generatingProcessType  == arg2%generatingProcessType)  .AND. &
        &     (arg1%iconVarNameLength      == arg2%iconVarNameLength)      .AND. &
        &     (arg1%dataDateTimeLength     == arg2%dataDateTimeLength)     .AND. &
        &     (arg1%validityDateTimeLength == arg2%validityDateTimeLength) .AND. &
        &     ((arg1%maxTileId < 1)     .EQV. (arg2%maxTileId < 1))        .AND. &
        &     arg1%gotRecord            .AND. arg2%gotRecord               .AND. &
        &     arg1%isValid              .AND. arg2%isValid

      ! Finally, do the expensive comparisons
      IF (.NOT. equal) THEN
        RETURN
      ELSEIF (ANY([arg1%iconVarNameLength, arg1%dataDateTimeLength, arg1%validityDateTimeLength] < 1)) THEN
        equal = .FALSE.
        RETURN
      ELSEIF (arg1%iconVarName(1:arg1%iconVarNameLength) /= arg2%iconVarName(1:arg2%iconVarNameLength)) THEN
        equal = .FALSE.
        RETURN
      ELSEIF (arg1%dataDateTime(1:arg1%dataDateTimeLength) /= arg2%dataDateTime(1:arg2%dataDateTimeLength)) THEN
        equal = .FALSE.
        RETURN
      ELSEIF (arg1%validityDateTime(1:arg1%validityDateTimeLength) /= arg2%validityDateTime(1:arg2%validityDateTimeLength)) THEN
        equal = .FALSE.
        RETURN
      ENDIF

      equal = .TRUE.

    END FUNCTION FileInventoryElement_equal

    !>
    !! @brief Integrate element into file inventory
    !!
    SUBROUTINE FileInventory_add(me, inventory_element, my_rank, root_rank, comm, comm_size, reset, successful)

      !-----------
      ! Arguments
      !-----------

      !> Passed-object dummy argument
      CLASS(t_FileInventory), INTENT(INOUT) :: me

      !> Inventory element
      TYPE(t_FileInventoryElement), INTENT(INOUT)  :: inventory_element

      !> For MPI communication
      INTEGER, INTENT(IN)  :: my_rank, root_rank, comm, comm_size

      !> Reset content of inventory_element to default values
      !> after integration into file inventory?
      LOGICAL,  INTENT(IN) :: reset

      !> Status flag
      LOGICAL,  INTENT(OUT) :: successful

      !-----------------
      ! Local variables
      !-----------------

      !> Local inventory element for caching the buffer content
      TYPE(t_FileInventoryElement) :: inventory_local

      !> Buffer for metadata of type integer
      INTEGER, ALLOCATABLE :: recv_buffer_int(:,:), send_buffer_int(:)

      !> Buffer for statistics of type real
      REAL(wp), ALLOCATABLE :: recv_buffer_real(:,:), send_buffer_real(:)

      !> Buffer for metadata of type character
      CHARACTER(LEN=vname_len + 2*max_datetime_str_len + uuid_string_length) :: send_buffer_str
      CHARACTER(LEN=vname_len + 2*max_datetime_str_len + uuid_string_length), &
        & ALLOCATABLE :: recv_buffer_str(:)

      !> Local comm size for allocation
      INTEGER :: comm_size_local

      !> Flag for root PE
      LOGICAL :: i_am_root_rank

      !> Position in string
      INTEGER :: position

      !> Loop index for mpi ranks
      INTEGER :: jrank

      !> String for gridUuid
      CHARACTER(LEN=uuid_string_length) :: uuid_str

      !> Inventory-element pointer
      TYPE(t_FileInventoryElement), POINTER :: elm => NULL()

      !> Status identifier
      INTEGER :: status

      !> Flag to indicate matching inventory elements
      LOGICAL :: found_match

      !> Status flag
      LOGICAL :: successful_local

      !> Loop counter
      INTEGER :: loop_counter

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":FileInventory_add"

      !----------------------------

      ! It is far from ideal to directly access the elements of t_FileInventoryElement instances
      ! from within this t_FileInventory procedure, but MPI_GATHER leaves us little choice.

      successful = .FALSE.

      ! In order to avoid a deadlock or something like that, we use 'finish' instead of 'return' here
      IF (inventory_element%isValid .AND. &
        & ANY([inventory_element%iconVarNameLength, inventory_element%dataDateTimeLength, &
        &      inventory_element%validityDateTimeLength] < 1)) THEN
        CALL finish(routine, "One of iconVarName, dataDateTime or validityDateTime is empty")
      ELSEIF (inventory_element%isValid .AND. &
        & ANY([inventory_element%statistics_counter, inventory_element%level_counter] < 1)) THEN
        CALL finish(routine, "One of statistics_counter or level_counter is < 1")
      ENDIF

      ! Is this the root PE?
      i_am_root_rank = (my_rank == root_rank)

      ! Allocate buffers:
      comm_size_local = MERGE(comm_size, 1, i_am_root_rank)
      ALLOCATE(recv_buffer_int(19,comm_size_local), send_buffer_int(19), &
        &      recv_buffer_real(3,comm_size_local), send_buffer_real(3), &
        &      recv_buffer_str(comm_size_local), STAT=status)
      IF (status /= SUCCESS) CALL finish(routine, "Allocation of buffers failed")

      send_buffer_int(:)  = -999
      send_buffer_real(:) = -999.0_wp
      send_buffer_str     = " "
      uuid_str            = " "

      ! Each PE, which got a GRIB record for decoding, fills the send buffers:

      IF (inventory_element%gotRecord .AND. inventory_element%isValid) THEN

        send_buffer_int(1)  = inventory_element%discipline
        send_buffer_int(2)  = inventory_element%parameterCategory
        send_buffer_int(3)  = inventory_element%parameterNumber
        send_buffer_int(4)  = inventory_element%levelType
        send_buffer_int(5)  = inventory_element%gridNumber
        send_buffer_int(6)  = inventory_element%gridPosition
        send_buffer_int(7)  = inventory_element%runClass
        send_buffer_int(8)  = inventory_element%experimentId
        send_buffer_int(9)  = inventory_element%generatingProcessType
        send_buffer_int(10) = inventory_element%iconVarNameLength
        send_buffer_int(11) = inventory_element%dataDateTimeLength
        send_buffer_int(12) = inventory_element%validityDateTimeLength
        send_buffer_int(13) = inventory_element%statistics_counter
        send_buffer_int(14) = inventory_element%maxTileID
        send_buffer_int(15) = inventory_element%level_counter
        send_buffer_int(16) = inventory_element%uniform_level_counter
        send_buffer_int(17) = MERGE(10, -10, inventory_element%missingValuesPresent)
        send_buffer_int(18) = MERGE(10, -10, inventory_element%gotRecord)
        send_buffer_int(19) = MERGE(10, -10, inventory_element%isValid)

        send_buffer_real(1) = inventory_element%min
        send_buffer_real(2) = inventory_element%max
        send_buffer_real(3) = inventory_element%mean_sum

        ! 1st iconVarName
        send_buffer_str(1:inventory_element%iconVarNameLength) = &
          & inventory_element%iconVarName(1:inventory_element%iconVarNameLength)
        position = inventory_element%iconVarNameLength
        ! 2nd dataDateTime
        send_buffer_str(position+1:position+inventory_element%dataDateTimeLength) = &
          & inventory_element%dataDateTime(1:inventory_element%dataDateTimeLength)
        position = position + inventory_element%dataDateTimeLength
        ! 3rd validityDateTime
        send_buffer_str(position+1:position+inventory_element%validityDateTimeLength) = &
          & inventory_element%validityDateTime(1:inventory_element%validityDateTimeLength)
        position = position + inventory_element%validityDateTimeLength
        ! 4th gridUuid
        CALL uuid_unparse(inventory_element%gridUuid, uuid_str)
        send_buffer_str(position+1:position+uuid_string_length) = uuid_str(1:uuid_string_length)

      ELSEIF (inventory_element%gotRecord) THEN

        ! Just to make sure (18 => gotRecord, 19 => isValid)
        send_buffer_int(18) = 10
        send_buffer_int(19) = -10

      ELSE

        ! Just to make sure
        send_buffer_int(18) = -10
        send_buffer_int(19) = -10

      ENDIF ! IF (inventory_element%gotRecord)

      ! The root PE gathers all buffer content
      CALL p_gather(sendbuf = send_buffer_int, & ! in
        &           recvbuf = recv_buffer_int, & ! inout
        &           p_dest  = root_rank,       & ! in
        &           comm    = comm             ) ! optin

      CALL p_gather(sendbuf = send_buffer_real, &
        &           recvbuf = recv_buffer_real, &
        &           p_dest  = root_rank,        &
        &           comm    = comm              )

      CALL p_gather(sbuf    = send_buffer_str, &
        &           recvbuf = recv_buffer_str, &
        &           p_dest  = root_rank,       &
        &           comm    = comm             )

      successful_local = .TRUE.

      IF (i_am_root_rank) THEN

        RANK_LOOP: DO jrank = 1, comm_size

          ! Does current buffer hold content from a valid GRIB record?
          IF ((recv_buffer_int(18, jrank) > 0) .AND. (recv_buffer_int(19, jrank) > 0)) THEN

            inventory_local%discipline             = recv_buffer_int(1, jrank)
            inventory_local%parameterCategory      = recv_buffer_int(2, jrank)
            inventory_local%parameterNumber        = recv_buffer_int(3, jrank)
            inventory_local%levelType              = recv_buffer_int(4, jrank)
            inventory_local%gridNumber             = recv_buffer_int(5, jrank)
            inventory_local%gridPosition           = recv_buffer_int(6, jrank)
            inventory_local%runClass               = recv_buffer_int(7, jrank)
            inventory_local%experimentId           = recv_buffer_int(8, jrank)
            inventory_local%generatingProcessType  = recv_buffer_int(9, jrank)
            inventory_local%iconVarNameLength      = recv_buffer_int(10, jrank)
            inventory_local%dataDateTimeLength     = recv_buffer_int(11, jrank)
            inventory_local%validityDateTimeLength = recv_buffer_int(12, jrank)
            inventory_local%statistics_counter     = recv_buffer_int(13, jrank)
            inventory_local%maxTileID              = recv_buffer_int(14, jrank)
            inventory_local%level_counter          = recv_buffer_int(15, jrank)
            inventory_local%uniform_level_counter  = recv_buffer_int(16, jrank)
            inventory_local%missingValuesPresent   = (recv_buffer_int(17, jrank) > 0)
            inventory_local%gotRecord              = (recv_buffer_int(18, jrank) > 0)
            inventory_local%isValid                = (recv_buffer_int(19, jrank) > 0)

            inventory_local%min      = recv_buffer_real(1, jrank)
            inventory_local%max      = recv_buffer_real(2, jrank)
            inventory_local%mean_sum = recv_buffer_real(3, jrank)

            ! 1st iconVarName
            inventory_local%iconVarName(1:inventory_local%iconVarNameLength) = &
              & recv_buffer_str(jrank)(1:inventory_local%iconVarNameLength)
            position = inventory_local%iconVarNameLength
            ! 2nd dataDateTime
            inventory_local%dataDateTime(1:inventory_local%dataDateTimeLength) = &
              & recv_buffer_str(jrank)(position+1:position+inventory_local%dataDateTimeLength)
            position = position + inventory_local%dataDateTimeLength
            ! 3rd validityDateTime
            inventory_local%validityDateTime(1:inventory_local%validityDateTimeLength) = &
              & recv_buffer_str(jrank)(position+1:position+inventory_local%validityDateTimeLength)
            position = position + inventory_local%validityDateTimeLength
            ! 4th gridUuid
            uuid_str(1:uuid_string_length) = recv_buffer_str(jrank)(position+1:position+uuid_string_length)
            CALL uuid_parse(uuid_str, inventory_local%gridUuid)

            ! Add current inventory element to file inventory:

            ! Update linked list
            IF (.NOT. ASSOCIATED(me%last_element)) THEN

              ! Initialization of file inventory
              ALLOCATE(me%first_element, STAT=status)
              IF (status /= SUCCESS) CALL finish(routine, "Allocation of first file inventory element failed")
              ! Assignment
              me%first_element = inventory_local
              me%last_element => me%first_element
              ! Initialize inventory-element counter
              me%element_counter = 1

            ELSE

              ! Traverse the linked list to search for a matching node
              elm              => me%first_element
              found_match      = .FALSE.
              loop_counter     = 0
              SEARCH_LOOP: DO WHILE (ASSOCIATED(elm))
                ! In case we find a match, we update the field-value statistics and exit the loop
                found_match = (elm == inventory_local)
                IF (found_match) THEN
                  ! Integrate new inventory element into
                  ! matching element of file inventory
                  CALL elm%update(inventory_element=inventory_local, successful=successful_local)
                  EXIT SEARCH_LOOP
                ENDIF ! IF (found_match)
                ! Just for safety reasons
                loop_counter = loop_counter + 1
                IF (loop_counter == me%element_counter) EXIT SEARCH_LOOP
                ! Set pointer to next list element
                elm => elm%next
              END DO SEARCH_LOOP
              elm => NULL()
              IF (.NOT. successful_local) THEN
                CALL warning(routine, "Integration of inventory element into file inventory failed")
                EXIT RANK_LOOP
              ENDIF

              ! If we found no match, we extend the inventory by a new element
              IF (.NOT. found_match) THEN
                ALLOCATE(me%last_element%next, STAT=status)
                IF (status /= SUCCESS) CALL finish(routine, "Allocation of next file inventory element failed")
                ! Assignment
                me%last_element%next = inventory_local
                me%last_element      => me%last_element%next
                ! Increment counter of inventory elements
                me%element_counter = me%element_counter + 1
              ENDIF ! IF (.NOT. found_match)

            ENDIF ! IF (.NOT. ASSOCIATED(me%last_element))

            ! Update counter of records in GRIB file
            me%record_counter = me%record_counter + 1

          ELSEIF (recv_buffer_int(18, jrank) > 0) THEN

            ! In this case, the GRIB record is not required and ignored
            me%record_counter    = me%record_counter + 1
            me%rejection_counter = me%rejection_counter + 1

          ENDIF ! IF ((recv_buffer_int(18, jrank) > 0) .AND. (recv_buffer_int(19, jrank) > 0)) ...

        END DO RANK_LOOP

      ENDIF ! IF (i_am_root_rank)

      DEALLOCATE(recv_buffer_int, send_buffer_int, recv_buffer_real, send_buffer_real, recv_buffer_str, STAT=status)
      IF (status /= SUCCESS) CALL finish(routine, "Deallocation of buffers failed")

      ! Reset content of invetory_element to default values, if required
      IF (reset) CALL inventory_element%reset()

      successful = successful_local

    END SUBROUTINE FileInventory_add

    !>
    !! @brief Clear the file inventory
    !!
    SUBROUTINE FileInventory_destruct(me, successful)

      !-----------
      ! Arguments
      !-----------

      !> Passed-object dummy argument
      CLASS(t_FileInventory), INTENT(INOUT) :: me

      !> Status flag
      LOGICAL,  INTENT(OUT) :: successful

      !-----------------
      ! Local variables
      !-----------------

      !> Inventory-element pointers
      TYPE(t_FileInventoryElement), POINTER :: elm => NULL(), &
        &                                      tmp => NULL()

      !> Status identifier
      INTEGER :: status

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":FileInventory_destruct"

      !----------------------------

      successful = .FALSE.

      ! Initialize auxiliary pointer elm
      elm => me%first_element

      CLEAR_LOOP: DO WHILE (ASSOCIATED(elm))
        tmp => elm
        elm => elm%next
        DEALLOCATE(tmp, STAT=status)
        IF (status /= SUCCESS) CALL finish(routine, "Destruction of file inventory failed")
      END DO CLEAR_LOOP

      me%last_element  => NULL()
      me%first_element => NULL()
      elm              => NULL()
      tmp              => NULL()

      me%element_counter   = 0
      me%record_counter    = 0
      me%rejection_counter = 0

      successful = .TRUE.

    END SUBROUTINE FileInventory_destruct

    !>
    !! @brief Print the file inventory
    !!
    SUBROUTINE FileInventory_print(me, grib_file_path, lIsFg, jg, legend, clear, successful)

      !-----------
      ! Arguments
      !-----------

      !> Passed-object dummy argument
      CLASS(t_FileInventory), INTENT(INOUT) :: me

      !> Path of GRIB file
      CHARACTER(LEN=*), INTENT(IN)  :: grib_file_path

      !> Type of products:
      LOGICAL,          INTENT(IN)  :: lIsFg

      !> Patch index
      INTEGER,          INTENT(IN)  :: jg

      !> Print legend?
      LOGICAL,          INTENT(IN)  :: legend

      !> Clear the file inventory after printing?
      LOGICAL,          INTENT(IN)  :: clear

      !> Status flag
      LOGICAL,          INTENT(OUT) :: successful

      !-----------------
      ! Local variables
      !-----------------

      !> Inventory-element pointer
      TYPE(t_FileInventoryElement), POINTER :: elm => NULL()

      !> File-inventory table
      TYPE(t_table) :: inventory_table

      !> Table row counter
      INTEGER :: row_counter

      !> ICON-internal variable name
      CHARACTER(LEN=vname_len) :: iconVarName

      !> Actual length of variable name
      INTEGER :: iconVarNameLength

      !> Reference and validity date and time in format:'YYYY-MM-DDThh:mm:ss'
      CHARACTER(LEN=max_datetime_str_len) :: dataDateTime, validityDateTime

      !> Actual lengths of reference and validity date and time
      INTEGER :: dataDateTimeLength, validityDateTimeLength

      !> Parameter triple
      INTEGER :: discipline, parameterCategory, parameterNumber

       !> GRIB level type
      INTEGER :: levelType

      !> Consecutive number of horizontal grid (GRIB key: numberOfGridUsed)
      INTEGER :: gridNumber

      !> Identifier of horizontal-grid element (cell, edge or vertex, GRIB key: numberOfGridInReference)
      INTEGER :: gridPosition

      !> Identifier of background generating process (GRIB key: backgroundProcess)
      INTEGER :: runClass

      !> DWD-specific experiment identifier (GRIB key: localNumberOfExperiment)
      INTEGER :: experimentId

      !> Identifier of process used to create the product (GRIB key: typeOfGeneratingProcess)
      INTEGER :: generatingProcessType

      !> Fingerprint of horizontal grid (GRIB key: uuidOfHGrid)
      TYPE(t_uuid) :: gridUuid

      !> Minimum, maximum and average of field values
      REAL(wp) :: min, max, mean

      !> Number of layers or levels
      INTEGER :: nlev

      !> Number of tiles
      INTEGER :: ntiles

      !> Number of levels/layers with uniform field values
      INTEGER  :: nunival

      !> Missing-value flag
      LOGICAL  :: missval

      !> General purpose string
      CHARACTER(LEN=30) :: string

      !> Status flag
      LOGICAL :: successful_local

      !> Table column headings
      CHARACTER(LEN=*), PARAMETER ::  col_iconVarName           = "variable",  &
        &                             col_triple                = "triple",    &
        &                             col_dataDateTime          = "ref. time", &
        &                             col_validityDateTime      = "val. time", &
        &                             col_levelType             = "levtype",   &
        &                             col_nlev                  = "nlev",      &
        &                             col_ntiles                = "ntiles",    &
        &                             col_generatingProcessType = "runtype",   &
        &                             col_runClass              = "class",     &
        &                             col_experimentId          = "expid",     &
        &                             col_gridNumber            = "gnum",      &
        &                             col_gridPosition          = "gpos",      &
        &                             col_missval               = "missval",   &
        &                             col_nunival               = "nunival",   &
        &                             col_min                   = "min",       &
        &                             col_mean                  = "mean",      &
        &                             col_max                   = "max"

      !> Procedure name
      CHARACTER(LEN=*), PARAMETER :: routine = modname//":FileInventory_print"

      !----------------------------

      successful = .FALSE.

      ! Is the inventory empty?
      IF (me%element_counter < 1) THEN
        CALL warning(routine, "File inventory is empty")
        RETURN
      ENDIF

      string = " "

      CALL initialize_table(table=inventory_table)

      ! Set 17 table columns
      CALL add_table_column(table=inventory_table, column_title=col_iconVarName)
      CALL add_table_column(table=inventory_table, column_title=col_triple)
      CALL add_table_column(table=inventory_table, column_title=col_dataDateTime)
      CALL add_table_column(table=inventory_table, column_title=col_validityDateTime)
      CALL add_table_column(table=inventory_table, column_title=col_levelType)
      CALL add_table_column(table=inventory_table, column_title=col_nlev)
      CALL add_table_column(table=inventory_table, column_title=col_ntiles)
      CALL add_table_column(table=inventory_table, column_title=col_generatingProcessType)
      CALL add_table_column(table=inventory_table, column_title=col_runClass)
      CALL add_table_column(table=inventory_table, column_title=col_experimentId)
      CALL add_table_column(table=inventory_table, column_title=col_gridNumber)
      CALL add_table_column(table=inventory_table, column_title=col_gridPosition)
      CALL add_table_column(table=inventory_table, column_title=col_missval)
      CALL add_table_column(table=inventory_table, column_title=col_nunival)
      CALL add_table_column(table=inventory_table, column_title=col_min)
      CALL add_table_column(table=inventory_table, column_title=col_mean)
      CALL add_table_column(table=inventory_table, column_title=col_max)

      ! Fill the table
      elm              => me%first_element
      row_counter      = 0
      successful_local = .TRUE.
      INVENTORY_LOOP: DO WHILE (ASSOCIATED(elm))

        ! Update row counter
        row_counter = row_counter + 1

        ! Get metadata of inventory element
        CALL elm%get(iconVarName            = iconVarName,            & ! out
          &          iconVarNameLength      = iconVarNameLength,      & ! out
          &          dataDateTime           = dataDateTime,           & ! out
          &          dataDateTimeLength     = dataDateTimeLength,     & ! out
          &          validityDateTime       = validityDateTime,       & ! out
          &          validityDateTimeLength = validityDateTimeLength, & ! out
          &          discipline             = discipline,             & ! out
          &          parameterCategory      = parameterCategory,      & ! out
          &          parameterNumber        = parameterNumber,        & ! out
          &          levelType              = levelType,              & ! out
          &          gridNumber             = gridNumber,             & ! out
          &          gridPosition           = gridPosition,           & ! out
          &          runClass               = runClass,               & ! out
          &          experimentId           = experimentId,           & ! out
          &          generatingProcessType  = generatingProcessType,  & ! out
          &          gridUuid               = gridUuid,               & ! out
          &          successful             = successful_local        ) ! out

        IF (.NOT. successful_local) THEN
          CALL warning(routine, "Getting metadata of inventory element failed")
          EXIT INVENTORY_LOOP
        ENDIF

        ! Get statistics of inventory element
        CALL elm%get(min        = min,            & ! out
          &          max        = max,            & ! out
          &          mean       = mean,           & ! out
          &          nlev       = nlev,           & ! out
          &          ntiles     = ntiles,         & ! out
          &          nunival    = nunival,        & ! out
          &          missval    = missval,        & ! out
          &          successful = successful_local) ! out

        IF (.NOT. successful_local) THEN
          CALL warning(routine, "Getting statistics of inventory element failed")
          EXIT INVENTORY_LOOP
        ENDIF

        ! Set table elements:

        ! iconVarName
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_iconVarName, &
          &                  entry_str    = iconVarName(1:iconVarNameLength))

        ! triple
        WRITE(string, fmt="(I3,1X,I3,1X,I3)") discipline, parameterCategory, parameterNumber
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_triple,      &
          &                  entry_str    = TRIM(string))

        ! Reference date and time
        CALL set_table_entry(table        = inventory_table,  &
          &                  irow         = row_counter,      &
          &                  column_title = col_dataDateTime, &
          &                  entry_str    = dataDateTime(1:dataDateTimeLength))

        ! Validity date and time
        CALL set_table_entry(table        = inventory_table,      &
          &                  irow         = row_counter,          &
          &                  column_title = col_validityDateTime, &
          &                  entry_str    = validityDateTime(1:validityDateTimeLength))

        ! Level/layer type
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_levelType,   &
          &                  entry_str    = TRIM(int2string(levelType)))

        ! Number of levels/layers
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_nlev,        &
          &                  entry_str    = TRIM(int2string(nlev)))

        ! Number of tiles
        IF (ntiles > 0) THEN
          string = TRIM(int2string(ntiles))
        ELSE
          string = " "
        ENDIF
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_ntiles,      &
          &                  entry_str    = TRIM(string))

        ! Generating process type (run type)
        CALL set_table_entry(table        = inventory_table,           &
          &                  irow         = row_counter,               &
          &                  column_title = col_generatingProcessType, &
          &                  entry_str    = TRIM(int2string(generatingProcessType)))

        ! Identifier of background generating process (run class)
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_runClass,    &
          &                  entry_str    = TRIM(int2string(runClass)))

        ! DWD-specific experiment identifier (experiment number)
        IF (experimentId < 0) THEN
          string = " "
        ELSE
          string = TRIM(int2string(experimentId))
        ENDIF
        CALL set_table_entry(table        = inventory_table,  &
          &                  irow         = row_counter,      &
          &                  column_title = col_experimentId, &
          &                  entry_str    = TRIM(string))

        ! Consecutive number of horizontal grid (grid number)
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_gridNumber,  &
          &                  entry_str    = TRIM(int2string(gridNumber)))

        ! Identifier of horizontal-grid element (cell, edge or vertex)
        SELECT CASE(gridPosition)
        CASE(1)
          string = "cell"
        CASE(2)
          string = "vert"
        CASE(3)
          string = "edge"
        CASE DEFAULT
          successful_local = .FALSE.
          CALL warning(routine, "Invalid grid position")
          EXIT INVENTORY_LOOP
        END SELECT
        CALL set_table_entry(table        = inventory_table,  &
          &                  irow         = row_counter,      &
          &                  column_title = col_gridPosition, &
          &                  entry_str    = TRIM(string))


        ! Are there missing values within the field?
        string = MERGE("X", " ", missval)
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_missval,     &
          &                  entry_str    = TRIM(string))

        ! Number of levels/layers with uniform field values
        IF (nunival > 0) THEN
          string = TRIM(int2string(nunival))
        ELSE
          string = " "
        ENDIF
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_nunival,     &
          &                  entry_str    = TRIM(string))

        ! Statistics: minimum of field values
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_min,         &
          &                  entry_str    = TRIM(real2string(min)))

        ! Statistics: average of field values
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_mean,        &
          &                  entry_str    = TRIM(real2string(mean)))

        ! Statistics: maximum of field values
        CALL set_table_entry(table        = inventory_table, &
          &                  irow         = row_counter,     &
          &                  column_title = col_max,         &
          &                  entry_str    = TRIM(real2string(max)))

        ! Just for safety reasons
        IF (row_counter == me%element_counter) EXIT INVENTORY_LOOP

        ! Set pointer to next list element
        elm => elm%next

      END DO INVENTORY_LOOP
      elm => NULL()

      IF (.NOT. successful_local) RETURN

      ! Print the inventory table:
      IF (lIsFg) THEN
        string = "first-guess"
      ELSE
        string = "analysis"
      ENDIF
      WRITE(0, '(A)')    " "
      WRITE(0, '(A)')    "-------------------------------------------------------------------------------------"
      WRITE(0, '(A)')    " Inventory of "//TRIM(string)//" file: "//TRIM(grib_file_path)
      WRITE(0, '(A,I0)') " - For domain:.................. ", jg
      WRITE(0, '(A,I0)') " - Total number of GRIB records: ", me%record_counter
      WRITE(0, '(A,I0)') " - Of them are ignored:......... ", me%rejection_counter
      WRITE(0, '(A)')    "-------------------------------------------------------------------------------------"

      CALL print_table(table=inventory_table, opt_delimiter=" | ")

      WRITE(0, '(A)')    "-------------------------------------------------------------------------------------"
      IF (legend) THEN
        ! Print a legend for the column headings
        WRITE(0, '(A)')  "Legend:"
        WRITE(0, '(A)')  "-------"
        WRITE(0, '(A)')  " - variable:.......Internal variable name (without runtime suffixes)"
        WRITE(0, '(A)')  " - triple:.........Parameter triple: discipline parameterCategory parameterNumber"
        WRITE(0, '(A)')  " - ref. time:......Reference date and time of data"
        WRITE(0, '(A)')  " - val. time:......Validity date and time of data (> ref. time, in general)"
        WRITE(0, '(A)')  " - levtype:........GRIB level identifier: typeOfFirstFixedSurface"
        WRITE(0, '(A)')  " - nlev:...........Number of levels or layers (> 1 for 3d fields)"
        WRITE(0, '(A)')  " - ntiles:.........Number of (surface/soil) tiles (if applicable, otherwise a blank)"
        WRITE(0, '(A)')  " - runtype:........Integer identifier of type of data-generating process"
        WRITE(0, '(A)')  "                   (typeOfGeneratingProcess: 0 'Analysis', 2 'Forecast', ...)"
        WRITE(0, '(A)')  " - class:..........Integer identifier for generating-centre-specific background process"
        WRITE(0, '(A)')  "                   (backgroundProcess: 0 'Main run', 1 'Pre-assimilation', 2 'Assimilation')"
        WRITE(0, '(A)')  " - expid:..........Generating-centre-specific experiment number (localNumberOfExperiment,"
        WRITE(0, '(A)')  "                   if applicable, otherwise a blank)"
        WRITE(0, '(A)')  " - gnum:...........Consecutive number of horizontal grid"
        WRITE(0, '(A)')  " - gpos:...........Position/element of horizontal grid: cell, edge, vertex"
        WRITE(0, '(A)')  " - missval:........'X'/' ': Missing values are/are not present in the field"
        WRITE(0, '(A)')  " - nunival:........Number of levels or layers where field values are uniform"
        WRITE(0, '(A)')  "                   (if applicable, otherwise a blank)"
        WRITE(0, '(A)')  " - min, mean, max: Min. mean and max. values of the entire (2d or 3d) field"
        WRITE(0, '(A)')  "-------------------------------------------------------------------------------------"
      ENDIF ! IF (legend)
      WRITE(0, '(A)')    " "

      ! Destruct table
      CALL finalize_table(table=inventory_table)

      ! Clear inventory if required
      IF (clear) CALL me%destruct(successful=successful_local)
      IF (.NOT. successful_local) THEN
        CALL warning(routine, "Destruction of file inventory failed")
        RETURN
      ENDIF

      successful = .TRUE.

    END SUBROUTINE FileInventory_print

END MODULE mo_input_request_list

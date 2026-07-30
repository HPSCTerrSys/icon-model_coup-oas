!
! MODULE: mo_oasis_land_coupling
! Subroutines for OASIS3-MCT land-atmosphere coupling
!
#ifdef COUP_OAS_ECLM

MODULE mo_oasis_land_coupling
  !
  ! Subroutines for OASIS3-MCT land-atmosphere coupling
  ! This module handles OASIS construction, mask creation, field sending, and field receiving
  !

  USE mo_oasis_define
  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_parallel_config,   ONLY: nproma, idx_1d
  USE mo_run_config,        ONLY: msg_level, iqv, &
                                   config_loascplscheme_exchcoef => loascplscheme_exchcoef
  USE mo_dynamics_config,   ONLY: nnow, nnow_rcf
  USE mo_nonhydro_state,    ONLY: p_nh_state
  USE mo_nonhydro_types,    ONLY: t_nh_state
  USE mo_impl_constants,    ONLY: min_rlcell_int, grf_bdywidth_c
  USE mo_model_domain,      ONLY: p_patch, t_patch
  USE mo_ext_data_state,    ONLY: ext_data
  USE mo_ext_data_types,    ONLY: t_external_data
  USE mo_nwp_phy_state,     ONLY: prm_diag
  USE mo_nwp_phy_types,     ONLY: t_nwp_phy_diag
  USE mo_nwp_lnd_types,     ONLY: t_lnd_state
  USE mo_loopindices,       ONLY: get_indices_c 
  USE mo_grid_config,       ONLY: n_dom
  USE mo_mpi,               ONLY: p_pe_work, p_comm_work_only, mpi_comm_null, p_sum
  USE mod_oasis,            ONLY: oasis_abort, oasis_set_couplcomm, oasis_enddef,      &
                                  oasis_def_partition, oasis_def_var, OASIS_In,        &
                                  OASIS_Out, OASIS_Real

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: construct_oasis_coupling
  PUBLIC :: oasis_create_mask
  PUBLIC :: oasis_send_fields
  PUBLIC :: oasis_receive_fields
  
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_oasis_land_coupling'

CONTAINS

  !-----------------------------------------------------------------------------
  !
  ! OASIS3-MCT subroutines
  !
  !-----------------------------------------------------------------------------

  SUBROUTINE construct_oasis_coupling(p_patch, ext_data)
    !
    ! Constructs and initializes the OASIS3-MCT coupling infrastructure
    ! This subroutine:
    ! - Sets up the coupling communicator
    ! - Defines partitions and grid decomposition
    ! - Defines send and receive variables
    ! - Allocates coupling field buffers
    ! - Creates the coupling mask
    ! - Finalizes OASIS definitions
    !
    ! Arguments
    TYPE(t_patch), INTENT(IN) :: p_patch(:)   !< Unstructured triangular grid and patch info
    TYPE(t_external_data), INTENT(IN) :: ext_data(:)  !< External data (land-sea mask, lake mask)


    ! Local variables
    INTEGER :: i_startblk, i_endblk, jb, jc, i_startidx, i_endidx
    INTEGER :: rl_start, rl_end, c, cc, cl, jg, ic
    CHARACTER(len=1024) :: oas_message
    INTEGER :: oas_prt = 6
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':construct_oasis_coupling'

    IF(msg_level >= 30 ) THEN
      WRITE(oas_prt,*) 'ICONOAS: ',TRIM(routine),' OASIS coupling initialization started'
      FLUSH(oas_prt)
    END IF

    ! take the (namelist-configurable) coupling scheme choice from run_nml
    loascplscheme_exchcoef = config_loascplscheme_exchcoef

    IF(msg_level >= 10 ) THEN
      WRITE(oas_message,*) 'ICONOAS: loascplscheme_exchcoef = ', loascplscheme_exchcoef
      CALL message(routine, oas_message)
    END IF

    ! local coupling communicator for work processors
    IF ( p_pe_work < p_patch(1)%proc0) THEN
      lpe_cpl = .FALSE.
      CALL oasis_set_couplcomm(mpi_comm_null)
    ELSE
      lpe_cpl = .TRUE.
      CALL oasis_set_couplcomm(p_comm_work_only)
    ENDIF

    IF(msg_level >= 30 ) THEN
      WRITE(oas_prt,*) 'ICONOAS: ',TRIM(routine),' oasis set couplcomm'
      FLUSH(oas_prt)
    END IF

    !  Check domain count
    IF (n_dom > 1) THEN
      CALL oasis_abort(oas_comp_id, oas_comp_name, &
        'Number of ICON domains > 1 when coupled to CLM.')
    END IF

    !  Process OASIS partition setup on coupling PEs
    IF ( lpe_cpl ) THEN

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int
      i_startblk = p_patch(1)%cells%start_blk(rl_start, 1)
      i_endblk   = p_patch(1)%cells%end_blk(rl_end, MAX(1,p_patch(1)%n_childdom))

      ! (1) Calculate length of oas_part on current process
      c = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch(1), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        c = c + i_endidx - i_startidx + 1
      END DO
      ALLOCATE( oas_part(c+2) )
      oas_part(1) = 4 ! 4: POINTS partitioning
      oas_part(2) = c ! number of points in the partition

      cc = p_sum(oas_part(2), comm=p_comm_work_only) ! number of points active in coupling
      cl = p_patch(1)%n_patch_cells_g - cc           ! number of points lateral boundaries

      IF(msg_level >= 10 ) THEN
        WRITE(oas_prt,*) 'ICONOAS: total grid points in coupling: ',cc,' lateral boundary points: ',cl
        FLUSH(oas_prt)
      END IF

      ! (2) Fill the partition with the cell indices
      c = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch(1), jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
          c = c + 1
          ic = idx_1d(jc,jb)
          oas_part(c+2) = p_patch(1)%cells%decomp_info%glb_index(ic)-cl ! shift global index for coupling
        END DO
      END DO

      oas_nlat = cc
      oas_nlon = 1

      IF(msg_level >= 30 ) THEN
        WRITE(oas_prt,*) 'ICONOAS: ',TRIM(routine),' partition information: c=',c,' oas_nlat=',oas_nlat
        FLUSH(oas_prt)
      END IF

      ! Define partition in OASIS
      CALL oasis_def_partition(oas_part_id, oas_part, oas_error, oas_nlat*oas_nlon)
      IF (oas_error /= 0) THEN
        CALL oasis_abort(oas_comp_id, oas_comp_name, 'Failure in oasis_def_partition')
      END IF

      IF(msg_level >= 30 ) THEN
        WRITE(oas_prt,*) 'ICONOAS: ',TRIM(routine),' oasis_def_partition completed'
        FLUSH(oas_prt)
      END IF

      ! Setup partition metadata for OASIS3-MCT
      oas_var_nodims(1) = 1   ! rang of field array
      oas_var_nodims(2) = 1   ! 'bundle' always 1 in OASIS3-MCT
      oas_vshape(1)     = 1
      oas_vshape(2)     = oas_part(2)

      ! Define send variables for OASIS3-MCT
      oas_snd_meta(1)%clpname = "ICOTEMPE"
      oas_snd_meta(2)%clpname = "ICOUWIND"
      oas_snd_meta(3)%clpname = "ICOVWIND"
      oas_snd_meta(4)%clpname = "ICOSPWAT"
      oas_snd_meta(5)%clpname = "ICOTHICK"
      oas_snd_meta(6)%clpname = "ICOPRESS"
      oas_snd_meta(7)%clpname = "ICODIRSW"
      oas_snd_meta(8)%clpname = "ICODIFSW"
      oas_snd_meta(9)%clpname = "ICOLONGW"
      oas_snd_meta(10)%clpname = "ICOCVPRE"
      oas_snd_meta(11)%clpname = "ICOGSPRE"

      ! Define receive variables for OASIS3-MCT
      IF (ALLOCATED(oas_rcv_meta)) THEN
        DEALLOCATE(oas_rcv_meta)
      END IF
      IF (loascplscheme_exchcoef) THEN
        oas_rcv_nvar = 8
      ELSE
        oas_rcv_nvar = 8
      END IF
      ALLOCATE(oas_rcv_meta(oas_rcv_nvar), stat = oas_error)
      IF (oas_error > 0) THEN
        CALL oasis_abort(oas_comp_id, oas_comp_name, &
          'Failure in allocating ICON receive metadata')
      END IF
      oas_rcv_meta(1)%clpname = "ICOINFRA"
      oas_rcv_meta(2)%clpname = "ICOALBED"
      oas_rcv_meta(3)%clpname = "ICOALBEI"
      IF (loascplscheme_exchcoef) THEN
        oas_rcv_meta(4)%clpname = "ICO_RAM1"
        oas_rcv_meta(5)%clpname = "ICO_RAH1"
        oas_rcv_meta(6)%clpname = "ICO_RAW1"
        oas_rcv_meta(7)%clpname = "ICOQVSFC"
        oas_rcv_meta(8)%clpname = "ICOTGRND"
      ELSE
        oas_rcv_meta(4)%clpname = "ICO_TAUX"
        oas_rcv_meta(5)%clpname = "ICO_TAUY"
        oas_rcv_meta(6)%clpname = "ICOSHFLX"
        oas_rcv_meta(7)%clpname = "ICOLHFLX"
        oas_rcv_meta(8)%clpname = "ICOTGRND"
      END IF

      ! Register send variables with OASIS
      DO jg = 1, SIZE(oas_snd_meta)
        CALL oasis_def_var(oas_snd_meta(jg)%vid, oas_snd_meta(jg)%clpname, oas_part_id, &
          oas_var_nodims, OASIS_Out, oas_vshape, OASIS_Real, oas_error)
        IF (oas_error /= 0) THEN
          WRITE(oas_message,*) 'Failure in oasis_def_var for ', oas_snd_meta(jg)%clpname
          CALL message(routine, oas_message)
          CALL oasis_abort(oas_comp_id, oas_comp_name, '')
        END IF
      END DO

      ! Register receive variables with OASIS
      DO jg = 1, SIZE(oas_rcv_meta)
        CALL oasis_def_var(oas_rcv_meta(jg)%vid, oas_rcv_meta(jg)%clpname, oas_part_id, &
          oas_var_nodims, OASIS_In, oas_vshape, OASIS_Real, oas_error)
        IF (oas_error /= 0) THEN
          WRITE(oas_message,*) 'Failure in oasis_def_var for ', oas_rcv_meta(jg)%clpname
          CALL message(routine, oas_message)
          CALL oasis_abort(oas_comp_id, oas_comp_name, '')
        END IF
      END DO

      IF(msg_level >= 30 ) THEN
        WRITE(oas_prt,*) 'ICONOAS: ',TRIM(routine),' oasis_def_partition completed'
        FLUSH(oas_prt)
      END IF

      ! Allocate memory for data exchange
      ALLOCATE( oas_snd_field(oas_vshape(1):oas_vshape(2),11), stat=oas_error )
      IF (oas_error > 0) CALL oasis_abort(oas_comp_id, oas_comp_name, &
          'Failure in allocating ICON send buffers')

      ALLOCATE( oas_rcv_field(oas_vshape(1):oas_vshape(2),oas_rcv_nvar), stat=oas_error )
      IF (oas_error > 0) CALL oasis_abort(oas_comp_id, oas_comp_name, &
          'Failure in allocating ICON receive buffers')

      ALLOCATE( oas_rcv_field_icon(nproma, p_patch(1)%nblks_c,oas_rcv_nvar), stat=oas_error )
      IF (oas_error > 0) CALL oasis_abort(oas_comp_id, oas_comp_name, &
          'Failure in allocating ICON receive buffer arrays')

      ! Allocate and initialize coupling mask
      ALLOCATE( loas_clm_c(nproma, p_patch(1)%nblks_c), stat=oas_error )
      IF (oas_error > 0) CALL oasis_abort(oas_comp_id, oas_comp_name, &
          'Failure in allocating ICON coupling mask')
      loas_clm_c(:,:) = .TRUE.

      ! Create the OASIS coupling mask
      CALL oasis_create_mask(p_patch(1), ext_data(1))

      IF(msg_level >= 20 ) THEN
        WRITE(oas_message,*) 'ICONOAS: Coupling mask created'
        CALL message(routine, oas_message)
      END IF

      ! Finalize OASIS definitions
      CALL oasis_enddef(oas_error)

      IF(msg_level >= 30 ) THEN
        WRITE(oas_prt,*) 'ICONOAS: ',TRIM(routine),' oasis_enddef completed on coupling processor'
        FLUSH(oas_prt)
      END IF

    ELSE ! lpe_cpl

      ! For non-coupling PEs, still call oasis_enddef
      CALL oasis_enddef(oas_error)

      IF(msg_level >= 30 ) THEN
        WRITE(oas_prt,*) 'ICONOAS: ',TRIM(routine),' oasis_enddef completed on non-coupling processor'
        FLUSH(oas_prt)
      END IF

    END IF ! lpe_cpl

    ! Initialize buffers with dummy values
    oas_snd_field = -1000._wp
    oas_rcv_field = -1000._wp
    oas_rcv_field_icon = -1000._wp

    IF(msg_level >= 30 ) THEN
      WRITE(oas_prt,*) 'ICONOAS: ',TRIM(routine),' OASIS coupling initialization completed'
      FLUSH(oas_prt)
    END IF

  END SUBROUTINE construct_oasis_coupling

  !-----------------------------------------------------------------------------

  SUBROUTINE oasis_create_mask(p_patch, ext_data)
    !
    ! Creates OASIS mask for land-atmosphere coupling
    ! Sets loas_clm_c based on land and lake information
    !
    ! Arguments
    TYPE(t_patch), INTENT(IN) :: p_patch   !< Unstructured triangular grid and patch info
    TYPE(t_external_data), INTENT(IN) :: ext_data  !< External data (land-sea mask, lake mask)

    ! Local variables
    INTEGER :: jb, jc, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    CHARACTER(len=1024) :: oas_message
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':oasis_create_mask'

    IF(msg_level >= 30 ) then
      WRITE(oas_message,*) 'iconoas setup oasis mask'
      CALL message(routine, oas_message)
    END IF

    ! set oasis mask
    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start, 1)
    i_endblk   = p_patch%cells%end_blk(rl_end, MAX(1,p_patch%n_childdom))
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx,&
        i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        loas_clm_c(jc,jb) = ext_data%atm%llsm_atm_c(jc,jb) .OR. ext_data%atm%llake_c(jc,jb)
!        loas_clm_c(jc,jb) = ext_data%atm%llsm_atm_c(jc,jb)
      END DO
    END DO

  END SUBROUTINE oasis_create_mask

  !-----------------------------------------------------------------------------

  SUBROUTINE oasis_send_fields(sim_time_oas, p_patch, p_nh_state, prm_diag)
    !
    ! Sends OASIS fields to the coupler
    !
    ! Arguments
    INTEGER, INTENT(IN) :: sim_time_oas  !< OASIS simulation time (0 for init, sim_time_oas for periodic)
    TYPE(t_patch), INTENT(IN) :: p_patch   !< Unstructured triangular grid and patch info
    TYPE(t_nh_state), INTENT(IN) :: p_nh_state   !< Non-hydrostatic state
    TYPE(t_nwp_phy_diag), INTENT(IN) :: prm_diag   !< Parametrized physics diagnostics

    ! Local variables
    INTEGER :: nlev, ii, jb, jc, isnd, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    CHARACTER(len=1024) :: oas_message
    INTEGER :: oas_prt = 6
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':oasis_send_fields'

    ! Initialize nlev
    nlev = p_patch%nlev

    IF(msg_level >= 30 ) then
      WRITE(oas_prt,*) 'iconoas: ',TRIM(routine),' oasis send fields'
      FLUSH(oas_prt)
    END IF

    oas_snd_field(:,7) = 0._wp
    oas_snd_field(:,8) = 0._wp

    rl_start = grf_bdywidth_c+1
    rl_end   = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start, 1)
    i_endblk   = p_patch%cells%end_blk(rl_end, MAX(1,p_patch%n_childdom))
    ii = 0
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx,&
        i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        ii = ii+1
!        WRITE(oas_message,*) 'DEBUG SPo swflx_dn_sfc_diff  ',ii, ' at 0: ', MAX(0._wp, prm_diag(1)%swflxsfc(jc,jb) - prm_diag(1)%swflx_dn_sfc_diff(jc,jb))
!        CALL message(routine, oas_message)
        oas_snd_field(ii,1)  = p_nh_state%diag%temp(jc,nlev,jb)
        oas_snd_field(ii,2)  = p_nh_state%diag%u(jc,nlev,jb)
        oas_snd_field(ii,3)  = p_nh_state%diag%v(jc,nlev,jb)
        oas_snd_field(ii,4)  = p_nh_state%prog(nnow(1))%tracer(jc,nlev,jb,iqv)
        oas_snd_field(ii,5)  = &
          p_nh_state%metrics%z_mc(jc,nlev,jb) - p_nh_state%metrics%z_ifc(jc,nlev+1,jb)
        oas_snd_field(ii,6)  = p_nh_state%diag%pres_sfc(jc,jb)
        oas_snd_field(ii,7)  = MAX(0._wp, prm_diag%swflxsfc(jc,jb) - prm_diag%swflx_dn_sfc_diff(jc,jb))
        oas_snd_field(ii,8)  = prm_diag%swflx_dn_sfc_diff(jc,jb)
        oas_snd_field(ii,9)  = prm_diag%lwflxsfc(jc,jb) + prm_diag%lwflx_up_sfc(jc,jb)
        oas_snd_field(ii,10) = prm_diag%rain_con_rate(jc,jb) + prm_diag%rain_gsp_rate(jc,jb)
        oas_snd_field(ii,11) = prm_diag%snow_con_rate(jc,jb) + prm_diag%snow_gsp_rate(jc,jb)! + &
        !  prm_diag(1)%ice_gsp_rate(jc,jb) + prm_diag(1)%graupel_gsp_rate(jc,jb) + &
        ! prm_diag(1)%hail_gsp_rate(jc,jb)
      END DO
    END DO

    ! Send and log fields
    DO isnd = 1, SIZE(oas_snd_meta)
      ! Send fields to OASIS
      CALL oasis_put(oas_snd_meta(isnd)%vid, sim_time_oas, oas_snd_field(:,isnd), oas_error)
      IF(msg_level >= 14 ) THEN
        WRITE(oas_message,*) 'ICONOAS: Sending  ', oas_snd_meta(isnd)%clpname, ' at ',sim_time_oas,' ', &
                MINVAL(oas_snd_field(:,isnd)), MAXVAL(oas_snd_field(:,isnd))
        CALL message(routine, oas_message)
      END IF
      IF (oas_error .NE. OASIS_Ok .AND. oas_error .LT. OASIS_Sent) THEN
        WRITE(oas_message,*) 'Failure in oasis_put of ', oas_snd_meta(isnd)%clpname
        CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)
      END IF
    END DO

  END SUBROUTINE oasis_send_fields

  !-----------------------------------------------------------------------------

!  SUBROUTINE oasis_receive_fields(sim_time_oas, p_patch, p_lnd_state)
   SUBROUTINE oasis_receive_fields(sim_time_oas, p_patch)
    !
    ! Receives OASIS fields from the coupler and transforms them into ICON format
    !
    ! Arguments
    INTEGER, INTENT(IN) :: sim_time_oas  !< OASIS simulation time (0 for init, sim_time_oas for periodic)
    TYPE(t_patch), INTENT(IN) :: p_patch   !< Unstructured triangular grid and patch info
!    TYPE(t_lnd_state), INTENT(INOUT) :: p_lnd_state  !< Land state (receives data from coupler)
 
    ! Local variables
    INTEGER :: ii, jb, jc, oas_i, rl_start, rl_end
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
    CHARACTER(len=1024) :: oas_message
    INTEGER :: oas_prt = 6
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':oasis_receive_fields'

    IF(msg_level >= 30 ) then
      WRITE(oas_prt,*) 'iconoas: ',TRIM(routine),' oasis receive fields'
      FLUSH(oas_prt)
    END IF

    DO oas_i = 1, SIZE(oas_rcv_meta)
      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) 'Receiving  ', oas_rcv_meta(oas_i)%clpname,&
          " at t=", sim_time_oas
        CALL message(routine, oas_message)
      END IF

      ! Receive field from OASIS
      CALL oasis_get(oas_rcv_meta(oas_i)%vid, sim_time_oas, oas_rcv_field(:,oas_i), oas_error)

      ! Logfiles min and max of received field
      IF(msg_level >= 14 ) THEN
        WRITE(oas_message,*) 'ICONOAS: Received  ', oas_rcv_meta(oas_i)%clpname, ' at t=', sim_time_oas,' ',&
                MINVAL(oas_rcv_field(:,oas_i)), MAXVAL(oas_rcv_field(:,oas_i))
        CALL message(routine, oas_message)
      END IF

      IF (oas_error .NE. OASIS_Ok .AND. oas_error .LT. OASIS_Recvd) THEN
        WRITE(oas_message,*) 'Failure in oasis_get of ', oas_rcv_meta(oas_i)%clpname
        CALL oasis_abort(oas_comp_id, oas_comp_name, oas_message)
      END IF

      ! transform thru oasis recieved variable into something icon can use
      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) "ICONOAS: icon transforming received vars"
        CALL message(routine, oas_message)
      END IF

      rl_start = grf_bdywidth_c+1
      rl_end   = min_rlcell_int
      i_startblk = p_patch%cells%start_blk(rl_start, 1)
      i_endblk   = p_patch%cells%end_blk(rl_end, MAX(1,p_patch%n_childdom))
      ii = 0
      DO jb = i_startblk, i_endblk
        CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, i_startidx,&
          i_endidx, rl_start, rl_end)
        DO jc = i_startidx, i_endidx
!          ii = idx_1d(jc,jb)
          ii = ii + 1
          oas_rcv_field_icon(jc,jb,oas_i) = oas_rcv_field(ii,oas_i)
        END DO
      END DO

      IF(msg_level > 30 ) THEN
        WRITE(oas_message,*) "ICONOAS: icon transformed received vars"
        CALL message(routine, oas_message)

        ! check:
        WRITE(oas_message,*) "ICONOAS: for ", oas_i, " got rank0-local min, max=", &
          MINVAL(oas_rcv_field(:,oas_i)), MAXVAL(oas_rcv_field(:,oas_i))
        CALL message(routine, oas_message)
      END IF
    END DO

  END SUBROUTINE oasis_receive_fields

  !-----------------------------------------------------------------------------

END MODULE mo_oasis_land_coupling

#endif

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

! Unit test to ensure the following functions from `mo_slice_array` work as expected:
!   get_3d_general
!   get_2d_general

PROGRAM test_slice_array

  USE mo_kind, ONLY: wp, vp, dp, sp
  USE mo_slice_array, ONLY: get_2d_general, get_3d_general
  USE mo_test_common, ONLY: test_fail, test_pass

  IMPLICIT NONE

  LOGICAL, PARAMETER :: debug = .FALSE. ! Not yet used
  CHARACTER(*), PARAMETER :: modname = 'test_slice_array'

  REAL(dp), TARGET :: r_data(5, 4, 3, 2, 1)
  REAL(sp), TARGET :: s_data(5, 4, 3, 2, 1)
  INTEGER, TARGET :: i_data(5, 4, 3, 2, 1)
  REAL(dp), POINTER :: r_ptr(:, :, :, :, :)
  REAL(sp), POINTER :: s_ptr(:, :, :, :, :)
  INTEGER, POINTER :: i_ptr(:, :, :, :, :)
  INTEGER, PARAMETER :: data_shape(5) = (/5, 4, 3, 2, 1/)
  INTEGER, PARAMETER :: data_size = PRODUCT(data_shape)
  INTEGER :: init_data(data_size)
  INTEGER :: i

  DO i = 1, SIZE(init_data)
    init_data(i) = i
  END DO

  r_data = RESHAPE(REAL(init_data, KIND=dp), SHAPE=data_shape)
  s_data = RESHAPE(REAL(init_data + data_size, sp), SHAPE=data_shape)
  i_data = RESHAPE(init_data + 2*data_size, SHAPE=data_shape)

  r_ptr => r_data
  s_ptr => s_data
  i_ptr => i_data

  IF (debug) THEN
    WRITE (*, *) "data_shape: ", data_shape
    !WRITE (*,*) "data: ", init_data
    !WRITE (*,*) "r_ptr: ", r_data
  END IF

  CALL check_get_3d_general(r_ptr, s_ptr, i_ptr)
  CALL check_get_2d_general(r_ptr, s_ptr, i_ptr)

  CALL test_pass

CONTAINS
  SUBROUTINE check_get_3d_general(r_data, s_data, i_data)
    REAL(dp), POINTER, INTENT(IN) :: r_data(:, :, :, :, :)
    REAL(sp), POINTER, INTENT(IN) :: s_data(:, :, :, :, :)
    INTEGER, POINTER, INTENT(IN) :: i_data(:, :, :, :, :)

    ! Local vars
    REAL(dp), POINTER :: r_ptr_3d(:, :, :)
    REAL(sp), POINTER :: s_ptr_3d(:, :, :)
    INTEGER, POINTER :: i_ptr_3d(:, :, :)

    REAL(wp), POINTER :: wp_ptr_3d(:, :, :)
    REAL(vp), POINTER :: vp_ptr_3d(:, :, :)

    r_ptr_3d => get_3d_general(r_data, (/3, 2/), (/.TRUE., .FALSE., .TRUE., .FALSE., .FALSE./))
    IF (ANY(r_data(3, :, 2, :, :) /= r_ptr_3d)) THEN
      CALL test_fail(modname, "get_3d_general(r_ptr)")
    END IF

    s_ptr_3d => get_3d_general(s_data, (/3, 2/), (/.TRUE., .FALSE., .TRUE., .FALSE., .FALSE./))
    IF (ANY(s_data(3, :, 2, :, :) /= s_ptr_3d)) THEN
      CALL test_fail(modname, "get_3d_general(s_ptr)")
    END IF

    i_ptr_3d => get_3d_general(i_data, (/3, 2/), (/.TRUE., .FALSE., .TRUE., .FALSE., .FALSE./))
    IF (ANY(i_data(3, :, 2, :, :) /= i_ptr_3d)) THEN
      CALL test_fail(modname, "get_3d_general(i_ptr)")
    END IF

#ifdef __SINGLE_PRECISION
    wp_ptr_3d => get_3d_general(s_data, (/3, 2/), (/.TRUE., .FALSE., .TRUE., .FALSE., .FALSE./))
    IF (ANY(s_data(3, :, 2, :, :) /= wp_ptr_3d)) THEN
#else
      wp_ptr_3d => get_3d_general(r_data, (/3, 2/), (/.TRUE., .FALSE., .TRUE., .FALSE., .FALSE./))
      IF (ANY(r_data(3, :, 2, :, :) /= wp_ptr_3d)) THEN
#endif
        CALL test_fail(modname, "get_3d_general(wp_ptr)")
      END IF

#if defined __SINGLE_PRECISION || defined __MIXED_PRECISION
      vp_ptr_3d => get_3d_general(s_data, (/3, 2/), (/.TRUE., .FALSE., .TRUE., .FALSE., .FALSE./))
      IF (ANY(s_data(3, :, 2, :, :) /= vp_ptr_3d)) THEN
#else
        vp_ptr_3d => get_3d_general(r_data, (/3, 2/), (/.TRUE., .FALSE., .TRUE., .FALSE., .FALSE./))
        IF (ANY(r_data(3, :, 2, :, :) /= vp_ptr_3d)) THEN
#endif
          CALL test_fail(modname, "get_3d_general(vp_ptr)")
        END IF

        END SUBROUTINE check_get_3d_general

        SUBROUTINE check_get_2d_general(r_data, s_data, i_data)
          REAL(dp), POINTER, INTENT(IN) :: r_data(:, :, :, :, :)
          REAL(sp), POINTER, INTENT(IN) :: s_data(:, :, :, :, :)
          INTEGER, POINTER, INTENT(IN) :: i_data(:, :, :, :, :)

          ! Local vars
          REAL(dp), POINTER :: r_ptr_2d(:, :)
          REAL(sp), POINTER :: s_ptr_2d(:, :)
          INTEGER, POINTER :: i_ptr_2d(:, :)

          REAL(wp), POINTER :: wp_ptr_2d(:, :)
          REAL(vp), POINTER :: vp_ptr_2d(:, :)

          r_ptr_2d => get_2d_general(r_data, (/3, 2, 2/), (/.TRUE., .FALSE., .TRUE., .TRUE., .FALSE./))
          IF (ANY(r_data(3, :, 2, 2, :) /= r_ptr_2d)) THEN
            CALL test_fail(modname, "get_2d_general(r_ptr)")
          END IF

          s_ptr_2d => get_2d_general(s_data, (/3, 2, 2/), (/.TRUE., .FALSE., .TRUE., .TRUE., .FALSE./))
          IF (ANY(s_data(3, :, 2, 2, :) /= s_ptr_2d)) THEN
            CALL test_fail(modname, "get_2d_general(s_ptr)")
          END IF

          i_ptr_2d => get_2d_general(i_data, (/3, 2, 2/), (/.TRUE., .FALSE., .TRUE., .TRUE., .FALSE./))
          IF (ANY(i_data(3, :, 2, 2, :) /= i_ptr_2d)) THEN
            CALL test_fail(modname, "get_2d_general(i_ptr)")
          END IF

#ifdef __SINGLE_PRECISION
          wp_ptr_2d => get_2d_general(s_data, (/3, 2, 2/), (/.TRUE., .FALSE., .TRUE., .TRUE., .FALSE./))
          IF (ANY(s_data(3, :, 2, 2, :) /= wp_ptr_2d)) THEN
#else
            wp_ptr_2d => get_2d_general(r_data, (/3, 2, 2/), (/.TRUE., .FALSE., .TRUE., .TRUE., .FALSE./))
            IF (ANY(r_data(3, :, 2, 2, :) /= wp_ptr_2d)) THEN
#endif
              CALL test_fail(modname, "get_2d_general(wp_ptr)")
            END IF

#if defined __SINGLE_PRECISION || defined __MIXED_PRECISION
            vp_ptr_2d => get_2d_general(s_data, (/3, 2, 2/), (/.TRUE., .FALSE., .TRUE., .TRUE., .FALSE./))
            IF (ANY(s_data(3, :, 2, 2, :) /= vp_ptr_2d)) THEN
#else
              vp_ptr_2d => get_2d_general(r_data, (/3, 2, 2/), (/.TRUE., .FALSE., .TRUE., .TRUE., .FALSE./))
              IF (ANY(r_data(3, :, 2, 2, :) /= vp_ptr_2d)) THEN
#endif
                CALL test_fail(modname, "get_2d_general(vp_ptr)")
              END IF

              END SUBROUTINE check_get_2d_general

              END PROGRAM test_slice_array

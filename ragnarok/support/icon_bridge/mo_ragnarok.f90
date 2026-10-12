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

MODULE mo_ragnarok
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_ragnarok, finalize_ragnarok, getKokkosVersion

  INTERFACE

    SUBROUTINE init_ragnarok() BIND(c, name="init_ragnarok")
    END SUBROUTINE init_ragnarok


    SUBROUTINE finalize_ragnarok() BIND(c, name="finalize_ragnarok")
    END SUBROUTINE finalize_ragnarok


    FUNCTION get_kokkos_version_c(len) RESULT(ret) BIND(c, name="retrieve_kokkos_version_c")
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(c_int), INTENT(INOUT) :: len
      TYPE(c_ptr) :: ret
    END FUNCTION get_kokkos_version_c

  END INTERFACE

  PUBLIC :: get_kokkos_version_c

CONTAINS

  SUBROUTINE getKokkosVersion(str)
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE
    CHARACTER(kind=c_char), POINTER, INTENT(OUT) :: str(:)

    TYPE(c_ptr) :: str_c
    INTEGER(c_int) :: len

    str_c = get_kokkos_version_c(len)
    CALL C_F_POINTER(str_c, str, (/len/))
  END SUBROUTINE getKokkosVersion

END MODULE mo_ragnarok

#ifdef COUP_OAS_ECLM
MODULE mo_oasis_define

USE mo_kind,  ONLY : wp
USE netcdf

USE mod_oasis

IMPLICIT NONE

INTEGER               :: oas_comp_id
CHARACTER(len=4)      :: oas_comp_name="icon"
INTEGER               :: kl_comm              ! local communicator 
INTEGER               :: oas_error            ! return error code
LOGICAL               :: lpe_cpl = .FALSE.
LOGICAL               :: loascplscheme_exchcoef = .FALSE.
INTEGER               :: oas_nlat
INTEGER               :: oas_nlon
INTEGER               :: oas_var_nodims(2)
INTEGER               :: oas_vshape(2)
INTEGER               :: oas_part_id
INTEGER               :: oas_rcv_nvar
INTEGER, ALLOCATABLE  :: oas_part(:)
TYPE :: t_oas_field
  CHARACTER(len = 8)  :: clpname
  INTEGER             :: vid
END TYPE t_oas_field
TYPE(t_oas_field), DIMENSION(12)  :: oas_snd_meta
TYPE(t_oas_field), ALLOCATABLE :: oas_rcv_meta(:)
REAL, ALLOCATABLE     :: oas_snd_field(:,:), oas_rcv_field(:,:)
REAL(wp), ALLOCATABLE :: oas_rcv_field_icon(:,:,:)
LOGICAL, ALLOCATABLE  :: loas_clm_c(:,:)      ! coupling mask
! SPo: Albedo (rcv vars 2, 3) and downward shortwave (snd vars 7, 8) are exchanged as
! bundled 3D fields carrying both radiation bands. Band 1 stays in the variable's own
! slot; band 2 goes to an appended slot (oas_nvar+2).
INTEGER, PARAMETER    :: oas_nbnd = 2       ! 1 = visible, 2 = near-infrared
REAL, ALLOCATABLE     :: oas_field_bnd(:,:) ! staging buffer for bundled put/get

PUBLIC :: &
 oas_comp_id, oas_comp_name, kl_comm, oas_error, oas_nlat, &
 oas_nlon, oas_var_nodims, oas_vshape, oas_part, &
 oas_snd_field, oas_rcv_field, oas_snd_meta, oas_rcv_meta, &
 loas_clm_c, loascplscheme_exchcoef

END MODULE mo_oasis_define
#endif

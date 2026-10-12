! Contains: Three routines of EVP dynamics. The driving routine is EVPdynamics.
! 2D indices are used to distinguish between boundary and internal nodes.
!
! Vladimir: added omp parallelization; rhs_a, rhs_m, rhs_mis are precalculated now
! The "original" AWI version of this solver is available in mo_ice_fem_evp_old!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2026, AWI
!
! Authors: Sergey Danilov, Qiang Wang
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
!
! This file has been modified for the use in ICON
! ---------------------------------------------------------------

!----------------------------
#include "omp_definitions.inc"
!----------------------------
!============================================================================
module mo_ice_fem_evp
  !
  ! Arrays defined here are used to keep mesh information
  !
  use mo_ice_fem_mesh
  use mo_ice_fem_types

  USE mo_kind,    ONLY: wp
#ifdef _OPENACC
    USE openacc, ONLY : acc_is_present
#endif
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PUBLIC :: init_evp_solver_coeffs
  PUBLIC :: EVPdynamics

  PRIVATE :: stress_tensor
  PRIVATE :: stress2rhs
  PRIVATE :: precalc4rhs
  PRIVATE :: index_si_elements

  PRIVATE
  ! some aggregated parameters used in the EVP solver
  REAL(wp):: val3=1.0_wp/3.0_wp
  REAL(wp):: vale, dte, det1, det2
  REAL(wp):: ax, ay

CONTAINS

!===================================================================

subroutine index_si_elements(lacc)
! Replaces "if" checking of whether sea ice is actually present in a given cell/node
#ifdef _OPENACC
    USE openacc, ONLY : acc_is_present
#endif
    IMPLICIT NONE

    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: elem,row, elnodes(3), i
    REAL(wp):: aa
    LOGICAL :: lzacc
    ! Temporary variables/buffers
    INTEGER :: buffy_array(myDim_elem2D)

    CALL set_acc_host_or_device(lzacc, lacc)

    ! count nodes with ice
    buffy_array(:) = 0._wp
    si_nod2D = 0

    !$ACC UPDATE SELF(m_ice, a_ice) IF(lzacc)
    DO i=1, myDim_nod2D
        row=myList_nod2D(i)

         aa = m_ice(row)*a_ice(row)

         IF (aa > 0._wp) THEN
            si_nod2D = si_nod2D + 1
            buffy_array(si_nod2D)=row
         ENDIF
    ENDDO

    if (allocated(si_idx_nodes)) then
      !$ACC EXIT DATA DELETE(si_idx_nodes) IF(lzacc .and. acc_is_present(si_idx_nodes))
      deallocate(si_idx_nodes)
    end if
    allocate(si_idx_nodes(si_nod2D))
    si_idx_nodes=buffy_array(1:si_nod2D)
    !$ACC ENTER DATA COPYIN(si_idx_nodes) IF(lzacc)

    ! count elements with ice
    buffy_array = 0._wp
    si_elem2D = 0

    DO i=1,myDim_elem2D
         elem=myList_elem2D(i)
         elnodes=elem2D_nodes(:,elem)

         aa=product(m_ice(elnodes))*product(a_ice(elnodes))

         IF (aa > 0._wp) THEN
            si_elem2D = si_elem2D + 1
            buffy_array(si_elem2D)=elem
         ENDIF
    ENDDO

    if (allocated(si_idx_elem)) then
      !$ACC EXIT DATA DELETE(si_idx_elem) IF(lzacc .and. acc_is_present(si_idx_elem))
      deallocate(si_idx_elem)
    end if
    allocate(si_idx_elem(si_elem2D))
    si_idx_elem=buffy_array(1:si_elem2D)
    !$ACC ENTER DATA COPYIN(si_idx_elem) IF(lzacc)

end subroutine index_si_elements
!===================================================================

subroutine precalc4rhs(lacc)
! Some of the quantities used in the solver do not change
! with the subcycles. Hence, can be precalculated and stored.
! Those are rhs_a, rhs_m, mass

  use mo_physical_constants,  ONLY: rhoi, rhos, grav

IMPLICIT NONE

LOGICAL, INTENT(IN), OPTIONAL :: lacc

INTEGER      :: row, elem, elnodes(3), nodels(6), k, i
REAL(wp) :: mass, aa
REAL(wp) :: cluster_area,elevation_elem(3),dx(3),dy(3)
REAL(wp) :: da(myDim_elem2D), dm(myDim_elem2D) ! temp storage for element-wise contributions from SSH
LOGICAL  :: lzacc

CALL set_acc_host_or_device(lzacc, lacc)

!$ACC DATA CREATE(da, dm) IF(lzacc)

!ICON_OMP_PARALLEL

!ICON_OMP_DO PRIVATE(i,row) SCHEDULE(static,4)
 !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
 DO i=1, myDim_nod2D
     row=myList_nod2D(i)
     rhs_a(row)=0.0_wp    ! these are used as temporal storage here
     rhs_m(row)=0.0_wp    ! for the contribution due to ssh
     rhs_mis(row)=0.0_wp

     rhs_u(row)=0.0_wp    ! these will be reinitialized at every subcycling iteration for non-ice-free nodes
     rhs_v(row)=0.0_wp    ! but fill all nodal vals with zeros here as well
 END DO
 !$ACC END PARALLEL LOOP
!ICON_OMP_END_DO

!ICON_OMP_DO PRIVATE(i,elem) SCHEDULE(static,4)
 !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
 DO i=1, myDim_elem2D
     elem=myList_elem2D(i)
     da(elem)=0.0_wp    ! initialize
     dm(elem)=0.0_wp    !
 END DO
 !$ACC END PARALLEL LOOP
!ICON_OMP_END_DO

!ICON_OMP_DO PRIVATE(i,elem,elnodes,aa,dx,dy,elevation_elem) SCHEDULE(static,4)
 !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(elnodes, dx, dy, elevation_elem) ASYNC(1) IF(lzacc)
 DO i=1,si_elem2D
     elem=si_idx_elem(i)
     elnodes(:)=elem2D_nodes(:,elem)

     dx(:)=bafux(:,elem)
     dy(:)=bafuy(:,elem)
     elevation_elem(:)=elevation(elnodes(:))

     ! use rhs_m and rhs_a for storing the contribution from elevation:
     aa=grav*voltriangle(elem)/3.0_wp
     da(elem)=-aa*sum(dx(:)*elevation_elem(:))
     dm(elem)=-aa*sum(dy(:)*elevation_elem(:))
 END DO
 !$ACC END PARALLEL LOOP
!ICON_OMP_END_DO

!ICON_OMP_DO  PRIVATE(i,row,nodels,k,elem)  SCHEDULE(static,4)
 !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) PRIVATE(nodels) ASYNC(1) IF(lzacc)
 DO i=1,si_nod2D
    row=si_idx_nodes(i)
    nodels=nod2D_elems(:,row)

     DO k=1,6
        elem=nodels(k)
        IF (elem > 0) THEN
         ! use rhs_m and rhs_a for storing the contribution from elevation:
        rhs_a(row)=rhs_a(row)+da(elem)
        rhs_m(row)=rhs_m(row)+dm(elem)
        END IF
     END DO
 END DO
 !$ACC END PARALLEL LOOP
!ICON_OMP_END_DO

!ICON_OMP_DO  PRIVATE(i,row,cluster_area,mass) SCHEDULE(static,4)
 !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
 DO i=1,si_nod2D
     row=si_idx_nodes(i)
     cluster_area=lmass_matrix(row)
     mass=cluster_area*(m_ice(row)*rhoi+m_snow(row)*rhos)

     rhs_a(row) = rhs_a(row)/cluster_area
     rhs_m(row) = rhs_m(row)/cluster_area
     rhs_mis(row) = mass
 END DO
 !$ACC END PARALLEL LOOP
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

!$ACC WAIT(1)
!$ACC END DATA

end subroutine precalc4rhs
!===================================================================

subroutine init_evp_solver_coeffs
! Calculates coefficients which are used for calculations in stress_tensor
! Called once during the initialization step at ice_init_fem

  USE mo_run_config,          ONLY: dtime
  USE mo_sea_ice_nml,         ONLY: evp_rheol_steps, Tevp_inv, theta_io, ellipse

  val3=1.0_wp/3.0_wp
  vale=1.0_wp/(ellipse**2)

  dte=dtime/(1.0_wp*REAL(evp_rheol_steps,wp))
  det1=1.0_wp+0.5_wp*Tevp_inv*dte
  det2=1.0_wp+0.5_wp*Tevp_inv*dte*ellipse**2    !RTSD corrected 8.3.2006
                                              ! There is error in CICE
                          ! manual.
  det1=1.0_wp/det1
  det2=1.0_wp/det2

  ! theta_io should be zero - set in ice_main.f90
  ax=cos(theta_io)
  ay=sin(theta_io)

end subroutine init_evp_solver_coeffs
!===================================================================
subroutine stress_tensor(si_elem2D, si_idx_elem, lacc)
!NEC$ always_inline
! EVP rheology implementation. Computes stress tensor components based on ice
! velocity field. They are stored as elemental arrays (sigma11, sigma22 and
! sigma12).

  USE mo_sea_ice_nml,         ONLY: delta_min, Tevp_inv, Pstar, c_pressure, luse_replacement_pressure, &
                                    seaice_stress_diag

implicit none

    INTEGER, INTENT(IN)               :: si_elem2D
    INTEGER, INTENT(IN), DIMENSION(:) :: si_idx_elem
    LOGICAL, INTENT(IN), OPTIONAL     :: lacc

    REAL(wp)   :: pressure, P, delta_inv!, aa
    INTEGER    :: elnodes(3), elem, i
    REAL(wp)   :: asum, msum, dx(3), dy(3)
    REAL(wp)   :: r1, r2, r3
    REAL(wp)   :: zeta, usum, vsum
    LOGICAL    :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA CREATE(elnodes, dx, dy) IF(lzacc)

! !ICON_OMP_DO        PRIVATE(i,elem,elnodes,dx,dy,vsum,usum,msum,asum,pressure, &
! !ICON_OMP                   delta_inv,zeta,r1,r2,r3)  SCHEDULE(static,4)
!NEC$ ivdep
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    DO i=1,si_elem2D
     elem = si_idx_elem(i)

! ATTENTION: the rows commented with !metrics contain terms due to
! differentiation of metrics.

     elnodes=elem2D_nodes(:,elem)

     dx=bafux(:,elem)
     dy=bafuy(:,elem)

!    meancos=sin_elem2D(elem)/cos_elem2D(elem)/earth_radius !is precalculated and stored as metrics_elem2D(elem)
     vsum=sum(v_ice(elnodes))                           !metrics
     usum=sum(u_ice(elnodes))                           !metrics

      ! ===== Deformation rate tensor on element elem:
     eps11(elem)=sum(dx*u_ice(elnodes))
     eps11(elem)=eps11(elem)-val3*vsum*metrics_elem2D(elem)                !metrics
     eps22(elem)=sum(dy*v_ice(elnodes))
     eps12(elem)=0.5_wp*sum(dy*u_ice(elnodes) + dx*v_ice(elnodes))
     eps12(elem)=eps12(elem)+0.5_wp*val3*usum*metrics_elem2D(elem)          !metrics
      ! ===== moduli:
     delta(elem)=(eps11(elem)**2+eps22(elem)**2)*(1.0_wp+vale)+4.0_wp*vale*eps12(elem)**2 + &
            2.0_wp*eps11(elem)*eps22(elem)*(1.0_wp-vale)
     delta(elem)=sqrt(delta(elem))
     msum=sum(m_ice(elnodes))*val3
     asum=sum(a_ice(elnodes))*val3

      ! ===== Hunke and Dukowicz c*h*p*
     pressure=pstar*(msum)*exp(-c_pressure*(1.0_wp-asum))
      ! =======================================
      ! ===== Here the EVP rheology piece starts
      ! =======================================
     pressure=0.5_wp*pressure
      ! ===== viscosity zeta should exceed zeta_min
      ! (done via limiting delta from above)
     !if(delta>pressure/zeta_min) delta=pressure/zeta_min
      ! ===== if viscosity is too big, it is limited too
      ! (done via correcting delta_inv)
     delta_inv=1.0_wp/max(delta(elem),delta_min)

     IF (luse_replacement_pressure) THEN
       !pressure=pressure*delta*delta_inv    ! Limiting pressure --- may not
       !                                     ! be needed. Should be tested
       P=pressure*delta(elem)/(delta(elem)+delta_min)
     ELSE
       P=pressure
     END IF

     IF (seaice_stress_diag) THEN
       ! ice strength / replacement pressure diagnostics
       p0_diag(elem) = pressure*2.0_wp ! compensate for 0.5*pressure above
       p_diag(elem) = p*2.0_wp         ! same for replacement pressure
     END IF

      ! ===== Limiting pressure/Delta  (zeta): still it may happen that zeta is too
      ! large in regions with fine mesh so that CFL criterion is violated.

      zeta=pressure*delta_inv
      ! This place was introduced by Hunke, but seemingly is not used
      ! in the current CICE. We artificially increase Clim_evp so
      ! that almost no limiting happens here
      !if (zeta>Clim_evp*voltriangle(elem)) then
      !zeta=Clim_evp*voltriangle(elem)
      !end if

      P=P*Tevp_inv
      zeta=zeta*Tevp_inv

     r1=zeta*(eps11(elem)+eps22(elem)) - P
     r2=zeta*(eps11(elem)-eps22(elem))
     r3=zeta*eps12(elem)
     si1(elem)=sigma11(elem)+sigma22(elem)
     si2(elem)=sigma11(elem)-sigma22(elem)

     si1(elem)=det1*(si1(elem)+dte*r1)
     si2(elem)=det2*(si2(elem)+dte*r2)
     sigma12(elem)=det2*(sigma12(elem)+dte*r3)
     sigma11(elem)=0.5_wp*(si1(elem)+si2(elem))
     sigma22(elem)=0.5_wp*(si1(elem)-si2(elem))

     ! Diagnostics only, not used in EVP time stepping
     IF (seaice_stress_diag) THEN
       ! =========================================
       ! Stress invariants (diagnostic)
       ! =========================================

       ! Normal stress invariant (mean normal stress)
       sigma_i(elem) = 0.5_wp * ( sigma11(elem) + sigma22(elem) )

       ! Shear stress invariant
       sigma_ii(elem) = sqrt( &
         (0.5_wp * (sigma11(elem) - sigma22(elem)))**2 &
         + sigma12(elem)**2 )

       ! ==========================================================
       ! Principal stresses (diagnostic)
       ! ==========================================================

       ! Maximum principal stress
       sigma_1(elem) = sigma_i(elem) + sigma_ii(elem)

       ! Minimum principal stress
       sigma_2(elem) = sigma_i(elem) - sigma_ii(elem)
     END IF

    END DO
    !$ACC END PARALLEL LOOP
! !ICON_OMP_END_DO

    !$ACC WAIT(1)
    !$ACC END DATA
end subroutine stress_tensor
!===================================================================
subroutine stress2rhs(si_nod2D, si_idx_nodes, lacc)
! EVP implementation:
! Computes the divergence of stress tensor and puts the result into the
! rhs vectors

IMPLICIT NONE

INTEGER, INTENT(IN) :: si_nod2D
INTEGER, INTENT(IN), DIMENSION(:) :: si_idx_nodes
LOGICAL, INTENT(IN), OPTIONAL     :: lacc

INTEGER  :: row, elem, k, i
REAL(wp) :: dx, dy
LOGICAL  :: lzacc

  CALL set_acc_host_or_device(lzacc, lacc)

! !ICON_OMP_DO        PRIVATE(i,k,row,elem,dx,dy) ICON_OMP_DEFAULT_SCHEDULE
!NEC$ ivdep
  !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
  DO i=1,si_nod2D

     row = si_idx_nodes(i)

     ! initialize
     rhs_u(row)=0.0_wp
     rhs_v(row)=0.0_wp

     DO k=1,6
        elem=nod2D_elems(k,row)
        dx=bafux_nod(k,row)
        dy=bafuy_nod(k,row)

        IF (elem > 0) THEN

        rhs_u(row)=rhs_u(row) - voltriangle(elem) * &
             (sigma11(elem)*dx+sigma12(elem)*dy &
             +sigma12(elem)*val3*metrics_elem2D(elem))                          !metrics
        rhs_v(row)=rhs_v(row) - voltriangle(elem) * &
             (sigma12(elem)*dx+sigma22(elem)*dy &
             -sigma11(elem)*val3*metrics_elem2D(elem))
        END IF
     END DO

     !  Vladimir: rhs_a and rhs_m are calculated in precalc4rhs
     rhs_u(row)=rhs_u(row)/rhs_mis(row) + rhs_a(row)
     rhs_v(row)=rhs_v(row)/rhs_mis(row) + rhs_m(row)

  END DO
  !$ACC END PARALLEL LOOP
! !ICON_OMP_END_DO

  !$ACC WAIT(1)
end subroutine stress2rhs
!===================================================================

subroutine EVPdynamics(lacc)
! EVP implementation. Does cybcycling and boundary conditions.
  USE mo_sea_ice_nml,           ONLY: evp_rheol_steps, Cd_io
  USE mo_physical_constants,    ONLY: rhoi, rhos, rho_ref

  USE mo_ice_fem_icon_init,     ONLY: exchange_nod2D

IMPLICIT NONE
LOGICAL, INTENT(IN), OPTIONAL :: lacc

integer     :: shortstep
REAL(wp)    ::  drag, inv_mass, det, umod, rhsu, rhsv
integer     ::  i,j
logical     :: lzacc

 CALL set_acc_host_or_device(lzacc, lacc)

 !$ACC DATA COPY(elem2D_nodes, bafux, bafuy, metrics_elem2D) &
 !$ACC   COPY(sigma11, sigma12, sigma22) &
 !$ACC   COPY(eps11, eps12, eps22) &
 !$ACC   COPY(si1, si2) &
 !$ACC   COPY(delta) &
 !$ACC   COPY(u_ice, v_ice, m_ice, a_ice, m_snow) &
 !$ACC   COPY(elevation, u_w, v_w, stress_atmice_x, stress_atmice_y) &
 !$ACC   COPY(myList_nod2D) &
 !$ACC   COPY(rhs_u, rhs_v, rhs_mis, rhs_a, rhs_m) &
 !$ACC   COPY(voltriangle, nod2D_elems, bafux_nod, bafuy_nod, index_nod2D) &
 !$ACC   COPY(coriolis_nod2D) &
 !$ACC   COPY(lmass_matrix, myList_elem2D) &
 !$ACC   IF(lzacc)

! We start evp iteration with stresses (sigma11,sigma12,sigma22) from the
! previous timestep for better convergence.

! index elements/nodes where sea ice is present for faster loops
    CALL index_si_elements(lacc=lzacc)
! precalculate several arrays that do not change during subcycling
    CALL precalc4rhs(lacc=lzacc)

 DO shortstep=1, evp_rheol_steps

     ! ===== Boundary conditions
     !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO j=1, myDim_nod2D+eDim_nod2D
        i=myList_nod2D(j)
        IF(index_nod2D(i)==1) THEN
          u_ice(i)=0.0_wp
          v_ice(i)=0.0_wp
        END IF
     END DO
     !$ACC END PARALLEL LOOP
     !$ACC WAIT(1)

     call stress_tensor(si_elem2D, si_idx_elem, lacc=lzacc)

     call stress2rhs(si_nod2D, si_idx_nodes, lacc=lzacc)

!ICON_OMP_PARALLEL
!ICON_OMP_DO        PRIVATE(j,i,inv_mass,umod,drag,rhsu,rhsv,det) ICON_OMP_DEFAULT_SCHEDULE
!NEC$ ivdep
     !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
     DO j=1,si_nod2D
       i=si_idx_nodes(j)
       IF (index_nod2D(i)>0) CYCLE          ! Skip boundary nodes
       IF (a_ice(i) > 0.01_wp) THEN             ! If ice is present, update velocities
         inv_mass=(rhoi*m_ice(i)+rhos*m_snow(i))/a_ice(i)
         inv_mass=max(inv_mass, 9.0_wp)        ! Limit the weighted mass
                                           ! if it is too small
         inv_mass=1.0_wp/inv_mass

         umod=sqrt((u_ice(i)-u_w(i))**2+(v_ice(i)-v_w(i))**2)
         drag=Cd_io*umod*rho_ref*inv_mass

         rhsu=u_ice(i)+dte*(drag*(ax*u_w(i)-ay*v_w(i))+inv_mass*stress_atmice_x(i)+rhs_u(i))
         rhsv=v_ice(i)+dte*(drag*(ax*v_w(i)+ay*u_w(i))+inv_mass*stress_atmice_y(i)+rhs_v(i))

         det=(1._wp+ax*drag*dte)**2+(dte*coriolis_nod2D(i)+dte*ay*drag)**2
         det=1.0_wp/det
         u_ice(i)=det*((1.0_wp+ax*drag*dte)*rhsu+dte*(coriolis_nod2D(i)+ay*drag)*rhsv)
         v_ice(i)=det*((1.0_wp+ax*drag*dte)*rhsv-dte*(coriolis_nod2D(i)+ay*drag)*rhsu)

       ELSE !Set velocity to 0 if ice absent
         u_ice(i) = 0.0_wp
         v_ice(i) = 0.0_wp
      ! better set ice velocity equal to water velocity
      ! u_ice(i)=u_w(i)
      ! v_ice(i)=v_w(i)

       END IF

     END DO
     !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

     CALL exchange_nod2D(u_ice, v_ice, lacc=lzacc)

 END DO

 !$ACC END DATA

end subroutine EVPdynamics
!===================================================================

end module mo_ice_fem_evp

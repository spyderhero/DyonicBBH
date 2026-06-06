! NPEM
! NP_calc_phi_gfs_EM_direct_magnetic_corrected.F90
!
! Computes electromagnetic Newman-Penrose scalars directly from the
! Einstein-Maxwell evolution variables E^i and B^i.
!
! This version is for the EM-direct formulation used in the present
! electric/magnetic-charge project. It does not reconstruct B^i from a
! vector potential A_i, and therefore is compatible with Dirac-monopole
! type magnetic charge initial data where a global smooth vector potential
! is not available.
!
! Evolution variables expected from ProcaBase / EinsteinMaxwell setup:
!   Ex,Ey,Ez  : contravariant electric field E^i
!   Bx,By,Bz  : contravariant magnetic field B^i
!
! Output grid functions:
!   phi0re, phi0im, phi1re, phi1im, phi2re, phi2im
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NPEM_calcPhiGF( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL gd(3,3), detgd
  CCTK_REAL Eu(3), Bu(3), Ed(3), Bd(3)
  CCTK_REAL xx(3)
  CCTK_REAL u_vec(3), v_vec(3), w_vec(3)
  CCTK_REAL dotp1, dotp2
  CCTK_REAL norm_u, norm_v, norm_w
  CCTK_REAL eps_axis

  CCTK_INT  i, j, k, m

  eps_axis = 1.0d-10

  !=== initialize NP scalar grid functions as zero ===
  if (calc_phi /= 0) then
     phi0re = 0.0d0
     phi0im = 0.0d0
     phi1re = 0.0d0
     phi1im = 0.0d0
     phi2re = 0.0d0
     phi2im = 0.0d0
  end if

  if (calc_phi == 0) return

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(i,j,k,m,gd,detgd,Eu,Bu,Ed,Bd,xx,u_vec,v_vec,w_vec, &
  !$OMP         dotp1,dotp2,norm_u,norm_v,norm_w)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !-------------- Get local metric ----------
    gd(1,1) = gxx(i,j,k)
    gd(1,2) = gxy(i,j,k)
    gd(1,3) = gxz(i,j,k)
    gd(2,2) = gyy(i,j,k)
    gd(2,3) = gyz(i,j,k)
    gd(3,3) = gzz(i,j,k)
    gd(2,1) = gd(1,2)
    gd(3,1) = gd(1,3)
    gd(3,2) = gd(2,3)

    detgd =       gd(1,1) * gd(2,2) * gd(3,3)                                &
            + 2.0d0 * gd(1,2) * gd(1,3) * gd(2,3)                            &
            -        gd(1,1) * gd(2,3) ** 2                                  &
            -        gd(2,2) * gd(1,3) ** 2                                  &
            -        gd(3,3) * gd(1,2) ** 2

    if (detgd <= 0.0d0) then
       phi0re(i,j,k) = 0.0d0
       phi0im(i,j,k) = 0.0d0
       phi1re(i,j,k) = 0.0d0
       phi1im(i,j,k) = 0.0d0
       phi2re(i,j,k) = 0.0d0
       phi2im(i,j,k) = 0.0d0
       cycle
    end if

    !-------------- EM variables --------------
    ! E^i and B^i are already evolved directly.  Do not compute B^i
    ! from curl(A), because magnetic monopole data generally do not admit
    ! a single smooth global vector potential.
    Eu(1) = Ex(i,j,k)
    Eu(2) = Ey(i,j,k)
    Eu(3) = Ez(i,j,k)

    Bu(1) = Bx(i,j,k)
    Bu(2) = By(i,j,k)
    Bu(3) = Bz(i,j,k)

    Ed(1) = gd(1,1) * Eu(1) + gd(1,2) * Eu(2) + gd(1,3) * Eu(3)
    Ed(2) = gd(2,1) * Eu(1) + gd(2,2) * Eu(2) + gd(2,3) * Eu(3)
    Ed(3) = gd(3,1) * Eu(1) + gd(3,2) * Eu(2) + gd(3,3) * Eu(3)

    Bd(1) = gd(1,1) * Bu(1) + gd(1,2) * Bu(2) + gd(1,3) * Bu(3)
    Bd(2) = gd(2,1) * Bu(1) + gd(2,2) * Bu(2) + gd(2,3) * Bu(3)
    Bd(3) = gd(3,1) * Bu(1) + gd(3,2) * Bu(2) + gd(3,3) * Bu(3)

    !------------ Orthonormal basis -----------
    ! u_vec: radial direction, v_vec/w_vec: angular directions.
    ! This follows the original NPEM construction and the Multipole-style
    ! z_orientation convention.
    xx(:) = (/ x(i,j,k), y(i,j,k), z(i,j,k) /)

    ! Points exactly on the z-axis have ill-defined angular basis vectors.
    if (xx(1)**2 + xx(2)**2 < 1.0d-12) xx(1) = xx(1) + eps_axis

    u_vec(:) = xx(:)
    v_vec(:) = (/ xx(1)*xx(3), xx(2)*xx(3), -xx(1)**2 - xx(2)**2 /)
    w_vec(:) = (/ -xx(2), xx(1), 0.0d0 /)

    ! Normalize u_vec
    norm_u = 0.0d0
    do m = 1, 3
       norm_u = norm_u + u_vec(m) * ( gd(m,1)*u_vec(1) + gd(m,2)*u_vec(2) + gd(m,3)*u_vec(3) )
    end do
    if (norm_u <= 0.0d0) cycle
    u_vec = u_vec / sqrt(norm_u)

    ! Orthogonalize and normalize v_vec
    dotp1 = 0.0d0
    do m = 1, 3
       dotp1 = dotp1 + u_vec(m) * ( gd(m,1)*v_vec(1) + gd(m,2)*v_vec(2) + gd(m,3)*v_vec(3) )
    end do
    v_vec = v_vec - dotp1 * u_vec

    norm_v = 0.0d0
    do m = 1, 3
       norm_v = norm_v + v_vec(m) * ( gd(m,1)*v_vec(1) + gd(m,2)*v_vec(2) + gd(m,3)*v_vec(3) )
    end do
    if (norm_v <= 0.0d0) cycle
    v_vec = v_vec / sqrt(norm_v)

    ! Orthogonalize and normalize w_vec
    dotp1 = 0.0d0
    dotp2 = 0.0d0
    do m = 1, 3
       dotp1 = dotp1 + u_vec(m) * ( gd(m,1)*w_vec(1) + gd(m,2)*w_vec(2) + gd(m,3)*w_vec(3) )
       dotp2 = dotp2 + v_vec(m) * ( gd(m,1)*w_vec(1) + gd(m,2)*w_vec(2) + gd(m,3)*w_vec(3) )
    end do
    w_vec = w_vec - dotp1 * u_vec - dotp2 * v_vec

    norm_w = 0.0d0
    do m = 1, 3
       norm_w = norm_w + w_vec(m) * ( gd(m,1)*w_vec(1) + gd(m,2)*w_vec(2) + gd(m,3)*w_vec(3) )
    end do
    if (norm_w <= 0.0d0) cycle
    w_vec = w_vec / sqrt(norm_w)

    !------------ Electromagnetic NP scalars ------------
    ! With the original tetrad convention used by this thorn:
    !   Phi0 = -1/2 [ E_v - B_w ] - i/2 [ E_w + B_v ]
    !   Phi2 = +1/2 [ E_v + B_w ] - i/2 [ E_w - B_v ]
    !   Phi1 = +1/2 [ E_u ]       + i/2 [ B_u ]
    ! This makes a magnetic monopole contribute directly to Im(Phi1),
    ! which is essential for the electric/magnetic-charge project.
    do m = 1, 3
       phi0re(i,j,k) = phi0re(i,j,k) - 0.5d0 * ( Ed(m) * v_vec(m) - Bd(m) * w_vec(m) )
       phi0im(i,j,k) = phi0im(i,j,k) - 0.5d0 * ( Ed(m) * w_vec(m) + Bd(m) * v_vec(m) )

       phi2re(i,j,k) = phi2re(i,j,k) + 0.5d0 * ( Ed(m) * v_vec(m) + Bd(m) * w_vec(m) )
       phi2im(i,j,k) = phi2im(i,j,k) - 0.5d0 * ( Ed(m) * w_vec(m) - Bd(m) * v_vec(m) )

       phi1re(i,j,k) = phi1re(i,j,k) + 0.5d0 * Ed(m) * u_vec(m)
       phi1im(i,j,k) = phi1im(i,j,k) + 0.5d0 * Bd(m) * u_vec(m)
    end do

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine NPEM_calcPhiGF

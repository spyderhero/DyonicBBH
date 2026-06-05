! Proca_adm_vs_bssn.F90 : Functions for converting between ADM and BSSN variables
!
! EM-direct version consistent with the Einstein-Maxwell evolution system
! based on arXiv:1205.1063.
!
! This file does not contain any electromagnetic evolution variables.  It only
! constructs the BSSN geometric variables used by the EM-direct RHS file:
!
!   chi, hxx, hxy, hxz, hyy, hyz, hzz, tracek
!
! The electromagnetic variables Ex,Ey,Ez,Bx,By,Bz,Psi,Phi are evolved in the
! corresponding RHS routine and are not modified here.
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Proca_adm2bssn( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL gg(3,3), gu(3,3), kk(3,3)
  CCTK_REAL detg, ch, hh(3,3), trk
  CCTK_INT  i, j, k, a, b

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(gg, gu, kk, detg, ch, hh, trk, i, j, k, a, b)
  do k = 1, cctk_lsh(3)
  do j = 1, cctk_lsh(2)
  do i = 1, cctk_lsh(1)

    !------------ Get local ADM variables -----------------
    gg(1,1) = gxx(i,j,k)
    gg(1,2) = gxy(i,j,k)
    gg(1,3) = gxz(i,j,k)
    gg(2,2) = gyy(i,j,k)
    gg(2,3) = gyz(i,j,k)
    gg(3,3) = gzz(i,j,k)
    gg(2,1) = gg(1,2)
    gg(3,1) = gg(1,3)
    gg(3,2) = gg(2,3)

    kk(1,1) = kxx(i,j,k)
    kk(1,2) = kxy(i,j,k)
    kk(1,3) = kxz(i,j,k)
    kk(2,2) = kyy(i,j,k)
    kk(2,3) = kyz(i,j,k)
    kk(3,3) = kzz(i,j,k)
    kk(2,1) = kk(1,2)
    kk(3,1) = kk(1,3)
    kk(3,2) = kk(2,3)
    !-------------------------------------------------------

    !------------- Calculate det(gamma_ij) ----------------
    detg =       gg(1,1) * gg(2,2) * gg(3,3)                         &
           + 2.0d0 * gg(1,2) * gg(1,3) * gg(2,3)                     &
           -       gg(1,1) * gg(2,3) ** 2                            &
           -       gg(2,2) * gg(1,3) ** 2                            &
           -       gg(3,3) * gg(1,2) ** 2

    if (detg <= 0.0d0) then
       call CCTK_WARN(0, "ADM 3-metric determinant is non-positive in Proca_adm2bssn")
    end if
    !-------------------------------------------------------

    !------------- Invert gamma_ij -------------------------
    gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / detg
    gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / detg
    gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / detg
    gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / detg
    gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / detg
    gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / detg
    gu(2,1) = gu(1,2)
    gu(3,1) = gu(1,3)
    gu(3,2) = gu(2,3)
    !-------------------------------------------------------

    !------------ Convert ADM geometry to BSSN geometry ----
    ! K = gamma^ij K_ij
    trk = 0.0d0
    do a = 1, 3
      do b = 1, 3
        trk = trk + gu(a,b) * kk(a,b)
      end do
    end do

    ! Convention used by the RHS routine:
    !   chi = det(gamma_ij)^[-1/(3*p)]
    !   h_ij = chi^p gamma_ij
    ! where p = conf_fac_exponent.  For p=1 this is the usual
    !   chi = gamma^(-1/3), h_ij = chi gamma_ij.
    ch = detg ** (-1.0d0 / (3.0d0 * conf_fac_exponent))
    if (ch < chi_floor) ch = chi_floor

    hh = ch ** conf_fac_exponent * gg
    !-------------------------------------------------------

    !------------ Write BSSN variables to grid functions ---
    chi(i,j,k) = ch

    hxx(i,j,k) = hh(1,1)
    hxy(i,j,k) = hh(1,2)
    hxz(i,j,k) = hh(1,3)
    hyy(i,j,k) = hh(2,2)
    hyz(i,j,k) = hh(2,3)
    hzz(i,j,k) = hh(3,3)

    tracek(i,j,k) = trk
    !-------------------------------------------------------

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine Proca_adm2bssn

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#define D1X4(q) ((-q(i+2,j,k) + 8.0d0*q(i+1,j,k) - 8.0d0*q(i-1,j,k) + q(i-2,j,k)) / dx12)
#define D1Y4(q) ((-q(i,j+2,k) + 8.0d0*q(i,j+1,k) - 8.0d0*q(i,j-1,k) + q(i,j-2,k)) / dy12)
#define D1Z4(q) ((-q(i,j,k+2) + 8.0d0*q(i,j,k+1) - 8.0d0*q(i,j,k-1) + q(i,j,k-2)) / dz12)
#define D1X6(q) ((q(i+3,j,k) - 9.0d0*q(i+2,j,k) + 45.0d0*q(i+1,j,k) - q(i-3,j,k) + 9.0d0*q(i-2,j,k) - 45.0d0*q(i-1,j,k)) * odx60)
#define D1Y6(q) ((q(i,j+3,k) - 9.0d0*q(i,j+2,k) + 45.0d0*q(i,j+1,k) - q(i,j-3,k) + 9.0d0*q(i,j-2,k) - 45.0d0*q(i,j-1,k)) * ody60)
#define D1Z6(q) ((q(i,j,k+3) - 9.0d0*q(i,j,k+2) + 45.0d0*q(i,j,k+1) - q(i,j,k-3) + 9.0d0*q(i,j,k-2) - 45.0d0*q(i,j,k-1)) * odz60)

subroutine EinsteinMaxwell_calc_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! This RHS implements the Einstein-Maxwell cleaning system of arXiv:1205.1063,
  ! Eq. (2.8), using the direct variables E^i, B^i, Psi, Phi.
  ! E^i and B^i are treated as contravariant spatial vectors.

  CCTK_REAL                alph, beta(3)
  CCTK_REAL                hh(3,3), hu(3,3), trk, dethh, ch
  CCTK_REAL                lE(3), lB(3), lPsi, lPhi

  CCTK_REAL                d1_alph(3), d1_beta(3,3)
  CCTK_REAL                d1_hh(3,3,3), d1_ch(3)
  CCTK_REAL                d1_lE(3,3), d1_lB(3,3), d1_lPsi(3), d1_lPhi(3)

  CCTK_REAL                ad1_lE(3), ad1_lB(3), ad1_lPsi, ad1_lPhi
  CCTK_REAL                d1_f(3)

  CCTK_REAL                cf1(3,3,3), cf2(3,3,3)
  CCTK_REAL                divE, divB
  CCTK_REAL                curlE, curlB

  CCTK_REAL                rhs_lE(3), rhs_lB(3), rhs_lPsi, rhs_lPhi

  CCTK_REAL                dx12, dy12, dz12
  CCTK_REAL                odx60, ody60, odz60
  CCTK_INT                 i, j, k
  CCTK_INT                 di, dj, dk
  CCTK_REAL, parameter ::  one = 1.0d0
  CCTK_INT                 a, b, c, m
  CCTK_REAL                levi_civita(3,3,3)

  integer                  istat
  logical                  use_jacobian
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ11, lJ12, lJ13
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ21, lJ22, lJ23
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: lJ31, lJ32, lJ33
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ111, ldJ112, ldJ113, ldJ122, ldJ123, ldJ133
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ211, ldJ212, ldJ213, ldJ222, ldJ223, ldJ233
  CCTK_REAL, dimension(cctk_ash(1),cctk_ash(2),cctk_ash(3)) :: ldJ311, ldJ312, ldJ313, ldJ322, ldJ323, ldJ333

  CCTK_POINTER             lJ11_ptr, lJ12_ptr, lJ13_ptr
  CCTK_POINTER             lJ21_ptr, lJ22_ptr, lJ23_ptr
  CCTK_POINTER             lJ31_ptr, lJ32_ptr, lJ33_ptr
  CCTK_POINTER             ldJ111_ptr, ldJ112_ptr, ldJ113_ptr, ldJ122_ptr, ldJ123_ptr, ldJ133_ptr
  CCTK_POINTER             ldJ211_ptr, ldJ212_ptr, ldJ213_ptr, ldJ222_ptr, ldJ223_ptr, ldJ233_ptr
  CCTK_POINTER             ldJ311_ptr, ldJ312_ptr, ldJ313_ptr, ldJ322_ptr, ldJ323_ptr, ldJ333_ptr

  CCTK_REAL                jac(3,3), hes(3,3,3)

  pointer (lJ11_ptr, lJ11), (lJ12_ptr, lJ12), (lJ13_ptr, lJ13)
  pointer (lJ21_ptr, lJ21), (lJ22_ptr, lJ22), (lJ23_ptr, lJ23)
  pointer (lJ31_ptr, lJ31), (lJ32_ptr, lJ32), (lJ33_ptr, lJ33)

  pointer (ldJ111_ptr, ldJ111), (ldJ112_ptr, ldJ112), (ldJ113_ptr, ldJ113), (ldJ122_ptr, ldJ122), (ldJ123_ptr, ldJ123), (ldJ133_ptr, ldJ133)
  pointer (ldJ211_ptr, ldJ211), (ldJ212_ptr, ldJ212), (ldJ213_ptr, ldJ213), (ldJ222_ptr, ldJ222), (ldJ223_ptr, ldJ223), (ldJ233_ptr, ldJ233)
  pointer (ldJ311_ptr, ldJ311), (ldJ312_ptr, ldJ312), (ldJ313_ptr, ldJ313), (ldJ322_ptr, ldJ322), (ldJ323_ptr, ldJ323), (ldJ333_ptr, ldJ333)

  call CCTK_IsFunctionAliased(istat, "MultiPatch_GetDomainSpecification")
  if (istat == 0) then
     use_jacobian = .false.
  else
     use_jacobian = .true.
  end if

  if (use_jacobian) then
     call CCTK_VarDataPtr(lJ11_ptr, cctkGH, 0, "Coordinates::J11")
     call CCTK_VarDataPtr(lJ12_ptr, cctkGH, 0, "Coordinates::J12")
     call CCTK_VarDataPtr(lJ13_ptr, cctkGH, 0, "Coordinates::J13")
     call CCTK_VarDataPtr(lJ21_ptr, cctkGH, 0, "Coordinates::J21")
     call CCTK_VarDataPtr(lJ22_ptr, cctkGH, 0, "Coordinates::J22")
     call CCTK_VarDataPtr(lJ23_ptr, cctkGH, 0, "Coordinates::J23")
     call CCTK_VarDataPtr(lJ31_ptr, cctkGH, 0, "Coordinates::J31")
     call CCTK_VarDataPtr(lJ32_ptr, cctkGH, 0, "Coordinates::J32")
     call CCTK_VarDataPtr(lJ33_ptr, cctkGH, 0, "Coordinates::J33")

     call CCTK_VarDataPtr(ldJ111_ptr, cctkGH, 0, "Coordinates::dJ111")
     call CCTK_VarDataPtr(ldJ112_ptr, cctkGH, 0, "Coordinates::dJ112")
     call CCTK_VarDataPtr(ldJ113_ptr, cctkGH, 0, "Coordinates::dJ113")
     call CCTK_VarDataPtr(ldJ122_ptr, cctkGH, 0, "Coordinates::dJ122")
     call CCTK_VarDataPtr(ldJ123_ptr, cctkGH, 0, "Coordinates::dJ123")
     call CCTK_VarDataPtr(ldJ133_ptr, cctkGH, 0, "Coordinates::dJ133")
     call CCTK_VarDataPtr(ldJ211_ptr, cctkGH, 0, "Coordinates::dJ211")
     call CCTK_VarDataPtr(ldJ212_ptr, cctkGH, 0, "Coordinates::dJ212")
     call CCTK_VarDataPtr(ldJ213_ptr, cctkGH, 0, "Coordinates::dJ213")
     call CCTK_VarDataPtr(ldJ222_ptr, cctkGH, 0, "Coordinates::dJ222")
     call CCTK_VarDataPtr(ldJ223_ptr, cctkGH, 0, "Coordinates::dJ223")
     call CCTK_VarDataPtr(ldJ233_ptr, cctkGH, 0, "Coordinates::dJ233")
     call CCTK_VarDataPtr(ldJ311_ptr, cctkGH, 0, "Coordinates::dJ311")
     call CCTK_VarDataPtr(ldJ312_ptr, cctkGH, 0, "Coordinates::dJ312")
     call CCTK_VarDataPtr(ldJ313_ptr, cctkGH, 0, "Coordinates::dJ313")
     call CCTK_VarDataPtr(ldJ322_ptr, cctkGH, 0, "Coordinates::dJ322")
     call CCTK_VarDataPtr(ldJ323_ptr, cctkGH, 0, "Coordinates::dJ323")
     call CCTK_VarDataPtr(ldJ333_ptr, cctkGH, 0, "Coordinates::dJ333")
  end if

  dx12 = 12.0d0*CCTK_DELTA_SPACE(1)
  dy12 = 12.0d0*CCTK_DELTA_SPACE(2)
  dz12 = 12.0d0*CCTK_DELTA_SPACE(3)

  odx60 = 1.0d0 / (60.0d0 * CCTK_DELTA_SPACE(1))
  ody60 = 1.0d0 / (60.0d0 * CCTK_DELTA_SPACE(2))
  odz60 = 1.0d0 / (60.0d0 * CCTK_DELTA_SPACE(3))

  call EinsteinMaxwell_adm2bssn(CCTK_PASS_FTOF)

  levi_civita = 0.0d0
  levi_civita(1,2,3) =  1.0d0
  levi_civita(2,3,1) =  1.0d0
  levi_civita(3,1,2) =  1.0d0
  levi_civita(1,3,2) = -1.0d0
  levi_civita(3,2,1) = -1.0d0
  levi_civita(2,1,3) = -1.0d0

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(alph, beta, hh, hu, trk, dethh, ch, &
  !$OMP lE, lB, lPsi, lPhi, &
  !$OMP d1_alph, d1_beta, d1_hh, d1_ch, &
  !$OMP d1_lE, d1_lB, d1_lPsi, d1_lPhi, &
  !$OMP ad1_lE, ad1_lB, ad1_lPsi, ad1_lPhi, d1_f, &
  !$OMP cf1, cf2, divE, divB, curlE, curlB, &
  !$OMP rhs_lE, rhs_lB, rhs_lPsi, rhs_lPhi, &
  !$OMP jac, hes, i, j, k, di, dj, dk, a, b, c, m)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    ch        = chi(i,j,k)
    hh(1,1)   = hxx(i,j,k)
    hh(1,2)   = hxy(i,j,k)
    hh(1,3)   = hxz(i,j,k)
    hh(2,2)   = hyy(i,j,k)
    hh(2,3)   = hyz(i,j,k)
    hh(3,3)   = hzz(i,j,k)
    hh(2,1)   = hh(1,2)
    hh(3,1)   = hh(1,3)
    hh(3,2)   = hh(2,3)

    trk       = tracek(i,j,k)
    alph      = alp(i,j,k)
    beta(1)   = betax(i,j,k)
    beta(2)   = betay(i,j,k)
    beta(3)   = betaz(i,j,k)

    lE(1)     = Ex(i,j,k)
    lE(2)     = Ey(i,j,k)
    lE(3)     = Ez(i,j,k)
    lB(1)     = Bx(i,j,k)
    lB(2)     = By(i,j,k)
    lB(3)     = Bz(i,j,k)
    lPsi      = Psi(i,j,k)
    lPhi      = Phi(i,j,k)

    if (use_jacobian) then
       jac(1,1) = lJ11(i,j,k)
       jac(1,2) = lJ12(i,j,k)
       jac(1,3) = lJ13(i,j,k)
       jac(2,1) = lJ21(i,j,k)
       jac(2,2) = lJ22(i,j,k)
       jac(2,3) = lJ23(i,j,k)
       jac(3,1) = lJ31(i,j,k)
       jac(3,2) = lJ32(i,j,k)
       jac(3,3) = lJ33(i,j,k)

       hes(1,1,1) = ldJ111(i,j,k)
       hes(1,1,2) = ldJ112(i,j,k)
       hes(1,1,3) = ldJ113(i,j,k)
       hes(1,2,1) = ldJ112(i,j,k)
       hes(1,2,2) = ldJ122(i,j,k)
       hes(1,2,3) = ldJ123(i,j,k)
       hes(1,3,1) = ldJ113(i,j,k)
       hes(1,3,2) = ldJ123(i,j,k)
       hes(1,3,3) = ldJ133(i,j,k)
       hes(2,1,1) = ldJ211(i,j,k)
       hes(2,1,2) = ldJ212(i,j,k)
       hes(2,1,3) = ldJ213(i,j,k)
       hes(2,2,1) = ldJ212(i,j,k)
       hes(2,2,2) = ldJ222(i,j,k)
       hes(2,2,3) = ldJ223(i,j,k)
       hes(2,3,1) = ldJ213(i,j,k)
       hes(2,3,2) = ldJ223(i,j,k)
       hes(2,3,3) = ldJ233(i,j,k)
       hes(3,1,1) = ldJ311(i,j,k)
       hes(3,1,2) = ldJ312(i,j,k)
       hes(3,1,3) = ldJ313(i,j,k)
       hes(3,2,1) = ldJ312(i,j,k)
       hes(3,2,2) = ldJ322(i,j,k)
       hes(3,2,3) = ldJ323(i,j,k)
       hes(3,3,1) = ldJ313(i,j,k)
       hes(3,3,2) = ldJ323(i,j,k)
       hes(3,3,3) = ldJ333(i,j,k)
    else
       jac      = 0.0d0
       jac(1,1) = 1.0d0
       jac(2,2) = 1.0d0
       jac(3,3) = 1.0d0
       hes      = 0.0d0
    end if

    dethh =       hh(1,1) * hh(2,2) * hh(3,3)                              &
            + 2.0d0 * hh(1,2) * hh(1,3) * hh(2,3)                          &
            -        hh(1,1) * hh(2,3) ** 2                                 &
            -        hh(2,2) * hh(1,3) ** 2                                 &
            -        hh(3,3) * hh(1,2) ** 2
    hu(1,1) = (hh(2,2) * hh(3,3) - hh(2,3) ** 2     ) / dethh
    hu(2,2) = (hh(1,1) * hh(3,3) - hh(1,3) ** 2     ) / dethh
    hu(3,3) = (hh(1,1) * hh(2,2) - hh(1,2) ** 2     ) / dethh
    hu(1,2) = (hh(1,3) * hh(2,3) - hh(1,2) * hh(3,3)) / dethh
    hu(1,3) = (hh(1,2) * hh(2,3) - hh(1,3) * hh(2,2)) / dethh
    hu(2,3) = (hh(1,3) * hh(1,2) - hh(2,3) * hh(1,1)) / dethh
    hu(2,1) = hu(1,2)
    hu(3,1) = hu(1,3)
    hu(3,2) = hu(2,3)

    if (derivs_order == 4) then
      d1_ch(1) = D1X4(chi)
      d1_ch(2) = D1Y4(chi)
      d1_ch(3) = D1Z4(chi)

      d1_alph(1) = D1X4(alp)
      d1_alph(2) = D1Y4(alp)
      d1_alph(3) = D1Z4(alp)

      d1_lPsi(1) = D1X4(Psi)
      d1_lPsi(2) = D1Y4(Psi)
      d1_lPsi(3) = D1Z4(Psi)

      d1_lPhi(1) = D1X4(Phi)
      d1_lPhi(2) = D1Y4(Phi)
      d1_lPhi(3) = D1Z4(Phi)

      d1_hh(1,1,1) = D1X4(hxx)
      d1_hh(1,1,2) = D1Y4(hxx)
      d1_hh(1,1,3) = D1Z4(hxx)
      d1_hh(1,2,1) = D1X4(hxy)
      d1_hh(1,2,2) = D1Y4(hxy)
      d1_hh(1,2,3) = D1Z4(hxy)
      d1_hh(1,3,1) = D1X4(hxz)
      d1_hh(1,3,2) = D1Y4(hxz)
      d1_hh(1,3,3) = D1Z4(hxz)
      d1_hh(2,2,1) = D1X4(hyy)
      d1_hh(2,2,2) = D1Y4(hyy)
      d1_hh(2,2,3) = D1Z4(hyy)
      d1_hh(2,3,1) = D1X4(hyz)
      d1_hh(2,3,2) = D1Y4(hyz)
      d1_hh(2,3,3) = D1Z4(hyz)
      d1_hh(3,3,1) = D1X4(hzz)
      d1_hh(3,3,2) = D1Y4(hzz)
      d1_hh(3,3,3) = D1Z4(hzz)
      d1_hh(2,1,:) = d1_hh(1,2,:)
      d1_hh(3,1,:) = d1_hh(1,3,:)
      d1_hh(3,2,:) = d1_hh(2,3,:)

      d1_beta(1,1) = D1X4(betax)
      d1_beta(1,2) = D1Y4(betax)
      d1_beta(1,3) = D1Z4(betax)
      d1_beta(2,1) = D1X4(betay)
      d1_beta(2,2) = D1Y4(betay)
      d1_beta(2,3) = D1Z4(betay)
      d1_beta(3,1) = D1X4(betaz)
      d1_beta(3,2) = D1Y4(betaz)
      d1_beta(3,3) = D1Z4(betaz)

      d1_lE(1,1) = D1X4(Ex)
      d1_lE(1,2) = D1Y4(Ex)
      d1_lE(1,3) = D1Z4(Ex)
      d1_lE(2,1) = D1X4(Ey)
      d1_lE(2,2) = D1Y4(Ey)
      d1_lE(2,3) = D1Z4(Ey)
      d1_lE(3,1) = D1X4(Ez)
      d1_lE(3,2) = D1Y4(Ez)
      d1_lE(3,3) = D1Z4(Ez)

      d1_lB(1,1) = D1X4(Bx)
      d1_lB(1,2) = D1Y4(Bx)
      d1_lB(1,3) = D1Z4(Bx)
      d1_lB(2,1) = D1X4(By)
      d1_lB(2,2) = D1Y4(By)
      d1_lB(2,3) = D1Z4(By)
      d1_lB(3,1) = D1X4(Bz)
      d1_lB(3,2) = D1Y4(Bz)
      d1_lB(3,3) = D1Z4(Bz)

      if( use_advection_stencils /= 0 ) then
        di = int( sign( one, beta(1) ) )
        dj = int( sign( one, beta(2) ) )
        dk = int( sign( one, beta(3) ) )

        d1_f(1) = di * ( -3.0d0*Ex(i-di,j,k) - 10.0d0*Ex(i,j,k) + 18.0d0*Ex(i+di,j,k) - 6.0d0*Ex(i+2*di,j,k) + Ex(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3.0d0*Ex(i,j-dj,k) - 10.0d0*Ex(i,j,k) + 18.0d0*Ex(i,j+dj,k) - 6.0d0*Ex(i,j+2*dj,k) + Ex(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3.0d0*Ex(i,j,k-dk) - 10.0d0*Ex(i,j,k) + 18.0d0*Ex(i,j,k+dk) - 6.0d0*Ex(i,j,k+2*dk) + Ex(i,j,k+3*dk)) / dz12
        ad1_lE(1) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3.0d0*Ey(i-di,j,k) - 10.0d0*Ey(i,j,k) + 18.0d0*Ey(i+di,j,k) - 6.0d0*Ey(i+2*di,j,k) + Ey(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3.0d0*Ey(i,j-dj,k) - 10.0d0*Ey(i,j,k) + 18.0d0*Ey(i,j+dj,k) - 6.0d0*Ey(i,j+2*dj,k) + Ey(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3.0d0*Ey(i,j,k-dk) - 10.0d0*Ey(i,j,k) + 18.0d0*Ey(i,j,k+dk) - 6.0d0*Ey(i,j,k+2*dk) + Ey(i,j,k+3*dk)) / dz12
        ad1_lE(2) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3.0d0*Ez(i-di,j,k) - 10.0d0*Ez(i,j,k) + 18.0d0*Ez(i+di,j,k) - 6.0d0*Ez(i+2*di,j,k) + Ez(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3.0d0*Ez(i,j-dj,k) - 10.0d0*Ez(i,j,k) + 18.0d0*Ez(i,j+dj,k) - 6.0d0*Ez(i,j+2*dj,k) + Ez(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3.0d0*Ez(i,j,k-dk) - 10.0d0*Ez(i,j,k) + 18.0d0*Ez(i,j,k+dk) - 6.0d0*Ez(i,j,k+2*dk) + Ez(i,j,k+3*dk)) / dz12
        ad1_lE(3) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3.0d0*Bx(i-di,j,k) - 10.0d0*Bx(i,j,k) + 18.0d0*Bx(i+di,j,k) - 6.0d0*Bx(i+2*di,j,k) + Bx(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3.0d0*Bx(i,j-dj,k) - 10.0d0*Bx(i,j,k) + 18.0d0*Bx(i,j+dj,k) - 6.0d0*Bx(i,j+2*dj,k) + Bx(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3.0d0*Bx(i,j,k-dk) - 10.0d0*Bx(i,j,k) + 18.0d0*Bx(i,j,k+dk) - 6.0d0*Bx(i,j,k+2*dk) + Bx(i,j,k+3*dk)) / dz12
        ad1_lB(1) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3.0d0*By(i-di,j,k) - 10.0d0*By(i,j,k) + 18.0d0*By(i+di,j,k) - 6.0d0*By(i+2*di,j,k) + By(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3.0d0*By(i,j-dj,k) - 10.0d0*By(i,j,k) + 18.0d0*By(i,j+dj,k) - 6.0d0*By(i,j+2*dj,k) + By(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3.0d0*By(i,j,k-dk) - 10.0d0*By(i,j,k) + 18.0d0*By(i,j,k+dk) - 6.0d0*By(i,j,k+2*dk) + By(i,j,k+3*dk)) / dz12
        ad1_lB(2) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3.0d0*Bz(i-di,j,k) - 10.0d0*Bz(i,j,k) + 18.0d0*Bz(i+di,j,k) - 6.0d0*Bz(i+2*di,j,k) + Bz(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3.0d0*Bz(i,j-dj,k) - 10.0d0*Bz(i,j,k) + 18.0d0*Bz(i,j+dj,k) - 6.0d0*Bz(i,j+2*dj,k) + Bz(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3.0d0*Bz(i,j,k-dk) - 10.0d0*Bz(i,j,k) + 18.0d0*Bz(i,j,k+dk) - 6.0d0*Bz(i,j,k+2*dk) + Bz(i,j,k+3*dk)) / dz12
        ad1_lB(3) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3.0d0*Psi(i-di,j,k) - 10.0d0*Psi(i,j,k) + 18.0d0*Psi(i+di,j,k) - 6.0d0*Psi(i+2*di,j,k) + Psi(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3.0d0*Psi(i,j-dj,k) - 10.0d0*Psi(i,j,k) + 18.0d0*Psi(i,j+dj,k) - 6.0d0*Psi(i,j+2*dj,k) + Psi(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3.0d0*Psi(i,j,k-dk) - 10.0d0*Psi(i,j,k) + 18.0d0*Psi(i,j,k+dk) - 6.0d0*Psi(i,j,k+2*dk) + Psi(i,j,k+3*dk)) / dz12
        ad1_lPsi = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( -3.0d0*Phi(i-di,j,k) - 10.0d0*Phi(i,j,k) + 18.0d0*Phi(i+di,j,k) - 6.0d0*Phi(i+2*di,j,k) + Phi(i+3*di,j,k)) / dx12
        d1_f(2) = dj * ( -3.0d0*Phi(i,j-dj,k) - 10.0d0*Phi(i,j,k) + 18.0d0*Phi(i,j+dj,k) - 6.0d0*Phi(i,j+2*dj,k) + Phi(i,j+3*dj,k)) / dy12
        d1_f(3) = dk * ( -3.0d0*Phi(i,j,k-dk) - 10.0d0*Phi(i,j,k) + 18.0d0*Phi(i,j,k+dk) - 6.0d0*Phi(i,j,k+2*dk) + Phi(i,j,k+3*dk)) / dz12
        ad1_lPhi = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

      else
        ad1_lE(1) = beta(1)*d1_lE(1,1) + beta(2)*d1_lE(1,2) + beta(3)*d1_lE(1,3)
        ad1_lE(2) = beta(1)*d1_lE(2,1) + beta(2)*d1_lE(2,2) + beta(3)*d1_lE(2,3)
        ad1_lE(3) = beta(1)*d1_lE(3,1) + beta(2)*d1_lE(3,2) + beta(3)*d1_lE(3,3)
        ad1_lB(1) = beta(1)*d1_lB(1,1) + beta(2)*d1_lB(1,2) + beta(3)*d1_lB(1,3)
        ad1_lB(2) = beta(1)*d1_lB(2,1) + beta(2)*d1_lB(2,2) + beta(3)*d1_lB(2,3)
        ad1_lB(3) = beta(1)*d1_lB(3,1) + beta(2)*d1_lB(3,2) + beta(3)*d1_lB(3,3)
        ad1_lPsi = beta(1)*d1_lPsi(1) + beta(2)*d1_lPsi(2) + beta(3)*d1_lPsi(3)
        ad1_lPhi = beta(1)*d1_lPhi(1) + beta(2)*d1_lPhi(2) + beta(3)*d1_lPhi(3)
      end if

    else if (derivs_order == 6) then
      d1_ch(1) = D1X6(chi)
      d1_ch(2) = D1Y6(chi)
      d1_ch(3) = D1Z6(chi)

      d1_alph(1) = D1X6(alp)
      d1_alph(2) = D1Y6(alp)
      d1_alph(3) = D1Z6(alp)

      d1_lPsi(1) = D1X6(Psi)
      d1_lPsi(2) = D1Y6(Psi)
      d1_lPsi(3) = D1Z6(Psi)

      d1_lPhi(1) = D1X6(Phi)
      d1_lPhi(2) = D1Y6(Phi)
      d1_lPhi(3) = D1Z6(Phi)

      d1_hh(1,1,1) = D1X6(hxx)
      d1_hh(1,1,2) = D1Y6(hxx)
      d1_hh(1,1,3) = D1Z6(hxx)
      d1_hh(1,2,1) = D1X6(hxy)
      d1_hh(1,2,2) = D1Y6(hxy)
      d1_hh(1,2,3) = D1Z6(hxy)
      d1_hh(1,3,1) = D1X6(hxz)
      d1_hh(1,3,2) = D1Y6(hxz)
      d1_hh(1,3,3) = D1Z6(hxz)
      d1_hh(2,2,1) = D1X6(hyy)
      d1_hh(2,2,2) = D1Y6(hyy)
      d1_hh(2,2,3) = D1Z6(hyy)
      d1_hh(2,3,1) = D1X6(hyz)
      d1_hh(2,3,2) = D1Y6(hyz)
      d1_hh(2,3,3) = D1Z6(hyz)
      d1_hh(3,3,1) = D1X6(hzz)
      d1_hh(3,3,2) = D1Y6(hzz)
      d1_hh(3,3,3) = D1Z6(hzz)
      d1_hh(2,1,:) = d1_hh(1,2,:)
      d1_hh(3,1,:) = d1_hh(1,3,:)
      d1_hh(3,2,:) = d1_hh(2,3,:)

      d1_beta(1,1) = D1X6(betax)
      d1_beta(1,2) = D1Y6(betax)
      d1_beta(1,3) = D1Z6(betax)
      d1_beta(2,1) = D1X6(betay)
      d1_beta(2,2) = D1Y6(betay)
      d1_beta(2,3) = D1Z6(betay)
      d1_beta(3,1) = D1X6(betaz)
      d1_beta(3,2) = D1Y6(betaz)
      d1_beta(3,3) = D1Z6(betaz)

      d1_lE(1,1) = D1X6(Ex)
      d1_lE(1,2) = D1Y6(Ex)
      d1_lE(1,3) = D1Z6(Ex)
      d1_lE(2,1) = D1X6(Ey)
      d1_lE(2,2) = D1Y6(Ey)
      d1_lE(2,3) = D1Z6(Ey)
      d1_lE(3,1) = D1X6(Ez)
      d1_lE(3,2) = D1Y6(Ez)
      d1_lE(3,3) = D1Z6(Ez)

      d1_lB(1,1) = D1X6(Bx)
      d1_lB(1,2) = D1Y6(Bx)
      d1_lB(1,3) = D1Z6(Bx)
      d1_lB(2,1) = D1X6(By)
      d1_lB(2,2) = D1Y6(By)
      d1_lB(2,3) = D1Z6(By)
      d1_lB(3,1) = D1X6(Bz)
      d1_lB(3,2) = D1Y6(Bz)
      d1_lB(3,3) = D1Z6(Bz)

      if( use_advection_stencils /= 0 ) then
        di = int( sign( one, beta(1) ) )
        dj = int( sign( one, beta(2) ) )
        dk = int( sign( one, beta(3) ) )

        d1_f(1) = di * ( 2.0d0*Ex(i-2*di,j,k) - 24.0d0*Ex(i-di,j,k) - 35.0d0*Ex(i,j,k) + 80.0d0*Ex(i+di,j,k) - 30.0d0*Ex(i+2*di,j,k) + 8.0d0*Ex(i+3*di,j,k) - Ex(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * ( 2.0d0*Ex(i,j-2*dj,k) - 24.0d0*Ex(i,j-dj,k) - 35.0d0*Ex(i,j,k) + 80.0d0*Ex(i,j+dj,k) - 30.0d0*Ex(i,j+2*dj,k) + 8.0d0*Ex(i,j+3*dj,k) - Ex(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * ( 2.0d0*Ex(i,j,k-2*dk) - 24.0d0*Ex(i,j,k-dk) - 35.0d0*Ex(i,j,k) + 80.0d0*Ex(i,j,k+dk) - 30.0d0*Ex(i,j,k+2*dk) + 8.0d0*Ex(i,j,k+3*dk) - Ex(i,j,k+4*dk) ) * odz60
        ad1_lE(1) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( 2.0d0*Ey(i-2*di,j,k) - 24.0d0*Ey(i-di,j,k) - 35.0d0*Ey(i,j,k) + 80.0d0*Ey(i+di,j,k) - 30.0d0*Ey(i+2*di,j,k) + 8.0d0*Ey(i+3*di,j,k) - Ey(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * ( 2.0d0*Ey(i,j-2*dj,k) - 24.0d0*Ey(i,j-dj,k) - 35.0d0*Ey(i,j,k) + 80.0d0*Ey(i,j+dj,k) - 30.0d0*Ey(i,j+2*dj,k) + 8.0d0*Ey(i,j+3*dj,k) - Ey(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * ( 2.0d0*Ey(i,j,k-2*dk) - 24.0d0*Ey(i,j,k-dk) - 35.0d0*Ey(i,j,k) + 80.0d0*Ey(i,j,k+dk) - 30.0d0*Ey(i,j,k+2*dk) + 8.0d0*Ey(i,j,k+3*dk) - Ey(i,j,k+4*dk) ) * odz60
        ad1_lE(2) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( 2.0d0*Ez(i-2*di,j,k) - 24.0d0*Ez(i-di,j,k) - 35.0d0*Ez(i,j,k) + 80.0d0*Ez(i+di,j,k) - 30.0d0*Ez(i+2*di,j,k) + 8.0d0*Ez(i+3*di,j,k) - Ez(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * ( 2.0d0*Ez(i,j-2*dj,k) - 24.0d0*Ez(i,j-dj,k) - 35.0d0*Ez(i,j,k) + 80.0d0*Ez(i,j+dj,k) - 30.0d0*Ez(i,j+2*dj,k) + 8.0d0*Ez(i,j+3*dj,k) - Ez(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * ( 2.0d0*Ez(i,j,k-2*dk) - 24.0d0*Ez(i,j,k-dk) - 35.0d0*Ez(i,j,k) + 80.0d0*Ez(i,j,k+dk) - 30.0d0*Ez(i,j,k+2*dk) + 8.0d0*Ez(i,j,k+3*dk) - Ez(i,j,k+4*dk) ) * odz60
        ad1_lE(3) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( 2.0d0*Bx(i-2*di,j,k) - 24.0d0*Bx(i-di,j,k) - 35.0d0*Bx(i,j,k) + 80.0d0*Bx(i+di,j,k) - 30.0d0*Bx(i+2*di,j,k) + 8.0d0*Bx(i+3*di,j,k) - Bx(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * ( 2.0d0*Bx(i,j-2*dj,k) - 24.0d0*Bx(i,j-dj,k) - 35.0d0*Bx(i,j,k) + 80.0d0*Bx(i,j+dj,k) - 30.0d0*Bx(i,j+2*dj,k) + 8.0d0*Bx(i,j+3*dj,k) - Bx(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * ( 2.0d0*Bx(i,j,k-2*dk) - 24.0d0*Bx(i,j,k-dk) - 35.0d0*Bx(i,j,k) + 80.0d0*Bx(i,j,k+dk) - 30.0d0*Bx(i,j,k+2*dk) + 8.0d0*Bx(i,j,k+3*dk) - Bx(i,j,k+4*dk) ) * odz60
        ad1_lB(1) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( 2.0d0*By(i-2*di,j,k) - 24.0d0*By(i-di,j,k) - 35.0d0*By(i,j,k) + 80.0d0*By(i+di,j,k) - 30.0d0*By(i+2*di,j,k) + 8.0d0*By(i+3*di,j,k) - By(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * ( 2.0d0*By(i,j-2*dj,k) - 24.0d0*By(i,j-dj,k) - 35.0d0*By(i,j,k) + 80.0d0*By(i,j+dj,k) - 30.0d0*By(i,j+2*dj,k) + 8.0d0*By(i,j+3*dj,k) - By(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * ( 2.0d0*By(i,j,k-2*dk) - 24.0d0*By(i,j,k-dk) - 35.0d0*By(i,j,k) + 80.0d0*By(i,j,k+dk) - 30.0d0*By(i,j,k+2*dk) + 8.0d0*By(i,j,k+3*dk) - By(i,j,k+4*dk) ) * odz60
        ad1_lB(2) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( 2.0d0*Bz(i-2*di,j,k) - 24.0d0*Bz(i-di,j,k) - 35.0d0*Bz(i,j,k) + 80.0d0*Bz(i+di,j,k) - 30.0d0*Bz(i+2*di,j,k) + 8.0d0*Bz(i+3*di,j,k) - Bz(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * ( 2.0d0*Bz(i,j-2*dj,k) - 24.0d0*Bz(i,j-dj,k) - 35.0d0*Bz(i,j,k) + 80.0d0*Bz(i,j+dj,k) - 30.0d0*Bz(i,j+2*dj,k) + 8.0d0*Bz(i,j+3*dj,k) - Bz(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * ( 2.0d0*Bz(i,j,k-2*dk) - 24.0d0*Bz(i,j,k-dk) - 35.0d0*Bz(i,j,k) + 80.0d0*Bz(i,j,k+dk) - 30.0d0*Bz(i,j,k+2*dk) + 8.0d0*Bz(i,j,k+3*dk) - Bz(i,j,k+4*dk) ) * odz60
        ad1_lB(3) = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( 2.0d0*Psi(i-2*di,j,k) - 24.0d0*Psi(i-di,j,k) - 35.0d0*Psi(i,j,k) + 80.0d0*Psi(i+di,j,k) - 30.0d0*Psi(i+2*di,j,k) + 8.0d0*Psi(i+3*di,j,k) - Psi(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * ( 2.0d0*Psi(i,j-2*dj,k) - 24.0d0*Psi(i,j-dj,k) - 35.0d0*Psi(i,j,k) + 80.0d0*Psi(i,j+dj,k) - 30.0d0*Psi(i,j+2*dj,k) + 8.0d0*Psi(i,j+3*dj,k) - Psi(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * ( 2.0d0*Psi(i,j,k-2*dk) - 24.0d0*Psi(i,j,k-dk) - 35.0d0*Psi(i,j,k) + 80.0d0*Psi(i,j,k+dk) - 30.0d0*Psi(i,j,k+2*dk) + 8.0d0*Psi(i,j,k+3*dk) - Psi(i,j,k+4*dk) ) * odz60
        ad1_lPsi = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

        d1_f(1) = di * ( 2.0d0*Phi(i-2*di,j,k) - 24.0d0*Phi(i-di,j,k) - 35.0d0*Phi(i,j,k) + 80.0d0*Phi(i+di,j,k) - 30.0d0*Phi(i+2*di,j,k) + 8.0d0*Phi(i+3*di,j,k) - Phi(i+4*di,j,k) ) * odx60
        d1_f(2) = dj * ( 2.0d0*Phi(i,j-2*dj,k) - 24.0d0*Phi(i,j-dj,k) - 35.0d0*Phi(i,j,k) + 80.0d0*Phi(i,j+dj,k) - 30.0d0*Phi(i,j+2*dj,k) + 8.0d0*Phi(i,j+3*dj,k) - Phi(i,j+4*dj,k) ) * ody60
        d1_f(3) = dk * ( 2.0d0*Phi(i,j,k-2*dk) - 24.0d0*Phi(i,j,k-dk) - 35.0d0*Phi(i,j,k) + 80.0d0*Phi(i,j,k+dk) - 30.0d0*Phi(i,j,k+2*dk) + 8.0d0*Phi(i,j,k+3*dk) - Phi(i,j,k+4*dk) ) * odz60
        ad1_lPhi = beta(1)*d1_f(1) + beta(2)*d1_f(2) + beta(3)*d1_f(3)

      else
        ad1_lE(1) = beta(1)*d1_lE(1,1) + beta(2)*d1_lE(1,2) + beta(3)*d1_lE(1,3)
        ad1_lE(2) = beta(1)*d1_lE(2,1) + beta(2)*d1_lE(2,2) + beta(3)*d1_lE(2,3)
        ad1_lE(3) = beta(1)*d1_lE(3,1) + beta(2)*d1_lE(3,2) + beta(3)*d1_lE(3,3)
        ad1_lB(1) = beta(1)*d1_lB(1,1) + beta(2)*d1_lB(1,2) + beta(3)*d1_lB(1,3)
        ad1_lB(2) = beta(1)*d1_lB(2,1) + beta(2)*d1_lB(2,2) + beta(3)*d1_lB(2,3)
        ad1_lB(3) = beta(1)*d1_lB(3,1) + beta(2)*d1_lB(3,2) + beta(3)*d1_lB(3,3)
        ad1_lPsi = beta(1)*d1_lPsi(1) + beta(2)*d1_lPsi(2) + beta(3)*d1_lPsi(3)
        ad1_lPhi = beta(1)*d1_lPhi(1) + beta(2)*d1_lPhi(2) + beta(3)*d1_lPhi(3)
      end if

    else
      call CCTK_WARN(0, "derivs_order not implemented in EinsteinMaxwell_calc_rhs_EM_direct_corrected.F90")
    end if

    if (use_jacobian) then
      call EinsteinMaxwell_d1_Scalar_apply_jacobian(d1_ch, jac)
      call EinsteinMaxwell_d1_Scalar_apply_jacobian(d1_alph, jac)
      call EinsteinMaxwell_d1_Scalar_apply_jacobian(d1_lPsi, jac)
      call EinsteinMaxwell_d1_Scalar_apply_jacobian(d1_lPhi, jac)
      call EinsteinMaxwell_d1_Vector_apply_jacobian(d1_beta, jac)
      call EinsteinMaxwell_d1_Vector_apply_jacobian(d1_lE, jac)
      call EinsteinMaxwell_d1_Vector_apply_jacobian(d1_lB, jac)
      call EinsteinMaxwell_d1_2nd_Tensor_apply_jacobian(d1_hh, jac, 1)
    end if

    ! Conformal Christoffel symbols of h_ij = tilde gamma_ij.
    cf1 = 0.0d0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          cf1(a,b,c) = 0.5d0 * (d1_hh(a,b,c) + d1_hh(a,c,b) - d1_hh(b,c,a))
        end do
      end do
    end do
    cf1(:,2,1) = cf1(:,1,2)
    cf1(:,3,1) = cf1(:,1,3)
    cf1(:,3,2) = cf1(:,2,3)

    cf2 = 0.0d0
    do a = 1, 3
      do b = 1, 3
        do c = b, 3
          do m = 1, 3
            cf2(a,b,c) = cf2(a,b,c) + hu(a,m) * cf1(m,b,c)
          end do
        end do
      end do
    end do
    cf2(:,2,1) = cf2(:,1,2)
    cf2(:,3,1) = cf2(:,1,3)
    cf2(:,3,2) = cf2(:,2,3)

    ! Physical divergence D_i V^i for gamma_ij = chi^{-p} h_ij.
    ! D_i V^i = partial_i V^i + tildeGamma^i_{m i} V^m - (3p/2) V^i partial_i chi / chi.
    divE = 0.0d0
    divB = 0.0d0
    do a = 1, 3
      divE = divE + d1_lE(a,a) - 1.5d0*conf_fac_exponent*lE(a)*d1_ch(a)/ch
      divB = divB + d1_lB(a,a) - 1.5d0*conf_fac_exponent*lB(a)*d1_ch(a)/ch
      do m = 1, 3
        divE = divE + cf2(a,m,a) * lE(m)
        divB = divB + cf2(a,m,a) * lB(m)
      end do
    end do

    rhs_lE = ad1_lE
    do a = 1, 3
      do m = 1, 3
        rhs_lE(a) = rhs_lE(a) - lE(m) * d1_beta(a,m)
      end do
    end do
    rhs_lE = rhs_lE + alph * trk * lE

    do a = 1, 3
      do b = 1, 3
        rhs_lE(a) = rhs_lE(a) - alph * ch**conf_fac_exponent * hu(a,b) * d1_lPsi(b)
      end do
      curlB = 0.0d0
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            curlB = curlB + levi_civita(a,b,c) * ch**(-conf_fac_exponent) *          &
                  ( hh(c,m) * lB(m) * d1_alph(b)                                    &
                  + alph * ( lB(m) * d1_hh(c,m,b) + hh(c,m) * d1_lB(m,b)             &
                  - ch**(-conf_fac_exponent) * hh(c,m) * lB(m) * d1_ch(b) ) )
          end do
        end do
      end do
      rhs_lE(a) = rhs_lE(a) + curlB
    end do

    rhs_lB = ad1_lB
    do a = 1, 3
      do m = 1, 3
        rhs_lB(a) = rhs_lB(a) - lB(m) * d1_beta(a,m)
      end do
    end do
    rhs_lB = rhs_lB + alph * trk * lB

    do a = 1, 3
      do b = 1, 3
        rhs_lB(a) = rhs_lB(a) - alph * ch**conf_fac_exponent * hu(a,b) * d1_lPhi(b)
      end do
      curlE = 0.0d0
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            curlE = curlE + levi_civita(a,b,c) * ch**(-conf_fac_exponent) *          &
                  ( hh(c,m) * lE(m) * d1_alph(b)                                    &
                  + alph * ( lE(m) * d1_hh(c,m,b) + hh(c,m) * d1_lE(m,b)             &
                  - ch**(-conf_fac_exponent) * hh(c,m) * lE(m) * d1_ch(b) ) )
          end do
        end do
      end do
      rhs_lB(a) = rhs_lB(a) - curlE
    end do

    rhs_lPsi = ad1_lPsi - alph * divE + alph * kappa * lPsi
    rhs_lPhi = ad1_lPhi - alph * divB + alph * kappa * lPhi

    rhs_Ex(i,j,k) = rhs_lE(1)
    rhs_Ey(i,j,k) = rhs_lE(2)
    rhs_Ez(i,j,k) = rhs_lE(3)
    rhs_Bx(i,j,k) = rhs_lB(1)
    rhs_By(i,j,k) = rhs_lB(2)
    rhs_Bz(i,j,k) = rhs_lB(3)
    rhs_Psi(i,j,k) = rhs_lPsi
    rhs_Phi(i,j,k) = rhs_lPhi

  end do
  end do
  end do
  !$OMP END PARALLEL DO

end subroutine EinsteinMaxwell_calc_rhs

subroutine EinsteinMaxwell_calc_rhs_bdry( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_REAL, parameter :: one  = 1.0d0
  CCTK_REAL, parameter :: zero = 0.0d0
  CCTK_INT ierr

  ierr = NewRad_Apply(cctkGH, Ex, rhs_Ex, E0(1), one, n_E(1))
  ierr = NewRad_Apply(cctkGH, Ey, rhs_Ey, E0(2), one, n_E(2))
  ierr = NewRad_Apply(cctkGH, Ez, rhs_Ez, E0(3), one, n_E(3))

  ierr = NewRad_Apply(cctkGH, Bx, rhs_Bx, B0(1), one, n_B(1))
  ierr = NewRad_Apply(cctkGH, By, rhs_By, B0(2), one, n_B(2))
  ierr = NewRad_Apply(cctkGH, Bz, rhs_Bz, B0(3), one, n_B(3))

  ierr = NewRad_Apply(cctkGH, Psi, rhs_Psi, zero, one, n_Psi)
  ierr = NewRad_Apply(cctkGH, Phi, rhs_Phi, zero, one, n_Phi)

end subroutine EinsteinMaxwell_calc_rhs_bdry

#undef D1X4
#undef D1Y4
#undef D1Z4
#undef D1X6
#undef D1Y6
#undef D1Z6

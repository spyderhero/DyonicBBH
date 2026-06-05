! NPScalars
! NP_calc_psi_gfs.F90 : Actual calculation of the NP scalar grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine NPProca_calcPsiGF( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS

  CCTK_REAL gd(3,3), gu(3,3), kd(3,3), detgd
  CCTK_REAL Ed(3), Bd(3), Bu(3)
  CCTK_REAL d1_gd(3,3,3), d1_gu(3,3,3), d2_gd(3,3,3,3), d1_kd(3,3,3),      &
            cd1_kd(3,3,3)
  CCTK_REAL cf1(3,3,3), cf2(3,3,3), ri(3,3), d1_cf2(3,3,3,3)
  CCTK_REAL dx(3), xx(3)
  CCTK_REAL dx12, dy12, dz12, dxsq12, dysq12, dzsq12, dxdy144, dxdz144,    &
            dydz144
  CCTK_REAL odx60, ody60, odz60, odxsq180, odysq180, odzsq180,             &
            odxdy3600, odxdz3600, odydz3600
  CCTK_REAL dx2, dy2, dz2, dxsq, dysq, dzsq, dxdy4, dxdz4, dydz4
  CCTK_REAL u_vec(3), v_vec(3), w_vec(3), dotp1, dotp2
  CCTK_REAL eps_axis
  CCTK_REAL eps_lc_u(3,3,3), elec(3,3), mag(3,3), electr, magtr

  CCTK_INT  i, j, k, m, n, p, a, b, c, d

  eps_axis = 1.0d-10

  if (calculate_NP_every .le. 0) then
     return
  end if

  if (MOD(cctk_iteration, calculate_NP_every) .ne. 0 ) then
     return
  endif

  dx(:)   = CCTK_DELTA_SPACE(:)

  !--- coefficients for sixth order
  odx60 = 1 / (60 * CCTK_DELTA_SPACE(1))
  ody60 = 1 / (60 * CCTK_DELTA_SPACE(2))
  odz60 = 1 / (60 * CCTK_DELTA_SPACE(3))

  odxsq180 = 1 / (180*CCTK_DELTA_SPACE(1)**2)
  odysq180 = 1 / (180*CCTK_DELTA_SPACE(2)**2)
  odzsq180 = 1 / (180*CCTK_DELTA_SPACE(3)**2)

  odxdy3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2))
  odxdz3600 = 1 / (3600*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(3))
  odydz3600 = 1 / (3600*CCTK_DELTA_SPACE(2)*CCTK_DELTA_SPACE(3))

  !--- coefficients for fourth order
  dx12    = 12 * dx(1)
  dy12    = 12 * dx(2)
  dz12    = 12 * dx(3)
  dxsq12  = 12 * dx(1)**2
  dysq12  = 12 * dx(2)**2
  dzsq12  = 12 * dx(3)**2
  dxdy144 = dx12 * dy12
  dxdz144 = dx12 * dz12
  dydz144 = dy12 * dz12

  !--- coefficients for second order
  dx2   = 2 * dx(1)
  dy2   = 2 * dx(2)
  dz2   = 2 * dx(3)
  dxsq  = dx(1)**2
  dysq  = dx(2)**2
  dzsq  = dx(3)**2
  dxdy4 = 4*dx(1) * dx(2)
  dxdz4 = 4*dx(1) * dx(3)
  dydz4 = 4*dx(2) * dx(3)

  !=== initialize source terms (gfs) as zero ===

  if (calc_psi /= 0) then
     psi0re = 0.0
     psi0im = 0.0
     psi4re = 0.0
     psi4im = 0.0
  end if

  if (calc_phi /= 0) then
     phi0re = 0.0
     phi0im = 0.0
     phi1re = 0.0
     phi1im = 0.0
     phi2re = 0.0
     phi2im = 0.0
  end if

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(gd, kd, &
  !$OMP Ed, detgd, gu, &
  !$OMP d1_gd, d1_kd, d2_gd, &
  !$OMP cf1, cf2, cd1_kd, d1_gu, d1_cf2, &
  !$OMP ri, eps_lc_u, Bu, Bd, &
  !$OMP xx, u_vec, v_vec, w_vec, dotp1, dotp2, &
  !$OMP elec, mag, electr, magtr, &
  !$OMP i, j, k, &
  !$OMP a, b, c, m, n)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
  do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
  do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

    !--------------Get local variables ----------
    gd(1,1) = gxx(i,j,k)
    gd(1,2) = gxy(i,j,k)
    gd(1,3) = gxz(i,j,k)
    gd(2,2) = gyy(i,j,k)
    gd(2,3) = gyz(i,j,k)
    gd(3,3) = gzz(i,j,k)
    gd(2,1) = gd(1,2)
    gd(3,1) = gd(1,3)
    gd(3,2) = gd(2,3)

    kd(1,1) = kxx(i,j,k)
    kd(1,2) = kxy(i,j,k)
    kd(1,3) = kxz(i,j,k)
    kd(2,2) = kyy(i,j,k)
    kd(2,3) = kyz(i,j,k)
    kd(3,3) = kzz(i,j,k)
    kd(2,1) = kd(1,2)
    kd(3,1) = kd(1,3)
    kd(3,2) = kd(2,3)

    ! Einstein-Maxwell direct field variables: no vector potential is used.

    Ed(1)   = gd(1,1) * Ex(i,j,k) + gd(1,2) * Ey(i,j,k) + gd(1,3) * Ez(i,j,k)
    Ed(2)   = gd(2,1) * Ex(i,j,k) + gd(2,2) * Ey(i,j,k) + gd(2,3) * Ez(i,j,k)
    Ed(3)   = gd(3,1) * Ex(i,j,k) + gd(3,2) * Ey(i,j,k) + gd(3,3) * Ez(i,j,k)

    !--------------------------------------------


    !-------------- Invert metric ---------------
    detgd =       gd(1,1) * gd(2,2) * gd(3,3)                                &
            + 2.0d0 * gd(1,2) * gd(1,3) * gd(2,3)                            &
            -        gd(1,1) * gd(2,3) ** 2                                  &
            -        gd(2,2) * gd(1,3) ** 2                                  &
            -        gd(3,3) * gd(1,2) ** 2

    if (detgd <= 0.0d0) then
       if (calc_psi /= 0) then
          psi0re(i,j,k) = 0.0d0
          psi0im(i,j,k) = 0.0d0
          psi4re(i,j,k) = 0.0d0
          psi4im(i,j,k) = 0.0d0
       end if
       if (calc_phi /= 0) then
          phi0re(i,j,k) = 0.0d0
          phi0im(i,j,k) = 0.0d0
          phi1re(i,j,k) = 0.0d0
          phi1im(i,j,k) = 0.0d0
          phi2re(i,j,k) = 0.0d0
          phi2im(i,j,k) = 0.0d0
       end if
       cycle
    end if

    gu(1,1) = (gd(2,2) * gd(3,3) - gd(2,3) ** 2     ) / detgd
    gu(2,2) = (gd(1,1) * gd(3,3) - gd(1,3) ** 2     ) / detgd
    gu(3,3) = (gd(1,1) * gd(2,2) - gd(1,2) ** 2     ) / detgd
    gu(1,2) = (gd(1,3) * gd(2,3) - gd(1,2) * gd(3,3)) / detgd
    gu(1,3) = (gd(1,2) * gd(2,3) - gd(1,3) * gd(2,2)) / detgd
    gu(2,3) = (gd(1,3) * gd(1,2) - gd(2,3) * gd(1,1)) / detgd
    gu(2,1) = gu(1,2)
    gu(3,1) = gu(1,3)
    gu(3,2) = gu(2,3)
    !--------------------------------------------


    !--------------------------------------------
    if ( NP_order == 6 ) then
      ! Sixth order derivatives

      !-------------- Centered 1st derivatives ----

      ! d1_gd(3,3,3)
      d1_gd(1,1,1) = (  gxx(i+3,j,k) - 9*gxx(i+2,j,k) + 45*gxx(i+1,j,k)      &
                      - gxx(i-3,j,k) + 9*gxx(i-2,j,k) - 45*gxx(i-1,j,k) ) * odx60
      d1_gd(1,2,1) = (  gxy(i+3,j,k) - 9*gxy(i+2,j,k) + 45*gxy(i+1,j,k)      &
                      - gxy(i-3,j,k) + 9*gxy(i-2,j,k) - 45*gxy(i-1,j,k) ) * odx60
      d1_gd(1,3,1) = (  gxz(i+3,j,k) - 9*gxz(i+2,j,k) + 45*gxz(i+1,j,k)      &
                      - gxz(i-3,j,k) + 9*gxz(i-2,j,k) - 45*gxz(i-1,j,k) ) * odx60
      d1_gd(2,2,1) = (  gyy(i+3,j,k) - 9*gyy(i+2,j,k) + 45*gyy(i+1,j,k)      &
                      - gyy(i-3,j,k) + 9*gyy(i-2,j,k) - 45*gyy(i-1,j,k) ) * odx60
      d1_gd(2,3,1) = (  gyz(i+3,j,k) - 9*gyz(i+2,j,k) + 45*gyz(i+1,j,k)      &
                      - gyz(i-3,j,k) + 9*gyz(i-2,j,k) - 45*gyz(i-1,j,k) ) * odx60
      d1_gd(3,3,1) = (  gzz(i+3,j,k) - 9*gzz(i+2,j,k) + 45*gzz(i+1,j,k)      &
                      - gzz(i-3,j,k) + 9*gzz(i-2,j,k) - 45*gzz(i-1,j,k) ) * odx60

      d1_gd(1,1,2) = (  gxx(i,j+3,k) - 9*gxx(i,j+2,k) + 45*gxx(i,j+1,k)      &
                      - gxx(i,j-3,k) + 9*gxx(i,j-2,k) - 45*gxx(i,j-1,k) ) * ody60
      d1_gd(1,2,2) = (  gxy(i,j+3,k) - 9*gxy(i,j+2,k) + 45*gxy(i,j+1,k)      &
                      - gxy(i,j-3,k) + 9*gxy(i,j-2,k) - 45*gxy(i,j-1,k) ) * ody60
      d1_gd(1,3,2) = (  gxz(i,j+3,k) - 9*gxz(i,j+2,k) + 45*gxz(i,j+1,k)      &
                      - gxz(i,j-3,k) + 9*gxz(i,j-2,k) - 45*gxz(i,j-1,k) ) * ody60
      d1_gd(2,2,2) = (  gyy(i,j+3,k) - 9*gyy(i,j+2,k) + 45*gyy(i,j+1,k)      &
                      - gyy(i,j-3,k) + 9*gyy(i,j-2,k) - 45*gyy(i,j-1,k) ) * ody60
      d1_gd(2,3,2) = (  gyz(i,j+3,k) - 9*gyz(i,j+2,k) + 45*gyz(i,j+1,k)      &
                      - gyz(i,j-3,k) + 9*gyz(i,j-2,k) - 45*gyz(i,j-1,k) ) * ody60
      d1_gd(3,3,2) = (  gzz(i,j+3,k) - 9*gzz(i,j+2,k) + 45*gzz(i,j+1,k)      &
                      - gzz(i,j-3,k) + 9*gzz(i,j-2,k) - 45*gzz(i,j-1,k) ) * ody60

      d1_gd(1,1,3) = (  gxx(i,j,k+3) - 9*gxx(i,j,k+2) + 45*gxx(i,j,k+1)      &
                      - gxx(i,j,k-3) + 9*gxx(i,j,k-2) - 45*gxx(i,j,k-1) ) * odz60
      d1_gd(1,2,3) = (  gxy(i,j,k+3) - 9*gxy(i,j,k+2) + 45*gxy(i,j,k+1)      &
                      - gxy(i,j,k-3) + 9*gxy(i,j,k-2) - 45*gxy(i,j,k-1) ) * odz60
      d1_gd(1,3,3) = (  gxz(i,j,k+3) - 9*gxz(i,j,k+2) + 45*gxz(i,j,k+1)      &
                      - gxz(i,j,k-3) + 9*gxz(i,j,k-2) - 45*gxz(i,j,k-1) ) * odz60
      d1_gd(2,2,3) = (  gyy(i,j,k+3) - 9*gyy(i,j,k+2) + 45*gyy(i,j,k+1)      &
                      - gyy(i,j,k-3) + 9*gyy(i,j,k-2) - 45*gyy(i,j,k-1) ) * odz60
      d1_gd(2,3,3) = (  gyz(i,j,k+3) - 9*gyz(i,j,k+2) + 45*gyz(i,j,k+1)      &
                      - gyz(i,j,k-3) + 9*gyz(i,j,k-2) - 45*gyz(i,j,k-1) ) * odz60
      d1_gd(3,3,3) = (  gzz(i,j,k+3) - 9*gzz(i,j,k+2) + 45*gzz(i,j,k+1)      &
                      - gzz(i,j,k-3) + 9*gzz(i,j,k-2) - 45*gzz(i,j,k-1) ) * odz60

      d1_gd(2,1,:) = d1_gd(1,2,:)
      d1_gd(3,1,:) = d1_gd(1,3,:)
      d1_gd(3,2,:) = d1_gd(2,3,:)


      ! d1_kd(3,3,3)
      d1_kd(1,1,1) = (  kxx(i+3,j,k) - 9*kxx(i+2,j,k) + 45*kxx(i+1,j,k)      &
                      - kxx(i-3,j,k) + 9*kxx(i-2,j,k) - 45*kxx(i-1,j,k) ) * odx60
      d1_kd(1,2,1) = (  kxy(i+3,j,k) - 9*kxy(i+2,j,k) + 45*kxy(i+1,j,k)      &
                      - kxy(i-3,j,k) + 9*kxy(i-2,j,k) - 45*kxy(i-1,j,k) ) * odx60
      d1_kd(1,3,1) = (  kxz(i+3,j,k) - 9*kxz(i+2,j,k) + 45*kxz(i+1,j,k)      &
                      - kxz(i-3,j,k) + 9*kxz(i-2,j,k) - 45*kxz(i-1,j,k) ) * odx60
      d1_kd(2,2,1) = (  kyy(i+3,j,k) - 9*kyy(i+2,j,k) + 45*kyy(i+1,j,k)      &
                      - kyy(i-3,j,k) + 9*kyy(i-2,j,k) - 45*kyy(i-1,j,k) ) * odx60
      d1_kd(2,3,1) = (  kyz(i+3,j,k) - 9*kyz(i+2,j,k) + 45*kyz(i+1,j,k)      &
                      - kyz(i-3,j,k) + 9*kyz(i-2,j,k) - 45*kyz(i-1,j,k) ) * odx60
      d1_kd(3,3,1) = (  kzz(i+3,j,k) - 9*kzz(i+2,j,k) + 45*kzz(i+1,j,k)      &
                      - kzz(i-3,j,k) + 9*kzz(i-2,j,k) - 45*kzz(i-1,j,k) ) * odx60

      d1_kd(1,1,2) = (  kxx(i,j+3,k) - 9*kxx(i,j+2,k) + 45*kxx(i,j+1,k)      &
                      - kxx(i,j-3,k) + 9*kxx(i,j-2,k) - 45*kxx(i,j-1,k) ) * ody60
      d1_kd(1,2,2) = (  kxy(i,j+3,k) - 9*kxy(i,j+2,k) + 45*kxy(i,j+1,k)      &
                      - kxy(i,j-3,k) + 9*kxy(i,j-2,k) - 45*kxy(i,j-1,k) ) * ody60
      d1_kd(1,3,2) = (  kxz(i,j+3,k) - 9*kxz(i,j+2,k) + 45*kxz(i,j+1,k)      &
                      - kxz(i,j-3,k) + 9*kxz(i,j-2,k) - 45*kxz(i,j-1,k) ) * ody60
      d1_kd(2,2,2) = (  kyy(i,j+3,k) - 9*kyy(i,j+2,k) + 45*kyy(i,j+1,k)      &
                      - kyy(i,j-3,k) + 9*kyy(i,j-2,k) - 45*kyy(i,j-1,k) ) * ody60
      d1_kd(2,3,2) = (  kyz(i,j+3,k) - 9*kyz(i,j+2,k) + 45*kyz(i,j+1,k)      &
                      - kyz(i,j-3,k) + 9*kyz(i,j-2,k) - 45*kyz(i,j-1,k) ) * ody60
      d1_kd(3,3,2) = (  kzz(i,j+3,k) - 9*kzz(i,j+2,k) + 45*kzz(i,j+1,k)      &
                      - kzz(i,j-3,k) + 9*kzz(i,j-2,k) - 45*kzz(i,j-1,k) ) * ody60

      d1_kd(1,1,3) = (  kxx(i,j,k+3) - 9*kxx(i,j,k+2) + 45*kxx(i,j,k+1)      &
                      - kxx(i,j,k-3) + 9*kxx(i,j,k-2) - 45*kxx(i,j,k-1) ) * odz60
      d1_kd(1,2,3) = (  kxy(i,j,k+3) - 9*kxy(i,j,k+2) + 45*kxy(i,j,k+1)      &
                      - kxy(i,j,k-3) + 9*kxy(i,j,k-2) - 45*kxy(i,j,k-1) ) * odz60
      d1_kd(1,3,3) = (  kxz(i,j,k+3) - 9*kxz(i,j,k+2) + 45*kxz(i,j,k+1)      &
                      - kxz(i,j,k-3) + 9*kxz(i,j,k-2) - 45*kxz(i,j,k-1) ) * odz60
      d1_kd(2,2,3) = (  kyy(i,j,k+3) - 9*kyy(i,j,k+2) + 45*kyy(i,j,k+1)      &
                      - kyy(i,j,k-3) + 9*kyy(i,j,k-2) - 45*kyy(i,j,k-1) ) * odz60
      d1_kd(2,3,3) = (  kyz(i,j,k+3) - 9*kyz(i,j,k+2) + 45*kyz(i,j,k+1)      &
                      - kyz(i,j,k-3) + 9*kyz(i,j,k-2) - 45*kyz(i,j,k-1) ) * odz60
      d1_kd(3,3,3) = (  kzz(i,j,k+3) - 9*kzz(i,j,k+2) + 45*kzz(i,j,k+1)      &
                      - kzz(i,j,k-3) + 9*kzz(i,j,k-2) - 45*kzz(i,j,k-1) ) * odz60

      d1_kd(2,1,:) = d1_kd(1,2,:)
      d1_kd(3,1,:) = d1_kd(1,3,:)
      d1_kd(3,2,:) = d1_kd(2,3,:)

      ! d1_lA removed: B^i is evolved directly in Einstein-Maxwell.


      !------------ Centered 2nd derivatives -----

      ! d2_gd(3,3,3,3)
      d2_gd(1,1,1,1) = (  2*gxx(i+3,j,k) - 27*gxx(i+2,j,k) + 270*gxx(i+1,j,k) - 490*gxx(i,j,k)&
                        + 2*gxx(i-3,j,k) - 27*gxx(i-2,j,k) + 270*gxx(i-1,j,k) ) * odxsq180
      d2_gd(1,2,1,1) = (  2*gxy(i+3,j,k) - 27*gxy(i+2,j,k) + 270*gxy(i+1,j,k) - 490*gxy(i,j,k)&
                        + 2*gxy(i-3,j,k) - 27*gxy(i-2,j,k) + 270*gxy(i-1,j,k) ) * odxsq180
      d2_gd(1,3,1,1) = (  2*gxz(i+3,j,k) - 27*gxz(i+2,j,k) + 270*gxz(i+1,j,k) - 490*gxz(i,j,k)&
                        + 2*gxz(i-3,j,k) - 27*gxz(i-2,j,k) + 270*gxz(i-1,j,k) ) * odxsq180
      d2_gd(2,2,1,1) = (  2*gyy(i+3,j,k) - 27*gyy(i+2,j,k) + 270*gyy(i+1,j,k) - 490*gyy(i,j,k)&
                        + 2*gyy(i-3,j,k) - 27*gyy(i-2,j,k) + 270*gyy(i-1,j,k) ) * odxsq180
      d2_gd(2,3,1,1) = (  2*gyz(i+3,j,k) - 27*gyz(i+2,j,k) + 270*gyz(i+1,j,k) - 490*gyz(i,j,k)&
                        + 2*gyz(i-3,j,k) - 27*gyz(i-2,j,k) + 270*gyz(i-1,j,k) ) * odxsq180
      d2_gd(3,3,1,1) = (  2*gzz(i+3,j,k) - 27*gzz(i+2,j,k) + 270*gzz(i+1,j,k) - 490*gzz(i,j,k)&
                        + 2*gzz(i-3,j,k) - 27*gzz(i-2,j,k) + 270*gzz(i-1,j,k) ) * odxsq180

      d2_gd(1,1,2,2) = (  2*gxx(i,j+3,k) - 27*gxx(i,j+2,k) + 270*gxx(i,j+1,k) - 490*gxx(i,j,k)&
                        + 2*gxx(i,j-3,k) - 27*gxx(i,j-2,k) + 270*gxx(i,j-1,k) ) * odysq180
      d2_gd(1,2,2,2) = (  2*gxy(i,j+3,k) - 27*gxy(i,j+2,k) + 270*gxy(i,j+1,k) - 490*gxy(i,j,k)&
                        + 2*gxy(i,j-3,k) - 27*gxy(i,j-2,k) + 270*gxy(i,j-1,k) ) * odysq180
      d2_gd(1,3,2,2) = (  2*gxz(i,j+3,k) - 27*gxz(i,j+2,k) + 270*gxz(i,j+1,k) - 490*gxz(i,j,k)&
                        + 2*gxz(i,j-3,k) - 27*gxz(i,j-2,k) + 270*gxz(i,j-1,k) ) * odysq180
      d2_gd(2,2,2,2) = (  2*gyy(i,j+3,k) - 27*gyy(i,j+2,k) + 270*gyy(i,j+1,k) - 490*gyy(i,j,k)&
                        + 2*gyy(i,j-3,k) - 27*gyy(i,j-2,k) + 270*gyy(i,j-1,k) ) * odysq180
      d2_gd(2,3,2,2) = (  2*gyz(i,j+3,k) - 27*gyz(i,j+2,k) + 270*gyz(i,j+1,k) - 490*gyz(i,j,k)&
                        + 2*gyz(i,j-3,k) - 27*gyz(i,j-2,k) + 270*gyz(i,j-1,k) ) * odysq180
      d2_gd(3,3,2,2) = (  2*gzz(i,j+3,k) - 27*gzz(i,j+2,k) + 270*gzz(i,j+1,k) - 490*gzz(i,j,k)&
                        + 2*gzz(i,j-3,k) - 27*gzz(i,j-2,k) + 270*gzz(i,j-1,k) ) * odysq180

      d2_gd(1,1,3,3) = (  2*gxx(i,j,k+3) - 27*gxx(i,j,k+2) + 270*gxx(i,j,k+1) - 490*gxx(i,j,k)&
                        + 2*gxx(i,j,k-3) - 27*gxx(i,j,k-2) + 270*gxx(i,j,k-1) ) * odzsq180
      d2_gd(1,2,3,3) = (  2*gxy(i,j,k+3) - 27*gxy(i,j,k+2) + 270*gxy(i,j,k+1) - 490*gxy(i,j,k)&
                        + 2*gxy(i,j,k-3) - 27*gxy(i,j,k-2) + 270*gxy(i,j,k-1) ) * odzsq180
      d2_gd(1,3,3,3) = (  2*gxz(i,j,k+3) - 27*gxz(i,j,k+2) + 270*gxz(i,j,k+1) - 490*gxz(i,j,k)&
                        + 2*gxz(i,j,k-3) - 27*gxz(i,j,k-2) + 270*gxz(i,j,k-1) ) * odzsq180
      d2_gd(2,2,3,3) = (  2*gyy(i,j,k+3) - 27*gyy(i,j,k+2) + 270*gyy(i,j,k+1) - 490*gyy(i,j,k)&
                        + 2*gyy(i,j,k-3) - 27*gyy(i,j,k-2) + 270*gyy(i,j,k-1) ) * odzsq180
      d2_gd(2,3,3,3) = (  2*gyz(i,j,k+3) - 27*gyz(i,j,k+2) + 270*gyz(i,j,k+1) - 490*gyz(i,j,k)&
                        + 2*gyz(i,j,k-3) - 27*gyz(i,j,k-2) + 270*gyz(i,j,k-1) ) * odzsq180
      d2_gd(3,3,3,3) = (  2*gzz(i,j,k+3) - 27*gzz(i,j,k+2) + 270*gzz(i,j,k+1) - 490*gzz(i,j,k)&
                        + 2*gzz(i,j,k-3) - 27*gzz(i,j,k-2) + 270*gzz(i,j,k-1) ) * odzsq180

      d2_gd(1,1,1,2) = (    -gxx(i-3,j+3,k) +   9*gxx(i-2,j+3,k) -   45*gxx(i-1,j+3,k) +   45*gxx(i+1,j+3,k) -   9*gxx(i+2,j+3,k) +    gxx(i+3,j+3,k) &
                        +  9*gxx(i-3,j+2,k) -  81*gxx(i-2,j+2,k) +  405*gxx(i-1,j+2,k) -  405*gxx(i+1,j+2,k) +  81*gxx(i+2,j+2,k) -  9*gxx(i+3,j+2,k) &
                        - 45*gxx(i-3,j+1,k) + 405*gxx(i-2,j+1,k) - 2025*gxx(i-1,j+1,k) + 2025*gxx(i+1,j+1,k) - 405*gxx(i+2,j+1,k) + 45*gxx(i+3,j+1,k) &
                        + 45*gxx(i-3,j-1,k) - 405*gxx(i-2,j-1,k) + 2025*gxx(i-1,j-1,k) - 2025*gxx(i+1,j-1,k) + 405*gxx(i+2,j-1,k) - 45*gxx(i+3,j-1,k) &
                        -  9*gxx(i-3,j-2,k) +  81*gxx(i-2,j-2,k) -  405*gxx(i-1,j-2,k) +  405*gxx(i+1,j-2,k) -  81*gxx(i+2,j-2,k) +  9*gxx(i+3,j-2,k) &
                        +    gxx(i-3,j-3,k) -   9*gxx(i-2,j-3,k) +   45*gxx(i-1,j-3,k) -   45*gxx(i+1,j-3,k) +   9*gxx(i+2,j-3,k) -    gxx(i+3,j-3,k) ) * odxdy3600
      d2_gd(1,2,1,2) = (    -gxy(i-3,j+3,k) +   9*gxy(i-2,j+3,k) -   45*gxy(i-1,j+3,k) +   45*gxy(i+1,j+3,k) -   9*gxy(i+2,j+3,k) +    gxy(i+3,j+3,k) &
                        +  9*gxy(i-3,j+2,k) -  81*gxy(i-2,j+2,k) +  405*gxy(i-1,j+2,k) -  405*gxy(i+1,j+2,k) +  81*gxy(i+2,j+2,k) -  9*gxy(i+3,j+2,k) &
                        - 45*gxy(i-3,j+1,k) + 405*gxy(i-2,j+1,k) - 2025*gxy(i-1,j+1,k) + 2025*gxy(i+1,j+1,k) - 405*gxy(i+2,j+1,k) + 45*gxy(i+3,j+1,k) &
                        + 45*gxy(i-3,j-1,k) - 405*gxy(i-2,j-1,k) + 2025*gxy(i-1,j-1,k) - 2025*gxy(i+1,j-1,k) + 405*gxy(i+2,j-1,k) - 45*gxy(i+3,j-1,k) &
                        -  9*gxy(i-3,j-2,k) +  81*gxy(i-2,j-2,k) -  405*gxy(i-1,j-2,k) +  405*gxy(i+1,j-2,k) -  81*gxy(i+2,j-2,k) +  9*gxy(i+3,j-2,k) &
                        +    gxy(i-3,j-3,k) -   9*gxy(i-2,j-3,k) +   45*gxy(i-1,j-3,k) -   45*gxy(i+1,j-3,k) +   9*gxy(i+2,j-3,k) -    gxy(i+3,j-3,k) ) * odxdy3600
      d2_gd(1,3,1,2) = (    -gxz(i-3,j+3,k) +   9*gxz(i-2,j+3,k) -   45*gxz(i-1,j+3,k) +   45*gxz(i+1,j+3,k) -   9*gxz(i+2,j+3,k) +    gxz(i+3,j+3,k) &
                        +  9*gxz(i-3,j+2,k) -  81*gxz(i-2,j+2,k) +  405*gxz(i-1,j+2,k) -  405*gxz(i+1,j+2,k) +  81*gxz(i+2,j+2,k) -  9*gxz(i+3,j+2,k) &
                        - 45*gxz(i-3,j+1,k) + 405*gxz(i-2,j+1,k) - 2025*gxz(i-1,j+1,k) + 2025*gxz(i+1,j+1,k) - 405*gxz(i+2,j+1,k) + 45*gxz(i+3,j+1,k) &
                        + 45*gxz(i-3,j-1,k) - 405*gxz(i-2,j-1,k) + 2025*gxz(i-1,j-1,k) - 2025*gxz(i+1,j-1,k) + 405*gxz(i+2,j-1,k) - 45*gxz(i+3,j-1,k) &
                        -  9*gxz(i-3,j-2,k) +  81*gxz(i-2,j-2,k) -  405*gxz(i-1,j-2,k) +  405*gxz(i+1,j-2,k) -  81*gxz(i+2,j-2,k) +  9*gxz(i+3,j-2,k) &
                        +    gxz(i-3,j-3,k) -   9*gxz(i-2,j-3,k) +   45*gxz(i-1,j-3,k) -   45*gxz(i+1,j-3,k) +   9*gxz(i+2,j-3,k) -    gxz(i+3,j-3,k) ) * odxdy3600
      d2_gd(2,2,1,2) = (    -gyy(i-3,j+3,k) +   9*gyy(i-2,j+3,k) -   45*gyy(i-1,j+3,k) +   45*gyy(i+1,j+3,k) -   9*gyy(i+2,j+3,k) +    gyy(i+3,j+3,k) &
                        +  9*gyy(i-3,j+2,k) -  81*gyy(i-2,j+2,k) +  405*gyy(i-1,j+2,k) -  405*gyy(i+1,j+2,k) +  81*gyy(i+2,j+2,k) -  9*gyy(i+3,j+2,k) &
                        - 45*gyy(i-3,j+1,k) + 405*gyy(i-2,j+1,k) - 2025*gyy(i-1,j+1,k) + 2025*gyy(i+1,j+1,k) - 405*gyy(i+2,j+1,k) + 45*gyy(i+3,j+1,k) &
                        + 45*gyy(i-3,j-1,k) - 405*gyy(i-2,j-1,k) + 2025*gyy(i-1,j-1,k) - 2025*gyy(i+1,j-1,k) + 405*gyy(i+2,j-1,k) - 45*gyy(i+3,j-1,k) &
                        -  9*gyy(i-3,j-2,k) +  81*gyy(i-2,j-2,k) -  405*gyy(i-1,j-2,k) +  405*gyy(i+1,j-2,k) -  81*gyy(i+2,j-2,k) +  9*gyy(i+3,j-2,k) &
                        +    gyy(i-3,j-3,k) -   9*gyy(i-2,j-3,k) +   45*gyy(i-1,j-3,k) -   45*gyy(i+1,j-3,k) +   9*gyy(i+2,j-3,k) -    gyy(i+3,j-3,k) ) * odxdy3600
      d2_gd(2,3,1,2) = (    -gyz(i-3,j+3,k) +   9*gyz(i-2,j+3,k) -   45*gyz(i-1,j+3,k) +   45*gyz(i+1,j+3,k) -   9*gyz(i+2,j+3,k) +    gyz(i+3,j+3,k) &
                        +  9*gyz(i-3,j+2,k) -  81*gyz(i-2,j+2,k) +  405*gyz(i-1,j+2,k) -  405*gyz(i+1,j+2,k) +  81*gyz(i+2,j+2,k) -  9*gyz(i+3,j+2,k) &
                        - 45*gyz(i-3,j+1,k) + 405*gyz(i-2,j+1,k) - 2025*gyz(i-1,j+1,k) + 2025*gyz(i+1,j+1,k) - 405*gyz(i+2,j+1,k) + 45*gyz(i+3,j+1,k) &
                        + 45*gyz(i-3,j-1,k) - 405*gyz(i-2,j-1,k) + 2025*gyz(i-1,j-1,k) - 2025*gyz(i+1,j-1,k) + 405*gyz(i+2,j-1,k) - 45*gyz(i+3,j-1,k) &
                        -  9*gyz(i-3,j-2,k) +  81*gyz(i-2,j-2,k) -  405*gyz(i-1,j-2,k) +  405*gyz(i+1,j-2,k) -  81*gyz(i+2,j-2,k) +  9*gyz(i+3,j-2,k) &
                        +    gyz(i-3,j-3,k) -   9*gyz(i-2,j-3,k) +   45*gyz(i-1,j-3,k) -   45*gyz(i+1,j-3,k) +   9*gyz(i+2,j-3,k) -    gyz(i+3,j-3,k) ) * odxdy3600
      d2_gd(3,3,1,2) = (    -gzz(i-3,j+3,k) +   9*gzz(i-2,j+3,k) -   45*gzz(i-1,j+3,k) +   45*gzz(i+1,j+3,k) -   9*gzz(i+2,j+3,k) +    gzz(i+3,j+3,k) &
                        +  9*gzz(i-3,j+2,k) -  81*gzz(i-2,j+2,k) +  405*gzz(i-1,j+2,k) -  405*gzz(i+1,j+2,k) +  81*gzz(i+2,j+2,k) -  9*gzz(i+3,j+2,k) &
                        - 45*gzz(i-3,j+1,k) + 405*gzz(i-2,j+1,k) - 2025*gzz(i-1,j+1,k) + 2025*gzz(i+1,j+1,k) - 405*gzz(i+2,j+1,k) + 45*gzz(i+3,j+1,k) &
                        + 45*gzz(i-3,j-1,k) - 405*gzz(i-2,j-1,k) + 2025*gzz(i-1,j-1,k) - 2025*gzz(i+1,j-1,k) + 405*gzz(i+2,j-1,k) - 45*gzz(i+3,j-1,k) &
                        -  9*gzz(i-3,j-2,k) +  81*gzz(i-2,j-2,k) -  405*gzz(i-1,j-2,k) +  405*gzz(i+1,j-2,k) -  81*gzz(i+2,j-2,k) +  9*gzz(i+3,j-2,k) &
                        +    gzz(i-3,j-3,k) -   9*gzz(i-2,j-3,k) +   45*gzz(i-1,j-3,k) -   45*gzz(i+1,j-3,k) +   9*gzz(i+2,j-3,k) -    gzz(i+3,j-3,k) ) * odxdy3600

      d2_gd(1,1,1,3) = (    -gxx(i-3,j,k+3) +   9*gxx(i-2,j,k+3) -   45*gxx(i-1,j,k+3) +   45*gxx(i+1,j,k+3) -   9*gxx(i+2,j,k+3) +    gxx(i+3,j,k+3) &
                        +  9*gxx(i-3,j,k+2) -  81*gxx(i-2,j,k+2) +  405*gxx(i-1,j,k+2) -  405*gxx(i+1,j,k+2) +  81*gxx(i+2,j,k+2) -  9*gxx(i+3,j,k+2) &
                        - 45*gxx(i-3,j,k+1) + 405*gxx(i-2,j,k+1) - 2025*gxx(i-1,j,k+1) + 2025*gxx(i+1,j,k+1) - 405*gxx(i+2,j,k+1) + 45*gxx(i+3,j,k+1) &
                        + 45*gxx(i-3,j,k-1) - 405*gxx(i-2,j,k-1) + 2025*gxx(i-1,j,k-1) - 2025*gxx(i+1,j,k-1) + 405*gxx(i+2,j,k-1) - 45*gxx(i+3,j,k-1) &
                        -  9*gxx(i-3,j,k-2) +  81*gxx(i-2,j,k-2) -  405*gxx(i-1,j,k-2) +  405*gxx(i+1,j,k-2) -  81*gxx(i+2,j,k-2) +  9*gxx(i+3,j,k-2) &
                        +    gxx(i-3,j,k-3) -   9*gxx(i-2,j,k-3) +   45*gxx(i-1,j,k-3) -   45*gxx(i+1,j,k-3) +   9*gxx(i+2,j,k-3) -    gxx(i+3,j,k-3) ) * odxdz3600
      d2_gd(1,2,1,3) = (    -gxy(i-3,j,k+3) +   9*gxy(i-2,j,k+3) -   45*gxy(i-1,j,k+3) +   45*gxy(i+1,j,k+3) -   9*gxy(i+2,j,k+3) +    gxy(i+3,j,k+3) &
                        +  9*gxy(i-3,j,k+2) -  81*gxy(i-2,j,k+2) +  405*gxy(i-1,j,k+2) -  405*gxy(i+1,j,k+2) +  81*gxy(i+2,j,k+2) -  9*gxy(i+3,j,k+2) &
                          - 45*gxy(i-3,j,k+1) + 405*gxy(i-2,j,k+1) - 2025*gxy(i-1,j,k+1) + 2025*gxy(i+1,j,k+1) - 405*gxy(i+2,j,k+1) + 45*gxy(i+3,j,k+1) &
                        + 45*gxy(i-3,j,k-1) - 405*gxy(i-2,j,k-1) + 2025*gxy(i-1,j,k-1) - 2025*gxy(i+1,j,k-1) + 405*gxy(i+2,j,k-1) - 45*gxy(i+3,j,k-1) &
                        -  9*gxy(i-3,j,k-2) +  81*gxy(i-2,j,k-2) -  405*gxy(i-1,j,k-2) +  405*gxy(i+1,j,k-2) -  81*gxy(i+2,j,k-2) +  9*gxy(i+3,j,k-2) &
                        +    gxy(i-3,j,k-3) -   9*gxy(i-2,j,k-3) +   45*gxy(i-1,j,k-3) -   45*gxy(i+1,j,k-3) +   9*gxy(i+2,j,k-3) -    gxy(i+3,j,k-3) ) * odxdz3600
      d2_gd(1,3,1,3) = (    -gxz(i-3,j,k+3) +   9*gxz(i-2,j,k+3) -   45*gxz(i-1,j,k+3) +   45*gxz(i+1,j,k+3) -   9*gxz(i+2,j,k+3) +    gxz(i+3,j,k+3) &
                        +  9*gxz(i-3,j,k+2) -  81*gxz(i-2,j,k+2) +  405*gxz(i-1,j,k+2) -  405*gxz(i+1,j,k+2) +  81*gxz(i+2,j,k+2) -  9*gxz(i+3,j,k+2) &
                        - 45*gxz(i-3,j,k+1) + 405*gxz(i-2,j,k+1) - 2025*gxz(i-1,j,k+1) + 2025*gxz(i+1,j,k+1) - 405*gxz(i+2,j,k+1) + 45*gxz(i+3,j,k+1) &
                        + 45*gxz(i-3,j,k-1) - 405*gxz(i-2,j,k-1) + 2025*gxz(i-1,j,k-1) - 2025*gxz(i+1,j,k-1) + 405*gxz(i+2,j,k-1) - 45*gxz(i+3,j,k-1) &
                        -  9*gxz(i-3,j,k-2) +  81*gxz(i-2,j,k-2) -  405*gxz(i-1,j,k-2) +  405*gxz(i+1,j,k-2) -  81*gxz(i+2,j,k-2) +  9*gxz(i+3,j,k-2) &
                        +    gxz(i-3,j,k-3) -   9*gxz(i-2,j,k-3) +   45*gxz(i-1,j,k-3) -   45*gxz(i+1,j,k-3) +   9*gxz(i+2,j,k-3) -    gxz(i+3,j,k-3) ) * odxdz3600
      d2_gd(2,2,1,3) = (    -gyy(i-3,j,k+3) +   9*gyy(i-2,j,k+3) -   45*gyy(i-1,j,k+3) +   45*gyy(i+1,j,k+3) -   9*gyy(i+2,j,k+3) +    gyy(i+3,j,k+3) &
                        +  9*gyy(i-3,j,k+2) -  81*gyy(i-2,j,k+2) +  405*gyy(i-1,j,k+2) -  405*gyy(i+1,j,k+2) +  81*gyy(i+2,j,k+2) -  9*gyy(i+3,j,k+2) &
                        - 45*gyy(i-3,j,k+1) + 405*gyy(i-2,j,k+1) - 2025*gyy(i-1,j,k+1) + 2025*gyy(i+1,j,k+1) - 405*gyy(i+2,j,k+1) + 45*gyy(i+3,j,k+1) &
                        + 45*gyy(i-3,j,k-1) - 405*gyy(i-2,j,k-1) + 2025*gyy(i-1,j,k-1) - 2025*gyy(i+1,j,k-1) + 405*gyy(i+2,j,k-1) - 45*gyy(i+3,j,k-1) &
                        -  9*gyy(i-3,j,k-2) +  81*gyy(i-2,j,k-2) -  405*gyy(i-1,j,k-2) +  405*gyy(i+1,j,k-2) -  81*gyy(i+2,j,k-2) +  9*gyy(i+3,j,k-2) &
                        +    gyy(i-3,j,k-3) -   9*gyy(i-2,j,k-3) +   45*gyy(i-1,j,k-3) -   45*gyy(i+1,j,k-3) +   9*gyy(i+2,j,k-3) -    gyy(i+3,j,k-3) ) * odxdz3600
      d2_gd(2,3,1,3) = (    -gyz(i-3,j,k+3) +   9*gyz(i-2,j,k+3) -   45*gyz(i-1,j,k+3) +   45*gyz(i+1,j,k+3) -   9*gyz(i+2,j,k+3) +    gyz(i+3,j,k+3) &
                        +  9*gyz(i-3,j,k+2) -  81*gyz(i-2,j,k+2) +  405*gyz(i-1,j,k+2) -  405*gyz(i+1,j,k+2) +  81*gyz(i+2,j,k+2) -  9*gyz(i+3,j,k+2) &
                        - 45*gyz(i-3,j,k+1) + 405*gyz(i-2,j,k+1) - 2025*gyz(i-1,j,k+1) + 2025*gyz(i+1,j,k+1) - 405*gyz(i+2,j,k+1) + 45*gyz(i+3,j,k+1) &
                        + 45*gyz(i-3,j,k-1) - 405*gyz(i-2,j,k-1) + 2025*gyz(i-1,j,k-1) - 2025*gyz(i+1,j,k-1) + 405*gyz(i+2,j,k-1) - 45*gyz(i+3,j,k-1) &
                        -  9*gyz(i-3,j,k-2) +  81*gyz(i-2,j,k-2) -  405*gyz(i-1,j,k-2) +  405*gyz(i+1,j,k-2) -  81*gyz(i+2,j,k-2) +  9*gyz(i+3,j,k-2) &
                        +    gyz(i-3,j,k-3) -   9*gyz(i-2,j,k-3) +   45*gyz(i-1,j,k-3) -   45*gyz(i+1,j,k-3) +   9*gyz(i+2,j,k-3) -    gyz(i+3,j,k-3) ) * odxdz3600
      d2_gd(3,3,1,3) = (    -gzz(i-3,j,k+3) +   9*gzz(i-2,j,k+3) -   45*gzz(i-1,j,k+3) +   45*gzz(i+1,j,k+3) -   9*gzz(i+2,j,k+3) +    gzz(i+3,j,k+3) &
                        +  9*gzz(i-3,j,k+2) -  81*gzz(i-2,j,k+2) +  405*gzz(i-1,j,k+2) -  405*gzz(i+1,j,k+2) +  81*gzz(i+2,j,k+2) -  9*gzz(i+3,j,k+2) &
                        - 45*gzz(i-3,j,k+1) + 405*gzz(i-2,j,k+1) - 2025*gzz(i-1,j,k+1) + 2025*gzz(i+1,j,k+1) - 405*gzz(i+2,j,k+1) + 45*gzz(i+3,j,k+1) &
                        + 45*gzz(i-3,j,k-1) - 405*gzz(i-2,j,k-1) + 2025*gzz(i-1,j,k-1) - 2025*gzz(i+1,j,k-1) + 405*gzz(i+2,j,k-1) - 45*gzz(i+3,j,k-1) &
                        -  9*gzz(i-3,j,k-2) +  81*gzz(i-2,j,k-2) -  405*gzz(i-1,j,k-2) +  405*gzz(i+1,j,k-2) -  81*gzz(i+2,j,k-2) +  9*gzz(i+3,j,k-2) &
                        +    gzz(i-3,j,k-3) -   9*gzz(i-2,j,k-3) +   45*gzz(i-1,j,k-3) -   45*gzz(i+1,j,k-3) +   9*gzz(i+2,j,k-3) -    gzz(i+3,j,k-3) ) * odxdz3600

      d2_gd(1,1,2,3) = (    -gxx(i,j-3,k+3) +   9*gxx(i,j-2,k+3) -   45*gxx(i,j-1,k+3) +   45*gxx(i,j+1,k+3) -   9*gxx(i,j+2,k+3) +    gxx(i,j+3,k+3) &
                        +  9*gxx(i,j-3,k+2) -  81*gxx(i,j-2,k+2) +  405*gxx(i,j-1,k+2) -  405*gxx(i,j+1,k+2) +  81*gxx(i,j+2,k+2) -  9*gxx(i,j+3,k+2) &
                        - 45*gxx(i,j-3,k+1) + 405*gxx(i,j-2,k+1) - 2025*gxx(i,j-1,k+1) + 2025*gxx(i,j+1,k+1) - 405*gxx(i,j+2,k+1) + 45*gxx(i,j+3,k+1) &
                        + 45*gxx(i,j-3,k-1) - 405*gxx(i,j-2,k-1) + 2025*gxx(i,j-1,k-1) - 2025*gxx(i,j+1,k-1) + 405*gxx(i,j+2,k-1) - 45*gxx(i,j+3,k-1) &
                        -  9*gxx(i,j-3,k-2) +  81*gxx(i,j-2,k-2) -  405*gxx(i,j-1,k-2) +  405*gxx(i,j+1,k-2) -  81*gxx(i,j+2,k-2) +  9*gxx(i,j+3,k-2) &
                        +    gxx(i,j-3,k-3) -   9*gxx(i,j-2,k-3) +   45*gxx(i,j-1,k-3) -   45*gxx(i,j+1,k-3) +   9*gxx(i,j+2,k-3) -    gxx(i,j+3,k-3) ) * odydz3600
      d2_gd(1,2,2,3) = (    -gxy(i,j-3,k+3) +   9*gxy(i,j-2,k+3) -   45*gxy(i,j-1,k+3) +   45*gxy(i,j+1,k+3) -   9*gxy(i,j+2,k+3) +    gxy(i,j+3,k+3) &
                        +  9*gxy(i,j-3,k+2) -  81*gxy(i,j-2,k+2) +  405*gxy(i,j-1,k+2) -  405*gxy(i,j+1,k+2) +  81*gxy(i,j+2,k+2) -  9*gxy(i,j+3,k+2) &
                        - 45*gxy(i,j-3,k+1) + 405*gxy(i,j-2,k+1) - 2025*gxy(i,j-1,k+1) + 2025*gxy(i,j+1,k+1) - 405*gxy(i,j+2,k+1) + 45*gxy(i,j+3,k+1) &
                        + 45*gxy(i,j-3,k-1) - 405*gxy(i,j-2,k-1) + 2025*gxy(i,j-1,k-1) - 2025*gxy(i,j+1,k-1) + 405*gxy(i,j+2,k-1) - 45*gxy(i,j+3,k-1) &
                        -  9*gxy(i,j-3,k-2) +  81*gxy(i,j-2,k-2) -  405*gxy(i,j-1,k-2) +  405*gxy(i,j+1,k-2) -  81*gxy(i,j+2,k-2) +  9*gxy(i,j+3,k-2) &
                        +    gxy(i,j-3,k-3) -   9*gxy(i,j-2,k-3) +   45*gxy(i,j-1,k-3) -   45*gxy(i,j+1,k-3) +   9*gxy(i,j+2,k-3) -    gxy(i,j+3,k-3) ) * odydz3600
      d2_gd(1,3,2,3) = (    -gxz(i,j-3,k+3) +   9*gxz(i,j-2,k+3) -   45*gxz(i,j-1,k+3) +   45*gxz(i,j+1,k+3) -   9*gxz(i,j+2,k+3) +    gxz(i,j+3,k+3) &
                        +  9*gxz(i,j-3,k+2) -  81*gxz(i,j-2,k+2) +  405*gxz(i,j-1,k+2) -  405*gxz(i,j+1,k+2) +  81*gxz(i,j+2,k+2) -  9*gxz(i,j+3,k+2) &
                        - 45*gxz(i,j-3,k+1) + 405*gxz(i,j-2,k+1) - 2025*gxz(i,j-1,k+1) + 2025*gxz(i,j+1,k+1) - 405*gxz(i,j+2,k+1) + 45*gxz(i,j+3,k+1) &
                        + 45*gxz(i,j-3,k-1) - 405*gxz(i,j-2,k-1) + 2025*gxz(i,j-1,k-1) - 2025*gxz(i,j+1,k-1) + 405*gxz(i,j+2,k-1) - 45*gxz(i,j+3,k-1) &
                        -  9*gxz(i,j-3,k-2) +  81*gxz(i,j-2,k-2) -  405*gxz(i,j-1,k-2) +  405*gxz(i,j+1,k-2) -  81*gxz(i,j+2,k-2) +  9*gxz(i,j+3,k-2) &
                        +    gxz(i,j-3,k-3) -   9*gxz(i,j-2,k-3) +   45*gxz(i,j-1,k-3) -   45*gxz(i,j+1,k-3) +   9*gxz(i,j+2,k-3) -    gxz(i,j+3,k-3) ) * odydz3600
      d2_gd(2,2,2,3) = (    -gyy(i,j-3,k+3) +   9*gyy(i,j-2,k+3) -   45*gyy(i,j-1,k+3) +   45*gyy(i,j+1,k+3) -   9*gyy(i,j+2,k+3) +    gyy(i,j+3,k+3) &
                        +  9*gyy(i,j-3,k+2) -  81*gyy(i,j-2,k+2) +  405*gyy(i,j-1,k+2) -  405*gyy(i,j+1,k+2) +  81*gyy(i,j+2,k+2) -  9*gyy(i,j+3,k+2) &
                        - 45*gyy(i,j-3,k+1) + 405*gyy(i,j-2,k+1) - 2025*gyy(i,j-1,k+1) + 2025*gyy(i,j+1,k+1) - 405*gyy(i,j+2,k+1) + 45*gyy(i,j+3,k+1) &
                        + 45*gyy(i,j-3,k-1) - 405*gyy(i,j-2,k-1) + 2025*gyy(i,j-1,k-1) - 2025*gyy(i,j+1,k-1) + 405*gyy(i,j+2,k-1) - 45*gyy(i,j+3,k-1) &
                        -  9*gyy(i,j-3,k-2) +  81*gyy(i,j-2,k-2) -  405*gyy(i,j-1,k-2) +  405*gyy(i,j+1,k-2) -  81*gyy(i,j+2,k-2) +  9*gyy(i,j+3,k-2) &
                        +    gyy(i,j-3,k-3) -   9*gyy(i,j-2,k-3) +   45*gyy(i,j-1,k-3) -   45*gyy(i,j+1,k-3) +   9*gyy(i,j+2,k-3) -    gyy(i,j+3,k-3) ) * odydz3600
      d2_gd(2,3,2,3) = (    -gyz(i,j-3,k+3) +   9*gyz(i,j-2,k+3) -   45*gyz(i,j-1,k+3) +   45*gyz(i,j+1,k+3) -   9*gyz(i,j+2,k+3) +    gyz(i,j+3,k+3) &
                        +  9*gyz(i,j-3,k+2) -  81*gyz(i,j-2,k+2) +  405*gyz(i,j-1,k+2) -  405*gyz(i,j+1,k+2) +  81*gyz(i,j+2,k+2) -  9*gyz(i,j+3,k+2) &
                        - 45*gyz(i,j-3,k+1) + 405*gyz(i,j-2,k+1) - 2025*gyz(i,j-1,k+1) + 2025*gyz(i,j+1,k+1) - 405*gyz(i,j+2,k+1) + 45*gyz(i,j+3,k+1) &
                        + 45*gyz(i,j-3,k-1) - 405*gyz(i,j-2,k-1) + 2025*gyz(i,j-1,k-1) - 2025*gyz(i,j+1,k-1) + 405*gyz(i,j+2,k-1) - 45*gyz(i,j+3,k-1) &
                        -  9*gyz(i,j-3,k-2) +  81*gyz(i,j-2,k-2) -  405*gyz(i,j-1,k-2) +  405*gyz(i,j+1,k-2) -  81*gyz(i,j+2,k-2) +  9*gyz(i,j+3,k-2) &
                        +    gyz(i,j-3,k-3) -   9*gyz(i,j-2,k-3) +   45*gyz(i,j-1,k-3) -   45*gyz(i,j+1,k-3) +   9*gyz(i,j+2,k-3) -    gyz(i,j+3,k-3) ) * odydz3600
      d2_gd(3,3,2,3) = (    -gzz(i,j-3,k+3) +   9*gzz(i,j-2,k+3) -   45*gzz(i,j-1,k+3) +   45*gzz(i,j+1,k+3) -   9*gzz(i,j+2,k+3) +    gzz(i,j+3,k+3) &
                        +  9*gzz(i,j-3,k+2) -  81*gzz(i,j-2,k+2) +  405*gzz(i,j-1,k+2) -  405*gzz(i,j+1,k+2) +  81*gzz(i,j+2,k+2) -  9*gzz(i,j+3,k+2) &
                        - 45*gzz(i,j-3,k+1) + 405*gzz(i,j-2,k+1) - 2025*gzz(i,j-1,k+1) + 2025*gzz(i,j+1,k+1) - 405*gzz(i,j+2,k+1) + 45*gzz(i,j+3,k+1) &
                        + 45*gzz(i,j-3,k-1) - 405*gzz(i,j-2,k-1) + 2025*gzz(i,j-1,k-1) - 2025*gzz(i,j+1,k-1) + 405*gzz(i,j+2,k-1) - 45*gzz(i,j+3,k-1) &
                        -  9*gzz(i,j-3,k-2) +  81*gzz(i,j-2,k-2) -  405*gzz(i,j-1,k-2) +  405*gzz(i,j+1,k-2) -  81*gzz(i,j+2,k-2) +  9*gzz(i,j+3,k-2) &
                        +    gzz(i,j-3,k-3) -   9*gzz(i,j-2,k-3) +   45*gzz(i,j-1,k-3) -   45*gzz(i,j+1,k-3) +   9*gzz(i,j+2,k-3) -    gzz(i,j+3,k-3) ) * odydz3600

      d2_gd(1,1,2,1) = d2_gd(1,1,1,2)
      d2_gd(1,2,2,1) = d2_gd(1,2,1,2)
      d2_gd(1,3,2,1) = d2_gd(1,3,1,2)
      d2_gd(2,2,2,1) = d2_gd(2,2,1,2)
      d2_gd(2,3,2,1) = d2_gd(2,3,1,2)
      d2_gd(3,3,2,1) = d2_gd(3,3,1,2)

      d2_gd(1,1,3,1) = d2_gd(1,1,1,3)
      d2_gd(1,2,3,1) = d2_gd(1,2,1,3)
      d2_gd(1,3,3,1) = d2_gd(1,3,1,3)
      d2_gd(2,2,3,1) = d2_gd(2,2,1,3)
      d2_gd(2,3,3,1) = d2_gd(2,3,1,3)
      d2_gd(3,3,3,1) = d2_gd(3,3,1,3)

      d2_gd(1,1,3,2) = d2_gd(1,1,2,3)
      d2_gd(1,2,3,2) = d2_gd(1,2,2,3)
      d2_gd(1,3,3,2) = d2_gd(1,3,2,3)
      d2_gd(2,2,3,2) = d2_gd(2,2,2,3)
      d2_gd(2,3,3,2) = d2_gd(2,3,2,3)
      d2_gd(3,3,3,2) = d2_gd(3,3,2,3)

      d2_gd(2,1,:,:) = d2_gd(1,2,:,:)
      d2_gd(3,1,:,:) = d2_gd(1,3,:,:)
      d2_gd(3,2,:,:) = d2_gd(2,3,:,:)

      !--------------------------------------------

    else if ( NP_order == 4 ) then

      ! Fourth order derivatives

      !-------------- Centered 1st derivatives ----

      ! d1_gd(3,3,3)
      d1_gd(1,1,1) = (   -gxx(i+2,j,k) + 8*gxx(i+1,j,k)                      &
                      - 8*gxx(i-1,j,k) +   gxx(i-2,j,k) ) / dx12
      d1_gd(1,2,1) = (   -gxy(i+2,j,k) + 8*gxy(i+1,j,k)                      &
                      - 8*gxy(i-1,j,k) +   gxy(i-2,j,k) ) / dx12
      d1_gd(1,3,1) = (   -gxz(i+2,j,k) + 8*gxz(i+1,j,k)                      &
                      - 8*gxz(i-1,j,k) +   gxz(i-2,j,k) ) / dx12
      d1_gd(2,2,1) = (   -gyy(i+2,j,k) + 8*gyy(i+1,j,k)                      &
                      - 8*gyy(i-1,j,k) +   gyy(i-2,j,k) ) / dx12
      d1_gd(2,3,1) = (   -gyz(i+2,j,k) + 8*gyz(i+1,j,k)                      &
                      - 8*gyz(i-1,j,k) +   gyz(i-2,j,k) ) / dx12
      d1_gd(3,3,1) = (   -gzz(i+2,j,k) + 8*gzz(i+1,j,k)                      &
                      - 8*gzz(i-1,j,k) +   gzz(i-2,j,k) ) / dx12

      d1_gd(1,1,2) = (   -gxx(i,j+2,k) + 8*gxx(i,j+1,k)                      &
                      - 8*gxx(i,j-1,k) +   gxx(i,j-2,k) ) / dy12
      d1_gd(1,2,2) = (   -gxy(i,j+2,k) + 8*gxy(i,j+1,k)                      &
                      - 8*gxy(i,j-1,k) +   gxy(i,j-2,k) ) / dy12
      d1_gd(1,3,2) = (   -gxz(i,j+2,k) + 8*gxz(i,j+1,k)                      &
                      - 8*gxz(i,j-1,k) +   gxz(i,j-2,k) ) / dy12
      d1_gd(2,2,2) = (   -gyy(i,j+2,k) + 8*gyy(i,j+1,k)                      &
                      - 8*gyy(i,j-1,k) +   gyy(i,j-2,k) ) / dy12
      d1_gd(2,3,2) = (   -gyz(i,j+2,k) + 8*gyz(i,j+1,k)                      &
                      - 8*gyz(i,j-1,k) +   gyz(i,j-2,k) ) / dy12
      d1_gd(3,3,2) = (   -gzz(i,j+2,k) + 8*gzz(i,j+1,k)                      &
                      - 8*gzz(i,j-1,k) +   gzz(i,j-2,k) ) / dy12

      d1_gd(1,1,3) = (   -gxx(i,j,k+2) + 8*gxx(i,j,k+1)                      &
                      - 8*gxx(i,j,k-1) +   gxx(i,j,k-2) ) / dz12
      d1_gd(1,2,3) = (   -gxy(i,j,k+2) + 8*gxy(i,j,k+1)                      &
                      - 8*gxy(i,j,k-1) +   gxy(i,j,k-2) ) / dz12
      d1_gd(1,3,3) = (   -gxz(i,j,k+2) + 8*gxz(i,j,k+1)                      &
                      - 8*gxz(i,j,k-1) +   gxz(i,j,k-2) ) / dz12
      d1_gd(2,2,3) = (   -gyy(i,j,k+2) + 8*gyy(i,j,k+1)                      &
                      - 8*gyy(i,j,k-1) +   gyy(i,j,k-2) ) / dz12
      d1_gd(2,3,3) = (   -gyz(i,j,k+2) + 8*gyz(i,j,k+1)                      &
                      - 8*gyz(i,j,k-1) +   gyz(i,j,k-2) ) / dz12
      d1_gd(3,3,3) = (   -gzz(i,j,k+2) + 8*gzz(i,j,k+1)                      &
                      - 8*gzz(i,j,k-1) +   gzz(i,j,k-2) ) / dz12

      d1_gd(2,1,:) = d1_gd(1,2,:)
      d1_gd(3,1,:) = d1_gd(1,3,:)
      d1_gd(3,2,:) = d1_gd(2,3,:)


      ! d1_kd(3,3,3)
      d1_kd(1,1,1) = (   -kxx(i+2,j,k) + 8*kxx(i+1,j,k)                      &
                      - 8*kxx(i-1,j,k) +   kxx(i-2,j,k) ) / dx12
      d1_kd(1,2,1) = (   -kxy(i+2,j,k) + 8*kxy(i+1,j,k)                      &
                      - 8*kxy(i-1,j,k) +   kxy(i-2,j,k) ) / dx12
      d1_kd(1,3,1) = (   -kxz(i+2,j,k) + 8*kxz(i+1,j,k)                      &
                      - 8*kxz(i-1,j,k) +   kxz(i-2,j,k) ) / dx12
      d1_kd(2,2,1) = (   -kyy(i+2,j,k) + 8*kyy(i+1,j,k)                      &
                      - 8*kyy(i-1,j,k) +   kyy(i-2,j,k) ) / dx12
      d1_kd(2,3,1) = (   -kyz(i+2,j,k) + 8*kyz(i+1,j,k)                      &
                      - 8*kyz(i-1,j,k) +   kyz(i-2,j,k) ) / dx12
      d1_kd(3,3,1) = (   -kzz(i+2,j,k) + 8*kzz(i+1,j,k)                      &
                      - 8*kzz(i-1,j,k) +   kzz(i-2,j,k) ) / dx12

      d1_kd(1,1,2) = (   -kxx(i,j+2,k) + 8*kxx(i,j+1,k)                      &
                      - 8*kxx(i,j-1,k) +   kxx(i,j-2,k) ) / dy12
      d1_kd(1,2,2) = (   -kxy(i,j+2,k) + 8*kxy(i,j+1,k)                      &
                      - 8*kxy(i,j-1,k) +   kxy(i,j-2,k) ) / dy12
      d1_kd(1,3,2) = (   -kxz(i,j+2,k) + 8*kxz(i,j+1,k)                      &
                      - 8*kxz(i,j-1,k) +   kxz(i,j-2,k) ) / dy12
      d1_kd(2,2,2) = (   -kyy(i,j+2,k) + 8*kyy(i,j+1,k)                      &
                      - 8*kyy(i,j-1,k) +   kyy(i,j-2,k) ) / dy12
      d1_kd(2,3,2) = (   -kyz(i,j+2,k) + 8*kyz(i,j+1,k)                      &
                      - 8*kyz(i,j-1,k) +   kyz(i,j-2,k) ) / dy12
      d1_kd(3,3,2) = (   -kzz(i,j+2,k) + 8*kzz(i,j+1,k)                      &
                      - 8*kzz(i,j-1,k) +   kzz(i,j-2,k) ) / dy12

      d1_kd(1,1,3) = (   -kxx(i,j,k+2) + 8*kxx(i,j,k+1)                      &
                      - 8*kxx(i,j,k-1) +   kxx(i,j,k-2) ) / dz12
      d1_kd(1,2,3) = (   -kxy(i,j,k+2) + 8*kxy(i,j,k+1)                      &
                      - 8*kxy(i,j,k-1) +   kxy(i,j,k-2) ) / dz12
      d1_kd(1,3,3) = (   -kxz(i,j,k+2) + 8*kxz(i,j,k+1)                      &
                      - 8*kxz(i,j,k-1) +   kxz(i,j,k-2) ) / dz12
      d1_kd(2,2,3) = (   -kyy(i,j,k+2) + 8*kyy(i,j,k+1)                      &
                      - 8*kyy(i,j,k-1) +   kyy(i,j,k-2) ) / dz12
      d1_kd(2,3,3) = (   -kyz(i,j,k+2) + 8*kyz(i,j,k+1)                      &
                      - 8*kyz(i,j,k-1) +   kyz(i,j,k-2) ) / dz12
      d1_kd(3,3,3) = (   -kzz(i,j,k+2) + 8*kzz(i,j,k+1)                      &
                      - 8*kzz(i,j,k-1) +   kzz(i,j,k-2) ) / dz12

      d1_kd(2,1,:) = d1_kd(1,2,:)
      d1_kd(3,1,:) = d1_kd(1,3,:)
      d1_kd(3,2,:) = d1_kd(2,3,:)


      ! d1_lA removed: B^i is evolved directly in Einstein-Maxwell.

      !--------------------------------------------


      !------------ Centered 2nd derivatives -----

      ! d2_gd(3,3,3,3)
      d2_gd(1,1,1,1) = (   -gxx(i+2,j,k) + 16*gxx(i+1,j,k) - 30*gxx(i,j,k)   &
                       + 16*gxx(i-1,j,k) -    gxx(i-2,j,k) ) / dxsq12
      d2_gd(1,2,1,1) = (   -gxy(i+2,j,k) + 16*gxy(i+1,j,k) - 30*gxy(i,j,k)   &
                       + 16*gxy(i-1,j,k) -    gxy(i-2,j,k) ) / dxsq12
      d2_gd(1,3,1,1) = (   -gxz(i+2,j,k) + 16*gxz(i+1,j,k) - 30*gxz(i,j,k)   &
                       + 16*gxz(i-1,j,k) -    gxz(i-2,j,k) ) / dxsq12
      d2_gd(2,2,1,1) = (   -gyy(i+2,j,k) + 16*gyy(i+1,j,k) - 30*gyy(i,j,k)   &
                       + 16*gyy(i-1,j,k) -    gyy(i-2,j,k) ) / dxsq12
      d2_gd(2,3,1,1) = (   -gyz(i+2,j,k) + 16*gyz(i+1,j,k) - 30*gyz(i,j,k)   &
                       + 16*gyz(i-1,j,k) -    gyz(i-2,j,k) ) / dxsq12
      d2_gd(3,3,1,1) = (   -gzz(i+2,j,k) + 16*gzz(i+1,j,k) - 30*gzz(i,j,k)   &
                       + 16*gzz(i-1,j,k) -    gzz(i-2,j,k) ) / dxsq12

      d2_gd(1,1,2,2) = (   -gxx(i,j+2,k) + 16*gxx(i,j+1,k) - 30*gxx(i,j,k)   &
                       + 16*gxx(i,j-1,k) -    gxx(i,j-2,k) ) / dysq12
      d2_gd(1,2,2,2) = (   -gxy(i,j+2,k) + 16*gxy(i,j+1,k) - 30*gxy(i,j,k)   &
                       + 16*gxy(i,j-1,k) -    gxy(i,j-2,k) ) / dysq12
      d2_gd(1,3,2,2) = (   -gxz(i,j+2,k) + 16*gxz(i,j+1,k) - 30*gxz(i,j,k)   &
                       + 16*gxz(i,j-1,k) -    gxz(i,j-2,k) ) / dysq12
      d2_gd(2,2,2,2) = (   -gyy(i,j+2,k) + 16*gyy(i,j+1,k) - 30*gyy(i,j,k)   &
                       + 16*gyy(i,j-1,k) -    gyy(i,j-2,k) ) / dysq12
      d2_gd(2,3,2,2) = (   -gyz(i,j+2,k) + 16*gyz(i,j+1,k) - 30*gyz(i,j,k)   &
                       + 16*gyz(i,j-1,k) -    gyz(i,j-2,k) ) / dysq12
      d2_gd(3,3,2,2) = (   -gzz(i,j+2,k) + 16*gzz(i,j+1,k) - 30*gzz(i,j,k)   &
                       + 16*gzz(i,j-1,k) -    gzz(i,j-2,k) ) / dysq12

      d2_gd(1,1,3,3) = (   -gxx(i,j,k+2) + 16*gxx(i,j,k+1) - 30*gxx(i,j,k)   &
                       + 16*gxx(i,j,k-1) -    gxx(i,j,k-2) ) / dzsq12
      d2_gd(1,2,3,3) = (   -gxy(i,j,k+2) + 16*gxy(i,j,k+1) - 30*gxy(i,j,k)   &
                       + 16*gxy(i,j,k-1) -    gxy(i,j,k-2) ) / dzsq12
      d2_gd(1,3,3,3) = (   -gxz(i,j,k+2) + 16*gxz(i,j,k+1) - 30*gxz(i,j,k)   &
                       + 16*gxz(i,j,k-1) -    gxz(i,j,k-2) ) / dzsq12
      d2_gd(2,2,3,3) = (   -gyy(i,j,k+2) + 16*gyy(i,j,k+1) - 30*gyy(i,j,k)   &
                       + 16*gyy(i,j,k-1) -    gyy(i,j,k-2) ) / dzsq12
      d2_gd(2,3,3,3) = (   -gyz(i,j,k+2) + 16*gyz(i,j,k+1) - 30*gyz(i,j,k)   &
                       + 16*gyz(i,j,k-1) -    gyz(i,j,k-2) ) / dzsq12
      d2_gd(3,3,3,3) = (   -gzz(i,j,k+2) + 16*gzz(i,j,k+1) - 30*gzz(i,j,k)   &
                       + 16*gzz(i,j,k-1) -    gzz(i,j,k-2) ) / dzsq12

      d2_gd(1,1,1,2) = (   -gxx(i-2,j+2,k) +  8*gxx(i-1,j+2,k) -  8*gxx(i+1,j+2,k) +   gxx(i+2,j+2,k)   &
                        + 8*gxx(i-2,j+1,k) - 64*gxx(i-1,j+1,k) + 64*gxx(i+1,j+1,k) - 8*gxx(i+2,j+1,k)   &
                        - 8*gxx(i-2,j-1,k) + 64*gxx(i-1,j-1,k) - 64*gxx(i+1,j-1,k) + 8*gxx(i+2,j-1,k)   &
                        +   gxx(i-2,j-2,k) -  8*gxx(i-1,j-2,k) +  8*gxx(i+1,j-2,k) -   gxx(i+2,j-2,k) ) / dxdy144
      d2_gd(1,2,1,2) = (   -gxy(i-2,j+2,k) +  8*gxy(i-1,j+2,k) -  8*gxy(i+1,j+2,k) +   gxy(i+2,j+2,k)   &
                        + 8*gxy(i-2,j+1,k) - 64*gxy(i-1,j+1,k) + 64*gxy(i+1,j+1,k) - 8*gxy(i+2,j+1,k)   &
                        - 8*gxy(i-2,j-1,k) + 64*gxy(i-1,j-1,k) - 64*gxy(i+1,j-1,k) + 8*gxy(i+2,j-1,k)   &
                        +   gxy(i-2,j-2,k) -  8*gxy(i-1,j-2,k) +  8*gxy(i+1,j-2,k) -   gxy(i+2,j-2,k) ) / dxdy144
      d2_gd(1,3,1,2) = (   -gxz(i-2,j+2,k) +  8*gxz(i-1,j+2,k) -  8*gxz(i+1,j+2,k) +   gxz(i+2,j+2,k)   &
                        + 8*gxz(i-2,j+1,k) - 64*gxz(i-1,j+1,k) + 64*gxz(i+1,j+1,k) - 8*gxz(i+2,j+1,k)   &
                        - 8*gxz(i-2,j-1,k) + 64*gxz(i-1,j-1,k) - 64*gxz(i+1,j-1,k) + 8*gxz(i+2,j-1,k)   &
                        +   gxz(i-2,j-2,k) -  8*gxz(i-1,j-2,k) +  8*gxz(i+1,j-2,k) -   gxz(i+2,j-2,k) ) / dxdy144
      d2_gd(2,2,1,2) = (   -gyy(i-2,j+2,k) +  8*gyy(i-1,j+2,k) -  8*gyy(i+1,j+2,k) +   gyy(i+2,j+2,k)   &
                        + 8*gyy(i-2,j+1,k) - 64*gyy(i-1,j+1,k) + 64*gyy(i+1,j+1,k) - 8*gyy(i+2,j+1,k)   &
                        - 8*gyy(i-2,j-1,k) + 64*gyy(i-1,j-1,k) - 64*gyy(i+1,j-1,k) + 8*gyy(i+2,j-1,k)   &
                        +   gyy(i-2,j-2,k) -  8*gyy(i-1,j-2,k) +  8*gyy(i+1,j-2,k) -   gyy(i+2,j-2,k) ) / dxdy144
      d2_gd(2,3,1,2) = (   -gyz(i-2,j+2,k) +  8*gyz(i-1,j+2,k) -  8*gyz(i+1,j+2,k) +   gyz(i+2,j+2,k)   &
                        + 8*gyz(i-2,j+1,k) - 64*gyz(i-1,j+1,k) + 64*gyz(i+1,j+1,k) - 8*gyz(i+2,j+1,k)   &
                        - 8*gyz(i-2,j-1,k) + 64*gyz(i-1,j-1,k) - 64*gyz(i+1,j-1,k) + 8*gyz(i+2,j-1,k)   &
                        +   gyz(i-2,j-2,k) -  8*gyz(i-1,j-2,k) +  8*gyz(i+1,j-2,k) -   gyz(i+2,j-2,k) ) / dxdy144
      d2_gd(3,3,1,2) = (   -gzz(i-2,j+2,k) +  8*gzz(i-1,j+2,k) -  8*gzz(i+1,j+2,k) +   gzz(i+2,j+2,k)   &
                        + 8*gzz(i-2,j+1,k) - 64*gzz(i-1,j+1,k) + 64*gzz(i+1,j+1,k) - 8*gzz(i+2,j+1,k)   &
                        - 8*gzz(i-2,j-1,k) + 64*gzz(i-1,j-1,k) - 64*gzz(i+1,j-1,k) + 8*gzz(i+2,j-1,k)   &
                        +   gzz(i-2,j-2,k) -  8*gzz(i-1,j-2,k) +  8*gzz(i+1,j-2,k) -   gzz(i+2,j-2,k) ) / dxdy144

      d2_gd(1,1,1,3) = (   -gxx(i-2,j,k+2) +  8*gxx(i-1,j,k+2) -  8*gxx(i+1,j,k+2) +   gxx(i+2,j,k+2)   &
                        + 8*gxx(i-2,j,k+1) - 64*gxx(i-1,j,k+1) + 64*gxx(i+1,j,k+1) - 8*gxx(i+2,j,k+1)   &
                        - 8*gxx(i-2,j,k-1) + 64*gxx(i-1,j,k-1) - 64*gxx(i+1,j,k-1) + 8*gxx(i+2,j,k-1)   &
                        +   gxx(i-2,j,k-2) -  8*gxx(i-1,j,k-2) +  8*gxx(i+1,j,k-2) -   gxx(i+2,j,k-2) ) / dxdz144
      d2_gd(1,2,1,3) = (   -gxy(i-2,j,k+2) +  8*gxy(i-1,j,k+2) -  8*gxy(i+1,j,k+2) +   gxy(i+2,j,k+2)   &
                        + 8*gxy(i-2,j,k+1) - 64*gxy(i-1,j,k+1) + 64*gxy(i+1,j,k+1) - 8*gxy(i+2,j,k+1)   &
                        - 8*gxy(i-2,j,k-1) + 64*gxy(i-1,j,k-1) - 64*gxy(i+1,j,k-1) + 8*gxy(i+2,j,k-1)   &
                        +   gxy(i-2,j,k-2) -  8*gxy(i-1,j,k-2) +  8*gxy(i+1,j,k-2) -   gxy(i+2,j,k-2) ) / dxdz144
      d2_gd(1,3,1,3) = (   -gxz(i-2,j,k+2) +  8*gxz(i-1,j,k+2) -  8*gxz(i+1,j,k+2) +   gxz(i+2,j,k+2)   &
                        + 8*gxz(i-2,j,k+1) - 64*gxz(i-1,j,k+1) + 64*gxz(i+1,j,k+1) - 8*gxz(i+2,j,k+1)   &
                        - 8*gxz(i-2,j,k-1) + 64*gxz(i-1,j,k-1) - 64*gxz(i+1,j,k-1) + 8*gxz(i+2,j,k-1)   &
                        +   gxz(i-2,j,k-2) -  8*gxz(i-1,j,k-2) +  8*gxz(i+1,j,k-2) -   gxz(i+2,j,k-2) ) / dxdz144
      d2_gd(2,2,1,3) = (   -gyy(i-2,j,k+2) +  8*gyy(i-1,j,k+2) -  8*gyy(i+1,j,k+2) +   gyy(i+2,j,k+2)   &
                        + 8*gyy(i-2,j,k+1) - 64*gyy(i-1,j,k+1) + 64*gyy(i+1,j,k+1) - 8*gyy(i+2,j,k+1)   &
                        - 8*gyy(i-2,j,k-1) + 64*gyy(i-1,j,k-1) - 64*gyy(i+1,j,k-1) + 8*gyy(i+2,j,k-1)   &
                        +   gyy(i-2,j,k-2) -  8*gyy(i-1,j,k-2) +  8*gyy(i+1,j,k-2) -   gyy(i+2,j,k-2) ) / dxdz144
      d2_gd(2,3,1,3) = (   -gyz(i-2,j,k+2) +  8*gyz(i-1,j,k+2) -  8*gyz(i+1,j,k+2) +   gyz(i+2,j,k+2)   &
                        + 8*gyz(i-2,j,k+1) - 64*gyz(i-1,j,k+1) + 64*gyz(i+1,j,k+1) - 8*gyz(i+2,j,k+1)   &
                        - 8*gyz(i-2,j,k-1) + 64*gyz(i-1,j,k-1) - 64*gyz(i+1,j,k-1) + 8*gyz(i+2,j,k-1)   &
                        +   gyz(i-2,j,k-2) -  8*gyz(i-1,j,k-2) +  8*gyz(i+1,j,k-2) -   gyz(i+2,j,k-2) ) / dxdz144
      d2_gd(3,3,1,3) = (   -gzz(i-2,j,k+2) +  8*gzz(i-1,j,k+2) -  8*gzz(i+1,j,k+2) +   gzz(i+2,j,k+2)   &
                        + 8*gzz(i-2,j,k+1) - 64*gzz(i-1,j,k+1) + 64*gzz(i+1,j,k+1) - 8*gzz(i+2,j,k+1)   &
                        - 8*gzz(i-2,j,k-1) + 64*gzz(i-1,j,k-1) - 64*gzz(i+1,j,k-1) + 8*gzz(i+2,j,k-1)   &
                        +   gzz(i-2,j,k-2) -  8*gzz(i-1,j,k-2) +  8*gzz(i+1,j,k-2) -   gzz(i+2,j,k-2) ) / dxdz144

      d2_gd(1,1,2,3) = (   -gxx(i,j-2,k+2) +  8*gxx(i,j-1,k+2) -  8*gxx(i,j+1,k+2) +   gxx(i,j+2,k+2)   &
                        + 8*gxx(i,j-2,k+1) - 64*gxx(i,j-1,k+1) + 64*gxx(i,j+1,k+1) - 8*gxx(i,j+2,k+1)   &
                        - 8*gxx(i,j-2,k-1) + 64*gxx(i,j-1,k-1) - 64*gxx(i,j+1,k-1) + 8*gxx(i,j+2,k-1)   &
                        +   gxx(i,j-2,k-2) -  8*gxx(i,j-1,k-2) +  8*gxx(i,j+1,k-2) -   gxx(i,j+2,k-2) ) / dydz144
      d2_gd(1,2,2,3) = (   -gxy(i,j-2,k+2) +  8*gxy(i,j-1,k+2) -  8*gxy(i,j+1,k+2) +   gxy(i,j+2,k+2)   &
                        + 8*gxy(i,j-2,k+1) - 64*gxy(i,j-1,k+1) + 64*gxy(i,j+1,k+1) - 8*gxy(i,j+2,k+1)   &
                        - 8*gxy(i,j-2,k-1) + 64*gxy(i,j-1,k-1) - 64*gxy(i,j+1,k-1) + 8*gxy(i,j+2,k-1)   &
                        +   gxy(i,j-2,k-2) -  8*gxy(i,j-1,k-2) +  8*gxy(i,j+1,k-2) -   gxy(i,j+2,k-2) ) / dydz144
      d2_gd(1,3,2,3) = (   -gxz(i,j-2,k+2) +  8*gxz(i,j-1,k+2) -  8*gxz(i,j+1,k+2) +   gxz(i,j+2,k+2)   &
                        + 8*gxz(i,j-2,k+1) - 64*gxz(i,j-1,k+1) + 64*gxz(i,j+1,k+1) - 8*gxz(i,j+2,k+1)   &
                        - 8*gxz(i,j-2,k-1) + 64*gxz(i,j-1,k-1) - 64*gxz(i,j+1,k-1) + 8*gxz(i,j+2,k-1)   &
                        +   gxz(i,j-2,k-2) -  8*gxz(i,j-1,k-2) +  8*gxz(i,j+1,k-2) -   gxz(i,j+2,k-2) ) / dydz144
      d2_gd(2,2,2,3) = (   -gyy(i,j-2,k+2) +  8*gyy(i,j-1,k+2) -  8*gyy(i,j+1,k+2) +   gyy(i,j+2,k+2)   &
                        + 8*gyy(i,j-2,k+1) - 64*gyy(i,j-1,k+1) + 64*gyy(i,j+1,k+1) - 8*gyy(i,j+2,k+1)   &
                        - 8*gyy(i,j-2,k-1) + 64*gyy(i,j-1,k-1) - 64*gyy(i,j+1,k-1) + 8*gyy(i,j+2,k-1)   &
                        +   gyy(i,j-2,k-2) -  8*gyy(i,j-1,k-2) +  8*gyy(i,j+1,k-2) -   gyy(i,j+2,k-2) ) / dydz144
      d2_gd(2,3,2,3) = (   -gyz(i,j-2,k+2) +  8*gyz(i,j-1,k+2) -  8*gyz(i,j+1,k+2) +   gyz(i,j+2,k+2)   &
                        + 8*gyz(i,j-2,k+1) - 64*gyz(i,j-1,k+1) + 64*gyz(i,j+1,k+1) - 8*gyz(i,j+2,k+1)   &
                        - 8*gyz(i,j-2,k-1) + 64*gyz(i,j-1,k-1) - 64*gyz(i,j+1,k-1) + 8*gyz(i,j+2,k-1)   &
                        +   gyz(i,j-2,k-2) -  8*gyz(i,j-1,k-2) +  8*gyz(i,j+1,k-2) -   gyz(i,j+2,k-2) ) / dydz144
      d2_gd(3,3,2,3) = (   -gzz(i,j-2,k+2) +  8*gzz(i,j-1,k+2) -  8*gzz(i,j+1,k+2) +   gzz(i,j+2,k+2)   &
                        + 8*gzz(i,j-2,k+1) - 64*gzz(i,j-1,k+1) + 64*gzz(i,j+1,k+1) - 8*gzz(i,j+2,k+1)   &
                        - 8*gzz(i,j-2,k-1) + 64*gzz(i,j-1,k-1) - 64*gzz(i,j+1,k-1) + 8*gzz(i,j+2,k-1)   &
                        +   gzz(i,j-2,k-2) -  8*gzz(i,j-1,k-2) +  8*gzz(i,j+1,k-2) -   gzz(i,j+2,k-2) ) / dydz144

      d2_gd(1,1,2,1) = d2_gd(1,1,1,2)
      d2_gd(1,2,2,1) = d2_gd(1,2,1,2)
      d2_gd(1,3,2,1) = d2_gd(1,3,1,2)
      d2_gd(2,2,2,1) = d2_gd(2,2,1,2)
      d2_gd(2,3,2,1) = d2_gd(2,3,1,2)
      d2_gd(3,3,2,1) = d2_gd(3,3,1,2)

      d2_gd(1,1,3,1) = d2_gd(1,1,1,3)
      d2_gd(1,2,3,1) = d2_gd(1,2,1,3)
      d2_gd(1,3,3,1) = d2_gd(1,3,1,3)
      d2_gd(2,2,3,1) = d2_gd(2,2,1,3)
      d2_gd(2,3,3,1) = d2_gd(2,3,1,3)
      d2_gd(3,3,3,1) = d2_gd(3,3,1,3)

      d2_gd(1,1,3,2) = d2_gd(1,1,2,3)
      d2_gd(1,2,3,2) = d2_gd(1,2,2,3)
      d2_gd(1,3,3,2) = d2_gd(1,3,2,3)
      d2_gd(2,2,3,2) = d2_gd(2,2,2,3)
      d2_gd(2,3,3,2) = d2_gd(2,3,2,3)
      d2_gd(3,3,3,2) = d2_gd(3,3,2,3)

      d2_gd(2,1,:,:) = d2_gd(1,2,:,:)
      d2_gd(3,1,:,:) = d2_gd(1,3,:,:)
      d2_gd(3,2,:,:) = d2_gd(2,3,:,:)
      !------------------------------------------

    else

      ! second order derivatives as default

      !-------------- Centered 1st derivatives ----

      ! d1_gd(3,3,3)
      d1_gd(1,1,1) = (gxx(i+1,j,k) - gxx(i-1,j,k)) / dx2
      d1_gd(1,2,1) = (gxy(i+1,j,k) - gxy(i-1,j,k)) / dx2
      d1_gd(1,3,1) = (gxz(i+1,j,k) - gxz(i-1,j,k)) / dx2
      d1_gd(2,2,1) = (gyy(i+1,j,k) - gyy(i-1,j,k)) / dx2
      d1_gd(2,3,1) = (gyz(i+1,j,k) - gyz(i-1,j,k)) / dx2
      d1_gd(3,3,1) = (gzz(i+1,j,k) - gzz(i-1,j,k)) / dx2

      d1_gd(1,1,2) = (gxx(i,j+1,k) - gxx(i,j-1,k)) / dy2
      d1_gd(1,2,2) = (gxy(i,j+1,k) - gxy(i,j-1,k)) / dy2
      d1_gd(1,3,2) = (gxz(i,j+1,k) - gxz(i,j-1,k)) / dy2
      d1_gd(2,2,2) = (gyy(i,j+1,k) - gyy(i,j-1,k)) / dy2
      d1_gd(2,3,2) = (gyz(i,j+1,k) - gyz(i,j-1,k)) / dy2
      d1_gd(3,3,2) = (gzz(i,j+1,k) - gzz(i,j-1,k)) / dy2

      d1_gd(1,1,3) = (gxx(i,j,k+1) - gxx(i,j,k-1)) / dz2
      d1_gd(1,2,3) = (gxy(i,j,k+1) - gxy(i,j,k-1)) / dz2
      d1_gd(1,3,3) = (gxz(i,j,k+1) - gxz(i,j,k-1)) / dz2
      d1_gd(2,2,3) = (gyy(i,j,k+1) - gyy(i,j,k-1)) / dz2
      d1_gd(2,3,3) = (gyz(i,j,k+1) - gyz(i,j,k-1)) / dz2
      d1_gd(3,3,3) = (gzz(i,j,k+1) - gzz(i,j,k-1)) / dz2

      d1_gd(2,1,:) = d1_gd(1,2,:)
      d1_gd(3,1,:) = d1_gd(1,3,:)
      d1_gd(3,2,:) = d1_gd(2,3,:)

      ! d1_kd(3,3,3)
      d1_kd(1,1,1) = (kxx(i+1,j,k) - kxx(i-1,j,k)) / dx2
      d1_kd(1,2,1) = (kxy(i+1,j,k) - kxy(i-1,j,k)) / dx2
      d1_kd(1,3,1) = (kxz(i+1,j,k) - kxz(i-1,j,k)) / dx2
      d1_kd(2,2,1) = (kyy(i+1,j,k) - kyy(i-1,j,k)) / dx2
      d1_kd(2,3,1) = (kyz(i+1,j,k) - kyz(i-1,j,k)) / dx2
      d1_kd(3,3,1) = (kzz(i+1,j,k) - kzz(i-1,j,k)) / dx2

      d1_kd(1,1,2) = (kxx(i,j+1,k) - kxx(i,j-1,k)) / dy2
      d1_kd(1,2,2) = (kxy(i,j+1,k) - kxy(i,j-1,k)) / dy2
      d1_kd(1,3,2) = (kxz(i,j+1,k) - kxz(i,j-1,k)) / dy2
      d1_kd(2,2,2) = (kyy(i,j+1,k) - kyy(i,j-1,k)) / dy2
      d1_kd(2,3,2) = (kyz(i,j+1,k) - kyz(i,j-1,k)) / dy2
      d1_kd(3,3,2) = (kzz(i,j+1,k) - kzz(i,j-1,k)) / dy2

      d1_kd(1,1,3) = (kxx(i,j,k+1) - kxx(i,j,k-1)) / dz2
      d1_kd(1,2,3) = (kxy(i,j,k+1) - kxy(i,j,k-1)) / dz2
      d1_kd(1,3,3) = (kxz(i,j,k+1) - kxz(i,j,k-1)) / dz2
      d1_kd(2,2,3) = (kyy(i,j,k+1) - kyy(i,j,k-1)) / dz2
      d1_kd(2,3,3) = (kyz(i,j,k+1) - kyz(i,j,k-1)) / dz2
      d1_kd(3,3,3) = (kzz(i,j,k+1) - kzz(i,j,k-1)) / dz2

      d1_kd(2,1,:) = d1_kd(1,2,:)
      d1_kd(3,1,:) = d1_kd(1,3,:)
      d1_kd(3,2,:) = d1_kd(2,3,:)

      ! d1_lA removed: B^i is evolved directly in Einstein-Maxwell.

      !------------ Centered 2nd derivatives -----

      ! d2_gd(3,3,3,3)
      d2_gd(1,1,1,1) = (gxx(i+1,j,k) - 2*gxx(i,j,k) + gxx(i-1,j,k)) / dxsq
      d2_gd(1,2,1,1) = (gxy(i+1,j,k) - 2*gxy(i,j,k) + gxy(i-1,j,k)) / dxsq
      d2_gd(1,3,1,1) = (gxz(i+1,j,k) - 2*gxz(i,j,k) + gxz(i-1,j,k)) / dxsq
      d2_gd(2,2,1,1) = (gyy(i+1,j,k) - 2*gyy(i,j,k) + gyy(i-1,j,k)) / dxsq
      d2_gd(2,3,1,1) = (gyz(i+1,j,k) - 2*gyz(i,j,k) + gyz(i-1,j,k)) / dxsq
      d2_gd(3,3,1,1) = (gzz(i+1,j,k) - 2*gzz(i,j,k) + gzz(i-1,j,k)) / dxsq

      d2_gd(1,1,2,2) = (gxx(i,j+1,k) - 2*gxx(i,j,k) + gxx(i,j-1,k)) / dysq
      d2_gd(1,2,2,2) = (gxy(i,j+1,k) - 2*gxy(i,j,k) + gxy(i,j-1,k)) / dysq
      d2_gd(1,3,2,2) = (gxz(i,j+1,k) - 2*gxz(i,j,k) + gxz(i,j-1,k)) / dysq
      d2_gd(2,2,2,2) = (gyy(i,j+1,k) - 2*gyy(i,j,k) + gyy(i,j-1,k)) / dysq
      d2_gd(2,3,2,2) = (gyz(i,j+1,k) - 2*gyz(i,j,k) + gyz(i,j-1,k)) / dysq
      d2_gd(3,3,2,2) = (gzz(i,j+1,k) - 2*gzz(i,j,k) + gzz(i,j-1,k)) / dysq

      d2_gd(1,1,3,3) = (gxx(i,j,k+1) - 2*gxx(i,j,k) + gxx(i,j,k-1)) / dzsq
      d2_gd(1,2,3,3) = (gxy(i,j,k+1) - 2*gxy(i,j,k) + gxy(i,j,k-1)) / dzsq
      d2_gd(1,3,3,3) = (gxz(i,j,k+1) - 2*gxz(i,j,k) + gxz(i,j,k-1)) / dzsq
      d2_gd(2,2,3,3) = (gyy(i,j,k+1) - 2*gyy(i,j,k) + gyy(i,j,k-1)) / dzsq
      d2_gd(2,3,3,3) = (gyz(i,j,k+1) - 2*gyz(i,j,k) + gyz(i,j,k-1)) / dzsq
      d2_gd(3,3,3,3) = (gzz(i,j,k+1) - 2*gzz(i,j,k) + gzz(i,j,k-1)) / dzsq

      d2_gd(1,1,1,2) = ( gxx(i+1,j+1,k) + gxx(i-1,j-1,k)                     &
                       - gxx(i+1,j-1,k) - gxx(i-1,j+1,k) ) / dxdy4
      d2_gd(1,2,1,2) = ( gxy(i+1,j+1,k) + gxy(i-1,j-1,k)                     &
                       - gxy(i+1,j-1,k) - gxy(i-1,j+1,k) ) / dxdy4
      d2_gd(1,3,1,2) = ( gxz(i+1,j+1,k) + gxz(i-1,j-1,k)                     &
                       - gxz(i+1,j-1,k) - gxz(i-1,j+1,k) ) / dxdy4
      d2_gd(2,2,1,2) = ( gyy(i+1,j+1,k) + gyy(i-1,j-1,k)                     &
                       - gyy(i+1,j-1,k) - gyy(i-1,j+1,k) ) / dxdy4
      d2_gd(2,3,1,2) = ( gyz(i+1,j+1,k) + gyz(i-1,j-1,k)                     &
                       - gyz(i+1,j-1,k) - gyz(i-1,j+1,k) ) / dxdy4
      d2_gd(3,3,1,2) = ( gzz(i+1,j+1,k) + gzz(i-1,j-1,k)                     &
                       - gzz(i+1,j-1,k) - gzz(i-1,j+1,k) ) / dxdy4

      d2_gd(1,1,1,3) = ( gxx(i+1,j,k+1) + gxx(i-1,j,k-1)                     &
                       - gxx(i+1,j,k-1) - gxx(i-1,j,k+1) ) / dxdz4
      d2_gd(1,2,1,3) = ( gxy(i+1,j,k+1) + gxy(i-1,j,k-1)                     &
                       - gxy(i+1,j,k-1) - gxy(i-1,j,k+1) ) / dxdz4
      d2_gd(1,3,1,3) = ( gxz(i+1,j,k+1) + gxz(i-1,j,k-1)                     &
                       - gxz(i+1,j,k-1) - gxz(i-1,j,k+1) ) / dxdz4
      d2_gd(2,2,1,3) = ( gyy(i+1,j,k+1) + gyy(i-1,j,k-1)                     &
                       - gyy(i+1,j,k-1) - gyy(i-1,j,k+1) ) / dxdz4
      d2_gd(2,3,1,3) = ( gyz(i+1,j,k+1) + gyz(i-1,j,k-1)                     &
                       - gyz(i+1,j,k-1) - gyz(i-1,j,k+1) ) / dxdz4
      d2_gd(3,3,1,3) = ( gzz(i+1,j,k+1) + gzz(i-1,j,k-1)                     &
                       - gzz(i+1,j,k-1) - gzz(i-1,j,k+1) ) / dxdz4

      d2_gd(1,1,2,3) = ( gxx(i,j+1,k+1) + gxx(i,j-1,k-1)                     &
                       - gxx(i,j+1,k-1) - gxx(i,j-1,k+1) ) / dydz4
      d2_gd(1,2,2,3) = ( gxy(i,j+1,k+1) + gxy(i,j-1,k-1)                     &
                       - gxy(i,j+1,k-1) - gxy(i,j-1,k+1) ) / dydz4
      d2_gd(1,3,2,3) = ( gxz(i,j+1,k+1) + gxz(i,j-1,k-1)                     &
                       - gxz(i,j+1,k-1) - gxz(i,j-1,k+1) ) / dydz4
      d2_gd(2,2,2,3) = ( gyy(i,j+1,k+1) + gyy(i,j-1,k-1)                     &
                       - gyy(i,j+1,k-1) - gyy(i,j-1,k+1) ) / dydz4
      d2_gd(2,3,2,3) = ( gyz(i,j+1,k+1) + gyz(i,j-1,k-1)                     &
                       - gyz(i,j+1,k-1) - gyz(i,j-1,k+1) ) / dydz4
      d2_gd(3,3,2,3) = ( gzz(i,j+1,k+1) + gzz(i,j-1,k-1)                     &
                       - gzz(i,j+1,k-1) - gzz(i,j-1,k+1) ) / dydz4

      d2_gd(1,1,2,1) = d2_gd(1,1,1,2)
      d2_gd(1,2,2,1) = d2_gd(1,2,1,2)
      d2_gd(1,3,2,1) = d2_gd(1,3,1,2)
      d2_gd(2,2,2,1) = d2_gd(2,2,1,2)
      d2_gd(2,3,2,1) = d2_gd(2,3,1,2)
      d2_gd(3,3,2,1) = d2_gd(3,3,1,2)

      d2_gd(1,1,3,1) = d2_gd(1,1,1,3)
      d2_gd(1,2,3,1) = d2_gd(1,2,1,3)
      d2_gd(1,3,3,1) = d2_gd(1,3,1,3)
      d2_gd(2,2,3,1) = d2_gd(2,2,1,3)
      d2_gd(2,3,3,1) = d2_gd(2,3,1,3)
      d2_gd(3,3,3,1) = d2_gd(3,3,1,3)

      d2_gd(1,1,3,2) = d2_gd(1,1,2,3)
      d2_gd(1,2,3,2) = d2_gd(1,2,2,3)
      d2_gd(1,3,3,2) = d2_gd(1,3,2,3)
      d2_gd(2,2,3,2) = d2_gd(2,2,2,3)
      d2_gd(2,3,3,2) = d2_gd(2,3,2,3)
      d2_gd(3,3,3,2) = d2_gd(3,3,2,3)

      d2_gd(2,1,:,:) = d2_gd(1,2,:,:)
      d2_gd(3,1,:,:) = d2_gd(1,3,:,:)
      d2_gd(3,2,:,:) = d2_gd(2,3,:,:)

    end if
    !------------------------------------------

    !------------ Christoffel symbols ---------
    cf1 = 0
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          cf1(a,b,c) = 0.5d0 * (d1_gd(a,b,c) + d1_gd(a,c,b) - d1_gd(b,c,a))
        end do
      end do
    end do

    cf2 = 0
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            cf2(a,b,c) = cf2(a,b,c) + gu(a,m) * cf1(m,b,c)
          end do
        end do
      end do
    end do
    !------------------------------------------


    !------------ covariant derivs ------------
    cd1_kd = d1_kd
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            cd1_kd(a,b,c) = cd1_kd(a,b,c) - cf2(m,a,c) * kd(m,b)           &
                                          - cf2(m,b,c) * kd(a,m)
          end do
        end do
      end do
    end do

    !------------------------------------------


    !----- d1 of inverse metric and cf2 -------
    d1_gu = 0
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          do m = 1, 3
            do n = 1, 3
              d1_gu(a,b,c) = d1_gu(a,b,c) - gu(a,m) * gu(b,n) * d1_gd(m,n,c)
            end do
          end do
        end do
      end do
    end do

    d1_cf2 = 0
    do a = 1, 3
      do b = 1, 3
        do c = 1, 3
          do d = 1, 3
            do m = 1, 3
              d1_cf2(a,b,c,d) = d1_cf2(a,b,c,d) + gu(a,m) * (d2_gd(c,m,b,d)&
                                + d2_gd(m,b,c,d) - d2_gd(b,c,m,d)) / 2     &
                                + d1_gu(a,m,d) * cf1(m,b,c)
            end do
          end do
        end do
      end do
    end do
    !------------------------------------------


    !------------ Ricci tensor ----------------
    ri = 0
    do a = 1, 3
      do b = 1, 3
        do m = 1, 3
          ri(a,b) = ri(a,b) + d1_cf2(m,b,a,m) - d1_cf2(m,m,a,b)
          do n = 1, 3
            ri(a,b) = ri(a,b) + cf2(m,m,n) * cf2(n,b,a)                    &
                              - cf2(m,b,n) * cf2(n,m,a)
          end do
        end do
      end do
    end do
    !------------------------------------------


    !------------ Levi-Civita tensor ----------
    eps_lc_u        = 0
    eps_lc_u(1,2,3) = 1
    eps_lc_u(2,3,1) = 1
    eps_lc_u(3,1,2) = 1
    eps_lc_u(3,2,1) = -1
    eps_lc_u(2,1,3) = -1
    eps_lc_u(1,3,2) = -1
    eps_lc_u = eps_lc_u / sqrt(detgd)

    ! eps_lc_d = eps_lc_u * detgd
    !------------------------------------------


    ! magnetic field
    ! In the magnetic-charge Einstein-Maxwell system B^i is an evolved
    ! variable. Do not reconstruct it from a vector potential.
    Bu(1) = Bx(i,j,k)
    Bu(2) = By(i,j,k)
    Bu(3) = Bz(i,j,k)

    Bd(1)   = gd(1,1) * Bu(1) + gd(1,2) * Bu(2) + gd(1,3) * Bu(3)
    Bd(2)   = gd(2,1) * Bu(1) + gd(2,2) * Bu(2) + gd(2,3) * Bu(3)
    Bd(3)   = gd(3,1) * Bu(1) + gd(3,2) * Bu(2) + gd(3,3) * Bu(3)


    !------------ Orthonormal basis -----------
    ! Starting vectors
    xx(:) = (/ x(i,j,k), y(i,j,k), z(i,j,k) /)

    ! All points on the z-axis are pathological, since the triad vectors
    ! in the phi and theta direction are not well-defined. Take points
    ! just a little off, say at x = +epsilon.
    !
    ! we always assume z_orientation, in agreement with the Multipole thorn

    if( xx(1)**2 + xx(2)**2 < 1.0d-12 ) xx(1) = xx(1) + 1.0d-10

    ! Build an orthonormal triad.  The angular basis is ill-defined
    ! exactly on the z-axis; perturb the auxiliary construction slightly there
    ! to avoid 0/0 while leaving tensor fields unchanged.
    if (xx(1)**2 + xx(2)**2 < 1.0d-12) xx(1) = xx(1) + eps_axis

    u_vec(:) = xx(:)
    v_vec(:) = (/ xx(1)*xx(3), xx(2)*xx(3), -xx(1)**2 - xx(2)**2 /)
    w_vec(:) = (/ -xx(2), xx(1), 0.0d0 /)

    ! Normalize u_vec.
    dotp1 = 0.0d0
    do m = 1, 3
       dotp1 = dotp1 + u_vec(m) * (gd(m,1)*u_vec(1) + gd(m,2)*u_vec(2) + gd(m,3)*u_vec(3))
    end do
    if (dotp1 <= 0.0d0) cycle
    u_vec = u_vec / sqrt(dotp1)

    ! Orthogonalize and normalize v_vec.
    dotp1 = 0.0d0
    do m = 1, 3
       dotp1 = dotp1 + u_vec(m) * (gd(m,1)*v_vec(1) + gd(m,2)*v_vec(2) + gd(m,3)*v_vec(3))
    end do
    v_vec = v_vec - dotp1 * u_vec

    dotp1 = 0.0d0
    do m = 1, 3
       dotp1 = dotp1 + v_vec(m) * (gd(m,1)*v_vec(1) + gd(m,2)*v_vec(2) + gd(m,3)*v_vec(3))
    end do
    if (dotp1 <= 0.0d0) cycle
    v_vec = v_vec / sqrt(dotp1)

    ! Orthogonalize and normalize w_vec.
    dotp1 = 0.0d0
    dotp2 = 0.0d0
    do m = 1, 3
       dotp1 = dotp1 + u_vec(m) * (gd(m,1)*w_vec(1) + gd(m,2)*w_vec(2) + gd(m,3)*w_vec(3))
       dotp2 = dotp2 + v_vec(m) * (gd(m,1)*w_vec(1) + gd(m,2)*w_vec(2) + gd(m,3)*w_vec(3))
    end do
    w_vec = w_vec - dotp1 * u_vec - dotp2 * v_vec

    dotp1 = 0.0d0
    do m = 1, 3
       dotp1 = dotp1 + w_vec(m) * (gd(m,1)*w_vec(1) + gd(m,2)*w_vec(2) + gd(m,3)*w_vec(3))
    end do
    if (dotp1 <= 0.0d0) cycle
    w_vec = w_vec / sqrt(dotp1)

    ! ud_vec = matmul( gd, u_vec )
    !------------------------------------------


    !-------------------------------------------


    ! Note: what follows is not exactly the electric and magnetic part of the
    ! Weyl... these specific forms are only adequate when computing Psi4 (or
    ! Psi0), which is what we are interested in here. should we ever need any
    ! of the other Weyl scalars, we must remember to check these formulae. in
    ! any case, for our current Einstein-Maxwell or Einstein-Proca evolutions,
    ! the matter terms decay faster than 1/r, so it should make no difference
    ! whether we add them here or not.

    !------------ EB part of Weyl -------------
    elec = ri
    mag  = 0
    do a = 1, 3
      do b = 1, 3
        elec(a,b) =   elec(a,b) + Ed(a) * Ed(b) + Bd(a) * Bd(b)
        do m = 1, 3
          do n = 1, 3
            elec(a,b) = elec(a,b) + gu(m,n) * (kd(a,b) * kd(m,n)           &
                      - kd(a,m) * kd(b,n))
            do p = 1, 3
              mag(a,b) = mag(a,b) + gd(b,p) * eps_lc_u(p,m,n) * cd1_kd(n,a,m)   &
                                  + gd(a,p) * eps_lc_u(p,m,n) * cd1_kd(n,b,m)
            end do
          end do
        end do
      end do
    end do
    mag = 0.5 * mag

    ! construct tracefree part of the electric and magnetic parts of the Weyl tensor
    electr = 0.0
    magtr  = 0.0
    do m = 1, 3
      do n = 1, 3
        electr = electr + gu(m,n) * elec(m,n)
        magtr  = magtr  + gu(m,n) * mag(m,n)
      end do
    end do

    elec = elec - gd * electr / 3.0
    mag  = mag  - gd * magtr  / 3.0


    !------------------------------------------

    ! NPScalars -- GW: Psi0, Psi4
    if (calc_psi /= 0) then
       do m = 1, 3
         do n = 1, 3
           psi0re(i,j,k) = psi0re(i,j,k)                                      &
                           + 0.5 * (   elec(m,n) * (   v_vec(m) * v_vec(n)    &
                                                     - w_vec(m) * w_vec(n) )  &
                                     + mag(m,n)  * (   v_vec(m) * w_vec(n)    &
                                                     + w_vec(m) * v_vec(n) ) )
           psi0im(i,j,k) = psi0im(i,j,k)                                      &
                           + 0.5 * (   elec(m,n) * (   v_vec(m) * w_vec(n)    &
                                                     + w_vec(m) * v_vec(n) )  &
                                     + mag(m,n)  * (   w_vec(m) * w_vec(n)    &
                                                     - v_vec(m) * v_vec(n) ) )
   
           psi4re(i,j,k) = psi4re(i,j,k)                                      &
                           + 0.5 * (   elec(m,n) * (   v_vec(m) * v_vec(n)    &
                                                     - w_vec(m) * w_vec(n) )  &
                                     - mag(m,n)  * (   v_vec(m) * w_vec(n)    &
                                                     + w_vec(m) * v_vec(n) ) )
           psi4im(i,j,k) = psi4im(i,j,k)                                      &
                           + 0.5 * (  -elec(m,n) * (   v_vec(m) * w_vec(n)    &
                                                     + w_vec(m) * v_vec(n) )  &
                                     + mag(m,n)  * (   w_vec(m) * w_vec(n)    &
                                                     - v_vec(m) * v_vec(n) ) )
         end do
       end do
    end if
    !------------------------------------------

    !------------ Phi0, Phi2, Phi1 ------------
    if (calc_phi /= 0) then
       do m = 1, 3
         phi0re(i,j,k) = phi0re(i,j,k) - 0.5d0 * (   Ed(m) * v_vec(m)          &
                                                 - Bd(m) * w_vec(m) )
         phi0im(i,j,k) = phi0im(i,j,k) - 0.5d0 * (   Ed(m) * w_vec(m)          &
                                                 + Bd(m) * v_vec(m) )
       
         phi2re(i,j,k) = phi2re(i,j,k) + 0.5d0 * (   Ed(m) * v_vec(m)          &
                                                 + Bd(m) * w_vec(m) )
         phi2im(i,j,k) = phi2im(i,j,k) - 0.5d0 * (   Ed(m) * w_vec(m)          &
                                                 - Bd(m) * v_vec(m) )
       
         phi1re(i,j,k) = phi1re(i,j,k) + 0.5d0 * Ed(m) * u_vec(m)
         phi1im(i,j,k) = phi1im(i,j,k) + 0.5d0 * Bd(m) * u_vec(m)
       end do
    end if
    !------------------------------------------


    ! if( abs(y(i,j,k)) < 1.0d-05 .and. abs(z(i,j,k)) < 1.0d-05 ) then

    !    write(*,*) 'i, j, k, x = ', i, j, k, x(i,j,k)
    !    write(*,*) 'psi4re     = ', psi4re(i,j,k)
    !    write(*,*) 'psi4im     = ', psi4im(i,j,k)
    !    write(*,*) 'Ed         = ', Ed
    !    write(*,*) 'Bd         = ', Bd
    !    write(*,*) 'detgd      = ', detgd

    ! end if

  end do
  end do
  end do

end subroutine NPProca_CalcPsiGF


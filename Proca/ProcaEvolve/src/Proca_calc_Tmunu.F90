#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine Proca_calc_Tmunu( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  ! This routine implements the Einstein-Maxwell source terms used in
  ! arXiv:1205.1063, Eq. (2.7), using E^i and B^i directly.
  ! No vector potential A_i, scalar fields, or EMD coupling terms are used.

  ! Fundamental variables
  CCTK_REAL                alph, beta(3)
  CCTK_REAL                gg(3,3), gu(3,3), detgg, sqrt_detgg
  CCTK_REAL                lE(3), lB(3)
  CCTK_REAL                Ed(3), Bd(3)
  CCTK_REAL                Tab(4,4)

  ! Auxiliary variables
  CCTK_REAL                eps_lc_d(3,3,3)

  ! Matter variables
  CCTK_REAL                srcE, srcjdi(3), srcSij(3,3)
  CCTK_REAL                srcE_F, srcjdi_F(3), srcSij_F(3,3)
  CCTK_REAL                jr_F, Sir_F(3)

  ! Misc variables
  CCTK_REAL                aux_F
  CCTK_REAL                xx(3), rr

  CCTK_REAL, parameter ::  one  = 1.0d0
  CCTK_REAL, parameter ::  half = 0.5d0
  CCTK_REAL, parameter ::  pi   = acos(-one)
  CCTK_REAL, parameter ::  pi4  = 4.0d0*pi
  CCTK_REAL, parameter ::  pi8  = 8.0d0*pi
  CCTK_INT                 i, j, k
  CCTK_INT                 a, b, c

  !$OMP PARALLEL DO COLLAPSE(3) &
  !$OMP PRIVATE(k,j,i,a,b,c,alph,beta,gg,gu,detgg,sqrt_detgg, &
  !$OMP         lE,lB,Ed,Bd,Tab,eps_lc_d,srcE,srcjdi,srcSij,  &
  !$OMP         srcE_F,srcjdi_F,srcSij_F,jr_F,Sir_F,aux_F,xx,rr)
  do k = 1+cctk_nghostzones(3), cctk_lsh(3)-cctk_nghostzones(3)
     do j = 1+cctk_nghostzones(2), cctk_lsh(2)-cctk_nghostzones(2)
        do i = 1+cctk_nghostzones(1), cctk_lsh(1)-cctk_nghostzones(1)

           !------------ Get local ADM variables ----------

           alph      = alp(i,j,k)

           beta(1)   = betax(i,j,k)
           beta(2)   = betay(i,j,k)
           beta(3)   = betaz(i,j,k)

           gg(1,1)   = gxx(i,j,k)
           gg(1,2)   = gxy(i,j,k)
           gg(1,3)   = gxz(i,j,k)
           gg(2,2)   = gyy(i,j,k)
           gg(2,3)   = gyz(i,j,k)
           gg(3,3)   = gzz(i,j,k)
           gg(2,1)   = gg(1,2)
           gg(3,1)   = gg(1,3)
           gg(3,2)   = gg(2,3)

           lE(1)     = Ex(i,j,k)
           lE(2)     = Ey(i,j,k)
           lE(3)     = Ez(i,j,k)

           lB(1)     = Bx(i,j,k)
           lB(2)     = By(i,j,k)
           lB(3)     = Bz(i,j,k)

           !------------ Invert 3-metric ----------------

           detgg   =     gg(1,1) * gg(2,2) * gg(3,3)                              &
                   + 2.0d0 * gg(1,2) * gg(1,3) * gg(2,3)                          &
                   -     gg(1,1) * gg(2,3) ** 2                                   &
                   -     gg(2,2) * gg(1,3) ** 2                                   &
                   -     gg(3,3) * gg(1,2) ** 2

           sqrt_detgg = sqrt(detgg)

           gu(1,1) = (gg(2,2) * gg(3,3) - gg(2,3) ** 2     ) / detgg
           gu(2,2) = (gg(1,1) * gg(3,3) - gg(1,3) ** 2     ) / detgg
           gu(3,3) = (gg(1,1) * gg(2,2) - gg(1,2) ** 2     ) / detgg
           gu(1,2) = (gg(1,3) * gg(2,3) - gg(1,2) * gg(3,3)) / detgg
           gu(1,3) = (gg(1,2) * gg(2,3) - gg(1,3) * gg(2,2)) / detgg
           gu(2,3) = (gg(1,3) * gg(1,2) - gg(2,3) * gg(1,1)) / detgg
           gu(2,1) = gu(1,2)
           gu(3,1) = gu(1,3)
           gu(3,2) = gu(2,3)

           ! Lowered electric and magnetic fields: E_i = gamma_ij E^j,
           ! B_i = gamma_ij B^j.
           Ed(1) = gg(1,1) * lE(1) + gg(1,2) * lE(2) + gg(1,3) * lE(3)
           Ed(2) = gg(2,1) * lE(1) + gg(2,2) * lE(2) + gg(2,3) * lE(3)
           Ed(3) = gg(3,1) * lE(1) + gg(3,2) * lE(2) + gg(3,3) * lE(3)

           Bd(1) = gg(1,1) * lB(1) + gg(1,2) * lB(2) + gg(1,3) * lB(3)
           Bd(2) = gg(2,1) * lB(1) + gg(2,2) * lB(2) + gg(2,3) * lB(3)
           Bd(3) = gg(3,1) * lB(1) + gg(3,2) * lB(2) + gg(3,3) * lB(3)

           !------------ Spatial Levi-Civita tensor epsilon_ijk ----------
           ! Eq. (2.4) of arXiv:1205.1063 uses epsilon_123 = sqrt(gamma).

           eps_lc_d        =  0.0d0
           eps_lc_d(1,2,3) =  sqrt_detgg
           eps_lc_d(2,3,1) =  sqrt_detgg
           eps_lc_d(3,1,2) =  sqrt_detgg
           eps_lc_d(3,2,1) = -sqrt_detgg
           eps_lc_d(2,1,3) = -sqrt_detgg
           eps_lc_d(1,3,2) = -sqrt_detgg

           !------------ Matter terms -----------------
           !
           ! Eq. (2.7) of arXiv:1205.1063:
           !
           ! rho = (E^2 + B^2)/(8*pi)
           ! j_i = epsilon_ijk E^j B^k/(4*pi)
           ! S_ij = [-E_i E_j - B_i B_j
           !         + 1/2 gamma_ij (E^2 + B^2)]/(4*pi)
           !
           ! Here E^2 = gamma_ij E^i E^j and similarly for B.

           srcE_F = 0.0d0
           do a = 1, 3
              do b = 1, 3
                 srcE_F = srcE_F + (lE(a) * lE(b) + lB(a) * lB(b)) * gg(a,b)
              end do
           end do
           srcE = srcE_F / pi8

           srcjdi_F = 0.0d0
           do a = 1, 3
              do b = 1, 3
                 do c = 1, 3
                    srcjdi_F(a) = srcjdi_F(a) + eps_lc_d(a,b,c) * lE(b) * lB(c)
                 end do
              end do
           end do
           srcjdi = srcjdi_F / pi4

           aux_F = srcE_F
           srcSij_F = half * aux_F * gg
           do a = 1, 3
              do b = 1, 3
                 srcSij_F(a,b) = srcSij_F(a,b) - Ed(a) * Ed(b) - Bd(a) * Bd(b)
              end do
           end do
           srcSij = srcSij_F / pi4

           !------------------------------------------
           ! Optional diagnostics in the coordinate radial direction.
           ! The stored values follow the physical EM source normalization,
           ! i.e. they include the 1/(4*pi) factor.

           xx(1) = x(i,j,k)
           xx(2) = y(i,j,k)
           xx(3) = z(i,j,k)

           rr = sqrt(xx(1)**2 + xx(2)**2 + xx(3)**2)
           if (rr < eps_r) rr = eps_r

           jr_F = 0.0d0
           do a = 1, 3
              jr_F = jr_F + srcjdi(a) * xx(a) / rr
           end do

           Sir_F = 0.0d0
           do a = 1, 3
              do b = 1, 3
                 Sir_F(a) = Sir_F(a) + srcSij(a,b) * xx(b) / rr
              end do
           end do

           jrF_gf(i,j,k)    = jr_F

           SxrF_gf(i,j,k)   = Sir_F(1)
           SyrF_gf(i,j,k)   = Sir_F(2)
           SzrF_gf(i,j,k)   = Sir_F(3)

           !------------------------------------------
           ! Fill T_{mu nu}.  The indexing convention here follows the
           ! original thorn: spatial indices are 1:3 and time is 4.
           !
           ! T_ab = S_ab
           ! T_0a = -alpha j_a + beta^b S_ab
           ! T_00 = alpha^2 rho - 2 alpha beta^a j_a
           !        + beta^a beta^b S_ab

           Tab = 0.0d0
           Tab(1:3,1:3) = srcSij(1:3,1:3)

           Tab(1:3,4) = -alph * srcjdi(1:3)
           do b = 1, 3
              Tab(1:3,4) = Tab(1:3,4) + beta(b) * srcSij(1:3,b)
           end do
           Tab(4,1:3) = Tab(1:3,4)

           Tab(4,4) = alph**2 * srcE
           do a = 1, 3
              Tab(4,4) = Tab(4,4) - 2.0d0 * alph * beta(a) * srcjdi(a)
              do b = 1, 3
                 Tab(4,4) = Tab(4,4) + beta(a) * beta(b) * srcSij(a,b)
              end do
           end do

           ! Store it in the Tmunu variables.  The original routine adds to
           ! eT**, so this replacement preserves that behavior.
           eTtt(i,j,k) = eTtt(i,j,k) + Tab(4,4)
           eTtx(i,j,k) = eTtx(i,j,k) + Tab(4,1)
           eTty(i,j,k) = eTty(i,j,k) + Tab(4,2)
           eTtz(i,j,k) = eTtz(i,j,k) + Tab(4,3)
           eTxx(i,j,k) = eTxx(i,j,k) + Tab(1,1)
           eTxy(i,j,k) = eTxy(i,j,k) + Tab(1,2)
           eTxz(i,j,k) = eTxz(i,j,k) + Tab(1,3)
           eTyy(i,j,k) = eTyy(i,j,k) + Tab(2,2)
           eTyz(i,j,k) = eTyz(i,j,k) + Tab(2,3)
           eTzz(i,j,k) = eTzz(i,j,k) + Tab(3,3)

        end do
     end do
  end do
  !$OMP END PARALLEL DO

end subroutine Proca_calc_Tmunu

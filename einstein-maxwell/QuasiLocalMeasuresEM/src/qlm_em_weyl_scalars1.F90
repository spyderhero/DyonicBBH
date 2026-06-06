#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_calc_weyl_scalars1 (CCTK_ARGUMENTS, hn)
  use adm_metric_simple
  use cctk
  use constants
  use qlm_em_derivs
  use qlm_em_variables
  use ricci
  use ricci4
  use tensor
  use tensor4
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL, parameter :: two=2, four=4
  CCTK_REAL    :: gg(3,3), dgg(3,3,3), ddgg(3,3,3,3), gg_dot(3,3), gg_dot2(3,3), dgg_dot(3,3,3)
  CCTK_REAL    :: kk(3,3), dkk(3,3,3), kk_dot(3,3)
  CCTK_REAL    :: tt(3,3)
  CCTK_REAL    :: dtg, gu(3,3), dgu(3,3,3), gamma(3,3,3), dgamma(3,3,3,3), ri(3,3), rsc
  CCTK_REAL    :: g4(0:3,0:3), dg4(0:3,0:3,0:3), ddg4(0:3,0:3,0:3,0:3)
  CCTK_REAL    :: gu4(0:3,0:3), dgu4(0:3,0:3,0:3)
  CCTK_REAL    :: gamma4(0:3,0:3,0:3), dgamma4(0:3,0:3,0:3,0:3)
  CCTK_REAL    :: ri4(0:3,0:3), rsc4
  CCTK_REAL    :: rm4(0:3,0:3,0:3,0:3), we4(0:3,0:3,0:3,0:3)
  CCTK_REAL    :: ll(0:3), nn(0:3)
  CCTK_COMPLEX :: mm(0:3)
  CCTK_REAL    :: ss(0:3)
  CCTK_REAL    :: nabla_ll(0:3,0:3), nabla_nn(0:3,0:3)
  CCTK_REAL    :: nabla_ss(0:3,0:3)
  CCTK_REAL    :: trkAH, kk_kk, tmpR
  
  integer      :: i, j
  integer      :: a, b, c, d
  CCTK_REAL    :: theta, phi
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Calculating Weyl scalars")
  end if
  
  ! Calculate the coordinates
  do j = 1+qlm_em_nghostsphi(hn), qlm_em_nphi(hn)-qlm_em_nghostsphi(hn)
     do i = 1+qlm_em_nghoststheta(hn), qlm_em_ntheta(hn)-qlm_em_nghoststheta(hn)
        theta = qlm_em_origin_theta(hn) + (i-1)*qlm_em_delta_theta(hn)
        phi   = qlm_em_origin_phi(hn)   + (j-1)*qlm_em_delta_phi(hn)
        
        ! Get the stuff from the arrays
        gg(1,1) = qlm_em_gxx(i,j)
        gg(1,2) = qlm_em_gxy(i,j)
        gg(1,3) = qlm_em_gxz(i,j)
        gg(2,2) = qlm_em_gyy(i,j)
        gg(2,3) = qlm_em_gyz(i,j)
        gg(3,3) = qlm_em_gzz(i,j)
        gg(2,1) = gg(1,2)
        gg(3,1) = gg(1,3)
        gg(3,2) = gg(2,3)
        
        dgg(1,1,1) = qlm_em_dgxxx(i,j)
        dgg(1,2,1) = qlm_em_dgxyx(i,j)
        dgg(1,3,1) = qlm_em_dgxzx(i,j)
        dgg(2,2,1) = qlm_em_dgyyx(i,j)
        dgg(2,3,1) = qlm_em_dgyzx(i,j)
        dgg(3,3,1) = qlm_em_dgzzx(i,j)
        dgg(1,1,2) = qlm_em_dgxxy(i,j)
        dgg(1,2,2) = qlm_em_dgxyy(i,j)
        dgg(1,3,2) = qlm_em_dgxzy(i,j)
        dgg(2,2,2) = qlm_em_dgyyy(i,j)
        dgg(2,3,2) = qlm_em_dgyzy(i,j)
        dgg(3,3,2) = qlm_em_dgzzy(i,j)
        dgg(1,1,3) = qlm_em_dgxxz(i,j)
        dgg(1,2,3) = qlm_em_dgxyz(i,j)
        dgg(1,3,3) = qlm_em_dgxzz(i,j)
        dgg(2,2,3) = qlm_em_dgyyz(i,j)
        dgg(2,3,3) = qlm_em_dgyzz(i,j)
        dgg(3,3,3) = qlm_em_dgzzz(i,j)
        dgg(2,1,:) = dgg(1,2,:)
        dgg(3,1,:) = dgg(1,3,:)
        dgg(3,2,:) = dgg(2,3,:)
        
        ddgg(1,1,1,1) = qlm_em_ddgxxxx(i,j)
        ddgg(1,2,1,1) = qlm_em_ddgxyxx(i,j)
        ddgg(1,3,1,1) = qlm_em_ddgxzxx(i,j)
        ddgg(2,2,1,1) = qlm_em_ddgyyxx(i,j)
        ddgg(2,3,1,1) = qlm_em_ddgyzxx(i,j)
        ddgg(3,3,1,1) = qlm_em_ddgzzxx(i,j)
        ddgg(1,1,1,2) = qlm_em_ddgxxxy(i,j)
        ddgg(1,2,1,2) = qlm_em_ddgxyxy(i,j)
        ddgg(1,3,1,2) = qlm_em_ddgxzxy(i,j)
        ddgg(2,2,1,2) = qlm_em_ddgyyxy(i,j)
        ddgg(2,3,1,2) = qlm_em_ddgyzxy(i,j)
        ddgg(3,3,1,2) = qlm_em_ddgzzxy(i,j)
        ddgg(1,1,1,3) = qlm_em_ddgxxxz(i,j)
        ddgg(1,2,1,3) = qlm_em_ddgxyxz(i,j)
        ddgg(1,3,1,3) = qlm_em_ddgxzxz(i,j)
        ddgg(2,2,1,3) = qlm_em_ddgyyxz(i,j)
        ddgg(2,3,1,3) = qlm_em_ddgyzxz(i,j)
        ddgg(3,3,1,3) = qlm_em_ddgzzxz(i,j)
        ddgg(1,1,2,2) = qlm_em_ddgxxyy(i,j)
        ddgg(1,2,2,2) = qlm_em_ddgxyyy(i,j)
        ddgg(1,3,2,2) = qlm_em_ddgxzyy(i,j)
        ddgg(2,2,2,2) = qlm_em_ddgyyyy(i,j)
        ddgg(2,3,2,2) = qlm_em_ddgyzyy(i,j)
        ddgg(3,3,2,2) = qlm_em_ddgzzyy(i,j)
        ddgg(1,1,2,3) = qlm_em_ddgxxyz(i,j)
        ddgg(1,2,2,3) = qlm_em_ddgxyyz(i,j)
        ddgg(1,3,2,3) = qlm_em_ddgxzyz(i,j)
        ddgg(2,2,2,3) = qlm_em_ddgyyyz(i,j)
        ddgg(2,3,2,3) = qlm_em_ddgyzyz(i,j)
        ddgg(3,3,2,3) = qlm_em_ddgzzyz(i,j)
        ddgg(1,1,3,3) = qlm_em_ddgxxzz(i,j)
        ddgg(1,2,3,3) = qlm_em_ddgxyzz(i,j)
        ddgg(1,3,3,3) = qlm_em_ddgxzzz(i,j)
        ddgg(2,2,3,3) = qlm_em_ddgyyzz(i,j)
        ddgg(2,3,3,3) = qlm_em_ddgyzzz(i,j)
        ddgg(3,3,3,3) = qlm_em_ddgzzzz(i,j)
        ddgg(2,1,:,:) = ddgg(1,2,:,:)
        ddgg(3,1,:,:) = ddgg(1,3,:,:)
        ddgg(3,2,:,:) = ddgg(2,3,:,:)
        ddgg(:,:,2,1) = ddgg(:,:,1,2)
        ddgg(:,:,3,1) = ddgg(:,:,1,3)
        ddgg(:,:,3,2) = ddgg(:,:,2,3)
        
        kk(1,1) = qlm_em_kxx(i,j)
        kk(1,2) = qlm_em_kxy(i,j)
        kk(1,3) = qlm_em_kxz(i,j)
        kk(2,2) = qlm_em_kyy(i,j)
        kk(2,3) = qlm_em_kyz(i,j)
        kk(3,3) = qlm_em_kzz(i,j)
        kk(2,1) = kk(1,2)
        kk(3,1) = kk(1,3)
        kk(3,2) = kk(2,3)
        
        dkk(1,1,1) = qlm_em_dkxxx(i,j)
        dkk(1,2,1) = qlm_em_dkxyx(i,j)
        dkk(1,3,1) = qlm_em_dkxzx(i,j)
        dkk(2,2,1) = qlm_em_dkyyx(i,j)
        dkk(2,3,1) = qlm_em_dkyzx(i,j)
        dkk(3,3,1) = qlm_em_dkzzx(i,j)
        dkk(1,1,2) = qlm_em_dkxxy(i,j)
        dkk(1,2,2) = qlm_em_dkxyy(i,j)
        dkk(1,3,2) = qlm_em_dkxzy(i,j)
        dkk(2,2,2) = qlm_em_dkyyy(i,j)
        dkk(2,3,2) = qlm_em_dkyzy(i,j)
        dkk(3,3,2) = qlm_em_dkzzy(i,j)
        dkk(1,1,3) = qlm_em_dkxxz(i,j)
        dkk(1,2,3) = qlm_em_dkxyz(i,j)
        dkk(1,3,3) = qlm_em_dkxzz(i,j)
        dkk(2,2,3) = qlm_em_dkyyz(i,j)
        dkk(2,3,3) = qlm_em_dkyzz(i,j)
        dkk(3,3,3) = qlm_em_dkzzz(i,j)
        dkk(2,1,:) = dkk(1,2,:)
        dkk(3,1,:) = dkk(1,3,:)
        dkk(3,2,:) = dkk(2,3,:)
        
        tt(1,1) = qlm_em_txx(i,j)
        tt(1,2) = qlm_em_txy(i,j)
        tt(1,3) = qlm_em_txz(i,j)
        tt(2,2) = qlm_em_tyy(i,j)
        tt(2,3) = qlm_em_tyz(i,j)
        tt(3,3) = qlm_em_tzz(i,j)
        tt(2,1) = tt(1,2)
        tt(3,1) = tt(1,3)
        tt(3,2) = tt(2,3)
        
        ll(0) = qlm_em_l0(i,j,hn)
        ll(1) = qlm_em_l1(i,j,hn)
        ll(2) = qlm_em_l2(i,j,hn)
        ll(3) = qlm_em_l3(i,j,hn)
        
        nn(0) = qlm_em_n0(i,j,hn)
        nn(1) = qlm_em_n1(i,j,hn)
        nn(2) = qlm_em_n2(i,j,hn)
        nn(3) = qlm_em_n3(i,j,hn)
        
        mm(0) = qlm_em_m0(i,j,hn)
        mm(1) = qlm_em_m1(i,j,hn)
        mm(2) = qlm_em_m2(i,j,hn)
        mm(3) = qlm_em_m3(i,j,hn)
        
        ss = (ll - nn) / sqrt(two)
        
        nabla_ll = qlm_em_tetrad_derivs(i,j)%nabla_ll
        nabla_nn = qlm_em_tetrad_derivs(i,j)%nabla_nn
        
        nabla_ss = (nabla_ll - nabla_nn) / sqrt(two)
        
        
        
        ! Calculate 4-metric
        call calc_det (gg, dtg)
        call calc_inv (gg, dtg, gu)
        call calc_invderiv (gu, dgg, dgu)
        call calc_connections (gu, dgg, gamma)
        call calc_connectionderivs (gu, dgg, dgu, ddgg, dgamma)
        call calc_ricci (gamma, dgamma, ri)
        call calc_trace (gu, ri, rsc)
        
        call calc_3metricdot_simple (kk, gg_dot)
        call calc_3metricderivdot_simple (dkk, dgg_dot)
        call calc_extcurvdot_simple (gg,gu,ri, kk, tt, kk_dot)
        call calc_3metricdot2_simple (kk_dot, gg_dot2)
        
        call calc_4metricderivs2_simple (gg, dgg, &
             ddgg, gg_dot, gg_dot2, dgg_dot, g4,dg4,ddg4)
        call calc_4inv (g4, gu4)
        call calc_4invderiv (gu4, dg4, dgu4)
        call calc_4connections (gu4,dg4, gamma4)
        call calc_4connectionderivs (gu4, dg4, dgu4, ddg4, dgamma4)
        call calc_4ricci (gamma4, dgamma4, ri4)
        call calc_4riemann (g4, gamma4, dgamma4, rm4)
        call calc_4trace (ri4, gu4, rsc4)
        call calc_4weyl (g4, rm4, ri4, rsc4, we4)
        
        ! debugging       
!        qlm_em_rsc4(i,j,hn) = rsc4 
        
        qlm_em_psi0(i,j,hn) = 0     ! transverse radiation along n
        qlm_em_psi1(i,j,hn) = 0     ! longitudinal radiation along n
        qlm_em_psi2(i,j,hn) = 0     ! Coulomb field and spin
        qlm_em_psi3(i,j,hn) = 0     ! longitudinal radiation along l
        qlm_em_psi4(i,j,hn) = 0     ! transverse radiation along l
        
        do a=0,3
           do b=0,3
              do c=0,3
                 do d=0,3
                    qlm_em_psi0(i,j,hn) = qlm_em_psi0(i,j,hn) + we4(a,b,c,d) * ll(a) * mm(b) * ll(c) * mm(d)
                    qlm_em_psi1(i,j,hn) = qlm_em_psi1(i,j,hn) + we4(a,b,c,d) * ll(a) * mm(b) * ll(c) * nn(d)
                    qlm_em_psi2(i,j,hn) = qlm_em_psi2(i,j,hn) + we4(a,b,c,d) * ll(a) * mm(b) * conjg(mm(c)) * nn(d)
                    qlm_em_psi3(i,j,hn) = qlm_em_psi3(i,j,hn) + we4(a,b,c,d) * ll(a) * nn(b) * conjg(mm(c)) * nn(d)
                    qlm_em_psi4(i,j,hn) = qlm_em_psi4(i,j,hn) + we4(a,b,c,d) * conjg(mm(a)) * nn(b) * conjg(mm(c)) * nn(d)
                 end do
              end do
           end do
        end do
        
        ! gr-qc/0104063, (3.3)
        qlm_em_i(i,j,hn) = + 3 * qlm_em_psi2(i,j,hn)**2 &
             &         - 4 * qlm_em_psi1(i,j,hn) * qlm_em_psi3(i,j,hn) &
             &         +     qlm_em_psi0(i,j,hn) * qlm_em_psi4(i,j,hn)
        qlm_em_j(i,j,hn) = -     qlm_em_psi2(i,j,hn)**3 &
             &         +     qlm_em_psi0(i,j,hn) * qlm_em_psi2(i,j,hn) * qlm_em_psi4(i,j,hn) &
             &         + 2 * qlm_em_psi1(i,j,hn) * qlm_em_psi2(i,j,hn) * qlm_em_psi3(i,j,hn) &
             &         -     qlm_em_psi1(i,j,hn)**2 * qlm_em_psi4(i,j,hn) &
             &         -     qlm_em_psi0(i,j,hn) * qlm_em_psi3(i,j,hn)**2
        
        ! gr-qc/0104063, (3.1)
        qlm_em_s(i,j,hn) = 27 * qlm_em_j(i,j,hn)**2 / qlm_em_i(i,j,hn)**3
        qlm_em_sdiff(i,j,hn) = (27 * qlm_em_j(i,j,hn)**2 - qlm_em_i(i,j,hn)**3) / sqrt((abs2(qlm_em_psi0(i,j,hn)) + abs2(qlm_em_psi1(i,j,hn)) + abs2(qlm_em_psi2(i,j,hn)) + abs2(qlm_em_psi3(i,j,hn)) + abs2(qlm_em_psi4(i,j,hn))) / 5)**3
        
        qlm_em_phi00(i,j,hn) = 0
        qlm_em_phi11(i,j,hn) = 0
        qlm_em_phi01(i,j,hn) = 0
        qlm_em_phi12(i,j,hn) = 0
        qlm_em_phi10(i,j,hn) = 0
        qlm_em_phi21(i,j,hn) = 0
        qlm_em_phi02(i,j,hn) = 0
        qlm_em_phi22(i,j,hn) = 0
        qlm_em_phi20(i,j,hn) = 0
        
        do a=0,3
           do b=0,3
              
              qlm_em_phi00(i,j,hn) = qlm_em_phi00(i,j,hn) - 1/two * ri4(a,b) * ll(a) * ll(b)
              qlm_em_phi11(i,j,hn) = qlm_em_phi11(i,j,hn) - 1/two * ri4(a,b) * (ll(a) * nn(b) + mm(a) * conjg(mm(b))) / 2
              qlm_em_phi01(i,j,hn) = qlm_em_phi01(i,j,hn) - 1/two * ri4(a,b) * ll(a) * mm(b)
              qlm_em_phi12(i,j,hn) = qlm_em_phi12(i,j,hn) - 1/two * ri4(a,b) * nn(a) * mm(b)
              qlm_em_phi10(i,j,hn) = qlm_em_phi10(i,j,hn) - 1/two * ri4(a,b) * ll(a) * conjg(mm(b))
              qlm_em_phi21(i,j,hn) = qlm_em_phi21(i,j,hn) - 1/two * ri4(a,b) * nn(a) * conjg(mm(b))
              qlm_em_phi02(i,j,hn) = qlm_em_phi02(i,j,hn) - 1/two * ri4(a,b) * mm(a) * mm(b)
              qlm_em_phi22(i,j,hn) = qlm_em_phi22(i,j,hn) - 1/two * ri4(a,b) * nn(a) * nn(b)
              qlm_em_phi20(i,j,hn) = qlm_em_phi20(i,j,hn) - 1/two * ri4(a,b) * conjg(mm(a)) * conjg(mm(b))
              
           end do
        end do
        
        qlm_em_lambda(i,j,hn) = rsc4 / 24
        
        qlm_em_lie_n_theta_l(i,j,hn) = &
             & + 2 * real (qlm_em_npsigma(i,j,hn) * qlm_em_nplambda(i,j,hn)) &
             & + 2 * real (qlm_em_psi2(i,j,hn)) &
             & + 4 * qlm_em_lambda(i,j,hn)
        
        

        trkAH = 0
        do a=1,3
           do b=1,3
              trkAH = trkAH + (gu(a,b) - ss(a) * ss(b)) * nabla_ss(a,b) 
           end do
        end do
        
        kk_kk = 0
        do a=1,3
           do b=1,3
              do c=1,3
                 do d=1,3
                    kk_kk = kk_kk + &
                         (gu(a,c) - ss(a) * ss(c)) * nabla_ss(a,b) * &
                         (gu(b,d) - ss(b) * ss(d)) * nabla_ss(c,d)
                 end do
              end do
           end do
        end do
        
        tmpR = 0
        do a=1,3
           do b=1,3
              tmpR = tmpR + ri(a,b) * ss(a) * ss(b)
           end do
        end do
        
        qlm_em_rsc(i,j,hn) = rsc - 2*tmpR + trkAH**2 - kk_kk
        
     end do
  end do
end subroutine qlm_em_calc_weyl_scalars1

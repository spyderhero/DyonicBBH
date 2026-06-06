#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_calc_twometric (CCTK_ARGUMENTS, hn)
  use adm_metric
  use cctk
  use qlm_em_boundary
  use qlm_em_derivs
  use qlm_em_variables
  use ricci
  use ricci2
  use tensor
  use tensor2
  use tensor4
  
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL :: gg(3,3), dgg(3,3,3)
  CCTK_REAL :: ee(3,2), dee(3,2,2)
  CCTK_REAL :: qq(2,2), dqq(2,2,2)
 
  CCTK_REAL :: delta_space(2)
  
  integer   :: i, j
  integer   :: a, b, c, d, e, f
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Calculating two-metric")
  end if
  
  delta_space(:) = (/ qlm_em_delta_theta(hn), qlm_em_delta_phi(hn) /)
  
  ! Calculate the two-metric
  do j = 1+qlm_em_nghostsphi(hn), qlm_em_nphi(hn)-qlm_em_nghostsphi(hn)
     do i = 1+qlm_em_nghoststheta(hn), qlm_em_ntheta(hn)-qlm_em_nghoststheta(hn)
        
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
        
        ee(1,1:2) = deriv (qlm_em_x(:,:,hn), i, j, delta_space)
        ee(2,1:2) = deriv (qlm_em_y(:,:,hn), i, j, delta_space)
        ee(3,1:2) = deriv (qlm_em_z(:,:,hn), i, j, delta_space)
        
        dee(1,1:2,1:2) = deriv2 (qlm_em_x(:,:,hn), i, j, delta_space)
        dee(2,1:2,1:2) = deriv2 (qlm_em_y(:,:,hn), i, j, delta_space)
        dee(3,1:2,1:2) = deriv2 (qlm_em_z(:,:,hn), i, j, delta_space)
        
        do a=1,2
           do b=1,2
              qq(a,b) = 0
              do c=1,3
                 do d=1,3
                    qq(a,b) = qq(a,b) + gg(c,d) * ee(c,a) * ee(d,b)
                 end do
              end do
           end do
        end do
        
        do a=1,2
           do b=1,2
              do c=1,2
                 dqq(a,b,c) = 0
                 do d=1,3
                    do e=1,3
                       do f=1,3
                          dqq(a,b,c) = dqq(a,b,c) + dgg(d,e,f) * ee(d,a) * ee(e,b) * ee(f,c)
                       end do
                       dqq(a,b,c) = dqq(a,b,c) + gg(d,e) * dee(d,a,c) * ee(e,b)
                       dqq(a,b,c) = dqq(a,b,c) + gg(d,e) * ee(d,a) * dee(e,b,c)
                    end do
                 end do
              end do
           end do
        end do
        
        ! Could also calculate this as:
        !    q^ab = m^a mbar^b + mbar^a m^b
        qlm_em_qtt(i,j,hn) = qq(1,1)
        qlm_em_qtp(i,j,hn) = qq(1,2)
        qlm_em_qpp(i,j,hn) = qq(2,2)
        
        qlm_em_dqttt(i,j) = dqq(1,1,1)
        qlm_em_dqtpt(i,j) = dqq(1,2,1)
        qlm_em_dqppt(i,j) = dqq(2,2,1)
        qlm_em_dqttp(i,j) = dqq(1,1,2)
        qlm_em_dqtpp(i,j) = dqq(1,2,2)
        qlm_em_dqppp(i,j) = dqq(2,2,2)
        
     end do
  end do
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_qtt(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_qtp(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_qpp(:,:,hn), +1)
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_dqttt(:,:), -1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_dqtpt(:,:), -1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_dqppt(:,:), -1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_dqttp(:,:), -1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_dqtpp(:,:), -1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_dqppp(:,:), -1)
  
end subroutine qlm_em_calc_twometric

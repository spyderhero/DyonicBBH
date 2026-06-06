#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_killing_test (CCTK_ARGUMENTS, hn)
  use cctk
  use qlm_em_boundary
  use qlm_em_derivs
  use qlm_em_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL :: qq(2,2), dqq(2,2,2)
  CCTK_REAL :: xi(2), dxi(2,2)
  CCTK_REAL :: lqq(2,2)
  
  CCTK_REAL    :: delta_space(2)
  
  integer   :: i,j
  integer   :: a, b, c
  CCTK_REAL :: theta, phi
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Testing Killing vector field")
  end if
  
  delta_space(:) = (/ qlm_em_delta_theta(hn), qlm_em_delta_phi(hn) /)
  
  ! Calculate the coordinates
  do j = 1+qlm_em_nghostsphi(hn), qlm_em_nphi(hn)-qlm_em_nghostsphi(hn)
     do i = 1+qlm_em_nghoststheta(hn), qlm_em_ntheta(hn)-qlm_em_nghoststheta(hn)
        theta = qlm_em_origin_theta(hn) + (i-1)*qlm_em_delta_theta(hn)
        phi   = qlm_em_origin_phi(hn)   + (j-1)*qlm_em_delta_phi(hn)
        
        ! 2-metric on the horizon
        qq(1,1) = qlm_em_qtt(i,j,hn)
        qq(1,2) = qlm_em_qtp(i,j,hn)
        qq(2,2) = qlm_em_qpp(i,j,hn)
        qq(2,1) = qq(1,2)
        
        dqq(1,1,1) = qlm_em_dqttt(i,j)
        dqq(1,2,1) = qlm_em_dqtpt(i,j)
        dqq(2,2,1) = qlm_em_dqppt(i,j)
        dqq(1,1,2) = qlm_em_dqttp(i,j)
        dqq(1,2,2) = qlm_em_dqtpp(i,j)
        dqq(2,2,2) = qlm_em_dqppp(i,j)
        dqq(2,1,:) = dqq(1,2,:)
        
        xi(1) = qlm_em_xi_t(i,j,hn)
        xi(2) = qlm_em_xi_p(i,j,hn)
        
        dxi(1,1:2) = deriv (qlm_em_xi_t(:,:,hn), i, j, delta_space)
        dxi(2,1:2) = deriv (qlm_em_xi_p(:,:,hn), i, j, delta_space)
        
        ! L_xi q_ab = q_cb d_a xi^c + q_ac d_b xi^c + xi^c d_c q_ab
        do a=1,2
           do b=1,2
              lqq(a,b) = 0
              do c=1,2
                 lqq(a,b) = lqq(a,b) + qq(c,b) * dxi(c,a) &
                      &              + qq(a,c) * dxi(c,b) &
                      &              + xi(c) * dqq(a,b,c)
              end do
           end do
        end do
        
        qlm_em_lqtt(i,j,hn) = lqq(1,1)
        qlm_em_lqtp(i,j,hn) = lqq(1,2)
        qlm_em_lqpp(i,j,hn) = lqq(2,2)
        
     end do
  end do
  
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_lqtt(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_lqtp(:,:,hn), +1)
  call em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_lqpp(:,:,hn), +1)
  
end subroutine qlm_em_killing_test

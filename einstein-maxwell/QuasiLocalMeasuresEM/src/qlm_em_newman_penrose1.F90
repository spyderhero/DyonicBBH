#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_calc_newman_penrose1 (CCTK_ARGUMENTS, hn)
  use adm_metric
  use cctk
  use classify
  use constants
  use qlm_em_derivs
  use qlm_em_variables
  use ricci4
  use tensor
  use tensor4
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL, parameter :: zero=0, one=1, two=2, half=1/two
  CCTK_REAL    :: ll(0:3), nn(0:3)
  CCTK_COMPLEX :: mm(0:3)
  CCTK_REAL    :: nabla_ll(0:3,0:3), nabla_nn(0:3,0:3)
  CCTK_COMPLEX :: nabla_mm(0:3,0:3)
  
  CCTK_REAL    :: t0, t1, t2
  logical      :: ce0, ce1, ce2
  CCTK_REAL    :: delta_space(2)
  
  integer      :: i, j
  integer      :: a, b
  CCTK_REAL    :: theta, phi
  
  if (veryverbose/=0) then
     call CCTK_INFO ("Calculating Newman-Penrose quantities")
  end if
  
  t0 = qlm_em_time(hn)
  t1 = qlm_em_time_p(hn)
  t2 = qlm_em_time_p_p(hn)
  
  ce0 = qlm_em_have_valid_data(hn) == 0
  ce1 = qlm_em_have_valid_data_p(hn) == 0
  ce2 = qlm_em_have_valid_data_p_p(hn) == 0
  
  delta_space(:) = (/ qlm_em_delta_theta(hn), qlm_em_delta_phi(hn) /)
  
  if (qlm_em_nghoststheta(hn)<2 .or. qlm_em_nghostsphi(hn)<2) call CCTK_WARN (0, "internal error")

  ! Calculate the coordinates
  do j = 1+qlm_em_nghostsphi(hn), qlm_em_nphi(hn)-qlm_em_nghostsphi(hn)
     do i = 1+qlm_em_nghoststheta(hn), qlm_em_ntheta(hn)-qlm_em_nghoststheta(hn)
        theta = qlm_em_origin_theta(hn) + (i-1)*qlm_em_delta_theta(hn)
        phi   = qlm_em_origin_phi(hn)   + (j-1)*qlm_em_delta_phi(hn)
        
        ! Get the stuff from the arrays
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
        
        nabla_ll(:,:) = qlm_em_tetrad_derivs(i,j)%nabla_ll(:,:)
        nabla_nn(:,:) = qlm_em_tetrad_derivs(i,j)%nabla_nn(:,:)
        nabla_mm(:,:) = qlm_em_tetrad_derivs(i,j)%nabla_mm(:,:)
        
        
        
        ! kappa   = - m^a  l^b D_b l_a
        ! tau     = - m^a  n^b D_b l_a
        ! sigma   = - m^a  m^b D_b l_a
        ! rho     = - m^a ~m^b D_b l_a   [ ~m = congj(m) ]
        qlm_em_npkappa  (i,j,hn) = 0
        qlm_em_nptau    (i,j,hn) = 0
        qlm_em_npsigma  (i,j,hn) = 0
        qlm_em_nprho    (i,j,hn) = 0
        ! epsilon = - 1/2 ( n^a  l^b D_b l_a - ~m^a  l^b D_b m_a )
        ! gamma   = - 1/2 ( n^a  n^b D_b l_a - ~m^a  n^b D_b m_a )
        ! beta    = - 1/2 ( n^a  m^b D_b l_a - ~m^a  m^b D_b m_a )
        ! alpha   = - 1/2 ( n^a ~m^b D_b l_a - ~m^a ~m^b D_b m_a )
        qlm_em_npepsilon(i,j,hn) = 0
        qlm_em_npgamma  (i,j,hn) = 0
        qlm_em_npbeta   (i,j,hn) = 0
        qlm_em_npalpha  (i,j,hn) = 0
        ! pi      = ~m^a  l^b D_b n_a
        ! nu      = ~m^a  n^b D_b n_a
        ! mu      = ~m^a  m^b D_b n_a
        ! lambda  = ~m^a ~m^b D_b n_a
        qlm_em_nppi     (i,j,hn) = 0
        qlm_em_npnu     (i,j,hn) = 0
        qlm_em_npmu     (i,j,hn) = 0
        qlm_em_nplambda (i,j,hn) = 0
        
!!$        qlm_em_lie_l_npsigma(i,j,hn) = 0
!!$        qlm_em_lie_n_npsigma(i,j,hn) = 0
        
!!$        ! Theta is the expansion
!!$        ! Theta_(l) is (- 2 Re rho)
!!$        ! Theta_(n) is (+ 2 Re mu)
!!$        ! Theta_(T) is [Theta_(l) + Theta_(n)] / sqrt(2)
!!$        ! Theta_(S) is [Theta_(l) - Theta_(n)] / sqrt(2)
        
!!$        qlm_em_lie_l_theta_l(i,j,hn) = 0
!!$        qlm_em_lie_l_theta_n(i,j,hn) = 0
!!$        qlm_em_lie_n_theta_l(i,j,hn) = 0
!!$        qlm_em_lie_n_theta_n(i,j,hn) = 0
        
        do a=0,3
           do b=0,3
              qlm_em_npkappa  (i,j,hn) = qlm_em_npkappa  (i,j,hn) - mm(a) *       ll(b)  * nabla_ll(a,b)
              qlm_em_nptau    (i,j,hn) = qlm_em_nptau    (i,j,hn) - mm(a) *       nn(b)  * nabla_ll(a,b)
              qlm_em_npsigma  (i,j,hn) = qlm_em_npsigma  (i,j,hn) - mm(a) *       mm(b)  * nabla_ll(a,b)
              qlm_em_nprho    (i,j,hn) = qlm_em_nprho    (i,j,hn) - mm(a) * conjg(mm(b)) * nabla_ll(a,b)
              qlm_em_npepsilon(i,j,hn) = qlm_em_npepsilon(i,j,hn) - half * (nn(a) *       ll(b)  * nabla_ll(a,b) - conjg(mm(a)) *       ll(b)  * nabla_mm(a,b))
              qlm_em_npgamma  (i,j,hn) = qlm_em_npgamma  (i,j,hn) - half * (nn(a) *       nn(b)  * nabla_ll(a,b) - conjg(mm(a)) *       nn(b)  * nabla_mm(a,b))
              qlm_em_npbeta   (i,j,hn) = qlm_em_npbeta   (i,j,hn) - half * (nn(a) *       mm(b)  * nabla_ll(a,b) - conjg(mm(a)) *       mm(b)  * nabla_mm(a,b))
              qlm_em_npalpha  (i,j,hn) = qlm_em_npalpha  (i,j,hn) - half * (nn(a) * conjg(mm(b)) * nabla_ll(a,b) - conjg(mm(a)) * conjg(mm(b)) * nabla_mm(a,b))
              qlm_em_nppi     (i,j,hn) = qlm_em_nppi     (i,j,hn) + conjg(mm(a)) *       ll(b)  * nabla_nn(a,b)
              qlm_em_npnu     (i,j,hn) = qlm_em_npnu     (i,j,hn) + conjg(mm(a)) *       nn(b)  * nabla_nn(a,b)
              qlm_em_npmu     (i,j,hn) = qlm_em_npmu     (i,j,hn) + conjg(mm(a)) *       mm(b)  * nabla_nn(a,b)
              qlm_em_nplambda (i,j,hn) = qlm_em_nplambda (i,j,hn) + conjg(mm(a)) * conjg(mm(b)) * nabla_nn(a,b)
              
           end do
        end do
        
     end do
  end do
end subroutine qlm_em_calc_newman_penrose1

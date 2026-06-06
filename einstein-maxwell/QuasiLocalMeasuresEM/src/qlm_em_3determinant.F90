#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



SUBROUTINE qlm_em_calc_3determinant (CCTK_ARGUMENTS, hn)
  USE adm_metric
  USE cctk
  USE qlm_em_boundary
  USE qlm_em_derivs
  USE qlm_em_variables
  USE tensor
  
  IMPLICIT NONE
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  INTEGER :: hn
  
  CCTK_REAL    :: gg(3,3)
  CCTK_REAL    :: alfa   
  CCTK_REAL    :: beta(3)
  CCTK_REAL    :: g4(0:3,0:3)
  CCTK_COMPLEX :: mm(0:3)
  CCTK_REAL    :: rr(3), rr_p(3), rr_p_p(3)
  CCTK_REAL    :: Ttilde(0:3)
  CCTK_REAL    :: a
  CCTK_COMPLEX :: d

  CCTK_REAL    :: t0, t1, t2
  logical      :: ce0, ce1, ce2
  CCTK_REAL    :: delta_space(2)

  INTEGER      :: i, j
  CCTK_REAL    :: theta, phi
  
  IF (veryverbose/=0) THEN
     CALL CCTK_INFO ("Computing 3-Volume Element of Surface")
  END IF
  
  t0 = qlm_em_time(hn)
  t1 = qlm_em_time_p(hn)
  t2 = qlm_em_time_p_p(hn)
  
  ce0 = qlm_em_have_valid_data(hn) == 0
  ce1 = qlm_em_have_valid_data_p(hn) == 0
  ce2 = qlm_em_have_valid_data_p_p(hn) == 0
  
  delta_space(:) = (/ qlm_em_delta_theta(hn), qlm_em_delta_phi(hn) /)
  
  ! Calculate the coordinates
  DO j = 1+qlm_em_nghostsphi(hn), qlm_em_nphi(hn)-qlm_em_nghostsphi(hn)
     DO i = 1+qlm_em_nghoststheta(hn), qlm_em_ntheta(hn)-qlm_em_nghoststheta(hn)
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
        
        alfa = qlm_em_alpha(i,j)
        
        beta(1) = qlm_em_betax(i,j)
        beta(2) = qlm_em_betay(i,j)
        beta(3) = qlm_em_betaz(i,j)
        
        mm(0) = qlm_em_m0(i,j,hn)
        mm(1) = qlm_em_m1(i,j,hn)
        mm(2) = qlm_em_m2(i,j,hn)
        mm(3) = qlm_em_m3(i,j,hn)
        
        ! Build the four-metric
        CALL calc_4metric (gg,alfa,beta, g4)
        
        ! Find a 3rd vector of triad
        rr(1) = qlm_em_x(i,j,hn)
        rr(2) = qlm_em_y(i,j,hn)
        rr(3) = qlm_em_z(i,j,hn)
        
        rr_p(1) = qlm_em_x_p(i,j,hn)
        rr_p(2) = qlm_em_y_p(i,j,hn)
        rr_p(3) = qlm_em_z_p(i,j,hn)
        
        rr_p_p(1) = qlm_em_x_p_p(i,j,hn)
        rr_p_p(2) = qlm_em_y_p_p(i,j,hn)
        rr_p_p(3) = qlm_em_z_p_p(i,j,hn)
        
        Ttilde(0)   = 1
        Ttilde(1:3) = timederiv (rr, rr_p, rr_p_p, t0,t1,t2, ce0,ce1,ce2)
        
        ! Compute some scalar products
        a=SUM(MATMUL(g4, Ttilde)*Ttilde)
        d=SUM(MATMUL(g4, mm)*Ttilde)
        
        ! This is the determinant
        qlm_em_3det(i,j,hn)=a-2*CONJG(d)*d
        
     END DO
  END DO

  CALL em_set_boundary (CCTK_PASS_FTOF, hn, qlm_em_3det(:,:,hn), +1)

END SUBROUTINE qlm_em_calc_3determinant

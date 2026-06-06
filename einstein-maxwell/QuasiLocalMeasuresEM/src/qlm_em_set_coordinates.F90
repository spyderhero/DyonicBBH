#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"



subroutine qlm_em_set_coordinates (CCTK_ARGUMENTS, hn)
  use cctk
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  
  CCTK_REAL :: theta, phi
  CCTK_REAL :: sin_theta, cos_theta, sin_phi, cos_phi
  
  integer :: i, j
  
  ! Calculate the coordinates
  do j = 1, qlm_em_nphi(hn)
     do i = 1, qlm_em_ntheta(hn)
        theta = qlm_em_origin_theta(hn) + (i-1)*qlm_em_delta_theta(hn)
        phi   = qlm_em_origin_phi(hn)   + (j-1)*qlm_em_delta_phi(hn)
        
        sin_theta = sin(theta)
        cos_theta = cos(theta)
        sin_phi   = sin(phi)
        cos_phi   = cos(phi)
        
        qlm_em_x(i,j,hn) = qlm_em_origin_x(hn) + qlm_em_shape(i,j,hn) * sin_theta * cos_phi
        qlm_em_y(i,j,hn) = qlm_em_origin_y(hn) + qlm_em_shape(i,j,hn) * sin_theta * sin_phi
        qlm_em_z(i,j,hn) = qlm_em_origin_z(hn) + qlm_em_shape(i,j,hn) * cos_theta
        
        qlm_em_x_p(i,j,hn) = qlm_em_origin_x_p(hn) + qlm_em_shape_p(i,j,hn) * sin_theta * cos_phi
        qlm_em_y_p(i,j,hn) = qlm_em_origin_y_p(hn) + qlm_em_shape_p(i,j,hn) * sin_theta * sin_phi
        qlm_em_z_p(i,j,hn) = qlm_em_origin_z_p(hn) + qlm_em_shape_p(i,j,hn) * cos_theta
        
        qlm_em_x_p_p(i,j,hn) = qlm_em_origin_x_p_p(hn) + qlm_em_shape_p_p(i,j,hn) * sin_theta * cos_phi
        qlm_em_y_p_p(i,j,hn) = qlm_em_origin_y_p_p(hn) + qlm_em_shape_p_p(i,j,hn) * sin_theta * sin_phi
        qlm_em_z_p_p(i,j,hn) = qlm_em_origin_z_p_p(hn) + qlm_em_shape_p_p(i,j,hn) * cos_theta
     end do
  end do
  
end subroutine qlm_em_set_coordinates

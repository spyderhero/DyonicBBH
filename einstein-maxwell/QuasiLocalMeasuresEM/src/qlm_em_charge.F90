#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"


subroutine qlm_em_compute_charge (CCTK_ARGUMENTS, hn)
  use cctk
  use constants
  use qlm_em_derivs
  use qlm_em_variables
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: hn
  integer :: i, j
  CCTK_REAL :: alpha, gg(3,3)
  CCTK_REAL :: E_f(3), B_f(3)
  CCTK_REAL :: dX_dtheta(3), dX_dphi(3), dS(3)
  CCTK_REAL :: sqrtgamma, gamma_det

  CCTK_REAL charge_electric_local, charge_magnetic_local

  character :: msg*1000

  charge_electric_local = 0.0d0
  charge_magnetic_local = 0.0d0

  do j = 1, qlm_em_nphi(hn)-1
    do i = 1, qlm_em_ntheta(hn)-1

      ! Tangent vectors in 3-space
      dX_dtheta(1) = qlm_em_x(i+1,j,hn) - qlm_em_x(i,j,hn)
      dX_dtheta(2) = qlm_em_y(i+1,j,hn) - qlm_em_y(i,j,hn)
      dX_dtheta(3) = qlm_em_z(i+1,j,hn) - qlm_em_z(i,j,hn)

      dX_dphi(1) = qlm_em_x(i,j+1,hn) - qlm_em_x(i,j,hn)
      dX_dphi(2) = qlm_em_y(i,j+1,hn) - qlm_em_y(i,j,hn)
      dX_dphi(3) = qlm_em_z(i,j+1,hn) - qlm_em_z(i,j,hn)

      ! Cross product → surface element vector
      dS(1) = dX_dtheta(2)*dX_dphi(3) - dX_dtheta(3)*dX_dphi(2)
      dS(2) = dX_dtheta(3)*dX_dphi(1) - dX_dtheta(1)*dX_dphi(3)
      dS(3) = dX_dtheta(1)*dX_dphi(2) - dX_dtheta(2)*dX_dphi(1)

      ! Compute sqrt(gamma) from qlm_em_gxx etc.
      alpha = qlm_em_alpha(i,j)
      gg(1,1) = qlm_em_gxx(i,j)
      gg(1,2) = qlm_em_gxy(i,j)
      gg(1,3) = qlm_em_gxz(i,j)
      gg(2,2) = qlm_em_gyy(i,j)
      gg(2,3) = qlm_em_gyz(i,j)
      gg(3,3) = qlm_em_gzz(i,j)
      gg(2,1) = gg(1,2)
      gg(3,1) = gg(1,3)
      gg(3,2) = gg(2,3)
      gamma_det = gg(1,1)*gg(2,2)*gg(3,3) &
                  + 2.0*gg(1,2)*gg(2,3)*gg(1,3) &
                  - gg(1,1)*gg(2,3)**2 &
                  - gg(2,2)*gg(1,3)**2 &
                  - gg(3,3)*gg(1,2)**2

      sqrtgamma = sqrt(gamma_det)

      ! Multiply surface element by sqrt(gamma)
      dS = dS * sqrtgamma

      ! Electric field at this point
      E_f(1) = qlm_em_ex(i,j)
      E_f(2) = qlm_em_ey(i,j)
      E_f(3) = qlm_em_ez(i,j)

      ! Magnetic field at this point
      B_f(1) = qlm_em_bx(i,j)
      B_f(2) = qlm_em_by(i,j)
      B_f(3) = qlm_em_bz(i,j)

      ! Flux contribution
      charge_electric_local = charge_electric_local + (E_f(1)*dS(1) + E_f(2)*dS(2) + E_f(3)*dS(3))
      charge_magnetic_local = charge_magnetic_local + (B_f(1)*dS(1) + B_f(2)*dS(2) + B_f(3)*dS(3))

    end do
  end do

  
  ! Divide by 4π
  qlm_em_electric_charge(hn) = charge_electric_local / (4.0*pi)
  qlm_em_magnetic_charge(hn) = charge_magnetic_local / (4.0*pi)

  write (msg, '("   Electric charge Qe:            ",g14.6)') qlm_em_electric_charge(hn)
  call CCTK_INFO (msg)
  write (msg, '("   Magnetic charge Qm:            ",g14.6)') qlm_em_magnetic_charge(hn)
  call CCTK_INFO (msg)

end subroutine qlm_em_compute_charge

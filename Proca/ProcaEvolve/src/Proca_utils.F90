! EinsteinMaxwell_utils_EM_direct_corrected.F90
!
! Utility routines for derivative transformations used by the EM-direct
! Einstein-Maxwell thorn.  These routines are intentionally kept with the
! original Proca_* names because the corrected RHS file calls these
! symbols directly.
!
! The routines transform derivative indices under a MultiPatch Jacobian.
! They do not transform the tensor/component index of vector grid functions.
! This matches the convention used by the original Proca utilities and
! by Proca_calc_rhs_EM_direct_corrected.F90, where Ex,Ey,Ez and
! Bx,By,Bz are evolved as the stored spatial components and only their
! coordinate derivatives are remapped.
!
! For the direct Einstein-Maxwell system based on arXiv:1205.1063, the RHS
! requires only first derivatives of scalars, vectors, and the conformal metric.
! The second-derivative utilities are retained for compatibility with existing
! build/link configurations, but the EM-direct RHS does not rely on them.

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Proca_d1_Scalar_apply_jacobian(dvar, jac)
  implicit none

  CCTK_REAL, intent(inout) :: dvar(3)
  CCTK_REAL, intent(in)    :: jac(3,3)
  CCTK_REAL                :: xdvar(3)
  CCTK_INT                 :: a, b

  xdvar = 0.0d0
  do a = 1, 3
     do b = 1, 3
        xdvar(a) = xdvar(a) + dvar(b) * jac(b,a)
     end do
  end do

  dvar = xdvar

end subroutine Proca_d1_Scalar_apply_jacobian
!
!=============================================================================
!
subroutine Proca_d2_Scalar_apply_jacobian(dvar, ddvar, jac, hes)
  implicit none

  CCTK_REAL, intent(inout) :: ddvar(3,3), dvar(3)
  CCTK_REAL, intent(in)    :: jac(3,3), hes(3,3,3)
  CCTK_REAL                :: xddvar(3,3), xdvar(3)
  CCTK_INT                 :: a, b, c, d

  xdvar  = 0.0d0
  xddvar = 0.0d0

  do a = 1, 3
     do b = 1, 3
        xdvar(a) = xdvar(a) + dvar(b) * jac(b,a)
     end do
  end do

  do b = 1, 3
     do c = 1, 3
        xddvar(1,b) = xddvar(1,b) + dvar(c) * hes(c,1,b)
        do d = 1, 3
           xddvar(1,b) = xddvar(1,b) + ddvar(c,d) * jac(c,1) * jac(d,b)
        end do
     end do
  end do

  do b = 2, 3
     do c = 1, 3
        xddvar(2,b) = xddvar(2,b) + dvar(c) * hes(c,2,b)
        do d = 1, 3
           xddvar(2,b) = xddvar(2,b) + ddvar(c,d) * jac(c,2) * jac(d,b)
        end do
     end do
  end do

  do c = 1, 3
     xddvar(3,3) = xddvar(3,3) + dvar(c) * hes(c,3,3)
     do d = 1, 3
        xddvar(3,3) = xddvar(3,3) + ddvar(c,d) * jac(c,3) * jac(d,3)
     end do
  end do

  xddvar(2,1) = xddvar(1,2)
  xddvar(3,1) = xddvar(1,3)
  xddvar(3,2) = xddvar(2,3)

  dvar  = xdvar
  ddvar = xddvar

end subroutine Proca_d2_Scalar_apply_jacobian
!
!=============================================================================
!
subroutine Proca_d1_Vector_apply_jacobian(dvar, jac)
  implicit none

  CCTK_REAL, intent(inout) :: dvar(3,3)
  CCTK_REAL, intent(in)    :: jac(3,3)
  CCTK_REAL                :: xdvar(3,3)
  CCTK_INT                 :: a, b, c

  xdvar = 0.0d0

  do a = 1, 3
     do b = 1, 3
        do c = 1, 3
           xdvar(a,b) = xdvar(a,b) + dvar(a,c) * jac(c,b)
        end do
     end do
  end do

  dvar = xdvar

end subroutine Proca_d1_Vector_apply_jacobian
!
!=============================================================================
!
subroutine Proca_d2_Vector_apply_jacobian(dvar, ddvar, jac, hes)
  implicit none

  CCTK_REAL, intent(inout) :: ddvar(3,3,3), dvar(3,3)
  CCTK_REAL, intent(in)    :: jac(3,3), hes(3,3,3)
  CCTK_REAL                :: xddvar(3,3,3), xdvar(3,3)
  CCTK_INT                 :: a, b, c, d, e

  xdvar  = 0.0d0
  xddvar = 0.0d0

  do a = 1, 3
     do b = 1, 3
        do c = 1, 3
           xdvar(a,b) = xdvar(a,b) + dvar(a,c) * jac(c,b)
        end do
     end do
  end do

  do a = 1, 3
     do c = 1, 3
        do d = 1, 3
           xddvar(a,1,c) = xddvar(a,1,c) + dvar(a,d) * hes(d,1,c)
           do e = 1, 3
              xddvar(a,1,c) = xddvar(a,1,c) + ddvar(a,d,e) * jac(d,1) * jac(e,c)
           end do
        end do
     end do

     do c = 2, 3
        do d = 1, 3
           xddvar(a,2,c) = xddvar(a,2,c) + dvar(a,d) * hes(d,2,c)
           do e = 1, 3
              xddvar(a,2,c) = xddvar(a,2,c) + ddvar(a,d,e) * jac(d,2) * jac(e,c)
           end do
        end do
     end do

     do d = 1, 3
        xddvar(a,3,3) = xddvar(a,3,3) + dvar(a,d) * hes(d,3,3)
        do e = 1, 3
           xddvar(a,3,3) = xddvar(a,3,3) + ddvar(a,d,e) * jac(d,3) * jac(e,3)
        end do
     end do
  end do

  xddvar(:,2,1) = xddvar(:,1,2)
  xddvar(:,3,1) = xddvar(:,1,3)
  xddvar(:,3,2) = xddvar(:,2,3)

  dvar  = xdvar
  ddvar = xddvar

end subroutine Proca_d2_Vector_apply_jacobian
!
!=============================================================================
!
subroutine Proca_d1_2nd_Tensor_apply_jacobian(dvar, jac, sym)
  implicit none

  CCTK_REAL, intent(inout) :: dvar(3,3,3)
  CCTK_REAL, intent(in)    :: jac(3,3)
  CCTK_INT,  intent(in)    :: sym
  CCTK_REAL                :: xdvar(3,3,3)
  CCTK_INT                 :: a, b, c, d

  xdvar = 0.0d0

  if (sym == 1) then
     ! Symmetric rank-2 object in the first two indices, e.g. d_k h_ij.
     do b = 1, 3
        do c = 1, 3
           do d = 1, 3
              xdvar(1,b,c) = xdvar(1,b,c) + dvar(1,b,d) * jac(d,c)
           end do
        end do
     end do

     do b = 2, 3
        do c = 1, 3
           do d = 1, 3
              xdvar(2,b,c) = xdvar(2,b,c) + dvar(2,b,d) * jac(d,c)
           end do
        end do
     end do

     do c = 1, 3
        do d = 1, 3
           xdvar(3,3,c) = xdvar(3,3,c) + dvar(3,3,d) * jac(d,c)
        end do
     end do

     xdvar(2,1,:) = xdvar(1,2,:)
     xdvar(3,1,:) = xdvar(1,3,:)
     xdvar(3,2,:) = xdvar(2,3,:)

  else if (sym == 0) then
     do a = 1, 3
        do b = 1, 3
           do c = 1, 3
              do d = 1, 3
                 xdvar(a,b,c) = xdvar(a,b,c) + dvar(a,b,d) * jac(d,c)
              end do
           end do
        end do
     end do

  else
     call CCTK_WARN(0, "invalid parameters for symmetry.")
  end if

  dvar = xdvar

end subroutine Proca_d1_2nd_Tensor_apply_jacobian

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

!=============================================================================
! Boundary registration for the Einstein-Maxwell direct E/B evolution system.
!
! This thorn evolves the electromagnetic variables directly as
!
!   E^i  = (Ex, Ey, Ez),
!   B^i  = (Bx, By, Bz),
!   Psi  = electric Gauss-constraint damping field,
!   Phi  = magnetic Gauss-constraint damping field.
!
! The actual outgoing radiative boundary update is applied in
! Proca_calc_rhs_bdry through NewRad_Apply.  Here we register the
! variables with the Boundary thorn using "none", so that symmetry boundary
! conditions are still enforced and no additional physical boundary condition
! is imposed on top of NewRad.
!=============================================================================

subroutine Proca_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr
  CCTK_INT, parameter :: one = 1

  !--------------------------------------------------------------------------
  ! Dynamical electromagnetic variables.
  !
  ! These group/variable names must agree with interface.ccl.  They are chosen
  ! to match the corrected direct E/B RHS file, where the variables are
  ! Ex,Ey,Ez,Bx,By,Bz,Psi,Phi and the RHS boundary routine applies NewRad to
  ! exactly these fields.
  !--------------------------------------------------------------------------

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,       &
       "ProcaBase::Ei", "none")
  if (ierr < 0)                                                             &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Ei!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,       &
       "ProcaBase::Bi", "none")
  if (ierr < 0)                                                             &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Bi!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,         &
       "ProcaBase::Psi", "none")
  if (ierr < 0)                                                             &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Psi!")

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, one, -one,         &
       "ProcaBase::Phi", "none")
  if (ierr < 0)                                                             &
       call CCTK_WARN(0, "Failed to register BC for ProcaBase::Phi!")

  !--------------------------------------------------------------------------
  ! Electromagnetic diagnostic variables written by the corrected Tmunu file.
  ! They do not need physical radiative BCs; registering "none" is enough to
  ! keep symmetry BC handling consistent.
  !--------------------------------------------------------------------------

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,       &
       "ProcaEvolve::EM_jr_F_gfs", "none")
  if (ierr < 0)                                                             &
       call CCTK_WARN(0, "Failed to register BC for ProcaEvolve::EM_jr_F_gfs!")

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,       &
       "ProcaEvolve::EM_Sir_F_gfs", "none")
  if (ierr < 0)                                                             &
       call CCTK_WARN(0, "Failed to register BC for ProcaEvolve::EM_Sir_F_gfs!")

end subroutine Proca_Boundaries

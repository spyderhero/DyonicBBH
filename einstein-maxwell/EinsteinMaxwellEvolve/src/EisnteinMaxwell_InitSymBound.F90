#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

!=============================================================================
! Symmetry registration for the Einstein-Maxwell direct E/B evolution system.
!
! This file is the direct-E/B analogue of the original EisnteinMaxwell
! InitSymBound routine.  It is consistent with the corrected RHS and Tmunu
! routines that evolve/use
!
!   E^i  = (Ex, Ey, Ez),
!   B^i  = (Bx, By, Bz),
!   Psi  = electric Gauss-constraint damping field,
!   Phi  = magnetic Gauss-constraint damping field.
!
! No vector potential A_i, scalar potential Aphi, Zeta field, dilaton, or
! scalar-field variables are registered here.
!
! Parity convention under Cartesian reflections:
!   - components of ordinary spatial vectors E^i and B^i have vector parity,
!     e.g. Vx is odd under x-reflection and even under y,z-reflections.
!   - scalar damping fields Psi and Phi are even under all reflections.
!   - radial scalar diagnostics jrF_gf are even.
!   - SirF_gf components transform as ordinary spatial covector/vector
!     components in Cartesian coordinates.
!=============================================================================

subroutine EisnteinMaxwell_InitSymBound( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  !--------------------------------------------------------------------------
  ! RHS variables for the direct electromagnetic evolution.
  ! These names must agree with the corrected RHS routine:
  !   EisnteinMaxwell_calc_rhs_EM_direct_corrected.F90
  !--------------------------------------------------------------------------

  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "EisnteinMaxwellEvolve::rhs_Ex"  )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "EisnteinMaxwellEvolve::rhs_Ey"  )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "EisnteinMaxwellEvolve::rhs_Ez"  )

  call SetCartSymVN( ierr, cctkGH, (/-1,  1,  1/), "EisnteinMaxwellEvolve::rhs_Bx"  )
  call SetCartSymVN( ierr, cctkGH, (/ 1, -1,  1/), "EisnteinMaxwellEvolve::rhs_By"  )
  call SetCartSymVN( ierr, cctkGH, (/ 1,  1, -1/), "EisnteinMaxwellEvolve::rhs_Bz"  )
  
  call SetCartSymVN( ierr, cctkGH, (/ 1,  1,  1/), "EisnteinMaxwellEvolve::rhs_Psi" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,  1,  1/), "EisnteinMaxwellEvolve::rhs_Phi" )

  !--------------------------------------------------------------------------
  ! Optional EM diagnostics written by the corrected Tmunu routine:
  !   EisnteinMaxwell_calc_Tmunu_EM_direct_corrected.F90
  ! Keep only the electromagnetic F-sector diagnostics.  The old EMD
  ! old scalar-sector diagnostics are intentionally not registered.
  !--------------------------------------------------------------------------

  call SetCartSymVN( ierr, cctkGH, (/ 1,  1,  1/), "EisnteinMaxwellEvolve::jrF_gf"  )

  call SetCartSymVN( ierr, cctkGH, (/-1,  1,  1/), "EisnteinMaxwellEvolve::SxrF_gf" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, -1,  1/), "EisnteinMaxwellEvolve::SyrF_gf" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,  1, -1/), "EisnteinMaxwellEvolve::SzrF_gf" )

end subroutine EisnteinMaxwell_InitSymBound

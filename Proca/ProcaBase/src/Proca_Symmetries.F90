! Proca_Symmetries.F90
!
! Register Cartesian symmetries for the Einstein-Maxwell base variables.
! The thorn name ProcaBase is kept for compatibility, but the variables are
! E^i, B^i, Psi and Phi, not a Proca/vector-potential system.
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine einsteinmaxwellbase_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

  ! E^i: contravariant Cartesian vector components.
  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaBase::Ex" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaBase::Ey" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaBase::Ez" )

  ! B^i is evolved directly.  For the radial monopole-type initial data used
  ! in the present project, the component parities are the same Cartesian
  ! vector parities as for E^i.
  call SetCartSymVN( ierr, cctkGH, (/-1, 1, 1/), "ProcaBase::Bx" )
  call SetCartSymVN( ierr, cctkGH, (/ 1,-1, 1/), "ProcaBase::By" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "ProcaBase::Bz" )

  ! Constraint-damping scalars.
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaBase::Psi" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaBase::Phi" )

end subroutine einsteinmaxwellbase_symmetries

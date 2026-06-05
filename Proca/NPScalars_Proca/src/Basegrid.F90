! NPScalars_Proca -- Einstein-Maxwell direct E/B version
! Basegrid.F90 : Register symmetries
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine NPProca_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT  ierr

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "NPScalars_Proca::psi0re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1,-1/), "NPScalars_Proca::psi0im" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "NPScalars_Proca::psi4re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1,-1/), "NPScalars_Proca::psi4im" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "NPScalars_Proca::phi0re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "NPScalars_Proca::phi0im" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "NPScalars_Proca::phi1re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1,-1/), "NPScalars_Proca::phi1im" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "NPScalars_Proca::phi2re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "NPScalars_Proca::phi2im" )

end subroutine NPProca_symmetries

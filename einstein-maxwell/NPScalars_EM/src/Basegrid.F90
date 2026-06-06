! NPScalars_EM
! Basegrid.F90 : Register symmetries
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

subroutine NPEM_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT  ierr

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "NPScalars_EM::phi0re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "NPScalars_EM::phi0im" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "NPScalars_EM::phi1re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1,-1/), "NPScalars_EM::phi1im" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,-1/), "NPScalars_EM::phi2re" )
  call SetCartSymVN( ierr, cctkGH, (/-1,-1, 1/), "NPScalars_EM::phi2im" )

end subroutine NPEM_symmetries

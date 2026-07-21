
#include "cctk.h"
#include "cctk_Arguments.h"

subroutine EinsteinMaxwell_zero_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS

  rhs_Ex    = 0
  rhs_Ey    = 0
  rhs_Ez    = 0

  rhs_Bx    = 0
  rhs_By    = 0
  rhs_Bz    = 0

  rhs_Psi  = 0
  rhs_Phi  = 0

end subroutine EinsteinMaxwell_zero_rhs

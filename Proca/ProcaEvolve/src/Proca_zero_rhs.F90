#include "cctk.h"
#include "cctk_Arguments.h"

subroutine Proca_zero_rhs( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS

  rhs_Ex  = 0.0d0
  rhs_Ey  = 0.0d0
  rhs_Ez  = 0.0d0

  rhs_Bx  = 0.0d0
  rhs_By  = 0.0d0
  rhs_Bz  = 0.0d0

  rhs_Psi = 0.0d0
  rhs_Phi = 0.0d0

end subroutine Proca_zero_rhs

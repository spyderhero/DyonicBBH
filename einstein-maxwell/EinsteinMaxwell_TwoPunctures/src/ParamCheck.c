/* TwoPunctures:  File  "TwoPunctures.c"*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

/* -------------------------------------------------------------------*/
void
EinsteinMaxwell_TwoPunctures_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (use_sources)
  {
    if (! CCTK_IsFunctionAliased ("Set_Rho_ADM"))
      CCTK_WARN (0, "Matter sources have been enabled, "
                 "but there is no aliased function for matter sources.");
  }

  /*
   * For a purely electric binary, a purely magnetic binary, or a globally
   * duality-rotated dyonic binary with q_plus/p_plus = q_minus/p_minus,
   * E and B are parallel everywhere in the conformal flat approximation and
   * the electromagnetic momentum density j_i ~ epsilon_ijk E^j B^k vanishes.
   *
   * For the mixed electric/magnetic binaries targeted by the new research
   * programme, q_plus*p_minus - q_minus*p_plus is generically nonzero.  Then
   * E x B is nonzero and the electromagnetic momentum constraint source is
   * nonzero.  This thorn still follows the original TwoPunctures workflow and
   * solves only the Hamiltonian constraint with the EM energy density source.
   * Therefore such data are useful for exploratory runs only if the user
   * knowingly accepts the residual momentum-constraint error.
   */
  if (! solve_momentum_constraint)
  {
    const CCTK_REAL duality_mismatch = par_q_plus * par_p_minus
                                    - par_q_minus * par_p_plus;
    if (fabs(duality_mismatch) > 1.0e-14)
    {
      if (! allow_nonzero_em_momentum_source)
      {
        CCTK_WARN(0, "Generic mixed electric/magnetic charges give nonzero E cross B momentum density, but this thorn does not solve the corresponding momentum constraint.  Set allow_nonzero_em_momentum_source=yes only if you intentionally accept this limitation.");
      }
      else if (verbose)
      {
        CCTK_WARN(1, "Proceeding with nonzero E cross B momentum density without solving the EM momentum constraint, because allow_nonzero_em_momentum_source=yes.");
      }
    }
  }
}

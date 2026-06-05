/* TwoPunctures_KerrEinsteinMaxwell:  File  "TwoPunctures_KerrEinsteinMaxwell.c"*/

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
TPKP_TwoPunctures_ParamCheck (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (use_sources)
  {
    if (! CCTK_IsFunctionAliased ("Set_Rho_ADM"))
      CCTK_WARN (0, "Matter sources have been enabled, "
                 "but there is no aliased function for matter sources.");
  }


  if( par_S_plus[1] != 0 || par_S_plus[2] != 0 || par_S_minus[0] != 0 || par_S_minus[1] != 0 || par_S_minus[2] != 0 )
    { 
        CCTK_WARN( 0, "Initialized wrong spin components. TwoPunctures_KerrEinsteinMaxwell expects spin in x direction, i.e., only par_S_plus[0] /= 0 ");
    }


  /* Mixed electric/magnetic charges generally produce a nonzero EM momentum
   * density j_i ~ epsilon_ijk E^j B^k.  This initial-data ansatz solves the
   * Hamiltonian correction but does not fully solve a coupled EM momentum
   * constraint with that source.
   */
  if (! allow_nonzero_em_momentum_source)
  {
    CCTK_REAL mismatch = fabs(par_q_plus * par_p_minus - par_q_minus * par_p_plus);
    if (mismatch > 1.0e-14)
      CCTK_WARN(1, "Mixed electric/magnetic charges imply nonzero E cross B momentum density; set allow_nonzero_em_momentum_source=yes after verifying the initial-data assumptions.");
  }


}

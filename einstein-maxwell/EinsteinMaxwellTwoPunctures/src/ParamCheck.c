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
EinsteinMaxwell_TwoPunctures_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (use_sources)
  {
    if (!CCTK_IsFunctionAliased("Set_Rho_ADM"))
    {
      CCTK_WARN(
        0,
        "Matter sources have been enabled, but there is no "
        "aliased function for matter sources."
      );
    }
  }

  {
    const CCTK_REAL Qp = par_electric_charge_plus;
    const CCTK_REAL Qm = par_electric_charge_minus;
    const CCTK_REAL Pp = par_magnetic_charge_plus;
    const CCTK_REAL Pm = par_magnetic_charge_minus;

    const CCTK_REAL mass_product =
        target_M_plus * target_M_minus;

    CCTK_REAL C_em;
    CCTK_REAL D_em;
    CCTK_REAL em_coulomb_factor;
    CCTK_REAL em_cross_factor;

    C_em = Qp*Qm + Pp*Pm;
    D_em = Qp*Pm - Pp*Qm;

    em_coulomb_factor = C_em / mass_product;
    em_cross_factor   = D_em / mass_product;

    if (output_duality_diagnostics)
    {
      CCTK_VInfo(
        CCTK_THORNSTRING,
        "EM charge diagnostics: "
        "C/(M_plus M_minus)=%.16e, "
        "D/(M_plus M_minus)=%.16e",
        (double)em_coulomb_factor,
        (double)em_cross_factor
      );
    }

    if (!solve_momentum_constraint
        && fabs(em_cross_factor) > 1.0e-14)
    {
      if (!allow_nonzero_em_momentum_source)
      {
        CCTK_WARN(
          0,
          "Nonzero em_cross_factor implies a nonzero "
          "electromagnetic momentum source, but the momentum "
          "constraint is not solved."
        );
      }
      else if (verbose)
      {
        CCTK_WARN(
          1,
          "Proceeding with nonzero electromagnetic momentum "
          "density without solving the momentum constraint."
        );
      }
    }
  }
}
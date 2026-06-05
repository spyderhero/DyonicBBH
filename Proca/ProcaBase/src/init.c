#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

void ProcaBase_Zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const np = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];

#pragma omp parallel for
  for (int i = 0; i < np; ++i) {
    Ex [i] = 0.0;
    Ey [i] = 0.0;
    Ez [i] = 0.0;
    Bx [i] = 0.0;
    By [i] = 0.0;
    Bz [i] = 0.0;
    Psi[i] = 0.0;
    Phi[i] = 0.0;
  }

  if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
    /* do nothing */
  } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {
      if (CCTK_ActiveTimeLevels(cctkGH, "ProcaBase::Ei") >= 2) {
#pragma omp parallel for
          for (int i = 0; i < np; ++i) {
              Ex_p [i] = 0.0;
              Ey_p [i] = 0.0;
              Ez_p [i] = 0.0;
              Bx_p [i] = 0.0;
              By_p [i] = 0.0;
              Bz_p [i] = 0.0;
              Psi_p[i] = 0.0;
              Phi_p[i] = 0.0;
          }
      }

      if (CCTK_ActiveTimeLevels(cctkGH, "ProcaBase::Ei") >= 3) {
#pragma omp parallel for
          for (int i = 0; i < np; ++i) {
              Ex_p_p [i] = 0.0;
              Ey_p_p [i] = 0.0;
              Ez_p_p [i] = 0.0;
              Bx_p_p [i] = 0.0;
              By_p_p [i] = 0.0;
              Bz_p_p [i] = 0.0;
              Psi_p_p[i] = 0.0;
              Phi_p_p[i] = 0.0;
          }
      }
  } else {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "Unsupported parameter value for InitBase::initial_data_setup_method");
  }
}

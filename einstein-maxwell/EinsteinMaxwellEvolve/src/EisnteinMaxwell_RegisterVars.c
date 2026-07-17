#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void EinsteinMaxwell_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0;
  CCTK_INT group = -1;
  CCTK_INT rhs   = -1;
  CCTK_INT var   = -1;

  /*
   * MoL registration for the Einstein-Maxwell system evolved with
   * E^i, B^i, Psi and Phi directly, without using a vector potential.
   *
   * This is consistent with EinsteinMaxwell_calc_rhs_EM_direct_corrected.F90:
   *   ProcaBase::Ei              <-> EinsteinMaxwellEvolve::rhs_Ei
   *   ProcaBase::Bi              <-> EinsteinMaxwellEvolve::rhs_Bi
   *   ProcaBase::Psi             <-> EinsteinMaxwellEvolve::rhs_Psi
   *   ProcaBase::Phi             <-> EinsteinMaxwellEvolve::rhs_Phi
   *
   * The old EMD/vector-potential variables Ai, Aphi, Zeta, phi1, phi2,
   * Kphi1 and Kphi2 are intentionally not registered here.
   */

  /* Save-and-restore ADMBase variables that are evolved by spacetime thorns. */
  group = CCTK_GroupIndex("ADMBase::lapse");
  if (group < 0) CCTK_ERROR("Could not find ADMBase::lapse");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  group = CCTK_GroupIndex("ADMBase::shift");
  if (group < 0) CCTK_ERROR("Could not find ADMBase::shift");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  group = CCTK_GroupIndex("ADMBase::metric");
  if (group < 0) CCTK_ERROR("Could not find ADMBase::metric");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  group = CCTK_GroupIndex("ADMBase::curv");
  if (group < 0) CCTK_ERROR("Could not find ADMBase::curv");
  ierr += MoLRegisterSaveAndRestoreGroup(group);

  /* Electric field E^i and its RHS. */
  group = CCTK_GroupIndex("ProcaBase::Ei");
  rhs   = CCTK_GroupIndex("EinsteinMaxwellEvolve::rhs_Ei");
  if (group < 0) CCTK_ERROR("Could not find ProcaBase::Ei");
  if (rhs   < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::rhs_Ei");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Magnetic field B^i and its RHS. */
  group = CCTK_GroupIndex("ProcaBase::Bi");
  rhs   = CCTK_GroupIndex("EinsteinMaxwellEvolve::rhs_Bi");
  if (group < 0) CCTK_ERROR("Could not find ProcaBase::Bi");
  if (rhs   < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::rhs_Bi");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Electric Gauss-constraint damping field Psi and its RHS. */
  var = CCTK_VarIndex("ProcaBase::Psi");
  rhs = CCTK_VarIndex("EinsteinMaxwellEvolve::rhs_Psi");
  if (var < 0) CCTK_ERROR("Could not find ProcaBase::Psi");
  if (rhs < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::rhs_Psi");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Magnetic Gauss-constraint damping field Phi and its RHS. */
  var = CCTK_VarIndex("ProcaBase::Phi");
  rhs = CCTK_VarIndex("EinsteinMaxwellEvolve::rhs_Phi");
  if (var < 0) CCTK_ERROR("Could not find ProcaBase::Phi");
  if (rhs < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::rhs_Phi");
  ierr += MoLRegisterEvolved(var, rhs);

  if (ierr) CCTK_ERROR("Problems registering Einstein-Maxwell variables with MoL");
}

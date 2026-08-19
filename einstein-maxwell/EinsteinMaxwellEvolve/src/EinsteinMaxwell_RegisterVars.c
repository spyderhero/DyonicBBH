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
   * E^i, B^i, EM_Psi and EM_Phi directly, without using a vector potential.
   *
   * Variable ownership:
   *   ProcaBase::Ei                <-> EinsteinMaxwellEvolve::rhs_Ei
   *   EinsteinMaxwellEvolve::Bi    <-> EinsteinMaxwellEvolve::rhs_Bi
   *   EinsteinMaxwellEvolve::EM_Psi   <-> EinsteinMaxwellEvolve::rhs_EM_Psi
   *   EinsteinMaxwellEvolve::EM_Phi   <-> EinsteinMaxwellEvolve::rhs_EM_Phi
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
  group = CCTK_GroupIndex("EinsteinMaxwellEvolve::Bi");
  rhs   = CCTK_GroupIndex("EinsteinMaxwellEvolve::rhs_Bi");
  if (group < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::Bi");
  if (rhs   < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::rhs_Bi");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* Electric Gauss-constraint damping field EM_Psi and its RHS. */
  var = CCTK_VarIndex("EinsteinMaxwellEvolve::EM_Psi");
  rhs = CCTK_VarIndex("EinsteinMaxwellEvolve::rhs_EM_Psi");
  if (var < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::EM_Psi");
  if (rhs < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::rhs_EM_Psi");
  ierr += MoLRegisterEvolved(var, rhs);

  /* Magnetic Gauss-constraint damping field EM_Phi and its RHS. */
  var = CCTK_VarIndex("EinsteinMaxwellEvolve::EM_Phi");
  rhs = CCTK_VarIndex("EinsteinMaxwellEvolve::rhs_EM_Phi");
  if (var < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::EM_Phi");
  if (rhs < 0) CCTK_ERROR("Could not find EinsteinMaxwellEvolve::rhs_EM_Phi");
  ierr += MoLRegisterEvolved(var, rhs);

  if (ierr)
    CCTK_ERROR("Problems registering Einstein-Maxwell variables with MoL");
}
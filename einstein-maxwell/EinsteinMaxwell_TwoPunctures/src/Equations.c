/* TwoPunctures:  File  "Equations.c"*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

/* U.d0[ivar]   = U[ivar];  (ivar = 0..nvar-1) */
/* U.d1[ivar]   = U[ivar]_x;  */
/* U.d2[ivar]   = U[ivar]_y;  */
/* U.d3[ivar]   = U[ivar]_z;  */
/* U.d11[ivar]  = U[ivar]_xx; */
/* U.d12[ivar]  = U[ivar]_xy; */
/* U.d13[ivar]  = U[ivar]_xz;*/
/* U.d22[ivar]  = U[ivar]_yy;*/
/* U.d23[ivar]  = U[ivar]_yz;*/
/* U.d33[ivar]  = U[ivar]_zz;*/

CCTK_REAL
EinsteinMaxwell_BY_KKofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  CCTK_REAL r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, nplus_parPplus, nminus_parPminus,
    Aij, AijAij, n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r_plus = sqrt (r2_plus);
  r_minus = sqrt (r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: nplus_parPplus = (n_+).(P_+); nminus_parPminus = (n_-).(P_-) */
  nplus_parPplus = 0;
  nminus_parPminus = 0;
  for (i = 0; i < 3; i++)
  {
    nplus_parPplus += n_plus[i] * par_P_plus[i];
    nminus_parPminus += n_minus[i] * par_P_minus[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  AijAij = 0;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				/* Bowen-York-Curvature :*/
      Aij =
	+ 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
                 + nplus_parPplus * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
		 + nminus_parPminus * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij -= +1.5 * (nplus_parPplus / r2_plus + nminus_parPminus / r2_minus);
      AijAij += Aij * Aij;
    }
  }

  return AijAij;
}

void
EinsteinMaxwell_BY_Aijofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL Aij[3][3])
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  CCTK_REAL r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, nplus_parPplus, nminus_parPminus,
    n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par_b) * (x - par_b) + y * y + z * z;
  r2_minus = (x + par_b) * (x + par_b) + y * y + z * z;
  r2_plus = sqrt (pow (r2_plus, 2) + pow (TP_epsilon, 4));
  r2_minus = sqrt (pow (r2_minus, 2) + pow (TP_epsilon, 4));
  if (r2_plus < pow(TP_Tiny,2))
    r2_plus = pow(TP_Tiny,2);
  if (r2_minus < pow(TP_Tiny,2))
    r2_minus = pow(TP_Tiny,2);
  r_plus = sqrt (r2_plus);
  r_minus = sqrt (r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par_b) / r_plus;
  n_minus[0] = (x + par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: nplus_parPplus = (n_+).(P_+); nminus_parPminus = (n_-).(P_-) */
  nplus_parPplus = 0;
  nminus_parPminus = 0;
  for (i = 0; i < 3; i++)
  {
    nplus_parPplus += n_plus[i] * par_P_plus[i];
    nminus_parPminus += n_minus[i] * par_P_minus[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * par_S_plus[2] - n_plus[2] * par_S_plus[1];
  np_Sp[1] = n_plus[2] * par_S_plus[0] - n_plus[0] * par_S_plus[2];
  np_Sp[2] = n_plus[0] * par_S_plus[1] - n_plus[1] * par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par_S_minus[2] - n_minus[2] * par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par_S_minus[0] - n_minus[0] * par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par_S_minus[1] - n_minus[1] * par_S_minus[0];
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {				/* Bowen-York-Curvature :*/
      Aij[i][j] =
        + 1.5 * (par_P_plus[i] * n_plus[j] + par_P_plus[j] * n_plus[i]
		 + nplus_parPplus * n_plus[i] * n_plus[j]) / r2_plus
	+ 1.5 * (par_P_minus[i] * n_minus[j] + par_P_minus[j] * n_minus[i]
		 + nminus_parPminus * n_minus[i] * n_minus[j]) / r2_minus
	- 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus
	- 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
	Aij[i][j] -= +1.5 * (nplus_parPplus / r2_plus + nminus_parPminus / r2_minus);
    }
  }
}

/*-----------------------------------------------------------*/
/********       electric and magnetic monopole fields       **********/
/*-----------------------------------------------------------*/
void
EinsteinMaxwell_Calc_EM_Fields(CCTK_REAL x1, CCTK_REAL y1, CCTK_REAL z1,
                    CCTK_REAL Ei[3], CCTK_REAL Bi[3], CCTK_REAL *F2)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL r_plus, r_minus;
  CCTK_REAL Er_plus, Er_minus;
  CCTK_REAL Br_plus, Br_minus;

  r_plus = sqrt ((x1 - par_b) * (x1 - par_b) + y1 * y1 + z1 * z1);
  r_minus = sqrt ((x1 + par_b) * (x1 + par_b) + y1 * y1 + z1 * z1);
  r_plus = pow (pow (r_plus, 4) + pow (TP_epsilon, 4), 0.25);
  r_minus = pow (pow (r_minus, 4) + pow (TP_epsilon, 4), 0.25);
  if (r_plus < TP_Tiny)
      r_plus = TP_Tiny;
  if (r_minus < TP_Tiny)
      r_minus = TP_Tiny;

  Er_plus  = par_Qp / (r_plus*r_plus);
  Er_minus = par_Qm / (r_minus*r_minus);

  Ei[0] = Er_plus * (x1 - par_b) / r_plus + Er_minus * (x1 + par_b) / r_minus;
  Ei[1] = Er_plus * y1 / r_plus           + Er_minus * y1 / r_minus;
  Ei[2] = Er_plus * z1 / r_plus           + Er_minus * z1 / r_minus;

  Br_plus  = parlus_parPplus / (r_plus*r_plus);
  Br_minus = par_Pm / (r_minus*r_minus);

  Bi[0] = Br_plus  * (x1 - par_b) / r_plus + Br_minus * (x1 + par_b) / r_minus;
  Bi[1] = Br_plus  * y1 / r_plus           + Br_minus * y1 / r_minus;
  Bi[2] = Br_plus  * z1 / r_plus           + Br_minus * z1 / r_minus;

  /* Flat-space norm of the conformal electromagnetic fields used in the
     Hamiltonian constraint. */
  *F2 = Ei[0]*Ei[0] + Ei[1]*Ei[1] + Ei[2]*Ei[2]
      + Bi[0]*Bi[0] + Bi[1]*Bi[1] + Bi[2]*Bi[2];
}

void
EinsteinMaxwell_Calc_E(CCTK_REAL x1, CCTK_REAL y1, CCTK_REAL z1,
                    CCTK_REAL Ei[3], CCTK_REAL *E2)
{
  CCTK_REAL Bi[3], F2;
  EinsteinMaxwell_Calc_EM_Fields(x1, y1, z1, Ei, Bi, &F2);
  *E2 = Ei[0]*Ei[0] + Ei[1]*Ei[1] + Ei[2]*Ei[2];
}

void
EinsteinMaxwell_Calc_B(CCTK_REAL x1, CCTK_REAL y1, CCTK_REAL z1,
                    CCTK_REAL Bi[3], CCTK_REAL *B2)
{
  CCTK_REAL Ei[3], F2;
  EinsteinMaxwell_Calc_EM_Fields(x1, y1, z1, Ei, Bi, &F2);
  *B2 = Bi[0]*Bi[0] + Bi[1]*Bi[1] + Bi[2]*Bi[2];
}

/*---------------------------------------------------------------------------------------*/
/********           Nonlinear Equations  (Hamiltonian Constraint)              ***********/
/*---------------------------------------------------------------------------------------*/
void
EinsteinMaxwell_NonLinEquations (CCTK_REAL rho_adm,
     CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
		 CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
		 CCTK_REAL y, CCTK_REAL z, derivs U, CCTK_REAL *values)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL r_plus, r_minus, psi, psi2, psi4, psi7, psi_ana, psi_ana2;
  CCTK_REAL Ei[3], Bi[3], F2;

  r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);
  r_plus = pow (pow (r_plus, 4) + pow (TP_epsilon, 4), 0.25);
  r_minus = pow (pow (r_minus, 4) + pow (TP_epsilon, 4), 0.25);
  if (r_plus < TP_Tiny)
      r_plus = TP_Tiny;
  if (r_minus < TP_Tiny)
      r_minus = TP_Tiny;

  psi = sqrt( pow(1 + 0.5 * par_m_plus / r_plus
              + 0.5 * par_m_minus/ r_minus , 2)
              - 0.25 * (pow( par_Qp/r_plus + par_Qm/r_minus, 2) + pow( parlus_parPplus/r_plus + par_Pm/r_minus, 2)) ) + U.d0[0];
  psi_ana = sqrt( pow(1 + 0.5 * par_m_plus / r_plus
              + 0.5 * par_m_minus/ r_minus , 2)
              - 0.25 * (pow( par_Qp/r_plus + par_Qm/r_minus, 2) + pow( parlus_parPplus/r_plus + par_Pm/r_minus, 2)) );
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi7 = psi * psi2 * psi4;
  psi_ana2 = psi_ana * psi_ana;

  EinsteinMaxwell_Calc_EM_Fields(x, y, z, Ei, Bi, &F2);

  values[0] =
    U.d11[0] + U.d22[0] + U.d33[0] + 0.125 * EinsteinMaxwell_BY_KKofxyz (x, y, z) / psi7 +
    2.0 * Pi / psi2/psi * rho_adm + 0.25 * F2 / psi2/psi
    - 0.25 * F2 / psi_ana2/psi_ana;

}

/*-----------------------------------------------------------*/
/********               Linear Equations                ***********/
/*-----------------------------------------------------------*/
void
EinsteinMaxwell_LinEquations (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
	      CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
	      CCTK_REAL y, CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL r_plus, r_minus, psi, psi2, psi4, psi8;
  CCTK_REAL Ei[3], Bi[3], F2;

  r_plus = sqrt ((x - par_b) * (x - par_b) + y * y + z * z);
  r_minus = sqrt ((x + par_b) * (x + par_b) + y * y + z * z);
  r_plus = pow (pow (r_plus, 4) + pow (TP_epsilon, 4), 0.25);
  r_minus = pow (pow (r_minus, 4) + pow (TP_epsilon, 4), 0.25);
  if (r_plus < TP_Tiny)
      r_plus = TP_Tiny;
  if (r_minus < TP_Tiny)
      r_minus = TP_Tiny;

  psi = sqrt( pow(1 + 0.5 * par_m_plus / r_plus
              + 0.5 * par_m_minus/ r_minus , 2)
              - 0.25 * (pow( par_Qp/r_plus + par_Qm/r_minus, 2) + pow( parlus_parPplus/r_plus + par_Pm/r_minus, 2)) ) + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi8 = psi4 * psi4;
  
  EinsteinMaxwell_Calc_EM_Fields(x, y, z, Ei, Bi, &F2);

  values[0] = dU.d11[0] + dU.d22[0] + dU.d33[0]
    - 0.875 * EinsteinMaxwell_BY_KKofxyz (x, y, z) / psi8 * dU.d0[0]
    - 0.75 * F2 / psi4 * dU.d0[0];
}

/*-----------------------------------------------------------*/

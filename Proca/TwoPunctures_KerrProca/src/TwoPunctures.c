/* TwoPunctures_KerrEinsteinMaxwell:  File  "TwoPunctures.c"*/

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

/* Swap two variables */
static inline
void swap (CCTK_REAL * restrict const a, CCTK_REAL * restrict const b)
{
  CCTK_REAL const t = *a; *a=*b; *b=t;
}
#undef SWAP
#define SWAP(a,b) (swap(&(a),&(b)))



static
void set_initial_guess(CCTK_POINTER_TO_CONST cctkGH,
                       derivs v)
{
  DECLARE_CCTK_PARAMETERS;

  int nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  CCTK_REAL *s_x, *s_y, *s_z;
  CCTK_REAL al, A, Am1, be, B, phi, R, r, X;
  CCTK_INT ivar, i, j, k, i3D, indx;
  derivs U;
  FILE *debug_file;

  if (solve_momentum_constraint)
    nvar = 4;

  s_x    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  s_y    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  s_z    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  TPKP_allocate_derivs (&U, nvar);
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++)
        {
          i3D = TPKP_Index(ivar,i,j,k,1,n1,n2,n3);

          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 2. * Pi * k / n3;

          /* Calculation of (X,R)*/
          TPKP_AB_To_XR (nvar, A, B, &X, &R, U);
          /* Calculation of (x,r)*/
          TPKP_C_To_c (nvar, X, R, &(s_x[i3D]), &r, U);
          /* Calculation of (y,z)*/
          TPKP_rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[i3D]), &(s_z[i3D]), U);
        }
  Set_Initial_Guess_for_u(cctkGH, n1*n2*n3, v.d0, s_x, s_y, s_z);
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++)
        {
          indx = TPKP_Index(ivar,i,j,k,1,n1,n2,n3);
          v.d0[indx]/=(-cos(Pih * (2 * i + 1) / n1)-1.0);
        }
  TPKP_Derivatives_AB3 (nvar, n1, n2, n3, v);
  if (do_initial_debug_output && CCTK_MyProc(cctkGH) == 0)
  {
    debug_file=fopen("initial.dat", "w");
    assert(debug_file);
    for (ivar = 0; ivar < nvar; ivar++)
      for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
        {
          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          Am1 = A -1.0;
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 0.0;
          indx = TPKP_Index(ivar,i,j,0,1,n1,n2,n3);
          U.d0[0] = Am1 * v.d0[indx];        /* U*/
          U.d1[0] = v.d0[indx] + Am1 * v.d1[indx];        /* U_A*/
          U.d2[0] = Am1 * v.d2[indx];        /* U_B*/
          U.d3[0] = Am1 * v.d3[indx];        /* U_3*/
          U.d11[0] = 2 * v.d1[indx] + Am1 * v.d11[indx];        /* U_AA*/
          U.d12[0] = v.d2[indx] + Am1 * v.d12[indx];        /* U_AB*/
          U.d13[0] = v.d3[indx] + Am1 * v.d13[indx];        /* U_AB*/
          U.d22[0] = Am1 * v.d22[indx];        /* U_BB*/
          U.d23[0] = Am1 * v.d23[indx];        /* U_B3*/
          U.d33[0] = Am1 * v.d33[indx];        /* U_33*/
        /* Calculation of (X,R)*/
        TPKP_AB_To_XR (nvar, A, B, &X, &R, U);
        /* Calculation of (x,r)*/
        TPKP_C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
        /* Calculation of (y,z)*/
        TPKP_rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[indx]), &(s_z[indx]), U);
        fprintf(debug_file,
                "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g "
                "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                (double)s_x[indx], (double)s_y[indx],
                (double)A,(double)B,
                (double)U.d0[0],
                (double)(-cos(Pih * (2 * i + 1) / n1)-1.0),
                (double)U.d1[0],
                (double)U.d2[0],
                (double)U.d3[0],
                (double)U.d11[0],
                (double)U.d22[0],
                (double)U.d33[0],
                (double)v.d0[indx],
                (double)v.d1[indx],
                (double)v.d2[indx],
                (double)v.d3[indx],
                (double)v.d11[indx],
                (double)v.d22[indx],
                (double)v.d33[indx]
                );
        }
    fprintf(debug_file, "\n\n");
    for (i=n2-10; i<n2; i++)
    {
      CCTK_REAL d;
      indx = TPKP_Index(0,0,i,0,1,n1,n2,n3);
      d = TPKP_PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,
              s_x[indx], 0.0, 0.0);
      fprintf(debug_file, "%.16g %.16g\n",
                (double)s_x[indx], (double)d);
    }
    fprintf(debug_file, "\n\n");
    for (i=n2-10; i<n2-1; i++)
    {
      CCTK_REAL d;
      int ip= TPKP_Index(0,0,i+1,0,1,n1,n2,n3);
      indx = TPKP_Index(0,0,i,0,1,n1,n2,n3);
      for (j=-10; j<10; j++)
      {
        d = TPKP_PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,
                s_x[indx]+(s_x[ip]-s_x[indx])*j/10,
                0.0, 0.0);
        fprintf(debug_file, "%.16g %.16g\n",
                (double)(s_x[indx]+(s_x[ip]-s_x[indx])*j/10), (double)d);
      }
    }
    fprintf(debug_file, "\n\n");
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
      {
        X = 2*(2.0*i/n1-1.0);
        R = 2*(1.0*j/n2);
        if (X*X+R*R > 1.0)
        {
          TPKP_C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
          TPKP_rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[indx]), &(s_z[indx]), U);
          *U.d0  = s_x[indx]*s_x[indx];
          *U.d1  = 2*s_x[indx];
          *U.d2  = 0.0;
          *U.d3 = 0.0;
          *U.d11 = 2.0;
          *U.d22 = 0.0;
          *U.d33 = *U.d12 = *U.d23 = *U.d13 = 0.0;
          TPKP_C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
          fprintf(debug_file,
                  "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                  (double)s_x[indx], (double)r, (double)X, (double)R, (double)U.d0[0],
                  (double)U.d1[0],
                  (double)U.d2[0],
                  (double)U.d3[0],
                  (double)U.d11[0],
                  (double)U.d22[0],
                  (double)U.d33[0]);
        }
      }
    fclose(debug_file);
  }
  free(s_z);
  free(s_y);
  free(s_x);
  TPKP_free_derivs (&U, nvar);
  /*exit(0);*/
}

/* -------------------------------------------------------------------*/
void TPKP_TwoPunctures(CCTK_ARGUMENTS);
void TPKP_TwoPunctures (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  * mp = par_m_plus;
  * mm = par_m_minus;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;

  int antisymmetric_lapse, averaged_lapse, pmn_lapse, brownsville_lapse;

  int const nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  int imin[3], imax[3];
  int const ntotal = n1 * n2 * n3 * nvar;
#if 0
  int percent10 = 0;
#endif
  static CCTK_REAL *F = NULL;
  static derivs u, v, cf_v;
  CCTK_REAL admMass;

  if (! F) {
    CCTK_REAL up, um;
    /* Solve only when called for the first time */
    F = dvector (0, ntotal - 1);
    TPKP_allocate_derivs (&u, ntotal);
    TPKP_allocate_derivs (&v, ntotal);
    TPKP_allocate_derivs (&cf_v, ntotal);

    if (use_sources) {
      CCTK_INFO ("Solving puncture equation for BH-NS/NS-NS system");
    } else {
      CCTK_INFO ("Solving puncture equation for BH-BH system");
    }
    CCTK_VInfo (CCTK_THORNSTRING, "b = %g", par_b);
    
    /* initialise to 0 */
    for (int j = 0; j < ntotal; j++)
    {
      cf_v.d0[j] = 0.0;
      cf_v.d1[j] = 0.0;
      cf_v.d2[j] = 0.0;
      cf_v.d3[j] = 0.0;
      cf_v.d11[j] = 0.0;
      cf_v.d12[j] = 0.0;
      cf_v.d13[j] = 0.0;
      cf_v.d22[j] = 0.0;
      cf_v.d23[j] = 0.0;
      cf_v.d33[j] = 0.0;
      v.d0[j] = 0.0;
      v.d1[j] = 0.0;
      v.d2[j] = 0.0;
      v.d3[j] = 0.0;
      v.d11[j] = 0.0;
      v.d12[j] = 0.0;
      v.d13[j] = 0.0;
      v.d22[j] = 0.0;
      v.d23[j] = 0.0;
      v.d33[j] = 0.0;
    }
    /* call for external initial guess */
    if (use_external_initial_guess)
    {
      set_initial_guess(cctkGH, v);
    }

    /* If bare masses are not given, iteratively solve for them given the 
       target ADM masses target_M_plus and target_M_minus and with initial 
       guesses given by par_m_plus and par_m_minus. */
    if(!(give_bare_mass)) {
      CCTK_REAL tmp, mp_adm_err, mm_adm_err;
      char valbuf[100];

      CCTK_REAL M_p = target_M_plus;
      CCTK_REAL M_m = target_M_minus;

      CCTK_VInfo (CCTK_THORNSTRING, "Attempting to find bare masses.");
      CCTK_VInfo (CCTK_THORNSTRING, "Target ADM masses: M_p=%g and M_m=%g",
                  (double) M_p, (double) M_m);
      CCTK_VInfo (CCTK_THORNSTRING, "ADM mass tolerance: %g", (double) adm_tol);

      /* Loop until both ADM masses are within adm_tol of their target */
      do {
        CCTK_VInfo (CCTK_THORNSTRING, "Bare masses: mp=%.15g, mm=%.15g",
                    (double)*mp, (double)*mm);
        TPKP_Newton (cctkGH, nvar, n1, n2, n3, v, Newton_tol, 1);

        TPKP_F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);

        up = TPKP_PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b, 0., 0.);
        um = TPKP_PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-par_b, 0., 0.);

        /* Calculate the ADM masses from the current bare mass guess */
        *mp_adm = (1 + up) * *mp + *mp * *mm / (4. * par_b);
        *mm_adm = (1 + um) * *mm + *mp * *mm / (4. * par_b);

        /* Check how far the current ADM masses are from the target */
        mp_adm_err = fabs(M_p-*mp_adm);
        mm_adm_err = fabs(M_m-*mm_adm);
        CCTK_VInfo (CCTK_THORNSTRING, "ADM mass error: M_p_err=%.15g, M_m_err=%.15g",
                    (double)mp_adm_err, (double)mm_adm_err);
        
        /* Invert the ADM mass equation and update the bare mass guess so that
           it gives the correct target ADM masses */
        tmp = -4*par_b*( 1 + um + up + um*up ) + 
                sqrt(16*par_b*M_m*(1 + um)*(1 + up) + 
                  pow(-M_m + M_p + 4*par_b*(1 + um)*(1 + up),2));
        *mp = (tmp + M_p - M_m)/(2.*(1 + up));
        *mm = (tmp - M_p + M_m)/(2.*(1 + um));
        
        /* Set the par_m_plus and par_m_minus parameters */
        sprintf (valbuf,"%.17g", (double) *mp);
        CCTK_ParameterSet ("par_m_plus", "TwoPunctures", valbuf);
        
        sprintf (valbuf,"%.17g", (double) *mm);
        CCTK_ParameterSet ("par_m_minus", "TwoPunctures", valbuf);
        
      } while ( (mp_adm_err > adm_tol) ||
                (mm_adm_err > adm_tol) );
                
      CCTK_VInfo (CCTK_THORNSTRING, "Found bare masses.");
    }

    TPKP_Newton (cctkGH, nvar, n1, n2, n3, v, Newton_tol, Newton_maxit);
 
    TPKP_F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);

    TPKP_SpecCoef(n1, n2, n3, 0, v.d0, cf_v.d0);

    CCTK_VInfo (CCTK_THORNSTRING,
		  "The two puncture masses are mp=%.17g and mm=%.17g",
                (double) *mp, (double) *mm);

    up = TPKP_PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b, 0., 0.);
    um = TPKP_PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-par_b, 0., 0.);

    /* Calculate the ADM masses from the current bare mass guess */
    *mp_adm = (1 + up) * *mp + *mp * *mm / (4. * par_b);
    *mm_adm = (1 + um) * *mm + *mp * *mm / (4. * par_b);

    CCTK_VInfo (CCTK_THORNSTRING, "Puncture 1 ADM mass is %g", (double) *mp_adm);
    CCTK_VInfo (CCTK_THORNSTRING, "Puncture 2 ADM mass is %g", (double) *mm_adm);

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    admMass = (*mp + *mm
               - 4*par_b*TPKP_PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, nvar, n1, n2, n3));
    CCTK_VInfo (CCTK_THORNSTRING, "The total ADM mass is %g", (double) admMass);
    *E = admMass;

    /*
      Run this in Mathematica (version 8 or later) with
        math -script <file>

      Needs["SymbolicC`"];
      co = Table["center_offset[" <> ToString[i] <> "]", {i, 0, 2}];
      r1 = co + {"par_b", 0, 0};
      r2 = co + {-"par_b", 0, 0};
      {p1, p2} = Table["par_P_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];
      {s1, s2} = Table["par_S_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];

      J = Cross[r1, p1] + Cross[r2, p2] + s1 + s2;

      JVar = Table["*J" <> ToString[i], {i, 1, 3}];
      Print[OutputForm@StringReplace[
        ToCCodeString@MapThread[CAssign[#1, CExpression[#2]] &, {JVar, J}], 
        "\"" -> ""]];
     */

    *J1 = -(center_offset[2]*par_P_minus[1]) + center_offset[1]*par_P_minus[2] - center_offset[2]*par_P_plus[1] + center_offset[1]*par_P_plus[2] + par_S_minus[0] + par_S_plus[0];
    *J2 = center_offset[2]*par_P_minus[0] - center_offset[0]*par_P_minus[2] + par_b*par_P_minus[2] + center_offset[2]*par_P_plus[0] - center_offset[0]*par_P_plus[2] - par_b*par_P_plus[2] + par_S_minus[1] + par_S_plus[1];
    *J3 = -(center_offset[1]*par_P_minus[0]) + center_offset[0]*par_P_minus[1] - par_b*par_P_minus[1] - center_offset[1]*par_P_plus[0] + center_offset[0]*par_P_plus[1] + par_b*par_P_plus[1] + par_S_minus[2] + par_S_plus[2];
  }

  if (CCTK_EQUALS(grid_setup_method, "Taylor expansion"))
  {
    gsm = GSM_Taylor_expansion;
  }
  else if (CCTK_EQUALS(grid_setup_method, "evaluation"))
  {
    gsm = GSM_evaluation;
  }
  else
  {
    CCTK_WARN (0, "internal error");
  }

  antisymmetric_lapse = CCTK_EQUALS(initial_lapse, "twopunctures-antisymmetric");
  averaged_lapse = CCTK_EQUALS(initial_lapse, "twopunctures-averaged");
	pmn_lapse = CCTK_EQUALS(initial_lapse, "psi^n");
  if (pmn_lapse)
		CCTK_VInfo(CCTK_THORNSTRING, "Setting initial lapse to psi^%f profile.",
               (double)initial_lapse_psi_exponent);
  brownsville_lapse = CCTK_EQUALS(initial_lapse, "brownsville");
  if (brownsville_lapse)
    CCTK_VInfo(CCTK_THORNSTRING, 
               "Setting initial lapse to a Brownsville-style profile "
               "with exp %f.",
               (double)initial_lapse_psi_exponent);

  CCTK_INFO ("Interpolating result");
  if (CCTK_EQUALS(metric_type, "static conformal")) {
    if (CCTK_EQUALS(conformal_storage, "factor")) {
      *conformal_state = 1;
    } else if (CCTK_EQUALS(conformal_storage, "factor+derivs")) {
      *conformal_state = 2;
    } else if (CCTK_EQUALS(conformal_storage, "factor+derivs+2nd derivs")) {
      *conformal_state = 3;
    }
  } else {
    *conformal_state = 0;
  }

  for (int d = 0; d < 3; ++ d)
  {
    /*
    imin[d] = 0           + (cctk_bbox[2*d  ] ? 0 : cctk_nghostzones[d]);
    imax[d] = cctk_lsh[d] - (cctk_bbox[2*d+1] ? 0 : cctk_nghostzones[d]);
    */
    imin[d] = 0;
    imax[d] = cctk_lsh[d];
  }

#pragma omp parallel for
  for (int k = imin[2]; k < imax[2]; ++k)
  {
    for (int j = imin[1]; j < imax[1]; ++j)
    {
      for (int i = imin[0]; i < imax[0]; ++i)
      {
#if 0
        /* We can't output this when running in parallel */
        if (percent10 != 10*(i+j*cctk_lsh[0]+k*cctk_lsh[0]*cctk_lsh[1]) /
                            (cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]))
        {
            percent10 = 10*(i+j*cctk_lsh[0]+k*cctk_lsh[0]*cctk_lsh[1]) /
                           (cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2]);
            CCTK_VInfo(CCTK_THORNSTRING, "%3d%% done", percent10*10);
        }
#endif

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

        CCTK_REAL xx, yy, zz;
        xx = x[ind] - center_offset[0];
        yy = y[ind] - center_offset[1];
        zz = z[ind] - center_offset[2];

        /* We implement swapping the x and z coordinates as follows.
           The bulk of the code that performs the actual calculations
           is unchanged.  This code looks only at local variables.
           Before the bulk --i.e., here-- we swap all x and z tensor
           components, and after the code --i.e., at the end of this
           main loop-- we swap everything back.  */
        if (swap_xz) {
          /* Swap the x and z coordinates */
          SWAP (xx, zz);
        }

//      Define radius and inverse radius here
        CCTK_REAL r_plus
          = sqrt(pow(xx - par_b, 2) + pow(yy, 2) + pow(zz, 2));
        CCTK_REAL r_minus
          = sqrt(pow(xx + par_b, 2) + pow(yy, 2) + pow(zz, 2));

        r_plus = pow (pow (r_plus, 4) + pow (TP_epsilon, 4), 0.25);
        r_minus = pow (pow (r_minus, 4) + pow (TP_epsilon, 4), 0.25);
        if (r_plus < TP_Tiny)
            r_plus = TP_Tiny;
        if (r_minus < TP_Tiny)
            r_minus = TP_Tiny;

#define EXTEND(M,r) \
          ( M * (3./8 * pow(r, 4) / pow(TP_Extend_Radius, 5) - \
                 5./4 * pow(r, 2) / pow(TP_Extend_Radius, 3) + \
                 15./8 / TP_Extend_Radius))


//      Define metric, curv, gauge functions and background conformal factor in quasi-isotropic coordinates
//      see Eqs. (70)-(73) in notes
        CCTK_REAL gij[3][3], detgK, Aij[3][3], alpK, betaK[3], psi04;
        TPKP_KQI_gijAijofxyz( xx, yy, zz, gij, &detgK, Aij, &alpK, betaK, &psi04 );


////DEBUG
//if( r_plus < 0.01 ){
//printf("=== debug in TwoPunctures.c -- after KQI_gij ===\n");
//printf("TwoPunctures.c: x, y, z, r_plus = %g \t %g \t %g \t %g\n", xx - par_b, yy, zz, r_plus);
//printf("TwoPunctures.c: psi04 = %g \n", psi04);
//printf("TwoPunctures.c: gij, detgK = %g\t %g\t %g\t %g\t %g\t %g\t %g\n", gij[0][0], gij[0][1], gij[0][2], gij[1][1], gij[1][2], gij[2][2], detgK);
//printf("TwoPunctures.c: Aij = %g\t %g\t %g\t %g\t %g\t %g\n", Aij[0][0], Aij[0][1], Aij[0][2], Aij[1][1], Aij[1][2], Aij[2][2]);
//printf("TwoPunctures.c: alph, beta = %g\t %g\t %g\t %g\n", alpK, betaK[0], betaK[1], betaK[2]);
//printf("===============================\n");
//fflush(stdout);
//}


//      Define direct Einstein-Maxwell electric/magnetic fields
        CCTK_REAL Ei[3], Bi[3], EM2;
        TPKP_Compute_EM_fields( xx, yy, zz, Ei, Bi, &EM2 );

//      Calculate correction function u
        CCTK_REAL U;
        switch (gsm)
        {
        case GSM_Taylor_expansion:
          U = TPKP_PunctTaylorExpandAtArbitPosition
            (0, nvar, n1, n2, n3, v, xx, yy, zz);
          break;
        case GSM_evaluation:
          U = TPKP_TPKP_PunctIntPolAtArbitPositionFast(0, nvar, n1, n2, n3, cf_v, xx, yy, zz);
          break;
        default:
          assert (0);
        }

//      Calculate conformal factor
        CCTK_REAL psi01, psi1;
        psi01 = pow( psi04, 0.25 );
        psi1  = psi01 * ( 1.0 + U );

////DEBUG
//if( r_plus < 0.01 ){
//printf("=== debug in TwoPunctures.c -- after calc conformal factor ===\n");
//printf("TwoPunctures.c: x, y, z, r_plus = %g \t %g \t %g \t %g\n", xx - par_b, yy, zz, r_plus);
//printf("TwoPunctures.c: U, psi04, psi0, psi = %g \t %g\t %g \t %g \n", puncture_u[ind], psi04, psi01, psi1);
//printf("===============================\n");
//fflush(stdout);
//}

//      just keep for backwards compatibility
        CCTK_REAL static_psi = 1;

        CCTK_REAL old_alp=1.0;
        if (multiply_old_lapse)
            old_alp = alp[ind];

        if ((*conformal_state > 0) || (pmn_lapse) || (brownsville_lapse)) {

          CCTK_REAL xp, yp, zp, rp, ir;
          CCTK_REAL s1, s3, s5;
          CCTK_REAL p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
          p = 1.0;
          px = py = pz = 0.0;
          pxx = pxy = pxz = 0.0;
          pyy = pyz = pzz = 0.0;

          /* first puncture */
          xp = xx - par_b;
          yp = yy;
          zp = zz;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TP_epsilon, 4), 0.25);
          if (rp < TP_Tiny)
              rp = TP_Tiny;
          ir = 1.0/rp;

          if (rp < TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* *mp *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          /* second puncture */
          xp = xx + par_b;
          yp = yy;
          zp = zz;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TP_epsilon, 4), 0.25);
          if (rp < TP_Tiny)
              rp = TP_Tiny;
          ir = 1.0/rp;

          if (rp < TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* *mm *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          if (*conformal_state >= 1) {
            static_psi = p;
            psi[ind] = static_psi;
          }
          if (*conformal_state >= 2) {
            psix[ind] = px / static_psi;
            psiy[ind] = py / static_psi;
            psiz[ind] = pz / static_psi;
          }
          if (*conformal_state >= 3) {
            psixx[ind] = pxx / static_psi;
            psixy[ind] = pxy / static_psi;
            psixz[ind] = pxz / static_psi;
            psiyy[ind] = pyy / static_psi;
            psiyz[ind] = pyz / static_psi;
            psizz[ind] = pzz / static_psi;
          }

          if (pmn_lapse)
            alp[ind] = pow(p, initial_lapse_psi_exponent);
          if (brownsville_lapse)
            alp[ind] = 2.0/(1.0+pow(p, initial_lapse_psi_exponent));

        } /* if conformal-state > 0 */

// Write to grid functions
//      /* correction function u */
        puncture_u[ind] = U;

//      /* metric components */
        gxx[ind] = gij[0][0] * pow(psi1, 4);
        gxy[ind] = gij[0][1] * pow(psi1, 4);
        gxz[ind] = gij[0][2] * pow(psi1, 4);
        gyy[ind] = gij[1][1] * pow(psi1, 4);
        gyz[ind] = gij[1][2] * pow(psi1, 4);
        gzz[ind] = gij[2][2] * pow(psi1, 4);

//      /* curvature components */
        kxx[ind] = Aij[0][0] / pow(psi1, 2);
        kxy[ind] = Aij[0][1] / pow(psi1, 2);
        kxz[ind] = Aij[0][2] / pow(psi1, 2);
        kyy[ind] = Aij[1][1] / pow(psi1, 2);
        kyz[ind] = Aij[1][2] / pow(psi1, 2);
        kzz[ind] = Aij[2][2] / pow(psi1, 2);

//      /* Lapse function */
        if (antisymmetric_lapse || averaged_lapse) {
          alp[ind] =
            ((1.0 -0.5* *mp /r_plus -0.5* *mm/r_minus)
            /(1.0 +0.5* *mp /r_plus +0.5* *mm/r_minus));

          if (r_plus < TP_Extend_Radius) {
            alp[ind] =
              ((1.0 -0.5*EXTEND(*mp, r_plus) -0.5* *mm/r_minus)
              /(1.0 +0.5*EXTEND(*mp, r_plus) +0.5* *mm/r_minus));
          }
          if (r_minus < TP_Extend_Radius) {
            alp[ind] =
              ((1.0 -0.5*EXTEND(*mm, r_minus) -0.5* *mp/r_plus)
              /(1.0 +0.5*EXTEND(*mp, r_minus) +0.5* *mp/r_plus));
          }
          
          if (averaged_lapse) {
            alp[ind] = 0.5 * (1.0 + alp[ind]);
          }
        }
        if (multiply_old_lapse)
          alp[ind] *= old_alp;

//      lapse for Kerr in quasi-isotropic coordinates
        if( CCTK_EQUALS(initial_lapse, "KQI_lapse") ){
            alp[ind] = alpK;
//            CCTK_VInfo(CCTK_THORNSTRING, 
//                "Setting initial lapse to that of Kerr in quasi-isotropic coordinates");
        }

//      pre-collapsed lapse
        if( CCTK_EQUALS(initial_lapse, "psi^n") ) {
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        }

//      /* Shift vector */
//      Set shift vector to its value in Kerr in quasi-isotropic coordinates
        if( CCTK_EQUALS(initial_shift, "KQI_shift") ){
          betax[ind] = betaK[0];
          betay[ind] = betaK[1];
          betaz[ind] = betaK[2];
        }


//      /* Direct Einstein-Maxwell variables */
//      Store contravariant E^i and B^i.  The psi^{-6}/sqrt(det(gtilde_KQI))
//      factor follows the original densitized-field convention of this thorn.
        Ex[ind]    = Ei[0] / pow(psi1, 6) / sqrt(fabs( detgK ));
        Ey[ind]    = Ei[1] / pow(psi1, 6) / sqrt(fabs( detgK ));
        Ez[ind]    = Ei[2] / pow(psi1, 6) / sqrt(fabs( detgK ));

        Bx[ind]    = Bi[0] / pow(psi1, 6) / sqrt(fabs( detgK ));
        By[ind]    = Bi[1] / pow(psi1, 6) / sqrt(fabs( detgK ));
        Bz[ind]    = Bi[2] / pow(psi1, 6) / sqrt(fabs( detgK ));

        Psi[ind]   = 0.0;
        Phi[ind]   = 0.0;

// DEBUG output
//if( psi04 < 0.1 ){
//if( r_plus < 0.01 ){
//printf("=== debug in TwoPunctures.c -- after writing to grid functions ===\n");
//printf("TwoPunctures.c: x, y, z, r_plus = %g \t %g \t %g \t %g\n", xx - par_b, yy, zz, r_plus);
//printf("TwoPunctures.c: U, psi04, psi0, psi = %g \t %g\t %g \t %g \n", puncture_u[ind], psi04, psi01, psi1);
//printf("TwoPunctures.c: gijC = %g \t %g \t %g \t %g \t %g \t %g \n", gxx[ind], gxy[ind], gxz[ind], gyy[ind], gyz[ind], gzz[ind]); 
//printf("TwoPunctures.c: KijC = %g \t %g \t %g \t %g \t %g \t %g \n", kxx[ind], kxy[ind], kxz[ind], kyy[ind], kyz[ind], kzz[ind]); 
//printf("TwoPunctures.c: alp = %g \t, beta = %g \t %g \t %g \n", alp[ind], betax[ind], betay[ind], betaz[ind]);
//printf("===============================\n");
//fflush(stdout);
//}

        if (swap_xz) {
          /* Swap the x and z components of all tensors */
          if (*conformal_state >= 2) {
            SWAP (psix[ind], psiz[ind]);
          }
          if (*conformal_state >= 3) {
            SWAP (psixx[ind], psizz[ind]);
            SWAP (psixy[ind], psiyz[ind]);
          }
          SWAP (gxx[ind], gzz[ind]);
          SWAP (gxy[ind], gyz[ind]);
          SWAP (kxx[ind], kzz[ind]);
          SWAP (kxy[ind], kyz[ind]);

          SWAP (betax[ind], betaz[ind]);

          SWAP (Bx[ind], Bz[ind]);
          SWAP (Ex[ind], Ez[ind]);

        } /* if swap_xz */

      } /* for i */
    }   /* for j */
  }     /* for k */

  if (use_sources && rescale_sources)
  {
    Rescale_Sources(cctkGH,
                    cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2],
                    x, y, z,
                    (*conformal_state > 0) ? psi : NULL,
                    gxx, gyy, gzz,
                    gxy, gxz, gyz);
  }

  if (0) {
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal - 1);
    TPKP_free_derivs (&u, ntotal);
    TPKP_free_derivs (&v, ntotal);
    TPKP_free_derivs (&cf_v, ntotal);
  }
}


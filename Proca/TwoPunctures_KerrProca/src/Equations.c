/* TwoPunctures_KerrEinsteinMaxwell:  File  "Equations.c"*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

// NOTE: We use x as the rotation axis!!!
//       That is ctheta = x / R
// FIXME: mass parameter as pointer

/* U.d0[ivar]   = U[ivar];  (ivar = 0..nvar-1) */
/* U.d1[ivar]   = U[ivar]_x;  */
/* U.d2[ivar]   = U[ivar]_y;  */
/* U.d3[ivar]   = U[ivar]_z;  */
/* U.d11[ivar]  = U[ivar]_xx; */
/* U.d12[ivar]  = U[ivar]_xy; */
/* U.d13[ivar]  = U[ivar]_xz; */
/* U.d22[ivar]  = U[ivar]_yy; */
/* U.d23[ivar]  = U[ivar]_yz; */
/* U.d33[ivar]  = U[ivar]_zz; */



/*-----------------------------------------------------------*/
/*--- Calculate A_{ij} A^{ij} -------------------------------*/
/*-----------------------------------------------------------*/
CCTK_REAL
TPKP_KQI_KKofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
  DECLARE_CCTK_PARAMETERS;

    //Structure of this function:
    //---------------------------------------------------------
    // 1) Define coordinates
    // 2) Define parameters
    // 3) Define metric functions 
    // 4) Define A2

    // Define coordinates
    CCTK_REAL xp = x - par_b;

    CCTK_REAL RR, RR2, RR6;
    CCTK_REAL rho, rho2;

    RR2 = xp * xp + y * y + z * z;
    RR2 = sqrt( pow( RR2, 2 ) + pow(TP_epsilon, 4));
    if( RR2 < pow(TP_Tiny,2))
        RR2 = pow(TP_Tiny,2);
    RR  = sqrt( RR2 );
    RR6 = pow( RR, 6 );

    rho2 = y*y + z*z;
//FIXME: no division by rho -> keep w/o eps_rho
//    rho2 = sqrt( pow( rho2, 2 ) + pow( TP_epsilon, 4 ));
//    if( rho2 < pow( TP_Tiny, 2) )
//        rho2 = pow( TP_Tiny, 2);
    rho  = sqrt( rho2 );

    CCTK_REAL ctheta, ctheta2, stheta, stheta2;
    ctheta  = xp / RR;
    ctheta2 = ctheta * ctheta;
    stheta  = rho / RR;
    stheta2 = stheta * stheta;


    // Define parameters
    // Note: in the original TP they use par_m_plus
    CCTK_REAL mass, mass2, spin, spin2, spin4, rBLp, rBLm;
    mass  = par_m_plus;
    mass2 = mass * mass;
    spin  = par_S_plus[0];
    spin2 = spin  * spin;
    spin4 = spin2 * spin2;
    rBLp  = mass + sqrt( mass2 - spin2 );
    rBLm  = mass - sqrt( mass2 - spin2 );


    // Define metric functions
    CCTK_REAL rBL, rBL2;
    rBL  = RR * ( 1.0 + 0.25 * rBLp / RR ) * ( 1.0 + 0.25 * rBLp / RR );
    rBL2 = rBL * rBL;

    CCTK_REAL Delt, Sigm, fctFF, fctFF2;
    Delt   = rBL2 + spin2 - 2 * mass * rBL;
    Sigm   = rBL2 + spin2 * ctheta2;
    fctFF  = ( rBL2 + spin2 ) * ( rBL2 + spin2 ) - Delt * spin2 * stheta2;
    fctFF2 = fctFF * fctFF;

    // Calculate AijAij
    CCTK_REAL AijAij, auxAijAij;
    auxAijAij = 2.0 * rBL2 * ( rBL2 + spin2 ) + Sigm * ( rBL2 - spin2 );

    AijAij = 2.0 * spin2 * mass2 * stheta2 / ( RR6 * fctFF2 )
                 * ( 4.0 * spin4 * rBL2 * Delt * ctheta2 * stheta2 + auxAijAij * auxAijAij );

    return AijAij;
}

/*----------------------------------------------------------------*/
/*--- Calculate gij, detg, Aij, \alpha, \beta^{i}, psi^{4}_{0} ---*/
/*----------------------------------------------------------------*/
void 
TPKP_KQI_gijAijofxyz( CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, 
                 CCTK_REAL gij[3][3], CCTK_REAL *detgK, CCTK_REAL Aij[3][3], 
                 CCTK_REAL *alpK, CCTK_REAL betaK[3], CCTK_REAL *psi04)
{
  DECLARE_CCTK_PARAMETERS;
  int ii, jj;

    //Structure of this function:
    //---------------------------------------------------------
    // 1) Define parameters
    // 2) Define coordinates
    // 3) Define metric functions 
    // 4) Initialize gauge functions
    // 5) Define metric, curvature etc in Cartesian coordinates (x,y,z) 

    // Define parameters
    // Note: in the original TP they use par_m_plus
    CCTK_REAL mass, mass2, spin, spin2, rBLp, rBLm;
    mass  = par_m_plus;
    mass2 = mass * mass;
    spin  = par_S_plus[0];
    spin2 = spin  * spin;
    rBLp  = mass + sqrt( mass2 - spin2 );
    rBLm  = mass - sqrt( mass2 - spin2 );
    /*----------------------------------*/

    // Define coordinates
    CCTK_REAL xp = x - par_b;

    CCTK_REAL RR, RR2, RR3;
    RR2 = xp * xp + y * y + z * z;
    RR2 = sqrt( pow( RR2, 2 ) + pow(TP_epsilon, 4));
    if( RR2 < pow(TP_Tiny,2))
        RR2 = pow(TP_Tiny,2);
    RR  = sqrt( RR2 );
    RR3 = RR2 * RR;

    CCTK_REAL rho, rho2, rho3;
    rho2 = y*y + z*z;
    rho2 = sqrt( pow( rho2, 2 ) + pow( TP_epsilon, 4 ));
    if( rho2 < pow( TP_Tiny, 2) )
        rho2 = pow( TP_Tiny, 2);
    rho  = sqrt( rho2 );
    rho3 = rho2 * rho;

    CCTK_REAL ctheta, ctheta2;
    ctheta  = xp / RR;
    ctheta2 = ctheta  * ctheta;

    CCTK_REAL stheta, stheta2;
    stheta  = rho / RR;
    stheta2 = stheta  * stheta;
    /*----------------------------------*/

    // Define metric functions
    CCTK_REAL rBL, rBL2;
    rBL  = RR * ( 1.0 + 0.25 * rBLp / RR ) * ( 1.0 + 0.25 * rBLp / RR );
    rBL2 = rBL * rBL;

    CCTK_REAL drBLdR;
    drBLdR = 1.0 - rBLp*rBLp / ( 16.0 * RR2 );

    CCTK_REAL Delt, Sigm, Sigm2, fctFF;
    Delt  = rBL2 + spin2 - 2 * mass * rBL;
    Sigm  = rBL2 + spin2 * ctheta2;
    Sigm2 = Sigm * Sigm;
    fctFF = ( rBL2 + spin2 ) * ( rBL2 + spin2 ) - Delt * spin2 * stheta2;

    CCTK_REAL fctGG, fctHH;
    fctGG = rBLm / ( RR2 * ( rBL - rBLm ) );
    fctHH = ( 2.0 * mass * rBL + Sigm ) / ( RR2 * Sigm2 );
    /*----------------------------------*/

    // Initialize gauge functions
    CCTK_REAL alpKQI;
    alpKQI = 0.25 * sqrt( Sigm * ( rBL - rBLm ) ) * ( 4.0 * RR - rBLp ) 
            / sqrt( RR * ( 2.0*mass*rBL * ( rBL2 + spin2 ) + Delt * Sigm ) );

    *alpK = alpKQI;

    CCTK_REAL bphi;
    bphi = 2.0 * spin * mass * rBL / fctFF;

    betaK[0] = 0;
    betaK[1] =   z * bphi;
    betaK[2] = - y * bphi;
    /*----------------------------------*/

    // Conformal factor and detg
    *psi04 = Sigm / RR2;

    *detgK = ( 1.0 + RR2 * fctGG ) * ( 1.0 + spin2 * rho2 * fctHH );

    // Metric
    gij[0][0] = 1.0 + xp * xp * fctGG;
    gij[0][1] =       xp * y  * fctGG;
    gij[0][2] =       xp * z  * fctGG;
    gij[1][1] = 1.0 + y  * y  * fctGG + spin2 * z * z * fctHH;
    gij[1][2] =       y  * z  * fctGG - spin2 * y * z * fctHH;
    gij[2][2] = 1.0 + z  * z  * fctGG + spin2 * y * y * fctHH;

    // Curvature
    CCTK_REAL auxAij, facAij, facAijX, facAijRho;
    auxAij    = 2.0 * rBL2 * ( rBL2 + spin2 ) + Sigm * ( rBL2 - spin2 );
    facAij    = alpKQI * spin * mass * stheta2 / ( RR3 * rho3 * Delt * sqrt( pow( Sigm, 3)) );
    facAijX   = 2.0 * rho * spin2 * rBL * Delt * ctheta * stheta + xp  * RR * drBLdR * auxAij;
    facAijRho = 2.0 * xp  * spin2 * rBL * Delt * ctheta * stheta - rho * RR * drBLdR * auxAij;

    // Aij[ii][jj]
    Aij[0][0] = 0;
    Aij[0][1] = - z * rho   * facAij * facAijX;
    Aij[0][2] =   y * rho   * facAij * facAijX;
    Aij[1][1] =  2.0 * y*z  * facAij * facAijRho;
    Aij[1][2] = (z*z - y*y) * facAij * facAijRho;
    Aij[2][2] = -2.0 * y*z  * facAij * facAijRho;
    /*---------------------------------------*/

//// DEBUG output
///printf("=== debug in Equations.c ===\n");
//printf("in Equations.c: x, y, z, RR = %g\t %g\t %g\t %g\n", xp, y, z, RR);
//printf("in Equations.c: gij = %g \t %g \t %g \t %g \t %g \t %g \n", (double) gij[0][0], gij[0][1], gij[0][2], gij[1][1], gij[1][2], gij[2][2] );
//printf("in Equations.c: Aij = %g \t %g \t %g \t %g \t %g \t %g \n", (double) Aij[0][0], Aij[0][1], Aij[0][2], Aij[1][1], Aij[1][2], Aij[2][2] );
//printf("in Equations.c: detgK, psi04, alpK = %g \t %g \t %g\n", detgK, psi04, alpK);
//printf("============================\n");
//fflush(stdout);

#if(0)
  if (isnan(Er)) {
    fprintf(stdout,
            "nan at x, y, z =  %.16g %.16g %.16g\n",
            xp, y, z ) ;
  }
#endif


}


/*-----------------------------------------------------------*/
/*--- Compute direct Einstein-Maxwell field quantities -------*/
/*-----------------------------------------------------------*/
/* conformal E^i in cart coordinates, plus A2 and EM2 */
//TODO: Recall that now we do not use the flat metric --> 
// adjust for different, NON-flat conformal metric!!!
void TPKP_Compute_EM_fields( CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
                             CCTK_REAL Ei[3], CCTK_REAL Bi[3], CCTK_REAL *EM2 )
{
  DECLARE_CCTK_PARAMETERS;

  /*
   * Direct Einstein-Maxwell initial data for electric and magnetic charges.
   * The fields returned here are conformal Cartesian vector components.
   * They are later rescaled by psi^{-6}/sqrt(det(gtilde_KQI)) before being
   * stored in ProcaBase::Ei and ProcaBase::Bi.
   *
   * The two punctures are placed at x=+par_b and x=-par_b.  par_q_plus/minus
   * are electric charges and par_p_plus/minus are magnetic charges.
   */
  CCTK_REAL dxp = x - par_b;
  CCTK_REAL dyp = y;
  CCTK_REAL dzp = z;
  CCTK_REAL dxm = x + par_b;
  CCTK_REAL dym = y;
  CCTK_REAL dzm = z;

  CCTK_REAL rp2 = dxp*dxp + dyp*dyp + dzp*dzp;
  CCTK_REAL rm2 = dxm*dxm + dym*dym + dzm*dzm;

  rp2 = pow(pow(rp2, 4) + pow(TP_epsilon, 4), 0.25);
  rm2 = pow(pow(rm2, 4) + pow(TP_epsilon, 4), 0.25);

  CCTK_REAL rp = sqrt(rp2);
  CCTK_REAL rm = sqrt(rm2);
  if (rp < TP_Tiny) rp = TP_Tiny;
  if (rm < TP_Tiny) rm = TP_Tiny;

  CCTK_REAL inv_rp3 = 1.0 / (rp*rp*rp);
  CCTK_REAL inv_rm3 = 1.0 / (rm*rm*rm);

  Ei[0] = par_q_plus  * dxp * inv_rp3 + par_q_minus  * dxm * inv_rm3;
  Ei[1] = par_q_plus  * dyp * inv_rp3 + par_q_minus  * dym * inv_rm3;
  Ei[2] = par_q_plus  * dzp * inv_rp3 + par_q_minus  * dzm * inv_rm3;

  Bi[0] = par_p_plus  * dxp * inv_rp3 + par_p_minus  * dxm * inv_rm3;
  Bi[1] = par_p_plus  * dyp * inv_rp3 + par_p_minus  * dym * inv_rm3;
  Bi[2] = par_p_plus  * dzp * inv_rp3 + par_p_minus  * dzm * inv_rm3;

  *EM2 = Ei[0]*Ei[0] + Ei[1]*Ei[1] + Ei[2]*Ei[2]
       + Bi[0]*Bi[0] + Bi[1]*Bi[1] + Bi[2]*Bi[2];
}


/*-----------------------------------------------------------*/
/*--- Nonlinear Equations -----------------------------------*/
/*-----------------------------------------------------------*/
void
TPKP_NonTPKP_LinEquations (CCTK_REAL rho_adm,
                 CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
                 CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
                 CCTK_REAL y, CCTK_REAL z, derivs U, CCTK_REAL *values)
{
  DECLARE_CCTK_PARAMETERS;
  int ii, jj;

    //Structure of this function:
    //---------------------------------------------------------
    // 1) Define coordinates
    // 2) Define parameters
    // 3) Define metric functions 
    // 4) Define auxiliary and derivatives of metric functions
    // 5) Define determinant and inverse metric
    // 6) Define direct Einstein-Maxwell fields
    // 7) Define terms of Hamiltonian
    // 8) Define Hamiltonian
    //---------------------------------------------------------

    // 1) Define coordinates
    CCTK_REAL xp = x - par_b;
    CCTK_REAL RR, RR2, RR3;
    RR2 = xp * xp + y * y + z * z;
    RR2 = sqrt( pow( RR2, 2 ) + pow( TP_epsilon, 4 ));
    if( RR2 < pow(TP_Tiny,2)) 
        RR2 = pow(TP_Tiny,2);
    RR  = sqrt( RR2 );
    RR3 = RR * RR2;

    CCTK_REAL rho, rho2;
    rho2 = y*y + z*z;
    rho2 = sqrt( pow( rho2, 2 ) + pow( TP_epsilon, 4 ));
    if( rho2 < pow( TP_Tiny, 2) )
        rho2 = pow( TP_Tiny, 2);
    rho  = sqrt( rho2 );

    CCTK_REAL ctheta, ctheta2;
    ctheta  = xp / RR;
    ctheta2 = ctheta  * ctheta;

    CCTK_REAL stheta, stheta2;
    stheta  = rho / RR;
    stheta2 = stheta  * stheta;

    // 2) Define parameters
    // Note: in the original TP they use par_m_plus
    CCTK_REAL mass, mass2, spin, spin2, rBLp, rBLp2, rBLm;
    mass  = par_m_plus;
    mass2 = mass * mass;
    spin  = par_S_plus[0];
    spin2 = spin  * spin;
    rBLp  = mass + sqrt( mass2 - spin2 );
    rBLp2 = rBLp * rBLp;
    rBLm  = mass - sqrt( mass2 - spin2 );

    // 3) Define metric functions 
    CCTK_REAL rBL, rBL2;
    rBL  = RR * ( 1.0 + 0.25 * rBLp / RR ) * ( 1.0 + 0.25 * rBLp / RR );
    rBL2 = rBL * rBL;

    CCTK_REAL Sigm, Sigm2, Sigm3;
    Sigm  = rBL2 + spin2 * ctheta2;
    Sigm2 = Sigm * Sigm;
    Sigm3 = Sigm * Sigm2;

    CCTK_REAL psi04, psi08;
    psi04 = Sigm / RR2;
    psi08 = psi04 * psi04;

    CCTK_REAL fctGG, fctHH;
    fctGG = rBLm / ( RR2 * ( rBL - rBLm ) );
    fctHH = ( 2.0 * mass * rBL + Sigm ) / ( RR2 * Sigm2 );

    // 4) Define auxiliary and derivatives of metric functions
    CCTK_REAL auxGG, auxGG2, auxGGx, auxHH;
    auxGG  = 1.0 + RR2 * fctGG;
    auxGG2 = auxGG * auxGG;
    auxGGx = 1.0 + xp*xp * fctGG;
    auxHH  = 1.0 + spin2 * rho2 * fctHH;

    CCTK_REAL drBLdR;
    drBLdR = 1.0 - 0.0625 * rBLp2 / RR2;

    CCTK_REAL dGdR, dHdR, dHdT;
    dGdR = - fctGG * ( 2.0 / RR + drBLdR / ( rBL - rBLm ) );
    dHdR = - 2.0 * ( 2.0*mass*rBL + Sigm + RR * drBLdR * ( rBL - mass + 4.0*mass*rBL2 / Sigm ) ) / ( RR3 * Sigm2 );
    dHdT = 2.0 * spin2 * ctheta * stheta * ( 4.0*mass*rBL + Sigm ) / ( RR2 * Sigm3 );

    // 5) Define determinant and inverse metric
    CCTK_REAL detg, guu[3][3];

    detg = auxGG * auxHH;

    guu[0][0] = ( 1.0 + rho2            * fctGG ) / auxGG;
    guu[0][1] =   - xp * y              * fctGG   / auxGG;
    guu[0][2] =   - xp * z              * fctGG   / auxGG;
    guu[1][1] = ( 1.0 + ( xp*xp + z*z ) * fctGG + spin2 * y*y * fctHH * auxGGx ) / ( auxGG * auxHH );
    guu[1][2] = ( - y  * z              * fctGG + spin2 * y*z * fctHH * auxGGx ) / ( auxGG * auxHH );
    guu[2][2] = ( 1.0 + ( xp*xp + y*y ) * fctGG + spin2 * z*z * fctHH * auxGGx ) / ( auxGG * auxHH );
    //Symmetries
    guu[1][0] = guu[0][1];
    guu[2][0] = guu[0][2];
    guu[2][1] = guu[1][2];

    // 6) Define direct Einstein-Maxwell fields
    CCTK_REAL Ei[3], Bi[3], EM2;
    TPKP_Compute_EM_fields(x, y, z, Ei, Bi, &EM2);

    // 7) Define terms of Hamiltonian
    // 7.1) Hamiltonian term1 = \gamma^{ij} D_{i} D_{j} u
    CCTK_REAL HamTerm1;
    CCTK_REAL HamTerm1a, HamTerm1b;
    // HamTerm1a = \gamma^{ij} \p_{i}\p_{j} u
    HamTerm1a =         guu[0][0] * U.d11[0]
                + 2.0 * guu[0][1] * U.d12[0]
                + 2.0 * guu[0][2] * U.d13[0]
                +       guu[1][1] * U.d22[0]
                + 2.0 * guu[1][2] * U.d23[0]
                +       guu[2][2] * U.d33[0];

    // HamTerm1b = \Gamma^{i} \p_{i} u
    CCTK_REAL CoefGam;
    CoefGam = 0.5 / ( RR2 * auxGG2 * auxHH );

    CCTK_REAL termGam1, termGam2, termGam3, termGam4, termGam5, termGam6x;
    termGam6x = 0.5 * spin2 * ( 4.0*xp * fctHH + rho * dHdT  ) / auxHH;
    termGam1  = RR3 * auxHH * dGdR;
    termGam2  = - spin2 * RR * rho2 * auxGG * dHdR;
    termGam3  = - spin2 * xp * rho * auxGG2 * dHdT;
    termGam4  = - 4.0 * spin2 * RR2 * fctHH * auxGG2;
    termGam5  = 2.0 * RR2 * fctGG * ( auxHH + 3.0 * auxGG * auxHH - auxGG );

    CCTK_REAL tGamx, tGamy, tGamz;
    tGamx = xp * CoefGam * ( termGam1 + termGam2 + termGam3 + termGam4 + termGam5 ) + termGam6x;
    tGamy = y  * CoefGam * ( termGam1 + termGam2 + termGam3 + termGam4 + termGam5 );
    tGamz = z  * CoefGam * ( termGam1 + termGam2 + termGam3 + termGam4 + termGam5 );

    HamTerm1b = tGamx * U.d1[0] + tGamy * U.d2[0] + tGamz * U.d3[0];

    HamTerm1 = HamTerm1a - HamTerm1b;

    // 7.2) Hamiltonian term2 = 1/(2psi04) \gamma^{ij} D_{i}u D_{j}\psi04
    CCTK_REAL HamTerm2;
    CCTK_REAL d1psi04[3], CoefDpsi04;

    CoefDpsi04 = 0.03125 / ( RR3 * Sigm );

    d1psi04[0] = CoefDpsi04 * xp * ( pow( 4.0*RR + rBLp, 3 ) - 64.0 * RR * Sigm + 32.0 * spin2 * RR );
    d1psi04[1] = CoefDpsi04 * y  * ( pow( 4.0*RR + rBLp, 3 ) - 64.0 * RR * Sigm );
    d1psi04[2] = CoefDpsi04 * z  * ( pow( 4.0*RR + rBLp, 3 ) - 64.0 * RR * Sigm );

    HamTerm2 =   guu[0][0] *   U.d1[0] * d1psi04[0]
               + guu[0][1] * ( U.d1[0] * d1psi04[1] + U.d2[0] * d1psi04[0] )
               + guu[0][2] * ( U.d1[0] * d1psi04[2] + U.d3[0] * d1psi04[0] )
               + guu[1][1] *   U.d2[0] * d1psi04[1] 
               + guu[1][2] * ( U.d2[0] * d1psi04[2] + U.d3[0] * d1psi04[1] )
               + guu[2][2] *   U.d3[0] * d1psi04[2];

    // 7.3) Hamiltonian term3 \sim A^{ij} A_{ij}
    CCTK_REAL HamTerm3;
    HamTerm3 = 0.125 * ( 1.0 - pow( 1.0 + U.d0[0], 8 ) ) * TPKP_KQI_KKofxyz (x, y, z)
                / ( psi08 * pow( 1.0 + U.d0[0], 7 ) );

    // 7.4) Hamiltonian electromagnetic term, proportional to E^2 + B^2
    CCTK_REAL HamTerm4;
    HamTerm4 = 0.25 * EM2 / ( detg * psi04 * pow( 1.0 + U.d0[0], 3 ) );

    // 8) Define Hamiltonian
    values[0] = HamTerm1 + HamTerm2 + HamTerm3 + HamTerm4;






//DEBUG output
//if( RR < 0.001 ){
//printf("=== debug in Equations.c -- after calculating the Hamiltonian ===\n");
//printf("Equations.c: x, y, z, RR, rho = %g \t %g \t %g \t %g \t %g \n", xp, y, z, RR, rho);
//printf("Equations.c: rBL, rBL2, Sigm, Sigm2, Sigm3, psi04, psi08 = %g\t %g\t %g\t %g\t %g\t %g\t %g\n", rBL, rBL2, Sigm, Sigm2, Sigm3, psi04, psi08);
//printf("Equations.c: fctGG, fctHH = %g\t %g\n", fctGG, fctHH);
//printf("Equations.c: auxGG, auxGG2, auxGGx, auxHH = %g\t %g\t %g\t %g\n", auxGG, auxGG2, auxGGx, auxHH);
//printf("Equations.c: drBLdR, dGdR, dHdR, dHdT = %g\t %g\t %g\t %g\n", drBLdR, dGdR, dHdR, dHdT);
//printf("Equations.c: detg, gxx, gxy, gxz = %g\t %g\t %g\t %g\n", detg, guu[0][0], guu[0][1], guu[0][2]);
//printf("Equations.c: gyx, gyy, gyz, gxz, gyz, gzz = %g\t %g\t %g\t %g\t %g\t %g\n", guu[1][0], guu[1][1], guu[1][2], guu[2][0], guu[2][1], guu[2][2]);
//printf("Equations.c: CoefDpsi04, d1psi04_x, d1psi04_y, d1psi04_z = %g\t %g\t %g\t %g\n", CoefDpsi04, d1psi04[0], d1psi04[1], d1psi04[2]);
//printf("Equations.c: U.d1, U.d2, U.d3 = %g\t %g\t %g\n", U.d1[0], U.d2[0], U.d3[0]);
//printf("Equations.c: HamTerm2 = %g\n", HamTerm2); 
//printf("Equations.c: U.d11, U.d12, U.d13, U.d22, U.d23, U.d33 = %g\t %g\t %g\t %g\t %g\t %g\n", U.d11[0], U.d12[0], U.d13[0], U.d22[0], U.d23[0], U.d33[0]);
//printf("Equations.c: CoefGam, termGam6x = %g\t %g\n", CoefGam, termGam6x);
//printf("Equations.c: termGam1, termGam2, termGam3, termGam4, termGam5 = %g\t %g\t %g\t %g\t %g\n", termGam1, termGam2, termGam3, termGam4, termGam5);
//printf("Equations.c: tGamx, tGamy, tGamz = %g\t %g\t %g\n", tGamx, tGamy, tGamz); 
//printf("Equations.c: U.d1, U.d2, U.d3 = %g\t %g\t %g\n", U.d1[0], U.d2[0], U.d3[0]);
//printf("Equations.c: HamTerm1a, HamTerm1b, HamTerm1 = %g\t %g\t %g\n", HamTerm1a, HamTerm1b, HamTerm1);
//printf("Equations.c: Terms of Hamiltonian: %g \t %g \t %g \t %g \n", HamTerm1, HamTerm2, HamTerm3, HamTerm4);
//printf("Equations.c: Hamiltonian = %g\n", values[0]);
//printf("===============================\n");
//fflush(stdout);
//}


// debug output
#if(0)
  if (fabs(values[0] ) > 1.e-9 ) {
    fprintf(stdout, "r_plus, r_minus =   %.16g  %.16g \n",
            r_plus, r_minus);
    fprintf(stdout,
            "x, y, z =  %.16g %.16g %.16g \nres = %.16g \n\n",
            x, y, z, values[0] ) ;
  }
#endif
}

/*-----------------------------------------------------------*/
/*--- Linear Equations --------------------------------------*/
/*-----------------------------------------------------------*/
void
TPKP_LinEquations (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R,
	      CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi,
	      CCTK_REAL y, CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values)
{

  DECLARE_CCTK_PARAMETERS;
  int ii, jj;

    //Structure of this function:
    //---------------------------------------------------------
    // 1) Define coordinates
    // 2) Define parameters
    // 3) Define metric functions 
    // 4) Define auxiliary and derivatives of metric functions
    // 5) Define determinant and inverse metric
    // 6) Define direct Einstein-Maxwell fields
    // 7) Define terms of linearized Hamiltonian
    // 8) Define linearized Hamiltonian
    //---------------------------------------------------------

    // 1) Define coordinates
    CCTK_REAL xp = x - par_b;
    CCTK_REAL RR, RR2, RR3;
    RR2 = xp * xp + y * y + z * z;
    RR2 = sqrt( pow( RR2, 2 ) + pow( TP_epsilon, 4 ));
    if( RR2 < pow(TP_Tiny,2)) 
        RR2 = pow(TP_Tiny,2);
    RR  = sqrt( RR2 );
    RR3 = RR * RR2;

    CCTK_REAL rho, rho2;
    rho2 = y*y + z*z;
    rho2 = sqrt( pow( rho2, 2 ) + pow( TP_epsilon, 4 ));
    if( rho2 < pow( TP_Tiny, 2) )
        rho2 = pow( TP_Tiny, 2);
    rho  = sqrt( rho2 );

    CCTK_REAL ctheta, ctheta2;
    ctheta  = xp / RR;
    ctheta2 = ctheta  * ctheta;

    CCTK_REAL stheta, stheta2;
    stheta  = rho / RR;
    stheta2 = stheta  * stheta;

    // 2) Define parameters
    // Note: in the original TP they use par_m_plus
    CCTK_REAL mass, mass2, spin, spin2, rBLp, rBLp2, rBLm;
    mass  = par_m_plus;
    mass2 = mass * mass;
    spin  = par_S_plus[0];
    spin2 = spin  * spin;
    rBLp  = mass + sqrt( mass2 - spin2 );
    rBLp2 = rBLp * rBLp;
    rBLm  = mass - sqrt( mass2 - spin2 );

    // 3) Define metric functions 
    CCTK_REAL rBL, rBL2;
    rBL  = RR * ( 1.0 + 0.25 * rBLp / RR ) * ( 1.0 + 0.25 * rBLp / RR );
    rBL2 = rBL * rBL;

    CCTK_REAL Sigm, Sigm2, Sigm3;
    Sigm  = rBL2 + spin2 * ctheta2;
    Sigm2 = Sigm * Sigm;
    Sigm3 = Sigm * Sigm2;

    CCTK_REAL psi04, psi08;
    psi04 = Sigm / RR2;
    psi08 = psi04 * psi04;

    CCTK_REAL fctGG, fctHH;
    fctGG = rBLm / ( RR2 * ( rBL - rBLm ) );
    fctHH = ( 2.0 * mass * rBL + Sigm ) / ( RR2 * Sigm2 );

    // 4) Define auxiliary and derivatives of metric functions
    CCTK_REAL auxGG, auxGG2, auxGGx, auxHH;
    auxGG  = 1.0 + RR2 * fctGG;
    auxGG2 = auxGG * auxGG;
    auxGGx = 1.0 + xp*xp * fctGG;
    auxHH  = 1.0 + spin2 * rho2 * fctHH;

    CCTK_REAL drBLdR;
    drBLdR = 1.0 - 0.0625 * rBLp2 / RR2;

    CCTK_REAL dGdR, dHdR, dHdT;
    dGdR = - fctGG * ( 2.0 / RR + drBLdR / ( rBL - rBLm ) );
    dHdR = - 2.0 * ( 2.0*mass*rBL + Sigm + RR * drBLdR * ( rBL - mass + 4.0*mass*rBL2 / Sigm ) ) / ( RR3 * Sigm2 );
    dHdT = 2.0 * spin2 * ctheta * stheta * ( 4.0*mass*rBL + Sigm ) / ( RR2 * Sigm3 );

    // 5) Define determinant and inverse metric
    CCTK_REAL detg, guu[3][3];

    detg = auxGG * auxHH;

    guu[0][0] = ( 1.0 + rho2            * fctGG ) / auxGG;
    guu[0][1] =   - xp * y              * fctGG   / auxGG;
    guu[0][2] =   - xp * z              * fctGG   / auxGG;
    guu[1][1] = ( 1.0 + ( xp*xp + z*z ) * fctGG + spin2 * y*y * fctHH * auxGGx ) / ( auxGG * auxHH );
    guu[1][2] = ( - y  * z              * fctGG + spin2 * y*z * fctHH * auxGGx ) / ( auxGG * auxHH );
    guu[2][2] = ( 1.0 + ( xp*xp + y*y ) * fctGG + spin2 * z*z * fctHH * auxGGx ) / ( auxGG * auxHH );
    //Symmetries
    guu[1][0] = guu[0][1];
    guu[2][0] = guu[0][2];
    guu[2][1] = guu[1][2];

    // 6) Define direct Einstein-Maxwell fields
    CCTK_REAL Ei[3], Bi[3], EM2;
    TPKP_Compute_EM_fields(x, y, z, Ei, Bi, &EM2);

    // 7) Define terms of Hamiltonian
    // 7.1) Hamiltonian term1 = \gamma^{ij} D_{i} D_{j} u
    CCTK_REAL HamTerm1;
    CCTK_REAL HamTerm1a, HamTerm1b;
    // HamTerm1a = \gamma^{ij} \p_{i}\p_{j} u
    HamTerm1a =         guu[0][0] * dU.d11[0]
                + 2.0 * guu[0][1] * dU.d12[0]
                + 2.0 * guu[0][2] * dU.d13[0]
                +       guu[1][1] * dU.d22[0]
                + 2.0 * guu[1][2] * dU.d23[0]
                +       guu[2][2] * dU.d33[0];

    // HamTerm1b = \Gamma^{i} \p_{i} u
    CCTK_REAL CoefGam;
    CoefGam = 0.5 / ( RR2 * auxGG2 * auxHH );

    CCTK_REAL termGam1, termGam2, termGam3, termGam4, termGam5, termGam6x;
    termGam6x = 0.5 * spin2 * ( 4.0*xp * fctHH + rho * dHdT  ) / auxHH;
    termGam1  = RR3 * auxHH * dGdR;
    termGam2  = - spin2 * RR * rho2 * auxGG * dHdR;
    termGam3  = - spin2 * xp * rho * auxGG2 * dHdT;
    termGam4  = - 4.0 * spin2 * RR2 * fctHH * auxGG2;
    termGam5  = 2.0 * RR2 * fctGG * ( auxHH + 3.0 * auxGG * auxHH - auxGG );

    CCTK_REAL tGamx, tGamy, tGamz;
    tGamx = xp * CoefGam * ( termGam1 + termGam2 + termGam3 + termGam4 + termGam5 ) + termGam6x;
    tGamy = y  * CoefGam * ( termGam1 + termGam2 + termGam3 + termGam4 + termGam5 );
    tGamz = z  * CoefGam * ( termGam1 + termGam2 + termGam3 + termGam4 + termGam5 );

    HamTerm1b = tGamx * dU.d1[0] + tGamy * dU.d2[0] + tGamz * dU.d3[0];

    HamTerm1 = HamTerm1a - HamTerm1b;

    // 7.2) Hamiltonian term2 = 1/(2psi04) \gamma^{ij} D_{i}u D_{j}\psi04
    CCTK_REAL HamTerm2;
    CCTK_REAL d1psi04[3], CoefDpsi04;

    CoefDpsi04 = 0.03125 / ( RR3 * Sigm );

    d1psi04[0] = CoefDpsi04 * xp * ( pow( 4.0*RR + rBLp, 3 ) - 64.0 * RR * Sigm + 32.0 * spin2 * RR );
    d1psi04[1] = CoefDpsi04 * y  * ( pow( 4.0*RR + rBLp, 3 ) - 64.0 * RR * Sigm );
    d1psi04[2] = CoefDpsi04 * z  * ( pow( 4.0*RR + rBLp, 3 ) - 64.0 * RR * Sigm );

    HamTerm2 =   guu[0][0] *   dU.d1[0] * d1psi04[0]
               + guu[0][1] * ( dU.d1[0] * d1psi04[1] + dU.d2[0] * d1psi04[0] )
               + guu[0][2] * ( dU.d1[0] * d1psi04[2] + dU.d3[0] * d1psi04[0] )
               + guu[1][1] *   dU.d2[0] * d1psi04[1] 
               + guu[1][2] * ( dU.d2[0] * d1psi04[2] + dU.d3[0] * d1psi04[1] )
               + guu[2][2] *   dU.d3[0] * d1psi04[2];

    // 7.3) Hamiltonian term3 \sim A^{ij} A_{ij}
    CCTK_REAL HamTerm3;
    HamTerm3 = 0.125 * ( 7.0 + pow( 1.0 + U.d0[0], 8 ) ) * TPKP_KQI_KKofxyz (x, y, z)
                / ( psi08 * pow( 1.0 + U.d0[0], 8 ) );

    // 7.4) Hamiltonian electromagnetic term, proportional to E^2 + B^2
    CCTK_REAL HamTerm4;
    HamTerm4 = 0.75 * EM2 / ( detg * psi04 * pow( 1.0 + U.d0[0], 4 ) );


    // 8) Define linearized Hamiltonian
    values[0] = HamTerm1 + HamTerm2 - dU.d0[0] * ( HamTerm3 + HamTerm4 );

}

/*-----------------------------------------------------------*/

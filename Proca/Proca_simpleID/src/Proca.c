
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"



/* Swap two variables */
static inline
void swap (CCTK_REAL * restrict const a, CCTK_REAL * restrict const b)
{
  CCTK_REAL const t = *a; *a=*b; *b=t;
}
#undef SWAP
#define SWAP(a,b) (swap(&(a),&(b)))

/* -------------------------------------------------------------------*/
void Proca (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int imin[3], imax[3];

  for (int d = 0; d < 3; ++ d)
  {
    /*
    imin[d] = 0           + (cctk_bbox[2*d  ] ? 0 : cctk_nghostzones[d]);
    imax[d] = cctk_lsh[d] - (cctk_bbox[2*d+1] ? 0 : cctk_nghostzones[d]);
    */
    imin[d] = 0;
    imax[d] = cctk_lsh[d];
  }


  for (int i = imin[0]; i < imax[0]; ++i) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int k = imin[2]; k < imax[2]; ++k) {

        const int ind = CCTK_GFINDEX3D (cctkGH, i, j, k);

        CCTK_REAL x1, y1, z1;
        x1 = x[ind];
        y1 = y[ind];
        z1 = z[ind];

        /* We implement swapping the x and z coordinates as follows.
           The bulk of the code that performs the actual calculations
           is unchanged.  This code looks only at local variables.
           Before the bulk --i.e., here-- we swap all x and z tensor
           components, and after the code --i.e., at the end of this
           main loop-- we swap everything back.  */
        if (swap_xz) {
          /* Swap the x and z coordinates */
          SWAP (x1, z1);
        }

        CCTK_REAL r_plus
          = sqrt(pow(x1 - par_b, 2) + pow(y1, 2) + pow(z1, 2));
        CCTK_REAL r_minus
          = sqrt(pow(x1 + par_b, 2) + pow(y1, 2) + pow(z1, 2));

        CCTK_REAL psi1 = sqrt( pow(1
                                   + 0.5 * par_m_plus / r_plus
                                   + 0.5 * par_m_minus/ r_minus , 2)
                               - 0.25 * pow( par_q_plus/r_plus
                                             + par_q_minus/r_minus, 2) ) ;

        gxx[ind] = pow (psi1, 4);
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = pow (psi1, 4);
        gyz[ind] = 0;
        gzz[ind] = pow (psi1, 4);

        kxx[ind] = 0;
        kxy[ind] = 0;
        kxz[ind] = 0;
        kyy[ind] = 0;
        kyz[ind] = 0;
        kzz[ind] = 0;


        /* EMG terms */

        Zeta[ind]  = 0;

        Ax[ind]    = 0;
        Ay[ind]    = 0;
        Az[ind]    = 0;

        Aphi[ind]  = 0;

        Ex[ind]    = (  par_q_plus * (x1-par_b)/(r_plus*r_plus*r_plus)
                      + par_q_minus* (x1+par_b)/(r_minus*r_minus*r_minus) )
                        / pow(psi1, 6) ;

        Ey[ind]    = (  par_q_plus * y1/(r_plus*r_plus*r_plus)
                      + par_q_minus* y1/(r_minus*r_minus*r_minus) )
                        / pow(psi1, 6) ;

        Ez[ind]    = (  par_q_plus * z1/(r_plus*r_plus*r_plus)
                      + par_q_minus* z1/(r_minus*r_minus*r_minus) )
                        / pow(psi1, 6) ;

        // lapse
        if ( CCTK_EQUALS(initial_lapse, "psi^n") ) {
          alp[ind] = pow(psi1, initial_lapse_psi_exponent);
        }

        if (swap_xz) {
          /* Swap the x and z components of all tensors */
          SWAP (gxx[ind], gzz[ind]);
          SWAP (gxy[ind], gyz[ind]);
          SWAP (kxx[ind], kzz[ind]);
          SWAP (kxy[ind], kyz[ind]);

          SWAP (Ax[ind], Az[ind]);
          SWAP (Ex[ind], Ez[ind]);

        } /* if swap_xz */


      } /* for k */
    }   /* for j */
  }     /* for i */

}


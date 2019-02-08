/* 
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include <stdint.h>
#include <stdlib.h>

#include "coefficients.h"

// Fortran dimensions: f1(NI,NJ,NK)
// NI      : length of a line
// NINJ    : length of a plane
// zlinear : zero if cubic interpolation along z, non zero if linear interpolation along z
// f1      : address of lower left corner of the 4 x 4 x 4 or 4 x 4 x 2 interpolation box
// interpolation is done along z, then along y, then along x
// values are interpolated using the cubic/linear polynomial coefficients
// pxyz(8) and pxyz(11) both = 0.0 when linear interpolation along z
// pxyz(9) = 1.0 - dz, pxyz(10) = dz are expected in the linear case
// in the z linear case, planes 0 and 1 are the same as are planes 2 and 3
// static inline void Tricublin_zyxf_beta(float *d, float *f1, double *px, double *py, double *pz, int NI, int NINJ){
static inline void Tricublin_zyxf1_inline(float *d, float *f1, double *pxyz, int NI, int NINJ, int zlinear){
  int ni = NI;
  int ninj = NINJ;
  int ninjl;    // ninj (cubic along z) or 0 (linear along z)
  int ni2, ni3;
  double *px = pxyz;
  double *py = pxyz+4;
  double *pz = pxyz+8;

  double va4[4], vb4[4], vc4[4], vd4[4], dst[4];
  int i, ninj2, ninj3;
  
  ni2 = ni + ni;
  ni3 = ni2 + ni;
  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be set to 0
  if(zlinear) ninjl = 0;

  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  // field 1
  for (i=0 ; i<4 ; i++){
    va4[i] = f1[i    ]*pz[0] + f1[i    +ninjl]*pz[1] +  f1[i    +ninj2]*pz[2] + f1[i    +ninj3]*pz[3];
    vb4[i] = f1[i+ni ]*pz[0] + f1[i+ni +ninjl]*pz[1] +  f1[i+ni +ninj2]*pz[2] + f1[i+ni +ninj3]*pz[3];
    vc4[i] = f1[i+ni2]*pz[0] + f1[i+ni2+ninjl]*pz[1] +  f1[i+ni2+ninj2]*pz[2] + f1[i+ni2+ninj3]*pz[3];
    vd4[i] = f1[i+ni3]*pz[0] + f1[i+ni3+ninjl]*pz[1] +  f1[i+ni3+ninj2]*pz[2] + f1[i+ni3+ninj3]*pz[3];
    dst[i] = va4[i]*py[0] + vb4[i]*py[1] + vc4[i]*py[2] + vd4[i]*py[3];
  }
  d[0] = dst[0]*px[0] + dst[1]*px[1] + dst[2]*px[2] + dst[3]*px[3];
}

// process n points
// for each point 1 value from ixyz, 3 values from pxyz are used (1 element)
// it is ASSUMED that along X and Y interpolation will always be CUBIC 
// ( 2 <= px < "ni"-1 ) ( 2 <= py < "nj"-1 )
// pz < 2 and pz >= nk - 1 will induce linear interpolation or extrapolation
void Tricublin_zyx1_n(float *d, float *f1, pxpypz *pxyz,  ztab *lv, int n){
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
/*printf*/("%12.7f %12.7f %12.7f\n",pxyz->px, pxyz->py, pxyz->pz);
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    Tricublin_zyxf1_inline(d, f1 + ixyz, cxyz, lv->ni, lv->nij, zlinear);         // interpolate
    d++;         // next result
    pxyz += 1;   // next set of positions
  }
}


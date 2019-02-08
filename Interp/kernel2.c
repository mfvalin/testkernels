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
// Fortran dimensions: d(3) , f(3,NI,NJ,NK)
// NI   : size of a 1D line
// NINJ : size of a 2D plane
// f : address of lower left corner of the 4 x 4 x 4 or 4 x 4 x 2 box
// interpolation is done along z, then along y, then along x
// 3 values are interpolated using the same cubic/linear polynomial coefficients
// pz(2) = 1.0 - dz, pz(3) = dz are expected in the linear along z case (0.0 <= dz <= 1.0)
// in the z linear case, planes 0 and 1 are the same as are planes 2 and 3
// note for future expansion:
// should both interpolations have to be done, planes 1 and 2 can be used for the z linear case,
// and planes 0, 1, 2, 3 for the z cubic case
// in that case another mechanism will have to be used to signal the z linear case
static inline void Tricublin_zyxf3_inline(float *d, float *f, double *pxyz, int NI, int NINJ, int zlinear){
// static inline void Tricublin_zyx3f_beta(float *d, float *f, double *px, double *py, double *pz, int NI, int NINJ){
  int ni = 3*NI;
  int ninj = 3*NINJ;
  double *px = pxyz;
  double *py = pxyz+4;
  double *pz = pxyz+8;
  int ninjl;    // ninj (cubic along z) or 0 (linear along z)
  float *s = f;
  double dst[13];
//   int64_t *L = (int64_t *)d;
  int ni2, ni3;

  double va4[12], vb4[12], vc4[12], vd4[12];
  int ninj2, ninj3, i;

  ni2 = ni + ni;
  ni3 = ni2 + ni;
  dst[12] = 0;
  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be set to 0
  if(zlinear) ninjl = 0;

  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  for (i=0 ; i<12 ; i++){
    va4[i] = s[i    ]*pz[0] + s[i    +ninjl]*pz[1] +  s[i    +ninj2]*pz[2] + s[i    +ninj3]*pz[3];
    vb4[i] = s[i+ni ]*pz[0] + s[i+ni +ninjl]*pz[1] +  s[i+ni +ninj2]*pz[2] + s[i+ni +ninj3]*pz[3];
    vc4[i] = s[i+ni2]*pz[0] + s[i+ni2+ninjl]*pz[1] +  s[i+ni2+ninj2]*pz[2] + s[i+ni2+ninj3]*pz[3];
    vd4[i] = s[i+ni3]*pz[0] + s[i+ni3+ninjl]*pz[1] +  s[i+ni3+ninj2]*pz[2] + s[i+ni3+ninj3]*pz[3];
    dst[i] = va4[i]*py[0] + vb4[i]*py[1] + vc4[i]*py[2] + vd4[i]*py[3];
  }
  d[0] = dst[0]*px[0] + dst[3]*px[1] + dst[6]*px[2] + dst[ 9]*px[3];
  d[1] = dst[1]*px[0] + dst[4]*px[1] + dst[7]*px[2] + dst[10]*px[3];
  d[2] = dst[2]*px[0] + dst[5]*px[1] + dst[8]*px[2] + dst[11]*px[3];
}

// process n triplets
// for each triplet 1 value from ixyz, 3 values from pxyz are used (1 element)
// it is ASSUMED that along X and Y interpolation will always be CUBIC 
// ( 2 <= px < "ni"-1 ) ( 2 <= py < "nj"-1 )
// pz < 2 and pz >= nk - 1 will induce linear interpolation or extrapolation
void Tricublin_zyx3_n(float *d, float *f3, pxpypz *pxyz,  ztab *lv, int n){
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    Tricublin_zyxf3_inline(d, f3 + ixyz*3, cxyz, lv->ni, lv->nij, zlinear);         // interpolate
    d+=3;         // next result
    pxyz += 1;   // next set of positions
  }
}

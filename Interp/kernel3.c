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

static inline void Tricublin_zyx_mm_d_inline(float *d, float *lin, float *min, float *max, float *f1, double *cxyz, int NI, int NINJ, int zlinear){
  int ni = NI;
  int ninj = NINJ;
  int ninjl;    // ninj (cubic along z) or 0 (linear along z)
  int ni2, ni3;
  double *px, *py, *pz;
  double va4[4], vb4[4], vc4[4], vd4[4], dst[4], dsl[4];
  int i, ninj2, ninj3;
  double ma, mi;

  px = cxyz;
  py = cxyz + 4;
  pz = cxyz + 8;
  ni2 = ni + ni;
  ni3 = ni2 + ni;
  ninjl = ninj ; // assuming cubic case. in the linear case, ninjl will be set to 0
  if(zlinear) ninjl = 0;

  ninj2 = ninj + ninjl;
  ninj3 = ninj2 + ninjl;
  // field 1
  // TODO: add min max code
  for (i=0 ; i<4 ; i++){   // tricubic or bicubic/linear interpolation
    va4[i] = f1[i    ]*pz[0] + f1[i    +ninjl]*pz[1] +  f1[i    +ninj2]*pz[2] + f1[i    +ninj3]*pz[3];
    vb4[i] = f1[i+ni ]*pz[0] + f1[i+ni +ninjl]*pz[1] +  f1[i+ni +ninj2]*pz[2] + f1[i+ni +ninj3]*pz[3];
    vc4[i] = f1[i+ni2]*pz[0] + f1[i+ni2+ninjl]*pz[1] +  f1[i+ni2+ninj2]*pz[2] + f1[i+ni2+ninj3]*pz[3];
    vd4[i] = f1[i+ni3]*pz[0] + f1[i+ni3+ninjl]*pz[1] +  f1[i+ni3+ninj2]*pz[2] + f1[i+ni3+ninj3]*pz[3];
    dst[i] = va4[i]*py[0] + vb4[i]*py[1] + vc4[i]*py[2] + vd4[i]*py[3];
  }
  d[0] = dst[0]*px[0] + dst[1]*px[1] + dst[2]*px[2] + dst[3]*px[3];
  for (i=1 ; i<3 ; i++){   // linear interpolation
    vb4[i] = f1[i+ni +ninjl]*cxyz[18] + f1[i+ni +ninj2]*cxyz[19];
    vc4[i] = f1[i+ni2+ninjl]*cxyz[18] + f1[i+ni2+ninj2]*cxyz[19];
    dsl[i] = vb4[i]*cxyz[16] + vc4[i]*cxyz[17];
  }
  lin[0] = dsl[1]*cxyz[13] + dsl[2]*cxyz[14];
  ma = f1[1 + ni + ninjl] ; mi = ma;   // point [1,1,1] of 2x2x2 inner box
#define MAX(a,b) ((a) > (b)) ? (a) : (b)
#define MIN(a,b) ((a) < (b)) ? (a) : (b)
  for (i=1 ; i<3 ; i++){              // min max of 2x2x2 inner box
    ma = MAX(ma , f1[i + ni  + ninjl]); ma = MAX(ma , f1[i + ni  + ninj2]);
    mi = MIN(mi , f1[i + ni  + ninjl]); mi = MIN(mi , f1[i + ni  + ninj2]);

    ma = MAX(ma , f1[i + ni2 + ninjl]); ma = MAX(ma , f1[i + ni2 + ninj2]);
    mi = MIN(mi , f1[i + ni2 + ninjl]); mi = MIN(mi , f1[i + ni2 + ninj2]);
  }
  *max = ma ; *min = mi;
}

void Tricublin_mono_zyx_n(float *d, float *l, float *mi, float *ma, float *f, pxpypz *pxyz,  ztab *lv, int n){
  int i;
  double cxyz[24];   // interpolation coefficients 4 for each dimension (x, y, z)
  int ixyz;          // unilinear index into array f1 (collapsed dimensions)
  int zlinear;       // non zero if linear interpolation
                     // all above computed in Vcoef_pxyz4, used in Tricublin_zyxf1
  while(n--){
    zlinear = Vcoef_pxyz4_inline(cxyz, &ixyz, pxyz->px, pxyz->py, pxyz->pz, lv);  // compute coefficients
    Tricublin_zyx_mm_d_inline(d, l, mi, ma, f + ixyz, cxyz, lv->ni, lv->nij, zlinear);         // interpolate
    d++;         // next result
    l++;
    mi++;
    ma++;
    pxyz += 1;   // next set of positions
  }
}

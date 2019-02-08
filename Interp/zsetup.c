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
#include <stdlib.h>
#include <stdint.h>
#include "coefficients.h"

// allocate lookup table set and
// return pointer to filled table set
// targets are expected to be positive, and monotonically increasing
ztab *Vsearch_setup(double *targets, int nk, int ni, int nj){
  int i;
  ztab *lv;
  double pad;

  lv = malloc(sizeof(ztab));
  if(lv == NULL) return NULL ;

  lv->z = malloc(nk * sizeof(double));
  if(NULL == lv->z){      // malloc failed
    free(lv);             // deallocate lv
    return NULL;
  }
  for(i=0 ; i<nk  ; i++) { lv->z[i] = targets[i] ; }   // z coordinate table
  
  lv->ocz = malloc(4 * nk * sizeof(double));
  if(NULL == lv->ocz){    // malloc failed
    free(lv->z);          // deallocate z
    free(lv);             // deallocate lv
    return NULL;
  }
  denominators( &(lv->ocz[0]) , targets[0], targets[1], targets[2], targets[3]);                 // level 0 coeffs are normally not used
  for(i=1 ; i<nk - 2   ; i++) {
    denominators( &(lv->ocz[4*i]) , lv->z[i-1], lv->z[i  ], lv->z[i+1], lv->z[i+2]);
  }
  denominators( &(lv->ocz[4*(nk-2)]) , targets[nk-4], targets[nk-3], targets[nk-2], targets[nk-1]);  // level nk-2 coeffs are normally not used
  denominators( &(lv->ocz[4*(nk-1)]) , targets[nk-4], targets[nk-3], targets[nk-2], targets[nk-1]);  // level nk-1 coeffs are normally not used

  lv->ni = ni;           // nb of points along x
  lv->nj = nj;           // nb of points along y
  lv->nk = nk;           // nb of points along z
  lv->nij = ni*nj;       // ni * nj
  lv->offi = 0;
  lv->offj = 0;
  return lv;             // return pointer to filled table
}
ztab *Vsearch_setup_plus(double *targets, int nk, int ni, int nj, int offseti, int offsetj){
  ztab *lv;
  lv = Vsearch_setup(targets, nk, ni, nj);
  if(lv != NULL) {
    lv->offi = offseti;
    lv->offj = offsetj;
  }
  return lv;             // return pointer to filled table
}

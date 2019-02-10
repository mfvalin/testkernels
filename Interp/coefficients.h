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
#if defined(USE_FLOAT)
#define double float
#endif

typedef struct{
  float px;    // position along x in index space
  float py;    // position along y in index space
  float pz;    // position along z in index space
//   float z;     // absolute position along z 
} pxpypz;

typedef struct{
  double *z;       // table for levels (nk doubles)
  double *ocz;     // pointer to inverse of coefficient denominators [4*nk doubles]
  uint32_t ni;     // distance between rows
  uint32_t nj;     // number of rows
  uint32_t nk;     // number of levels (max 255)
  uint32_t nij;    // ni*nj, distance between levels
  uint32_t offi;   // offset to add to i position (e.g. global to local grid remapping)
  uint32_t offj;   // offset to add to j position (e.g. global to local grid remapping)
}ztab;

static double cp167 =  1.0/6.0;
static double cm167 = -1.0/6.0;
static double cm5 = -0.5;
static double cp5 = 0.5;
static double one = 1.0;
static double two = 2.0;

// triple product used for Lagrange cubic polynomials coefficients in the non constant case
#define TRIPRD(x,a,b,c) ((x-a)*(x-b)*(x-c))

// inverse denominators r from positions a b c d
static inline void denominators(double *r, double a, double b, double c, double d){
  r[0] = 1.0 / TRIPRD(a,b,c,d);
  r[1] = 1.0 / TRIPRD(b,a,c,d);
  r[2] = 1.0 / TRIPRD(c,a,b,d);
  r[3] = 1.0 / TRIPRD(d,a,b,c);
}

static inline int Vcoef_pxyz4_inline(double *cxyz, int *offset, float px8, float py8, float pz8, ztab *lv){
  int ix, iy, iz, ijk, zlinear;
  double pxy[2], *base, *pos;
  double zza, zzb, zzc, zzd, zzab, zzcd, dz, px, py, pz;
  int i, j;

    px = px8 ;            // fractional index positions along x, y, z (float to double)
    py = py8 ;
    pz = pz8 ;

    ix = px ;
    px = px - ix;                   // px is now deltax (fractional part of px)
    ix = ix + lv->offi;             // global to local grid remapping
    ijk = ix - 2;                   // x displacement (elements), ix assumed to always be >1 and < ni-1

    cxyz[12] = 0.0;
    cxyz[13] = 1.0 - px;             // linear interpolation coefficients along x
    cxyz[14] = px;
    cxyz[15] = 0.0;

    iy = py ;
    py = py - iy;                   // py is now deltay (fractional part of py)
    iy = iy + lv->offj;             // global to local grid remapping
    ijk = ijk + (iy - 2) * lv->ni;  // add y displacement (rows), ix assumed to always be >1 and < nj-1

    cxyz[16] = 1.0 - py;             // linear interpolation coefficients along y
    cxyz[17] = py;

    iz = pz ; 
    if(iz<1) iz = 1; 
    if(iz>lv->nk-1) iz = lv->nk-1;  // iz < 1 or iz > nk-1 will result in linear extrapolation
    dz = pz - iz;                   // dz is now "fractional" part of pz  (may be <0 or >1 if extrapolating)
    ijk = ijk + (iz -1) * lv->nij;  // add z displacement (2D planes)
    cxyz[18] = 1.0 - dz;             // linear interpolation coefficients along z
    cxyz[19] = dz;

    iz--;                           // iz needs to be in "origin 0" (C index from Fortran index)
    zlinear = (iz - 1) | (lv->nk - 3 - iz); 
    zlinear >>= 31;                  // nonzero only if iz < 1 or iz > nk -3 (top and bottom intervals)
    if(! zlinear) ijk = ijk - lv->nij;  // not the linear case, go down one 2D plane to get lower left corner of 4x4x4 cube
    *offset = ijk;

    // now we can compute the coefficients along z using iz and dz
    if(zlinear){
      cxyz[ 8] = 0.0;                    // coefficients for linear interpolation along z
      cxyz[ 9] = 1.0 - dz;
      cxyz[10] = dz;
      cxyz[11] = 0.0;
    }else{
      base  = &(lv->ocz[4*iz]);  // precomputed inverses of denominators
      pos   = &(lv->z[iz]);
      pz  = dz * pos[1] + (1.0 - dz) * pos[0];   // pz is now an absolute position
      zza = pz - pos[-1] ; zzb = pz - pos[0] ; zzc = pz - pos[1] ; zzd = pz - pos[2] ; 
      zzab = zza * zzb ; zzcd = zzc * zzd;
      cxyz[ 8] = zzb * zzcd * base[0];   //   cxyz[16] = TRIPRD(pz,pos[1],pos[2],pos[3]) * base[0];
      cxyz[ 9] = zza * zzcd * base[1];   //   cxyz[17] = TRIPRD(pz,pos[0],pos[2],pos[3]) * base[1];
      cxyz[10] = zzd * zzab * base[2];   //   cxyz[18] = TRIPRD(pz,pos[0],pos[1],pos[3]) * base[2];
      cxyz[11] = zzc * zzab * base[3];   //   cxyz[19] = TRIPRD(pz,pos[0],pos[1],pos[2]) * base[3];
    }

//     cxyz[ 0] = (px*px  - px)   * (-(cp167*px) + cp133);
//     cxyz[ 1] = (cp5*px - one)  * (px*px       - one);
//     cxyz[ 2] = (px*px  + px)   * (-(cp5*px)   + one);
//     cxyz[ 3] = (px*px  + px)   * (cp167*px    - cp167);
// 
//     cxyz[ 8] = (py*py  - py)   * (-(cp167*py) + cp133);
//     cxyz[ 9] = (cp5*py - one)  * (py*py       - one);
//     cxyz[10] = (py*py  + py)   * (-(cp5*py)   + one);
//     cxyz[11] = (py*py  + py)   * (cp167*py    - cp167);

    cxyz[ 0] = cm167*px*(px-one)*(px-two);        // coefficients for cubic interpolation along x
    cxyz[ 1] = cp5*(px+one)*(px-one)*(px-two);
    cxyz[ 2] = cm5*px*(px+one)*(px-two);
    cxyz[ 3] = cp167*px*(px+one)*(px-one);

    cxyz[ 4] = cm167*py*(py-one)*(py-two);        // coefficients for cubic interpolation along y
    cxyz[ 5] = cp5*(py+one)*(py-one)*(py-two);
    cxyz[ 6] = cm5*py*(py+one)*(py-two);
    cxyz[ 7] = cp167*py*(py+one)*(py-one);

  return zlinear;  // linear / cubic flag
}

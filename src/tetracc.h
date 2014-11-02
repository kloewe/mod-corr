/*----------------------------------------------------------------------
  File    : tetracc.h
  Contents: compute pairwise tetrachoric correlation coefficients
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#ifndef __TETRACC__
#define __TETRACC__

#ifdef _MSC_VER
#define uint32_t   unsigned __int32
#define uint64_t   unsigned __int64
#define inline     __inline
#else                           /* MSC still does not support C99 */
#include <stdint.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif
#ifdef __SSSE3__
#include <tmmintrin.h>
#endif
#ifdef __SSE4_1__
#include <smmintrin.h>
#endif
#ifdef __POPCNT__
#include <popcntintrin.h>
#endif

#include <assert.h>

/*----------------------------------------------------------------------
  Data Type Definition / Recursion Handling
----------------------------------------------------------------------*/
#ifdef REAL                     /* if REAL is defined */
#  undef  _TCC_PASS             /* ensure _TCC_PASS is undefined */
#  define _TCC_PASS 0           /* define macro for single pass */
#  ifndef SUFFIX                /* function names get no suffix */
#  define SUFFIX                /* (only single set of functions) */
#  endif
#elif !defined _TCC_PASS        /* if in first pass of two */
#  undef  _TCC_PASS             /* ensure _TCC_PASS is undefined */
#  define _TCC_PASS 1           /* define macro for first pass */
#  define REAL      float       /* first pass: single precision */
#  define SUFFIX    Flt         /* function name suffix is 'Flt' */
#else                           /* if in second pass of two */
#  undef  _TCC_PASS             /* ensure _TCC_PASS is undefined */
#  define _TCC_PASS 2           /* define macro for second pass */
#  define REAL      double      /* second pass: double precision */
#  define SUFFIX    Dbl         /* function name suffix is 'Dbl' */
#endif

#define float  1                /* to check the definition of REAL */
#define double 2
#if REAL==float                 /* if single precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  0       /* clear indicator for double */
#else                           /* if double precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  1       /* set   indicator for double */
#endif
#undef float                    /* delete definitions */
#undef double                   /* used for type checking */

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#ifndef TCC_AUTO
#define TCC_AUTO        0       /* use automatic choice */
#define TCC_LUT16       1       /* use a lookup table for 16 bit */
#define TCC_SSE2        2       /* use SSE2  computations */
#define TCC_SSSE3       3       /* use SSSE3 computations */
#define TCC_POP32       4       /* use _mm_popcnt_u32() */
#define TCC_POP64       5       /* use _mm_popcnt_u64() */
#define TCC_M128I       6       /* use 128 bit ints & popcnt64 */

#define TCC_VARIANT  0x0f       /* mask for variant */
#define TCC_TILED    0x10       /* flag for tiled    version */
#define TCC_THREAD   0x20       /* flag for threaded version */
#endif

#ifndef SFXNAME                 /* macros to generate function names */
#define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#define SFXNAME_2(n,s)  n##s    /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
extern int SFXNAME(tetracc)  (REAL *data, REAL *res, int N, int T);
extern int SFXNAME(tetraccx) (REAL *data, REAL *res, int N, int T,
                              int var, ...);

extern REAL* SFXNAME(make_cmap) (int T);
extern void init_popcnt (void);

/*----------------------------------------------------------------------
  Variables
----------------------------------------------------------------------*/
extern unsigned char popcnt[1 << 16];

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#if   _TCC_PASS <= 0
#define tetracc(d,r,N,T)     tetraccx(d,r,N,T,0)
#elif _TCC_PASS <= 1
#define tetraccFlt(d,r,N,T)  tetraccxFlt(d,r,N,T,0)
#elif _TCC_PASS <= 2
#define tetraccDbl(d,r,N,T)  tetraccxDbl(d,r,N,T,0)
#endif

/*----------------------------------------------------------------------
  Inline Functions
----------------------------------------------------------------------*/

inline int pcand_lut16 (uint32_t *a, uint32_t *b, int n)
{                               /* --- pop. count of conjunction */
  int k, s;                     /* loop variable, population count */

  assert(a && b && (n > 0));    /* check the function arguments */
  for (k = s = 0; k < n; k++) { /* traverse the two series and */
    uint32_t x = a[k] & b[k];   /* compute elementwise conjunction */
    s += (int)popcnt[x & 0xffff] +(int)popcnt[x >> 16];
  }                             /* count set bits in the conjunction */
  return s;                     /* return the population count */
}  /* pcand_lut16() */

/*--------------------------------------------------------------------*/
#if defined __POPCNT__ && defined __SSE4_1__
inline int pcand_m128i (uint32_t *a, uint32_t *b, int n)
{                               /* --- pop. count of conjunction */
  int k;                        /* loop variable */
  int s = 0;                    /* sum of population counts */

  for (k = 0; k < n; k += 4) {  /* traverse the binarized data */
    __m128i x = _mm_and_si128(_mm_load_si128((__m128i*)(a+k)),
                              _mm_load_si128((__m128i*)(b+k)));
    s += (int)_mm_popcnt_u64((uint64_t)_mm_extract_epi64(x, 0))
      +  (int)_mm_popcnt_u64((uint64_t)_mm_extract_epi64(x, 1));
  }                             /* count set bits in the conjunction */
  return s;                     /* and return this number */
}  /* pcand_m128i() */

#endif  /* #if defined __POPCNT__ && defined __SSE4_1__ */
/*----------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------*/
#if _TCC_PASS == 1              /* if if in first of two passes */
#undef REAL
#undef SUFFIX
#include "tetracc.h"            /* process header recursively */
#elif _TCC_PASS == 2
#undef REAL
#endif

#undef SUFFIX
#undef SFXNAME
#undef SFXNAME_1
#undef SFXNAME_2
#undef REAL_IS_DOUBLE

#undef  _TCC_PASS
#define __TETRACC__
#endif

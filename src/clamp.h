/*----------------------------------------------------------------------
  File    : clamp.h
  Contents: clamp values to the interval [min,max]
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#ifndef CLAMP_H

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "assert.h"
#ifdef __SSE2__
#include "emmintrin.h"
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions for Data Type / Recursion Handling
----------------------------------------------------------------------*/
#ifdef REAL                     /* if REAL is defined, */
#  undef  _CLAMP_PASS           /* ensure _CLAMP_PASS is undefined */
#  define _CLAMP_PASS 0         /* define macro for single pass */
#  ifndef SUFFIX                /* function names get no suffix */
#  define SUFFIX                /* (only single set of functions) */
#  endif
#elif !defined _CLAMP_PASS      /* if in first pass of two */
#  undef  _CLAMP_PASS           /* ensure _CLAMP_PASS is undefined */
#  define _CLAMP_PASS 1         /* define macro for first pass */
#  define REAL      float       /* first pass: single precision */
#  define SUFFIX    _flt        /* function name suffix is '_flt' */
#else                           /* if in second pass of two */
#  undef  _CLAMP_PASS           /* ensure _CLAMP_PASS is undefined */
#  define _CLAMP_PASS 2         /* define macro for second pass */
#  define REAL      double      /* second pass: double precision */
#  define SUFFIX    _dbl        /* function name suffix is '_dbl' */
#endif

/*--------------------------------------------------------------------*/
#define float  1                /* to check the definition of REAL */
#define double 2

#if REAL==float                 /* if single precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  0       /* clear indicator for double */
#elif REAL==double              /* if double precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  1       /* set   indicator for double */
#else
#error "REAL must be either 'float' or 'double'"
#endif

#undef float                    /* delete definitions */
#undef double                   /* used for type checking */
/*--------------------------------------------------------------------*/

#ifndef SFXNAME                 /* macros to generate function names */
#define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#define SFXNAME_2(n,s)  n##s    /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

/*----------------------------------------------------------------------
  Inline Functions
----------------------------------------------------------------------*/

inline REAL SFXNAME(clamp1) (REAL r, REAL minval, REAL maxval)
{
  if (r < minval) return minval;
  if (r > maxval) return maxval;
  return r;
}  /* clamp1() */

/*--------------------------------------------------------------------*/

inline REAL SFXNAME(clamp2) (REAL r, REAL minval, REAL maxval)
{
  if (r < 0) return (r > minval) ? r : minval;
  else       return (r < maxval) ? r : maxval;
}  /* clamp2() */

/*--------------------------------------------------------------------*/

inline REAL SFXNAME(clamp3) (REAL r, REAL minval, REAL maxval)
{
  return (fabs(r) <= maxval) ? r : (r < 0) ? minval : maxval;
}  /* clamp3() */

/*--------------------------------------------------------------------*/

#ifdef __SSE2__
inline REAL SFXNAME(clamp4) (REAL val, REAL minval, REAL maxval)
{ /* min(max(r,-1),1) */
  #if REAL_IS_DOUBLE            /* data is double precision */
  _mm_store_sd(&val, _mm_min_sd(
            _mm_max_sd(_mm_set_sd(val), _mm_set_sd(minval)),
            _mm_set_sd(maxval)));
  #else                         /* data is single precision */
  _mm_store_ss(&val, _mm_min_ss(
            _mm_max_ss(_mm_set_ss(val), _mm_set_ss(minval)),
            _mm_set_ss(maxval)));
  #endif
  return val;
}  /* clamp4() */
#endif

/*----------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------*/
#if   _CLAMP_PASS == 1            /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "clamp.h"                /* process header recursively */
#elif _CLAMP_PASS == 2
#undef REAL
#endif

#undef SUFFIX
#undef SFXNAME
#undef SFXNAME_1
#undef SFXNAME_2
#undef REAL_IS_DOUBLE

#undef  _CLAMP_PASS

#define CLAMP_H
#endif  /* #ifndef CLAMP_H */

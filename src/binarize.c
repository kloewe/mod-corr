/*----------------------------------------------------------------------
  File    : binarize.c
  Contents: functions to binarize real-valued series
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifdef _MSC_VER
#ifndef uint32_t
#define uint16_t   unsigned __int16
#define uint32_t   unsigned __int32
#define uint64_t   unsigned __int64
#endif                          /* MSC still does not support C99 */
#else
#include <stdint.h>
#endif

#include "binarize.h"

/*----------------------------------------------------------------------
  Data Type Definition / Recursion Handling
----------------------------------------------------------------------*/
#ifndef REAL                    /* if to compile as float and double */
#  define REAL    float         /* first pass: single precision */
#  define SUFFIX  _flt          /* function name suffix is '_flt' */
#else                           /* if REAL defined or in second pass, */
#  ifndef __BINARIZE_2__        /* prevent (further) recursion by */
#  define __BINARIZE_2__        /* defining indicator for second pass */
#  endif
#  ifndef SUFFIX                /* unless in recursion (second pass), */
#  define SUFFIX                /* function names get no suffix */
#  endif
#endif

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#ifndef SFXNAME                 /* macros to generate function names */
#  define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#  define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#  define SFXNAME_2(n,s)  n##s  /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/

static REAL SFXNAME(mean) (REAL *data, size_t n)
{                               /* --- compute mean of a data array */
  size_t i;                     /* loop variable */
  REAL   sum = 0;               /* sum of data values */

  assert(data && (n > 0));      /* check the function arguments */
  for (i = 0; i < n; i++)       /* traverse the data array */
    sum += data[i];             /* and sum the values */
  return sum /(REAL)n;          /* compute and return the mean value */
} /* mean() */

/*--------------------------------------------------------------------*/

static REAL SFXNAME(ksmall) (REAL *data, size_t n, size_t k, REAL *buf)
{                               /* --- k-smallest by quick selection */
  size_t a, b, l, r;            /* loop variables, section boundaries */
  REAL   p, x;                  /* pivot element, exchange buffer */

  assert(data && buf && (n > 0) && (k < n));
  if (n <= 2) {                 /* if at most two elements */
    if (n <= 1) return data[0]; /* return the only element or */
    if (k > 0) return (data[0] > data[1]) ? data[0] : data[1];
    else       return (data[0] > data[1]) ? data[1] : data[0];
  }                             /* return the greater/smaller element */
  if (buf != data)              /* copy data to the buffer if needed */
    memcpy(buf, data, n *sizeof(REAL));
  a = 0; b = n-1;               /* indices of first and last element */
  while (a < b) {               /* partition loop */
    l = a; r = b;               /* get the current index range */
    p = buf[(l+r)/2];           /* and the middle element as pivot */
    if (buf[l] > buf[r]) { x = buf[l]; buf[l] = buf[r]; buf[r] = x; }
    if      (p < buf[l]) p = buf[l];    /* compute median of three */
    else if (p > buf[r]) p = buf[r];    /* to find a better pivot */
    while (1) {                 /* split and exchange loop */
      while (buf[++l] < p);     /* skip smaller elems. on the left */
      while (buf[--r] > p);     /* skip greater elems. on the right */
      if (l >= r) break;        /* if less than two elements, abort */
      x = buf[l]; buf[l] = buf[r]; buf[r] = x;
    }                           /* otherwise exchange elements */
    if (l > k) b = l-1;         /* adapt the boundaries of */
    else       a = r+1;         /* the array section in which */
  }                             /* to search for k-smallest element */
  return buf[k];                /* return the k-smallest element */
}  /* ksmall() */

/*--------------------------------------------------------------------*/

void* SFXNAME(binarize) (REAL *data, int N, int T, REAL *thhs, int bpi)
{                               /* --- binarize data with thresholds */
  int  i, k;                    /* loop variables */
  int  X;                       /* size of (padded) data arrays */
  REAL t;                       /* threshold for binarization */
  void *bits;                   /* binarized data */
  REAL *buf = NULL;             /* buffer for median finding */

  assert(data && (N > 0) && (T > 0));
  /* default is bpi = 32; if bpi is not in {16,64,128,256,512}, */
  /* this function behaves as if bpi = 32. */
  if (thhs == BIN_MEDIAN) {     /* if to binarize with median */
    buf = malloc((size_t)T *sizeof(REAL));
    if (!buf) return NULL;      /* allocate buffer for median finding */
  }                             /* (avoid destroying the data) */

  if ((bpi ==  64) || (bpi == 128)    /* if 64/128/256/512 */
  ||  (bpi == 256) || (bpi == 512)) { /* bits per integer */
    uint64_t *u64;              /* compute binarized array size and */
    X = (bpi == 512) ? 8 *((T+511) >> 9)         /* allocate memory */
      : (bpi == 256) ? 4 *((T+255) >> 8)
      : (bpi == 128) ? 2 *((T+127) >> 7) : (T+63) >> 6;
    bits = u64 = calloc((size_t)N *(size_t)X, sizeof(uint64_t));
    if (!bits) { free(buf); return NULL; }
    for (i = 0; i < N; i++) {   /* traverse the data arrays */
      if      (thhs == BIN_MEAN)
        t = SFXNAME(mean)  (data+i*T, (size_t)T);
      else if (thhs == BIN_MEDIAN)
        t = SFXNAME(ksmall)(data+i*T, (size_t)T, (size_t)T/2, buf);
      else t = thhs[i];         /* get threshold for binarization */
      u64[i*X+X-1] = 0;         /* clear last entry (for 128 bits) */
      for (k = 0; k < T; k++)   /* set bits for all values */
        if (data[i*T+k] >= t)   /* at least as large as the threshold */
          u64[i*X+(k >> 6)] |= (uint64_t)1 << (k & 63);
    } }
  else if (bpi == 16) {         /* if 16 bits per integer */
    uint16_t *u16;              /* compute binarized array size and */
    X = (T+15) >> 4;            /* allocate memory for binarized data */
    bits = u16 = calloc((size_t)N *(size_t)X, sizeof(uint16_t));
    if (!bits) { free(buf); return NULL; }
    for (i = 0; i < N; i++) {   /* traverse the data arrays */
      if      (thhs == BIN_MEAN)
        t = SFXNAME(mean)  (data+i*T, (size_t)T);
      else if (thhs == BIN_MEDIAN)
        t = SFXNAME(ksmall)(data+i*T, (size_t)T, (size_t)T/2, buf);
      else t = thhs[i];         /* get threshold for binarization */
      for (k = 0; k < T; k++)   /* set bits for all values */
        if (data[i*T+k] >= t)   /* at least as large as the threshold */
          u16[i*X+(k >> 4)] = (uint16_t)(u16[i*X+(k >> 4)] | (1 << (k & 15)));
    } }
  else {                        /* if 32 bits per integer (default) */
    uint32_t *u32;              /* compute binarized array size and */
    X = (T+31) >> 5;            /* allocate memory for binarized data */
    bits = u32 = calloc((size_t)N *(size_t)X, sizeof(uint32_t));
    if (!bits) { free(buf); return NULL; }
    for (i = 0; i < N; i++) {   /* traverse the data arrays */
      if      (thhs == BIN_MEAN)
        t = SFXNAME(mean)  (data+i*T, (size_t)T);
      else if (thhs == BIN_MEDIAN)
        t = SFXNAME(ksmall)(data+i*T, (size_t)T, (size_t)T/2, buf);
      else t = thhs[i];         /* get threshold for binarization */
      for (k = 0; k < T; k++)   /* set bits for all values */
        if (data[i*T+k] >= t)   /* at least as large as the threshold */
          u32[i*X+(k >> 5)] |= (uint32_t)1 << (k & 31);
    }
  }
  if (thhs == BIN_MEDIAN)       /* if to binarize with median, */
    free(buf);                  /* delete the allocated buffer */
  return bits;                  /* return created binary data */
}  /* binarize() */

/*--------------------------------------------------------------------*/

#ifndef __BINARIZE_2__
#define __BINARIZE_2__
#undef REAL
#undef SUFFIX
#define REAL    double          /* second pass: double precision */
#define SUFFIX  _dbl            /* function name suffix is '_dbl' */
#include "binarize.c"           /* process source recursively */
#endif

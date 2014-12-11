/*----------------------------------------------------------------------
  File    : pcc.c
  Contents: compute pairwise Pearson correlation coefficients
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#ifndef _WIN32                  /* if Linux/Unix system */
#define _POSIX_C_SOURCE 200809L /* needed for clock_gettime() */
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <pthread.h>
#endif

#include "cpuinfo.h"
#include "clamp.h"
#include "pcc.h"

/*----------------------------------------------------------------------
  Data Type Definition / Recursion Handling
----------------------------------------------------------------------*/
#ifdef REAL                     /* if REAL is defined, */
#  undef  _PCC_PASS             /* ensure _PCC_PASS is undefined */
#  define _PCC_PASS 0           /* define macro for single pass */
#  ifndef SUFFIX                /* function names get no suffix */
#  define SUFFIX                /* (only single set of functions) */
#  endif
#elif !defined _PCC_PASS        /* if in first pass of two */
#  undef  _PCC_PASS             /* ensure _PCC_PASS is undefined */
#  define _PCC_PASS 1           /* define macro for first pass */
#  define REAL      float       /* first pass: single precision */
#  define SUFFIX    _flt        /* function name suffix is '_flt' */
#else                           /* if in second pass of two */
#  undef  _PCC_PASS             /* ensure _PCC_PASS is undefined */
#  define _PCC_PASS 2           /* define macro for second pass */
#  define REAL      double      /* second pass: double precision */
#  define SUFFIX    _dbl        /* function name suffix is '_dbl' */
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

#ifndef SFXNAME                 /* macros to generate function names */
#define SFXNAME(n)      SFXNAME_1(n,SUFFIX)
#define SFXNAME_1(n,s)  SFXNAME_2(n,s)
#define SFXNAME_2(n,s)  n##s    /* the two step recursion is needed */
#endif                          /* to ensure proper expansion */

/*----------------------------------------------------------------------
  Processing Options
----------------------------------------------------------------------*/
#define LOWER       0           /* compute lower triangular matrix */
                                /* default: upper */
#define ROWS        0           /* use row pointers for tiling */
                                /* default: index computation */
#define USERTHD     0           /* copy user flag for threading */
                                /* default: choose automatically */
#define clamp       clamp1      /* choose clamping implementation */

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#ifndef M_PI                    /* if not defined (no Posix source) */
#define M_PI            3.14159265358979323846
#endif                          /* pi in double precision */

#ifndef DELROWS                 /* if not yet defined */
#if ROWS                        /* if to use an array of row starts */
#define DELROWS         free(rows)
#else                           /* if to use index computation */
#define DELROWS                 /* there are no rows to delete */
#endif
#endif

#if !ROWS && !defined INDEX     /* if index function not defined yet */
#if LOWER                       /* if lower triangular matrix */
#define INDEX(i,j,N)    ((size_t)(i)*(size_t)((i)-1)/2+(size_t)(j))
#else                           /* if upper triangular matrix */
#define INDEX(i,j,N)    ((size_t)(i)*((size_t)(N)+(size_t)(N) \
                        -(size_t)(i)-3)/2-1+(size_t)(j))
#endif                          /* index computation for result */
#endif

#ifndef GET_THREAD              /* if not yet defined */
#if USERTHD                     /* if to respect user flag */
#define GET_THREAD(var) ((var) & PCC_THREAD)
#else                           /* copy the threading flag */
#define GET_THREAD(var) PCC_THREAD
#endif                          /* otherwise force threading flag */
#endif

#ifndef THREAD                  /* if not yet defined */
#ifdef _WIN32                   /* if Microsoft Windows system */
#define THREAD          HANDLE  /* threads are identified by handles */
#define THREAD_OK       0       /* return value is DWORD */
#define WORKERDEF(n,p)  DWORD WINAPI SFXNAME(n) (LPVOID p)
#else                           /* if Linux/Unix system */
#define THREAD          pthread_t  /* use the POSIX thread type */
#define THREAD_OK       NULL    /* return value is void* */
#define WORKERDEF(n,p)  void*        SFXNAME(n) (void* p)
#endif                          /* definition of a worker function */
#endif

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- thread worker data --- */
  REAL *diff;                   /* differences to mean value */
  REAL *rssd;                   /* roots of sum of squared deviats. */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* array of row starts */
  #else                         /* if to use index computation */
  REAL *res;                    /* result array */
  #endif
  int  N;                       /* number of series */
  int  T;                       /* original length of each serie */
  int  X;                       /* padded   length of each serie */
  int  tile;                    /* tile size (columns to group) */
  int  s, e;                    /* index of start and end serie */
} SFXNAME(WORK);                /* (thread worker data) */

#ifndef WORKERTYPE              /* if not yet defined */
#define WORKERTYPE              /* define worker function type */
#ifdef _WIN32                   /* if Microsoft Windows system */
typedef DWORD WINAPI WORKER (LPVOID);
#else                           /* if Linux/Unix system */
typedef void*        WORKER (void*);
#endif                          /* worker for parallel execution */
#endif

/*----------------------------------------------------------------------
  Auxiliary Functions
----------------------------------------------------------------------*/
#if ROWS                        /* if to use an array of row starts */

static REAL** SFXNAME(make_rows) (REAL *res, int N)
{                               /* --- create array of row starts */
  int  i;                       /* loop variable */
  REAL **rows;                  /* starts of output rows */

  assert(res && (N > 0));       /* check the function arguments */
  rows = malloc((size_t)N *sizeof(REAL*));
  if (!rows) return NULL;       /* allocate array of row starts */
  for (i = 0; i < N; i++) {     /* traverse the matrix rows */
    #if LOWER                   /* if lower triangular matrix */
    rows[i] = res;       res += i; }
    #else                       /* if upper triangular matrix */
    rows[i] = res-(i+1); res += N-(i+1); }
    #endif                      /* compute addresses of row starts */
  return rows;                  /* return the created array */
}  /* make_rows() */

#endif  /* #if ROWS */
/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/

void SFXNAME(init_naive) (REAL *data, int N, int T,
                                 REAL *diff, int X, REAL *rssd)
{                               /* --- init. differences and rssds */
  int  i, k;                    /* loop variables */
  REAL d, sum, sqr;             /* difference, sum (of squares) */

  assert(diff && rssd && (T > 0) && (X >= T));
  for (i = 0; i < N; i++) {     /* traverse the data arrays */
    for (sum = 0, k = 0; k < T; k++)
      sum += data[i*T+k];       /* sum the data values and */
    sum /= (REAL)T;             /* compute their arithmetic mean */
    for (sqr = 0, k = 0; k < T; k++) {
      d = diff[i*X+k] = data[i*T+k] -sum;
      sqr += d*d;               /* compute differences to mean and */
    }                           /* sum of squared deviations from it */
    for ( ; k < X; k++) diff[i*X+k] = 0;
    rssd[i] = (REAL)sqrt(sqr);  /* pad array with differences to mean */
  }                             /* root of sum of squared deviations */
}  /* init_naive() */

/*--------------------------------------------------------------------*/

static int SFXNAME(pcc_naive) (REAL *diff, REAL *rssd, REAL *res,
                               int N, int T, int X)
{                               /* --- compute Pearson's corr. coeff. */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(diff && rssd && res && (N > 0) && (T > 0));
  for (i = 0; i < N; i++) {     /* traverse the data arrays */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* traverse the pair of series */
      sum = SFXNAME(pair_naive)(diff+i*X, diff+j*X, T);
      *res++ = SFXNAME(clamp)(sum /(rssd[i] *rssd[j]), -1.0, +1.0);
    }                           /* compute correlation coefficient */
  }                             /* (Pearson's r) */
  return 0;                     /* return 'ok' */
}  /* pcc_naive() */

/*--------------------------------------------------------------------*/

static int SFXNAME(pcc_naive_tiled) (REAL *diff, REAL *rssd, REAL *res,
                                     int N, int T, int X, int tile)
{                               /* --- compute Pearson's corr. coeff. */
  int  i, j, m, e;              /* loop variables */
  REAL sum;                     /* sum of products */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(diff && rssd && res && (N > 0) && (T > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  for (m = 0; m < N; m += tile) {  /* traverse the stripes/tiles */
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* traverse the pair of series */
        sum = SFXNAME(pair_naive)(diff+i*X, diff+j*X, T);
        #if ROWS                /* if row starts are available */
        rows[i][j]
        #else                   /* if to use index computation */
        res[INDEX(i,j,N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(rssd[i] *rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
  }                             /* with a linear index */

  DELROWS;                      /* deallocate row starts */
  return 0;                     /* return 'ok' */
}  /* pcc_naive_tiled() */

/*--------------------------------------------------------------------*/

static void SFXNAME(rct_naive) (SFXNAME(WORK) *w,
                                int ra, int rb, int ca, int cb)
{                               /* --- cache-oblivious, rectangle */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(w                      /* check the funktion arguments */
  &&    (ra >= 0) && (ra < w->N) && (rb > ra) && (rb <= w->N)
  &&    (ca >= 0) && (ca < w->N) && (cb > ca) && (cb <= w->N));
  j = (ra == ca) ? 2*w->tile : w->tile;
  if (rb-ra > j) {              /* if larger than minimal tile size */
    i = (ra+rb)/2;              /* halven the tile size and */
    j = (ca+cb)/2;              /* process parts recursively */
    #if LOWER                   /* if lower triangular matrix */
    SFXNAME(rct_naive)(w, ra, i, ca, j);
    SFXNAME(rct_naive)(w, i, rb, ca, j);
    SFXNAME(rct_naive)(w, i, rb, j, cb);
    SFXNAME(rct_naive)(w, ra, i, j, cb); }
    #else                       /* if upper triangular matrix */
    SFXNAME(rct_naive)(w, ra, i, ca, j);
    SFXNAME(rct_naive)(w, ra, i, j, cb);
    SFXNAME(rct_naive)(w, i, rb, j, cb);
    SFXNAME(rct_naive)(w, i, rb, ca, j); }
    #endif                      /* split into four rectangles */
  else {                        /* if no larger than min. tile size */
    for (i = ra; i < rb; i++) { /* traverse the rows */
      for (j = ca; j < cb; j++){/* traverse the columns */
        sum = SFXNAME(pair_naive)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
  }
}  /* rct_naive() */

/*--------------------------------------------------------------------*/

static void SFXNAME(trg_naive) (SFXNAME(WORK) *w, int a, int b)
{                               /* --- cache-oblivious, triangle */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(w                      /* check the funktion arguments */
  &&    (a >= 0) && (a < w->N) && (b > a) && (b <= w->N));
  if (b-a > 2*w->tile) {        /* if larger than minimal tile size */
    i = (a+b)/2;                /* halven the tile size */
    SFXNAME(trg_naive)(w, a, i);
    #if LOWER                   /* if lower triangular matrix */
    SFXNAME(rct_naive)(w, i, b, a, i);
    #else                       /* if upper triangular matrix */
    SFXNAME(rct_naive)(w, a, i, i, b);
    #endif                      /* compute corresponding square */
    SFXNAME(trg_naive)(w, i, b); }
  else {                        /* if no larger than min. tile size */
    for (i = a; i < b; i++) {   /* traverse the data arrays */
      #if LOWER                 /* if lower triangular matrix */
      for (j = a;   j < i; j++){/* traverse the smaller indices */
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < b; j++){/* traverse the greater indices */
      #endif                    /* traverse the pair of series */
        sum = SFXNAME(pair_naive)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
  }
}  /* trg_naive() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_naive, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {   /* traverse row indices */
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* traverse column indices */
        sum = SFXNAME(pair_naive)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_naive() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_naive_tiled, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int  i, j, m, e;              /* loop variables */
  REAL sum;                     /* sum of products */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif                    /* traverse the pair of series */
          sum = SFXNAME(pair_naive)(w->diff+i*w->X,w->diff+j*w->X,w->T);
          #if ROWS              /* if row starts are available */
          w->rows[i][j]
          #else                 /* if to use index computation */
          w->res[INDEX(i,j,w->N)]
          #endif                /* compute correlation coefficient */
            = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
        }                       /* (Pearson's r) and */
      }                         /* store it in the result */
    }                           /* with a linear index */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_naive_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_naive_cobl, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int s, e, k;                  /* loop variables */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    SFXNAME(trg_naive)(w, w->s, w->e);
    k = w->e -w->s;             /* compute leading triangle */
    #if LOWER                   /* if lower triangular matrix */
    for (s = w->e; s < w->N; s += k) {
      e = (s+k < w->N) ? s+k : w->N;
    #else                       /* if upper triangular matrix */
    for (s = 0;    s < w->s; s += k) {
      e = (s+k < w->s) ? s+k : w->s;
    #endif                      /* traverse the row indices */
      SFXNAME(rct_naive)(w, s, e, w->s, w->e);
    }                           /* split the strip into squares */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    s    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > s)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = s;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_naive_cobl() */

/*--------------------------------------------------------------------*/
#ifdef __SSE2__

void SFXNAME(init_sse2) (REAL *data, int N, int T,
                                REAL *diff, int X, REAL *rssd)
{                               /* --- init. differences and rssds */
  int     i, k, n;              /* loop variables */
  REAL    d, sum, sqr;          /* difference, sum (of squares) */
  #if REAL_IS_DOUBLE            /* data is double precision */
  __m128d s, e, m;              /* registers for SSE2 computations */
  #else                         /* data is single precision */
  __m128  s, e, m;              /* registers for SSE2 computations */
  #endif

  assert(diff && rssd && (T > 0) && (X >= T));
  for (i = 0; i < N; i++) {     /* traverse the data arrays */

    /* --- compute mean value --- */
    n =    (int)(((size_t)(data+i*T)   & 15) /sizeof(REAL));
    if      (n > T) n = T;      /* compute offset to alignment */
    else if (n > 0) n = (int)(16/sizeof(REAL)) -n;
    for (sum = 0, k = 0; k < n; k++)
      sum += data[i*T+k];       /* sum the initial data values */
    n = T -(int)(((size_t)(data+i*T+T) & 15) /sizeof(REAL));
    #if REAL_IS_DOUBLE          /* data is double precision */
    s = _mm_setzero_pd();       /* initialize the sum (2 values) */
    for ( ; k < n; k += 2)      /* sum two values in parallel */
      s = _mm_add_pd(s, _mm_load_pd(data+i*T+k));
    s = _mm_add_pd(s, _mm_shuffle_pd(s, s, 1));
    sum += _mm_cvtsd_f64(s);    /* sum two sums horizontally */
    #else                       /* data is single precision */
    s = _mm_setzero_ps();       /* initialize the sum (4 values) */
    for ( ; k < n; k += 4)      /* sum four values in parallel */
      s = _mm_add_ps(s, _mm_load_ps(data+i*T+k));
    s = _mm_add_ps(s, _mm_movehl_ps(s, s));
    s = _mm_add_ss(s, _mm_shuffle_ps(s, s, 1));
    sum += _mm_cvtss_f32(s);    /* sum four sums horizontally */
    #endif
    for ( ; k < T; k++)         /* traverse the remaining values */
      sum += data[i*T+k];       /* sum the data values and */
    sum /= (REAL)T;             /* compute their arithmetic mean */

    /* --- compute differences to mean and rssds --- */
    #if REAL_IS_DOUBLE          /* data is double precision */
    m = _mm_set1_pd(sum);       /* initialize the mean and */
    s = _mm_setzero_pd();       /* the sum of squares (2 values) */
    for (k = 0; k < (T & ~1); k += 2) {
      e = _mm_sub_pd(_mm_loadu_pd(data+i*T+k), m);
      _mm_store_pd(diff+i*X+k, e);
      s = _mm_add_pd(s, _mm_mul_pd(e, e));
    }                           /* compute two diffs/squares in par. */
    s = _mm_add_pd(s, _mm_shuffle_pd(s, s, 1));
    sqr = _mm_cvtsd_f64(s);     /* sum two sums horizontally */
    #else                       /* data is single precision */
    m = _mm_set1_ps(sum);       /* initialize the mean and */
    s = _mm_setzero_ps();       /* initialize the sum (4 values) */
    for (k = 0; k < (T & ~3); k += 4) {
      e = _mm_sub_ps(_mm_loadu_ps(data+i*T+k), m);
      _mm_store_ps(diff+i*X+k, e);
      s = _mm_add_ps(s, _mm_mul_ps(e, e));
    }                           /* compute four diffs/squares in par. */
    s = _mm_add_ps(s, _mm_movehl_ps(s, s));
    s = _mm_add_ss(s, _mm_shuffle_ps(s, s, 1));
    sqr = _mm_cvtss_f32(s);     /* sum four sums horizontally */
    #endif
    for ( ; k < T; k++) {       /* traverse the remaining values */
      d = diff[i*X+k] = data[i*T+k] -sum;
      sqr += d*d;               /* compute differences to mean and */
    }                           /* sum of squared deviations from it */
    for ( ; k < X; k++) diff[i*X+k] = 0;
    rssd[i] = (REAL)sqrt(sqr);  /* pad array with differences to mean */
  }                             /* root of sum of squared deviations */
}  /* init_sse2() */

/*--------------------------------------------------------------------*/

static int SFXNAME(pcc_sse2) (REAL *diff, REAL *rssd, REAL *res,
                              int N, int T, int X)
{                               /* --- compute Pearson's corr. coeff. */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(diff && rssd && res && (N > 0) && (T > 0));
  for (i = 0; i < N; i++) {     /* traverse the data arrays */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* traverse the pair of series */
      sum = SFXNAME(pair_sse2)(diff+i*X, diff+j*X, T);
      *res++ = SFXNAME(clamp)(sum /(rssd[i] *rssd[j]), -1.0, +1.0);
    }                           /* sum the values of the sum and */
  }                             /* compute correlation coefficient */
  return 0;                     /* return 'ok' */
}  /* pcc_sse2() */

/*--------------------------------------------------------------------*/

static int SFXNAME(pcc_sse2_tiled) (REAL *diff, REAL *rssd, REAL *res,
                                    int N, int T, int X, int tile)
{                               /* --- compute Pearson's corr. coeff. */
  int  i, j, m, e;              /* loop variables */
  REAL sum;                     /* sum of products */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(diff && rssd && res && (N > 0) && (T > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  for (m = 0; m < N; m += tile) {
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* traverse the pair of series */
        sum = SFXNAME(pair_sse2)(diff+i*X, diff+j*X, T);
        #if ROWS                /* if row starts are available */
        rows[i][j]
        #else                   /* if to use index computation */
        res[INDEX(i,j,N)]
        #endif
          = SFXNAME(clamp)(sum /(rssd[i] *rssd[j]), -1.0, +1.0);
      }                         /* sum the four values of the sum and */
    }                           /* compute correlation coefficient */
  }

  DELROWS;                      /* deallocate row starts */
  return 0;                     /* return 'ok' */
}  /* pcc_sse2_tiled() */

/*--------------------------------------------------------------------*/

static void SFXNAME(rct_sse2) (SFXNAME(WORK) *w,
                               int ra, int rb, int ca, int cb)
{                               /* --- cache-oblivious, rectangle */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(w                      /* check the funktion arguments */
  &&    (ra >= 0) && (ra < w->N) && (rb > ra) && (rb <= w->N)
  &&    (ca >= 0) && (ca < w->N) && (cb > ca) && (cb <= w->N));
  j = (ra == ca) ? 2*w->tile : w->tile;
  if (rb-ra > j) {              /* if larger than minimal tile size */
    i = (ra+rb)/2;              /* halven the tile size and */
    j = (ca+cb)/2;              /* process parts recursively */
    #if LOWER                   /* if lower triangular matrix */
    SFXNAME(rct_sse2)(w, ra, i, ca, j);
    SFXNAME(rct_sse2)(w, i, rb, ca, j);
    SFXNAME(rct_sse2)(w, i, rb, j, cb);
    SFXNAME(rct_sse2)(w, ra, i, j, cb); }
    #else                       /* if upper triangular matrix */
    SFXNAME(rct_sse2)(w, ra, i, ca, j);
    SFXNAME(rct_sse2)(w, ra, i, j, cb);
    SFXNAME(rct_sse2)(w, i, rb, j, cb);
    SFXNAME(rct_sse2)(w, i, rb, ca, j); }
    #endif                      /* split into four rectangles */
  else {                        /* if no larger than min. tile size */
    for (i = ra; i < rb; i++) { /* traverse the rows */
      for (j = ca; j < cb; j++){/* traverse the columns */
        sum = SFXNAME(pair_sse2)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
  }
}  /* rct_sse2() */

/*--------------------------------------------------------------------*/

static void SFXNAME(trg_sse2) (SFXNAME(WORK) *w, int a, int b)
{                               /* --- cache-oblivious, triangle */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(w                      /* check the funktion arguments */
  &&    (a >= 0) && (a < w->N) && (b > a) && (b <= w->N));
  if (b-a > 2*w->tile) {        /* if larger than minimal tile size */
    i = (a+b)/2;                /* halven the tile size */
    SFXNAME(trg_sse2)(w, a, i);
    #if LOWER                   /* if lower triangular matrix */
    SFXNAME(rct_sse2)(w, i, b, a, i);
    #else                       /* if upper triangular matrix */
    SFXNAME(rct_sse2)(w, a, i, i, b);
    #endif                      /* compute corresponding square */
    SFXNAME(trg_sse2)(w, i, b); }
  else {                        /* if no larger than min. tile size */
    for (i = a; i < b; i++) {   /* traverse the data arrays */
      #if LOWER                 /* if lower triangular matrix */
      for (j = a;   j < i; j++){/* traverse the smaller indices */
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < b; j++){/* traverse the greater indices */
      #endif                    /* traverse the pair of series */
        sum = SFXNAME(pair_sse2)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
  }
}  /* trg_sse2() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_sse2, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* traverse the row indices */
        sum = SFXNAME(pair_sse2)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_sse2() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_sse2_tiled, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int  i, j, m, e;              /* loop variables */
  REAL sum;                     /* sum of products */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif
          sum = SFXNAME(pair_sse2)(w->diff+i*w->X,w->diff+j*w->X, w->T);
          #if ROWS              /* if row starts are available */
          w->rows[i][j]
          #else                 /* if to use index computation */
          w->res[INDEX(i,j,w->N)]
          #endif
            = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
        }                       /* compute correlation coefficient */
      }                         /* (Pearson's r) */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_sse2_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_sse2_cobl, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int s, e, k;                  /* loop variables */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    SFXNAME(trg_sse2)(w, w->s, w->e);
    k = w->e -w->s;             /* compute leading triangle */
    #if LOWER                   /* if lower triangular matrix */
    for (s = w->e; s < w->N; s += k) {
      e = (s+k < w->N) ? s+k : w->N;
    #else                       /* if upper triangular matrix */
    for (s = 0;    s < w->s; s += k) {
      e = (s+k < w->s) ? s+k : w->s;
    #endif                      /* traverse the row indices */
      SFXNAME(rct_sse2)(w, s, e, w->s, w->e);
    }                           /* split the strip into squares */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    s    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > s)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = s;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_sse2_cobl() */

#endif  /* #ifdef __SSE2__ */
/*--------------------------------------------------------------------*/
#ifdef __AVX__

void SFXNAME(init_avx) (REAL *data, int N, int T,
                               REAL *diff, int X, REAL *rssd)
{                               /* --- init. differences and rssds */
  int     i, k, n;              /* loop variables */
  REAL    d, sum, sqr;          /* difference, sum (of squares) */
  #if REAL_IS_DOUBLE            /* data is double precision */
  __m256d s, e, m;              /* registers for AVX  computations */
  __m128d x;                    /* register  for SSE2 computations */
  #else                         /* data is single precision */
  __m256  s, e, m;              /* registers for AVX  computations */
  __m128  x;                    /* register  for SSE2 computations */
  #endif

  assert(diff && rssd && (T > 0) && (X >= T));
  for (i = 0; i < N; i++) {     /* traverse the data arrays */

    /* --- compute mean value --- */
    n =    (int)(((size_t)(data+i*T)   & 31) /sizeof(REAL));
    if      (n > T) n = T;      /* compute offset to alignment */
    else if (n > 0) n = (int)(32/sizeof(REAL)) -n;
    for (sum = 0, k = 0; k < n; k++)
      sum += data[i*T+k];       /* sum the initial data values */
    n = T -(int)(((size_t)(data+i*T+T) & 31) /sizeof(REAL));
    #if REAL_IS_DOUBLE          /* data is double precision */
    s = _mm256_setzero_pd();    /* initialize the sum (4 values) */
    for ( ; k < n; k += 4)      /* sum four values in parallel */
      s = _mm256_add_pd(s, _mm256_load_pd(data+i*T+k));
    x = _mm_add_pd(_mm256_extractf128_pd(s, 0),
                   _mm256_extractf128_pd(s, 1));
    x = _mm_add_pd(x, _mm_shuffle_pd(x, x, 1));
    sum += _mm_cvtsd_f64(x);    /* sum four sums horizontally */
    #else                       /* data is single precision */
    s = _mm256_setzero_ps();    /* initialize the sum (8 values) */
    for ( ; k < n; k += 8)      /* sum eight values in parallel */
      s = _mm256_add_ps(s, _mm256_load_ps(data+i*T+k));
    s = _mm256_hadd_ps(s, s);   /* do horizontal sums in upper */
    s = _mm256_hadd_ps(s, s);   /* and lower half of the register */
    x = _mm_add_ss(_mm256_castps256_ps128(s),
                   _mm256_extractf128_ps(s, 1));
    sum += _mm_cvtss_f32(x);    /* sum eight sums horizontally */
    #endif
    for ( ; k < T; k++)         /* traverse the remaining values */
      sum += data[i*T+k];       /* sum the data values and */
    sum /= (REAL)T;             /* compute their arithmetic mean */

    /* --- compute differences to mean and rssds --- */
    #if REAL_IS_DOUBLE          /* data is double precision */
    m = _mm256_set1_pd(sum);    /* initialize the mean and */
    s = _mm256_setzero_pd();    /* the sum of squares (4 values) */
    for (k = 0; k < (T & ~3); k += 4) {
      e = _mm256_sub_pd(_mm256_loadu_pd(data+i*T+k), m);
      _mm256_store_pd(diff+i*X+k, e);
      s = _mm256_add_pd(s, _mm256_mul_pd(e, e));
    }                           /* compute four values in parallel */
    x = _mm_add_pd(_mm256_extractf128_pd(s, 0),
                   _mm256_extractf128_pd(s, 1));
    x = _mm_add_pd(x, _mm_shuffle_pd(x, x, 1));
    sqr = _mm_cvtsd_f64(x);     /* sum four sums horizontally */
    #else                       /* data is single precision */
    m = _mm256_set1_ps(sum);    /* initialize the mean and */
    s = _mm256_setzero_ps();    /* the sum of squares (8 values) */
    for (k = 0; k < (T & ~7); k += 8) {
      e = _mm256_sub_ps(_mm256_loadu_ps(data+i*T+k), m);
      _mm256_store_ps(diff+i*X+k, e);
      s = _mm256_add_ps(s, _mm256_mul_ps(e, e));
    }                           /* compute eight values in parallel */
    s = _mm256_hadd_ps(s, s);   /* do horizontal sums in upper */
    s = _mm256_hadd_ps(s, s);   /* and lower half of the register */
    x = _mm_add_ss(_mm256_castps256_ps128(s),
                   _mm256_extractf128_ps(s, 1));
    sqr = _mm_cvtss_f32(x);     /* sum eight sums horizontally */
    #endif
    for ( ; k < T; k++) {       /* traverse the remaining values */
      d = diff[i*X+k] = data[i*T+k] -sum;
      sqr += d*d;               /* compute differences to mean and */
    }                           /* sum of squared deviations from it */
    for ( ; k < X; k++) diff[i*X+k] = 0;
    rssd[i] = (REAL)sqrt(sqr);  /* pad array with differences to mean */
  }                             /* root of sum of squared deviations */
}  /* init_avx() */

/*--------------------------------------------------------------------*/

static int SFXNAME(pcc_avx) (REAL *diff, REAL *rssd, REAL *res,
                             int N, int T, int X)
{                               /* --- compute Pearson's corr. coeff. */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(diff && rssd && res && (N > 0) && (T > 0));
  for (i = 0; i < N; i++) {     /* traverse the data arrays */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* traverse the pair of series */
      sum = SFXNAME(pair_avx)(diff+i*X, diff+j*X, T);
      *res++ = SFXNAME(clamp)(sum /(rssd[i] *rssd[j]), -1.0, +1.0);
    }                           /* sum the values of the sum and */
  }                             /* compute correlation coefficient */
  return 0;                     /* return 'ok' */
}  /* pcc_avx() */

/*--------------------------------------------------------------------*/

static int SFXNAME(pcc_avx_tiled) (REAL *diff, REAL *rssd, REAL *res,
                                   int N, int T, int X, int tile)
{                               /* --- compute Pearson's corr. coeff. */
  int  i, j, m, e;              /* loop variables */
  REAL sum;                     /* sum of products */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(diff && rssd && res && (N > 0) && (T > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  for (m = 0; m < N; m += tile) {
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* traverse the pair of series */
        sum = SFXNAME(pair_avx)(diff+i*X, diff+j*X, T);
        #if ROWS                /* if row starts are available */
        rows[i][j]
        #else                   /* if to use index computation */
        res[INDEX(i,j,N)]
        #endif
          = SFXNAME(clamp)(sum /(rssd[i] *rssd[j]), -1.0, +1.0);
      }                         /* sum the four values of the sum and */
    }                           /* compute correlation coefficient */
  }

  DELROWS;                      /* deallocate row starts */
  return 0;                     /* return 'ok' */
}  /* pcc_avx_tiled() */

/*--------------------------------------------------------------------*/

static void SFXNAME(rct_avx) (SFXNAME(WORK) *w,
                              int ra, int rb, int ca, int cb)
{                               /* --- cache-oblivious, rectangle */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(w                      /* check the funktion arguments */
  &&    (ra >= 0) && (ra < w->N) && (rb > ra) && (rb <= w->N)
  &&    (ca >= 0) && (ca < w->N) && (cb > ca) && (cb <= w->N));
  j = (ra == ca) ? 2*w->tile : w->tile;
  if (rb-ra > j) {              /* if larger than minimal tile size */
    i = (ra+rb)/2;              /* halven the tile size and */
    j = (ca+cb)/2;              /* process parts recursively */
    #if LOWER                   /* if lower triangular matrix */
    SFXNAME(rct_avx)(w, ra, i, ca, j);
    SFXNAME(rct_avx)(w, i, rb, ca, j);
    SFXNAME(rct_avx)(w, i, rb, j, cb);
    SFXNAME(rct_avx)(w, ra, i, j, cb); }
    #else                       /* if upper triangular matrix */
    SFXNAME(rct_avx)(w, ra, i, ca, j);
    SFXNAME(rct_avx)(w, ra, i, j, cb);
    SFXNAME(rct_avx)(w, i, rb, j, cb);
    SFXNAME(rct_avx)(w, i, rb, ca, j); }
    #endif                      /* split into four rectangles */
  else {                        /* if no larger than min. tile size */
    for (i = ra; i < rb; i++) { /* traverse the rows */
      for (j = ca; j < cb; j++){/* traverse the columns */
        sum = SFXNAME(pair_avx)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
  }
}  /* rct_avx() */

/*--------------------------------------------------------------------*/

static void SFXNAME(trg_avx) (SFXNAME(WORK) *w, int a, int b)
{                               /* --- cache-oblivious, triangle */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(w                      /* check the funktion arguments */
  &&    (a >= 0) && (a < w->N) && (b > a) && (b <= w->N));
  if (b-a > 2*w->tile) {        /* if larger than minimal tile size */
    i = (a+b)/2;                /* halven the tile size */
    SFXNAME(trg_avx)(w, a, i);
    #if LOWER                   /* if lower triangular matrix */
    SFXNAME(rct_avx)(w, i, b, a, i);
    #else                       /* if upper triangular matrix */
    SFXNAME(rct_avx)(w, a, i, i, b);
    #endif                      /* compute corresponding square */
    SFXNAME(trg_avx)(w, i, b); }
  else {                        /* if no larger than min. tile size */
    for (i = a; i < b; i++) {   /* traverse the data arrays */
      #if LOWER                 /* if lower triangular matrix */
      for (j = a;   j < i; j++){/* traverse the smaller indices */
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < b; j++){/* traverse the greater indices */
      #endif                    /* traverse the pair of series */
        sum = SFXNAME(pair_avx)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
  }
}  /* trg_avx() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_avx, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int  i, j;                    /* loop variables */
  REAL sum;                     /* sum of products */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* traverse the row indices */
        sum = SFXNAME(pair_avx)(w->diff+i*w->X, w->diff+j*w->X, w->T);
        #if ROWS                /* if row starts are available */
        w->rows[i][j]
        #else                   /* if to use index computation */
        w->res[INDEX(i,j,w->N)]
        #endif                  /* compute correlation coefficient */
          = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
      }                         /* (Pearson's r) and */
    }                           /* store it in the result */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_avx() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_avx_tiled, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int  i, j, m, e;              /* loop variables */
  REAL sum;                     /* sum of products */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif
          sum = SFXNAME(pair_avx)(w->diff+i*w->X,w->diff+j*w->X, w->T);
          #if ROWS              /* if row starts are available */
          w->rows[i][j]
          #else                 /* if to use index computation */
          w->res[INDEX(i,j,w->N)]
          #endif
            = SFXNAME(clamp)(sum /(w->rssd[i] *w->rssd[j]), -1.0, +1.0);
        }                       /* compute correlation coefficient */
      }                         /* (Pearson's r) */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_avx_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_avx_cobl, p)
{                               /* --- compute pcc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int s, e, k;                  /* loop variables */

  assert(p);                    /* check the function argument */
  while (1) {                   /* process two strip parts */
    SFXNAME(trg_avx)(w, w->s, w->e);
    k = w->e -w->s;             /* compute leading triangle */
    #if LOWER                   /* if lower triangular matrix */
    for (s = w->e; s < w->N; s += k) {
      e = (s+k < w->N) ? s+k : w->N;
    #else                       /* if upper triangular matrix */
    for (s = 0;    s < w->s; s += k) {
      e = (s+k < w->s) ? s+k : w->s;
    #endif                      /* traverse the row indices */
      SFXNAME(rct_avx)(w, s, e, w->s, w->e);
    }                           /* split the strip into squares */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    s    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > s)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = s;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_avx_cobl() */

#endif  /* #ifdef __AVX__ */
/*--------------------------------------------------------------------*/

static int SFXNAME(pcc_single) (REAL *diff, REAL *rssd, REAL *res,
                                int N, int T, int X,
                                int var, int tile)
{                               /* --- single thread version */
  SFXNAME(WORK) w;              /* cache-oblivious recursion data */

  assert(diff && rssd && res && (N > 0) && (T > 0) && (X >= T));
  if (var & PCC_COBL) {         /* if to use cache-oblivious version */
    w.diff = diff; w.rssd = rssd;
    w.N = N; w.T = T; w.X = X;  /* note the data parameters */
    w.tile = tile;              /* and the minimal tile size */
    #if ROWS                    /* if to use an array of row starts */
    w.rows = SFXNAME(make_rows)(res, N);
    if (!w.rows) return -1;     /* allocate and initialize the array */
    #else                       /* if to use index computation */
    w.res  = res;               /* note the result array */
    #endif
    switch (var & PCC_VARIANT){ /* evaluate the variant */
      #ifdef __AVX__            /* if AVX instructions available */
      case PCC_AVX:             /* if to use AVX computation */
        SFXNAME(trg_avx)  (&w, 0, N); break;
      #endif                    /* compute correlation coefficient */
      #ifdef __SSE2__           /* if SSE2 instructions available */
      case PCC_SSE2:            /* if to use SSE2 computation */
        SFXNAME(trg_sse2) (&w, 0, N); break;
      #endif                    /* compute correlation coefficient */
      default:                  /* if to use naive computation */
        SFXNAME(trg_naive)(&w, 0, N); break;
    }                           /* compute correlation coefficient */
    #if ROWS                    /* if to use an array of row starts */
    free(w.rows);               /* deallocate row starts */
    #endif
    return 0;                   /* return 'ok' */
  }

  switch (var & PCC_VARIANT) {  /* evaluate the variant */
    #ifdef __AVX__              /* if AVX instructions available */
    case PCC_AVX:               /* if to use AVX computation */
      if (var & PCC_TILED)      /* distinguished tiled version */
           SFXNAME(pcc_avx_tiled)  (diff, rssd, res, N, T, X, tile);
      else SFXNAME(pcc_avx)        (diff, rssd, res, N, T, X);
      break;                    /* compute correlation coefficient */
    #endif
    #ifdef __SSE2__             /* if SSE2 instructions available */
    case PCC_SSE2:              /* if to use SSE2 computation */
      if (var & PCC_TILED)      /* distinguished tiled version */
           SFXNAME(pcc_sse2_tiled) (diff, rssd, res, N, T, X, tile);
      else SFXNAME(pcc_sse2)       (diff, rssd, res, N, T, X);
      break;                    /* compute correlation coefficient */
    #endif
    default:                    /* if to use naive computation */
      if (var & PCC_TILED)      /* distinguished tiled version */
           SFXNAME(pcc_naive_tiled)(diff, rssd, res, N, T, X, tile);
      else SFXNAME(pcc_naive)      (diff, rssd, res, N, T, X);
      break;                    /* compute correlation coefficient */
  }  /* switch (var & PCC_VARIANT) */

  return 0;                     /* return 'ok' */
}  /* pcc_single() */

/*--------------------------------------------------------------------*/

static int SFXNAME(pcc_multi) (REAL *diff, REAL *rssd, REAL *res,
                               int N, int T, int X,
                               int var, int tile, int nthd)
{                               /* --- multi-thread version */
  int  i, k, r = 0;             /* loop variables, error status */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif
  THREAD *threads;              /* thread handles */
  WORKER *worker;               /* worker for parallel execution */
  SFXNAME(WORK) *w;             /* data   for worker thread */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  DWORD  thid;                  /* dummy for storing the thread id */
  #endif                        /* (not really needed here) */

  assert(diff && rssd && res && (N > 0) && (T > 0) && (nthd > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  /* --- get thread data and worker --- */
  k = (N/2 +nthd-1) /nthd;      /* compute the number of series */
  if (k <= 0) k = 1;            /* to be processed per thread */
  if (k <= tile)                /* if only one tile per thread, */
    var &= ~PCC_TILED;          /* use the untiled version */
  threads = malloc((size_t)nthd *sizeof(THREAD));
  if (!threads) { DELROWS; return -1; }
  w = malloc((size_t)nthd *sizeof(SFXNAME(WORK)) );
  if (!w) { free(threads); DELROWS; return -1; }
  switch (var & PCC_VARIANT) {  /* evaluate the variant */
    #ifdef __AVX__              /* if AVX instructions available */
    case PCC_AVX:               /* if to use AVX computation */
      worker = (var & PCC_COBL)  ? SFXNAME(wrk_avx_cobl)
             : (var & PCC_TILED) ? SFXNAME(wrk_avx_tiled)
                                 : SFXNAME(wrk_avx);   break;
    #endif                      /* get the AVX worker */
    #ifdef __SSE2__             /* if SSE2 instructions available */
    case PCC_SSE2:              /* if to use SSE2 computation */
      worker = (var & PCC_COBL)  ? SFXNAME(wrk_sse2_cobl)
             : (var & PCC_TILED) ? SFXNAME(wrk_sse2_tiled)
                                 : SFXNAME(wrk_sse2);  break;
    #endif                      /* get the SSE2 worker */
    default:                    /* if to use naive computation */
      worker = (var & PCC_COBL)  ? SFXNAME(wrk_naive_cobl)
             : (var & PCC_TILED) ? SFXNAME(wrk_naive_tiled)
                                 : SFXNAME(wrk_naive); break;
  }                             /* get the default worker (naive) */

  /* --- execute the threads --- */
  for (i = 0; i < nthd; i++) {  /* traverse the threads */
    w[i].diff = diff;           /* store differences to mean and */
    w[i].rssd = rssd;           /* roots of sums of squared deviats. */
    #if ROWS                    /* if row starts are available */
    w[i].rows = rows;           /* note the row starts */
    #else                       /* if to use index computation */
    w[i].res  = res;            /* note the result array */
    #endif
    w[i].N    = N;              /* note the number of series, */
    w[i].T    = T;              /* the length of each serie */
    w[i].X    = X;              /* and the padded length, */
    w[i].tile = tile;           /* and the tile size */
    w[i].s    = i*k;            /* compute and store start index */
    if (w[i].s >= N/2) break;   /* if beyond half, already done */
    w[i].e    = w[i].s +k;      /* compute and store end index */
    if (w[i].e >= N/2) w[i].e = N -w[i].s;
    #ifdef _WIN32               /* if Microsoft Windows system */
    threads[i] = CreateThread(NULL, 0, worker, w+i, 0, &thid);
    if (!threads[i]) { r = -1; break; }
    #else                       /* if Linux/Unix system */
    if (pthread_create(threads+i, NULL, worker, w+i) != 0) {
      r = -1; break; }          /* create a thread for each strip */
    #endif                      /* to compute the strips in parallel */
  }
  #ifdef _WIN32                 /* if Microsoft Windows system */
  WaitForMultipleObjects(i, threads, TRUE, INFINITE);
  while (--i >= 0)              /* wait for threads to finish, */
    CloseHandle(threads[i]);    /* then close all thread handles */
  #else                         /* if Linux/Unix system */
  while (--i >= 0)              /* wait for threads to finish */
    pthread_join(threads[i], NULL);
  #endif                        /* (join threads with this one) */

  free(threads);                /* deallocate thread handles */
  free(w);                      /* deallocate parallelization data */
  DELROWS;                      /* deallocate row starts */
  return r;                     /* return error status */
}  /* pcc_multi() */

/*--------------------------------------------------------------------*/

int SFXNAME(pccx) (REAL *data, REAL *res, int N, int T, int var, ...)
{                               /* --- compute Pearson's corr. coeff. */
  int     r;                    /* buffer, return value */
  int     tile = 0;             /* size of the tiles */
  int     nthd = proccnt();     /* get the number of processors */
  int     X;                    /* size of padded data arrays */
  REAL    *rssd;                /* roots of sums of squared deviats. */
  REAL    *diff;                /* (differences to) mean values */
  va_list args;                 /* list of variable arguments */

  assert(data && res && (N > 0) && (T > 0));

  /* --- get variable arguments --- */
  va_start(args, var);          /* start variable arguments */
  if (var & (PCC_TILED|PCC_COBL))  /* if to use a tiled version */
    tile = va_arg(args, int);   /* get the tile size */
  if (var &  PCC_THREAD)        /* if to use a threaded version */
    nthd = va_arg(args, int);   /* get the number of threads */
  va_end(args);                 /* end variable arguments */

  /* --- choose and check variant --- */
  while (!(var & PCC_VARIANT)){ /* if automatic choice of variant */
    #ifdef __AVX__              /* if AVX instructions available */
    if (hasAVX()) {             /* and actually supported by CPU */
      var = PCC_AVX  |PCC_COBL|GET_THREAD(var); break; }
    #endif                      /* otherwise (next best choice) */
    #ifdef __SSE2__             /* if SSE2 instructions available */
    if (hasSSE2()) {            /* and actually supported by CPU */
      var = PCC_SSE2 |PCC_COBL|GET_THREAD(var); break; }
    #endif                      /* otherwise (fallback) */
      var = PCC_NAIVE|PCC_COBL|GET_THREAD(var); break;
  }                             /* use naive computation */
  #ifdef __AVX__                /* if AVX instructions available */
  if (((var & PCC_VARIANT) == PCC_AVX)  && !hasAVX())
  #else                         /* check processor capabilities */
  if  ((var & PCC_VARIANT) == PCC_AVX)
  #endif                        /* otherwise AVX is impossible */
    { fprintf(stderr,  "AVX not supported!\n"); return -1; }
  #ifdef __SSE2__               /* if SSE2 instructions available */
  if (((var & PCC_VARIANT) == PCC_SSE2) && !hasSSE2())
  #else                         /* check processor capabilities */
  if  ((var & PCC_VARIANT) == PCC_SSE2)
  #endif                        /* otherwise SSE2 is impossible */
    { fprintf(stderr, "SSE2 not supported!\n"); return -1; }

  /* --- handle tiling and choose tile size --- */
  if      ((var & PCC_VARIANT) == PCC_SSE2)
    #if REAL_IS_DOUBLE          /* if data is double precision */
    X = (T+1) & ~1;             /* compute padded data array size */
    #else                       /* if data is single precision */
    X = (T+3) & ~3;             /* compute padded data array size */
    #endif                      /* (two or four numbers) */
  else if ((var & PCC_VARIANT) == PCC_AVX)
    #if REAL_IS_DOUBLE          /* if data is double precision */
    X = (T+3) & ~3;             /* compute padded data array size */
    #else                       /* if data is single precision */
    X = (T+7) & ~7;             /* compute padded data array size */
    #endif                      /* (four or eight numbers) */
  else                          /* if naive computation, */
    X = T;                      /* no padding is needed */
  if (var & PCC_COBL) {         /* if to use cache-oblivious tiling */
    var &= ~PCC_TILED;          /* disable standard tiling */
    if (tile <= 0) {            /* if to choose minimal tile size */
      r = (int)((size_t)(16*1024)/sizeof(REAL)/(size_t)X) -1;
      tile = (r < 16) ? 16 : r; /* compute a suitable value */
    }                           /* for the minimal tile size */
  }                             /* assuming 16kB first level cache */
  if ((var & PCC_TILED)         /* if to use tiling and to choose */
  &&  (tile <= 0)) {            /* the tile size automatically */
    #if 1
    r = (int)((size_t)(1024*1024)/sizeof(REAL)/(size_t)X);
    #else
    r = (int)((size_t)(  32*1024)/sizeof(REAL)/(size_t)X);
    #endif
    for (tile = 1; tile < r; tile <<= 1);
    tile = (tile >> 1) -1;      /* find the tile size from the size */
  }                             /* of the arrays and a cache estimate */
  if ((tile <= 1) || (tile >= N))
    var &= ~PCC_TILED;          /* check for a useful tile size */

  /* --- compute differences and rssds --- */
  rssd = malloc(((size_t)N*(size_t)X +(size_t)N) *sizeof(REAL) +31);
  diff = (REAL*)(((ptrdiff_t)(rssd +N) +31) & ~31);
  if (!rssd) return -1;         /* allocate some work memory */
  switch (var & PCC_VARIANT) {  /* evaluate the variant */
    #ifdef __AVX__              /* if AVX instructions available */
    case PCC_AVX:               /* if to use AVX computation */
      SFXNAME(init_avx)  (data, N, T, diff, X, rssd); break;
    #endif                      /* compute differences and rssds */
    #ifdef __SSE2__             /* if SSE2 instructions available */
    case PCC_SSE2:              /* if to use SSE2 computation */
      SFXNAME(init_sse2) (data, N, T, diff, X, rssd); break;
    #endif                      /* compute differences and rssds */
    default:                    /* if to use naive computation */
      SFXNAME(init_naive)(data, N, T, diff, X, rssd); break;
  }                             /* compute differences and rssds */

  /* --- compute correlation coefficients --- */
  if (nthd <= 1)                /* if to use only one thread, */
    var &= ~PCC_THREAD;         /* do not use multi-thread version */
  if (var & PCC_THREAD)         /* if multi-thread  version */
    r = SFXNAME(pcc_multi) (diff, rssd, res, N, T, X, var, tile, nthd);
  else                          /* if single-thread version */
    r = SFXNAME(pcc_single)(diff, rssd, res, N, T, X, var, tile);

  free(rssd);                   /* deallocate work memory */
  return r;                     /* return the error status */
}  /* pccx() */

/*----------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------*/
#if _PCC_PASS == 1              /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "pcc.c"                /* process file recursively */
#endif

/*----------------------------------------------------------------------
  Main Function
----------------------------------------------------------------------*/
#ifdef PCC_MAIN
#undef REAL                     /* delete definitions */
#undef SUFFIX                   /* of precision selectors */
#if 1
#define REAL    float           /* test case: single precision */
#define SUFFIX  _flt            /* function name suffix is '_flt' */
#else
#define REAL    double          /* test case: double precision */
#define SUFFIX  _dbl            /* function name suffix is '_dbl' */
#endif

/*--------------------------------------------------------------------*/

static double timer (void)
{                               /* --- get current time */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  return 0.001 *(double)timeGetTime();
  #else                         /* if Linux/Unix system */
  struct timespec tp;           /* POSIX time specification */
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return (double)tp.tv_sec +1e-9 *(double)tp.tv_nsec;
  #endif                        /* return time in seconds */
}  /* timer() */

/*--------------------------------------------------------------------*/

static void check (REAL *res1, REAL *res2, size_t len)
{                               /* --- compare result arrays */
  size_t k;                     /* loop variable */
  REAL   d, s, m;               /* difference, sum, maximum */

  assert(res1 && res2);         /* check the function arguments */
  for (s = m = 0, k = 0; k < len; k++) {
    d = (REAL)fabs(res1[k]-res2[k]); /* traverse the result arrays */
    s += d; if (d > m) m = d;   /* compute difference and sum them */
  }                             /* and compute their maximum */
  printf("   max: %-13.8g sum: %-13.8g\n", m, s);
}  /* check() */

/*--------------------------------------------------------------------*/

int main (int argc, char* argv[])
{                               /* --- main function for testing */
  int      i, n = 1;            /* loop variable */
  int      N, T;                /* number of voxels and time points */
  size_t   z;                   /* size of the result array */
  REAL     *data, *res1, *res2; /* data and result arrays */
  double   t;                   /* timer for measurements */

  if (argc < 3) {               /* check the function arguments */
    fprintf(stderr, "usage: %s N T [runs]\n", argv[0]); return 0; }
  if (argc > 3) n = atoi(argv[3]);
  N = atoi(argv[1]);            /* get the number of voxels */
  T = atoi(argv[2]);            /* and the number of time points */
  z = (size_t)N*(size_t)(N-1)/2;/* compute size of result arrays */
#if 1
  srand((unsigned)time(NULL));  /* seed the random number generator */
#else
  srand(1);                     /* seed the random number generator */
#endif

  data = malloc((size_t)N *(size_t)T *sizeof(REAL));
  res1 = malloc((z+z) *sizeof(REAL));
  if (!data || !res1) {         /* allocate data and result arrays */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }
  res2 = res1 +z;               /* get second result array */
  for (i = 0; i < T*N; i++)     /* fill data with random numbers */
    data[i] = (REAL)(rand()/((double)RAND_MAX+1));

  /* --- naive computation --- */
  memset(res1, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res1, N, T, PCC_NAIVE);
  printf("naive      : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  printf("   (%d repetition(s))\n", n);

  /* --- naive computation, tiled --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_NAIVE|PCC_TILED, 0);
  printf("naive tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- naive computation, cache-oblivious --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_NAIVE|PCC_COBL, 0);
  printf("naive cobl : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- naive computation, threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_NAIVE|PCC_THREAD, 4);
  printf("naive thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- naive computation, tiled & threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_NAIVE|PCC_TILED|PCC_THREAD,0,4);
  printf("naive t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- naive computation, cache-oblivious & threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_NAIVE|PCC_COBL|PCC_THREAD, 0,4);
  printf("naive cot'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  if (z <= 10)
    for (i = 0; i < (int)z; i++)
      printf("%.16g\n", res1[i]);

#ifdef __SSE2__
  /* --- SSE2 computation --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_SSE2);
  printf("sse2       : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSE2 computation, tiled --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_SSE2|PCC_TILED, 0);
  printf("sse2  tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSE2 computation, cache-oblivious --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_SSE2|PCC_COBL, 0);
  printf("sse2  cobl : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSE2 computation, threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_SSE2|PCC_THREAD, 4);
  printf("sse2  thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSE2 computation, tiled & threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_SSE2|PCC_TILED|PCC_THREAD, 0,4);
  printf("sse2  t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSE2 computation, cache-oblivious & threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_SSE2|PCC_COBL|PCC_THREAD, 0, 4);
  printf("sse2  cot'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  if (z <= 10)
    for (i = 0; i < (int)z; i++)
      printf("%.16g\n", res2[i]);
#endif

#ifdef __AVX__
  /* --- AVX computation --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_AVX);
  printf("avx        : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- AVX computation, tiled --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_AVX|PCC_TILED, 0);
  printf("avx   tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- AVX computation, cache-oblivious --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_AVX|PCC_COBL, 0);
  printf("avx   cobl : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- AVX computation, threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_AVX|PCC_THREAD, 4);
  printf("avx   thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- AVX computation, tiled & threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_AVX|PCC_TILED|PCC_THREAD, 0, 4);
  printf("avx   t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- AVX computation, cache-oblivious & threaded --- */
  memset(res2, 0, z *sizeof(REAL));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(pccx)(data, res2, N, T, PCC_AVX|PCC_COBL|PCC_THREAD, 0, 4);
  printf("avx   cot'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  if (z <= 10)
    for (i = 0; i < (int)z; i++)
      printf("%.16g\n", res2[i]);
#endif

  free(res1);                   /* delete the result arrays */
  free(data);                   /* and the data array */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif  /* #ifdef PCC_MAIN */

#undef PCC_MAIN

/*----------------------------------------------------------------------
  File    : tetracc.c
  Contents: compute pairwise tetrachoric correlation coefficients
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#ifndef _WIN32                  /* if Linux/Unix system */
#define _POSIX_C_SOURCE 200809L /* needed for clock_gettime() */
#endif
#include <stdio.h>
#include <stdlib.h>
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
#include "binarize.h"
#include "tetracc.h"

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
#  define SUFFIX    _flt        /* function name suffix is '_flt' */
#else                           /* if in second pass of two */
#  undef  _TCC_PASS             /* ensure _TCC_PASS is undefined */
#  define _TCC_PASS 2           /* define macro for second pass */
#  define REAL      double      /* second pass: double precision */
#  define SUFFIX    _dbl        /* function name suffix is '_dbl' */
#endif

/*--------------------------------------------------------------------*/
#define float  1                /* to check the definition of REAL */
#define double 2

#if   REAL == float             /* if single precision data */
#undef  REAL_IS_DOUBLE
#define REAL_IS_DOUBLE  0       /* clear indicator for double */
#elif REAL == double            /* if double precision data */
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
  Processing Options
----------------------------------------------------------------------*/
#define LOWER       0           /* compute lower triangular matrix */
                                /* default: upper triangular matrix */
#define ROWS        0           /* use row pointers for tiling */
                                /* default: index computation */
#define SSE2_16     0           /* use reduction to 16 bits in SSE2 */
                                /* default: reduction to 8 bits */
#define USERTHD     0           /* copy user flag for threading */
                                /* default: choose automatically */

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

#if ROWS                        /* if row starts are available */
#define RESULT(i,j,N)   rows[i][j]
#else                           /* if to use index computation */
#define RESULT(i,j,N)   res[INDEX(i,j,N)]
#endif                          /* unify left value for result */

#ifndef GET_THREAD              /* if not yet defined */
#if USERTHD                     /* if to respect user flag */
#define GET_THREAD(var) ((var) & TCC_THREAD)
#else                           /* copy the threading flag */
#define GET_THREAD(var) TCC_THREAD
#endif                          /* otherwise force threading flag */
#endif

#ifndef THREAD                  /* if not yet defined */
#ifdef _WIN32                   /* if Microsoft Windows system */
#define THREAD          HANDLE  /* threads are identified by handles */
#define THREAD_OK       0       /* return value is DWORD */
#define WORKERDEF(n,p)  DWORD WINAPI SFXNAME(n) (LPVOID p)
#else                           /* if Linux/Unix system */
#define THREAD          pthread_t
#define THREAD_OK       NULL    /* return value is void* */
#define WORKERDEF(n,p)  void*        SFXNAME(n) (void* p)
#endif                          /* use the POSIX thread type */
#endif

/*----------------------------------------------------------------------
  Type Definitions
----------------------------------------------------------------------*/
typedef struct {                /* --- thread worker data --- */
  void *bits;                   /* binarized data */
  REAL *cmap;                   /* map from n_11 to correl. coeff. */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* array of row starts */
  #else                         /* if to use index computation */
  REAL *res;                    /* result array */
  #endif
  int  N;                       /* number of series */
  int  X;                       /* length of each binarized serie */
  int  tile;                    /* tile size (columns to group) */
  int  s, e;                    /* index of start and end serie */
} SFXNAME(WORK);                /* (thread worker data) */

#ifndef WORKERTYPE              /* if not yet defined */
#define WORKERTYPE              /* define worker data types */
#ifdef _WIN32                   /* if Microsoft Windows system */
typedef DWORD WINAPI WORKER (LPVOID);
#else                           /* if Linux/Unix system */
typedef void*        WORKER (void*);
#endif                          /* worker for parallel execution */
#endif

/*----------------------------------------------------------------------
  Auxiliary Functions
----------------------------------------------------------------------*/
#ifndef TAB_POPCNT              /* if table is not yet defined */
#define TAB_POPCNT              /* table of bit counts in 16 bits */
unsigned char popcnt[1 << 16];

void init_popcnt (void)
{                               /* --- init. bit count lookup table */
  int i, k;                     /* loop variables */

  if (popcnt[65535]) return;    /* if already initialized, abort */
  for (i = 0; ++i < 65536; ) {  /* traverse all word values */
    #if 1                       /* magic computation (~ as in SSE2) */
    k = i - ((i >> 1) & 0x55555555);
    k = (k & 0x33333333) +((k >> 2)  & 0x33333333);
    k = (((k +(k >> 4)) & 0xf0f0f0f) * 0x01010101) >> 24;
    #else                       /* fallback (foolproof method) */
    int j;                      /* loop variable */
    for (k = 0, j = i; j; j >>= 1)
      k += j & 1;               /* count the bits in the value */
    #endif
    popcnt[i] = (unsigned char)k;
  }                             /* store this number in the table */
}  /* init_popcnt() */

#endif  /* #ifndef TAB_POPCNT */
/*--------------------------------------------------------------------*/
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
/*--------------------------------------------------------------------*/

REAL* SFXNAME(make_cmap) (int T)
{                               /* --- create map from n11 to cc */
  int  i;                       /* loop variable */
  REAL *cmap;                   /* map for cosine function */

  assert(T > 0);                /* check the function argument */
  cmap = malloc(((size_t)T+1) *sizeof(REAL));
  if (!cmap) return NULL;       /* allocate a lookup table and */
  for (i = T; i >= 0; i--)      /* initialize map from n_11 to cc */
    cmap[i] = (REAL)-cos(2*M_PI*(double)i/(double)T);
  return cmap;                  /* return the created map */
}  /* make_cmap() */

/*----------------------------------------------------------------------
  Main Functions
----------------------------------------------------------------------*/

int SFXNAME(tcc_lut16) (uint32_t *bits, REAL *res, int N, int T)
{                               /* --- lookup table to count bits */
  int  i, j;                    /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */

  assert(bits && res && (N > 0) && (T > 0));
  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) return -1;         /* (avoid cosine computations) */
  init_popcnt();                /* initialize the bit count table */

  X = (T+31) >> 5;              /* compute binarized array size */
  for (i = 0; i < N; i++) {     /* traverse the matrix rows */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* count set bits in the conjunction */
      s = pcand_lut16(bits+i*X, bits+j*X, X);
      *res++ = cmap[s];         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  return 0;                     /* return 'ok' */
}  /* tcc_lut16() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_lut16, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j;                /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {   /* traverse row indices */
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* count set bits in the conjunction */
        s = pcand_lut16(bits+i*w->X, bits+j*w->X, w->X);
        w->RESULT(i,j,w->N) = w->cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_lut16() */

/*--------------------------------------------------------------------*/

int SFXNAME(tcc_lut16_tiled) (uint32_t *bits, REAL *res,
                              int N, int T, int tile)
{                               /* --- lookup table to count bits */
  int  i, j, m, e;              /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */
  #if ROWS
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(bits && res && (N > 0) && (T > 0) && (tile > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) { DELROWS; return -1; }
  init_popcnt();                /* initialize the bit count table */

  X = (T+31) >> 5;              /* compute binarized array size */
  for (m = 0; m < N; m += tile) {  /* traverse the stripes/tiles */
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* count set bits in the conjunction */
        s = pcand_lut16(bits+i*X, bits+j*X, X);
        RESULT(i,j,N) = cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  DELROWS;                      /* deallocate the row starts */
  return 0;                     /* return 'ok' */
}  /* tcc_lut16_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_lut16_tiled, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j, m, e;          /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif                    /* count set bits in the conjunction */
          s = pcand_lut16(bits+i*w->X, bits+j*w->X, w->X);
          w->RESULT(i,j,w->N) = w->cmap[s];
        }                       /* compute correlation coefficient */
      }                         /* (Pearson's r) and store it */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_lut16_tiled() */

/*----------------------------------------------------------------------
Reference for the SSE2 popcnt algorithms:
  I.S. Haque, V.S. Pande, and W.P. Walters
  Anatomy of High-Performance 2D Similarity Calculations
  Journal of Chemical Information and Modeling 51(9):2345-351
  American Chemical Society, Washington, DC, USA 2011
  dx.doi.org/10.1021/ci200235e
In the original source the function that processes a block uses 4 times
loop unrolling. This was removed here, since it increases the alignment
requirements to 256 bytes, which is too large for this application.
With this version alignment needs to be only 16 bytes = 4 x uint32_t.
----------------------------------------------------------------------*/
#ifdef __SSE2__                 /* if SSE2 instructions available */
#ifndef PCAND_SSE2              /* if not defined yet */
#define PCAND_SSE2              /* define array bit counting function */
#if SSE2_16                     /* if reduction to 16 bit integers */

static inline __m128i pca_sse2 (uint32_t *a, uint32_t *b, int n)
{                               /* --- popcnt of conj. of a block */
  const __m128i zero  = _mm_setzero_si128();
  const __m128i mask0 = _mm_set1_epi32(0x55555555);
  const __m128i mask1 = _mm_set1_epi32(0x33333333);
  const __m128i mask2 = _mm_set1_epi32(0x0f0f0f0f);
  const __m128i mask3 = _mm_set1_epi32(0x00FF00FF);
  int     i;                    /* loop variable */
  __m128i c;                    /* 8 parallel popcnt sums of 16 bit */

  assert(a && b && (n <= 64) && ((n & 3) == 0));
  c = _mm_setzero_si128();      /* initialize the sums */
  for (i = 0; i < n; i += 4) {  /* traverse the arrays */
    /* compute conjunction of array elements (4 uint32_t at a time) */
    __m128i v = _mm_and_si128(_mm_load_si128((__m128i*)(a+i)),
                              _mm_load_si128((__m128i*)(b+i)));
    v = _mm_sub_epi16(v, _mm_and_si128(mask0, _mm_srli_epi16(v, 1)));
    v = _mm_add_epi16(_mm_and_si128(mask1, v),
                      _mm_and_si128(mask1, _mm_srli_epi16(v, 2)));
    v = _mm_and_si128(mask2, _mm_add_epi16(v, _mm_srli_epi16(v, 4)));
    v = _mm_and_si128(mask3, _mm_add_epi16(v, _mm_srli_epi16(v, 8)));
    c = _mm_add_epi16(c, v);    /* reduce consecutively to pairs, */
  }                             /* then nibbles, then bytes, then words */
  /* reduce {i16,..6x..,i16} to {i32,i32,i32,i32} */
  return _mm_add_epi32(_mm_unpacklo_epi16(c, zero),
                       _mm_unpackhi_epi16(c, zero));
}  /* pca_sse2 */

/*----------------------------------------------------------------------*/

static inline int pcand_sse2 (uint32_t *a, uint32_t *b, int n)
{                               /* --- pop. count of conjunction */
  int     s;                    /* sum of population counts */
  __m128i c;                    /* two popcnt sums of 32 bit each */

  assert(a && b && (n > 0));    /* check the function arguments */
  c = _mm_setzero_si128();      /* clear the two 16bit sums */
  for ( ; n > 2048; n -= 2048){ /* process in blocks of 64 values */
    c = _mm_add_epi32(c, pca_sse2(a, b, 64));
    a += 2048; b += 2048;       /* process a block and */
  }                             /* skip the processed elements */
  if (n > 0)                    /* process the remaining elements */
    c = _mm_add_epi32(c, pca_sse2(a, b, n));
  c = _mm_add_epi32(c, _mm_shuffle_epi32(c, _MM_SHUFFLE(0,1,2,3)));
  c = _mm_add_epi32(c, _mm_shuffle_epi32(c, 2));
  _mm_store_ss((float*)&s, (__m128)c);
  return s;                     /* combine and return sums */
}  /* pcand_sse2() */

/*----------------------------------------------------------------------*/
#else  /* #if SSE2_16 .. */     /* if reduction to 8 bit integers */

static inline __m128i _mm_srli1_epi8 (__m128i v)
{ return _mm_avg_epu8(v, _mm_set1_epi8((char)-1)); }

/*----------------------------------------------------------------------*/

static inline __m128i _mm_srli2_epi8 (__m128i v)
{ return _mm_srli1_epi8(_mm_srli1_epi8(v)); }

/*----------------------------------------------------------------------*/

static inline __m128i pca_sse2 (uint32_t *a, uint32_t *b, int n)
{                               /* --- popcnt of conj. of a block */
  const __m128i zero  = _mm_setzero_si128();
  const __m128i mask0 = _mm_set1_epi32(0x55555555);
  const __m128i mask1 = _mm_set1_epi32(0x33333333);
  const __m128i mask2 = _mm_set1_epi32(0x0f0f0f0f);
  int     i;                    /* loop variable */
  __m128i c;                    /* 8 parallel popcnt sums of 16 bit */

  assert(a && b && (n <= 64) && ((n & 3) == 0));
  c = _mm_setzero_si128();      /* initialize the sums */
  for (i = 0; i < n; i += 4) {  /* traverse the arrays */
    /* compute conjunction of array elements (4 uint32_t at a time) */
    __m128i v = _mm_and_si128(_mm_load_si128((__m128i*)(a+i)),
                              _mm_load_si128((__m128i*)(b+i)));
    v = _mm_sub_epi8(v, _mm_and_si128(mask0, _mm_srli1_epi8(v)));
    v = _mm_add_epi8(_mm_and_si128(mask1, v),
                     _mm_and_si128(mask1,_mm_srli2_epi8(v)));
    v = _mm_and_si128(mask2, _mm_add_epi8(v,
          _mm_srli_epi16(_mm_andnot_si128(mask2, v), 4)));
    c = _mm_add_epi8(c, v);     /* reduce consecutively to pairs, */
  }                             /* then nibbles, then bytes */
  /* reduce {i8,..14x..,i8} to {-,-,-,i16,-,-,-,i16} */
  return _mm_sad_epu8(c, zero); /* combine sums into two and return */
}  /* pca_sse2 */

/*----------------------------------------------------------------------*/

static inline int pcand_sse2 (uint32_t *a, uint32_t *b, int n)
{                               /* --- pop. count of conjunction */
  int     s;                    /* sum of population counts */
  __m128i c;                    /* two popcnt sums of 32 bit each */

  assert(a && b && ((n & 3) == 0));
  c = _mm_setzero_si128();      /* clear the two 16bit sums */
  for ( ; n > 64; n -= 64) {    /* process in blocks of 64 values */
    c = _mm_add_epi32(c, pca_sse2(a, b, 64));
    a += 64; b += 64;           /* process a block and */
  }                             /* skip the processed elements */
  if (n > 0)                    /* process the remaining elements */
    c = _mm_add_epi32(c, pca_sse2(a, b, n));
  c = _mm_add_epi32(c, _mm_shuffle_epi32(c, 2));
  _mm_store_ss((float*)&s, _mm_castsi128_ps(c));
  return s;                     /* combine the two 32 bit sums */
}  /* pcand_sse2() */           /* and store and return the result */

#endif  /* #if SSE2_16 .. else .. */
#endif  /* #ifndef PCAND_SSE2 .. */
/*--------------------------------------------------------------------*/

int SFXNAME(tcc_sse2) (uint32_t *bits, REAL *res, int N, int T)
{                               /* --- SSE2 function to count bits */
  int  i, j;                    /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */

  assert(bits && res && (N > 0) && (T > 0));
  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) return -1;         /* (avoid cosine computations) */

  X = 4 *((T+127) >> 7);        /* compute binarized array size */
  for (i = 0; i < N; i++) {     /* traverse the matrix rows */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* count set bits in the conjunction */
      s = pcand_sse2(bits+i*X, bits+j*X, X);
      *res++ = cmap[s];         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  return 0;                     /* return 'ok' */
}  /* tcc_sse2() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_sse2, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j;                /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {   /* traverse row indices */
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* count set bits in the conjunction */
        s = pcand_sse2(bits+i*w->X, bits+j*w->X, w->X);
        w->RESULT(i,j,w->N) = w->cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_sse2() */

/*--------------------------------------------------------------------*/

int SFXNAME(tcc_sse2_tiled) (uint32_t *bits, REAL *res,
                             int N, int T, int tile)
{                               /* --- SSE2 function to count bits */
  int  i, j, m, e;              /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(bits && res && (N > 0) && (T > 0) && (tile > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) { DELROWS; return -1; }

  X = 4 *((T+127) >> 7);        /* compute binarized array size */
  for (m = 0; m < N; m += tile) {  /* traverse the stripes/tiles */
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* count set bits in the conjunction */
        s = pcand_sse2(bits+i*X, bits+j*X, X);
        RESULT(i,j,N) = cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  DELROWS;                      /* deallocate the row starts */
  return 0;                     /* return 'ok' */
}  /* tcc_sse2_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_sse2_tiled, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j, m, e;          /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif                    /* count set bits in the conjunction */
          s = pcand_sse2(bits+i*w->X, bits+j*w->X, w->X);
          w->RESULT(i,j,w->N) = w->cmap[s];
        }                       /* compute correlation coefficient */
      }                         /* (Pearson's r) and store it */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_sse2_tiled() */

#endif  /* #if __SSE2__ */
/*----------------------------------------------------------------------
Reference for the SSSE3 popcnt algorithm:
  I.S. Haque, V.S. Pande, and W.P. Walters
  Anatomy of High-Performance 2D Similarity Calculations
  Journal of Chemical Information and Modeling 51(9):2345-351
  American Chemical Society, Washington, DC, USA 2011
  dx.doi.org/10.1021/ci200235e
In the original source the function that processes a block uses 4 times
loop unrolling. This was removed here, since it increases the alignment
requirements to 256 bytes, which is too large for this application.
With this version alignment needs to be only 16 bytes = 4 x uint32_t.
----------------------------------------------------------------------*/
#if __SSSE3__                   /* if SSSE3 instructions available */
#ifndef PCAND_SSSE3             /* if not defined yet */
#define PCAND_SSSE3             /* define array bit counting function */

static __m128i pca_ssse3 (uint32_t *a, uint32_t *b, int n)
{                               /* --- popcnt of conj. of a block */
  const __m128i zero = _mm_setzero_si128();
  const __m128i mask = _mm_set1_epi32(0x0f0f0f0f);
  const __m128i lut4 = _mm_set_epi32 (0x04030302, 0x03020201,
                                      0x03020201, 0x02010100);
  /* A look-up table of population counts in each possible nibble,    */
  /* from least to most significant: 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4. */
  int     i;                    /* loop variable */
  __m128i c;                    /* 16 parallel popcnt sums of 8 bit */

  assert(a && b && (n <= 64) && ((n & 3) == 0));
  c = _mm_setzero_si128();      /* initialize the sums */
  for (i = 0; i < n; i += 4) {  /* traverse the arrays */
    /* compute conjunction of array elements (4 uint32_t at a time) */
    __m128i v = _mm_and_si128(_mm_load_si128((__m128i*)(a+i)),
                              _mm_load_si128((__m128i*)(b+i)));
    /* split each byte into low and high nibble */
    __m128i l = _mm_and_si128(mask, v);
    __m128i h = _mm_and_si128(mask, _mm_srli_epi16(v, 4));
    /* compute popcnt of each byte in two halves with pshufb and lut4 */
    c = _mm_add_epi8(c, _mm_shuffle_epi8(lut4, l));
    c = _mm_add_epi8(c, _mm_shuffle_epi8(lut4, h));
  }                             /* sum the 16 population counts */
  /* reduce {i8,..14x..,i8} -> {-,-,-,i16,-,-,-,i16} */
  return _mm_sad_epu8(c, zero); /* combine sums into two and return */
}  /* pca_ssse3() */

/*--------------------------------------------------------------------*/

static int pcand_ssse3 (uint32_t *a, uint32_t *b, int n)
{                               /* --- pop. count of conjunction */
  int     r;                    /* to store the result */
  __m128i c;                    /* two popcnt sums of 32 bit each */

  assert(a && b && ((n & 3) == 0));
  c = _mm_setzero_si128();      /* clear the two 16bit sums */
  for ( ; n > 64; n -= 64) {    /* process in blocks of 64 values */
    c = _mm_add_epi32(c, pca_ssse3(a, b, 64));
    a += 64; b += 64;           /* process a block and */
  }                             /* skip the processed elements */
  if (n > 0)                    /* process the remaining elements */
    c = _mm_add_epi32(c, pca_ssse3(a, b, n));
  c = _mm_add_epi32(c, _mm_shuffle_epi32(c, 2));
  _mm_store_ss((float*)&r, _mm_castsi128_ps(c));
  return r;                     /* combine and return sums */
}  /* pcand_ssse3() */

#endif  /* #ifndef PCAND_SSSE3 */
/*--------------------------------------------------------------------*/

int SFXNAME(tcc_ssse3) (uint32_t *bits, REAL *res, int N, int T)
{                               /* --- SSSE3 function to count bits */
  int  i, j;                    /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */

  assert(bits && res && (N > 0) && (T > 0));
  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) return -1;         /* (avoid cosine computations) */

  X = 4 *((T+127) >> 7);        /* compute binarized array size */
  for (i = 0; i < N; i++) {     /* traverse the matrix rows */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* count set bits in the conjunction */
      s = pcand_ssse3(bits+i*X, bits+j*X, X);
      *res++ = cmap[s];         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  return 0;                     /* return 'ok' */
}  /* tcc_ssse3() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_ssse3, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j;                /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function arguments */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {   /* traverse row indices */
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* count set bits in the conjunction */
        s = pcand_ssse3(bits+i*w->X, bits+j*w->X, w->X);
        w->RESULT(i,j,w->N) = w->cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_ssse3() */

/*--------------------------------------------------------------------*/

int SFXNAME(tcc_ssse3_tiled) (uint32_t *bits, REAL *res,
                              int N, int T, int tile)
{                               /* --- SSSE3 function to count bits */
  int  i, j, m, e;              /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(bits && res && (N > 0) && (T > 0) && (tile > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) { DELROWS; return -1; }

  X = 4 *((T+127) >> 7);        /* compute binarized array size */
  for (m = 0; m < N; m += tile) {  /* traverse the stripes/tiles */
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* count set bits in the conjunction */
        s = pcand_ssse3(bits+i*X, bits+j*X, X);
        RESULT(i,j,N) = cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  DELROWS;                      /* deallocate the row starts */
  return 0;                     /* return 'ok' */
}  /* tcc_ssse3_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_ssse3_tiled, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j, m, e;          /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif                    /* count set bits in the conjunction */
          s = pcand_ssse3(bits+i*w->X, bits+j*w->X, w->X);
          w->RESULT(i,j,w->N) = w->cmap[s];
        }                       /* compute correlation coefficient */
      }                         /* (Pearson's r) and store it */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_ssse3_tiled() */

#endif  /* #if __SSSE3__ */
/*--------------------------------------------------------------------*/
#if __POPCNT__                  /* if popcnt instructions available */
#ifndef PCAND_POP32             /* if not defined yet */
#define PCAND_POP32             /* define array bit counting function */

static inline int pcand_pop32 (uint32_t *a, uint32_t *b, int n)
{                               /* --- pop. count of conjunction */
  int k, s;                     /* loop variable, population count */

  assert(a && b && (n > 0));    /* check the function arguments */
  for (k = s = 0; k < n; k++)   /* traverse conjunction of the series */
    s += (int)_mm_popcnt_u32(a[k] & b[k]);
  return s;                     /* return the population count */
}  /* pcand_pop32() */

#endif  /* ifndef PCAND_POP32 */
/*--------------------------------------------------------------------*/

int SFXNAME(tcc_pop32) (uint32_t *bits, REAL *res, int N, int T)
{                               /* --- _mm_popcnt_u32() to count bits */
  int  i, j;                    /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */

  assert(bits && res && (N > 0) && (T > 0));
  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) return -1;         /* (avoid cosine computations) */

  X = (T+31) >> 5;              /* compute binarized array size */
  for (i = 0; i < N; i++) {     /* traverse the matrix rows */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* count set bits in the conjunction */
      s = pcand_pop32(bits+i*X, bits+j*X, X);
      *res++ = cmap[s];         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  return 0;                     /* return 'ok' */
}  /* tcc_pop32() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_pop32, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j;                /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {   /* traverse row indices */
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* count set bits in the conjunction */
        s = pcand_pop32(bits+i*w->X, bits+j*w->X, w->X);
        w->RESULT(i,j,w->N) = w->cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_pop32() */

/*--------------------------------------------------------------------*/

int SFXNAME(tcc_pop32_tiled) (uint32_t *bits, REAL *res,
                              int N, int T, int tile)
{                               /* --- _mm_popcnt_u32() to count bits */
  int  i, j, m, e;              /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(bits && res && (N > 0) && (T > 0) && (tile > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) { DELROWS; return -1; }

  X = (T+31) >> 5;              /* compute binarized array size */
  for (m = 0; m < N; m += tile) {  /* traverse the stripes/tiles */
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* count set bits in the conjunction */
        s = pcand_pop32(bits+i*X, bits+j*X, X);
        RESULT(i,j,N) = cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  DELROWS;                      /* deallocate the row starts */
  return 0;                     /* return 'ok' */
}  /* tcc_pop32_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_pop32_tiled, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j, m, e;          /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function arguments */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif                    /* count set bits in the conjunction */
          s = pcand_pop32(bits+i*w->X, bits+j*w->X, w->X);
          w->RESULT(i,j,w->N) = w->cmap[s];
        }                       /* compute correlation coefficient */
      }                         /* (Pearson's r) and store it */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_pop32_tiled() */

#endif  /* #if __POPCNT__ */
/*--------------------------------------------------------------------*/
#if __POPCNT__                  /* if popcnt instructions available */
#ifndef PCAND_POP64             /* if not defined yet */
#define PCAND_POP64             /* define array bit counting function */

static inline int pcand_pop64 (uint64_t *a, uint64_t *b, int n)
{                               /* --- pop. count of conjunction */
  int k, s;                     /* loop variable, population count */

  assert(a && b && (n > 0));    /* check the function arguments */
  for (k = s = 0; k < n; k++)   /* traverse conjunction of the series */
    s += (int)_mm_popcnt_u64(a[k] & b[k]);
  return s;                     /* return the population count */
}  /* pcand_pop64() */

#endif  /* #ifndef PCAND_POP64 */
/*--------------------------------------------------------------------*/

int SFXNAME(tcc_pop64) (uint64_t *bits, REAL *res, int N, int T)
{                               /* --- _mm_popcnt_u64() to count bits */
  int  i, j;                    /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */

  assert(bits && res && (N > 0) && (T > 0));
  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) return -1;         /* (avoid cosine computations) */

  X = (T+63) >> 6;              /* compute binarized array size */
  for (i = 0; i < N; i++) {     /* traverse the matrix rows */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* count set bits in the conjunction */
      s = pcand_pop64(bits+i*X, bits+j*X, X);
      *res++ = cmap[s];         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  return 0;                     /* return 'ok' */
}  /* tcc_pop64() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_pop64, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j;                /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint64_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {   /* traverse row indices */
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* count set bits in the conjunction */
        s = pcand_pop64(bits+i*w->X, bits+j*w->X, w->X);
        w->RESULT(i,j,w->N) = w->cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_pop64() */

/*--------------------------------------------------------------------*/

int SFXNAME(tcc_pop64_tiled) (uint64_t *bits, REAL *res,
                              int N, int T, int tile)
{                               /* --- _mm_popcnt_u64() to count bits */
  int  i, j, m, e;              /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(bits && res && (N > 0) && (T > 0) && (tile > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) { DELROWS; return -1; }

  X = (T+63) >> 6;              /* compute binarized array size */
  for (m = 0; m < N; m += tile) {  /* traverse the stripes/tiles */
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* count set bits in the conjunction */
        s = pcand_pop64(bits+i*X, bits+j*X, X);
        RESULT(i,j,N) = cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  DELROWS;                      /* deallocate the row starts */
  return 0;                     /* return 'ok' */
}  /* tcc_pop64_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_pop64_tiled, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j, m, e;          /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint64_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif                    /* count set bits in the conjunction */
          s = pcand_pop64(bits+i*w->X, bits+j*w->X, w->X);
          w->RESULT(i,j,w->N) = w->cmap[s];
        }                       /* compute correlation coefficient */
      }                         /* (Pearson's r) and store it */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_pop64_tiled() */

#endif  /* #ifdef __POPCNT__ */
/*--------------------------------------------------------------------*/
#if defined __POPCNT__ && defined __SSE4_1__

int SFXNAME(tcc_m128i) (uint32_t *bits, REAL *res, int N, int T)
{                               /* --- _m128i integers and popcnt64 */
  int  i, j;                    /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */

  assert(bits && res && (N > 0) && (T > 0));
  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) return -1;         /* (avoid cosine computations) */

  X = 4 *((T+127) >> 7);        /* compute binarized array size */
  for (i = 0; i < N; i++) {     /* traverse the matrix rows */
    #if LOWER                   /* if lower triangular matrix */
    for (j = 0;   j < i; j++) { /* traverse the smaller indices */
    #else                       /* if upper triangular matrix */
    for (j = i+1; j < N; j++) { /* traverse the greater indices */
    #endif                      /* count set bits in the conjunction */
      s = pcand_m128i(bits+i*X, bits+j*X, X);
      *res++ = cmap[s];         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  return 0;                     /* return 'ok' */
}  /* tcc_m128i() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_m128i, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j;                /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (i = w->s; i < w->e; i++) {   /* traverse row indices */
      #if LOWER                 /* if lower triangular matrix */
      for (j = 0;   j < i;    j++) {
      #else                     /* if upper triangular matrix */
      for (j = i+1; j < w->N; j++) {
      #endif                    /* count set bits in the conjunction */
        s = pcand_m128i(bits+i*w->X, bits+j*w->X, w->X);
        w->RESULT(i,j,w->N) = w->cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_m128i() */

/*--------------------------------------------------------------------*/

int SFXNAME(tcc_m128i_tiled) (uint32_t *bits, REAL *res,
                              int N, int T, int tile)
{                               /* --- _m128i integers and popcnt64 */
  int  i, j, m, e;              /* loop variables */
  int  X;                       /* size of padded data arrays */
  int  s;                       /* number of set bits in conjunction */
  REAL *cmap;                   /* map for cosine function */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif

  assert(bits && res && (N > 0) && (T > 0) && (tile > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */

  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) { DELROWS; return -1; }

  X = 4 *((T+127) >> 7);        /* compute binarized array size */
  for (m = 0; m < N; m += tile) {  /* traverse the stripes/tiles */
    e = (m+tile < N) ? m+tile : N; /* get end column of stripe/tile */
    #if LOWER                   /* if lower triangular matrix */
    for (i = m+1; i < N; i++) { /* traverse the greater indices */
      for (j = (i < e) ? i : e; --j >= m; ) {
    #else                       /* if upper triangular matrix */
    for (i = 0; i < e; i++) {   /* traverse the smaller indices */
      for (j = (i >= m) ? i+1 : m; j < e; j++) {
    #endif                      /* count set bits in the conjunction */
        s = pcand_m128i(bits+i*X, bits+j*X, X);
        RESULT(i,j,N) = cmap[s];
      }                         /* map the resulting n_11 value */
    }                           /* with the tabulated cosine function */
  }                             /* to the correlation coefficient */

  free(cmap);                   /* deallocate the n_11 to cc map */
  DELROWS;                      /* deallocate the row starts */
  return 0;                     /* return 'ok' */
}  /* tcc_m128i_tiled() */

/*--------------------------------------------------------------------*/

static WORKERDEF(wrk_m128i_tiled, p)
{                               /* --- compute tetracc of one strip */
  SFXNAME(WORK) *w = p;         /* type the argument pointer */
  int      i, j, m, e;          /* loop variables */
  int      s;                   /* number of set bits in conjunction */
  uint32_t *bits;               /* binarized data */

  assert(p);                    /* check the function argument */
  bits = w->bits;               /* type the binarized data */
  while (1) {                   /* process two strip parts */
    for (m = w->s; m < w->e; m += w->tile) {
      e = (m +w->tile < w->e) ? m +w->tile : w->e;
      #if LOWER                 /* if lower triangular matrix */
      for (i = m+1; i < w->N; i++) { /* traverse greater indices */
        for (j = (i < e) ? i : e; --j >= m; ) {
      #else                     /* if upper triangular matrix */
      for (i = 0; i < e; i++) { /* traverse the smaller indices */
        for (j = (i >= m) ? i+1 : m; j < e; j++) {
      #endif                    /* count set bits in the conjunction */
          s = pcand_m128i(bits+i*w->X, bits+j*w->X, w->X);
          w->RESULT(i,j,w->N) = w->cmap[s];
        }                       /* compute correlation coefficient */
      }                         /* (Pearson's r) and store it */
    }
    if (w->s > w->N/2) break;   /* if second strip done, abort */
    i    = w->N -w->e;          /* get start of opposite stripe */
    if (w->e > i)      break;   /* if no opposite strip, abort */
    w->e = w->N -w->s;          /* get start and end index */
    w->s = i;                   /* of the opposite strip */
  }
  return THREAD_OK;             /* return a dummy result */
}  /* wrk_m128i_tiled() */

#endif  /* #if defined __POPCNT__ && defined __SSE4_1__ */
/*--------------------------------------------------------------------*/

static int SFXNAME(tcc_thread) (void *bits, REAL *res, int N, int T,
                                int var, int tile, int nthd)
{                               /* --- parallel version */
  int  i, k, r = 0;             /* loop variables, error status */
  int  X;                       /* size of padded data arrays */
  REAL *cmap;                   /* map for cosine function */
  #if ROWS                      /* if to use an array of row starts */
  REAL **rows;                  /* starts of output rows */
  #endif
  THREAD *threads;              /* thread handles */
  WORKER *worker;               /* worker for parallel execution */
  SFXNAME(WORK) *w;             /* data   for worker thread */
  #ifdef _WIN32                 /* if Microsoft Windows system */
  DWORD  thid;                  /* dummy for storing the thread id */
  #endif                        /* (not really needed here) */

  assert(bits && res && (N > 0) && (T > 0) && (nthd > 0));
  #if ROWS                      /* if to use an array of row starts */
  rows = SFXNAME(make_rows)(res, N);
  if (!rows) return -1;         /* allocate and initialize */
  #endif                        /* an array of row starts */
  cmap = SFXNAME(make_cmap)(T); /* create map from n_11 to cc */
  if (!cmap) { DELROWS; return -1; }

  /* get thread data and worker */
  k = (N/2 +nthd-1) /nthd;      /* compute the number of series */
  if (k <= 0) k = 1;            /* to be processed per thread */
  if (k <= tile)                /* if only one tile per thread, */
    var &= ~TCC_TILED;          /* use the untiled version */
  threads = malloc((size_t)nthd *sizeof(THREAD));
  if (!threads) { free(cmap); DELROWS; return -1; }
  w = malloc((size_t)nthd *sizeof(SFXNAME(WORK)));
  if (!w) { free(threads); free(cmap); DELROWS; return -1; }
  switch (var & TCC_VARIANT) {  /* evaluate the variant */
    #ifdef __SSE2__             /* if SSE2 instructions available */
    case TCC_SSE2:              /* if to use SSE2 computation */
      X = 4 *((T+127) >> 7);    /* compute binarized array size */
      worker = (var & TCC_TILED) ? SFXNAME(wrk_sse2_tiled)
                                 : SFXNAME(wrk_sse2);  break;
    #endif                      /* get the SSE2 worker */
    #ifdef __SSSE3__            /* if SSSE3 instructions available */
    case TCC_SSSE3:             /* if to use SSSE3 computation */
      X = 4 *((T+127) >> 7);    /* compute binarized array size */
      worker = (var & TCC_TILED) ? SFXNAME(wrk_ssse3_tiled)
                                 : SFXNAME(wrk_ssse3); break;
    #endif                      /* get the SSSE3 worker */
    #ifdef __POPCNT__           /* if popcnt instructions available */
    case TCC_POP32:             /* if to use _mm_popcnt_u32() */
      X = (T+31) >> 5;          /* compute binarized array size */
      worker = (var & TCC_TILED) ? SFXNAME(wrk_pop32_tiled)
                                 : SFXNAME(wrk_pop32); break;
    case TCC_POP64:             /* if to use _mm_popcnt_u32() */
      X = (T+63) >> 6;          /* compute binarized array size */
      worker = (var & TCC_TILED) ? SFXNAME(wrk_pop64_tiled)
                                 : SFXNAME(wrk_pop64); break;
    #endif                      /* get the popcnt32 worker */
    #if defined __POPCNT__ && defined __SSE4_1__
    case TCC_M128I:             /* if to use 128 bit integers */
      X = 4 *((T+127) >> 7);    /* compute binarized array size */
      worker = (var & TCC_TILED) ? SFXNAME(wrk_m128i_tiled)
                                 : SFXNAME(wrk_m128i); break;
    #endif                      /* get the 128 bit integer worker */
    default:                    /* if to use a 16 bit lookup table */
      X = (T+31) >> 5;          /* compute binarized array size */
      init_popcnt();            /* initialize the bit count table */
      worker = (var & TCC_TILED) ? SFXNAME(wrk_lut16_tiled)
                                 : SFXNAME(wrk_lut16); break;
  }                             /* get the default worker (lut16) */

  /* create the threads */
  for (i = 0; i < nthd; i++) {  /* traverse the threads */
    w[i].bits = bits;           /* store the binarized data and */
    w[i].cmap = cmap;           /* the map from n_11 to cc */
    #if ROWS                    /* if row starts are available */
    w[i].rows = rows;           /* note the row starts */
    #else                       /* if to use index computation */
    w[i].res  = res;            /* note the result array */
    #endif
    w[i].N    = N;              /* note the number of series, */
    w[i].X    = X;              /* the length of each serie, */
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
  free(cmap);                   /* deallocate n_11 to cc map */
  DELROWS;                      /* deallocate row starts */
  return r;                     /* return error status */
}  /* tcc_thread() */

/*--------------------------------------------------------------------*/

int SFXNAME(tetraccx) (REAL *data, REAL *res, int N, int T, int var,...)
{                               /* --- compute tetrachoric corr. cf. */
  va_list args;                 /* list of variable arguments */
  int     r;                    /* return value, buffer */
  int     tile = 0;             /* size of the tiles */
  int     nthd = proccnt();     /* get the number of processors */
  void    *bits;                /* binarized data */
  static int bpi[] = {  32,     /* TCC_AUTO  */
                        32,     /* TCC_LUT16 */
                       128,     /* TCC_SSE2  */
                       128,     /* TCC_SSSE3 */
                        32,     /* TCC_POP32 */
                        64,     /* TCC_POP64 */
                       128,     /* TCC_M182I */
                        32 };   /* <reserved> */

  assert(data && res && (N > 0) && (T > 0));
  va_start(args, var);          /* start variable arguments */
  if (var & TCC_TILED)          /* if to use a tiled version */
    tile = va_arg(args, int);   /* get the tile size */
  if (var & TCC_THREAD)         /* if to use a threaded version */
    nthd = va_arg(args, int);   /* get the number of threads */
  va_end(args);                 /* end variable arguments */

  while (!(var & TCC_VARIANT)){ /* if automatic choice of variant */
    #ifdef __POPCNT__           /* if popcnt instruction available */
    if (hasPOPCNT()) {          /* and actually supported by CPU */
      var = TCC_POP64|GET_THREAD(var); break; }
    #endif                      /* otherwise (next best choice) */
    #ifdef __SSSE3__            /* if SSSE3 instructions available */
    if (hasSSSE3()) {           /* and actually supported by CPU */
      var = TCC_SSSE3|GET_THREAD(var); break; }
    #endif                      /* otherwise (fallback) */
      var = TCC_LUT16|GET_THREAD(var); break;
  }                             /* use 16 bit lookup table */
  #ifdef __SSE2__               /* if SSE2 instructions available */
  if (((var & TCC_VARIANT) == TCC_SSE2)  && !hasSSE2())
  #else                         /* check processor capabilities */
  if  ((var & TCC_VARIANT) == TCC_SSE2)
  #endif                        /* otherwise SSE2 is impossible */
    { fprintf(stderr,   "SSE2 not supported!\n"); return -1; }
  #ifdef __SSSE3__              /* if SSSE3 instructions available */
  if (((var & TCC_VARIANT) == TCC_SSSE3) && !hasSSSE3())
  #else                         /* check processor capabilities */
  if  ((var & TCC_VARIANT) == TCC_SSSE3)
  #endif                        /* otherwise SSSE3 is impossible */
    { fprintf(stderr,  "SSSE3 not supported!\n"); return -1; }
  #ifdef __POPCNT__             /* if POPCNT instructions available */
  if (((var & TCC_VARIANT) == TCC_POP32) && !hasPOPCNT())
  #else                         /* check processor capabilities */
  if  ((var & TCC_VARIANT) == TCC_POP32)
  #endif                        /* otherwise POP32 is impossible */
    { fprintf(stderr, "POPCNT not supported!\n"); return -1; }
  #ifdef __POPCNT__             /* if POPCNT instructions available */
  if (((var & TCC_VARIANT) == TCC_POP64) && !hasPOPCNT())
  #else                         /* check processor capabilities */
  if  ((var & TCC_VARIANT) == TCC_POP64)
  #endif                        /* otherwise POP64 is impossible */
    { fprintf(stderr, "POPCNT not supported!\n"); return -1; }
  #if defined __POPCNT__ && defined __SSE4_1__
  if (((var & TCC_VARIANT) == TCC_M128I)/* if POPCNT and SSE4.1 */
  &&  (!hasSSE41() || !hasPOPCNT()))  /* instructions available */
  #else                         /* check processor capabilities */
  if  ((var & TCC_VARIANT) == TCC_M128I)
  #endif                        /* otherwise M128I is impossible */
    { fprintf(stderr, "POPCNT or SSE4.1 not supported!\n"); return -1; }

  if ((var & TCC_TILED)         /* if to use tiling and to choose */
  &&  (tile <= 0)) {            /* the tile size automatically */
    r = (int)((size_t)(8*1024*1024)/sizeof(REAL)/(size_t)T);
    for (tile = 1; tile < r; tile <<= 1);
    tile = (tile >> 1) -1;      /* find the tile size from the size */
  }                             /* of the arrays and a cache estimate */
  if ((tile <= 1) || (tile >= N))
    var &= ~TCC_TILED;          /* check for a useful tile size */
  if (nthd <= 1)                /* if to use only one thread, */
    var &= ~TCC_THREAD;         /* do not use the threaded version */
  bits = SFXNAME(binarize)(data, N, T, BIN_MEDIAN,
                           bpi[var & TCC_VARIANT]);
  if (!bits) return -1;         /* binarize the data w.r.t. median */
  if (var & TCC_THREAD) {       /* if to use multiple threads */
    r = SFXNAME(tcc_thread)(bits, res, N, T, var, tile, nthd);
    free(bits); return r;       /* run thread to compute correlation */
  }                             /* free memory, return error status */

  switch (var & TCC_VARIANT) {  /* evaluate the variant */
    #ifdef __SSE2__             /* if SSE2 instructions avaliable */
    case TCC_SSE2:              /* if to use SSE2 computation */
      r = (var & TCC_TILED)     /* distinguish tiled version */
        ? SFXNAME(tcc_sse2_tiled)(bits, res, N, T, tile)
        : SFXNAME(tcc_sse2)      (bits, res, N, T);
      break;                    /* compute correlation coefficient */
    #endif
    #ifdef __SSSE3__            /* if SSSE3 instructions available */
    case TCC_SSSE3:             /* if to use SSSE3 computation */
      r = (var & TCC_TILED)     /* distinguished tiled version */
        ? SFXNAME(tcc_ssse3_tiled)(bits, res, N, T, tile)
        : SFXNAME(tcc_ssse3)      (bits, res, N, T);
      break;                    /* compute correlation coefficient */
    #endif
    #ifdef __POPCNT__           /* if popcnt instructions available */
    case TCC_POP32:             /* if to use _mm_popcnt_u32() */
      r = (var & TCC_TILED)     /* distinguished tiled version */
        ? SFXNAME(tcc_pop32_tiled)(bits, res, N, T, tile)
        : SFXNAME(tcc_pop32)      (bits, res, N, T);
      break;                    /* compute correlation coefficient */
    case TCC_POP64:             /* if to use use _mm_popcnt_u64() */
      r = (var & TCC_TILED)     /* distinguished tiled version */
        ? SFXNAME(tcc_pop64_tiled)(bits, res, N, T, tile)
        : SFXNAME(tcc_pop64)      (bits, res, N, T);
      break;                    /* compute correlation coefficient */
    #endif
    #if defined __POPCNT__ && defined __SSE4_1__
    case TCC_M128I:             /* if to use 128 bit integers */
      r = (var & TCC_TILED)     /* distinguished tiled version */
        ? SFXNAME(tcc_m128i_tiled)(bits, res, N, T, tile)
        : SFXNAME(tcc_m128i)      (bits, res, N, T);
      break;                    /* compute correlation coefficient */
    #endif
    default:                    /* lookup table for 16 bit integers */
      r = (var & TCC_TILED)     /* distinguished tiled version */
        ? SFXNAME(tcc_lut16_tiled)(bits, res, N, T, tile)
        : SFXNAME(tcc_lut16)      (bits, res, N, T);
      break;                    /* compute correlation coefficient */
  }  /* switch (var & TCC_VARIANT) */

  free(bits);                   /* deallocate work memory */
  return r;                     /* return the error status */
}  /* tetraccx() */

/*----------------------------------------------------------------------
  Recursion Handling
----------------------------------------------------------------------*/
#if _TCC_PASS == 1              /* if in first of two passes */
#undef REAL
#undef SUFFIX
#include "tetracc.c"            /* process header recursively */
#endif

/*----------------------------------------------------------------------
  Main Function
----------------------------------------------------------------------*/
#ifdef TCC_MAIN
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

static void check (float *res1, float *res2, size_t len)
{                               /* --- compare result arrays */
  size_t k;                     /* loop variable */
  float  d, s, m;               /* difference, sum, maximum */

  assert(res1 && res2);         /* check the function arguments */
  for (s = m = 0, k = 0; k < len; k++) {
    d = (float)fabs(res1[k]-res2[k]);/* traverse the result arrays */
    s += d; if (d > m) m = d;   /* compute difference and sum them */
  }                             /* and compute their maximum */
  printf("   max: %-13.8g sum: %-13.8g\n", m, s);
}  /* check() */

/*--------------------------------------------------------------------*/

int main (int argc, char* argv[])
{                               /* --- main function for testing */
  int    i, n = 1;              /* loop variable */
  int    N, T;                  /* number of voxels and time points */
  size_t z;                     /* size of the result array */
  float  *data, *res1, *res2;   /* data and result arrays */
  double t;                     /* timer for measurements */

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

  data = malloc((size_t)N *(size_t)T *sizeof(float));
  res1 = malloc((z+z) *sizeof(float));
  if (!data || !res1) {         /* allocate data and result arrays */
    fprintf(stderr, "%s: not enough memory\n", argv[0]); return -1; }
  res2 = res1 +z;               /* get second result array */
  for (i = 0; i < T*N; i++)     /* fill data with random numbers */
    data[i] = (float)(rand()/((double)RAND_MAX+1));

  /* --- 16 bit lookup table -- */
  memset(res1, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res1, N, T, TCC_LUT16);
  printf("lut16      : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  printf("   (%d repetition(s))\n", n);

  /* --- 16 bit lookup table, tiled -- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_LUT16|TCC_TILED, 0);
  printf("lut16 tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- 16 bit lookup table, threaded -- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_LUT16|TCC_THREAD, 4);
  printf("lut16 thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- 16 bit lookup table, tiled & threaded -- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T,
                      TCC_LUT16|TCC_TILED|TCC_THREAD, 0, 4);
  printf("lut16 t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

#ifdef __SSE2__
  /* --- SSE2 instructions --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_SSE2);
  printf("sse2       : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSE2 instructions, tiled --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_SSE2|TCC_TILED, 0);
  printf("sse2  tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSE2 instructions, threaded --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_SSE2|TCC_THREAD, 4);
  printf("sse2  thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSE2 instructions, tiled & threaded --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T,
                      TCC_SSE2|TCC_TILED|TCC_THREAD, 0, 4);
  printf("sse2  t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);
#endif

#ifdef __SSSE3__
  /* --- SSSE3 instructions --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_SSSE3);
  printf("ssse3      : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSSE3 instructions, tiled --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_SSSE3|TCC_TILED, 0);
  printf("ssse3 tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSSE3 instructions, threaded --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_SSSE3|TCC_THREAD, 4);
  printf("ssse3 thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- SSSE3 instructions, tiled --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T,
                      TCC_SSSE3|TCC_TILED|TCC_THREAD, 0, 4);
  printf("ssse3 t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);
#endif

#ifdef __POPCNT__
  /* --- _mm_popcnt_u32() --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_POP32);
  printf("pop32      : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- _mm_popcnt_u32(), tiled --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_POP32|TCC_TILED, 0);
  printf("pop32 tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- _mm_popcnt_u32(), threaded --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_POP32|TCC_THREAD, 4);
  printf("pop32 thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- _mm_popcnt_u32(), tiled & threaded --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T,
                      TCC_POP32|TCC_TILED|TCC_THREAD, 0, 4);
  printf("pop32 t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);
#endif

#ifdef __POPCNT__
  /* --- _mm_popcnt_u64() --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_POP64);
  printf("pop64      : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- _mm_popcnt_u64(), tiled --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_POP64|TCC_TILED, 0);
  printf("pop64 tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- _mm_popcnt_u64(), threaded --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_POP64|TCC_THREAD, 4);
  printf("pop64 thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- _mm_popcnt_u64(), tiled & threaded --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T,
                      TCC_POP64|TCC_TILED|TCC_THREAD, 0, 4);
  printf("pop64 t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);
#endif

#if defined __POPCNT__ && defined __SSE4_1__
  /* --- __m128i and _mm_popcnt_u64() --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_M128I);
  printf("m128i      : %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- __m128i and _mm_popcnt_u64(), tiled --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_M128I|TCC_TILED, 0);
  printf("m128i tiled: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- __m128i and _mm_popcnt_u64(), thread --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T, TCC_M128I|TCC_THREAD, 4);
  printf("m128i thd'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);

  /* --- __m128i and _mm_popcnt_u64(), tiled & threaded --- */
  memset(res2, 0, z *sizeof(float));
  t = timer();
  for (i = 0; i < n; i++)
    SFXNAME(tetraccx)(data, res2, N, T,
                      TCC_M128I|TCC_TILED|TCC_THREAD, 0, 4);
  printf("m128i t&t'd: %7d %7d %8.2f", N, T, (timer()-t)/(double)n);
  check(res1, res2, z);
#endif

  free(res1);                   /* delete the result arrays */
  free(data);                   /* and the data array */
  return 0;                     /* return 'ok' */
}  /* main() */

#endif  /* #ifdef TCC_MAIN */

#undef TCC_MAIN

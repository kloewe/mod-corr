/*----------------------------------------------------------------------
  File    : binarize.h
  Contents: functions to binarize real-valued series
  Author  : Kristian Loewe, Christian Borgelt
----------------------------------------------------------------------*/
#ifndef __BINARIZE__
#define __BINARIZE__
#include <stddef.h>

/*----------------------------------------------------------------------
  Preprocessor Definitions
----------------------------------------------------------------------*/
#define BIN_MEDIAN  ((void*)1)  /* >= median is set to 1 */
#define BIN_MEAN    ((void*)2)  /* >= mean   is set to 1 */

/*----------------------------------------------------------------------
  Functions
----------------------------------------------------------------------*/
#ifdef REAL
extern void* binarize    (REAL   *data, int N, int T,
                          REAL   *thhs, int bpi);
#else
extern void* binarizeFlt (float  *data, int N, int T,
                          float  *thhs, int bpi);
extern void* binarizeDbl (double *data, int N, int T,
                          double *thhs, int bpi);
#endif

#endif  /* #ifndef __BINARIZE__ */

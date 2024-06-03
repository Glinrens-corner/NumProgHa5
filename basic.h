
#ifndef BASIC_H
#define BASIC_H

#include <math.h>
#include <stdbool.h>
#include <assert.h>

#define BLACK "\033[0;30m"
#define BBLACK "\033[1;30m"
#define RED "\033[0;31m"
#define BRED "\033[1;31m"
#define GREEN "\033[0;32m"
#define BGREEN "\033[1;32m"
#define YELLOW "\033[0;33m"
#define BYELLOW "\033[1;33m"
#define BLUE "\033[0;34m"
#define BBLUE "\033[1;34m"
#define PURPLE "\033[0;35m"
#define BPURPLE "\033[1;35m"
#define CYAN "\033[0;36m"
#define BCYAN "\033[1;36m"
#define WHITE "\033[0;37m"
#define BWHITE "\033[1;37m"
#define NORMAL "\033[0;0m"

/* ------------------------------------------------------------
 * Float type
 * ------------------------------------------------------------ */

#ifdef USE_FLOAT
typedef float real;
#define ABS(x) fabsf(x)
#define SQRT(x) sqrtf(x)
#else
typedef double real;
#define ABS(x) fabs(x)
#define SQRT(x) sqrt(x)
#endif

/* ------------------------------------------------------------
 * Prefix for C function declarations
 * ------------------------------------------------------------ */

#ifdef __cplusplus
#define HEADER_PREFIX extern "C"
#else
#define HEADER_PREFIX
#endif

/* ------------------------------------------------------------
 * Stopwatch
 * ------------------------------------------------------------ */

typedef struct _stopwatch stopwatch;
typedef stopwatch *pstopwatch;

HEADER_PREFIX pstopwatch
new_stopwatch();

HEADER_PREFIX void
del_stopwatch(pstopwatch sw);

HEADER_PREFIX void
start_stopwatch(pstopwatch sw);

HEADER_PREFIX real
stop_stopwatch(pstopwatch sw);

#endif

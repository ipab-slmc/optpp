
#ifndef PDS_INC
#define PDS_INC

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#if defined __cplusplus && defined HAVE_STD
#include <cstdio>
#else
#include <stdio.h>
#endif

#if defined(__cplusplus)
#include<algorithm>
#else
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

#include "proto.h"

#define INT_SIZE        sizeof (int)
#define DBL_SIZE        sizeof (double)

#if (defined(__SUNPRO_C) || defined(__SUNPRO_CC))
#define READ_TYPE (char *)
#endif
#ifndef READ_TYPE
#define READ_TYPE (void *)
#endif

#if (defined(__SUNPRO_C) || defined(__SUNPRO_CC))
#define WRITE_TYPE (char *)
#endif
#ifndef WRITE_TYPE
#define WRITE_TYPE (void *)
#endif

#if (defined(__SUNPRO_C) || defined(__SUNPRO_CC))
#define FREE_TYPE (char *)
#endif
#ifndef FREE_TYPE
#define FREE_TYPE (void *)
#endif

#endif


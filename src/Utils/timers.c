/*------------------------------------------------------------------------
// Copyright (C) 1993,1994: 
// J.C. Meza
// Sandia National Laboratories
// meza@california.sandia.gov
//----------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include "OPT++_config.h"
#endif

#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#if defined(HAVE_SYS_TIMES_H) && !(defined(_MSC_VER) || defined(__MINGW32__))
#include <sys/times.h>
#endif
#if defined(HAVE_SYS_RESOURCE_H) && !(defined(_MSC_VER) || defined(__MINGW32__))
#include <sys/resource.h>
#endif
#include <stddef.h>

#if !defined(HAVE_TIMES) && (defined(_MSC_VER) || defined(__MINGW32__))
#include <windows.h>
#endif

#ifndef HZ
#define HZ 100
#endif

double get_cpu_time()
{
/* ********************************************************************
**
**  Name: get_cpu_time
**
**  Purpose: general purpose CPU timing routine.
**
**  Arguments: none.
**
**  Return Value: user CPU time in (double) seconds.
**
**  Revision History:
**
**  10-May-94 -- initial development of get_cpu_time ().
**
** *******************************************************************/
#if defined(HAVE_SYS_TIMES_H) && !(defined(_MSC_VER) || defined(__MINGW32__))
    struct tms tms;
#endif
    double time;
/*
**  Begin get_cpu_time.
*/
#if defined(HAVE_TIMES) && !(defined(_MSC_VER) || defined(__MINGW32__))
    times (&tms);
    time = (double) tms.tms_utime / (double) HZ;
#else
    FILETIME creationTime,kernTime,userTime,exitTime;
    HANDLE process = GetCurrentProcess();
    GetProcessTimes(process,&creationTime,&exitTime,&kernTime,&userTime);
    SYSTEMTIME sysTime;
    FileTimeToSystemTime(&userTime,&sysTime);
    time = (double)sysTime.wSecond + (double)sysTime.wMilliseconds * 0.001;
#endif
    return time;
/*
**  End get_cpu_time.
*/
}

double get_wall_clock_time()
{
/* ********************************************************************
**
**  Name: get_wall_clock_time
**
**  Purpose: general purpose wall-clock timing routine.
**
**  Arguments: none.
**
**  Return Value: time in (double) seconds since the Epoch.
**
**  Notes: The Paragon specific dclock() routine is used to avoid
**         unnecessary references from each node back to the boot
**         node as is required for system calls like gettimeofday(),
**         getrusage() or times().
**
**  Revision History:
**
**  10-May-94 -- TXF; initial development of get_wall_clock_time ().
**
** *******************************************************************/

    double time;
#ifdef HAVE_SYS_TIME_H
    struct timeval tp;
#endif

/*
**  Begin get_wall_clock_time.
*/

    void* tzp = 0;
    gettimeofday (&tp, tzp);
    time = (double) tp.tv_sec + ((double) tp.tv_usec / (double) 1.0e06);

    return(time);
/*
**  End get_wall_clock_time.
*/
}

/* Modified from http://mywebpage.netscape.com/yongweiwu/timeval.h.txt */
#if !defined(HAVE_GETTIMEOFDAY) && (defined(_MSC_VER) || defined(__MINGW32__))
int gettimeofday (struct timeval *tv, void* tz)
{
  union {
    __int64 ns100; /*time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } now;

  GetSystemTimeAsFileTime (&now.ft);
  tv->tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
  tv->tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL);
  return (0);
}
#endif

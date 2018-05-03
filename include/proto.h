
#ifndef PROTO_INC
#define PROTO_INC

/* ********************************************************************
**
**  Name: proto.h.
**
**  Purpose: prototype definitions for PDS.
**
**  Revision History:
**
**  24-Mar-94 -- initial development of proto.h ().
**  12-May-94 -- changed name of search module to make_search.
**  02-Jun-94 -- alphabetized the list and removed out of date entries.
**
** *******************************************************************/
/*
**  Begin proto.h.
*/

#ifdef __cplusplus
extern "C" {
#endif

int pdscom(char *);

int pdsdone(int, int, int, double, double *, double *, double *,
	    double, double, double, double, int, int, char *, int);

int depth(int, int, int);

int pdseql(int, double, double *);
     
int pdsget(int, FILE *, int *, double *, int *, char *);

int pdsglb(int, double *, double *, char *);

int pdsgop(double *, int, double *, char *);

int pdshrk(int, int, int *, int *);

int make_search(int, int, int *, int *, int *, int *, int *, int *,
		int *);

double pdslen(int, int, double *, double, double *);

int order(int, int *, int *);

int quick(int, int *, int *);

int pdsrgt(int, double, double *);

int pdscld(int, double, double *);

int sort(int, int *, int *, int *, int *);

int writes(int, int, int, int, int, int *, int *);

int pdsdgn(int, double *, double *, double *, double *, int *, double *);

int pdsupd(long int, int, int, int *, double *, double *, double);

#ifdef __cplusplus
}
#endif
/*
**  End proto.h.
*/
#endif

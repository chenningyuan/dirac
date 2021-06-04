/*
 *
 * FILE    : dirac.h
 * DESCR.  : Header file for Dirac subroutine
 * AUTHOR  : J. Thyssen - 980709
 *           portability to AIX and HPUX simplified by P. Salek, 020608
 *
 */

#define MACRO_DIRAC dirac_
#define MACRO_DIRNOD dirnod_

#if defined (INT_STAR8) || defined (IA64)
void MACRO_DIRAC(long *iparcal, long *imytid,   \
                 long *imparid , long *inumnod, \
                 long *mempointer, long *lwork);

void MACRO_DIRNOD(long *iparcal, long *imytid,   \
                  long *imparid , long *inumnod, \
                  long *mempointer, long *lwork);
#else
void MACRO_DIRAC(int *iparcal, int *imytid,   \
                 int *imparid , int *inumnod, \
                 int *mempointer, int *lwork);

void MACRO_DIRNOD(int *iparcal, int *imytid,   \
                  int *imparid , int *inumnod, \
                  int *mempointer, int *lwork);
#endif

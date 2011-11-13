/*
 *
 * VERSION : $Revision: 5419 $
 * DATE    : $Date: 2007-01-19 13:41:33 +0100 (Fri, 19 Jan 2007) $
 * FILE    : dirac.h
 * DESCR.  : Header file for Dirac subroutine
 * AUTHOR  : J. Thyssen - 980709
 *           portability to AIX and HPUX simplified by P. Salek, 020608
 *
 */

#if defined (SYS_CRAY) || defined (SYS_T3D)
#define MACRO_DIRAC DIRAC
#define MACRO_DIRNOD DIRNOD
#else
#define MACRO_DIRAC dirac_
#define MACRO_DIRNOD dirnod_
#endif

#if defined (INT_STAR8) || defined (IA64)
void MACRO_DIRAC(long *iparcal, long *imytid,   \
                 long *imparid , long *inumnod, \
                 long *mempointer, long *lwork);

void MACRO_DIRNOD(long *iparcal, long *imytid,  \
                 long *imparid , long *inumnod, \
                 long *mempointer, long *lwork);
#else
void MACRO_DIRAC(int *iparcal, int *imytid,   \
                 int *imparid , int *inumnod, \
                 int *mempointer, int *lwork);

void MACRO_DIRNOD(int *iparcal, int *imytid,  \
                 int *imparid , int *inumnod, \
                 int *mempointer, int *lwork);
#endif

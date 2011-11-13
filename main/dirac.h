/*
 *
 * VERSION : $Revision: 1.4 $
 * DATE    : $Date: 2002/11/20 15:55:07 $
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
void MACRO_DIRAC(int *iparcal, int *imytid,   \
                 int *imparid , int *inumnod, \
                 int *mempointer, int *lwork);

void MACRO_DIRNOD(int *iparcal, int *imytid,  \
                 int *imparid , int *inumnod, \
                 int *mempointer, int *lwork);

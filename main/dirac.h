/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Header file for Dirac subroutine
 * AUTHOR  : J. Thyssen - 980709
 *
 */

#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_HPUX)
#define MACRO_DIRAC DIRAC
#define MACRO_DIRNOD DIRNOD
#elif defined (SYS_AIX)
#define MACRO_DIRAC dirac
#define MACRO_DIRNOD dirnod
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

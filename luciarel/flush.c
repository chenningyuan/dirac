/* RCS *******************************************************************/
/* $Id$ */
/* $Log$
 * Revision 1.2  2002/07/09 18:22:50  tow
 * hjaaj:
 * 1) merge of krmc-lucia-branch into main trunk as new krmc-esr-branch
 *    ALL FURTHER MCSCF development should be on this branch!
 * 1a) conflict between files in lucita/ and luciarel/ with same name
 *    needs to be resolved !!!
 * 2) some mcscf modifications as preparation for new multistate ESR module
 *    - compiles and passes ftnchek, but has not been tested yet with
 *      mcscf test runs
 *
/* Revision 1.1.2.1  2001/01/10 17:36:44  fleig
/* Merged LUCIAREL into DIRAC.
/*
 * Revision 1.1  1992/09/17 10:56:23  marco
 * Initial revision
 * */
/* RCS *******************************************************************/

#include <stdio.h>

#if !defined(IRIS) && !defined(CONVEX)
int flush(iunit)
   long int *iunit;
{
   if (*long iunit == 6) fflush(stdout);
}
#endif

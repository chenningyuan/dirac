/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Execute Dirac or Dirnod
 * AUTHOR  : J. Thyssen - 980713 
 *
 */

#include <stdio.h>

#include "exec_dirac.head"
 ;
#include "dirac.h"


#include "exec_dirac.head"

{
   if ( (*mytid) > 0 )
   {
#if defined (VAR_MPI)
     /*
      * Execute slave
      */

      MACRO_DIRNOD(parcal,mytid,mparid,numnod,mempointer,size);
#else
      printf("Trying to execute DIRNOD in a non-MPI version\n");
      exit(1);
#endif
   }
   else
   {

     /*
      * Execute master
      */

      MACRO_DIRAC(parcal,mytid,mparid,numnod,mempointer,size);
   }
}

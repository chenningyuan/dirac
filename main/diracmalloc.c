/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Allocate memory and return pointer
 *           This subroutine contains all the system-dependent stuff
 * AUTHOR  : J. Thyssen - 980709
 *
 */

#include <stdlib.h>


#include "diracmalloc.head"
 ;

#if defined (SYS_CRAY) || defined (SYS_T3D)
#define CALL_QUIT QUIT
#elif defined (SYS_AIX || defined (SYS_HPUX)
#define CALL_QUIT quit
#else 
#define CALL_QUIT quit_
#endif /* SYS_CRAY */

#include "diracmalloc.head"

{

void *pointer;

#if defined (SYS_AIX)  || defined (SYS_IRIX)  || defined (SYS_SUN) \
 || defined (SYS_HPUX) || defined (SYS_LINUX) || defined (SYS_CRAY) \
 || defined (SYS_T3D)  
   pointer = malloc(size);
   if ( pointer == NULL )
   {
      printf("DIRAC could not allocate memory (%i bytes).\n",size);
      CALL_QUIT("Dirac could not allocate memory",31); 
   }
#elif defined (SYS_CONVEX)
   pointer = nalloc(size);
   if ( pointer == NULL )
   {
      printf("DIRAC could not allocate memory (%i bytes).\n",size);
      CALL_QUIT("Dirac could not allocate memory",31); 
   }
#else
   Malloc????
#endif
   return pointer;
}





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


#include "diracmalloc.head"

{

void *pointer;

#if defined (SYS_AIX)  || defined (SYS_IRIX)  || defined (SYS_SUN) \
 || defined (SYS_PHUX) || defined (SYS_LINUX) 
   pointer = malloc(size);
   if ( pointer == NULL )
   {
      printf("DIRAC could not allocate memory.\n");
      exit(1);
   }
#elif defined (SYS_CRAY) || defined (SYS_T3D)  
   pointer = hpalloc(size);
   if ( pointer == NULL )
   {
      printf("DIRAC could not allocate memory (%i bytes).\n",size);
      exit(1);
   }
#elif defined (SYS_CONVEX)
   pointer = nalloc(size);
   if ( pointer == NULL )
   {
      printf("DIRAC could not allocate memory (%i bytes).\n",size);
      exit(1);
   }
#else
   Malloc????
#endif
   return pointer;
}





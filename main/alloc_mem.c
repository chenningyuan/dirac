/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Allocate memory for dirac
 *           Get memory size either from enviroment variables or
 *           default values.
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdlib.h>
#include <stdio.h>


#include "alloc_mem.head"
 ;
#include "diracmalloc.head"
 ;
#include "get_memsize.head"
 ;

#define SIZE_REAL (sizeof(double))

#include "alloc_mem.head"

{
   if (! get_memsize(mytid,size))
   {
      *size = INSTALL_WRKMEM;
   }
   if ( ! mytid )
      printf("DIRAC allocating %i words of memory\n",*size);
   return diracmalloc(SIZE_REAL*(*size));
}





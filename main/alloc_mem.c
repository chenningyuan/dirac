/*
       Copyright (c) 2004 by the authors of Dirac (see below).
       All Rights Reserved.
    
       The source code in this file is part of
       Dirac, a relativistic ab initio electronic structure program,
       release DIRAC04.0 (2004),
       written by H. J. Aa. Jensen, T. Saue, and L. Visscher
       with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
       O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
       P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
       P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
       (http://dirac.chem.sdu.dk).
    
       For a suitable BibTEX entry, see:
          http://dirac.chem.sdu.dk/doc/reference.shtml
    
       This source code is provided under a written licence and may be
       used, copied, transmitted, or stored only in accord with that
       written licence.
    
       In particular, no part of the source code or compiled modules may
       be distributed outside the research group of the licence holder.
       This means also that persons (e.g. post-docs) leaving the research
       group of the licence holder may not take any part of Dirac,
       including modified files, with him/her, unless that person has
       obtained his/her own licence.
    
       For questions concerning this copyright write to:
          dirac-admin@dirac.chem.sdu.dk
    
       For information on how to get a licence see:
          http://dirac.chem.sdu.dk
*/
/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * FILE    : alloc_mem.c
 * DESCR.  : Allocate memory for dirac
 *           Get memory size either from enviroment variables or
 *           default values.
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <sys/utsname.h>



#include "alloc_mem.head"
 ;
#include "diracmalloc.head"
 ;
#include "get_memsize.head"
 ;

#define SIZE_REAL (sizeof(double))

#include "alloc_mem.head"

{
   struct utsname *un;

   un = (struct utsname *) malloc(sizeof(struct utsname));

   if ( un != NULL )
   {
      uname(un);
   }

   if (! get_memsize(mytid,size))
   {
      *size = INSTALL_WRKMEM;
   }

   if ( un != NULL )
   {
     if ( ! mytid )
     {
#if defined (INT_STAR8)
      printf("Master   ( %s ): DIRAC allocating %li words of memory\n",un->nodename,*size);
#else
      printf("Master   ( %s ): DIRAC allocating %i words of memory\n",un->nodename,*size);
#endif
      fflush( stdout );
     }
     else
     {
#if defined (INT_STAR8)
      printf("Node %3i ( %s ): DIRAC allocating %li words of memory\n",mytid,un->nodename,*size);
#else
      printf("Node %3i ( %s ): DIRAC allocating %i words of memory\n",mytid,un->nodename,*size);
#endif
      fflush( stdout );
     }
     free(un);
   }
   else
   {
     if ( ! mytid )
     {
#if defined (INT_STAR8)
      printf("Master : DIRAC allocating %li words of memory\n",*size);
#else
      printf("Master : DIRAC allocating %i words of memory\n",*size);
#endif
      fflush( stdout );
     }
     else
     {
#if defined (INT_STAR8)
      printf("Node %3i : DIRAC allocating %li words of memory\n",mytid,*size);
#else
      printf("Node %3i : DIRAC allocating %i words of memory\n",mytid,*size);
#endif
      fflush( stdout );
     }
   }
     
   return diracmalloc((size_t) SIZE_REAL*(*size));
}

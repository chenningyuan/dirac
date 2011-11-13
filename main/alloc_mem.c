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
 * VERSION : $Revision: 1.11.12.1 $
 * DATE    : $Date: 2004/09/01 08:13:49 $
 * FILE    : alloc_mem.c
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
      printf("Master : DIRAC allocating %i words of memory\n",*size);
   else
   {
      printf("Node %2i: DIRAC allocating %i words of memory\n",mytid,*size);
      fflush( stdout );
   }
   return diracmalloc((size_t) SIZE_REAL*(*size));
}





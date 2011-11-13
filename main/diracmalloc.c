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
 * VERSION : $Revision: 1.13.12.1 $
 * DATE    : $Date: 2004/09/01 08:13:51 $
 * FILE    : diracmalloc.c
 * DESCR.  : Allocate memory and return pointer
 *           This subroutine contains all the system-dependent stuff
 * AUTHOR  : J. Thyssen - 980709
 *           portability to AIX and HPUX simplified by P. Salek, 020608
 *
 */

#include <stdlib.h>


#include "diracmalloc.head"
 ;

#if defined (SYS_CRAY) || defined (SYS_T3D)
#define CALL_QUIT QUIT
#else 
#define CALL_QUIT quit_
#endif /* SYS_CRAY */

#include "diracmalloc.head"

{

void *pointer;

#if defined (SYS_AIX)  || defined (SYS_IRIX)  || defined (SYS_SUN) \
 || defined (SYS_HPUX) || defined (SYS_LINUX) || defined (SYS_CRAY) \
 || defined (SYS_T3D)  || defined (SYS_DEC)   || defined (SYS_SX) \
 || defined (SYS_DARWIN)
   pointer = malloc(size);
   if ( pointer == NULL )
   {
      printf("DIRAC could not allocate memory (%li bytes).\n",size);
      CALL_QUIT("Dirac could not allocate memory",31); 
   }
#elif defined (SYS_CONVEX)
   pointer = nalloc(size);
   if ( pointer == NULL )
   {
      printf("DIRAC could not allocate memory (%li bytes).\n",size);
      CALL_QUIT("Dirac could not allocate memory",31); 
   }
#else
   printf("\n MALLOC not defined for this architecture !\n");
   exit(1);
#endif
   return pointer;
}

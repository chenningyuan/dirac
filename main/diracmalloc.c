 /*
  *
  *...   Copyright (c) 1998 by the authors of Dirac (see below).
  *...   All Rights Reserved.
  *...
  *...   The source code in this file is part of 
  *...   "Dirac, a relativistic ab initio electronic structure program, 
  *...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
  *...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
  *...
  *...   This source code is provided under a written licence and may be 
  *...   used, copied, transmitted, or stored only in accord with that
  *...   written licence.
  *...
  *...   In particular, no part of the source code or compiled modules may 
  *...   be distributed outside the research group of the licence holder. 
  *...   This means also that persons (e.g. post-docs) leaving the research
  *...   group of the licence holder may not take any part of Dirac,
  *...   including modified files, with him/her, unless that person has
  *...   obtained his/her own licence. 
  *...
  *...   For questions concerning this copyright write to:
  *...      dirac-admin@dirac.chem.ou.dk
  *...
  *...   For information on how to get a licence see:
  *...      http://dirac.chem.ou.dk/Dirac
  *
  *
  */
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
#elif defined (SYS_AIX)|| defined (SYS_HPUX)
#define CALL_QUIT quit
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
   Malloc????
#endif
   return pointer;
}





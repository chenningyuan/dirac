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
  *...   group of the licence holder may not take any part of Dalton,
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
 * DESCR.  : Print logo if master node
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>

#include "print_logo.head"
 ;


#include "print_logo.head"

{
#include "title.h"
#include "logo.h"
#include "subtitle.h"

   if (mytid != 0) return;

   printf("%s",strtitle);
   printf("%s",strlogo);
   printf("%s",strsubtitle);

   sysinf(parcal);
   
}

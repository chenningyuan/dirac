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
 * DESCR.  : Print timing statistic to stdout
 * AUTHOR  : J. Thyssen - 980713 
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "pr_timing.head"
 ;
#include "csectid.head"
 ;
#include "diracmalloc.head"
 ;


#include "pr_timing.head"

{
   char *sectid;

   if ( ! mytid )
      printf("*****************************************************\n");
   sectid = (char *) malloc(30);
   csectid(gt_last->cpu - gt_first->cpu,sectid);
   if ( numnod > 0 )
      printf(">>>> Total CPU  time used in DIRAC (Node %i): %s\n", \
         mytid,sectid);
   else
      printf(">>>> Total CPU  time used in DIRAC: %s\n",sectid);

   csectid(gt_last->wall - gt_first->wall,sectid);
   if ( numnod > 0 )
      printf(">>>> Total WALL time used in DIRAC (Node %i): %s\n", \
         mytid,sectid);
   else
      printf(">>>> Total WALL time used in DIRAC: %s\n",sectid);
   free(sectid);

   if ( ! mytid )
      printf("********** E N D   of   D I R A C  output  **********\n");

}

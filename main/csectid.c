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
 * DESCR.  : Return nicely formatted character string
 *           C-version of gp/gpsaue.F:SECTID
 * AUTHOR  : J. Thyssen - 980713 (T. Saue)
 *
 */

#include <stdio.h>
#include <string.h>

#include "csectid.head"
 ;


#include "csectid.head"

{
   int itid, iday, ihour, imin;

   itid = (int) tid;

   if (itid >= 86400)
   {
      iday  = itid / 86400;
      itid  = itid % 86400;
      ihour = itid / 3600;
      itid  = itid % 3600;
      imin  = itid / 60;
      if ( ! sprintf(sectid,"%3id%2ih%2imin",iday,ihour,imin) )
         strcpy(sectid,"undef");
   }
   else if ( itid >= 3600 )
   {
      ihour = itid / 3600;
      itid  = itid % 3600;
      imin  = itid / 60;
      itid  = itid % 60;
      if ( ! sprintf(sectid,"%2ih%2imin%2is",ihour,imin,itid) )
         strcpy(sectid,"undef");
   }
   else if ( itid >= 60 )
   {
      imin  = itid / 60;
      tid   = tid - imin*60;
      if ( ! sprintf(sectid,"%2imin%6.3fs",imin,tid) )
         strcpy(sectid,"undef");
   }
   else
   {
      if ( ! sprintf(sectid,"%11.8fs",tid) )
         strcpy(sectid,"undef");
   }
}

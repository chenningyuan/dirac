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
 * VERSION : $Revision: 1.4.12.1 $
 * DATE    : $Date: 2004/09/01 08:13:50 $
 * FILE    : csectid.c
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

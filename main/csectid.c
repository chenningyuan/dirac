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

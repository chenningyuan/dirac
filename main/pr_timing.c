/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Print timing statistic to stdout
 * AUTHOR  : J. Thyssen - 980713 
 *
 */

#include "pr_timing.head"
 ;
#include "csectid.head"
 ;
#include "diracmalloc.head"
 ;


#include "pr_timing.head"

{
   char *sectid;

   sectid = (char *) malloc(30);
   csectid(cpu_used,sectid);
   if ( numnod > 0 )
      printf(">>>> Total CPU  time used in DIRAC (Node %i): %s\n", \
         mytid,sectid);
   else
      printf(">>>> Total CPU  time used in DIRAC: %s\n",sectid);

   csectid(cpu_used,sectid);
   if ( numnod > 0 )
      printf(">>>> Total WALL time used in DIRAC (Node %i): %s\n", \
         mytid,sectid);
   else
      printf(">>>> Total WALL time used in DIRAC: %s\n",sectid);
   free(sectid);

}

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

#include "pr_timing.head"
 ;
#include "csectid.head"
 ;
#include "diracmalloc.head"
 ;


#include "pr_timing.head"

{
   char *sectid;

#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_HPUX)
    FLSHFO(stdout);
#elif defined (SYS_AIX)
   flshfo(stdout);
#else
   flshfo_(stdout);
#endif

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

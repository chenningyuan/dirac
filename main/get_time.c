/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Get timing
 * AUTHOR  : J. Thyssen - 980716
 *
 */

#include "get_time.head"
 ;

#include "get_time.struct"


#include "get_time.head"

{
#if defined (SYS_CRAY) || defined (SYS_T3D)
   GETTIM(&(gt->cpu),&(gt->wall));
#elif defined (SYS_AIX) || defined (SYS_HPUX)
   gettim(&(gt->cpu),&(gt->wall));
#else
   gettim_(&(gt->cpu),&(gt->wall));
#endif
}

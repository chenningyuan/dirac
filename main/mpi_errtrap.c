/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Trap MPI errors
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>
#include <string.h>
#if defined (VAR_MPI)
#include <mpi.h>
#endif

#include "mpi_errtrap.head"
 ;


#include "mpi_errtrap.head"

{
   int len;

#if defined (VAR_MPI)

#if defined (SYS_AIX)

   flshfo(stdout);

   fprintf(stdout,"******************************************\n");
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"Comm                 %d\n",*comm);
   fprintf(stdout,"Mesg                 %s\n",errmesg);
   fprintf(stdout,"????                 %d\n",*oeh1);
   fprintf(stdout,"????                 %d\n",*oeh2);
   fprintf(stdout,"******************************************\n");

   len = strlen(errmesg);
   quit(errmesg,len);
   
#else /* SYS_AIX */

#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_HPUX)
    FLSHFO(stdout);
#else
   flshfo_(stdout);
#endif

   fprintf(stdout,"******************************************\n");
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"Comm                 %d\n",*comm);
   fprintf(stdout,"The PID is %d\n",getpid());
   fprintf(stdout,"******************************************\n");

   len = 0;
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_HPUX)
   QUIT(NULL,len);
#else /* SYS_CRAY etc. */
   quit_(NULL,len);
#endif /* SYS_CRAY etc. */


#endif /* SYS_AIX */

#endif /* VAR_MPI */
}

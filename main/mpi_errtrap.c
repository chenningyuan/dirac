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
#if defined (VAR_MPI)
#include <mpi.h>
#endif

#include "mpi_errtrap.head"
 ;


#include "mpi_errtrap.head"

{
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
   while(1);
   
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
   while(1);

#endif /* SYS_AIX */

#endif /* VAR_MPI */
}

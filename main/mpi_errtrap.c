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
   char *errorstring;

#if defined (VAR_MPI)

#if defined (SYS_AIX)

/*   flshfo(stdout); */

   fprintf(stdout,"******************************************\n");
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"Comm                 %d\n",*comm);
   fprintf(stdout,"Function             %s\n",errmesg);
   fprintf(stdout,"PID                  %d\n",getpid());
   len = 100;
   errorstring = (char *) diracmalloc(len);
   MPI_Error_string(*err,errorstring,&len);
   fprintf(stdout,"Error string         %s\n",errorstring);
   fprintf(stdout,"******************************************\n");
   quit(errorstring,len);
   free(errorstring);

#else /* SYS_AIX */

#if defined (SYS_CRAY) || defined (SYS_T3D)
    FLSHFO(stdout);
#elif defined (SYS_AIX)  || defined (SYS_HPUX)
   flshfo(stdout);
#else
   flshfo_(stdout);
#endif

   fprintf(stdout,"******************************************\n");
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"Comm                 %d\n",*comm);
   fprintf(stdout,"The PID              %d\n",getpid());
   len = 100;
   errorstring = (char *) diracmalloc(len);
   MPI_Error_string(*err,errorstring,&len);
   fprintf(stdout,"Error string         %s\n",errorstring);
   fprintf(stdout,"******************************************\n");

#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_HPUX)
   QUIT(errorstring,len);
#else /* SYS_CRAY etc */
   quit_(errorstring,len);
#endif /* SYS_CRAY etc */
   free(errorstring);

#endif /* SYS_AIX */

#endif /* VAR_MPI */
}

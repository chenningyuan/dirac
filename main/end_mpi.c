/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Kill MPI-slaves if parallel run
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>

#if defined (VAR_MPI)
#include <infmpi.h>
#include <mpi.h>
#endif /* VAR_MPI */

#include "end_mpi.head"
 ;


#include "end_mpi.head"

{

#if defined (VAR_MPI)
   int ierr;

   if ( mytid == 0 )
   {

      /*
       *  Send send signal to slaves
       */

      int nfinis = -1;
      if ( ! ierr = MPI_Bcast(&nfinis,1,mpi_integer,mparid,mpi_comm_world) )
         printf("MPI_Bcast error %i in end_mpi()\n",ierr);
   }
   else
   {

      /*
       *  Finalize MPI
       */
      
#if defined (SYS_IRIX)
      sleep(1);
#endif

      if ( ! ierr = MPI_Finalize() )
         printf("MPI_Finalize error %i in end_mpi()\n",ierr);
   }

   /*
    * Free MPI error handler
    */

   if ( ! ierr = MPI_Errhandler_Free(&errhandler) )
      printf("MPI_Errhandler_Free error %i in end_mpi()\n",ierr);

#endif  /* VAR_MPI */
}

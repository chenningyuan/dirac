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
#include <infmpi.h>
#include <mpi.h>
#endif

#include "mpi_errtrap.head"
 ;


#include "mpi_errtrap.head"

{
#if defined (VAR_MPI)
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"The PID is %d\n",getpid());
   while(1);
#endif
}

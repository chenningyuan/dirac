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
#include <infmpi.h>
#include <mpi.h>

#include "mpi_errtrap.head"
 ;


#include "mpi_errtrap.head"

{
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"The PID is %d\n",getpid());
   while(1);
}

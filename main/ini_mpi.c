/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Initialize MPI if parallel run
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>

#if defined (VAR_MPI)
#define schedmode 13
#define sgs_gang   4
#include <infmpi.h>
#include <mpi.h>
#endif /* VAR_MPI */

#include "ini_mpi.head"
 ;
#include "mpi_errtrap.head"
 ;


#include "ini_mpi.head"

{

#if defined (SYS_IBM)
 /*
  * IBM VM/CMS: 212 : disable end of record error for formatted I/O'
  * (such that variable length input records may be read)
  *
  */
  errset(212,256,-1,-1);
#endif /* SYS_IBM */

#if defined (VAR_MPI)

   int par,lwrkmm,ierr;

   /*
    * Check if this is a parallel calculation or not;
    * check for DIRPAR in environment
    */

   par = 0
   dirpar = diracmalloc(2);
   dirpar = getenv("DIRPAR");

   if (dirpar != NULL)
   {
      par = atoi(dirpar);
   }
   else
   {
      printf("DIRPAR error; DIRPAR = %s\n",dirpar);
      par = 0;
   }

   free(dirpar);

   /*
    * MPI start
    */

   if ( ! ierr = MPI_Init(NULL,NULL) )
      printf(stderr,"MPI_Init error %i in ini_mpi()\n",ierr);

   /*
    * Set up MPI error handler
    */

   if ( ! ierr = MPI_Errhandler_create(&mpi_errtrap,errhandler) )
      printf(stderr,"MPI_Errhandler_create error %i in ini_mpi()\n",ierr);

   if ( ! ierr = MPI_Errhandler_set(MPI_Comm_world,errhandler) )
      printf(stderr,"MPI_Errhandler_create error %i in ini_mpi()\n",ierr);
   
   /*
    * Get number of nodes and the number of this node
    */

   if ( ! ierr = MPI_Comm_rank(mpi_comm_world,mytid) )
      printf(stderr,"MPI_Comm_rank error %i in ini_mpi()\n",ierr);

   if ( ! ierr = MPI_Comm_size(mpi_comm_world,numnod) )
      printf(stderr,"MPI_Comm_size error %i in ini_mpi()\n",ierr);

#if defined (SYS_IRIX)
   /*
    * Impose gang scheduling to avoid problems with idle nodes
    * when high system load
    *
    */
   if ( (mytid == 0) || (schedctl(schedmode,sgs_gang,0) < 0) )
   {
      printf("*** ERROR gang scheduling\n");
      quit("Error in gang scheduling!");
   }
#endif /* SYS_IRIX */
   (*numnod)--;
   *mparid = 0;
    
#else  /* VAR_MPI */
   *mytid =  0;
   *parcal = 0;
   *mparid = 0;
   *numnod = 0;
#endif  /* VAR_MPI */
}

 /*
  *
  *...   Copyright (c) 1998 by the authors of Dirac (see below).
  *...   All Rights Reserved.
  *...
  *...   The source code in this file is part of 
  *...   "Dirac, a relativistic ab initio electronic structure program, 
  *...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
  *...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
  *...
  *...   This source code is provided under a written licence and may be 
  *...   used, copied, transmitted, or stored only in accord with that
  *...   written licence.
  *...
  *...   In particular, no part of the source code or compiled modules may 
  *...   be distributed outside the research group of the licence holder. 
  *...   This means also that persons (e.g. post-docs) leaving the research
  *...   group of the licence holder may not take any part of Dirac,
  *...   including modified files, with him/her, unless that person has
  *...   obtained his/her own licence. 
  *...
  *...   For questions concerning this copyright write to:
  *...      dirac-admin@dirac.chem.ou.dk
  *...
  *...   For information on how to get a licence see:
  *...      http://dirac.chem.ou.dk/Dirac
  *
  *
  */
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
#include <stdlib.h>

#if defined (VAR_MPI)
#define schedmode 13
#define sgs_gang   4
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

   int par,ierr;
   char *dirpar;

   /*
    * Check if this is a parallel calculation or not;
    * check for DIRPAR in environment
    */

   par = 0;
   dirpar = getenv("DIRPAR");

   if (dirpar != NULL)
   {
      par = atoi(dirpar);
   }
   else
   {
      par = 0;
   }

   if ( ! par )
   {
      *mytid =  0;
      *parcal = 0;
      *mparid = 0;
      *numnod = 0;
   }
   else
   {

      /*
       * MPI start
       */

      if ( (ierr = MPI_Init(NULL,NULL)) )
         printf("MPI_Init error %i in ini_mpi()\n",ierr);

      /*
       * Set up MPI error handler
       */

      if ( (ierr = MPI_Errhandler_create(&mpi_errtrap,errhandler)) )
         printf("MPI_Errhandler_create error %i in ini_mpi()\n",ierr);

      if ( (ierr = MPI_Errhandler_set(MPI_COMM_WORLD,*errhandler)) )
         printf("MPI_Errhandler_create error %i in ini_mpi()\n",ierr);

      /*
       * Get number of nodes and the number of this node
       */

      if ( (ierr = MPI_Comm_rank(MPI_COMM_WORLD,mytid)) )
         printf("MPI_Comm_rank error %i in ini_mpi()\n",ierr);

      if ( (ierr = MPI_Comm_size(MPI_COMM_WORLD,numnod)) )
         printf("MPI_Comm_size error %i in ini_mpi()\n",ierr);

#if defined (SYS_IRIX)
      /*
       * Impose gang scheduling to avoid problems with idle nodes
       * when high system load
       *
       */
      if ( (mytid == 0) && (schedctl(schedmode,sgs_gang,0) < 0) )
      {
         printf("*** ERROR gang scheduling\n");
         exit(1);
      }
#endif /* SYS_IRIX */
      (*numnod)--;
      *mparid = 0;
      *parcal = 1;
   }
    
#else  /* VAR_MPI */
   *mytid =  0;
   *parcal = 0;
   *mparid = 0;
   *numnod = 0;
#endif  /* VAR_MPI */
}

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
 * FILE    : end_mpi.c
 * DESCR.  : Kill MPI-slaves if parallel run
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>

#if defined (VAR_MPI)
#include <mpi.h>
#endif /* VAR_MPI */

#include "end_mpi.head"
 ;


#include "end_mpi.head"

{

#if defined (VAR_MPI)
   int ierr;


  /*
   *  Send finish signal to slaves from master
   */

/* hjaaj: signal to slaves is now sent in SUBROUTINE DIRAC /July 2000
   if ( mytid == 0 )
   {
      int nfinis = -1;
      if ( (ierr = MPI_Bcast(&nfinis,1,MPI_INT,mparid,MPI_COMM_WORLD)) )
         printf("MPI_Bcast error %i in end_mpi()\n",ierr);
   }
*/

   /*
    * Free MPI error handler
    */

   if (  ( ierr = MPI_Errhandler_free(&errhandler)) )
      printf("MPI_Errhandler_Free error %i in end_mpi()\n",ierr);

   /*
    *  Finalize MPI
    */
      
#if defined (SYS_IRIX)
      sleep(1);
#endif

#if defined (MPE)
   MPE_Finish_log("mpilog");
#endif

   if ( ( ierr = MPI_Finalize()) )
      printf("MPI_Finalize error %i in end_mpi()\n",ierr);


#endif  /* VAR_MPI */
}

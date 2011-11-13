/*
       Copyright (c) 2004 by the authors of Dirac (see below).
       All Rights Reserved.
    
       The source code in this file is part of
       Dirac, a relativistic ab initio electronic structure program,
       release DIRAC04.0 (2004),
       written by H. J. Aa. Jensen, T. Saue, and L. Visscher
       with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
       O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
       P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
       P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
       (http://dirac.chem.sdu.dk).
    
       For a suitable BibTEX entry, see:
          http://dirac.chem.sdu.dk/doc/reference.shtml
    
       This source code is provided under a written licence and may be
       used, copied, transmitted, or stored only in accord with that
       written licence.
    
       In particular, no part of the source code or compiled modules may
       be distributed outside the research group of the licence holder.
       This means also that persons (e.g. post-docs) leaving the research
       group of the licence holder may not take any part of Dirac,
       including modified files, with him/her, unless that person has
       obtained his/her own licence.
    
       For questions concerning this copyright write to:
          dirac-admin@dirac.chem.sdu.dk
    
       For information on how to get a licence see:
          http://dirac.chem.sdu.dk
*/
/*
 *
 * VERSION : $Revision: 1.18.4.1 $
 * DATE    : $Date: 2004/09/01 08:13:52 $
 * FILE    : ini_mpi.c
 * DESCR.  : Initialize MPI if parallel run
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>
#include <stdlib.h>

#if defined (VAR_MPI)
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
    * New default Nov 2002: parallel if DIRPAR not defined
    */

   dirpar = getenv("DIRPAR");

   if (dirpar != NULL)
      par = atoi(dirpar);
   else
      par = 1;

   if ( ! par )
   {
      *mytid =  0;
      *parcal = 0;
      *mparid = 0;
      *numnod = 1;
   }
   else
   {

      /*
       * MPI start
       */

      if ( (ierr = MPI_Init( &argc, &argv )) )
         printf("MPI_Init error %i in ini_mpi()\n",ierr);

      /*
       * Set up MPI error handler
       */

/* MI: use correct type casting for the error function
     'int MPI_Errhandler_create(MPI_Handler_function *function,MPI_Errhandler *errhandler)'   */
     /* if ( (ierr = MPI_Errhandler_create(&mpi_errtrap, errhandler)) )  */
      if ( (ierr = MPI_Errhandler_create((MPI_Handler_function *) &mpi_errtrap, errhandler)) )  
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

      /* printf(" mytid and numnod %i %i\n",*mytid,*numnod); */

      if ( *mytid == 0 ) sleep(1);
      /* if master then sleep 1 second
       * so the nodes can print to output before the master
       * (a synchronize would have been better) /hjaaj
       */

      *mparid = 0;
      *parcal = 1;

#if defined (MPE)
      /* initialize MPE if precompiler option MPE is set */
      MPE_Init_log();
      if ( *mytid == 0 ) {
	MPE_Describe_state(1, 2, "Init", "black");
	MPE_Describe_state(3, 4, "Wait", "red");
	MPE_Describe_state(5, 6, "Send/recv task", "yellow");
	MPE_Describe_state(7, 8, "Send/recv data", "green");
	MPE_Describe_state(9,10, "Calc", "blue");
      }
      MPE_Start_log();
#endif  /* MPE */

   }
    
#else  /* VAR_MPI */
   *mytid =  0;
   *parcal = 0;
   *mparid = 0;
   *numnod = 1;
#endif  /* VAR_MPI */
}

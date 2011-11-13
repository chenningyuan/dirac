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
 * VERSION : $Revision: 1.9.4.1 $
 * DATE    : $Date: 2004/09/01 08:13:51 $
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

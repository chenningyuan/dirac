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
 * VERSION : $Revision: 1.21.4.1 $
 * DATE    : $Date: 2004/09/01 08:13:52 $
 * FILE    : main.c
 * DESCR.  : Main program for dirac
 * AUTHOR  : J. Thyssen - 980709
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "alloc_mem.head"
 ;
#include "print_logo.head"
 ;
#include "ini_mpi.head"
 ;
#include "end_mpi.head"
 ;
#include "exec_dirac.head"
 ;
#include "pr_timing.head"
 ;
#include "get_time.struct"
#include "get_time.head"
 ;
#include "mpi_errtrap.h"

#include "diracmalloc.head"
 ;

#if defined (SYS_DEC)
extern int MAIN__( int argc, char** argv )
#elif defined (SYS_SX)
extern int c_main_( int argc, char** argv )
#elif defined (VAR_PGCC)
extern int MAIN_( int argc, char** argv )
#elif defined (LF95)
extern int MAIN__( int argc, char** argv )
#else
extern int main( int argc, char** argv )
#endif

{

   int mytid, size, parcal, mparid, numnod;
   errhandler_type errhandler;
   int *mempointer;
   get_time_struct *gt_first, *gt_last;

   /*
    * Start timing
    */


#if defined (SYS_SX)
   gt_first = (get_time_struct *) malloc(sizeof(get_time_struct));
   gt_last  = (get_time_struct *) malloc(sizeof(get_time_struct));
#else
   gt_first = (get_time_struct *) diracmalloc(sizeof(get_time_struct));
   gt_last  = (get_time_struct *) diracmalloc(sizeof(get_time_struct));
#endif

   get_time(gt_first);

   /*
    * Initialize MPI stuff
    */

   ini_mpi( &parcal, &mytid, &mparid, &numnod, &errhandler, 
            argc, argv );

   /*
    * Distribute environment variables
    */

   if ( parcal )
   {
      distribute_env( "BASDIR", mytid );
      distribute_env( "DIRWRK", mytid );
      distribute_env( "DIRNOD", mytid );
   }


   /*
    * Allocate memory
    */

   mempointer = alloc_mem(mytid,&size);

   /*
    * Print logo, sysinfo etc.
    */

   /* if master then sleep 1 second
    * so that the slaves can finish their output
    * before the master starts writing
    */
   if ( mytid == 0 ) sleep(1);

   if ( mytid == 0 ) print_logo(parcal);

   /*
    * Execute dirac
    */

   exec_dirac(&parcal,&mytid,&mparid,&numnod,mempointer,&size);

   /*
    * End timing
    */

   get_time(gt_last);

   /*
    * Write out timing
    */

   pr_timing(gt_first, gt_last, numnod, mytid);

   /*
    * End MPI
    */

   if ( parcal ) end_mpi(mparid,mytid,errhandler);

   /*
    * Free memory
    */

   free(mempointer);
   free(gt_first);
   free(gt_last);

   /*
    * End
    */

   exit(0);
}

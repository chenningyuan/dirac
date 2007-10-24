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
 * VERSION : $Revision$
 * DATE    : $Date$
 * FILE    : mpi_errtrap.c
 * DESCR.  : Trap MPI errors
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>
#include <string.h>

/* MI: this is already in "mpi_errtrap.head":
---------------------------------------------
#if defined (VAR_MPI)
#include <mpi.h>
#endif
*/

#include "mpi_errtrap.head"


{

#if defined (VAR_MPI)

/*MI: to satisfy pgcc  */
#include "diracmalloc.head"
 ; 

#if defined (INT_STAR8)
   long len;
   char *errorstring;
#else
   int len;
   char *errorstring;
#endif

#if defined (SYS_CRAY) || defined (SYS_T3D)
   FLSHFO(stdout);
#else
   flshfo_(stdout);
#endif

/* MI: "*err" declared in mpi_errtrap header  */   

   fprintf(stdout,"******************************************\n");
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"Comm                 %d\n",*comm);
   fprintf(stdout,"The PID              %d\n",getpid());
/* SK: KRMC MO-transformations produce file-close errors on IBM systems 
   in Juelich resulting in an endless MPI error loop --> this is a quick 
   and dirty fix until the real bug is found (introduced between 
   JULY 2007 and OCTOBER 2007.*/ 
#if defined (SYS_AIX)
   exit(0);
#endif
   len = 100;
   errorstring = (char *) diracmalloc(len);

/* MI: defined as
#include "mpi.h"
int MPI_Error_string( int errorcode, char *string, int *resultlen )
*/
   /* MPI_Error_string(*err,errorstring,&len); */
   MPI_Error_string((int) *err,  errorstring,  (int *)len);

   fprintf(stdout,"Error string         %s\n",errorstring);
   fprintf(stdout,"******************************************\n");

#if defined (SYS_CRAY) || defined (SYS_T3D)
   QUIT(errorstring,len);
#else
   quit_(errorstring);
#endif

   free(errorstring);

#endif /* VAR_MPI */
}


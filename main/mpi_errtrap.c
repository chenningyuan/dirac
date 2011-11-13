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
 * VERSION : $Revision: 1.15.12.1 $
 * DATE    : $Date: 2004/09/01 08:13:53 $
 * FILE    : mpi_errtrap.c
 * DESCR.  : Trap MPI errors
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdio.h>
#include <string.h>
#if defined (VAR_MPI)
#include <mpi.h>
#endif

#include "mpi_errtrap.head"
 ;


#include "mpi_errtrap.head"

{

#if defined (VAR_MPI)

   int len;
   char *errorstring;

#if defined (SYS_CRAY) || defined (SYS_T3D)
    FLSHFO(stdout);
#else
   flshfo_(stdout);
#endif

   fprintf(stdout,"******************************************\n");
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"Comm                 %d\n",*comm);
   fprintf(stdout,"The PID              %d\n",getpid());
   len = 100;
   errorstring = (char *) diracmalloc(len);
   MPI_Error_string(*err,errorstring,&len);
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

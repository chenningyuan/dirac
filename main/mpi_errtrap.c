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

   int len;
   char *errorstring;

#if defined (VAR_MPI)

#if defined (SYS_AIX)

/*   flshfo(stdout); */

   fprintf(stdout,"******************************************\n");
   fprintf(stdout,"MPI ERROR: error no. %d\n",*err);
   fprintf(stdout,"Comm                 %d\n",*comm);
   fprintf(stdout,"Function             %s\n",errmesg);
   fprintf(stdout,"PID                  %d\n",getpid());
   len = 100;
   errorstring = (char *) diracmalloc(len);
   MPI_Error_string(*err,errorstring,&len);
   fprintf(stdout,"Error string         %s\n",errorstring);
   fprintf(stdout,"******************************************\n");
   quit(errorstring,len);
   free(errorstring);

#else /* SYS_AIX */

#if defined (SYS_CRAY) || defined (SYS_T3D)
    FLSHFO(stdout);
#elif defined (SYS_AIX)  || defined (SYS_HPUX)
   flshfo(stdout);
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
#elif defined (SYS_AIX)  || defined (SYS_HPUX)
   quit(errorstring,len);
#else
   quit_(errorstring);
#endif

   free(errorstring);

#endif /* SYS_AIX */

#endif /* VAR_MPI */
}

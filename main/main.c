/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Main program for dirac
 * AUTHOR  : J. Thyssen - 980709
 *
 */

#include <stdio.h>

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

void main()

{

   int mytid, size, parcal, mparid, numnod, errhandler;
   int *mempointer;
#if defined (SYS_CRAY) || defined (SYS_T3D)
   real cpu_start, cpu_end, wall_start, wall_end;
#else
   double cpu_start, cpu_end, wall_start, wall_end;
#endif
   /*
    * Start timing
    */

   gettim_(&cpu_start,&wall_start);

   /*
    * Initialize MPI stuff
    */

   ini_mpi(&parcal, &mytid, &mparid, &numnod, &errhandler);

   /*
    * Allocate memory
    */

   mempointer = alloc_mem(mytid,&size);

   /*
    * Print logo, sysinfo etc.
    */

   print_logo(mytid);
  
   /* sysinf(); */
 

   /*
    * Execute dirac
    */

   exec_dirac(&parcal,&mytid,&mparid,&numnod,mempointer,&size);

   /*
    * Free memory
    */

   free(mempointer);

   /*
    * End timing
    */

   gettim_(&cpu_end,&wall_end);

   /*
    * Write out timing
    */

   pr_timing(cpu_end - cpu_start, wall_end - wall_start, numnod, mytid);

   /*
    * End MPI
    */

   end_mpi(mparid,mytid,errhandler);
}

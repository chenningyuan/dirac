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
#include "get_time.struct"
#include "get_time.head"
 ;

void main()

{

   int mytid, size, parcal, mparid, numnod, errhandler;
   int *mempointer;
   get_time_struct *gt_first, *gt_last;

   /*
    * Start timing
    */

   gt_first = (get_time_struct *) diracmalloc(sizeof(get_time_struct));
   gt_last  = (get_time_struct *) diracmalloc(sizeof(get_time_struct));
   get_time(gt_first);

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
  
   sysinf();
 

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

   end_mpi(mparid,mytid,errhandler);

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

#ifdef VAR_G77
/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * FILE    : dirac_linux.c 
 * SOURCE  : $Source$
 * PURPOSE : Allocate memory for Dirac
 *           and execute Dirac.
 * CPRGHT  : (c) Jan 8 1998 by Joern Thyssen
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void main()

{
   double *work;
   int lwork;
   char *dirwrk,*temp;
   int str_dirwrk[2];

   temp = malloc(80);
   dirwrk = temp;
   printf("DIRAC: \n\tDefault work memory:        %10i words\n",
      INSTALL_WRKMEM);
   dirwrk = getenv("DIRWRK");
   if ( dirwrk == NULL ) 
   {
      dirwrk = temp;
      strcpy(dirwrk,"undefined");
      lwork = INSTALL_WRKMEM;
   }
   else
   {
      lwork = atoi(dirwrk);
      printf("DIRAC: \n\tDIRWRK enviroment variable: %10i words\n",
         lwork);
   }
   printf("DIRAC: \n\tAllocating:                 %10i words\n",
      lwork);
   (int) work = malloc(8*lwork);
   printf("work %i %i \n",work,lwork);
   if ( work == NULL )
   {
      fprintf(stderr,"\n\n  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---\n");
      fprintf(stderr," MEMORY ALLOCATION ERROR.\n\n");
      return;
   }

   str_dirwrk[0] = (int) dirwrk;
   str_dirwrk[1] = strlen(dirwrk);
   *work = 10101.0;
   dirac_(work,&lwork,str_dirwrk[0],str_dirwrk[1]);
}
#endif

#ifdef VAR_G77
/*
 *
 * VERSION : $Version$
 * DATE    : $Date$
 * FILE    : dirac_linux.c ($RCSfile$)
 * PURPOSE : Allocate memory for Dirac
 *          and execute Dirac.
 * CPRGHT  : (c) Jan 8 1998 by Joern Thyssen
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/utsname.h>
#include <string.h>

void main()

{
   double *work;
   int lwork;
   char *dirwrk,*temp;
   struct utsname *un;
   int str_dirwrk[2];

   un = malloc(sizeof(struct utsname));
   if (uname(un) != 0)  {
      printf("uname error...(not important)\n");
      strcpy(un->nodename,"unknown");
      strcpy(un->nodename,"unknown");
      strcpy(un->release,"unknown");
      strcpy(un->sysname,"unknown");
      strcpy(un->machine,"unknown");
   }
   temp = malloc(80);
   dirwrk = temp;
   printf("LAUNCHING DIRAC on %s (Linux release %s)\n",
      un->nodename,un->release);
   printf("DIRAC(%s %s): \n\tDefault work memory:        %10i words\n",
      un->sysname,un->machine,INSTALL_WRKMEM);
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
      printf("DIRAC(%s %s): \n\tDIRWRK enviroment variable: %10i words\n",
         un->sysname,un->machine,lwork);
   }
   printf("DIRAC(%s %s): \n\tAllocating:                 %10i words\n",
      un->sysname,un->machine,lwork);
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

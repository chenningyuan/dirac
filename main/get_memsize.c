/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * SOURCE  : $Source$
 * DESCR.  : Get memory size from enviroment variables if existing
 * AUTHOR  : J. Thyssen - 980713
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "get_memsize.head"
 ;
#include "diracmalloc.head"
 ;


#include "get_memsize.head"

{
   char *strwrk,*str;
   int lwrkmm,ierr,returncode;

   strwrk = (char *) diracmalloc(20);
   str    = (char *) diracmalloc(20);
   if (mytid == 0)
   {
      strcpy(str,"DIRWRK");
   }
   else
   {
      strcpy(str,"DIRWRK");
   }

   strwrk = getenv(str);
   returncode = 0;
   if (strwrk != NULL)
   {
      *size = atoi(strwrk);
      if ( *size != 0 )
      {
         returncode = 1;
      }
   }
   free(strwrk);
   free(str);
   return returncode;
}





/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * FILE    : gpjth_c.c
 * SOURCE  : $Source$
 *
 * frac2i:
 *
 * Take string with a fraction (e.g. either 1/2 or 0.5)
 * Written by J. Thyssen - Jan 28 1998
 *
 *
 * sysinf:
 *
 * Print info about system to stdout
 * Written by J. Thyssen - Feb 9 1998
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/utsname.h>
#include "patchlevel.h"

#if defined (SYS_CRAY) || defined (SYS_T3D)
#include <fortran.h>
#define frac2i FRAC2I
#define sysinf SYSINF
/*
 * These HPUX flags does not work
 *
 * #elif defined (SYS_HPUX)
 * #define frac2i FRAC2I
 * #define sysinf SYSINF
 *
 */
#elif defined (SYS_AIX)
#define frac2i frac2i
#define sysinf sysinf
#else
#define frac2i frac2i_
#define sysinf sysinf_
#endif /*CRAY*/

void c_frac2i();
void frac2i();
void sysinf();


/* 
 * Take string s and return the fraction n/d
 */

void c_frac2i(s,f) 
char *s;
double *f;

{
   char *s1;
   int d,n;

   s1 = strstr(s,"/");
   /*
    * If no '/' present in string then it is a float
    * else it is a 'n/d'
    */
   if ( s1 == (char *) NULL ) 
   {
      *f = atof(s);
   }
   else
   {
      s1++;
      d = atoi(s1);
      n = atoi(s);
      *f = ((double) n) / ((double) d);
   }
}

void
#if defined(SYS_CRAY) || defined (SYS_T3D)
frac2i(f,fcd)
double *f;
_fcd fcd;
#else
frac2i(f,s,l)
double *f;
char *s;
int l;
#endif

{
   char *temp;
#if defined(SYS_CRAY) || defined (SYS_T3D)
   int l;
   char *s;

   l = _fcdlen(fcd);
   s = _fcdtocp(fcd);
#endif
   temp = (char *) malloc(l+1);
   strncpy(temp, s, l);
   temp[l] = '\0';
   c_frac2i(temp,f);
   free(temp);
}


void sysinf()

{
   struct utsname *un;

   un = (struct utsname *) malloc(sizeof(struct utsname));

   if ( un != NULL )
   {
      uname(un);
   }

   printf( " Version           : %s %i.%i (patchlevel %i) \n" , \
      VER_STR,MAJOR_REV,MINOR_REV,PATCHLEVEL);
   if ( un != NULL )
   {
      printf( " Hostname          : %s\n" , un->nodename);
      printf( " Operating system  : %s %s\n" , un->sysname, un->release);
      printf( " Machine           : %s\n" , un->machine);
   }

   free(un);
}


/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * FILE    : gpjth_c.c
 * SOURCE  : $Source$
 *
 * Take string with a fraction (e.g. either 1/2 or 0.5)
 * Written by J. Thyssen - Jan 28 1998
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined (SYS_CRAY)
#include <fortran.h>
#define frac2i FRAC2I
#elif defined (SYS_AIX) || defined (SYS_HPUX)
#define frac2i frac2i
#else
#define frac2i frac2i_
#endif /*CRAY*/

void c_frac2i();
void frac2i();


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
#if defined(SYS_CRAY)
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
#if defined(_CRAY)
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

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
  *...   group of the licence holder may not take any part of Dalton,
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
#define c_readreorder C_READREORDER
#elif defined (SYS_HPUX)
/* nothing */
#elif defined (SYS_AIX)
#define frac2i frac2i
#define c_readreorder c_readreorder
#else
#define frac2i frac2i_
#define c_readreorder c_readreorder_
#endif /*CRAY*/

void c_frac2i();
void frac2i();
void sysinf();
void c_readreorder();

/*
 * Take string of tokens and return new ordering
 */

void c_readreorder(char *string, int *ireord, int imoord[], int *mxreord)

{
   char *token;
   int itemp;

   *ireord = 0;
   
   /*
   printf("string: %s\n",string);
    */
   token = strtok(string,",");

   while ( token != NULL )
   {
      /*
      printf("Token %2i: %s\n",*ireord,token);
       */
      itemp = atoi(token);
      if ( itemp <= 0 ) return;
      if ( (*ireord = *ireord + 1) > (*mxreord) ) 
      {
         *ireord = -1;
         return;
      }
      imoord[(*ireord)-1] = atoi(token);
      token = strtok(NULL,",");
   }
}

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


void sysinf(int parcal)

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
   if ( ! parcal )
      printf( " Run               : Sequential\n");
   else
      printf( " Run               : Parallel\n");
   free(un);
}


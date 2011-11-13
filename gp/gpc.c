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
 * VERSION : $Revision: 1.9.2.1 $
 * DATE    : $Date: 2004/09/01 08:10:54 $
 * FILE    : gpc.c
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
#include "compinfo.h"

#if defined (SYS_CRAY) || defined (SYS_T3D)
#include <fortran.h>
#define frac2i FRAC2I
#define creadreorder CREADREORDER
#elif defined (SYS_HPUX)
/* nothing */
#elif defined (SYS_AIX)
#define frac2i frac2i
#define creadreorder creadreorder
#else
#define frac2i frac2i_
#define creadreorder creadreorder_
#define compoff compoff_
#endif /*CRAY*/

void c_frac2i();
void frac2i();
/* void sysinf(); */
void sysinf(int );
/* void creadreorder(); */
void creadreorder(char *, int *, int *, int *);
/* void distribute_env(); */
void distribute_env(const char*, int);

/*
 * Take string of tokens and return new ordering
 */

void creadreorder(char *string, int *ireord, int imoord[], int *mxreord)

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

   /* Information about compilation */

   printf(" Last compilation  : %s %s (by %s@%s)\n", 
      __TIME__,__DATE__,WHO_COMPILED,WHERE_COMPILED);

   free(un);
}

#if defined (IA64)
void compoff_(double *wrk,double *work, long *kbase) {
  *kbase = (long) (wrk-work+1);
/*
 *printf( "compoff: wrk=%p work=%p diff=%p\n", wrk, work, wrk-work+1);
 *printf( "kbase=%d\n", *kbase );
 *printf( "int %d long %d short int %d\n", sizeof(int), sizeof(long),
 *sizeof(short int));
*/
}
#else
void compoff_(double *wrk,double *work, int *kbase) {
  *kbase = (int) (wrk-work+1);
}
#endif

void distribute_env(const char* name, int mytid)
{

 /* Routine to set the value of environment variables on the slaves
    equal to that of the master if it was not defined on the slave
    Written by Luuk Visscher, october 2002                          */

#if !defined (VAR_MPI)
  return;
#else
#include <mpi.h>
  char *env_local, *env_master, *newname;
  int len_local, len_master, len_name;

  len_local = 0;
  if( (env_local=getenv(name)) != NULL) len_local=strlen(env_local)+1;

  if ( mytid == 0 ) len_master = len_local; 
  MPI_Bcast(&len_master, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if ( mytid == 0 ) 
     env_master = env_local; 
  else
     env_master = (char *) malloc(len_master);

  MPI_Bcast(env_master, len_master, MPI_CHAR, 0, MPI_COMM_WORLD);

  if ( mytid != 0 && len_local == 0 && len_master != 0) 
  { 
    len_name = strlen(name) + len_master + 1;
    newname = (char *) malloc(len_name);

    /*
    printf("Environment variable %s not set on node %i taking value %s from master\n",
            name,mytid,env_master);
    */
    sprintf(newname,"%s=%s",name,env_master);
    putenv(newname);
    free(env_master);
  }

#endif
}

int icommsize_(const int *istart, const int *iend)
{
  return iend-istart+1;
}


#ifdef VAR_IFC
 /* extra routines needed by intel's ifc compiler which does not
 * provide them itself.
 * Pawel Salek, pawsa@theochem.kth.se, 2002.01.31
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <time.h>
#include <unistd.h>

void unlink_(char * fname, int len)
{
  int i;
  for(i=len; i>0 && !isspace(fname[i]); i--)
    ;
  fname[i] = '\0';
  unlink(fname);
}

void getenv_(const char* str, char* val, int strl, int vall)
{
  char name[32];
  char *envv;
  int i;

  if(strl>sizeof(name)-1) {
    printf("oops\n");
    return;
  }
  strncpy(name, str, strl);
  name[strl] = '\0';

  /*printf("getenv(%s) called with strl=%d, vall=%d '%s'\n",
    name, strl, vall, getenv(name));*/
  if( (envv=getenv(name)) != NULL)
    strncpy(val, envv, vall);
  else val[0] = '\0';
  /*printf("getenv(%s) returned '%s' (%d, %d)\n", name, val, strl, vall);*/
  for(i=strlen(val); i< vall; i++)
    val[i] = ' ';
}

void system_(char* str, int strlen)
{
  printf("Bloody crayio. I ignore your commands.\n");
}

void flush_(int* lu)
{
  /* no flushing */
}

float etime_(float* et)
{
  static int ticksclk = -1;
  struct tms tm;

  times(&tm);
  if(ticksclk<0) ticksclk = sysconf(_SC_CLK_TCK);  
  et[0] = tm.tms_utime/(double)ticksclk;
  et[1] = tm.tms_stime/(double)ticksclk;
  return et[0] + et[1];
}

void fdate_(char* dt, int dtlen)
{
  time_t curr_time;
  time(&curr_time);
  strncpy(dt, ctime(&curr_time), 24);
}

int time_(void)
{ 
  return (int)time(NULL);
}

void hostnm_(char* buf, int buflen)
{
  int i;
  gethostname(buf, buflen);
  /* printf("Hostname: '%s' (%d)\n", buf, buflen); */
  for(i=strlen(buf); i<buflen; i++)
    buf[i] = ' ';
}

double* malloc_(int cnt)
{
  return malloc(cnt);
}
#endif

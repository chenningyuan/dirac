/*
!
!      Copyright (c) 2008 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
*/
/*
 *
 * VERSION : $Revision: 7750 $
 * DATE    : $Date: 2008-09-17 13:58:49 +0200 (Wed, 17 Sep 2008) $
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
#include <sys/time.h>
#include <sys/resource.h>
#include "patchlevel.h"
#include "compinfo.h"


#ifndef  NO_C99
#include <stdint.h>
typedef int64_t int_star_8_t;
#else
typedef long int_star_8_t; /* This should be long long on many machines. */
#endif


/*
 * Take string of tokens and return new ordering
 */

void creadreorder(char *string, int *ireord, int imoord[], int *mxreord)

{
   char *token;
   int itemp;

   *ireord = 0;
   
   token = strtok(string,",");

   while ( token != NULL )
   {
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

static void print_time(time_t secs, FILE *dst)
{
    if (secs > 60*60*24)
    {
	int days = secs / (60*60*24);
	fprintf(dst,"%idays",days);
	secs -= days*60*60*24;
    }
    if (secs > 60*60)
    {
	int hrs = secs / (60*60);
	fprintf(dst,"%ih",hrs);
	secs -= hrs*60*60;
    }
    if (secs > 60)
    {
	int mins = secs / (60);
	fprintf(dst,"%imin",mins);
	secs -= mins*60;
    }
    fprintf(dst,"%is",(int)secs);
}

/* Start timers, or stop and print cpu and wall time*/
void stopwatch_(int *start, int *node_nr)
{
    static struct timeval start_wall;
    if (*start)
    {
	gettimeofday(&start_wall, NULL);
    }
    else
    {
	struct rusage rusage;
	if (getrusage(RUSAGE_SELF,&rusage) == 0)
	{
	    /* These are the fields supported by linux, we may add others on other systems. See the man page of getrusage*/
	    printf(">>>> Node %i, utime: %lld, stime: %lld, minflt: %ld, majflt: %ld, nvcsw: %ld, nivcsw: %ld\n",
		   *node_nr, (long long int)(rusage.ru_utime.tv_sec), (long long int)(rusage.ru_stime.tv_sec), 
		   rusage.ru_minflt, rusage.ru_majflt, rusage.ru_nvcsw, rusage.ru_nivcsw, rusage.ru_maxrss);
	}
	if (*node_nr == 0) /* print wall time only on the master */
	{
	    struct timeval end_wall;
	    gettimeofday(&end_wall, NULL);
	    printf(">>>> Total WALL time used in DIRAC: ");
	    print_time(end_wall.tv_sec - start_wall.tv_sec,stdout);
	    printf("\n");
	}
    }
}

void sysinf_(int *parcal)

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
   if ( ! *parcal )
      printf( " Run               : Sequential\n");
   else
      printf( " Run               : Parallel\n");

   /* Information about compilation */

   printf(" Last compilation  : %s %s (by %s@%s)\n", 
      __TIME__,__DATE__,WHO_COMPILED,WHERE_COMPILED);

   printf(" FORTRAN compiler  : %s \n", FCOMP);
   printf(" C compiler        : %s \n", CCOMP);
   printf(" Changeset         : %s \n", VERSION);
   fflush(stdout);
   free(un);
}

#define quit        quit_
    /* NOTE: if this does not link with the FORTRAN subroutine QUIT,
     * then try quit without the underscore */
extern void quit(const char* str, int len);

void compoff_(double *wrk, double *work, int_star_8_t *k_offset) {
  *k_offset = wrk-work+1;
/* Aug 2008: compoff_ shall always return k_offset as a FORTRAN INTEGER*8
 * so that we don't get nonsense if address difference is bigger than 31
 * bit. /HJAaJ
 *
*/
  if ( sizeof(*k_offset) != 8 )
  {
    printf( "ERROR in COMPOFF in gpc.c: long is not I*8");
    printf( "compoff: wrk=%p work=%p diff=%d\n", wrk, work, wrk-work+1);
    printf( "k_offset=%d\n", *k_offset );
    printf( "int %d long %d short int %d int_star_8_t %d\n", sizeof(int), sizeof(long), sizeof(short int),sizeof(int_star_8_t));
    quit( "ERROR in COMPOFF in gpc.c: long is not I*8",42);
    exit(1);
  }
/*
 *printf( "compoff: wrk=%p work=%p diff=%p\n", wrk, work, wrk-work+1);
 *printf( "k_offset=%d\n", *k_offset );
 *printf( "int %d long %d short int %d\n", sizeof(int), sizeof(long),
 *sizeof(short int));
*/
}

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


/*
  Workaround for buggy fortran SYSTEM() on some
  platforms (intel fortran 8.x-9.x)
*/
int csystem_(const char *cmd, int cmdlen)
{
  int i;
  char *s;
  s = malloc(cmdlen+1);
  for (i=0;i<cmdlen;i++)
    s[i] = cmd[i];
  s[cmdlen] = 0;
  i = system(cmd);
  free(s);
  return i;
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
    printf("gpc.c:getenv(): input string too long (I suck!)\n");
    return;
  }
  strncpy(name, str, strl);
  name[strl] = '\0';

  if( (envv=getenv(name)) != NULL)
    strncpy(val, envv, vall);
  else val[0] = '\0';

  for(i=strlen(val); i< vall; i++)
    val[i] = ' ';
}

void system_(char* str, int len)
{    
    char *buf = alloca(len+1);
    buf[len] = 0;
    strncpy(buf,str,len);
    printf("VAR_IFC defined. Ignoring call system(%s)\n",buf);
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
#endif /* VAR_IFC */

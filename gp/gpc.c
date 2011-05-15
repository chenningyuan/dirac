/* dirac_copyright_start */
/*
 *     Copyright (c) 2010 by the authors of DIRAC.
 *     All Rights Reserved.
 *
 *     This source code is part of the DIRAC program package.
 *     It is provided under a written license and may be used,
 *     copied, transmitted, or stored only in accordance to the
 *     conditions of that written license.
 *
 *     In particular, no part of the source code or compiled modules may
 *     be distributed outside the research group of the license holder.
 *     This means also that persons (e.g. post-docs) leaving the research
 *     group of the license holder may not take any part of Dirac,
 *     including modified files, with him/her, unless that person has
 *     obtained his/her own license.
 *
 *     For information on how to get a license, as well as the
 *     author list and the complete list of contributors to the
 *     DIRAC program, see: http://dirac.chem.vu.nl
 */
/* dirac_copyright_end */

/*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
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

#if !defined (SYS_WINDOWS)
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <math.h>
#include "patchlevel.h"
#include "compinfo.h"
#include "fortran_interface.h"

void printaddr_(void *p)
{
   fprintf(stderr,"Address is %p\n",p);
}


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
void stopwatch_(fortran_int_t *start, fortran_int_t *node_nr)
{
#if !defined (SYS_WINDOWS)
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
         /* miro: deactivated, because often it's beeing printed among master data, hurting thus output

	    printf(">>>> Node %i, utime: %lld, stime: %lld, minflt: %ld, majflt: %ld, nvcsw: %ld, nivcsw: %ld\n",
		   *node_nr, (long long int)(rusage.ru_utime.tv_sec), (long long int)(rusage.ru_stime.tv_sec), 
		   rusage.ru_minflt, rusage.ru_majflt, rusage.ru_nvcsw, rusage.ru_nivcsw, rusage.ru_maxrss);
          */
	}

	if (*node_nr == 0) /* print wall time only on the master */
	{
	    struct timeval end_wall;
	    gettimeofday(&end_wall, NULL);

	    printf(">>>> Node %i, utime: %lld, stime: %lld, minflt: %ld, majflt: %ld, nvcsw: %ld, nivcsw: %ld, maxrss: %ld\n",
		   (int)*node_nr, (long long int)(rusage.ru_utime.tv_sec), (long long int)(rusage.ru_stime.tv_sec), 
		   rusage.ru_minflt, rusage.ru_majflt, rusage.ru_nvcsw, rusage.ru_nivcsw, rusage.ru_maxrss);

	    printf(">>>> Total WALL time used in DIRAC: ");
	    print_time(end_wall.tv_sec - start_wall.tv_sec,stdout);
	    printf("\n");
	}
    }
#endif
}

void sysinf_(fortran_int_t *parcal)
{

#ifdef THIS_IS_CMAKE_BUILD
   printf("\n");
   printf("\n");
   printf("Compilation settings\n");
   printf("--------------------\n");
   printf("\n");
   printf("CMAKE_SYSTEM_NAME           : %s \n", CMAKE_SYSTEM_NAME           );
   printf("\n");
   printf("ENABLE_EXTERNAL_MATH        : %s \n", ENABLE_EXTERNAL_MATH        );
   printf("ENABLE_64BIT_INTEGERS       : %s \n", ENABLE_64BIT_INTEGERS       );
   printf("ENABLE_MPI                  : %s \n", ENABLE_MPI                  );
   printf("ENABLE_XCFUN                : %s \n", ENABLE_XCFUN                );
   printf("ENABLE_LOCAL_DISKS          : %s \n", ENABLE_LOCAL_DISKS          );
   printf("NO_SSH_ACCESS_TO_NODES      : %s \n", NO_SSH_ACCESS_TO_NODES      );
   printf("\n");
   printf("CMAKE_BUILD_TYPE            : %s \n", CMAKE_BUILD_TYPE            );
   printf("\n");
   printf("CMAKE_Fortran_COMPILER      : %s \n", CMAKE_Fortran_COMPILER      );
   printf("FORTRAN_COMPILER_VERSION    : %s \n", FORTRAN_COMPILER_VERSION    );
   printf("CMAKE_Fortran_FLAGS         : %s \n", CMAKE_Fortran_FLAGS         );
   printf("CMAKE_Fortran_FLAGS_DEBUG   : %s \n", CMAKE_Fortran_FLAGS_DEBUG   );
   printf("CMAKE_Fortran_FLAGS_RELEASE : %s \n", CMAKE_Fortran_FLAGS_RELEASE );
   printf("\n");
   printf("CMAKE_C_COMPILER            : %s \n", CMAKE_C_COMPILER            );
   printf("C_COMPILER_VERSION          : %s \n", C_COMPILER_VERSION          );
   printf("CMAKE_C_FLAGS               : %s \n", CMAKE_C_FLAGS               );
   printf("CMAKE_C_FLAGS_DEBUG         : %s \n", CMAKE_C_FLAGS_DEBUG         );
   printf("CMAKE_C_FLAGS_RELEASE       : %s \n", CMAKE_C_FLAGS_RELEASE       );
   printf("\n");
   printf("BLAS_LIBRARIES              : %s \n", BLAS_LIBRARIES              );
   printf("LAPACK_LIBRARIES            : %s \n", LAPACK_LIBRARIES            );
   printf("\n");
   printf("SVN_REVISION                : %s \n", SVN_REVISION                );
   printf("GIT_REVISION                : %s \n", GIT_REVISION                );

#else

#if !defined (SYS_WINDOWS)
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
#endif /* if !defined (SYS_WINDOWS) */
   if ( ! *parcal )
      printf( " Run               : Sequential\n");
   else
      printf( " Run               : Parallel\n");

   printf(" Last compilation  : %s %s (by %s@%s:%s)\n", 
      __TIME__,__DATE__,WHO_COMPILED,WHERE_COMPILED, COMPILED_DIRAC_DIR);

   printf(" FORTRAN compiler  : %s \n", FCOMP);
   printf(" C compiler        : %s \n", CCOMP);
   printf(" External libraries: %s \n", LIBRARIES);
   printf(" Internal libraries: DIRAC_BLAS=%s  DIRAC_LAPACK=%s \n", DIRAC_BLAS, DIRAC_LAPACK);
   printf(" ILP data model    : %s \n", ILP_DATA_MODEL);
   printf(" SVN URL           : %s \n", DIRAC_SVN_URL);
   printf(" SVN Revision      : %s \n", DIRAC_SVN_REVISION);

#endif /* #ifdef THIS_IS_CMAKE_BUILD */
   fflush(stdout);
}

#define quit        quit_
    /* NOTE: if this does not link with the FORTRAN subroutine QUIT,
     * then try quit without the underscore */
extern void quit(const char* str, fortran_int_t len);

void compoff_(double *wrk, double *work, fortran_int8_t *k_offset) {
  *k_offset = wrk-work+1;
/* Aug 2008: compoff_ shall always return k_offset as a FORTRAN INTEGER*8
 * so that we don't get nonsense if address difference is bigger than 31
 * bit. /HJAaJ
 *
*/
  if ( sizeof(*k_offset) != 8 )
  {
    printf( "ERROR in COMPOFF in gpc.c: long is not I*8");
    printf( "compoff: wrk=%p work=%p diff=%li\n", wrk, work, (long)(wrk-work+1));
    printf( "k_offset=%li\n", (long)*k_offset );
    printf( "int %i long %i short int %i fortran_int8_t %i\n", (int)sizeof(int), (int)sizeof(long), (int)sizeof(short int),(int)sizeof(fortran_int8_t));
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

int icommsize_(const fortran_int_t *istart, const fortran_int_t *iend)
{
  return iend-istart+1;
}


/*
  Workaround for buggy fortran SYSTEM() on some
  platforms (intel fortran 8.x-9.x)
*/
int csystem_(const char *cmd, fortran_int_t cmdlen)
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

void unlink_(char * fname, fortran_int_t len)
{
  int i;
  for(i=len; i>0 && !isspace(fname[i]); i--)
    ;
  fname[i] = '\0';
  unlink(fname);
}

void getenv_(const char* str, char* val, fortran_int_t strl, fortran_int_t vall)
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

void system_(char* str, fortran_int_t len)
{    
    char *buf = alloca(len+1);
    buf[len] = 0;
    strncpy(buf,str,len);
    printf("VAR_IFC defined. Ignoring call system(%s)\n",buf);
}

void flush_(fortran_int_t* lu)
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

void fdate_(char* dt, fortran_int_t dtlen)
{
  time_t curr_time;
  time(&curr_time);
  strncpy(dt, ctime(&curr_time), 24);
}

int time_(void)
{ 
  return (int)time(NULL);
}

void hostnm_(char* buf, fortran_int_t buflen)
{
  int i;
  gethostname(buf, buflen);
  /* printf("Hostname: '%s' (%d)\n", buf, buflen); */
  for(i=strlen(buf); i<buflen; i++)
    buf[i] = ' ';
}

double* malloc_(fortran_int_t cnt)
{
  return malloc(cnt);
}
#endif /* VAR_IFC */

/* Try to cause a segmentation fault, to enable debuggers to
   trap this and inspect the program.
*/
void crash_(void)
{
  *(int *)0 = 0;
}

/*
 Often  we need to detect "ill" real numbers - NaN, +Inf, -Inf.
 For that we use C <math.h> functions isnan(x), isinf(x).
 since not all Fortran compilers have these, IMHO(MI) functions.
*/
int is_ill_number_(double x)
{
if (isnan(x) != 0 || isinf(x) ==-1 || isinf(x) ==+1)
{
 return 1;
}
else{return 0;}
}


static double *mem_start = 0;
static double *mem_high = 0;

void set_mem_start_(double *here)
{
  mem_start = here;
  mem_high = here;
}

void push_water_mark_(double *here)
{
  if (mem_high < here)
    mem_high = here;
}

int high_water_mark_(void)
{
  return mem_high - mem_start;
}

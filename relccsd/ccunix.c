/* Utility routines for molpro88 unix version
	Coded by U. Welz, Universitaet Bielefeld, January 1988
	Further information from Peter Knowles or Joachim Werner

I/O routines:	openc(unit,fname,size,status)
		character*(*) fname
		integer unit,size,status
			open file fname with unit number unit
			integer status codes as defined below

		closec(unit)

		rdabsf(unit,a,l,p)
		wrabsf(unit,a,l,p)
		integer unit,l,p
		double precision a
			read,write respectively l words on unit with buffer a at
			offset p words relative to beginning of file
			all counting done in double words

Timing routines:second()
			return cpu time in seconds
		timing(cpu,sys,io)
		double precision cpu,sys,io
			return cpu, system, i/o time (all seconds)
			the i/o time is computed by couning disk transfers in
			rdabsf, wrabsf assuming SEEK and SPEED parameters given
			below.

Bit routines:	popcnt(pat)
		leadz(pat)
			as documented for cray, except operating on 32 bit word
	*/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <sys/times.h>
#include <sys/param.h>

extern char  *getenv(); 
#define SEEK		0.015		/* average seektime in seconds */
#define SPEED		307250		/* speed in words per second */

#define MAXLENGTH	15
#define MAXUNIT		99

#if defined (SYS_CRAY)
#define WORT	        8	
#else
#define WORT	        4	
#endif

#define PATHMAX         512 
#define UNKNOWN		1
#define SCRATCH		0
#define NEW		2
#define OLD		3

typedef struct
{ long  fd, addr, size;
  char fname[15];
} FILE_DEFINITION;

FILE_DEFINITION files[MAXUNIT];
int nwords=0, nio=0;

#if defined (SYS_IRIX)
void openc_(unit,fname,size,status) long *unit, *size, *status; char *fname;
#endif
#if defined (SYS_CRAY) || defined (SYS_AIX)
void openc(unit,fname,size,status) long *unit, *size, *status; char *fname;
#endif

{ char *cp, *env ;
  char  name[PATHMAX]; 


  if (*unit>MAXUNIT || *unit<0)
	{ fprintf(stderr,"openc: Unit number out of range (%d)\n",*unit);
	  return;
	}
  cp=fname+strlen(fname)-1;
  if (*cp==' ') { while (*cp--==' '); *(cp+2)=NULL; }
  if (strlen(fname)-1>MAXLENGTH) fprintf(stderr,"openc: filename too long >15\n");
       /*  Get alternate file name from environment, 
           if it exists, use it, if not, use program supplied name 
        */  
  strcpy(name,fname); 
  if((env=getenv(fname)) !=NULL) 
  strcpy(name,env); 

  switch(*status)
{ case SCRATCH:	sprintf(files[*unit].fname,"Tmp%d",getpid());
		files[*unit].fd=open(name,O_RDWR|O_CREAT,0666);
		unlink(files[*unit].fname);
		break;
  case UNKNOWN:	if ((files[*unit].fd=open(name,O_RDWR|O_CREAT,0666))==-1)
		{ fprintf(stderr,"openc: Error in opening file %s\n",fname);
				  exit(1);
				}
		break;
  case NEW:	if ((files[*unit].fd=open(name,O_RDWR|O_CREAT|O_TRUNC,0666))==-1)
			{ fprintf(stderr,"openc: Error in opening file %s\n",fname);
			  exit(1);
			}
		break;
  case OLD:	if ((files[*unit].fd=open(name,O_RDWR))==-1)
			{ fprintf(stderr,"openc: Error in opening file %s\n",fname);
			  exit(1);
			}
		break;
  default:	fprintf(stderr,"openc: Unknown status\n"); exit(1);
}
  files[*unit].size= *size;
  strncpy(files[*unit].fname,fname,14); *(files[*unit].fname+14)=NULL;
  *size=(lseek(files[*unit].fd,0L,2)+511)/512;
  files[*unit].addr= -1;
}

#if defined (SYS_IRIX)
void wrabsf_(unit,a,l,p) long *unit, *l, *p; char *a;
#endif
#if defined (SYS_CRAY) || (SYS_AIX)
void wrabsf(unit,a,l,p) long *unit, *l, *p; char *a;
#endif

{ long addr, m, n;

  if (*unit>MAXUNIT || *unit<0)
	{ fprintf(stderr,"wrabs: Unit number out of range (%d)\n",*unit);
	  return;
	}
  if (!*files[*unit].fname)
	{ fprintf(stderr,"wrabs: write without open file unit=%d\n",*unit);
	  return;
	}
  addr= *p * WORT;
  m= *l * WORT;
  if (addr!=files[*unit].addr)
	if (lseek(files[*unit].fd,addr,0)==-1)
	{ fprintf(stderr,"wrabs: Error in lseek of (%d:%s) -> pointer=%d\n",*unit,files[*unit].fname,*p);
	  files[*unit].addr= -1;
	  return;
	}
  if ((n=write (files[*unit].fd,a,m))!=m)
	fprintf(stderr,"wrabs: Error in writing %d words in file %s with unit %d\n",*p,files[*unit].fname,*unit);
  nio++; nwords+= *l;
  files[*unit].addr=addr+n;
}

#if defined (SYS_IRIX)
void rdabsf_(unit,a,l,p) long *unit, *l, *p; char *a;
#endif
#if defined (SYS_CRAY) || (SYS_AIX)
void rdabsf(unit,a,l,p) long *unit, *l, *p; char *a;
#endif

{ long addr, m, n;

  if (*unit>MAXUNIT || *unit<0)
	{ fprintf(stderr,"rdabs: Unit number out of range (%d)\n",*unit);
	  return;
	}
  if (!*files[*unit].fname)
	{ fprintf(stderr,"rdabs: read without open file unit=%d\n",*unit);
	  return;
	}
  addr= *p * WORT;
  m= *l *WORT;
  if (addr!=files[*unit].addr)
	if (lseek(files[*unit].fd,addr,0)==-1)
	{ fprintf(stderr,"rdabs: Error in lseek of (%d:%s) -> pointer=%d\n",*unit,files[*unit].fname,*p);
	  files[*unit].addr= -1;
	  return;
	}
  if ((n=read (files[*unit].fd,a,m))!=m)
	fprintf(stderr,"rdabs: Error in reading %d words in file %s with unit %d\n",*p,files[*unit].fname,*unit);
  nio++; nwords+= *l;
  files[*unit].addr=addr+n;
}

#if defined (SYS_IRIX)
void closec_(unit) long *unit;
#endif
#if defined (SYS_CRAY) || (SYS_AIX)
void closec(unit) long *unit;
#endif

{ if (*unit>MAXUNIT || *unit<0)
	{ fprintf(stderr,"closec: unit out of range\n"); return; }
  if (files[*unit].fd)
	{ close(files[*unit].fd);
	  files[*unit].fd=files[*unit].addr=files[*unit].size=0;
	  *files[*unit].fname=NULL;
	}
  return;
}


 /*
 *
 * VERSION : $Revision$
 * DATE    : $Date$
 * FILE    : ccunix.c
 * SOURCE  : $Source$
 *
        Coded by U. Welz, Universitaet Bielefeld, January 1988
        Further information from Peter Knowles or Joachim Werner
        Stripped and revised for RELCCSD Luuk Visscher

I/O routines:   openc(unit,fname,size,status)
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

        */
#include <unistd.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#if defined (SYS_T3D) || defined (SYS_CRAY)
#include <fortran.h>
#endif

#include <sys/times.h>
#include <sys/param.h>

extern char  *getenv(); 

#define MAXLENGTH       15
#define MAXUNIT         99

#if defined (SYS_CRAY) || defined (SYS_T3D)
#define WORT            8       
#else
#define WORT            4       
#endif

#define PATHMAX         512 
#define SCRATCH         0
#define UNKNOWN         1
#define NEW             2
#define OLD             3

/* Check that the large file option is enabled, if not use normal lseek */

#if defined SMALL_CFILES
typedef off_t off64_t;
#define lseek64 lseek
#endif

typedef struct
{ int  fd, size;
  off64_t addr;
  char fname[15];
} FILE_DEFINITION;

FILE_DEFINITION files[MAXUNIT];
int nwords=0, nio=0;

#if defined (SYS_CRAY) || defined (SYS_T3D)
#define MACRO_QUIT QUITC
#elif defined (SYS_AIX) || defined (SYS_HPUX)
#define MACRO_QUIT quitc
#else
#define MACRO_QUIT quitc_
#endif

void MACRO_QUIT();

#if defined (SYS_IRIX) || defined (SYS_LINUX) || defined (SYS_DEC) || \
    defined (SYS_SUN) || defined (SYS_DARWIN)
void openc_();
void wrabsf_();
void rdabsf_();
void closec_();
#elif defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_T3D) || defined (SYS_CRAY)
void openc();
void wrabsf();
void rdabsf();
void closec();
#else
????
#endif



#if defined (SYS_IRIX) || defined (SYS_LINUX) || defined (SYS_DEC) || \
    defined (SYS_SUN) || defined (SYS_DARWIN)
void openc_(unit,fname,size,status) off64_t *size; int *unit, *status; char *fname;
#elif defined (SYS_CRAY) || defined (SYS_T3D)
void openc(unit,fcd,size,status) off64_t *size; int *unit, *status; _fcd fcd;
#elif defined (SYS_AIX) || defined (SYS_HPUX)
void openc(unit,fname,size,status) off64_t *size; int *unit, *status; char *fname;
#else
????
#endif

{ char *cp, *env ;
  char  name[PATHMAX]; 
#if defined (SYS_CRAY) || defined (SYS_T3D)
  char *fname;

  fname = _fcdtocp(fcd);
#endif

  if (*unit>MAXUNIT || *unit<0)
        { fprintf(stderr,"openc: Unit number out of range (%d)\n",*unit);
          return;
        }
  cp=fname+strlen(fname)-1;
  if (*cp==' ') { while (*cp--==' '); *(cp+2)=NULL; }
  if (strlen(fname)-1>MAXLENGTH) fprintf(stderr,"openc: filename too int >15\n");
       /*  Get alternate file name from environment, 
           if it exists, use it, if not, use program supplied name 
        */  
  strcpy(name,fname); 
  if((env=getenv(fname)) !=NULL) 
  strcpy(name,env); 

  switch(*status)
{ case SCRATCH: sprintf(files[*unit].fname,"Tmp%d",getpid());
                files[*unit].fd=open(name,O_RDWR|O_CREAT,0666);
                unlink(files[*unit].fname);
                break;
  case UNKNOWN: if ((files[*unit].fd=open(name,O_RDWR|O_CREAT,0666))==-1)
                { fprintf(stderr,"openc: Error in opening file %s\n",fname);
                          MACRO_QUIT();
                          exit(1);
                                }
                break;
  case NEW:     if ((files[*unit].fd=open(name,O_RDWR|O_CREAT|O_TRUNC,0666))==-1)
                        { fprintf(stderr,"openc: Error in opening file %s\n",fname);
                          MACRO_QUIT();
                          exit(1);
                        }
                break;
  case OLD:     if ((files[*unit].fd=open(name,O_RDWR))==-1)
                        { fprintf(stderr,"openc: Error in opening file %s\n",fname);
                          MACRO_QUIT();
                          exit(1);
                        }
                break;
  default:      fprintf(stderr,"openc: Unknown status\n"); 
                          MACRO_QUIT();
                          exit(1);
}
  files[*unit].size= *size;
  strncpy(files[*unit].fname,fname,14); *(files[*unit].fname+14)=NULL;
  *size=(lseek64(files[*unit].fd,0L,2)+511)/512;
  files[*unit].addr= -1;

}

#if defined (SYS_IRIX) || defined (SYS_LINUX) || defined (SYS_DEC) || \
    defined (SYS_SUN) || defined (SYS_DARWIN)
void wrabsf_(unit,a,l,p) off64_t *p; int *unit, *l; char *a;
#elif defined (SYS_CRAY) || defined (SYS_T3D)
void wrabsf(unit,fcd,l,p) off64_t *p; int *unit, *l; _fcd fcd;
#elif defined (SYS_AIX) || defined (SYS_HPUX) 
void wrabsf(unit,a,l,p) off64_t *p; int *unit, *l; char *a;
#else
???
#endif

{ off64_t addr; int  m, n;
#if defined (SYS_CRAY) || defined (SYS_T3D)
  char *a;

  a = _fcdtocp(fcd);
#endif

  if (*unit>MAXUNIT || *unit<0)
        { fprintf(stderr,"wrabs: Unit number out of range (%d)\n",*unit);
          MACRO_QUIT();
          return;
        }
  if (!*files[*unit].fname)
        { fprintf(stderr,"wrabs: write without open file unit=%d\n",*unit);
          MACRO_QUIT();
          return;
        }
  addr= *p * WORT;
  m= *l * WORT;

  if (addr!=files[*unit].addr)
        if (lseek64(files[*unit].fd,addr,0)==-1)
        { fprintf(stderr,"wrabs: Error in lseek64 of (%d:%s) -> pointer=%d\n",*unit,files[*unit].fname,*p);
          files[*unit].addr= -1;
          MACRO_QUIT();
          return;
        }

  if ((n=write (files[*unit].fd,a,m))!=m)
      { fprintf(stderr,"wrabs: Error in writing %d words in file %s with unit %d\n",*p,files[*unit].fname,*unit);
        MACRO_QUIT();
      }
  nio++; nwords+= *l;
  files[*unit].addr=addr+n;
}

#if defined (SYS_IRIX) || defined (SYS_LINUX) || defined (SYS_DEC) || \
    defined (SYS_SUN) || defined (SYS_DARWIN)
void rdabsf_(unit,a,l,p) off64_t *p; int *unit, *l; char *a;
#elif defined (SYS_CRAY) || defined (SYS_T3D)
void rdabsf(unit,fcd,l,p) off64_t *p; int *unit, *l; _fcd fcd;
#elif defined (SYS_AIX) || defined (SYS_HPUX)
void rdabsf(unit,a,l,p) off64_t *p; int *unit, *l; char *a;
#else
????
#endif

{ off64_t addr; int m, n;
#if defined (SYS_CRAY) || defined (SYS_T3D)
  char *a;

  a = _fcdtocp(fcd);
#endif

  if (*unit>MAXUNIT || *unit<0)
        { fprintf(stderr,"rdabs: Unit number out of range (%d)\n",*unit);
          MACRO_QUIT();
          return;
        }
  if (!*files[*unit].fname)
        { fprintf(stderr,"rdabs: read without open file unit=%d\n",*unit);
          MACRO_QUIT();
          return;
        }
  addr= *p * WORT;
  m= *l *WORT;

  if (addr!=files[*unit].addr)
        if (lseek64(files[*unit].fd,addr,0)==-1)
        { fprintf(stderr,"rdabs: Error in lseek64 of (%d:%s) -> pointer=%d\n",*unit,files[*unit].fname,*p);
          files[*unit].addr= -1;
          MACRO_QUIT();
          return;
        }

  if ((n=read (files[*unit].fd,a,m))!=m)
      { fprintf(stderr,"rdabs: Error in reading %d words in file %s with unit %d\n",*p,files[*unit].fname,*unit);
        MACRO_QUIT();
      }
  nio++; nwords+= *l;
  files[*unit].addr=addr+n;
}

#if defined (SYS_IRIX) || defined (SYS_LINUX) || defined (SYS_DEC) || \
    defined (SYS_SUN) || defined (SYS_DARWIN)
void closec_(unit) int *unit;
#elif defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_AIX) || defined (SYS_HPUX)
void closec(unit) int *unit;
#else
????
#endif

{ if (*unit>MAXUNIT || *unit<0)
        { fprintf(stderr,"closec: unit out of range\n");
           MACRO_QUIT();
           return;
        }
  if (files[*unit].fd)
        { close(files[*unit].fd);
          files[*unit].fd=files[*unit].addr=files[*unit].size=0;
          *files[*unit].fname=NULL;
        }
  return;
}



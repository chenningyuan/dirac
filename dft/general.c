/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* general.c:
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02
   NOTES: Adding new functionals:
   a. use fun-slater.c as template.
   b. add 'extern Functional MyFunctional;' to functionals.h
   c. add '&MyFunctional' to available_functionals below.
   d. have a beer. Or some crackers, if you prefer.
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif

/* Use BSD's strncasecmp(); if there is a platform that has no strncasecmp()
 * ask pawsa@theochem.kth.se for replacement */
#define _BSD_SOURCE 1

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#define __CVERSION__

#include "general.h"
#include "functionals.h"
#include "inforb.h"
#include "dcbham.h"

#ifdef VAR_MPI
#include <mpi.h>
#include "infpar.h"
#endif

/* C-wide constants */
const int  ZEROI = 0,   ONEI = 1, THREEI = 3, FOURI = 4;
const real ZEROR = 0.0, ONER = 1.0, TWOR = 2.0, FOURR = 4.0;

/* =================================================================== */
/* dftinput:

   read DFT functional from given line. The calling convention assumes
   Sun-style passing parameters. ATLAS linear algebra package
   http://www.netlib.org/atlas/ or http://math-atlas.sourceforge.net/
   contains an elaborate discuttion of character type variable passing
   conventions, in particular little bit below
   http://math-atlas.sourceforge.net/errata.html#RH7.0
*/
static char* DftConfString = NULL;
static char* EmbXcConfString = NULL;
static char* EmbKinConfString = NULL;

extern void dftsetcam_(real *w, real *b);

static void
dft_set_hf_weight(real hfxfac)
{
  dcrham_.hfxfac = hfxfac;
}

static real
dft_get_hf_weight(void)
{
    return dcrham_.hfxfac;
}

static void
dft_set_cam_param(real w, real be) { dftsetcam_(&w, &be);}

static int
dft_func_input (const char* line, int* inperr, int len, char** ConfString,
                Functional* avail[], Functional** func);

int
FSYM(dftinput)(const char* line, int * inperr, int len)
{
    fun_printf        = fort_print;
    fun_set_hf_weight = dft_set_hf_weight;
    fun_get_hf_weight = dft_get_hf_weight;
    fun_set_cam_param = dft_set_cam_param;
    return dft_func_input (line, inperr, len, &DftConfString,
                           available_functionals, &selected_func); 
}

int
FSYM(dft_emb_set_xcfunc)(const char* line, int * inperr, int len)
{
    return dft_func_input (line, inperr, len, &EmbXcConfString,
                           available_functionals, &emb_xc_func); 
}

int
FSYM(dft_emb_set_kinfunc)(const char* line, int * inperr, int len)
{
    return dft_func_input (line, inperr, len, &EmbKinConfString,
                           available_kin_functionals, &emb_kin_func); 
}

static int
dft_func_input (const char* line, int* inperr, int len, char** ConfString,
                Functional* avail[], Functional** func) 
{
    char func_name[20];
    int i, off;

    for(i=len-1; i>=0 && isspace((int)line[i]); i--)
        ;
    if(*ConfString) free(*ConfString);
    i++;
    for(off=0; line[off] && isspace((int)line[off]); off++)
        ;
    *ConfString = malloc(i+1-off);
    strncpy(*ConfString, line+off, i-off); 
    (*ConfString)[i-off] = '\0';
    sscanf(*ConfString,"%20s", func_name);
    for(i=0; avail[i]; i++)
        if(strcasecmp(avail[i]->name, func_name)==0) {
            int ok;
            *func = avail[i];
            ok = (*func)->read ?
                (*func)->read((*ConfString)+strlen(func_name)) : 1;
            if(!ok) (*inperr)++;
            return ok;
        }

    fort_print("Unknown functional '%s'. Aborting.\n", func_name);
    funlistfuncs(avail);
    (*inperr)++;

    return 0; /* failed to find */
}
static int do_alda  = 0; /* do not do ALDA approximation by default */

void
doalda_(const int *alda)
{
    do_alda = *alda;
}


void
dftpot0_(FirstDrv *ds, const real* weight, const real* rho, const real* grad)
{
    FunFirstFuncDrv drvs;
    FunDensProp dp;
    /* HP does not grok C99's { *rhoa, *rhob, *grada, *gradb }; */
    dp.rhoa  = dp.rhob  = *rho  *0.5;
    dp.grada = dp.gradb = *grad *0.5;
    dp.gradab = dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;
    drv1_clear(&drvs);
    selected_func->first(&drvs, *weight, &dp);
    ds->fR = drvs.df1000;
    ds->fZ = drvs.df0010 + 0.5*drvs.df00001* *grad;
}

void
dftpot1_(SecondDrv *ds, const real* w, const real* rho, const real* grad,
         const int* triplet)
{
    FunSecondFuncDrv drvs;
    FunDensProp dp;
    drv2_clear(&drvs);
    /* HP does not grok C99's { *rho*0.5, *rho*0.5, *grad*0.5, *grad*0.5 }; */
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5; 
    dp.gradab = dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;
    if(do_alda)
        LDAFunctional.second(&drvs, *w, &dp);
    else
        selected_func->second(&drvs, *w, &dp);
         
    ds->fR  = 0; /* NOT USED */
    ds->fZ  = drvs.df0010 + 0.5*drvs.df00001* *grad;
    ds->fRR = 0.5*(drvs.df2000 + drvs.df1100);
    ds->fRZ = 0.5*(drvs.df1010 + drvs.df1001 + drvs.df10001* *grad);
    ds->fZZ = 0.5*(drvs.df0020 + drvs.df0011 + drvs.df00001) + drvs.df00101* * grad + 0.25*drvs.df00002* *grad* *grad;
    ds->fRG = 0; /* NOT USED */
    ds->fZG = 0; /* NOT USED */
    ds->fGG = 0; /* NOT USED */
    ds->fG  = 0; /* NOT USED */
}

void
dftpot2_(ThirdDrv *ds, const real* w, const real* rho, const real* grad,
         const int* triplet)
{
    FunThirdFuncDrv drvs;
    FunDensProp dp;
    drv3_clear(&drvs);
    /* HP does not grok C99's { *rho*0.5, *rho*0.5, *grad*0.5, *grad*0.5 }; */
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5; 
    dp.gradab = dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;
    if(do_alda)
        LDAFunctional.third(&drvs, *w, &dp);
    else
        selected_func->third(&drvs, *w, &dp);
    
    ds->fR  = 0; /* NOT USED */
    ds->fZ  = drvs.df0010 + 0.5*drvs.df00001* *grad;
    ds->fRR = 0.5*(drvs.df2000 + drvs.df1100);
//    ds->fRZ = 0.5*(drvs.df1010 + drvs.df1001 + drvs.df10001* *grad);
//    ds->fZZ = 0.5*(drvs.df0020 + drvs.df0011 + drvs.df00001) + drvs.df00101* * grad + 0.25*drvs.df00002* *grad* *grad;
//    ds->fRG = 0; /* NOT USED */
//    ds->fZG = 0; /* NOT USED */
//    ds->fGG = 0; /* NOT USED */
//    ds->fG  = 0; /* NOT USED */
    ds->fRRR = 0.25*(drvs.df3000 + 3*drvs.df2100);
//    ds->fRRZ = 0.25*(drvs.df2010 + drvs.df2001 + 2*drvs.df1110 + 
//	            (drvs.df20001 + drvs.df11001)* *grad);
//    ds->fRZZ = 0.25*(drvs.df1020 + *drvs.df1002) + 0.5*drvs.df1011 +
//	       0.5*(drvs.df10101 + drvs.df10011)* *grad +
//	       0.25*drvs.df10002* *grad* *grad + 0.5*drvs.df10001;
//    ds->fZZZ = 0.25*drvs.df0030 + 0.75*drvs.df0021 + 0.75*drvs.df00201* *grad +
//	       0.75*drvs.df00111* *grad + 0.75*drvs.df00102* *grad* *grad +
//	       0.125*drvs.df00003* *grad* *grad* *grad + 1.5*drvs.df00101 + 0.75*drvs.df00002* *grad;
}

void
sdftpot1_(SecondDrv *ds, const real* w, const real* rho, const real* grad)
{
  FunSecondFuncDrv drvs;
  FunDensProp dp;

  drv2_clear(&drvs);

  dp.rhoa   = dp.rhob  = *rho *0.5;
  dp.grada  = dp.gradb = *grad*0.5;
  dp.gradab = dp.grada * dp.gradb;

  dp.subsystem[0] = NULL;
  dp.subsystem[1] = NULL;

  if(do_alda)
    LDAFunctional.second(&drvs, *w, &dp);
  else
    selected_func->second(&drvs, *w, &dp);
       
  ds->fR  = 0; /* NOT USED */
  ds->fZ  = 0; /* NOT USED */
  ds->fRR = 0.5*(drvs.df2000 - drvs.df1100);
  ds->fRZ = 0; /* NOT USED */
  ds->fZZ = 0.5*(drvs.df0020 - drvs.df0011 - drvs.df00001);
  ds->fRG = 0; /* NOT USED */
  ds->fZG = 0; /* NOT USED */
  ds->fGG = 0; /* NOT USED */
  ds->fG  = 0; /* NOT USED */
}

extern void quit_(const char* str, int len);
void
dalton_quit(const char* format, ...)
{
    char line[128];
    int len;
    va_list a;
 
    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = strlen(line);
    quit_(line, len);
}

/* Helper functions. Could be bracketed with #ifdef DEBUG or something */
extern void fortwrt_(const char *line, const int *linelen);
int
fort_print(const char* format, ...)
{
    char line[128];
    int len;
    va_list a;

/* DEC does not recognize the vsnprintf command */
/* unless when linked with -ldb */
#ifndef SYS_DEC
    va_start(a, format);
    vsnprintf(line, sizeof(line), format, a);
    va_end(a);
    len = strlen(line);
    fortwrt_(line, &len);
    return len;
#endif
}

/* ################################################################### */
/* MPI related functions                                               */
/* ################################################################### */
/* dft_sync_func:
   synchronizes functional data between nodes.
   Needs to be called only once (after functional initialization/change).
   NO-OP for serial code.
   In any case, do it only once, otherwise the functional becomes
   a sum of functionals.
*/
void
dftsyncfunc_(int* amImaster)
{
#ifdef VAR_MPI
    static int Initialized = 0;
    int success; /* jvs: this is not used !!! */
    int dft_conf_len =  DftConfString ? strlen(DftConfString)+1 : 0;
    
    if(infpar_.nodes==0||Initialized) return;

    MPI_Bcast(&do_alda,      1, MPI_INT, infpar_.mparid, MPI_COMM_WORLD);
    MPI_Bcast(&dft_conf_len, 1, MPI_INT, infpar_.mparid, MPI_COMM_WORLD);
    if(dft_conf_len>0) {
        char* line = malloc(dft_conf_len);
        if(*amImaster)
            strcpy(line, DftConfString);

        MPI_Bcast(line, dft_conf_len, MPI_CHAR, infpar_.mparid, MPI_COMM_WORLD);
        if(!*amImaster)
            dftinput_(line, &success, dft_conf_len);

        free(line);
        Initialized = 1;
    }
#endif /* VAR_MPI */
}

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

static int
dft_func_input (const char* line, int* inperr, int len, char** ConfString,
                Functional* avail[], Functional** func);

int
FSYM(dftinput)(const char* line, int * inperr, int len)
{
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
    dftlistfuncs_(avail);
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
    /* we should really use a separate structure here
     * but this time, we cheat */
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
         
    ds->fR  = 0.5*(drvs.df1000 + drvs.df0100);
    ds->fZ  = drvs.df0010;
    ds->fRR = 0.5*(drvs.df2000 + drvs.df1100);
    ds->fRZ = 0.5*(drvs.df1010 + drvs.df1001);
    ds->fZZ = 0.5*(drvs.df0020 + drvs.df0011);
    ds->fRG = 0.5*drvs.df10001;
    ds->fZG = 0.5*drvs.df00101;
    ds->fGG = 0.25*drvs.df00002;
    ds->fG  = 0.5*drvs.df00001;
}

#if defined(VAR_MPI)
#include <mpi.h>
#define MASTER_NO 0
/* =================================================================== */
/* General parallel routines.
 * dft_lifesupport_sync - one-time sync of basic common block data
 *                   that is crucial for initialization.
 * dft_wake_slaves - wakes slaves for a specified DFT evaluator (perhaps this
 *                   can be generalized?). It gets a registered property
 *                   evaluator address as an argument and broadcasts a
 *                   message that will call corresponding slave code
 *                   via dft_cslave().
 * FIXME: use rather a register approach that does not require modification
 * of this file every time new property evaluator is added.
 */

struct {
    DFTPropEvalMaster master_func;
    DFTPropEvalSlave  slave_func;
} PropEvaluatorList[] = {
#if 0
    { (DFTPropEvalMaster)dft_kohn_sham_,   dft_kohn_sham_slave   },
    { (DFTPropEvalMaster)dft_kohn_shamab_, dft_kohn_shamab_slave },
    { (DFTPropEvalMaster)dft_lin_resp_,    dft_lin_resp_slave    },
#endif
    { (DFTPropEvalMaster)FSYM2(dft_lin_respf),   dft_lin_respf_slave   },
#if 0
    { (DFTPropEvalMaster)dft_lin_respab_,  dft_lin_respab_slave  },
    { (DFTPropEvalMaster)dft_mol_grad_,    dft_mol_grad_slave    },
#endif
    { (DFTPropEvalMaster)dftqrcf_,         dft_qr_resp_slave     }
};

/* mpi_sync_data:
   sync given set of data with slaves.
*/
void
mpi_sync_data(const SyncData* data, int count)
{
    int i;
    for(i=0; i<count; i++) {
        MPI_Bcast(data[i].data, data[i].count, data[i].type,
                  0, MPI_COMM_WORLD);
    }
}


void
dft_wake_slaves(DFTPropEvalMaster evaluator)
{
    static int iprtyp = 5; /* magic DFT/C number */
    static int iprint = 0;
    int id, mynum;

    MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
    if(mynum != 0)
        return; /* slaves do not wake up other slaves */

    for(id=0; 
        id<ELEMENTS(PropEvaluatorList) && 
            PropEvaluatorList[id].master_func != evaluator;
        id++)
        ;
        
    if(id>=ELEMENTS(PropEvaluatorList)) {
        /* this would really be an programming error.... */
        fprintf(stderr, "Evaluator not registered. No slaves activated.\n");
        return;
    }
    /* ignore MPI errors */
    MPI_Bcast(&iprtyp,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&iprint,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&id,    1, MPI_INT, 0, MPI_COMM_WORLD);
    FSYM(dftintbcast)();
}

/* dft_cslave:
   a slave task handler. Receives the task ID and calls apropriate
   property evaluator.
   We broadcast also some basic properties and common block data
   that slaves should absolutely know but for some reason do not.
*/
void
FSYM2(dft_cslave)(real* work,int*lwork,int*iprint)
{
    int rank, size;
    if(MPI_Comm_rank(MPI_COMM_WORLD, &rank) ||
       MPI_Comm_size(MPI_COMM_WORLD, &size)) printf("MPI error\n");
    else {
        int id;
        MPI_Bcast(&id,1,MPI_INT, MASTER_NO, MPI_COMM_WORLD);
        FSYM(dftintbcast)();
        (PropEvaluatorList[id].slave_func)(work, lwork, iprint);
    }
}
#else
void
FSYM2(dft_cslave)(real* work,int*lwork,int*iprint)
{
   fort_print("DFT slave called but does nothing now.");
}
#endif /* VAR_MPI */


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

/* ===================================================================
 * Parallel support.
 * =================================================================== */
#ifdef VAR_MPI
#include <mpi.h>

/* dftfuncsync synchronises the selected functional between all
 * nodes. */
void
FSYM(dftfuncsync)(int *mynum, int *nodes)
{
    static int done = 0;
    int len = DftConfString ? strlen(DftConfString) + 1: 0;
    if(done) return;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(len>0) {
        int res;
        char *line = malloc(len);
        if(*mynum == 0)
            strcpy(line, DftConfString);
        MPI_Bcast(line, len, MPI_CHAR, 0, MPI_COMM_WORLD);
        if(*mynum != 0) {
            int res;
            FSYM(dftsetfunc)(line, &res, len);
        }
        free(line);
    }
    done = 1;
}

#endif
 

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
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* functionals.c:
   Program for computing and testing functional routines in the DFT module.
   (c) Pawel Salek, pawsa@theochem.kth.se, 2001-08-02

   The test build can be done by:
   g77 -O -DTEST_BUILD functionals.c -o functionals -L. -ldft -lm
   or 
   cc -O -DTEST_BUILD functionals.c -o functionals -L. -ldft -lm -lg2c

   NOTES: this file is short but in a separate file to reduce the
   number of dependences and be able to easily compile the code for
   the TEST_BUILD. 

   NOTES:
   the closed shell assumptions make following derivatives equivalent:

   3000 = 0300
   2100 = 1200
   2010 = 0201
   1110 = 1101
   0210 = 2001
   1020 = 0102
   0120 = 1002
   0030 = 0003
   1011 = 0111
   0021 = 0012
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define __CVERSION__

#include "functionals.h"
#include "dcbham.h"
#include "inforb.h"
#ifdef VAR_MPI
#include <mpi.h>
#include "infpar.h"
#endif

Functional* available_functionals[] = {
    /* generic functionals */
    &SlaterFunctional,
    &VWNFunctional,
    &B88CorrFunctional,
    &BeckeFunctional,
    &LYPFunctional,
    &PW86xFunctional,
    &P86cFunctional,
    /* mixed functionals */
    &XAlphaFunctional,
    &LDAFunctional,
    &GGAKeyFunctional,
    &BLYPFunctional,
    &B3LYPFunctional,
    &WAHFunctional,

    /* generic embedding functional */
    &EmbeddingFunctional,
    NULL
};    

/* kinetic energy functional */
Functional* available_kin_functionals[] = {
    &TF_KinFunctional,
    &VW_KinFunctional,
    NULL
};
Functional* selected_func = &LDAFunctional;

/* =================================================================== */

void
drv1_clear(FirstFuncDrv* gga)
{
    gga->df1000 = gga->df0100 = gga->df0010 = gga->df0001 = gga->df00001 = 0;
}
                                                                                void
drv2_clear(SecondFuncDrv* gga)
{
    gga->df1000 = gga->df0100 = gga->df0010 = gga->df0001 = gga->df00001 = 0;
    gga->df2000 = gga->df0200 = gga->df0020 = gga->df0002 = gga->df00002 = 0;
    gga->df1010 = gga->df0101 = gga->df1001 = gga->df0110 = 0;
    gga->df1100 = gga->df0011 = gga->df10001 = gga->df01001 = 0;
}
                                                                                
void
drv3_clear(ThirdFuncDrv* gga)
{
    gga->df1000 = gga->df0100 = gga->df0010 = gga->df0001 = gga->df00001 = 0;
    gga->df2000 = gga->df0200 = gga->df0020 = gga->df0002 = 0;
    gga->df1100 = gga->df0011 = 0;
    gga->df1010 = gga->df1001 = gga->df0101 = gga->df0110 = 0;
    gga->df10001 = gga->df01001 = 0;
                                                                                
    gga->df3000 = gga->df0300 = gga->df0030 = gga->df0003 = 0;
    gga->df2100 = gga->df1200 = gga->df0012 = gga->df0021 = 0;
    gga->df2010 = gga->df2001 = gga->df0201 = gga->df0210 = 0;
    gga->df1020 = gga->df1002 = gga->df0102 = gga->df0120 = 0;
    gga->df1110 = gga->df1101 = gga->df1011 = gga->df0111 = 0;
    gga->df20001 = gga->df02001 = gga->df11001 = 0;
}

/* fortran (and not only) functional stub routines */
int
dftisgga_(void)
{
    return selected_func->is_gga();
}

int
dftisemb_(void)
{
    return !strcmp(selected_func->name, "embedding");
}

void
dft_set_hf_weight(real hfxfac)
{
#ifndef NO_BACKWARD_COMP
  dcrham_.hfxfac = hfxfac;
#endif
}

real
dft_get_hf_weight(void)
{
#ifndef NO_BACKWARD_COMP
    return dcrham_.hfxfac;
#else
    return 0;
#endif
}

real
dftenergy_(const real* rho, const real* grad)
{
    DftDensProp dp;
    /* HP does not grok { *rho*0.5, *rho*0.5, *grad*0.5, *grad*0.5 }; */
    dp.rhoa  = dp.rhob  = *rho *0.5;
    dp.grada = dp.gradb = *grad*0.5; 
    dp.gradab= dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;
    return selected_func->func(&dp);
}

void
dft_check_consistency()
{
    static const char* OPEN_ERR =
	"\nERROR: currently, DFT does not handle open shells "
	"and active orbitals.";
    if(inforb_.nasht>0) {
	fort_print(OPEN_ERR); fprintf(stderr, "%s\n\n", OPEN_ERR);
	exit(1);
    }
}

/* dftreport:
   report the data and check consistency of the input (FIXME: this should
   be done by a separate routine.  
*/
void
dftreport_(void)
{
    fort_print("\n     This is a DFT calculation of type: %s",
               selected_func->name);
    if(selected_func->report)
        selected_func->report();

    dft_check_consistency();
}

/* Helper functions. Could be bracketed with #ifdef DEBUG or something */
void
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
#ifndef NO_BACKWARD_COMP
    fortwrt_(line, &len);
#endif
#endif
}

void
dftlistfuncs_(void)
{
    int i;
    fort_print("\nAvailable functionals:");
    for(i=0; available_functionals[i]; i++)
        fort_print(available_functionals[i]->name);
}


/* =================================================================== */
/* dftinput:

   read DFT input from given line. The calling convention assumes
   Sun-style passing parameters. ATLAS linear algebra package
   http://www.netlib.org/atlas/ or http://math-atlas.sourceforge.net/
   contains an elaborate discuttion of character type variable passing
   conventions,  in particular little bit below
   http://math-atlas.sourceforge.net/errata.html#RH7.0
*/
static char* DftConfString = NULL;
static int do_alda = 0;
int
dftinput_(char* line, int* success, int len)
{
    char func_name[20];
    int i;
    for(i=len-1; i>=0 && isspace((int)line[i]); i--)
        ;
    line[i+1] = '\0';
    if(DftConfString) free(DftConfString);
    DftConfString = strdup(line);
    /* fort_print("Reading DFT input from '%s'...", line); */
    sscanf(line,"%20s", func_name);

    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            selected_func = available_functionals[i];
            *success = selected_func->read ?
                selected_func->read(line+strlen(func_name)) : 1;
            return;
        }
    fort_print("Unknown functional '%s'. Aborting.\n", func_name);
    *success = 0; /* failed to find */
}

void
doalda_(const int *alda)
{
    do_alda = *alda;
}

void
dftpot0_(FirstDrv *ds, const real* weight, const real* rho, const real* grad)
{
    FirstFuncDrv drvs;
    DftDensProp dp;
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
    SecondFuncDrv drvs;
    DftDensProp dp;
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
    ds->fZZ = 0.5*(drvs.df0020 + drvs.df0011 + drvs.df00001);
}

/* dftpot2_:
   computes third order derivatives of selected functional with respect
   to rho and zeta=|\nabla\rho|
*/
#ifdef TESTED_THE_CODE_BELOW
void
dftpot2_(ThirdDrv *ds, real factor, real rho, real grad, int isgga)
{
    ThirdFuncDrv drvs;
    DftDensProp dp;
    /* HP does not grok { *rho*0.5, *rho*0.5, *grad*0.5, *grad*0.5 }; */
    dp.rhoa  = dp.rhob  = rho *0.5;
    dp.grada = dp.gradb = grad*0.5; 
    dp.gradab = dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;

    drv3_clear(&drvs);
    selected_func->third(&drvs, factor, &dp);
    /* transform to density from density_alpha derivatives below */
    /* This could be a separate function. */
    /* we treat singlet here */
    ds->fR   = (drvs.df1000);
    ds->fRR  = (drvs.df2000 + drvs.df1100);
    ds->fRRR = (drvs.df3000 + 3*drvs.df2100);

    if(isgga) { /* FORMULAE seem ok but were never really tested */
        real grada = 0.5*grad;
        real grada2= grada*grada;
        real grada3= grada2*grada;
	/* transform to closed shell, second time. Oh, I love this mess. */
        ds->fZ  = drvs.df0010/(2*grada)
            + drvs.df00001;
        ds->fZZ  = (drvs.df0020 + drvs.df0011)/(4*grada2) 
	    -drvs.df0010/(4*grada3);
        ds->fRZ  = (drvs.df1010 + drvs.df1001)/(2*grada);
        ds->fRG  = 2*drvs.df10001;
        ds->fRRZ = (drvs.df2010+drvs.df2001+2*drvs.df1110)/(2*grada);
        ds->fRRG = drvs.df20001+drvs.df11001;
        ds->fRRGX = 2*(drvs.df20001+drvs.df11001);
         
        ds->fRZZ = (drvs.df1020+drvs.df0120+2*drvs.df1011)/(4*grada2) 
                      - (drvs.df1010+drvs.df1001)/(4*grada3);
        ds->fZZZ = ((drvs.df0030 + 3*drvs.df0021)/grada3 
                   -3*(drvs.df0020 + drvs.df0011)/(grada2*grada2)
                   +3*drvs.df0010/(grada3*grada2))/8.0; 
    } else {
        ds->fZ = ds->fZZ = ds->fRZ = 0;
        ds->fRRZ  = 0;
        ds->fZZZ  = 0; 
        ds->fG    = ds->fRG = ds->fRRG = 0;
        ds->fRRGX = 0;  
    }
}
#endif
/* MIXED FUNCTIONALS */
typedef struct FuncList_ FuncList;

struct FuncList_ {
    Functional* func;
    real        weight;
    FuncList*   next;
};

FuncList* gga_fun_list = NULL;

static int
gga_isgga(void)
{ 
    int res = 0;

    FuncList* lst;
    for(lst=gga_fun_list; lst && !res; lst=lst->next) 
        res |= lst->func->is_gga();
    return res;
}

static FuncList*
add_functional(FuncList* lst, Functional* f, float weight)
{
    FuncList* n = malloc(sizeof(FuncList));
    n->func = f; n->weight = weight; n->next = lst;
    return n;
}
    
int
fun_false() { return 0; }
int
fun_true()  { return 1; }


static int
xalpha_read(const char* conf_line)
{
    float weight;
    int res = (sscanf(conf_line, "%g", &weight)==1);
    if(res) 
        gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 
                                      1.5*weight);
    dft_set_hf_weight(0);
    return res;
}

static int
lda_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,   1.0);
    dft_set_hf_weight(0);
    return 1;
}

static int
blyp_read(const char* conf_line)
{
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, 1.0);
    gga_fun_list = add_functional(gga_fun_list, &B88CorrFunctional,  1.0);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    1.0);
    dft_set_hf_weight(0);
    return 1;
}

static int
b3lyp_read(const char* conf_line)
{
    static const real lypw = 0.81, dirw = 0.8;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, dirw);
    gga_fun_list = add_functional(gga_fun_list, &B88CorrFunctional,0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lypw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1-lypw);
    dft_set_hf_weight(1-dirw);
    return 1;
}

static int
wah_read(const char* conf_line)
{
    static const real sw = 0.95, lw = 0.81;
    gga_fun_list = add_functional(gga_fun_list, &SlaterFunctional, sw);
    gga_fun_list = add_functional(gga_fun_list, &B88CorrFunctional,0.72);
    gga_fun_list = add_functional(gga_fun_list, &LYPFunctional,    lw);
    gga_fun_list = add_functional(gga_fun_list, &VWNFunctional,    1-lw);
    dft_set_hf_weight(1 - sw);
    return 1;
}

/* gga_key_read:
   read general GGA input: 
*/
static int
gga_key_read(const char* conf_line)
{
  int res = 1, i;
  float f;
  const char* str = conf_line;

  dft_set_hf_weight(0);

  while(*str) {
    while(*str && isspace((int)*str)) str++; /* skip whitespace */
    if(*str =='\0') break; /* line ended by whitespace */
    if(strncasecmp("HF=", str, 3)==0) {
        if(sscanf(str+3,"%g", &f) != 1) {
            fort_print("GGAKey: HF not followed by the weight: ", conf_line);
            res = 0;
        } else     dft_set_hf_weight(f);
    } else {
        for(i=0; available_functionals[i]; i++) {
            int len = strlen(available_functionals[i]->name);
            if(strncasecmp(available_functionals[i]->name, str, len)==0 &&
               str[len] == '=') {
                if(sscanf(str+len+1,"%g", &f) != 1) {
                    fort_print("GGAKey: keyword '%s' not followed by "
                               "weight: %s", available_functionals[i]->name, 
                               conf_line);
                    res = 0;
                } else {
                    gga_fun_list = add_functional(gga_fun_list, 
                                                  available_functionals[i], f);
                    break; /* weight properly read, break the 'for' loop */
                }
            }
        }  
        if(available_functionals[i] == NULL) {
            fort_print("GGAKey: functional '%s' not recognised: ", str);
            res = 0;
        }
    }
    while(*str && !isspace((int)*str)) str++; /* skip nonws */
  }
    
  return res;
}


static void
gga_report(void)
{
    FuncList* lst;
    fort_print("Weighted mixed functional:");
    if(dft_get_hf_weight()>0)
        fort_print("%30s: %10.5f", "HF exchange", dft_get_hf_weight());
    for(lst=gga_fun_list; lst; lst=lst->next) 
        fort_print("%30s: %10.5f", lst->func->name, lst->weight);
}

static real
gga_energy(DftDensProp* dp)
{
    real res = 0;
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        res += lst->weight*lst->func->func(dp);
    return res;
}

static void
gga_first(FirstFuncDrv *ds, real factor, DftDensProp* dp) 
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->first(ds, factor*lst->weight, dp);
}

static void
gga_second(SecondFuncDrv *ds, real factor, DftDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->second(ds, factor*lst->weight, dp);
}

static void
gga_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp)
{
    FuncList* lst;
    for(lst=gga_fun_list; lst; lst=lst->next) 
        lst->func->third(ds, factor*lst->weight, dp);
}

static real
lda_energy(DftDensProp* dp)
{
    return SlaterFunctional.func(dp) + VWNFunctional.func(dp);
}

static void
lda_first(FirstFuncDrv *ds, real factor, DftDensProp* dp) 
{
    SlaterFunctional.first(ds, factor, dp);
    VWNFunctional.first   (ds, factor, dp);
}

static void
lda_second(SecondFuncDrv *ds, real factor, DftDensProp* dp)
{
    SlaterFunctional.second(ds, factor, dp);
    VWNFunctional.second   (ds, factor, dp);
}

static void
lda_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp)
{
    SlaterFunctional.third(ds, factor, dp);
    VWNFunctional.third   (ds, factor, dp);
}

Functional XAlphaFunctional = {
  "XAlpha",       /* name */
  fun_false,   /* gga-corrected */
  xalpha_read, 
  NULL,
  gga_energy, 
  gga_first,
  gga_second,
  gga_third
};

Functional LDAFunctional = {
  "LDA",         /* name */
  fun_false,     /* not gga-corrected */
  lda_read,
  NULL,
  lda_energy,
  lda_first,
  lda_second,
  lda_third
};

Functional GGAKeyFunctional = {
  "GGAKey",      /* name */
  gga_isgga,     /* gga-corrected */
  gga_key_read,
  gga_report,
  gga_energy, 
  gga_first,
  gga_second,
  gga_third
};

Functional BLYPFunctional = {
  "BLYP",      /* name */
  gga_isgga,     /* gga-corrected */
  blyp_read,
  gga_report,
  gga_energy, 
  gga_first,
  gga_second,
  gga_third
};



Functional B3LYPFunctional = {
  "B3LYP",      /* name */
  gga_isgga,     /* gga-corrected */
  b3lyp_read,
  gga_report,
  gga_energy, 
  gga_first,
  gga_second,
  gga_third
};

Functional WAHFunctional = {
  "WAH",      /* name */
  gga_isgga,     /* gga-corrected */
  wah_read,
  gga_report,
  gga_energy, 
  gga_first,
  gga_second,
  gga_third
};

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

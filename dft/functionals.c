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

 */

#define _BSD_SOURCE 1

#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define __CVERSION__

#include "functionals.h"

Functional* available_functionals[] = {
    /* generic functionals */
    &BeckeFunctional,
    &KTFunctional,
    &LB94Functional,
    &LYPFunctional,
    &OPTXFunctional,
    &P86cFunctional,
    &PW86xFunctional,
    &Pw91cFunctional,
    &PZ81Functional,
    &PbecFunctional,
    &PbexFunctional,
    &SlaterFunctional,
    &VWN3Functional,
    &VWN5Functional,
    &VWNIFunctional,
    &VWNFunctional,
    &XAlphaFunctional,
    /* mixed functionals */
    &B3LYPFunctional,
    &B3LYPGaussFunctional,
    &B3P86Functional,
    &B3P86GFunctional,
    &BLYPFunctional,
    &BP86Functional,
    &BPW91Functional,
    &Camb3lypFunctional,
    &GGAKeyFunctional,
    &KT1Functional,
    &KT2Functional,
    &KT3Functional,
    &LDAFunctional,
    &OLYPFunctional,
    &PBE0Functional,
    &PBEFunctional,
    &SVWN3Functional,
    &SVWN5Functional,
    NULL
};

/* kinetic energy functionals */
Functional* available_kin_functionals[] = {
    &TF_KinFunctional,
    &VW_KinFunctional,
    NULL
};

static int my_printf(const char *fmt, ...)
{
    int i;va_list ap; va_start(ap, fmt); i= vprintf(fmt, ap); va_end(ap);
    puts("");
    return i;
}
 
static void set_hf_weight(real w)         {}
static real get_hf_weight(void)           {return 0;}
static void set_cam_param(real w, real b) {}

Functional* selected_func = &LDAFunctional;
Functional* emb_xc_func   = &LDAFunctional;
Functional* emb_kin_func  = NULL;
int (*fun_printf)(const char *fmt, ...) = my_printf;
void (*fun_set_hf_weight)(real w)         = set_hf_weight;
real (*fun_get_hf_weight)(void)           = get_hf_weight;
void (*fun_set_cam_param)(real w, real b) = set_cam_param;

/* =================================================================== */
void
drv1_clear(FunFirstFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

void
drv2_clear(FunSecondFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

void
drv3_clear(FunThirdFuncDrv* gga)
{
    memset(gga, 0, sizeof(*gga));
}

enum FunError
fun_select_by_name(const char *conf_string)
{
    int ok, i;
    char func_name[20];

    sscanf(conf_string,"%20s", func_name);
    for(i=0; available_functionals[i]; i++)
        if(strcasecmp(available_functionals[i]->name, func_name)==0) {
            selected_func = available_functionals[i];
            ok = selected_func->read ?
                selected_func->read(conf_string+strlen(func_name)) : 1;
            return ok ? FUN_OK : FUN_CONF_ERROR;
        }
    return FUN_UNKNOWN;
}

int fun_true(void)  { return 1; }
int fun_false(void) { return 0; }

/* Fortran interface. We specify different names suffixes so that
 * library can be linked with code compiled with different compilers
 * or different compilation options. */

void
funset(const char *str, int *info, int len)
{
    switch(fun_select_by_name(str)) {
    case FUN_OK:         *info = 0; break;
    case FUN_UNKNOWN   : *info = 1; break;
    case FUN_CONF_ERROR: *info = 2; break;
    }
}
void
funset_(const char *str, int *info, int len)
{
    switch(fun_select_by_name(str)) {
    case FUN_OK:         *info = 0; break;
    case FUN_UNKNOWN   : *info = 1; break;
    case FUN_CONF_ERROR: *info = 2; break;
    }
}

int
funisgga(void)
{ return selected_func->is_gga(); }
int
funisgga_(void)
{ return selected_func->is_gga(); }

/* end of backward compatiblity section */

int
dftisemb_(void)
{
    return !strcmp(selected_func->name, "embedding");
    /* or return selected_func == EmbFunctional */
}

real
funenergy(const FunDensProp *dp)
{ return selected_func->func(dp); }
real
funenergy_(const FunDensProp *dp)
{ return selected_func->func(dp); }

void
funfirst(FunFirstFuncDrv *df, real *factor, const FunDensProp *dp)
{ selected_func->first(df, *factor, dp); }
void
funfirst_(FunFirstFuncDrv *df, real *factor, const FunDensProp *dp)
{ selected_func->first(df, *factor, dp); }

void
funlistfuncs(Functional *l[])
{
    int i;
    
    fun_printf("\nAvailable functionals:");
    for(i=0; l[i]; i++)
        fun_printf(l[i]->name);
}

/* =================================================================== */
/*           fortran (and not only) functional stub routines           */
/* =================================================================== */

/* dftreport:
   report the selected functional and its configuration.
*/
void
dftreport_(void)
{
    fun_printf("\n     This is a DFT calculation of type: %s",
               selected_func->name);
    if(selected_func->report)
        selected_func->report();
}

void
dftlistfuncs_(void)
{
    funlistfuncs(available_functionals);
}

/* declare both known fortran name-mangled variants */
int
dftisgga_(void)
{ return selected_func->is_gga(); }

int
dftisgga__(void)
{ return selected_func->is_gga(); }

real
dftenergy_(const real* rho, const real* grad)
{
    FunDensProp dp;
    /* HP does not grok { *rho*0.5, *rho*0.5, *grad*0.5, *grad*0.5 }; */
    if (*rho > 1e-20)	
        dp.rhoa  = dp.rhob  = *rho *0.5;
    else dp.rhoa = dp.rhob = 1e-20;
    if (*grad > 1e-20) 
        dp.grada = dp.gradb = *grad*0.5;
    else dp.grada = dp.gradb = 1e-20; 
    dp.gradab= dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;
    return selected_func->func(&dp);
}

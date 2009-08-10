#include <iostream>
#include <cstdlib>
#include "taylor.h"
#include "xc_fortran_interface.h"
#include "auto_functionals.h"
#include "dispatch.h"

using namespace std;

/*
  Use the variable ordering R,S,Z,X,Y; which are
  defined as

      R = \rho_\alpha + \rho_\beta
      S = \rho_\alpha - \rho_\beta

      Z = \nabla R \cdot \nabla R
      X = \nabla S \cdot \nabla S
      Y = \nabla R \cdot \nabla S

      and conversely

      na = (R + S)/2
      nb = (R - S)/2

      Dna^2 = (X+Z+2*Y)/4
      Dnb^2 = (X+Z-2*Y)/4
      DnaDnb = (Z - X)/4

   It it assumed below that all derivatives are taken at S = 0,
   and then also X = Y = 0.

 */

// If you change MAX_ORDER you also have to update TEMPLATE_TAB!
#define MAX_ORDER 8

template<int Ndeg>
static void lda_helper(int functional,
		       int order,
		       double R,
		       double *res)
{
  taylor<2,Ndeg> *t = reinterpret_cast< taylor<2,Ndeg> * >(res);
  switch(functional)
    {
    case LDA_VWN:
      *t = vwn_correlation<Ndeg>(R);
      break;
    case LDA_SLATER:
      *t = slater_exchange<Ndeg>(R);
      break;
    default:
      assert(0 && "Unknown LDA functional");
      break;
    }
}

extern "C" void
lda_eval_(const int *functional,
	  const int *order,
	  const double* R,
	  double* result)
{
  static void (*lda_helper_tab[])(int, int, double, double *) =
    {FUNLIST(lda_helper,MAX_ORDER)};
  if (*order >= 0 && *order < MAX_ORDER)
    {
      lda_helper_tab[*order](*functional,*order,*R,result);
    }
  else
    {
      cerr << "lda_eval(): Invalid functional order " << *order << 
	", MAX_ORDER = " << MAX_ORDER << endl;
      exit(-1);
    }
}

template<int Nvar, int Ndeg>
static void contract_helper(double *dxc_dD,
			    const double* dxc_drho,
			    const double* drho_dD)
{
  taylor<1, Ndeg> *pout = reinterpret_cast< taylor<1,Ndeg>* >(dxc_dD);
  const taylor<Nvar,Ndeg> *pxc = 
    reinterpret_cast< const taylor<Nvar,Ndeg>* >(dxc_drho);
  taylor<1, Ndeg> at[Nvar];
  // Insert the relations dExc/dD = dExc/drho * drho/dD, by
  // inserting drho_i = drho_dD[i]*dD
  for (int i=0;i<Nvar;i++)
    {
      at[0] = 0;
      at[0][1] = drho_dD[i];
    }
  *pout = pxc->eval(at); //Result is a Ndeg taylor expansion in D
}

/*
  Given drho_dD (an array of length 2), and the xc
  derivatives up to some order (from lda_eval),
  calculate dxc_dD up to order.
 */
extern "C" void
lda_contract_(double *dxc_dD,
	      const int *order,
	      const double* dxc_drho,
	      const double* drho_dD)
{
  switch(*order)
    {
    case 0:
      contract_helper<2,0>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 1:
      contract_helper<2,1>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 2:
      contract_helper<2,2>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 3:
      contract_helper<2,3>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 4:
      contract_helper<2,4>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 5:
      contract_helper<2,5>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 6:
      contract_helper<2,6>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 7:
      contract_helper<2,7>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 8:
      contract_helper<2,8>(dxc_dD,dxc_drho,drho_dD);
      break;
    default:
      assert(0 && "Unimplemented order in lda_contract");
    };
}


/*
  Given drho_dD (an array of length 5), and the xc
  derivatives up to some order (from gga_eval),
  calculate dxc_dD up to order.
 */
extern "C" void
gga_contract_(double *dxc_dD,
	      const int *order,
	      const double* dxc_drho,
	      const double* drho_dD)
{
  switch(*order)
    {
    case 0:
      contract_helper<5,0>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 1:
      contract_helper<5,1>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 2:
      contract_helper<5,2>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 3:
      contract_helper<5,3>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 4:
      contract_helper<5,4>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 5:
      contract_helper<5,5>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 6:
      contract_helper<5,6>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 7:
      contract_helper<5,7>(dxc_dD,dxc_drho,drho_dD);
      break;
    case 8:
      contract_helper<5,8>(dxc_dD,dxc_drho,drho_dD);
      break;
    default:
      assert(0 && "Unimplemented order in gga_contract");
    };
}

template<int Ndeg>
static void gga_helper(int functional,
		       int order,
		       double R, double Z,
		       double *res)
{
#ifndef TAYLOR_QD
  taylor<5,Ndeg> *t = reinterpret_cast< taylor<5,Ndeg> * >(res);
#else
  taylor<5,Ndeg> tmp;
  taylor<5,Ndeg> *t = &tmp;
  unsigned int cw;
  fpu_fix_start(&cw);
#endif
  switch(functional)
    {
    case GGA_LYP:
      *t = lyp_correlation<Ndeg>(R,Z);
      break;
    case GGA_BECKE:
      *t = becke_exchange<Ndeg>(R,Z);
      break;
    case GGA_PBEX:
      *t = pbe_exchange<Ndeg>(R,Z);
      break;
    case GGA_PBEC:
      *t = pbe_correlation<Ndeg>(R,Z);
      break;
    default:
      assert(0 && "Unknown GGA functional");
      break;
    }
#ifdef TAYLOR_QD
  for (int i=0;i<tmp.size;i++)
    res[i] = to_double(tmp[i]);
  fpu_fix_end(&cw);
#endif
}

extern "C" void
gga_eval_(const int *functional,
	  const int *order,
	  const double *R, const double *Z,
	  double *result)
{
  static void (*gga_helper_tab[])(int, int, double, double,
					      double *) =
    {FUNLIST(gga_helper,MAX_ORDER)};
  if (*order >= 0 && *order < MAX_ORDER)
    {
      gga_helper_tab[*order](*functional,*order,*R,*Z,result);
    }
  else
    {
      cerr << "gga_eval(): Invalid functional order " << *order << 
	", MAX_ORDER = " << MAX_ORDER << endl;
      exit(-1);
    }
}

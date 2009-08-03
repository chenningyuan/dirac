#include "taylor.h"
#include "xc_fortran_interface.h"
#include "auto_functionals.h"

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

template<int Ndeg>
static void lda_helper(int functional,
		       int order,
		       double R,
		       double *res)
{
#ifndef TAYLOR_QD
  taylor<2,Ndeg> *t = reinterpret_cast< taylor<2,Ndeg> * >(res);
#else
  taylor<2,Ndeg> tmp;
  taylor<2,Ndeg> *t = &tmp;
  unsigned int cw;
  fpu_fix_start(&cw);
#endif
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
#ifdef TAYLOR_QD
  for (int i=0;i<tmp.size;i++)
    res[i] = to_double(tmp[i]);
  fpu_fix_end(&cw);
#endif
}

extern "C" void
lda_eval_(const int *functional,
	  const int *order,
	  const double* R,
	  double* result)
{
  switch(*order)
    {
    case 0:
      lda_helper<0>(*functional,*order,*R,result);
      break;
    case 1:
      lda_helper<1>(*functional,*order,*R,result);
      break;
    case 2:
      lda_helper<2>(*functional,*order,*R,result);
      break;
    case 3:
      lda_helper<3>(*functional,*order,*R,result);
      break;
    case 4:
      lda_helper<4>(*functional,*order,*R,result);
      break;
    case 5:
      lda_helper<5>(*functional,*order,*R,result);
      break;
    case 6:
      lda_helper<6>(*functional,*order,*R,result);
      break;
    case 7:
      lda_helper<7>(*functional,*order,*R,result);
      break;
    case 8:
      lda_helper<8>(*functional,*order,*R,result);
      break;
    default:
      assert(0 && "Unimplemented order of LDA functional");
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
  switch(*order)
    {
    case 0:
      gga_helper<0>(*functional,*order,*R, *Z,result);
      break;
    case 1:
      gga_helper<1>(*functional,*order,*R, *Z,result);
      break;
    case 2:
      gga_helper<2>(*functional,*order,*R, *Z,result);
      break;
    case 3:
      gga_helper<3>(*functional,*order,*R, *Z,result);
      break;
    case 4:
      gga_helper<4>(*functional,*order,*R, *Z,result);
      break;
    case 5:
      gga_helper<5>(*functional,*order,*R, *Z,result);
      break;
    default:
      assert(0 && "Unimplemented order of GGA functional");
    };
}

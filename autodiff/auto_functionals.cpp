#include "taylor.h"
#include "auto_functionals.h"

// VWN:

const tnum_t  para[] = {-0.10498,   0.0621814,  3.72744, 12.9352};
const tnum_t ferro[] = {-0.325,     0.0310907,  7.06042, 18.0578};
const tnum_t inter[] = {-0.0047584, -0.0337737, 1.13107, 13.0045};

static tnum_t a(const tnum_t p[])
{
  return p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1;
}

static tnum_t b(const tnum_t p[])
{
  return 2*(p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1) + 2;
}

static tnum_t c(const tnum_t p[])
{
  return 2*p[2]*(1/sqrt(4*p[3] - p[2]*p[2]) - p[0]/
		 ((p[0]*p[0] + p[0]*p[2] + p[3])*
		  sqrt(4*p[3] - p[2]*p[2])/(p[2] + 2*p[0])));
}

template<int N>
static taylor<2,N> x(const taylor<2,N>& s,const tnum_t p[])
{
  return s*s + p[2]*s + p[3];
}

template<int N>
static taylor<2,N> y(const taylor<2,N>& s,const tnum_t p[])
{
  return s - p[0];
}

template<int N>
static taylor<2,N> z(const taylor<2,N>& s,const tnum_t p[])
{
  return sqrt(4*p[3] - p[2]*p[2])/(2*s + p[2]);
}

template<int N>
static taylor<2,N> f(const taylor<2,N>& s,const tnum_t p[])
{
  return 0.5*p[1]*(2*log(s) + a(p)*log(x(s, p)) - b(p)*log(y(s, p)) + c(p)*atan(z(s, p)));
}

// Return a N:th order Taylor expansion of the VWN functional around
// (R0,S0)
template<int N>
static taylor<2,N> vwn_correlation(tnum_t R0)
{
  assert(R0 > 1e-20);
  taylor<2,N> R(R0,0,AUTO_FUNCTIONAL_SCALING); // R = R0 + x
  taylor<2,N> S(0,1,AUTO_FUNCTIONAL_SCALING); // S = S0 + y
  taylor<2,N> s = pow(3/(4*TAYLOR_PI),tnum_t(1.0)/6)/pow(R,tnum_t(1.0)/6);
  taylor<2,N> t = S/R;
  taylor<2,N> g = 1.92366105093154*(pow(1 + t,tnum_t(4.0)/3) + pow(1 - t,tnum_t(4.0)/3) - 2);
  taylor<2,N> d = g*((f(s, ferro) - f(s, para))*pow(t,4) +
    f(s, inter)*(1 - pow(t,4))*0.584822305543806);
  return R*(f(s, para) + d);
}


/*
  Because of a design "oddity" in C++ we have to force
  instantiation here to all orders we need of the functionals.
  - ulf
 */
template
taylor<2,0> vwn_correlation(tnum_t R0);
template
taylor<2,1> vwn_correlation(tnum_t R0);
template
taylor<2,2> vwn_correlation(tnum_t R0);
template
taylor<2,3> vwn_correlation(tnum_t R0);
template
taylor<2,4> vwn_correlation(tnum_t R0);
template
taylor<2,5> vwn_correlation(tnum_t R0);
template
taylor<2,6> vwn_correlation(tnum_t R0);
template
taylor<2,7> vwn_correlation(tnum_t R0);
template
taylor<2,8> vwn_correlation(tnum_t R0);



template<int N>
static taylor<2,N> slater_exchange(tnum_t R0)
{
  taylor<2,N> R(R0,0,AUTO_FUNCTIONAL_SCALING), S(0,1,AUTO_FUNCTIONAL_SCALING);
  const tnum_t prefactor = -3.0/4.0*pow(6/TAYLOR_PI, tnum_t(1.0)/3.0);
  return prefactor*(pow((R+S)/2,tnum_t(4.0)/3) + pow((R-S)/2,tnum_t(4.0)/3));
}

template
taylor<2,0> slater_exchange(tnum_t R0);
template
taylor<2,1> slater_exchange(tnum_t R0);
template
taylor<2,2> slater_exchange(tnum_t R0);
template
taylor<2,3> slater_exchange(tnum_t R0);
template
taylor<2,4> slater_exchange(tnum_t R0);
template
taylor<2,5> slater_exchange(tnum_t R0);
template
taylor<2,6> slater_exchange(tnum_t R0);
template
taylor<2,7> slater_exchange(tnum_t R0);
template
taylor<2,8> slater_exchange(tnum_t R0);



/*  Implementation of LYP functional and its derivatives
   (c) Pawel Salek, pawsa@theochem.kth.se, oct 2001
    Z. Rinkevicius modification for open-shell, general 5 variables formalism.
    autodiff version by Ulf Ekstrom 2009 */
template<int N>
static taylor<5,N> lyp_correlation(tnum_t R0, tnum_t Z0)
{
    const tnum_t A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    const tnum_t CF = 0.3*pow(3*TAYLOR_PI*TAYLOR_PI,2.0/3.0);

    taylor<5,N> rho(R0,0,AUTO_FUNCTIONAL_SCALING);
    taylor<5,N> rho2 = rho*rho;
    taylor<5,N> rhom13 = pow(rho,-1.0/3.0);
    taylor<5,N> S(0,1,AUTO_FUNCTIONAL_SCALING), 
      Z(Z0,2,AUTO_FUNCTIONAL_SCALING), 
      X(0,3,AUTO_FUNCTIONAL_SCALING), 
      Y(0,4,AUTO_FUNCTIONAL_SCALING);
    taylor<5,N> ngrada2 = (X+Z+2*Y)/4;
    taylor<5,N> ngradb2 = (X+Z-2*Y)/4;
    taylor<5,N> rhoa = (rho + S)/2;
    taylor<5,N> rhob = (rho - S)/2;

    taylor<5,N> denom = 1+D*rhom13;
    taylor<5,N> omega = exp(-C*rhom13)/denom*pow(rho,-11.0/3.0);
    taylor<5,N> delta = rhom13*(C + D/denom);

    taylor<5,N> t1 =  pow(2.0,11.0/3.0)*CF*(pow(rhoa,8.0/3.0) +
					    pow(rhob,8.0/3.0));
    taylor<5,N> t2 =  ((47.0 - 7.0*delta)/18.0)*Z;
    taylor<5,N> t3 = -(2.5 -delta/18.0)*(ngrada2+ngradb2);
    taylor<5,N> t4 =  (11.0-delta)/9.0*(rhoa*ngrada2 + rhob*ngradb2)/rho;
    taylor<5,N> t5 = -2.0/3.0*rho2*Z;
    taylor<5,N> t6 = ((2.0/3.0*rho2-rhoa*rhoa)*ngradb2 +
          (2.0/3.0*rho2-rhob*rhob)*ngrada2);
    return -A*(4*rhoa*rhob/(denom*rho)
	      +B*omega*(rhoa*rhob*(t1+t2+t3+t4)+t5+t6));
}


template
taylor<5,0> lyp_correlation(tnum_t R0, tnum_t Z0);
template
taylor<5,1> lyp_correlation(tnum_t R0, tnum_t Z0);
template
taylor<5,2> lyp_correlation(tnum_t R0, tnum_t Z0);
template
taylor<5,3> lyp_correlation(tnum_t R0, tnum_t Z0);
template
taylor<5,4> lyp_correlation(tnum_t R0, tnum_t Z0);
template
taylor<5,5> lyp_correlation(tnum_t R0, tnum_t Z0);


static const tnum_t BECKE_THRESHOLD = 1e-14;
static const tnum_t BETA = 0.0042;
template<int N>
static taylor<5,N> becke_exchange(tnum_t R0,tnum_t Z0)
{
  
  taylor<5,N> ea(0),eb(0);
  taylor<5,N> R(R0,0,AUTO_FUNCTIONAL_SCALING), S(0,1,AUTO_FUNCTIONAL_SCALING);
  taylor<5,N> Z(Z0,2,AUTO_FUNCTIONAL_SCALING), 
    X(0,3,AUTO_FUNCTIONAL_SCALING), Y(0,4,AUTO_FUNCTIONAL_SCALING);
  taylor<5,N> rhoa = (R + S)/2;
  taylor<5,N> rhob = (R - S)/2;

  if (rhob[0]>BECKE_THRESHOLD)
    {
      taylor<5,N> gb = sqrt(X+Z-2*Y)/2;  
      taylor<5,N> xb = gb*pow(rhob,-4.0/3.0);
      taylor<5,N> rb = pow(rhob,4.0/3.0);
      taylor<5,N> denomb = 1.0 +6.0*xb*BETA*asinh(xb);
       eb = rb*xb*xb/denomb; 
    } 
  if (rhoa[0]>BECKE_THRESHOLD) 
    {
      taylor<5,N> ga = sqrt(X+Z+2*Y)/2;  
      taylor<5,N> xa = ga*pow(rhoa,-4.0/3.0);
      taylor<5,N> ra = pow(rhoa,4.0/3.0);
      taylor<5,N> denoma = 1.0 +6.0*BETA*xa*asinh(xa);
      ea = ra*xa*xa/denoma;
    }
  return -BETA*(ea+eb);
}


template
taylor<5,0> becke_exchange(tnum_t R0, tnum_t Z0);
template
taylor<5,1> becke_exchange(tnum_t R0, tnum_t Z0);
template
taylor<5,2> becke_exchange(tnum_t R0, tnum_t Z0);
template
taylor<5,3> becke_exchange(tnum_t R0, tnum_t Z0);
template
taylor<5,4> becke_exchange(tnum_t R0, tnum_t Z0);
template
taylor<5,5> becke_exchange(tnum_t R0, tnum_t Z0);

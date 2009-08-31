#include "taylor.h"
#include "auto_functionals.h"


// VWN:

const double  para[] = {-0.10498,   0.0621814,  3.72744, 12.9352};
const double ferro[] = {-0.325,     0.0310907,  7.06042, 18.0578};
const double inter[] = {-0.0047584, -0.0337737, 1.13107, 13.0045};

static double a(const double p[])
{
  return p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1;
}

static double b(const double p[])
{
  return 2*(p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1) + 2;
}

static double c(const double p[])
{
  return 2*p[2]*(1/sqrt(4*p[3] - p[2]*p[2]) - p[0]/
		 ((p[0]*p[0] + p[0]*p[2] + p[3])*
		  sqrt(4*p[3] - p[2]*p[2])/(p[2] + 2*p[0])));
}

template<int N>
static taylor<2,N> x(const taylor<2,N>& s,const double p[])
{
  return s*s + p[2]*s + p[3];
}

template<int N>
static taylor<2,N> y(const taylor<2,N>& s,const double p[])
{
  return s - p[0];
}

template<int N>
static taylor<2,N> z(const taylor<2,N>& s,const double p[])
{
  return sqrt(4*p[3] - p[2]*p[2])/(2*s + p[2]);
}

template<int N>
static taylor<2,N> f(const taylor<2,N>& s,const double p[])
{
  return 0.5*p[1]*(2*log(s) + a(p)*log(x(s, p)) - b(p)*log(y(s, p)) + 
		   c(p)*atan(z(s, p)));
}

// Return a N:th order Taylor expansion of the VWN functional around
// (R0,S0)
template<int N>
static taylor<2,N> vwn_correlation(const double d[2])
{
  taylor<2,N> R,S;
  get_R(R,d);
  assert(R[0] > 1e-20);
  get_S(S,d);
  taylor<2,N> t,g,dd,s = pow(3/(4*TAYLOR_PI),double(1.0)/6)/pow(R,double(1.0)/6);
  t = S/R;
  g = 1.92366105093154*(pow(1 + t,double(4.0)/3) + pow(1 - t,double(4.0)/3) - 2);
  dd = g*((f(s, ferro) - f(s, para))*pow(t,4) +
	 f(s, inter)*(1 - pow(t,4))*0.584822305543806);
  return R*(f(s, para) + dd);
}


/*
  Because of a design "oddity" in C++ we have to force
  instantiation here to all orders we need of the functionals.
  - ulf
 */
template
taylor<2,0> vwn_correlation(const double d[2]);
template
taylor<2,1> vwn_correlation(const double d[2]);
template
taylor<2,2> vwn_correlation(const double d[2]);
template
taylor<2,3> vwn_correlation(const double d[2]);
template
taylor<2,4> vwn_correlation(const double d[2]);
template
taylor<2,5> vwn_correlation(const double d[2]);
template
taylor<2,6> vwn_correlation(const double d[2]);
template
taylor<2,7> vwn_correlation(const double d[2]);
template
taylor<2,8> vwn_correlation(const double d[2]);



template<int N>
static taylor<2,N> slater_exchange(const double d[2])
{
  
  taylor<2,N> R,S;
  get_R(R,d);
  get_S(S,d);
  const double prefactor = -3.0/4.0*pow(6/TAYLOR_PI, double(1.0)/3.0);
  return prefactor*(pow((R+S)/2,double(4.0)/3) + pow((R-S)/2,double(4.0)/3));
}

template
taylor<2,0> slater_exchange(const double d[2]);
template
taylor<2,1> slater_exchange(const double d[2]);
template
taylor<2,2> slater_exchange(const double d[2]);
template
taylor<2,3> slater_exchange(const double d[2]);
template
taylor<2,4> slater_exchange(const double d[2]);
template
taylor<2,5> slater_exchange(const double d[2]);
template
taylor<2,6> slater_exchange(const double d[2]);
template
taylor<2,7> slater_exchange(const double d[2]);
template
taylor<2,8> slater_exchange(const double d[2]);



/*  Implementation of LYP functional and its derivatives
   (c) Pawel Salek, pawsa@theochem.kth.se, oct 2001
    Z. Rinkevicius modification for open-shell, general 5 variables formalism.
    autodiff version by Ulf Ekstrom 2009 */
template<int N>
static taylor<5,N> lyp_correlation(const double d[5])
{
    const double A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    const double CF = 0.3*pow(3*TAYLOR_PI*TAYLOR_PI,2.0/3.0);

    taylor<5,N> rho;
    get_R(rho,d);
    taylor<5,N> rho2 = rho*rho;
    taylor<5,N> rhom13 = pow(rho,-1.0/3.0);
    taylor<5,N> Z, ngrada2, ngradb2, rhoa, rhob;
    get_Z(Z,d);
    get_rhoa(rhoa,d);
    get_rhob(rhob,d);
    get_grada2(ngrada2,d);
    get_gradb2(ngradb2,d);

    taylor<5,N> denom = 1+D*rhom13;
    taylor<5,N> omega = exp(-C*rhom13)/denom*pow(rho,-11.0/3.0);
    taylor<5,N> delta = rhom13*(C + D/denom);

    taylor<5,N> t1 =  pow(2.0,11.0/3.0)*CF*(pow(rhoa,8.0/3.0) +
					    pow(rhob,8.0/3.0));
    t1 +=  ((47.0 - 7.0*delta)/18.0)*Z;
    t1 += -(2.5 -delta/18.0)*(ngrada2+ngradb2);
    t1 +=  (11.0-delta)/9.0*(rhoa*ngrada2 + rhob*ngradb2)/rho;
    taylor<5,N> t5 = -2.0/3.0*rho2*Z;
    t5 += ((2.0/3.0*rho2-rhoa*rhoa)*ngradb2 +
	   (2.0/3.0*rho2-rhob*rhob)*ngrada2);
    return -A*(4*rhoa*rhob/(denom*rho)
	      +B*omega*(rhoa*rhob*t1+t5));
}


template
taylor<5,0> lyp_correlation(const double d[5]);
template
taylor<5,1> lyp_correlation(const double d[5]);
template
taylor<5,2> lyp_correlation(const double d[5]);
template
taylor<5,3> lyp_correlation(const double d[5]);
template
taylor<5,4> lyp_correlation(const double d[5]);
template
taylor<5,5> lyp_correlation(const double d[5]);


static const double BECKE_THRESHOLD = 1e-14;
static const double BETA = 0.0042;
template<int N>
static taylor<5,N> becke_exchange(const double d[5])
{
  
  taylor<5,N> ea = 0,eb = 0;
  taylor<5,N> rhoa,rhob;
  get_rhoa(rhoa,d);
  get_rhob(rhob,d);

  if (rhob[0]>BECKE_THRESHOLD)
    {
      taylor<5,N> gb2;
      get_gradb2(gb2,d);
      taylor<5,N> xb = sqrt(gb2)*pow(rhob,-4.0/3.0);
      taylor<5,N> rb = pow(rhob,4.0/3.0);
      taylor<5,N> denomb = 1.0 +6.0*xb*BETA*asinh(xb);
       eb = rb*xb*xb/denomb; 
    } 
  if (rhoa[0]>BECKE_THRESHOLD) 
    {
      taylor<5,N> ga2;
      get_grada2(ga2,d);
      taylor<5,N> xa = sqrt(ga2)*pow(rhoa,-4.0/3.0);
      taylor<5,N> ra = pow(rhoa,4.0/3.0);
      taylor<5,N> denoma = 1.0 +6.0*BETA*xa*asinh(xa);
      ea = ra*xa*xa/denoma;
    }
  return -BETA*(ea+eb);
}


template
taylor<5,0> becke_exchange(const double d[5]);
template
taylor<5,1> becke_exchange(const double d[5]);
template
taylor<5,2> becke_exchange(const double d[5]);
template
taylor<5,3> becke_exchange(const double d[5]);
template
taylor<5,4> becke_exchange(const double d[5]);
template
taylor<5,5> becke_exchange(const double d[5]);

// This file is ment to be included from xc_fun.cpp 

/*  Implementation of LYP functional and its derivatives
   (c) Pawel Salek, pawsa@theochem.kth.se, oct 2001
    Z. Rinkevicius modification for open-shell, general 5 variables formalism.
    autodiff version by Ulf Ekstrom 2009 */
template<class num>
static num lyp_correlation(const densvars<num> &d)
{
    const double A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    using xc_constants::CF;

    num rho2 = d.n*d.n;
    num rhom13 = pow(d.n,-1.0/3.0);
    num denom = 1+D*rhom13;
    num omega = exp(-C*rhom13)/denom*pow(d.n,-11.0/3.0);
    num delta = rhom13*(C + D/denom);

    num t1 = pow(2.0,11.0/3.0)*CF*(pow(d.a,8.0/3.0) +
				   pow(d.b,8.0/3.0));
    t1 +=  ((47.0 - 7.0*delta)/18.0)*d.gnn;
    t1 += -(2.5 -delta/18.0)*(d.gaa+d.gbb);
    t1 +=  (11.0-delta)/9.0*(d.a*d.gaa + d.b*d.gbb)/d.n;
    num t5 = -2.0/3.0*rho2*d.gnn + 
      ((2.0/3.0*rho2-d.a*d.a)*d.gbb +
       (2.0/3.0*rho2-d.b*d.b)*d.gaa);
    return -A*(4*d.a*d.b/(denom*d.n)
	      +B*omega*(d.a*d.b*t1+t5));
}


template<class num>
static num becke_exchange(const densvars<num> &d)
{
  const double BECKE_THRESHOLD = 1e-14;
  const double BETA = 0.0042;  
  num ea = 0, eb = 0;

  if (d.b>BECKE_THRESHOLD)
    {
      num xb = sqrt(d.gbb)*pow(d.b,-4.0/3.0);
      num rb = pow(d.b,4.0/3.0);
      num denomb = 1.0 +6.0*xb*BETA*asinh(xb);
      eb = rb*xb*xb/denomb; 
    } 
  if (d.a>BECKE_THRESHOLD) 
    {
      num xa = sqrt(d.gaa)*pow(d.a,-4.0/3.0);
      num ra = pow(d.a,4.0/3.0);
      num denoma = 1.0 +6.0*BETA*xa*asinh(xa);
      ea = ra*xa*xa/denoma;
    }
  return -BETA*(ea+eb);
}

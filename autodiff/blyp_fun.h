// This file is ment to be included from xc_fun.cpp 

/*  Implementation of LYP functional and its derivatives
   (c) Pawel Salek, pawsa@theochem.kth.se, oct 2001
    Z. Rinkevicius modification for open-shell, general 5 variables formalism.
    autodiff version by Ulf Ekstrom 2009 */
template<int N>
static taylor<5,N> lyp_correlation(const taylor<5,N> &R,
				   const taylor<5,N> &S,
				   const taylor<5,N> &Z,
				   const taylor<5,N> &X,
				   const taylor<5,N> &Y)
{
    const double A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    const double CF = 0.3*pow(3*TAYLOR_PI*TAYLOR_PI,2.0/3.0);

    taylor<5,N> rho2 = R*R;
    taylor<5,N> rhom13 = pow(R,-1.0/3.0);
    taylor<5,N> ngrada2 = get_grada2(Z,X,Y), 
      ngradb2 = get_gradb2(Z,X,Y), 
      rhoa = get_rhoa(R,S), 
      rhob = get_rhob(R,S);

    taylor<5,N> denom = 1+D*rhom13;
    taylor<5,N> omega = exp(-C*rhom13)/denom*pow(R,-11.0/3.0);
    taylor<5,N> delta = rhom13*(C + D/denom);

    taylor<5,N> t1 = pow(2.0,11.0/3.0)*CF*(pow(rhoa,8.0/3.0) +
					   pow(rhob,8.0/3.0));
    t1 +=  ((47.0 - 7.0*delta)/18.0)*Z;
    t1 += -(2.5 -delta/18.0)*(ngrada2+ngradb2);
    t1 +=  (11.0-delta)/9.0*(rhoa*ngrada2 + rhob*ngradb2)/R;
    taylor<5,N> t5 = -2.0/3.0*rho2*Z;
    t5 += ((2.0/3.0*rho2-rhoa*rhoa)*ngradb2 +
	   (2.0/3.0*rho2-rhob*rhob)*ngrada2);
    return -A*(4*rhoa*rhob/(denom*R)
	      +B*omega*(rhoa*rhob*t1+t5));
}


template<int N>
static taylor<5,N> becke_exchange(const taylor<5,N> &R,
				  const taylor<5,N> &S,
				  const taylor<5,N> &Z,
				  const taylor<5,N> &X,
				  const taylor<5,N> &Y)
{
  static const double BECKE_THRESHOLD = 1e-14;
  static const double BETA = 0.0042;  
  taylor<5,N> ea = 0,eb = 0;
  taylor<5,N> rhoa = get_rhoa(R,S),rhob = get_rhob(R,S);

  if (rhob[0]>BECKE_THRESHOLD)
    {
      taylor<5,N> gb2 = get_gradb2(Z,X,Y);
      taylor<5,N> xb = sqrt(gb2)*pow(rhob,-4.0/3.0);
      taylor<5,N> rb = pow(rhob,4.0/3.0);
      taylor<5,N> denomb = 1.0 +6.0*xb*BETA*asinh(xb);
       eb = rb*xb*xb/denomb; 
    } 
  if (rhoa[0]>BECKE_THRESHOLD) 
    {
      taylor<5,N> ga2 = get_grada2(Z,X,Y);;
      taylor<5,N> xa = sqrt(ga2)*pow(rhoa,-4.0/3.0);
      taylor<5,N> ra = pow(rhoa,4.0/3.0);
      taylor<5,N> denoma = 1.0 +6.0*BETA*xa*asinh(xa);
      ea = ra*xa*xa/denoma;
    }
  return -BETA*(ea+eb);
}

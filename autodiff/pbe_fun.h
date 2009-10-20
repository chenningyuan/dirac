/*
  PBE Exchange
  Imported from Dalton fun-pbex.c, originally by Pawel Salek(?).
 */

template<class num>
static num pbe_exchange(const num &rhoa,
			const num &rhob,
			const num &grada2,
			const num &gradb2)
{
  double t1 = pow(2.0,0.33333333333333);
  return  0.5*(-1.477117532764045*t1*pow(rhob,1.333333333333333)*
             (1.804-0.804/(0.00449279948023*gradb2/pow(rhob,2.666666666666667)+
			   1.0))
	       -1.477117532764045*t1*pow(rhoa,1.333333333333333)*
	     (1.804-0.804/(0.00449279948023*grada2/pow(rhoa,2.666666666666667)+
	     1.0)));
}


template<class num>
static num pbe_correlation_unpolarized(const num &R,
				       const num &Z)
{
  num t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

  t1 = R;
  t2 = 1/pow(t1,0.33333333333333);
  t3 = 0.1325688999052*t2+1.0;
  t4 = log(16.0818243221511/(2.225569421150687*t2+0.18970043257476/
      pow(t1,0.66666666666667)+5.98255043577108/pow(t1,0.16666666666667)+
        0.80042863499936/sqrt(t1))+1.0);
  t5 = Z;
  t6 = 1/pow(t1,2.333333333333334);
  t7 = pow(t5,2.0);
  t8 = 1/pow(t1,4.666666666666667);
  t9 = exp(1.999923215200099*t3*t4)-1.0;
  t10 = 1/t9;

  return t1*(0.03109219370394*
          log(2.146079601698033*(0.00864486614818*
				 t7*t8*t10+0.0634682060977*t5*t6)/
	      (0.01855257090002*t7*t8/pow(t9,2.0)
	       +0.13620782246265*t5*t6*t10+1.0)+1.0)-0.062182*t3*t4);
}

// PBEc implementation from MOLPRO manual, by ulfek.
namespace PBEc_internal
{

  const double iota = 0.0716;
  const double kappa = 0.004235;
  const double nu = 16*pow(3*M_PI*M_PI,1.0/3)/M_PI;
  const double lambda = nu*kappa;
  
  template<class num>
  static num zeta(const num &a, const num &b)
  {
    return (a-b)/(a+b);
  }

  template<class num>
  static num r(const num &a, const num &b)
  {
    return 0.25*pow(3,1.0/3)*pow(4,2.0/3)*pow(1/(M_PI*(a+b)),1.0/3);
  }
  
  template<class num>
  static num u(const num &zeta)
  {
    return 0.5*(pow(1+zeta,2.0/3)+(pow(1-zeta,2.0/3)));
  }

  template<class num>
  static num omega(const num &z)
  {
    return (pow(1+z,4.0/3)+pow(1-z,4.0/3)-2)/(2*pow(2,1.0/3)-2);
  }
  
  
  template<class num>
  static num M(const num &d, const num &a, const num &b,
	       const num &rho)
  {
    const double Z = -0.001667;
    return 0.5*nu*(phi(r(a,b))-kappa-3.0/7*Z)*d*d*
      exp(-335.9789467*pow(9/(pow(M_PI,5)*rho),1.0/3)*d*d);
  }
  
  template<class num>
  static num K(const num &d, const num &a, const num &b)
  {
    0.25*lambda*lambda/iota*
      log(1+2*iota*
	  (d*d+N(a,b)*pow(d,4))/
	  (lambda*(1+N(a,b)*d*d+pow(N(a,b),2)*pow(d,4))));
  }
  
  
  template<class num>
  static num e(const num &r, const double t[7])
  {
    return -2*t[0]*(1+t[1]*r)*log(1+0.5/
				  (t[0]*(t[2]*sqrt(r)+t[3]*r+t[4]*pow(r,3.0/2)+t[5]*pow(r,t[6]+1))));
  }

  template<class num>
  static num phi(const num &r)
  {
    const double 
      Xi = 23.266,
      Phi = 0.007389,
      Lambda = 8.723,
      Upsilon = 0.472;
    return 1.0/1000*(2.568+r*(Xi+Phi*r))/(1+r*(Lambda+r*(Upsilon+10*Phi*r)));
  }

  template<class num>
  static num epsilon(const num &a, const num &b)
  {
    // The T..P parameters, but transposed for easier use with the e function
    const double TUVWXYP[3][7] = {{0.031091,0.21370,7.5957 ,3.5876,1.6382 ,0.49294,1},
				  {0.01554 ,0.20548,14.1189,6.1977,3.3662 ,0.62517,1},
				  {0.016887,0.11125,10.357 ,3.6231,0.88026,0.49671,1}};
    const double c = 1.709921;
    num rval = r(a,b);
    num zetaval = zeta(a,b);
    num zeta4 = pow(zetaval,4);
    num omegaval = omega(zetaval);
    return e(rval,TUVWXYP[0])
      - e(rval,TUVWXYP[2])*omegaval*(1-zeta4)/c 
      + (e(rval,TUVWXYP[1]) - e(rval,TUVWXYP[0]))*omegaval*zeta4;
  }
  
  template<class num>
  static num A(const num &a, const num &b)
  {
    return 2*iota/(lambda*(exp(-2*iota*epsilon(a,b)/
			       (lambda*lambda*pow(u(zeta(a,b)),3))) - 1));
  }

  template<class num>
  static num H(const num &d, const num &a, const num &b)
  {
  num d2 = d*d;
  num d2A = d2*A(a,b);
  return 0.5*lambda*lambda/iota*pow(u(zeta(a,b)),3)*
    log(1+2*(iota*d2*(1 + d2A))/(lambda*(1+d2A*(1+d2A))));
  }

  template<class num>
  static num pbe_correlation(const num &R,
			     const num &S,
			     const num &Z,
			     const num &rhoa,
			     const num &rhob)
  {
    const num &sigma = Z;
    num d = 1.0/12*(sqrt(sigma)*pow(3,5.0/6))/
      (u(zeta(rhoa,rhob))*pow(M_PI,-1.0/6)*pow(R,7.0/6));
    return R*(epsilon(rhoa,rhob) + H(d,rhoa,rhob));
  }

}

using PBEc_internal::pbe_correlation;

template<class num>
static num rpbe_exchange(const num &rhoa,
			 const num &rhob,
			 const num &grada,
			 const num &gradab,
			 const num &gradb)
{
  double t1 = pow(2.0,.3333333333333333);

   /* code */
  return 0.5*(-1.477117532764045*(0.804*(1.0-1.0/exp(.004492799480231091*pow(gradb,2.0)/pow(rhob,2.666666666666667)))+
        1.0)*t1*pow(rhob,1.333333333333333)-1.477117532764045*(0.804*
        (1.0-1.0/exp(.004492799480231091*pow(grada,
        2.0)/pow(rhoa,2.666666666666667)))+1.0)*t1*pow(rhoa,1.333333333333333));
}

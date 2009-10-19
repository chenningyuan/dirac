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
static num pbe_correlation(const num &R,
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
#if 0
template<class num>
static zeta(const num &a, const num &b)
{
  return (a-b)/(a+b);
}

template<class num>
static r(const num &a, const num &b)
{
  return 0.25*pow(3,1.0/3)*pow(4,2.0/3)*pow(1/(M_PI*(a+b)),1.0/3);
}

template<class num>
static u(const num &a, const num &b)
{
  return 0.5*(pow(1+zeta(a,b),2.0/3)+(pow(1-zeta(a,b),2.0/3)))
}

template<class num>
static H(const num &d, const num &a, const num &b)
{
  return 0.5*(pow(1+zeta(a,b),2.0/3)+(pow(1-zeta(a,b),2.0/3)))
}



template<class num>
static A(const num &a, const num &b)
{
  return 2*iota/(lambda*(exp(-2*iota*epsilon(a,b)/
			     (lambda*lambda)) - 1));
}

template<class num>
static K(const num &d, const num &a, const num &b)
{
  0.25*lambda*lambda/iota*
    log(1+2*iota*
	(d*d+N(a,b)*pow(d,4))/
	(lambda*(1+N(a,b)*d*d+pow(N(a,b),2)*pow(d,4))));
}

template<class num>
static e(const num &r,
	 double t,
	 double u,
	 double v,
	 double w,
	 double x,
	 double y,
	 double p)
{
  return -2*t*(1+u*r)*log(1+0.5/
	    (t*(v*sqrt(r)+w*r+x*pow(r,3.0/2)+y*pow(r,p+1))));
}

template<class num>
static omega(const num &z)
{

template<class num>
static num pbe_correlation_new(const num &R,
			       const num &Z)
{
  
}
#endif
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

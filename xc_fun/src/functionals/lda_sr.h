


template<class num>
num esrc_ldaerfspin(const densvars<num> &d, double mu)
{
/*
C     Short-range spin-dependent LDA correlation functional from 
C       Paziani, Moroni, Gori-Giorgi and Bachelet, PRB 73, 155111 (2006)
C
C     subroutine adapted from Molpro
*/
  return d.n*(pw92_eps(d) - ecorrlr(d,mu));
}

template<class num>
num Qrpa(const num &x)
{
  const double Acoul = 2.0*(log(2.0)-1.0)/(M_PI*M_PI);
  const double a2    = 5.84605;
  const double c2    = 3.91744;
  const double d2    = 3.44851;
  const double b2 = d2 - 3.0/(2*M_PI*Acoul)*pow(4.0/(9.0*M_PI),1.0/3.0);
  return Acoul*log( (1+x*(a2+x*(b2+c2*x)))/(1+x*(a2+d2*x)) );
}

template<class num>
num dpol(const num &rs)
{
  const double cf  = pow(9.0*M_PI/4.0,1.0/3.0);
  const double p2p = 0.04;
  const double p3p = 0.4319;
  num rs2 = rs*rs;
  return pow(2.0,5.0/3.0)/5.0*pow(cf,2)/rs2*
    (1.0 + (p3p - 0.454555)*rs)/(1.0 + p3p*rs + p2p*rs2);
}

template<class num>
num ecorrlr(const densvars<num> &d, double mu)
{
  const double alpha=pow(4.0/9.0/M_PI,1.0/3.0);
  double cf=1/alpha; //TODO: The normal CF?
  num phi=(pow(1.0+d.zeta,2.0/3.0)+pow(1.0-d.zeta,2.0/3.0))/2.0;
  // cc parameters from the fit
  const double adib = 0.784949;
  const double q1a = -0.388;
  const double q2a = 0.676;
  const double q3a = 0.547;
  const double t1a = -4.95;
  const double t2a = 1.0;
  const double t3a = 0.31;
  
  num b0 = adib*d.r_s;
  num rs2 = d.r_s*d.r_s;
  num rs3 = rs2*d.r_s;

  num d2anti=(q1a*d.r_s+q2a*rs2)*exp(-q3a*d.r_s)/rs2;
  num d3anti=(t1a*d.r_s+t2a*rs2)*exp(-t3a*d.r_s)/rs3;

  const num &z = d.zeta;
  num z2 = d.zeta*d.zeta;
  num coe2 = -3.0/8.0/rs3*(1.0-z2)*(g0f(d.r_s)-0.5);

  num coe3 = -(1.0-z2)*g0f(d.r_s)/(sqrt(2.0*M_PI)*rs3);

  num coe4 = -9.0/64.0/rs3*(pow((1.0+z)/2.0,2)*
			    dpol(d.r_s*pow(2/(1.0+z),1.0/3.0))+
			    pow((1.0-z)/2.0,2)*
			    dpol(d.r_s*pow(2.0/(1.0-z),1.0/3.0))+
			    (1-z*z)*d2anti-
			    pow(cf,2)/10.0*(pow(1.0+z,8.0/3.0)
					    +pow(1-z,8.0/3.0))/rs2);

  num coe5 = -9.0/40.0/(sqrt(2.0*M_PI)*rs3)*
    (pow((1.0+z)/2.0,2)*
     dpol(d.r_s*pow(2.0/(1.0+z),1.0/3.0))+
     pow((1.0-z)/2.0,2)*
     dpol(d.r_s*pow(2.0/(1.0-z),1.0/3.0))+(1.0-z2)*
     d3anti);
    
  num ec = pw92_eps(d);
  //call ecPW(d.r_s,z,ec,ecd,ecz)
  
  num b06 = pow(b0,6);
  num b08 = pow(b0,8);
  num a1 = 4*b06*coe3 + b08*coe5;
  num a2 = 4*b06*coe2 + b08*coe4 + 6*pow(b0,4)*ec;
  num a3 = b08*coe3;
  num a4 = b06*(pow(b0,2)*coe2 + 4*ec);
    
  return (pow(phi,3)*Qrpa(mu*sqrt(d.r_s)/phi)+a1*pow(mu,3)+a2*pow(mu,4)+a3*pow(mu,5)+
          a4*pow(mu,6)+pow(b0*mu,8)*ec)/pow(1.0+pow(b0,2)*pow(mu,2),4);
}

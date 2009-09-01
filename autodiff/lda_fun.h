// This file is ment to be included from xc_fun.cpp, not to be used
// on its own.

// VWN:

static double vwn_a(const double p[])
{
  return p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1;
}

static double vwn_b(const double p[])
{
  return 2*(p[0]*p[2]/(p[0]*p[0] + p[0]*p[2] + p[3]) - 1) + 2;
}

static double vwn_c(const double p[])
{
  return 2*p[2]*(1/sqrt(4*p[3] - p[2]*p[2]) - p[0]/
		 ((p[0]*p[0] + p[0]*p[2] + p[3])*
		  sqrt(4*p[3] - p[2]*p[2])/(p[2] + 2*p[0])));
}

template<int N>
static taylor<2,N> vwn_x(const taylor<2,N>& s,const double p[])
{
  return s*s + p[2]*s + p[3];
}

template<int N>
static taylor<2,N> vwn_y(const taylor<2,N>& s,const double p[])
{
  return s - p[0];
}

template<int N>
static taylor<2,N> vwn_z(const taylor<2,N>& s,const double p[])
{
  return sqrt(4*p[3] - p[2]*p[2])/(2*s + p[2]);
}

template<int N>
static taylor<2,N> vwn_f(const taylor<2,N>& s,const double p[])
{
  return 0.5*p[1]*(2*log(s) + vwn_a(p)*log(vwn_x(s, p)) - vwn_b(p)*log(vwn_y(s, p)) + 
		   vwn_c(p)*atan(vwn_z(s, p)));
}

template<int N>
static taylor<2,N> vwn_correlation(const double d[2])
{
  static const double  para[] = {-0.10498,   0.0621814,  3.72744, 12.9352};
  static const double ferro[] = {-0.325,     0.0310907,  7.06042, 18.0578};
  static const double inter[] = {-0.0047584, -0.0337737, 1.13107, 13.0045};
  taylor<2,N> R,S;
  get_R(R,d);
  assert(R[0] > 1e-20);
  get_S(S,d);
  taylor<2,N> t,g,dd,s = pow(3/(4*TAYLOR_PI),double(1.0)/6)/pow(R,double(1.0)/6);
  t = S/R;
  g = 1.92366105093154*(pow(1 + t,double(4.0)/3) + pow(1 - t,double(4.0)/3) - 2);
  dd = g*((vwn_f(s, ferro) - vwn_f(s, para))*pow(t,4) +
	  vwn_f(s, inter)*(1 - pow(t,4))*0.584822305543806);
  return R*(vwn_f(s, para) + dd);
}

template<int N>
static taylor<2,N> slater_exchange(const double d[2])
{
  
  taylor<2,N> R,S;
  get_R(R,d);
  get_S(S,d);
  const double prefactor = -3.0/4.0*pow(6/TAYLOR_PI, double(1.0)/3.0);
  return prefactor*(pow((R+S)/2,double(4.0)/3) + pow((R-S)/2,double(4.0)/3));
}

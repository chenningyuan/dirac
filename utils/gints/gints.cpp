#include <iostream>
#include <cmath>
#include "taylor.h"
#include "matrix.h"
#include "boys.h"

/*
  Internally use hermite-gaussians D^n g_s(x).
  A shell is a set of hermite gaussians with a given
  exponent and location.

  Functions needed:
  Calculating Cartesian/Spherical gaussian shell coefficients
  x Overlap given two shells and coefficients on those shells
  x Kinetic energy - " " -
  x Shell products -> third shell coefficients
  Local potentials: Product shell -> Nuclear attraction, multipole
  x Shell/Shell electrostatic
 */


using namespace std;

#define MAX_L 3

class shell
{
public:
  double alpha;
  double R[3];
};

matrix<int,polymul_internal::binomial<3+4*MAX_L,4*MAX_L>::value,3> poly3_exponents;
int poly3_term_index[4*MAX_L+1][4*MAX_L+1][4*MAX_L+1]; // Index of the term with the given exponents


taylor<3,4*MAX_L> fac3(0);
taylor<6,2*MAX_L> fac6(0);

int fac(int n)
{
  assert(n>=0);
  assert(n<13);
  int f = 1;
  while(n>1)
    {
      f *= n;
      n--;
    }
  return f;
}

template<class numtype, int Nvar, int Ndeg>
void polydfac(polynomial<numtype, Nvar,Ndeg> &p)
{
  int exponents[Nvar];
  for (int i = 0;i<Nvar;i++)
    exponents[i] = 0;
  for (int i = 0;i<p.size();i++)
    {
      p.exponents(i,exponents);
      numtype f = 1;
      for (int j=0;j<Nvar;j++)
	f *= fac(exponents[j]);
      p[i] = f;
    }
}

void setup_tables(void)
{
  polynomial<int,3,4*MAX_L> p;
  int e[3];
  for (int i=0;i<p.size();i++)
    {
      p.exponents(i,e);
      for (int j=0;j<3;j++)
	poly3_exponents(i,j) = e[j];
      poly3_term_index[e[0]][e[1]][e[2]] = i;
    }
  for (int i=0;i<4*MAX_L+1;i++)
    for (int j=0;j<4*MAX_L+1;j++)
      for (int k=0;k<4*MAX_L+1;k++)
	{
	  int ijk[3] = {i,j,k};
	  poly3_term_index[i][j][k] = p.term_index(ijk);
	}						  
  polydfac(fac6);
  polydfac(fac3);
}

template<class num, int Nvar, int Ndeg>
void printpoly(ostream &dst, const polynomial<num,Nvar,Ndeg> &p)
{
  int exps[Nvar] = {0};
  for (int i=0;i<p.size();i++)
    {
      p.exponents(i,exps);
      for (int j=0;j<Nvar;j++)
	cout << exps[j] << " ";
      cout << "  " << p[i] << endl;
      assert(i == p.term_index(exps));
    }
}


// Multiply by  i!j!.. to each term, to get derivatives from
// Taylor coefficients
template<int N>
void dfactors(taylor<3,N> &t)
{
  assert(N <= 4*MAX_L); 
  for (int i=0;i<t.size();i++)
    t[i] *= fac3[i];
}



// Calculate the taylor expansion of (x + d)^2
template<int N>
void d2_taylor(taylor<3,N> &d2, const double d[3])
{
  d2[0] = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
  if (N>0)
    {
      d2[1] = 2*d[0];
      d2[2] = 2*d[1];
      d2[3] = 2*d[2];
    }
  if (N>1)
    {
      d2[4] = 1; //x^2
      d2[5] = 0;
      d2[6] = 0;
      d2[7] = 1; //y^2
      d2[8] = 0;
      d2[9] = 1; //z^2
    }
  for (int i=10;i<d2.size();i++)
    d2[i] = 0;
}


// Evaluate functions of this shell at R. This is the definition
// of the basis functions.
template<int N>
void funvalue(double f[],
	      const taylor<3,N> c[], int clen,
	      const shell &s, const double r[3])
{
  // f_ijk(r) = d^i/dR1 d^j/dR2 d^k/dR3 exp(-a(r-R)^2)
  double d[3];
  for (int i=0;i<3;i++)
    d[i] = r[i] - s.R[i];
  taylor<3,N> d2;
  d2_taylor(d2,d);
  taylor<3,N> fn = exp((-s.alpha)*d2);
  dfactors(fn);
  for (int i=0;i<clen;i++)
    f[i] = c[i].dot(fn);
}

/* 
   Represent the product D^n_R1 D^m_R2 g(a1,x-R1)*g(a2,x-R2) as a sum
   sum_k c_k D^k_RP g(a1+a2,x-RP)
   with RP(0) located at (a1*R1 + a2*R2)/(a1+a2).
   R1 is variable 0, R2 variable 1 and RP variable 2.
   Use a generating function for the expansion coefficients.

   All terms R1^n R2^m P^k with k > n + m are zero., so a bit
   tarded.

*/
template<int N>
void gauss_product_1d(taylor<3,N> &c, 
		      double R1, double R2,
		      double alpha1, double alpha2)
{
  double asum = alpha1 + alpha2;
  taylor<3,N> arg(0),r1mr2(R1-R2);
  if (N>0)
    {
      r1mr2[1] = 1;
      r1mr2[2] = -1;
      arg[1] = alpha1/asum;
      arg[2] = alpha2/asum;
    }
  arg *= taylor<3,N>(0,2);

  c = exp(-alpha1*alpha2/asum*pow(r1mr2,2))*exp(2*asum*arg);
  // Normalization
  double normfac[N+1];
  normfac[0] = 1;
  for (int i=1;i<=N;i++)
    normfac[i] = normfac[i-1]/(2*asum*i);
  for (int i=0;i<c.size();i++)
    c[i] *= normfac[poly3_exponents(i,2)];
  dfactors(c);
}


/* Given vectors on two shells s1 and s2, compute the new shell
   s12 of the product function and the expansion of the product
   at that shell. One expansion is computed for each pairs of vectors
   from c1 and c2. */
template<int N1, int N2, int Nc1, int Nc2>
void shell_products(matrix<taylor<3,N1+N2>,Nc1,Nc2> &prod, shell &s12,
		    const taylor<3,N1> c1[], const shell &s1,
		    const taylor<3,N2> c2[], const shell &s2)
{
  s12.alpha = s1.alpha + s2.alpha;
  for (int i=0;i<3;i++)
    s12.R[i] = (s1.alpha*s1.R[i] + s2.alpha*s2.R[i])/s12.alpha;
  taylor<3,2*(N1+N2)> c[3];
  gauss_product_1d(c[0],s1.R[0],s2.R[0],s1.alpha,s2.alpha);
  gauss_product_1d(c[1],s1.R[1],s2.R[1],s1.alpha,s2.alpha);
  gauss_product_1d(c[2],s1.R[2],s2.R[2],s1.alpha,s2.alpha);
  for (int i = 0; i < Nc1; i++)
    for (int j = 0; j < Nc2; j++)
      {
	prod(i,j).zero();
	int exponents[3];
	for (int k = 0; k < taylor<3,N1+N2>::size(); k++)
	  for (int m = 0;m < taylor<3,N1>::size(); m++)
	    for (int n = 0;n < taylor<3,N2>::size(); n++)
	      {
		double w = 1;
		for (int p=0;p<3;p++)
		  {		
		    exponents[0] = poly3_exponents(m,p);
		    exponents[1] = poly3_exponents(n,p);
		    exponents[2] = poly3_exponents(k,p);
		    //tabulated version of
		    //w *= c[p][c[p].term_index(exponents)];		    
		    w *= c[p][poly3_term_index[exponents[0]]
		      [exponents[1]][exponents[2]]];
		  }
		prod(i,j)[k] += w*c1[i][m]*c2[j][n];
	      }
      }
}


// Given a function that depends only on the distance
// between two shells, contract over vectors on these
// shells. f is the taylor expansion of the function,
// in the vector d[3], the distance between the shells.
// Application to overlap, kinetic energy and Coulomb
template<int N1, int N2, int Nc1, int Nc2>
void contract_d12(matrix<double,Nc1,Nc2> &res, 
		  const taylor<3,N1> c1[], 
		  const taylor<3,N2> c2[], 
		  const taylor<3,N1+N2> &f)
{
  taylor<3,N1> tmp1;
  taylor<3,N2> tmp2;
  for (int i=0;i<Nc2;i++)
    {
      tmp2 = c2[i];
      tmp2.invert_parity();
      polycontract(tmp1,tmp2,f);
      for (int j=0;j<Nc1;j++)
	res(j,i) = tmp1.dot(c1[j]);
    }
}

template<int N>
void overlap_kernel(taylor<3,N> &ovl,
		    double alpha1, double alpha2, const double *d)
{
  double asum = alpha1 + alpha2;
  taylor<3,N> d2;
  d2_taylor(d2,d);
  ovl = pow(M_PI/asum,3.0/2.0)*exp((-(alpha1*alpha2)/asum)*d2);
  dfactors(ovl);
}

template<int N>
void kinetic_kernel(taylor<3,N> &kin,
		    double alpha1, double alpha2, const double *d)
{
  double asum = alpha1 + alpha2;
  double h = alpha1*alpha2/asum;
  taylor<3,N> d2;
  d2_taylor(d2,d);
  kin = (pow(M_PI/asum,3.0/2.0)*h)*(3-(2*h)*d2)*exp((-h)*d2);
  dfactors(kin);
}

template<int N>
void coulomb_kernel(taylor<3,N> &k,
		    double alpha1, double alpha2, const double *d)
{
  double asum = alpha1 + alpha2;
  taylor<3,N> d2;
  d2_taylor(d2,d);
  k = (2*pow(M_PI,5.0/2.0)/(alpha1*alpha2*sqrt(asum)))*
    boys((alpha1*alpha2/asum)*d2);
  dfactors(k);
}

template<int N1, int N2, int Nc1, int Nc2>
void overlap(matrix<double,Nc1,Nc2> &S, 
	     const taylor<3,N1> c1[], 
	     const shell &s1,
	     const taylor<3,N2> c2[],
	     const shell &s2)
{
  taylor<3,N1+N2> ovl;
  double d[3];
  for (int i=0;i<3;i++)
    d[i] = s1.R[i] - s2.R[i];
  overlap_kernel(ovl, s1.alpha, s2.alpha, d);
  contract_d12(S,c1,c2,ovl);
}

template<int N1, int N2, int Nc1, int Nc2>
void kinetic(matrix<double,Nc1,Nc2> &T, 
	     const taylor<3,N1> c1[], 
	     const shell &s1,
	     const taylor<3,N2> c2[],
	     const shell &s2)
{
  taylor<3,N1+N2> kin;
  double d[3];
  for (int i=0;i<3;i++)
    d[i] = s1.R[i] - s2.R[i];
  kinetic_kernel(kin, s1.alpha, s2.alpha, d);
  contract_d12(T,c1,c2,kin);
}

template<int N1, int N2, int Nc1, int Nc2>
void coulomb(matrix<double,Nc1,Nc2> &C, 
	     const taylor<3,N1> c1[], 
	     const shell &s1,
	     const taylor<3,N2> c2[],
	     const shell &s2)
{
  taylor<3,N1+N2> kernel;
  double d[3];
  for (int i=0;i<3;i++)
    d[i] = s1.R[i] - s2.R[i];
  coulomb_kernel(kernel, s1.alpha, s2.alpha, d);
  contract_d12(C,c1,c2,kernel);
}

/* 
   Create shell coefficients for normalized cartesian or spherical
   gaussians.
*/
template<int Ndeg>
void cartesian_function(taylor<3,Ndeg> &fn, const shell &s, int ijk[3])
{
  matrix<double,1,1> S;
  fn.zero();
  fn[poly3_term_index[ijk[0]][ijk[1]][ijk[2]]] = 1;
  overlap(S,&fn,&fn,s,s);
  fn[poly3_term_index[ijk[0]][ijk[1]][ijk[2]]] /= sqrt(S(0,0));
}
template<int Ndeg>
void spherical_function(taylor<3,Ndeg> &fn, const shell &s, int L, int M)
{
  assert(L <= Ndeg);
  //  assert(abs(M) <= L);
}


template<int Ntot>
void overlap_6d(taylor<6,Ntot> &ovl, 
		double alpha1, double alpha2, const double *d)
{
  double asum = alpha1 + alpha2;
  // Inefficient like hell:
  taylor<6,Ntot> d1(d[0],0), d2(d[1],1), d3(d[2],2);
  d1[4] = -1;
  d2[5] = -1;
  d3[6] = -1;
  taylor<6,Ntot> dsq = pow(d1,2) + pow(d2,2) + pow(d3,2);
  ovl = pow(M_PI/asum,3.0/2.0)*exp(-(alpha1*alpha2)/asum*dsq);
  for (int i=0;i<ovl.size();i++)
    ovl[i] *= fac6[i];
}


template<int Ntot>
void kinetic_6d(taylor<6,Ntot> &kin, 
		double alpha1, double alpha2, const double *d)
{
  double asum = alpha1 + alpha2;
  // Inefficient like hell:
  taylor<6,Ntot> d1(d[0],0), d2(d[1],1), d3(d[2],2);
  d1[4] = -1;
  d2[5] = -1;
  d3[6] = -1;
  taylor<6,Ntot> dsq = pow(d1,2) + pow(d2,2) + pow(d3,2);
  double h = alpha1*alpha2/asum;

  kin = (h*(3-2*h*dsq))*pow(M_PI/asum,3.0/2.0)*exp(-(alpha1*alpha2)/asum*dsq);
  for (int i=0;i<kin.size();i++)
    kin[i] *= fac6[i];
}

void check_funvalue(void)
{
  shell s;
  s.alpha = 7;
  s.R[0] = 0.1;
  s.R[1] = -0.2;
  s.R[2] = 0.3;
  double r[3] = {0.4,0.5,0.6};
  taylor<3,0> c[1];
  c[0].zero();
  c[0][0] = 1;
  double f;
  // s function test
  funvalue(&f,c,1,s,r);
  double fcorr = exp(-s.alpha*(pow(r[0] - s.R[0],2)+
			       pow(r[1] - s.R[1],2)+
			       pow(r[2] - s.R[2],2)));
  if (fabs(f-fcorr)>1e-15)
    {
      cerr << "check_funvalue, ERROR in s function. expected " <<
	fcorr << " got " << f << endl;
    }
  taylor<3,1> cp[1];
  cp[0].zero();
  cp[0][1] = 1;
  // p function test
  funvalue(&f,cp,1,s,r);
  fcorr = -2*s.alpha*(r[0] - s.R[0])*exp(-s.alpha*(pow(r[0] - s.R[0],2)+
						   pow(r[1] - s.R[1],2)+
						   pow(r[2] - s.R[2],2)));
  if (fabs(f-fcorr)>1e-15)
    {
      cerr << "check_funvalue, ERROR in p function. expected " <<
	fcorr << " got " << f << endl;
    }
  cerr << "End of check_funval." << endl;
}

void check_overlap(void)
{
  // First slow method
  taylor<6,5> ovl6d;
  shell s1,s2;
  double d[3];
  s1.alpha = 2;
  s1.R[0] = 1.2;
  s1.R[1] = -0.2;
  s1.R[2] = 0.5;
  s2.alpha = 0.3;
  s2.R[0] = 0.8;
  s2.R[1] = -0.5;
  s2.R[2] = 0.1;
  for (int i=0;i<3;i++)
    d[i] = s1.R[i] - s2.R[i];
  overlap_6d(ovl6d,s1.alpha,s2.alpha,d);
  
  // compare to "fast" method
  int ncheck = 0;
  int e6[6];
  taylor<3,2> c1[1];
  taylor<3,3> c2[1];
  matrix<double,1,1> S;
  for (int i=0;i<c1[0].size();i++)
    for (int j=0;j<c2[0].size();j++)
      {
	c1[0].zero();
	c1[0][i] = 1;
	c2[0].zero();
	c2[0][j] = 1;
	overlap(S,c1,s1,c2,s2);
	c1[0].exponents(i,e6);
	c2[0].exponents(j,e6+3);
	if (fabs(S(0,0) - ovl6d[ovl6d.term_index(e6)]) > 1e-15)
	  {
	    cerr << "Error in overlap, expected " << 
	      ovl6d[ovl6d.term_index(e6)] <<
	      " got " << S(0,0) << endl;
	  }
	ncheck++;
      }
  cout << "End of check overlap. Checked " << ncheck << endl;
}


void check_kinetic(void)
{
  // First slow method
  taylor<6,5> kin6d;
  shell s1,s2;
  double d[3];
  s1.alpha = 2;
  s1.R[0] = 1.2;
  s1.R[1] = -0.2;
  s1.R[2] = 0.5;
  s2.alpha = 0.3;
  s2.R[0] = 0.8;
  s2.R[1] = -0.5;
  s2.R[2] = 0.1;
  for (int i=0;i<3;i++)
    d[i] = s1.R[i] - s2.R[i];
  kinetic_6d(kin6d,s1.alpha,s2.alpha,d);
  
  // compare to "fast" method
  int ncheck = 0;
  int e6[6];
  taylor<3,2> c1[1];
  taylor<3,3> c2[1];
  matrix<double,1,1> S;
  for (int i=0;i<c1[0].size();i++)
    for (int j=0;j<c2[0].size();j++)
      {
	c1[0].zero();
	c1[0][i] = 1;
	c2[0].zero();
	c2[0][j] = 1;
	kinetic(S,c1,s1,c2,s2);
	c1[0].exponents(i,e6);
	c2[0].exponents(j,e6+3);
	if (fabs(S(0,0) - kin6d[kin6d.term_index(e6)]) > 1e-15)
	  {
	    cerr << "Error in kinetic, expected " << 
	      kin6d[kin6d.term_index(e6)] <<
	      " got " << S(0,0) << endl;
	  }
	ncheck++;
      }
  cout << "End of check kinetic. Checked " << ncheck << endl;
}


void check_shell_product(void)
{
  shell s1,s2;
  s1.alpha = 2;
  s1.R[0] = 0;//1.2;
  s1.R[1] = 0;//-0.2;
  s1.R[2] = 0;//0.5;
  s2.alpha = 0.3;
  /*  s2.R[0] = 0.8;
  s2.R[1] = -0.5;
  s2.R[2] = 0.1;*/
  s2.R[0] = s1.R[0];
  s2.R[1] = s1.R[1];
  s2.R[2] = s1.R[2];
  taylor<3,2> c1[polymul_internal::binomial<3+2,2>::value];
  taylor<3,3> c2[polymul_internal::binomial<3+3,3>::value];
  matrix<taylor<3,5>,polymul_internal::binomial<3+2,2>::value,
    polymul_internal::binomial<3+3,3>::value> c12;
  shell s12;
  for (int i=0;i<c1[0].size();i++)
    {
      c1[i].zero();
      c1[i][i] = 1;
    }
  for (int i=0;i<c2[0].size();i++)
    {
      c2[i].zero();
      c2[i][i] = 1;
    }
  shell_products(c12,s12,c1,s1,c2,s2);
  /* Evaluate and compare at some points */
  for (int i=0;i<c1[0].size();i++)
    for (int j=0;j<c2[0].size();j++)
  for (double x = -1; x < 1; x+= 0.21)
    for (double y = -1; y < 1; y+= 0.22)
      for (double z = -1; z < 1; z+= 0.23)
	{
	  double r[] = {x,y,z};
	  double f1,f2,f12;
	  funvalue(&f1,&c1[i],1,s1,r);
	  funvalue(&f2,&c2[j],1,s2,r);
	  funvalue(&f12,&c12(i,j),1,s12,r);
	  if (fabs(f1*f2 - f12) > 1e-15)
	    cerr << "Checking shell products (" << 0 << ", " << 0 << ") expected " <<
	      f1*f2 << ", got " << f12 << endl;
	}
  cerr << "End of check_shell_products\n";
}


void evalbench(void)
{
  double fval[5],fsum = 0;
  taylor<3,3> forbs[5];
  shell s;
  double r[3] = {4,5,6};
  s.alpha = 3;
  s.R[0] = 1;
  s.R[1] = 2;
  s.R[2] = 3;

  // fill with randomness
  for (int i=0;i<5;i++)
    for (int j=0;j<forbs[i].size();j++)
      forbs[i][j] = cos(i*1.5+j);
  for (int i=0;i<1e6;i++)
    {
      funvalue(fval,forbs,5,s,r);
      fsum += fval[0]+fval[1];
    }
  cout << fsum;
}

void coulombbench(void)
{
  taylor<3,4> ddorbs[25],ddrhs;
  shell s1,s2;
  matrix<double,25,1> J;
  s1.alpha = 3;
  s1.R[0] = 1;
  s1.R[1] = 2;
  s1.R[2] = 3;
  s2.alpha = 2;
  s2.R[0] = 1;
  s2.R[1] = 2;
  s2.R[2] = 3;

  // fill with randomness
  for (int i=0;i<25;i++)
    for (int j=0;j<ddorbs[i].size();j++)
      ddorbs[i][j] = cos(i*1.5+j);
  for (int j=0;j<ddrhs.size();j++)
    ddrhs[j] = j;
  double fsum = 0;
  for (int i=0;i<1e6;i++)
    {
      coulomb(J,ddorbs,s1,&ddrhs,s2);
      fsum += J(0,0);
    }
  cout << fsum;
}


int main(void)
{
  coulombbench();
  return 0;
}

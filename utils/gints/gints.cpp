#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include "taylor.h"
#include "matrix.h"
#include "boys.h"

/*
  One and two electron program, by Ulf Ekstrom 2009
 */

using namespace std;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define MAX_L 3

class shell
{
public:
  double alpha;
  double R[3];
};

// Number of coefficients for s,p,d,.. shells 
const int shell_size[2*MAX_L+1] = {taylor<3,0>::size,taylor<3,1>::size,
				   taylor<3,2>::size,taylor<3,3>::size,
                                   taylor<3,4>::size,taylor<3,5>::size,
				   taylor<3,6>::size};

tmatrix<int,taylor<3,4*MAX_L>::size,3> poly3_exponents;
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
  for (int i = 0;i<p.size;i++)
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
  for (int i=0;i<p.size;i++)
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

// Take a column from matrix C and interpret it as a taylor series.
template<int N>
inline taylor<3,N> &column_as_taylor(matrix<double> &C, int col)
{
  assert(C.rows == (taylor<3,N>::size));
  return *reinterpret_cast<taylor<3,N>*>(C.column_ptr(col));
}

template<int N>
inline const taylor<3,N> &column_as_taylor(const matrix<double> &C, int col)
{
  assert(C.rows == (taylor<3,N>::size));
  return *reinterpret_cast<const taylor<3,N>*>(C.column_ptr(col));
}

// Multiply by  i!j!.. to each term, to get derivatives from
// Taylor coefficients
template<int N>
void dfactors(taylor<3,N> &t)
{
  assert(N <= 4*MAX_L); 
  for (int i=0;i<t.size;i++)
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
  for (int i=10;i<d2.size;i++)
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

// Given a function that depends only on the distance
// between two shells, contract over vectors on these
// shells. f is the taylor expansion of the function,
// in the vector d[3], the distance between the shells.
// Application to overlap, kinetic energy and Coulomb
template<int N1, int N2>
void contract_d12(matrix<double> &res,
		  const matrix<double> &c1,
		  const matrix<double> &c2,
		  const taylor<3,N1+N2> &f)
{
  taylor<3,N1> tmp1;
  taylor<3,N2> tmp2;
  for (int i=0;i<res.cols;i++)
    {
      tmp2 = column_as_taylor<N2>(c2,i);
      tmp2.invert_parity();
      polycontract(tmp1,tmp2,f);
      for (int j=0;j<res.rows;j++)
	res(j,i) = tmp1.dot(column_as_taylor<N1>(c1,j));
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
  taylor<3,MIN(N,2)> d2;
  d2_taylor(d2,d);
  d2 *= alpha1*alpha2/asum;
  boys_opt(k,d2);
  k *= 2*pow(M_PI,5.0/2.0)/(alpha1*alpha2*sqrt(asum));
  dfactors(k);
}

template<int N1, int N2>
void overlap(matrix<double> &S, 
	     const matrix<double> &c1, 
	     const shell &s1,
	     const matrix<double> &c2,
	     const shell &s2)
{
  taylor<3,N1+N2> ovl;
  double d[3];
  for (int i=0;i<3;i++)
    d[i] = s1.R[i] - s2.R[i];
  overlap_kernel(ovl, s1.alpha, s2.alpha, d);
  contract_d12<N1,N2>(S,c1,c2,ovl);
}

template<int N1, int N2>
void kinetic(matrix<double> &T, 
	     const matrix<double> &c1, 
	     const shell &s1,
	     const matrix<double> &c2,
	     const shell &s2)
{
  taylor<3,N1+N2> kin;
  double d[3];
  for (int i=0;i<3;i++)
    d[i] = s1.R[i] - s2.R[i];
  kinetic_kernel(kin, s1.alpha, s2.alpha, d);
  contract_d12<N1,N2>(T,c1,c2,kin);
}

template<int N1, int N2>
void coulomb(matrix<double> &C, 
	     const matrix<double> &c1, 
	     const shell &s1,
	     const matrix<double> &c2,
	     const shell &s2)
{
  taylor<3,N1+N2> kernel;
  double d[3];
  for (int i=0;i<3;i++)
    d[i] = s1.R[i] - s2.R[i];
  coulomb_kernel(kernel, s1.alpha, s2.alpha, d);
  contract_d12<N1,N2>(C,c1,c2,kernel);
}

template<int N>
void pointcharge(matrix<double> &C,
		 const matrix<double> &c1, 
		 const shell &s,
		 const double r[3],
		 double Z)
{
  assert(C.rows == c1.cols);
  assert(c1.rows == (taylor<3,N>::size));
  double d[3] = {s.R[0]-r[0],
		 s.R[1]-r[1],
		 s.R[2]-r[2]};
  taylor<3,MIN(N,2)> d2;
  taylor<3,N> k;
  d2_taylor(d2,d);
  d2 *= s.alpha;
  boys_opt(k,d2);
  k *= -2*Z*M_PI/s.alpha;
  dfactors(k);
  for (int i=0;i<c1.cols;i++)
    C(i,0) = k.dot(column_as_taylor<N>(c1,i));
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
  for (int i=0;i<c.size;i++)
    c[i] *= normfac[poly3_exponents(i,2)];
  dfactors(c);
}

/* Given vectors on two shells s1 and s2, compute the new shell
   s12 of the product function and the expansion of the product
   at that shell. One expansion is computed for each pairs of vectors
   from c1 and c2. */
template<int N1, int N2>
void shell_product(matrix<double> &prod, shell &s12,
		   const matrix<double> &c1, const shell &s1,
		   const matrix<double> &c2, const shell &s2)
{
  //  cerr << "This is shell_product<" << N1 << ", " << N2 << ">\n";
  assert(prod.rows == (taylor<3,N1+N2>::size));
  assert(c1.rows == (taylor<3,N1>::size));
  assert(c2.rows == (taylor<3,N2>::size));
  assert(prod.cols == c1.cols*c2.cols);

  s12.alpha = s1.alpha + s2.alpha;
  for (int i=0;i<3;i++)
    s12.R[i] = (s1.alpha*s1.R[i] + s2.alpha*s2.R[i])/s12.alpha;
  taylor<3,2*(N1+N2)> c[3];
  gauss_product_1d(c[0],s1.R[0],s2.R[0],s1.alpha,s2.alpha);
  gauss_product_1d(c[1],s1.R[1],s2.R[1],s1.alpha,s2.alpha);
  gauss_product_1d(c[2],s1.R[2],s2.R[2],s1.alpha,s2.alpha);

  prod.zero();
  int dcol = 0;
  for (int i = 0; i < c1.cols; i++)
    for (int j = 0; j < c2.cols; j++)
      {
	int exponents[3];
	for (int k = 0; k < taylor<3,N1+N2>::size; k++)
	  for (int m = 0;m < taylor<3,N1>::size; m++)
	    for (int n = 0;n < taylor<3,N2>::size; n++)
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
		prod(k,dcol) += w*c1(m,i)*c2(n,j);
	      }
	dcol++;
      }
}

/* 
   Create shell coefficients for normalized cartesian or spherical
   gaussians. In principle not difficult to make this generic. See
   inverse Weierstrass transform.
*/
template<int Ndeg>
void cartesian_function_unnormalized(taylor<3,Ndeg> &fn,
				     const shell &s, int ijk[3])
{
  assert(ijk[0] + ijk[1] + ijk[2] <= Ndeg);  
  double c1d[3][Ndeg+1];
  for (int i=0;i<3;i++)
    for (int j=0;j<=Ndeg;j++)
      c1d[i][j] = 0;
  fn = 0;
  for (int i=0;i<3;i++)
    switch(ijk[i])
      {
      case 0: // s
	c1d[i][0] = 1;
	break;
      case 1: // p
	c1d[i][1] = 1/(2*s.alpha);
	break;
      case 2: // d
	c1d[i][0] = 1/(2*s.alpha);
	c1d[i][2] = 1/pow(2*s.alpha,2);
	break;
      case 3: // f
	if (Ndeg>2) //shuts up the compiler about array bounds
	  {
	    c1d[i][1] = 3/pow(2*s.alpha,2);
	    c1d[i][3] = 1/pow(2*s.alpha,3);
	  }
	break;
      default:
	assert(0 && "Unimplemented cartesian -> hermite order :-(");
	break;
      }
  for (int i=0;i<=Ndeg;i++)
    for (int j=0;j<=Ndeg-i;j++)
      for (int k=0;k<=Ndeg-i-j;k++)
	fn[poly3_term_index[i][j][k]] = c1d[0][i]*c1d[1][j]*c1d[2][k];
}

template<int Ndeg>
void cartesian_function(taylor<3,Ndeg> &fn, const shell &s, int ijk[3])
{
  cartesian_function_unnormalized(fn, s, ijk);
  matrix<double> S(1,1);
  matrix<double> c(&fn[0],fn.size,1,0); // Wrap fn in a matrix
  overlap<Ndeg,Ndeg>(S,c,s,c,s);
  fn *= 1/sqrt(S(0,0));
}

#include "call_tables.inc"

class shell_basis : public matrix<double>
{
public:
  shell_basis(const double R[3], double alpha, int order, int nvec) 
    : matrix<double>(shell_size[order],nvec)
  {
    assert(order <= 2*MAX_L);
    this->order = order;
    s.alpha = alpha;
    s.R[0] = R[0];
    s.R[1] = R[1];
    s.R[2] = R[2];
  }
  // Make cartesian functions for now
  shell_basis(const double R[3], double alpha, int order) 
    : matrix<double>(shell_size[order],
	     order > 0 ? shell_size[order] - shell_size[order-1] : 1)
  {
    assert(order <= MAX_L);
    this->order = order;
    s.alpha = alpha;
    s.R[0] = R[0];
    s.R[1] = R[1];
    s.R[2] = R[2];    
    int ijk[3] = {order,0,0};
    for (int i=0;i<cols;i++)
      {
	switch(order)
	  {
	  case 0:
	    cartesian_function(column_as_taylor<0>(*this,i), s, ijk);
	    break;
	  case 1:
	    cartesian_function(column_as_taylor<1>(*this,i), s, ijk);
	    break;
	  case 2:
	    cartesian_function(column_as_taylor<2>(*this,i), s, ijk);
	    break;
	  case 3:
	    cartesian_function(column_as_taylor<3>(*this,i), s, ijk);
	    break;
	  default:
	    assert(0 && "Unimplemented order");
	  };
	taylor<3,MAX_L>::next_exponents(3,ijk);
      }
  }
  // construct a product basis for all functions on b1 and b2
  shell_basis(const shell_basis &b1, const shell_basis &b2) 
    : matrix<double>(shell_size[b1.order+b2.order],b1.cols*b2.cols)
  {
    order = b1.order+b2.order;
    assert(order <= 2*MAX_L);
    shell_product_tab[b1.order][b2.order](*this,s,b1,b1.s,b2,b2.s);
  }
  int nr_functions(void) const
  {
    return cols;
  }
  double exponent(void) const
  {
    return s.alpha;
  }
  void overlap(matrix<double> &S, const shell_basis &b2)
  {
    overlap_tab[order][b2.order](S,*this,s,b2,b2.s);
  }
  void kinetic(matrix<double> &S, const shell_basis &b2)
  {
    kinetic_tab[order][b2.order](S,*this,s,b2,b2.s);
  }
  void coulomb(matrix<double> &S, const shell_basis &b2)
  {
    coulomb_tab[order][b2.order](S,*this,s,b2,b2.s);
  }
  void point_charge(matrix<double> &S, const double R[3], double Z)
  {
    pointcharge_tab[order](S,*this,s,R,Z);
  }
protected:
  int order;
  shell s;
};

class basis
{
public:
  basis(void) 
  { 
    nr_functions = 0; 
  }
  void add_shell(const double R[3], double alpha, int L)
  {
    shell_basis *b = new shell_basis(R,alpha,L);
    shells.push_back(b);
    nr_functions += b->cols;
  }
  void overlap(matrix<double> &S, const basis &b2)
  {
    assert(S.rows == nr_functions);
    assert(S.cols == b2.nr_functions);
    int ic = 0;
    for (unsigned int i=0;i<shells.size();i++)
      {
	int jc = 0;
	for (unsigned int j=0;j<b2.shells.size();j++)
	  {
	    matrix<double> S12 = 
	      S.block(ic,jc,shells[i]->cols,b2.shells[j]->cols);
	    shells[i]->overlap(S12,*b2.shells[j]);
	    jc += b2.shells[j]->cols;
	  }
	ic += shells[i]->cols;
      }
  }
  void kinetic(matrix<double> &S, const basis &b2)
  {
    assert(S.rows == nr_functions);
    assert(S.cols == b2.nr_functions);
    int ic = 0;
    for (unsigned int i=0;i<shells.size();i++)
      {
	int jc = 0;
	for (unsigned int j=0;j<b2.shells.size();j++)
	  {
	    matrix<double> S12 = 
	      S.block(ic,jc,shells[i]->cols,b2.shells[j]->cols);
	    shells[i]->kinetic(S12,*b2.shells[j]);
	    jc += b2.shells[j]->cols;
	  }
	ic += shells[i]->cols;
      }
  }  
  void read_from(istream &src)
  {
    double R[3],alpha;
    int L;
    while (!src.fail())
      {
	src >> L >> R[0] >> R[1] >> R[2] >> alpha;
	if (!src.fail())
	  add_shell(R,alpha,L);
      }
  }
  int nr_funs(void)
  {
    return nr_functions;
  }
  //protected:
  vector<shell_basis *> shells;
  int nr_functions;
};


class product_basis : public basis
{
public:
  product_basis(const basis &b1, const basis &b2) : basis()
  { 
    nr_functions = 0;
    if (&b1 == &b2)
      {
	int i1 = 0;
	for (unsigned int i=0;i<b1.shells.size();i++)
	  {
	    int i2 = 0;
	    for (unsigned int j=0;j<=i;j++)
	      {
		shell_basis *b = new shell_basis(*b1.shells[i],*b2.shells[j]);
		shells.push_back(b);
		nr_functions += b->cols;
		for (int ii=i1;ii<i1+b1.shells[i]->cols;ii++)
		  for (int jj=i2;jj<i2+b2.shells[j]->cols;jj++)
		    origin.push_back(make_pair(ii,jj));
		i2 += b2.shells[j]->cols;
	      }
	    i1 += b1.shells[i]->cols;
	  }
      }
    else
      {
	int i1 = 0;
	for (unsigned int i=0;i<b1.shells.size();i++)
	  {
	    int i2 = 0;
	    for (unsigned int j=0;j<b2.shells.size();j++)
	      {
		shell_basis *b = new shell_basis(*b1.shells[i],*b2.shells[j]);
		shells.push_back(b);
		nr_functions += b->cols;
		for (int ii=i1;ii<i1+b1.shells[i]->cols;ii++)
		  for (int jj=i2;jj<i2+b2.shells[j]->cols;jj++)
		    origin.push_back(make_pair(ii,jj));
		i2 += b2.shells[j]->cols;
	      }
	    i1 += b1.shells[i]->cols;
	  }
      }
  }
  // Given a list of point charges and positions
  // Q1 x1 y1 z1
  // Q2 x2 y2 z2 
  // ...
  // Calculate the matrix between the elements of the originating 
  // basis sets.
  void point_charges(matrix<double> &F, 
		     const matrix<double> &charge_and_pos)
  {
    assert(charge_and_pos.cols == 4);
    F.zero();
    int ic = 0;   
    for (unsigned int i=0;i<shells.size();i++)
      {
	for (int j=0;j<charge_and_pos.rows;j++)
	  {
	    double R[3] = {charge_and_pos(j,1),
			   charge_and_pos(j,2),
			   charge_and_pos(j,3)};
	    matrix<double> A(shells[i]->cols, 1);
	    shells[i]->point_charge(A, R, charge_and_pos(j,0));
	    for (int k=0;k<shells[i]->cols;k++)
	      F(origin[ic+k]) += A(k,0);
	  }
	ic += shells[i]->cols;
      }
  }
  // Computes integral_psi(1) r_12^-1 psi(2). If psi1
  // and psi2 are overlaps distributions is this two electron
  // integrals (and a lot of them!)
  // b2 is the basis the density matrix lives on.
  void coulomb(matrix<double> &J, const basis &b2)
  {
    assert(J.rows == nr_functions);
    assert(J.cols == b2.nr_functions);
    int ic = 0;
    for (unsigned int i=0;i<shells.size();i++)
      {
	int jc = 0;
	unsigned int j;
#if 0 //skip this optimization for now
	if (&b2 == this)
	  j = i;
	else
#endif
	  j = 0;
	for (;j<b2.shells.size();j++)
	  {
	    matrix<double> J12 = 
	      J.block(ic,jc,shells[i]->cols,b2.shells[j]->cols);
	    shells[i]->coulomb(J12,*b2.shells[j]);
	    jc += b2.shells[j]->cols;
	  }
	ic += shells[i]->cols;
      }
  }
  // Computes integral_psi(1) r_12^-1 psi(2). If psi1
  // and psi2 are overlaps distributions is this two electron
  // integrals (and a lot of them!)
  // b2 is the basis the density matrix lives on.
  void coulomb_contract(matrix<double> &J,
			const matrix<double> &D, 
			const product_basis &b2)
  {
    //TODO: optimize if this == &b2
    J.zero();
    int ic = 0;   
    for (unsigned int i=0;i<shells.size();i++)
      {
	int jc = 0;
	for (unsigned int j=0;j<b2.shells.size();j++)
	  {
	    matrix<double> J12(shells[i]->cols, b2.shells[j]->cols);
	    shells[i]->coulomb(J12, *b2.shells[j]);
	    for (int k=0;k<shells[i]->cols;k++)
	      for (int m=0;m<b2.shells[j]->cols;m++)
		J(origin[ic+k]) += J12(k,m)*D(b2.origin[jc+m]);
	    jc += b2.shells[j]->cols;
	  }
	ic += shells[i]->cols;
      }
  }
  pair<int,int> get_origin(int i)
  {
    return origin.at(i);
  }
  //protected:
  // For each function, tell where the product came from
  vector< pair<int,int> > origin;
};



void print_usage(ostream &dst)
{
  dst << "Usage: gints -b BASIS [-bb BASIS] [-o|-k|-n CHARGES]\n";
  dst << "Options:\n";
  dst << " -o    Overlap\n";
  dst << " -k    Kinetic energy\n";
  dst << " -n    Nuclear attraction (point charges)\n";
  dst << "       Charges are read from CHARGES, in format Q X Y Z\n";
  dst << " -c    Two electron Coulomb (and exchange) integrals\n";
  dst << "If a second basis set file is give use this as\n";
  dst << "the bra ('right') basis.\n";
  dst << "The basis set file has a format\n";
  dst << "L X Y Z EXP\n";
  dst << "for each shell. L is the angular momentum nr.\n";
}

int main(int argc, char *argv[])
{
  // Parse command line options
  string type("none");
  string point_charge_path;
  string basis_path;
  string second_basis_path("");
  for (int i=1;i<argc;i++)
    {
      if (argv[i] == string("-h"))
	{
	  print_usage(cout);
	  return EXIT_SUCCESS;
	}
      else if (argv[i] == string("-o"))
	type = "overlap";
      else if (argv[i] == string("-k"))
	type = "kinetic";
      else if (argv[i] == string("-c"))
	type = "coulomb";
      else if (argv[i] == string("-n"))
	{
	  type = "pointcharge";
	  if (argc>i+1)
	    {
	      point_charge_path = string(argv[i+1]);
	      i++;
	    }
	  else
	    {
	      cerr << "Expected FILE argument to -n, quitting.\n";
	      return EXIT_FAILURE;
	    }
	}
      else if (argv[i] == string("-b"))
	{
	  if (argc>i+1)
	    {
	      basis_path = string(argv[i+1]);
	      i++;
	    }
	  else
	    {
	      cerr << "Expected FILE argument to -b, quitting.\n";
	      return EXIT_FAILURE;
	    }
	}
      else if (argv[i] == string("-bb"))
	{
	  if (argc>i+1)
	    {
	      second_basis_path = string(argv[i+1]);
	      i++;
	    }
	  else
	    {
	      cerr << "Expected FILE argument to -bb, quitting.\n";
	      return EXIT_FAILURE;
	    }
	}
      else
	{
	  cerr << "Unknown option '" << argv[i] << "', try -h.\n";
	  return EXIT_FAILURE;
	}
    }
  if (type == "none")
    {
      cerr << "Nothing to calculate, quitting.\n";
      return EXIT_FAILURE;
    }
  if (basis_path == "")
    {
      cerr << "No basis given, quitting.\n";
      return EXIT_FAILURE;
    }

  setup_tables();

  basis b;
  {
    ifstream basfile(basis_path.c_str());
    if (basfile.fail())
      {
	cerr << "Error opening basis file " << basis_path << ", quitting.\n";
	return EXIT_FAILURE;
      }
    b.read_from(basfile);
    cerr << "Read basis: " << b.nr_funs() << " functions, (" << 
      b.shells.size() << " shells).\n";
  }
  basis *b2;
  if (second_basis_path != "")
    {
      b2 = new basis();
      ifstream basfile(second_basis_path.c_str());
      if (basfile.fail())
	{
	  cerr << "Error opening basis file " << second_basis_path 
	       << ", quitting.\n";
	  return EXIT_FAILURE;
	}
      b2->read_from(basfile);
      cerr << "Read second basis: " << b2->nr_funs() << " functions, (" 
	   << b2->shells.size() << " shells).\n";
    }
  else
    {
      b2 = &b;
    }
  basis &bb = *b2;

  cout.precision(16);
  cout << scientific;
  
  if (type == "overlap")
    {
      matrix<double> S(b.nr_funs(),bb.nr_funs());
      b.overlap(S,bb);
      S.print(cout);
    }
  else if (type == "kinetic")
    {
      matrix<double> S(b.nr_funs(),bb.nr_funs());
      b.kinetic(S,bb);
      S.print(cout);
    }
  else if (type == "pointcharge")
    {
      vector<double> q;
      ifstream chargef(point_charge_path.c_str());
      if (chargef.fail())
	{
	  cerr << "gints: Error opening point charge file '" 
	       << point_charge_path << "', quitting.\n";
	  return EXIT_FAILURE;
	}
      double x;
      chargef >> x;
      while (!chargef.fail())
	{
	  q.push_back(x);
	  chargef >> x;
	}
      if (q.size() % 4 != 0 || !chargef.eof())
	{
	  cerr << "gints: Invalid data in point charge file '" 
	       << point_charge_path << "', quitting.\n";
	  return EXIT_FAILURE;
	}
      matrix<double> charges(q.size()/4,4);
      int k = 0;
      for (unsigned int i=0;i<q.size()/4;i++)
	for (int j=0;j<4;j++)
	  charges(i,j) = q[k++];
      product_basis prodbas(b,bb);
      matrix<double> F(b.nr_funs(),bb.nr_funs());
      prodbas.point_charges(F,charges);
      if (&b == &bb)
	F.symmetrize_from_lower();
      F.print(cout);
    }
  else if (type == "coulomb")
    {
      product_basis pbas(b,bb);
      matrix<double> J(pbas.nr_funs(),pbas.nr_funs());
      pbas.coulomb(J,pbas);
      for (int i=0;i<J.rows;i++)
	for (int j=0;j<J.cols;j++)
	  {
	    pair<int,int> ij1 = pbas.get_origin(i), ij2 = pbas.get_origin(j);
	    // 1-based indexing in the output.
	    cout << ij1.first+1 << " " << ij1.second+1 << " "
		 << ij2.first +1<< " " << ij2.second+1 << " "
		 << J(i,j) << endl;
	  }
    }
  return EXIT_SUCCESS;
}

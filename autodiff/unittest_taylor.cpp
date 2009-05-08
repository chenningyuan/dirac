#include <iostream>
#include "taylor.h"

using namespace std;

// Some taylor coefficients around x = 2.5, calculated using
// Maxima bfloat(taylor(...))
const double exp_good[] = {1.218249396070347e1,
			   1.218249396070347e1,
			   6.091246980351737e0,
			   2.030415660117246e0,
			   5.076039150293114e-1,
			   1.015207830058623e-1,
			   1.692013050097705e-2};

const double log_good[] = {9.162907318741551e-1,
			   4.0e-1,
			   -8.0e-2,
			   2.133333333333333e-2,
			   -6.4e-3,
			   2.048e-3,
			   -6.826666666666667e-4};

const double atan_good[] = {1.1902899496825317,
			    1.379310344827586e-1,
			    -4.756242568370987e-2,
			    1.552612516571679e-2,
			    -4.750587107528691e-3,
			    1.336092873197889e-3,
			    -3.310338712605161e-4};


//exp(1/(2.5 + 3x + 7y))
const double composed_good[] = { 1.49182469764127,
				-7.160758548678098e-1, //x
				-1.670843661358223, //y
				 1.031149231009646, //x^2
				 4.812029744711682, //xy
				 5.614034702163628}; //y^2

#define NR_COEFF_CHECK 6

template<int Nvar, int Ndeg>
int taylor_check(const char *label,
		  const taylor<Nvar,Ndeg>& t,
		  const double correct[],
		  double thres)
{
  int nfail = 0;
  for (int i=0;i<NR_COEFF_CHECK;i++)
    {
      double err = 2*fabs(t[i] - correct[i])/(1 + fabs(correct[i]));
      if (err > thres)
	{
	  cout << "Error in coefficient " << i << " of " << label << ". Correct:" << correct[i] << ", taylor:" << t[i] << " error: " <<  correct[i] - t[i] << endl;
	  nfail++;
	}
    }
  return nfail;
}

template<int Nvar, int Ndeg>
int taylor_compare(const taylor<Nvar,Ndeg>& t1,
		   const taylor<Nvar,Ndeg>& t2,
		   double thres)
{
  int nfail = 0;
  for (int i=0;i<t1.size();i++)
    {
      double err = 2*fabs(t1[i] - t2[i])/(1 + 0.5*fabs(t1[i]+t2[i]));
      if (err > thres)
	{
	  cout << "Difference in coefficient " << i << ". Correct:" << t1[i] << ", taylor:" << t2[i] << " difference: " <<  t1[i] - t2[i] << endl;
	  nfail++;
	}
    }
  return nfail;
}


template<int Nvar, int Ndeg>
double taylor_sumsq(const taylor<Nvar,Ndeg>& t)
{
  double s = 0;
  for (int i=0;i<t.size();i++)
    s += t[i]*t[i];
  return s;
}

double fac(int n)
{
  double f = 1;
  while (n > 1)
    f *= n--;
  return f;
}

int main(void)
{
  double x0 = 2.5;
  taylor<1,6> tin(x0,0);
  int res = 0;

  cout.precision(16);

  // Test some elementary functions (also test composition etc)
  res += taylor_check("exp",exp(tin),exp_good,1e-15);
  res += taylor_check("log",log(tin),log_good,1e-15);
  res += taylor_check("atan",atan(tin),atan_good,1e-15);

  // Test multidimensional multiply consistency
  taylor<4,5> tmul(0),tacc(1),t1,t2,t3;
  for (int i=0;i<t1.size();i++)
    {
      t1[i] = i+1;
      t2[i] = -3+i*i;
      t3[i] = i + i*i + 1;
    }
  // multiplication should commute
  tacc=t1;
  tacc*=t2;
  tacc*=t3;

  tmul=t2;
  tmul*=t1;
  tmul*=t3;
  res += taylor_compare(tacc,tmul,1e-16);
  // add and subtract a bit
  if (taylor_sumsq(tmul-tacc) > 1e-16)
    {
      cout << "Subtraction error\n";
      res++;
    }
  // Evaluate derivatives at 0 and compare to coefficients
  for (int i=0;i<tin.size();i++)
    {
      tin[i] = 2 + i + i*i;
      if (fabs(fac(i)*tin[i] - tin.eval(0,i)) > 1e-15)
	{
	  cout << "Error evaluating derivative " << i << endl;
	  res ++;
	}
    }
  // Evaluate a simple second order polynomial and check the results
  tin.zero();
  tin[0] = 7;
  tin[1] = 3;
  tin[2] = 23;
  if (fabs(tin.eval(0.0)-7) > 1e-15)
    {
      cout << "Error evaluating polynomial (i)\n";
      res++;
    }
  if (fabs(tin.eval(8.0)-1503) > 1e-15)
    {
      cout << "Error evaluating polynomial (ii), I got" << tin.eval(8) << endl;
      res++;
    }
  if (fabs(tin.eval(6,1)-279) > 1e-15)
    {
      cout << "Error evaluating polynomial (iii)\n";
      res++;
    }
  if (fabs(tin.eval(-2,1)+89) > 1e-15)
    {
      cout << "Error evaluating polynomial (iv), got " << tin.eval(-2,1) << endl;
      res++;
    }
  if (fabs(tin.eval(-2,2)-46) > 1e-15)
    {
      cout << "Error evaluating polynomial (v), got "  << tin.eval(-2,2) << endl;
      res++;
    }
#ifdef TODO
  // Evaluate multivariate polynomial and check the results
  {
    taylor<2,1> p;
    for (int i=0;i<p.size();i++)
      p[i] = i+1;
    double x[2] = {3,7.1};
    if (fabs(p.eval(x)-473.26) > 1e-15)
      {
	cout << "Error evaluating multivariate, got " << p.eval(x) << endl;
	res++;
      }
  }
#endif
  
  // composition
  taylor<2,2> p(2.5);
  p[1] = 3;
  p[2] = 7; // Now p = 2.5 + 3*x + 7*y
  taylor<2,2> ret = exp(1/p);
  res += taylor_check("exp(1/(2.5+3*x+7*y))",ret,composed_good,1e-15);

  // Term indexing
  int exponents1[] = {2,0,1};
  int exponents2[] = {0,1,1};
  /* 1 x y z x^2 xy xz y^2 yz z^2 x^3 x^2y x^2z xy^2 xyz xz^2 y^3 y^2z yz^2 z^3
     0 1 2 3 4   5  6  7   8  9   10  11   12   13   14  15   16  17   18  19*/
  if (taylor_index(3,exponents1) != 12 ||
      taylor_index(3,exponents2) != 8)
    {
      cout << "Error in taylor::index_of()\n";
      res++;
    }


  return res;
}

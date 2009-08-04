#include <iostream>
#include <fstream>
#include <cmath>
#include <qd/qd_real.h>

using namespace std;

static qd_real fac(int n)
{
  qd_real f = 1;
  while (n > 1)
    f *= n--;
  return f;
}

static qd_real boys_from_taylor(qd_real x, int n)
{
  qd_real xn = 1, sum = 0;
  int k = 0;
  while (1)
    {
      qd_real term = xn/(2*n+1+2*k);
      sum += term;
      if (fabs(term/sum) < 1e-53)
	return pow(-1,n)*sum;
      k++;
      xn = xn*(-x)/k;
    }
}

int main(void)
{
  unsigned int cw;
  fpu_fix_start(&cw); 
  int maxdeg = 20;
  int nsteps = 100;
  qd_real delta = 0.1, x = 0;
  ofstream header("boys_table.h");
  header << "#ifndef BOYS_TABLE_H"  << endl;
  header << "#define BOYS_TABLE_H"  << endl << endl;
  header << "#define BOYS_MAXDEG " << maxdeg << endl;
  header << "#define BOYS_STEPS " << nsteps << endl;
  header << "#define BOYS_DELTA " << delta << endl << endl;
  header << "extern const double boystab[BOYS_STEPS][BOYS_MAXDEG];\n\n";
  header << "#endif"  << endl;
  header.close();
  ofstream body("boys_table.cpp");
  body << "#include \"boys_table.h\"\n";
  body << "const double boystab[BOYS_STEPS][BOYS_MAXDEG] = {\n";
  body.precision(40);
  for (int i=0;i<nsteps;i++)
    {
      body << "{\n";
      for (int j=0;j<maxdeg;j++)
	body << boys_from_taylor(x,j) << ",\n";
      body << "},\n";
      x += delta;
    }
  body << "};\n";
  fpu_fix_end(&cw);
  cout << "Generated boys_table.cpp and boys_table.h\n";
  return 0;
}

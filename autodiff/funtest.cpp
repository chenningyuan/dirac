#include <iostream>
#include "taylor.h"
#include "auto_functionals.h"
using namespace std;

// Program to print functional derivatives in high 
// precision and normal precision.
//
// Compile with -DTAYLOR_QD to use quadratic double 
// precision. You need libqd for this to work, see taylor.h
//
// I built this program like
// c++ -Wall funtest.cpp taylor.cpp auto_functionals.cpp -lqd
// -Ulf

int main(void)
{
  taylor<2,5> res;
#ifdef TAYLOR_QD
  unsigned int cw;
  fpu_fix_start(&cw); 
  cout << "Using quadratic double precision (62-64 digits)" << endl;
#else
  cout << "Using standard double precision (15-16 digits)" << endl;
#endif
  // How many digits to print
  cout.precision(30);

  res = vwn_correlation<5>(0.13);
  for (int i=0;i<res.size();i++)
    cout << res[i] << endl;

#ifdef TAYLOR_QD
  fpu_fix_end(&cw);
#endif
  return 0;
}

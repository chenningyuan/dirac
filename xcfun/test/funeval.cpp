#include <iostream>
#include <cstdlib>
#include "xcfun.h"

using namespace std;

#define TEST_ORDER 2

int main(int argc, char *argv[])
{
  xc_functional fun;
  fun.set_parameter("slaterx",1.0);
  fun.set_mode(XC_VARS_AB);
  cout << "Type is " << fun.get_type() << endl;
  cout << "Input length " << fun.input_length() << endl;
  cout << "Output length " << fun.output_length(TEST_ORDER) << endl;
  cout << "Maximal supported order " << fun.get_max_order() << endl;
  if (fun.get_max_order() >= TEST_ORDER)
    {
      double dv[fun.input_length()];
      for (int i=0;i<fun.input_length();i++)
	cin >> dv[i];
      double res[fun.output_length(TEST_ORDER)];
      fun.eval(res,TEST_ORDER,dv);
      for (int i=0;i<fun.output_length(TEST_ORDER);i++)
	cout << res[i] << endl;
    }
  return 0;
}

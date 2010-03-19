#include <iostream>
#include <cstdlib>
#include "xcfun.h"

using namespace std;

int main(int argc, char *argv[])
{
  xc_functional fun;

  fun.set_parameter("slaterx",1.0);
  fun.set_parameter("vwn5c",1.0);
  
  fun.set_mode(XC_VARS_AB);
  cout << "Maximal supported order " << fun.get_max_order() << endl;

  cout << endl << "Functional:" << endl;
  cout << fun.describe_parameter("slaterx") << endl << endl;
  cout << fun.describe_parameter("vwn5c") << endl << endl; 

  int order = 2;

  // Allocate and read input
  cout << "Reading " << fun.input_length() << " density variables" << endl;
  double density[fun.input_length()];
  for (int i=0;i<fun.input_length();i++)
    cin >> density[i];

  // Allocate space for the output
  double result[fun.output_length(order)];

  // Calculate energy and derivatives
  fun.eval(result,order,density);

  // Print result
  cout << "Output (" << fun.output_length(order) << ") values" << endl;
  for (int i=0;i<fun.output_length(order);i++)
    cout << result[i] << endl;

  return 0;
}

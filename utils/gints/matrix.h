#ifndef MATRIX_H
#define MATRIX_H

#include <fstream>

template<class numtype, int nrows, int ncols>
class matrix
{
 public:
  matrix() {}
  void zero(void)
  {
    for (int i=0;i<ncols;i++)
      for (int j=0;j<nrows;j++)
	c[i][j] = 0;
  }
  numtype &operator()(int row, int col) 
    { 
      assert(col>=0 && col < ncols);
      assert(row>=0 && row < nrows);
      return c[col][row]; 
    }
  const numtype &operator()(int row, int col) const 
  { 
    assert(col>=0 && col < ncols);
    assert(row>=0 && row < nrows);
    return c[col][row]; 
  }
  void print(std::ostream &dst) const
  {
    for (int i=0;i<nrows;i++)
      {
	for (int j=0;j<ncols;j++)
	  dst << c[j][i] << " ";
	dst << "\n";
      }
  }
 protected:
  numtype c[ncols][nrows];
};

#endif

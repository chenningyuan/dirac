#ifndef MATRIX_H
#define MATRIX_H

#include <fstream>
#include <utility>

template<class numtype, int nrows, int ncols>
class tmatrix
{
 public:
  tmatrix() {}
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


template<class numtype>
class matrix
{
 public:
  matrix(int nrows, int ncols) 
  {
    assert(nrows > 0);
    assert(ncols > 0);
    rows = nrows;
    cols = ncols;
    pitch = nrows;
    owns_c = true;
    c = new numtype[rows*cols];
  }
  matrix(numtype *wrapped_data, int nrows, int ncols, int pitch) 
  {
    rows = nrows;
    cols = ncols;
    this->pitch = pitch;
    this->c = wrapped_data;
    owns_c = false;
  }
  ~matrix(void)
  {
    if (owns_c)
      delete[] c;
  }
  void zero(void)
  {
    for (int i=0;i<cols;i++)
      for (int j=0;j<rows;j++)
	c[i*pitch+j] = 0;
  }
  void eye(void)
  {
    for (int i=0;i<cols;i++)
      for (int j=0;j<rows;j++)
	c[i*pitch+j] = (i == j) ? 1 : 0;
  }
  void symmetrize_from_lower(void)
  {
    assert(rows == cols);
    for (int i=0;i<cols;i++)
      for (int j=i+1;j<cols;j++)
	c[j*pitch+i] = c[i*pitch+j];
  }
  numtype *column_ptr(int col) 
  { 
    assert(col>=0 && col < cols);
    return &c[col*pitch];
  }
  const numtype *column_ptr(int col) const
  { 
    assert(col>=0 && col < cols);
    return &c[col*pitch];
  }
  inline const numtype &get(int row, int col) const
  { 
    assert(col>=0 && col < cols);
    assert(row>=0 && row < rows);
    return c[col*pitch + row];
  }
  inline numtype &operator()(int row, int col) 
  { 
    assert(col>=0 && col < cols);
    assert(row>=0 && row < rows);
    return c[col*pitch + row];
  }
  inline const numtype &operator()(int row, int col) const 
  { 
    assert(col>=0 && col < cols);
    assert(row>=0 && row < rows);
    return c[col*pitch + row];
  }
  inline numtype &operator()(const std::pair<int,int> &element)
  { 
    assert(element.second>=0 && element.second < cols);
    assert(element.first>=0 && element.first < rows);
    return c[element.second*pitch + element.first];
  }
  inline const numtype &operator()(const std::pair<int,int> &element) const
  { 
    assert(element.second>=0 && element.second < cols);
    assert(element.first>=0 && element.first < rows);
    return c[element.second*pitch + element.first];
  }
  void print(std::ostream &dst) const
  {
    for (int i=0;i<rows;i++)
      {
	for (int j=0;j<cols;j++)
	  dst << (*this)(i,j) << " ";
	dst << "\n";
      }
  }
  matrix<numtype> block(int firstrow, int firstcol, int nrows, int ncols)
  {
    assert(firstcol>=0 && firstcol + ncols <= cols);
    assert(firstrow>=0 && firstrow + nrows <= rows);
    return matrix(c + firstcol*pitch + firstrow,nrows,ncols,pitch);
  }
  int rows, cols;
 protected:
  int pitch;
  numtype *c;
  bool owns_c;
};

#endif

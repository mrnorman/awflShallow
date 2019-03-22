
#ifndef _SARRAY_H_
#define _SARRAY_H_

#include <iostream>
#include <iomanip>

/*
  This is intended to be a simple, low-overhead class to do multi-dimensional arrays
  without pointer dereferencing. It supports indexing and cout only up to 3-D.

  It templates based on array dimension sizes, which conveniently allows overloaded
  functions in the TransformMatrices class.
*/

template <class T, unsigned long D0, unsigned long D1=1, unsigned long D2=1> class SArray {

  public :

  typedef unsigned long ulong;

  T data[D0*D1*D2];

  SArray() { }
  SArray(SArray &&in) {
    for (int i=0; i < D0*D1*D2; i++) { data[i] = in.data[i]; }
  }
  SArray &operator=(SArray &&in) {
    for (int i=0; i < D0*D1*D2; i++) { data[i] = in.data[i]; }
  }
  ~SArray() { }

  inline T &operator()(ulong const i0)       {
    return data[i0];
  }
  inline T &operator()(ulong const i0, ulong const i1)       {
    return data[i0*D1 + i1];
  }
  inline T &operator()(ulong const i0, ulong const i1, ulong const i2)       {
    return data[i0*D1*D2 + i1*D2 + i2];
  }

  inline T  operator()(ulong const i0) const {
    return data[i0];
  }
  inline T  operator()(ulong const i0, ulong const i1) const {
    return data[i0*D1 + i1];
  }
  inline T  operator()(ulong const i0, ulong const i1, ulong const i2) const {
    return data[i0*D1*D2 + i1*D2 + i2];
  }

  inline friend std::ostream &operator<<(std::ostream& os, SArray const &v) {
    for (ulong i=0; i<D0*D1*D2; i++) {
      os << std::setw(12) << v.data[i] << "\n";
    }
    os << "\n";
    return os;
  }

};

#endif

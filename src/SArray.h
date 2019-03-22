
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

  template <class I, ulong E0> inline SArray<T,E0> operator*(SArray<I,D0> const &rhs) {
    //This template could match either vector-vector or matrix-vector multiplication
    if ( (D1*D2 == 1) ) {
      // Both 1-D Arrays --> Element-wise multiplication
      SArray<T,D0> ret;
      for (ulong i=0; i<D0; i++) {
        ret.data[i] = data[i] * rhs.data[i];
      }
      return ret;
    } else {
      // Matrix-Vector multiplication
      SArray<T,D1> ret;
      for (ulong j=0; j<D1; j++) {
        T tot = 0;
        for (ulong i=0; i<D0; i++) {
          tot += (*this)(i,j) * rhs(i);
        }
        ret(j) = tot;
      }
      return ret;
    }
  }

  template <class I, ulong E0> inline SArray<T,E0,D1> operator*(SArray<I,E0,D0> const &rhs) {
    //This template matches Matrix-Matrix multiplication
    SArray<T,E0,D1> ret;
    for (ulong j=0; j<E0; j++) {
      for (ulong i=0; i<D1; i++) {
        T tot = 0;
        for (ulong k=0; k<D0; k++) {
          tot += (*this)(k,i) * rhs(j,k);
        }
        ret(j,i) = tot;
      }
    }
    return ret;
  }

  inline void operator=(double rhs) {
    //Scalar assignment
    for (ulong i=0; i<D0*D1*D2; i++) {
      data[i] = rhs;
    }
  }

  inline T sum() {
    //Scalar division
    T sum = 0.;
    for (ulong i=0; i<D0*D1*D2; i++) {
      sum += data[i];
    }
    return sum;
  }

  inline void operator/=(double rhs) {
    //Scalar division
    for (ulong i=0; i<D0*D1*D2; i++) {
      data[i] = data[i] / rhs;
    }
  }

  inline SArray<T,D0,D1,D2> operator*(double rhs) {
    //Scalar multiplication
    SArray<T,D0,D1,D2> ret;
    for (ulong i=0; i<D0*D1*D2; i++) {
      ret.data[i] = data[i] * rhs;
    }
    return ret;
  }

  inline SArray<T,D0,D1,D2> operator/(double rhs) {
    //Scalar division
    SArray<T,D0,D1,D2> ret;
    for (ulong i=0; i<D0*D1*D2; i++) {
      ret.data[i] = data[i] / rhs;
    }
    return ret;
  }

  inline friend std::ostream &operator<<(std::ostream& os, SArray const &v) {
    if (D1*D2 == 1) {
      for (ulong i=0; i<D0; i++) {
        os << std::setw(12) << v(i) << "\n";
      }
    } else if (D2 == 1) {
      for (ulong j=0; j<D1; j++) {
        for (ulong i=0; i<D0; i++) {
          os << std::setw(12) << v(i,j) << " ";
        }
        os << "\n";
      }
    }
    return os;
  }

};

#endif

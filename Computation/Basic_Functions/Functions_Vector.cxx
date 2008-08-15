// Copyright (C) 2001-2008 Vivien Mallet
// File authors: Vivien Mallet (main part), Marc Duruflé.
//
// This file is part of Seldon library.
// Seldon library provides matrices and vectors structures for
// linear algebra.
//
// Seldon is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Seldon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the Seldon home page:
//     http://spacetown.free.fr/lib/seldon/

#ifndef SELDON_FILE_FUNCTIONS_VECTOR_CXX

/*
  Functions defined in this file:

  alpha.X -> X
  Mlt(alpha, X)

  alpha.X + Y -> Y
  Add(alpha, X, Y)
  
  X.Y
  DotProd(X, Y)
  DotProdConj(X, Y)
  
  Omega*X
  GenRot(x, y, cos, sin)
  ApplyRot(x, y, cos, sin)
  
*/

namespace Seldon
{


  /////////
  // Mlt //


  template <class T0,
	    class T1, class Storage1, class Allocator1>
  void Mlt(const T0 alpha,
	   Vector<T1, Storage1, Allocator1>& X)  throw()
  {
    T1 alpha_ = alpha;

    typename Vector<T1, Storage1, Allocator1>::pointer data = X.GetData();

    for (int i = 0; i < X.GetDataSize(); i++)
      data[i] = alpha_ * data[i];
  }


  // Mlt //
  /////////



  /////////
  // Add //


  template <class T0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Add(const T0 alpha,
	   const Vector<T1, Storage1, Allocator1>& X,
	   Vector<T2, Storage2, Allocator2>& Y)  throw(WrongDim, NoMemory)
  {
    if (alpha != T0(0))
      {
	T1 alpha_ = alpha;

	int ma = X.GetM();
	
#ifdef SELDON_CHECK_BOUNDARIES
	CheckDim(X, Y, "Add(alpha, X, Y)");
#endif

	for (int i = 0; i < ma; i++)
	  Y(i) += alpha_ * X(i);
      }
  }


  // Add //
  /////////
  
  
  
  /////////////
  // DotProd //
  
  
  //! Scalar product between two vectors
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  T1 DotProd(const Vector<T1, Storage1, Allocator1>& X,
	     const Vector<T2, Storage2, Allocator2>& Y)
  {
    T1 value(0);

#ifdef SELDON_CHECK_BOUNDARIES
    CheckDim(X, Y, "DotProd(X, Y)");
#endif

    for (int i = 0; i < X.GetM(); i++)
      value += X(i)*Y(i);
    
    return value;
  }
  
  
  //! Scalar product between two vectors
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  T1 DotProdConj(const Vector<T1, Storage1, Allocator1>& X,
		 const Vector<T2, Storage2, Allocator2>& Y)
  {
    T1 value(0);

#ifdef SELDON_CHECK_BOUNDARIES
    CheckDim(X, Y, "DotProdConj(X, Y)");
#endif

    for (int i = 0; i < X.GetM(); i++)
      value += X(i)*Y(i);
    
    return value;
  }
  
  
  //! Scalar product between two vectors
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  complex<T1> DotProdConj(const Vector<complex<T1>, Storage1, Allocator1>& X,
			  const Vector<T2, Storage2, Allocator2>& Y)
  {
    complex<T1> value(0);

#ifdef SELDON_CHECK_BOUNDARIES
    CheckDim(X, Y, "DotProdConj(X, Y)");
#endif

    for (int i = 0; i < X.GetM(); i++)
      value += conj(X(i))*Y(i);
    
    return value;
  }
  
  
  // DotProd //
  /////////////
  
  
  
  //////////////
  // ApplyRot //
  
  
  //! Computation of rotation between two points
  template<class T>
  void GenRot(T& a_in, T& b_in, T& c_, T& s_)
  {
    // old BLAS version
    T roe;
    if (abs(a_in) > abs(b_in))
      roe = a_in;
    else
      roe = b_in;
      
    T scal = abs(a_in) + abs(b_in);
    T r, z;
    if (scal != T(0))
      {
	T a_scl = a_in / scal;
	T b_scl = b_in / scal;
	r = scal * sqrt(a_scl * a_scl + b_scl * b_scl);
	if (roe < T(0))
	  r *= -1;
	
	c_ = a_in / r;
	s_ = b_in / r;
	z = 1;
	if (abs(a_in) > abs(b_in))
	  z = s_;
	else if ((abs(b_in) >= abs(a_in)) && (c_ != T(0)))
	  z = T(1) / c_;
      }
    else
      {
	c_ = 1;	s_ = 0; r = 0; z = 0;
      }
    a_in = r;
    b_in = z;
  }
  
  
  //! Computation of rotation between two points
  template<class T>
  void GenRot(complex<T>& a_in, complex<T>& b_in, T& c_, complex<T>& s_)
  {
      
    T a = abs(a_in), b = abs(b_in);
    if ( a == T(0) )
      {
	c_ = T(0);
	s_ = complex<T>(1,0);
	a_in = b_in;
      }
    else
      {
	T scale = a + b;
	T a_scal = abs(a_in/scale);
	T b_scal = abs(b_in/scale);
	T norm = sqrt(a_scal*a_scal+b_scal*b_scal) * scale;
	
	c_ = a / norm;
	complex<T> alpha = a_in/a;
	s_ =  alpha * conj(b_in)/norm;
	a_in = alpha * norm;
      }
    b_in = complex<T>(0,0);
  }
  
  
  //! Rotation of a point in 2-D
  template<class T>
  void ApplyRot(T& x, T& y, const T c_, const T s_)
  {
    T temp = c_*x + s_*y;
    y = c_*y - s_*x;
    x = temp;
  }
  
  
  //! Rotation of a complex point in 2-D
  template<class T>
  void ApplyRot(complex<T>& x, complex<T>& y,
		const T& c_, const complex<T>& s_)
  {
    complex<T> temp = s_*y + c_*x;
    y = -conj(s_)*x + c_*y;
    x = temp;
  }
  
  
  // ApplyRot //
  //////////////
  
  
  
  //////////////
  // CHECKDIM //


  //! Checks the compatibility of the dimensions.
  /*! Checks that X + Y is possible according to the dimensions of
    the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "X + Y".
  */
  template <class T0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void CheckDim(const Vector<T0, Storage0, Allocator0>& X,
		const Vector<T1, Storage1, Allocator1>& Y,
		string function = "", string op = "X + Y")
  {
    if (X.GetLength() != Y.GetLength())
      throw WrongDim(function, string("Operation ") + op
		     + string(" not permitted:")
		     + string("\n     X (") + to_str(&X) + string(") is a ")
		     + string("vector of length ") + to_str(X.GetLength())
		     + string(";\n     Y (") + to_str(&Y) + string(") is a ")
		     + string("vector of length ") + to_str(Y.GetLength())
		     + string("."));
  }


  // CHECKDIM //
  //////////////


  
  ///////////////
  // CONJUGATE //
  
  
  //! Sets a vector to its conjugate.
  template<class T, class Prop, class Allocator>
  void Conjugate(Vector<T, Prop, Allocator>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) = conj(X(i));
  }
  
  
  // CONJUGATE //
  ///////////////
  
  
} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_VECTOR_CXX
#endif

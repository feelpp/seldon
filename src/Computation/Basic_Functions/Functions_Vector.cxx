// Copyright (C) 2001-2005 Vivien Mallet
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
	CheckDim(X, Y, "Mlt(X, Y)");
#endif

	for (int i = 0; i < ma; i++)
	  Y(i) += alpha_ * X(i);
      }
  }


  // Add //
  /////////


      
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


} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_VECTOR_CXX
#endif

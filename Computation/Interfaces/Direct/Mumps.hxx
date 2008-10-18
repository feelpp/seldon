// Copyright (C) 2001-2008 Vivien Mallet
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

#ifndef SELDON_FILE_MUMPS_HXX

// including Mumps headers
extern "C"
{
#include "dmumps_c.h"
#include "zmumps_c.h"
#include "mpi.h"
}

namespace Seldon
{
  template<class T>
  class TypeMumps
  {
  };
 
  
  //! class containing MUMPS data structure
  template<>
  class TypeMumps<double>
  {
  public :
    typedef DMUMPS_STRUC_C data;
    typedef double* pointer;
  };
  
  
  //! class containing MUMPS data structure
  template<>
  class TypeMumps<complex<double> >
  {
  public :
    typedef ZMUMPS_STRUC_C data;
    typedef mumps_double_complex* pointer;
  };
  
  
  //! object used to solve linear system by calling mumps subroutines
  template<class T>
  class MatrixMumps
  {
  protected :
    int rank; //!< rank of processor
    int type_ordering; //!< ordering scheme (AMD, Metis, etc)
    //! object containing Mumps data structure
    typename TypeMumps<T>::data struct_mumps;
    //! double* or complex<double>*
    typedef typename TypeMumps<T>::pointer pointer;
    int print_level;
    
    // internal method
    void CallMumps();
    
  public :
    MatrixMumps();
    ~MatrixMumps();
    
    void Clear();
    
    void InitSymmetricMatrix();
    void InitUnSymmetricMatrix();
    void SelectOrdering(int num_ordering);
    void HideMessages();
    void ShowMessages();
    int GetInfoFactorization();
    
    template<class Prop,class Storage,class Allocator>
    void FindOrdering(Matrix<T, Prop, Storage, Allocator> & mat,
		      IVect& numbers, bool keep_matrix = false);
    
    template<class Prop,class Storage,class Allocator>
    void FactorizeMatrix(Matrix<T,Prop,Storage,Allocator> & mat,
			 bool keep_matrix = false);
    
    template<class Prop1, class Storage1, class Allocator1,
	     class Prop2, class Storage2, class Allocator2>
    void GetSchurMatrix(Matrix<T,Prop1,Storage1,Allocator1>& mat,
			const IVect& num,
			Matrix<T,Prop2,Storage2,Allocator2> & mat_schur,
			bool keep_matrix = false);
    
    template<class Allocator2>
    void Solve(Vector<T,Vect_Full,Allocator2>& x);
    
  };
  
}

#define SELDON_FILE_MUMPS_HXX
#endif




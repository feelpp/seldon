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

#ifndef SELDON_FILE_SUPERLU_HXX

extern "C"
{
#include "Computation/Interfaces/Direct/superlu_interface.h"
}

namespace Seldon
{
  //! class interfacing SuperLU functions
  template<class T>
  class MatrixSuperLU_Base
  {
  public :
    //! objects of SuperLU
    SuperMatrix A, L, U, B;
    SCformat *Lstore;  //! object of SuperLU
    NCformat *Ustore;  //! object of SuperLU
    //! permutation array
    IVect perm_r, perm_c;
    
    int permc_spec; //!< ordering scheme
    int n; //!< number of rows
    bool display_info; //!< display information about factorization ?
    
  public :
    MatrixSuperLU_Base();
    ~MatrixSuperLU_Base();
    
    void Clear();
    void HideMessages();
    void ShowMessages();
    
  };
  
  //! empty matrix
  template<class T>
  class MatrixSuperLU : public MatrixSuperLU_Base<T>
  {
  };
  
  //! class interfacing SuperLU functions in double precision
  template<>
  class MatrixSuperLU<double> : public MatrixSuperLU_Base<double>
  {
  public:
    MatrixSuperLU() : MatrixSuperLU_Base<double>() {}
    ~MatrixSuperLU() {}
    
    
    template<class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<double,Prop,Storage,Allocator> & mat,
			 bool keep_matrix = false);
    
    template<class Allocator2>
    void Solve(Vector<double,Vect_Full,Allocator2>& x);
    
  };
  
  //! class interfacing SuperLU functions in complex double precision
  template<>
  class MatrixSuperLU<complex<double> >
    : public MatrixSuperLU_Base<complex<double> >
  {
  public:
    MatrixSuperLU() : MatrixSuperLU_Base<complex<double> >() {}
    ~MatrixSuperLU() {}
    
    template<class Prop, class Storage,class Allocator>
    void FactorizeMatrix(Matrix<complex<double>,Prop,Storage,Allocator> & mat,
			 bool keep_matrix = false);
    
    template<class Allocator2>
    void Solve(Vector<complex<double>,Vect_Full,Allocator2>& x);
  };
}

#define SELDON_FILE_SUPERLU_HXX
#endif

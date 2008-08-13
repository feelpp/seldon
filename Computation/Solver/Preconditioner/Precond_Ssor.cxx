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

#ifndef SELDON_FILE_PRECOND_SSOR_CXX

namespace Seldon
{
  template<class T, class Prop, class Storage, class Allocator>
  class SsorPreconditioner
  {
  protected :
    bool symmetric_precond; //!< true for Symmetric relaxation
    int nb_iter; //!< number of iterations
    T omega; //!< relaxation parameter
    Matrix<T, Prop, Storage, Allocator>* A_; //!< pointer to matrix
    
  public :
    SsorPreconditioner();
    ~SsorPreconditioner(){}
    
    void Init(Matrix<T, Prop, Storage, Allocator>& B);
    void InitSymmetricPreconditioning() { symmetric_precond = true; }
    void InitUnSymmetricPreconditioning() { symmetric_precond = false; }
    void SetParameterRelaxation(const T& param) { omega = param; }
    void SetNumberIterations(int nb_iterations) { nb_iter = nb_iterations; }
    
    template<class Vector, class Matrix>
    void Solve(const Matrix& A, const Vector& r, Vector& z);
    template<class Vector, class Matrix>
    void TransSolve(const Matrix& A, const Vector& r, Vector& z);
    
  };
  
  
  //! Default constructor
  template<class T, class Prop, class Storage, class Allocator>
  SsorPreconditioner<T, Prop, Storage, Allocator>::SsorPreconditioner()
  {
    nb_iter = 1; omega = T(1);
    symmetric_precond = true;
  }
  
  
  //! Initialization with a given matrix
  template<class T, class Prop, class Storage, class Allocator>
  void SsorPreconditioner<T, Prop, Storage, Allocator>::
  Init(Matrix<T, Prop, Storage, Allocator>& B)
  {
    A_ = &B; 
  }
  
  
  //! Solves M z = r
  template<class T, class Prop, class Storage, class Allocator>
  template<class Vector, class Matrix>
  void SsorPreconditioner<T, Prop, Storage, Allocator>::
  Solve(const Matrix& A, const Vector& r, Vector& z)
  {    
    z.Zero();
    if (symmetric_precond)
      Seldon::SOR(*A_, z, r, omega, nb_iter, 0);
    else
      Seldon::SOR(*A_, z, r, omega, nb_iter, 2);
  }
  
  
  //! Solves M^t z = r
  template<class T, class Prop, class Storage, class Allocator>
  template<class Vector, class Matrix>
  void SsorPreconditioner<T, Prop, Storage, Allocator>::
  TransSolve(const Matrix& A, const Vector& r, Vector& z)
  {    
    z.Zero();
    if (symmetric_precond)
      Seldon::SOR(*A_, z, r, omega, nb_iter, 0);
    else
      Seldon::SOR(*A_, z, r, omega, nb_iter, 3);
  }
  
}

#define SELDON_FILE_PRECOND_SSOR_CXX
#endif

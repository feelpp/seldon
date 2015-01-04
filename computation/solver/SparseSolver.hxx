// Copyright (C) 2010 Vivien Mallet
// Copyright (C) 2010 Marc Durufle
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_FILE_COMPUTATION_SPARSESOLVER_HXX

#ifdef SELDON_WITH_PRECONDITIONING
#include "computation/solver/preconditioner/IlutPreconditioning.hxx"
#endif

#include "Ordering.hxx"

namespace Seldon
{
  
  //! Default solver in Seldon
  template<class T, class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class SparseSeldonSolver
  {
  protected :
    //! Verbosity level
    int print_level;
    //! Threshold for pivoting
    typename ClassComplexType<T>::Treal permtol;
    //! Symmetric matrix.
    Matrix<T, Symmetric, ArrayRowSymSparse, Allocator> mat_sym;
    //! Unsymmetric matrix.
    Matrix<T, General, ArrayRowSparse, Allocator> mat_unsym;
    //! Permutation arrays.
    IVect permutation_row, permutation_col;
    //! if true the factorisation is contained in mat_sym
    bool symmetric_matrix;
    
  public :
    
    SparseSeldonSolver();
    
    void Clear();
    
    void HideMessages();
    void ShowMessages();
    
    int64_t GetMemorySize() const;
    
    double GetPivotThreshold() const;
    void SetPivotThreshold(const double&);
    
    template<class T0, class Storage0, class Allocator0>
    void FactorizeMatrix(const IVect& perm,
                         Matrix<T0, General, Storage0, Allocator0>& mat,
                         bool keep_matrix = false);

    template<class T0, class Storage0, class Allocator0>
    void FactorizeMatrix(const IVect& perm,
                         Matrix<T0, Symmetric, Storage0, Allocator0>& mat,
                         bool keep_matrix = false);

    template<class Vector1>
    void Solve(Vector1& z);
    
    template<class TransStatus, class Vector1>
    void Solve(const TransStatus& TransA, Vector1& z);
    
  };

  template<class T, class Treal, class Allocator>
  void GetLU(Matrix<T, General, ArrayRowSparse, Allocator>& A,
	     IVect& iperm, IVect& rperm, 
	     const Treal& permtol, int print_level);

  template<class real, class cplx,
	   class Allocator, class Storage2, class Allocator2>
  void SolveLU(const Matrix<real, General, ArrayRowSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x);

  template<class real, class cplx, class TransStatus,
	   class Allocator, class Storage2, class Allocator2>
  void SolveLU(const TransStatus& transA,
               const Matrix<real, General, ArrayRowSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x);

  template<class T, class Allocator>
  void GetLU(Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A, int print_level);

  template<class real, class cplx, class Allocator,
           class Storage2, class Allocator2>
  void SolveLU(const Matrix<real, Symmetric, ArrayRowSymSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x);

  template<class T0, class Prop, class Storage, class Allocator,
	   class T, class Alloc2>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A,
	     SparseSeldonSolver<T, Alloc2>& mat_lu,
	     bool keep_matrix = false);

  template<class T0, class Prop, class Storage, class Allocator,
	   class T, class Alloc2>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A,
	     SparseSeldonSolver<T, Alloc2>& mat_lu,
	     IVect& permut, bool keep_matrix = false);

  template<class T, class Alloc2, class T1, class Allocator>
  void SolveLU(SparseSeldonSolver<T, Alloc2>& mat_lu,
	       Vector<T1, VectFull, Allocator>& x);

  template<class T, class Alloc2, class T1, class Allocator, class Transpose_status>
  void SolveLU(const Transpose_status& TransA,
	       SparseSeldonSolver<T, Alloc2>& mat_lu,
	       Vector<T1, VectFull, Allocator>& x);
    
}  // namespace Seldon.


#define SELDON_FILE_COMPUTATION_SPARSESOLVER_HXX
#endif

// Copyright (C) 2010 Marc Duruflé
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

#ifndef SELDON_FILE_SPARSE_CHOLESKY_FACTORISATION_HXX

namespace Seldon
{
  
  //! Class grouping different Cholesky solvers
  template<class T>
  class SparseCholeskySolver
  {
  protected :
    //! Verbose level
    int print_level;
    //! ordering to use
    int type_ordering;
    //! permutation array
    IVect permutation;
    //! solver to use
    int type_solver;
    //! size of factorized linear system 
    int n;
    //! Cholesky factors
    Matrix<T, Symmetric, ArrayRowSymSparse> mat_sym;
    //! temporary vector
    Vector<T> xtmp;
    
#ifdef SELDON_WITH_CHOLMOD
    MatrixCholmod mat_chol;
#endif
        
  public :
    // available solvers
    enum {SELDON_SOLVER, CHOLMOD};
    
    SparseCholeskySolver();
    
    void HideMessages();
    void ShowMessages();
    void ShowFullHistory();
    
    void Clear();
    
    int GetM() const;
    int GetN() const;
    
    int GetTypeOrdering() const;
    void SetOrdering(const IVect&);
    void SetTypeOrdering(int);
    
    void SelectDirectSolver(int);
    int GetDirectSolver();
    
    template<class MatrixSparse>
    void Factorize(MatrixSparse& A, bool keep_matrix = false);
    
    template<class TransStatus, class Vector1>
    void Solve(const TransStatus& TransA, Vector1& x);

    template<class TransStatus, class Vector1>
    void Mlt(const TransStatus& TransA, Vector1& x);
    
  };


  template<class T, class Prop, class Allocator>
  void GetCholesky(Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
                   int print_level = 0);

  template<class classTrans, 
           class T0, class T1, class Prop, class Storage,
           class Allocator1, class Allocator2>
  void SolveCholesky(const classTrans& TransA,
                     const Matrix<T0, Prop, ArrayRowSymSparse, Allocator1>& A,
                     Vector<T1, Storage, Allocator2>& x);

  template<class classTrans,
	   class T0, class T1, class Prop, class Allocator>
  void SolveCholesky(const classTrans& TransA,
		     const Matrix<T0, Prop, RowSymSparse>& A,
                     Vector<T1, VectFull, Allocator>& X);

  template<class classTrans, 
           class T0, class T1, class Prop, class Storage,
           class Allocator1, class Allocator2>
  void MltCholesky(const classTrans& TransA,
                   const Matrix<T0, Prop, ArrayRowSymSparse, Allocator1>& A,
                   Vector<T1, Storage, Allocator2>& x);

  template<class classTrans, 
           class T0, class T1, class Prop, class Storage,
           class Allocator1, class Allocator2>
  void MltCholesky(const classTrans& TransA,
                   const Matrix<T0, Prop, RowSymSparse, Allocator1>& A,
                   Vector<T1, Storage, Allocator2>& x);
  
}

#define SELDON_FILE_SPARSE_CHOLESKY_FACTORISATION_HXX
#endif

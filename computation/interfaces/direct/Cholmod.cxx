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


#ifndef SELDON_FILE_CHOLMOD_CXX

#include "Cholmod.hxx"

namespace Seldon
{

  MatrixCholmod::MatrixCholmod()
  {
    // Sets default parameters.
    cholmod_start(&param_chol);

    n = 0;
  }



  MatrixCholmod::~MatrixCholmod()
  {
    if (n > 0)
      {
        Clear();
        cholmod_finish(&param_chol);
      }
  }


  void MatrixCholmod::Clear()
  {
    if (n > 0)
      {
        n = 0;
        cholmod_free_factor(&L, &param_chol);
      }
  }


  template<class Prop, class Storage, class Allocator>
  void MatrixCholmod::
  FactorizeMatrix(Matrix<double, Prop, Storage, Allocator> & mat,
                  bool keep_matrix)
  {
    Clear();

    n = mat.GetM();
    Matrix<double, Symmetric, RowSymSparse, MallocAlloc<double> > Acsc;
    Copy(mat, Acsc);
    if (!keep_matrix)
      mat.Clear();

    // Initialization of sparse matrix.
    cholmod_sparse A;

    A.nrow = n;
    A.ncol = n;
    A.nzmax = Acsc.GetDataSize();
    A.nz = NULL;
    A.p = Acsc.GetPtr();
    A.i = Acsc.GetInd();
    A.x = Acsc.GetData();
    A.z = NULL;
    A.stype = -1;
    A.xtype = CHOLMOD_REAL;
    A.dtype = CHOLMOD_DOUBLE;
    A.sorted = true;
    A.packed = true;

    L = cholmod_analyze(&A, &param_chol);

    cholmod_factorize(&A, L, &param_chol);

    cholmod_change_factor(CHOLMOD_REAL, true, L->is_super,
                          true, L->is_monotonic, L, &param_chol);
  }


  template<class Transpose_status, class Allocator>
  void MatrixCholmod::Solve(const Transpose_status& TransA,
                            Vector<double, VectFull, Allocator>& x)
  {
    // and dense right hand side
    cholmod_dense b_rhs;
    b_rhs.nrow = x.GetM();
    b_rhs.ncol = 1;
    b_rhs.nzmax = b_rhs.nrow;
    b_rhs.d = b_rhs.nrow;
    b_rhs.x = x.GetData();
    b_rhs.z = NULL;
    b_rhs.xtype = CHOLMOD_REAL;
    b_rhs.dtype = CHOLMOD_DOUBLE;

    cholmod_dense* x_sol, *y;
    if (TransA.Trans())
      {
        y = cholmod_solve(CHOLMOD_Lt, L, &b_rhs, &param_chol);
        x_sol = cholmod_solve(CHOLMOD_Pt, L, y, &param_chol);
      }
    else
      {
        y = cholmod_solve(CHOLMOD_P, L, &b_rhs, &param_chol);
        x_sol = cholmod_solve(CHOLMOD_L, L, y, &param_chol);
      }

    double* data = reinterpret_cast<double*>(x_sol->x);
    for (int i = 0; i < x.GetM(); i++)
      x(i) = data[i];

    cholmod_free_dense(&x_sol, &param_chol);
    cholmod_free_dense(&y, &param_chol);
  }


  template<class T, class Prop, class Storage, class Allocator>
  void GetCholesky(Matrix<T, Prop, Storage, Allocator>& A,
                   MatrixCholmod& mat_chol, bool keep_matrix = false)
  {
    mat_chol.FactorizeMatrix(A, keep_matrix);
  }

  template<class T, class Allocator, class Transpose_status>
  void
  SolveCholesky(const Transpose_status& TransA,
                MatrixCholmod& mat_chol, Vector<T, VectFull, Allocator>& x)
  {
    mat_chol.Solve(TransA, x);
  }

}

#define SELDON_FILE_CHOLMOD_CXX
#endif

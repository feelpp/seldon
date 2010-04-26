// Copyright (C) 2010 Vivien Mallet
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


#ifndef SELDON_FILE_COMPUTATION_SPARSESOLVER_CXX


#include "SparseSolver.hxx"


namespace Seldon
{


  /////////////
  // SOLVELU //


  //! Solves a sparse linear system using LU factorization.
  /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
    \f$ X \f$ and \f$ Y \f$ are vectors.
    \param[in] M the sparse matrix of the linear system, to be factorized in
    LU form by UMFPACK, SuperLU or Mumps. On exit, \a M is cleared.
    \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
    solution \f$ X \f$ of the system.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                   Vector<T1, Storage1, Allocator1>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, Y, "SolveLU(M, Y)");
#endif

#if defined(SELDON_WITH_UMFPACK)
    MatrixUmfPack<T0> matrix_lu;
#elif defined(SELDON_WITH_SUPERLU)
    MatrixSuperLU<T0> matrix_lu;
#elif defined(SELDON_WITH_MUMPS)
    MatrixMumps<T0> matrix_lu;
#endif

#if !defined(SELDON_WITH_UMFPACK) && !defined(SELDON_WITH_SUPERLU)      \
  && !defined(SELDON_WITH_MUMPS)
    throw Undefined("SolveLU(M, Y)", "No sparse solver (UMFPACK, SuperLU "
                    "or Mumps) has been activated.");
#else
    GetLU(M, matrix_lu);
    SolveLU(matrix_lu, Y);
#endif
  }


  //! \copydoc SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M, Vector<T1, Storage1, Allocator1>& Y)
  template <class T, class Prop0, class Allocator0, class Allocator1>
  void Solve(Matrix<T, Prop0, ColSparse, Allocator0>& M,
             Vector<T, VectFull, Allocator1>& Y)
  {
    SparseSolve(M, Y);
  }


  //! \copydoc SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M, Vector<T1, Storage1, Allocator1>& Y)
  template <class T, class Prop0, class Allocator0, class Allocator1>
  void Solve(Matrix<T, Prop0, RowSparse, Allocator0>& M,
             Vector<T, VectFull, Allocator1>& Y)
  {
    SparseSolve(M, Y);
  }


  // SOLVELU //
  /////////////


}  // namespace Seldon.


#define SELDON_FILE_COMPUTATION_SPARSESOLVER_CXX
#endif

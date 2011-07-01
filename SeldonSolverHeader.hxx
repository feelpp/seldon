// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_SELDON_SOLVER_HEADER_HXX

// additional classes and functions for sparse matrices
#include "matrix_sparse/Matrix_Conversions.hxx"
#include "matrix_sparse/Matrix_ArraySparse.hxx"
#include "matrix_sparse/Permutation_ScalingMatrix.hxx"
#include "matrix_sparse/Relaxation_MatVect.hxx"
#include "matrix_sparse/Functions_MatrixArray.hxx"


// interfaces with direct solvers
#ifdef SELDON_WITH_MUMPS
#include "computation/interfaces/direct/Mumps.hxx"
#endif

#ifdef SELDON_WITH_UMFPACK
#include "computation/interfaces/direct/UmfPack.hxx"
#endif

#ifdef SELDON_WITH_SUPERLU
#include "computation/interfaces/direct/SuperLU.hxx"
#endif

#ifdef SELDON_WITH_PASTIX
#include "computation/interfaces/direct/Pastix.hxx"
#endif

#ifdef SELDON_WITH_PRECONDITIONING
#include "computation/solver/preconditioner/IlutPreconditioning.hxx"
#endif

#include "computation/solver/SparseSolver.hxx"

// iterative solvers and preconditioning
#include "computation/solver/iterative/Iterative.hxx"
#include "computation/solver/preconditioner/Precond_Ssor.cxx"

// Cholesky Solver
#ifdef SELDON_WITH_CHOLMOD
#include "computation/interfaces/direct/Cholmod.hxx"
#endif

#include "computation/solver/SparseCholeskyFactorisation.hxx"

// eigenvalue stuff
#ifdef SELDON_WITH_ARPACK
#include "computation/interfaces/eigenvalue/Arpack.hxx"
#include "computation/interfaces/eigenvalue/EigenvalueSolver.hxx"
#endif

#define SELDON_FILE_SELDON_SOLVER_HEADER_HXX
#endif

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

#ifndef SELDON_FILE_SELDON_SOLVER_HXX

// additional classes and functions for sparse matrices
#include "matrix_sparse/Matrix_Conversions.cxx"
#include "matrix_sparse/Matrix_ArraySparse.cxx"
#include "matrix_sparse/Matrix_ArrayComplexSparse.cxx"
#include "matrix_sparse/Permutation_ScalingMatrix.cxx"
#include "matrix_sparse/Relaxation_MatVect.cxx"
#include "matrix_sparse/Functions_MatrixArray.cxx"


// interfaces with direct solvers
#ifdef SELDON_WITH_MUMPS
#include "computation/interfaces/direct/Mumps.cxx"
#endif

#ifdef SELDON_WITH_UMFPACK
#include "computation/interfaces/direct/UmfPack.cxx"
#endif

#ifdef SELDON_WITH_SUPERLU
#include "computation/interfaces/direct/SuperLU.cxx"
#endif

// iterative solvers and preconditioning

#include "computation/solver/iterative/Iterative.cxx"
#include "computation/solver/preconditioner/Precond_Ssor.cxx"

#define SELDON_FILE_SELDON_SOLVER_HXX
#endif

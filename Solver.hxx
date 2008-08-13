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

#ifndef SELDON_FILE_SOLVER_HXX

// additional classes and functions for sparse matrices
#include "MatrixSparse/Functions_Arrays.cxx"
#include "MatrixSparse/Matrix_ArraySparse.cxx"
#include "MatrixSparse/Matrix_ArrayComplexSparse.cxx"
#include "MatrixSparse/Matrix_Conversions.cxx"
#include "MatrixSparse/Permutation_ScalingMatrix.cxx"
#include "MatrixSparse/Relaxation_MatVect.cxx"
#include "MatrixSparse/Functions_MatrixArray.cxx"


// interfaces with direct solvers
#ifdef SELDON_WITH_MUMPS
#include "Computation/Interfaces/Direct/Mumps.cxx"
#endif

#ifdef SELDON_WITH_UMFPACK
#include "Computation/Interfaces/Direct/UmfPack.cxx"
#endif

#ifdef SELDON_WITH_SUPERLU
#include "Computation/Interfaces/Direct/SuperLU.cxx"
#endif

// iterative solvers and preconditioning

#include "Computation/Solver/Iterative/Iterative.cxx"
#include "Computation/Solver/Preconditioner/Precond_Ssor.cxx"

#define SELDON_FILE_SOLVER_HXX
#endif

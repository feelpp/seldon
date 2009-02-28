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

#ifndef SELDON_FILE_MUMPS_CXX

#include "Mumps.hxx"

namespace Seldon
{
  
  //! Mumps is called in double precision
  template<>
  inline void MatrixMumps<double>::CallMumps()
  {
    dmumps_c(&struct_mumps);
  }
  
  
  //! Mumps is called in complex double precision
  template<>
  inline void MatrixMumps<complex<double> >::CallMumps()
  {
    zmumps_c(&struct_mumps);
  }
  
  
  //! initialization
  template<class T>
  inline MatrixMumps<T>::MatrixMumps()
  {
    // MPI initialization for sequential execution
    int ierr = MPI_Init(0, 0);
    ierr = MPI_Comm_rank(-987654, &rank);
    
    // parameters for mumps
    struct_mumps.job = -1;
    struct_mumps.par = 1;
    struct_mumps.sym = 0; // 0 -> unsymmetric matrix
    struct_mumps.comm_fortran = -987654;
    
    // mumps is called
    CallMumps();
    
    // other parameters
    struct_mumps.n = 0;
    struct_mumps.icntl[6] = 7;
    type_ordering = -1; // default : we let Mumps choose the ordering
    print_level = 0;
  }
  
  
  //! informs mumps that the matrix is symmetric
  template<class T>
  inline void MatrixMumps<T>::InitSymmetricMatrix()
  {
    // the user has to infom the symmetry during the initialization stage
    struct_mumps.job = -1;
    struct_mumps.par = 1;
    struct_mumps.sym = 2; // general symmetric matrix
    
    // mumps is called
    CallMumps();
    
    struct_mumps.icntl[13] = 20;
    
    // the print level is set in mumps
    if (print_level >= 0)
      {
	struct_mumps.icntl[0] = 6;
	struct_mumps.icntl[1] = 0;
	struct_mumps.icntl[2] = 6;
	struct_mumps.icntl[3] = 2;
      }
    else
      {
	struct_mumps.icntl[0] = -1;
	struct_mumps.icntl[1] = -1;
	struct_mumps.icntl[2] = -1;
	struct_mumps.icntl[3] = 0;
      }
  }
  
  
  //! informs mumps that the matrix is unsymmetric
  template<class T>
  inline void MatrixMumps<T>::InitUnSymmetricMatrix()
  {
    // the user has to infom the symmetry during the initialization stage
    struct_mumps.job = -1;
    struct_mumps.par = 1;
    struct_mumps.sym = 0; // general unsymmetric matrix
    
    // mumps is called
    CallMumps();
    
    struct_mumps.icntl[13] = 20;
    if (print_level >= 0)
      {
	struct_mumps.icntl[0] = 6;
	struct_mumps.icntl[1] = 0;
	struct_mumps.icntl[2] = 6;
	struct_mumps.icntl[3] = 2;
      }
    else
      {
	struct_mumps.icntl[0] = -1;
	struct_mumps.icntl[1] = -1;
	struct_mumps.icntl[2] = -1;
	struct_mumps.icntl[3] = 0;
      }
  }
  
  
  //! selects another ordering scheme
  template<class T>
  inline void MatrixMumps<T>::SelectOrdering(int num_ordering)
  {
    struct_mumps.icntl[6] = num_ordering;
  }
  
  
  //! clears factorization
  template<class T>
  MatrixMumps<T>::~MatrixMumps()
  {
    Clear();
  }
  
  
  //! clears factorization
  template<class T>
  inline void MatrixMumps<T>::Clear()
  {
    if (struct_mumps.n > 0)
      {
	struct_mumps.job = -2;
	CallMumps(); /* Terminate instance */
	MPI_Finalize();
	struct_mumps.n = 0;
      }
  }
  
  
  //! no display from Mumps
  template<class T>
  inline void MatrixMumps<T>::HideMessages()
  {
    print_level = -1;
    
    struct_mumps.icntl[0] = -1;
    struct_mumps.icntl[1] = -1;
    struct_mumps.icntl[2] = -1;
    struct_mumps.icntl[3] = 0;
    
  }
  
  
  //! standard display
  template<class T>
  inline void MatrixMumps<T>::ShowMessages()
  {
    print_level = 0;
    
    struct_mumps.icntl[0] = 6;
    struct_mumps.icntl[1] = 0;
    struct_mumps.icntl[2] = 6;
    struct_mumps.icntl[3] = 2;
    
  }
  
  
  //! computes row numbers
  /*!
    \param[in,out] mat matrix whose we want to find the ordering
    \param[out] numbers new row numbers
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class Prop,class Storage,class Allocator>
  void MatrixMumps<T>::FindOrdering(Matrix<T, Prop, Storage, Allocator> & mat,
				    IVect& numbers, bool keep_matrix)
  {
    int n = mat.GetM(), nnz = mat.GetNonZeros();
    // conversion in coordinate format
    IVect num_row, num_col; Vector<T, Vect_Full, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);
    // no values needed to renumber
    values.Clear();
    if (!keep_matrix)
      mat.Clear();
    
    /* Define the problem on the host */
    if (rank == 0)
      {
	struct_mumps.n = n; struct_mumps.nz = nnz;
	struct_mumps.irn = num_row.GetData();
	struct_mumps.jcn = num_col.GetData();
      }
    
    /* Call the MUMPS package. */
    struct_mumps.job = 1; // we analyse the system
    CallMumps();
    
    numbers.Reallocate(n);
    for (int i = 0; i < n; i++)
      numbers(i) = struct_mumps.sym_perm[i]-1;
  }
  
  
  //! factorization of a given matrix
  /*!
    \param[in,out] mat matrix to factorize
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class Prop, class Storage, class Allocator>
  void MatrixMumps<T>::FactorizeMatrix(Matrix<T,Prop,Storage,Allocator> & mat,
				       bool keep_matrix)
  {
    int n = mat.GetM(), nnz = mat.GetNonZeros();
    // conversion in coordinate format with fortran convention (1-index)
    IVect num_row, num_col; Vector<T, Vect_Full, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);
    if (!keep_matrix)
      mat.Clear();
    
    /* Define the problem on the host */
    if (rank == 0)
      {
	struct_mumps.n = n; struct_mumps.nz = nnz;
	struct_mumps.irn = num_row.GetData(); struct_mumps.jcn = num_col.GetData();
	struct_mumps.a = reinterpret_cast<pointer>(values.GetData());
      }
    
    /* Call the MUMPS package. */
    struct_mumps.job = 4; // we analyse and factorize the system
    CallMumps();
  }
  
  
  //! returns information about factorization performed
  template<class T>
  int MatrixMumps<T>::GetInfoFactorization() const
  {
    return struct_mumps.info[0];
  }
  
  
  //! Computation of Schur complement.
  /*!
    \param[in,out] mat initial matrix.
    \param[in] num numbers to keep in Schur complement.
    \param[out] mat_schur Schur matrix.
    \param[in] keep_matrix if false, \a mat is cleared.
  */
  template<class T> template<class Prop1, class Storage1, class Allocator,
			     class Prop2, class Storage2, class Allocator2>
  void MatrixMumps<T>::
  GetSchurMatrix(Matrix<T, Prop1, Storage1, Allocator>& mat, const IVect& num,
		 Matrix<T, Prop2, Storage2, Allocator2> & mat_schur,
		 bool keep_matrix)
  {
    int n_schur = num.GetM(), n = mat.GetM();
    // Subscripts are changed to respect fortran convention
    IVect index_schur(n_schur);
    for (int i = 0; i < n_schur; i++)
      index_schur(i) = num(i)+1;
    
    // array that will contain values of Schur matrix
    Vector<T, Vect_Full, Allocator2> vec_schur(n_schur*n_schur);
    
    struct_mumps.icntl[18] = n_schur;
    struct_mumps.size_schur = n_schur;
    struct_mumps.listvar_schur = index_schur.GetData();
    struct_mumps.schur = reinterpret_cast<pointer>(vec_schur.GetData());
    
    // factorization of the matrix
    FactorizeMatrix(mat, keep_matrix);
    
    // resetting parameters related to Schur complement
    struct_mumps.icntl[18] = 0;
    struct_mumps.size_schur = 0;
    struct_mumps.listvar_schur = NULL;
    struct_mumps.schur = NULL;
    
    // schur complement stored by rows
    int nb = 0;
    mat_schur.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n ;j++)
	mat_schur(i,j) = vec_schur(nb++);
    
    vec_schur.Clear(); index_schur.Clear();
  }
  
  
  //! resolution of a linear system using the computed factorization
  /*!
    \param[in,out] x right-hand-side on input, solution on output
    It is assumed that a call to FactorizeMatrix has been done before
  */
  template<class T> template<class Allocator2>
  void MatrixMumps<T>::Solve(Vector<T,Vect_Full,Allocator2>& x)
  {
    struct_mumps.rhs = reinterpret_cast<pointer>(x.GetData());
    struct_mumps.job = 3; // we solve system
    CallMumps();
  }
  
  
  template<class T, class Storage, class Allocator>
  void GetLU(Matrix<T,Symmetric,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.InitSymmetricMatrix();
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  
  
  template<class T, class Storage, class Allocator>
  void GetLU(Matrix<T,General,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.InitUnSymmetricMatrix();
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  
  
  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T,Symmetric,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
		      const IVect& num, MatrixFull& schur_matrix, bool keep_matrix = false)
  {
    mat_lu.InitSymmetricMatrix();
    mat_lu.GetSchurMatrix(A, num, schur_matrix, keep_matrix);
  }
  
  
  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T,General,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
		      const IVect& num, MatrixFull& schur_matrix, bool keep_matrix = false)
  {
    mat_lu.InitUnSymmetricMatrix();
    mat_lu.GetSchurMatrix(A, num, schur_matrix, keep_matrix);
  }
  
  
  template<class T, class Allocator>
  void SolveLU(MatrixMumps<T>& mat_lu, Vector<T, Vect_Full, Allocator>& x)
  {
    mat_lu.Solve(x);
  }
}

#define SELDON_FILE_MUMPS_CXX
#endif

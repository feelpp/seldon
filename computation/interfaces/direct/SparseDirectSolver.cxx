// Copyright (C) 2015 Marc Duruflé
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

#ifndef SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_CXX

#include "SparseDirectSolver.hxx"

namespace Seldon
{

  //! Default constructor
  template<class T>
  SparseDirectSolver<T>::SparseDirectSolver()
  {
    n = 0;
    type_ordering = SparseMatrixOrdering::AUTO;

    type_solver = SELDON_SOLVER;
    
    // we try to use an other available solver
    // The order of preference is Pastix, Mumps, Pardiso, UmfPack and SuperLU
#ifdef SELDON_WITH_SUPERLU
    type_solver = SUPERLU;
#endif
#ifdef SELDON_WITH_UMFPACK
    type_solver = UMFPACK;
#endif
#ifdef SELDON_WITH_PARDISO
    type_solver = PARDISO;
#endif
#ifdef SELDON_WITH_WSMP
    type_solver = WSMP;
#endif
#ifdef SELDON_WITH_MUMPS
    type_solver = MUMPS;
#endif
#ifdef SELDON_WITH_PASTIX
    type_solver = PASTIX;
#endif
    
    nb_threads_per_node = 1;
    threshold_matrix = 0;
    enforce_unsym_ilut = false;
  }
  
  
  //! hiding all messages
  template<class T>
  void SparseDirectSolver<T>::HideMessages()
  {
    mat_seldon.HideMessages();
    
#ifdef SELDON_WITH_MUMPS
    mat_mumps.HideMessages();
#endif

#ifdef SELDON_WITH_PARDISO
    mat_pardiso.HideMessages();
#endif

#ifdef SELDON_WITH_SUPERLU
    mat_superlu.HideMessages();
#endif

#ifdef SELDON_WITH_UMFPACK
    mat_umf.HideMessages();
#endif

#ifdef SELDON_WITH_PASTIX
    mat_pastix.HideMessages();
#endif

#ifdef SELDON_WITH_PRECONDITIONING
    mat_ilut.SetPrintLevel(0);
#endif

  }
  
  
  //! displaying basic messages
  template<class T>
  void SparseDirectSolver<T>::ShowMessages()
  {
    mat_seldon.ShowMessages();
    
#ifdef SELDON_WITH_MUMPS
    mat_mumps.ShowMessages();
#endif

#ifdef SELDON_WITH_PARDISO
    mat_pardiso.ShowMessages();
#endif

#ifdef SELDON_WITH_SUPERLU
    mat_superlu.ShowMessages();
#endif

#ifdef SELDON_WITH_UMFPACK
    mat_umf.ShowMessages();
#endif

#ifdef SELDON_WITH_PASTIX
    mat_pastix.ShowMessages();
#endif

#ifdef SELDON_WITH_PRECONDITIONING
    mat_ilut.SetPrintLevel(1);
#endif
    
  }
  
  
  //! displaying all the messages
  template<class T>
  void SparseDirectSolver<T>::ShowFullHistory()
  {
    mat_seldon.ShowMessages();
    
#ifdef SELDON_WITH_MUMPS
    mat_mumps.ShowMessages();
#endif

#ifdef SELDON_WITH_PARDISO
    mat_pardiso.ShowMessages();
#endif

#ifdef SELDON_WITH_SUPERLU
    mat_superlu.ShowMessages();
#endif

#ifdef SELDON_WITH_UMFPACK
    mat_umf.ShowMessages();
#endif

#ifdef SELDON_WITH_PASTIX
    mat_pastix.ShowFullHistory();
#endif

#ifdef SELDON_WITH_PRECONDITIONING
    mat_ilut.SetPrintLevel(1);
#endif
    
  }
    
  
  //! clearing factorisation
  template<class T>
  void SparseDirectSolver<T>::Clear()
  {
    if (n > 0)
      {
	n = 0;
	mat_seldon.Clear();
	
#ifdef SELDON_WITH_MUMPS
	mat_mumps.Clear();
#endif

#ifdef SELDON_WITH_PARDISO
	mat_pardiso.Clear();
#endif
	
#ifdef SELDON_WITH_SUPERLU
	mat_superlu.Clear();
#endif
	
#ifdef SELDON_WITH_UMFPACK
	mat_umf.Clear();
#endif
	
#ifdef SELDON_WITH_PASTIX
	mat_pastix.Clear();
#endif

#ifdef SELDON_WITH_WSMP
	mat_wsmp.Clear();
#endif

#ifdef SELDON_WITH_PRECONDITIONING
        mat_ilut.Clear();
#endif
    
      }    
  }
    
  template<class T> 
  bool SparseDirectSolver<T>::AffectOrdering()
  {
    bool user_ordering = false;
    // we set the ordering for each direct solver interfaced
    switch (type_ordering)
      {
      case SparseMatrixOrdering::AUTO :
	{
	  // we choose the default strategy
	  // proposed by the direct solver that will be called
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
	      mat_mumps.SelectOrdering(7);
#endif
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_PASTIX
	      mat_pastix.SelectOrdering(API_ORDER_SCOTCH);
#endif
	    }
	  else if (type_solver == PARDISO)
	    {
#ifdef SELDON_WITH_PARDISO
	      mat_pardiso.SelectOrdering(2);
#endif
	    }
	  else if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_UMFPACK
	      mat_umf.SelectOrdering(UMFPACK_ORDERING_AMD);
#endif
	    }
	  else if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      mat_superlu.SelectOrdering(superlu::COLAMD);
#endif 
	    }
	  else if (type_solver == WSMP)
	    {
	    }
	  else
	    {
              
              type_ordering = SparseMatrixOrdering::IDENTITY;
              
#ifdef SELDON_WITH_UMFPACK
              type_ordering = SparseMatrixOrdering::AMD;
#endif

#ifdef SELDON_WITH_MUMPS
              type_ordering = SparseMatrixOrdering::METIS;
#endif

#ifdef SELDON_WITH_PARDISO
              type_ordering = SparseMatrixOrdering::METIS;
#endif

#ifdef SELDON_WITH_PASTIX
              type_ordering = SparseMatrixOrdering::SCOTCH;
#endif
              
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::IDENTITY :
      case SparseMatrixOrdering::REVERSE_CUTHILL_MCKEE :
      case SparseMatrixOrdering::USER :
	{
	  user_ordering = true;
	}
	break;
      case SparseMatrixOrdering::PORD :
      case SparseMatrixOrdering::AMF :
      case SparseMatrixOrdering::QAMD :
	{
	  // Mumps orderings
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
	      if (type_ordering == SparseMatrixOrdering::PORD)
		mat_mumps.SelectOrdering(4);
	      else if (type_ordering == SparseMatrixOrdering::AMF)
		mat_mumps.SelectOrdering(2);
	      else
		mat_mumps.SelectOrdering(6);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::SCOTCH :
      case SparseMatrixOrdering::PTSCOTCH :
	{
	  // available for Mumps and Pastix
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
              if (type_ordering==SparseMatrixOrdering::PTSCOTCH)
                mat_mumps.SelectParallelOrdering(1);
              else
                mat_mumps.SelectOrdering(3);
#endif
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_SCOTCH
	      mat_pastix.SelectOrdering(API_ORDER_SCOTCH);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::METIS :
      case SparseMatrixOrdering::PARMETIS :
	{
	  // available for Mumps and Pastix
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
              if (type_ordering==SparseMatrixOrdering::PARMETIS)
                mat_mumps.SelectParallelOrdering(2);
              else
                mat_mumps.SelectOrdering(5);
#endif
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_PASTIX
	      mat_pastix.SelectOrdering(API_ORDER_METIS);
#endif
	    }
	  else if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_UMFPACK
	      mat_umf.SelectOrdering(UMFPACK_ORDERING_METIS);
#endif
	    }
          /*
            currently not implemented in SuperLU
	  else if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
              if (type_ordering==SparseMatrixOrdering::PARMETIS)
                mat_superlu.SelectOrdering(superlu::PARMETIS);
              else
                mat_superlu.SelectOrdering(superlu::METIS_AT_PLUS_A);
#endif
	    }
          */
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::AMD :
      case SparseMatrixOrdering::COLAMD :
	{
	  // default ordering for UmfPack
	  if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_UMFPACK
	      mat_umf.SelectOrdering(UMFPACK_ORDERING_AMD);
#endif
	    }
	  else if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
              mat_mumps.SelectOrdering(0);
#endif
            }
	  else if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      mat_superlu.SelectOrdering(superlu::COLAMD);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
      case SparseMatrixOrdering::MMD_AT_PLUS_A:
        {
          if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      mat_superlu.SelectOrdering(superlu::MMD_AT_PLUS_A);
#endif
	    }
        }
        break;
      case SparseMatrixOrdering::MMD_ATA:
        {
          if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      mat_superlu.SelectOrdering(superlu::MMD_ATA);
#endif
	    }
        }
	break;
      }
    
    return user_ordering;
  }
  
  //! computation of the permutation vector in order to reduce fill-in
  template<class T> template<class MatrixSparse>
  void SparseDirectSolver<T>::ComputeOrdering(MatrixSparse& A)
  {
    bool user_ordering = AffectOrdering();
    
    if (user_ordering)
      {
	// case where the ordering is not natively available in the direct solver
	// computing separetely the ordering
	FindSparseOrdering(A, permut, type_ordering);
        
        if (type_solver == MUMPS)
	  {
#ifdef SELDON_WITH_MUMPS
	    mat_mumps.SetPermutation(permut);
#endif
	  }
	else if (type_solver == PASTIX)
	  {
#ifdef SELDON_WITH_PASTIX
	    mat_pastix.SetPermutation(permut);
#endif
	  }
        else if (type_solver == PARDISO)
	  {
#ifdef SELDON_WITH_PARDISO
	    mat_pardiso.SetPermutation(permut);
#endif
	  }
	else if (type_solver == UMFPACK)
	  {
#ifdef SELDON_WITH_UMFPACK
	    mat_umf.SetPermutation(permut);
#endif
	  }
	else if (type_solver == SUPERLU)
	  {
#ifdef SELDON_WITH_SUPERLU
	    mat_superlu.SetPermutation(permut);
#endif 
	  }
	else
	  {
	  }	  
      }
	
  }
  
  
  //! factorisation of matrix A
  /*!
    LU factorisation is stored in the current object.
    You can ask to clear the matrix given on input (to spare memory)
   */
  template<class T> template<class MatrixSparse>
  void SparseDirectSolver<T>::Factorize(MatrixSparse& A, bool keep_matrix)
  {
    ComputeOrdering(A);
    
    n = A.GetM();
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
	GetLU(A, mat_umf, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with UmfPack support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
        mat_superlu.SetNumberOfThreadPerNode(nb_threads_per_node);
	GetLU(A, mat_superlu, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else if (type_solver == PARDISO)
      {
#ifdef SELDON_WITH_PARDISO
	GetLU(A, mat_pardiso, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Pardiso support.");
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	GetLU(A, mat_mumps, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
        mat_pastix.SetNumberOfThreadPerNode(nb_threads_per_node);
	GetLU(A, mat_pastix, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == WSMP)
      {
#ifdef SELDON_WITH_WSMP
        mat_wsmp.SetNumberOfThreadPerNode(nb_threads_per_node);
	GetLU(A, mat_wsmp, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Wsmp support.");
#endif
      }
    else if (type_solver == ILUT)
      {
#ifdef SELDON_WITH_PRECONDITIONING        
        // setting some parameters
        if (enforce_unsym_ilut || (!IsSymmetricMatrix(A)))
          mat_ilut.SetUnsymmetricAlgorithm();
        else
          mat_ilut.SetSymmetricAlgorithm();
        
        // then performing factorization
	GetLU(A, mat_ilut, permut, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with the preconditioners.");
#endif
      }
    else
      {
	GetLU(A, mat_seldon, permut, keep_matrix);
      }

  }
   

  //! Returns error code of the direct solver
  template <class T>
  int SparseDirectSolver<T>::GetInfoFactorization(int& ierr) const
  {
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
        ierr = mat_umf.GetInfoFactorization();
        switch (ierr)
          {
          case UMFPACK_OK :
            return FACTO_OK;
          case UMFPACK_WARNING_singular_matrix :
            return NUMERICALLY_SINGULAR_MATRIX;
          case UMFPACK_ERROR_out_of_memory :
            return OUT_OF_MEMORY;
          case UMFPACK_ERROR_invalid_Numeric_object :
          case UMFPACK_ERROR_invalid_Symbolic_object :
          case UMFPACK_ERROR_argument_missing :
          case UMFPACK_ERROR_different_pattern :
          case UMFPACK_ERROR_invalid_system :
            return INVALID_ARGUMENT;
          case UMFPACK_ERROR_n_nonpositive :
            return INCORRECT_NUMBER_OF_ROWS;
          case UMFPACK_ERROR_invalid_matrix :
            return MATRIX_INDICES_INCORRECT;
          case UMFPACK_ERROR_invalid_permutation :
            return INVALID_PERMUTATION;
          case UMFPACK_ERROR_ordering_failed :
            return ORDERING_FAILED;
          default :
            return INTERNAL_ERROR;
          }
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
        ierr = mat_superlu.GetInfoFactorization();        
        if (ierr > 0)
          return INTERNAL_ERROR;
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
        ierr = mat_mumps.GetInfoFactorization();
        switch (ierr)
          {
          case -2 :
            // nz out of range
            return MATRIX_INDICES_INCORRECT;
          case -3 :
            // invalid job number
            return INVALID_ARGUMENT;
          case -4 :
            // invalid permutation
            return INVALID_PERMUTATION;
          case -5 :
            // problem of real workspace allocation
            return OUT_OF_MEMORY;
          case -6 :
            // structurally singular matrix
            return STRUCTURALLY_SINGULAR_MATRIX;
          case -7 :
            // problem of integer workspace allocation
            return OUT_OF_MEMORY;
          case -10 :
            // numerically singular matrix
            return NUMERICALLY_SINGULAR_MATRIX;
          case -13 :
            // allocate failed
            return OUT_OF_MEMORY;
          case -16 :
            // N out of range
            return INCORRECT_NUMBER_OF_ROWS;
          case -22 :
            // invalid pointer
            return INVALID_ARGUMENT;
          case 1 :
            // index out of range
            return MATRIX_INDICES_INCORRECT;
	  case 0 :
	    return FACTO_OK;
          default :
            return INTERNAL_ERROR;
          }
#endif
      }
    else if (type_solver == PARDISO)
      {
#ifdef SELDON_WITH_PARDISO
        ierr = mat_pardiso.GetInfoFactorization();
        switch (ierr)
          {
          case -1 :
            return INVALID_ARGUMENT;
          case -2 :
          case -9 :
            return OUT_OF_MEMORY;
          case -3 :
            return INVALID_PERMUTATION;
          case -4 :
            return NUMERICALLY_SINGULAR_MATRIX;
          case -6 :
            return ORDERING_FAILED;
          case -7 :
            return STRUCTURALLY_SINGULAR_MATRIX;
          case -8 :
            return OVERFLOW_32BIT;            
	  case 0 :
	    return FACTO_OK;
          default :
            return INTERNAL_ERROR;
          }
#endif
      }
    
    return FACTO_OK;
  }
  
    
  //! x_solution is overwritten by solution of A x = b
  /*!
    We assume that Factorize has been called previously
  */
  template<class T> template<class Vector1>
  void SparseDirectSolver<T>::Solve(Vector1& x_solution)
  {
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
	Seldon::SolveLU(mat_umf, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with UmfPack support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
	Seldon::SolveLU(mat_superlu, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else if (type_solver == PARDISO)
      {
#ifdef SELDON_WITH_PARDISO
	Seldon::SolveLU(mat_pardiso, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with Pardiso support.");
#endif
      } 
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	Seldon::SolveLU(mat_mumps, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	Seldon::SolveLU(mat_pastix, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == WSMP)
      {
#ifdef SELDON_WITH_WSMP
	Seldon::SolveLU(mat_wsmp, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with Wsmp support.");
#endif
      }
    else if (type_solver == ILUT)
      {
#ifdef SELDON_WITH_PRECONDITIONING
	mat_ilut.Solve(x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with the preconditioners.");
#endif
      }
    else
      {
	Seldon::SolveLU(mat_seldon, x_solution);
      }
  }
  
  
  //! x_solution is overwritten with solution of A x = b or A^T x = b
  template<class T> template<class Vector1>
  void SparseDirectSolver<T>
  ::Solve(const SeldonTranspose& TransA, Vector1& x_solution)
  {
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
	Seldon::SolveLU(TransA, mat_umf, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with UmpfPack support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
	Seldon::SolveLU(TransA, mat_superlu, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      } 
    else if (type_solver == PARDISO)
      {
#ifdef SELDON_WITH_PARDISO
	Seldon::SolveLU(TransA, mat_pardiso, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with Pardiso support.");
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	Seldon::SolveLU(TransA, mat_mumps, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	Seldon::SolveLU(TransA, mat_pastix, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == WSMP)
      {
#ifdef SELDON_WITH_WSMP
	Seldon::SolveLU(TransA, mat_wsmp, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with Wsmp support.");
#endif
      }
    else if (type_solver == ILUT)
      {
#ifdef SELDON_WITH_PRECONDITIONING
	mat_ilut.Solve(TransA, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with the preconditioners.");
#endif
      }
    else
      {
	Seldon::SolveLU(TransA, mat_seldon, x_solution);
      }
  }
  

  //! x_solution is overwritten by solution of A x = b
  /*!
    Multiple right hand sides
    We assume that Factorize has been called previously
  */
  template<class T> template<class T1, class Alloc1>
  void SparseDirectSolver<T>
  ::Solve(Matrix<T1, General, ColMajor, Alloc1>& x_solution)
  {
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
        throw Undefined("SparseDirectSolver::Solve(Matrix&)",
                        "Multiple right hand sides not available in UmfPack.");
#else
        throw Undefined("SparseDirectSolver::Solve(Matrix&)",
                        "Seldon was not compiled with UmfPack support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
	throw Undefined("SparseDirectSolver::Solve(Matrix&)",
                        "Multiple right hand sides not available in SuperLU.");
#else
        throw Undefined("SparseDirectSolver::Solve(Matrix&)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else if (type_solver == PARDISO)
      {
#ifdef SELDON_WITH_PARDISO
        throw Undefined("SparseDirectSolver::Solve(Matrix&)",
                        "Multiple right hand sides not available in Pardiso.");
#else
        throw Undefined("SparseDirectSolver::Solve(Matrix&)",
                        "Seldon was not compiled with Pardiso support.");
#endif
      } 
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	Seldon::SolveLU(mat_mumps, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
        throw Undefined("SparseDirectSolver::Solve(Matrix&)",
	                 "Multiple right hand sides not available in Pastix.");
#else
        throw Undefined("SparseDirectSolver::Solve(Matrix&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == WSMP)
      {
#ifdef SELDON_WITH_WSMP
	Seldon::SolveLU(mat_wsmp, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with Wsmp support.");
#endif
      }
    else if (type_solver == ILUT)
      {
#ifdef SELDON_WITH_PRECONDITIONING
	throw Undefined("SparseDirectSolver::Solve(Matrix&)",
	                "Multiple right hand sides not available in ILUT.");
#else
        throw Undefined("SparseDirectSolver::Solve(Matrix&)",
                        "Seldon was not compiled with the preconditioners.");
#endif
      }
    else
      {
	throw Undefined("SparseDirectSolver::Solve(Matrix&)",
	                "Multiple right hand sides not available in Seldon.");	
      }
  }

  
#ifdef SELDON_WITH_MPI
  //! Factorization of a matrix
  /*!
    The matrix is given on each processor of the communicator in CSC
    format. If the matrix is assumed to be symmetric, you provide only
    the lower part of the matrix.
    \param[in] comm_facto communicator grouping processors involved in the
    factorisation.
    \param[inout] Ptr start indices
    \param[inout] Row column indices
    \param[inout] Val values of non-zero entries
    \param[in] glob_num global column numbers, each column of the global
    matrix is associated with one processor and only one
    \param[in] sym if true, the matrix is assumed to be symmetric
    \param[in] keep_matrix if false the input matrix is erased
   */
  template<class T> template<class Tint>
  void SparseDirectSolver<T>::
  FactorizeDistributed(MPI::Comm& comm_facto,
                       Vector<Tint>& Ptr, Vector<Tint>& Row, Vector<T>& Val,
                       const IVect& glob_num, bool sym,
                       bool reorder, bool keep_matrix)
  {
    bool user_ordering = AffectOrdering();
    if (user_ordering)
      {
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "user orderings not available in parallel");
      }
    
    n = Ptr.GetM()-1;
    if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS        
	mat_mumps.FactorizeDistributedMatrix(comm_facto, Ptr, Row,
					     Val, glob_num, sym, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
        mat_pastix.SetNumberOfThreadPerNode(nb_threads_per_node);
        mat_pastix.FactorizeDistributedMatrix(comm_facto, Ptr, Row,
					      Val, glob_num, sym, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == WSMP)
      {
#ifdef SELDON_WITH_WSMP
        mat_wsmp.SetNumberOfThreadPerNode(nb_threads_per_node);
        mat_wsmp.FactorizeDistributedMatrix(comm_facto, Ptr, Row,
                                            Val, glob_num, sym, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "Seldon was not compiled with Wsmp support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU_DIST
        mat_superlu.SetNumberOfThreadPerNode(nb_threads_per_node);
        mat_superlu.FactorizeDistributedMatrix(comm_facto, Ptr, Row,
                                               Val, glob_num, sym, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else
      {
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "The method is defined for Mumps and Pastix only.");
      }


  }
  
  
  //! solution of linear system Ax = b by using LU factorization
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T> template<class Vector1>
  void SparseDirectSolver<T>::
  SolveDistributed(MPI::Comm& comm_facto, Vector1& x_solution,
		   const IVect& glob_number)
  {
    if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	mat_mumps.SolveDistributed(comm_facto, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(MPI::Comm&,"
                        " Vector&, IVect&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
        mat_pastix.SolveDistributed(comm_facto, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(MPI::Comm&,"
                        " Vector&, IVect&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == WSMP)
      {
#ifdef SELDON_WITH_WSMP
        mat_wsmp.SolveDistributed(comm_facto, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(MPI::Comm&,"
                        " Vector&, IVect&)",
                        "Seldon was not compiled with Wsmp support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU_DIST
        mat_superlu.SolveDistributed(comm_facto, x_solution);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(MPI::Comm&,"
                        " Vector&, IVect&)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else
      {
        throw Undefined("SparseDirectSolver::SolveDistributed(MPI::Comm&,"
                        " Vector&, IVect&)",
                        "The method is defined for Mumps and Pastix only.");
      }

  }
  
  
  //! solution of linear system A^T x = b by using LU factorization (without scaling)
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T> template<class Vector1>
  void SparseDirectSolver<T>::
  SolveDistributed(MPI::Comm& comm_facto, const SeldonTranspose& TransA,
                   Vector1& x_solution, const IVect& glob_number)
  {
    if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	mat_mumps.SolveDistributed(comm_facto, TransA, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(TransStatus, "
                        "MPI::Comm&, Vector&, IVect&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	mat_pastix.SolveDistributed(comm_facto, TransA, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(TransStatus, "
                        "MPI::Comm&, Vector&, IVect&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU_DIST
	mat_superlu.SolveDistributed(comm_facto, TransA, x_solution);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(TransStatus, "
                        "MPI::Comm&, Vector&, IVect&)",
                        "Seldon was not compiled with SuperLU_DIST support.");
#endif
      }
    else if (type_solver == WSMP)
      {
#ifdef SELDON_WITH_WSMP
	mat_wsmp.SolveDistributed(comm_facto, TransA, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(TransStatus, "
                        "MPI::Comm&, Vector&, IVect&)",
                        "Seldon was not compiled with Wsmp support.");
#endif
      }
    else
      {
        throw Undefined("SparseDirectSolver::SolveDistributed(TransStatus, "
                        "MPI::Comm&, Vector&, IVect&)",
                        "The method is not defined for this solver.");
      }
  }
#endif
  
}

#define SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_CXX
#endif


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

#ifndef SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_HXX

namespace Seldon
{
  
  //! Class grouping different direct solvers
  template<class T>
  class SparseDirectSolver
  {
  protected :
    //! ordering to use
    int type_ordering;
    //! solver to use
    int type_solver;
    //! number of threads (for Pastix)
    int nb_threads_per_node;
    //! ordering (if supplied by the user)
    IVect permut;
    //! size of factorized linear system 
    int n;
    
#ifdef SELDON_WITH_UMFPACK
    MatrixUmfPack<T> mat_umf; //!< Umfpack solver
#endif
#ifdef SELDON_WITH_SUPERLU
    MatrixSuperLU<T> mat_superlu; //!< SuperLU solver
#endif
#ifdef SELDON_WITH_PARDISO
    MatrixPardiso<T> mat_pardiso; //!< Pardiso solver
#endif
#ifdef SELDON_WITH_MUMPS
    MatrixMumps<T> mat_mumps; //!< Mumps solver
#endif
#ifdef SELDON_WITH_PASTIX
    MatrixPastix<T> mat_pastix; //!< Pastix solver
#endif
#ifdef SELDON_WITH_WSMP
    MatrixWsmp<T> mat_wsmp; //!< Wsmp solver
#endif
    
#ifdef SELDON_WITH_PRECONDITIONING
    //! ILUT solver
    IlutPreconditioning<T> mat_ilut;
#endif
    
    //! threshold for ilut solver
    double threshold_matrix;
    //! use of non-symmetric ilut ?
    bool enforce_unsym_ilut;
    
    //! default solver
    SparseSeldonSolver<T> mat_seldon;
        
  public :
    // available solvers
    enum {SELDON_SOLVER, UMFPACK, SUPERLU, MUMPS, PASTIX, ILUT, PARDISO, WSMP};
    
    // error codes
    enum {FACTO_OK, STRUCTURALLY_SINGULAR_MATRIX,
          NUMERICALLY_SINGULAR_MATRIX, OUT_OF_MEMORY, INVALID_ARGUMENT,
          INCORRECT_NUMBER_OF_ROWS, MATRIX_INDICES_INCORRECT,
          INVALID_PERMUTATION, ORDERING_FAILED, INTERNAL_ERROR,
          OVERFLOW_32BIT};
    
    SparseDirectSolver();
    
    void HideMessages();
    void ShowMessages();
    void ShowFullHistory();
    
    void Clear();
    
    int GetM() const;
    int GetN() const;
    
    int GetTypeOrdering() const;
    void SetPermutation(const IVect&);
    bool AffectOrdering();
    void SelectOrdering(int);
    
    void SetNumberOfThreadPerNode(int m);
    int GetNumberOfThreadPerNode() const;

    void SetPivotThreshold(const double&);
    
    template<class MatrixSparse>
    void ComputeOrdering(MatrixSparse& A);
    
    void SelectDirectSolver(int);
    void SetNonSymmetricIlut();

    int GetDirectSolver();

    void RefineSolution();
    void DoNotRefineSolution();

    void SetCoefficientEstimationNeededMemory(double);
    void SetMaximumCoefficientEstimationNeededMemory(double);
    void SetIncreaseCoefficientEstimationNeededMemory(double);
    
    double GetThresholdMatrix() const;
    void SetThresholdMatrix(const double&);

    template<class MatrixSparse>
    void Factorize(MatrixSparse& A, bool keep_matrix = false);
    
    int GetInfoFactorization(int& ierr) const;
    
    template<class Vector1>
    void Solve(Vector1& x);
    
    template<class Vector1>
    void Solve(const SeldonTranspose& TransA, Vector1& x);

    template<class T1, class Alloc1>
    void Solve(Matrix<T1, General, ColMajor, Alloc1>& x);

#ifdef SELDON_WITH_MPI
    template<class Tint>
    void FactorizeDistributed(MPI::Comm& comm_facto,
                              Vector<Tint>& Ptr, Vector<Tint>& IndRow,
                              Vector<T>& Val, const IVect& glob_num,
                              bool sym, bool keep_matrix = false);
    
    template<class Vector1>
    void SolveDistributed(MPI::Comm& comm_facto, Vector1& x_solution,
                          const IVect& glob_number);
    
    template<class Vector1>
    void SolveDistributed(MPI::Comm& comm_facto,
			  const SeldonTranspose& TransA, Vector1& x_solution,
                          const IVect& glob_number);
#endif
    
    
  };

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                   Vector<T1, Storage1, Allocator1>& Y);


  template <class T, class Prop0, class Allocator0, class Allocator1>
  void GetAndSolveLU(Matrix<T, Prop0, ColSparse, Allocator0>& M,
		     Vector<T, VectFull, Allocator1>& Y);

  
  template <class T, class Prop0, class Allocator0, class Allocator1>
  void GetAndSolveLU(Matrix<T, Prop0, RowSparse, Allocator0>& M,
		     Vector<T, VectFull, Allocator1>& Y);
   
}

#define SELDON_FILE_COMPUTATION_SPARSE_DIRECT_SOLVER_HXX
#endif

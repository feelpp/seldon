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
    double permtol;
    //! Symmetric matrix.
    Matrix<T, Symmetric, ArrayRowSymSparse, Allocator> mat_sym;
    //! Unsymmetric matrix.
    Matrix<T, General, ArrayRowSparse, Allocator> mat_unsym;
    //! Permutation arrays.
    IVect permutation_row, permutation_col;
    //! Temporary vector.
    Vector<T, VectFull, Allocator> xtmp;
    //! if true the factorisation is contained in mat_sym
    bool symmetric_matrix;
    
  public :
    
    SparseSeldonSolver();
    
    void Clear();
    
    void HideMessages();
    void ShowMessages();
    
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
  
  
  //! basic class grouping different ordering strategies
  class SparseMatrixOrdering
  {
  public :
    // supported orderings
    enum {IDENTITY, REVERSE_CUTHILL_MCKEE, PORD,
	  SCOTCH, METIS, AMD, COLAMD, QAMD, USER, AUTO};    
  };
  
  
  template<class T, class Prop, class Storage, class Allocator,
	   class Tint, class Alloc>
  void FindSparseOrdering(Matrix<T, Prop, Storage, Allocator>& A,
			  Vector<Tint, VectFull, Alloc>& num, int type);
  
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
#ifdef SELDON_WITH_MUMPS
    MatrixMumps<T> mat_mumps; //!< Mumps solver
#endif
#ifdef SELDON_WITH_PASTIX
    MatrixPastix<T> mat_pastix; //!< Pastix solver
#endif
    
#ifdef SELDON_WITH_PRECONDITIONING
    //! ILUT solver
    IlutPreconditioning<double, T, NewAlloc<T> > mat_ilut;
#endif
    
    //! threshold for ilut solver
    double threshold_matrix;
    //! use of non-symmetric ilut ?
    bool enforce_unsym_ilut;
    
    //! default solver
    SparseSeldonSolver<T> mat_seldon;
        
  public :
    // available solvers
    enum {SELDON_SOLVER, UMFPACK, SUPERLU, MUMPS, PASTIX, ILUT};
    // error codes
    enum {FACTO_OK, STRUCTURALLY_SINGULAR_MATRIX,
          NUMERICALLY_SINGULAR_MATRIX, OUT_OF_MEMORY, INVALID_ARGUMENT,
          INCORRECT_NUMBER_OF_ROWS, MATRIX_INDICES_INCORRECT,
          INVALID_PERMUTATION, ORDERING_FAILED, INTERNAL_ERROR};
    
    SparseDirectSolver();
    
    void HideMessages();
    void ShowMessages();
    void ShowFullHistory();
    
    void Clear();
    
    int GetM() const;
    int GetN() const;
    
    int GetTypeOrdering() const;
    void SetPermutation(const IVect&);
    void SelectOrdering(int);
    
    void SetNbThreadPerNode(int m);
    
    template<class MatrixSparse>
    void ComputeOrdering(MatrixSparse& A);
    
    void SelectDirectSolver(int);
    void SetNonSymmetricIlut();

    int GetDirectSolver();

    void RefineSolution();
    void DoNotRefineSolution();
    
    double GetThresholdMatrix() const;
    
    template<class MatrixSparse>
    void Factorize(MatrixSparse& A, bool keep_matrix = false);
    
    int GetInfoFactorization(int& ierr) const;
    
    template<class Vector1>
    void Solve(Vector1& x);
    
    template<class TransStatus, class Vector1>
    void Solve(const TransStatus& TransA, Vector1& x);

#ifdef SELDON_WITH_MPI
    template<class Tint>
    void FactorizeDistributed(MPI::Comm& comm_facto,
                              Vector<Tint>& Ptr, Vector<Tint>& IndRow,
                              Vector<T>& Val, const IVect& glob_num,
                              bool sym, bool keep_matrix = false);
    
    template<class Vector1>
    void SolveDistributed(MPI::Comm& comm_facto, Vector1& x_solution,
                          const IVect& glob_number);
    
    template<class TransStatus, class Vector1>
    void SolveDistributed(MPI::Comm& comm_facto,
			  const TransStatus& TransA, Vector1& x_solution,
                          const IVect& glob_number);
#endif
    
    
  };

  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                   Vector<T1, Storage1, Allocator1>& Y);


  template <class T, class Prop0, class Allocator0, class Allocator1>
  void Solve(Matrix<T, Prop0, ColSparse, Allocator0>& M,
             Vector<T, VectFull, Allocator1>& Y);


  template <class T, class Prop0, class Allocator0, class Allocator1>
  void Solve(Matrix<T, Prop0, RowSparse, Allocator0>& M,
             Vector<T, VectFull, Allocator1>& Y);


}  // namespace Seldon.


#define SELDON_FILE_COMPUTATION_SPARSESOLVER_HXX
#endif

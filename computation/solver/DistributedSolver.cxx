#ifndef SELDON_FILE_DISTRIBUTED_SOLVER_CXX

#include "DistributedSolver.hxx"

namespace Seldon
{

  //! Default constructor
  template <class T>
  SparseDistributedSolver<T>::SparseDistributedSolver()
    : SparseDirectSolver<T>()
  {
    diagonal_scaling_left = false;
    diagonal_scaling_right = false;
    
#ifdef SELDON_WITH_PASTIX
    this->mat_pastix.RefineSolution();
#endif

#ifdef SELDON_WITH_MPI
    nodl_scalar_ = 1;
    nb_unknowns_scal_ = 1;
    comm_ = &MPI::COMM_SELF;
    ProcSharingRows_ = NULL;
    SharingRowNumbers_ = NULL;
#endif
  }


  //! Changes the level of displayed messages
  template <class T>
  void SparseDistributedSolver<T>::SetPrintLevel(int print_level)
  { 
    if (print_level >= 2)
      this->ShowMessages(); 
    else if (print_level >= 6)
      this->ShowFullHistory();
    else
      this->HideMessages();    
  }
  
  
  //! Clears LU matrices
  template <class T>
  void SparseDistributedSolver<T>::Clear()
  {
    SparseDirectSolver<T>::Clear();
    
    diagonal_scale_left.Clear();
    diagonal_scale_right.Clear();
  }
  
    
  //! scales matrix with 1 / \sum a_{i, j} (rows and columns)
  template<class T> template<class MatrixSparse>
  void SparseDistributedSolver<T>::ScaleMatrixRowCol(MatrixSparse& A)
  {
    int n = A.GetM();
    diagonal_scaling_left = true;
    diagonal_scaling_right = true;
    
    GetRowColSum(diagonal_scale_left, diagonal_scale_right, A);
    
    // forming D^-1 A C^-1 where d_i = \sum_j | a_ij |
    // and c_j = \sum_i | a_ij |
    for (int i = 0; i < n; i++)
      {
        if (diagonal_scale_left(i) != 0)
          diagonal_scale_left(i) = abs(1.0/diagonal_scale_left(i));
        else
          diagonal_scale_left(i) = 1.0;
        
        if (diagonal_scale_right(i) != 0)
          diagonal_scale_right(i) = abs(1.0/diagonal_scale_right(i));
        else
          diagonal_scale_right(i) = 1.0;
      }
    
    ScaleMatrix(A, diagonal_scale_left, diagonal_scale_right);
  }
  
  
  //! factorization of matrix
  template<class T> template<class MatrixSparse>
  void SparseDistributedSolver<T>::Factorize(MatrixSparse& A,
					     bool keep_matrix,
					     bool scale_matrix)
  {
    if (scale_matrix)
      ScaleMatrixRowCol(A);
    
    // then factorizing the modified matrix
    SparseDirectSolver<T>::Factorize(A, keep_matrix);
  }
  
  
  //! factorization of distributed matrix
  template<class T>
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void SparseDistributedSolver<T>
  ::Factorize(DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      bool keep_matrix, bool scale_matrix)
  {
    if (scale_matrix)
      ScaleMatrixRowCol(A);

#ifdef SELDON_WITH_MPI  
    MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      {
        comm_ = &comm;
        SparseDirectSolver<T>::Factorize(A, keep_matrix);
      }
    else
      {
        // factorisation of a distributed matrix
        comm_ = &comm;
        ProcSharingRows_ = &A.GetProcessorSharingRows();
        SharingRowNumbers_ = &A.GetSharingRowNumbers();
        nodl_scalar_ = A.GetNodlScalar();
        nb_unknowns_scal_ = A.GetNbScalarUnknowns();
        bool sym_matrix = IsSymmetricMatrix(A);
        
        // assembles distributed matrix
        if ((this->mat_seldon.GetPrintLevel() >= 1) &&(comm.Get_rank() == 0))
          cout << "Assembling the distributed matrix..." << endl;
        
        Vector<Int_wp> Ptr, IndRow;
        Vector<T> Val; Prop0 sym;
        bool sym_pattern = false;
        if (this->type_solver == this->PASTIX)
          sym_pattern = true;
        
        AssembleDistributed(A, sym, comm, global_col_numbers,
                            local_col_numbers,
			    Ptr, IndRow, Val, sym_pattern);
        
        if ((this->mat_seldon.GetPrintLevel() >= 1) &&(comm.Get_rank() == 0))
          cout << "Factorizing the distributed matrix..." << endl;
        
        // factorizes the matrix
        this->FactorizeDistributed(comm, Ptr, IndRow, Val,
                                   global_col_numbers, sym_matrix);
      }
#else
    SparseDirectSolver<T>::Factorize(A, keep_matrix);
#endif
  }
  
  
#ifdef SELDON_WITH_MPI
  //! assembles the vector (adds values of shared rows)
  template<class T> template<class T2>
  void SparseDistributedSolver<T>::AssembleVec(Vector<T2>& X) const
  {
    AssembleVector(X, MPI::SUM, *ProcSharingRows_, *SharingRowNumbers_,
                   *comm_, nodl_scalar_, nb_unknowns_scal_, 20);
  }

  //! assembles the vector (adds values of shared rows)
  template<class T> template<class T2>
  void SparseDistributedSolver<T>::AssembleVec(Matrix<T2, General, ColMajor>& A) const
  {
    if (nb_unknowns_scal_ > 1)
      {
	cout << "not implemented " << endl;
	abort();
      }

    int nrhs = A.GetN();
    Vector<T2> X;
    X.SetData(nrhs*A.GetM(), A.GetData());
    AssembleVector(X, MPI::SUM, *ProcSharingRows_, *SharingRowNumbers_,
                   *comm_, A.GetM(), nrhs, 21);
    
    X.Nullify();
  }
#endif
  
  
  //! solution of linear system Ax = b by using LU factorization
  //! (with scaling)
  /*!
    \param[out] x_solution solution
    \param[in] b_rhs right hand side
   */
  template<class T> template<class Vector1>
  void SparseDistributedSolver<T>::Solve(Vector1& x_solution,
					 const Vector1& b_rhs)
  {
    Copy(b_rhs, x_solution);
    Solve(x_solution);
  }
  
  
  //! solution of linear system Ax = b by using LU factorization
  //! (with scaling)
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T> template<class T1>
  void SparseDistributedSolver<T>::Solve(Vector<T1>& x_solution)
  {
    if (diagonal_scaling_left)
      for (int i = 0; i < x_solution.GetM(); i++)
        x_solution(i) *= diagonal_scale_left(i);
    
#ifdef SELDON_WITH_MPI
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      SparseDirectSolver<T>::Solve(x_solution);
    else
      {
        // extracting right hand side (we remove overlapped dofs)
        int n = local_col_numbers.GetM();
        Vector<T1> x_sol_extract(n);
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_sol_extract(i) = x_solution(local_col_numbers(i));
        
        SolveLU_Distributed(comm, SeldonNoTrans, *this,
                            x_sol_extract, global_col_numbers);
        
        x_solution.Fill(0);
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_solution(local_col_numbers(i)) = x_sol_extract(i);
        
        // adding overlapped components
        this->AssembleVec(x_solution);
      }
#else
    SparseDirectSolver<T>::Solve(x_solution);
#endif
    
    if (diagonal_scaling_right)
      for (int i = 0; i < x_solution.GetM(); i++)
        x_solution(i) *= diagonal_scale_right(i);
  }
  
  
  //! solution of linear system A^T x = b by using LU factorization
  //! (with scaling)
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T> template<class T1>
  void SparseDistributedSolver<T>::TransSolve(Vector<T1>& x_solution)
  {
    if (diagonal_scaling_right)
      for (int i = 0; i < x_solution.GetM(); i++)
        x_solution(i) *= diagonal_scale_right(i);

#ifdef SELDON_WITH_MPI
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      SparseDirectSolver<T>::Solve(SeldonTrans, x_solution);
    else
      {
        // extracting right hand side (we remove overlapped dofs)
        int n = local_col_numbers.GetM();
        Vector<T1> x_sol_extract(n);
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_sol_extract(i) = x_solution(local_col_numbers(i));
        
        SolveLU_Distributed(comm, SeldonTrans, *this,
                            x_sol_extract, global_col_numbers);
        
        x_solution.Fill(0);
        for (int i = 0; i < local_col_numbers.GetM(); i++)
          x_solution(local_col_numbers(i)) = x_sol_extract(i);
        
        // adding overlapped components
        this->AssembleVec(x_solution);
      }
#else
    SparseDirectSolver<T>::Solve(SeldonTrans, x_solution);
#endif

    if (diagonal_scaling_left)
      for (int i = 0; i < x_solution.GetM(); i++)
        x_solution(i) *= diagonal_scale_left(i);
  }
  
  
  //! solution of linear system Ax = b by using LU factorization (with scaling)
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
   */
  template<class T> template<class T1>
  void SparseDistributedSolver<T>::Solve(Matrix<T1, General, ColMajor>& x_solution)
  {
    if (diagonal_scaling_left)
      ScaleLeftMatrix(x_solution, diagonal_scale_left);
    
#ifdef SELDON_WITH_MPI
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      SparseDirectSolver<T>::Solve(x_solution);
    else
      {
        // extracting right hand side (we remove overlapped dofs)
	int n = local_col_numbers.GetM();
	int N = x_solution.GetM(), nrhs = x_solution.GetN();
        Matrix<T1, General, ColMajor> x_sol_extract(n, nrhs);
	for (int k = 0; k < nrhs; k++)
	  for (int i = 0; i < local_col_numbers.GetM(); i++)
	    x_sol_extract(i, k) = x_solution(local_col_numbers(i), k);
        
	x_solution.Clear();
        SolveLU_Distributed(comm, SeldonNoTrans, *this,
                            x_sol_extract, global_col_numbers);
        
	x_solution.Reallocate(N, nrhs);
        x_solution.Fill(0);
	for (int k = 0; k < nrhs; k++)
	  for (int i = 0; i < local_col_numbers.GetM(); i++)
	    x_solution(local_col_numbers(i), k) = x_sol_extract(i, k);
        
        // adding overlapped components
        this->AssembleVec(x_solution);
      }
#else
    SparseDirectSolver<T>::Solve(x_solution);
#endif
    
    if (diagonal_scaling_right)
      ScaleLeftMatrix(x_solution, diagonal_scale_right);
  }
  

#ifdef SELDON_WITH_MPI
  template<class T>
  void SolveLU_Distributed(MPI::Comm& comm, const SeldonTranspose& transA,
			   SparseDistributedSolver<T>& mat_lu,
                           Vector<T>& x, Vector<int>& global_col)
  {
    mat_lu.SolveDistributed(comm, transA, x, global_col);
  }

  
  template<class T>
  void SolveLU_Distributed(MPI::Comm& comm, const SeldonTranspose& transA,
                           SparseDistributedSolver<complex<T> >& mat_lu,
                           Vector<T>& x, Vector<int>& global_col)
  {
    throw WrongArgument("SolveLU(SparseDistributedSolver, Vector, Vector)",
			"the result must be complex");
  }
  
  
  template<class T>
  void SolveLU_Distributed(MPI::Comm& comm, const SeldonTranspose& transA,
			   SparseDistributedSolver<T>& mat_lu,
                           Vector<complex<T> >& x, Vector<int>& global_col)
  {
    Vector<T> xreal(x.GetM()), ximag(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
        xreal(i) = real(x(i));
        ximag(i) = imag(x(i));
      }
    
    mat_lu.SolveDistributed(comm, transA, xreal, global_col);
    mat_lu.SolveDistributed(comm, transA, ximag, global_col);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T>(xreal(i), ximag(i));
  }
  
  
  template<class T>
  void SolveLU_Distributed(MPI::Comm& comm, const SeldonTranspose& transA,
			   SparseDistributedSolver<T>& mat_lu,
                           Matrix<T, General, ColMajor>& x, Vector<int>& global_col)
  {
    mat_lu.SolveDistributed(comm, transA, x, global_col);
  }

  
  template<class T>
  void SolveLU_Distributed(MPI::Comm& comm, const SeldonTranspose& transA,
                           SparseDistributedSolver<complex<T> >& mat_lu,
                           Matrix<T, General, ColMajor>& x, Vector<int>& global_col)
  {
    throw WrongArgument("SolveLU(SparseDistributedSolver, Vector, Vector)",
			"the result must be complex");
  }
  
  
  template<class T>
  void SolveLU_Distributed(MPI::Comm& comm, const SeldonTranspose& transA,
			   SparseDistributedSolver<T>& mat_lu,
                           Matrix<complex<T>, General, ColMajor>& x,
			   Vector<int>& global_col)
  {
    throw WrongArgument("SolveLU(SparseDistributedSolver, Vector, Vector)",
			"the matrix must be real");
  }
#endif  


  //! Computation of a Schur Complement
  /*!
    \param[in] mat_direct initial matrix
    \param[in] num row numbers of the Schur Complement
    \param[out] mat_schur Schur Complement
   */
  template<class T> template<class MatrixSparse, class MatrixFull>
  void SparseDistributedSolver<T>::GetSchurComplement(MatrixSparse& mat_direct, 
						      const IVect& num,
						      MatrixFull& mat_schur)
  {
    if (this->type_solver == this->MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	GetSchurMatrix(mat_direct, this->mat_mumps, num, mat_schur);
#else
	cout<<"Recompile Montjoie with Mumps."
	    << " Schur complement can't be performed otherwise"<<endl;
	abort();
#endif
      }
    else
      {
        cout << "Try to use Mumps, not implemented with other solvers " <<endl;
        abort();
      }
  }

  
  //! returns in byte memory used by the direct solver
  template<class T>
  int64_t SparseDistributedSolver<T>::GetMemorySize() const
  {
    int64_t taille = sizeof(Treal)*(diagonal_scale_left.GetM()
				    +diagonal_scale_right.GetM());
    
    taille += this->permut.GetM()*sizeof(int);
#ifdef SELDON_WITH_UMFPACK
    taille += this->mat_umf.GetMemorySize();
#endif
#ifdef SELDON_WITH_SUPERLU
    taille += this->mat_superlu.GetMemorySize();
#endif
#ifdef SELDON_WITH_PARDISO
    taille += this->mat_pardiso.GetMemorySize();
#endif
#ifdef SELDON_WITH_MUMPS
    taille += this->mat_mumps.GetMemorySize();
#endif
#ifdef SELDON_WITH_PASTIX
    taille += this->mat_pastix.GetMemorySize();
#endif
    
#ifdef SELDON_WITH_PRECONDITIONING
    taille += this->mat_ilut.GetMemorySize();
#endif
    
    taille += this->mat_seldon.GetMemorySize();

    return taille;
  }
  
  
  //! factorizes a matrix
  template<class T, class MatrixSparse>
  void GetLU(SparseDistributedSolver<T>& mat_lu, MatrixSparse& A,
             bool keep_matrix, bool scale_matrix, T& x_test)
  {
    mat_lu.Factorize(A, keep_matrix, scale_matrix);
  }


  //! factorizes a matrix
  template<class T, class MatrixSparse, class T0>
  void GetLU(SparseDistributedSolver<T>& mat_lu, MatrixSparse& A,
             bool keep_matrix, bool scale_matrix, T0& x_test)
  {
    throw WrongArgument("GetLU(SparseDistributedSolver)",
			"the types are not compatible");
  }


  //! factorizes a matrix
  template<class T, class MatrixSparse, class T0>
  void GetLU(SparseDistributedSolver<T>& mat_lu, MatrixSparse& A,
             bool keep_matrix, bool scale_matrix)
  {
    typename MatrixSparse::entry_type x;
    GetLU(mat_lu, A, keep_matrix, scale_matrix, x);
  }
  
  
  //! solves A x = b
  template<class T0, class T1>
  void SolveLU(SparseDistributedSolver<T0>& mat_lu,
	       Vector<T1>& x, const Vector<T1>& b)
  {
    mat_lu.Solve(x, b);
  }


  //! solves A x = b
  template<class T>
  void SolveLU(SparseDistributedSolver<complex<T> >& mat_lu,
	       Vector<T>& x, const Vector<T>& b)
  {
    throw WrongArgument("SolveLU(SparseDistributedSolver, Vector, Vector)",
			"the result must be complex");
  }


  //! solves A x = b
  template<class T0, class T1>
  void SolveLU(SparseDistributedSolver<T0>& mat_lu,
	       Matrix<T1, General, ColMajor>& x)
  {
    mat_lu.Solve(x);
  }


  //! solves A x = b
  template<class T>
  void SolveLU(SparseDistributedSolver<complex<T> >& mat_lu,
	       Matrix<T, General, ColMajor>& x)
  {
    throw WrongArgument("SolveLU(SparseDistributedSolver, Vector, Vector)",
			"the result must be complex");
  }
  
} // namespace Seldon

#define SELDON_FILE_DISTRIBUTED_SOLVER_CXX
#endif

#ifndef SELDON_FILE_DISTRIBUTED_SOLVER_HXX

namespace Seldon
{

  // int used for pastix
#ifdef INTSIZE64
  typedef int64_t Int_wp;
#else
  typedef int32_t Int_wp;
#endif

  //! general class for direct solver
  template<class T>
  class SparseDistributedSolver : public SparseDirectSolver<T>
  {
  protected :
    typedef typename ClassComplexType<T>::Treal Treal;
    
    bool diagonal_scaling_left; //!< left scaling ?
    bool diagonal_scaling_right; //!< right scaling ?
    Vector<Treal> diagonal_scale_left; //!< left scaling
    Vector<Treal> diagonal_scale_right; //!< right scaling

#ifdef SELDON_WITH_MPI
    int nodl_scalar_, nb_unknowns_scal_;
    MPI::Comm* comm_;
    IVect* ProcSharingRows_;
    Vector<IVect>* SharingRowNumbers_;
    IVect global_col_numbers, local_col_numbers;
    
    template<class T2>
    void AssembleVec(Vector<T2>& X) const;

    template<class T2>
    void AssembleVec(Matrix<T2, General, ColMajor>& A) const;
#endif

    template<class MatrixSparse>
    void ScaleMatrixRowCol(MatrixSparse& A);

  public :
    
    SparseDistributedSolver();
    
    void SetPrintLevel(int print);
    
    void Clear();

    template<class MatrixSparse>
    void Factorize(MatrixSparse& A, bool keep_matrix = false,
                   bool scale_matrix = false);
    
    template<class T0, class Prop0, class Storage0, class Allocator0>
    void Factorize(DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
                   bool keep_matrix = false, bool scale_matrix = false);
    
    template<class Vector1>
    void Solve(Vector1& x_solution, const Vector1& b_rhs);
    
    template<class T1>
    void Solve(Vector<T1>& x_solution);

    template<class T1>
    void TransSolve(Vector<T1>& x_solution);

    template<class T1>
    void Solve(Matrix<T1, General, ColMajor>& x_solution);    
    
    template<class MatrixSparse, class MatrixFull>
    void GetSchurComplement(MatrixSparse& mat_direct,
			    const IVect& num, MatrixFull& mat_schur);
    
    int64_t GetMemorySize() const;
        
  };


#ifdef SELDON_WITH_MPI
  template<class TransA, class T>
  void SolveLU_Distributed(MPI::Comm& comm, const TransA& transA,
			   SparseDistributedSolver<T>& mat_lu,
                           Vector<T>& x, Vector<int>& global_col);

  template<class TransA, class T>
  void SolveLU_Distributed(MPI::Comm& comm, const TransA& transA,
                           SparseDistributedSolver<complex<T> >& mat_lu,
                           Vector<T>& x, Vector<int>& global_col);

  template<class TransA, class T>
  void SolveLU_Distributed(MPI::Comm& comm, const TransA& transA,
			   SparseDistributedSolver<T>& mat_lu,
                           Vector<complex<T> >& x, Vector<int>& global_col);
#endif    

  template<class T, class MatrixSparse>
  void GetLU(SparseDistributedSolver<T>& mat_lu, MatrixSparse& A,
             bool keep_matrix, bool scale_matrix, T& x_test);

  template<class T, class MatrixSparse, class T0>
  void GetLU(SparseDistributedSolver<T>& mat_lu, MatrixSparse& A,
             bool keep_matrix, bool scale_matrix, T0& x_test);

  template<class T, class MatrixSparse, class T0>
  void GetLU(SparseDistributedSolver<T>& mat_lu, MatrixSparse& A,
             bool keep_matrix = false, bool scale_matrix = false);

  template<class T0, class T1>
  void SolveLU(SparseDistributedSolver<T0>& mat_lu,
	       Vector<T1>& x, const Vector<T1>& b);

  template<class T>
  void SolveLU(SparseDistributedSolver<complex<T> >& mat_lu,
	       Vector<T>& x, const Vector<T>& b);

  template<class T0, class T1>
  void SolveLU(SparseDistributedSolver<T0>& mat_lu,
	       Matrix<T1, General, ColMajor>& x);

  template<class T>
  void SolveLU(SparseDistributedSolver<complex<T> >& mat_lu,
	       Matrix<T, General, ColMajor>& x);
  
}

#define SELDON_FILE_DISTRIBUTED_SOLVER_HXX
#endif
  

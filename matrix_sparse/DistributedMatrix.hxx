#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_HXX

namespace Seldon
{
  
  //! matrix distributed over all the processors
  template<class T, class Prop, class Storage,
           class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class DistributedMatrix : public Matrix<T, Prop, Storage, Allocator>
  {
    template<class T0, class Prop0, class Storage0, class Allocator0>
    friend class DistributedMatrix;
    
    template<class T0, class Prop0, class Storage0, class Allocator0>
    friend class DistributedMatrix_BlockDiag;
    
  public :
    //! type of non-zero entries (double, complex<double>, etc)
    typedef typename Matrix<T, Prop,
			    Storage, Allocator>::entry_type entry_type;
    
  protected :
    //! row numbers shared with other processors
    /*!
      It is the same array required in a distributed vector.
      Some rows can be shared by several processors.
    */
    IVect* OverlapRowNumbers;
    
    //! processor where each shared row should be assembled
    /*!
      For rows which have already been counted in other processors,
      you specify on which processor these rows are assembled
    */
    IVect* OverlapProcNumbers;
    
    //! global row numbers 
    IVect* GlobalRowNumbers;
    
    //! list of processors sharing rows with the current one
    IVect* ProcSharingRows;
    
    //! for each processor sharing rows, list of "local numbers" shared
    /*!
      It is assumed that these local numbers are "sorted", such that
      local numbers shared with processor j on processor i are corresponding
      with local numbers shared with processor i on processor j
    */
    Vector<IVect>* SharingRowNumbers;
    
    //! number of "scalar" unknowns
    int nodl_scalar_, nb_unknowns_scal_;
    
    //! total number of rows (on all processors)
    int nglob_;
    
    //! MPI communicator
    MPI::Comm* comm_;
    
    //! additional values on rows with non-local columns
    Vector<Vector<entry_type, VectSparse>, VectFull,
           NewAlloc<Vector<entry_type, VectSparse> > > dist_col;
    
    //! additional values on columns with non-local rows
    Vector<Vector<entry_type, VectSparse>, VectFull,
                  NewAlloc<Vector<entry_type, VectSparse> > > dist_row;
    
    //! distant processor for additional values
    Vector<IVect> proc_col, proc_row;
    
    //! global row/col numbers (needed for MltAdd)
    IVect global_row_to_recv, global_col_to_recv;
    IVect ptr_global_row_to_recv, ptr_global_col_to_recv;
    
    //! local row/col numbers (needed for MltAdd)
    Vector<IVect> local_row_to_send, local_col_to_send;
    
    //! processor numbers (needed for MltAdd)
    IVect proc_col_to_recv, proc_col_to_send,
      proc_row_to_recv, proc_row_to_send;
    
    //! if true local numbers are present in dist_row/dist_col 
    //! instead of global numbers
    bool local_number_distant_values;
    
    //! number of distant non-zero entries
    int size_max_distant_row, size_max_distant_col;
    
    // internal functions
    void EraseArrayForMltAdd();
    void SwitchToGlobalNumbers();
    
    template<class TypeDist>
    void SortAndAssembleDistantInteractions(TypeDist& dist_val,
					    Vector<IVect>& dist_proc,
                                            IVect& glob_num,
					    IVect& ptr_glob_num,
					    IVect& proc_glob,
                                            Vector<IVect>& local_num,
					    IVect& proc_local);
    
    template<class T2>
    void ScatterValues(const Vector<T2>& X, const IVect& num_recv,
                       const IVect&, const IVect& proc_recv,
                       const Vector<IVect>& num_send,
                       const IVect& proc_send, Vector<T2>& Xcol) const;
    
    template<class T2>
    void AssembleValues(const Vector<T2>& Xcol, const IVect& num_recv,
                        const IVect&, const IVect& proc_recv,
                        const Vector<IVect>& num_send,
                        const IVect& proc_send, Vector<T2>& X) const;
    
    void AssembleValuesMin(const IVect& Xcol, const IVect& Xcol_proc,
                           const IVect& num_recv, const IVect& ptr_num_recv,
                           const IVect& proc_recv,
                           const Vector<IVect>& num_send,
			   const IVect& proc_send,
                           IVect& Y, IVect& Yproc) const;
    
    void AssembleVecMin(Vector<int>& X, Vector<int>& Xproc) const;
    
    template<class T0, class TypeDist>
    void RemoveSmallEntryDistant(const T0&, TypeDist&, Vector<IVect>&);
    
    template<class T0> void GetRowSumDistantRow(Vector<T0>& vec_sum) const;
    template<class T0> void GetRowSumDistantCol(Vector<T0>& vec_sum) const;

    template<class T0> void GetColSumDistantRow(Vector<T0>& vec_sum) const;
    template<class T0> void GetColSumDistantCol(Vector<T0>& vec_sum) const;
    
  public :
    // constructors
    DistributedMatrix();
    DistributedMatrix(int m, int n);
    
    MPI::Comm& GetCommunicator();
    const MPI::Comm& GetCommunicator() const;
  
    int64_t GetMemorySize() const;

    // initialisation of pointers
    void Init(int n, IVect*, IVect*, IVect*,
              int, int, IVect*, Vector<IVect>*, MPI::Comm&);
  
    template<class T0, class Prop0, class Storage0, class Allocator0>
    void Init(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>&);
    
    void Init(Vector<IVect>&, IVect&, IVect&, IVect&,
              int, int, IVect&, Vector<IVect>&, MPI::Comm&,
              bool distribute_row = true);

    void Init(IVect&, IVect&, IVect&,
              int, int, IVect&, Vector<IVect>&, MPI::Comm&);
    
    // memory management
    void Reallocate(int m, int n);
    void Resize(int m, int n);
    void Clear();
    
    // basic functions
    DistributedMatrix<T, Prop, Storage, Allocator>&
    operator=(const DistributedMatrix<T, Prop, Storage, Allocator>& X);
    
    template<class T0>
    DistributedMatrix<T, Prop, Storage, Allocator>& operator *=(const T0& x);
    
    int GetGlobalM() const;
    IVect& GetGlobalRowNumber();
    const IVect& GetGlobalRowNumber() const;
    IVect& GetOverlapRowNumber();
    const IVect& GetOverlapRowNumber() const;
    IVect& GetOverlapProcNumber();
    const IVect& GetOverlapProcNumber() const;
    int GetNodlScalar() const;
    int GetNbScalarUnknowns() const;
    IVect& GetProcessorSharingRows();
    Vector<IVect>& GetSharingRowNumbers();

    // convenient functions overloaded    
    int GetNonZeros() const;
    int GetDataSize() const;    
    
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);
    
    void SetIdentity();
    void Zero();
    void Fill();
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName, bool cplx = false) const;
    void WriteText(ostream& FileStream, bool cplx = false) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName, bool cplx = false);
    void ReadText(istream& FileStream, bool cplx = false);
    
    // adding values to matrix
    void AddDistantInteraction(int i, int jglob, int proc,
                               const entry_type& val);
    
    void AddRowDistantInteraction(int iglob, int j, int proc,
                                  const entry_type& val);

    // functions for matrix-vector product
    int GetMaxDataSizeDistantCol() const;
    int GetMaxDataSizeDistantRow() const;
    bool IsReadyForMltAdd() const;
    
    void PrepareMltAdd();
    
    template<class T2>
    void ScatterRowValues(const Vector<T2>& X, Vector<T2>& Xcol) const;

    template<class T2>
    void ScatterColValues(const Vector<T2>& X, Vector<T2>& Xcol) const;
    
    template<class T2>
    void AssembleRowValues(const Vector<T2>& Xrow, Vector<T2>& X) const;

    template<class T2>
    void AssembleColValues(const Vector<T2>& Xrow, Vector<T2>& X) const;
    
    template<class T2>
    void AssembleVec(Vector<T2>&) const;
    
    template<class T2, class Storage2, class Allocator2,
             class T4, class Storage4, class Allocator4>
    void MltAddCol(const class_SeldonNoTrans& Trans,
                   const Vector<T2, Storage2, Allocator2>& X,
                   Vector<T4, Storage4, Allocator4>& Y) const;

    template<class T2, class Storage2, class Allocator2,
             class T4, class Storage4, class Allocator4>
    void MltAddCol(const class_SeldonTrans& Trans,
                   const Vector<T2, Storage2, Allocator2>& X,
                   Vector<T4, Storage4, Allocator4>& Y) const;

    template<class T2, class Storage2, class Allocator2,
             class T4, class Storage4, class Allocator4>
    void MltAddCol(const class_SeldonConjTrans& Trans,
                   const Vector<T2, Storage2, Allocator2>& X,
                   Vector<T4, Storage4, Allocator4>& Y) const;
    
    template<class T2, class Storage2, class Allocator2,
             class T4, class Storage4, class Allocator4>
    void MltAddRow(const class_SeldonNoTrans& Trans,
                   const Vector<T2, Storage2, Allocator2>& X,
                   Vector<T4, Storage4, Allocator4>& Y) const;

    template<class T2, class Storage2, class Allocator2,
             class T4, class Storage4, class Allocator4>
    void MltAddRow(const class_SeldonTrans& Trans,
                   const Vector<T2, Storage2, Allocator2>& X,
                   Vector<T4, Storage4, Allocator4>& Y) const;

    template<class T2, class Storage2, class Allocator2,
             class T4, class Storage4, class Allocator4>
    void MltAddRow(const class_SeldonConjTrans& Trans,
                   const Vector<T2, Storage2, Allocator2>& X,
                   Vector<T4, Storage4, Allocator4>& Y) const;
    
    // functions for assembling matrix
    template<class T0, class Allocator0>
    void GetDistributedRows(Matrix<T0, General,
			    ArrayRowSparse, Allocator0>& rows) const;
    
    template<class T0, class Allocator0>
    void GetDistributedColumns(Matrix<T0, General,
			       ArrayColSparse, Allocator0>& rows,
			       bool sym_pattern) const;
    
    template<class MatrixSparse, class Tint, class T0> friend void
    AssembleDistributed(MatrixSparse& A,
			Symmetric& sym, const MPI::Comm& comm,
                        IVect& row_numbers, IVect& local_row_numbers,
                        Vector<Tint>& PtrA, Vector<Tint>& IndA,
                        Vector<T0>& ValA, bool sym_pattern);
    
    template<class MatrixSparse, class Tint, class T0> friend void
    AssembleDistributed(MatrixSparse& A,
			General& prop, const MPI::Comm& comm,
                        IVect& col_numbers, IVect& local_col_numbers,
                        Vector<Tint>& PtrA, Vector<Tint>& IndA,
                        Vector<T0>& ValA, bool sym_pattern);

    // friend functions
    template<class T1, class Prop1, class Storage1, class Allocator1> friend
    void MltMin(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
		IVect& Y, IVect& Yproc);
    
    template<class T0, class T1, class Prop1, class Storage1,
	     class Allocator1, class T2, class Prop2, class Storage2,
	     class Allocator2> friend void
    Add(const T0& alpha,
	const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A, 
        DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B);
    
    template<class T1, class Prop1, class Storage1, class Allocator1>
    friend typename ClassComplexType<T1>::Treal
    MaxAbs(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

    template<class T0, class T1, class Prop1,
	     class Storage1, class Allocator1>
    friend void
    GetRowSum(Vector<T0>& vec_sum,
	      const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);
    
    template<class T0, class T1, class Prop1,
	     class Storage1, class Allocator1>
    friend void
    GetColSum(Vector<T0>& vec_sum,
	      const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);
 
    template<class T0, class T1, class Prop1,
	     class Storage1, class Allocator1>
    friend void 
    GetRowColSum(Vector<T0>& row_sum, Vector<T0>& col_sum,
		 const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

    template<class T1, class Prop1, class Storage1, class Allocator1>
    friend void
    Transpose(DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

    template<class T1, class Prop1, class Storage1, class Allocator1>
    friend void 
    Transpose(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
	      DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B);
    
    template<class T1, class Prop1, class Storage1, class Allocator1>
    friend void 
    Conjugate(DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);
    
    template<class T0, class Storage0, class Allocator0,
             class T1, class Allocator1>
    friend void 
    ScaleRightMatrix(DistributedMatrix<T0, General, Storage0, Allocator0>& A,
		     const Vector<T1, VectFull, Allocator1>& coef);
    
    template<class T0, class Storage0, class Allocator0,
             class T1, class Allocator1>
    friend void 
    ScaleLeftMatrix(DistributedMatrix<T0, General, Storage0, Allocator0>& A,
		    const Vector<T1, VectFull, Allocator1>& coef);

    template<class T0, class Prop0, class Storage0, class Allocator0,
             class T1, class Allocator1, class T2, class Allocator2>
    friend void 
    ScaleMatrix(DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		const Vector<T1, VectFull, Allocator1>& Drow,
		const Vector<T2, VectFull, Allocator2>& Dcol);
    
    template<class T1, class Prop0, class Storage0, class Allocator0>
    friend void
    EraseCol(const IVect& col_number,
             DistributedMatrix<T1, Prop0, Storage0, Allocator0>& A);
    
    template<class T1, class Prop0, class Storage0, class Allocator0>
    friend void
    EraseRow(const IVect& col_number,
             DistributedMatrix<T1, Prop0, Storage0, Allocator0>& A);
    
    template<class T1, class Prop1, class Storage1, class Allocator1,
             class T2, class Prop2, class Storage2, class Allocator2>
    friend void 
    Copy(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
	 DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B);
    
    template<class T0, class Prop0, class Storage0, class Allocator0,
             class T1, class Prop1, class Storage1, class Allocator1>
    friend void
    CopySubMatrix(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>&,
                  const IVect& row, const IVect& col,
                  DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B);
    
    template<class T0, class Prop0, class Storage0, class Allocator0,
             class T1, class Prop1, class Storage1, class Allocator1>
    friend void
    ConvertToBlockDiagonal(const DistributedMatrix<T0, Prop0,
			   Storage0, Allocator0>& A,
                           DistributedMatrix<T1, Prop1,
			   Storage1, Allocator1>& B);
    
  };


  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3& beta,
              Vector<T4, Storage4, Allocator4>& Yres, bool assemble = true);

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const class_SeldonNoTrans& Trans,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3& beta,
              Vector<T4, Storage4, Allocator4>& Y, bool assemble = true);

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const class_SeldonTrans& Trans,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3& beta,
              Vector<T4, Storage4, Allocator4>& Yres, bool assemble = true);

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const class_SeldonConjTrans& Trans,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3& beta,
              Vector<T4, Storage4, Allocator4>& Yres, bool assemble = true);
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& M,
	   const Vector<T1, Storage1, Allocator1>& X,
	   Vector<T2, Storage2, Allocator2>& Y, bool assemble = true);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T3& alpha,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble = true);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T3& alpha,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<complex<T3>, Storage3, Allocator3>& Y,
	   bool assemble = true);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(int alpha,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble = true);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonNoTrans& Trans,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble = true);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonTrans& Trans,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble = true);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonConjTrans& Trans,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble = true);  

  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  Norm1(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  NormInf(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A); 

  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(const DistributedMatrix<T, Prop,
		     Storage, Allocator>& A,
                     DistributedMatrix<T, Prop, Storage, Allocator>& B);
  
  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(DistributedMatrix<T, Prop, Storage, Allocator>& A);

  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T4, class Prop4, class Storage4, class Allocator4>
  void Mlt(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
           const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
           DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C);

  template<class T0,
           class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T3,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
              const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
              const T3& beta,
              DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C);

  template<class T0, class TransA,
           class T1, class Prop1, class Storage1, class Allocator1,
           class TransB, class T2, class Prop2,
	   class Storage2, class Allocator2, class T3,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha, const TransA& transA,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
              const TransB& transB,
              const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
              const T3& beta,
              DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C);

  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void GetRow(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
              int i, Vector<T1, VectSparse, Allocator1>& X);

  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void GetCol(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
              int i, Vector<T1, VectSparse, Allocator1>& X);

  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
              int i, DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A);

  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
              int i, DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A);

  template<class T, class Prop, class Storage, class Allocator>
  void ApplyPermutation(DistributedMatrix<T, Prop, Storage, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);

  template<class T, class Prop, class Storage, class Allocator>
  void 
  ApplyInversePermutation(DistributedMatrix<T, Prop, Storage, Allocator>& A,
			  const Vector<int>& row_perm,
			  const Vector<int>& col_perm);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SOR(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   Vector<T2, Storage2, Allocator2>& X,
	   const Vector<T1, Storage1, Allocator1>& B,
	   const T3& omega, int iter, int type_ssor = 2);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SOR(const class_SeldonTrans& transM,
	   const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   Vector<T2, Storage2, Allocator2>& X,
	   const Vector<T1, Storage1, Allocator1>& B,
	   const T3& omega, int iter, int type_ssor = 3);

  template<class T1, class Prop, class Storage, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const DistributedMatrix<T1, Prop, Storage, Allocator>& A,
              const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>,
	      VectSparse, Allocator3>& V);

  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void Copy(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
            DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B);
  
  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void
  CopyReal(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
	   DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B);
  
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                    int m, int n,
		    DistributedMatrix<T, Prop, Storage, Allocator>& B);
  
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const Matrix<T, Prop, Storage, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    Matrix<T, Prop, Storage, Allocator>& B);

  template<class T, class Prop, class Storage, class Allocator>
  typename ClassComplexType<T>::Treal
  NormFro(const DistributedMatrix<T, Prop, Storage, Allocator>& A);

}

#define SELDON_FILE_DISTRIBUTED_MATRIX_HXX
#endif


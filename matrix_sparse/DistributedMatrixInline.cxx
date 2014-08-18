#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_INLINE_CXX

#include "DistributedMatrix.hxx"

namespace Seldon
{
  
  //! default constructor
  template<class T, class Prop, class Storage, class Allocator>
  inline DistributedMatrix<T, Prop, Storage, Allocator>::DistributedMatrix()
    : Matrix<T, Prop, Storage, Allocator>()
  {
    GlobalRowNumbers = NULL;
    OverlapProcNumbers = NULL;
    OverlapRowNumbers = NULL;
    ProcSharingRows = NULL;
    SharingRowNumbers = NULL;
    nodl_scalar_ = 0;
    nb_unknowns_scal_ = 1;
    nglob_ = 0;
    comm_ = &MPI::COMM_SELF;
    
    local_number_distant_values = false;
    size_max_distant_row = 0;
    size_max_distant_col = 0;
  }
  
  
  //! construction of an m by n matrix
  /*!
    Here m and n are the number of rows and columns of the local matrix
  */
  template<class T, class Prop, class Storage, class Allocator>
  inline DistributedMatrix<T, Prop, Storage, Allocator>::
  DistributedMatrix(int m, int n)
    : Matrix<T, Prop, Storage, Allocator>(m, n)
  {
    GlobalRowNumbers = NULL;
    OverlapProcNumbers = NULL;
    OverlapRowNumbers = NULL;
    ProcSharingRows = NULL;
    SharingRowNumbers = NULL;
    nglob_ = m;
    nodl_scalar_ = 0;
    nb_unknowns_scal_ = 1;
    comm_ = &MPI::COMM_SELF;
    
    dist_col.Reallocate(m);
    dist_row.Reallocate(n);
    proc_col.Reallocate(m);
    proc_row.Reallocate(n);

    local_number_distant_values = false;
    size_max_distant_row = 0;
    size_max_distant_col = 0;
  }

  
  //! returns MPI communicator (processors that will share the matrix)
  template<class T, class Prop, class Storage, class Allocator>
  inline MPI::Comm& DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetCommunicator()
  {
    return *comm_;
  }
  
  
  //! returns MPI communicator (processors that will share the matrix)
  template<class T, class Prop, class Storage, class Allocator>
  inline const MPI::Comm& DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetCommunicator() const
  {
    return *comm_;
  }
  
  
  //! Initialisation of pointers
  /*!
    This method is mandatory, otherwise pointers are set to NULL
    \param[in] n global number of rows (as if it was on a single processor)
    \param[in] row_num local to global numbering
    \param[in] overlap_num rows already counted in another processor
    \param[in] proc_num original processor for overlapped rows
  */
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>::
  Init(int n, IVect* row_num, IVect* overlap_num, IVect* proc_num,
       int Nvol, int nb_u, IVect* MatchingProc,
       Vector<IVect>* MatchingDofNumber, MPI::Comm& comm)
  {
    nglob_ = n;
    GlobalRowNumbers = row_num;
    OverlapRowNumbers = overlap_num;
    OverlapProcNumbers = proc_num;
    
    nodl_scalar_ = Nvol;
    nb_unknowns_scal_ = nb_u;
    ProcSharingRows = MatchingProc;
    SharingRowNumbers = MatchingDofNumber;
    
    comm_ = &comm;    
  }
  
  
  //! inits pointers with those of A
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Init(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    OverlapRowNumbers = A.OverlapRowNumbers;
    OverlapProcNumbers = A.OverlapProcNumbers;
    GlobalRowNumbers = A.GlobalRowNumbers;
    ProcSharingRows = A.ProcSharingRows;
    SharingRowNumbers = A.SharingRowNumbers;
    nodl_scalar_ = A.nodl_scalar_;
    nb_unknowns_scal_ = A.nb_unknowns_scal_;
    nglob_ = A.nglob_;
    comm_ = A.comm_;
  }
  
  
  //! changing the size of the local matrix, previous values are lost
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Reallocate(int m, int n)
  {
    // previous values are erased for simplicity
    Clear();
    
    Matrix<T, Prop, Storage, Allocator>::Reallocate(m, n);
    
    dist_col.Reallocate(m);
    dist_row.Reallocate(n);
    proc_col.Reallocate(m);
    proc_row.Reallocate(n);
  }
  
  
  //! changing the size of the local matrix, previous values are kept
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Resize(int m, int n)
  {
    if (this->m_ != m)
      {
        dist_col.Resize(m);
        proc_col.Resize(m);
      }

    if (this->n_ != n)
      {
        dist_row.Resize(n);
        proc_row.Resize(n);
      }
    
    Matrix<T, Prop, Storage, Allocator>::Resize(m, n);    
  }
  
  
  //! matrix is cleared
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>::Clear()
  {
    Matrix<T, Prop, Storage, Allocator>::Clear();
    
    dist_col.Clear();
    proc_col.Clear();
    dist_row.Clear();
    proc_row.Clear();
    
    global_row_to_recv.Clear(); global_col_to_recv.Clear();
    ptr_global_row_to_recv.Clear(); ptr_global_col_to_recv.Clear();
    local_row_to_send.Clear(); local_col_to_send.Clear();
    proc_col_to_recv.Clear(); proc_col_to_send.Clear();
    proc_row_to_recv.Clear(); proc_row_to_send.Clear();
    size_max_distant_row = 0;
    size_max_distant_col = 0;
    local_number_distant_values = false;
  }
  
  
  //! equality *this = X 
  template<class T, class Prop, class Storage, class Allocator>
  inline DistributedMatrix<T, Prop, Storage, Allocator>&
  DistributedMatrix<T, Prop, Storage, Allocator>::
  operator=(const DistributedMatrix<T, Prop, Storage, Allocator>& X)
  {
    Seldon::Copy(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(X),
		 *this);
    
    dist_col = X.dist_col;
    dist_row = X.dist_row;

    proc_col = X.proc_col;
    proc_row = X.proc_row;
    
    GlobalRowNumbers = X.GlobalRowNumbers;
    OverlapProcNumbers = X.OverlapProcNumbers;
    OverlapRowNumbers = X.OverlapRowNumbers;

    ProcSharingRows = X.ProcSharingRows;
    SharingRowNumbers = X.SharingRowNumbers;
    
    nodl_scalar_ = X.nodl_scalar_;
    nb_unknowns_scal_ = X.nb_unknowns_scal_;
    nglob_ = X.nglob_;
    
    comm_ = X.comm_;
    
    global_row_to_recv = X.global_row_to_recv;
    global_col_to_recv = X.global_col_to_recv;
    ptr_global_row_to_recv = X.ptr_global_row_to_recv;
    ptr_global_col_to_recv = X.ptr_global_col_to_recv;
    
    local_row_to_send = X.local_row_to_send;
    local_col_to_send = X.local_col_to_send;
    proc_col_to_recv = X.proc_col_to_recv;
    proc_col_to_send = X.proc_col_to_send;
    proc_row_to_recv = X.proc_row_to_recv;
    proc_row_to_send = X.proc_row_to_send;
    local_number_distant_values = X.local_number_distant_values;
    
    size_max_distant_row = X.size_max_distant_row;
    size_max_distant_col = X.size_max_distant_col;
    
    return *this;
  }
  
  
  //! returns the global number of rows
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetGlobalM() const
  {
    return nglob_;
  }
  
  
  //! returns local to global numbering
  template<class T, class Prop, class Storage, class Allocator>
  inline const IVect& DistributedMatrix<T, Prop, Storage, Allocator>::
  GetGlobalRowNumber() const
  {
    if (this->GlobalRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *GlobalRowNumbers;
  }
  
  
  //! returns local to global numbering
  template<class T, class Prop, class Storage, class Allocator>
  inline IVect& DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetGlobalRowNumber()
  {
    if (this->GlobalRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *GlobalRowNumbers;
  }


  //! returns rows already counted in another processor
  template<class T, class Prop, class Storage, class Allocator>
  inline const IVect& DistributedMatrix<T, Prop, Storage, Allocator>::
  GetOverlapRowNumber() const
  {
    if (this->OverlapRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapRowNumbers;
  }
  

  //! returns rows already counted in another processor
  template<class T, class Prop, class Storage, class Allocator>
  inline IVect& DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetOverlapRowNumber()
  {
    if (this->OverlapRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapRowNumbers;
  }

  
  //! returns processor numbers of the original rows
  template<class T, class Prop, class Storage, class Allocator>
  inline const IVect& DistributedMatrix<T, Prop, Storage, Allocator>::
  GetOverlapProcNumber() const
  {
    if (this->OverlapProcNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapProcNumbers;
  }
  

  //! returns processor numbers of the original rows
  template<class T, class Prop, class Storage, class Allocator>
  inline IVect& DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetOverlapProcNumber()
  {
    if (this->OverlapProcNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *OverlapProcNumbers;
  }


  //! returns the number of scalar unknowns
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetNodlScalar() const
  {
    return nodl_scalar_;
  }
  

  //! returns the number of scalar unknowns
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetNbScalarUnknowns() const
  {
    return nb_unknowns_scal_;
  }

  
  //! returns processor numbers for each set of shared rows
  template<class T, class Prop, class Storage, class Allocator>
  inline IVect& DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetProcessorSharingRows()
  {
    if (this->ProcSharingRows == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *ProcSharingRows;
  }
  
  
  //! returns row numbers for each set of shared rows
  template<class T, class Prop, class Storage, class Allocator>
  inline Vector<IVect>& DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetSharingRowNumbers()
  {
    if (this->SharingRowNumbers == NULL)
      {
        cout << "You should call Init of DistributedMatrix" << endl;
        abort();
      }    

    return *SharingRowNumbers;
  }
  
  
  //! adding a non-zero entry, between a local row/column
  //! and a non-local row/column
  /*!
    This method is used only internally
    \param[in] dist_col array to which entry is added
    \param[in] proc_col processors associated with dist_col
    \param[in] jglob global row/column number
    \param[in] proc2 distant processor containing the global row/column
    \param[in] val value of the non-zero entry
  */
  template<class T, class Allocator>
  inline void AddDistantValue(Vector<T, VectSparse, Allocator>& dist_col,
                              IVect& proc_col,
			      int jglob, int proc2, const T& val)
  {
    int pos = 0;
    int size_row = dist_col.GetM();
    while ((pos < size_row) && (dist_col.Index(pos) < jglob))
      pos++;
    
    if ((pos < size_row) && (dist_col.Index(pos) == jglob))
      {
        // already existing entry
        dist_col.Value(pos) += val;
      }
    else
      {
        // new entry
        Vector<T> value(size_row);
        IVect index(size_row), proc(size_row);
        for (int k = 0; k < size_row; k++)
          {
            index(k) = dist_col.Index(k);
            value(k) = dist_col.Value(k);
            proc(k) = proc_col(k);
          }
        
        dist_col.Reallocate(size_row+1);
        proc_col.Reallocate(size_row+1);
        for (int k = 0; k < pos; k++)
          {
            dist_col.Index(k) = index(k);
            dist_col.Value(k) = value(k);
            proc_col(k) = proc(k);
          }
        
        dist_col.Index(pos) = jglob;
        dist_col.Value(pos) = val;
        proc_col(pos) = proc2;
        for (int k = pos+1; k <= size_row; k++)
          {
            dist_col.Index(k) = index(k-1);
            dist_col.Value(k) = value(k-1);
            proc_col(k) = proc(k-1);
          }
        
      }
  }
  
  
  //! adding a non-zero entry, between a local column
  //! and a non-local column
  /*!
    We assume here that the row number is local to the current processor
    \param[in] i local row number
    \param[in] jglob global column number
    \param[in] proc2 distant processor containing the global column
    \param[in] val value of the non-zero entry
  */
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>::
  AddDistantInteraction(int i, int jglob, int proc2, const entry_type& val)
  {
    if (local_number_distant_values)
      SwitchToGlobalNumbers();
    
    AddDistantValue(dist_col(i), proc_col(i), jglob, proc2, val);
  }
  
  
  //! adding a non-zero entry, between a local row and a non-local row
  /*!
    We assume here that the column number is local to the current processor
    \param[in] iglob global row number
    \param[in] j local column number
    \param[in] proc2 distant processor containing the global row
    \param[in] val value of the non-zero entry
  */
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AddRowDistantInteraction(int iglob, int j,
			     int proc2, const entry_type& val)
  {
    if (local_number_distant_values)
      SwitchToGlobalNumbers();
    
    AddDistantValue(dist_row(j), proc_row(j), iglob, proc2, val);
  }
  
  
  //! returns the maximum number of values in dist_col to exchange
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetMaxDataSizeDistantCol() const
  {
    return size_max_distant_col;
  }
  
  
  //! returns the maximum number of values in dist_row to exchange
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetMaxDataSizeDistantRow() const
  { 
    return size_max_distant_row;
  }
  
  
  //! returns true if the matrix is ready to perform a matrix-vector product
  template<class T, class Prop, class Storage, class Allocator>
  inline bool DistributedMatrix<T, Prop, Storage, Allocator>
  ::IsReadyForMltAdd() const
  {
    return local_number_distant_values;
  }

  
  //! grouping all the local rows of the matrix into CSR form
  /*!
    Creation of a sparse matrix B, that will contain the local rows
    of the current matrix. The column numbers are global.
    Values placed on non-local rows are ignored
  */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Allocator0>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetDistributedRows(Matrix<T0, General, ArrayRowSparse,
		       Allocator0>& B) const
  {    
    int m = this->GetM();
    int n = this->GetGlobalM();
    
    Copy(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), B);
    
    // now, we are using global numbers
    // and removing lower part of the matrix if symmetric
    B.Resize(m, n);
    const IVect& RowNumber = this->GetGlobalRowNumber();
    for (int i = 0; i < m; i++)
      {
	int size_row = B.GetRowSize(i);
        size_row += dist_col(i).GetM();
	IVect index(size_row);
	Vector<T0> value(size_row);
	int nb = 0;
	int num_row = RowNumber(i);
	if (IsSymmetricMatrix(*this))
          {
            // local values
            for (int j = 0; j < B.GetRowSize(i); j++)
              {
                int num_col = RowNumber(B.Index(i, j));
                if (num_row <= num_col)
                  {
                    index(nb) = num_col;
                    value(nb) = B.Value(i, j);
                    nb++;
                  }
              }
            
            // distant values
            for (int j = 0; j < dist_col(i).GetM(); j++)
              {
                int num_col = dist_col(i).Index(j);
                if (this->local_number_distant_values)
                  num_col = global_col_to_recv(num_col);
                
                if (num_row <= num_col)
                  {
                    index(nb) = num_col;
                    value(nb) = dist_col(i).Value(j);
                    nb++;
                  }
              }
          }
        else
          {
            // local values
            for (int j = 0; j < B.GetRowSize(i); j++)
              {
                int num_col = RowNumber(B.Index(i, j));
                index(nb) = num_col;
                value(nb) = B.Value(i, j);
                nb++;
              }
            
            // distant values
            for (int j = 0; j < dist_col(i).GetM(); j++)
              {
                int num_col = dist_col(i).Index(j);
                if (this->local_number_distant_values)
                  num_col = global_col_to_recv(num_col);
                
                index(nb) = num_col;
                value(nb) = dist_col(i).Value(j);
                nb++;
              }
          }
        
        Seldon::Assemble(nb, index, value);
        
	B.ReallocateRow(i, nb);
        for (int j = 0; j < nb; j++)
          {
            B.Index(i, j) = index(j);
            B.Value(i, j) = value(j);
          }
      }
  }
  
  
  //! grouping all the local columns of the matrix into CSC form
  /*!
    Creation of a sparse matrix B, that will contain the local columns
    of the current matrix. The row numbers are global.
    Values placed on non-local columns are ignored
  */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Allocator0>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>::
  GetDistributedColumns(Matrix<T0, General, ArrayColSparse, Allocator0>& B,
                        bool sym_pattern) const
  {
    int m = this->GetGlobalM();
    int n = this->GetN();
    
    // conversion to CSC format of local part and symmetrisation of pattern
    IVect Ptr, Ind; Vector<T0> Val;
    Prop sym;
    ConvertToCSC(*this, sym, Ptr, Ind, Val, sym_pattern);
    
    // for row numbers, we put global numbers and we add some distant entries
    // (i.e entries with local columns,
    // and null values by symmetry of local rows )
    B.Clear(); B.Reallocate(m, n);
    const IVect& RowNumber = this->GetGlobalRowNumber();
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
        size_col += dist_row(i).GetM();
        if (sym_pattern)
          size_col += dist_col(i).GetM();
        
	IVect index(size_col);
	Vector<T0> value(size_col);
	int nb = 0;
        // local values
        for (int j = Ptr(i); j < Ptr(i+1); j++)
          {
            index(nb) = RowNumber(Ind(j));
            value(nb) = Val(j);
            nb++;
          }
        
        // distant values
        for (int j = 0; j < dist_row(i).GetM(); j++)
          {
            index(nb) = dist_row(i).Index(j);
            if (this->local_number_distant_values)    
              index(nb) = global_row_to_recv(index(nb));
            
            value(nb) = dist_row(i).Value(j);
            nb++;
          }
        
        // values due to symmetrisation of pattern
        if (sym_pattern)
          for (int j = 0; j < dist_col(i).GetM(); j++)
            {
              index(nb) = dist_col(i).Index(j);
              if (this->local_number_distant_values)    
                index(nb) = global_col_to_recv(index(nb));
            
              value(nb) = 0;
              nb++;
            }
        
        Seldon::Assemble(nb, index, value);
        
	B.ReallocateColumn(i, nb);
        for (int j = 0; j < nb; j++)
          {
            B.Index(i, j) = index(j);
            B.Value(i, j) = value(j);
          }
      }
  }
   
}

#define SELDON_FILE_DISTRIBUTED_MATRIX_INLINE_CXX
#endif


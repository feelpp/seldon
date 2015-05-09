// Copyright (C) 2014 INRIA
// Author(s): Marc Duruflé
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
  
  
  //! adds non-zero entries in a row of the matrix
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AddInteractionRow(int i, int num_val, const Vector<int>& col,
		      const Vector<entry_type>& val)
  {
    Matrix<T, Prop, Storage, Allocator>::
      AddInteractionRow(i, num_val, col, val);
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


#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ApplySor(Vector<T>& x, const Vector<T>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(*this, x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ApplySor(const class_SeldonTrans& trans, Vector<T>& x, const Vector<T>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans, *this, x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAddComplex(alpha, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAddComplex(alpha, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAddComplex(alpha, trans, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAddComplex(alpha, trans, *this, x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    MltComplex(*this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    MltComplex(*this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    MltComplex(trans, *this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    MltComplex(trans, *this, x, y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltComplex(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y)
  {
    Mlt(A, X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltComplex(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y)
  {
    Mlt(A, X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltComplex(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y)
  {
    Mlt(A, X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltComplex(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y)
  {
    throw WrongArgument("MltComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltComplex(const SeldonTranspose& trans,
		  const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y)
  {
    Mlt(trans, A, X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltComplex(const SeldonTranspose& trans,
		  const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y)
  {
    Mlt(trans, A, X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltComplex(const SeldonTranspose& trans,
		  const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y)
  {
    Mlt(trans, A, X, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltComplex(const SeldonTranspose& trans,
		  const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y)
  {
    throw WrongArgument("MltComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAddComplex(const T0& alpha,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAdd(alpha, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAddComplex(const complex<T0>& alpha,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAddComplex(const complex<T0>& alpha,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAddComplex(const T0& alpha,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltAddComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAddComplex(const T0& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAdd(alpha, trans, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAddComplex(const complex<T0>& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, trans, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAddComplex(const complex<T0>& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAdd(alpha, trans, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAddComplex(const T0& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltAddComplex", "Incompatible matrix-vector product");			
  }

#endif
  
}

#define SELDON_FILE_DISTRIBUTED_MATRIX_INLINE_CXX
#endif


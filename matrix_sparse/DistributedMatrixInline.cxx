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
  
  
  //! returns the global number of rows
  template<class T, class Prop, class Storage, class Allocator>
  inline int DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetGlobalM() const
  {
    return nglob_;
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
  AddDistantInteraction(int i, int jglob, int proc2, const T& val)
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
			     int proc2, const T& val)
  {
    if (local_number_distant_values)
      SwitchToGlobalNumbers();
    
    AddDistantValue(dist_row(j), proc_row(j), iglob, proc2, val);
  }
  
  
  //! adds non-zero entries in a row of the matrix
  template<class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AddInteractionRow(int i, int num_val, const Vector<int>& col,
		      const Vector<T>& val)
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
    MltAdd(alpha, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Treal& alpha, const SeldonTranspose& trans,
		 const Vector<Treal>& x,
		 const Treal& beta, Vector<Treal>& y) const
  {
    MltAdd(alpha, trans, *this, x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddVector(const Tcplx& alpha, const SeldonTranspose& trans,
		 const Vector<Tcplx>& x,
		 const Tcplx& beta, Vector<Tcplx>& y) const
  {
    MltAdd(alpha, trans, *this, x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(*this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(*this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Treal>& x, Vector<Treal>& y) const
  {
    Mlt(trans, *this, x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void DistributedMatrix<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<Tcplx>& x, Vector<Tcplx>& y) const
  {
    Mlt(trans, *this, x, y);
  }
#endif


  /****************************************
   * Mlt, MltAdd for distributed matrices *
   ****************************************/


  template<class T, class Prop, class Storage, class Allocator>
  inline void Mlt(const T& alpha,
                  DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }

  template<class T, class Prop, class Storage, class Allocator>
  inline void Mlt(const T& alpha,
                  DistributedMatrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    MltScalar(alpha, A);
  }


  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
		  DistributedMatrix<T, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const complex<T>& alpha,
		  const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
		  DistributedMatrix<complex<T>, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }


  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const DistributedMatrix<complex<T>, Prop1, Storage1, Allocator1>& A,
		  DistributedMatrix<complex<T>, Prop2, Storage2, Allocator2>& B)
  {
    AddMatrix(alpha, A, B);
  }

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  inline void Add(const T& alpha,
		  const DistributedMatrix<complex<T>, Prop1, Storage1, Allocator1>& A,
		  DistributedMatrix<T, Prop2, Storage2, Allocator2>& B)
  {
    throw WrongArgument("Add", "incompatible types");    
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    MltVector(A, X, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y,
		  bool assemble)
  {
    MltVector(A, X, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    MltVector(trans, A, X, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		  const Vector<complex<T0> >& X, Vector<complex<T0> >& Y,
		  bool assemble)
  {
    MltVector(trans, A, X, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void Mlt(const SeldonTranspose& trans,
		  const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		  const Vector<T0>& X, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAddVector(alpha, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const complex<T0>& alpha,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAddVector(alpha, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltAddComplex", "Incompatible matrix-vector product");			
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    MltAddVector(alpha, trans, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		     const Vector<complex<T0> >& X, const complex<T0>& beta,
		     Vector<complex<T0> >& Y, bool assemble)
  {
    MltAddVector(alpha, trans, A, X, beta, Y, assemble);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  inline void MltAdd(const T0& alpha, const SeldonTranspose& trans,
		     const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
		     const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble)
  {
    throw WrongArgument("MltAddComplex", "Incompatible matrix-vector product");			
  }

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const T& omega, int iter, int type_ssor)
  {
    SorVector(M, Y, X, omega, iter, type_ssor);
  }


  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  inline void SOR(const class_SeldonTrans&,
		  const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
		  Vector<T, Storage2, Allocator2>& Y,
		  const Vector<T, Storage1, Allocator1>& X,
		  const T& omega, int iter, int type_ssor)
  {
    SorVector(SeldonTrans, M, Y, X, omega, iter, type_ssor);
  }
  
}

#define SELDON_FILE_DISTRIBUTED_MATRIX_INLINE_CXX
#endif


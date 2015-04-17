// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_MATRIX_BASE_CXX

#include "Matrix_Base.hxx"

namespace Seldon
{

  // Matrix allocator.
  template <class T, class Allocator>
  Allocator Matrix_Base<T, Allocator>::allocator_;


#ifdef SELDON_WITH_VIRTUAL
  //! Applies S.O.R method to solve A x = r.
  template<class T> void VirtualMatrix<T>
  ::ApplySor(Vector<T>& x, const Vector<T>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    throw Undefined("ApplySOR", "Not implemented");
  }
  

  //! Applies S.O.R method to solve A x = r or A^T x = r.
  template<class T> void VirtualMatrix<T>
  ::ApplySor(const class_SeldonTrans&, Vector<T>& x, const Vector<T>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    throw Undefined("ApplySOR", "Not implemented");
  }


  //! Computes y = beta + alpha A x.
  template<class T> void VirtualMatrix<T>
  ::MltAddVector(const T& alpha, const Vector<T>& x,
		 const T& beta, Vector<T>& y) const
  {
    throw Undefined("MltAddVector", "Not implemented");
  }


  //! Computes y = A x
  template<class T> void VirtualMatrix<T>
  ::MltVector(const Vector<T>& x, Vector<T>& y) const
  {
    throw Undefined("MltVector", "Not implemented");
  }
  

  //! Computes y = A^T x
  template<class T> void VirtualMatrix<T>
  ::MltVector(const class_SeldonTrans&, const Vector<T>& x, Vector<T>& y) const
  {
    throw Undefined("MltVector", "Not implemented");
  }
#endif

} // namespace Seldon.

#define SELDON_FILE_MATRIX_BASE_CXX
#endif

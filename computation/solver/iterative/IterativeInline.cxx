// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_ITERATIVE_INLINE_CXX

namespace Seldon
{

  //! Default constructor
  inline Preconditioner_Base::Preconditioner_Base()
  {
  }


  //! Solves M z = r
  /*!
    Identity preconditioner M = I
  */
  template<class Matrix1, class Vector1>
  inline void Preconditioner_Base::Solve(const Matrix1& A, const Vector1& r,
                                         Vector1& z)
  {
    Copy(r,z);
  }


  //! Solves M^t z = r
  /*!
    Identity preconditioner M = I
  */
  template<class Matrix1, class Vector1>
  inline void Preconditioner_Base::
  TransSolve(const Matrix1& A, const Vector1 & r, Vector1 & z)
  {
    Solve(A, r, z);
  }

} // end namespace

#define SELDON_FILE_ITERATIVE_INLINE_CXX
#endif


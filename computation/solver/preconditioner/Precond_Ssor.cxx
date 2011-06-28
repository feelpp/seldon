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


#ifndef SELDON_FILE_PRECOND_SSOR_CXX

namespace Seldon
{
  
  template <class T>
  class SorPreconditioner
  {
  protected :
    bool symmetric_precond; //!< true for Symmetric relaxation
    int nb_iter; //!< number of iterations
    T omega; //!< relaxation parameter

  public :
    SorPreconditioner();
    inline ~SorPreconditioner(){}

    inline void InitSymmetricPreconditioning() { symmetric_precond = true; }
    inline void InitUnSymmetricPreconditioning() { symmetric_precond = false; }
    inline void SetParameterRelaxation(const T& param) { omega = param; }
    inline void SetNumberIterations(int nb_iterations) { nb_iter = nb_iterations; }

    template<class Vector1, class Matrix1>
    void Solve(const Matrix1& A, const Vector1& r, Vector1& z,
	       bool init_guess_null = true);

    template<class Vector1, class Matrix1>
    void TransSolve(const Matrix1& A, const Vector1& r, Vector1& z,
		    bool init_guess_null = true);

  };


  //! Default constructor
  template<class T>
  inline SorPreconditioner<T>::SorPreconditioner()
  {
    nb_iter = 1; omega = T(1);
    symmetric_precond = true;
  }


  //! Solves M z = r
  template<class T> template<class Vector1, class Matrix1>
  inline void SorPreconditioner<T>::
  Solve(const Matrix1& A, const Vector1& r, Vector1& z, bool init_guess_null)
  {
    if (init_guess_null)
      z.Zero();

    if (symmetric_precond)
      Seldon::SOR(A, z, r, omega, nb_iter, 0);
    else
      Seldon::SOR(A, z, r, omega, nb_iter, 2);
  }


  //! Solves M^t z = r
  template<class T> template<class Vector1, class Matrix1>
  inline void SorPreconditioner<T>::
  TransSolve(const Matrix1& A, const Vector1& r,
	     Vector1& z, bool init_guess_null)
  {
    if (init_guess_null)
      z.Zero();
    
    if (symmetric_precond)
      Seldon::SOR(SeldonTrans, A, z, r, omega, nb_iter, 0);
    else
      Seldon::SOR(SeldonTrans, A, z, r, omega, nb_iter, 3);
    
  }

}

#define SELDON_FILE_PRECOND_SSOR_CXX
#endif

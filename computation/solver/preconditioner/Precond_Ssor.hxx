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


#ifndef SELDON_FILE_PRECOND_SSOR_HXX

namespace Seldon
{
  
  template <class T>
  class SorPreconditioner
  {
  protected :
    bool symmetric_precond; //!< true for Symmetric relaxation
    int nb_iter; //!< number of iterations
    typename ClassComplexType<T>::Treal omega; //!< relaxation parameter

  public :
    SorPreconditioner();

    void InitSymmetricPreconditioning();
    void InitUnSymmetricPreconditioning();
    
    void SetParameterRelaxation(const typename ClassComplexType<T>::Treal& param);
    
    void SetNumberIterations(int nb_iterations);

    template<class Vector1, class Matrix1>
    void Solve(const Matrix1& A, const Vector1& r, Vector1& z,
	       bool init_guess_null = true);

    template<class Vector1, class Matrix1>
    void TransSolve(const Matrix1& A, const Vector1& r, Vector1& z,
		    bool init_guess_null = true);

  };

}

#define SELDON_FILE_PRECOND_SSOR_HXX
#endif

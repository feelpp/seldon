// Copyright (C) 2001-2008 Vivien Mallet
//
// This file is part of Seldon library.
// Seldon library provides matrices and vectors structures for
// linear algebra.
//
// Seldon is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Seldon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the Seldon home page:
//     http://spacetown.free.fr/lib/seldon/

#ifndef SELDON_FILE_ITERATIVE_CGNE_CXX

namespace Seldon
{

  //! Solves a linear system using Conjugate Gradient Normal Equation (CGNE)
  /*!
    Solves the unsymmetric linear system A x = b.
    
    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.
    
    \param[in] A  Complex General Matrix
    \param[inout] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Right hand side of the linear system
    \param[in] M Left preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix, class Vector, class Preconditioner>
  int Cgne(Matrix& A, Vector& x, const Vector& b,
	   Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;
    
    typedef typename Vector::value_type Complexe;
    Complexe rho(1), rho_1(0), alpha, beta, delta;
    Vector p(N), q(N), r(N), z(N);
    Titer dp;
    
    // x should be equal to 0
    // see Cg to understand implementation
    // we solve A^t A x = A^t b
    // left-preconditioner is equal to M M^t
    
    // q = A^t b
    Mlt(SeldonTrans, A, b, q);
    // we initialize iter
    int success_init = iter.Init(q);
    if (success_init != 0)
      return iter.ErrorCode();
    
    if (!iter.IsInitGuess_Null())
      {
	// r = A^t b - A^t A x
	Mlt(A, x, p);
	Mlt(SeldonTrans, A, p, r); Mlt(Complexe(-1), r);
	Seldon::Add(Complexe(1), q, r);
      }
    else
      {
	Seldon::Copy(q, r);
	x.Zero();
      }
    
    dp = Norm2(q);
    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (! iter.Finished(dp))
      {
	// Preconditioning
	M.TransSolve(A, r, q);
	M.Solve(A, q, z);
	
	rho = DotProd(r,z);
	if (rho == Complexe(0))
	  {
	    iter.Fail(1, "Cgne breakdown #1");
	    break;
	  }
	
	if (iter.First())
	  Seldon::Copy(z, p);
	else
	  {
	    beta = rho / rho_1;
	    Mlt(beta, p);
	    Seldon::Add(Complexe(1), z, p);
	  }
	
	// instead of q = A*p, we compute q = A^t A *p
	Mlt(A, p, q);
	Mlt(SeldonTrans, A, q, z);
	
	delta = DotProd(p, z);
	if (delta == Complexe(0))
	  {
	    iter.Fail(2, "Cgne breakdown #2");
	    break;
	  }
	alpha = rho / delta;
	
	Seldon::Add(alpha, p, x);
	Seldon::Add(-alpha, z, r);
	
	rho_1 = rho;
	dp = Norm2(r);
	
	// two iterations, because of two multiplications with A
	++iter;
	++iter;
      }
    
    return iter.ErrorCode();
  }
  
  
} // end namespace

#define SELDON_FILE_ITERATIVE_CGNE_CXX
#endif

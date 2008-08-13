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

#ifndef SELDON_FILE_ITERATIVE_CGS_CXX

namespace Seldon 
{
  
  //! Solves linear system using Conjugate Gradient Squared (CGS)
  /*!
    Solves the unsymmetric linear system Ax = b 
    using the Conjugate Gradient Squared method.
    
    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.
        
    See: P. Sonneveld, CGS, a fast Lanczos-type solver for nonsymmetric linear
    systems, SIAM, J.Sci. Statist. Comput., 10(1989), pp. 36-52
    
    \param[in] A Complex General Matrix 
    \param[inout] x Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner   
    \param[in] iter Iteration parameters
  */  
  template <class Titer, class Matrix, class Vector, class Preconditioner>
  int Cgs(Matrix& A, Vector& x, const Vector& b,
	  Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;
    
    typedef typename Vector::value_type Complexe;
    Complexe rho_1, rho_2(0), alpha, beta, delta;
    Vector p(N), phat(N), q(N), qhat(N), vhat(N), u(N), uhat(N),
      r(N), rtilde(N);
    
    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();
    
    // we compute the initial residual r = b - Ax
    Seldon::Copy(b,r);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r);
    else
      x.Zero();
    
    Seldon::Copy(r, rtilde);
    
    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are reached
    while (! iter.Finished(r)) 
      {
	rho_1 = DotProd(rtilde, r);
	
	if (rho_1 == Complexe(0)) 
	  {
	    iter.Fail(1, "Cgs breakdown #1");
	    break;
	  }
	  
	if (iter.First()) 
	  {
	    Seldon::Copy(r, u);
	    Seldon::Copy(u, p);
	  } 
	else 
	  {
	    // u = r + beta*q
	    // p = beta*(beta*p +q) + u  where beta = rho_i/rho_{i-1}
	    beta = rho_1 / rho_2;
	    Seldon::Copy(r, u);
	    Seldon::Add(beta, q, u);
	    Mlt(beta, p);
	    Seldon::Add(Complexe(1), q, p);
	    Mlt(beta, p);
	    Seldon::Add(Complexe(1), u, p);
	  }
	
	// preconditioning phat = M^{-1} p
	M.Solve(A, p, phat);
	
	// matrix vector product vhat = A*phat
	Mlt(A, phat, vhat); ++iter;
	delta = DotProd(rtilde, vhat);
	if (delta == Complexe(0)) 
	  {
	    iter.Fail(2, "Cgs breakdown #2");
	    break;
	  }
	// q = u-alpha*vhat  where alpha = rho_i/(rtilde,vhat)
	alpha = rho_1 /delta; 
	Seldon::Copy(u,q);
	Seldon::Add(-alpha, vhat, q);
	
	//  u =u+q
	Seldon::Add(Complexe(1), q, u);
	M.Solve(A, u, uhat);
	
	Seldon::Add(alpha, uhat, x);
	Mlt(A, uhat, qhat);
	Seldon::Add(-alpha, qhat, r);
	
	rho_2 = rho_1;
	
	++iter;
      }
    
    return iter.ErrorCode();
  }
  
} // end namespace

#define SELDON_FILE_ITERATIVE_CGS_CXX
#endif


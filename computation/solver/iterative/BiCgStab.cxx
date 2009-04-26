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

#ifndef SELDON_FILE_ITERATIVE_BICGSTAB_CXX

namespace Seldon
{

  //! Implements  BiConjugate Gradient Stabilized (BICG-STAB)
  /*!
    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.
    
    See: H. Van der Vorst, Bi-CGSTAB: A fast and smoothly converging variant
    of BiCG for the solution of nonsysmmetric linear systems, SIAM J. Sci.
    Statist. Comput. 13(1992), pp. 631-644
    
    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix, class Vector, class Preconditioner>
  int BiCgStab(Matrix& A, Vector& x, const Vector& b,
	       Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;
    
    typedef typename Vector::value_type Complexe;
    Complexe rho_1, rho_2(0), alpha(0), beta, omega(0), sigma;
    Vector p(N), phat(N), s(N), shat(N), t(N), v(N), r(N), rtilde(N);
    
    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();
    
    // we compute the residual r = b - Ax
    Seldon::Copy(b, r);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r);
    else
      x.Zero();
    
    Seldon::Copy(r, rtilde);
    
    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (! iter.Finished(r))
      {
	
	rho_1 = DotProdConj(rtilde, r);
	if (rho_1 == Complexe(0))
	  {
	    iter.Fail(1, "Bicgstab breakdown #1");
	    break;
	  }
	
	if (iter.First())
	  Seldon::Copy(r, p);
	else
	  {
	    if (omega == Complexe(0))
	      {
		iter.Fail(2, "Bicgstab breakdown #2");
		break;
	      }
	    // p= r + beta*(p-omega*v)
	    // beta = rho_i/rho_{i-1} * alpha/omega
	    beta = (rho_1 / rho_2) * (alpha / omega);
	    Seldon::Add(-omega, v, p);
	    Mlt(beta, p);
	    Seldon::Add(Complexe(1), r, p);
	  }
	// preconditioning phat = M^{-1} p
	M.Solve(A, p, phat);
	
	// product matrix vector  v = A*phat
	Mlt(A, phat, v);
	
	// s=r-alpha*v  where alpha = rho_i / (v,rtilde)
	sigma = DotProdConj(rtilde, v);
	if (sigma == Complexe(0))
	  {
	    iter.Fail(3, "Bicgstab breakdown #3");
	    break;
	  }
	alpha = rho_1 / sigma;
	Seldon::Copy(r, s);
	Seldon::Add(-alpha, v, s);
	
	// we increment iter, bicgstab has two products matrix vector
	++iter;
	if (iter.Finished(s))
	  {
	    // x=x+alpha*phat
	    Seldon::Add(alpha, phat, x);
	    break;
	  }
	
	// preconditioning shat = M^{-1} s
	M.Solve(A, s, shat);
	
	// product matrix vector t = A*shat
	Mlt(A, shat, t);
	
	omega = DotProdConj(t, s) / DotProdConj(t, t);
	
	// new iterate x=x+alpha*phat+omega*shat
	Seldon::Add(alpha, phat, x);
	Seldon::Add(omega, shat, x);
	
	// new residual r=s-omega*t
	Seldon::Copy(s, r);
	Seldon::Add(-omega, t, r);
	
	rho_2 = rho_1;
	
	++iter;
      }
    
    return iter.ErrorCode();
  }
  
} // end namespace

#define SELDON_FILE_ITERATIVE_BICGSTAB_CXX
#endif

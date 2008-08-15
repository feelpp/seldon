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

#ifndef SELDON_FILE_ITERATIVE_BICGCR_CXX

namespace Seldon
{
  
  //! Solves a linear system by using BiCgCr
  /*!
    Solves the symmetric linear system Ax = b using the
    BiConjugate Gradient Conjugate Residual method.
    
    See Iterative Methods for the Solution of Very Large Complex Symmetric Linear Systems of
    Equations in Electrodynamics,
    Markus Clemens, Thomas Weiland
    Fachbereich 18 Elektrische Nachrichtentechnik, Fachgebiet Theorie Elektromagnetischer Felder,
    Technische Hochschule Darmstadt, Schlossgartenstr. 8, D-64289 Darmstadt, Germany
    
    \param[in] A  Complex Symmetric Matrix
    \param[inout] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix, class Vector, class Preconditioner>
  int BiCgcr(Matrix& A, Vector& x, const Vector& b,
	     Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;
    
    typedef typename Vector::value_type Complexe;
    Complexe rho, mu, alpha, beta, tau;
    Vector v(N), w(N), s(N), z(N), p(N), a(N);
    v.Zero(); w.Zero(); s.Zero(); z.Zero();  p.Zero(); a.Zero();
    
    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();
    
    // we compute the residual v = b - Ax
    Seldon::Copy(b, v);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), v);
    else
      x.Zero();
    
    iter.SetNumberIteration(0);
    // s = M*v   p = s
    M.Solve(A, v, s); Seldon::Copy(s, p);
    // a = A*p   w = M*a
    Mlt(A, p, a); M.Solve(A, a, w);
    // we made one product matrix vector
    ++iter;
    
    while (! iter.Finished(v))
      {
	rho = DotProd(w, v);
	mu = DotProd(w, a);
	
	// new iterate x = x + alpha*p0  where alpha = (r1,r0)/(p1,p1)
	if (mu == Complexe(0))
	  {
	    iter.Fail(1, "Bicgcr breakdown #1");
	    break;
	  }
	alpha = rho/mu;
	Seldon::Add(alpha, p, x);
	
	// new residual r0 = r0 - alpha * p1
	// r1 = r1 - alpha*p2
	Seldon::Add(-alpha, a, v);
	Seldon::Add(-alpha, w, s);
	
	Mlt(A, s, z);
	tau = DotProd(w, z);
	
	if (tau == Complexe(0))
	  {
	    iter.Fail(2, "Bicgcr breakdown #2");
	    break;
	  }
	
	beta = tau/mu;
	
	Mlt(Complexe(-beta), p);
	Seldon::Add(Complexe(1), s, p);
	Mlt(Complexe(-beta), a);
	Seldon::Add(Complexe(1), z, a);
	
	M.Solve(A, a, w);

	++iter;
      }
    
    return iter.ErrorCode();
  }
  
} // end namespace

#define SELDON_FILE_ITERATIVE_BICGCR_CXX
#endif

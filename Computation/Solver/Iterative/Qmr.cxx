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

#ifndef SELDON_FILE_ITERATIVE_QMR_CXX

namespace Seldon
{
    
  //! Solves a linear system by using Quasi-Minimal Residual (QMR)
  /*!
    Solves the unsymmetric linear system Ax = b using the
    Quasi-Minimal Residual method.
    
    See: R. W. Freund and N. M. Nachtigal, A quasi-minimal residual method for
    non-Hermitian linear systems, Numerical Math., 60(1991), pp. 315-339
    
    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix, class Vector, class Preconditioner>
  int Qmr(Matrix& A, Vector& x, const Vector& b,
	  Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;
    
    typedef typename Vector::value_type Complexe;
    Complexe delta, ep(0), beta;
    Titer  rho, rho_1, xi;
    Complexe theta_1, gamma_1;
    Complexe theta(0), gamma(1), eta(-1);
    
    Vector r(N), y(N), z_tld(N); r.Zero();
    Vector v(N), w(N), p_tld(N);
    Vector p(N), q(N), d(N), s(N);
    
    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();
    
    // r = b - Ax
    Seldon::Copy(b, r);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r);
    else
      x.Zero();
    
    Seldon::Copy(r, v);
    
    M.Solve(A, v, y);
    rho = Norm2(y);
    
    Seldon::Copy(r, w);
    xi = Norm2(w);
    
    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are reached
    while (! iter.Finished(r))
      {
	if (rho == Titer(0))
	  {
	    iter.Fail(1, "Qmr breakdown #1");
	    break;
	  }
	if (xi == Titer(0))
	  {
	    iter.Fail(2, "Qmr breakdown #2");
	    break;
	  }
	  
	// v = v / rho
	// y = y / rho
	Mlt(Complexe(1./rho), v);
	Mlt(Complexe(1./rho), y);
	
	// w = w / xi
	Mlt(Complexe(1./xi), w);
	
	delta = DotProd(w, y);
	if (delta == Complexe(0))
	  {
	    iter.Fail(3, "Qmr breakdown #3");
	    break;
	  }
	
	M.TransSolve(A, w, z_tld);
	
	if (iter.First())
	  {
	    Seldon::Copy(y, p);
	    Seldon::Copy(z_tld, q);
	  }
	else
	  {
	    // p = y - (xi delta / ep) p
	    // q = z_tld - (rho delta / ep) q
	    Mlt(Complexe(-(xi  * delta / ep)), p);
	    Seldon::Add(Complexe(1), y, p);
	    Mlt(Complexe(-(rho  * delta / ep)), q);
	    Seldon::Add(Complexe(1), z_tld, q);
	  }
	
	// product matrix vector p_tld = A*p
	Mlt(A, p, p_tld);
	
	ep = DotProd(q, p_tld);
	if (ep == Complexe(0))
	  {
	    iter.Fail(4, "Qmr breakdown #4");
	    break;
	  }
	
	beta = ep / delta;
	if (beta == Complexe(0))
	  {
	    iter.Fail(5, "Qmr breakdown #5");
	    break;
	  }
	  
	// v = -beta v + p_tld
	Mlt(Complexe(-beta), v);
	Seldon::Add(Complexe(1), p_tld, v);
	M.Solve(A, v, y);
	
	rho_1 = rho;
	rho = Norm2(y);
	
	// product matrix vector z_tld = A q
	Mlt(SeldonTrans, A, q, z_tld);
	// w = z_tld - beta*w
	Mlt(Complexe(-beta), w); Seldon::Add(Complexe(1), z_tld, w);
	
	xi = Norm2(w);
	
	gamma_1 = gamma;
	theta_1 = theta;
	
	++iter;
	theta = rho / (gamma_1 * beta);
	gamma = Complexe(1) / sqrt(1.0 + theta * theta);
	
	if (gamma == Titer(0))
	  {
	    iter.Fail(6, "Qmr breakdown #6");
	    break;
	  }
	
	eta = -eta * rho_1 * gamma * gamma / (beta * gamma_1 * gamma_1);
	
	if (iter.First())
	  {
	    Seldon::Copy(p, d);
	    Mlt(eta, d);
	    Seldon::Copy(p_tld, s);
	    Mlt(eta, s);
	  }
	else
	  {
	    Complexe tmp = (theta_1 * theta_1 * gamma * gamma);
	    Mlt(tmp, d);
	    Seldon::Add(eta, p, d);
	    Mlt(tmp, s);
	    Seldon::Add(eta, p_tld, s);
	  }
	Seldon::Add(Complexe(1), d, x);
	Seldon::Add(-Complexe(1), s, r);
	
	++iter;
      }
    
    return iter.ErrorCode();
  }
    
} // end namespace

#define SELDON_FILE_ITERATIVE_QMR_CXX
#endif

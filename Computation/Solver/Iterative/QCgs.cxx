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

#ifndef SELDON_FILE_ITERATIVE_QCGS_CXX

namespace Seldon
{
  
  //! Solves linear system using Quasi-minimized Conjugate Gradient Squared
  /*!
    Solves the unsymmetric linear system Ax = b
    using the Quasi-minimized Conjugate Gradient Squared method.
    
    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.
    
    see A Comparative Study of Preconditioned Lanczos Methods
    for Nonsymmetric Linear Systems
    C. H. Tong, Sandia Report
    
    \param[in] A Complex General Matrix
    \param[inout] x Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] iter Iteration parameters
  */
  template <class Titer, class Matrix, class Vector, class Preconditioner>
  int QCgs(Matrix& A, Vector& x, const Vector& b,
	   Preconditioner& M, Iteration<Titer> & iter)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;
    
    typedef typename Vector::value_type Complexe;
    Complexe rho_1, rho_2, mu,nu,alpha,beta,sigma,delta;
    Vector p(N), q(N), r(N), rtilde(N), u(N),
      phat(N), r_qcgs(N), x_qcgs(N), v(N);
    
    // we initialize iter
    int success_init = iter.Init(b);
    if (success_init != 0)
      return iter.ErrorCode();
    
    // we compute the residual r = b - Ax
    Seldon::Copy(b,r);
    if (!iter.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), r);
    else
      x.Zero();
    
    Seldon::Copy(r, rtilde);
    rho_1 = Complexe(1);
    q.Zero(); p.Zero();
    Seldon::Copy(r, r_qcgs); Seldon::Copy(x, x_qcgs);
    Seldon::Matrix<Complexe, Symmetric, RowSymPacked> bt_b(2,2),bt_b_m1(2,2);
    Seldon::Vector<Complexe> bt_rn(2);
    
    iter.SetNumberIteration(0);
    // Loop until the stopping criteria are reached
    while (! iter.Finished(r_qcgs))
      {
	rho_2 = DotProd(rtilde, r);
	
	if (rho_1 == Complexe(0))
	  {
	    iter.Fail(1, "QCgs breakdown #1");
	    break;
	  }
	beta = rho_2/rho_1;
	
	// u = r + beta*q
	// p = beta*(beta*p +q) + u  where beta = rho_i/rho_{i-1}
	Seldon::Copy(r, u);
	Seldon::Add(beta, q, u);
	Mlt(beta, p);
	Seldon::Add(Complexe(1), q, p);
	Mlt(beta, p);
	Seldon::Add(Complexe(1), u, p);
	
	// preconditioning phat = M^{-1} p
	M.Solve(A, p, phat);
	
	// product matrix vector vhat = A*phat
	Mlt(A, phat, v); ++iter;
	sigma = DotProd(rtilde, v);
	if (sigma == Complexe(0))
	  {
	    iter.Fail(2, "Qcgs breakdown #2");
	    break;
	  }
	// q = u-alpha*v  where alpha = rho_i/(rtilde,v)
	alpha = rho_2 /sigma;
	Seldon::Copy(u, q);
	Seldon::Add(-alpha, v, q);
	
	// u = u +q
	Seldon::Add(Complexe(1), q, u);
	// x = x + alpha u
	Seldon::Add(alpha, u, x);
	// preconditioning phat = M^{-1} u
	M.Solve(A, u, phat);
	// product matrix vector q = A*phat
	Mlt(A, phat, u);
	
	// r = r - alpha u
	Seldon::Add(-alpha, u, r);
	// u = r_qcgs - r_n
	Seldon::Copy(r_qcgs, u);
	Seldon::Add(-Complexe(1), r, u);
	bt_b(0,0) = DotProd(u,u);
	bt_b(1,1) = DotProd(v,v);
	bt_b(1,0) = DotProd(u,v);
	
	// we compute inverse of bt_b
	delta = bt_b(0,0)*bt_b(1,1) - bt_b(1,0)*bt_b(0,1);
	if (delta == Complexe(0))
	  {
	    iter.Fail(3, "Qcgs breakdown #3");
	    break;
	  }
	bt_b_m1(0,0) = bt_b(1,1)/delta;
	bt_b_m1(1,1) = bt_b(0,0)/delta;
	bt_b_m1(1,0) = -bt_b(1,0)/delta;
	
	bt_rn(0) = -DotProd(u, r); bt_rn(1) = -DotProd(v, r);
	mu = bt_b_m1(0,0)*bt_rn(0)+bt_b_m1(0,1)*bt_rn(1);
	nu = bt_b_m1(1,0)*bt_rn(0)+bt_b_m1(1,1)*bt_rn(1);
	
	// smoothing step
	// r_qcgs = r + mu (r_qcgs - r) + nu v
	Seldon::Copy(r, r_qcgs); Seldon::Add(mu, u, r_qcgs);
	Seldon::Add(nu, v, r_qcgs);
	// x_qcgs = x + mu (x_qcgs - x) - nu p
	Mlt(mu,x_qcgs);
	Seldon::Add(Complexe(1)-mu, x, x_qcgs);
	Seldon::Add(-nu, p, x_qcgs);
	
	rho_1 = rho_2;
	
	++iter;
      }
    
    M.Solve(A, x_qcgs, x);
    return iter.ErrorCode();
  }

} // end namespace

#define SELDON_FILE_ITERATIVE_QCGS_CXX
#endif


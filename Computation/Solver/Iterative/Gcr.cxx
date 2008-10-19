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

#ifndef SELDON_FILE_ITERATIVE_GCR_CXX

namespace Seldon
{

  //! Solves a linear system by using Generalized Conjugate Residual (GCR)
  /*!
    Solves the linear system Ax = b with restarted Preconditioned
    Generalized Conjugate Residual Algorithm.
    
    return value of 0 indicates convergence within the
    maximum number of iterations (determined by the iter object).
    return value of 1 indicates a failure to converge.
    
    See: Y. Saad, Iterative Methods for Sparse Linear System, PWS Publishing
    Company, 1996

    \param[in] A  Complex General Matrix
    \param[in,out] x  Vector on input it is the initial guess
    on output it is the solution
    \param[in] b  Vector right hand side of the linear system
    \param[in] M Right preconditioner
    \param[in] outer Iteration parameters
  */
  template <class Titer, class Matrix, class Vector, class Preconditioner>
  int Gcr(Matrix& A, Vector& x, const Vector& b,
	  Preconditioner& M, Iteration<Titer> & outer)
  {
    const int N = A.GetM();
    if (N <= 0)
      return 0;
    
    typedef typename Vector::value_type Complexe;
    int m = outer.GetRestart();
    // we initialize outer
    int success_init = outer.Init(b);
    if (success_init != 0)
      return outer.ErrorCode();
    
    Seldon::Vector<Vector,Vect_Full,NewAlloc<Vector> > p(m+1), w(m+1);
    
    Vector beta(m+1);
    
    Vector r(N), q(N), u(N);
    
    for (int i = 0; i < (m+1); i++)
      {
	p(i).Reallocate(N);
	p(i).Zero();
	w(i).Reallocate(N);
	w(i).Zero();
      }

    // we compute initial residual
    Seldon::Copy(b,u);
    if (!outer.IsInitGuess_Null())
      MltAdd(Complexe(-1), A, x, Complexe(1), u);
    else
      x.Zero();
    
    M.Solve(A, u, r);
    
    Complexe alpha,delta;
    
    Titer normr = Norm2(r);
    outer.SetNumberIteration(0);
    // Loop until the stopping criteria are satisfied
    while (! outer.Finished(r))
      {
	// m is the maximum number of inner iterations
	Iteration<Titer> inner(outer);
	inner.SetNumberIteration(outer.GetNumberIteration());
	inner.SetMaxNumberIteration(outer.GetNumberIteration()+m);
	Seldon::Copy(r, p(0));
	Mlt(Titer(1)/normr, p(0));
	
	int j = 0;
	
	while (! inner.Finished(r) )
	  {
	    // product matrix vector u=A*p(j)
	    Mlt(A, p(j), u);
	    
	    // preconditioning
	    M.Solve(A, u, w(j));
	    
	    beta(j) = DotProdConj(w(j), w(j));
	    if (beta(j) == Complexe(0))
	      {
		outer.Fail(1, "Gcr breakdown #1");
		break;
	      }
	    
	    // new iterate x = x + alpha*p(j) new residual r = r - alpha*w(j)
	    // where alpha = (conj(r_j),A*p_j)/(A*p_j,A*p_j)
	    alpha = DotProdConj(w(j), r) / beta(j);
	    Seldon::Add(alpha, p(j), x);
	    Seldon::Add(-alpha, w(j), r);
	    
	    ++inner;
	    ++outer;
	    // product Matrix vector u = A*r
	    Mlt(A, r, u);
	    M.Solve(A, u, q);
	    
	    Seldon::Copy(r, p(j+1));
	    // we compute direction p(j+1) = r(j+1) +
	    // \sum_{i=0..j} ( -(A*r_j+1,A*p_i)/(A*p_i,A*p_i) p(i))
	    for (int i = 0; i <= j; i++)
	      {
		delta = -DotProdConj(w(i), q)/beta(i);
		Seldon::Add(delta, p(i), p(j+1));
	      }
	    
	    ++inner;
	    ++outer;
	    ++j;
	  }
	normr = Norm2(r);
      }
    
    return outer.ErrorCode();
  }
  
} // end namespace

#define SELDON_FILE_ITERATIVE_GCR_CXX
#endif

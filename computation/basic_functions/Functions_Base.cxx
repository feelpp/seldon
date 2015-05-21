#ifndef SELDON_FILE_FUNCTIONS_BASE_CXX

// in this file, we put functions such as Add, Mlt, MltAdd
// with a reduced number of templates in order
// to forbid undesired instantiations of generic functions

/*
  Functions defined in this file:
  
  M X -> Y
  Mlt(M, X, Y)

  alpha M X -> Y
  Mlt(alpha, M, X, Y)

  M X -> Y
  M^T X -> Y
  Mlt(Trans, M, X, Y)

*/

namespace Seldon
{

  /**********************************
   * Functions in Functions_MatVect *
   **********************************/
  
  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  void Mlt(const SeldonTranspose& Trans,
	   const Matrix<T, Prop1, Storage1, Allocator1>& M,
	   const Vector<T, Storage2, Allocator2>& X,
	   Vector<T, Storage3, Allocator3>& Y)
  {
    if (Trans.NoTrans())
      MltVector(M, X, Y);
    else if (Trans.Trans())
      MltVector(SeldonTrans, M, X, Y);
    else
      MltVector(SeldonConjTrans, M, X, Y);
  }


  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  void Mlt(const SeldonTranspose& Trans,
	   const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	   const Vector<T, Storage2, Allocator2>& X,
	   Vector<T, Storage3, Allocator3>& Y)
  {
    throw WrongArgument("Mlt", "Incompatible matrix-vector product");
  }


  template<class T, class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage3, class Allocator3>
  void Mlt(const SeldonTranspose& Trans,
	   const Matrix<T, Prop1, Storage1, Allocator1>& M,
	   const Vector<complex<T>, Storage2, Allocator2>& X,
	   Vector<complex<T>, Storage3, Allocator3>& Y)
  {
    if (Trans.NoTrans())
      MltVector(M, X, Y);
    else if (Trans.Trans())
      MltVector(SeldonTrans, M, X, Y);
    else
      MltVector(SeldonConjTrans, M, X, Y);
  }
  
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha, const SeldonTranspose& Trans,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<T, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<T, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      MltAddVector(alpha, M, X, beta, Y);
    else if (Trans.Trans())
      MltAddVector(alpha, SeldonTrans, M, X, beta, Y);
    else
      MltAddVector(alpha, SeldonConjTrans, M, X, beta, Y);
  }


  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha, const SeldonTranspose& Trans,
	      const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      MltAddVector(alpha, M, X, beta, Y);
    else if (Trans.Trans())
      MltAddVector(alpha, SeldonTrans, M, X, beta, Y);
    else
      MltAddVector(alpha, SeldonConjTrans, M, X, beta, Y);
  }
  
  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha, const SeldonTranspose& Trans,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      MltAddVector(alpha, M, X, beta, Y);
    else if (Trans.Trans())
      MltAddVector(alpha, SeldonTrans, M, X, beta, Y);
    else
      MltAddVector(alpha, SeldonConjTrans, M, X, beta, Y);
  }

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const T& alpha, const SeldonTranspose& Trans,
	      const Matrix<complex<T>, Prop1, Storage1, Allocator1>& M,
	      const Vector<T, Storage2, Allocator2>& X,
	      const T& beta,
	      Vector<T, Storage4, Allocator4>& Y)
  {
    throw WrongArgument("MltAdd", "Incompatible matrix-vector product");
  }

  template<class T,
	   class Prop1, class Storage1, class Allocator1,
	   class Storage2, class Allocator2,
	   class Storage4, class Allocator4>
  void MltAdd(const complex<T>& alpha, const SeldonTranspose& Trans,
	      const Matrix<T, Prop1, Storage1, Allocator1>& M,
	      const Vector<complex<T>, Storage2, Allocator2>& X,
	      const complex<T>& beta,
	      Vector<complex<T>, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      MltAddVector(alpha, M, X, beta, Y);
    else if (Trans.Trans())
      MltAddVector(alpha, SeldonTrans, M, X, beta, Y);
    else
      MltAddVector(alpha, SeldonConjTrans, M, X, beta, Y);
  }

  template<class T0, class Prop0, class Storage0, class Allocator0>
  void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x)
  {
    Vector<T0> xr(x.GetM()), xi(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
	xr = real(x(i));
	xi = imag(x(i));
      }
    
    SolveLuVector(A, pivot, xr);
    SolveLuVector(A, pivot, xi);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T0>(xr(i), xi(i));
  }

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<T0, General, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x)
  {
    Vector<T0> xr(x.GetM()), xi(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
	xr = real(x(i));
	xi = imag(x(i));
      }
    
    SolveLuVector(trans, A, pivot, xr);
    SolveLuVector(trans, A, pivot, xi);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T0>(xr(i), xi(i));
  }

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<T0, Symmetric, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x)
  {
    Vector<T0> xr(x.GetM()), xi(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
	xr = real(x(i));
	xi = imag(x(i));
      }
    
    SolveLuVector(A, pivot, xr);
    SolveLuVector(A, pivot, xi);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T0>(xr(i), xi(i));
  }
  
}

#define SELDON_FILE_FUNCTIONS_BASE_CXX
#endif

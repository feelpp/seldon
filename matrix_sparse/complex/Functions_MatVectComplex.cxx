// Copyright (C) 2001-2011 Marc Duruflé
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


#ifndef SELDON_FILE_FUNCTIONS_MATVECT_COMPLEX_CXX

namespace Seldon
{

  
  /*** Complex sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i, j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    complex<T1> zero(0);
    complex<T1> temp;

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * X(real_ind[j]);
	Y(i) += alpha * temp;
      }

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T1>(T1(0), imag_data[j]) * X(imag_ind[j]);
	Y(i) += alpha * temp;
      }
  }
  
  
  /*** Symmetric complex sparse matrices ***/
  
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i, j;
    complex<T1> zero(0);
    complex<T1> temp;

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i<ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * X(real_ind[j]);
	Y(i) += alpha * temp;
      }
    for (i = 0; i<ma-1; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	if (real_ind[j] != i)
	  Y(real_ind[j]) += alpha * real_data[j] * X(i);

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T1>(T1(0), imag_data[j]) * X(imag_ind[j]);
	Y(i) += alpha * temp;
      }
    for (i = 0; i<ma-1; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	if (imag_ind[j] != i)
	  Y(imag_ind[j]) += alpha * complex<T1>(T1(0), imag_data[j]) * X(i);
  }

  
  /*** Complex sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i, j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	Y(real_ind[j]) += alpha * real_data[j] * X(i);

    for (i = 0; i < ma; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	Y(imag_ind[j]) += alpha * complex<T1>(T1(0), imag_data[j]) * X(i);
  }


  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i, j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	Y(real_ind[j]) += alpha * real_data[j] * X(i);

    for (i = 0; i < ma; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	Y(imag_ind[j]) += alpha * complex<T1>(T1(0), - imag_data[j]) * X(i);
  }
  
  
  /*** Symmetric complex sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i, j;
    complex<T1> zero(0);
    complex<T1> temp;

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * X(real_ind[j]);
	Y(i) += temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	if (real_ind[j] != i)
	  Y(real_ind[j]) += real_data[j] * X(i);

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T1>(T1(0), - imag_data[j]) * X(imag_ind[j]);
	Y(i) += temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	if (imag_ind[j] != i)
	  Y(imag_ind[j]) += complex<T1>(T1(0), - imag_data[j]) * X(i);
  }

  
    /*** ArrayRowSymComplexSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const Matrix<T1, Symmetric,
	      ArrayRowSymComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		if (p == i)
		  C(i) += complex<T1>(0, val) * B(i);
		else
		  {
		    C(i) += complex<T1>(0, val) * B(p);
		    C(p) += complex<T1>(0, val) * B(i);
		  }
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * complex<T1>(0, A.ValueImag(i, k));
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonConjTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(),n,p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		if (p == i)
		  C(i) -= complex<T1>(0, val) * B(i);
		else
		  {
		    C(i) -= complex<T1>(0, val) * B(p);
		    C(p) -= complex<T1>(0, val) * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * complex<T1>(0, A.ValueImag(i, k));
		if (p == i)
		  C(i) -= val_cplx * B(i);
		else
		  {
		    C(i) -= val_cplx * B(p);
		    C(p) -= val_cplx * B(i);
		  }
	      }
	  }
      }
  }


  /*** ArrayRowComplexSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const Matrix<T1, General,
	      ArrayRowComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(), n, p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1, 0))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(i) += val * B(p);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(i) += complex<T1>(0, val) * B(p);
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		C(i) += val_cplx * B(p);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * complex<T1>(0, A.ValueImag(i, k));
		C(i) += val_cplx * B(p);
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans, const Matrix<T1, General,
	      ArrayRowComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonTrans& Trans, const Matrix<T1, General,
	      ArrayRowComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(),n,p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(p) += val * B(i);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(p) += complex<T1>(0, val) * B(i);
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		C(p) += val_cplx * B(i);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * complex<T1>(0, A.ValueImag(i, k));
		C(p) += val_cplx * B(i);
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonConjTrans& Trans, const Matrix<T1, General,
	      ArrayRowComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(),n,p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(p) += val * B(i);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(p) -= complex<T1>(0, val) * B(i);
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		C(p) += val_cplx * B(i);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx -= alpha * complex<T1>(0, A.ValueImag(i, k));
		C(p) += val_cplx * B(i);
	      }
	  }
      }
  }
  
  
  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SOR(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
	   Vector<complex<T2>, Storage2, Allocator2>& X,
	   const Vector<complex<T1>, Storage1, Allocator1>& B,
	   const T3& omega, int iter, int type_ssor = 2)
  {
    complex<T1> temp(0);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetLength() || ma != B.GetLength())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    int* ptr_real = A.GetRealPtr();
    int* ptr_imag = A.GetImagPtr();
    int* ind_real = A.GetRealInd();
    int* ind_imag = A.GetImagInd();
    typename Matrix<T0, Prop0, RowComplexSparse, Allocator0>::pointer
      data_real = A.GetRealData();
    typename Matrix<T0, Prop0, RowComplexSparse, Allocator0>::pointer
      data_imag = A.GetImagData();
    complex<T0> ajj(1);

    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = 0; j < ma; j++)
	  {
	    temp = complex<T1>(0);
	    ajj = A(j,j);
	    for (int k = ptr_real[j]; k < ptr_real[j+1]; k++)
	      temp += data_real[k] * X(ind_real[k]);

	    for (int k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
	      temp += complex<T1>(0, data_imag[k]) * X(ind_imag[k]);

	    temp = B(j) - temp + ajj*X(j);
	    X(j) = (T2(1) - omega) * X(j) + omega * temp / ajj;
	  }

    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = ma-1; j >= 0; j--)
	  {
	    temp = complex<T1>(0);
	    ajj = A(j,j);
	    for (int k = ptr_real[j]; k < ptr_real[j+1]; k++)
	      temp += data_real[k] * X(ind_real[k]);

	    for (int k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
	      temp += complex<T1>(0, data_imag[k]) * X(ind_imag[k]);

	    temp = B(j) - temp + ajj * X(j);
	    X(j) = (T2(1) - omega) * X(j) + omega * temp / ajj;
	  }
  }


  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SOR(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
	   Vector<complex<T2>, Storage2, Allocator2>& X,
	   const Vector<complex<T1>, Storage1, Allocator1>& B,
	   const T3& omega, int iter, int type_ssor = 2)
  {
    complex<T1> temp(0);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetLength() || ma != B.GetLength())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    complex<T0> ajj(1);

    if (type_ssor%2 == 0)
      for (int i = 0; i < iter; i++)
	{
	  for (int j = 0; j < ma; j++)
	    {
	      temp = complex<T1>(0);
	      ajj = complex<T0>(0);
	      for (int k = 0; k < A.GetRealRowSize(j); k++)
		{
		  temp += A.ValueReal(j,k) * X(A.IndexReal(j,k));
		  if (A.IndexReal(j,k) == j)
		    ajj += complex<T0>(A.ValueReal(j,k), 0);
		}
	      for (int k = 0; k < A.GetImagRowSize(j); k++)
		{
		  temp += complex<T0>(0,A.ValueImag(j,k))
		    * X(A.IndexImag(j,k));
		  if (A.IndexImag(j,k) == j)
		    ajj += complex<T0>(0, A.ValueImag(j,k));
		}

	      temp = B(j) - temp + ajj * X(j);
	      X(j) = (T2(1) - omega) * X(j) + omega * temp / ajj;
	    }
	}

    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	for (int j = ma-1; j >= 0; j--)
	  {
	    temp = complex<T1>(0);
	    ajj = complex<T0>(0);
	    for (int k = 0; k < A.GetRealRowSize(j); k++)
	      {
		temp += A.ValueReal(j,k) * X(A.IndexReal(j,k));
		if (A.IndexReal(j,k) == j)
		  ajj += complex<T0>(A.ValueReal(j,k), 0);
	      }
	    for (int k = 0; k < A.GetImagRowSize(j); k++)
	      {
		temp += complex<T0>(0, A.ValueImag(j,k))
		  * X(A.IndexImag(j,k));
		if (A.IndexImag(j,k) == j)
		  ajj += complex<T0>(0, A.ValueImag(j,k));
	      }

	    temp = B(j) - temp + ajj * X(j);
	    X(j) = (T2(1) - omega) * X(j) + omega * temp / ajj;
	  }
  }


  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SOR(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
	   Vector<complex<T2>, Storage2, Allocator2>& X,
	   const Vector<complex<T1>, Storage1, Allocator1>& B,
	   const T3& omega, int iter, int type_ssor = 2)
  {
    complex<T1> temp(0);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetLength() || ma != B.GetLength())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    int* ptr_real = A.GetRealPtr();
    int* ptr_imag = A.GetImagPtr();
    int* ind_real = A.GetRealInd();
    int* ind_imag = A.GetImagInd();
    T0* data_real = A.GetRealData();
    T0* data_imag = A.GetImagData();

    Vector<complex<T2>, Storage2, Allocator2> Y(ma);
    Y.Zero();
    complex<T0> ajj(1);
    int p;
    complex<T0> val(0);

    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A, A is symmetric, so L = U^t
    // forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
	  for (int j = 0; j < ma; j++)
	    {
	      // first we do X = (U + (1-omega)/omega D) X + B
	      temp = complex<T1>(0);
	      ajj = complex<T0>(0);
	      for (int k = ptr_real[j]; k < ptr_real[j+1]; k++)
		{
		  p = ind_real[k];
		  val = complex<T0>(data_real[k], 0);
		  if (p == j)
		    ajj += val;
		  else
		    temp += val * X(p);
		}
	      for (int k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
		{
		  p = ind_imag[k];
		  val = complex<T0>(0, data_imag[k]);
		  if (p == j)
		    ajj += val;
		  else
		    temp += val * X(p);
		}

	      temp = B(j) - temp;
	      X(j) = (T2(1) - omega) / omega * ajj * X(j) + temp;
	    }

	  for (int j = 0; j < ma; j++)
	    {
	      ajj = complex<T0>(0);
	      // Then we solve (D/omega - L) X = X.
	      for (int k = ptr_real[j]; k < ptr_real[j+1]; k++)
		{
		  p = ind_real[k];
		  val = complex<T0>(data_real[k], 0);
		  if (p == j)
		    ajj += val;
		}

	      for (int k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
		{
		  p = ind_imag[k];
		  val = complex<T0>(0, data_imag[k]);
		  if (p == j)
		    ajj += val;
		}
	      X(j) *= omega / ajj;

	      for (int k = ptr_real[j]; k < ptr_real[j+1]; k++)
		{
		  p = ind_real[k];
		  val = complex<T0>(data_real[k], 0);
		  if (p != j)
		    X(p) -= val*X(j);
		}

	      for (int k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
		{
		  p = ind_imag[k];
		  val = complex<T0>(0, data_imag[k]);
		  if (p != j)
		    X(p) -= val*X(j);
		}
	    }
	}

    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
	  Y.Zero();

	  for (int j = 0; j < ma; j++)
	    {
	      ajj = complex<T0>(0);
	      // Then we compute X = (L + (1-omega)/omega D) X + B.
	      for (int k = ptr_real[j]; k < ptr_real[j+1]; k++)
		{
		  p = ind_real[k];
		  val = complex<T0>(data_real[k], 0);
		  if (p == j)
		    ajj += val;
		  else
		    Y(p) += val * X(j);
		}

	      for (int k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
		{
		  p = ind_imag[k];
		  val = complex<T0>(0, data_imag[k]);
		  if (p == j)
		    ajj += val;
		  else
		    Y(p) += val * X(j);
		}
	      X(j) = (T2(1) - omega) / omega * ajj * X(j) + B(j) - Y(j);
	    }

	  for (int j = (ma-1); j >= 0; j--)
	    {
	      temp = complex<T1>(0);
	      ajj = complex<T0>(0);
	      // Then we solve (D/omega - U) X = X.
	      for (int k = ptr_real[j]; k < ptr_real[j+1]; k++)
		{
		  p = ind_real[k];
		  val = complex<T0>(data_real[k], 0);
		  if (p == j)
		    ajj += val;
		  else
		    temp += val * X(p);
		}

	      for (int k = ptr_imag[j]; k < ptr_imag[j+1]; k++)
		{
		  p = ind_imag[k];
		  val = complex<T0>(0, data_imag[k]);
		  if (p == j)
		    ajj += val;
		  else
		    temp += val * X(p);
		}
	      X(j) = (X(j) - temp) * omega / ajj;
	    }
	}
  }


  //! Successive overrelaxation.
  /*!
    Solving A X = B by using S.O.R algorithm.
    omega is the relaxation parameter, iter the number of iterations.
    type_ssor = 2 forward sweep
    type_ssor = 3 backward sweep
    type_ssor = 0 forward and backward sweep
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SOR(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
	   Vector<complex<T2>, Storage2, Allocator2>& X,
	   const Vector<complex<T1>, Storage1, Allocator1>& B,
	   const T3& omega, int iter, int type_ssor = 2)
  {
    complex<T1> temp(0);

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("SOR", "Matrix must be squared.");

    if (ma != X.GetLength() || ma != B.GetLength())
      throw WrongDim("SOR", "Matrix and vector dimensions are incompatible.");
#endif

    Vector<complex<T2>, Storage2, Allocator2> Y(ma);
    Y.Zero();
    complex<T0> ajj(1);
    int p;
    complex<T0> val(0);
    // Let us consider the following splitting : A = D - L - U
    // D diagonal of A
    // L lower part of A
    // U upper part of A, A is symmetric, so L = U^t
    // forward sweep
    // (D/omega - L) X^{n+1/2} = (U + (1-omega)/omega D) X^n + B
    if (type_ssor % 2 == 0)
      for (int i = 0; i < iter; i++)
	{
	  for (int j = 0; j < ma; j++)
	    {
	      // First we do X = (U + (1-omega)/omega D) X + B
	      temp = complex<T1>(0);
	      ajj = complex<T0>(0);
	      for (int k = 0; k < A.GetRealRowSize(j); k++)
		{
		  p = A.IndexReal(j,k);
		  val = complex<T0>(A.ValueReal(j,k), 0);
		  if (p == j)
		    ajj += val;
		  else
		    temp += val * X(p);
		}
	      for (int k = 0; k < A.GetImagRowSize(j); k++)
		{
		  p = A.IndexImag(j,k);
		  val = complex<T0>(0, A.ValueImag(j,k));
		  if (p == j)
		    ajj += val;
		  else
		    temp += val * X(p);
		}

	      temp = B(j) - temp;
	      X(j) = (T2(1) - omega) / omega * ajj * X(j) + temp;
	    }

	  for (int j = 0; j < ma; j++)
	    {
	      ajj = complex<T0>(0);
	      // Then we solve (D/omega - L) X = X.
	      for (int k = 0; k < A.GetRealRowSize(j); k++)
		{
		  p = A.IndexReal(j,k);
		  val = complex<T0>(A.ValueReal(j,k), 0);
		  if (p == j)
		    ajj += val;
		}
	      for (int k = 0; k < A.GetImagRowSize(j); k++)
		{
		  p = A.IndexImag(j,k);
		  val = complex<T0>(0, A.ValueImag(j,k));
		  if (p == j)
		    ajj += val;
		}
	      X(j) *= omega / ajj;
	      for (int k = 0; k < A.GetRealRowSize(j); k++)
		{
		  p = A.IndexReal(j,k);
		  val = complex<T0>(A.ValueReal(j,k), 0);
		  if (p != j)
		    X(p) -= val * X(j);
		}
	      for (int k = 0; k < A.GetImagRowSize(j); k++)
		{
		  p = A.IndexImag(j,k);
		  val = complex<T0>(0, A.ValueImag(j,k));
		  if (p != j)
		    X(p) -= val*X(j);
		}
	    }
	}

    // Backward sweep.
    // (D/omega - U) X^{n+1} = (L + (1-omega)/omega D) X^{n+1/2} + B
    if (type_ssor % 3 == 0)
      for (int i = 0; i < iter; i++)
	{
	  Y.Zero();

	  for (int j = 0; j < ma; j++)
	    {
	      ajj = complex<T0>(0);
	      // Then we compute X = (L + (1-omega)/omega D) X + B.
	      for (int k = 0; k < A.GetRealRowSize(j); k++)
		{
		  p = A.IndexReal(j,k);
		  val = complex<T0>(A.ValueReal(j,k), 0);
		  if (p == j)
		    ajj += val;
		  else
		    Y(p) += val * X(j);
		}
	      for (int k = 0; k < A.GetImagRowSize(j); k++)
		{
		  p = A.IndexImag(j,k);
		  val = complex<T0>(0, A.ValueImag(j,k));
		  if (p == j)
		    ajj += val;
		  else
		    Y(p) += val * X(j);
		}
	      X(j) = (T2(1) - omega) / omega * ajj * X(j) + B(j) - Y(j);
	    }

	  for (int j = ma-1; j >= 0; j--)
	    {
	      temp = complex<T1>(0);
	      ajj = complex<T0>(0);
	      // Then we solve (D/omega - U) X = X.
	      for (int k = 0; k < A.GetRealRowSize(j); k++)
		{
		  p = A.IndexReal(j,k);
		  val = complex<T0>(A.ValueReal(j,k), 0);
		  if (p == j)
		    ajj += val;
		  else
		    temp += val * X(p);
		}
	      for (int k = 0; k < A.GetImagRowSize(j); k++)
		{
		  p = A.IndexImag(j,k);
		  val = complex<T0>(0, A.ValueImag(j,k));
		  if (p == j)
		    ajj += val;
		  else
		    temp += val * X(p);
		}
	      X(j) = (X(j) - temp) * omega / ajj;
	    }
	}
  }

}

#define SELDON_FILE_FUNCTIONS_MATVECT_COMPLEX_CXX
#endif


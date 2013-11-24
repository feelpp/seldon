// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_COMPLEX_CXX

/*
  Functions defined in this file:
  (storage RowComplexSparse, ArrayRowComplexSparse, etc)
  
  alpha.A + B -> B
  Add(alpha, A, B)

  alpha.M -> M
  Mlt(alpha, M)

  A = A(row_perm, col_perm)
  ApplyPermutation(A, row_perm, col_perm)

  A(row_perm, col_perm) = A
  ApplyInversePermutation(A, row_perm, col_perm)
  
  A = Drow * A * Dcol
  ScaleMatrix(A, Drow, Dcol)
  
  A = Drow * A
  ScaleLeftMatrix(A, Drow)

  A = A * Dcol
  ScaleRightMatrix(A, Dcol)
  
  IsComplexMatrix(A)
*/

namespace Seldon
{
  
  //! B = B + alpha A
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric,
           ArrayRowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymComplexSparse, Allocator2>& B)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T2> > value(2*B.GetN());
    IVect index(2*B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	ni = A.GetImagRowSize(i);
        for (int j = 0; j < n; j++)
	  {
	    value(j) = alpha*complex<T2>(A.ValueReal(i, j), 0);
	    index(j) = A.IndexReal(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(j+n) = alpha*complex<T2>(0, A.ValueImag(i, j));
	    index(j+n) = A.IndexImag(i, j);
	  }
	
	B.AddInteractionRow(i, n+ni, index, value);
      }
  }
  

  //! B = B + alpha A  
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, Symmetric,
	   ArrayRowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B)
  {
    int m = B.GetM(), n;
    Vector<T2, VectFull, Allocator2> value(B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	for (int j = 0; j < n; j++)
	  value(j) = alpha*complex<T1>(A.ValueReal(i, j), 0);

	B.AddInteractionRow(i, n, A.GetRealInd(i), value.GetData());
	
        n = A.GetImagRowSize(i);
        for (int j = 0; j < n; j++)
	  value(j) = alpha*complex<T1>(0, A.ValueImag(i, j));

	B.AddInteractionRow(i, n, A.GetImagInd(i), value.GetData());
      }
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, General,
	   ArrayRowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowSparse, Allocator2>& B)
  {
    int m = B.GetM(), n;
    Vector<T2, VectFull, Allocator2> value(B.GetN());
    for (int i = 0; i < m; i++)
      {
	n = A.GetRealRowSize(i);
        for (int j = 0; j < n; j++)
          value(j) = alpha*complex<T1>(A.ValueReal(i, j), 0);
        
        B.AddInteractionRow(i, n, A.GetRealInd(i), value.GetData());
              
        n = A.GetImagRowSize(i);
        for (int j = 0; j < n; j++)
          value(j) = alpha*complex<T1>(0, A.ValueImag(i, j));
            
        B.AddInteractionRow(i, n, A.GetImagInd(i), value.GetData());
      }
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, General,
	   ArrayRowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowComplexSparse, Allocator2>& B)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T2> > value(2*B.GetN()); IVect index(2*B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	ni = A.GetImagRowSize(i);
        for (int j = 0; j < n; j++)
          {
            value(j) = alpha*complex<T2>(A.ValueReal(i, j), 0);
            index(j) = A.IndexReal(i, j);
          }
        
        for (int j = 0; j < ni; j++)
          {
            value(n+j) = alpha*complex<T2>(0, A.ValueImag(i, j));
            index(n+j) = A.IndexImag(i, j);
          }
            
        B.AddInteractionRow(i, n+ni, index, value);
      }
  }
  
  
  //! C = C + complex(A,B)
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2,
           class T3, class Prop3, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& A,
	   const Matrix<T2, Prop2, ArrayRowSymSparse, Allocator2>& B,
	   Matrix<T3, Prop3, ArrayRowSymComplexSparse, Allocator3>& C)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T3> > value(2*B.GetN()); IVect index(2*B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	ni = B.GetRowSize(i);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = alpha*complex<T3>(A.Value(i, j), 0);
	    index(j) = A.Index(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = alpha*complex<T3>(0, B.Value(i, j));
	    index(n+j) = B.Index(i, j);
	  }

	C.AddInteractionRow(i, n+ni, index, value);
      }
  }


  //! C = C + complex(A,B)
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2,
           class T3, class Prop3, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
	   const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
	   Matrix<T3, Prop3, ArrayRowComplexSparse, Allocator3>& C)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T3> > value(2*B.GetN());
    IVect index(2*B.GetN());
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	ni = B.GetRowSize(i);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = alpha*complex<T3>(A.Value(i, j), 0);
	    index(j) = A.Index(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = alpha*complex<T3>(0, B.Value(i, j));
	    index(n+j) = B.Index(i, j);
	  }

	C.AddInteractionRow(i, n+ni, index, value);
      }
  }
  

  template<class T, class Allocator>
  void Add_csr_ptr(const T& alpha, int* ptr_A, int* ind_A, T* data_A,
                   int* ptr_B, int* ind_B, T* data_B, int m,
                   Vector<int, VectFull, CallocAlloc<int> >& Ptr,
                   Vector<int, VectFull, CallocAlloc<int> >& Ind,
                   Vector<T, VectFull, Allocator>& Val)
  {
    int i = 0;
    int j = 0;
    int k;
    
    // A and B might have the same structure
    // Loop over all non-zeros. If the structures of A and B differ at any
    // time, the loop is broken and a different strategy is undertaken.
    for (i = 0; i < m; i++)
      if (ptr_A[i + 1] == ptr_B[i + 1])
        {
          for (j = ptr_A[i]; j < ptr_A[i + 1]; j++)
            if (ind_A[j] == ind_B[j])
              data_B[j] += alpha * data_A[j];
            else
              break;
          if (j != ptr_A[i + 1])
            break;
        }
      else
        break;
    
    // Success: A and B have the same structure.
    if (i == m)
      return;
    
    // The addition is performed row by row in the following lines. Thus the
    // additions already performed in the current line, if started, should be
    // canceled.
    for (k = ptr_A[i]; k < j; k++)
      if (ind_A[k] == ind_B[k])
        data_B[k] -= alpha * data_A[k];

    // Number of non zero entries currently found.
    int Nnonzero = ptr_A[i];
    
    // counting the number of non-zero entries
    int kb, jb(0), ka, ja(0);
    for (int i2 = i; i2 < m; i2++)
      {
        kb = ptr_B[i2];
        
        for (ka = ptr_A[i2]; ka < ptr_A[i2 + 1]; ka++)
          {
            ja = ind_A[ka];
            while (kb < ptr_B[i2 + 1] && ind_B[kb] < ja)
              {
                kb++;
                Nnonzero++;
              }
            
            if (kb < ptr_B[i2 + 1] && ja == ind_B[kb])
              kb++;
            
            Nnonzero++;
          }

        while (kb < ptr_B[i2 + 1])
          {
            kb++;
            Nnonzero++;
          }
      }
    
    // A and B do not have the same structure. An intermediate matrix will be
    // needed. The first i rows have already been added. These computations
    // are preserved in arrays Ptr, Ind Val.
    Ptr.Reallocate(m+1); Ind.Reallocate(Nnonzero);
    Val.Reallocate(Nnonzero);
    for (int i2 = 0; i2 <= i; i2++)
      Ptr(i2) = ptr_B[i2];
    
    for (j = 0; j < ptr_B[i]; j++)
      {
        Ind(j) = ind_B[j];
        Val(j) = data_B[j];
      }

    // Now deals with the remaining lines.
    Nnonzero = ptr_A[i];
    for (; i < m; i++)
      {
        kb = ptr_B[i];
        if (kb < ptr_B[i + 1])
          jb = ind_B[kb];
        for (ka = ptr_A[i]; ka < ptr_A[i + 1]; ka++)
          {
            ja = ind_A[ka];
            while (kb < ptr_B[i + 1] && jb < ja)
              // For all elements in B that are before the ka-th element of A.
              {
                Ind(Nnonzero) = jb;
                Val(Nnonzero) = data_B[kb];
                kb++;
                if (kb < ptr_B[i + 1])
                  jb = ind_B[kb];
                Nnonzero++;
              }

            if (kb < ptr_B[i + 1] && ja == jb)
              // The element in A is also in B.
              {
                Ind(Nnonzero) = jb;
                Val(Nnonzero) = data_B[kb] + alpha * data_A[ka];
                kb++;
                if (kb < ptr_B[i + 1])
                  jb = ind_B[kb];
              }
            else
              {
                Ind(Nnonzero) = ja;
                Val(Nnonzero) = alpha * data_A[ka];
              }
            Nnonzero++;
          }

        // The remaining elements from B.
        while (kb < ptr_B[i + 1])
          {
            Ind(Nnonzero) = jb;
            Val(Nnonzero) = data_B[kb];
            kb++;
            if (kb < ptr_B[i + 1])
              jb = ind_B[kb];
            Nnonzero++;
          }

        Ptr(i + 1) = Nnonzero;
      }    
  }
  
  
  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, RowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, RowSymComplexSparse, Allocator2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > PtrReal, IndReal, PtrImag, IndImag;
    Vector<T2, VectFull, Allocator2> DataReal, DataImag;
    Add_csr_ptr(alpha, A.GetRealPtr(), A.GetRealInd(), A.GetRealData(),
                B.GetRealPtr(), B.GetRealInd(), B.GetRealData(), B.GetM(),
                PtrReal, IndReal, DataReal);

    Add_csr_ptr(alpha, A.GetImagPtr(), A.GetImagInd(), A.GetImagData(),
                B.GetImagPtr(), B.GetImagInd(), B.GetImagData(), B.GetM(),
                PtrImag, IndImag, DataImag);

    B.SetData(B.GetM(), B.GetN(), DataReal, PtrReal, IndReal,
              DataImag, PtrImag, IndImag);
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, General, RowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, RowComplexSparse, Allocator2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > PtrReal, IndReal, PtrImag, IndImag;
    Vector<T2, VectFull, Allocator2> DataReal, DataImag;
    Add_csr_ptr(alpha, A.GetRealPtr(), A.GetRealInd(), A.GetRealData(),
                B.GetRealPtr(), B.GetRealInd(), B.GetRealData(), B.GetM(),
                PtrReal, IndReal, DataReal);

    Add_csr_ptr(alpha, A.GetImagPtr(), A.GetImagInd(), A.GetImagData(),
                B.GetImagPtr(), B.GetImagInd(), B.GetImagData(), B.GetM(),
                PtrImag, IndImag, DataImag);

    B.SetData(B.GetM(), B.GetN(), DataReal, PtrReal, IndReal,
              DataImag, PtrImag, IndImag);
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const complex<T0>& alpha,
	   const Matrix<T1, Symmetric, RowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, RowSymComplexSparse, Allocator2>& B)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Add(Matrix<RowSymComplexSparse>)",
                      "Function not implemented for complex scalars");

    Vector<int, VectFull, CallocAlloc<int> > PtrReal, IndReal, PtrImag, IndImag;
    Vector<T2, VectFull, Allocator2> DataReal, DataImag;
    Add_csr_ptr(real(alpha), A.GetRealPtr(), A.GetRealInd(), A.GetRealData(),
                B.GetRealPtr(), B.GetRealInd(), B.GetRealData(), B.GetM(),
                PtrReal, IndReal, DataReal);

    Add_csr_ptr(real(alpha), A.GetImagPtr(), A.GetImagInd(), A.GetImagData(),
                B.GetImagPtr(), B.GetImagInd(), B.GetImagData(), B.GetM(),
                PtrImag, IndImag, DataImag);

    B.SetData(B.GetM(), B.GetN(), DataReal, PtrReal, IndReal,
              DataImag, PtrImag, IndImag);
  }


  //! B = B + alpha A
  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const complex<T0>& alpha,
	   const Matrix<T1, General, RowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, RowComplexSparse, Allocator2>& B)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Add(Matrix<RowComplexSparse>)",
                      "Function not implemented for complex scalars");

    Vector<int, VectFull, CallocAlloc<int> > PtrReal, IndReal, PtrImag, IndImag;
    Vector<T2, VectFull, Allocator2> DataReal, DataImag;
    Add_csr_ptr(real(alpha), A.GetRealPtr(), A.GetRealInd(), A.GetRealData(),
                B.GetRealPtr(), B.GetRealInd(), B.GetRealData(), B.GetM(),
                PtrReal, IndReal, DataReal);

    Add_csr_ptr(real(alpha), A.GetImagPtr(), A.GetImagInd(), A.GetImagData(),
                B.GetImagPtr(), B.GetImagInd(), B.GetImagData(), B.GetM(),
                PtrImag, IndImag, DataImag);

    B.SetData(B.GetM(), B.GetN(), DataReal, PtrReal, IndReal,
              DataImag, PtrImag, IndImag);
  }
  

  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i, j) *= alpha;
        
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i, j) *= alpha;
        
      }
  }


  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const complex<T0>& alpha,
	   Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Mlt(Matrix<ArrayRowComplexSparse>)",
                      "Function not implemented for complex scalars");
    
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i, j) *= real(alpha);
        
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i, j) *= real(alpha);        
      }
  }


  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const T0& alpha,
           Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i,j) *= alpha;

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i,j) *= alpha;
      }
  }
  

  //! multiplication by a scalar  
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const complex<T0>& alpha,
           Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Mlt(Matrix<ArrayRowComplexSparse>)",
                      "Function not implemented for complex scalars");
    
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i,j) *= real(alpha);

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i,j) *= real(alpha);
      }
  }


  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    T* data_A = A.GetRealData();
    for (int i = 0; i < A.GetRealDataSize(); i++)
      data_A[i] *= alpha;

    data_A = A.GetImagData();
    for (int i = 0; i < A.GetImagDataSize(); i++)
      data_A[i] *= alpha;
  }
  

  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const complex<T0>& alpha,
	   Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Mlt(Matrix<RowComplexSparse>)",
                      "Function not implemented for complex scalars");
    
    T* data_A = A.GetRealData();
    T0 alpha_r = real(alpha);
    for (int i = 0; i < A.GetRealDataSize(); i++)
      data_A[i] *= alpha_r;

    data_A = A.GetImagData();
    for (int i = 0; i < A.GetImagDataSize(); i++)
      data_A[i] *= alpha_r;
  }


  //! multiplication by a scalar
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const T0& alpha, Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    T* data_A = A.GetRealData();
    for (int i = 0; i < A.GetRealDataSize(); i++)
      data_A[i] *= alpha;

    data_A = A.GetImagData();
    for (int i = 0; i < A.GetImagDataSize(); i++)
      data_A[i] *= alpha;
  }
  

  //! multiplication by a scalar  
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const complex<T0>& alpha,
           Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    if (imag(alpha) != T0(0))
      throw Undefined("Mlt(Matrix<RowSymComplexSparse>)",
                      "Function not implemented for complex scalars");
    
    T* data_A = A.GetRealData();
    T0 alpha_r = real(alpha);
    for (int i = 0; i < A.GetRealDataSize(); i++)
      data_A[i] *= alpha_r;

    data_A = A.GetImagData();
    for (int i = 0; i < A.GetImagDataSize(); i++)
      data_A[i] *= alpha_r;
  }

  
  //! Permutation of a general matrix stored by rows.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm)
  {
    int m = A.GetM();
    IVect ind_tmp, iperm(m), rperm(m);
    for (int i = 0; i < m; i++)
      {
	iperm(i) = i;
	rperm(i) = i;
      }
    
    // A(rperm(i),:) will be the place where is the initial row i.

    // Algorithm avoiding the allocation of another matrix.
    for (int i = 0; i < m; i++)
      {
	// We get the index of row where the row initially placed on row i is.
	int i2 = rperm(i);
	// We get the new index of this row.
	int i_ = row_perm(i);

	// We fill ind_tmp of the permuted indices of columns of row i.
	int nr = A.GetRealRowSize(i2);
	ind_tmp.Reallocate(nr);
	for (int j = 0; j < nr; j++)
	  ind_tmp(j) = col_perm(A.IndexReal(i2, j));

	// We swap the two rows i and its destination row_perm(i).
	A.SwapRealRow(i2, i_);
	A.ReplaceRealIndexRow(i_, ind_tmp);

	int ni = A.GetImagRowSize(i2);
	ind_tmp.Reallocate(ni);
	for (int j = 0; j < ni; j++)
	  ind_tmp(j) = col_perm(A.IndexImag(i2, j));

	A.SwapImagRow(i2, i_);
	A.ReplaceImagIndexRow(i_, ind_tmp);
        
	// We update the indices iperm and rperm in order to keep in memory
	// the place where the row row_perm(i) is.
	int i_tmp = iperm(i_);
	iperm(i_) = iperm(i2);
	iperm(i2) = i_tmp;
	rperm(iperm(i_)) = i_;
	rperm(iperm(i2)) = i2;

	// We assemble the row i.
	A.AssembleRealRow(i_);
        A.AssembleImagRow(i_);
      }
  }

  
  //! Permutation of a symmetric matrix stored by rows.
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop,
                               ArrayRowSymComplexSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm)
  {
    // It is assumed that the permuted matrix is still symmetric! For example,
    // the user can provide row_perm = col_perm.
    int m = A.GetM();
    int nnz_real = A.GetRealDataSize(), nnz_imag = A.GetImagDataSize();
    IVect IndRow(nnz_real), IndCol(nnz_real);
    Vector<T, VectFull, Allocator> Val(nnz_real);

    // First we convert the matrix in coordinate format and we permute the
    // indices.
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    int k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  {
	    IndRow(k) = row_perm(i);
	    Val(k) = A.ValueReal(i,j);
	    IndCol(k) = col_perm(A.IndexReal(i, j));
	    if (IndCol(k) <= IndRow(k))
	      {
		// We store only the superior part of the symmetric matrix.
		int ind_tmp = IndRow(k);
		IndRow(k) = IndCol(k);
		IndCol(k) = ind_tmp;
	      }
	    k++;
	  }
      }

    // We sort by row number.
    Sort(nnz_real, IndRow, IndCol, Val);

    // A is filled.
    k = 0;
    for (int i = 0; i < m; i++)
      {
	int first_index = k;
	// We get the size of the row i.
	while (k < nnz_real && IndRow(k) <= i)
	  k++;

	int size_row = k - first_index;
	// If row not empty.
	if (size_row > 0)
	  {
	    A.ReallocateRealRow(i, size_row);
	    k = first_index;
	    Sort(k, k+size_row-1, IndCol, Val);
	    for (int j = 0; j < size_row; j++)
	      {
		A.IndexReal(i,j) = IndCol(k);
		A.ValueReal(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearRealRow(i);
      }

    // Same procedure for imaginary part.

    IndRow.Reallocate(nnz_imag);
    IndCol.Reallocate(nnz_imag);
    Val.Reallocate(nnz_imag);

    // First we convert the matrix in coordinate format and we permute the
    // indices.
    // IndRow -> indices of the permuted rows
    // IndCol -> indices of the permuted columns
    k = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  {
	    IndRow(k) = row_perm(i);
	    Val(k) = A.ValueImag(i,j);
	    IndCol(k) = col_perm(A.IndexImag(i,j));
	    if (IndCol(k) <= IndRow(k))
	      {
		// We store only the superior part of the symmetric matrix.
		int ind_tmp = IndRow(k);
		IndRow(k) = IndCol(k);
		IndCol(k) = ind_tmp;
	      }
	    k++;
	  }
      }
    // We sort by row number.
    Sort(nnz_imag, IndRow, IndCol, Val);

    // A is filled
    k = 0;
    for (int i = 0; i < m; i++)
      {
	int first_index = k;
	// We get the size of the row i.
	while (k < nnz_imag && IndRow(k) <= i)
	  k++;
	int size_row = k - first_index;
	// If row not empty.
	if (size_row > 0)
	  {
	    A.ReallocateImagRow(i, size_row);
	    k = first_index;
	    Sort(k, k+size_row-1, IndCol, Val);
	    for (int j = 0; j < size_row; j++)
	      {
		A.IndexImag(i,j) = IndCol(k);
		A.ValueImag(i,j) = Val(k);
		k++;
	      }
	  }
	else
	  A.ClearImagRow(i);
      }
  }

  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, col_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    Vector<int> inv_col_perm(col_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    for (int i = 0; i < col_perm.GetM(); i++)
      inv_col_perm(col_perm(i)) = i;
    
    ApplyInversePermutation(A, inv_row_perm, inv_col_perm);
  }

  
  //! Permutation of rows and columns of a matrix
  /*!
    B(i, j) = A(row_perm(i), row_perm(j)) and A = B.
    Equivalent Matlab operation: A = A(row_perm, row_perm)
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    Vector<int> inv_row_perm(row_perm.GetM());
    for (int i = 0; i < row_perm.GetM(); i++)
      inv_row_perm(row_perm(i)) = i;

    ApplyInversePermutation(A, inv_row_perm, inv_row_perm);
  }
  
  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowSymComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale_left(i) * scale_right(A.IndexReal(i, j));

	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale_left(i) * scale_right(A.IndexImag(i, j));
      }
  }


  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale_left(i) * scale_right(A.IndexReal(i, j));

	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale_left(i) * scale_right(A.IndexImag(i, j));
      }
  }

  
  //! Each row is scaled.
  /*!
    We compute diag(S)*A where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ArrayRowComplexSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale(i);

	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale(i);
      }
  }
  
  
  //! Each column is scaled.
  /*!
    We compute A*diag(S) where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, ArrayRowComplexSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    for (int i = 0; i < m; i++ )
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++ )
	  A.ValueReal(i,j) *= scale(A.IndexReal(i, j));

	for (int j = 0; j < A.GetImagRowSize(i); j++ )
	  A.ValueImag(i,j) *= scale(A.IndexImag(i, j));
      }
  }

  
  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, RowSymComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T1* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T1* data_imag = A.GetImagData();
    for (int i = 0; i < m; i++ )
      {
	for (int j = ptr_real[i]; j < ptr_real[i+1]; j++)
	  data_real[j] *= scale_left(i) * scale_right(ind_real[j]);

	for (int j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
	  data_imag[j] *= scale_left(i) * scale_right(ind_imag[j]);
      }
  }


  //! Each row and column are scaled.
  /*!
    We compute diag(scale_left)*A*diag(scale_right).
  */
  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, RowComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right)
  {
    int m = A.GetM();
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T1* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T1* data_imag = A.GetImagData();
    for (int i = 0; i < m; i++ )
      {
	for (int j = ptr_real[i]; j < ptr_real[i+1]; j++ )
	  data_real[j] *= scale_left(i) * scale_right(ind_real[j]);

	for (int j = ptr_imag[i]; j < ptr_imag[i+1]; j++ )
	  data_imag[j] *= scale_left(i) * scale_right(ind_imag[j]);
      }
  }

  
  //! Each row is scaled.
  /*!
    We compute diag(S)*A where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, RowComplexSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    int* ptr_real = A.GetRealPtr();
    T1* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    T1* data_imag = A.GetImagData();
    for (int i = 0; i < m; i++ )
      {
	for (int j = ptr_real[i]; j < ptr_real[i+1]; j++ )
	  data_real[j] *= scale(i);

	for (int j = ptr_imag[i]; j < ptr_imag[i+1]; j++ )
	  data_imag[j] *= scale(i);
      }
  }
  
  
  //! Each column is scaled.
  /*!
    We compute A*diag(S) where S = scale.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, RowComplexSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale)
  {
    int m = A.GetM();
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T1* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T1* data_imag = A.GetImagData();
    for (int i = 0; i < m; i++ )
      {
	for (int j = ptr_real[i]; j < ptr_real[i+1]; j++ )
	  data_real[j] *= scale(ind_real[j]);

	for (int j = ptr_imag[i]; j < ptr_imag[i+1]; j++ )
	  data_imag[j] *= scale(ind_imag[j]);
      }
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayRowComplexSparse, Allocator>& A,
                 Matrix<T, General, ArrayRowComplexSparse, Allocator>& B)
  {
    B.Clear();
    
    int m = A.GetM();
    int n = A.GetN();
    Vector<int, VectFull, CallocAlloc<int> > ptr_r(n), ptr_i(n);
    
    B.Reallocate(n, m);

    // For each column j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    ptr_r.Zero();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRealRowSize(i); j++)
        ptr_r(A.IndexReal(i, j))++;

    ptr_i.Zero();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetImagRowSize(i); j++)
        ptr_i(A.IndexImag(i, j))++;
    
    for (int i = 0; i < n; i++)
      {
        B.ReallocateRealRow(i, ptr_r(i));
        B.ReallocateImagRow(i, ptr_i(i));
      }
    
    // filling matrix B
    ptr_r.Zero();
    ptr_i.Zero();
    for (int i = 0; i < m; i++)
      {
        for (int jp = 0; jp < A.GetRealRowSize(i); jp++)
          {
            int j = A.IndexReal(i, jp);
            int k = ptr_r(j);
            ++ptr_r(j);
            B.ValueReal(j, k) = A.ValueReal(i, jp);
            B.IndexReal(j, k) = i;
          }

        for (int jp = 0; jp < A.GetImagRowSize(i); jp++)
          {
            int j = A.IndexImag(i, jp);
            int k = ptr_i(j);
            ++ptr_i(j);
            B.ValueImag(j, k) = A.ValueImag(i, jp);
            B.IndexImag(j, k) = i;
          }
      }
    
    // sorting numbers 
    B.Assemble();
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, RowComplexSparse, Allocator>& A,
                 Matrix<T, General, RowComplexSparse, Allocator>& B)
  {
    B.Clear();
    
    int m = A.GetM();
    int n = A.GetN();
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T* data_imag = A.GetImagData();
    Vector<int, VectFull, CallocAlloc<int> > ptr_r(n), ptr_i(n);
    
    // For each column j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    ptr_r.Zero();
    for (int i = 0; i < m; i++)
      for (int j = ptr_real[i]; j < ptr_real[i+1]; j++)
        ptr_r(ind_real[j])++;

    ptr_i.Zero();
    for (int i = 0; i < m; i++)
      for (int j = ptr_imag[i]; j < ptr_imag[i+1]; j++)
        ptr_i(ind_imag[j])++;
    
    Vector<T, VectFull, Allocator> ValReal(A.GetRealDataSize());
    Vector<T, VectFull, Allocator> ValImag(A.GetImagDataSize());
    Vector<int, VectFull, CallocAlloc<int> > PtrReal(n+1), PtrImag(n+1),
      IndReal(A.GetRealDataSize()), IndImag(A.GetImagDataSize());
    
    PtrReal(0) = 0; PtrImag(0) = 0;
    for (int i = 0; i < n; i++)
      {
        PtrReal(i+1) = PtrReal(i) + ptr_r(i);
        PtrImag(i+1) = PtrImag(i) + ptr_i(i);
      }

    // filling matrix B
    ptr_r.Zero();
    ptr_i.Zero();
    for (int i = 0; i < m; i++)
      {
        for (int jp = ptr_real[i]; jp < ptr_real[i+1]; jp++)
          {
            int j = ind_real[jp];
            int k = ptr_r(j);
            ValReal(PtrReal(j) + k) = data_real[jp];
            IndReal(PtrReal(j) + k) = i;
            ++ptr_r(j);
          }

        for (int jp = ptr_imag[i]; jp < ptr_imag[i+1]; jp++)
          {
            int j = ind_imag[jp];
            int k = ptr_i(j);
            ValImag(PtrImag(j) + k) = data_imag[jp];
            IndImag(PtrImag(j) + k) = i;
            ++ptr_i(j);
          }
      }
    
    // sorting numbers
    for (int i = 0; i < n; i++)
      {
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
      }
    
    B.SetData(n, m, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  T MaxAbs(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    T res(0);
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                res = max(res, abs(A.ValueImag(i, ji)));
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                res = max(res, ComplexAbs(complex<T>(A.ValueReal(i, j), 
						     A.ValueImag(i, ji))));
                ji++;
              }
            else
              res = max(res, abs(A.ValueReal(i, j)));
          }
        
        while (ji < size_imag)
          {
            res = max(res, abs(A.ValueImag(i, ji)));
            ji++;
          }
      }
    
    return res;
  }


  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  T MaxAbs(const Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    T res(0);
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T* data_imag = A.GetImagData();
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = ptr_imag[i];
        for (int j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                res = max(res, abs(data_imag[ji]));
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                res = max(res, ComplexAbs(complex<T>(data_real[j], 
						     data_imag[ji])));
                ji++;
              }
            else
              res = max(res, abs(data_real[j]));
          }
        
        while (ji < ptr_imag[i+1])
          {
            res = max(res, abs(data_imag[ji]));
            ji++;
          }
      }
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  T Norm1(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    Vector<T> sum(A.GetN());
    sum.Fill(T(0));
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                sum(A.IndexImag(i, ji)) += abs(A.ValueImag(i, ji));
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                sum(k) += ComplexAbs(complex<T>(A.ValueReal(i, j), 
						A.ValueImag(i, ji)));
                ji++;
              }
            else
              sum(k) += abs(A.ValueReal(i, j));
          }
        
        while (ji < size_imag)
          {
            sum(A.IndexImag(i, ji)) += abs(A.ValueImag(i, ji));
            ji++;
          }
      }
    
    return sum.GetNormInf();
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  T Norm1(const Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    Vector<T> sum(A.GetN());
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T* data_imag = A.GetImagData();
    sum.Fill(T(0));
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = ptr_imag[i];
        for (int j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                sum(ind_imag[ji]) += abs(data_imag[ji]);
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                sum(k) += ComplexAbs(complex<T>(data_real[j], 
						data_imag[ji]));
                ji++;
              }
            else
              sum(k) += abs(data_real[j]);
          }
        
        while (ji < ptr_imag[i+1])
          {
            sum(ind_imag[ji]) += abs(data_imag[ji]);
            ji++;
          }
      }
    
    return sum.GetNormInf();
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  T NormInf(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    T res(0), sum;
    for (int i = 0; i < A.GetM(); i++)
      {
        sum = T(0);
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                sum += abs(A.ValueImag(i, ji));
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                sum += ComplexAbs(complex<T>(A.ValueReal(i, j), 
					     A.ValueImag(i, ji)));
                ji++;
              }
            else
              sum += abs(A.ValueReal(i, j));
          }
        
        while (ji < size_imag)
          {
            sum += abs(A.ValueImag(i, ji));
            ji++;
          }
        
        res = max(res, sum);
      }
    
    return res;
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  T NormInf(const Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    T res(0), sum;
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T* data_imag = A.GetImagData();
    for (int i = 0; i < A.GetM(); i++)
      {
        sum = T(0);
        int ji = ptr_imag[i];
        for (int j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                sum += abs(data_imag[ji]);
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                sum += ComplexAbs(complex<T>(data_real[j], 
					     data_imag[ji]));
                ji++;
              }
            else
              sum += abs(data_real[j]);
          }
        
        while (ji < ptr_imag[i+1])
          {
            sum += abs(data_imag[ji]);
            ji++;
          }
        
        res = max(res, sum);
      }
    
    return res;
  }

  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  T MaxAbs(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    T res(0);
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                res = max(res, abs(A.ValueImag(i, ji)));
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                res = max(res, ComplexAbs(complex<T>(A.ValueReal(i, j), 
						     A.ValueImag(i, ji))));
                ji++;
              }
            else
              res = max(res, abs(A.ValueReal(i, j)));
          }
        
        while (ji < size_imag)
          {
            res = max(res, abs(A.ValueImag(i, ji)));
            ji++;
          }
      }
    
    return res;
  }


  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  T MaxAbs(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    T res(0);
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T* data_imag = A.GetImagData();
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = ptr_imag[i];
        for (int j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                res = max(res, abs(data_imag[ji]));
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                res = max(res, ComplexAbs(complex<T>(data_real[j], 
						     data_imag[ji])));
                ji++;
              }
            else
              res = max(res, abs(data_real[j]));
          }
        
        while (ji < ptr_imag[i+1])
          {
            res = max(res, abs(data_imag[ji]));
            ji++;
          }
      }
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  T Norm1(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    T val;
    Vector<T> sum(A.GetN());
    sum.Fill(T(0));
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = 0;
        int size_imag = A.GetImagRowSize(i);
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int k = A.IndexReal(i, j);
            while ( ji < size_imag && A.IndexImag(i, ji) < k)
              {
                val = abs(A.ValueImag(i, ji));
                sum(A.IndexImag(i, ji)) += val;
                if (A.IndexImag(i, ji) != i)
                  sum(i) += val;
                
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                val = ComplexAbs(complex<T>(A.ValueReal(i, j), 
                                            A.ValueImag(i, ji)));
                sum(k) += val;
                if (k != i)
                  sum(i) += val;
                
                ji++;
              }
            else
              {
                val = abs(A.ValueReal(i, j));
                sum(k) += val;
                if (k != i)
                  sum(i) += val;
              }
          }
        
        while (ji < size_imag)
          {
            val = abs(A.ValueImag(i, ji));
            sum(A.IndexImag(i, ji)) += val;
            if (A.IndexImag(i, ji) != i)
              sum(i) += val;
            
            ji++;
          }
      }
    
    return sum.GetNormInf();
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  T Norm1(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    Vector<T> sum(A.GetN());
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    T* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    T* data_imag = A.GetImagData();
    sum.Fill(T(0));
    T val;
    for (int i = 0; i < A.GetM(); i++)
      {
        int ji = ptr_imag[i];
        for (int j = ptr_real[i]; j < ptr_real[i+1]; j++)
          {
            int k = ind_real[j];
            while ( ji < ptr_imag[i+1] && ind_imag[ji] < k)
              {
                val = abs(data_imag[ji]);
                sum(ind_imag[ji]) += val;
                if (ind_imag[ji] != i)
                  sum(i) += val;
                
                ji++;
              }
            
            if ( ji < ptr_imag[i+1] && (ind_imag[ji] == k))              
              {
                val = ComplexAbs(complex<T>(data_real[j], 
                                            data_imag[ji]));
                sum(k) += val;
                if (k != i)
                  sum(i) += val;
                
                ji++;
              }
            else
              {
                val = abs(data_real[j]);
                sum(k) += val;
                if (k != i)
                  sum(i) += val;
              }
          }
        
        while (ji < ptr_imag[i+1])
          {
            val = abs(data_imag[ji]);
            sum(ind_imag[ji]) += val;
            if (ind_imag[ji] != i)
              sum(i) += val;
            
            ji++;
          }
      }
    
    return sum.GetNormInf();
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  T NormInf(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    return Norm1(A);
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  T NormInf(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    return Norm1(A);
  }
  
  
  //! A is replaced by its conjugate
  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetImagRowSize(i); j++)
        A.ValueImag(i, j) = -A.ValueImag(i, j);
  }


  //! A is replaced by its conjugate
  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetImagRowSize(i); j++)
        A.ValueImag(i, j) = -A.ValueImag(i, j); 
  }


  //! A is replaced by its conjugate
  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    T* data = A.GetImagData();
    for (int i = 0; i < A.GetImagDataSize(); i++)
      data[i] = -data[i];
  }


  //! A is replaced by its conjugate
  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    T* data = A.GetImagData();
    for (int i = 0; i < A.GetImagDataSize(); i++)
      data[i] = -data[i];
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, Symmetric,
                 ArrayRowSymComplexSparse, Allocator>& A,
                 Matrix<T, Symmetric,
                 ArrayRowSymComplexSparse, Allocator>& B)
  {
    B = A;
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayRowComplexSparse, Allocator>& A)
  {
    Matrix<T, General, ArrayRowComplexSparse, Allocator> B(A);
    Transpose(B, A);
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, Symmetric, ArrayRowSymComplexSparse, Allocator>& A)
  {
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, Symmetric,
                 RowSymComplexSparse, Allocator>& A,
                 Matrix<T, Symmetric,
                 RowSymComplexSparse, Allocator>& B)
  {
    B = A;
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, RowComplexSparse, Allocator>& A)
  {
    Matrix<T, General, RowComplexSparse, Allocator> B(A);
    Transpose(B, A);
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, Symmetric, RowSymComplexSparse, Allocator>& A)
  {
  }
  
  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, RowComplexSparse, Allocator>& A)
  {
    return true;
  }


  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ColComplexSparse, Allocator>& A)
  {
    return true;
  }

  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A)
  {
    return true;
  }

  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ColSymComplexSparse, Allocator>& A)
  {
    return true;
  }

  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A)
  {
    return true;
  }

  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ArrayColComplexSparse, Allocator>& A)
  {
    return true;
  }

  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A)
  {
    return true;
  }

  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>& A)
  {
    return true;
  }

}

#define SELDON_FILE_FUNCTIONS_MATRIX_COMPLEX_CXX
#endif


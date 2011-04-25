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
  
  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric,
	   ArrayRowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymComplexSparse, Allocator2>& B)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T2> > value;
    IVect index;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	ni = A.GetImagRowSize(i);
	value.Reallocate(n + ni);
	index.Reallocate(n + ni);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = A.ValueReal(i, j);
	    index(j) = A.IndexReal(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(j+n) = complex<T2>(0, 1) * A.ValueImag(i, j);
	    index(j+n) = A.IndexImag(i, j);
	  }

	Mlt(alpha, value);
	B.AddInteractionRow(i, n+ni, index, value);
      }
  }
  
  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, Symmetric,
	   ArrayRowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B)
  {
    int m = B.GetM(), n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.ValueReal(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetRealInd(i), value.GetData());
	n = A.GetImagRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = complex<T1>(0, 1) * A.ValueImag(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetImagInd(i), value.GetData());
      }
  }


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, General,
	   ArrayRowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowSparse, Allocator2>& B)
  {
    int m = B.GetM(),n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0; i < m; i++)
      {
	n = A.GetRealRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.ValueReal(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetRealInd(i), value.GetData());
	n = A.GetImagRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = complex<T1>(0, 1) * A.ValueImag(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetImagInd(i), value.GetData());
      }
  }


  template<class T0, class T1, class T2, class Allocator1,class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, General,
	   ArrayRowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowComplexSparse, Allocator2>& B)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T2> > value; IVect index;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	ni = A.GetImagRowSize(i);
	value.Reallocate(n + ni);
	index.Reallocate(n + ni);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = A.ValueReal(i, j);
	    index(j) = A.IndexReal(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = complex<T2>(0, 1) * A.ValueImag(i, j);
	    index(n+j) = A.IndexImag(i, j);
	  }

	Mlt(alpha, value);
	B.AddInteractionRow(i, n+ni, index, value);
      }
  }
  
  
  // C = C + complex(A,B)
  template<class T0, class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
	   const Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B,
	   Matrix<T3, Symmetric, ArrayRowSymComplexSparse, Allocator3>& C)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T3> > value; IVect index;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	ni = B.GetRowSize(i);
	value.Reallocate(n+ni);
	index.Reallocate(n+ni);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = A.Value(i, j);
	    index(j) = A.Index(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = complex<T3>(0, 1) * B.Value(i, j);
	    index(n+j) = B.Index(i, j);
	  }

	Mlt(alpha, value);
	C.AddInteractionRow(i, n+ni, index, value);
      }
  }


  // C = C + complex(A,B)
  template<class T0, class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	   const Matrix<T2, General, ArrayRowSparse, Allocator2>& B,
	   Matrix<T3, General, ArrayRowComplexSparse, Allocator3>& C)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T3> > value;
    IVect index;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	ni = B.GetRowSize(i);
	value.Reallocate(n + ni);
	index.Reallocate(n + ni);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = A.Value(i, j);
	    index(j) = A.Index(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = complex<T3>(0, 1) * B.Value(i, j);
	    index(n+j) = B.Index(i, j);
	  }

	Mlt(alpha, value);
	C.AddInteractionRow(i, n+ni, index, value);
      }
  }
  
  
  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, General, ArrayRowComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i, j) *= alpha;
        
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i, j) *= alpha;
        
      }
  }


  template<class T0, class T, class Allocator>
  void Mlt(const complex<T0>& alpha,
	   Matrix<T, General, ArrayRowComplexSparse, Allocator>& A)
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


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, Matrix<T, Symmetric,
	   ArrayRowSymComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i,j) *= alpha;

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i,j) *= alpha;
      }
  }
  
  
  template<class T0, class T, class Allocator>
  void Mlt(const complex<T0>& alpha, Matrix<T, Symmetric,
	   ArrayRowSymComplexSparse, Allocator>& A)
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
                sum(A.IndexImag(i, ji)) += abs(A.ValueImag(i, ji));
                if (A.IndexImag(i, ji) != i)
                  sum(i) += abs(A.ValueImag(i, ji));
                
                ji++;
              }
            
            if ( ji < size_imag && (A.IndexImag(i, ji) == k))              
              {
                T val = ComplexAbs(complex<T>(A.ValueReal(i, j), 
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
  
  
  //! A is replaced by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, General, ArrayRowComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetImagRowSize(i); j++)
        A.ValueImag(i, j) = -A.ValueImag(i, j);
  }


  //! A is replaced by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, Symmetric, ArrayRowSymComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetImagRowSize(i); j++)
        A.ValueImag(i, j) = -A.ValueImag(i, j); 
  }

  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayRowSymComplexSparse, Allocator>& A,
                 Matrix<T, Symmetric, ArrayRowSymComplexSparse, Allocator>& B)
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


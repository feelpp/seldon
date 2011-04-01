// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2009 Vivien Mallet
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
	   Matrix<T2, Symmetric, ArrayRowSparse, Allocator2>& B)
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
	  A.ValueReal(i,j) *= alpha;
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i,j) *= alpha;
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
  
  
  //! Permutation of a symmetric matrix stored by rows.
  /*!
    B(row_perm(i),col_perm(j)) = A(i,j) and A = B.
    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop,
                               ArrayRowSymComplexSparse, Allocator>& A,
                               const IVect& row_perm,const IVect& col_perm)
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
    We compute diag(S)*A where S = scale.  In order to keep symmetry, the
    operation is performed on upper part of the matrix, considering that lower
    part is affected by this operation.
  */
  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ArrayRowSymComplexSparse, Allocator1>& A,
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
  
}

#define SELDON_FILE_FUNCTIONS_MATRIX_COMPLEX_CXX
#endif


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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_ARRAY_CXX

/*
  Functions defined in this file:
  (storage ArrayRowSparse, ArrayColSparse, etc)

  X = A(i, :)
  GetRow(A, i, X)

  X = A(:, j)
  GetCol(A, j, X)

  A(i, :) = X
  SetRow(X, i, A)

  A(:, j) = X
  SetCol(X, j, A)

  alpha.M*X + beta.Y -> Y
  MltAdd(alpha, M, X, beta, Y)

  alpha.A + B -> B
  Add(alpha, A, B)

  alpha.M -> M
  Mlt(alpha, M)
  
  Highest absolute value of A.
  MaxAbs(A)

  1-norm of matrix A.
  Norm1(A)

  infinity norm of matrix A.
  NormInf(A)

  transpose of matrix A
  Transpose(A)

  B = transpose(A)
  Transpose(A, B)
  
  conjugate of transpose of matrix A
  TransposeConj(A)

  conjugate of matrix A
  Conjugate(A)

*/

namespace Seldon
{

  
  ////////////////////
  // GetRow, SetRow //
  
  
  //! Extracts a row from a sparse matrix
  /*!
    \param A sparse matrix
    \param i row index
    \param X row extracted
    X = A(i, :)
  */
  template<class T0, class Allocator0,
	   class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, ArrayRowSparse, Allocator0>& A,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = A.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif
    
    int size_row = A.GetRowSize(i);
    X.Reallocate(size_row);
    for (int j = 0; j < size_row; j++)
      {
	X.Index(j) = A.Index(i, j);
	X.Value(j) = A.Value(i, j);
      }
  }
  

  //! Extracts a row from a sparse matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, ArrayColSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    list<pair<int, T0> > vec;
    for (int j = 0; j < M.GetN(); j++)
      for (int k = 0; k < M.GetColumnSize(j); k++)
	if (M.Index(j, k) == i)
	  vec.push_back(make_pair(j, M.Value(j, k)));
    
    typename list<pair<int, T0> >::iterator it;
    X.Reallocate(vec.size());
    int j = 0;
    for (it = vec.begin(); it != vec.end(); it++)
      {
	X.Index(j) = it->first;
	X.Value(j) = it->second;
	j++;
      }
  }

  
  //! Extracts a row from a sparse matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    list<pair<int, T0> > vec;
    // beginning of row i
    for (int j = 0; j < i; j++)
      for (int k = 0; k < M.GetRowSize(j); k++)
	if (M.Index(j, k) == i)
	  vec.push_back(make_pair(j, M.Value(j, k)));
    
    // end of row i
    for (int k = 0; k < M.GetRowSize(i); k++)
      vec.push_back(make_pair(M.Index(i, k), M.Value(i, k)));
    
    typename list<pair<int, T0> >::iterator it;
    X.Reallocate(vec.size());
    int j = 0;
    for (it = vec.begin(); it != vec.end(); it++)
      {
	X.Index(j) = it->first;
	X.Value(j) = it->second;
	j++;
      }
  }
  
  
  //! Extracts a row from a sparse matrix
  /*!
    \param M sparse matrix
    \param i row index
    \param X row extracted
    X = M(i, :)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    list<pair<int, T0> > vec;
    // beginning of row i
    for (int k = 0; k < M.GetColumnSize(i); k++)
      vec.push_back(make_pair(M.Index(i, k), M.Value(i, k)));
    
    // end of row i
    for (int j = i+1; j < M.GetN(); j++)
      for (int k = 0; k < M.GetColumnSize(j); k++)
	if (M.Index(j, k) == i)
	  vec.push_back(make_pair(j, M.Value(j, k)));
    
    typename list<pair<int, T0> >::iterator it;
    X.Reallocate(vec.size());
    int j = 0;
    for (it = vec.begin(); it != vec.end(); it++)
      {
	X.Index(j) = it->first;
	X.Value(j) = it->second;
	j++;
      }
  }
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, ArrayRowSparse, Allocator0>& M)
  {
    M.ClearRow(i);
    M.ReallocateRow(i, X.GetM());    
    for (int k = 0; k < X.GetM(); k++)
      {
	M.Index(i, k) = X.Index(k);
	M.Value(i, k) = X.Value(k);
      }
  }
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, ArrayColSparse, Allocator0>& M)
  {
    T0 val;
    int p = 0;
    for (int j = 0; j < M.GetN(); j++)
      {
	while ( (p < X.GetM()) && (X.Index(p) < j))
	  p++;
	
	bool present_X = false;
	if ( (p < X.GetM()) && (X.Index(p) == j))	  
	  {
	    present_X = true;
	    val = X.Value(p);
	  }
	
	int size_col = M.GetColumnSize(j);
	bool present_val = false;
	int k = 0;
	while ( (k < size_col) && (M.Index(j, k) < i))
	  k++;
	
	if ( (k < size_col) && (M.Index(j, k) == i))
	  {
	    if (!present_X)
	      {
		// reducing size of column
		if (size_col > 1)
		  {
		    if (k == size_col-1)
		      M.ResizeColumn(j, size_col-1);
		    else
		      {
			int last_col = M.Index(j, size_col-1);
			val = M.Value(j, size_col-1);
			M.ResizeColumn(j, size_col-1);
			for (int q = k; q < size_col-2; q++)
			  {
			    M.Index(j, q) = M.Index(j, q+1);
			    M.Value(j, q) = M.Value(j, q+1);
			  }
			
			M.Index(j, size_col-2) = last_col;
			M.Value(j, size_col-2) = val;
		      }
		  }
		else
		  M.ClearColumn(j);
	      }
	    else
	      M.Value(j, k) = val;
	    
	    present_val = true;
	  }
	
	if (!present_val && present_X)
	  {
	    // increasing size of column
	    M.ResizeColumn(j, size_col+1);
	    for (int q = size_col; q > k; q--)
	      {
		M.Index(j, q) = M.Index(j, q-1);
		M.Value(j, q) = M.Value(j, q-1);
	      }
	    
	    M.Index(j, k) = i;
	    M.Value(j, k) = val;
	  }
      }
  }

  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M)
  {
    T0 val;
    int p = 0;
    for (int j = 0; j < i; j++)
      {
	while ( (p < X.GetM()) && (X.Index(p) < j))
	  p++;
	
	bool present_X = false;
	if ( (p < X.GetM()) && (X.Index(p) == j))	  
	  {
	    present_X = true;
	    val = X.Value(p);
	    p++;
	  }
	
	int size_row = M.GetRowSize(j);
	bool present_val = false;
	int k = 0;
	while ( (k < size_row) && (M.Index(j, k) < i))
	  k++;
	
	if ( (k < size_row) && (M.Index(j, k) == i))
	  {
	    if (!present_X)
	      {
		// reducing size of row
		if (size_row > 1)
		  {
		    if (k == size_row-1)
		      M.ResizeRow(j, size_row-1);
		    else
		      {
			int last_col = M.Index(j, size_row-1);
			val = M.Value(j, size_row-1);
			M.ResizeRow(j, size_row-1);
			for (int q = k; q < size_row-2; q++)
			  {
			    M.Index(j, q) = M.Index(j, q+1);
			    M.Value(j, q) = M.Value(j, q+1);
			  }
			
			M.Index(j, size_row-2) = last_col;
			M.Value(j, size_row-2) = val;
		      }
		  }
		else
		  M.ClearRow(j);
	      }
	    else
	      M.Value(j, k) = val;
	    
	    present_val = true;
	  }
	
	if (!present_val && present_X)
	  {
	    // increasing size of row
	    M.ResizeRow(j, size_row+1);
	    for (int q = size_row; q > k; q--)
	      {
		M.Index(j, q) = M.Index(j, q-1);
		M.Value(j, q) = M.Value(j, q-1);
	      }
	    
	    M.Index(j, k) = i;
	    M.Value(j, k) = val;
	  }
      }
    
    // and changing row i
    M.ClearRow(i);
    if (p < X.GetM())
      {
	M.ReallocateRow(i, X.GetM() - p);
	int k = 0;
	while (p < X.GetM())
	  {
	    M.Index(i, k) = X.Index(p);
	    M.Value(i, k) = X.Value(p);
	    k++;
	    p++;
	  }
      }
  }
  
  
  //! Sets a row of a matrix
  /*!
    \param M matrix
    \param i row index
    \param X new row of M
    M(i, :) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M)
  {
    T1 val;
    int p = 0;
    // setting column i
    M.ClearColumn(i);
    while ( (p < X.GetM()) && (X.Index(p) <= i))
      p++;
    
    if (p > 0)
      {
	M.ReallocateColumn(i, p);
	for (int k = 0; k < p; k++)
	  {
	    M.Index(i, k) = X.Index(k);
	    M.Value(i, k) = X.Value(k);
	  }	
      }
    
    // then modifying last columns
    for (int j = i+1; j < M.GetN(); j++)
      {
	while ( (p < X.GetM()) && (X.Index(p) < j))
	  p++;
	
	bool present_X = false;
	if ( (p < X.GetM()) && (X.Index(p) == j))	  
	  {
	    present_X = true;
	    val = X.Value(p);
	  }
	
	int size_col = M.GetColumnSize(j);
	bool present_val = false;
	int k = 0;
	while ( (k < size_col) && (M.Index(j, k) < i))
	  k++;
	
	if ( (k < size_col) && (M.Index(j, k) == i))
	  {
	    if (!present_X)
	      {
		// reducing size of column
		if (size_col > 1)
		  {
		    if (k == size_col-1)
		      M.ResizeColumn(j, size_col-1);
		    else
		      {
			int last_col = M.Index(j, size_col-1);
			val = M.Value(j, size_col-1);
			M.ResizeColumn(j, size_col-1);
			for (int q = k; q < size_col-2; q++)
			  {
			    M.Index(j, q) = M.Index(j, q+1);
			    M.Value(j, q) = M.Value(j, q+1);
			  }
			
			M.Index(j, size_col-2) = last_col;
			M.Value(j, size_col-2) = val;
		      }
		  }
		else
		  M.ClearColumn(j);
	      }
	    else
	      M.Value(j, k) = val;
	    
	    present_val = true;
	  }
	
	if (!present_val && present_X)
	  {
	    // increasing size of column
	    M.ResizeColumn(j, size_col+1);
	    for (int q = size_col; q > k; q--)
	      {
		M.Index(j, q) = M.Index(j, q-1);
		M.Value(j, q) = M.Value(j, q-1);
	      }
	    
	    M.Index(j, k) = i;
	    M.Value(j, k) = val;
	  }

      }
  }
  
  
  // GetRow, SetRow //
  ////////////////////


  ////////////////////
  // GetCol, SetCol //

  
  //! Extracts a column from a matrix
  /*!
    \param M matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, ArrayRowSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int n = M.GetN();
    if (j < 0 || j >= n)
      throw WrongIndex("GetCol()",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif

    int m = M.GetM();

    list<pair<int, T0> > vec;
    for (int i = 0; i < m; i++)
      for (int k = 0; k < M.GetRowSize(i); k++)
	if (M.Index(i, k) == j)
	  vec.push_back(make_pair(i, M.Value(i, k)));
    
    typename list<pair<int, T0> >::iterator it;
    X.Reallocate(vec.size());
    int i = 0;
    for (it = vec.begin(); it != vec.end(); it++)
      {
	X.Index(i) = it->first;
	X.Value(i) = it->second;
	i++;
      }
  }
  
  
  //! Extracts a column from a matrix
  /*!
    \param M matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template<class T0, class Allocator0,
	   class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, ArrayColSparse, Allocator0>& A,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int n = A.GetN();
    if (j < 0 || j >= n)
      throw WrongIndex("GetCol()",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif
    
    int size_col = A.GetColumnSize(j);
    X.Reallocate(size_col);
    for (int k = 0; k < size_col; k++)
      {
	X.Index(k) = A.Index(j, k);
	X.Value(k) = A.Value(j, k);
      }
  }

  
  //! Extracts a column from a sparse matrix
  /*!
    \param M sparse matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
    // symmetric matrix row = col
    GetRow(M, j, X);
  }

  
  //! Extracts a column from a sparse matrix
  /*!
    \param M sparse matrix
    \param j column index
    \param X column extracted
    X = M(:, j)
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
    // symmetric matrix row = col
    GetRow(M, j, X);
  }

  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, ArrayRowSparse, Allocator0>& M)
  {
    T0 val;
    int p = 0;
    for (int i = 0; i < M.GetM(); i++)
      {
	while ( (p < X.GetM()) && (X.Index(p) < i))
	  p++;
	
	bool present_X = false;
	if ( (p < X.GetM()) && (X.Index(p) == i))	  
	  {
	    present_X = true;
	    val = X.Value(p);
	  }
	
	int size_row = M.GetRowSize(i);
	bool present_val = false;
	int k = 0;
	while ( (k < size_row) && (M.Index(i, k) < j))
	  k++;
	

	if ( (k < size_row) && (M.Index(i, k) == j))
	  {
	    if (!present_X)
	      {
		// reducing size of row
		if (size_row > 1)
		  {
		    if (k == size_row-1)
		      M.ResizeRow(j, size_row-1);
		    else
		      {
			int last_row = M.Index(i, size_row-1);
			val = M.Value(i, size_row-1);
			M.ResizeRow(i, size_row-1);
			for (int q = k; q < size_row-2; q++)
			  {
			    M.Index(i, q) = M.Index(i, q+1);
			    M.Value(i, q) = M.Value(i, q+1);
			  }
			
			M.Index(i, size_row-2) = last_row;
			M.Value(i, size_row-2) = val;
		      }
		  }
		else
		  M.ClearRow(i);
	      }
	    else
	      M.Value(i, k) = val;
	    
	    present_val = true;
	  }
	
	if (!present_val && present_X)
	  {
	    // increasing size of row
	    M.ResizeRow(i, size_row+1);
	    for (int q = size_row; q > k; q--)
	      {
		M.Index(i, q) = M.Index(i, q-1);
		M.Value(i, q) = M.Value(i, q-1);
	      }
	    
	    M.Index(i, k) = j;
	    M.Value(i, k) = val;
	  }	
      }
  }
  
  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, ArrayColSparse, Allocator0>& M)
  {
    M.ClearColumn(j);
    M.ReallocateColumn(j, X.GetM());    
    for (int k = 0; k < X.GetM(); k++)
      {
	M.Index(j, k) = X.Index(k);
	M.Value(j, k) = X.Value(k);
      }
  }
  
  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, ArrayRowSymSparse, Allocator0>& M)
  {
    // symmetric matrix, row = column
    SetRow(X, j, M);
  }
   
  
  //! Sets a column of a matrix
  /*!
    \param M matrix
    \param j column index
    \param X new column of M
    M(:, j) = X
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, Symmetric, ArrayColSymSparse, Allocator0>& M)
  {
    // symmetric matrix, row = column
    SetRow(X, j, M);
  }
  
  
  // GetCol, SetCol //
  ////////////////////
  
  
  ////////////
  // MltAdd //


  /*** ArrayRowSymSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const Matrix<T1, Symmetric,
	      ArrayRowSymSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (B.GetM() <= 0)
      return;

    T4 zero;
    SetComplexZero(zero);

    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T3 val;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);

		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = alpha * A.Value(i, k);

		if (p==i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymSparse, Allocator1>& A,
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
	      ArrayRowSymSparse, Allocator1>& A,
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
	      ArrayRowSymSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (B.GetM() <= 0)
      return;

    T4 zero;
    SetComplexZero(zero);

    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T3 val;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = conj(A.Value(i, k));

		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = alpha * conj(A.Value(i, k));

		if (p==i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
  }

  
  /*** ArrayColSymSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const Matrix<T1, Symmetric,
	      ArrayColSymSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (B.GetM() <= 0)
      return;

    T4 zero;
    SetComplexZero(zero);

    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T3 val;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);

		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = alpha * A.Value(i, k);

		if (p==i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayColSymSparse, Allocator1>& A,
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
	      ArrayColSymSparse, Allocator1>& A,
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
	      ArrayColSymSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (B.GetM() <= 0)
      return;

    T4 zero;
    SetComplexZero(zero);

    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T3 val;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = conj(A.Value(i, k));

		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetColumnSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = alpha * conj(A.Value(i, k));

		if (p==i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
  }
  
  
  /*** ArrayRowSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
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
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(i) += val * B(p);
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(i) += alpha * val * B(p);
              }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T4, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
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
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(p) += val * B(i);
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(p) += alpha * val * B(i);
	      }
	  }
      }
  }

  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
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
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = conj(A.Value(i, k));
		C(p) += val * B(i);
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = conj(A.Value(i, k));
		C(p) += alpha * val * B(i);
	      }
	  }
      }
  }
  
  
  /*** ArrayColSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero;
    SetComplexZero(zero);
    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);

    if (alpha == T0(1))
      {
	for (int i = 0 ; i < A.GetN(); i++)
          for (int k = 0; k < A.GetColumnSize(i); k++)
            C(A.Index(i, k)) += A.Value(i, k) * B(i);
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < A.GetN(); i++)
          for (int k = 0; k < A.GetColumnSize(i); k++)
            C(A.Index(i, k)) += alpha * A.Value(i, k) * B(i);
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T4, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero;
    SetComplexZero(zero);
    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);

    if (alpha == T0(1))
      {
	for (int i = 0 ; i < A.GetN(); i++)
          for (int k = 0; k < A.GetColumnSize(i); k++)
            C(i) += A.Value(i, k) * B(A.Index(i, k));
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < A.GetN(); i++)
          for (int k = 0; k < A.GetColumnSize(i); k++)
            C(i) += alpha * A.Value(i, k) * B(A.Index(i, k));
      }
  }

  
  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    T4 zero;
    SetComplexZero(zero);
    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);

    if (alpha == T0(1))
      {
	for (int i = 0 ; i < A.GetN(); i++)
          for (int k = 0; k < A.GetColumnSize(i); k++)
            C(i) += conj(A.Value(i, k)) * B(A.Index(i, k));
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < A.GetN(); i++)
          for (int k = 0; k < A.GetColumnSize(i); k++)
            C(i) += alpha * conj(A.Value(i, k)) * B(A.Index(i, k));
      }
  }

  
  // MltAdd //
  ////////////



  /////////
  // Add //


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowSparse, Allocator2>& B)
  {
    int m = B.GetM(),n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.Value(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetIndex(i), value.GetData());
      }
  }

  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, General, ArrayColSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayColSparse, Allocator2>& B)
  {
    int m = B.GetN(), n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetColumnSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.Value(i, j);
        
	Mlt(alpha, value);
	B.AddInteractionColumn(i, n, A.GetIndex(i), value.GetData());
      }
  }

  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B)
  {
    int m = B.GetM(),n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.Value(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetIndex(i), value.GetData());
      }
  }
  
  
  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, ArrayColSymSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayColSymSparse, Allocator2>& B)
  {
    int m = B.GetN(), n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetColumnSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.Value(i, j);
        
	Mlt(alpha, value);
	B.AddInteractionColumn(i, n, A.GetIndex(i), value.GetData());
      }
  }


  // C = C + complex(A,B)
  template<class T0, class T1, class T2, class T3, class Allocator1,
	   class Allocator2, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	   const Matrix<T2, General, ArrayRowSparse, Allocator2>& B,
	   Matrix<complex<T3>, General, ArrayRowSparse, Allocator3>& C)
  {
    int m = B.GetM(),n1,n2,size_row;;
    Vector<complex<T3>, VectFull, Allocator3> val_row;
    IVect ind_row;
    for (int i = 0 ; i < m ; i++)
      {
	n1 = A.GetRowSize(i);
	n2 = B.GetRowSize(i);
	size_row = n1 + n2;
	val_row.Reallocate(size_row);
	ind_row.Reallocate(size_row);
	for (int j = 0 ; j < n1 ; j++)
	  {
	    ind_row(j) = A.Index(i, j);
	    val_row(j) = alpha*complex<T3>(A.Value(i, j), 0);
	  }

	for (int j = 0 ; j < n2 ; j++)
	  {
	    ind_row(j+n1) = B.Index(i, j);
	    val_row(j+n1) = alpha * complex<T3>(B.Value(i, j));
	  }
        
	C.AddInteractionRow(i, size_row, ind_row, val_row);
      }
  }


  // C = C + complex(A,B)
  template<class T0, class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
	   const Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B,
	   Matrix<complex<T3>, Symmetric, ArrayRowSymSparse, Allocator3>& C)
  {
    int m = B.GetM(), n1, n2, size_row;
    Vector<complex<T3>, VectFull, Allocator3> val_row;
    IVect ind_row;
    for (int i = 0 ; i < m ; i++)
      {
	n1 = A.GetRowSize(i);
	n2 = B.GetRowSize(i);
	size_row = n1 + n2;
	val_row.Reallocate(size_row);
	ind_row.Reallocate(size_row);
	for (int j = 0 ; j < n1 ; j++)
	  {
	    ind_row(j) = A.Index(i, j);
	    val_row(j) = alpha * complex<T3>(A.Value(i, j), 0);
	  }

	for (int j = 0 ; j < n2 ; j++)
	  {
	    ind_row(j+n1) = B.Index(i, j);
	    val_row(j+n1) = alpha * complex<T3>(B.Value(i, j));
	  }

	C.AddInteractionRow(i, size_row, ind_row, val_row);
      }
  }

  
  // Add //
  /////////



  /////////
  // Mlt //


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, Matrix<T, General, ArrayRowSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, Matrix<T, General, ArrayColSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, Symmetric, ArrayColSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  // Matrix-matrix product (sparse matrix against full matrix)
  template<class T1, class Allocator1, class T2, class Prop2,
	   class Storage2, class Allocator2, class T3, class Prop3,
	   class Storage3, class Allocator3>
  void Mlt(const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	   const Matrix<T2, Prop2, Storage2, Allocator2>& B,
	   Matrix<T3, Prop3, Storage3, Allocator3>& C)
  {
    int m = A.GetM();
    int n = B.GetN();
    C.Reallocate(m,n);
    T3 val;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	{
	  val = T3(0);
	  for (int ind = 0; ind < A.GetRowSize(i); ind++)
	    {
	      int k = A.Index(i, ind);
	      val += A.Value(i, ind) * B(k, j);
	    }
	  C(i, j) = val;
	}
  }


  // Mlt //
  /////////
  
  
  ///////////
  // Norms //
  
  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        res = max(res, abs(A.Value(i, j)));
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal> sum(A.GetN());
    sum.Fill(Treal(0));
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        sum(A.Index(i, j)) += abs( A.Value(i, j));
    
    return sum.GetNormInf();
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal res(0), sum;
    for (int i = 0; i < A.GetM(); i++)
      {
        sum = Treal(0);
        for (int j = 0; j < A.GetRowSize(i); j++)
          sum += abs(A.Value(i, j));
        
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
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayColSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        res = max(res, abs(A.Value(i, j)));
    
    return res;
  }
  
  
  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayColSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal res(0), sum;
    for (int i = 0; i < A.GetN(); i++)
      {
        sum = Treal(0);
        for (int j = 0; j < A.GetColumnSize(i); j++)
          sum += abs(A.Value(i, j));
        
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
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayColSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal> sum(A.GetM());
    sum.Fill(Treal(0));
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        sum(A.Index(i, j)) += abs(A.Value(i, j));
    
    return sum.GetNormInf();
  }
  
  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        res = max(res, abs(A.Value(i, j)));
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal> sum(A.GetN());
    sum.Fill(Treal(0));
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          sum(A.Index(i, j)) += abs( A.Value(i, j));
          if (A.Index(i, j) != i)
            sum(i) += abs(A.Value(i, j));
        }
    
    return sum.GetNormInf();
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    return Norm1(A);
  }
  
  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  MaxAbs(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A)
  {
    typename ClassComplexType<T>::Treal res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        res = max(res, abs(A.Value(i, j)));
    
    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm1(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Vector<Treal> sum(A.GetN());
    sum.Fill(Treal(0));
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        {
          sum(A.Index(i, j)) += abs( A.Value(i, j));
          if (A.Index(i, j) != i)
            sum(i) += abs(A.Value(i, j));
        }
    
    return sum.GetNormInf();
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Allocator>
  typename ClassComplexType<T>::Treal
  NormInf(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A)
  {
    return Norm1(A);
  }
  
  
  // Norms //
  ///////////
  
  
  ///////////////
  // Transpose //
  
  
  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayRowSparse, Allocator>& A,
                 Matrix<T, General, ArrayRowSparse, Allocator>& B)
  {
    B.Clear();
    
    int m = A.GetM();
    int n = A.GetN();
    Vector<int, VectFull, CallocAlloc<int> > ptr_T(n);
    
    B.Reallocate(n, m);

    // For each column j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    ptr_T.Zero();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        ptr_T(A.Index(i, j))++;
    
    for (int i = 0; i < n; i++)
      B.ReallocateRow(i, ptr_T(i));
    
    // filling matrix B
    ptr_T.Zero();
    for (int i = 0; i < m; i++)
      for (int jp = 0; jp < A.GetRowSize(i); jp++)
    	{
	  int j = A.Index(i, jp);
	  int k = ptr_T(j);
	  ++ptr_T(j);
          B.Value(j, k) = A.Value(i, jp);
          B.Index(j, k) = i;
    	}
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayRowSparse, Allocator>& A)
  {
    Matrix<T, General, ArrayRowSparse, Allocator> Acopy(A);
    Transpose(Acopy, A);
  }
  

  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayColSparse, Allocator>& A,
                 Matrix<T, General, ArrayColSparse, Allocator>& B)
  {
    B.Clear();
    
    int m = A.GetM();
    int n = A.GetN();
    Vector<int, VectFull, CallocAlloc<int> > ptr_T(m);
    
    B.Reallocate(n, m);

    // For each row j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    ptr_T.Zero();
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        ptr_T(A.Index(i, j))++;
    
    for (int i = 0; i < m; i++)
      B.ReallocateColumn(i, ptr_T(i));
    
    // filling matrix B
    ptr_T.Zero();
    for (int i = 0; i < n; i++)
      for (int jp = 0; jp < A.GetColumnSize(i); jp++)
    	{
	  int j = A.Index(i, jp);
	  int k = ptr_T(j);
	  ++ptr_T(j);
          B.Value(j, k) = A.Value(i, jp);
          B.Index(j, k) = i;
    	}
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayColSparse, Allocator>& A)
  {
    Matrix<T, General, ArrayColSparse, Allocator> Acopy(A);
    Transpose(Acopy, A);
  }


  //! Replaces A by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, General, ArrayRowSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        A.Value(i, j) = conj(A.Value(i, j));
  }


  //! Replaces A by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, General, ArrayColSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        A.Value(i, j) = conj(A.Value(i, j));
  }


  //! Replaces A by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        A.Value(i, j) = conj(A.Value(i, j));
  }


  //! Replaces A by its conjugate
  template<class T, class Allocator>
  void Conjugate(Matrix<T, Symmetric, ArrayColSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
        A.Value(i, j) = conj(A.Value(i, j));
  }

  
  // Transpose //
  ///////////////
  
  
} // namespace Seldon

#define SELDON_FILE_FUNCTIONS_MATRIX_ARRAY_CXX
#endif

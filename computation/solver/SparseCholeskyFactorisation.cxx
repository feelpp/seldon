// Copyright (C) 2010 Marc Duruflé
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


#ifndef SELDON_FILE_SPARSE_CHOLESKY_FACTORIZATION_CXX

namespace Seldon
{

  //! Implementation of Cholesky factorization for sparse symmetric matrix.
  /*! This method may be slow for large matrices. For large matrices, it is
    more efficient to use an external library (Cholmod for example)
    \warning The diagonal value is set to invert of diagonal value of true L
    so that no division is performed when SolveCholesky is called.
  */
  template<class T, class Prop, class Allocator>
  void GetCholesky(Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    int n = A.GetN();
    T t, s, fact;
    int j_col, jrow, index_lu, jpos;
    Vector<T> Row_Val(n);
    IVect Index(n), Row_Ind(n);
    Row_Val.Fill(0);
    Row_Ind.Fill(-1);

    Index.Fill(-1);

    // Conversion to unsymmetric matrix.
    Matrix<T, General, ArrayRowSparse, Allocator> B;
    Copy(A, B);

    // A is cleared.
    A.Clear();
    A.Reallocate(n, n);

    // Main loop over rows.
    for (int i_row = 0; i_row < n; i_row++)
      {
	int size_row = B.GetRowSize(i_row);

        // we are separating lower from upper part
	int length_lower = 0, length_upper = 1;
        Row_Ind(i_row) = i_row;
	Row_Val(i_row) = 0.0;
	Index(i_row) = i_row;

        for (int j = 0; j < size_row; j++)
          {
            int k = B.Index(i_row, j);
            t = B.Value(i_row, j);
            if (k < i_row)
              {
                Row_Ind(length_lower) = k;
                Row_Val(length_lower) = t;
                Index(k) = length_lower;
                length_lower++;
              }
            else if (k == i_row)
              Row_Val(i_row) = t;
            else
              {
                jpos = i_row + length_upper;
		Row_Ind(jpos) = k;
		Row_Val(jpos) = t;
		Index(k) = jpos;
		length_upper++;
              }
          }

        B.ClearRow(i_row);

	j_col = 0;
	int length = 0;

	// Previous rows are eliminated.
	while (j_col < length_lower)
	  {
	    jrow = Row_Ind(j_col);

            // We search first element in lower part.
            int k = j_col;
            for (int j = (j_col+1) ; j < length_lower; j++)
	      {
		if (Row_Ind(j) < jrow)
		  {
		    jrow = Row_Ind(j);
		    k = j;
		  }
	      }


	    if (k != j_col)
	      {
                // If k different from j_col, we are exchanging positions.
		int j = Row_Ind(j_col);
		Row_Ind(j_col) = Row_Ind(k);
		Row_Ind(k) = j;

		Index(jrow) = j_col;
		Index(j) = k;

		s = Row_Val(j_col);
		Row_Val(j_col) = Row_Val(k);
	        Row_Val(k) = s;
	      }

	    // Zero out element in row.
	    Index(jrow) = -1;
	    fact = Row_Val(j_col) * A.Value(jrow, 0);

	    // Combines current row and row jrow.
	    for (int k = 1; k < A.GetRowSize(jrow); k++)
	      {
		s = fact * A.Value(jrow, k);
		int j = A.Index(jrow, k);

		jpos = Index(j);
		if (j >= i_row)
		  {
                    // Dealing with upper part.
		    if (jpos == -1)
		      {
			// This is a fill-in element.
                        int i = i_row + length_upper;
			Row_Ind(i) = j;
			Index(j) = i;
			Row_Val(i) = -s;
			length_upper++;
		      }
		    else
		      {
                        // This is not a fill-in element.
			Row_Val(jpos) -= s;
		      }
		  }
		else
		  {
		    // Dealing  with lower part.
		    if (jpos == -1)
		      {
			// This is a fill-in element.
			Row_Ind(length_lower) = j;
			Index(j) = length_lower;
			Row_Val(length_lower) = -s;
			length_lower++;
		      }
		    else
		      {
			// This is not a fill-in element.
			Row_Val(jpos) -= s;
		      }
		  }
	      }

	    // Stores this pivot element
            // (from left to right -- no danger of overlap
            // with the working elements in L (pivots).
	    Row_Val(length) = fact;
	    Row_Ind(length) = jrow;
	    ++length;
	    j_col++;
	  }

        for (int k = 0; k < length_upper; k++)
	  Index(Row_Ind(i_row + k)) = -1;

	// Now we can store the uppert part of row.
	size_row = length_upper;
	A.ReallocateRow(i_row, length_upper);

	// We store inverse of square root of diagonal element of u.
        if (Row_Val(i_row) < 0)
          {
            cout << "Error during Cholesky factorization " << endl;
            cout << "Matrix must be definite positive " << endl;
            cout << "but diagonal element of row " << i_row
                 << "is equal to " << Row_Val(i_row) << endl;

#ifdef SELDON_WITH_ABORT
            abort();
#endif
          }

	A.Value(i_row, 0) = 1.0 / Row_Val(i_row);
        A.Index(i_row, 0) = i_row;
        index_lu = 1;

	// and extra-diagonal terms
	for (int k = i_row + 1; k < i_row + length_upper; k++)
	  {
	    A.Index(i_row, index_lu) = Row_Ind(k);
	    A.Value(i_row, index_lu) = Row_Val(k);
	    index_lu++;
	  }
      }

    // Diagonal of A is replaced by its square root.
    for (int i = 0; i < n; i++)
      {
        A.Value(i, 0) = sqrt(A.Value(i,0));
        // and other elements multiplied by this value.
        for (int k = 1; k < A.GetRowSize(i); k++)
          A.Value(i, k) *= A.Value(i, 0);
      }
  }


  // Resolution of L x = y.
  template<class classTrans,
           class T0, class T1, class Prop, class Storage,
           class Allocator1, class Allocator2>
  void SolveCholesky(const classTrans& TransA,
                     const Matrix<T0, Prop, ArrayRowSymSparse, Allocator1>& A,
                     Vector<T1, Storage, Allocator2>& x)
  {
    int n = A.GetM();
    if (n <= 0)
      return;

    if (TransA.Trans())
      {
        // We solve L^T x = x
        int j;
        for (int i = n - 1; i >= 0; i--)
          {
            T1 val = x(i);
            for (int k = 1; k < A.GetRowSize(i) ; k++)
              {
                j = A.Index(i, k);
                val -= A.Value(i, k) * x(j);
              }

            x(i) = val * A.Value(i, 0);
          }
      }
    else
      {
        // We solve L x = x
        int j;
        for (int i = 0; i < n; i++)
          {
            x(i) *= A.Value(i, 0);
            for (int k = 1; k < A.GetRowSize(i) ; k++)
              {
                j = A.Index(i, k);
                x(j) -= A.Value(i, k) * x(i);
              }
          }
      }
  }


  // Computation of y = L x.
  template<class classTrans,
           class T0, class T1, class Prop, class Storage,
           class Allocator1, class Allocator2>
  void MltCholesky(const classTrans& TransA,
                   const Matrix<T0, Prop, ArrayRowSymSparse, Allocator1>& A,
                   Vector<T1, Storage, Allocator2>& x)
  {
    int n = A.GetM();
    if (n <= 0)
      return;

    if (TransA.Trans())
      {
        // We overwrite x by L^T x
        int j;
        for (int i = 0; i < n; i++)
          {
            T1 val = x(i) / A.Value(i, 0);
            for (int k = 1; k < A.GetRowSize(i) ; k++)
              {
                j = A.Index(i, k);
                val += A.Value(i, k)*x(j);
              }

            x(i) = val;
          }
      }
    else
      {
        // We overwrite x with L x
        int j;
        for (int i = n - 1; i >= 0; i--)
          {
            for (int k = 1; k < A.GetRowSize(i) ; k++)
              {
                j = A.Index(i, k);
                x(j) += A.Value(i, k)*x(i);
              }

            x(i) /= A.Value(i, 0);
          }
      }
  }

}

#define SELDON_FILE_SPARSE_CHOLESKY_FACTORIZATION_CXX
#endif

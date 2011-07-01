// Copyright (C) 2010 Vivien Mallet
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


#ifndef SELDON_FILE_COMPUTATION_SPARSESOLVER_CXX


#include "SparseSolverInline.cxx"

namespace Seldon
{
  
  /*************************
   * Default Seldon solver *
   *************************/
  
  
  template<class T, class Allocator>
  template<class T0, class Storage0, class Allocator0>
  void SparseSeldonSolver<T, Allocator>::
  FactorizeMatrix(const IVect& perm,
		  Matrix<T0, General, Storage0, Allocator0>& mat,
		  bool keep_matrix)
  {
    IVect inv_permutation;
    
    // We convert matrix to unsymmetric format.
    Copy(mat, mat_unsym);

    // Old matrix is erased if needed.
    if (!keep_matrix)
      mat.Clear();
    
    // We keep permutation array in memory, and check it.
    int n = mat_unsym.GetM();
    if (perm.GetM() != n)
      {
        cout << "Numbering array should have the same size as matrix.";
        cout << endl;
        abort();
      }
    
    permutation_row.Reallocate(n);
    permutation_col.Reallocate(n);
    inv_permutation.Reallocate(n);
    inv_permutation.Fill(-1);
    for (int i = 0; i < n; i++)
      {
        permutation_row(i) = i;
        permutation_col(i) = i;
        inv_permutation(perm(i)) = i;
      }
    
    for (int i = 0; i < n; i++)
      if (inv_permutation(i) == -1)
        {
          cout << "Error in the numbering array." << endl;
          abort();
        }
    
    IVect iperm = inv_permutation;
    
    // Rows of matrix are permuted.
    ApplyInversePermutation(mat_unsym, perm, perm);
    
    // Temporary vector used for solving.
    xtmp.Reallocate(n);

    // Factorization is performed.
    // Columns are permuted during the factorization.
    symmetric_matrix = false;
    inv_permutation.Fill();
    GetLU(mat_unsym, permutation_col, inv_permutation, permtol, print_level);

    // Combining permutations.
    IVect itmp = permutation_col;
    for (int i = 0; i < n; i++)
      permutation_col(i) = iperm(itmp(i));

    permutation_row = perm;

  }
  
  
  template<class T, class Allocator> 
  template<class T0, class Storage0, class Allocator0>
  void SparseSeldonSolver<T, Allocator>::
  FactorizeMatrix(const IVect& perm,
		  Matrix<T0, Symmetric, Storage0, Allocator0>& mat,
		  bool keep_matrix)
  {
    IVect inv_permutation;

    // We convert matrix to symmetric format.
    Copy(mat, mat_sym);

    // Old matrix is erased if needed.
    if (!keep_matrix)
      mat.Clear();

    // We keep permutation array in memory, and check it.
    int n = mat_sym.GetM();
    if (perm.GetM() != n)
      {
        cout << "Numbering array should have the same size as matrix.";
        cout << endl;
        abort();
      }

    permutation_row.Reallocate(n);
    inv_permutation.Reallocate(n);
    inv_permutation.Fill(-1);
    for (int i = 0; i < n; i++)
      {
        permutation_row(i) = perm(i);
        inv_permutation(perm(i)) = i;
      }

    for (int i = 0; i < n; i++)
      if (inv_permutation(i) == -1)
        {
          cout << "Error in the numbering array." << endl;
          abort();
        }

    // Matrix is permuted.
    ApplyInversePermutation(mat_sym, perm, perm);

    // Temporary vector used for solving.
    xtmp.Reallocate(n);

    // Factorization is performed.
    symmetric_matrix = true;
    GetLU(mat_sym, print_level);
  }
  
  
  template<class T, class Allocator> template<class Vector1>
  void SparseSeldonSolver<T, Allocator>::Solve(Vector1& z)
  {
    if (symmetric_matrix)
      {
	for (int i = 0; i < z.GetM(); i++)
	  xtmp(permutation_row(i)) = z(i);
	
	SolveLU(mat_sym, xtmp);
	
	for (int i = 0; i < z.GetM(); i++)
	  z(i) = xtmp(permutation_row(i));
      }
    else
      {
	for (int i = 0; i < z.GetM(); i++)
	  xtmp(permutation_row(i)) = z(i);
	
	SolveLU(mat_unsym, xtmp);
	
	for (int i = 0; i < z.GetM(); i++)
	  z(permutation_col(i)) = xtmp(i);
      }
  }
  
  
  template<class T, class Allocator> template<class TransStatus, class Vector1>
  void SparseSeldonSolver<T, Allocator>::Solve(const TransStatus& TransA, Vector1& z)
  {
    if (symmetric_matrix)
      Solve(z);
    else
      {
	if (TransA.Trans())
	  {
	    for (int i = 0; i < z.GetM(); i++)
	      xtmp(i) = z(permutation_col(i));
	    
	    SolveLU(SeldonTrans, mat_unsym, xtmp);
	    
	    for (int i = 0; i < z.GetM(); i++)
	      z(i) = xtmp(permutation_row(i));
	  }
	else
	  Solve(z);
      }
  }
  
  
  /************************************************
   * GetLU and SolveLU used by SeldonSparseSolver *
   ************************************************/
  
  
  //! LU factorisation with pivot of unsymmetric matrix
  template<class T, class Allocator>
  void GetLU(Matrix<T, General, ArrayRowSparse, Allocator>& A,
	     IVect& iperm, IVect& rperm, 
	     double permtol, int print_level)
  {
    int n = A.GetN();
    
    T fact, s, t;
    double tnorm, zero = 0.0;
    int length_lower, length_upper, jpos, jrow, i_row, j_col;
    int i, j, k, length, size_row, index_lu;
    
    T czero, cone;
    SetComplexZero(czero);
    SetComplexOne(cone);
    Vector<T, VectFull, Allocator> Row_Val(n);
    IVect Index(n), Row_Ind(n), Index_Diag(n);
    Row_Val.Fill(czero);
    Row_Ind.Fill(-1);
    Index_Diag.Fill(-1);

    Index.Fill(-1);
    // main loop
    int new_percent = 0, old_percent = 0;
    for (i_row = 0; i_row < n; i_row++)
      {
        // Progress bar if print level is high enough.
        if (print_level > 0)
          {
            new_percent = int(double(i_row)/(n-1)*80);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#"; cout.flush();
              }

            old_percent = new_percent;
          }

	size_row = A.GetRowSize(i_row);
	tnorm = zero;
	
	// tnorm is the sum of absolute value of coefficients of row i_row.
	for (k = 0 ; k < size_row; k++)
          if (A.Value(i_row, k) != czero)
            tnorm += abs(A.Value(i_row, k));

	if (tnorm == zero)
	  {
            cout << "Structurally singular matrix." << endl;
            cout << "Norm of row " << i_row << " is equal to 0." << endl;
            abort();
          }

        // Unpack L-part and U-part of row of A.
	length_upper = 1;
	length_lower = 0;
	Row_Ind(i_row) = i_row;
	Row_Val(i_row) = czero;
	Index(i_row) = i_row;

	for (j = 0; j < size_row; j++)
	  {
            k = rperm(A.Index(i_row, j));

	    t = A.Value(i_row,j);
	    if (k < i_row)
	      {
		Row_Ind(length_lower) = k;
		Row_Val(length_lower) = t;
		Index(k) = length_lower;
		++length_lower;
	      }
	    else if (k == i_row)
	      {
		Row_Val(i_row) = t;
	      }
	    else
	      {
		jpos = i_row + length_upper;
		Row_Ind(jpos) = k;
		Row_Val(jpos) = t;
		Index(k) = jpos;
		length_upper++;
	      }
          }

	j_col = 0;
	length = 0;

	// Eliminates previous rows.
	while (j_col < length_lower)
	  {
	    // In order to do the elimination in the correct order, we must
            // select the smallest column index.
	    jrow = Row_Ind(j_col);
	    k = j_col;

	    // Determine smallest column index.
	    for (j = j_col + 1; j < length_lower; j++)
	      {
		if (Row_Ind(j) < jrow)
		  {
		    jrow = Row_Ind(j);
		    k = j;
		  }
              }

	    if (k != j_col)
	      {
		// Exchanging column coefficients.
		j = Row_Ind(j_col);
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

            // Gets the multiplier for row to be eliminated (jrow).
	    // first_index_upper points now on the diagonal coefficient.
	    fact = Row_Val(j_col) * A.Value(jrow, Index_Diag(jrow));
	    
	    // Combines current row and row jrow.
	    for (k = (Index_Diag(jrow)+1); k < A.GetRowSize(jrow); k++)
	      {
		s = fact * A.Value(jrow,k);
		j = rperm(A.Index(jrow,k));
		
		jpos = Index(j);
		
		if (j >= i_row)
		  {
		    
		    // Dealing with upper part.
		    if (jpos == -1)
		      {
			// This is a fill-in element.
			i = i_row + length_upper;
			Row_Ind(i) = j;
			Index(j) = i;
			Row_Val(i) = -s;
			++length_upper;
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
			// this is a fill-in element
			Row_Ind(length_lower) = j;
			Index(j) = length_lower;
			Row_Val(length_lower) = -s;
			++length_lower;
		      }
		    else
		      {
			// This is not a fill-in element.
			Row_Val(jpos) -= s;
		      }
		  }
	      }
	    
	    // Stores this pivot element from left to right -- no danger
	    // of overlap with the working elements in L (pivots).
	    Row_Val(length) = fact;
	    Row_Ind(length) = jrow;
	    ++length;	    
	    j_col++;
	  }
	
	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          Index(Row_Ind(i_row+k )) = -1;
	
	size_row = length;
	A.ReallocateRow(i_row,size_row);

        // store L-part
	index_lu = 0;
        for (k = 0 ; k < length ; k++)
          {
            A.Value(i_row,index_lu) = Row_Val(k);
            A.Index(i_row,index_lu) = iperm(Row_Ind(k));
            ++index_lu;
          }

	// Saves pointer to beginning of row i_row of U.
	Index_Diag(i_row) = index_lu;

        // Updates. U-matrix -- first apply dropping strategy.
	length = 0;
	for (k = 1; k <= (length_upper-1); k++)
	  {
	    ++length;
	    Row_Val(i_row+length) = Row_Val(i_row+k);
	    Row_Ind(i_row+length) = Row_Ind(i_row+k);
	  }

	length++;

	// Determines next pivot.
        int imax = i_row;
        double xmax = abs(Row_Val(imax));
        double xmax0 = xmax;
        for ( k = i_row + 1; k <= i_row + length - 1; k++)
          {
            tnorm = abs(Row_Val(k));
            if ((tnorm > xmax) && (tnorm*permtol > xmax0))
              {
                imax = k;
                xmax = tnorm;
              }
          }
	
        // Exchanges Row_Val.
        s = Row_Val(i_row);
        Row_Val(i_row) = Row_Val(imax);
        Row_Val(imax) = s;

        // Updates iperm and reverses iperm.
        j = Row_Ind(imax);
        i = iperm(i_row);
        iperm(i_row) = iperm(j);
        iperm(j) = i;

        // Reverses iperm.
        rperm(iperm(i_row)) = i_row;
        rperm(iperm(j)) = j;

	// Copies U-part in original coordinates.
        int index_diag = index_lu;
	A.ResizeRow(i_row, size_row+length);
	
        for (k = i_row ; k <= i_row + length - 1; k++)
          {
            A.Index(i_row,index_lu) = iperm(Row_Ind(k));
            A.Value(i_row,index_lu) = Row_Val(k);
            ++index_lu;
          }

        // Stores inverse of diagonal element of u.
	A.Value(i_row, index_diag) = cone / Row_Val(i_row);

      } // end main loop

    if (print_level > 0)
      cout << endl;

    if (print_level > 0)
      cout << "The matrix takes " <<
        int((A.GetDataSize()*(sizeof(T)+4))/(1024*1024)) << " MB" << endl;

    for (i = 0; i < n; i++ )
      for (j = 0; j < A.GetRowSize(i); j++)
        A.Index(i,j) = rperm(A.Index(i,j));
    
  }
  
  
  template<class cplx,
	   class Allocator, class Storage2, class Allocator2>
  void SolveLU(const Matrix<cplx, General, ArrayRowSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x)
  {
    SolveLU(SeldonNoTrans, A, x);
  }
  
  
  //! Resolution of LU y = x (x is overwritten with y)
  /*!  L and U are assumed to be stored in A. The diagonal of A contains the
    inverse of the diagonal of U.
  */
  template<class cplx, class TransStatus,
	   class Allocator, class Storage2, class Allocator2>
  void SolveLU(const TransStatus& transA,
               const Matrix<cplx, General, ArrayRowSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x)
  {
    int i, k, n, k_;
    cplx inv_diag;
    n = A.GetM();

    if (transA.Trans())
      {
	// Forward solve (with U^T).
	for (i = 0 ; i < n ; i++)
	  {
	    k_ = 0; k = A.Index(i,k_);
	    while ( k < i)
	      {
		k_++;
		k = A.Index(i,k_);
	      }
	    
	    x(i) *= A.Value(i,k_);
	    for ( k = (k_+1); k < A.GetRowSize(i) ; k++)
	      {
		x(A.Index(i,k)) -= A.Value(i,k) * x(i);
	      }
	  }
	
	// Backward solve (with L^T).
	for (i = n-1 ; i>=0  ; i--)
	  {
	    k_ = 0; k = A.Index(i,k_);
	    while ( k < i)
	      {
		x(k) -= A.Value(i,k_)*x(i);
		k_++;
		k = A.Index(i,k_);
	      }
	  }
      }
    else
      {
	// Forward solve.
	for (i = 0; i < n; i++)
	  {
	    k_ = 0; k = A.Index(i,k_);
	    while ( k < i)
	      {
		x(i) -= A.Value(i,k_)*x(k);
		k_++;
		k = A.Index(i,k_);
	      }
	  }
	
	// Backward solve.
	for (i = n-1; i >= 0; i--)
	  {
	    k_ = 0; k = A.Index(i,k_);
	    while ( k < i)
	      {
		k_++;
		k = A.Index(i,k_);
	      }
	    
	    inv_diag = A.Value(i,k_);
	    for ( k = (k_+1); k < A.GetRowSize(i) ; k++)
	      x(i) -= A.Value(i,k) * x(A.Index(i,k));
	    
	    x(i) *= inv_diag;
	  }
      }
  }

  
  //! LDLt factorisation without pivot for symmetric matrix.
  template<class T, class Allocator>
  void GetLU(Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A, int print_level)
  {
    int size_row;
    int n = A.GetN();
    double zero = 0.0;
    
    T fact, s, t; double tnorm;
    int length_lower, length_upper, jpos, jrow;
    int i_row, j_col, index_lu, length;
    int i, j, k;


    Vector<T, VectFull, Allocator> Row_Val(n);
    IVect Index(n), Row_Ind(n);
    Row_Val.Zero(); Row_Ind.Fill(-1);

    Index.Fill(-1);

    // We convert A into an unsymmetric matrix.
    Matrix<T, General, ArrayRowSparse, Allocator> B;
    Seldon::Copy(A, B);

    A.Clear();
    A.Reallocate(n, n);

    // Main loop.
    int new_percent = 0, old_percent = 0;
    for (i_row = 0; i_row < n; i_row++)
      {
        // Progress bar if print level is high enough.
        if (print_level > 0)
          {
            new_percent = int(double(i_row)/(n-1)*80);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#"; cout.flush();
              }

            old_percent = new_percent;
          }

        // 1-norm of the row of initial matrix.
	size_row = B.GetRowSize(i_row);
	tnorm = zero;
	for (k = 0 ; k < size_row; k++)
          tnorm += abs(B.Value(i_row, k));

	if (tnorm == zero)
	  {
            cout << "Structurally singular matrix." << endl;
            cout << "Norm of row " << i_row << " is equal to 0." << endl;
            abort();
          }

        // Separating lower part from upper part for this row.
	length_upper = 1;
	length_lower = 0;
	Row_Ind(i_row) = i_row;
	Row_Val(i_row) = 0.0;
	Index(i_row) = i_row;

	for (j = 0; j < size_row; j++)
	  {
	    k = B.Index(i_row,j);
            t = B.Value(i_row,j);
	    if (k < i_row)
	      {
		Row_Ind(length_lower) = k;
		Row_Val(length_lower) = t;
		Index(k) = length_lower;
		++length_lower;
	      }
	    else if (k == i_row)
	      {
		Row_Val(i_row) = t;
	      }
	    else
	      {
		jpos = i_row + length_upper;
		Row_Ind(jpos) = k;
		Row_Val(jpos) = t;
		Index(k) = jpos;
		length_upper++;
	      }
          }

        // This row of B is cleared.
        B.ClearRow(i_row);

	j_col = 0;
	length = 0;


        // We eliminate previous rows.
	while (j_col <length_lower)
	  {
	    // In order to do the elimination in the correct order, we must
            // select the smallest column index.
	    jrow = Row_Ind(j_col);
	    k = j_col;

	    // We determine smallest column index.
	    for (j = (j_col+1) ; j < length_lower; j++)
	      {
		if (Row_Ind(j) < jrow)
		  {
		    jrow = Row_Ind(j);
		    k = j;
		  }
	      }

            // If needed, we exchange positions of this element in
            // Row_Ind/Row_Val so that it appears first.
	    if (k != j_col)
	      {

		j = Row_Ind(j_col);
		Row_Ind(j_col) = Row_Ind(k);
		Row_Ind(k) = j;

		Index(jrow) = j_col;
		Index(j) = k;

		s = Row_Val(j_col);
		Row_Val(j_col) = Row_Val(k);
	        Row_Val(k) = s;
	      }

            // Zero out element in row by setting Index to -1.
	    Index(jrow) = -1;

	    // Gets the multiplier for row to be eliminated.
	    fact = Row_Val(j_col) * A.Value(jrow, 0);
	    
	    // Combines current row and row jrow.
	    for (k = 1; k < A.GetRowSize(jrow); k++)
	      {
		s = fact * A.Value(jrow,k);
		j = A.Index(jrow,k);
		
		jpos = Index(j);
		if (j >= i_row)
		  {
		    
		    // Dealing with upper part.
		    if (jpos == -1)
		      {
			
			// This is a fill-in element.
			i = i_row + length_upper;
			Row_Ind(i) = j;
			Index(j) = i;
			Row_Val(i) = -s;
			++length_upper;
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
			++length_lower;
		      }
		    else
		      {
			// This is not a fill-in element.
			Row_Val(jpos) -= s;
		      }
		  }
	      }
	    
	    // We store this pivot element (from left to right -- no
	    // danger of overlap with the working elements in L (pivots).
	    Row_Val(length) = fact;
	    Row_Ind(length) = jrow;
	    ++length;
	    j_col++;
	  }

	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          Index(Row_Ind(i_row+k )) = -1;

	// Updating U-matrix
	length = 0;
	for (k = 1; k <= (length_upper-1); k++)
	  {
	    ++length;
	    Row_Val(i_row+length) = Row_Val(i_row+k);
	    Row_Ind(i_row+length) = Row_Ind(i_row+k);
	  }

	length++;

	// Copies U-part in matrix A.
	A.ReallocateRow(i_row, length);
	index_lu = 1;
	for (k = (i_row+1) ; k <= (i_row+length-1) ; k++)
	  {
	    A.Index(i_row,index_lu) = Row_Ind(k);
	    A.Value(i_row,index_lu) = Row_Val(k);
	    ++index_lu;
	  }
		
	// Stores the inverse of the diagonal element of u.
	A.Value(i_row,0) = 1.0 / Row_Val(i_row);
	
      } // end main loop.

    if (print_level > 0)
      cout<<endl;

    // for each row of A, we divide by diagonal value
    for (int i = 0; i < n; i++)
      for (int j = 1; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= A.Value(i,0);

    if (print_level > 0)
      cout << "The matrix takes " <<
        int((A.GetDataSize()*(sizeof(T)+4))/(1024*1024)) << " MB" << endl;

  }
  
  
  //! Resolution of L D L^t y = x (result is overwritten in x)
  /*! The factor L^t is assumed to be stored in matrix A. The diagonal of A is
    equal to the inverse of diagonal D.
  */
  template<class real, class cplx, class Allocator,
           class Storage2, class Allocator2>
  void SolveLU(const Matrix<real, Symmetric, ArrayRowSymSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x)
  {
    int n = A.GetM(), j_row;
    cplx tmp;

    // We solve first L y = b.
    for (int i_col = 0; i_col < n ; i_col++)
      {
	for (int k = 1; k < A.GetRowSize(i_col) ; k++)
	  {
	    j_row = A.Index(i_col, k);
	    x(j_row) -= A.Value(i_col, k)*x(i_col);
	  }
      }

    // Inverting by diagonal D.
    for (int i_col = 0; i_col < n ; i_col++)
      x(i_col) *= A.Value(i_col,0);

    // Then we solve L^t x = y.
    for (int i_col = n-1; i_col >=0 ; i_col--)
      {
	tmp = x(i_col);
	for (int k = 1; k < A.GetRowSize(i_col) ; k++)
	  {
	    j_row = A.Index(i_col,k);
	    tmp -= A.Value(i_col,k)*x(j_row);
	  }
	
	x(i_col) = tmp;
      }
  }
  
  
  /********************************************
   * GetLU and SolveLU for SeldonSparseSolver *
   ********************************************/
  
  
  template<class T, class Storage, class Allocator, class Alloc2>
  void GetLU(Matrix<T, Symmetric, Storage, Allocator>& A,
	     SparseSeldonSolver<T, Alloc2>& mat_lu,
	     bool keep_matrix)
  {
    // identity ordering
    IVect perm(A.GetM());
    perm.Fill();
    	
    // factorisation
    mat_lu.FactorizeMatrix(perm, A, keep_matrix);
  }


  template<class T, class Storage, class Allocator, class Alloc2>
  void GetLU(Matrix<T, General, Storage, Allocator>& A,
	     SparseSeldonSolver<T, Alloc2>& mat_lu,
	     bool keep_matrix)
  {
    // identity ordering
    IVect perm(A.GetM());
    perm.Fill();
    
    // factorisation
    mat_lu.FactorizeMatrix(perm, A, keep_matrix);
  }
  

  template<class T, class Alloc2, class Allocator>
  void SolveLU(SparseSeldonSolver<T, Alloc2>& mat_lu,
	       Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  template<class T, class Alloc2, class Allocator, class Transpose_status>
  void SolveLU(const Transpose_status& TransA,
	       SparseSeldonSolver<T, Alloc2>& mat_lu,
	       Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }
  
  
  /******************
   * Ordering stuff *
   ******************/
  
  
  //! Construct reverse Cuthill-McKee ordering from a given matrix
  template<class T, class Prop, class Storage, class Allocator,
	   class Tint, class Alloc>
  void FindReverseCuthillMcKeeOrdering(const Matrix<T, Prop,
				       Storage, Allocator>& A,
				       Vector<Tint, VectFull, Alloc>& num)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	num.Clear();
	return;
      }
    
    // pattern of A+A' is retrieved in CSC format
    Vector<Tint, VectFull, Alloc> Ptr, Ind;
    Vector<T, VectFull, Allocator> Value;
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Value, true);
    Value.Clear();
    
    // allocation of arrays
    Vector<Tint, VectFull, Alloc> vertex(n), neighbor(n);
    Vector<bool> RowUsed(n);
    IVect nb_neighbor(n);
    
    num.Reallocate(n);
    num.Fill(0);
    RowUsed.Fill(false);
    vertex.Fill(0);
    neighbor.Fill(0);
    nb_neighbor.Fill(0);
    
    // beginning with first row
    int nb_row = 1;
    int nb_active_row = 1;
    vertex(0) = 0;
    num(0) = 0;
    RowUsed(0) = true;
    int first_free_row = 1;
    
    // loop until all rows are found
    while (nb_row < n)
      {
	// searching all neighboring rows
	int nb = 0;
	for (int i = 0; i < nb_active_row; i++)
	  {
	    int irow = vertex(i);
	    for (int j = Ptr(irow); j < Ptr(irow+1); j++)
	      {
		int icol = Ind(j);
		if (!RowUsed(icol))
		  {
		    RowUsed(icol) = true;
		    neighbor(nb) = icol;
		    nb_neighbor(nb) = Ptr(icol+1) - Ptr(icol);
		    nb++;
		  }
	      }
	  }
	
	if (nb == 0)
	  {
	    // no row has been found, therefore there are independent
	    // blocks in the matrix, searching next free row
	    while (RowUsed(first_free_row))
	      first_free_row++;
	    
	    RowUsed(first_free_row) = true;
	    num(nb_row) = first_free_row;
	    nb_active_row = 1;
	    vertex(0) = first_free_row;
	    nb_row++;
	  }
	else
	  {
	    // sorting neighboring vertices by the number of adjacent vertices
	    Sort(nb, nb_neighbor, neighbor);
	    
	    // then adding those rows
	    nb_active_row = nb;
	    for (int i = 0; i < nb; i++)
	      {
		vertex(i) = neighbor(i);
		num(nb_row) = neighbor(i);
		nb_row++;
	      }
	  }
      }
    
    // reverting final result
    for (int i = 0; i < n; i++)
      vertex(i) = num(i);
    
    for (int i = 0; i < n; i++)
      num(n-1-i) = vertex(i); 
  }
  
}  // namespace Seldon.


#define SELDON_FILE_COMPUTATION_SPARSESOLVER_CXX
#endif

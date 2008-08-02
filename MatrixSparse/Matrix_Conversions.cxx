// Copyright (C) 2001-2008 Vivien Mallet
// File authors: Vivien Mallet (main part), Marc Duruflé.
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

#ifndef SELDON_FILE_MATRIX_CONVERSIONS_CXX

namespace Seldon
{
  
  /*
    From CSR formats to "Matlab" coordinate format
   */
  
  
  //! conversion from RowSparse to coordinate format
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse, Allocator1>& A,
				    IVect& IndRow, IVect& IndCol, Vector<T, Storage, Allocator2>& Val, int index = 0)
  {
    int m = A.GetM(); int nnz = A.GetDataSize();
    IndRow.Reallocate(nnz); IndCol.Reallocate(nnz); Val.Reallocate(nnz);
    int* ptr = A.GetPtr(); int* ind = A.GetInd(); T* val = A.GetData();
    for (int i = 0 ; i < m ; i++)
      {
	for (int j = ptr[i] ; j< ptr[i+1] ; j++)
	  {
	    IndRow(j) = i + index;
	    IndCol(j) = ind[j] + index;
	    Val(j) = val[j];
	  }
      }
  }
  
  
  //! conversion from ColSparse to coordinate format
  template<class T,class Prop,class Storage,class Allocator1,class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse, Allocator1>& A, IVect& IndRow, IVect& IndCol,
				    Vector<T,Storage,Allocator2>& Val, int index = 0)
  {    
    int n = A.GetN(); int nnz = A.GetDataSize();
    IndCol.Reallocate(nnz); IndRow.Reallocate(nnz); Val.Reallocate(nnz);
    int* ptr = A.GetPtr(); int* ind = A.GetInd(); T* val = A.GetData();
    for (int i = 0 ; i< n ; i++)
      {
	for (int j = ptr[i] ; j< ptr[i+1] ; j++)
	  {
	    IndCol(j) = i + index;
	    IndRow(j) = ind[j] + index;
	    Val(j) = val[j];
	  }
      }
  }
  
  
  //! conversion from RowSymSparse to coordinate format
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymSparse, Allocator1>& A,
				    IVect& IndRow, IVect& IndCol, Vector<T, Storage, Allocator2>& Val, int index = 0)
  {
    int m = A.GetM(); int nnz = A.GetDataSize();
    IndRow.Reallocate(nnz); IndCol.Reallocate(nnz); Val.Reallocate(nnz);
    int* ptr = A.GetPtr(); int* ind = A.GetInd(); T* val = A.GetData();
    for (int i = 0 ; i< m ; i++)
      {
	for (int j = ptr[i] ; j< ptr[i+1] ; j++)
	  {
	    IndRow(j) = i + index;
	    IndCol(j) = ind[j] + index;
	    Val(j) = val[j];
	  }
      }
  }
  
  
  //! conversion from ColSymSparse to coordinate format
  template<class T,class Prop,class Storage,class Allocator1,class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymSparse, Allocator1>& A, IVect& IndRow, IVect& IndCol,
				    Vector<T,Storage,Allocator2>& Val, int index = 0)
  {    
    int n = A.GetN(); int nnz = A.GetDataSize();
    IndCol.Reallocate(nnz); IndRow.Reallocate(nnz); Val.Reallocate(nnz);
    int* ptr = A.GetPtr(); int* ind = A.GetInd(); T* val = A.GetData();
    for (int i = 0 ; i< n ; i++)
      {
	for (int j = ptr[i] ; j< ptr[i+1] ; j++)
	  {
	    IndCol(j) = i + index;
	    IndRow(j) = ind[j] + index;
	    Val(j) = val[j];
	  }
      }
  }
  
  
  /*
    From Sparse Array formats to "Matlab" coordinate format
   */
  
  
  //! conversion from ArrayRowSparse to coordinate format
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSparse, Allocator1>& A,
				    IVect& IndRow, IVect& IndCol, Vector<T, Storage, Allocator2>& Val, int index = 0)
  {
    int m = A.GetM(); int nnz = A.GetDataSize();
    // allocating arrays
    IndRow.Reallocate(nnz); IndCol.Reallocate(nnz); Val.Reallocate(nnz);
    int nb = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    IndRow(nb) = i + index;
	    IndCol(nb) = A.Index(i,j) + index;
	    Val(nb) = A.Value(i,j);
	    nb++;
	  }
      }
  }
  
  
  //! conversion from ArrayColSparse to coordinate format
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSparse, Allocator1>& A,
				    IVect& IndRow, IVect& IndCol, Vector<T, Storage, Allocator2>& Val, int index = 0)
  {
    int n = A.GetN(); int nnz = A.GetDataSize();
    // allocating arrays
    IndRow.Reallocate(nnz); IndCol.Reallocate(nnz); Val.Reallocate(nnz);
    int nb = 0;
    for (int i = 0; i < n; i++)
      {
	for (int j = 0; j < A.GetColumnSize(i); j++)
	  {
	    IndRow(nb) = A.Index(i,j) + index;
	    IndCol(nb) = i + index;
	    Val(nb) = A.Value(i,j);
	    nb++;
	  }
      }
  }
  
  
  //! conversion from ArrayRowSparse to coordinate format
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
				    IVect& IndRow, IVect& IndCol, Vector<T, Storage, Allocator2>& Val, int index = 0)
  {
    int m = A.GetM(); int nnz = A.GetDataSize();
    // allocating arrays
    IndRow.Reallocate(nnz); IndCol.Reallocate(nnz); Val.Reallocate(nnz);
    int nb = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    IndRow(nb) = i + index;
	    IndCol(nb) = A.Index(i,j) + index;
	    Val(nb) = A.Value(i,j);
	    nb++;
	  }
      }
  }
  
  
  //! conversion from ArrayColSparse to coordinate format
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymSparse, Allocator1>& A,
				    IVect& IndRow, IVect& IndCol, Vector<T, Storage, Allocator2>& Val, int index = 0)
  {
    int n = A.GetN(); int nnz = A.GetDataSize();
    // allocating arrays
    IndRow.Reallocate(nnz); IndCol.Reallocate(nnz); Val.Reallocate(nnz);
    int nb = 0;
    for (int i = 0; i < n; i++)
      {
	for (int j = 0; j < A.GetColumnSize(i); j++)
	  {
	    IndRow(nb) = A.Index(i,j) + index;
	    IndCol(nb) = i + index;
	    Val(nb) = A.Value(i,j);
	    nb++;
	  }
      }
  }
  
  
  /*
    From CSR to other CSR formats
   */
  
  
  //! B = A
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T,Prop,ColSparse,Alloc1>& A, Matrix<T,Prop,ColSparse,Alloc2>& B)
  {
    int m = A.GetM(), n = A.GetN(), nnz = A.GetDataSize();
    int* ptr_ = A.GetPtr(); int* ind_ = A.GetInd(); T* data_ = A.GetData();
    
    IVect Ptr(n+1), Ind(nnz); Vector<T,Vect_Full,Alloc2> Val(nnz);
    for (int i = 0; i <= n; i++)
      Ptr(i) = ptr_[i];
    
    for (int i = 0; i < nnz; i++)
      {
	Ind(i) = ind_[i];
	Val(i) = data_[i];
      }
    
    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  //! conversion from row-sparse to column-sparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T,Prop,RowSparse,Alloc1>& A, Matrix<T,Prop,ColSparse,Alloc2>& B)
  {
    int m = A.GetM(), n = A.GetN(), nnz = A.GetDataSize();
    int* ptr_ = A.GetPtr(); int* ind_ = A.GetInd(); T* data_ = A.GetData();
    
    // computation of indices of beginning of columns
    IVect Ptr(n+1); Ptr.Fill(0);
    // counting the number of entries per column
    for (int i = 0; i < nnz; i++)
      Ptr(ind_[i])++; 
    
    // incrementing in order to get indices
    int increment = 0, size, num_col;
    for (int i = 0; i < n; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }
    // last index
    Ptr(n) = increment;
    
    // we will use Offset to get current positions of new entries
    IVect Offset = Ptr;
    IVect Ind(nnz); Vector<T,Vect_Full,Alloc2> Val(nnz);
    
    // loop over rows
    for (int i = 0; i < m; i++)
      for (int j = ptr_[i]; j < ptr_[i+1]; j++)
	{
	  num_col = ind_[j];
	  Ind(Offset(num_col)) = i;
	  Val(Offset(num_col)) = data_[j];
	  Offset(num_col)++;
	}
    
    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  //! conversion from row-sparse to column-sparse
  template<class T, class Prop1, class Prop2, class Storage, class Alloc1, class Alloc2>
  void ConvertMatrixSymSparse_to_ColSparse(const Matrix<T,Prop1,Storage,Alloc1>& A, Matrix<T,Prop2,ColSparse,Alloc2>& B)
  {
    int m = A.GetM(), n = A.GetN();
    int* ptr_ = A.GetPtr(); int* ind_ = A.GetInd(); T* data_ = A.GetData();
    
    // computation of indices of beginning of columns
    IVect Ptr(n+1); Ptr.Fill(0);
    // couting the number of entries per column
    for (int i = 0; i < m; i++)
      for (int j = ptr_[i]; j < ptr_[i+1]; j++)
	{
	  Ptr(ind_[j])++;
	  if (i != ind_[j])
	    Ptr(i)++;
	}
    
    // incrementing in order to get indices
    int nnz = 0, size, num_col;
    for (int i = 0; i < n; i++)
      {
	size = Ptr(i);
	Ptr(i) = nnz;
	nnz += size;
      }
    // last index
    Ptr(n) = nnz;
    
    // we will use Offset to get current positions of new entries
    IVect Offset = Ptr;
    IVect Ind(nnz); Vector<T,Vect_Full,Alloc2> Val(nnz);
    
    // loop over rows
    for (int i = 0; i < m; i++)
      for (int j = ptr_[i]; j < ptr_[i+1]; j++)
	{
	  num_col = ind_[j];
	  Ind(Offset(num_col)) = i;
	  Val(Offset(num_col)) = data_[j];
	  Offset(num_col)++;
	  if (i != ind_[j])
	    {
	      Ind(Offset(i)) = num_col;
	      Val(Offset(i)) = data_[j];
	      Offset(i)++;
	    }
	}
    
    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  //! conversion from RowSymSparse to column-sparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T,Prop1,RowSymSparse,Alloc1>& A, Matrix<T,Prop2,ColSparse,Alloc2>& B)
  {
    ConvertMatrixSymSparse_to_ColSparse(A, B);
  }
  
  
  //! conversion from ColSymSparse to column-sparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T,Prop1,ColSymSparse,Alloc1>& A, Matrix<T,Prop2,ColSparse,Alloc2>& B)
  {
    ConvertMatrixSymSparse_to_ColSparse(A, B);
  }
  
  
  /*
    From ArraySparse matrices to CSR matrices
   */
  
  
  //! conversion from ArrayRowSparse to RowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    // matrix (m,n) with nnz entries
    int nnz = mat_array.GetDataSize();
    int m = mat_array.GetM(); int n = mat_array.GetN();
    
    // allocating arrays needed for CSR format
    Vector<T1> Val(nnz);
    IVect IndRow(m+1);
    IVect IndCol(nnz);
    
    // we fill arrays
    int ind = 0;
    IndRow(0) = 0; 
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k<mat_array.GetRowSize(i);k++)
	  {
	    IndCol(ind) = mat_array.Index(i,k);
	    Val(ind) = mat_array.Value(i,k);
	    ind++;
	  }
	IndRow(i+1) = ind;
      }
    
    mat_csr.SetData(m, n, Val, IndRow, IndCol);
  }
  
  
  //! conversion from ArrayRowSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csr)
  {
    // matrix (m,n) with nnz entries
    int nnz = mat_array.GetDataSize();
    int m = mat_array.GetM(); int n = mat_array.GetN();
   
    // conversion in coordinate format
    Vector<T1> Val; IVect IndRow, IndCol;
    ConvertMatrix_to_Coordinates(mat_array, IndRow, IndCol, Val);
    
    // sorting with respect to column numbers
    Sort(IndCol, IndRow, Val);
    
    // constructing pointer array Ptr
    IVect Ptr(n+1);
    
    // counting non-zero entries per column
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i)+1)++;
    
    // accumulation to get pointer array
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i+1) += Ptr(i); 
    
    mat_csr.SetData(m, n, Val, Ptr, IndRow);
  }
  
  
  //! conversion from ArrayRowSymSparse to RowSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    // number of rows and non-zero entries
    int nnz = mat_array.GetDataSize();
    int m = mat_array.GetM();
    
    // allocationr of arrays for CSR format
    Vector<T1> Val(nnz);
    Vector<int,Vect_Full,CallocAlloc<int> > IndRow(m+1);
    Vector<int,Vect_Full,CallocAlloc<int> > IndCol(nnz);
    
    int ind = 0;
    IndRow(0) = 0; 
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < mat_array.GetRowSize(i); k++)
	  {
	    IndCol(ind) = mat_array.Index(i,k);
	    Val(ind) = mat_array.Value(i,k);
	    ind++;
	  }
	IndRow(i+1)=ind;
      }
    
    mat_csr.SetData(m,m,Val,IndRow,IndCol);
  }
  
  
  //! conversion from ArrayRowComplexSparse to RowComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowComplexSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr)
  {
    // non-zero entries (real and imaginary part)
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();
    
    // allocation of arrays for CSR
    Vector<T0> Val_real(nnz_real),Val_imag(nnz_imag);
    IVect IndRow_real(m+1), IndRow_imag(m+1);
    IVect IndCol_real(nnz_real), IndCol_imag(nnz_imag);
    
    int ind_real = 0,ind_imag = 0;
    IndRow_real(0) = 0; IndRow_imag(0) = 0; 
    // loop over rows
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i,k);
	    Val_real(ind_real) = mat_array.ValueReal(i,k);
	    ind_real++;
	  }
	
	IndRow_real(i+1) = ind_real;
	for (int k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i,k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i,k);
	    ind_imag++;
	  }
	
	IndRow_imag(i+1) = ind_imag;
      }
    
    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real, Val_imag, IndRow_imag, IndCol_imag);
  }
  
  
  //! conversion from ArrayRowSymComplexSparse to RowSymComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymComplexSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr)
  {
    // non-zero entries (real and imaginary part)
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();
    
    // allocation of arrays for CSR
    Vector<T0> Val_real(nnz_real),Val_imag(nnz_imag);
    IVect IndRow_real(m+1), IndRow_imag(m+1);
    IVect IndCol_real(nnz_real), IndCol_imag(nnz_imag);
    
    int ind_real = 0,ind_imag = 0;
    IndRow_real(0) = 0; IndRow_imag(0) = 0; 
    // loop over rows
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i,k);
	    Val_real(ind_real) = mat_array.ValueReal(i,k);
	    ind_real++;
	  }
	
	IndRow_real(i+1) = ind_real;
	for (int k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i,k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i,k);
	    ind_imag++;
	  }
	
	IndRow_imag(i+1) = ind_imag;
      }
    
    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real, Val_imag, IndRow_imag, IndCol_imag);
  }  
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse,Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int nnz = A.GetDataSize();
    int n = A.GetM();
    IVect IndRow(nnz),IndCol(nnz);
    Vector<T1> Val(nnz);
    int ind = 0;
    for (int i = 0 ; i < n ; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    if (A.Index(i,j) != i)
	      {
		IndRow(ind) = i;
		IndCol(ind) = A.Index(i,j);
		Val(ind) = A.Value(i,j);
		ind++;
	      }
	  }
      }
    Sort(ind,IndCol,IndRow,Val);
    nnz = ind; ind = 0;
    
    B.Reallocate(n,n);
    for (int i=0; i<n ; i++)
      {
	int first_index = ind;
	while ((ind<nnz)&&(IndCol(ind)<=i))
	  ind++;
	int size_lower = ind-first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeRow(i,size_row);
	ind = first_index;
	for (int j = 0 ; j < size_lower ; j++)
	  {
	    B.Index(i,j) = IndRow(ind);
	    B.Value(i,j) = Val(ind);
	    ind++;
	  }
	for (int j = 0 ; j < size_upper ; j++)
	  {
	    B.Index(i,size_lower+j) = A.Index(i,j);
	    B.Value(i,size_lower+j) = A.Value(i,j);
	  }
	B.AssembleRow(i);
      }
  }

  
} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_MATRIX_CXX
#endif

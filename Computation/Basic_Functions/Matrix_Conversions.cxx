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
  template<class T, class Prop, class Storage, class Allocator1, class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse, Allocator1>& A,
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
  
  template<class T,class Prop,class Storage,class Allocator1,class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse, Allocator1>& A, IVect& IndRow, IVect& IndCol,
				    Vector<T,Storage,Allocator2>& Val, int index)
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
  
  template<class T,class Prop,class Storage,class Allocator1,class Allocator2>
  void ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymSparse, Allocator1>& A, IVect& IndRow, IVect& IndCol,
				    Vector<T,Storage,Allocator2>& Val, int index)
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
  
} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_MATRIX_CXX
#endif

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


#ifndef SELDON_FILE_MATRIX_CONVERSIONS_CXX


#include "Matrix_Conversions.hxx"


namespace Seldon
{

  /*
    From CSR formats to "Matlab" coordinate format.
  */


  //! Conversion from RowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i < m; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndRow(j) = i + index;
	  IndCol(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }


  //! Conversion from ColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int n = A.GetN();
    int nnz = A.GetDataSize();
    IndCol.Reallocate(nnz);
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i< n; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndCol(j) = i + index;
	  IndRow(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }


  //! Conversion from RowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  if (ind[ptr[i]] == i)
	    nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = ind[j] + index;
	      Val(nb) = val[j];
	      Ptr(ind[j])++;
	      nb++;

	      if (ind[j] != i)
		{
		  IndRow(nb) = ind[j] + index;
		  IndCol(nb) = i + index;
		  Val(nb) = val[j];
		  Ptr(i)++;
		  nb++;
		}
	    }

	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ... and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }

      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j< ptr[i + 1]; j++)
	    {
	      IndRow(j) = i + index;
	      IndCol(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }


  //! Conversion from ColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    if (ind[j] == i)
	      nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = ind[j] + index;
	      Val(nb) = val[j];
	      Ptr(ind[j])++;
	      nb++;

	      if (ind[j] != i)
		{
		  IndRow(nb) = ind[j] + index;
		  IndCol(nb) = i + index;
		  Val(nb) = val[j];
		  Ptr(i)++;
		  nb++;
		}
	    }

	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }

      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j< ptr[i + 1]; j++)
	    {
	      IndRow(j) = i + index;
	      IndCol(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }


  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from ArrayRowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int nb = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	{
	  IndRow(nb) = i + index;
	  IndCol(nb) = A.Index(i, j) + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int n = A.GetN();
    int nnz = A.GetDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int nb = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetColumnSize(i); j++)
	{
	  IndRow(nb) = A.Index(i, j) + index;
	  IndCol(nb) = i + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayRowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i, j) == i)
	      nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      Ptr(A.Index(i, j))++;
	      nb++;

	      if (A.Index(i, j) != i)
		{
		  IndRow(nb) = A.Index(i, j) + index;
		  IndCol(nb) = i + index;
		  Val(nb) = A.Value(i, j);
		  Ptr(i)++;
		  nb++;
		}
	    }

        // Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }


  //! Conversion from ArrayColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i, j) == i)
	      nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      Ptr(A.Index(i, j))++;
	      nb++;

	      if (A.Index(i, j) != i)
		{
		  IndRow(nb) = A.Index(i, j) + index;
		  IndCol(nb) = i + index;
		  Val(nb) = A.Value(i, j);
		  Ptr(i)++;
		  nb++;
		}
	    }

	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndRow(nb) = A.Index(i, j) + index;
	      IndCol(nb) = i + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }


  /*
    From "Matlab" coordinate format to CSR formats.
  */


  //! Conversion from coordinate format to RowSparse.
  /*! Contrary to the other conversion functions
    ConvertMatrix_from_Coordinates, this one accepts duplicates.
    \param[in] IndRow_ row indexes of the non-zero elements.
    \param[in] IndCol_ column indexes of the non-zero elements.
    \param[in] Val values of the non-zero elements.
    \param[out] A matrix defined by \a IndRow, \a IndCol and \a Val.
    \param[in] index index of the first column and of the first row.
  */
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSparse, Allocator3>& A,
				 int index = 0)
  {
    int i, j;

    int Nelement = IndRow_.GetLength();

    Vector<int, VectFull, CallocAlloc<int> > IndRow(Nelement),
      IndCol(Nelement);
    for (int i = 0; i < Nelement; i++)
      {
	IndRow(i) = IndRow_(i) - index;
	IndCol(i) = IndCol_(i) - index;
      }
    IndRow_.Clear();
    IndCol_.Clear();

    Sort(IndRow, IndCol, Val);

    Vector<int, VectFull, CallocAlloc<int> > ptr(A.GetM() + 1);
    for (i = 0; i <= IndRow(0); i++)
      ptr(i) = 0;

    int row = IndRow(0);
    int previous_i = 0;

    for (i = 1; i < Nelement; i++)
      if (IndRow(i) != row)
        {
          for (j = row + 1; j <= IndRow(i); j++)
            ptr(j) = i;

          Sort(previous_i, i-1, IndCol, Val);
          row = IndRow(i);
          previous_i = i;
        }

    Sort(previous_i, Nelement - 1, IndCol, Val);
    for (i = IndRow(Nelement - 1); i < A.GetM(); i++)
      ptr(i + 1) = Nelement;

    A.SetData(A.GetM(), A.GetN(), Val, ptr, IndCol);
  }


  //! Conversion from coordinate format to ColSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSparse, Allocator3>& A,
				 int index = 0)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(n + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i) + 1)++;
      }

    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndRow'
    for (int i = 0; i < n; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndRow, Val);

    A.SetData(m, n, Val, Ptr, IndRow);
  }


  //! Conversion from coordinate format to RowSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymSparse, Allocator3>& A,
				 int index = 0)
  {
    // Assuming there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i) + 1)++;
      }

    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndCol'.
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);

    A.SetData(m, n, Val, Ptr, IndCol);
  }


  //! Conversion from coordinate format to ColSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymSparse, Allocator3>& A,
				 int index = 0)
  {
    // Assuming there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i) + 1)++;
      }

    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndCol'.
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);

    A.SetData(m, n, Val, Ptr, IndCol);
  }


  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from coordinate format to ArrayRowSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSparse,
				 Allocator3>& A,
				 int index = 0)
  {
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Number of elements per row.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateRow(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayColSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSparse,
				 Allocator3>& A,
				 int index = 0)
  {
    // Assuming there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // Sorts array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(n);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < n; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateColumn(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndRow(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayRowSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymSparse,
				 Allocator3>& A,
				 int index = 0)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Clear(); A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
          // sorting column numbers
          Sort(offset, offset+Ptr(i)-1, IndCol, Val);

          // putting values in A
	  A.ReallocateRow(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
            }

	  offset += Ptr(i);
	}
  }


  //! Conversion from coordinate format to ArrayColSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymSparse,
				 Allocator3>& A,
				 int index = 0)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateColumn(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }


  /*
    From CSR to other CSR formats.
  */


  //! B = A.
  template<class T, class Prop, class Storage, class Allocator>
  void Copy(const Matrix<T, Prop, Storage, Allocator>& A,
	    Matrix<T, Prop, Storage, Allocator>& B)
  {
    B = A;
  }


  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    int i, j;

    int m = A.GetM(), n = A.GetN(), nnz = A.GetDataSize();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Computation of the indices of the beginning of columns.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);
    // Counting the number of entries per column.
    for (i = 0; i < nnz; i++)
      Ptr(ind_[i])++;

    // Incrementing in order to get indices.
    int increment = 0, size, num_col;
    for (i = 0; i < n; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }
    // Last index.
    Ptr(n) = increment;

    // 'Offset' will be used to get current positions of new entries.
    Vector<int, VectFull, CallocAlloc<int> > Offset = Ptr;
    Ind.Reallocate(nnz);
    Val.Reallocate(nnz);

    // Loop over rows.
    for (i = 0; i < m; i++)
      for (j = ptr_[i]; j < ptr_[i + 1]; j++)
	{
	  num_col = ind_[j];
	  Ind(Offset(num_col)) = i;
	  Val(Offset(num_col)) = data_[j];
	  Offset(num_col)++;
	}
  }


  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop, RowSparse, Alloc1>& A,
	    Matrix<T, Prop, ColSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr;
    Vector<int, VectFull, CallocAlloc<int> > Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    IVect IndRow, IndCol;
    Vector<T> Val;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);

    int m = A.GetM();
    int n = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();
    // counting number of non-zero entries
    int nnz = 0;
    for (int i = 0; i < IndCol.GetM(); i++)
      if (IndRow(i) <= IndCol(i))
        {
          Ptr(IndRow(i) + 1)++;
          nnz++;
        }

    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);

    // filling Ind and Value
    Ind.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0;
    for (int i = 0; i < IndCol.GetM(); i++)
      if (IndRow(i) <= IndCol(i))
        {
          Ind(nnz) = IndCol(i);
          Value(nnz) = Val(i);
          nnz++;
        }
  }


  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
	    Matrix<T, Prop, RowSymSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;
    Vector<T, VectFull, Alloc2> Value;

    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop1, class Storage, class Tint,
	   class Alloc1, class Alloc2, class Alloc3, class Alloc4>
  void
  ConvertSymmetricToCSC(const Matrix<T, Prop1, Storage, Alloc1>& A,
                        General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                        Vector<Tint, VectFull, Alloc3>& Ind,
                        Vector<T, VectFull, Alloc4>& Val)
  {
    int i, j;

    int m = A.GetM(), n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Computation of the indices of the beginning of columns.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);
    // Counting the number of entries per column.
    for (i = 0; i < m; i++)
      for (j = ptr_[i]; j < ptr_[i + 1]; j++)
	{
	  Ptr(ind_[j])++;
	  if (i != ind_[j])
	    Ptr(i)++;
	}

    // Incrementing in order to get indices.
    int nnz = 0, size, num_col;
    for (i = 0; i < n; i++)
      {
	size = Ptr(i);
	Ptr(i) = nnz;
	nnz += size;
      }
    // Last index.
    Ptr(n) = nnz;

    // 'Offset' will be used to get current positions of new entries.
    Vector<int, VectFull, CallocAlloc<int> > Offset = Ptr;
    Ind.Reallocate(nnz);
    Val.Reallocate(nnz);

    // Loop over rows.
    for (i = 0; i < m; i++)
      for (j = ptr_[i]; j < ptr_[i + 1]; j++)
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
  }


  template<class T, class Prop1, class Tint,
	   class Alloc1, class Alloc2, class Alloc3, class Alloc4>
  void
  ConvertToCSC(const Matrix<T, Prop1, RowSymSparse, Alloc1>& A,
               General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
               Vector<Tint, VectFull, Alloc3>& Ind,
               Vector<T, VectFull, Alloc4>& Val)
  {
    ConvertSymmetricToCSC(A, sym, Ptr, Ind, Val);
  }


  template<class T, class Prop1, class Tint,
	   class Alloc1, class Alloc2, class Alloc3, class Alloc4>
  void
  ConvertToCSC(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
               General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
               Vector<Tint, VectFull, Alloc3>& Ind,
               Vector<T, VectFull, Alloc4>& Val)
  {
    ConvertSymmetricToCSC(A, sym, Ptr, Ind, Val);
  }


  //! Conversion from column-oriented sparse to row-oriented sparse.
  /*!
    \param[in] A matrix to be converted.
    \param[out] B converted matrix.
  */
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void Copy(const Matrix<T1, Prop1, ColSparse, Alloc1>& A,
	    Matrix<T2, Prop2, RowSparse, Alloc2>& B)
  {
    int i, j;

    int m = A.GetM();
    int n = A.GetN();
    int  nnz = A.GetDataSize();

    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T1* data = A.GetData();

    // Computation of the indexes of the beginning of rows.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(n + 1);
    Ptr.Fill(0);
    // Counting the number of entries per row.
    for (i = 0; i < nnz; i++)
      Ptr(ind[i])++;

    // Incrementing in order to get the row indexes.
    int increment = 0, size, num_row;
    for (i = 0; i < n; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }
    // Last index.
    Ptr(n) = increment;

    // 'Offset' will be used to get current positions of new entries.
    Vector<int, VectFull, CallocAlloc<int> > Offset = Ptr;
    Vector<int, VectFull, CallocAlloc<int> > Ind(nnz);
    Vector<T2, VectFull, Alloc2> Val(nnz);

    // Loop over the columns.
    for (j = 0; j < n; j++)
      for (i = ptr[j]; i < ptr[j + 1]; i++)
	{
	  num_row = ind[i];
	  Ind(Offset(num_row)) = j;
	  Val(Offset(num_row)) = data[i];
	  Offset(num_row)++;
	}

    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! Conversion from RowSymSparse to column-sparse.
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, RowSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr;
    Vector<int, VectFull, CallocAlloc<int> > Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! Conversion from ColSymSparse to column-sparse.
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr;
    Vector<int, VectFull, CallocAlloc<int> > Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  /*
    From ArraySparse matrices to CSR matrices.
  */


  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& IndRow,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with 'nnz' entries.
    int nnz = A.GetDataSize();
    int m = A.GetM();

    // Allocating arrays needed for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    // Filling the arrays.
    int ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < A.GetRowSize(i); k++)
	  {
	    IndCol(ind) = A.Index(i, k);
	    Val(ind) = A.Value(i, k);
	    ind++;
	  }
	IndRow(i + 1) = ind;
      }
  }


  //! Conversion from ArrayRowSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int, VectFull, CallocAlloc<int> > IndRow;
    Vector<int, VectFull, CallocAlloc<int> > IndCol;

    General sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Val, IndRow, IndCol);
  }

  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    int i;

    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();

    // Conversion in coordinate format.
    Vector<Tint, VectFull, CallocAlloc<Tint> > IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);

    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);

    // Counting non-zero entries per column.
    for (i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
  }


  //! Conversion from ArrayRowSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csr)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Val;

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    General sym;
    ConvertToCSC(mat_array, sym, Ptr, IndRow, Val);

    mat_csr.SetData(m, n, Val, Ptr, IndRow);
  }


  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& IndRow,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Number of rows and non-zero entries.
    int nnz = A.GetDataSize();
    int m = A.GetM();

    // Allocation of arrays for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    int ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < A.GetRowSize(i); k++)
	  {
	    IndCol(ind) = A.Index(i, k);
	    Val(ind) = A.Value(i, k);
	    ind++;
	  }
	IndRow(i + 1) = ind;
      }
  }


  //! Conversion from ArrayRowSymSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int, VectFull, CallocAlloc<int> > IndRow;
    Vector<int, VectFull, CallocAlloc<int> > IndCol;

    Symmetric sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }


  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& AllVal)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    Vector<T> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }

    Sort(ind, IndCol, IndRow, Val);

    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz + ind);
    AllVal.Reallocate(nnz+ind);
    nnz = ind;
    ind = 0;

    int offset = 0; Ptr(0) = 0;
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;

        int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;

	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    Ind(offset+j) = IndRow(ind);
	    AllVal(offset+j) = Val(ind);
            ind++;
	  }

	for (j = 0; j < size_upper; j++)
	  {
	    Ind(offset + size_lower + j) = A.Index(i, j);
	    AllVal(offset + size_lower + j) = A.Value(i, j);
          }

        offset += size_row; Ptr(i+1) = offset;
      }
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;
    Vector<T1> AllVal;

    int n = A.GetM();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, AllVal);

    B.SetData(n, n, AllVal, Ptr, Ind);
  }


  //! Conversion from ArrayRowComplexSparse to RowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();

    // Allocation of arrays for CSR.
    Vector<T0, VectFull, Allocator1> Val_real(nnz_real), Val_imag(nnz_imag);
    Vector<int, VectFull, CallocAlloc<int> > IndRow_real(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndRow_imag(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndCol_real(nnz_real);
    Vector<int, VectFull, CallocAlloc<int> > IndCol_imag(nnz_imag);

    int ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }

	IndRow_real(i + 1) = ind_real;
	for (k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }

	IndRow_imag(i + 1) = ind_imag;
      }

    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }


  //! Conversion from ArrayRowSymComplexSparse to RowSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0,
       ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();

    // Allocation of arrays for CSR.
    Vector<T0, VectFull, Allocator1> Val_real(nnz_real), Val_imag(nnz_imag);
    Vector<int, VectFull, CallocAlloc<int> > IndRow_real(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndRow_imag(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndCol_real(nnz_real);
    Vector<int, VectFull, CallocAlloc<int> > IndCol_imag(nnz_imag);

    int ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }

	IndRow_real(i + 1) = ind_real;
	for (int k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }

	IndRow_imag(i + 1) = ind_imag;
      }

    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz),IndCol(nnz);
    Vector<T1, VectFull, Allocator1> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;

    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
	int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeRow(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = A.Index(i, j);
	    B.Value(i, size_lower + j) = A.Value(i, j);
	  }
	B.AssembleRow(i);
      }
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int, VectFull, CallocAlloc<int> > IndRow(nnz), IndCol(nnz);
    Vector<T1> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;

    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
	int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeColumn(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = A.Index(i, j);
	    B.Value(i, size_lower + j) = A.Value(i, j);
	  }
	B.AssembleColumn(i);
      }
  }


  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& IndCol,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    int i, k;

    // Matrix (m,n) with 'nnz' entries.
    int nnz = A.GetDataSize();
    int m = A.GetM();
    int n = A.GetN();

    // Allocating arrays needed for CSC format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(n+1);

    // Filling the arrays.
    int ind = 0;
    IndCol(0) = 0;
    for (i = 0; i < n; i++)
      {
	for (k = 0; k < A.GetColumnSize(i); k++)
	  {
	    IndRow(ind) = A.Index(i, k);
	    Val(ind) = A.Value(i, k);
	    ind++;
	  }

	IndCol(i + 1) = ind;
      }
  }


  //! Conversion from ArrayColSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csc)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int, VectFull, CallocAlloc<int> > IndRow;
    Vector<int, VectFull, CallocAlloc<int> > IndCol;

    General sym;
    ConvertToCSC(mat_array, sym, IndCol, IndRow, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();

    mat_csc.SetData(m, n, Val, IndCol, IndRow);
  }


  //! Conversion from ArrayRowSparse to ArrayColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& A,
	    Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i;

    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();

    // Conversion in coordinate format.
    Vector<T1> Val;
    Vector<int, VectFull, CallocAlloc<int> > IndRow, IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);

    // Constructing pointer array 'Ptr'.
    Vector<int, VectFull, CallocAlloc<int> > Ptr(n + 1);

    // Counting non-zero entries per column.
    for (i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    // we fill matrix B
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
	if (size_col > 0)
	  {
	    B.ReallocateColumn(i, size_col);
	    for (int j = Ptr(i); j < Ptr(i+1); j++)
	      {
		B.Index(i, j-Ptr(i)) = IndRow(j);
		B.Value(i, j-Ptr(i)) = Val(j);
	      }
	  }
      }
  }


  /***********************
   * GetSymmetricPattern *
   ***********************/


  template<class T, class Prop, class Storage, class Allocator,
           class Tint, class Allocator2, class Allocator3>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Vector<Tint, VectFull, Allocator2>& Ptr,
                           Vector<Tint, VectFull, Allocator3>& Ind)
  {
    int n = A.GetM();

    // Converting to coordinates.
    Vector<int, VectFull, CallocAlloc<int> > IndRow, IndCol;
    Vector<T> Value;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value);

    // clearing values
    Value.Clear();

    // Sorting columns too.
    Vector<int, VectFull, CallocAlloc<int> > IndRow2, IndCol2, Index(n);
    IndRow2 = IndRow;
    IndCol2 = IndCol;
    Sort(IndCol2.GetM(), IndCol2, IndRow2);

    Tint max_nnz = 0;
    for (int i = 0; i < IndRow.GetM(); i++)
      if (IndRow(i) <= IndCol(i))
        max_nnz++;

    for (int i = 0; i < IndRow.GetM(); i++)
      if (IndCol2(i) <= IndRow2(i))
        max_nnz++;

    // then symmetrization of pattern and conversion to csr.
    Ptr.Reallocate(n+1);
    Ind.Reallocate(max_nnz);
    Tint j_begin = 0, j_end = 0, size_row = 0;
    Tint j2_begin = 0, j2_end = 0;
    Ptr(0) = 0;
    for (int i = 0; i < A.GetM(); i++)
      {
        j_begin = j_end;
        size_row = 0;
        // We retrieve column numbers.
        while ( (j_end < IndRow.GetM()) && (IndRow(j_end) == i))
          {
            if (IndRow(j_end) <= IndCol(j_end))
              {
                Index(size_row) = IndCol(j_end);
                size_row++;
              }

            j_end++;
          }

        j2_begin = j2_end;
        while ( (j2_end < IndCol2.GetM()) && (IndCol2(j2_end) == i))
          {
            if (IndCol2(j2_end) <= IndRow2(j2_end))
              {
                Index(size_row) = IndRow2(j2_end);
                size_row++;
              }

            j2_end++;
          }

        // Sorting indexes.
        Assemble(size_row, Index);

        // Updating Ptr, Ind.
        for (int j = 0; j < size_row; j++)
          Ind(Ptr(i) + j) = Index(j);

        Ptr(i+1) = Ptr(i) + size_row;
      }

    IndRow2.Clear(); IndCol2.Clear();
    IndRow.Clear(); IndCol.Clear();
    Ind.Resize(Ptr(n));
  }


  template<class T, class Prop, class Storage, class Allocator, class AllocI>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Matrix<int, Symmetric, RowSymSparse, AllocI>& B)
  {
    Vector<int, VectFull, CallocAlloc<int> > Ptr, Ind;

    GetSymmetricPattern(A, Ptr, Ind);

    int n = A.GetM();
    Vector<int, VectFull, CallocAlloc<int> > Val(Ptr(n));
    // We put Ptr and Ind into the matrix B.
    B.SetData(n, n, Val, Ptr, Ind);
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_CONVERSIONS_CXX
#endif

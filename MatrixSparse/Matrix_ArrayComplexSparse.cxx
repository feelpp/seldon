// Copyright (C) 2001-2008 Vivien Mallet
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

#ifndef FILE_MATRIX_ARRAY_COMPLEX_SPARSE_CXX

#include "Matrix_ArrayComplexSparse.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/
  
  
  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_ArrayComplexSparse()
    : ind_real(),val_real(),ind_imag(),val_imag()
  {
    nz_real_ = 0;
    nz_imag_ = 0;
    this->m_ = 0;
    this->n_ = 0;
  }
  
  
  //! Constructor.
  /*!
    Builds an empty i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_ArrayComplexSparse(int i, int j):
    ind_real(Storage::GetFirst(i, j)), val_real(Storage::GetFirst(i, j)),
    ind_imag(Storage::GetFirst(i, j)), val_imag(Storage::GetFirst(i, j))
  {
    nz_real_ = 0;
    nz_imag_ = 0;
    this->m_ = Storage::GetFirst(i, j);
    this->n_ = Storage::GetSecond(i, j);
  }
  
  
  /**************
   * DESTRUCTOR *
   **************/
  
  
  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ~Matrix_ArrayComplexSparse()
  {
    this->m_ = 0;
    this->n_ = 0;
  }
  
  
  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix is
    empty (0 by 0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_ArrayComplexSparse();
  }
  
  
  /*********************
   * MEMORY MANAGEMENT *
   *********************/
  
  
  //! Returns used memory in bytes
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetMemorySize() const
  {
    int size = sizeof(int) * (12 + 2*ind_real.GetM() + 2*val_real.GetM()
			      + 2*ind_imag.GetM() + 2*val_imag.GetM()
			      + GetRealDataSize() + GetImagDataSize())
      + sizeof(T) * (GetRealDataSize() + GetImagDataSize());
    
    return size;
  }
  
  
  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Data is lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  Reallocate(int i, int j)
  {
    // Clears previous entries.
    Clear();
    
    this->m_ = i;
    this->n_ = j;
    
    int n = Storage::GetFirst(i,j);
    val_real.Reallocate(n);
    ind_real.Reallocate(n);
    val_imag.Reallocate(n);
    ind_imag.Reallocate(n);
    nz_real_ = 0;
    nz_imag_ = 0;
  }
  
  
  //! Reallocates additional memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    \note Data is kept.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  Resize(int i, int j)
  {
    if ((i != this->m_) || (j != this->n_))
      {
	int n = Storage::GetFirst(m_, n_);
	int new_n = Storage::GetFirst(i, j);
	Vector<IVect, Vect_Full, NewAlloc<IVect> > new_ind_real;
	Vector<vect_value, Vect_Full, NewAlloc<vect_value> > new_val_real;
	Vector<IVect, Vect_Full, NewAlloc<IVect> > new_ind_imag;
	Vector<vect_value, Vect_Full, NewAlloc<vect_value> > new_val_imag;
	new_ind_real.Reallocate(new_n);
	new_val_real.Reallocate(new_n);
	new_ind_imag.Reallocate(new_n);
	new_val_imag.Reallocate(new_n);
	
	nz_real_ = 0;
	nz_imag_ = 0;
	for (int k = 0; k < min(n, new_n); k++)
	  {
	    new_ind_real(k).SetData(ind_real(k).GetM(),
				    ind_real(k).GetData());
	    ind_real(k).Nullify();
	    new_val_real(k).SetData(val_real(k).GetM(),
				    val_real(k).GetData());
	    val_real(k).Nullify();
	    nz_real_ += new_ind_real(k).GetM();
	    
	    new_ind_imag(k).SetData(ind_imag(k).GetM(),
				    ind_imag(k).GetData());
	    ind_imag(k).Nullify();
	    new_val_imag(k).SetData(val_imag(k).GetM(),
				    val_imag(k).GetData());
	    val_imag(k).Nullify();
	    nz_imag_ += new_ind_imag(k).GetM();
	  }
	
	ind_real.Clear();
	ind_real.SetData(new_ind_real.GetM(), new_ind_real.GetData());
	new_ind_real.Nullify();
	val_real.Clear();
	val_real.SetData(new_val_real.GetM(), new_val_real.GetData());
	new_val_real.Nullify();
	
	ind_imag.Clear();
	ind_imag.SetData(new_ind_imag.GetM(), new_ind_imag.GetData());
	new_ind_imag.Nullify();
	val_imag.Clear();
	val_imag.SetData(new_val_imag.GetM(), new_val_imag.GetData());
	new_val_imag.Nullify();
	
	this->m_ = i;
	this->n_ = j;
      }
  }
  
  
  //! Changes the size of i-th row of the matrix (real part of the row).
  /*!
    \param[in] i row number.
    \param[in] j new size of row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ReallocateRealRow(int i, int j)
  {
    if (j != this->GetRealRowSize(i))
      {
	nz_real_ += j - GetRealRowSize(i);
	val_real(i).Reallocate(j);
	ind_real(i).Reallocate(j);
      }
  }
  
  
  //! Changes the size of i-th row of the matrix (imaginary part of the row).
  /*!
    \param[in] i row number.
    \param[in] j new size of row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ReallocateImagRow(int i, int j)
  {
    if (j != this->GetImagRowSize(i))
      {
	nz_imag_ += j - GetImagRowSize(i);
	val_imag(i).Reallocate(j);
	ind_imag(i).Reallocate(j);
      }
  }
  
  
  //! chAnges the size of i-th row of the matrix (real part of the row).
  /*!
    \param[in] i row number.
    \param[in] j new size of row.
    \note Data is kept.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ResizeRealRow(int i, int j)
  {
    int m = this->GetRealRowSize(i);
    if (j != m)
      {
	nz_real_ += j - m;
	pointer data_real_;
	data_real_ = allocator_.allocate(j, this);
	IVect Ind_New(j);
	
	for (int k = 0; k < min(m,j); k++)
	  {
	    data_real_[k] = val_real(i)(k);
	    Ind_New(k) = ind_real(i)(k);
	  }
	
	ind_real(i).SetData(j, Ind_New.GetData());
	Ind_New.Nullify();
	val_real(i).SetData(j, data_real_);
      }
  }
  
  
  //! Changes the size of i-th row of the matrix (imaginary part of the row).
  /*!
    \param[in] i row number.
    \param[in] j new size of row.
    \note Data is kept.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ResizeImagRow(int i, int j)
  {
    int m = this->GetImagRowSize(i);
    if (j != m)
      {
	nz_imag_ += j - m;
	pointer data_imag_;
	data_imag_ = allocator_.allocate(j, this);
	IVect Ind_New(j);
	
	for (int k = 0; k < min(m, j); k++)
	  {
	    data_imag_[k] = val_imag(i)(k);
	    Ind_New(k) = ind_imag(i)(k);
	  }
	
	ind_imag(i).SetData(j, Ind_New.GetData());
	Ind_New.Nullify();
	val_imag(i).SetData(j, data_imag_);
      }
  }
  
  
  //! Real part of row i is removed.
  /*!
    \param[in] i row number.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ClearRealRow(int i)
  {
    int l = this->GetRealRowSize(i);
    if (l > 0)
      {
	nz_real_ -= l;
	val_real(i).Clear();
	ind_real(i).Clear();
      }
  }
  
  
  //! Imaginary part of row i is removed.
  /*!
    \param[in] i row number.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ClearImagRow(int i)
  {
    int l = this->GetImagRowSize(i);
    if (l > 0)
      {
	nz_imag_ -= l;
	val_imag(i).Clear();
	ind_imag(i).Clear();
      }
  }
  
  
  //! Two rows/columns are swapped (real part).
  /*!
    \param[in] i first row number.
    \param[in] i_ second row number.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  SwapRealRow(int i,int i_)
  {
    if (i != i_)
      {
	int j_ = ind_real(i_).GetM();
	int j = ind_real(i).GetM();
	int* tmp = ind_real(i).GetData();
	ind_real(i).Nullify();
	pointer tmp_val = val_real(i).GetData();
	val_real(i).Nullify();
	ind_real(i).SetData(j_, ind_real(i_).GetData());
	val_real(i).SetData(j_, val_real(i_).GetData());
	ind_real(i_).Nullify();
	val_real(i_).Nullify();
	ind_real(i_).SetData(j, tmp);
	val_real(i_).SetData(j, tmp_val);
      }
  }
  
  
  //! Two rows/columns are swapped (imaginary part).
  /*!
    \param[in] i first row number.
    \param[in] i_ second row number.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  SwapImagRow(int i,int i_)
  {
    if (i != i_)
      {
	int j_ = ind_imag(i_).GetM();
	int j = ind_imag(i).GetM();
	int* tmp = ind_imag(i).GetData();
	ind_imag(i).Nullify();
	pointer tmp_val = val_imag(i).GetData();
	val_imag(i).Nullify();
	ind_imag(i).SetData(j_, ind_imag(i_).GetData());
	val_imag(i).SetData(j_, val_imag(i_).GetData());
	ind_imag(i_).Nullify();
	val_imag(i_).Nullify();
	ind_imag(i_).SetData(j, tmp);
	val_imag(i_).SetData(j, tmp_val);
      }
  }
  
  
  //! Column numbers of row i are modified (real part).
  /*!
    Useful method when a permutation is applied to a matrix.
    \param[in] i row number.
    \param[in] new_index column numbers.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ReplaceIndexRealRow(int i, IVect& new_index)
  {
    ind_real(i).Clear();
    ind_real(i).SetData(new_index.GetM(), new_index.GetData());
    new_index.Nullify();
  }
  
  
  //! Column numbers of row i are modified (imaginary part).
  /*!
    Useful method when a permutation is applied to a matrix.
    \param[in] i row number.
    \param[in] new_index column numbers.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ReplaceIndexImagRow(int i, IVect& new_index)
  {
    ind_imag(i).Clear();
    ind_imag(i).SetData(new_index.GetM(), new_index.GetData());
    new_index.Nullify();
  }
  
  
  /*******************
   * BASIC FUNCTIONS *
   *******************/
  
  
  //! Returns the number of rows.
  /*!
    \return the number of rows.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T,Prop,Storage,Allocator>::GetM() const
  {
    return m_;
  }


  //! Returns the number of columns.
  /*!
    \return the number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T,Prop,Storage,Allocator>::GetN() const
  {
    return n_;
  }
  
  
  //! Returns the number of non-zero elements (real part).
  /*!
    \return The number of non-zero elements for real part of matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealNonZeros() const
  {
    return nz_real_;
  }
  
  
  //! Returns the number of non-zero elements (imaginary part).
  /*!
    \return The number of non-zero elements for imaginary part of matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagNonZeros() const
  {
    return nz_imag_;
  }
  
  
  //! Returns the number of elements stored in memory (real part).
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealDataSize() const
  {
    return nz_real_;
  }
  
  
  //! Returns the number of elements stored in memory (imaginary part).
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagDataSize() const
  {
    return nz_imag_;
  }
  
  
  //! Returns the number of elements stored in memory (real+imaginary part).
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetDataSize() const
  {
    return nz_real_ + nz_imag_;
  }
  
  
  //! Returns the size of row i (real part).
  /*!
    \param[in] i row number.
    \return The number of non-zero entries in row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealRowSize(int i) const
  {
    return ind_real(i).GetM();
  }
  
  
  //! Returns the size of row i (imaginary part).
  /*!
    \param[in] i row number.
    \return The number of non-zero entries in row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagRowSize(int i) const
  {
    return ind_imag(i).GetM();
  }
  
  
  //! Returns column indices of non-zero entries in row (real part).
  /*!
    \param[in] i row number.
    \return The array of column indices of non-zero entries
    of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealInd(int i) const
  {
    return ind_real(i).GetData();
  }
  
  
  //! Returns values of non-zero entries of a row (real part).
  /*!
    \param[in] i row number.
    \return The array of values of non-zero entries of row i.
  */
  template <class T, class Prop, class Storage, class Allocator> inline
  typename Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::pointer
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::GetRealData(int i)
    const
  {
    return val_real(i).GetData();
  }
  
  
  //! Returns column indices of non-zero entries in row (imaginary part).
  /*!
    \param[in] i row number.
    \return the array of column indices of non-zero entries
    of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagInd(int i) const
  {
    return ind_imag(i).GetData();
  }
  
  
  //! Returns values of non-zero entries of a row (imaginary part).
  /*!
    \param[in] i row number.
    \return The array of values of non-zero entries of row i.
  */
  template <class T, class Prop, class Storage, class Allocator> inline
  typename Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::pointer
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::GetImagData(int i)
    const
  {
    return val_imag(i).GetData();
  }
  
  
  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline complex<T>
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::operator()
    (int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    int k;
    int a, b;
    complex<T> res(0);
    
    a = Storage::GetFirst(i, j);
    b = Storage::GetSecond(i, j);
    
    for (k = 0; k < ind_real(a).GetM(); k++)
      if (ind_real(a)(k) == b)
	res += complex<T>(val_real(a)(k), 0.0);
    
    for (k = 0; k < ind_imag(a).GetM(); k++)
      if (ind_imag(a)(k) == b)
	res += complex<T>(0.0, val_imag(a)(k));
    
    return res;
  }
  
  
  //! Returns j-th non-zero value of row i (real part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::const_reference
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ValueReal(int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::value", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->GetRealRowSize(i))
      throw WrongCol("Matrix_ArraySparse::value", "Index should be in [0, " +
		     to_str(this->GetRealRowSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return val_real(i)(j);
  }
  
  
  //! Returns j-th non-zero value of row i (real part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::reference
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ValueReal(int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::value", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->GetRealRowSize(i))
      throw WrongCol("Matrix_ArraySparse::value", "Index should be in [0, " +
		     to_str(this->GetRealRowSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return val_real(i)(j);
  }
  
  
  //! Returns column number of j-th non-zero value of row i (real part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return Column number of j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  IndexReal(int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::index", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->GetRealRowSize(i))
      throw WrongCol("Matrix_ArraySparse::index", "Index should be in [0, " +
		     to_str(this->GetRealRowSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return ind_real(i)(j);
  }
  
  
  //! Returns column number of j-th non-zero value of row i (real part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return column number of j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int& Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  IndexReal(int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::index", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->GetRealRowSize(i))
      throw WrongCol("Matrix_ArraySparse::index", "Index should be in [0, "
		     + to_str(this->GetRealRowSize(i)-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif
    
    return ind_real(i)(j);
  }
  
  
  //! Returns j-th non-zero value of row i (imaginary part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  const_reference Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ValueImag(int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::value", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->GetImagRowSize(i))
      throw WrongCol("Matrix_ArraySparse::value", "Index should be in [0, " +
		     to_str(this->GetImagRowSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return val_imag(i)(j);
  }
  
  
  //! Returns j-th non-zero value of row i (imaginary part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  reference Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  ValueImag (int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::value", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->GetImagRowSize(i))
      throw WrongCol("Matrix_ArraySparse::value", "Index should be in [0, " +
		     to_str(this->GetImagRowSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return val_imag(i)(j);
  }
  
  
  //! Returns column number of j-th non-zero value of row i (imaginary part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return Column number of j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  IndexImag(int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::index", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->GetImagRowSize(i))
      throw WrongCol("Matrix_ArraySparse::index", "Index should be in [0, " +
		     to_str(this->GetImagRowSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return ind_imag(i)(j);
  }
  
  
  //! Returns column number of j-th non-zero value of row i (imaginary part).
  /*!
    \param[in] i row number.
    \param[in] j local number.
    \return column number of j-th non-zero entry of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int& Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  IndexImag(int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArrayComplexSparse::index",
		     "Index should be in [0, " + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    
    if (j < 0 || j >= this->GetImagRowSize(i))
      throw WrongCol("Matrix_ArraySparse::index", "Index should be in [0, " +
		     to_str(this->GetImagRowSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return ind_imag(i)(j);
  }
  
  
  //! Returns column indices of non-zero entries in row (real part).
  /*!
    \param[in] i row number.
    \return The array of column indices of non-zero entries of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const IVect& Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealIndexRow(int i) const
  {
    return ind_real(i);
  }
  
  
  //! Returns column indices of non-zero entries in row (imaginary part).
  /*!
    \param[in] i row number.
    \return The array of column indices of non-zero entries of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const IVect& Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagIndexRow(int i) const
  {
    return ind_imag(i);
  }
  
  
  //! Returns values of non-zero entries of a row (real part).
  /*!
    \param[in] i row number.
    \return The array of values of non-zero entries of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const Vector<T,Vect_Full,Allocator>&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealValueRow(int i) const
  {
    return val_real(i);
  }
  
  
  //! Returns values of non-zero entries of a row (imaginary part).
  /*!
    \param[in] i row number.
    \return The array of values of non-zero entries of row i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const Vector<T,Vect_Full,Allocator>&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagValueRow(int i) const
  {
    return val_imag(i);
  }
  
  
  //! coefficient a is added in the real part of matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] a value to add.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  AddRealInteraction(int i, int j, const T& a)
  {
    int n = ind_real(i).GetM();
    for (int k = 0; k < n; k++)
      if (ind_real(i)(k) == j)
	{
	  val_real(i)(k) += a;
	  return;
	}
    
    // The interaction doesn't exist, we add it on the correct position.
    vect_value new_val(n+1); IVect new_ind(n+1);
    int k = 0;
    while (k < n && ind_real(i)(k) < j)
      {
	new_ind(k) = ind_real(i)(k);
	new_val(k) = val_real(i)(k);
	k++;
      }
    
    // The new entry.
    new_ind(k) = j;
    new_val(k) = a;
    k++;
    
    while (k <= n)
      {
	new_ind(k) = ind_real(i)(k-1);
	new_val(k) = val_real(i)(k-1);
	k++;
      }
    
    ind_real(i).Clear();
    val_real(i).Clear();
    n++;
    ind_real(i).SetData(n, new_ind.GetData());
    val_real(i).SetData(n, new_val.GetData());
    new_ind.Nullify();
    new_val.Nullify();
    nz_real_++;
  }
  
  
  //! Coefficients are added in the row of a matrix (real part).
  /*!
    The method sorts given coefficients and adds them
    in the correct positions.
    \param[in] i row number.
    \param[in] nb_interac number of coefficients to add.
    \param[in] col_interac column numbers.
    \param[in] val_interac coefficients.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template<class Storage1,class Allocator1>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  AddRealInteraction(int i, int nb_interac, IVect col_interac,
		     Vector<T,Storage1,Allocator1> val_interac)
  {
    int n = ind_real(i).GetM();
    // Column numbers are sorted and we add values for same numbers.
    Seldon::Assemble(nb_interac, col_interac, val_interac);
    
    // Gets already existing interactions.
    int nb_new = 0;
    Vector<bool> new_interac(nb_interac);
    new_interac.Fill(true);
    int k = 0;
    for (int j = 0; j < nb_interac; j++)
      {
	while (k < n && ind_real(i)(k) < col_interac(j))
	  k++;
	
	if (k < n && col_interac(j) == ind_real(i)(k))
	  {
	    new_interac(j) = false;
	    val_real(i)(k) += val_interac(j);
	  }
	else
	  nb_new++;
      }
    
    if (nb_new > 0)
      {
	// Some interactions doesn't exist, we add them on the correct
	// position.
	vect_value new_val(n+nb_new); IVect new_ind(n+nb_new);
	int nb = 0, k = 0;
	for (int j = 0; j < nb_interac; j++)
	  if (new_interac(j))
	    {
	      while (k < n && ind_real(i)(k) < col_interac(j))
		{
		  new_ind(nb) = ind_real(i)(k);
		  new_val(nb) = val_real(i)(k);
		  k++;
		  nb++;
		}
	      
	      // The new entry.
	      new_ind(nb) = col_interac(j);
	      new_val(nb) = val_interac(j);
	      nb++;
	    }
	
	while (k < n)
	  {
	    new_ind(nb) = ind_real(i)(k);
	    new_val(nb) = val_real(i)(k);
	    k++;
	    nb++;
	  }
	
	n += nb_new;
	ind_real(i).Clear();
	val_real(i).Clear();
	ind_real(i).SetData(n, new_ind.GetData());
	val_real(i).SetData(n, new_val.GetData());
	new_ind.Nullify();
	new_val.Nullify();
	nz_real_ += nb_new;
      }
  }
  
  
  //! Coefficient a is added in the imaginary part of matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] a value to add.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  AddImagInteraction(int i, int j, const T& a)
  {
    int n = ind_imag(i).GetM();
    for (int k = 0; k < n; k++)
      if (ind_imag(i)(k) == j)
	{
	  val_imag(i)(k) += a;
	  return;
	}
    
    // The interaction doesn't exist, we add it on the correct position.
    vect_value new_val(n+1); IVect new_ind(n+1);
    int k = 0;
    while (k < n && ind_imag(i)(k) < j)
      {
	new_ind(k) = ind_imag(i)(k);
	new_val(k) = val_imag(i)(k);
	k++;
      }
    
    // The new entry.
    new_ind(k) = j;
    new_val(k) = a;
    k++;
    
    while (k <= n)
      {
	new_ind(k) = ind_imag(i)(k-1);
	new_val(k) = val_imag(i)(k-1);
	k++;
      }
    
    ind_imag(i).Clear();
    val_imag(i).Clear();
    n++;
    ind_imag(i).SetData(n, new_ind.GetData());
    val_imag(i).SetData(n, new_val.GetData());
    new_ind.Nullify();
    new_val.Nullify();
    nz_imag_++;
  }
  
  
  //! Coefficients are added in the row of a matrix (imaginary part).
  /*!
    The method sorts given coefficients and adds them
    in the correct positions.
    \param[in] i row number.
    \param[in] nb_interac number of coefficients to add.
    \param[in] col_interac column numbers.
    \param[in] val_interac coefficients.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template<class Storage1,class Allocator1>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  AddImagInteraction(int i, int nb_interac, IVect col_interac,
		     Vector<T,Storage1,Allocator1> val_interac)
  {
    int n = ind_imag(i).GetM();
    // Column numbers are sorted and we add values for same numbers.
    Seldon::Assemble(nb_interac, col_interac, val_interac);
    
    // Gets already existing interactions.
    int nb_new = 0;
    Vector<bool> new_interac(nb_interac);
    new_interac.Fill(true);
    int k = 0;
    for (int j = 0; j < nb_interac; j++)
      {
	while (k < n && ind_imag(i)(k) < col_interac(j))
	  k++;
	
	if (k < n && col_interac(j) == ind_imag(i)(k))
	  {
	    new_interac(j) = false;
	    val_imag(i)(k) += val_interac(j);
	  }
	else
	  nb_new++;
      }
    
    if (nb_new > 0)
      {
	// Some interactions doesn't exist,we add them on the correct
	// position.
	vect_value new_val(n+nb_new);
	IVect new_ind(n+nb_new);
	int nb = 0, k = 0;
	for (int j = 0; j < nb_interac; j++)
	  if (new_interac(j))
	    {
	      while (k < n && ind_imag(i)(k) < col_interac(j))
		{
		  new_ind(nb) = ind_imag(i)(k);
		  new_val(nb) = val_imag(i)(k);
		  k++; nb++;
		}
	      
	      // The new entry.
	      new_ind(nb) = col_interac(j);
	      new_val(nb) = val_interac(j);
	      nb++;
	    }
	
	while (k < n)
	  {
	    new_ind(nb) = ind_imag(i)(k);
	    new_val(nb) = val_imag(i)(k);
	    k++;
	    nb++;
	  }
	
	n += nb_new;
	ind_imag(i).Clear();
	val_imag(i).Clear();
	ind_imag(i).SetData(n, new_ind.GetData());
	val_imag(i).SetData(n, new_val.GetData());
	new_ind.Nullify();
	new_val.Nullify();
	nz_imag_ += nb_new;
      }
  }
  
  
  /************************
   * CONVENIENT FUNCTIONS *
   ************************/
  
  
  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }
  
  
  //! Displays real part of row i.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  PrintRealRow(int i) const
  {
    cout << " Row " << i << endl;
    cout << " Index " << endl << ind_real(i) << endl;
    cout << " Value " << endl << val_real(i) << endl;
  }
  
  
  //! Displays imaginary part of row i.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  PrintImagRow(int i) const
  {
    cout << " Row " << i << endl;
    cout << " Index " << endl << ind_imag(i) << endl;
    cout << " Value " << endl << val_imag(i) << endl;
  }
  
  
  //! Assembles the matrix.
  /*!
    All the row numbers are sorted.
    If same row numbers exist, values are added.
    \warning If you are using the methods AddInteraction/AddInteractions,
    you don't need to call that method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Assemble()
  {
    for (int i = 0; i < m_ ; i++)
      {
	AssembleRealRow(i);
	AssembleImagRow(i);
      }
  }
  
  
  //! Assembles a row (real part).
  /*!
    \param[in] i row number.
    \warning If you are using the methods AddInteraction/AddInteractions,
    you don't need to call that method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  AssembleRealRow(int i)
  {
    int new_size = ind_real(i).GetM();
    Seldon::Assemble(new_size, ind_real(i), val_real(i));
    ResizeRealRow(i, new_size);
  }
  
  
  //! Assembles a row (imaginary part).
  /*!
    \param[in] i row number.
    \warning If you are using the methods AddInteraction/AddInteractions,
    you don't need to call that method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  AssembleImagRow(int i)
  {
    int new_size = ind_imag(i).GetM();
    Seldon::Assemble(new_size, ind_imag(i), val_imag(i));
    ResizeImagRow(i, new_size);
  }
  
  
  //! Matrix is initialized to the identity matrix.
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  SetIdentity()
  {
    this->n_ = this->m_;
    this->nz_real_ = this->m_;
    this->nz_imag_ = 0;
    for (int i = 0; i < this->m_; i++)
      {
	ind_real(i).Reallocate(1);
	ind_real(0) = i;
	val_real(i).Reallocate(1);
	val_real(0) = T(1);
      }
  }
  
  
  //! Non-zero entries are set to 0 (but not removed).
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Zero()
  {
    for (int i = 0; i < val_real.GetM(); i++)
      {
	val_real(i).Zero();
	val_imag(i).Zero();
      }
  }
  
  
  //! Non-zero entries are filled with values 0, 1, 2, 3 ...
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Fill()
  {
    int value = 0;
    for (int i = 0; i < val_real.GetM(); i++)
      {
	for (int j = 0; j < val_real(i).GetM(); j++)
	  val_real(i)(j) = value++;
	
	for (int j = 0; j < val_imag(i).GetM(); j++)
	  val_imag(i)(j) = value++;
      }
  }
  
  
  //! Non-zero entries are set to a given value x.
  template <class T, class Prop, class Storage, class Allo> template<class T0>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allo>::
  Fill(const complex<T0>& x)
  {
    for (int i = 0; i < val_real.GetM(); i++)
      {
	val_real(i).Fill(real(x));
	val_imag(i).Fill(imag(x));
      }
  }
  
  
  //! Non-zero entries are set to a given value x.
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  inline Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::operator=
  (const complex<T0>& x)
  {
    this->Fill(x);
  }
  
  
  //! Non-zero entries take a random value.
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  FillRand()
  {
    for (int i = 0; i < val_real.GetM(); i++)
      {
	val_real(i).FillRand();
	val_imag(i).FillRand();
      }
  }
  
  
  ///////////////////////
  // MATRIX<ROWSPARSE> //
  ///////////////////////
  
  
  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::Matrix()  throw():
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>()
  {
  }
  
  
  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>(i, j)
  {
  }
  
  
  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  AddInteraction(int i, int j, const complex<T>& val)
  {
    if (real(val) != T(0))
      Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>::
	AddRealInteraction(i, j, real(val));
	
    if (imag(val) != T(0))
      Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>::
	AddImagInteraction(i, j, imag(val));
  }
  
  
  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<complex<T>, Vect_Full, Alloc1>& val)
  {
    IVect new_col_real(nb), new_col_imag(nb);
    Vector<T> new_val_real(nb), new_val_imag(nb);
    int nb_real = 0, nb_imag = 0;
    for (int j = 0; j < nb; j++)
      {
	if (real(val(j)) != T(0))
	  {
	    new_col_real(nb_real) = col(j);
	    new_val_real(nb_real) = real(val(j));
	    nb_real++;
	  }
	
	if (imag(val(j)) != T(0))
	  {
	    new_col_imag(nb_imag) = col(j);
	    new_val_imag(nb_imag) = imag(val(j));
	    nb_imag++;
	  }
      }
    
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>::
      AddRealInteraction(i, nb_real, new_col_real, new_val_real);
    
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse, Allocator>::
      AddImagInteraction(i, nb_imag, new_col_imag, new_val_imag);
  }
  
  
  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<complex<T>, Vect_Full, Alloc1>& val)
  {
    for (int j = 0; j < nb ; j++)
      {
        if (real(val(j)) != T(0))
	  Matrix_ArrayComplexSparse<T, Prop,
	    ArrayRowComplexSparse, Allocator>::
	    AddRealInteraction(row(j), i, real(val(j)));
	
	if (imag(val) != T(0))
	  Matrix_ArrayComplexSparse<T, Prop,
	    ArrayRowComplexSparse, Allocator>::
	    AddImagInteraction(row(j), j, imag(val(j)));
      }
  }
  
  
  //////////////////////////
  // MATRIX<ROWSYMSPARSE> //
  //////////////////////////
  
  
  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::Matrix()  throw():
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowSymComplexSparse, Allocator>()
  {
  }
  
  
  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArrayComplexSparse<T, Prop, ArrayRowSymComplexSparse,
			      Allocator>(i, j)
  {
  }
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop,  class Allocator>
  inline complex<T> Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  operator() (int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    int k;
    int a,b;
    complex<T> res(0);
    
    if (i <= j)
      {
	a = i; b = j;
      }
    else
      {
	a = j; b = i;
      }
    
    for (k = 0; k < this->ind_real(a).GetM(); k++)
      {
	if (this->ind_real(a)(k) == b)
	  res += complex<T>(this->val_real(a)(k), 0.0);
      }
    
    for (k = 0; k < this->ind_imag(a).GetM() ; k++)
      {
	if (this->ind_imag(a)(k) == b)
	  res += complex<T>(0.0, this->val_imag(a)(k));
      }
    
    return res;
  }
  
  
  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  AddInteraction(int i, int j, const complex<T>& val)
  {
    if (i <= j)
      {
	if (real(val) != T(0))
	  Matrix_ArrayComplexSparse<T, Prop,
	    ArrayRowSymComplexSparse, Allocator>::
	    AddRealInteraction(i, j, real(val));
	
	if (imag(val) != T(0))
	  Matrix_ArrayComplexSparse<T, Prop,
	    ArrayRowSymComplexSparse, Allocator>::
	    AddImagInteraction(i, j, imag(val));
      }
  }
  
  
  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<complex<T>, Vect_Full, Alloc1>& val)
  {
    IVect new_col_real(nb), new_col_imag(nb);
    Vector<T> new_val_real(nb), new_val_imag(nb);
    int nb_real = 0, nb_imag = 0;
    for (int j = 0; j < nb; j++)
      if (i <= col(j))
	{
	  if (real(val(j)) != T(0))
	    {
	      new_col_real(nb_real) = col(j);
	      new_val_real(nb_real) = real(val(j));
	      nb_real++;
	    }
	  
	  if (imag(val(j)) != T(0))
	    {
	      new_col_imag(nb_imag) = col(j);
	      new_val_imag(nb_imag) = imag(val(j));
	      nb_imag++;
	    }
	}
    
    Matrix_ArrayComplexSparse<T,Prop,ArrayRowSymComplexSparse,Allocator>::
      AddRealInteraction(i, nb_real, new_col_real, new_val_real);
    
    Matrix_ArrayComplexSparse<T,Prop,ArrayRowSymComplexSparse,Allocator>::
      AddImagInteraction(i, nb_imag, new_col_imag, new_val_imag);
  }
  
  
  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<complex<T>, Vect_Full, Alloc1>& val)
  {
    // Symmetric matrix, row = column.
    AddInteractionRow(i, nb, row, val);
  }

  
  // operator<< overloaded for matrices.
  template <class T, class Prop, class Storage, class Allocator>
  ostream&
  operator <<(ostream& out,
	      const Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>& A)
  {
    out << "Real value" << endl;
    out << A.val_real << endl;
    out << "Real index" << endl;
    out << A.ind_real << endl;
    out << "Imaginary value" << endl;
    out << A.val_imag << endl;
    out << "Imaginary index" << endl;
    out << A.ind_imag << endl;
    
    return out;
  }

  
} // namespace Seldon

#define FILE_MATRIX_ARRAY_COMPLEX_SPARSE_CXX
#endif

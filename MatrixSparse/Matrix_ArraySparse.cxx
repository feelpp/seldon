// Copyright (C) 2001-2003 Vivien Mallet
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

#ifndef FILE_MATRIX_ARRAY_SPARSE_CXX

#include "Matrix_ArraySparse.hxx"

namespace Seldon
{

  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocator>::Matrix_ArraySparse()
    : ind(),val()
  {
    nz_ = 0;
    this->m_ = 0; this->n_ = 0;
  }
  
  
  //! Constructor.
  /*!
    Builds an empty i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Matrix_ArraySparse(int i, int j) :
    ind(Storage::GetFirst(i, j)),val(Storage::GetFirst(i, j))
  {
    nz_=0;
    this->m_ = Storage::GetFirst(i,j); this->n_ = Storage::GetSecond(i,j);
  }
  
  
  //! Destructor.
  template <class T, class Prop, class Storage, class Allocat>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocat>::~Matrix_ArraySparse()
  {
    this->m_ = 0;
    this->n_ = 0;
    this->nz_ = 0;
  }
  
  
  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix is empty
    (0 by 0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_ArraySparse();
  }
  
  
  /*********************
   * MEMORY MANAGEMENT *
   *********************/
  
  
  //! returns size occupied by the matrix in bytes
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  GetMemorySize() const
  {
    return (sizeof(int)*(7+2*ind.GetM()+2*val.GetM()+GetDataSize())
	    +sizeof(T)*GetDataSize());
  }
  
  
  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Data is lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Reallocate(int i, int j)
  {
    // clears previous entries
    Clear();
    
    this->m_ = i;
    this->n_ = j;
    
    int n = Storage::GetFirst(i,j);
    val.Reallocate(n); ind.Reallocate(n);
    nz_ = 0;
  }
  
  
  //! Reallocates additional memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    Data is kept
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Resize(int i, int j)
  {
    if ((i != this->m_) || (j != this->n_))
      {
	int n = Storage::GetFirst(m_,n_);
	int new_n = Storage::GetFirst(i,j);
	Vector<IVect,Vect_Full,NewAlloc<IVect> > new_ind;
	Vector<vect_value,Vect_Full,NewAlloc<vect_value> > new_val;
	new_ind.Reallocate(new_n); new_val.Reallocate(new_n);
	
	nz_ = 0;
	for (int k=0 ; k<min(n,new_n) ; k++)
	  {
	    new_ind(k).SetData(ind(k).GetM(),ind(k).GetData());
	    ind(k).Nullify();
	    new_val(k).SetData(val(k).GetM(),val(k).GetData());
	    val(k).Nullify();
	    nz_ += new_ind(k).GetM();
	  }
	
	ind.Clear(); ind.SetData(new_ind.GetM(),new_ind.GetData());
	new_ind.Nullify();
	val.Clear(); val.SetData(new_val.GetM(),new_val.GetData());
	new_val.Nullify();
	
	this->m_=i;
	this->n_=j;
	
      }
  }
  
  
  //! changes the size of i-th row/column of the matrix
  /*!
    \param[in] i row/column number
    \param[in] j new size of row/column
    \warning data may be lost
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  ReallocateVector(int i, int j)
  {
    if (j != this->GetVectorSize(i))
      {
	nz_ += j - GetVectorSize(i);
	val(i).Reallocate(j); ind(i).Reallocate(j);
      }
  }
  
  
  //! changes the size of i-th row/column of the matrix
  /*!
    Data is kept
    \param[in] i row/column number
    \param[in] j new size of row/column
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  ResizeVector(int i, int j)
  {
    int m = this->GetVectorSize(i);
    if (j != m)
      {
	nz_ += j - m;
	pointer data_;
	data_ = allocator_.allocate(j, this);
	IVect Ind_New(j);
	
	for (int k = 0; k < min(m,j); k++)
	  {
	    data_[k] = val(i)(k);
	    Ind_New(k) = ind(i)(k);
	  }
	
	ind(i).SetData(j,Ind_New.GetData());
	Ind_New.Nullify();
	val(i).SetData(j,data_);
      }
  }
  
  
  //! row-column i is removed
  /*!
    \param[in] i row/column number
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  ClearVector(int i)
  {
    int l = this->GetVectorSize(i);
    if (l > 0)
      {
	nz_ -= l;
	val(i).Clear(); ind(i).Clear();
      }
  }
  
  
  //! two rows/columns are swapped
  /*!
    \param[in] i first row/column number
    \param[in] i_ second row/column number
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  SwapVector(int i, int i_)
  {
    if (i != i_)
      {
	int j_ = ind(i_).GetM(); int j = ind(i).GetM();
	int* tmp = ind(i).GetData(); ind(i).Nullify();
	pointer tmp_val = val(i).GetData(); val(i).Nullify();
	ind(i).SetData(j_,ind(i_).GetData());
	val(i).SetData(j_,val(i_).GetData());
	ind(i_).Nullify(); val(i_).Nullify();
	ind(i_).SetData(j,tmp);
	val(i_).SetData(j,tmp_val);
      }
  }
  
  
  //! column numbers of row i are modified
  /*!
    useful method when a permutation is applied to a matrix
    \param[in] i row number (or column number)
    \param[in] new_index column numbers (or row numbers)
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  ReplaceIndexVector(int i, IVect& new_index)
  {
    ind(i).Clear();
    ind(i).SetData(new_index.GetM(),new_index.GetData());
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
  inline int Matrix_ArraySparse<T,Prop,Storage,Allocator>::GetM() const
  {
    return m_;
  }
  
  
  //! Returns the number of columns.
  /*!
    \return the number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T,Prop,Storage,Allocator>::GetN() const
  {
    return n_;
  }
  
  
  //! Returns the number of non-zero entries.
  /*!
    \return The number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetNonZeros()
    const
  {
    return nz_;
  }
  
  
  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetDataSize()
    const
  {
    return nz_;
  }
  
  
  //! Returns the size of row/column i
  /*!
    \param[in] i row/column number
    \return the number of non-zero entries in row/column i
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  GetVectorSize(int i) const
  {
    return ind(i).GetM();
  }
  
  
  //! Returns (row or column) indices of non-zero entries in row
  /*!
    \param[in] i row (or column) number
    \return the array of column (or row) indices of non-zero entries
    of row (or column) i
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetInd(int i)
    const
  {
    return ind(i).GetData();
  }
  
  
  //! Returns values of non-zero entries of a row/column
  /*!
    \param[in] i row (or column) number
    \return The array of values of non-zero entries of row/column i
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_ArraySparse<T, Prop, Storage, Allocator>::pointer
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetData(int i) const
  {
    return val(i).GetData();
  }
  
  
  //! Returns row (or column) indices of non-zero entries
  /*!
    \return Arrays of row (or column) indices of non-zero entries
    There is a different array for each row/column
  */
  template <class T, class Prop, class Storage, class Allocat>
  inline IVect* Matrix_ArraySparse<T, Prop, Storage, Allocat>::GetInd() const
  {
    return ind.GetData();
  }
  
  
  //! Returns values of non-zero entries
  /*!
    \return Arrays of values of non-zero entries
    There is a different array for each row/column
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline
  typename Matrix_ArraySparse<T, Prop, Storage, Allocator>::vect_value_ptr
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetData() const
  {
    return val.GetData();
  }
  
  
  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_ArraySparse<T, Prop, Storage, Allocator>::value_type
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::operator() (int i, int j)
    const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->m_-1) +
		     "],but is equal to " + to_str(i) + ".");
    
    if ((j < 0)||(j >= this->n_))
      throw WrongCol("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->n_-1) +
		     "], but is equal to " + to_str(j) + ".");
#endif
    
    int k;
    int a, b;
    
    a = Storage::GetFirst(i,j);
    b = Storage::GetSecond(i,j);
    for (k = 0; k < ind(a).GetM(); k++)
      {
	if (ind(a)(k) == b)
	  return val(a)(k);
      }
    
    return T(0);
  }
  
  
  //! Access operator
  /*!
    Returns the value of element (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_ArraySparse<T, Prop, Storage, Allocator>::reference
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::operator() (int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->m_-1) +
		     "],but is equal to " + to_str(i) + ".");
    
    if ((j < 0)||(j >= this->n_))
      throw WrongCol("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->n_-1) +
		     "], but is equal to " + to_str(j) + ".");
#endif
    
    int a, b;
    a = Storage::GetFirst(i,j);
    b = Storage::GetSecond(i,j);
    int size_vec = ind(a).GetM();
    int k = 0;
    // we search the entry
    while ((k < size_vec)&&(ind(a)(k) < b))
      k++;
    
    if ((k >= size_vec)||(ind(a)(k) != b))
      {
	// the entry does not exist, we add a null value
	AddInteraction(a, b, T(0));
      }
    
    return val(a)(k);
  }
  
  
  //! Returns j-th non-zero value of row/column i
  /*!
    \param[in] i row/column number
    \param[in] j local number
    \return j-th non-zero entry of row/column i
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  const_reference Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Value (int i, int j) const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix_ArraySparse::value", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if ((j < 0)||(j >= this->GetVectorSize(i)))
      throw WrongCol("Matrix_ArraySparse::value", "Index should be in [0, " +
		     to_str(this->GetVectorSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return val(i)(j);
  }
  
  
  //! Returns j-th non-zero value of row/column i
  /*!
    \param[in] i row/column number
    \param[in] j local number
    \return j-th non-zero entry of row/column i
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_ArraySparse<T, Prop, Storage, Allocator>::reference
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Value (int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix_ArraySparse::value", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if ((j < 0)||(j >= this->GetVectorSize(i)))
      throw WrongCol("Matrix_ArraySparse::value", "Index should be in [0, " +
		     to_str(this->GetVectorSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return val(i)(j);
  }
  
  
  //! Returns column/row number of j-th non-zero value of row/column i
  /*!
    \param[in] i row/column number
    \param[in] j local number
    \return column/row number of j-th non-zero value of row/column i
  */
  template <class T, class Prop, class Storage, class Allocator> inline
  int Matrix_ArraySparse<T, Prop, Storage, Allocator>::Index(int i, int j)
    const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix_ArraySparse::index", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if ((j < 0)||(j >= this->GetVectorSize(i)))
      throw WrongCol("Matrix_ArraySparse::index", "Index should be in [0, " +
		     to_str(this->GetVectorSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return ind(i)(j);
  }
  
  
  //! Returns column/row number of j-th non-zero value of row/column i
  /*!
    \param[in] i row/column number
    \param[in] j local number
    \return column/row number of j-th non-zero value of row/column i
  */
  template <class T, class Prop, class Storage, class Allocator> inline
  int& Matrix_ArraySparse<T, Prop, Storage, Allocator>::Index(int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix_ArraySparse::index", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    
    if ((j < 0)||(j >= this->GetVectorSize(i)))
      throw WrongCol("Matrix_ArraySparse::index", "Index should be in [0, " +
		     to_str(this->GetVectorSize(i)-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    return ind(i)(j);
  }
  
  
  //! Redefines the matrix
  /*!
    \param[in] m number of rows
    \param[in] n number of columns
    \param[in] nnz number of non-zero entries
    \param[in] ind_ptr arrays of row/column numbers
    \param[in] val_ptr arrays of values
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  SetData(int m, int n, int nnz, IVect* ind_ptr,
	  typename Matrix_ArraySparse<T, Prop, Storage,
	  Allocator>::vect_value_ptr val_ptr)
  {
    m_ = m; n_ = n; nz_ = nnz;
    ind.Clear(); val.Clear();
    ind.SetData(n, ind_ptr); val.SetData(n, val_ptr);
  }
  
  
  //!  Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Nullify()
  {
    m_ = 0;
    n_ = 0;
    nz_ = 0;
    ind.Nullify();
    val.Nullify();
  }
  
  
  //! coefficient a is added in the matrix
  /*!
    \param[in] i row/column number
    \param[in] j column/row number
    \param[in] a value to add
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  AddInteraction(int i, int j, const T& a)
  {
    int n = ind(i).GetM();
    // we look for the position where the entry is
    int pos = 0;
    while ((pos < n)&&(ind(i)(pos) < j))
      pos++;
    
    // if the entry already exists, we add a
    if (pos < n)
      if (ind(i)(pos) == j)
	{
	  val(i)(pos) += a;
	  return;
	}
    
    // the interaction doesn't exist, the row/column is reallocated
    vect_value new_val(n+1); IVect new_ind(n+1);
    for (int k = 0; k < pos; k++)
      {
	new_ind(k) = ind(i)(k);
	new_val(k) = val(i)(k);
      }
    
    // the new entry
    new_ind(pos) = j; new_val(pos) = a;
    
    // end of the row/column
    for (int k = pos+1; k <= n; k++)
      {
	new_ind(k) = ind(i)(k-1);
	new_val(k) = val(i)(k-1);
      }
    
    ind(i).Clear(); val(i).Clear(); n++;
    ind(i).SetData(n, new_ind.GetData());
    val(i).SetData(n, new_val.GetData());
    new_ind.Nullify(); new_val.Nullify();
    
    // new entry -> we increase the number of non-zero entries
    nz_++;
  }
  
  
  //! coefficients are added in the row/column of a matrix
  /*!
    The method sorts given coefficients and adds them
    in the correct positions
    \param[in] i row/column number
    \param[in] nb_interac number of coefficients to add
    \param[in] col_interac column/row numbers
    \param[in] val_interac coefficients
  */
  template <class T, class Prop, class Storage, class Allocator>
  template<class Storage1, class Allocator1>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  AddInteractionVector(int i, int nb_interac, IVect col_interac,
		       Vector<T, Storage1, Allocator1> val_interac)
  {
    // There is no reference in arguments, because we want to
    // perform a local copy and keep the input arguments unchanged
        
    // number of elements in the row/column i
    int n = ind(i).GetM();
    
    Seldon::Assemble(nb_interac, col_interac, val_interac);
    
    // we add new entries which have already a slot in the matrix
    int nb_new = 0;
    Vector<bool> new_interac(nb_interac); new_interac.Fill(true);
    int k = 0;
    for (int j = 0; j < nb_interac; j++)
      {
	while ((k < n)&&(ind(i)(k) < col_interac(j)))
	  k++;
	
	if ((k < n)&&(col_interac(j) == ind(i)(k)))
	  {
	    new_interac(j) = false;
	    val(i)(k) += val_interac(j);
	  }
	else
	  nb_new++;
      }
    
    if (nb_new > 0)
      {
	// some new entries have no slot, we add them on the correct position
	vect_value new_val(n+nb_new); IVect new_ind(n+nb_new);
	int nb = 0, k = 0;
	for (int j = 0; j < nb_interac; j++)
	  if (new_interac(j))
	    {
	      while ((k < n)&&(ind(i)(k) < col_interac(j)))
		{
		  new_ind(nb) = ind(i)(k);
		  new_val(nb) = val(i)(k);
		  k++; nb++;
		}
	      
	      // the new entry
	      new_ind(nb) = col_interac(j);
	      new_val(nb) = val_interac(j); nb++;
	    }
	
	// last entries
	while (k < n)
	  {
	    new_ind(nb) = ind(i)(k);
	    new_val(nb) = val(i)(k);
	    k++; nb++;
	  }
	
	n += nb_new;
	ind(i).Clear(); val(i).Clear();
	ind(i).SetData(n, new_ind.GetData());
	val(i).SetData(n, new_val.GetData());
	new_ind.Nullify(); new_val.Nullify();
	nz_ += nb_new;
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
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cout << (*this)(i,j) << "\t";
	cout << endl;
      }
  }
  
  
  //! Displays a row/column
  /*!
    \param[in] i row/column number
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::PrintVector(int i)
    const
  {
    cout<<" Row "<<i<<endl;
    cout<<" Index "<<endl<<ind(i)<<endl;
    cout<<" Value "<<endl<<val(i)<<endl;
  }
  
  
  //! Assembles the matrix
  /*!
    All the column/row numbers are sorted
    If same column/row numbers exist, values are added
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Assemble()
  {
    for (int i = 0; i < m_; i++)
      AssembleVector(i);
  }
  
  
  //! Removes small coefficients from entries
  /*!
    \param[in] epsilon entries whose values are below epsilon are removed
  */
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  RemoveSmallEntry(const T0& epsilon)
  {
    for (int i = 0; i < ind.GetM(); i++)
      {
	int nb = 0;
	for (int j = 0; j < ind(i).GetM(); j++)
	  if (abs(val(i)(j)) > epsilon)
	    {
	      ind(i)(nb) = ind(i)(j);
	      val(i)(nb) = val(i)(j);
	      nb++;
	    }
	
	if (nb != ind(i).GetM())
	  ResizeVector(i,nb);
      }
  }
  
  
  //! Assembles a row/column
  /*!
    \param[in] i row/column number
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  AssembleVector(int i)
  {
    int new_size = ind(i).GetM();
    Seldon::Assemble(new_size, ind(i), val(i));
    ResizeVector(i,new_size);
  }
  
  
  //! Matrix is initialized to the identity matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::SetIdentity()
  {
    this->n_ = this->m_; this->nz_ = this->m_;
    for (int i = 0; i < this->m_; i++)
      {
	ind(i).Reallocate(1); ind(0) = i;
	val(i).Reallocate(1); val(0) = T(1);
      }
  }
  
  
  //! Non-zero entries are set to 0 (but not removed)
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Zero()
  {
    for (int i = 0; i < val.GetM(); i++)
      val(i).Zero();
  }
  
  
  //! Non-zero entries are filled with values 0, 1, 2, 3 ...
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Fill()
  {
    int value = 0;
    for (int i = 0; i < val.GetM(); i++)
      for (int j = 0; j < val(i).GetM(); j++)
	val(i)(j) = value++;
  }
  
  
  //! Non-zero entries are set to a given value x
  template <class T, class Prop, class Storage, class Allo> template<class T0>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allo>::Fill(const T0& x)
  {
    for (int i = 0; i < val.GetM(); i++)
      val(i).Fill(x);
  }
  
  
  //! Non-zero entries are set to a given value x
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocator>&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);
  }
  
  
  //! Non-zero entries take a random value
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::FillRand()
  {
    for (int i = 0; i < val.GetM(); i++)
      val(i).FillRand();
  }
  
  
  /////////////////////////////
  // MATRIX<ARRAY_COLSPARSE> //
  /////////////////////////////

  
  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSparse, Allocator>::Matrix()  throw():
    Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>()
  {
  }
  
  
  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>(i, j)
  {
  }
  
  
  //! Clears column i
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::ClearColumn(int i)
  {
    this->ClearVector(i);
  }
  
  
  //! Reallocates column i
  /*!
    \param[in] i column number
    \param[in] j new number of non-zero entries in the column
  */
  template <class T, class Prop, class Alloc> inline
  void Matrix<T, Prop, ArrayColSparse, Alloc>::ReallocateColumn(int i,int j)
  {
    this->ReallocateVector(i,j);
  }
  
  
  //! Reallocates column i
  /*!
    \param[in] i column number
    \param[in] j new number of non-zero entries in the column
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::ResizeColumn(int i,int j)
  {
    this->ResizeVector(i,j);
  }
  
  
  //! Swaps two columns
  /*!
    \param[in] i first column number
    \param[in] j second column number
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::SwapColumn(int i,int j)
  {
    this->SwapVector(i,j);
  }
  
  
  //! Sets row numbers of non-zero entries of a column
  /*!
    \param[in] i column number
    \param[in] new_index new row numbers
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  ReplaceIndexColumn(int i, IVect& new_index)
  {
    this->ReplaceIndexVector(i, new_index);
  }
  
  
  //! Returns the number of non-zero entries of a column
  /*!
    \param[in] i column number
    \return The number of non-zero entries of the column i
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColSparse, Allocator>::GetColumnSize(int i) const
  {
    return this->GetVectorSize(i);
  }
  
  
  //! Displays non-zero values of a column
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::PrintColumn(int i) const
  {
    cout<<" Column "<<i<<endl;
    cout<<" Index "<<endl<<this->ind(i)<<endl;
    cout<<" Value "<<endl<<this->val(i)<<endl;
  }
  
  
  //! Assembles a column
  /*!
    \param[in] i column number
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::AssembleColumn(int i)
  {
    this->AssembleVector(i);
  }
  
  
  //! Adds a coefficient in the matrix
  /*!
    \param[in] i row number
    \param[in] j column number
    \param[in] val coefficient to add
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    Matrix_ArraySparse<T,Prop,ArrayColSparse,Allocator>::
      AddInteraction(j, i, val);
  }
  
  
  //! Adds coefficients in a row
  /*!
    \param[in] i row number
    \param[in] nb number of coefficients to add
    \param[in] col_ column numbers of coefficients
    \param[in] val_ values of coefficients
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* val_)
  {
    IVect col; col.SetData(nb, col_);
    Vector<T> val; val.SetData(nb, val_);
    AddInteractionRow(i, nb, col, val);
    col.Nullify(); val.Nullify();
  }
  
  
  //! Adds coefficients in a column
  /*!
    \param[in] i column number
    \param[in] nb number of coefficients to add
    \param[in] row_ row numbers of coefficients
    \param[in] val_ values of coefficients
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* val_)
  {
    IVect row; row.SetData(nb, row_);
    Vector<T> val; val.SetData(nb, val_);
    AddInteractionColumn(i, nb, row, val);
    row.Nullify(); val.Nullify();
  }
  
  
  //! Adds coefficients in a row
  /*!
    \param[in] i row number
    \param[in] nb number of coefficients to add
    \param[in] col column numbers of coefficients
    \param[in] val values of coefficients
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T,Vect_Full,Alloc1>& val)
  {
    for (int j = 0; j < nb; j++)
      Matrix_ArraySparse<T,Prop,ArrayColSparse,Allocator>::
	AddInteraction(col(j), i, val(j));
  }
  
  
  //! Adds coefficients in a column
  /*!
    \param[in] i column number
    \param[in] nb number of coefficients to add
    \param[in] row row numbers of coefficients
    \param[in] val values of coefficients
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<T,Vect_Full,Alloc1>& val)
  {
    Matrix_ArraySparse<T,Prop,ArrayColSparse,Allocator>::
      AddInteractionVector(i, nb, row, val);
  }
  
  
  /////////////////////////////
  // MATRIX<ARRAY_ROWSPARSE> //
  /////////////////////////////
  
  
  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSparse, Allocator>::Matrix()  throw():
    Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>()
  {
  }
  
  
  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>(i, j)
  {
  }
  
  
  //! Clears a row
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::ClearRow(int i)
  {
    this->ClearVector(i);
  }
  
  
  //! Changes the size of a row
  /*!
    \param[in] i row number
    \param[in] j new number of non-zero entries of the row
    \warning data may be lost
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  ReallocateRow(int i, int j)
  {
    this->ReallocateVector(i,j);
  }
  
  
  //! Changes the size of a row
  /*!
    \param[in] i row number
    \param[in] j new number of non-zero entries of the row
    Data is kept
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::ResizeRow(int i, int j)
  {
    this->ResizeVector(i,j);
  }
  
  
  //! Swaps two rows
  /*!
    \param[in] i first row number
    \param[in] j second row number
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::SwapRow(int i,int j)
  {
    this->SwapVector(i,j);
  }
  
  
  //! Sets column numbers of non-zero entries of a row
  /*!
    \param[in] i column number
    \param[in] new_index new column numbers
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  ReplaceIndexRow(int i, IVect& new_index)
  {
    this->ReplaceIndexVector(i,new_index);
  }
  
  
  //! Returns the number of non-zero entries of a row
  /*!
    \param[in] i row number
    \return The number of non-zero entries of the row i
  */
  template <class T, class Prop, class Allocator> inline
  int Matrix<T, Prop, ArrayRowSparse, Allocator>::GetRowSize(int i) const
  {
    return this->GetVectorSize(i);
  }
  
  
  //! Displays non-zero values of a row
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::PrintRow(int i) const
  {
    this->PrintVector(i);
  }
  
  
  //! Assembles a row
  /*!
    \param[in] i row number
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::AssembleRow(int i)
  {
    this->AssembleVector(i);
  }
  
  //! Adds a coefficient in the matrix
  /*!
    \param[in] i row number
    \param[in] j column number
    \param[in] val coefficient to add
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    Matrix_ArraySparse<T,Prop,ArrayRowSparse,Allocator>::
      AddInteraction(i, j, val);
  }
  
  
  //! Adds coefficients in a row
  /*!
    \param[in] i row number
    \param[in] nb number of coefficients to add
    \param[in] col_ column numbers of coefficients
    \param[in] val_ values of coefficients
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* val_)
  {
    IVect col; col.SetData(nb, col_);
    Vector<T> val; val.SetData(nb, val_);
    AddInteractionRow(i, nb, col, val);
    col.Nullify(); val.Nullify();
  }
  
  
  //! Adds coefficients in a column
  /*!
    \param[in] i column number
    \param[in] nb number of coefficients to add
    \param[in] row_ row numbers of coefficients
    \param[in] val_ values of coefficients
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* val_)
  {
    IVect row; row.SetData(nb, row_);
    Vector<T> val; val.SetData(nb, val_);
    AddInteractionColumn(i, nb, row, val);
    row.Nullify(); val.Nullify();
  }
  
  
  //! Adds coefficients in a row
  /*!
    \param[in] i row number
    \param[in] nb number of coefficients to add
    \param[in] col column numbers of coefficients
    \param[in] val values of coefficients
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T,Vect_Full,Alloc1>& val)
  {
    Matrix_ArraySparse<T,Prop,ArrayRowSparse,Allocator>::
      AddInteractionVector(i, nb, col, val);
  }
  
  
  //! Adds coefficients in a column
  /*!
    \param[in] i column number
    \param[in] nb number of coefficients to add
    \param[in] row row numbers of coefficients
    \param[in] val values of coefficients
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<T,Vect_Full,Alloc1>& val)
  {
    for (int j = 0; j < nb; j++)
      Matrix_ArraySparse<T,Prop,ArrayRowSparse,Allocator>::
	AddInteraction(row(j), i, val(j));
  }
  
  
  ////////////////////////////////
  // MATRIX<ARRAY_COLSYMSPARSE> //
  ////////////////////////////////
  
  
  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSymSparse, Allocator>::Matrix()  throw():
    Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>()
  {
  }
  
  
  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>(i, j)
  {
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
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, ArrayColSymSparse, Allocator>::value_type
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::operator() (int i, int j)
    const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if ((j < 0)||(j >= this->n_))
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    int k;
    int a, b;
    
    a = i; b = j;
    Sort(a, b);
    
    for (k = 0; k < this->ind(b).GetM(); k++)
      if (this->ind(b)(k) == a)
	return this->val(b)(k);
    
    return T(0);
  }
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, ArrayColSymSparse, Allocator>::reference
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::operator() (int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if ((j < 0)||(j >= this->n_))
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    Sort(i, j);
    
    int size_vec = this->ind(j).GetM();
    int k = 0;
    // we search the entry
    while ((k < size_vec)&&(this->ind(j)(k) < i))
      k++;
    
    if ((k >= size_vec)||(this->ind(j)(k) != i))
      {
	// the entry does not exist, we add a null value
	AddInteraction(i, j, T(0));
      }
    
    return this->val(j)(k);
  }
  
  
  //! Clears a column
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSymSparse, Allocator>::ClearColumn(int i)
  {
    this->ClearVector(i);
  }
  
  
  //! Reallocates column i
  /*!
    \param[in] i column number
    \param[in] j new number of non-zero entries in the column
    \warning Data may be lost
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ReallocateColumn(int i, int j)
  {
    this->ReallocateVector(i,j);
  }
  
  
  //! Reallocates column i
  /*!
    \param[in] i column number
    \param[in] j new number of non-zero entries in the column
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ResizeColumn(int i, int j)
  {
    this->ResizeVector(i,j);
  }
  
  
  //! Swaps two columns
  /*!
    \param[in] i first column number
    \param[in] j second column number
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  SwapColumn(int i, int j)
  {
    this->SwapVector(i,j);
  }
  
  
  //! Sets row numbers of non-zero entries of a column
  /*!
    \param[in] i column number
    \param[in] new_index new row numbers
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ReplaceIndexColumn(int i, IVect& new_index)
  {
    this->ReplaceIndexVector(i,new_index);
  }
  
  
  //! Returns the number of non-zero entries of a column
  /*!
    \param[in] i column number
    \return The number of non-zero entries of the column i
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  GetColumnSize(int i) const
  {
    return this->GetVectorSize(i);
  }
  
  
  //! Displays non-zero values of a column
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  PrintColumn(int i) const
  {
    cout<<" Column "<<i<<endl;
    cout<<" Index "<<endl<<this->ind(i)<<endl;
    cout<<" Value "<<endl<<this->val(i)<<endl;
  }
  
  
  //! Assembles a column
  /*!
    \param[in] i column number
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AssembleColumn(int i)
  {
    this->AssembleVector(i);
  }
  
  
  //! Adds coefficients in a row
  /*!
    \param[in] i row number
    \param[in] nb number of coefficients to add
    \param[in] col_ column numbers of coefficients
    \param[in] val_ values of coefficients
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* val_)
  {
    IVect col; col.SetData(nb, col_);
    Vector<T> val; val.SetData(nb, val_);
    AddInteractionRow(i, nb, col, val);
    col.Nullify(); val.Nullify();
  }
  
  
  //! Adds coefficients in a column
  /*!
    \param[in] i column number
    \param[in] nb number of coefficients to add
    \param[in] row_ row numbers of coefficients
    \param[in] val_ values of coefficients
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* val_)
  {
    IVect row; row.SetData(nb, row_);
    Vector<T> val; val.SetData(nb, val_);
    AddInteractionColumn(i, nb, row, val);
    row.Nullify(); val.Nullify();
  }
  
  
  //! Adds a coefficient in the matrix
  /*!
    \param[in] i row number
    \param[in] j column number
    \param[in] val coefficient to add
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      Matrix_ArraySparse<T,Prop,ArrayColSymSparse,Allocator>::
	AddInteraction(j, i, val);
  }
  
  
  //! Adds coefficients in a row
  /*!
    \param[in] i row number
    \param[in] nb number of coefficients to add
    \param[in] col column numbers of coefficients
    \param[in] val values of coefficients
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T,Vect_Full,Alloc1>& val)
  {
    for (int j = 0; j < nb; j++)
      if (col(j) <= i)
	Matrix_ArraySparse<T,Prop,ArrayColSymSparse,Allocator>::
	  AddInteraction(col(j), i, val(j));
  }
  
  
  //! Adds coefficients in a column
  /*!
    \param[in] i column number
    \param[in] nb number of coefficients to add
    \param[in] row row numbers of coefficients
    \param[in] val values of coefficients
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<T,Vect_Full,Alloc1>& val)
  {
    IVect new_row(nb); Vector<T, Vect_Full, Alloc1> new_val(nb);
    nb = 0;
    for (int j = 0; j < new_row.GetM(); j++)
      if (row(j) <= i)
	{
	  new_row(nb) = row(j);
	  new_val(nb) = val(j); nb++;
	}
    
    Matrix_ArraySparse<T,Prop,ArrayColSymSparse,Allocator>::
      AddInteractionVector(i, nb, new_row, new_val);
  }
  
  
  ////////////////////////////////
  // MATRIX<ARRAY_ROWSYMSPARSE> //
  ////////////////////////////////
  
  
  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Matrix()  throw():
    Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>()
  {
  }
  
  
  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>(i, j)
  {
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
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, ArrayRowSymSparse, Allocator>::value_type
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::operator() (int i, int j)
    const
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix_ArraySparse::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if ((j < 0)||(j >= this->n_))
      throw WrongCol("Matrix_ArraySparse::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    int k;
    int a, b;
    
    a = i; b = j;
    Sort(a,b);
    
    for (k = 0; k < this->ind(a).GetM(); k++)
      {
	if (this->ind(a)(k) == b)
	  return this->val(a)(k);
      }
    
    return T(0);

  }
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline typename Matrix<T, Prop, ArrayRowSymSparse, Allocator>::reference
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::operator() (int i, int j)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if ((i < 0)||(i >= this->m_))
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if ((j < 0)||(j >= this->n_))
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif
    
    Sort(i, j);
    
    int size_vec = this->ind(i).GetM();
    int k = 0;
    // we search the entry
    while ((k < size_vec)&&(this->ind(i)(k) < j))
      k++;
    
    if ((k >= size_vec)||(this->ind(i)(k) != j))
      {
	// the entry does not exist, we add a null value
	AddInteraction(i, j, T(0));
      }
    
    return this->val(i)(k);
  }
  
  
  //! Clears a row
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::ClearRow(int i)
  {
    this->ClearVector(i);
  }
  
  
  //! Reallocates row i
  /*!
    \param[in] i row number
    \param[in] j new number of non-zero entries in the row
    \warning Data may be lost
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ReallocateRow(int i,int j)
  {
    this->ReallocateVector(i,j);
  }
  
  
  //! Reallocates column i
  /*!
    \param[in] i column number
    \param[in] j new number of non-zero entries in the row
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ResizeRow(int i,int j)
  {
    this->ResizeVector(i,j);
  }
  
  
  //! Swaps two rows
  /*!
    \param[in] i first row number
    \param[in] j second row number
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  SwapRow(int i,int j)
  {
    this->SwapVector(i,j);
  }
  
  
  //! Sets column numbers of non-zero entries of a row
  /*!
    \param[in] i row number
    \param[in] new_index new column numbers
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ReplaceIndexRow(int i,IVect& new_index)
  {
    this->ReplaceIndexVector(i,new_index);
  }
  
  
  //! Returns the number of non-zero entries of a row
  /*!
    \param[in] i row number
    \return The number of non-zero entries of the row i
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayRowSymSparse, Allocator>::GetRowSize(int i)
    const
  {
    return this->GetVectorSize(i);
  }
  
  
  //! Displays non-zero values of a column
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::PrintRow(int i)
    const
  {
    this->PrintVector(i);
  }
  
  
  //! Assembles a column
  /*!
    \param[in] i column number
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::AssembleRow(int i)
  {
    this->AssembleVector(i);
  }
  
  
  //! Adds a coefficient in the matrix
  /*!
    \param[in] i row number
    \param[in] j column number
    \param[in] val coefficient to add
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      Matrix_ArraySparse<T,Prop,ArrayRowSymSparse,Allocator>::
	AddInteraction(i, j, val);
  }
  
  
  //! Adds coefficients in a row
  /*!
    \param[in] i row number
    \param[in] nb number of coefficients to add
    \param[in] col_ column numbers of coefficients
    \param[in] val_ values of coefficients
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* val_)
  {
    IVect col; col.SetData(nb, col_);
    Vector<T> val; val.SetData(nb, val_);
    AddInteractionRow(i, nb, col, val);
    col.Nullify(); val.Nullify();
  }
  
  
  //! Adds coefficients in a column
  /*!
    \param[in] i column number
    \param[in] nb number of coefficients to add
    \param[in] row_ row numbers of coefficients
    \param[in] val_ values of coefficients
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* val_)
  {
    IVect row; row.SetData(nb, row_);
    Vector<T> val; val.SetData(nb, val_);
    AddInteractionColumn(i, nb, row, val);
    row.Nullify(); val.Nullify();
  }
  
  
  //! Adds coefficients in a row
  /*!
    \param[in] i row number
    \param[in] nb number of coefficients to add
    \param[in] col column numbers of coefficients
    \param[in] val values of coefficients
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T,Vect_Full,Alloc1>& val)
  {
    IVect new_col(nb); Vector<T, Vect_Full, Alloc1> new_val(nb);
    nb = 0;
    for (int j = 0; j < new_col.GetM(); j++)
      if (i <= col(j))
	{
	  new_col(nb) = col(j);
	  new_val(nb) = val(j); nb++;
	}
    
    Matrix_ArraySparse<T,Prop,ArrayRowSymSparse,Allocator>::
      AddInteractionVector(i, nb, new_col, new_val);
  }
  
  
  //! Adds coefficients in a column
  /*!
    \param[in] i column number
    \param[in] nb number of coefficients to add
    \param[in] row row numbers of coefficients
    \param[in] val values of coefficients
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<T,Vect_Full,Alloc1>& val)
  {
    // symmetric matrix, row = column
    AddInteractionRow(i, nb, row, val);
  }
  
  
  //! operator<< overloaded for sparse matrices.
  template <class T, class Prop, class Storage, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix_ArraySparse<T, Prop, Storage, Allocator>& A)
  {

    out<<"Value"<<endl;
    out<<A.val<<endl;
    out<<"Index"<<endl;
    out<<A.ind<<endl;
    
    return out;
    
  }

} // namespace Seldon

#define FILE_MATRIX_ARRAY_SPARSE_CXX
#endif

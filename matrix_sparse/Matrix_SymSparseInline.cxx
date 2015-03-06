// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2003-2011 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_SYMSPARSE_INLINE_CXX

#include "Matrix_SymSparse.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::Matrix_SymSparse():
    Matrix_Base<T, Allocator>()
  {
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymSparse(int i, int j): Matrix_Base<T, Allocator>()
  {
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;

    Reallocate(i, j);
  }


  //! Constructor.
  /*! Builds a sparse matrix of size i x j , with nz non-zero (stored)
    elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized. Indices of non-zero entries
    are not initialized either.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Matrix_SymSparse(int i, int j, int nz):
    Matrix_Base<T, Allocator>()
  {
    this->nz_ = 0;
    ind_ = NULL;
    ptr_ = NULL;

    Reallocate(i, j, nz);
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to i so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Matrix_SymSparse(int i, int j,
		   Vector<T, Storage0, Allocator0>& values,
		   Vector<int, Storage1, Allocator1>& ptr,
		   Vector<int, Storage2, Allocator2>& ind):
    Matrix_Base<T, Allocator>(i, j)
  {
    nz_ = values.GetLength();

#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.

    if (ind.GetLength() != nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(nz_) + " values but "
		       + to_str(ind.GetLength()) + " row or column indices.");
      }

    if (ptr.GetLength() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices contains ")
		       + to_str(ptr.GetLength()-1) + string(" row or column ")
		       + string("start  indices (plus the number of non-zero")
		       + " entries) but there are " + to_str(i)
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (static_cast<long int>(2 * nz_ - 2) / static_cast<long int>(i + 1)
	>= static_cast<long int>(i))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_SymSparse::Matrix_SymSparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(values.GetLength())
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

    this->ptr_ = ptr.GetData();
    this->ind_ = ind.GetData();
    this->data_ = values.GetData();

    ptr.Nullify();
    ind.Nullify();
    values.Nullify();
  }


  //! Copy constructor
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::
  Matrix_SymSparse(const Matrix_SymSparse<T, Prop, Storage, Allocator>& A)
  {
    this->m_ = 0;
    this->n_ = 0;
    this->nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>::~Matrix_SymSparse()
  {
    this->Clear();
  }


  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit,
    the matrix is empty (0x0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (ptr_ != NULL)
	  {
	    free(ptr_);
	    ptr_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (ind_ != NULL)
	  {
	    free(ind_);
	    ind_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    this->allocator_.deallocate(this->data_, nz_);
	    this->data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->nz_ = 0;
	this->data_ = NULL;
      }
#endif

    this->nz_ = 0;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory (the number
    of non-zero elements).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_SymSparse<T, Prop, Storage, Allocator>::GetNonZeros() const
  {
    return nz_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory (the number
    of non-zero elements).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_SymSparse<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return nz_;
  }


  //! returns size of matrix in bytes
  template<class T, class Prop, class Storage, class Allocator>
  inline int64_t Matrix_SymSparse<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    int64_t taille = this->GetPtrSize()*sizeof(int);
    int coef = sizeof(T) + sizeof(int); // for each non-zero entry
    taille += coef*int64_t(this->nz_);
    return taille;
  }
  

  //! Returns (row or column) start indices.
  /*!
    Returns the array ('ptr_') of start indices.
    \return The array of start indices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_SymSparse<T, Prop, Storage, Allocator>::GetPtr() const
  {
    return ptr_;
  }


  //! Returns (row or column) indices of non-zero entries.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries. This array defines non-zero entries
    indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_SymSparse<T, Prop, Storage, Allocator>::GetInd() const
  {
    return ind_;
  }


  //! Returns the length of the array of start indices.
  /*!
    \return The length of the array of start indices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_SymSparse<T, Prop, Storage, Allocator>::GetPtrSize() const
  {
    return (this->m_ + 1);
  }


  //! Returns the length of the array of (column or row) indices.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored). This array defines non-zero
    entries indices (that are stored) if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_SymSparse<T, Prop, Storage, Allocator>::GetIndSize() const
  {
    return nz_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access method.
  /*! Returns reference to element (\a i, \a j) 
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Matrix_SymSparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_SymSparse<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    return Val(i, j);
  }
  
  
  //! Add a value to a non-zero entry.
  /*! This function adds \a val to the element (\a i, \a j), provided that
    this element is already a non-zero entry. Otherwise 
    a non-zero entry is inserted equal to \a val.
    \param[in] i row index.
    \param[in] j column index.
    \param[in] val value to be added to the element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      Get(i, j) += val;
  }


  //! Sets an element (i, j) to a value
  /*! This function sets \a val to the element (\a i, \a j)
    \param[in] i row index.
    \param[in] j column index.
    \param[in] val A(i, j) = val
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymSparse<T, Prop, Storage, Allocator>&
  Matrix_SymSparse<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_SymSparse<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }

  
  //////////////////////////
  // MATRIX<COLSYMSPARSE> //
  //////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymSparse, Allocator>::Matrix():
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>(i, j, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with nz non-zero (stored) elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymSparse, Allocator>::Matrix(int i, int j, int nz):
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>(i, j, nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row start indices.
    \param ind column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix<T, Prop, ColSymSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& values,
	 Vector<int, Storage1, Allocator1>& ptr,
	 Vector<int, Storage2, Allocator2>& ind):
    Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>(i, j, values, ptr, ind)
  {
  }



  //////////////////////////
  // MATRIX<ROWSYMSPARSE> //
  //////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymSparse, Allocator>::Matrix():
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>(i, j, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with nz non-zero (stored) elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements that are stored.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymSparse, Allocator>::Matrix(int i, int j, int nz):
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>(i, j, nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr column start indices.
    \param ind row indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix<T, Prop, RowSymSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& values,
	 Vector<int, Storage1, Allocator1>& ptr,
	 Vector<int, Storage2, Allocator2>& ind):
    Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>(i, j, values, ptr, ind)
  {
  }
    
  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMSPARSE_INLINE_CXX
#endif

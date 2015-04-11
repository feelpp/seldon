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


#ifndef SELDON_FILE_MATRIX_COMPLEXSPARSE_INLINE_CXX

#include "Matrix_ComplexSparse.hxx"

namespace Seldon
{

  
  inline int ColComplexSparse::GetFirst(int i, int j)
  {
    return j;
  }
  inline int ColComplexSparse::GetSecond(int i, int j)
  {
    return i;
  }
  inline int ColComplexSparse::GetBeginLoop(int i)
  {
    return 0;
  }


  inline int RowComplexSparse::GetFirst(int i, int j)
  {
    return i;
  }
  inline int RowComplexSparse::GetSecond(int i, int j)
  {
    return j;
  }
  inline int RowComplexSparse::GetBeginLoop(int i)
  {
    return 0;
  }

  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_ComplexSparse(): Matrix_Base<T, Allocator>()
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_ComplexSparse(int i, int j): Matrix_Base<T, Allocator>()
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
    
    Reallocate(i, j);
  }


  //! Constructor.
  /*! Builds a sparse matrix of size i by j , with real_nz
    non-zero elements in the real part of the matrix and imag_nz
    non-zero elements in the imaginary part of the matrix.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements in the real part.
    \param imag_nz number of non-zero elements in the imaginary part.
    \note Matrix values are not initialized. Indices of non-zero entries
    are not initialized either.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_ComplexSparse(int i, int j, int real_nz, int imag_nz):
    Matrix_Base<T, Allocator>()
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
    
    Reallocate(i, j, real_nz, imag_nz);
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_ComplexSparse(int i, int j,
		       Vector<value_type, Storage0, Allocator0>& real_values,
		       Vector<int, Storage1, Allocator1>& real_ptr,
		       Vector<int, Storage2, Allocator2>& real_ind,
		       Vector<value_type, Storage0, Allocator0>& imag_values,
		       Vector<int, Storage1, Allocator1>& imag_ptr,
		       Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_Base<T, Allocator>(i, j)
  {

    real_nz_ = real_values.GetLength();
    imag_nz_ = imag_values.GetLength();

#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.

    if (real_ind.GetLength() != real_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(real_nz_)
		       + " values (real part) but "
		       + to_str(real_ind.GetLength())
		       + " row or column indices.");
      }

    if (imag_ind.GetLength() != imag_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(imag_nz_)
		       + " values (imaginary part) but "
		       + to_str(imag_ind.GetLength())
		       + " row or column indices.");
      }

    if (real_ptr.GetLength()-1 != Storage::GetFirst(i, j))
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (real part)")
		       + " contains " + to_str(real_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(Storage::GetFirst(i, j))
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(j) + " matrix).");
      }

    if (imag_ptr.GetLength()-1 != Storage::GetFirst(i, j))
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (imaginary part)")
		       + " contains " + to_str(imag_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(Storage::GetFirst(i, j))
		       + " rows or columns ("
		       + to_str(i) + " by " + to_str(j) + " matrix).");
      }

    if ((real_nz_ > 0
         && (j == 0
             || static_cast<long int>(real_nz_-1) / static_cast<long int>(j)
             >= static_cast<long int>(i)))
        ||
        (imag_nz_ > 0
         && (j == 0
             || static_cast<long int>(imag_nz_-1) / static_cast<long int>(j)
             >= static_cast<long int>(i))))
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_ComplexSparse::")
		       + string("Matrix_ComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(real_values.GetLength())
		       + " values for the real part and "
		       + to_str(real_values.GetLength()) + string(" values")
		       + string(" for the imaginary part) than elements")
		       + " in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

    this->real_ptr_ = real_ptr.GetData();
    this->imag_ptr_ = imag_ptr.GetData();
    this->real_ind_ = real_ind.GetData();
    this->imag_ind_ = imag_ind.GetData();
    this->real_data_ = real_values.GetData();
    this->imag_data_ = imag_values.GetData();

    real_ptr.Nullify();
    imag_ptr.Nullify();
    real_ind.Nullify();
    imag_ind.Nullify();
    real_values.Nullify();
    imag_values.Nullify();
  }


  //! Copy constructor
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_ComplexSparse(const Matrix_ComplexSparse<T, Prop,
			 Storage, Allocator>& A)
  {
    this->m_ = 0;
    this->n_ = 0;
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;

    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::~Matrix_ComplexSparse()
  {
    this->Clear();
  }


  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix
    is empty (0x0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ComplexSparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (real_ptr_ != NULL)
	  {
	    free(real_ptr_);
	    real_ptr_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	real_ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (imag_ptr_ != NULL)
	  {
	    free(imag_ptr_);
	    imag_ptr_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	imag_ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (real_ind_ != NULL)
	  {
	    free(real_ind_);
	    real_ind_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	real_ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (imag_ind_ != NULL)
	  {
	    free(imag_ind_);
	    imag_ind_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	imag_ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->real_data_ != NULL)
	  {
	    this->allocator_.deallocate(this->real_data_, real_nz_);
	    this->real_data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->real_nz_ = 0;
	this->real_data_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->imag_data_ != NULL)
	  {
	    this->allocator_.deallocate(this->imag_data_, imag_nz_);
	    this->imag_data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->imag_nz_ = 0;
	this->imag_data_ = NULL;
      }
#endif

    this->real_nz_ = 0;
    this->imag_nz_ = 0;
  }

  
  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the cumulated number of non-zero entries of both the real and
    the imaginary part.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return real_nz_ + imag_nz_;
  }


  //! returns size of A in bytes used to store the matrix
  template<class T, class Prop, class Storage, class Allocator>
  inline int64_t Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetMemorySize() const
  {
    int64_t taille = 2*this->GetRealPtrSize()*sizeof(int);
    int coef = sizeof(value_type) + sizeof(int); // for each non-zero entry
    taille += coef*int64_t(this->real_nz_ + this->imag_nz_);
    return taille;
  }
  
  
  //! Returns (row or column) start indices for the real part.
  /*!
    Returns the array ('ptr_') of start indices for the real part.
    \return The array of start indices for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetRealPtr() const
  {
    return real_ptr_;
  }


  //! Returns (row or column) start indices for the imaginary part.
  /*!
    Returns the array ('ptr_') of start indices for the imaginary part.
    \return The array of start indices for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetImagPtr() const
  {
    return imag_ptr_;
  }


  //! Returns (row or column) indices of non-zero entries for the real part.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries for the real part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetRealInd() const
  {
    return real_ind_;
  }


  //! Returns (row or column) indices of non-zero entries
  //! for the imaginary part.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries for the imaginary part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetImagInd() const
  {
    return imag_ind_;
  }


  //! Returns the length of the array of start indices for the real part.
  /*!
    \return The length of the array of start indices for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealPtrSize() const
  {
    return (Storage::GetFirst(this->m_, this->n_) + 1);
  }


  //! Returns the length of the array of start indices for the imaginary part.
  /*!
    \return The length of the array of start indices for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagPtrSize() const
  {
    return (Storage::GetFirst(this->m_, this->n_) + 1);
  }


  //! Returns the length of the array of (column or row) indices
  //! for the real part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries for the real part. This array defines non-zero entries
    indices if coupled with (column or row) start indices.
    \return The length of the array of (column or row) indices
    for the real part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealIndSize() const
  {
    return real_nz_;
  }


  //! Returns the length of the array of (column or row) indices
  //! for the imaginary part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries for the imaginary part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The length of the array of (column or row) indices
    for the imaginary part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagIndSize() const
  {
    return imag_nz_;
  }

  
  //! Returns the length of the array of (column or row) indices
  //! for the real part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries for the real part. This array defines non-zero entries
    indices if coupled with (column or row) start indices.
    \return The length of the array of (column or row) indices
    for the real part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealDataSize() const
  {
    return real_nz_;
  }


  //! Returns the length of the array of (column or row) indices
  //! for the imaginary part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries for the imaginary part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The length of the array of (column or row) indices
    for the imaginary part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagDataSize() const
  {
    return imag_nz_;
  }

  
  //! Returns the array of values of the real part.
  /*!
    \return The array 'real_data_' of values of the real part..
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type*
  Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetRealData() const
  {
    return real_data_;
  }


  //! Returns the array of values of the imaginary part.
  /*!
    \return The array 'imag_data_' of values of the imaginary part..
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Allocator::value_type*
  Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetImagData() const
  {
    return imag_data_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access method
  /*! Returns the real part of element (\a i, \a j)
    if it can be returned as a reference. If the non-zero entry
    does not exit, it is created
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetReal(int i, int j) const
  {
    return ValReal(i, j);
  }
  
  
  //! Access method
  /*! Returns the imaginary part of element (\a i, \a j)
    if it can be returned as a reference. If the non-zero entry
    does not exit, it is created
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::value_type&
  Matrix_ComplexSparse<T, Prop, Storage, Allocator>::GetImag(int i, int j) const
  {
    return ValImag(i, j);
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
  inline void Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const entry_type& val)
  {
    if (real(val) != value_type(0))
      GetReal(i, j) += real(val);

    if (imag(val) != value_type(0))
      GetImag(i, j) += imag(val);
  }


  //! Adds values to several non-zero entries on a given row
  template <class T, class Prop, class Storage, class Allocator>
  template<class Alloc1>
  inline void Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::AddInteractionRow(int i, int nb, const Vector<int>& col,
		      const Vector<entry_type, VectFull, Alloc1>& val)
  {
    throw Undefined("AddInteractionRow", "Not implemented");
  }


  //! Sets an element (i, j) to a value
  /*! This function sets \a val to the element (\a i, \a j)
    \param[in] i row index.
    \param[in] j column index.
    \param[in] val A(i, j) = val
  */  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const entry_type& val)
  {
    GetReal(i, j) = real(val);
    GetImag(i, j) = imag(val);
  }

  
  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ComplexSparse<T, Prop, Storage, Allocator>&
  Matrix_ComplexSparse<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_ComplexSparse<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }

  
  //////////////////////////////
  // MATRIX<COLCOMPLEXSPARSE> //
  //////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColComplexSparse, Allocator>::Matrix() :
    Matrix_ComplexSparse<T, Prop, ColComplexSparse, Allocator>()
  {
  }


  //! Builds a i by j matrix.
  /*!
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColComplexSparse, Allocator>::Matrix(int i, int j):
    Matrix_ComplexSparse<T, Prop, ColComplexSparse, Allocator>(i, j, 0, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with real_nz and imag_nz non-zero elements for
    the real part and the imaginary part respectively.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements for the real part.
    \param imag_nz number of non-zero elements for the imaginary part.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColComplexSparse, Allocator>::Matrix(int i, int j,
							      int real_nz,
							      int imag_nz):
    Matrix_ComplexSparse<T, Prop, ColComplexSparse, Allocator>(i, j,
							       real_nz,
							       imag_nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix<T, Prop, ColComplexSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<value_type, Storage0, Allocator0>& real_values,
	 Vector<int, Storage1, Allocator1>& real_ptr,
	 Vector<int, Storage2, Allocator2>& real_ind,
	 Vector<value_type, Storage0, Allocator0>& imag_values,
	 Vector<int, Storage1, Allocator1>& imag_ptr,
	 Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_ComplexSparse<T, Prop, ColComplexSparse, Allocator>(i, j,
							       real_values,
							       real_ptr,
							       real_ind,
							       imag_values,
							       imag_ptr,
							       imag_ind)
  {
  }



  //////////////////////////////
  // MATRIX<ROWCOMPLEXSPARSE> //
  //////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowComplexSparse, Allocator>::Matrix() :
    Matrix_ComplexSparse<T, Prop, RowComplexSparse, Allocator>()
  {
  }


  //! Builds a i by j matrix.
  /*!
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator> inline
  Matrix<T, Prop, RowComplexSparse, Allocator>::Matrix(int i, int j):
    Matrix_ComplexSparse<T, Prop, RowComplexSparse, Allocator>(i, j, 0, 0)
  {
  }


  /*! Builds a i by j matrix with real_nz and imag_nz non-zero elements for
    the real part and the imaginary part respectively.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements for the real part.
    \param imag_nz number of non-zero elements for the imaginary part.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator> inline
  Matrix<T, Prop, RowComplexSparse, Allocator>::Matrix(int i, int j,
						       int real_nz,
						       int imag_nz):
    Matrix_ComplexSparse<T, Prop, RowComplexSparse, Allocator>(i, j,
							       real_nz,
							       imag_nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix<T, Prop, RowComplexSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<value_type, Storage0, Allocator0>& real_values,
	 Vector<int, Storage1, Allocator1>& real_ptr,
	 Vector<int, Storage2, Allocator2>& real_ind,
	 Vector<value_type, Storage0, Allocator0>& imag_values,
	 Vector<int, Storage1, Allocator1>& imag_ptr,
	 Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_ComplexSparse<T, Prop, RowComplexSparse, Allocator>(i, j,
							       real_values,
							       real_ptr,
							       real_ind,
							       imag_values,
							       imag_ptr,
							       imag_ind)
  {
  }
  
  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_COMPLEXSPARSE_INLINE_CXX
#endif

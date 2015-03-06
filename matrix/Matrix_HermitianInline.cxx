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


#ifndef SELDON_FILE_MATRIX_HERMITIAN_INLINE_CXX

#include "Matrix_Hermitian.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>::Matrix_Hermitian():
    Matrix_Base<T, Allocator>()
  {
    me_ = NULL;
  }


  //! Main constructor.
  /*! Builds a i x j full matrix.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Matrix_Hermitian(int i, int j):
    Matrix_Base<T, Allocator>(i, i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	me_ = reinterpret_cast<pointer*>(calloc(i, sizeof(pointer)));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
	this->data_ = NULL;
      }
    if (me_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
      }
    if ((me_ == NULL) && (i != 0))
      throw NoMemory("Matrix_Hermitian::Matrix_Hermitian(int, int)",
		     string("Unable to allocate memory for a matrix of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(i)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(i)
		     + " elements).");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->allocator_.allocate(i * i, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(me_);
	me_ = NULL;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(me_);
	me_ = NULL;
      }
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Matrix_Hermitian::Matrix_Hermitian(int, int)",
		     string("Unable to allocate memory for a matrix of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(i)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(i)
		     + " elements).");
#endif

    pointer ptr = this->data_;
    int lgth = i;
    for (int k = 0; k < i; k++, ptr += lgth)
      me_[k] = ptr;
  }


  //! Copy constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Matrix_Hermitian(const Matrix_Hermitian<T, Prop, Storage, Allocator>& A):
    Matrix_Base<T, Allocator>()
  {
    this->m_ = 0;
    this->n_ = 0;
    this->data_ = NULL;
    this->me_ = NULL;

    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>::~Matrix_Hermitian()
  {
    Clear();
  }


  //! Clears the matrix.
  /*!
    Destructs the matrix.
    \warning On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>::Clear()
  {
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    this->allocator_.deallocate(this->data_, this->m_ * this->n_);
	    this->data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (me_ != NULL)
	  {
	    free(me_);
	    me_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
      }
#endif
    
    this->m_ = 0;
    this->n_ = 0;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of rows multiplied by the number of columns
    because the matrix is full.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_Hermitian<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return this->m_ * this->n_;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x i matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j)
  {

    if (i != this->m_)
      {
	this->m_ = i;
	this->n_ = i;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    me_ = reinterpret_cast<pointer*>( realloc(me_,
						      i * sizeof(pointer)) );

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    me_ = NULL;
	    this->data_ = NULL;
	  }
	if (me_ == NULL)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    this->data_ = NULL;
	  }
	if (me_ == NULL && i != 0)
	  throw NoMemory("Matrix_Hermitian::Reallocate(int, int)",
			 string("Unable to reallocate memory")
			 + " for a matrix of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(i)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(i)
			 + " elements).");
#endif

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(this->allocator_.reallocate(this->data_,
								    i * i,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    free(me_);
	    me_ = NULL;
	    this->data_ = NULL;
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    free(me_);
	    me_ = NULL;
	  }
	if (this->data_ == NULL && i != 0)
	  throw NoMemory("Matrix_Hermitian::Reallocate(int, int)",
			 string("Unable to reallocate memory")
			 + " for a matrix of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(i)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(i)
			 + " elements).");
#endif

	pointer ptr = this->data_;
	int lgth = Storage::GetSecond(i, i);
	for (int k = 0; k < Storage::GetFirst(i, i); k++, ptr += lgth)
	  me_[k] = ptr;
      }
  }


  //! Changes the size of the matrix and sets its data array
  //! (low level method).
  /*!
    The matrix is first cleared (memory is freed). The matrix is then resized
    to a i x j matrix, and the data array of the matrix is set to 'data'.
    'data' elements are not duplicated: the new data array of the matrix is
    the 'data' array. It is useful to create a matrix from pre-existing data.
    \param i new number of rows.
    \param j new number of columns.
    \param data new array storing elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The matrix
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::SetData(int i, int j,
	    typename Matrix_Hermitian<T, Prop, Storage, Allocator>
	    ::pointer data)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = i;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	me_ = reinterpret_cast<pointer*>(calloc(i, sizeof(pointer)));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
	this->data_ = NULL;
	return;
      }
    if (me_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
	return;
      }
#endif

    this->data_ = data;

    pointer ptr = this->data_;
    int lgth = i;
    for (int k = 0; k < i; k++, ptr += lgth)
      me_[k] = ptr;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
    \warning Memory is not released except for me_.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (me_ != NULL)
	  {
	    free(me_);
	    me_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
      }
#endif

    this->data_ = NULL;
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
  inline const typename
  Matrix_Hermitian<T, Prop, Storage, Allocator>::value_type
  Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::operator() const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::operator() const",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    if (i > j)
      return conj(me_[Storage::GetSecond(i, j)][Storage::GetFirst(i, j)]);
    else
      return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];

  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");    
    if (i > j)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }

  
  //! Returns access to an element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_Hermitian::Get(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Returns the element (i, j).
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
#endif
    
    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }
  
  
  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::operator[] (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_Hermitian::operator[] (int)",
		       string("Index should be in [0, ")
		       + to_str(this->GetDataSize()-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::operator[] (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_Hermitian::operator[] (int) const",
		       string("Index should be in [0, ")
		       + to_str(this->GetDataSize()-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif

    return this->data_[i];
  }

  
  //! Sets an element of the matrix
  /*!
    \param i row index
    \param j column index
    \param x sets a(i, j) = x
   */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& x)
  {
    if (i > j)
      this->Val(j, i) = conj(x);
    else
      this->Val(i, j) = x;
  }
  

  //! Duplicates a matrix (assignement operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>&
  Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::operator=(const Matrix_Hermitian<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Copy(const Matrix_Hermitian<T, Prop, Storage, Allocator>& A)
  {
    this->Reallocate(A.GetM(), A.GetN());

    this->allocator_.memorycpy(this->data_, A.GetData(), this->GetDataSize());
  }

  
  ////////////////////
  // MATRIX<COLHERM> //
  ////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColHerm, Allocator>::Matrix():
    Matrix_Hermitian<T, Prop, ColHerm, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j full column-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColHerm, Allocator>::Matrix(int i, int j):
    Matrix_Hermitian<T, Prop, ColHerm, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills a matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }

  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator= (const Matrix<T, Prop,
                                                  ColHerm,
                                                  Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
    \warning The imaginary part of x should be null to keep an Hermitian
    matrix.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator*= (const T0& x)
  {
    for (int i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }



  ////////////////////
  // MATRIX<ROWHERM> //
  ////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowHerm, Allocator>::Matrix():
    Matrix_Hermitian<T, Prop, RowHerm, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i by j full row-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowHerm, Allocator>::Matrix(int i, int j):
    Matrix_Hermitian<T, Prop, RowHerm, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills a matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }

  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator= (const Matrix<T, Prop,
                                                       RowHerm,
                                                       Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator*= (const T0& x)
  {
    for (int i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }

} // namespace Seldon.

#define SELDON_FILE_MATRIX_HERMITIAN_INLINE_CXX
#endif

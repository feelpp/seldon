// Copyright (C) 2001-2009 Vivien Mallet
// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_ARRAY3D_INLINE_CXX

#include "Array3D.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the array is an empty 0x0x0 3D array.
  */
  template <class T, class Allocator>
  inline Array3D<T, Allocator>::Array3D()
  {
    length1_ = 0;
    length2_ = 0;
    length3_ = 0;

    length23_ = 0;

    data_ = NULL;
  }


  //! Main constructor.
  /*! Builds a i x j x k 3D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
  */
  template <class T, class Allocator>
  inline Array3D<T, Allocator>::Array3D(int i, int j, int k)
  {
    length1_ = i;
    length2_ = j;
    length3_ = k;

    length23_ = length2_ * length3_;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	data_ = array3D_allocator_.allocate(i*j*k, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length1_ = 0;
	length2_ = 0;
	length3_ = 0;
	length23_ = 0;
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0)
      throw NoMemory("Array3D::Array3D(int, int, int)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " elements).");
#endif

  }


  //! Copy constructor.
  template <class T, class Allocator>
  inline Array3D<T, Allocator>::Array3D(const Array3D<T, Allocator>& A)
  {
    length1_ = 0;
    length2_ = 0;
    length3_ = 0;

    length23_ = 0;

    data_ = NULL;

    Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Array3D<T, Allocator>::~Array3D()
  {
    length1_ = 0;
    length2_ = 0;
    length3_ = 0;
    length23_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (data_ != NULL)
	  {
	    array3D_allocator_.deallocate(data_, length1_ * length23_);
	    data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	data_ = NULL;
      }
#endif

  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the 3D array.
  /*!
    On exit, the array is a i x j x k 3D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Allocator>
  inline void Array3D<T, Allocator>::Reallocate(int i, int j, int k)
  {
    if (i != length1_ || j != length2_ || k != length3_)
      {
	length1_ = i;
	length2_ = j;
	length3_ = k;

	length23_ = j * k;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif
            
	    data_ =
	      reinterpret_cast<pointer>(array3D_allocator_.reallocate(data_,
								      i*j*k,
								      this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    length1_ = 0;
	    length2_ = 0;
	    length3_ = 0;
	    length23_ = 0;
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0)
	  throw NoMemory("Array3D::Reallocate(int, int, int)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " elements).");
#endif

      }
  }


  //! Changes the size of the array and sets its data array
  //! (low level method).
  /*!
    The 3D array is first cleared (memory is freed). The 3D array is then
    resized to a i x j x k, and the data array of the 3D array is set to
    'data'.  'data' elements are not duplicated: the new data array of the 3D
    array is the 'data' array. It is useful to create a 3D array from
    pre-existing data.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param data new array storing elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The matrix
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline void Array3D<T, Allocator>::SetData(int i, int j, int k, 
					     typename Array3D<T, Allocator>
					     ::pointer data)
  {
    Clear();
    
    length1_ = i;
    length2_ = j;
    length3_ = k;
    length23_ = j * k;
    data_ = data;
  }


  //! Clears the 3D array without releasing memory.
  /*!
    On exit, the 3D array is empty and the memory has not been released.
    It is useful for low level manipulations on a 3D arrat instance.
  */
  template <class T, class Allocator>
  inline void Array3D<T, Allocator>::Nullify()
  {
    length1_ = 0;
    length2_ = 0;
    length3_ = 0;
    length23_ = 0;
    data_ = NULL;
  }


  //! Clears the array.
  /*!
    Destructs the array.
    \warning On exit, the 3D array is empty.
  */
  template <class T, class Allocator>
  inline void Array3D<T, Allocator>::Clear()
  {
    this->~Array3D();
    this->length1_ = 0;
    this->length2_ = 0;
    this->length3_ = 0;
  }


  /*****************
   * BASIC METHODS *
   *****************/


  //! Returns the length in dimension #1.
  /*!
    \return The length in dimension #1.
  */
  template <class T, class Allocator>
  inline int Array3D<T, Allocator>::GetLength1() const
  {
    return length1_;
  }


  //! Returns the length in dimension #2.
  /*!
    \return The length in dimension #2.
  */
  template <class T, class Allocator>
  inline int Array3D<T, Allocator>::GetLength2() const
  {
    return length2_;
  }


  //! Returns the length in dimension #3.
  /*!
    \return The length in dimension #3.
  */
  template <class T, class Allocator>
  inline int Array3D<T, Allocator>::GetLength3() const
  {
    return length3_;
  }


  //! Returns the number of elements in the 3D array.
  /*!
    Returns the number of elements stored by the 3D array, i.e.
    the product of the lengths in the three dimensions.
    \return The number of elements in the 3D array.
  */
  template <class T, class Allocator>
  inline int Array3D<T, Allocator>::GetSize() const
  {
    return length1_ * length23_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory by
    the array, i.e. the product of lengths in the three
    dimensions.
    \return The number of elements stored in the array.
  */
  template <class T, class Allocator>
  inline int Array3D<T, Allocator>::GetDataSize() const
  {
    return length1_ * length23_;
  }


  //! Returns a pointer to the data array.
  /*!
    Returns a pointer to data, i.e. the data array 'data_' which stores the
    values.
    \return A pointer to the data array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::pointer Array3D<T, Allocator>
  ::GetData() const
  {
    return data_;
  }


  //! Returns a pointer to an element of data array.
  /*!
    Returns a pointer to an element of data array.
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return A pointer to the data array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::pointer Array3D<T, Allocator>
  ::GetDataPointer(int i, int j, int k) const
  {
    return data_ + i * length23_ + j * length3_ + k;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::reference
  Array3D<T, Allocator>::operator() (int i, int j, int k)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length1_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length1_-1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length2_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length2_-1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length3_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length3_-1) + "], but is equal to "
		       + to_str(k) + ".");
#endif

    return data_[i * length23_ + j * length3_ + k];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, class Allocator>
  inline typename Array3D<T, Allocator>::const_reference
  Array3D<T, Allocator>::operator() (int i, int j, int k) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length1_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length1_-1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length2_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length2_-1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length3_)
      throw WrongIndex("Array3D::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length3_-1) + "], but is equal to "
		       + to_str(k) + ".");
#endif

    return data_[i*length23_ + j*length3_ + k];
  }

  //! Duplicates a 3D array (assignment operator).
  /*!
    \param A 3D array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Array3D<T, Allocator>& Array3D<T, Allocator>::operator=
  (const Array3D<T, Allocator>& A)
  {
    Copy(A);

    return *this;
  }

  //! Duplicates a 3D array.
  /*!
    \param A 3D array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Array3D<T, Allocator>::Copy(const Array3D<T, Allocator>& A)
  {
    Reallocate(A.GetLength1(), A.GetLength2(), A.GetLength3());

    array3D_allocator_.memorycpy(data_, A.GetData(), GetDataSize());
  }
  
} // namespace Seldon.

#define SELDON_FILE_ARRAY3D_INLINE_CXX
#endif

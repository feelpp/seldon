// Copyright (C) 2010 Lin Wu
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


#ifndef SELDON_FILE_ARRAY_INLINE_CXX


namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the array is an empty array.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>::Array()
  {
    if (N < ARRAY_MINRANK || N > ARRAY_MAXRANK)
      {
	string msg = string("Array dimension should be in [") +
	  to_str(ARRAY_MINRANK) + string(", ") + to_str(ARRAY_MAXRANK) + "].";
	throw WrongDim("Array<T, N, Allocator>::Array(int, ...)", msg);
      }
    length_ = NULL;
    offset_ = NULL;
    data_ = NULL;
  }


  //! Main constructor.
  /*! Builds a i x j x k 3D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>::Array(int i, int j, int k)
  {
    if (N != 3)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 3.");

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	length_ = length_allocator_.allocate(N, this);
	offset_ = length_allocator_.allocate(N, this);
	data_ = array_allocator_.allocate(i*j*k, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length_ = NULL;
	offset_ = NULL;
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    offset_[0] = j * k;
    offset_[1] = k;
    offset_[2] = i * j * k;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l 4D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>::Array(int i, int j, int k, int l)
  {
    if (N != 4)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 4.");

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	length_ = length_allocator_.allocate(N, this);
	offset_ = length_allocator_.allocate(N, this);
	data_ = array_allocator_.allocate(i*j*k*l, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length_ = NULL;
	offset_ = NULL;
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    offset_[0] = j * k * l;
    offset_[1] = k * l;
    offset_[2] = l;
    offset_[3] = i * j * k * l;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m 5D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>::Array(int i, int j, int k, int l, int m)
  {
    if (N != 5)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 5.");

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	length_ = length_allocator_.allocate(N, this);
	offset_ = length_allocator_.allocate(N, this);
	data_ = array_allocator_.allocate(i*j*k*l*m, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length_ = NULL;
	offset_ = NULL;
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x "
		     + to_str(m) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    offset_[0] = j * k * l * m;
    offset_[1] = k * l * m;
    offset_[2] = l * m;
    offset_[3] = m;
    offset_[4] = i * j * k * l * m;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m x n 6D array, but data is not initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>
  ::Array(int i, int j, int k, int l, int m, int n)
  {
    if (N != 6)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 6.");

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	length_ = length_allocator_.allocate(N, this);
	offset_ = length_allocator_.allocate(N, this);
	data_ = array_allocator_.allocate(i*j*k*l*m*n, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length_ = NULL;
	offset_ = NULL;
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	n != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(n)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x " + 
		     to_str(m) + " x " + to_str(n) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    length_[5] = n;
    offset_[0] = j * k * l * m * n;
    offset_[1] = k * l * m * n;
    offset_[2] = l * m * n;
    offset_[3] = m * n;
    offset_[4] = n;
    offset_[5] = i * j * k * l * m * n;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m x n x o 7D array, but data is not
    initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>
  ::Array(int i, int j, int k, int l, int m, int n, int o)
  {
    if (N != 7)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 7.");

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	length_ = length_allocator_.allocate(N, this);
	offset_ = length_allocator_.allocate(N, this);
	data_ = array_allocator_.allocate(i*j*k*l*m*n*o, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length_ = NULL;
	offset_ = NULL;
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	n != 0 && o != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(n)
			      * static_cast<long int>(o)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x " + 
		     to_str(m) + " x " + to_str(n) + " x " + to_str(o) +
		     " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    length_[5] = n;
    length_[6] = o;
    offset_[0] = j * k * l * m * n * o;
    offset_[1] = k * l * m * n * o;
    offset_[2] = l * m * n * o;
    offset_[3] = m * n * o;
    offset_[4] = n * o;
    offset_[5] = o;
    offset_[6] = i * j * k * l * m * n * o;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m x n x o x p 8D array, but data is not
    initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \param p length in dimension #8.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>
  ::Array(int i, int j, int k, int l, int m, int n, int o, int p)
  {
    if (N != 8)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 8.");

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	length_ = length_allocator_.allocate(N, this);
	offset_ = length_allocator_.allocate(N, this);
	data_ = array_allocator_.allocate(i*j*k*l*m*n*o*p, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length_ = NULL;
	offset_ = NULL;
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	n != 0 && o != 0 && p != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(n)
			      * static_cast<long int>(o)
			      * static_cast<long int>(p)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x " + 
		     to_str(m) + " x " + to_str(n) + " x " + to_str(o)
		     + " x " + to_str(p) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    length_[5] = n;
    length_[6] = o;
    length_[7] = p;
    offset_[0] = j * k * l * m * n * o * p;
    offset_[1] = k * l * m * n * o * p;
    offset_[2] = l * m * n * o * p;
    offset_[3] = m * n * o * p;
    offset_[4] = n * o * p;
    offset_[5] = o * p;
    offset_[6] = p;
    offset_[7] = i * j * k * l * m * n * o * p;

  }


  //! Main constructor.
  /*! Builds a i x j x k x l x m x n x o x p x q 9D array, but data is not
    initialized.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \param p length in dimension #8.
    \param q length in dimension #9.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>
  ::Array(int i, int j, int k, int l, int m, int n, int o, int p, int q)
  {
    if (N != 9)
      throw WrongDim("Array<T, N, Allocator>::Array(int, ...)",
		     "Array dimension should be 9.");

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	length_ = length_allocator_.allocate(N, this);
	offset_ = length_allocator_.allocate(N, this);
	data_ = array_allocator_.allocate(i*j*k*l*m*n*o*p*q, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length_ = NULL;
	offset_ = NULL;
	data_ = NULL;
      }
    if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	n != 0 && o != 0 && p != 0 && q != 0)
      throw NoMemory("Array::Array(int, ...)",
		     string("Unable to allocate memory for an array of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(k)
			      * static_cast<long int>(l)
			      * static_cast<long int>(m)
			      * static_cast<long int>(n)
			      * static_cast<long int>(o)
			      * static_cast<long int>(p)
			      * static_cast<long int>(q)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " x " + to_str(k) + " x " + to_str(l) + " x " + 
		     to_str(m) + " x " + to_str(n) + " x " + to_str(o) + " x "
		     + to_str(p) + " x " + to_str(q) + " elements).");
#endif

    length_[0] = i;
    length_[1] = j;
    length_[2] = k;
    length_[3] = l;
    length_[4] = m;
    length_[5] = n;
    length_[6] = o;
    length_[7] = p;
    length_[8] = q;
    offset_[0] = j * k * l * m * n * o * p * q;
    offset_[1] = k * l * m * n * o * p * q;
    offset_[2] = l * m * n * o * p * q;
    offset_[3] = m * n * o * p * q;
    offset_[4] = n * o * p * q;
    offset_[5] = o * p * q;
    offset_[6] = p * q;
    offset_[7] = q;
    offset_[8] = i * j * k * l * m * n * o * p * q;

  }


  //! Copy constructor.
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>::Array(const Array<T, N, Allocator>& A)
  {
    length_ = NULL;
    offset_ = NULL;
    data_ = NULL;

    Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>::~Array()
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (data_ != NULL)
	  {
	    array_allocator_.deallocate(data_, offset_[N-1]);
	    data_ = NULL;
	  }
	if (length_ != NULL)
	  {
	    length_allocator_.deallocate(length_, N);
	    length_ = NULL;
	  }
	if (offset_ != NULL)
	  {
	    length_allocator_.deallocate(offset_, N);
	    offset_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	length_ = NULL;
	offset_ = NULL;
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
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>::Reallocate(int i, int j, int k)
  {
    if (N != 3)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 3.");

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    length_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(length_,
								  N,
								  this));
	    offset_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(offset_,
								  N,
								  this));
	    data_ =
	      reinterpret_cast<pointer>(array_allocator_.reallocate(data_,
								    i*j*k,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    length_ = NULL;
	    offset_ = NULL;
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " elements).");
#endif

	length_[0] = i;
	length_[1] = j;
	length_[2] = k;
	offset_[0] = j * k;
	offset_[1] = k;
	offset_[2] = i * j * k;
      }
  }


  //! Reallocates memory to resize the 4D array.
  /*!
    On exit, the array is a i x j x k x l 4D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>::Reallocate(int i, int j, int k, int l)
  {
    if (N != 4)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 4.");

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    length_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(length_,
								  N,
								  this));
	    offset_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(offset_,
								  N,
								  this));
	    data_ =
	      reinterpret_cast<pointer>(array_allocator_.reallocate(data_,
								    i*j*k*l,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    length_ = NULL;
	    offset_ = NULL;
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " elements).");
#endif

	length_[0] = i;
	length_[1] = j;
	length_[2] = k;
	length_[3] = l;
	offset_[0] = j * k * l;
	offset_[1] = k * l;
	offset_[2] = l;
	offset_[3] = i * j * k * l;
      }
  }


  //! Reallocates memory to resize the 5D array.
  /*!
    On exit, the array is a i x j x k x l x m 5D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m)
  {
    if (N != 5)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 5.");

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    length_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(length_,
								  N,
								  this));
	    offset_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(offset_,
								  N,
								  this));
	    data_ =
	      reinterpret_cast<pointer>(array_allocator_.reallocate(data_,
								    i*j*k*l*m,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    length_ = NULL;
	    offset_ = NULL;
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " elements).");
#endif

	length_[0] = i;
	length_[1] = j;
	length_[2] = k;
	length_[3] = l;
	length_[4] = m;
	offset_[0] = j * k * l * m;
	offset_[1] = k * l * m;
	offset_[2] = l * m;
	offset_[3] = m;
	offset_[4] = i * j * k * l * m;
      }
  }


  //! Reallocates memory to resize the 6D array.
  /*!
    On exit, the array is a i x j x k x l x m x n 6D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m, int n)
  {
    if (N != 6)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 6.");

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4]
	|| n != length_[5])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    length_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(length_,
								  N,
								  this));
	    offset_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(offset_,
								  N,
								  this));
	    data_ =
	      reinterpret_cast<pointer>(array_allocator_.reallocate(data_,
								    i*j*k*l
								    *m*n,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    length_ = NULL;
	    offset_ = NULL;
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	    n != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(n)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " x " + to_str(n)
			 + " elements).");
#endif

	length_[0] = i;
	length_[1] = j;
	length_[2] = k;
	length_[3] = l;
	length_[4] = m;
	length_[5] = n;
	offset_[0] = j * k * l * m * n;
	offset_[1] = k * l * m * n;
	offset_[2] = l * m * n;
	offset_[3] = m * n;
	offset_[4] = n;
	offset_[5] = i * j * k * l * m * n;
      }
  }


  //! Reallocates memory to resize the 7D array.
  /*!
    On exit, the array is a i x j x k x l x m x n x o 7D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m, int n, int o)
  {
    if (N != 7)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 7.");

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4]
	|| n != length_[5] || o != length_[6])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    length_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(length_,
								  N,
								  this));
	    offset_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(offset_,
								  N,
								  this));
	    data_ =
	      reinterpret_cast<pointer>(array_allocator_.reallocate(data_,
								    i*j*k*l
								    *m*n*o,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    length_ = NULL;
	    offset_ = NULL;
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	    n != 0 && o != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(n)
				  * static_cast<long int>(o)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " x " + to_str(n)
			 + " x " + to_str(o) + " elements).");
#endif

	length_[0] = i;
	length_[1] = j;
	length_[2] = k;
	length_[3] = l;
	length_[4] = m;
	length_[5] = n;
	length_[6] = o;
	offset_[0] = j * k * l * m * n * o;
	offset_[1] = k * l * m * n * o;
	offset_[2] = l * m * n * o;
	offset_[3] = m * n * o;
	offset_[4] = n * o;
	offset_[5] = o;
	offset_[6] = i * j * k * l * m * n * o;
      }
  }


  //! Reallocates memory to resize the 8D array.
  /*!
    On exit, the array is a i x j x k x l x m x n x o x p 8D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \param p length in dimension #8.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m, int n, int o, int p)
  {
    if (N != 8)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 8.");

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4]
	|| n != length_[5] || o != length_[6] || p != length_[7])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    length_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(length_,
								  N,
								  this));
	    offset_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(offset_,
								  N,
								  this));
	    data_ =
	      reinterpret_cast<pointer>(array_allocator_.reallocate(data_,
								    i*j*k*l
								    *m*n*o*p,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    length_ = NULL;
	    offset_ = NULL;
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	    n != 0 && o != 0 && p != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(n)
				  * static_cast<long int>(o)
				  * static_cast<long int>(p)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " x " + to_str(n)
			 + " x " + to_str(o) + " x " + to_str(p)
			 + " elements).");
#endif

	length_[0] = i;
	length_[1] = j;
	length_[2] = k;
	length_[3] = l;
	length_[4] = m;
	length_[5] = n;
	length_[6] = o;
	length_[7] = p;
	offset_[0] = j * k * l * m * n * o * p;
	offset_[1] = k * l * m * n * o * p;
	offset_[2] = l * m * n * o * p;
	offset_[3] = m * n * o * p;
	offset_[4] = n * o * p;
	offset_[5] = o * p;
	offset_[6] = p;
	offset_[7] = i * j * k * l * m * n * o * p;
      }
  }


  //! Reallocates memory to resize the 9D array.
  /*!
    On exit, the array is a i x j x k x l x m x n x o x p x q 9D array.
    \param i length in dimension #1.
    \param j length in dimension #2.
    \param k length in dimension #3.
    \param l length in dimension #4.
    \param m length in dimension #5.
    \param n length in dimension #6.
    \param o length in dimension #7.
    \param p length in dimension #8.
    \param q length in dimension #9.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>
  ::Reallocate(int i, int j, int k, int l, int m, int n, int o, int p, int q)
  {
    if (N != 9)
      throw WrongDim("Array::Reallocate(int, ...)",
		     "Array dimension should be 9.");

    if (length_ == NULL || i != length_[0] || j != length_[1]
	|| k != length_[2] || l != length_[3] || m != length_[4]
	|| n != length_[5] || o != length_[6] || p != length_[7]
	|| q != length_[8])
      {

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    length_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(length_,
								  N,
								  this));
	    offset_ = reinterpret_cast<length_pointer>(length_allocator_.
						       reallocate(offset_,
								  N,
								  this));
	    data_ =
	      reinterpret_cast<pointer>(array_allocator_.reallocate(data_,
								    i*j*k*l
								    *m*n*o*
								    p*q,
								    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    length_ = NULL;
	    offset_ = NULL;
	    data_ = NULL;
	  }
	if (data_ == NULL && i != 0 && j != 0 && k != 0 && l != 0 && m != 0 &&
	    n != 0 && o != 0 && p != 0 && q != 0)
	  throw NoMemory("Array::Reallocate(int, ...)",
			 string("Unable to reallocate memory")
			 + " for an array of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(k)
				  * static_cast<long int>(l)
				  * static_cast<long int>(m)
				  * static_cast<long int>(n)
				  * static_cast<long int>(o)
				  * static_cast<long int>(p)
				  * static_cast<long int>(q)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " x " + to_str(k) + " x " + to_str(l)
			 + " x " + to_str(m) + " x " + to_str(n)
			 + " x " + to_str(o) + " x " + to_str(p)
			 + " x " + to_str(q) + " elements).");
#endif

	length_[0] = i;
	length_[1] = j;
	length_[2] = k;
	length_[3] = l;
	length_[4] = m;
	length_[5] = n;
	length_[6] = o;
	length_[7] = p;
	length_[8] = q;
	offset_[0] = j * k * l * m * n * o * p * q;
	offset_[1] = k * l * m * n * o * p * q;
	offset_[2] = l * m * n * o * p * q;
	offset_[3] = m * n * o * p * q;
	offset_[4] = n * o * p * q;
	offset_[5] = o * p * q;
	offset_[6] = p * q;
	offset_[7] = q;
	offset_[8] = i * j * k * l * m * n * o * p * q;
      }
  }


  //! Clears the array.
  /*!
    Destructs the array.
    \warning On exit, the array is empty.
  */
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>::Clear()
  {
    this->~Array();
  }


  /*****************
   * BASIC METHODS *
   *****************/

  //! Returns the length in dimension #1.
  /*!
    \param dimension index for dimension.
    \return The length in dimension #1.
  */
  template <class T, int N, class Allocator>
  inline int Array<T, N, Allocator>::GetLength(int dimension) const
  {
    return length_[dimension];
  }


  //! Returns the number of elements in the 3D array.
  /*!
    Returns the number of elements stored by the 3D array, i.e.
    the product of the lengths in the three dimensions.
    \return The number of elements in the 3D array.
  */
  template <class T, int N, class Allocator>
  inline int Array<T, N, Allocator>::GetSize() const
  {
    if (offset_ == NULL)
      return 0;
    else
      return offset_[N - 1];
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory by
    the array, i.e. the product of lengths in the three
    dimensions.
    \return The number of elements stored in the array.
  */
  template <class T, int N, class Allocator>
  inline int Array<T, N, Allocator>::GetDataSize() const
  {
    if (offset_ == NULL)
      return 0;
    else
      return offset_[N - 1];
  }


  //! Returns a pointer to the data array.
  /*!
    Returns a pointer to data, i.e. the data array 'data_' which stores the
    values.
    \return A pointer to the data array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::pointer Array<T, N, Allocator>
  ::GetData() const
  {
    return data_;
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
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>::operator() (int i, int j, int k)
  {
    if (N != 3)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 3.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \return Element (i, j, k) of the 3D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>::operator() (int i, int j, int k) const
  {
    if (N != 3)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 3.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \return Element (i, j, k, l) of the 4D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l)
  {
    if (N != 4)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 4.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] + l];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \return Element (i, j, k, l) of the 4D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l) const
  {
    if (N != 4)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 4.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] + l];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \return Element (i, j, k, l, m) of the 5D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l, int m)
  {
    if (N != 5)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 5.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \return Element (i, j, k, l, m) of the 5D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l, int m) const
  {
    if (N != 5)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 5.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \return Element (i, j, k, l, m, n) of the 6D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n)
  {
    if (N != 6)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 6.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length_[5])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length_[5] - 1) + "], but is equal to "
		       + to_str(n) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \return Element (i, j, k, l, m, n) of the 6D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n) const
  {
    if (N != 6)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 6.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length_[5])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length_[5] - 1) + "], but is equal to "
		       + to_str(n) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \return Element (i, j, k, l, m, n, o) of the 7D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n, int o)
  {
    if (N != 7)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 7.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length_[5])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length_[5] - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length_[6])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length_[6] - 1) + "], but is equal to "
		       + to_str(o) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] + o];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \return Element (i, j, k, l, m, n, o) of the 7D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n, int o) const
  {
    if (N != 7)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 7.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length_[5])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length_[5] - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length_[6])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length_[6] - 1) + "], but is equal to "
		       + to_str(o) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] + o];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o, p).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \param p index along dimension #8.
    \return Element (i, j, k, l, m, n, o, p) of the 8D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n, int o, int p)
  {
    if (N != 8)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 8.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length_[5])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length_[5] - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length_[6])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length_[6] - 1) + "], but is equal to "
		       + to_str(o) + ".");
    if (p < 0 || p >= length_[7])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #8 should be in [0, ")
		       + to_str(length_[7] - 1) + "], but is equal to "
		       + to_str(p) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] +
		 o * offset_[6] + p];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o, p).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \param p index along dimension #8.
    \return Element (i, j, k, l, m, n, o, p) of the 8D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>
  ::operator() (int i, int j, int k, int l, int m, int n, int o, int p) const
  {
    if (N != 8)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 8.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length_[5])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length_[5] - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length_[6])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length_[6] - 1) + "], but is equal to "
		       + to_str(o) + ".");
    if (p < 0 || p >= length_[7])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #8 should be in [0, ")
		       + to_str(length_[7] - 1) + "], but is equal to "
		       + to_str(p) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] +
		 o * offset_[6] + p];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o, p, q).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \param p index along dimension #8.
    \param q index along dimension #9.
    \return Element (i, j, k, l, m, n, o, p, q) of the 9D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l, int m,
				      int n, int o, int p, int q)
  {
    if (N != 9)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 9.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length_[5])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length_[5] - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length_[6])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length_[6] - 1) + "], but is equal to "
		       + to_str(o) + ".");
    if (p < 0 || p >= length_[7])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #8 should be in [0, ")
		       + to_str(length_[7] - 1) + "], but is equal to "
		       + to_str(p) + ".");
    if (q < 0 || q >= length_[8])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #9 should be in [0, ")
		       + to_str(length_[8] - 1) + "], but is equal to "
		       + to_str(q) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] +
		 o * offset_[6] + p * offset_[7] + q];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j, k, l, m, n, o, p, q).
    \param i index along dimension #1.
    \param j index along dimension #2.
    \param k index along dimension #3.
    \param l index along dimension #4.
    \param m index along dimension #5.
    \param n index along dimension #6.
    \param o index along dimension #7.
    \param p index along dimension #8.
    \param q index along dimension #9.
    \return Element (i, j, k, l, m, n, o, p, q) of the 9D array.
  */
  template <class T, int N, class Allocator>
  inline typename Array<T, N, Allocator>::const_reference
  Array<T, N, Allocator>::operator() (int i, int j, int k, int l, int m,
				      int n, int o, int p, int q) const
  {
    if (N != 9)
      throw WrongDim("Array<T, N, Allocator>::operator(int, ...)",
		     "Array dimension should be 9.");

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= length_[0])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #1 should be in [0, ")
		       + to_str(length_[0] - 1) + "], but is equal to "
		       + to_str(i) + ".");
    if (j < 0 || j >= length_[1])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #2 should be in [0, ")
		       + to_str(length_[1] - 1) + "], but is equal to "
		       + to_str(j) + ".");
    if (k < 0 || k >= length_[2])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #3 should be in [0, ")
		       + to_str(length_[2] - 1) + "], but is equal to "
		       + to_str(k) + ".");
    if (l < 0 || l >= length_[3])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #4 should be in [0, ")
		       + to_str(length_[3] - 1) + "], but is equal to "
		       + to_str(l) + ".");
    if (m < 0 || m >= length_[4])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #5 should be in [0, ")
		       + to_str(length_[4] - 1) + "], but is equal to "
		       + to_str(m) + ".");
    if (n < 0 || n >= length_[5])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #6 should be in [0, ")
		       + to_str(length_[5] - 1) + "], but is equal to "
		       + to_str(n) + ".");
    if (o < 0 || o >= length_[6])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #7 should be in [0, ")
		       + to_str(length_[6] - 1) + "], but is equal to "
		       + to_str(o) + ".");
    if (p < 0 || p >= length_[7])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #8 should be in [0, ")
		       + to_str(length_[7] - 1) + "], but is equal to "
		       + to_str(p) + ".");
    if (q < 0 || q >= length_[8])
      throw WrongIndex("Array::operator()",
		       string("Index along dimension #9 should be in [0, ")
		       + to_str(length_[8] - 1) + "], but is equal to "
		       + to_str(q) + ".");
#endif

    return data_[i * offset_[0] + j * offset_[1] + k * offset_[2] +
		 l * offset_[3] + m * offset_[4] + n * offset_[5] +
		 o * offset_[6] + p * offset_[7] + q];
  }


  //! Duplicates an array (assignment operator).
  /*!
    \param A array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, int N, class Allocator>
  inline Array<T, N, Allocator>& Array<T, N, Allocator>::operator=
  (const Array<T, N, Allocator>& A)
  {
    Copy(A);

    return *this;
  }
  
} // namespace Seldon.

#define SELDON_FILE_ARRAY_INLINE_CXX
#endif

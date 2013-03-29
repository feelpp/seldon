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


#ifndef SELDON_FILE_ARRAY_CXX


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
  Array<T, N, Allocator>::Array(const Array<T, N, Allocator>& A)
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
  int Array<T, N, Allocator>::GetLength(int dimension) const
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
  int Array<T, N, Allocator>::GetSize() const
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
  int Array<T, N, Allocator>::GetDataSize() const
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
  typename Array<T, N, Allocator>::pointer Array<T, N, Allocator>
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

  //! Duplicates an array.
  /*!
    \param A array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, int N, class Allocator>
  inline void Array<T, N, Allocator>::Copy(const Array<T, N, Allocator>& A)
  {
    if (N == 3)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2));
    if (N == 4)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3));
    if (N == 5)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4));
    if (N == 6)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5));
    if (N == 7)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6));
    if (N == 8)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6), A.GetLength(7));
    if (N == 9)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6), A.GetLength(7), A.GetLength(8));

    array_allocator_.memorycpy(data_, A.GetData(), GetDataSize());
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the array stores complex
    structures, use 'Fill' instead.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Zero()
  {
    array_allocator_.memoryset(data_, char(0),
			       GetDataSize()*sizeof(value_type));
  }


  //! Fills the array.
  /*!
    On exit, the array is filled with 1, 2, 3, 4, ... The order of
    those numbers depends on the storage.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Fill()
  {
    for (int i = 0; i < GetDataSize(); i++)
      data_[i] = i;
  }


  //! Fills the array with a given value.
  /*!
    On exit, the array is filled with 'x'.
    \param x the value to fill the array with.
  */
  template <class T, int N, class Allocator>
  template <class T0>
  void Array<T, N, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < GetDataSize(); i++)
      data_[i] = x;
  }


  //! Fills the array randomly.
  /*!
    On exit, the array is filled with random values.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::FillRand()
  {
    srand(time(NULL));
    for (int i = 0; i < GetDataSize(); i++)
      data_[i] = rand();
  }


  //! Displays the array on the standard output.
  /*!
    Displays elements on the standard output, in text format.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Print() const
  {
    if (N == 3)
      {
	int i, j, k;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  cout << (*this)(i, j, k) << '\t';
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 4)
      {
	int i, j, k, l;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      cout << (*this)(i, j, k, l) << '\t';
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 5)
      {
	int i, j, k, l, m;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  cout << (*this)(i, j, k, l, m) << '\t';
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 6)
      {
	int i, j, k, l, m, n;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      cout << (*this)(i, j, k, l, m, n) << '\t';
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 7)
      {
	int i, j, k, l, m, n, o;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  cout << (*this)(i, j, k, l, m, n, o)
				       << '\t';
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 8)
      {
	int i, j, k, l, m, n, o, p;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  {
				    for (p = 0; p < GetLength(7); p++)
				      cout << (*this)(i, j, k, l, m, n, o, p)
					   << '\t';
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 9)
      {
	int i, j, k, l, m, n, o, p, q;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  {
				    for (p = 0; p < GetLength(7); p++)
				      {
					for (q = 0; q < GetLength(8); q++)
					  cout << (*this)(i, j, k, l, m, n,
							  o, p, q)
					       << '\t';
					cout << endl;
				      }
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the array in a file.
  /*!
    Stores the array in a file in binary format. The dimensions (integer)
    are written, and array elements are then written in the same order as in
    memory.
    \param FileName output file name.
  */
  template <class T, int N, class Allocator> void Array<T, N, Allocator>
  ::Write(string FileName, bool with_size = true) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Array::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the array to an output stream.
  /*!
    Writes the array to an output stream in binary format. The dimensions
    (integerS) are written, and array elements are then written in the same
    order as in memory.
    \param FileStream output stream.
  */
  template <class T, int N, class Allocator> void Array<T, N, Allocator>
  ::Write(ofstream& FileStream, bool with_size = true) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Array::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    if (with_size)
      {
	for (int i = 0; i < N; i++)
	  FileStream.write(reinterpret_cast<char*>
			   (const_cast<int*>(&length_[i])),
			   sizeof(int));
      }

    FileStream.write(reinterpret_cast<char*>(data_),
		     offset_[N - 1] * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Array::Write(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string("  The output file may have been removed")
		    + " or there is no space left on device.");
#endif

  }


  //! Reads the array from a file.
  /*!
    Reads a array stored in binary format in a file.  The dimensions
    (integers) of the array are read, and array elements are then read in the
    same order as it should be in memory
    \param FileName input file name.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Read(string FileName, bool with_size = true)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Array::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    Read(FileStream, with_size);

    FileStream.close();
  }


  //! Reads the array from an input stream.
  /*!
    Reads a array in binary format from an input stream.  The dimensions
    (integers) of the array are read, and array elements are then read in the
    same order as it should be in memory.
    \param FileStream input stream.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>
  ::Read(ifstream& FileStream, bool with_size = true)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Read(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    if (with_size)
      {
	if (N == 3)
	  {
	    int new_l1, new_l2, new_l3;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3);
	  }
	if (N == 4)
	  {
	    int new_l1, new_l2, new_l3, new_l4;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4);
	  }
	if (N == 5)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5);
	  }
	if (N == 6)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6);
	  }
	if (N == 7)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7);
	  }
	if (N == 8)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7,
	      new_l8;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l8), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7, new_l8);
	  }
	if (N == 9)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7,
	      new_l8, new_l9;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l8), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l9), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7, new_l8, new_l9);
	  }
      }

    FileStream.read(reinterpret_cast<char*>(data_), 
		    offset_[N - 1] * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Array::Read(ifstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The intput file may have been removed")
		    + " or may not contain enough data.");
#endif

  }

  //! operator<< overloaded for a 3D array.
  /*!
    \param out output stream.
    \param A the 3D array.
    \return The updated stream.
  */
  template <class T, int N, class Allocator>
  ostream& operator << (ostream& out,
			const Array<T, N, Allocator>& A)
  {
    if (N == 3)
      {
	int i, j, k;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  out << A(i, j, k) << '\t';
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 4)
      {
	int i, j, k, l;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      out << A(i, j, k, l) << '\t';
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 5)
      {
	int i, j, k, l, m;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  out << A(i, j, k, l, m) << '\t';
			out << endl;
		      }
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 6)
      {
	int i, j, k, l, m, n;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      out << A(i, j, k, l, m, n) << '\t';
			    out << endl;
			  }
			out << endl;
		      }
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 7)
      {
	int i, j, k, l, m, n, o;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  cout << A(i, j, k, l, m, n, o)
				       << '\t';
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 8)
      {
	int i, j, k, l, m, n, o, p;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  {
				    for (p = 0; p < A.GetLength(7); p++)
				      cout << A(i, j, k, l, m, n, o, p)
					   << '\t';
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 9)
      {
	int i, j, k, l, m, n, o, p, q;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  {
				    for (p = 0; p < A.GetLength(7); p++)
				      {
					for (q = 0; q < A.GetLength(8); q++)
					  cout << A(i, j, k, l, m, n, o, p, q)
					       << '\t';
					cout << endl;
				      }
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    return out;
  }


} // namespace Seldon.

#define SELDON_FILE_ARRAY_CXX
#endif

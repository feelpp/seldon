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

#ifndef SELDON_FILE_SPARSE_VECTOR_CXX

#include "SparseVector.hxx"

namespace Seldon
{
  
  
  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  Vector<T, Vect_Sparse, Allocator>::Vector()  throw():
    Vector<T, Vect_Full, Allocator>()
  {
    index_ = NULL;
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param i length of the vector.
  */
  template <class T, class Allocator>
  Vector<T, Vect_Sparse, Allocator>::Vector(int i):
    Vector<T, Vect_Full, Allocator>(i)
  {
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	  
	this->index_ = index_allocator_.allocate(i, this);
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->index_ = NULL;
	this->data_ = NULL;
      }
    
    if (this->index_ == NULL)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Vector<Vect_Sparse>::Vector(int)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(i*sizeof(T)) + " bytes ("
		     + to_str(i) + " elements).");
#endif

  }
  
  
  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param V vector to be copied.
  */
  template <class T, class Allocator>
  Vector<T, Vect_Sparse, Allocator>::
  Vector(const Vector<T, Vect_Sparse, Allocator>& V) :
    Vector<T, Vect_Full, Allocator>()
  {
    this->index_ = NULL;
    Copy(V);
  }

  
  /**************
   * DESTRUCTOR *
   **************/
  
  
  //! Destructor.
  template <class T, class Allocator>
  Vector<T, Vect_Sparse, Allocator>::~Vector()
  {
    // data_ is released
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	if (this->data_ != NULL)
	  {
	    this->vect_allocator_.deallocate(this->data_, this->m_);
	    this->data_ = NULL;
	  }
	
	if (index_ != NULL)
	  {
	    index_allocator_.deallocate(index_, this->m_);
	    index_ = NULL;
	  }
	
	this->m_ = 0;
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
	index_ = NULL;
	this->m_ = 0;
	return;
      }
#endif

  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Clears the vector.
  /*!
    Destructs the vector.
    \warning On exit, the vector is an empty vector.
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Sparse, Allocator>::Clear()
  {
    this->~Vector();
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Previous non-zero entries are removed
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Sparse, Allocator>::Reallocate(int i)
  {
    if (i != this->m_)
      {

	this->m_ = i;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(this->vect_allocator_.reallocate(this->data_,
									 i, this) );
										      
	    index_
	      = reinterpret_cast<int*>(this->index_allocator_.reallocate(index_,
									 i, this) );

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->data_ = NULL;
	    this->index_ = NULL;
	    return;
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    this->index_ = NULL;
	    return;
	  }
#endif

      }
  }
  
  
  //! Changes the number of non-zero entries of the vector
  /*!
    Changes the number of non-zero entries to i. Previous values are kept.
    \param n new number of non-zero entries of the vector.
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Sparse, Allocator>::Resize(int n)
  {
    if (n == this->m_)
      return;
    
    Vector<T, Vect_Full, Allocator> new_value(n);
    Vector<int> new_index(n);
    int nmin = min(this->m_, n);
    for (int i = 0; i < nmin; i++)
      {
	new_value(i) = this->data_[i];
	new_index(i) = index_[i];
      }
    
    SetData(new_value, new_index);
  }
  
  
  //! Changes the length of the vector and sets its data array
  //! (low level method).
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param i new length of the vector.
    \param data the new data array. 'data' contains the new elements of the
    vector and must therefore contain 'i' elements.
    \param index the new index array. 'index' contains the new row numbers
    of the vector and must therefore contain 'i' elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The vector
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Sparse, Allocator>
  ::SetData(int i, T* data, int* index)
  {
    this->Clear();

    this->m_ = i;
    
    this->data_ = data;
    this->index_ = index;
  }
  
  
  //! Changes the length of the vector and sets its data array
  //! (low level method).
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param data the new data array. 'data' contains the values
    \param index the new index array. 'index' contains the new row numbers
    of the vector and must have the same size than data.
    \note vectors data and index are empty after the method
  */
  template <class T, class Allocator> template<class Allocator2>
  void Vector<T, Vect_Sparse, Allocator>::SetData(Vector<T, Vect_Full, Allocator2>& data,
						  Vector<int>& index)
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (data.GetM() != index.GetM())
      throw WrongDim(string("Vector<Vect_Sparse>::SetData ") +
		     " Vectors data and index should have the same size "
		     + " size of data : " + to_str(data.GetM())
		     + "size of index " + to_str(index.GetM()) );
#endif

    SetData(data.GetM(), data.GetData(), index.GetData());
    data.Nullify();
    index.Nullify();
  }

  //! Clears the vector without releasing memory.
  /*!
    On exit, the vector is empty and the memory has not been released.
    It is useful for low level manipulations on a Vector instance.
    \warning Memory is not released.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->data_ = NULL;
    this->index_ = NULL;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/

  
  //! Access operator.
  /*!
    \param i index.
    \return The value of the non-zero element i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, Vect_Sparse, Allocator>::reference
  Vector<T, Vect_Sparse, Allocator>::Value(int i)
  {
    
#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<Vect_Sparse>::Value(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif
    
    return this->data_[i];
  }
  
  
  //! Access operator.
  /*!
    \param i index.
    \return The value of the non-zero element i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, Vect_Sparse, Allocator>::const_reference
  Vector<T, Vect_Sparse, Allocator>::Value(int i) const
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<Vect_Sparse>::Value(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif
    
    return this->data_[i];
  }
  
  
  //! Access operator.
  /*!
    \param i index.
    \return The row number of the non-zero element i.
  */
  template <class T, class Allocator>
  inline int& Vector<T, Vect_Sparse, Allocator>::Index(int i)
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<Vect_Sparse>::Index(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif
    
    return this->index_[i];
  }
  
  
  //! Access operator.
  /*!
    \param i index.
    \return The row number of the non-zero element i.
  */
  template <class T, class Allocator>
  inline int Vector<T, Vect_Sparse, Allocator>::Index(int i) const
  {

#ifdef SELDON_CHECK_BOUNDARIES
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<Vect_Sparse>::Index(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif
    
    return this->index_[i];
  }
  
  
  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator>
  inline typename Vector<T, Vect_Sparse, Allocator>::reference
  Vector<T, Vect_Sparse, Allocator>::operator() (int i)
  {
    int k = 0;
    // we search the entry
    while (k < this->m_ && index_[k] < i)
      k++;
    
    if (k >= this->m_ || index_[k] != i)
      // The entry does not exist, we add a null value.
      AddInteraction(i, T(0));
    
    return this->data_[k];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator>
  inline typename Vector<T, Vect_Sparse, Allocator>::value_type
  Vector<T, Vect_Sparse, Allocator>::operator() (int i) const
  {
    int k = 0;
    // we search the entry
    while (k < this->m_ && index_[k] < i)
      k++;
    
    if (k >= this->m_ || index_[k] != i)
      // The entry does not exist, we return a null value
      return T(0);
    
    return this->data_[k];
  }
  
  
  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, Vect_Sparse, Allocator>& Vector<T, Vect_Sparse, Allocator>
  ::operator= (const Vector<T, Vect_Sparse, Allocator>& X)
  {
    this->Copy(X);

    return *this;
  }
  

  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, Vect_Sparse, Allocator>
  ::Copy(const Vector<T, Vect_Sparse, Allocator>& X)
  {
    this->Reallocate(X.GetLength());
    
    this->vect_allocator_.memorycpy(this->data_, X.GetData(), this->m_);
    this->index_allocator_.memorycpy(this->index_, X.GetIndex(), this->m_);
  }
  
  
  /*******************
   * BASIC FUNCTIONS *
   *******************/
  

  //! Returns a pointer to array containing row numbers
  template <class T, class Allocator>
  int* Vector<T, Vect_Sparse, Allocator>::GetIndex() const
  {
    return this->index_;
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/

  
  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  Vector<T, Vect_Sparse, Allocator>&
  Vector<T, Vect_Sparse, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }
  
  
  //! Displays the vector.
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::Print() const
  {
    for (int i = 0; i < this->GetLength(); i++)
      cout << (Index(i)+1) << ' ' <<Value(i) << '\n';
  }
  
  
  //! Assembles the vector
  /*!
    \warning If you are using the methods AddInteraction,
    you don't need to call that method
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::Assemble()
  {
    int new_size = this->m_;
    Vector<T, Vect_Full, Allocator> values(new_size);
    Vector<int> index(new_size);
    for (int i = 0; i < new_size; i++)
      {
	values(i) = this->data_[i];
	index(i) = index_[i];
      }
    
    Seldon::Assemble(new_size, index, values);
    index.Resize(new_size);
    values.Resize(new_size);
    SetData(values, index);
  }
  
  
  //! Removes small entries
  template <class T, class Allocator> template<class T0>
  void Vector<T, Vect_Sparse, Allocator>::RemoveSmallEntry(const T0& epsilon)
  {
    int new_size = this->m_;
    Vector<T, Vect_Full, Allocator> values(new_size);
    Vector<int> index(new_size);
    new_size = 0;
    for (int i = 0; i < this->m_; i++)
      if (abs(this->data_[i]) > epsilon)
	{
	  values(new_size) = this->data_[i];
	  index(new_size) = index_[i];
	  new_size++;
	}
    
    index.Resize(new_size);
    values.Resize(new_size);
    SetData(values, index);
  }
  
  
  //! Coefficient a is added in the vector.
  /*!
    \param[in] i row number.
    \param[in] a value to add.
  */
  template <class T, class Allocator> inline
  void Vector<T, Vect_Sparse, Allocator>::AddInteraction(int i, const T& val)
  {
    // We look for the position where the entry is.
    int pos = 0;
    while (pos < this->m_ && index_[pos] < i)
      pos++;
    
    // If the entry already exists, we add a.
    if (pos < this->m_)
      if (index_[pos] == i)
	{
	  this->data_[pos] += val;
	  return;
	}
    
    // The interaction doesn't exist, the vector is reallocated.
    Vector<T, Vect_Full, Allocator> new_val(this->m_+1);
    Vector<int> new_ind(this->m_+1);
    for (int k = 0; k < pos; k++)
      {
	new_ind(k) = index_[k];
	new_val(k) = this->data_[k];
      }
    
    // The new entry.
    new_ind(pos) = i;
    new_val(pos) = val;
    
    // End of the row/column.
    for (int k = pos+1; k <= this->m_; k++)
      {
	new_ind(k) = index_[k-1];
	new_val(k) = this->data_[k-1];
      }
    
    SetData(new_val, new_ind);
  }
  
  //! Coefficients are added in the vector
  /*!
    The method sorts given coefficients and adds them
    in the correct positions.
    \param[in] n number of coefficients to add.
    \param[in] row row numbers.
    \param[in] values coefficients.
    \param[in] already_sorted true if row numbers are already sorted
  */
  template <class T, class Allocator> inline
  void Vector<T, Vect_Sparse, Allocator>::
  AddInteractionRow(int n, int* row, T* values, bool already_sorted = false)
  {
    Vector<int> ind;
    Vector<T, Vect_Full, Allocator> val;
    ind.SetData(n, row);
    val.SetData(n, values);
    AddInteractionRow(n, ind, val, already_sorted);
    ind.Nullify();
    val.Nullify();
  }
  
  
  //! Coefficients are added in the vector
  /*!
    The method sorts given coefficients and adds them
    in the correct positions.
    \param[in] nb_interac number of coefficients to add.
    \param[in] col_interac row numbers.
    \param[in] val_interac coefficients.
    \param[in] already_sorted true if row numbers are already sorted
  */
  template <class T, class Allocator> template<class Alloc1>
  void Vector<T, Vect_Sparse, Allocator>::
  AddInteractionRow(int nb_interac, Vector<int> col_interac,
		    Vector<T, Vect_Full, Alloc1> val_interac,
		    bool already_sorted = false)
  {
    // There is no reference in arguments, because we want to perform a local
    // copy and keep the input arguments unchanged.
    // row numbers are sorted if required
    if (!already_sorted)
      Seldon::Assemble(nb_interac, col_interac, val_interac);
    
    // We add new entries which have already a slot in the matrix.
    int nb_new = 0;
    Vector<bool> new_interac(nb_interac);
    new_interac.Fill(true);
    int k = 0;
    for (int j = 0; j < nb_interac; j++)
      {
	while (k < this->m_ && index_[k] < col_interac(j))
	  k++;
	
	if (k < this->m_ && col_interac(j) == index_[k])
	  {
	    new_interac(j) = false;
	    this->data_[k] += val_interac(j);
	  }
	else
	  nb_new++;
      }
    
    if (nb_new > 0)
      {
	// Some new entries have no slot, we add them on the correct position.
	Vector<T> new_val(this->m_+nb_new);
	Vector<int> new_ind(this->m_+nb_new);
	int nb = 0, k = 0;
	for (int j = 0; j < nb_interac; j++)
	  if (new_interac(j))
	    {
	      while (k < this->m_ && index_[k] < col_interac(j))
		{
		  new_ind(nb) = index_[k];
		  new_val(nb) = this->data_[k];
		  k++;
		  nb++;
		}
	      
	      // The new entry.
	      new_ind(nb) = col_interac(j);
	      new_val(nb) = val_interac(j); nb++;
	    }
	
	// Last entries.
	while (k < this->m_)
	  {
	    new_ind(nb) = index_[k];
	    new_val(nb) = this->data_[k];
	    k++;
	    nb++;
	  }
	
	SetData(new_val, new_ind);
      }
  }
  
  
  /**************************
   * OUTPUT/INPUT FUNCTIONS *
   **************************/
  
  
  //! Writes the vector in a file.
  /*!
    The length and number of non-zero entries
    of the vector (two integers) and non-zero elements of the vector are
    stored in binary format. For non-zero elements, row numbers
    are first written, then values.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<Vect_Sparse>::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);
    
    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    The number of non-zero entries
    of the vector (two integers) and non-zero elements of the vector are
    stored in binary format. For non-zero elements, row numbers
    are first written, then values.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::Write(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Sparse>::Write(ofstream& FileStream)",
                    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));

    FileStream.write(reinterpret_cast<char*>(this->index_),
		     this->m_ * sizeof(int));

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->m_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Sparse>::Write(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The output file may have been removed")
		    + "or there is no space left on device.");
#endif

  }


  //! Writes the vector in a file.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<Vect_Sparse>::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif
    
    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Sparse>::WriteText(ofstream& FileStream)",
                    "Stream is not ready.");
#endif
    
    // 1-index
    for (int i = 0; i < this->m_ - 1; i++)
      FileStream<<(Index(i)+1)<<" "<<Value(i)<<'\n';
    
    // no empty line at the end of the file
    if (this->m_ > 0)
      FileStream<<(Index(this->m_-1)+1)<<" "<<Value(this->m_-1);
    
#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Sparse>::WriteText(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The output file may have been removed")
		    + "or there is no space left on device.");
#endif
    
  }


  //! Sets the vector from a file.
  /*!
    Sets the vector according to a binary file that stores the
    number of rows and non-zero entries of the vector (two integers).
    For non-zero elements, row numbers are first read, then values.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<Vect_Sparse>::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets the vector according to a binary file stream that stores the
    number of rows and non-zero entries of the vector (two integers).
    For non-zero elements, row numbers are first read, then values.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Sparse>::Read(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    int m;
    FileStream.read(reinterpret_cast<char*>(&m), sizeof(int));
    this->Reallocate(m);
    
    FileStream.read(reinterpret_cast<char*>(this->index_),
		    m * sizeof(int));

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    m * sizeof(value_type));
    
#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Sparse>::Read(ifstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The intput file may have been removed")
		    + " or may not contain enough data.");
#endif

  }

  
  //! Sets the vector from a file.
  /*!
    Sets all elements of the vector according to a text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<Vect_Sparse>::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream);

    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets all elements of the vector according to a text format. The length is not
    stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, Vect_Sparse, Allocator>::ReadText(istream& FileStream)
  {
    // previous vector is cleared
    Clear();
    
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<Vect_Sparse>::ReadText(ofstream& FileStream)",
                    "Stream is not ready.");
#endif
    
    Vector<T, Vect_Full, Allocator> values;
    Vector<int> row_numbers;
    T entry; int row = 0;
    int nb_elt = 0;
    while (!FileStream.eof())
      {
	// new entry is read (1-index)
	FileStream>>row>>entry;
	
	if (FileStream.fail())
	  break;
	else
	  {
#ifdef SELDON_CHECK_IO
	    if (row < 1)
	      throw IOError(string("Vector<Vect_Sparse>::ReadText")+
			    "(istream& FileStream)",
			    string("Error : Row number should be greater ")
			    + "than 0 but is equal to " + to_str(row));
#endif
	
	    nb_elt++; row--;
	    
	    // inserting new element
	    if (nb_elt > values.GetM())
	      {
		values.Resize(2*nb_elt);
		row_numbers.Resize(2*nb_elt);
	      }
	    
	    values(nb_elt-1) = entry;
	    row_numbers(nb_elt-1) = row;
	  }
      }
    
    if (nb_elt > 0)
      {
	// we allocate to the right size
	this->Reallocate(nb_elt);
	for (int i = 0; i < nb_elt; i++)
	  {
	    Index(i) = row_numbers(i);
	    Value(i) = values(i);
	  }
      }
  }
  
  
  //! operator<< overloaded for sparse vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, Vect_Sparse, Allocator>& V)
  {
    for (int i = 0; i < V.GetLength(); i++)
      out << (V.Index(i) + 1) << ' ' << V.Value(i) << '\n';
    
    return out;
  }
  
  
} // namespace Seldon.

#define SELDON_FILE_SPARSE_VECTOR_CXX
#endif

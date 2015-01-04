// Copyright (C) 2010, INRIA
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef SELDON_FILE_VECTOR_VECTORCOLLECTION_INLINE_CXX


#include "VectorCollection.hxx"


namespace Seldon
{


  //////////////////////
  // VECTORCOLLECTION //
  //////////////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Default constructor.
  /*!
    Nothing is allocated. The vector length is set to zero.
  */
  template <class T, class Allocator>
  inline Vector<T, Collection, Allocator>::Vector():
    Vector_Base<T, Allocator>(), label_map_(), label_vector_()
  {
    Nvector_ = 0;
  }


  //! Main constructor.
  /*! Builds a vector collection of a given size.
    \param[in] i length of the vector.
  */
  template <class T, class Allocator>
  inline Vector<T, Collection, Allocator>::Vector(int i):
    Vector_Base<T, Allocator>(0), length_(i), length_sum_(i), vector_(i),
    label_map_(), label_vector_()
  {
    Nvector_ = i;
    for (int k = 0; k < i; k++)
      length_(k) = length_sum_(k) = 0;
  }


  //! Copy constructor.
  /*! Builds a copy of a vector collection.
    \param[in] V vector collection to be copied.
  */
  template <class T, class Allocator>
  inline Vector<T, Collection, Allocator>::
  Vector(const Vector<T, Collection, Allocator>& V):
    Vector_Base<T, Allocator>(V), Nvector_(0)
  {
    Copy(V);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  /*! The inner vectors are nullified so that their memory blocks should not
    be deallocated.
  */
  template <class T, class Allocator >
  inline Vector<T, Collection, Allocator>::~Vector()
  {
    for (int i = 0; i < Nvector_; i++)
      vector_(i).Nullify();
    label_map_.clear();
    label_vector_.clear();
  }


  //! Clears the vector collection.
  /*! The inner vectors are nullified so that their memory blocks should not
    be deallocated.
  */
  template <class T, class Allocator >
  inline void Vector<T, Collection, Allocator>::Clear()
  {
    for (int i = 0; i < Nvector_; i++)
      vector_(i).Nullify();
    vector_.Clear();
    length_sum_.Clear();
    length_.Clear();
    Nvector_ = 0;
    this->m_ = 0;
    label_map_.clear();
    label_vector_.clear();
  }


  //! Reallocates the vector collection.
  /*! This method first clears the collection. Then it allocates a new vector
    of size \a i, and puts this vector in the collection. On exit, the
    collection is only composed of this vector of size \a i.
    \param[in] i new size.
  */
  template <class T, class Allocator >
  inline void Vector<T, Collection, Allocator>::Reallocate(int i)
  {
    Clear();
    vector_type v;
    v.Reallocate(i);
    AddVector(v);
    v.Nullify();
  }


  //! Clears the vector collection.
  /*! The inner vectors are cleared and the memory
    blocks are deallocated.
  */
  template <class T, class Allocator >
  inline void Vector<T, Collection, Allocator>::Deallocate()
  {
    for (int i = 0; i < Nvector_; i++)
      vector_(i).Clear();
    vector_.Clear();
    length_sum_.Clear();
    length_.Clear();
    Nvector_ = 0;
    this->m_ = 0;
    label_map_.clear();
    label_vector_.clear();
  }


  /*****************
   * BASIC METHODS *
   *****************/


  //! Returns the total number of elements.
  /*!
    \return The total length of the vector.
  */
  template <class T, class Allocator >
  inline int Vector<T, Collection, Allocator>::GetM() const
  {
    return this->m_;
  }


  //! Returns the total number of elements.
  /*!
    \return The total length of the vector.
  */
  template <class T, class Allocator >
  inline int Vector<T, Collection, Allocator>::GetLength() const
  {
    return this->m_;
  }


  //! Returns the number of aggregated vectors.
  /*!
    \return The total number of aggregated vectors.
  */
  template <class T, class Allocator >
  inline int Vector<T, Collection, Allocator>::GetNvector() const
  {
    return Nvector_;
  }


  //! Returns the length vector of the underlying vectors.
  /*!
    \return The lengths of the underlying vectors.
  */
  template <class T, class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<T, Collection, Allocator>::GetVectorLength() const
  {
    return length_;
  }


  //! Returns the cumulative sum of the lengths of the underlying vectors.
  /*!
    \return The cumulative sum of the lengths of the underlying vectors.
  */
  template <class T, class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<T, Collection, Allocator>::GetLengthSum() const
  {
    return length_sum_;
  }


  //! Returns the vector index of the aggregated vector named \a name.
  /*!
    \param[in] name name of the aggregated vector.
    \return The vector index of the aggregated vector.
  */
  template <class T, class Allocator >
  inline int Vector<T, Collection, Allocator>
  ::GetVectorIndex(string name) const
  {
    map<string,int>::const_iterator label_iterator = label_map_.find(name);
    if (label_iterator == label_map_.end())
      throw WrongArgument("VectorCollection::GetVectorIndex(string name)",
			  "Unknown vector name: \"" + name + "\".");
    return label_iterator->second;
  }


  /*! \brief Returns the index of the first element of the aggregated vector
    named \a name. */
  /*!
    \param[in] name name of the aggregated vector.
    \return The index of the first element of the aggregated vector.
  */
  template <class T, class Allocator >
  inline int Vector<T, Collection, Allocator>::GetIndex(string name) const
  {
    map<string,int>::const_iterator label_iterator = label_map_.find(name);
    if (label_iterator == label_map_.end())
      throw WrongArgument("VectorCollection::GetIndex(string name)",
			  string("Unknown vector name: \"") + name + "\".");
    return (label_iterator->second == 0) ?
      0 : length_sum_(label_iterator->second - 1);
  }


  //! Returns the list of vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::collection_reference
  Vector<T, Collection, Allocator>::GetVector()
  {
    return vector_;
  }


  //! Returns the list of vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::const_collection_reference
  Vector<T, Collection, Allocator>::GetVector() const
  {
    return vector_;
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::vector_reference
  Vector<T, Collection, Allocator>::GetVector(int i)
  {
    return vector_(i);
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, class Allocator >
  inline typename
  Vector<T, Collection, Allocator>::const_vector_reference
  Vector<T, Collection, Allocator>::GetVector(int i) const
  {
    return vector_(i);
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] name the name of the vector to be returned.
    \return The aggregated vector named \a name.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::vector_reference
  Vector<T, Collection, Allocator>::GetVector(string name)
  {
    map<string,int>::iterator label_iterator;
    label_iterator = label_map_.find(name);
    if (label_iterator == label_map_.end())
      throw WrongArgument("VectorCollection::SetVector(string name)",
			  string("Unknown vector name: \"") + name + "\".");
    return GetVector(label_iterator->second);
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] name the name of the vector to be returned.
    \return The aggregated vector named \a name.
  */
  template <class T, class Allocator >
  inline typename
  Vector<T, Collection, Allocator>::const_vector_reference
  Vector<T, Collection, Allocator>::GetVector(string name) const
  {
    map<string,int>::const_iterator label_iterator;
    label_iterator = label_map_.find(name);
    if (label_iterator == label_map_.end())
      throw WrongArgument("VectorCollection::SetVector(string name)",
			  string("Unknown vector name: \"") + name + "\".");
    return GetVector(label_iterator->second);
  }


  /*********************************
   * ELEMENT ACCESS AND ASSIGNMENT *
   *********************************/


  //! Access operator.
  /*!
    \param[in] i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::reference
  Vector<T, Collection, Allocator>::operator() (int i)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("VectorCollection::operator()",
		       string("Index should be in [0, ")
                       + to_str(this->m_ - 1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    int j = 0;
    while (i >= length_sum_(j))
      j++;
    return (j == 0) ? vector_(j)(i) : vector_(j)(i - length_sum_(j - 1));
  }


  //! Access operator.
  /*!
    \param[in] i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator >
  inline typename Vector<T, Collection, Allocator>::const_reference
  Vector<T, Collection, Allocator>::operator() (int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("VectorCollection::operator()",
		       string("Index should be in [0, ")
                       + to_str(this->m_ - 1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    int j = 0;
    while (i >= length_sum_(j))
      j++;
    return (j == 0) ? vector_(j)(i) : vector_(j)(i - length_sum_(j - 1));
  }


  //! Duplicates a vector collection (assignment operator).
  /*!
    \param[in] X vector collection to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator >
  inline Vector<T, Collection, Allocator>&
  Vector<T, Collection, Allocator>::operator=
  (const Vector<T, Collection, Allocator>& X)
  {
    this->Copy(X);
    return *this;
  }


  //! Duplicates a vector collection.
  /*!
    \param[in] X vector collection to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator >
  inline void Vector<T, Collection, Allocator>
  ::Copy(const Vector<T, Collection, Allocator>& X)
  {
    Clear();
    for (int i = 0; i < X.GetNvector(); i++)
      AddVector(X.GetVector(i));
    label_map_.insert(X.label_map_.begin(), X.label_map_.end());
    label_vector_.assign(X.label_vector_.begin(), X.label_vector_.end());
  }


  //! Copies the values of a full vector into the current vector.
  /*!
    \param[in] X full vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator >
  template <class T0, class Allocator0>
  inline void Vector<T, Collection, Allocator>
  ::Copy(const Vector<T0, VectFull, Allocator0>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (this->m_ != X.GetM())
      throw WrongIndex("VectorCollection::Copy(X)",
		       string("The size of X should be equal to ")
                       + to_str(this->m_ - 1)
		       + ", but is equal to " + to_str(X.GetM()) + ".");
#endif

    for (int i = 0; i < X.GetM(); i++)
      (*this)(i) = X(i);
  }


  //! Multiplies a vector collection by a scalar.
  /*!
    \param[in] alpha scalar.
  */
  template <class T, class Allocator>
  template<class T0>
  inline Vector<T, Collection, Allocator>&
  Vector<T, Collection, Allocator>::operator*= (const T0& alpha)
  {
    for (int i = 0; i < this->Nvector_; i++)
      this->vector_(i) *= alpha;

    return *this;
  }
  
  
} // namespace Seldon.


#define SELDON_FILE_VECTOR_VECTORCOLLECTION_INLINE_CXX
#endif

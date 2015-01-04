// Copyright (C) 2001-2009 INRIA
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


#ifndef SELDON_FILE_VECTOR_HETEROGENEOUSCOLLECTION_INLINE_CXX


#include "HeterogeneousCollection.hxx"


namespace Seldon
{


  ////////////////////////////////////
  // VECTOR HETEROGENEOUSCOLLECTION //
  ////////////////////////////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Default constructor.
  /*!
    Nothing is allocated. The number of vectors is set to zero.
  */
  template <class T, template <class U> class Allocator >
  inline Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::Vector(): Vector_Base<T, Allocator<T> >(), label_map_(),
	      label_vector_()
  {
    Nvector_ = 0;
  }


  //! Copy constructor.
  /*! Builds a copy of a vector collection.
    \param[in] V vector collection to be copied.
  */
  template <class T, template <class U> class Allocator >
  inline Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::Vector(const
	   Vector<FloatDouble, DenseSparseCollection, Allocator<T> >& V):
    Vector_Base<T, Allocator<T> >(V), label_map_(), label_vector_()
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
  template <class T, template <class U> class Allocator >
  inline Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::~Vector()
  {
    float_dense_c_.Clear();
    float_sparse_c_.Clear();
    double_dense_c_.Clear();
    double_sparse_c_.Clear();
    collection_.Clear();
    subvector_.Clear();
    length_sum_.Clear();
    length_.Clear();
    Nvector_ = 0;
    this->m_ = 0;
    label_map_.clear();
    label_vector_.clear();
  }


  //! Clears the vector collection.
  /*! The inner vectors are nullified so that their memory blocks should not
    be deallocated.
  */
  template <class T, template <class U> class Allocator >
  inline void Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::Clear()
  {
    float_dense_c_.Clear();
    float_sparse_c_.Clear();
    double_dense_c_.Clear();
    double_sparse_c_.Clear();
    collection_.Clear();
    subvector_.Clear();
    length_sum_.Clear();
    length_.Clear();
    Nvector_ = 0;
    this->m_ = 0;
    label_map_.clear();
    label_vector_.clear();
  }


  //! Clears the vector collection.
  /*! The inner vectors are cleared and the memory
    blocks are deallocated.
  */
  template <class T, template <class U> class Allocator >
  inline void Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::Deallocate()
  {
    float_dense_c_.Deallocate();
    float_sparse_c_.Deallocate();
    double_dense_c_.Deallocate();
    double_sparse_c_.Deallocate();
    collection_.Clear();
    subvector_.Clear();
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
  template <class T, template <class U> class Allocator >
  inline int Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetM()
    const
  {
    return this->m_;
  }


  //! Returns the total number of elements.
  /*!
    \return The total length of the vector.
  */
  template <class T, template <class U> class Allocator >
  inline int Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetLength() const
  {
    return this->m_;
  }


  //! Returns the number of aggregated vectors.
  /*!
    \return The total number of aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline int Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetNvector() const
  {
    return Nvector_;
  }


  //! Returns the length vector of the underlying vectors.
  /*!
    \return The lengths of the underlying vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetVectorLength() const
  {
    return length_;
  }


  //! Returns the cumulative sum of the lengths of the underlying vectors.
  /*!
    \return The cumulative sum of the lengths of the underlying vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetLengthSum() const
  {
    return length_sum_;
  }


  //! Returns the collection indexes of the underlying vectors.
  /*!
    \return The collection index vector of the underlying vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetCollectionIndex() const
  {
    return collection_;
  }


  //! Returns the index of the underlying vectors in the inner collection.
  /*!
    \return The index vector of the underlying vectors in the inner
    collection.
  */
  template <class T, template <class U> class Allocator >
  inline const Vector<int, VectFull, MallocAlloc<int> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetSubvectorIndex() const
  {
    return length_sum_;
  }


  //! Returns the list of float dense vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::float_dense_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetFloatDense()
  {
    return float_dense_c_;
  }


  //! Returns the list of float dense vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const
  typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::float_dense_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetFloatDense() const
  {
    return float_dense_c_;
  }


  //! Returns the list of float sparse vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline typename
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::float_sparse_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetFloatSparse()
  {
    return float_sparse_c_;
  }


  //! Returns the list of float sparse vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline const
  typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::float_sparse_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetFloatSparse() const
  {
    return float_sparse_c_;
  }


  //! Returns the list of double dense vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::double_dense_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetDoubleDense()
  {
    return double_dense_c_;
  }


  //! Returns the list of double dense vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline
  const typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::double_dense_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetDoubleDense() const
  {
    return double_dense_c_;
  }


  //! Returns the list of double sparse vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::double_sparse_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::GetDoubleSparse()
  {
    return double_sparse_c_;
  }


  //! Returns the list of double sparse vectors.
  /*!
    \return The list of the aggregated vectors.
  */
  template <class T, template <class U> class Allocator >
  inline
  const typename Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::double_sparse_c&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetDoubleSparse() const
  {
    return double_sparse_c_;
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, template <class U> class Allocator >
  inline void Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetVector(int i, typename
	      Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
	      ::float_dense_v& vector) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Nvector_)
      throw WrongIndex("Vector<FloatDouble, DenseSparseCollection>"
                       "::GetVector(int i, Vector<float, VectFull>&)",
		       string("Index should be in [0, ")
                       + to_str(Nvector_ - 1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    if (collection_(i) != 0)
      WrongArgument("Vector<FloatDouble, DenseSparseCollection>"
                    "::GetVector(int i, Vector<float, VectFull>&)",
                    string("The ") + to_str(i) + "-th inner vector "
                    "is of type " + GetType(i) + ".");

    vector.SetData(float_dense_c_.GetVector(subvector_(i)));
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, template <class U> class Allocator >
  inline void Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetVector(int i, typename
	      Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
	      ::float_sparse_v& vector) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Nvector_)
      throw WrongIndex("Vector<FloatDouble, DenseSparseCollection>"
                       "::GetVector(int i, Vector<float, VectSparse>&)",
		       string("Index should be in [0, ")
                       + to_str(Nvector_ - 1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    if (collection_(i) != 1)
      WrongArgument("Vector<FloatDouble, DenseSparseCollection>"
                    "::GetVector(int i, Vector<float, VectSparse>&)",
                    string("The ") + to_str(i) + "-th inner vector "
                    "is of type " + GetType(i) + ".");

    vector.SetData(float_sparse_c_.GetVector(subvector_(i)));
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, template <class U> class Allocator >
  inline void Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetVector(int i, typename
	      Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
	      ::double_dense_v& vector) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Nvector_)
      throw WrongIndex("Vector<FloatDouble, DenseSparseCollection>"
                       "::GetVector(int i, Vector<double, VectDense>&)",
		       string("Index should be in [0, ")
                       + to_str(Nvector_ - 1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    if (collection_(i) != 2)
      WrongArgument("Vector<FloatDouble, DenseSparseCollection>"
                    "::GetVector(int i, Vector<double, VectDense>&)",
                    string("The ") + to_str(i) + "-th inner vector "
                    "is of type " + GetType(i) + ".");

    vector.SetData(double_dense_c_.GetVector(subvector_(i)));
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, template <class U> class Allocator >
  inline void Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetVector(int i, typename
	      Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
	      ::double_sparse_v& vector) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Nvector_)
      throw WrongIndex("Vector<FloatDouble, DenseSparseCollection>"
                       "::GetVector(int i, Vector<double, VectSparse>&)",
		       string("Index should be in [0, ")
                       + to_str(Nvector_ - 1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    if (collection_(i) != 3)
      WrongArgument("Vector<FloatDouble, DenseSparseCollection>"
                    "::GetVector(int i, Vector<double, VectSparse>&)",
                    string("The ") + to_str(i) + "-th inner vector "
                    "is of type " + GetType(i) + ".");

    vector.SetData(double_sparse_c_.GetVector(subvector_(i)));
  }


  //! Returns one of the aggregated vectors.
  /*!
    \param[in] i the index of the vector to be returned.
    \return The \a i th aggregated vector.
  */
  template <class T, template <class U> class Allocator >
  template <class T0, class Storage0, class Allocator0>
  inline void Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::GetVector(string name,Vector<T0, Storage0, Allocator0>& vector) const
  {
    map<string,int>::const_iterator label_iterator;
    label_iterator = label_map_.find(name);
    if (label_iterator == label_map_.end())
      throw WrongArgument("Vector<FloatDouble, DenseSparseCollection>"
                          "::SetVector(string name)",
			  string("Unknown vector name ") + name + ".");
    GetVector(label_iterator->second, vector);
  }


  /*********************************
   * ELEMENT ACCESS AND ASSIGNMENT *
   *********************************/


  //! Access operator.
  /*!
    \param[in] i index.
    \return The value of the vector at 'i'.
  */
  template <class T, template <class U> class Allocator >
  inline double
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::operator() (int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<FloatDouble, DenseSparse>::operator()",
		       string("Index should be in [0, ")
                       + to_str(this->m_ - 1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    int j = 0;
    while (i >= length_sum_(j))
      j++;

    switch (collection_(j))
      {
      case 0:
	return (j == 0) ? double(float_dense_c_.GetVector(subvector_(j))(i)) :
	  double(float_dense_c_.
		 GetVector(subvector_(j))(i - length_sum_(j - 1)));
      case 1:
	return (j == 0) ? double(float_sparse_c_.GetVector(subvector_(j))(i)):
	  double(float_sparse_c_.
		 GetVector(subvector_(j))(i - length_sum_(j - 1)));
      case 2:
	return j == 0 ? double_dense_c_.GetVector(subvector_(j))(i) :
	  double_dense_c_.GetVector(subvector_(j))(i - length_sum_(j - 1));
      case 3:
	return j == 0 ? double_sparse_c_.GetVector(subvector_(j))(i) :
	  double_sparse_c_.GetVector(subvector_(j))(i - length_sum_(j - 1));
      default:
	return 0.;
      }
  }


  //! Duplicates a vector collection (assignment operator).
  /*!
    \param[in] X vector collection to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, template <class U> class Allocator >
  inline Vector<FloatDouble, DenseSparseCollection, Allocator<T> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >::operator=
  (const Vector<FloatDouble, DenseSparseCollection, Allocator<T> >& X)
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
  template <class T, template <class U> class Allocator >
  inline void Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::Copy(const Vector<FloatDouble, DenseSparseCollection, Allocator<T> >& X)
  {
    Clear();
    this->m_ = X.GetM();
    collection_.Copy(X.collection_);
    subvector_.Copy(X.subvector_);
    length_.Copy(X.length_);
    length_sum_.Copy(X.length_sum_);

    float_dense_c_.Copy(X.float_dense_c_);
    float_sparse_c_.Copy(X.float_sparse_c_);
    double_dense_c_.Copy(X.double_dense_c_);
    double_sparse_c_.Copy(X.double_sparse_c_);

    label_map_.insert(X.label_map_.begin(), X.label_map_.end());
    label_vector_.assign(X.label_vector_.begin(), X.label_vector_.end());
  }


  //! Multiplies a vector collection by a scalar.
  /*!
    \param[in] alpha scalar.
  */
  template <class T, template <class U> class Allocator >
  template<class T0>
  inline Vector<FloatDouble, DenseSparseCollection, Allocator<T> >&
  Vector<FloatDouble, DenseSparseCollection, Allocator<T> >
  ::operator*= (const T0& alpha)
  {
    float_dense_c_ *= alpha;
    float_sparse_c_ *= alpha;
    double_dense_c_ *= alpha;
    double_sparse_c_ *= alpha;
    return *this;
  }
  
  
} // namespace Seldon.


#define SELDON_FILE_VECTOR_HETEROGENEOUSCOLLECTION_INLINE_CXX
#endif

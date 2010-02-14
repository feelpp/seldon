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


#ifndef SELDON_FILE_VECTOR_VECTOR2_CXX


#include "Vector2.hxx"


namespace Seldon
{


  /////////////
  // VECTOR2 //
  /////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Default constructor.
  /*!
    Nothing is allocated.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector2<T, Allocator0, Allocator1>::Vector2()
  {
  }


  //! Constructor.
  /*! The vector of vectors is allocated with \a length empty vectors.
    \param[in] length the length of the vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector2<T, Allocator0, Allocator1>::Vector2(int length)
  {
    data_.Reallocate(length);
  }


  //! Constructor.
  /*! The vector of vectors and the inner vectors are allocated.
    \param[in] length the lengths of the inner vectors. The vector of vectors
    will obviously have as many elements as \a length has.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector2<T, Allocator0, Allocator1>::Vector2(const Vector<int>& length)
  {
    data_.Clear();
    int m = length.GetSize();
    data_.Reallocate(m);
    for(int i = 0; i < m; i++)
      data_(i).Reallocate(length(i));
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor. The vector of vectors and the inner vectors are deallocated.
  template <class T, class Allocator0, class Allocator1>
  Vector2<T, Allocator0, Allocator1>::~Vector2()
  {
  }


  /*****************************
   * MANAGEMENT OF THE VECTORS *
   *****************************/


  //! Returns the size along dimension 1.
  /*!
    \return The size along dimension 1.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetSize() const
  {
    return data_.GetSize();
  }


  //! Returns the size along dimension 1.
  /*!
    \return The size along dimension 1.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetLength() const
  {
    return data_.GetLength();
  }


  //! Returns the size of the inner vector #\a i.
  /*!
    \param[in] i index.
    \return The size of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetSize(int i) const
  {
    return data_(i).GetSize();
  }


  //! Returns the size of the inner vector #\a i.
  /*!
    \param[in] i index.
    \return The size of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetLength(int i) const
  {
    return data_(i).GetLength();
  }


  //! Reallocates the vector of vector.
  /*!
    \param[in] M the new size of the vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Reallocate(int M)
  {
    data_.Reallocate(M);
  }


  //! Reallocates the inner vector #\a i.
  /*!
    \param[in] i index of the inner vector to be reallocated.
    \param[in] N the new size of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Reallocate(int i, int N)
  {
    data_(i).Reallocate(N);
  }


  //! Reallocates the whole structure.
  /*!
    \param[in] length the new lengths of the inner vectors. The vector of
    vectors will obviously have as many elements as \a length has.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::Reallocate(const Vector<int>& length)
  {
    int m = length.GetSize();
    data_.Reallocate(m);
    for(int i = 0; i < m; i++)
      data_(i).Reallocate(length(i));
  }


  //! Appends an element at the end of the inner vector #\a i.
  /*!
    \param[in] i index of the inner vector to which \a x should be appended.
    \param[in] x element to be appended.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::PushBack(int i, const T& x)
  {
    data_(i).PushBack(x);
  }


  //! Appends an inner vector at the end of the vector.
  /*!
    \param[in] X vector to be appended.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::PushBack(const Vector<T, VectFull, Allocator0>& X)
  {
    data_.PushBack(X);
  }


  //! Appends a vector of vectors.
  /*! The inner vectors of \a V are appended to the current instance, in the
    same order as they appear in \a V.
    \param[in] V vector of vectors to be appended.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::PushBack(const Vector<Vector<T, VectFull, Allocator0>,
                          VectFull, Allocator1>& V)
  {
    for (int i = 0; i < V.GetLength(); i++)
      data_.PushBack(V(i));
  }


  //! Appends a vector of vectors.
  /*! The inner vectors of \a V are appended to the current instance, in the
    same order as they appear in \a V.
    \param[in] V vector of vectors to be appended.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::PushBack(const Vector2<T, Allocator0, Allocator1>& V)
  {
    PushBack(V.GetVector());
  }


  //! Clears the vector.
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Clear()
  {
    data_.Clear();
  }


  //! Clears a given vector.
  /*!
    \param[in] i index of the vector to be cleared.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Clear(int i)
  {
    data_(i).Clear();
  }


  //! Fills the vector with a given value.
  /*!
    \param[in] x value to fill the vector with.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Fill(const T& x)
  {
    for (int i = 0; i < data_.GetLength(); i++)
      data_(i).Fill(x);
  }


  //! Returns the vector of vectors.
  /*!
    \return The vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector<Vector<T, VectFull, Allocator0>, VectFull, Allocator1>&
  Vector2<T, Allocator0, Allocator1>::GetVector()
  {
    return data_;
  }


  //! Returns the vector of vectors.
  /*!
    \return The vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  const Vector<Vector<T, VectFull, Allocator0>, VectFull, Allocator1>
  Vector2<T, Allocator0, Allocator1>::GetVector() const
  {
    return data_;
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \return The inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector<T, VectFull, Allocator0>&
  Vector2<T, Allocator0, Allocator1>::GetVector(int i)
  {
    return data_(i);
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \return The inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  const Vector<T, VectFull, Allocator0>&
  Vector2<T, Allocator0, Allocator1>::GetVector(int i) const
  {
    return data_(i);
  }


  /*********************************
   * ELEMENT ACCESS AND ASSIGNMENT *
   *********************************/


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \return The inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  const Vector<T, VectFull, Allocator0>&
  Vector2<T, Allocator0, Allocator1>::operator() (int i) const
  {
    return data_(i);
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \return The inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector<T, VectFull, Allocator0>&
  Vector2<T, Allocator0, Allocator1>::operator() (int i)
  {
    return data_(i);
  }


  //! Returns an element of a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \param[in] j index of the element in the inner vector #\a i.
    \return The element #\a j of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  typename Vector2<T, Allocator0, Allocator1>::const_reference
  Vector2<T, Allocator0, Allocator1>::operator() (int i, int j) const
  {
    return data_(i)(j);
  }


  //! Returns an element of a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \param[in] j index of the element in the inner vector #\a i.
    \return The element #\a j of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  typename Vector2<T, Allocator0, Allocator1>::reference
  Vector2<T, Allocator0, Allocator1>::operator() (int i, int j)
  {
    return data_(i)(j);
  }


  /**********************
   * CONVENIENT METHODS *
   **********************/


  //! Checks whether another Vector2 instance has the same structure.
  /*! Checks whether another Vector2 instance has the same structure as the
    current instance. The structures are the same if both instances have the
    same number of inner vectors, and if the inner vectors have the same
    lengths.
    \param[in] V Vector2 instance whose structure is compared to that of the
    current instance.
    \return True if the current instance as the same structure as \a V, false
    otherwise.
  */
  template <class T, class Allocator0, class Allocator1>
  template <class V2>
  bool Vector2<T, Allocator0, Allocator1>::HasSameStructure(const V2& V) const
  {
    if (V.GetLength() != GetLength())
      return false;
    for (int i = 0; i < GetLength(); i++)
      if (V.GetLength(i) != GetLength(i))
	return false;
    return true;
  }


  //! Displays the vector.
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Print() const
  {
    for (int i = 0; i < data_.GetLength(); i++)
      {
        cout << "Vector " << i << ": ";
        data_(i).Print();
      }
  }


} // namespace Seldon.


#define SELDON_FILE_VECTOR_VECTOR2_CXX
#endif

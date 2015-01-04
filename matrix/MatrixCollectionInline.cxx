// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu
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


#ifndef SELDON_FILE_MATRIX_COLLECTION_INLINE_CXX

#include "MatrixCollection.hxx"


namespace Seldon
{


  //////////////////////
  // MATRIXCOLLECTION //
  //////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline MatrixCollection<T, Prop, Storage, Allocator>::MatrixCollection():
    Matrix_Base<T, Allocator>(), Mlocal_(), Mlocal_sum_(1),
    Nlocal_(), Nlocal_sum_(1), matrix_()
  {
    nz_ = 0;
    Mmatrix_ = 0;
    Nmatrix_ = 0;
    Mlocal_sum_.Fill(0);
    Nlocal_sum_.Fill(0);
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline MatrixCollection<T, Prop, Storage, Allocator>
  ::MatrixCollection(int i, int j): Matrix_Base<T, Allocator>(i, j),
				    Mlocal_(i), Mlocal_sum_(i + 1),
				    Nlocal_(j), Nlocal_sum_(j + 1),
				    matrix_(i, j)
  {
    nz_ = 0;
    Mmatrix_ = i;
    Nmatrix_ = j;
    Mlocal_.Fill(0);
    Nlocal_.Fill(0);
    Mlocal_sum_.Fill(0);
    Nlocal_sum_.Fill(0);
  }


  //! Copy constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline MatrixCollection<T, Prop, Storage, Allocator>
  ::MatrixCollection(const MatrixCollection<T, Prop, Storage, Allocator>& A)
    : Matrix_Base<T, Allocator>()
  {
    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline MatrixCollection<T, Prop, Storage, Allocator>::~MatrixCollection()
  {
    Clear();
  }


  //! Clears the matrix.
  /*!
    Destructs the matrix.
    \warning On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void MatrixCollection<T, Prop, Storage, Allocator>::Clear()
  {
    for (int i = 0; i < Mmatrix_; i++)
      for (int j = 0; j < Nmatrix_; j++)
	matrix_(i, j).Nullify();

    matrix_.Clear();

    nz_ = 0;
    Mmatrix_ = 0;
    Nmatrix_ = 0;
    Mlocal_.Clear();
    Nlocal_.Clear();
    Mlocal_sum_.Clear();
    Nlocal_sum_.Clear();
  }


  //! Clears the matrix.
  /*!
    Destructs the matrix.
    \warning On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void MatrixCollection<T, Prop, Storage, Allocator>::Nullify()
  {
    for (int i = 0; i < Mmatrix_; i++)
      for (int j = 0; j < Nmatrix_; j++)
	matrix_(i, j).Nullify();

    nz_ = 0;
    Mmatrix_ = 0;
    Nmatrix_ = 0;
    Mlocal_.Clear();
    Nlocal_.Clear();
    Mlocal_sum_.Clear();
    Nlocal_sum_.Clear();
  }


  //! Clears a given underlying matrix.
  /*!
    \param[in] i row of the underlying matrix to be nullified.
    \param[in] j column of the underlying matrix to be nullified.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void MatrixCollection<T, Prop, Storage, Allocator>
  ::Nullify(int i, int j)
  {
    nz_ -= matrix_(i, j).GetDataSize();
    matrix_(i, j).Nullify();
  }


  //! Deallocates underlying the matrices.
  template <class T, class Prop, class Storage, class Allocator>
  inline void MatrixCollection<T, Prop, Storage, Allocator>::Deallocate()
  {
    for (int i = 0; i < Mmatrix_; i++)
      for (int j = 0; j < Nmatrix_; j++)
	GetMatrix(i, j).Clear();
    this->~MatrixCollection();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of rows.
  /*!
    \return the total number of rows. It is the sum of the number of rows in
    the underlying matrices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetM() const
  {
    return this->m_;
  }


  //! Returns the number of rows.
  /*!
    \return the total number of rows. It is the sum of the number of rows
    in the underlying matrices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetMmatrix() const
  {
    return Mmatrix_;
  }


  //! Returns the number of rows in an underlying matrix.
  /*!
    \param[in] i row index of the underlying matrix.
    \return The number of rows in the underlying matrices with row index \a i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetM(int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::GetM()",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
		     + to_str(i) + ".");
#endif

    return Mlocal_(i);
  }


  //! Returns the number of columns.
  /*!
    \return the total number of columns. It is the sum of the number of
    columns in the underlying matrices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetN() const
  {
    return this->n_;
  }


  //! Returns the number of columns.
  /*!
    \return the total number of columns. It is the sum of the number of
    columns in the underlying matrices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetNmatrix() const
  {
    return Nmatrix_;
  }


  //! Returns the number of columns in an underlying matrix.
  /*!
    \param[in] j column index of the underlying matrix.
    \return The number of columns in the underlying matrices with column index
    \a j.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetN(int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::GetN()",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    return Nlocal_(j);
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetSize() const
  {
    return this->m_ * this->n_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int MatrixCollection<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return nz_;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix collection.
  /*! On exit, the matrix is a matrix collection with \a i x \a j underlying
    matrices.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void MatrixCollection<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j)
  {
    nz_ = 0;
    Mmatrix_ = i;
    Nmatrix_ = j;
    Mlocal_.Reallocate(i);
    Nlocal_.Reallocate(j);
    Mlocal_sum_.Reallocate(i + 1);
    Nlocal_sum_.Reallocate(j + 1);
    Mlocal_.Fill(0.);
    Nlocal_.Fill(0.);
    Mlocal_sum_.Fill(0.);
    Nlocal_sum_.Fill(0.);
    matrix_.Reallocate(i, j);
  }


  //! Sets an underlying  matrix in the matrix collection.
  /*!
    \param[in] i row of the underlying matrix to be set.
    \param[in] j column of the underlying matrix to be set.
    \param[in] matrix new value for the underlying matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0, class Prop0, class Storage0, class Allocator0>
  inline void MatrixCollection<T, Prop, Storage, Allocator>
  ::SetMatrix(int i, int j, const Matrix<T0, Prop0, Storage0, Allocator0>& A)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::SetMatrix()",
                     string("Line index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::SetMatrix()",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
    if ((Mlocal_(i) != 0) && (Mlocal_(i) != A.GetM()))
      throw WrongDim("MatrixCollection::SetMatrix()",
		     string("The matrix expected should have ")
		     + to_str(this->Mlocal_(i)) + " rows, but has "
		     + to_str(A.GetM()) + " rows.");
    if ((Nlocal_(j) != 0) && (Nlocal_(j) != A.GetN()))
      throw WrongDim("MatrixCollection::SetMatrix()",
		     string("The matrix expected should have ")
		     + to_str(this->Nlocal_(j)) + " columns, but has "
		     + to_str(A.GetN()) + " columns.");
#endif

    nz_ = A.GetDataSize() - matrix_(i, j).GetDataSize();

    int Mdiff = A.GetM() - Mlocal_(i);
    int Ndiff = A.GetN() - Nlocal_(j);

    Mlocal_(i) = A.GetM();
    Nlocal_(j) = A.GetN();

    for (int k = i + 1; k < Mmatrix_ + 1; k++)
      Mlocal_sum_(k) += Mdiff;

    for (int k = j + 1; k < Nmatrix_ + 1; k++)
      Nlocal_sum_(k) += Ndiff;

    this->m_ = Mlocal_sum_(Mmatrix_);
    this->n_ = Nlocal_sum_(Nmatrix_);

    matrix_(i, j).SetData(A.GetM(), A.GetN(), A.GetData());
  }


  //! Sets an underlying matrix in the matrix collection.
  /*!
    \param[in] i row of the underlying matrix to be set.
    \param[in] j column of the underlying matrix to be set.
    \param[in] matrix new value for the underlying matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0, class Prop0, class Allocator0>
  inline void MatrixCollection<T, Prop, Storage, Allocator>
  ::SetMatrix(int i, int j, const Matrix<T0, Prop0, RowSparse, Allocator0>& A)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::SetMatrix()",
                     string("Index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::SetMatrix()",
                     string("Index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
    if ((Mlocal_(i) != 0) && (Mlocal_(i) != A.GetM()))
      throw WrongDim("MatrixCollection::SetMatrix()",
		     string("The matrix expected should have ")
		     + to_str(this->Mlocal_(i)) + " rows, but has "
		     + to_str(A.GetM()) + " rows.");
    if ((Nlocal_(j) != 0) && (Nlocal_(j) != A.GetN()))
      throw WrongDim("MatrixCollection::SetMatrix()",
		     string("The matrix expected should have ")
		     + to_str(this->Nlocal_(j)) + " columns, but has "
		     + to_str(A.GetN()) + " columns.");
#endif

    nz_ = A.GetDataSize() - matrix_(i, j).GetDataSize();

    int Mdiff = A.GetM() - Mlocal_(i);
    int Ndiff = A.GetN() - Nlocal_(j);

    Mlocal_(i) = A.GetM();
    Nlocal_(j) = A.GetN();

    for (int k = i + 1; k < Mmatrix_ + 1; k++)
      Mlocal_sum_(k) += Mdiff;

    for (int k = j + 1; k < Nmatrix_ + 1; k++)
      Nlocal_sum_(k) += Ndiff;

    this->m_ = Mlocal_sum_(Mmatrix_);
    this->n_ = Nlocal_sum_(Nmatrix_);

    matrix_(i, j).SetData(A.GetM(), A.GetN(), A.GetNonZeros(), A.GetData(),
                          A.GetPtr(), A.GetInd());
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access to an underlying matrix.
  /*!
    Returns the underlying matrix (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return The matrix collection (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline
  typename MatrixCollection<T, Prop, Storage, Allocator>::matrix_reference
  MatrixCollection<T, Prop, Storage, Allocator>::GetMatrix(int i, int j)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::GetMatrix(int, int)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::GetMatrix(int, int)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    return matrix_(i, j);
  }


  //! Access to an underlying matrix.
  /*!
    Returns the underlying matrix (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return The matrix collection (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename MatrixCollection<T, Prop, Storage, Allocator>
  ::const_matrix_reference
  MatrixCollection<T, Prop, Storage, Allocator>::GetMatrix(int i,
							   int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::GetMatrix(int, int)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::GetMatrix(int, int)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    return matrix_(i, j);
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename MatrixCollection<T, Prop, Storage, Allocator>
  ::value_type
  MatrixCollection<T, Prop, Storage, Allocator>::operator() (int i,
							     int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->Mlocal_sum_(Mmatrix_))
      throw WrongRow("MatrixCollection::operator(int, int)",
                     string("Row index should be in [0, ")
                     + to_str(this->Mlocal_sum_(Mmatrix_) - 1)
                     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->Nlocal_sum_(Nmatrix_))
      throw WrongCol("MatrixCollection::operator(int, int)",
                     string("Column index should be in [0, ")
                     + to_str(this->Nlocal_sum_(Nmatrix_) - 1)
                     + "], but is equal to " + to_str(j) + ".");
#endif

    int i_global = 0;
    while (i >= Mlocal_sum_(i_global))
      i_global++;
    i_global--;

    int j_global = 0;
    while (j >= Nlocal_sum_(j_global))
      j_global++;
    j_global--;

    return matrix_(i_global, j_global)(i - Mlocal_sum_(i_global),
                                       j - Nlocal_sum_(j_global));
  }


  //! Duplicates a matrix collection (assignment operator).
  /*!
    \param[in] A matrix collection to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline MatrixCollection<T, Prop, Storage, Allocator>&
  MatrixCollection<T, Prop, Storage, Allocator>
  ::operator= (const MatrixCollection<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Duplicates a matrix collection.
  /*!
    \param[in] A matrix collection to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void MatrixCollection<T, Prop, Storage, Allocator>
  ::Copy(const MatrixCollection<T, Prop, Storage, Allocator>& A)
  {
    Clear();

    this->nz_ = A.nz_;
    this->m_ = A.GetM();
    this->n_ = A.GetN();
    Mmatrix_ = A.Mmatrix_;
    Nmatrix_ = A.Nmatrix_;
    this->Mlocal_ = A.Mlocal_;
    this->Mlocal_sum_ = A.Mlocal_sum_;
    this->Nlocal_ = A.Nlocal_;
    this->Nlocal_sum_ = A.Nlocal_sum_;

    matrix_.Reallocate(Mmatrix_, Nmatrix_);

    for (int i = 0; i < Mmatrix_; i++)
      for (int j = 0; j < Nmatrix_; j++)
	SetMatrix(i, j, A.GetMatrix(i, j));
  }
  
  
  ////////////////////////
  // COLMAJORCOLLECTION //
  ////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajorCollection, Allocator>::Matrix():
    MatrixCollection<T, Prop, ColMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajorCollection, Allocator>::Matrix(int i, int j):
    MatrixCollection<T, Prop, ColMajor, Allocator>(i, j)
  {
  }


  ////////////////////////
  // ROWMAJORCOLLECTION //
  ////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajorCollection, Allocator>::Matrix():
    MatrixCollection<T, Prop, RowMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajorCollection, Allocator>::Matrix(int i, int j):
    MatrixCollection<T, Prop, RowMajor, Allocator>(i, j)
  {
  }


  ////////////////////////////
  // COLSYMPACKEDCOLLECTION //
  ////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymPackedCollection, Allocator>::Matrix():
    MatrixCollection<T, Prop, ColSymPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColSymPackedCollection, Allocator>
  ::Matrix(int i, int j):
    MatrixCollection<T, Prop, ColSymPacked, Allocator>(i, j)
  {
  }


  ////////////////////////////
  // ROWSYMPACKEDCOLLECTION //
  ////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymPackedCollection, Allocator>::Matrix():
    MatrixCollection<T, Prop, RowSymPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j collection matrix.
    \param[in] i number of rows of matrices.
    \param[in] j number of columns of matrices.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowSymPackedCollection, Allocator>
  ::Matrix(int i, int j):
    MatrixCollection<T, Prop, RowSymPacked, Allocator>(i, j)
  {
  }


} // namespace Seldon.


#define SELDON_FILE_MATRIX_COLLECTION_INLINE_CXX
#endif

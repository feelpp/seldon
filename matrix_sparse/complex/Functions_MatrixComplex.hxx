// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_COMPLEX_HXX

namespace Seldon
{
    
  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric,
           ArrayRowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymComplexSparse, Allocator2>& B);
  
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, Symmetric,
	   ArrayRowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B);

  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, General,
	   ArrayRowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowSparse, Allocator2>& B);

  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, General,
	   ArrayRowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowComplexSparse, Allocator2>& B);

  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2,
           class T3, class Prop3, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& A,
	   const Matrix<T2, Prop2, ArrayRowSymSparse, Allocator2>& B,
	   Matrix<T3, Prop3, ArrayRowSymComplexSparse, Allocator3>& C);

  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2,
           class T3, class Prop3, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
	   const Matrix<T2, Prop2, ArrayRowSparse, Allocator2>& B,
	   Matrix<T3, Prop3, ArrayRowComplexSparse, Allocator3>& C);

  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, RowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, RowSymComplexSparse, Allocator2>& B);

  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, General, RowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, RowComplexSparse, Allocator2>& B);

  template<class T0, class T1, class Allocator1, class T2, class Allocator2>
  void Add(const complex<T0>& alpha,
	   const Matrix<T1, Symmetric, RowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, RowSymComplexSparse, Allocator2>& B);

  template<class T0, class T1, class Allocator1,
           class T2, class Allocator2>
  void Add(const complex<T0>& alpha,
	   const Matrix<T1, General, RowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, RowComplexSparse, Allocator2>& B);
  
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A);

  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const complex<T0>& alpha,
	   Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A);

  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const T0& alpha,
           Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A);
  
  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const complex<T0>& alpha,
           Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A);

  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, Prop, RowComplexSparse, Allocator>& A);

  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const complex<T0>& alpha,
	   Matrix<T, Prop, RowComplexSparse, Allocator>& A);

  template<class T0, class T, class Prop, class Allocator>
  void Mlt(const T0& alpha, Matrix<T, Prop, RowSymComplexSparse, Allocator>& A);

  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop,
                               ArrayRowSymComplexSparse, Allocator>& A,
                               const IVect& row_perm, const IVect& col_perm);

  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);
  
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);

  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowSymComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);

  template<class Prop, class T1, class Allocator1,
	   class T2, class Allocator2, class T3, class Allocator3>
  void ScaleMatrix(Matrix<T1, Prop, ArrayRowComplexSparse, Allocator1>& A,
		   const Vector<T2, VectFull, Allocator2>& scale_left,
		   const Vector<T3, VectFull, Allocator3>& scale_right);

  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleLeftMatrix(Matrix<T1, Prop, ArrayRowComplexSparse, Allocator1>& A,
		       const Vector<T2, VectFull, Allocator2>& scale);

  template<class T1, class Allocator1,
	   class Prop, class T2, class Allocator2>
  void ScaleRightMatrix(Matrix<T1, Prop, ArrayRowComplexSparse, Allocator1>& A,
			const Vector<T2, VectFull, Allocator2>& scale);

  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, ArrayRowComplexSparse, Allocator>& A,
                 Matrix<T, General, ArrayRowComplexSparse, Allocator>& B);

  template<class T, class Allocator>
  void Transpose(const Matrix<T, General, RowComplexSparse, Allocator>& A,
                 Matrix<T, General, RowComplexSparse, Allocator>& B);

  template <class T, class Prop, class Allocator>
  T MaxAbs(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T MaxAbs(const Matrix<T, Prop, RowComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T Norm1(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T Norm1(const Matrix<T, Prop, RowComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T NormInf(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T NormInf(const Matrix<T, Prop, RowComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T MaxAbs(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T MaxAbs(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T Norm1(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T Norm1(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T NormInf(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A);

  template <class T, class Prop, class Allocator>
  T NormInf(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A);

  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A);

  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A);

  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, RowComplexSparse, Allocator>& A);

  template<class T, class Prop, class Allocator>
  void Conjugate(Matrix<T, Prop, RowSymComplexSparse, Allocator>& A);

  template<class T, class Allocator>
  void Transpose(const Matrix<T, Symmetric,
                 ArrayRowSymComplexSparse, Allocator>& A,
                 Matrix<T, Symmetric, ArrayRowSymComplexSparse, Allocator>& B);

  template<class T, class Allocator>
  void Transpose(Matrix<T, General, ArrayRowComplexSparse, Allocator>& A);

  template<class T, class Allocator>
  void Transpose(Matrix<T, Symmetric,
                 ArrayRowSymComplexSparse, Allocator>& A);

  template<class T, class Allocator>
  void Transpose(const Matrix<T, Symmetric,
                 RowSymComplexSparse, Allocator>& A,
                 Matrix<T, Symmetric, RowSymComplexSparse, Allocator>& B);

  template<class T, class Allocator>
  void Transpose(Matrix<T, General, RowComplexSparse, Allocator>& A);

  template<class T, class Allocator>
  void Transpose(Matrix<T, Symmetric,
                 RowSymComplexSparse, Allocator>& A);

  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, RowComplexSparse, Allocator>& A);

  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ColComplexSparse, Allocator>& A);
  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A);
  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ColSymComplexSparse, Allocator>& A);
  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A);
  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ArrayColComplexSparse, Allocator>& A);

  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A);
  
  template<class T, class Prop, class Allocator>
  bool IsComplexMatrix(const Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>& A);
  
}

#define SELDON_FILE_FUNCTIONS_MATRIX_COMPLEX_HXX
#endif


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


#ifndef SELDON_FILE_LAPACK_LINEAREQUATIONS_HXX


namespace Seldon
{

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, ColMajor, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, ColMajor, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, RowMajor, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, RowMajor, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, ColSym, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, ColSym, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<complex<float>, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<complex<double>, Prop0, ColSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, RowSym, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, RowSym, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<complex<float>, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(SeldonUplo Uplo,
	     Matrix<complex<double>, Prop0, RowSymPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, ColHermPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, ColHermPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<float>, Prop0, RowHermPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void GetLU(Matrix<complex<double>, Prop0, RowHermPacked, Allocator0>& A,
	     Vector<int, VectFull, Allocator1>& P,
	     LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const SeldonTranspose& TransA,
	       const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const SeldonTranspose& TransA,
	       const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const SeldonTranspose& TransA,
	       const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const SeldonTranspose& TransA,
	       const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const SeldonTranspose& TransA,
	       const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const SeldonTranspose& TransA,
	       const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const SeldonTranspose& TransA,
	       const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const SeldonTranspose& TransA,
	       const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, ColSym, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, ColSym, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColSymPacked,
	       Allocator0>& A, const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColSymPacked,
	       Allocator0>& A, const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, RowSym, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, RowSym, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, RowSymPacked,
	       Allocator0>& A, const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, RowSymPacked,
	       Allocator0>& A, const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColHermPacked,
	       Allocator0>& A, const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColHermPacked,
	       Allocator0>& A, const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	       const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, RowHermPacked,
	       Allocator0>& A, const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, RowHermPacked,
	       Allocator0>& A, const Vector<int, VectFull, Allocator1>& P,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColUpTriang,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColUpTriang,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<float>, Prop0, ColUpTriang,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<double>, Prop0, ColUpTriang,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColLoTriang,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColLoTriang,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<float>, Prop0, ColLoTriang,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<double>, Prop0, ColLoTriang,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColUpTriangPacked,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColUpTriangPacked,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<float>, Prop0, ColUpTriangPacked,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<double>, Prop0, ColUpTriangPacked,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, ColLoTriangPacked,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, ColLoTriangPacked,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<float>, Prop0, ColLoTriangPacked,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<double>, Prop0, ColLoTriangPacked,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void
  SolveLU(const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	  Vector<complex<float>, VectFull, Allocator2>& b,
	  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void
  SolveLU(const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	  Vector<complex<double>, VectFull, Allocator2>& b,
	  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void
  SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	  const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	  Vector<complex<float>, VectFull, Allocator2>& b,
	  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void
  SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	  const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	  Vector<complex<double>, VectFull, Allocator2>& b,
	  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void
  SolveLU(const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	  Vector<complex<float>, VectFull, Allocator2>& b,
	  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void
  SolveLU(const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	  Vector<complex<double>, VectFull, Allocator2>& b,
	  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void
  SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	  const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	  Vector<complex<float>, VectFull, Allocator2>& b,
	  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void
  SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	  const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	  Vector<complex<double>, VectFull, Allocator2>& b,
	  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, RowUpTriangPacked,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, RowUpTriangPacked,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<float>, Prop0, RowUpTriangPacked,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<double>, Prop0, RowUpTriangPacked,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<float>, Prop0, RowLoTriangPacked,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const Matrix<complex<double>, Prop0, RowLoTriangPacked,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<float>, Prop0, RowLoTriangPacked,
	       Allocator0>& A,
	       Vector<complex<float>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0, class Allocator2>
  void SolveLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
	       const Matrix<complex<double>, Prop0, RowLoTriangPacked,
	       Allocator0>& A,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<float, Prop0, ColMajor,
				  Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, ColMajor,
				   Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColMajor, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColMajor, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<float, Prop0, RowMajor,
				  Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, RowMajor,
				   Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  RowMajor, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowMajor, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<float, Prop0, ColSym,
				  Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, ColSym,
				   Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColSym, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColSym, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<float, Prop0,
				  ColSymPacked, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, ColSymPacked,
				   Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColSymPacked, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColSymPacked, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<float, Prop0, RowSym,
				  Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, RowSym,
				   Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0, RowSym,
				  Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowSym, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<float, Prop0, RowSymPacked,
				  Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, RowSymPacked,
				   Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  RowSymPacked, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowSymPacked, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColHerm, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColHerm, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm,  double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColHermPacked, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColHermPacked, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  RowHerm, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowHerm, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);
  
  template<class Prop0, class Allocator0, class Allocator1>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  RowHermPacked, Allocator0>& A,
				  Vector<int, VectFull, Allocator1>& P,
				  SeldonNorm norm, float anorm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0, class Allocator1>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowHermPacked, Allocator0>& A,
				   Vector<int, VectFull, Allocator1>& P,
				   SeldonNorm norm, double anorm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<float, Prop0, ColUpTriang,
				  Allocator0>& A, SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, ColUpTriang,
				   Allocator0>& A, SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColUpTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColUpTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<float, Prop0, ColUpTriang,
				  Allocator0>& A, SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<double, Prop0, ColUpTriang,
				   Allocator0>& A, SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<complex<float>, Prop0,
				  ColUpTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<complex<double>, Prop0,
				   ColUpTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<float, Prop0,
				  ColLoTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, ColLoTriang,
				   Allocator0>& A, SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColLoTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColLoTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<float, Prop0,
				  ColLoTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<double, Prop0, ColLoTriang,
				   Allocator0>& A, SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<complex<float>, Prop0,
				  ColLoTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<complex<double>, Prop0,
				   ColLoTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<float, Prop0,
				  ColUpTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<double, Prop0,
				   ColUpTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColUpTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColUpTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<float, Prop0,
				  ColUpTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<double, Prop0,
				   ColUpTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<complex<float>, Prop0,
				  ColUpTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<complex<double>, Prop0,
				   ColUpTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<float, Prop0,
				  ColLoTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<double, Prop0,
				   ColLoTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  ColLoTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   ColLoTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<float, Prop0,
				  ColLoTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<double, Prop0,
				   ColLoTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<complex<float>, Prop0,
				  ColLoTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<complex<double>, Prop0,
				   ColLoTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<float, Prop0, RowUpTriang,
				  Allocator0>& A, SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<double, Prop0, RowUpTriang,
				   Allocator0>& A, SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  RowUpTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowUpTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<float, Prop0, RowUpTriang,
				  Allocator0>& A, SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<double, Prop0, RowUpTriang,
				   Allocator0>& A, SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<complex<float>, Prop0,
				  RowUpTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<complex<double>, Prop0,
				   RowUpTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<float, Prop0,
				  RowLoTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<double, Prop0,
				   RowLoTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  RowLoTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowLoTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<float, Prop0,
				  RowLoTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<double, Prop0,
				   RowLoTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<complex<float>, Prop0,
				  RowLoTriang, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<complex<double>, Prop0,
				   RowLoTriang, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<float, Prop0,
				  RowUpTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<double, Prop0,
				   RowUpTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  RowUpTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowUpTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<float, Prop0,
				  RowUpTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<double, Prop0,
				   RowUpTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<complex<float>, Prop0,
				  RowUpTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<complex<double>, Prop0,
				   RowUpTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<float, Prop0,
				  RowLoTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<double, Prop0,
				   RowLoTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const Matrix<complex<float>, Prop0,
				  RowLoTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const Matrix<complex<double>, Prop0,
				   RowLoTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<float, Prop0,
				  RowLoTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<double, Prop0,
				   RowLoTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  float ReciprocalConditionNumber(const SeldonDiag& DiagA,
				  const Matrix<complex<float>, Prop0,
				  RowLoTriangPacked, Allocator0>& A,
				  SeldonNorm norm,
				  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0>
  double ReciprocalConditionNumber(const SeldonDiag& DiagA,
				   const Matrix<complex<double>, Prop0,
				   RowLoTriangPacked, Allocator0>& A,
				   SeldonNorm norm,
				   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<float, Prop0, ColMajor, Allocator0>& A,
			const Matrix<float, Prop0, ColMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<float, VectFull, Allocator3>& x,
			const Vector<float, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
			const Matrix<double, Prop0, ColMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<double, VectFull, Allocator3>& x,
			const Vector<double, VectFull, Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			ColMajor, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			ColMajor, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			const Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			ColMajor, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			ColMajor, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const SeldonTranspose& TransA,
			const Matrix<float, Prop0, ColMajor, Allocator0>& A,
			const Matrix<float, Prop0, ColMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<float, VectFull, Allocator3>& x,
			const Vector<float, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const SeldonTranspose& TransA,
			const Matrix<double, Prop0, ColMajor, Allocator0>& A,
			const Matrix<double, Prop0, ColMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<double, VectFull, Allocator3>& x,
			const Vector<double, VectFull, Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void
  RefineSolutionLU(const SeldonTranspose& TransA,
		   const Matrix<complex<float>, Prop0, ColMajor,
		   Allocator0>& A,
		   const Matrix<complex<float>, Prop0, ColMajor,
		   Allocator1>& Alu,
		   const Vector<int, VectFull, Allocator2>& P,
		   Vector<complex<float>, VectFull, Allocator3>& x,
		   const Vector<complex<float>, VectFull, Allocator4>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const SeldonTranspose& TransA,
			const Matrix<complex<double>, Prop0, ColMajor,
			Allocator0>& A,
			const Matrix<complex<double>, Prop0, ColMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
			const Matrix<float, Prop0, RowMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<float, VectFull, Allocator3>& x,
			const Vector<float, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
			const Matrix<double, Prop0, RowMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<double, VectFull, Allocator3>& x,
			const Vector<double, VectFull, Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			RowMajor, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			RowMajor, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			const Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			RowMajor, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			RowMajor, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const SeldonTranspose& TransA,
			const Matrix<float, Prop0, RowMajor, Allocator0>& A,
			const Matrix<float, Prop0, RowMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<float, VectFull, Allocator3>& x,
			const Vector<float, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const SeldonTranspose& TransA,
			const Matrix<double, Prop0, RowMajor, Allocator0>& A,
			const Matrix<double, Prop0, RowMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<double, VectFull, Allocator3>& x,
			const Vector<double, VectFull, Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const SeldonTranspose& TransA,
			const Matrix<complex<float>, Prop0, RowMajor,
			Allocator0>& A,
			const Matrix<complex<float>, Prop0, RowMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			Vector<complex<float>, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const SeldonTranspose& TransA,
			const Matrix<complex<double>, Prop0, RowMajor,
			Allocator0>& A,
			const Matrix<complex<double>, Prop0, RowMajor,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<float, Prop0, ColSym, Allocator0>& A,
			const Matrix<float, Prop0, ColSym,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<float, VectFull, Allocator3>& x,
			const Vector<float, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<double, Prop0, ColSym, Allocator0>& A,
			const Matrix<double, Prop0, ColSym,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<double, VectFull, Allocator3>& x,
			const Vector<double, VectFull, Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			ColSym, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			ColSym, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			const Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			ColSym, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			ColSym, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<float, Prop0, ColSymPacked,
			Allocator0>& A,
			const Matrix<float, Prop0, ColSymPacked,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<float, VectFull, Allocator3>& x,
			const Vector<float, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<double, Prop0, ColSymPacked,
			Allocator0>& A,
			const Matrix<double, Prop0, ColSymPacked,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<double, VectFull, Allocator3>& x,
			const Vector<double, VectFull, Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			ColSymPacked, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			ColSymPacked, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			const Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			ColSymPacked, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			ColSymPacked, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<float, Prop0, RowSym, Allocator0>& A,
			const Matrix<float, Prop0, RowSym,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<float, VectFull, Allocator3>& x,
			const Vector<float, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<double, Prop0, RowSym, Allocator0>& A,
			const Matrix<double, Prop0, RowSym,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<double, VectFull, Allocator3>& x,
			const Vector<double, VectFull, Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			RowSym, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			RowSym, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			const Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			RowSym, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			RowSym, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<float, Prop0, RowSymPacked,
			Allocator0>& A,
			const Matrix<float, Prop0, RowSymPacked,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<float, VectFull, Allocator3>& x,
			const Vector<float, VectFull, Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<double, Prop0, RowSymPacked,
			Allocator0>& A,
			const Matrix<double, Prop0, RowSymPacked,
			Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<double, VectFull, Allocator3>& x,
			const Vector<double, VectFull, Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			RowSymPacked, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			RowSymPacked, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			const Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			RowSymPacked, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			RowSymPacked, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			ColHerm, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			ColHerm, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			const Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			ColHerm, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			ColHerm, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			ColHermPacked, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			ColHermPacked, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			const Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			ColHermPacked, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			ColHermPacked, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			const Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			RowHerm, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			RowHerm, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			RowHerm, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			RowHerm, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<float>, Prop0,
			RowHermPacked, Allocator0>& A,
			const Matrix<complex<float>, Prop0,
			RowHermPacked, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<float>, VectFull, Allocator3>& x,
			Vector<complex<float>, VectFull,
			Allocator4>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2,
	    class Allocator3, class Allocator4>
  void RefineSolutionLU(const Matrix<complex<double>, Prop0,
			RowHermPacked, Allocator0>& A,
			const Matrix<complex<double>, Prop0,
			RowHermPacked, Allocator1>& Alu,
			const Vector<int, VectFull, Allocator2>& P,
			Vector<complex<double>, VectFull, Allocator3>& x,
			Vector<complex<double>, VectFull,
			Allocator4>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<float, Prop0, ColUpTriang,
			Allocator0>& A,
			Vector<float, VectFull, Allocator1>& x,
			const Vector<float, VectFull, Allocator2>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<double, Prop0, ColUpTriang,
			Allocator0>& A,
			Vector<double, VectFull, Allocator1>& x,
			const Vector<double, VectFull, Allocator2>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<float>, Prop0, ColUpTriang,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<double>, Prop0, ColUpTriang,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<float, Prop0, ColUpTriang,
		   Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& x,
		   const Vector<float, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<double, Prop0, ColUpTriang,
		   Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& x,
		   const Vector<double, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<float>, Prop0, ColUpTriang,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<double>, Prop0, ColUpTriang,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<float, Prop0, ColLoTriang,
			Allocator0>& A,
			Vector<float, VectFull, Allocator1>& x,
			const Vector<float, VectFull, Allocator2>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<double, Prop0, ColLoTriang,
			Allocator0>& A,
			Vector<double, VectFull, Allocator1>& x,
			const Vector<double, VectFull, Allocator2>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<float>, Prop0, ColLoTriang,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<double>, Prop0, ColLoTriang,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<float, Prop0, ColLoTriang,
		   Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& x,
		   const Vector<float, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<double, Prop0, ColLoTriang,
		   Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& x,
		   const Vector<double, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<float>, Prop0, ColLoTriang,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<double>, Prop0, ColLoTriang,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<float, Prop0, ColUpTriangPacked,
			Allocator0>& A,
			Vector<float, VectFull, Allocator1>& x,
			const Vector<float, VectFull, Allocator2>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<double, Prop0, ColUpTriangPacked,
			Allocator0>& A,
			Vector<double, VectFull, Allocator1>& x,
			const Vector<double, VectFull, Allocator2>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<float>, Prop0, ColUpTriangPacked,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<double>, Prop0, ColUpTriangPacked,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<float, Prop0, ColUpTriangPacked,
		   Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& x,
		   const Vector<float, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<double, Prop0, ColUpTriangPacked,
		   Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& x,
		   const Vector<double, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<float>, Prop0, ColUpTriangPacked,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<double>, Prop0, ColUpTriangPacked,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<float, Prop0, ColLoTriangPacked,
			Allocator0>& A,
			Vector<float, VectFull, Allocator1>& x,
			const Vector<float, VectFull, Allocator2>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<double, Prop0, ColLoTriangPacked,
			Allocator0>& A,
			Vector<double, VectFull, Allocator1>& x,
			const Vector<double, VectFull, Allocator2>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<float>, Prop0, ColLoTriangPacked,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<double>, Prop0, ColLoTriangPacked,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<float, Prop0, ColLoTriangPacked,
		   Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& x,
		   const Vector<float, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<double, Prop0, ColLoTriangPacked,
		   Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& x,
		   const Vector<double, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<float>, Prop0, ColLoTriangPacked,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<double>, Prop0, ColLoTriangPacked,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<float, Prop0, RowUpTriang,
			Allocator0>& A,
			Vector<float, VectFull, Allocator1>& x,
			const Vector<float, VectFull, Allocator2>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<double, Prop0, RowUpTriang,
			Allocator0>& A,
			Vector<double, VectFull, Allocator1>& x,
			const Vector<double, VectFull, Allocator2>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<float>, Prop0, RowUpTriang,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<double>, Prop0, RowUpTriang,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<float, Prop0, RowUpTriang,
		   Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& x,
		   const Vector<float, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<double, Prop0, RowUpTriang,
		   Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& x,
		   const Vector<double, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<float>, Prop0, RowUpTriang,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<double>, Prop0, RowUpTriang,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<float, Prop0, RowLoTriang,
			Allocator0>& A,
			Vector<float, VectFull, Allocator1>& x,
			const Vector<float, VectFull, Allocator2>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<double, Prop0, RowLoTriang,
			Allocator0>& A,
			Vector<double, VectFull, Allocator1>& x,
			const Vector<double, VectFull, Allocator2>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<float>, Prop0, RowLoTriang,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<double>, Prop0, RowLoTriang,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<float, Prop0, RowLoTriang,
		   Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& x,
		   const Vector<float, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<double, Prop0, RowLoTriang,
		   Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& x,
		   const Vector<double, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<float>, Prop0, RowLoTriang,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<double>, Prop0, RowLoTriang,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<float, Prop0, RowUpTriangPacked,
			Allocator0>& A,
			Vector<float, VectFull, Allocator1>& x,
			const Vector<float, VectFull, Allocator2>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<double, Prop0, RowUpTriangPacked,
			Allocator0>& A,
			Vector<double, VectFull, Allocator1>& x,
			const Vector<double, VectFull, Allocator2>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<float>, Prop0, RowUpTriangPacked,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<double>, Prop0, RowUpTriangPacked,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<float, Prop0, RowUpTriangPacked,
		   Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& x,
		   const Vector<float, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<double, Prop0, RowUpTriangPacked,
		   Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& x,
		   const Vector<double, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<float>, Prop0, RowUpTriangPacked,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<double>, Prop0, RowUpTriangPacked,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<float, Prop0, RowLoTriangPacked,
			Allocator0>& A,
			Vector<float, VectFull, Allocator1>& x,
			const Vector<float, VectFull, Allocator2>& b,
			float& ferr, float& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void RefineSolutionLU(const Matrix<double, Prop0, RowLoTriangPacked,
			Allocator0>& A,
			Vector<double, VectFull, Allocator1>& x,
			const Vector<double, VectFull, Allocator2>& b,
			double& ferr, double& berr,
			LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<float>, Prop0, RowLoTriangPacked,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   const Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const Matrix<complex<double>, Prop0, RowLoTriangPacked,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   const Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<float, Prop0, RowLoTriangPacked,
		   Allocator0>& A,
		   Vector<float, VectFull, Allocator1>& x,
		   const Vector<float, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<double, Prop0, RowLoTriangPacked,
		   Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& x,
		   const Vector<double, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<float>, Prop0, RowLoTriangPacked,
		   Allocator0>& A,
		   Vector<complex<float>, VectFull, Allocator1>& x,
		   Vector<complex<float>, VectFull, Allocator2>& b,
		   float& ferr, float& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void
  RefineSolutionLU(const SeldonTranspose& TransA, const SeldonDiag& DiagA,
		   const Matrix<complex<double>, Prop0, RowLoTriangPacked,
		   Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& x,
		   Vector<complex<double>, VectFull, Allocator2>& b,
		   double& ferr, double& berr,
		   LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, ColMajor, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, RowMajor, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, RowMajor, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, ColSym, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, ColSym, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColSymPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColSymPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, RowSym, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, RowSym, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowSymPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowSymPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColHermPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColHermPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowHermPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowHermPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColUpTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColUpTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<float>, Prop0, ColUpTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<double>, Prop0, ColUpTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, ColLoTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, ColLoTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<float>, Prop0, ColLoTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<double>, Prop0, ColLoTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowUpTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowUpTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<float>, Prop0, RowUpTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<double>, Prop0, RowUpTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<float>, Prop0, RowLoTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(Matrix<complex<double>, Prop0, RowLoTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<float>, Prop0, RowLoTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template <class Prop0, class Allocator0>
  void GetInverse(const SeldonDiag& DiagA,
		  Matrix<complex<double>, Prop0, RowLoTriangPacked,
		  Allocator0>& A,
		  LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void GetScalingFactors(const Matrix<float, Prop0, ColMajor, Allocator0>& A,
			 Vector<float, VectFull, Allocator1>& row_scale,
			 Vector<float, VectFull, Allocator2>& col_scale,
			 float& row_condition_number,
			 float& col_condition_number, float& amax,
			 LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void GetScalingFactors(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
			 Vector<double, VectFull, Allocator1>& row_scale,
			 Vector<double, VectFull, Allocator2>& col_scale,
			 double& row_condition_number,
			 double& col_condition_number, double& amax,
			 LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void GetScalingFactors(const Matrix<complex<float>, Prop0, ColMajor,
			 Allocator0>& A,
			 Vector<float, VectFull, Allocator1>& row_scale,
			 Vector<float, VectFull, Allocator2>& col_scale,
			 float& row_condition_number,
			 float& col_condition_number, float& amax,
			 LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void GetScalingFactors(const Matrix<complex<double>, Prop0, ColMajor,
			 Allocator0>& A,
			 Vector<double, VectFull, Allocator1>& row_scale,
			 Vector<double, VectFull, Allocator2>& col_scale,
			 double& row_condition_number,
			 double& col_condition_number, double& amax,
			 LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void GetScalingFactors(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
			 Vector<float, VectFull, Allocator1>& row_scale,
			 Vector<float, VectFull, Allocator2>& col_scale,
			 float& row_condition_number,
			 float& col_condition_number, float& amax,
			 LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void GetScalingFactors(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
			 Vector<double, VectFull, Allocator1>& row_scale,
			 Vector<double, VectFull, Allocator2>& col_scale,
			 double& row_condition_number,
			 double& col_condition_number, double& amax,
			 LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void GetScalingFactors(const Matrix<complex<float>, Prop0, RowMajor,
			 Allocator0>& A,
			 Vector<float, VectFull, Allocator1>& row_scale,
			 Vector<float, VectFull, Allocator2>& col_scale,
			 float& row_condition_number,
			 float& col_condition_number, float& amax,
			 LapackInfo& info = lapack_info);

  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2>
  void GetScalingFactors(const Matrix<complex<double>, Prop0, RowMajor,
			 Allocator0>& A,
			 Vector<double, VectFull, Allocator1>& row_scale,
			 Vector<double, VectFull, Allocator2>& col_scale,
			 double& row_condition_number,
			 double& col_condition_number, double& amax,
			 LapackInfo& info = lapack_info);

  template<class Prop, class Allocator>
  void GetCholesky(Matrix<double, Prop, RowSymPacked, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class Prop, class Allocator>
  void GetCholesky(Matrix<double, Prop, ColSymPacked, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class Prop, class Allocator>
  void GetCholesky(Matrix<double, Prop, RowSym, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class Prop, class Allocator>
  void GetCholesky(Matrix<double, Prop, ColSym, Allocator>& A,
		   LapackInfo& info = lapack_info);
  
  template<class Prop, class Allocator>
  void GetCholesky(Matrix<complex<double>, Prop, RowHermPacked, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class Prop, class Allocator>
  void GetCholesky(Matrix<complex<double>, Prop, ColHermPacked, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class Prop, class Allocator>
  void GetCholesky(Matrix<complex<double>, Prop, RowHerm, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class Prop, class Allocator>
  void GetCholesky(Matrix<complex<double>, Prop, ColHerm, Allocator>& A,
		   LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void SolveCholesky(const Transp& TransA,
		     const Matrix<double, Prop, RowSymPacked, Allocator>& A,
		     Vector<double, VectFull, Allocator2>& X,
		     LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void SolveCholesky(const Transp& TransA,
		     const Matrix<double, Prop, ColSymPacked, Allocator>& A,
		     Vector<double, VectFull, Allocator2>& X,
		     LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void SolveCholesky(const Transp& TransA,
		     const Matrix<double, Prop, RowSym, Allocator>& A,
		     Vector<double, VectFull, Allocator2>& X,
		     LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void SolveCholesky(const Transp& TransA,
		     const Matrix<double, Prop, ColSym, Allocator>& A,
		     Vector<double, VectFull, Allocator2>& X,
		     LapackInfo& info = lapack_info);
  
  template<class Transp, class Prop, class Allocator, class Allocator2>
  void SolveCholesky(const Transp& TransA,
		     const Matrix<complex<double>, Prop, RowHermPacked, Allocator>& A,
		     Vector<complex<double>, VectFull, Allocator2>& X,
		     LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void SolveCholesky(const Transp& TransA,
		     const Matrix<complex<double>, Prop, ColHermPacked, Allocator>& A,
		     Vector<complex<double>, VectFull, Allocator2>& X,
		     LapackInfo& info = lapack_info);
  
  template<class Transp, class Prop, class Allocator, class Allocator2>
  void SolveCholesky(const Transp& TransA,
		     const Matrix<complex<double>, Prop, RowHerm, Allocator>& A,
		     Vector<complex<double>, VectFull, Allocator2>& X,
		     LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void SolveCholesky(const Transp& TransA,
		     const Matrix<complex<double>, Prop, ColHerm, Allocator>& A,
		     Vector<complex<double>, VectFull, Allocator2>& X,
		     LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void MltCholesky(const Transp& TransA,
                   const Matrix<double, Prop, RowSymPacked, Allocator>& A,
                   Vector<double, VectFull, Allocator2>& X,
                   LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void MltCholesky(const Transp& TransA,
                   const Matrix<double, Prop, ColSymPacked, Allocator>& A,
                   Vector<double, VectFull, Allocator2>& X,
                   LapackInfo& info = lapack_info);
  
  template<class Transp, class Prop, class Allocator, class Allocator2>
  void MltCholesky(const Transp& TransA,
                   const Matrix<double, Prop, RowSym, Allocator>& A,
                   Vector<double, VectFull, Allocator2>& X,
                   LapackInfo& info = lapack_info);

  template<class Transp, class Prop, class Allocator, class Allocator2>
  void MltCholesky(const Transp& TransA,
                   const Matrix<double, Prop, ColSym, Allocator>& A,
                   Vector<double, VectFull, Allocator2>& X,
                   LapackInfo& info = lapack_info);
  
  template<class Transp, class Prop, class Allocator, class Allocator2>
  void MltCholesky(const Transp& TransA,
                   const Matrix<complex<double>,
		   Prop, RowHermPacked, Allocator>& A,
                   Vector<complex<double>, VectFull, Allocator2>& X,
                   LapackInfo& info = lapack_info);
  
  template<class Transp, class Prop, class Allocator, class Allocator2>
  void MltCholesky(const Transp& TransA,
                   const Matrix<complex<double>,
		   Prop, ColHermPacked, Allocator>& A,
                   Vector<complex<double>, VectFull, Allocator2>& X,
                   LapackInfo& info = lapack_info);
  
  template<class Transp, class Prop, class Allocator, class Allocator2>
  void MltCholesky(const Transp& TransA,
                   const Matrix<complex<double>,
		   Prop, RowHerm, Allocator>& A,
                   Vector<complex<double>, VectFull, Allocator2>& X,
                   LapackInfo& info = lapack_info);
  
  template<class Transp, class Prop, class Allocator, class Allocator2>
  void MltCholesky(const Transp& TransA,
                   const Matrix<complex<double>,
		   Prop, ColHerm, Allocator>& A,
                   Vector<complex<double>, VectFull, Allocator2>& X,
                   LapackInfo& info = lapack_info);

  // Generic method, which factorizes a matrix and solve the linear system
  // b is overwritten by the solution
  template<class T, class Prop, class Storage, class Allocator,
	   class Allocator1, class Allocator2>
  void GetAndSolveLU(Matrix<T, Prop, Storage, Allocator>& A,
		     Vector<int, VectFull, Allocator1>& P,
		     Vector<T, VectFull, Allocator2>& b,
		     LapackInfo& info = lapack_info);

} // namespace Seldon.

#define SELDON_FILE_LAPACK_LINEAREQUATIONS_HXX
#endif

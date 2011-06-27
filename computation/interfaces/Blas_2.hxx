// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2001-2011 Marc Duruflé
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


#ifndef SELDON_FILE_BLAS_2_HXX

namespace Seldon
{

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, ColUpTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, ColUpTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, ColUpTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, ColUpTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, ColLoTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, ColLoTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, ColLoTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, ColLoTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, RowUpTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, RowUpTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, RowUpTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, RowUpTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, RowLoTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, RowLoTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, RowLoTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, RowLoTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0,
	      ColHermPacked, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0,
	      ColHermPacked, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0,
	      ColHermPacked, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0,
	      ColHermPacked, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0,
	      RowHermPacked, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0,
	      RowHermPacked, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0,
	      RowHermPacked, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0,
	      RowHermPacked, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator1>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop0, ColMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator1>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop0, ColMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<float> alpha,
		   const Vector<complex<float>, VectFull, Allocator1>& X,
		   const Vector<complex<float>, VectFull, Allocator2>& Y,
		   Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<double> alpha,
		   const Vector<complex<double>, VectFull, Allocator1>& X,
		   const Vector<complex<double>, VectFull, Allocator2>& Y,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<float> alpha,
		   const Vector<complex<float>, VectFull, Allocator1>& X,
		   const SeldonConjugate& ConjY,
		   const Vector<complex<float>, VectFull, Allocator2>& Y,
		   Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<double> alpha,
		   const Vector<complex<double>, VectFull, Allocator1>& X,
		   const SeldonConjugate& ConjY,
		   const Vector<complex<double>, VectFull, Allocator2>& Y,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator1>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop0, RowMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator1>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop0, RowMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<float> alpha,
		   const Vector<complex<float>, VectFull, Allocator1>& X,
		   const Vector<complex<float>, VectFull, Allocator2>& Y,
		   Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<double> alpha,
		   const Vector<complex<double>, VectFull, Allocator1>& X,
		   const Vector<complex<double>, VectFull, Allocator2>& Y,
		   Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<float> alpha,
		   const Vector<complex<float>, VectFull, Allocator1>& X,
		   const SeldonConjugate& ConjY,
		   const Vector<complex<float>, VectFull, Allocator2>& Y,
		   Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<double> alpha,
		   const Vector<complex<double>, VectFull, Allocator1>& X,
		   const SeldonConjugate& ConjY,
		   const Vector<complex<double>, VectFull, Allocator2>& Y,
		   Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   Matrix<float, Prop1, ColSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   Matrix<double, Prop1, ColSymPacked, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   Matrix<float, Prop1, ColSym, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   Matrix<double, Prop1, ColSym, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      Matrix<complex<float>, Prop1, ColHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      Matrix<complex<double>, Prop1, ColHermPacked, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      Matrix<complex<float>, Prop1, ColHerm, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      Matrix<complex<double>, Prop1, ColHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, ColSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, ColSymPacked, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, ColSym, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, ColSym, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, ColHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, ColHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, ColHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, ColHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   Matrix<float, Prop1, RowSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   Matrix<double, Prop1, RowSymPacked, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   Matrix<float, Prop1, RowSym, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   Matrix<double, Prop1, RowSym, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      Matrix<complex<float>, Prop1, RowHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      Matrix<complex<double>, Prop1, RowHermPacked, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      Matrix<complex<float>, Prop1, RowHerm, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      Matrix<complex<double>, Prop1, RowHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, RowSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, RowSymPacked, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, RowSym, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, RowSym, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, RowHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, RowHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, RowHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, RowHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop1, ColSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop1, ColSymPacked, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop1, ColSym, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop1, ColSym, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      Matrix<complex<float>, Prop1, ColHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      Matrix<complex<double>, Prop1, ColHermPacked, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      Matrix<complex<float>, Prop1, ColHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      Matrix<complex<double>, Prop1, ColHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, ColSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, ColSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, ColSym, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, ColSym, Allocator1>& A);
  
  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, ColHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, ColHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, ColHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, ColHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop1, RowSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop1, RowSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop1, RowSym, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop1, RowSym, Allocator1>& A);
 
  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      Matrix<complex<float>, Prop1, RowHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      Matrix<complex<double>, Prop1, RowHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      Matrix<complex<float>, Prop1, RowHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      Matrix<complex<double>, Prop1, RowHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, RowSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, RowSymPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, RowSym, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, RowSym, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, RowHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, RowHermPacked, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, RowHerm, Allocator1>& A);

  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, RowHerm, Allocator1>& A);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	ColUpTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	ColUpTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	ColLoTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	ColLoTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	RowUpTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	RowUpTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	RowLoTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	RowLoTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, ColUpTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	ColUpTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, ColUpTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	ColUpTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, ColLoTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	ColLoTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, ColLoTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	ColLoTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, RowUpTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	RowUpTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowUpTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	RowUpTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, RowLoTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	RowLoTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowLoTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X);

  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	RowLoTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X);

} // namespace Seldon.

#define SELDON_FILE_BLAS_2_HXX
#endif

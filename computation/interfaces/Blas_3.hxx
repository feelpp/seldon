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


#ifndef SELDON_FILE_BLAS_3_HXX

namespace Seldon
{

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const float& alpha,
	      const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const double& alpha,
	      const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const float& alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const double& alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<float>& alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<double>& alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const float& alpha,
	      const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const double& alpha,
	      const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const float& alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const double& alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<float>& alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const complex<double>& alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const float& alpha,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const double& alpha,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const float& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Matrix<float, Prop1, ColMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const double& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Matrix<double, Prop1, ColMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, ColSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const float& alpha,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const double& alpha,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const float& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Matrix<float, Prop1, RowMajor, Allocator1>& B,
	      const float& beta,
	      Matrix<float, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const double& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Matrix<double, Prop1, RowMajor, Allocator1>& B,
	      const double& beta,
	      Matrix<double, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, RowSym, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<float>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B,
	      const complex<float>& beta,
	      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1,
	    class Prop2, class Allocator2>
  void MltAdd(const SeldonSide& Side,
	      const complex<double>& alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B,
	      const complex<double>& beta,
	      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& C);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const float& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<float, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const double& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<double, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<float>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Mlt(const SeldonSide& Side,
	   const complex<double>& alpha,
	   const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	   Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, ColMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const float& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<float, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const double& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<double, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<float>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<float>, Prop1, RowMajor, Allocator1>& B);

  template <class Prop0, class Allocator0,
	    class Prop1, class Allocator1>
  void Solve(const SeldonSide& Side,
	     const complex<double>& alpha,
	     const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
	     Matrix<complex<double>, Prop1, RowMajor, Allocator1>& B);

} // namespace Seldon.

#define SELDON_FILE_BLAS_3_HXX
#endif

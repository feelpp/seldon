// Copyright (C) 2003-2011 Marc Duruflé
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

#ifndef SELDON_FILE_MATRIX_COMPLEX_CONVERSIONS_HXX

namespace Seldon
{
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);


  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<complex<T>, VectFull, Allocator4>& Val,
			       int index = 0, bool sym = false);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowComplexSparse,
				 Allocator4>& A,
				 int index = 0);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColComplexSparse,
				 Allocator4>& A,
				 int index = 0);


  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymComplexSparse,
				 Allocator4>& A,
				 int index = 0);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymComplexSparse,
				 Allocator4>& A,
				 int index = 0);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowComplexSparse,
				 Allocator4>& A,
				 int index = 0);
  
  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColComplexSparse,
				 Allocator4>& A,
				 int index = 0);

  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymComplexSparse,
				 Allocator4>& A,
				 int index = 0);

  
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<complex<T>, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymComplexSparse,
				 Allocator4>& A,
				 int index = 0);  
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr);
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0,
       ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr);
  
  

} 

#define SELDON_FILE_MATRIX_COMPLEX_CONVERSIONS_HXX
#endif 

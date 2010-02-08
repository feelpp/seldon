// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2010 Vivien Mallet
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


#ifndef SELDON_FILE_MATRIX_CONVERSIONS_HXX


namespace Seldon
{


  /*
    From CSR formats to "Matlab" coordinate format.
  */


  //! Conversion from RowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse,
			       Allocator1>& A,
			       Vector<int>& IndRow, Vector<int>& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false);


  //! Conversion from ColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse,
			       Allocator1>& A,
			       Vector<int>& IndRow, Vector<int>& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false);


  //! Conversion from RowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       RowSymSparse, Allocator1>& A,
			       Vector<int>& IndRow, Vector<int>& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false);


  //! Conversion from ColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ColSymSparse, Allocator1>& A,
			       Vector<int>& IndRow, Vector<int>& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false);


  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from ArrayRowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ArrayRowSparse, Allocator1>& A,
			       Vector<int>& IndRow, Vector<int>& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false);


  //! Conversion from ArrayColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ArrayColSparse, Allocator1>& A,
			       Vector<int>& IndRow, Vector<int>& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false);


  //! Conversion from ArrayRowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ArrayRowSymSparse, Allocator1>& A,
			       Vector<int>& IndRow, Vector<int>& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false);


  //! Conversion from ArrayColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ArrayColSymSparse, Allocator1>& A,
			       Vector<int>& IndRow, Vector<int>& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false);


  /*
    From "Matlab" coordinate format to CSR formats.
  */


  //! Conversion from coordinate format to RowSparse.
  /*! Contrary to the other conversion functions
    ConvertMatrix_from_Coordinates, this one accepts duplicates.
    \param[in] IndRow row indexes of the non-zero elements.
    \param[in] IndCol column indexes of the non-zero elements.
    \param[in] Val values of the non-zero elements.
    \param[out] A matrix defined by \a IndRow, \a IndCol and \a Val.
    \param[in] index index of the first column and of the first row.
  */
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, RowSparse, Allocator>& A,
				 int index = 0);


  //! Conversion from coordinate format to ColSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, ColSparse, Allocator>& A,
				 int index = 0);


  //! Conversion from coordinate format to RowSymSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, RowSymSparse, Allocator>& A,
				 int index = 0);


  //! Conversion from coordinate format to ColSymSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, ColSymSparse, Allocator>& A,
				 int index = 0);


  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from coordinate format to ArrayRowSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, ArrayRowSparse,
				 Allocator>& A,
				 int index = 0);


  //! Conversion from coordinate format to ArrayColSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, ArrayColSparse,
				 Allocator>& A,
				 int index = 0);


  //! Conversion from coordinate format to ArrayRowSymSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop,
				 ArrayRowSymSparse, Allocator>& A,
				 int index = 0);


  //! Conversion from coordinate format to ArrayColSymSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(Vector<int>& IndRow, Vector<int>& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop,
				 ArrayColSymSparse, Allocator>& A,
				 int index = 0);


  /*
    From CSR to other CSR formats.
  */


  //! B = A.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop, ColSparse, Alloc1>& A,
	    Matrix<T, Prop, ColSparse, Alloc2>& B);


  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop, RowSparse, Alloc1>& A,
	    Matrix<T, Prop, ColSparse, Alloc2>& B);


  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop1, class Prop2,
	   class Storage, class Alloc1, class Alloc2>
  void
  ConvertMatrixSymSparse_to_ColSparse(const Matrix<T, Prop1,
				      Storage, Alloc1>& A,
				      Matrix<T, Prop2, ColSparse, Alloc2>& B);


  //! Conversion from RowSymSparse to column-sparse.
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, RowSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B);


  //! Conversion from ColSymSparse to column-sparse.
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B);


  /*
    From ArraySparse matrices to CSR matrices.
  */


  //! Conversion from ArrayRowSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr);


  //! Conversion from ArrayRowSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csr);


  //! Conversion from ArrayRowSymSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr);


  //! Conversion from ArrayRowComplexSparse to RowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr);


  //! Conversion from ArrayRowSymComplexSparse to RowSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0,
       ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr);


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B);


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B);


  //! Conversion from ArrayColSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csc);


  //! Conversion from ArrayRowSparse to ArrayColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& A,
	    Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B);


} // namespace Seldon.


#define SELDON_FILE_MATRIX_CONVERSIONS_HXX
#endif

#include "SeldonHeader.hxx"
#include "matrix_sparse/Matrix_ArraySparse.hxx"
#include "matrix_sparse/complex/Matrix_ComplexSparse.hxx"
#include "matrix_sparse/complex/Matrix_SymComplexSparse.hxx"
#include "matrix_sparse/complex/Matrix_ArrayComplexSparse.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/basic_functions/Functions_Vector.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#endif

namespace Seldon
{
  SELDON_EXTERN template void 
  CheckDim(const Vector<double, VectFull>&,
	   const Vector<double, VectFull>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const Vector<complex<double>, VectFull>&,
	   const Vector<complex<double>, VectFull>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ColMajor>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, RowMajor>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, RowSparse>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, RowSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ArrayRowSparse>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ArrayRowSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, RowSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, RowSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ArrayRowSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ArrayRowSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ColSparse>&,
			 const Vector<double>&, const Vector<double>&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ColSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ArrayColSparse>&,
	   const Vector<double>&, const Vector<double>&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, General, ArrayColSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ColSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ColSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ArrayColSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string);
  
  SELDON_EXTERN template void 
  CheckDim(const Matrix<complex<double>, Symmetric, ArrayColSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, RowComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ArrayRowComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, RowSymComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ArrayRowSymComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ColComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, General, ArrayColComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ColSymComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const Matrix<double, Symmetric, ArrayColSymComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, RowMajor>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, RowMajor>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ColMajor>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ColMajor>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, RowSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, RowSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ArrayRowSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ArrayRowSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, RowSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, Symmetric, RowSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, ArrayRowSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, Symmetric, ArrayRowSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ColSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ColSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ArrayColSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, General, ArrayColSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&, const Matrix<double, Symmetric, ColSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, Symmetric, ColSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, ArrayColSymSparse>&,
	   const Vector<double>&, const Vector<double>&, string, string);
  
  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<complex<double>, Symmetric, ArrayColSymSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, RowComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ArrayRowComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, RowSymComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, ArrayRowSymComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ColComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, General, ArrayColComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, ColSymComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);

  SELDON_EXTERN template void 
  CheckDim(const SeldonTranspose&,
	   const Matrix<double, Symmetric, ArrayColSymComplexSparse>&,
	   const Vector<complex<double> >&, const Vector<complex<double> >&, string, string);
  
}



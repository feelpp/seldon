#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#include "SeldonHeader.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/Lapack_Eigenvalues.cxx"
#endif

#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;

namespace Seldon
{
  SELDON_EXTERN template void GetEigenvalues(Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@complex, General, @storage_blasGE>&, Vector<@complex>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, Matrix<@real, General, @storage_blasGE>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@complex, General, @storage_blasGE>&, Vector<@complex>&, Matrix<@complex, General, @storage_blasGE>&, LapackInfo& info);

}


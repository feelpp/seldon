#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#include "SeldonHeader.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "computation/interfaces/Blas_1.cxx"
#include "computation/basic_functions/Functions_Vector.cxx"
#endif

#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;

namespace Seldon
{
  SELDON_EXTERN template void ApplyRot(Vector<@real>&, Vector<@real>&, const @real, const @real);
  SELDON_EXTERN template void ApplyModifRot(Vector<@real>&, Vector<@real>&, const @real*);
  SELDON_EXTERN template void Swap(Vector<@real_complex>&, Vector<@real_complex>&);
  SELDON_EXTERN template void Mlt(const @scalar, Vector<@scalar>&);
  SELDON_EXTERN template void Copy(const Vector<@real_complex>&, Vector<@real_complex>&);
  SELDON_EXTERN template void Add(const @real_complex, const Vector<@real_complex>&, Vector<@real_complex>&);
  SELDON_EXTERN template @real_complex DotProd(const Vector<@real_complex>&, const Vector<@real_complex>&);
  SELDON_EXTERN template @complex DotProdConj(const Vector<@complex>&, const Vector<@complex>&);
  SELDON_EXTERN template float Norm1(const Vector<float>&);
  SELDON_EXTERN template double Norm1(const Vector<double>&);
  SELDON_EXTERN template float Norm1(const Vector<complexfloat>&);
  SELDON_EXTERN template double Norm1(const Vector<complexdouble>&);
  
  SELDON_EXTERN template float Norm2(const Vector<float>&);
  SELDON_EXTERN template double Norm2(const Vector<double>&);
  SELDON_EXTERN template float Norm2(const Vector<complexfloat>&);
  SELDON_EXTERN template double Norm2(const Vector<complexdouble>&);

  SELDON_EXTERN template size_t GetMaxAbsIndex(const Vector<@real_complex>&);
}


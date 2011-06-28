#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#include "SeldonHeader.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "vector/Vector.cxx"
#include "matrix/Matrix_Base.cxx"
#include "matrix/Matrix_Pointers.cxx"
#include "matrix/Matrix_Symmetric.cxx"
#include "matrix/Matrix_SymPacked.cxx"
#include "matrix/Matrix_Hermitian.cxx"
#include "matrix/Matrix_HermPacked.cxx"
#include "share/Allocator.cxx"
#include "share/Storage.cxx"
#include "share/MatrixFlagInline.cxx"
#include "computation/interfaces/Lapack_LinearEquations.cxx"
#include "computation/interfaces/Lapack_LeastSquares.cxx"
#include "computation/interfaces/Lapack_Eigenvalues.cxx"
#endif

#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;

namespace Seldon
{
  /* Linear equations */
     
  SELDON_EXTERN template void GetLU(Matrix<@real_complex, General, @storage_blasGE>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<@real_complex, Symmetric, @storage_blasS>&, Vector<int>&, LapackInfo&);
  SELDON_EXTERN template void GetLU(Matrix<@complex, Hermitian, @storage_blasH>&, Vector<int>&, LapackInfo&);
  
  SELDON_EXTERN template void SolveLU(const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<int>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, const Matrix<@real_complex, General, @storage_blasGE>&, const Vector<int>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const Matrix<@real_complex, Symmetric, @storage_blasS>&, const Vector<int>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const Matrix<@complex, Hermitian, @storage_blasH>&, const Vector<int>&, Vector<@complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const Matrix<@real_complex, General, @storage_blasT>&, Vector<@real_complex>&, LapackInfo&);
  SELDON_EXTERN template void SolveLU(const SeldonTranspose&, const SeldonDiag&, const Matrix<@real_complex, General, @storage_blasT>&, Vector<@real_complex>&, LapackInfo&);
  
  SELDON_EXTERN template @real ReciprocalConditionNumber(const Matrix<@real, General, @storage_blasGE>&, Vector<int>&, SeldonNorm, @real, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const Matrix<complexdouble, General, @storage_blasGE>&, Vector<int>&, SeldonNorm, double, LapackInfo&);
  SELDON_EXTERN template @real ReciprocalConditionNumber(const Matrix<@real, Symmetric, @storage_blasS>&, Vector<int>&, SeldonNorm, @real, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const Matrix<complexdouble, Symmetric, @storage_blasS>&, Vector<int>&, SeldonNorm, double, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const Matrix<complexdouble, Hermitian, @storage_blasH>&, Vector<int>&, SeldonNorm, double, LapackInfo&);
  SELDON_EXTERN template @real ReciprocalConditionNumber(const Matrix<@real, General, @storage_blasT>&, SeldonNorm, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const Matrix<complexdouble, General, @storage_blasT>&, SeldonNorm, LapackInfo&);
  SELDON_EXTERN template @real ReciprocalConditionNumber(const SeldonDiag&, const Matrix<@real, General, @storage_blasT>&, SeldonNorm, LapackInfo&);
  SELDON_EXTERN template double ReciprocalConditionNumber(const SeldonDiag&, const Matrix<complexdouble, General, @storage_blasT>&, SeldonNorm, LapackInfo&);
  
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<@real, General, @storage_blasGE>&, const Matrix<@real, General, @storage_blasGE>&, const Vector<int>&, Vector<@real>&, const Vector<@real>&, @real&, @real&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, General, RowMajor>&, const Matrix<complexdouble, General, RowMajor>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, General, ColMajor>&, const Matrix<complexdouble, General, ColMajor>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const Matrix<@real, General, @storage_blasGE>&, const Matrix<@real, General, @storage_blasGE>&, const Vector<int>&, Vector<@real>&, const Vector<@real>&, @real&, @real&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const Matrix<complexdouble, General, ColMajor>&, const Matrix<complexdouble, General, ColMajor>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const SeldonTranspose&, const Matrix<complexdouble, General, RowMajor>&, const Matrix<complexdouble, General, RowMajor>&, const Vector<int>&, Vector<complexdouble>&, Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<@real, Symmetric, @storage_blasS>&, const Matrix<@real, Symmetric, @storage_blasS>&, const Vector<int>&, Vector<@real>&, const Vector<@real>&, @real&, @real&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Symmetric, @storage_blasS>&, const Matrix<complexdouble, Symmetric, @storage_blasS>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Hermitian, RowHerm>&, const Matrix<complexdouble, Hermitian, RowHerm>&, const Vector<int>&, Vector<complexdouble>&, Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Hermitian, RowHermPacked>&, const Matrix<complexdouble, Hermitian, RowHermPacked>&, const Vector<int>&, Vector<complexdouble>&, Vector<complexdouble>&, double&, double&, LapackInfo&);
  SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Hermitian, ColHerm>&, const Matrix<complexdouble, Hermitian, ColHerm>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
SELDON_EXTERN template void RefineSolutionLU(const Matrix<complexdouble, Hermitian, ColHermPacked>&, const Matrix<complexdouble, Hermitian, ColHermPacked>&, const Vector<int>&, Vector<complexdouble>&, const Vector<complexdouble>&, double&, double&, LapackInfo&);
  
  SELDON_EXTERN template void GetInverse(Matrix<@real_complex, General, @storage_blasGE>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<@real_complex, Symmetric, @storage_blasS>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<@complex, Hermitian, @storage_blasH>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(Matrix<@real_complex, General, @storage_blasT>&, LapackInfo&);
  SELDON_EXTERN template void GetInverse(const SeldonDiag&, Matrix<@real_complex, General, @storage_blasT>&, LapackInfo&);
  
  SELDON_EXTERN template void GetScalingFactors(const Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, @real&, @real&, @real&, LapackInfo&);
  SELDON_EXTERN template void GetScalingFactors(const Matrix<complexdouble, General, @storage_blasGE>&, Vector<double>&, Vector<double>&, double&, double&, double&, LapackInfo&);
  
  SELDON_EXTERN template void GetCholesky(Matrix<double, Symmetric, @storage_blasS>&, LapackInfo&);
  SELDON_EXTERN template void GetCholesky(Matrix<complexdouble, Hermitian, @storage_blasH>&, LapackInfo&);
  
  SELDON_EXTERN template void SolveCholesky(const @trans&, const Matrix<double, Symmetric, @storage_blasS>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void SolveCholesky(const @trans&, const Matrix<complexdouble, Hermitian, @storage_blasH>&, Vector<complexdouble>&, LapackInfo&);
  
  SELDON_EXTERN template void MltCholesky(const @trans&, const Matrix<double, Symmetric, @storage_blasS>&, Vector<double>&, LapackInfo&);
  SELDON_EXTERN template void MltCholesky(const @trans&, const Matrix<complexdouble, Hermitian, @storage_blasH>&, Vector<complexdouble>&, LapackInfo&);
  
  /* Least-squares */
  
  

  /* Eigenvalue problems */
  
  SELDON_EXTERN template void GetEigenvalues(Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvalues(Matrix<@complex, General, @storage_blasGE>&, Vector<@complex>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@real, General, @storage_blasGE>&, Vector<@real>&, Vector<@real>&, Matrix<@real, General, @storage_blasGE>&, LapackInfo& info);
  SELDON_EXTERN template void GetEigenvaluesEigenvectors(Matrix<@complex, General, @storage_blasGE>&, Vector<@complex>&, Matrix<@complex, General, @storage_blasGE>&, LapackInfo& info);

}


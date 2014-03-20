#ifndef SELDON_FILE_MKL_SPARSE_HXX

extern "C"
{
  float  cblas_sdoti(const int N, const float *X, const int *indx,
		     const float *Y);
  
  double cblas_ddoti(const int N, const double *X, const int *indx,
		     const double *Y);
  
  void cblas_cdotui_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);
  
  void cblas_cdotci_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);

  void cblas_zdotui_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);

  void cblas_zdotci_sub(const int N, const void *X, const int *indx,
                        const void *Y, void *dotui);
  
  void cblas_saxpyi(const int N, const float alpha, const float *X,
		    const int *indx, float *Y);

  void cblas_daxpyi(const int N, const double alpha, const double *X,
		    const int *indx, double *Y);
  
  void cblas_caxpyi(const int N, const void *alpha, const void *X,
		    const int *indx, void *Y);
  
  void cblas_zaxpyi(const int N, const void *alpha, const void *X,
		    const int *indx, void *Y);  
  
  void cblas_sgthr(const int N, const float* y, float* x, const int* indx);

  void cblas_dgthr(const int N, const double* y, double* x, const int* indx);

  void cblas_cgthr(const int N, const void* y, void* x, const int* indx);

  void cblas_zgthr(const int N, const void* y, void* x, const int* indx);

  void cblas_sgthrz(const int N, float* y, float* x, const int* indx);

  void cblas_dgthrz(const int N, double* y, double* x, const int* indx);

  void cblas_cgthrz(const int N, void* y, void* x, const int* indx);

  void cblas_zgthrz(const int N, void* y, void* x, const int* indx);
  
  void cblas_sroti(const int N, float* x, const int* indx, float* y,
		   const float c, const float s);

  void cblas_droti(const int N, double* x, const int* indx, double* y,
		   const double c, const double s);
  
  void cblas_ssctr(const int N, float* x, const int* indx, float* y);
  
  void cblas_dsctr(const int N, double* x, const int* indx, double* y);
  
  void cblas_csctr(const int N, void* x, const int* indx, void* y);
  
  void cblas_zsctr(const int N, void* x, const int* indx, void* y);
  
}

namespace Seldon
{
  
  /***********************
   * Sparse Blas Level 1 *
   ***********************/
  
  
  template<class Alloc1, class Alloc2>
  void Add(const float& alpha, Vector<float, VectSparse, Alloc1>& X,
	   Vector<float, VectFull, Alloc2>& Y);

  template<class Alloc1, class Alloc2>
  void Add(const double& alpha, Vector<double, VectSparse, Alloc1>& X,
	   Vector<double, VectFull, Alloc2>& Y);
  
  template<class Alloc1, class Alloc2>
  void Add(const complex<float>& alpha, Vector<complex<float>, VectSparse, Alloc1>& X,
	   Vector<complex<float>, VectFull, Alloc2>& Y);

  template<class Alloc1, class Alloc2>
  void Add(const complex<double>& alpha, Vector<complex<double>, VectSparse, Alloc1>& X,
	   Vector<complex<double>, VectFull, Alloc2>& Y);

  template<class Alloc1, class Alloc2>
  float DotProd(const Vector<float, VectSparse, Alloc1>& x,
		const Vector<float, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  double DotProd(const Vector<double, VectSparse, Alloc1>& x,
		 const Vector<double, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  complex<float> DotProd(const Vector<complex<float>, VectSparse, Alloc1>& x,
			 const Vector<complex<float>, VectFull, Alloc2>& y);
  
  template<class Alloc1, class Alloc2>
  complex<double> DotProd(const Vector<complex<double>, VectSparse, Alloc1>& x,
			  const Vector<complex<double>, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  float DotProdConj(const Vector<float, VectSparse, Alloc1>& x,
		    const Vector<float, VectFull, Alloc2>& y);
  
  template<class Alloc1, class Alloc2>
  double DotProdConj(const Vector<double, VectSparse, Alloc1>& x,
		     const Vector<double, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  complex<float> DotProdConj(const Vector<complex<float>, VectSparse, Alloc1>& x,
			     const Vector<complex<float>, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  complex<double> DotProdConj(const Vector<complex<double>, VectSparse, Alloc1>& x,
			      const Vector<complex<double>, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<float, VectFull, Alloc1>& y,
			 Vector<float, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<double, VectFull, Alloc1>& y,
			 Vector<double, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<complex<float>, VectFull, Alloc1>& y,
			 Vector<complex<float>, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<complex<double>, VectFull, Alloc1>& y,
			 Vector<complex<double>, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<float, VectFull, Alloc1>& y,
			     Vector<float, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<double, VectFull, Alloc1>& y,
			     Vector<double, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<complex<float>, VectFull, Alloc1>& y,
			     Vector<complex<float>, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<complex<double>, VectFull, Alloc1>& y,
			     Vector<complex<double>, VectSparse, Alloc2>& x);

  template<class Alloc1, class Alloc2>
  void ApplyRot(Vector<float, VectSparse, Alloc1>& x,
		Vector<float, VectFull, Alloc2>& y,
		const float& c, const float& s);

  template<class Alloc1, class Alloc2>
  void ApplyRot(Vector<double, VectSparse, Alloc1>& x,
		Vector<double, VectFull, Alloc2>& y,
		const double& c, const double& s);

  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<float, VectSparse, Alloc1>& x,
			  Vector<float, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<double, VectSparse, Alloc1>& x,
			  Vector<double, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<complex<float>, VectSparse, Alloc1>& x,
			  Vector<complex<float>, VectFull, Alloc2>& y);

  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<complex<double>, VectSparse, Alloc1>& x,
			  Vector<complex<double>, VectFull, Alloc2>& y);
  
} // end namespace

#define SELDON_FILE_MKL_SPARSE_HXX
#endif

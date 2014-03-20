#ifndef SELDON_FILE_MKL_SPARSE_CXX

#include "Mkl_Sparse.hxx"

namespace Seldon
{

  //! computes Y = Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const float& alpha, Vector<float, VectSparse, Alloc1>& X,
	   Vector<float, VectFull, Alloc2>& Y)
  {
    cblas_saxpyi(X.GetM(), alpha, X.GetData(), X.GetIndex(), Y.GetData());
  }

  
  //! computes Y = Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const double& alpha, Vector<double, VectSparse, Alloc1>& X,
	   Vector<double, VectFull, Alloc2>& Y)
  {
    cblas_daxpyi(X.GetM(), alpha, X.GetData(), X.GetIndex(), Y.GetData());
  }


  //! computes Y = Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const complex<float>& alpha, Vector<complex<float>, VectSparse, Alloc1>& X,
	   Vector<complex<float>, VectFull, Alloc2>& Y)
  {
    cblas_caxpyi(X.GetM(), &alpha, X.GetDataVoid(), X.GetIndex(), Y.GetDataVoid());
  }


  //! computes Y = Y + alpha X
  template<class Alloc1, class Alloc2>
  void Add(const complex<double>& alpha, Vector<complex<double>, VectSparse, Alloc1>& X,
	   Vector<complex<double>, VectFull, Alloc2>& Y)
  {
    cblas_zaxpyi(X.GetM(), &alpha, X.GetDataVoid(), X.GetIndex(), Y.GetDataVoid());
  }

  
  //! returns x^T y
  template<class Alloc1, class Alloc2>
  float DotProd(const Vector<float, VectSparse, Alloc1>& x,
		const Vector<float, VectFull, Alloc2>& y)
  {
    return cblas_sdoti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  double DotProd(const Vector<double, VectSparse, Alloc1>& x,
		 const Vector<double, VectFull, Alloc2>& y)
  {
    return cblas_ddoti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  complex<float> DotProd(const Vector<complex<float>, VectSparse, Alloc1>& x,
			 const Vector<complex<float>, VectFull, Alloc2>& y)
  {
    complex<float> scal;
    cblas_cdotui_sub(x.GetM(), x.GetDataVoid(), x.GetIndex(),
		     y.GetDataVoid(), reinterpret_cast<void*>(&scal));
    
    return scal;
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  complex<double> DotProd(const Vector<complex<double>, VectSparse, Alloc1>& x,
			 const Vector<complex<double>, VectFull, Alloc2>& y)
  {
    complex<double> scal;
    cblas_zdotui_sub(x.GetM(), x.GetDataVoid(), x.GetIndex(),
		     y.GetDataVoid(), reinterpret_cast<void*>(&scal));
    
    return scal;
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  float DotProdConj(const Vector<float, VectSparse, Alloc1>& x,
		    const Vector<float, VectFull, Alloc2>& y)
  {
    return cblas_sdoti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! returns x^T y
  template<class Alloc1, class Alloc2>
  double DotProdConj(const Vector<double, VectSparse, Alloc1>& x,
		     const Vector<double, VectFull, Alloc2>& y)
  {
    return cblas_ddoti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! returns x^H y
  template<class Alloc1, class Alloc2>
  complex<float> DotProdConj(const Vector<complex<float>, VectSparse, Alloc1>& x,
			     const Vector<complex<float>, VectFull, Alloc2>& y)
  {
    complex<float> scal;
    cblas_cdotci_sub(x.GetM(), x.GetDataVoid(), x.GetIndex(),
		     y.GetDataVoid(), reinterpret_cast<void*>(&scal));
    
    return scal;
  }


  //! returns x^H y
  template<class Alloc1, class Alloc2>
  complex<double> DotProdConj(const Vector<complex<double>, VectSparse, Alloc1>& x,
			      const Vector<complex<double>, VectFull, Alloc2>& y)
  {
    complex<double> scal;
    cblas_zdotci_sub(x.GetM(), x.GetDataVoid(), x.GetIndex(),
		     y.GetDataVoid(), reinterpret_cast<void*>(&scal));
    
    return scal;
  }

  
  //! x = y(index)
  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<float, VectFull, Alloc1>& y,
			 Vector<float, VectSparse, Alloc2>& x)
  {
    cblas_sgthr(x.GetM(), y.GetData(), x.GetData(), x.GetIndex());
  }


  //! x = y(index)
  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<double, VectFull, Alloc1>& y,
			 Vector<double, VectSparse, Alloc2>& x)
  {
    cblas_dgthr(x.GetM(), y.GetData(), x.GetData(), x.GetIndex());
  }


  //! x = y(index)
  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<complex<float>, VectFull, Alloc1>& y,
			 Vector<complex<float>, VectSparse, Alloc2>& x)
  {
    cblas_cgthr(x.GetM(), y.GetDataVoid(), x.GetDataVoid(), x.GetIndex());
  }


  //! x = y(index)
  template<class Alloc1, class Alloc2>
  void GatherSparseEntry(const Vector<complex<double>, VectFull, Alloc1>& y,
			 Vector<complex<double>, VectSparse, Alloc2>& x)
  {
    cblas_zgthr(x.GetM(), y.GetDataVoid(), x.GetDataVoid(), x.GetIndex());
  }


  //! x = y(index) and y(index) = 0
  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<float, VectFull, Alloc1>& y,
			 Vector<float, VectSparse, Alloc2>& x)
  {
    cblas_sgthrz(x.GetM(), y.GetData(), x.GetData(), x.GetIndex());
  }


  //! x = y(index) and y(index) = 0
  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<double, VectFull, Alloc1>& y,
			     Vector<double, VectSparse, Alloc2>& x)
  {
    cblas_dgthrz(x.GetM(), y.GetData(), x.GetData(), x.GetIndex());
  }


  //! x = y(index) and y(index) = 0
  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<complex<float>, VectFull, Alloc1>& y,
			     Vector<complex<float>, VectSparse, Alloc2>& x)
  {
    cblas_cgthrz(x.GetM(), y.GetDataVoid(), x.GetDataVoid(), x.GetIndex());
  }


  //! x = y(index) and y(index) = 0
  template<class Alloc1, class Alloc2>
  void GatherSparseEntryZero(const Vector<complex<double>, VectFull, Alloc1>& y,
			     Vector<complex<double>, VectSparse, Alloc2>& x)
  {
    cblas_zgthrz(x.GetM(), y.GetDataVoid(), x.GetDataVoid(), x.GetIndex());
  }

  
  //! applies rotation to vectors x and y
  template<class Alloc1, class Alloc2>
  void ApplyRot(Vector<float, VectSparse, Alloc1>& x,
		Vector<float, VectFull, Alloc2>& y,
		const float& c, const float& s)
  {
    cblas_sroti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData(), c, s);
  }


  //! applies rotation to vectors x and y
  template<class Alloc1, class Alloc2>
  void ApplyRot(Vector<double, VectSparse, Alloc1>& x,
		Vector<double, VectFull, Alloc2>& y,
		const double& c, const double& s)
  {
    cblas_droti(x.GetM(), x.GetData(), x.GetIndex(), y.GetData(), c, s);
  }
  
  
  //! y(index) = x
  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<float, VectSparse, Alloc1>& x,
			  Vector<float, VectFull, Alloc2>& y)
  {
    cblas_ssctr(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! y(index) = x
  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<double, VectSparse, Alloc1>& x,
			  Vector<double, VectFull, Alloc2>& y)
  {
    cblas_dsctr(x.GetM(), x.GetData(), x.GetIndex(), y.GetData());
  }


  //! y(index) = x
  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<complex<float>, VectSparse, Alloc1>& x,
			  Vector<complex<float>, VectFull, Alloc2>& y)
  {
    cblas_csctr(x.GetM(), x.GetDataVoid(), x.GetIndex(), y.GetDataVoid());
  }


  //! y(index) = x
  template<class Alloc1, class Alloc2>
  void ScatterSparseEntry(const Vector<complex<double>, VectSparse, Alloc1>& x,
			  Vector<complex<double>, VectFull, Alloc2>& y)
  {
    cblas_zsctr(x.GetM(), x.GetDataVoid(), x.GetIndex(), y.GetDataVoid());
  }
  
} // end namespace

#define SELDON_FILE_MKL_SPARSE_CXX
#endif

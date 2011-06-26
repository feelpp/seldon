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


#ifndef SELDON_FILE_FUNCTIONS_MATVECT_HXX

namespace Seldon
{

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	   const Vector<T1, Storage1, Allocator1>& X,
	   Vector<T2, Storage2, Allocator2>& Y);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T3& alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y);
  
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T3& alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<complex<T3>, Storage3, Allocator3>& Y);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(int alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y);
  
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonNoTrans& Trans,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonTrans& Trans,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonConjTrans& Trans,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Collection, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, ColSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
            class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, ColSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Collection, Allocator2>& X,
	      const T3 beta, Vector<T4, Collection, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Collection, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<complex<T1>, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, ColSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, ColSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<complex<T1>, Prop1, ColSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<complex<T1>, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, ColSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, ColSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<complex<T1>, Prop1, ColSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Collection, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowMajorCollection, Allocator1>& M,
	      const Vector<T2, Collection, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Collection, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, ColMajorCollection, Allocator1>& M,
	      const Vector<T2, Collection, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Collection, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& Y);
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  inline void Gauss(Matrix<T0, Prop0, Storage0, Allocator0>& M,
		    Vector<T1, Storage1, Allocator1>& X);
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  inline void GaussSeidel(const Matrix<T0, Prop0, Storage0, Allocator0>& M,			  
			  Vector<T2, Storage2, Allocator2>& Y,
                          const Vector<T1, Storage1, Allocator1>& X,
			  int iter, int type_algo = 2);
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3>
  inline void SOR(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		  Vector<T2, Storage2, Allocator2>& Y,
		  const Vector<T1, Storage1, Allocator1>& X,
		  const T3& omega, int iter, int type_ssor = 2);
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3>
  inline void SOR(const class_SeldonTrans& transM,
		  const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		  Vector<T2, Storage2, Allocator2>& Y,
		  const Vector<T1, Storage1, Allocator1>& X,
		  const T3& omega, int iter, int type_ssor = 3);
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	       Vector<T1, Storage1, Allocator1>& Y);
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetAndSolveLU(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                     Vector<T1, Storage1, Allocator1>& Y);
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "");
  
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, RowMajorCollection, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Collection, Allocator2>& Y,
		string function = "");

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, ColMajorCollection, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Collection, Allocator2>& Y,
		string function = "");
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "");
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const SeldonTranspose& trans,
		const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "", string op = "M X + Y -> Y");
  
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		string function = "", string op = "M X");

}  // namespace Seldon.

#define SELDON_FUNCTIONS_MATVECT_HXX
#endif

// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_SUPERLU_HXX

#include "superlu_interface.h"

namespace Seldon
{

  void SetComplexOne(doublecomplex& one);
  
  //! class interfacing SuperLU functions
  template<class T>
  class MatrixSuperLU_Base
  {
  protected :
    //! objects of SuperLU
    SuperMatrix L, U, B;
    GlobalLU_t Glu; //!< object of SuperLU
    
#ifdef SELDON_WITH_SUPERLU_MT
    int_t nprocs;
    superlumt_options_t options; //!< options
    Gstat_t stat; //!< statistics

    double diag_pivot_thresh, drop_tol;
    yes_no_t usepr, refact;
    fact_t fact;
    
    SCPformat *Lstore;  //!< object of SuperLU
    NCPformat *Ustore;  //!< object of SuperLU

#else
    SCformat *Lstore;  //!< object of SuperLU
    NCformat *Ustore;  //!< object of SuperLU
    superlu_options_t options; //!< options
    SuperLUStat_t stat; //!< statistics
#endif

    //! permutation array
    Vector<int_t> perm_r, perm_c;

    colperm_t permc_spec; //!< ordering scheme
    int_t n; //!< number of rows
    bool display_info; //!< display information about factorization ?
    //! Error code returned by SuperLU.
    int_t info_facto;

  public :
    MatrixSuperLU_Base();
    ~MatrixSuperLU_Base();

    const Vector<int_t>& GetRowPermutation() const;
    const Vector<int_t>& GetColPermutation() const;

    void Init(int_t size, int_t& panel_size, int_t& relax);
    void SetNumberOfThreadPerNode(int p);
    
    void SelectOrdering(colperm_t type);
    void SetPermutation(const IVect&);

    void Clear();
    void HideMessages();
    void ShowMessages();

    int GetInfoFactorization() const;
  };


  //! empty matrix
  template<class T>
  class MatrixSuperLU : public MatrixSuperLU_Base<T>
  {
  };


  //! class interfacing SuperLU functions in double precision
  template<>
  class MatrixSuperLU<double> : public MatrixSuperLU_Base<double>
  {
  public:
    MatrixSuperLU() : MatrixSuperLU_Base<double>() {}
    
    int64_t GetMemorySize() const;
    
    template<class Prop, class Allocator>
    void GetLU(Matrix<double, Prop, ColSparse, Allocator>& Lmat,
               Matrix<double, Prop, ColSparse, Allocator>& Umat,
               bool permuted = true);

    template<class Prop, class Allocator>
    void GetLU(Matrix<double, Prop, RowSparse, Allocator>& Lmat,
               Matrix<double, Prop, RowSparse, Allocator>& Umat,
               bool permuted = true);

    template<class T0, class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, Prop, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Allocator2>
    void Solve(Vector<double, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
               Vector<double, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(Matrix<double, General, ColMajor, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
               Matrix<double, General, ColMajor, Allocator2>& x);
  };


  //! class interfacing SuperLU functions in complex double precision
  template<>
  class MatrixSuperLU<complex<double> >
    : public MatrixSuperLU_Base<complex<double> >
  {
  public:
    MatrixSuperLU() : MatrixSuperLU_Base<complex<double> >() {}

    int64_t GetMemorySize() const;
    
    template<class Prop, class Allocator>
    void GetLU(Matrix<complex<double>, Prop, ColSparse, Allocator>& Lmat,
               Matrix<complex<double>, Prop, ColSparse, Allocator>& Umat,
               bool permuted = true);

    template<class Prop, class Allocator>
    void GetLU(Matrix<complex<double>, Prop, RowSparse, Allocator>& Lmat,
               Matrix<complex<double>, Prop, RowSparse, Allocator>& Umat,
               bool permuted = true);

    template<class T0, class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T0, Prop,
			 Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Allocator2>
    void Solve(Vector<complex<double>, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
               Vector<complex<double>, VectFull, Allocator2>& x);

    template<class Allocator2>
    void Solve(Matrix<complex<double>, General, ColMajor, Allocator2>& x);

    template<class Allocator2>
    void Solve(const SeldonTranspose& TransA,
               Matrix<complex<double>, General, ColMajor, Allocator2>& x);

  };

  template<class T0, class Prop, class Storage, class Allocator, class T>
  void GetLU(Matrix<T0, Prop, Storage, Allocator>& A, MatrixSuperLU<T>& mat_lu,
	     bool keep_matrix = false);

  template<class T, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x);

  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x);
  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class T, class Prop, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<T>& mat_lu, Matrix<T, Prop, ColMajor, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixSuperLU<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<double>& mat_lu,
	       Vector<complex<double>, VectFull, Allocator>& x);

  template<class Allocator>
  void SolveLU(MatrixSuperLU<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);
  
  template<class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<complex<double> >& mat_lu,
	       Vector<double, VectFull, Allocator>& x);
  
}

#define SELDON_FILE_SUPERLU_HXX
#endif

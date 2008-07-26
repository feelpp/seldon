// Copyright (C) 2001-2008 Vivien Mallet
//
// This file is part of Seldon library.
// Seldon library provides matrices and vectors structures for
// linear algebra.
// 
// Seldon is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Seldon is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License (file "license") for more details.
//
// For more information, please see the Seldon home page:
//     http://spacetown.free.fr/lib/seldon/

#ifndef SELDON_FILE_SUPERLU_CXX

#include "Computation/Interfaces/Direct/SuperLU.hxx"

namespace Seldon
{
  //! default constructor
  template<class T>
  MatrixSuperLU_Base<T>::MatrixSuperLU_Base()
  {
    //   permc_spec = 0: use the natural ordering 
    //   permc_spec = 1: use minimum degree ordering on structure of A'*A
    //   permc_spec = 2: use minimum degree ordering on structure of A'+A
    //   permc_spec = 3: use approximate mininum degree column ordering
    n = 0;
    permc_spec = 2;
    ShowMessages();
  }
  
  //! destructor
  template<class T>
  MatrixSuperLU_Base<T>::~MatrixSuperLU_Base()
  {
    Clear();
  }
  
  //! same effect as a call to the destructor
  template<class T>
  void MatrixSuperLU_Base<T>::Clear()
  {
    if (n > 0)
      {
	// SuperLU objects are cleared
	Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
	n = 0;
      }
  }
  
  //! no message from SuperLU
  template<class T>
  void MatrixSuperLU_Base<T>::HideMessages()
  {
    display_info = false;
  }
  
  //! allows messages from SuperLU
  template<class T>
  void MatrixSuperLU_Base<T>::ShowMessages()
  {
    display_info = true;
  }
  
  //! factorization of matrix in double precision using SuperLU
  template<class Prop,class Storage,class Allocator>
  void MatrixSuperLU<double>::
  FactorizeMatrix(Matrix<double,Prop,Storage,Allocator> & mat,
		  bool keep_matrix)
  {
    // conversion in CSR format
    n = mat.GetN();
    Matrix<double,General,ColSparse> Acsr;
    Copy(mat, Acsr);
    if (!keep_matrix)
      mat.Clear();
    
    // we get renumbering vectors perm_r and perm_c
    int nnz = Acsr.GetDataSize();
    dCreate_CompCol_Matrix(&A, n, n, nnz, Acsr.GetData(), Acsr.GetInd(),
			   Acsr.GetPtr(), SLU_NC, SLU_D, SLU_GE);
    
    perm_r.Reallocate(n);
    perm_c.Reallocate(n);
    
    get_perm_c(permc_spec, &A, perm_c.GetData());
    
    // factorization -> no right hand side
    int nb_rhs = 0, info;
    dCreate_Dense_Matrix(&B, n, nb_rhs, NULL, n, SLU_DN, SLU_D, SLU_GE);
    
    dgssv(&A, perm_c.GetData(), perm_r.GetData(), &L, &U, &B, &info);
    
    if ((info==0)&&(display_info))
      {
	mem_usage_t mem_usage;
	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
	cout<<"No of nonzeros in factor L = "<<Lstore->nnz<<endl;
	cout<<"No of nonzeros in factor U = "<<Ustore->nnz<<endl;
	cout<<"No of nonzeros in L+U     = "<<(Lstore->nnz+Ustore->nnz)<<endl;
	int panel_size = sp_ienv(1);
	dQuerySpace(&L, &U, panel_size, &mem_usage);
	cout<<"Memory used for factorisation in Mo "
	    <<mem_usage.total_needed/(1024*1024)<<endl;
      }
    
    Acsr.Nullify();
  }
  
  //! resolution of linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(Vector<double,Vect_Full,Allocator2>& x)
  {
    char trans('N');
    int nb_rhs = 1, info;
    dCreate_Dense_Matrix(&B, x.GetM(), nb_rhs,
			 x.GetData(), x.GetM(), SLU_DN, SLU_D, SLU_GE);
    
    dgstrs(&trans, &L, &U, perm_r.GetData(), perm_c.GetData(), &B, &info);
  }
  
  //! factorization of matrix in complex double precision using SuperLU
  template<class Prop, class Storage,class Allocator>     
  void MatrixSuperLU<complex<double> >::
  FactorizeMatrix(Matrix<complex<double>,Prop,Storage,Allocator> & mat,
		  bool keep_matrix)
  {
    // conversion in CSR format
    n = mat.GetN();
    Matrix<complex<double>,General,ColSparse> Acsr;
    Copy(mat, Acsr);
    if (!keep_matrix)
      mat.Clear();
    
    // we get renumbering vectors perm_r and perm_c
    int nnz = Acsr.GetDataSize();
    zCreate_CompCol_Matrix(&A, n, n, nnz, Acsr.GetDataVoid(), Acsr.GetInd(),
			   Acsr.GetPtr(), SLU_NC, SLU_D, SLU_GE);
    perm_r.Reallocate(n);
    perm_c.Reallocate(n);
    
    get_perm_c(permc_spec, &A, perm_c.GetData());
    int nb_rhs = 0, info;
    zCreate_Dense_Matrix(&B, n, nb_rhs, NULL, n, SLU_DN, SLU_D, SLU_GE);
    zgssv(&A, perm_c.GetData(), perm_r.GetData(), &L, &U, &B, &info);
    
    if ((info==0)&&(display_info))
      {
	mem_usage_t mem_usage;
	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
	cout<<"No of nonzeros in factor L = "<<Lstore->nnz<<endl;
	cout<<"No of nonzeros in factor U = "<<Ustore->nnz<<endl;
	cout<<"No of nonzeros in L+U     = "<<(Lstore->nnz+Ustore->nnz)<<endl;
	int panel_size = sp_ienv(1);
	zQuerySpace(&L,&U,panel_size,&mem_usage);
	cout<<"Memory used for factorisation in Mo "
	    <<mem_usage.total_needed/1e6<<endl;
      }
    
    Acsr.Nullify();
  }
  
  //! resolution of linear system A x = b
  template<class Allocator2>     
  void MatrixSuperLU<complex<double> >::
  Solve(Vector<complex<double>,Vect_Full,Allocator2>& x)
  {
    char trans('N');
    int nb_rhs = 1, info;
    zCreate_Dense_Matrix(&B, x.GetM(), nb_rhs, x.GetDataVoid(),
			 x.GetM(), SLU_DN, SLU_D, SLU_GE);
    
    zgstrs (&trans, &L, &U, perm_r.GetData(), perm_c.GetData(), &B, &info);
  }
  
  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T,Prop,Storage,Allocator>& A, MatrixSuperLU<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  
  template<class T, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu, Vector<T, Vect_Full, Allocator>& x)
  {
    mat_lu.Solve(x);
  }
  
}

#define SELDON_FILE_SUPERLU_CXX
#endif

// Copyright (C) 2001-2010 Marc Duruflé
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

#ifndef SELDON_FILE_PASTIX_CXX

#include "Pastix.hxx"

namespace Seldon
{

  template<class T>
  MatrixPastix<T>::MatrixPastix()
  {
    pastix_data = NULL;
    n = 0;
    for (int i = 0; i < 64; i++)
      {
	iparm[i] = 0;
	dparm[i] = 0;
      }

    distributed = false;

    // initializing parameters
    pastix_initParam(iparm, dparm);

    iparm[IPARM_RHS_MAKING] = API_RHS_B;
    iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
    // iparm[IPARM_VERBOSE] = API_VERBOSE_NO;
  }


  template<class T>
  MatrixPastix<T>::~MatrixPastix()
  {
    Clear();
  }


  template<>
  void MatrixPastix<double>::
  CallPastix(const MPI_Comm& comm, int* colptr, int* row,
             double* val, double* b, int nrhs)
  {
    if (distributed)
      d_dpastix(&pastix_data, comm, n, colptr, row, val,
                col_num.GetData(), perm.GetData(), invp.GetData(),
                b, nrhs, iparm, dparm);
    else
      d_pastix(&pastix_data, comm, n, colptr, row, val,
               perm.GetData(), invp.GetData(), b, nrhs, iparm, dparm);
  }


  template<>
  void MatrixPastix<double>::
  CheckMatrix(const MPI_Comm& comm, int** ptr_, int** row_, double** val_)
  {
    if (distributed)
      d_pastix_checkMatrix(comm, iparm[IPARM_VERBOSE],
                           iparm[IPARM_SYM], API_YES, n, ptr_, row_,
                           val_, NULL, 1);
    else
      d_pastix_checkMatrix(comm, iparm[IPARM_VERBOSE],
                           iparm[IPARM_SYM], API_YES, n, ptr_, row_,
                           val_, NULL, 1);
  }


  template<>
  void MatrixPastix<complex<double> >::
  CheckMatrix(const MPI_Comm& comm, int** ptr_, int** row_,
              complex<double>** val_)
  {
    if (distributed)
      z_pastix_checkMatrix(comm, iparm[IPARM_VERBOSE],
                           iparm[IPARM_SYM], API_YES, n, ptr_, row_,
                           reinterpret_cast<DCOMPLEX**>(val_), NULL, 1);
    else
      z_pastix_checkMatrix(comm, iparm[IPARM_VERBOSE],
                           iparm[IPARM_SYM], API_YES, n, ptr_, row_,
                           reinterpret_cast<DCOMPLEX**>(val_), NULL, 1);
  }


  template<>
  void MatrixPastix<complex<double> >::
  CallPastix(const MPI_Comm& comm, int* colptr, int* row,
             complex<double>* val, complex<double>* b, int nrhs)
  {
    if (distributed)
      z_dpastix(&pastix_data, comm, n, colptr, row,
                reinterpret_cast<DCOMPLEX*>(val),
                col_num.GetData(), perm.GetData(), invp.GetData(),
                reinterpret_cast<DCOMPLEX*>(b), nrhs, iparm, dparm);
    else
      /* z_dpastix(&pastix_data, comm, n, colptr, row,
         reinterpret_cast<DCOMPLEX*>(val),
         col_num.GetData(), perm.GetData(), invp.GetData(),
         reinterpret_cast<DCOMPLEX*>(b), nrhs, iparm, dparm); */
      z_pastix(&pastix_data, comm, n, colptr, row,
               reinterpret_cast<DCOMPLEX*>(val),
               perm.GetData(), invp.GetData(),
               reinterpret_cast<DCOMPLEX*>(b), nrhs, iparm, dparm);
  }


  template<class T>
  void MatrixPastix<T>::Clear()
  {
    if (n > 0)
      {
	int nrhs = 1;
	iparm[IPARM_START_TASK] = API_TASK_CLEAN;
	iparm[IPARM_END_TASK] = API_TASK_CLEAN;

	CallPastix(MPI_COMM_WORLD, NULL, NULL, NULL, NULL, nrhs);

	perm.Clear(); invp.Clear(); col_num.Clear();
	n = 0; pastix_data = NULL; distributed = false;
      }
  }


  template<class T>
  void MatrixPastix<T>::HideMessages()
  {
    iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
  }


  template<class T>
  void MatrixPastix<T>::ShowMessages()
  {
    iparm[IPARM_VERBOSE] = API_VERBOSE_NO;
  }


  template<class T>
  template<class T0, class Prop, class Storage, class Allocator>
  void MatrixPastix<T>::
  FindOrdering(Matrix<T0, Prop, Storage, Allocator> & mat,
               IVect& numbers, bool keep_matrix)
  {
    // We clear the previous factorization, if any.
    Clear();

    n = mat.GetN();
    if (n <= 0)
      return;

    int nrhs = 1, nnz = 0;
    int* ptr_ = NULL; int* rows_ = NULL;
    Matrix<int, Symmetric, RowSymSparse> As;

    iparm[IPARM_SYM] = API_SYM_YES;
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    GetSymmetricPattern(mat, As);
    if (!keep_matrix)
      mat.Clear();

    ptr_ = As.GetPtr();
    // Changing to 1-index notation.
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    nnz = As.GetNonZeros();
    rows_ = As.GetInd();
    for (int i = 0; i < nnz; i++)
      rows_[i]++;

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // We get ordering only.
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_ORDERING;

    CallPastix(MPI_COMM_WORLD, ptr_, rows_, NULL, NULL, nrhs);

    numbers = perm;
  }


  template<class T> template<class Storage, class Allocator>
  void MatrixPastix<T>
  ::FactorizeMatrix(Matrix<T, General, Storage, Allocator> & mat,
                    bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    n = mat.GetN();
    if (n <= 0)
      return;

    int nrhs = 1, nnz = 0;
    int* ptr_ = NULL; int* rows_ = NULL;
    T* values_ = NULL;
    Matrix<T, General, ColSparse, MallocAlloc<T> > A;

    iparm[IPARM_SYM] = API_SYM_NO;
    iparm[IPARM_FACTORIZATION] = API_FACT_LU;

    Copy(mat, A);
    if (!keep_matrix)
      mat.Clear();

    ptr_ = A.GetPtr();
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    nnz = A.GetNonZeros();
    rows_ = A.GetInd();
    for (int i = 0; i < nnz; i++)
      rows_[i]++;

    values_ = A.GetData();

    // destruction of pointers will be performed at the end
    A.Nullify();

    // pattern is symmetrized if needed
    CheckMatrix(MPI_COMM_WORLD, &ptr_, &rows_, &values_);

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // factorization only
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;

    CallPastix(MPI_COMM_WORLD, ptr_, rows_, values_, NULL, nrhs);

    if (iparm[IPARM_VERBOSE] != API_VERBOSE_NOT)
      cout << "Factorization successful" << endl;

    // deallocation of pointers
    free(ptr_); free(rows_); free(values_);

  }


  template<class T> template<class Storage, class Allocator>
  void MatrixPastix<T>::
  FactorizeMatrix(Matrix<T, Symmetric, Storage, Allocator> & mat,
                  bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    n = mat.GetN();
    if (n <= 0)
      return;

    //col_num.Reallocate(n);
    //for (int i = 0; i < n; i++)
    //col_num(i) = i+1;

    int nrhs = 1, nnz = 0;
    int* ptr_ = NULL; int* rows_ = NULL;
    T* values_ = NULL;
    Matrix<T, Symmetric, RowSymSparse, Allocator> As;

    iparm[IPARM_SYM]           = API_SYM_YES;
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    Copy(mat, As);
    if (!keep_matrix)
      mat.Clear();

    ptr_ = As.GetPtr();
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    nnz = As.GetNonZeros();
    rows_ = As.GetInd();
    for (int i = 0; i < nnz; i++)
      rows_[i]++;

    values_ = As.GetData();

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // factorization only
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
    // iparm[IPARM_VERBOSE]          = 4;

    CallPastix(MPI_COMM_WORLD, ptr_, rows_, values_, NULL, nrhs);

  }


  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  template<class T> template<class Allocator2, class Transpose_status>
  void MatrixPastix<T>::Solve(const Transpose_status& TransA,
                              Vector<T, VectFull, Allocator2>& x)
  {
    int nrhs = 1;
    T* rhs_ = x.GetData();

    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    iparm[IPARM_END_TASK] = API_TASK_SOLVE;

    CallPastix(MPI_COMM_WORLD, NULL, NULL, NULL, rhs_, nrhs);

  }


#ifdef SELDON_WITH_MPI
  template<class T> template<class Prop, class Allocator>
  void MatrixPastix<T>
  ::FactorizeDistributedMatrix(Matrix<T, General, ColSparse, Allocator>& A,
                               const Prop& sym, const IVect& glob_number,
                               bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    n = A.GetN();
    if (n <= 0)
      return;

    distributed = true;

    iparm[IPARM_SYM] = API_SYM_YES;
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    iparm[IPARM_GRAPHDIST] = API_YES;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int* ptr_ = A.GetPtr();
    int nrhs = 1;
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    int nnz = A.GetNonZeros();
    int* rows_ = A.GetInd();
    for (int i = 0; i < nnz; i++)
      rows_[i]++;

    T* values_ = A.GetData();

    col_num.Reallocate(n);
    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();
    for (int i = 0; i < n; i++)
      col_num(i) = glob_number(i)+1;

    // factorization only
    // iparm[IPARM_VERBOSE] = API_VERBOSE_YES;
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
    // iparm[IPARM_VERBOSE] = 4;

    CallPastix(MPI_COMM_WORLD, ptr_, rows_, values_, NULL, nrhs);

  }


  template<class T> template<class Allocator2>
  void MatrixPastix<T>
  ::SolveDistributed(Vector<T, Vect_Full, Allocator2>& x,
                     const IVect& glob_num)
  {
    SolveDistributed(SeldonNoTrans, x, glob_num);
  }


  template<class T> template<class Allocator2, class Transpose_status>
  void MatrixPastix<T>::SolveDistributed(const Transpose_status& TransA,
                                         Vector<T, Vect_Full, Allocator2>& x,
                                         const IVect& glob_num)
  {
    int nrhs = 1;
    T* rhs_ = x.GetData();

    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    iparm[IPARM_END_TASK] = API_TASK_SOLVE;

    CallPastix(MPI_COMM_WORLD, NULL, NULL, NULL, rhs_, nrhs);

  }
#endif


  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T, Prop, Storage, Allocator>& A, MatrixPastix<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }


  template<class T, class Allocator>
  void SolveLU(MatrixPastix<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  template<class T, class Allocator, class Transpose_status>
  void SolveLU(const Transpose_status& TransA,
	       MatrixPastix<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }

} // end namespace

#define SELDON_FILE_PASTIX_CXX
#endif

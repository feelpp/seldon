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

  MatrixPastix::MatrixPastix()
  {
    pastix_data = NULL;
    n = 0;
    for (int i = 0; i < 64; i++)
      {
	iparm[i] = 0;
	dparm[i] = 0;
      }

    // initializing parameters
    pastix_initParam(iparm, dparm);
    iparm[IPARM_RHS_MAKING]          = API_RHS_B;
    iparm[IPARM_VERBOSE]          = API_VERBOSE_NOT;
  }


  MatrixPastix::~MatrixPastix()
  {
    Clear();
  }


  void MatrixPastix::Clear()
  {
    if (n > 0)
      {
	int nrhs = 0;
	iparm[IPARM_START_TASK] = API_TASK_CLEAN;
	iparm[IPARM_END_TASK] = API_TASK_CLEAN;
	dpastix(&pastix_data, MPI_COMM_WORLD, n,
		NULL, NULL, NULL, col_num.GetData(), perm.GetData(),
		invp.GetData(), NULL, nrhs, iparm, dparm);

	perm.Clear(); invp.Clear(); col_num.Clear();
	n = 0; pastix_data = NULL;
      }
  }


  void MatrixPastix::HideMessages()
  {
    iparm[IPARM_VERBOSE]          = API_VERBOSE_NOT;
  }


  void MatrixPastix::ShowMessages()
  {
    iparm[IPARM_VERBOSE]          = API_VERBOSE_NO;
  }


  template<class T, class Prop, class Storage, class Allocator>
  void MatrixPastix::FindOrdering(Matrix<T, Prop, Storage, Allocator> & mat,
				  IVect& numbers, bool keep_matrix)
  {

  }


  template<class T, class Storage, class Allocator>
  void MatrixPastix
  ::FactorizeMatrix(Matrix<T, General, Storage, Allocator> & mat,
                    bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    n = mat.GetN();
    int nrhs = 0, nnz = 0;
    int* ptr_ = NULL; int* rows_ = NULL;
    pastix_float_t* values_ = NULL;
    Matrix<T, General, ColSparse, Allocator> A;

    iparm[IPARM_SYM]           = API_SYM_NO;
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

    values_ = reinterpret_cast<pastix_float_t*>(A.GetData());

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // factorization only
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
    pastix(&pastix_data, MPI_COMM_WORLD, n, ptr_, rows_, values_,
	   perm.GetData(), invp.GetData(), NULL, nrhs, iparm, dparm);

  }


  template<class T, class Storage, class Allocator>
  void MatrixPastix
  ::FactorizeMatrix(Matrix<T, Symmetric, Storage, Allocator> & mat,
                    bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    n = mat.GetN();
    int nrhs = 0, nnz = 0;
    int* ptr_ = NULL; int* rows_ = NULL;
    pastix_float_t* values_ = NULL;
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

    values_ = reinterpret_cast<pastix_float_t*>(As.GetData());

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // factorization only
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
    pastix(&pastix_data, MPI_COMM_WORLD, n, ptr_, rows_, values_,
	   perm.GetData(), invp.GetData(), NULL, nrhs, iparm, dparm);

  }


  template<class T, class Allocator2>
  void MatrixPastix::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  template<class T, class Allocator2, class Transpose_status>
  void MatrixPastix::Solve(const Transpose_status& TransA,
			   Vector<T, VectFull, Allocator2>& x)
  {
    int nrhs = 1;
    pastix_float_t* rhs_ = reinterpret_cast<pastix_float_t*>(x.GetData());

    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    iparm[IPARM_END_TASK] = API_TASK_SOLVE;
    pastix(&pastix_data, MPI_COMM_WORLD, n, NULL, NULL, NULL,
	   perm.GetData(), invp.GetData(), rhs_, nrhs, iparm, dparm);

  }


#ifdef SELDON_WITH_MPI
  template<class T, class Prop, class Allocator>
  void MatrixPastix
  ::FactorizeDistributedMatrix(Matrix<T, General, ColSparse, Allocator>& A,
                               const Prop& sym, const IVect& glob_number,
                               bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    n = A.GetN();
    int* ptr_ = A.GetPtr();
    int nrhs = 0;
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    int nnz = A.GetNonZeros();
    int* rows_ = A.GetInd();
    for (int i = 0; i < nnz; i++)
      rows_[i]++;

    pastix_float_t* values_ = reinterpret_cast<pastix_float_t*>(A.GetData());

    col_num.Reallocate(n);
    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();
    for (int i = 0; i < n; i++)
      col_num(i) = glob_number(i)+1;

    int* num_loc_ = col_num.GetData();

    // check matrix
    // pastix_checkMatrix(MPI_COMM_WORLD, API_VERBOSE_NO,
    // API_SYM_YES,  API_YES,
    // n, &ptr_, &rows_, &values_, &num_loc_);

    // factorization only
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
    dpastix(&pastix_data, MPI_COMM_WORLD, n, ptr_, rows_, values_, num_loc_,
	    perm.GetData(), invp.GetData(), NULL, nrhs, iparm, dparm);

  }


  template<class T, class Allocator2>
  void MatrixPastix::SolveDistributed(Vector<T, Vect_Full, Allocator2>& x,
                                      const IVect& glob_num)
  {
    SolveDistributed(SeldonNoTrans, x, glob_num);
  }


  template<class T, class Allocator2, class Transpose_status>
  void MatrixPastix::SolveDistributed(const Transpose_status& TransA,
				      Vector<T, Vect_Full, Allocator2>& x,
                                      const IVect& glob_num)
  {
    int nrhs = 1;
    pastix_float_t* rhs_ = reinterpret_cast<pastix_float_t*>(x.GetData());

    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    iparm[IPARM_END_TASK] = API_TASK_SOLVE;
    dpastix(&pastix_data, MPI_COMM_WORLD, n, NULL, NULL, NULL,
            col_num.GetData(), perm.GetData(), invp.GetData(), rhs_, nrhs,
            iparm, dparm);
  }
#endif

  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T, Prop, Storage, Allocator>& A, MatrixPastix& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }

  template<class T, class Allocator>
  void SolveLU(MatrixPastix& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }

  template<class T, class Allocator, class Transpose_status>
  void SolveLU(const Transpose_status& TransA,
	       MatrixPastix& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }

} // end namespace

#define SELDON_FILE_PASTIX_CXX
#endif

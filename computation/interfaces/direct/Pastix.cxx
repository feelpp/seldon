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

  //! Default constructor.
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

    // Factorization of a matrix on a single processor.
    distributed = false;

    // Refinement by default.
    refine_solution = true;

    // initializing parameters
    pastix_initParam(iparm, dparm);

    iparm[IPARM_RHS_MAKING] = API_RHS_B;
    iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
    // iparm[IPARM_VERBOSE] = API_VERBOSE_NO;
  }


  //! destructor
  template<class T>
  MatrixPastix<T>::~MatrixPastix()
  {
    Clear();
  }


  //! Calling main C-function pastix.
  template<>
  void MatrixPastix<double>::
  CallPastix(const MPI_Comm& comm, pastix_int_t* colptr, pastix_int_t* row,
             double* val, double* b, pastix_int_t nrhs)
  {
    if (distributed)
      d_dpastix(&pastix_data, comm, n, colptr, row, val,
                col_num.GetData(), perm.GetData(), invp.GetData(),
                b, nrhs, iparm, dparm);
    else
      d_pastix(&pastix_data, comm, n, colptr, row, val,
               perm.GetData(), invp.GetData(), b, nrhs, iparm, dparm);
  }


  //! Calling main C-function pastix.
  template<>
  void MatrixPastix<complex<double> >::
  CallPastix(const MPI_Comm& comm, pastix_int_t* colptr, pastix_int_t* row,
             complex<double>* val, complex<double>* b, pastix_int_t nrhs)
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


  //! Clearing factorization.
  template<class T>
  void MatrixPastix<T>::Clear()
  {
    if (n > 0)
      {
	pastix_int_t nrhs = 1;
	iparm[IPARM_START_TASK] = API_TASK_CLEAN;
	iparm[IPARM_END_TASK] = API_TASK_CLEAN;

	CallPastix(MPI_COMM_WORLD, NULL, NULL, NULL, NULL, nrhs);

	perm.Clear();
        invp.Clear();
        col_num.Clear();
	n = 0;
        pastix_data = NULL;
        distributed = false;

        MPI_Comm_free(&comm_facto);
      }
  }


  //! no message will be displayed
  template<class T>
  void MatrixPastix<T>::HideMessages()
  {
    iparm[IPARM_VERBOSE] = API_VERBOSE_NOT;
  }


  //! Low level of display.
  template<class T>
  void MatrixPastix<T>::ShowMessages()
  {
    iparm[IPARM_VERBOSE] = API_VERBOSE_NO;
  }

  //! You can require that solution is refined after LU resolution.
  template<class T>
  void MatrixPastix<T>::RefineSolution()
  {
    refine_solution = true;
  }


  //! You can require that solution is not refined (faster).
  template<class T>
  void MatrixPastix<T>::DoNotRefineSolution()
  {
    refine_solution = false;
  }


  template<class T>
  void MatrixPastix<T>::CreateCommunicator()
  {
    MPI_Group single_group;
    int rank, nb_cpu;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_cpu);

    if (distributed)
      {
        IVect rank_(nb_cpu); rank_.Fill();
        MPI_Comm_group(MPI_COMM_WORLD, &single_group);
	MPI_Group_incl(single_group, nb_cpu, rank_.GetData(), &single_group);
	MPI_Comm_create(MPI_COMM_WORLD, single_group, &comm_facto);
      }
    else
      {
        MPI_Comm_group(MPI_COMM_WORLD, &single_group);
	MPI_Group_incl(single_group, 1, &rank, &single_group);
	MPI_Comm_create(MPI_COMM_WORLD, single_group, &comm_facto);
      }
  }


  //! Returning ordering found by Scotch.
  template<class T>
  template<class T0, class Prop, class Storage, class Allocator, class Tint>
  void MatrixPastix<T>::
  FindOrdering(Matrix<T0, Prop, Storage, Allocator> & mat,
               Vector<Tint>& numbers, bool keep_matrix)
  {
    // We clear the previous factorization, if any.
    Clear();

    distributed = false;
    CreateCommunicator();

    n = mat.GetN();
    if (n <= 0)
      return;

    pastix_int_t nrhs = 1, nnz = 0;
    pastix_int_t* ptr_ = NULL;
    pastix_int_t* ind_ = NULL;
    Vector<pastix_int_t> Ptr, Ind;

    iparm[IPARM_SYM] = API_SYM_YES;
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    GetSymmetricPattern(mat, Ptr, Ind);
    if (!keep_matrix)
      mat.Clear();

    ptr_ = Ptr.GetData();
    // Changing to 1-index notation.
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    nnz = Ind.GetM();
    ind_ = Ind.GetData();
    for (int i = 0; i < nnz; i++)
      ind_[i]++;

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // We get ordering only.
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_ORDERING;

    CallPastix(comm_facto, ptr_, ind_, NULL, NULL, nrhs);

    numbers = perm;
  }


  //! Factorization of unsymmetric matrix
  template<class T> template<class Storage, class Allocator>
  void MatrixPastix<T>
  ::FactorizeMatrix(Matrix<T, General, Storage, Allocator> & mat,
                    bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    distributed = false;
    CreateCommunicator();
    n = mat.GetN();
    if (n <= 0)
      return;

    pastix_int_t nrhs = 1, nnz = 0;
    pastix_int_t* ptr_ = NULL;
    pastix_int_t* ind_ = NULL;
    T* values_ = NULL;
    Vector<pastix_int_t, VectFull, MallocAlloc<pastix_int_t> > Ptr, IndRow;
    Vector<T, VectFull, MallocAlloc<T> > Val;

    iparm[IPARM_SYM] = API_SYM_NO;
    iparm[IPARM_FACTORIZATION] = API_FACT_LU;

    General prop;
    ConvertToCSC(mat, prop, Ptr, IndRow, Val, true);
    if (!keep_matrix)
      mat.Clear();

    ptr_ = Ptr.GetData();
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    nnz = IndRow.GetM();
    ind_ = IndRow.GetData();
    for (int i = 0; i < nnz; i++)
      ind_[i]++;

    values_ = Val.GetData();

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // factorization only
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;

    CallPastix(comm_facto, ptr_, ind_, values_, NULL, nrhs);

    if (iparm[IPARM_VERBOSE] != API_VERBOSE_NOT)
      cout << "Factorization successful" << endl;
  }


  //! Factorization of symmetric matrix.
  template<class T> template<class Storage, class Allocator>
  void MatrixPastix<T>::
  FactorizeMatrix(Matrix<T, Symmetric, Storage, Allocator> & mat,
                  bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    distributed = false;
    CreateCommunicator();
    n = mat.GetN();
    if (n <= 0)
      return;

    pastix_int_t nrhs = 1, nnz = 0;
    pastix_int_t* ptr_ = NULL;
    pastix_int_t* ind_ = NULL;

    T* values_ = NULL;
    Vector<pastix_int_t, VectFull, MallocAlloc<pastix_int_t> > Ptr, IndRow;
    Vector<T, VectFull, MallocAlloc<T> > Val;

    Symmetric prop;
    ConvertToCSR(mat, prop, Ptr, IndRow, Val);

    iparm[IPARM_SYM]           = API_SYM_YES;
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    if (!keep_matrix)
      mat.Clear();

    ptr_ = Ptr.GetData();
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    nnz = IndRow.GetM();
    ind_ = IndRow.GetData();
    for (int i = 0; i < nnz; i++)
      ind_[i]++;

    values_ = Val.GetData();

    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();

    // factorization only
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
    // iparm[IPARM_VERBOSE]          = 4;

    CallPastix(comm_facto, ptr_, ind_, values_, NULL, nrhs);
  }


  //! solving A x = b (A is already factorized)
  template<class T> template<class Allocator2>
  void MatrixPastix<T>::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! solving A x = b or A^T x = b (A is already factorized)
  template<class T> template<class Allocator2, class Transpose_status>
  void MatrixPastix<T>::Solve(const Transpose_status& TransA,
                              Vector<T, VectFull, Allocator2>& x)
  {
    pastix_int_t nrhs = 1;
    T* rhs_ = x.GetData();

    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    if (refine_solution)
      iparm[IPARM_END_TASK] = API_TASK_REFINE;
    else
      iparm[IPARM_END_TASK] = API_TASK_SOLVE;

    CallPastix(comm_facto, NULL, NULL, NULL, rhs_, nrhs);
  }


#ifdef SELDON_WITH_MPI

  //! Modifies the number of threads per node.
  template<class T>
  void MatrixPastix<T>::SetNbThreadPerNode(int nb_thread)
  {
    iparm[IPARM_THREAD_NBR] = nb_thread;
  }


  //! Distributed factorization (on several nodes).
  template<class T>
  template<class Alloc1, class Alloc2, class Alloc3, class Tint>
  void MatrixPastix<T>::
  FactorizeDistributedMatrix(Vector<pastix_int_t, VectFull, Alloc1>& Ptr,
                             Vector<pastix_int_t, VectFull, Alloc2>& IndRow,
                             Vector<T, VectFull, Alloc3>& Val,
                             const Vector<Tint>& glob_number,
                             bool sym, bool keep_matrix)
  {
    // we clear previous factorization if present
    Clear();

    n = Ptr.GetM() - 1;
    if (n <= 0)
      return;

    distributed = true;

    CreateCommunicator();

    iparm[IPARM_SYM] = API_SYM_YES;

    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    iparm[IPARM_GRAPHDIST] = API_YES;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    pastix_int_t* ptr_ = Ptr.GetData();
    pastix_int_t nrhs = 1;
    // changing to 1-index notation
    for (int i = 0; i <= n; i++)
      ptr_[i]++;

    pastix_int_t nnz = IndRow.GetM();
    pastix_int_t* ind_ = IndRow.GetData();
    for (int i = 0; i < nnz; i++)
      ind_[i]++;

    T* values_ = Val.GetData();

    col_num.Reallocate(n);
    perm.Reallocate(n); invp.Reallocate(n);
    perm.Fill(); invp.Fill();
    for (int i = 0; i < n; i++)
      col_num(i) = glob_number(i)+1;

    // factorization only
    // iparm[IPARM_VERBOSE] = API_VERBOSE_YES;
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK] = API_TASK_NUMFACT;

    CallPastix(comm_facto, ptr_, ind_, values_, NULL, nrhs);
  }


  template<class T> template<class Allocator2, class Tint>
  void MatrixPastix<T>::SolveDistributed(Vector<T, Vect_Full, Allocator2>& x,
                                         const Vector<Tint>& glob_num)
  {
    SolveDistributed(SeldonNoTrans, x, glob_num);
  }


  template<class T>
  template<class Allocator2, class Transpose_status, class Tint>
  void MatrixPastix<T>::SolveDistributed(const Transpose_status& TransA,
                                         Vector<T, Vect_Full, Allocator2>& x,
                                         const Vector<Tint>& glob_num)
  {
    pastix_int_t nrhs = 1;
    T* rhs_ = x.GetData();

    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    if (refine_solution)
      iparm[IPARM_END_TASK] = API_TASK_REFINE;
    else
      iparm[IPARM_END_TASK] = API_TASK_SOLVE;

    CallPastix(comm_facto, NULL, NULL, NULL, rhs_, nrhs);
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

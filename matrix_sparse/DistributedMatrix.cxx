// Copyright (C) 2014 INRIA
// Author(s): Marc Duruflé
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

#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_CXX

#include "DistributedMatrix.hxx"

namespace Seldon
{

  //! equality *this = X 
  template<class T, class Prop, class Storage, class Allocator>
  inline DistributedMatrix<T, Prop, Storage, Allocator>&
  DistributedMatrix<T, Prop, Storage, Allocator>::
  operator=(const DistributedMatrix<T, Prop, Storage, Allocator>& X)
  {
    Seldon::Copy(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(X),
		 *this);
    
    dist_col = X.dist_col;
    dist_row = X.dist_row;

    proc_col = X.proc_col;
    proc_row = X.proc_row;
    
    GlobalRowNumbers = X.GlobalRowNumbers;
    OverlapProcNumbers = X.OverlapProcNumbers;
    OverlapRowNumbers = X.OverlapRowNumbers;

    ProcSharingRows = X.ProcSharingRows;
    SharingRowNumbers = X.SharingRowNumbers;
    
    nodl_scalar_ = X.nodl_scalar_;
    nb_unknowns_scal_ = X.nb_unknowns_scal_;
    nglob_ = X.nglob_;
    
    comm_ = X.comm_;
    
    global_row_to_recv = X.global_row_to_recv;
    global_col_to_recv = X.global_col_to_recv;
    ptr_global_row_to_recv = X.ptr_global_row_to_recv;
    ptr_global_col_to_recv = X.ptr_global_col_to_recv;
    
    local_row_to_send = X.local_row_to_send;
    local_col_to_send = X.local_col_to_send;
    proc_col_to_recv = X.proc_col_to_recv;
    proc_col_to_send = X.proc_col_to_send;
    proc_row_to_recv = X.proc_row_to_recv;
    proc_row_to_send = X.proc_row_to_send;
    local_number_distant_values = X.local_number_distant_values;
    
    size_max_distant_row = X.size_max_distant_row;
    size_max_distant_col = X.size_max_distant_col;
    
    return *this;
  }
  
    
  //! initialisation of a distributed matrix with global row numbers
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Init(Vector<IVect>& all_rows,
         IVect& row_num, IVect& overlap_num, IVect& proc_num,
         int Nvol, int nb_u, IVect& MatchingProc,
         Vector<IVect>& MatchingDofNumber,
	 MPI::Comm& comm, bool distribute_row)
  {
    // throwing an error if the rows are not sorted
    if ((comm.Get_rank() == 0) && (distribute_row))
      for (int i = 0; i < comm.Get_size(); i++)
        {
          bool sorted = true;
          for (int j = 1; j < all_rows(i).GetM(); j++)
            if (all_rows(i)(j) < all_rows(i)(j-1))
              sorted = false;
          
          if (!sorted)
            {
              cout << "Row numbers must be sorted in this function" << endl;
              abort();
              //Sort(all_rows(i));
            }
        }
    
    // finding the global number of rows
    int nglob = 0;
    if (comm.Get_rank() == 0)
      for (int i = 0; i < all_rows.GetM(); i++)
        for (int j = 0; j < all_rows(i).GetM(); j++)
          if (all_rows(i)(j)+1 > nglob)
            nglob = all_rows(i)(j)+1;
    
    if (comm.Get_rank() == 0)
      {
        // sending row numbers to all processors
        for (int i = 1; i < comm.Get_size(); i++)
          {
            int n = all_rows(i).GetM();
            comm.Ssend(&nglob, 1, MPI::INTEGER, i, 101);
            if (distribute_row)
              {
                comm.Ssend(&n, 1, MPI::INTEGER, i, 102);
                comm.Ssend(all_rows(i).GetData(), all_rows(i).GetM(),
                           MPI::INTEGER, i, 103);
              }
          }
        
        if (distribute_row)
          row_num = all_rows(0);
        
        // minimal processor for each row
        Vector<int> MinProc(nglob);
        MinProc.Fill(comm.Get_size());
        for (int i = 0; i < all_rows.GetM(); i++)
          for (int j = 0; j < all_rows(i).GetM(); j++)
            if (MinProc(all_rows(i)(j)) > i)
              MinProc(all_rows(i)(j)) = i;
        
        // counting all the rows that are overlapped 
	// (i.e. shared with a processor of lower rank)
        for (int i = 0; i < all_rows.GetM(); i++)
          {
            int nb_overlap = 0;
            for (int j = 0; j < all_rows(i).GetM(); j++)
              if (MinProc(all_rows(i)(j)) < i)
                nb_overlap++;
            
            IVect num(nb_overlap), proc(nb_overlap);
            nb_overlap = 0;
            for (int j = 0; j < all_rows(i).GetM(); j++)
              if (MinProc(all_rows(i)(j)) < i)
                {
                  num(nb_overlap) = j;            
                  proc(nb_overlap) = MinProc(all_rows(i)(j));
                  nb_overlap++;
                }
            
            // sending overlapped rows
            if (i == 0)
              {
                overlap_num = num;
                proc_num = proc;
              }
            else
              {
                comm.Ssend(&nb_overlap, 1, MPI::INTEGER, i, 104);
                if (nb_overlap > 0)
                  {
                    comm.Ssend(num.GetData(), nb_overlap,
			       MPI::INTEGER, i, 105);
                    comm.Ssend(proc.GetData(), nb_overlap,
			       MPI::INTEGER, i, 106);
                  }
              }
          }
        
        // counting the number of processors sharing each global row
        MinProc.Fill(0);
        for (int i = 0; i < all_rows.GetM(); i++)
          for (int j = 0; j < all_rows(i).GetM(); j++)
            MinProc(all_rows(i)(j))++;
        
        int nb_shared_row = 0;
        for (int i = 0; i < nglob; i++)
          {
            if (MinProc(i) > 1)
              {
                MinProc(i) = nb_shared_row;
                nb_shared_row++;
              }
            else
              MinProc(i) = -1;
          }
        
        // for each shared row, storing the list of processors
        Vector<IVect> ListProcRow(nb_shared_row);
        for (int i = 0; i < all_rows.GetM(); i++)
          for (int j = 0; j < all_rows(i).GetM(); j++)
            {
              int n = all_rows(i)(j);
              if (MinProc(n) >= 0)
                {
                  int nloc = MinProc(n);
                  ListProcRow(nloc).PushBack(i);
                }
            }
        
        // then arrays MatchingDofNumber are constructed
        for (int i = 0; i < all_rows.GetM(); i++)
          {
            // searching all the processors interacting
            Vector<bool> ProcUsed(comm.Get_size());
            ProcUsed.Fill(false);
            for (int j = 0; j < all_rows(i).GetM(); j++)
              {
                int n = all_rows(i)(j);
                if (MinProc(n) >= 0)
                  {
                    int nloc = MinProc(n);
                    for (int k = 0; k < ListProcRow(nloc).GetM(); k++)
                      ProcUsed(ListProcRow(nloc)(k)) = true;
                  }
              }
            
            int nb_proc_interac = 0;
            for (int k = 0; k < ProcUsed.GetM(); k++)
              if ((k != i) && (ProcUsed(k)))
                nb_proc_interac++;
            
            Vector<int> matching_proc(nb_proc_interac);
            Vector<Vector<int> > matching_row(nb_proc_interac);
            nb_proc_interac = 0;
            for (int k = 0; k < ProcUsed.GetM(); k++)
              if ((k != i) && (ProcUsed(k)))
                {
                  // counting rows shared with processor k
                  int nb_row_interac = 0;
                  for (int j = 0; j < all_rows(i).GetM(); j++)
                    {
                      int n = all_rows(i)(j);
                      if (MinProc(n) >= 0)
                        {
                          int nloc = MinProc(n);
                          for (int k2 = 0;
			       k2 < ListProcRow(nloc).GetM(); k2++)
                            if (ListProcRow(nloc)(k2) == k)
                              nb_row_interac++;
                        }
                    }
                  
                  // filling arrays
                  matching_proc(nb_proc_interac) = k;
                  matching_row(nb_proc_interac).Reallocate(nb_row_interac);
                  nb_row_interac = 0;
                  for (int j = 0; j < all_rows(i).GetM(); j++)
                    {
                      int n = all_rows(i)(j);
                      if (MinProc(n) >= 0)
                        {
                          int nloc = MinProc(n);
                          for (int k2 = 0;
			       k2 < ListProcRow(nloc).GetM(); k2++)
                            if (ListProcRow(nloc)(k2) == k)
                              {
                                matching_row(nb_proc_interac)
				  (nb_row_interac) = j;
                                
				nb_row_interac++;
                              }
                        }
                    }
                                   
                  nb_proc_interac++;
                }
            
            // sending arrays to the processor i
            if (i == 0)
              {
                MatchingProc = matching_proc;
                MatchingDofNumber = matching_row;
              }
            else
              {
                comm.Ssend(&nb_proc_interac, 1, MPI::INTEGER, i, 107);
                if (nb_proc_interac > 0)
                  {
                    comm.Ssend(matching_proc.GetData(), nb_proc_interac,
                               MPI::INTEGER, i, 108);
                    
                    for (int k = 0; k < nb_proc_interac; k++)
                      {
                        int nb_row = matching_row(k).GetM();
                        comm.Ssend(&nb_row, 1, MPI::INTEGER, i, 109);
                        comm.Ssend(matching_row(k).GetData(), nb_row,
                                   MPI::INTEGER, i, 110);
                      }
                  }
              }
          }
      }
    else
      {
        // receiving row numbers
        MPI::Status status; int n;
        comm.Recv(&nglob, 1, MPI::INTEGER, 0, 101, status);
        if (distribute_row)
          {
            comm.Recv(&n, 1, MPI::INTEGER, 0, 102, status);
            row_num.Reallocate(n);
            comm.Recv(row_num.GetData(), n, MPI::INTEGER, 0, 103, status);
          }
        
        // receiving overlapped numbers
        comm.Recv(&n, 1, MPI::INTEGER, 0, 104, status);
        if (n > 0)
          {
            overlap_num.Reallocate(n);
            proc_num.Reallocate(n);
            comm.Recv(overlap_num.GetData(), n, MPI::INTEGER,0, 105, status);
            comm.Recv(proc_num.GetData(), n, MPI::INTEGER, 0, 106, status);
          }
        else
          {
            overlap_num.Clear();
            proc_num.Clear();
          }
        
        // receiving local numbers of rows shared with other processors
        comm.Recv(&n, 1, MPI::INTEGER, 0, 107, status);
        if (n > 0)
          {
            MatchingProc.Reallocate(n);
            MatchingDofNumber.Reallocate(n);
            comm.Recv(MatchingProc.GetData(), n,
		      MPI::INTEGER, 0, 108, status);
            
	    for (int k = 0; k < MatchingProc.GetM(); k++)
              {
                comm.Recv(&n, 1, MPI::INTEGER, 0, 109, status);
                MatchingDofNumber(k).Reallocate(n);
                comm.Recv(MatchingDofNumber(k).GetData(), n, MPI::INTEGER,
                          0, 110, status);
              }
          }
        else
          {
            MatchingProc.Clear();
            MatchingDofNumber.Clear();
          }
      }
    
    // once the arrays are constructed, we store the pointers
    nglob_ = nglob;
    GlobalRowNumbers = &row_num;
    OverlapRowNumbers = &overlap_num;
    OverlapProcNumbers = &proc_num;
    
    nodl_scalar_ = Nvol;
    nb_unknowns_scal_ = nb_u;
    ProcSharingRows = &MatchingProc;
    SharingRowNumbers = &MatchingDofNumber;
    
    comm_ = &comm;
  }
  

  //! initialisation of a distributed matrix with global row numbers
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Init(IVect& row_num, IVect& overlap_num, IVect& proc_num,
         int Nvol, int nb_u, IVect& MatchingProc,
         Vector<IVect>& MatchingDofNumber, MPI::Comm& comm)
  {
    // throwing an error if the rows are not sorted
    bool sorted = true;
    for (int j = 1; j < row_num.GetM(); j++)
      if (row_num(j) < row_num(j-1))
        sorted = false;
    
    if (!sorted)
      {
        cout << "Row numbers must be sorted in this function" << endl;
        abort();        
        //Sort(row_num);
      }
    
    // gathering row numbers to the processor 0
    Vector<IVect> all_rows(comm.Get_size());
    if (comm.Get_rank() == 0)
      {
        all_rows(0) = row_num;
        int nodl_par; MPI::Status status;
        for (int i = 1; i < comm.Get_size(); i++)
          {
            comm.Recv(&nodl_par, 1, MPI::INTEGER, i, 13, status);
            all_rows(i).Reallocate(nodl_par);
            comm.Recv(all_rows(i).GetData(), nodl_par,
		      MPI::INTEGER, i, 14, status);
          }
      }
    else
      {
        int nodl = row_num.GetM();
        comm.Send(&nodl, 1, MPI::INTEGER, 0, 13);
        comm.Send(row_num.GetData(), nodl, MPI::INTEGER, 0, 14);
      }
    
    // then calling Init with all_rows
    Init(all_rows, row_num, overlap_num, proc_num, Nvol, nb_u,
         MatchingProc, MatchingDofNumber, comm, false);
  }
  
  
  //! multiplication by a scalar
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  DistributedMatrix<T, Prop, Storage, Allocator>&
  DistributedMatrix<T, Prop, Storage, Allocator>::operator *=(const T0& x)
  {
    Mlt(x, static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this));
    
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i) *= x;

    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i) *= x;
    
    return *this;
  }
  

  //! erases any array used for MltAdd
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::EraseArrayForMltAdd()
  {
    global_row_to_recv.Clear(); global_col_to_recv.Clear();
    ptr_global_row_to_recv.Clear(); ptr_global_col_to_recv.Clear();
    local_row_to_send.Clear(); local_col_to_send.Clear();
    
    proc_col_to_recv.Clear(); proc_col_to_send.Clear();
    proc_row_to_recv.Clear(); proc_row_to_send.Clear();
    
    local_number_distant_values = false;
    size_max_distant_row = 0;
    size_max_distant_col = 0;
  }
  

  //! erases informations for matrix-vector product
  //! and reverts dist_row/dist_col to global numbers
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::SwitchToGlobalNumbers()
  {
    if (!local_number_distant_values)
      return;
    
    // changing row numbers
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        dist_row(i).Index(j) = global_row_to_recv(dist_row(i).Index(j));
    
    // then col numbers
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        dist_col(i).Index(j) = global_col_to_recv(dist_col(i).Index(j));
    
    // erasing datas needed for matrix-vector product
    EraseArrayForMltAdd();
  }
  
  
  //! returns the size of memory used to store the matrix
  template<class T, class Prop, class Storage, class Allocator>
  int64_t DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetMemorySize() const
  {
    int64_t taille = Matrix<T, Prop, Storage, Allocator>::GetMemorySize();
    taille += sizeof(int)*(3*dist_col.GetM()+3*dist_row.GetM());
    taille += sizeof(int)*(3*proc_col.GetM()+3*proc_row.GetM());
    taille += sizeof(int)*(3*local_row_to_send.GetM()
			   +3*local_col_to_send.GetM());
    taille += sizeof(int)*(global_row_to_recv.GetM()
			   +global_col_to_recv.GetM()+
                           ptr_global_row_to_recv.GetM()
			   +ptr_global_col_to_recv.GetM()+
                           proc_col_to_recv.GetM()+proc_col_to_send.GetM()+
                           proc_row_to_recv.GetM()+proc_row_to_send.GetM());

    for (int i = 0; i < proc_row.GetM(); i++)
      taille += sizeof(int)*proc_row(i).GetM();

    for (int i = 0; i < proc_col.GetM(); i++)
      taille += sizeof(int)*proc_col(i).GetM();

    for (int i = 0; i < local_row_to_send.GetM(); i++)
      taille += sizeof(int)*local_row_to_send(i).GetM();

    for (int i = 0; i < local_col_to_send.GetM(); i++)
      taille += sizeof(int)*local_col_to_send(i).GetM();
    
    for (int i = 0; i < dist_row.GetM(); i++)
      taille += (sizeof(entry_type)+sizeof(int))*dist_row(i).GetM();

    for (int i = 0; i < dist_col.GetM(); i++)
      taille += (sizeof(entry_type)+sizeof(int))*dist_col(i).GetM();
    
    return taille;
  }
  
  
  //! returns the number of non-zero entries stored in the matrix
  /*!
    It returns the number of non-zero entries stored on the current processor.
    To obtain an "estimation" of the global number of non-zero entries,
    MPI_Reduce should be called
   */
  template<class T, class Prop, class Storage, class Allocator>
  int DistributedMatrix<T, Prop, Storage, Allocator>::GetNonZeros() const
  {
    int nnz = Matrix<T, Prop, Storage, Allocator>::GetNonZeros();
    for (int i = 0; i < dist_col.GetM(); i++)
      nnz += dist_col(i).GetM();

    for (int i = 0; i < dist_row.GetM(); i++)
      nnz += dist_row(i).GetM();
    
    return nnz;
  }
  
  
  //! returns the number of elements stored
  /*!
    It returns the number of elements (of type T)
    stored on the current processor.
    To obtain an "estimation" of the global number of elements
    MPI_Reduce should be called
  */
  template<class T, class Prop, class Storage, class Allocator>
  int DistributedMatrix<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    int size = Matrix<T, Prop, Storage, Allocator>::GetDataSize();
    
    int size_int = 0;
    for (int i = 0; i < dist_col.GetM(); i++)
      {
        size_int += 2*dist_col(i).GetM();
        size += dist_col(i).GetM();
      }

    for (int i = 0; i < dist_row.GetM(); i++)
      {
        size_int += 2*dist_row(i).GetM();
        size += dist_row(i).GetM();
      }
    
    size_int += 2*(global_row_to_recv.GetM() + global_col_to_recv.GetM());
    
    for (int i = 0; i < local_row_to_send.GetM(); i++)
      size_int += local_row_to_send(i).GetM();

    for (int i = 0; i < local_col_to_send.GetM(); i++)
      size_int += local_col_to_send(i).GetM();
    
    size_int += 2*(dist_col.GetM() + dist_row.GetM()) +
      2*(local_row_to_send.GetM() + local_col_to_send.GetM()) + 25;
    
    return size + size_int/2;
  }
  

  //! adding a non-zero entry, between a local row/column
  //! and a non-local row/column
  /*!
    This method is used only internally
    \param[in] dist_col array to which entry is added
    \param[in] proc_col processors associated with dist_col
    \param[in] jglob global row/column number
    \param[in] proc2 distant processor containing the global row/column
    \param[in] val value of the non-zero entry
  */
  template<class T, class Allocator>
  void AddDistantValue(Vector<T, VectSparse, Allocator>& dist_col,
		       IVect& proc_col,
		       int jglob, int proc2, const T& val)
  {
    int pos = 0;
    int size_row = dist_col.GetM();
    while ((pos < size_row) && (dist_col.Index(pos) < jglob))
      pos++;
    
    if ((pos < size_row) && (dist_col.Index(pos) == jglob))
      {
        // already existing entry
        dist_col.Value(pos) += val;
      }
    else
      {
        // new entry
        Vector<T> value(size_row);
        IVect index(size_row), proc(size_row);
        for (int k = 0; k < size_row; k++)
          {
            index(k) = dist_col.Index(k);
            value(k) = dist_col.Value(k);
            proc(k) = proc_col(k);
          }
        
        dist_col.Reallocate(size_row+1);
        proc_col.Reallocate(size_row+1);
        for (int k = 0; k < pos; k++)
          {
            dist_col.Index(k) = index(k);
            dist_col.Value(k) = value(k);
            proc_col(k) = proc(k);
          }
        
        dist_col.Index(pos) = jglob;
        dist_col.Value(pos) = val;
        proc_col(pos) = proc2;
        for (int k = pos+1; k <= size_row; k++)
          {
            dist_col.Index(k) = index(k-1);
            dist_col.Value(k) = value(k-1);
            proc_col(k) = proc(k-1);
          }
        
      }
  }

  
  //! removes non-zero entries contained in dist_vec below epsilon
  //! dist_proc is modified in the same way
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class TypeDist>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::RemoveSmallEntryDistant(const T0& epsilon,
                            TypeDist& dist_vec, Vector<IVect>& dist_proc)
  {
    for (int i = 0; i < dist_vec.GetM(); i++)
      {
        int nb = 0, size = dist_vec(i).GetM();
        for (int j = 0; j < size; j++)
          if (abs(dist_vec(i).Value(j)) > epsilon)
            nb++;
        
        if (nb < size)
          {
            IVect num(size), proc(size); Vector<entry_type> val(size);
            for (int j = 0; j < size; j++)
              {
                num(j) = dist_vec(i).Index(j);
                val(j) = dist_vec(i).Value(j);
                proc(j) = dist_proc(i)(j);
              }
            
            dist_vec(i).Reallocate(nb);
            dist_proc(i).Reallocate(nb);
            nb = 0;
            for (int j = 0; j < size; j++)
              if (abs(val(j)) > epsilon)
                {
                  dist_vec(i).Index(nb) = num(j);
                  dist_vec(i).Value(nb) = val(j);
                  dist_proc(i)(nb) = proc(j);
                  nb++;
                }
          }
      }
  }
  

  //! removes small entries of the matrix
  /*!
    The result of this function will be different from a sequential
    execution since values stored on different processors
    may add up and be eliminated or not
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::RemoveSmallEntry(const T0& epsilon)
  {
    Seldon::RemoveSmallEntry(static_cast<Matrix<T, Prop, Storage,
			     Allocator>& >(*this),
                             epsilon);
    
    RemoveSmallEntryDistant(epsilon, dist_col, proc_col);
    RemoveSmallEntryDistant(epsilon, dist_row, proc_row);
  }
  
  
  //! sets matrix to identity matrix
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::SetIdentity()
  {
    Matrix<T, Prop, Storage, Allocator>::SetIdentity();

    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return;
    
    // setting to 0 diagonal of overlapped rows
    const IVect& overlap = this->GetOverlapRowNumber();
    entry_type zero; SetComplexZero(zero);
    for (int i = 0; i < overlap.GetM(); i++)
      this->Set(overlap(i), overlap(i), zero);

    for (int i = 0; i < dist_col.GetM(); i++)
      {
        dist_col(i).Clear();
        proc_col(i).Clear();
      }
    
    for (int i = 0; i < dist_row.GetM(); i++)
      {
        dist_row(i).Clear();
        proc_row(i).Clear();
      }      
    
    EraseArrayForMltAdd();
  }
  
    
  //! sets values of non-zero entries to 0
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Zero()
  {
    Matrix<T, Prop, Storage, Allocator>::Zero();
    
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i).Zero();
    
    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i).Zero();    
  }
  
  
  //! sets values of non-zero entries to 0, 1, 2, etc
  /*!
    The result will be different from a sequential execution
    since values can be stored in different processors
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Fill()
  {
    int value(0);
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          SetComplexReal(value, dist_col(i).Value(j));
          value++;
        }
    
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          SetComplexReal(value, dist_row(i).Value(j));
          value++;
        }    
  }
  
  
  //! sets values of non-zero entries to x
  /*!
    The result will be different from a sequential execution
    since values can be stored in different processors
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    Matrix<T, Prop, Storage, Allocator>::Fill(x);
        
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i).Fill(x);

    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i).Fill(x);
  }
  
  
  //! sets values of non-zero entries to random values
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::FillRand()
  {
    Matrix<T, Prop, Storage, Allocator>::FillRand();
    
    for (int i = 0; i < dist_col.GetM(); i++)
      dist_col(i).FillRand();

    for (int i = 0; i < dist_row.GetM(); i++)
      dist_row(i).FillRand();
  }
  
  
  //! writes the matrix in binary format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return Matrix<T, Prop, Storage, Allocator>::Write(FileName);
    
    abort();
  }
  

  //! writes the matrix in binary format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return Matrix<T, Prop, Storage, Allocator>::Write(FileStream);
    
    abort();
  }
  

  
  //! writes the matrix in text format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::WriteText(string FileName, bool cplx) const
  {
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return Matrix<T, Prop, Storage, Allocator>::WriteText(FileName, cplx);
    
    // a different file name for each processor
    string name = GetBaseString(FileName) + "_P" + to_str(comm.Get_rank())
      + "." + GetExtension(FileName);
    
    // opening the stream
    ofstream FileStream(name.c_str());
    FileStream.precision(15);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("DistributedMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + name + "\".");
#endif
    
    // then writing datas
    WriteText(FileStream, cplx);
    
    // closing files
    FileStream.close();
  }
  

  //! writes the matrix in text format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream, bool cplx) const
  {    
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return Matrix<T, Prop, Storage, Allocator>::WriteText(FileStream, cplx);
    
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("DistributedMatrix::WriteText(ofstream& FileStream)",
		    "Stream is not ready.");
#endif
    
    // converting local part into coordinate form
    Vector<int, VectFull, CallocAlloc<int> > IndRow, IndCol;
    Vector<entry_type> Value;
    ConvertMatrix_to_Coordinates(*this, IndRow, IndCol, Value, 0, true);
    
    // extending arrays in order to contain non-local part
    int N = 0;
    for (int i = 0; i < dist_col.GetM(); i++)
      N += dist_col(i).GetM();

    for (int i = 0; i < dist_row.GetM(); i++)
      N += dist_row(i).GetM();

    int old_size = IndRow.GetM();
    IndRow.Resize(old_size+N);
    IndCol.Resize(old_size+N);
    Value.Resize(old_size+N);
    
    // filling non-local part
    const IVect& global = this->GetGlobalRowNumber();
    N = old_size;
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          IndRow(N) = global(i) + 1;
          IndCol(N) = dist_col(i).Index(j) + 1;
          if (local_number_distant_values)
            IndCol(N) = global_col_to_recv(dist_col(i).Index(j)) + 1;
          
          Value(N) = dist_col(i).Value(j);
          N++;
        }

    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          IndCol(N) = global(i) + 1;
          IndRow(N) = dist_row(i).Index(j) + 1;
          if (local_number_distant_values)
            IndRow(N) = global_row_to_recv(dist_row(i).Index(j)) + 1;
          
          Value(N) = dist_row(i).Value(j);
          N++;
        }
    
    // changing numbers of local part
    for (int i = 0; i < old_size; i++)
      {
        IndRow(i) = global(IndRow(i)) + 1;
        IndCol(i) = global(IndCol(i)) + 1;
      }

    // writing values on the stream
    WriteCoordinateMatrix(FileStream, IndRow, IndCol, Value, cplx);
  }
  
  
  //! reads the matrix in binary format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::Read(string FileName)
  {
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return Matrix<T, Prop, Storage, Allocator>::Read(FileName);

    abort();
  }
  

  //! reads the matrix in binary format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream)
  {
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return Matrix<T, Prop, Storage, Allocator>::Read(FileStream);

    abort();
  }
  

  //! reads the matrix in text format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ReadText(string FileName, bool cplx)
  {
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return Matrix<T, Prop, Storage, Allocator>::ReadText(FileName, cplx);

    abort();
  }
  

  //! reads the matrix in text format
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ReadText(istream& FileStream, bool cplx)
  {
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return Matrix<T, Prop, Storage, Allocator>::ReadText(FileStream, cplx);

    abort();
  }
  
  
  //! changes global numbers in proc_row to local numbers
  /*!
    \param[in] dist_val list of distant non-zero entries. the distant row
    numbers are contained in dist_val(i).Index(j) for each column i
    \param[in] dist_proc processor rank for each distant non-zero entry
    \param[out] glob_num list of distant row numbers, sorted by processor rank
    it will contain row numbers for processor proc_glob(0),
    then processor proc_glob(1), etc
    \param[out] ptr_glob_num beginning of index for each processor
    in array glob_num the row numbers of processor proc_glob(i)
    in array glob_num are stored 
    in glob_num(ptr_glob_num(i):ptr_glob_num(i+1))
    \param[out] proc_glob list of distant processors involved
    in array dist_num
    \param[out] local_num for each processor proc_local(i)
    local_num(i) will store the local row numbers that should be sent
    to processor proc_local(i) (since these local
    rows will be distant in distant processor proc_local(i)),
    somehow glob_num stores the distants rows that will be received,
    whereas local_num stores the local rows that will be sent
    \param[out] proc_local list of distant processors that will receive
    local row numbers
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class TypeDist>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::SortAndAssembleDistantInteractions(TypeDist& dist_val,
				       Vector<IVect>& dist_proc,
                                       IVect& glob_num,
				       IVect& ptr_glob_num, IVect& proc_glob,
                                       Vector<IVect>& local_num,
				       IVect& proc_local)
  {
    MPI::Comm& comm = *comm_;
    // counting the number of distant interactions
    int N = 0;
    for (int i = 0; i < dist_proc.GetM(); i++)
      N += dist_proc(i).GetM();
    
    // sorting distant interactions by processor number
    IVect all_proc(N), all_num(N), permut(N);
    IVect nb_inter_per_proc(comm.Get_size());
    permut.Fill(); int nb = 0;
    nb_inter_per_proc.Fill(0);
    for (int i = 0; i < dist_proc.GetM(); i++)
      for (int j = 0; j < dist_proc(i).GetM(); j++)
        {
          all_proc(nb) = dist_proc(i)(j);
          all_num(nb) = dist_val(i).Index(j);
          nb_inter_per_proc(dist_proc(i)(j))++;
          nb++;
        }
    
    Sort(all_proc, all_num, permut);
    
    // number of processors involved ?
    int nb_global_proc = 0;
    for (int i = 0; i < nb_inter_per_proc.GetM(); i++)
      if (nb_inter_per_proc(i) > 0)
        nb_global_proc++;
    
    proc_glob.Reallocate(nb_global_proc);
    
    // for each processor, sorting numbers,
    // and counting how many different numbers there are
    int offset = 0;
    int nb_glob = 0;
    IVect nb_num_per_proc(comm.Get_size());
    nb_num_per_proc.Fill(0);
    nb_global_proc = 0;
    for (int p = 0; p < nb_inter_per_proc.GetM(); p++)
      if (nb_inter_per_proc(p) > 0)
        {
          int size = nb_inter_per_proc(p);
          Sort(offset, offset+size-1, all_num, permut);
          int prec = all_num(offset);
          int nb_glob_p = 1;
          for (int j = 1; j < size; j++)
            if (all_num(offset+j) != prec)
              {
                prec = all_num(offset+j);
                nb_glob_p++;
              }
          
          nb_num_per_proc(p) = nb_glob_p;
          nb_glob += nb_glob_p;
          proc_glob(nb_global_proc) = p;
          nb_global_proc++;
          offset += size;
        }
    
    // grouping global numbers    
    all_proc.Clear(); IVect local(N);
    offset = 0; glob_num.Reallocate(nb_glob);
    ptr_glob_num.Reallocate(nb_global_proc+1);
    nb_glob = 0; nb_global_proc = 0; ptr_glob_num(0) = 0;
    for (int p = 0; p < nb_inter_per_proc.GetM(); p++)
      if (nb_inter_per_proc(p) > 0)
        {
          int size = nb_inter_per_proc(p);
          int prec = all_num(offset);
          ptr_glob_num(nb_global_proc+1) = ptr_glob_num(nb_global_proc)
	    + nb_num_per_proc(p);
          
	  glob_num(nb_glob) = prec; nb_glob++;
          int nb_glob_p = 1;
          for (int j = 0; j < size; j++)
            {
              if (all_num(offset+j) != prec)
                {
                  prec = all_num(offset+j);
                  glob_num(nb_glob) = prec;
                  nb_glob_p++; nb_glob++;
                }
              
              local(offset+j) = nb_glob-1;
            }
          
          nb_global_proc++;
          offset += size;
        }
    
    // changing numbers in dist_val
    IVect inv_permut(N);
    for (int i = 0; i < N; i++)
      inv_permut(permut(i)) = i;
    
    nb_glob = 0;
    for (int i = 0; i < dist_val.GetM(); i++)
      for (int j = 0; j < dist_val(i).GetM(); j++)
        {
          int n = inv_permut(nb_glob);
          dist_val(i).Index(j) = local(n);
          nb_glob++;
        }
    
    // exchanging nb_num_per_proc
    IVect nb_num_send(comm.Get_size());
    comm.Alltoall(nb_num_per_proc.GetData(), 1, MPI::INTEGER,
                  nb_num_send.GetData(), 1, MPI::INTEGER);
    
    // sending numbers
    Vector<MPI::Request> request_send(comm.Get_size()),
      request_recv(comm.Get_size());
    
    nb_global_proc = 0;
    for (int p = 0; p < comm.Get_size(); p++)
      if (nb_num_per_proc(p) > 0)
        {
          int size = nb_num_per_proc(p);
          request_send(p)
	    = comm.Isend(&glob_num(ptr_glob_num(nb_global_proc)), size,
			 MPI::INTEGER, p, 17);
          
          nb_global_proc++;
        }
        
    int nb_local_proc = 0;
    for (int p = 0; p < comm.Get_size(); p++)
      if (nb_num_send(p) > 0)
        nb_local_proc++;
    
    local_num.Reallocate(nb_local_proc);
    proc_local.Reallocate(nb_local_proc);
    
    // receiving numbers
    MPI::Status status; nb_local_proc = 0;
    for (int p = 0; p < comm.Get_size(); p++)
      if (nb_num_send(p) > 0)
        {
          proc_local(nb_local_proc) = p;
          local_num(nb_local_proc).Reallocate(nb_num_send(p));
          comm.Recv(local_num(nb_local_proc).GetData(), nb_num_send(p),
                    MPI::INTEGER, p, 17, status);
          
          nb_local_proc++;
        }

    for (int i = 0; i < request_send.GetM(); i++)
      if (nb_num_send(i) > 0)
        request_send(i).Wait(status);
    
    // global to local conversion
    IVect Glob_to_local(this->GetGlobalM());
    const IVect& RowNumber = this->GetGlobalRowNumber();
    Glob_to_local.Fill(-1);
    for (int i = 0; i < RowNumber.GetM(); i++)
      Glob_to_local(RowNumber(i)) = i;
    
    // replacing global numbers with local numbers
    for (int i = 0; i < local_num.GetM(); i++)
      for (int j = 0; j < local_num(i).GetM(); j++)
        local_num(i)(j) = Glob_to_local(local_num(i)(j));    
  }
  
  
  //! prepares a matrix vector product with the distributed matrix
  /*!
    This method prepares the matrix-vector, such that a call to MltAdd
    with a distributed will perform exchanges of values between values
    in order to produce the correct result.
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::PrepareMltAdd()
  {
    if (local_number_distant_values)
      return;
    
    local_number_distant_values = true;
    MPI::Comm& comm = *comm_;
    if (comm.Get_size() == 1)
      return;

    int N = 0;
    for (int i = 0; i < dist_col.GetM(); i++)
      N += dist_col(i).GetM();
    
    comm.Allreduce(&N, &size_max_distant_col, 1, MPI::INTEGER, MPI::MAX);
    
    N = 0;
    for (int i = 0; i < dist_row.GetM(); i++)
      N += dist_row(i).GetM();
    
    comm.Allreduce(&N, &size_max_distant_row, 1, MPI::INTEGER, MPI::MAX);
    
    if (size_max_distant_col > 0)
      SortAndAssembleDistantInteractions(dist_col, proc_col,
					 global_col_to_recv,
                                         ptr_global_col_to_recv,
					 proc_col_to_recv,
                                         local_col_to_send, proc_col_to_send);
    
    if (size_max_distant_row > 0)
      SortAndAssembleDistantInteractions(dist_row, proc_row,
					 global_row_to_recv,
                                         ptr_global_row_to_recv,
					 proc_row_to_recv,
                                         local_row_to_send, proc_row_to_send);
  }
  
  
  //! internal function
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2>
  void DistributedMatrix<T, Prop, Storage, Allocator>::
  ScatterValues(const Vector<T2>& X, const IVect& num_recv,
                const IVect& ptr_num_recv, const IVect& proc_recv,
                const Vector<IVect>& num_send, const IVect& proc_send,
		Vector<T2>& Xcol) const
  {
    // sending datas
    MPI::Comm& comm = *comm_;
    Vector<Vector<T2> > xsend(proc_send.GetM()), xrecv(proc_recv.GetM());
    Vector<Vector<int64_t> > xsend_tmp(proc_send.GetM()),
      xrecv_tmp(proc_recv.GetM());
    
    int tag = 30;
    Vector<MPI::Request> request_send(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xsend(i).Reallocate(nb);
        for (int j = 0; j < nb; j++)
          xsend(i)(j) = X(num_send(i)(j));
        
        request_send(i) =
          MpiIsend(comm, xsend(i), xsend_tmp(i), nb, proc_send(i), tag);
      }
    
    // receiving datas
    Vector<MPI::Request> request_recv(proc_recv.GetM());
    int N = 0;
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i); N += nb;
        xrecv(i).Reallocate(nb);
        request_recv(i) =
          MpiIrecv(comm, xrecv(i), xrecv_tmp(i), nb, proc_recv(i), tag);
      }
    
    // waiting that transfers are effective
    MPI::Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      request_send(i).Wait(status);

    for (int i = 0; i < request_recv.GetM(); i++)
      request_recv(i).Wait(status);
    
    xsend.Clear();
    // completing receives
    for (int i = 0; i < request_recv.GetM(); i++)
      MpiCompleteIrecv(xrecv(i), xrecv_tmp(i), xrecv(i).GetM());
    
    // values are stored in Xcol
    Xcol.Reallocate(N); N = 0;
    for (int i = 0; i < proc_recv.GetM(); i++)
      for (int j = 0; j < ptr_num_recv(i+1) - ptr_num_recv(i); j++)
        Xcol(N++) = xrecv(i)(j);
  }


  //! internal function
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2>
  void DistributedMatrix<T, Prop, Storage, Allocator>::
  AssembleValues(const Vector<T2>& Xcol, const IVect& num_recv,
                 const IVect& ptr_num_recv, const IVect& proc_recv,
                 const Vector<IVect>& num_send, const IVect& proc_send,
		 Vector<T2>& X) const
  {
    // sending datas
    MPI::Comm& comm = *comm_;
    Vector<Vector<T2> > xsend(proc_recv.GetM()), xrecv(proc_send.GetM());    
    Vector<Vector<int64_t> > xsend_tmp(proc_recv.GetM()),
      xrecv_tmp(proc_send.GetM());
    
    int tag = 32, N = 0;
    Vector<MPI::Request> request_send(proc_recv.GetM());
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i);
        xsend(i).Reallocate(nb);
        for (int j = 0; j < nb; j++)
          xsend(i)(j) = Xcol(N++);
        
        request_send(i) =
          MpiIsend(comm, xsend(i), xsend_tmp(i), nb, proc_recv(i), tag);
      }
    
    // receiving datas
    Vector<MPI::Request> request_recv(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xrecv(i).Reallocate(nb);
        request_recv(i) =
          MpiIrecv(comm, xrecv(i), xrecv_tmp(i), nb, proc_send(i), tag);
      }
    
    // waiting that transfers are effective
    MPI::Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      request_send(i).Wait(status);

    for (int i = 0; i < request_recv.GetM(); i++)
      request_recv(i).Wait(status);
    
    xsend.Clear();
    // completing receives
    for (int i = 0; i < request_recv.GetM(); i++)
      MpiCompleteIrecv(xrecv(i), xrecv_tmp(i), xrecv(i).GetM());
    
    // values are added to X
    for (int i = 0; i < num_send.GetM(); i++)
      for (int j = 0; j < num_send(i).GetM(); j++)
        X(num_send(i)(j)) += xrecv(i)(j);
  }

  
  //! assembles the results for each row, by taking the minimum of Yproc
  //! then the minimum of Y
  /*!
    Instead of summing values as in function AssembleRowValues,
    the minimum is taken for Yproc, and if there is equality in Yproc,
    the minimum is taken for Y
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>::
  AssembleValuesMin(const IVect& Xcol, const IVect& Xcol_proc,
                    const IVect& num_recv, const IVect& ptr_num_recv,
		    const IVect& proc_recv,
                    const Vector<IVect>& num_send, const IVect& proc_send,
                    IVect& Y, IVect& Yproc) const
  {
    // sending datas
    MPI::Comm& comm = *comm_;
    Vector<Vector<int> > xsend(proc_recv.GetM()), xrecv(proc_send.GetM());    
    int tag = 35, N = 0;
    Vector<MPI::Request> request_send(proc_recv.GetM());
    for (int i = 0; i < proc_recv.GetM(); i++)
      {
        int nb = ptr_num_recv(i+1) - ptr_num_recv(i);
        xsend(i).Reallocate(2*nb);
        for (int j = 0; j < nb; j++)
          {
            xsend(i)(j) = Xcol(N);
            xsend(i)(nb+j) = Xcol_proc(N);
            N++;
          }
        
        request_send(i) = comm.Isend(xsend(i).GetDataVoid(), 2*nb,
                                     GetMpiDataType(Xcol), proc_recv(i), tag);
      }
    
    // receiving datas
    Vector<MPI::Request> request_recv(proc_send.GetM());
    for (int i = 0; i < proc_send.GetM(); i++)
      {
        int nb = num_send(i).GetM();
        xrecv(i).Reallocate(2*nb);
        request_recv(i) = comm.Irecv(xrecv(i).GetDataVoid(), 2*nb,
                                     GetMpiDataType(Xcol), proc_send(i), tag);
      }
    
    // waiting that transfers are effective
    MPI::Status status;
    for (int i = 0; i < request_send.GetM(); i++)
      request_send(i).Wait(status);

    for (int i = 0; i < request_recv.GetM(); i++)
      request_recv(i).Wait(status);
    
    xsend.Clear();
    // values are assembled in X
    for (int i = 0; i < num_send.GetM(); i++)
      for (int j = 0; j < num_send(i).GetM(); j++)
        {
          int nb = num_send(i).GetM();
          int proc = xrecv(i)(nb+j);
          int col = xrecv(i)(j);
          if (proc < Yproc(num_send(i)(j)))
            {
              Yproc(num_send(i)(j)) = proc;
              Y(num_send(i)(j)) = col;
            }
          else if (proc == Yproc(num_send(i)(j)))
            {
              if (col < Y(num_send(i)(j)))
                Y(num_send(i)(j)) = col;
            }
        }
  }


  //! Sends/receives values of X on distant rows (similar to ScatterColValues)
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ScatterRowValues(const Vector<T2>& X, Vector<T2>& Xcol) const
  {
    ScatterValues(X, global_row_to_recv, ptr_global_row_to_recv,
		  proc_row_to_recv,
                  local_row_to_send, proc_row_to_send, Xcol);
  }

  
  //! Sends/receives values of X on distant columns
  /*!
    this method exchanges values of X between processors so that Xcol will
    contain the values of X on distant columns dist_col
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::ScatterColValues(const Vector<T2>& X, Vector<T2>& Xcol) const
  {
    ScatterValues(X, global_col_to_recv, ptr_global_col_to_recv,
		  proc_col_to_recv,
                  local_col_to_send, proc_col_to_send, Xcol);
  }


  //! Sends/receives values of Xrow on distant rows and adds them to X
  /*!
    this method exchanges values of Xrow between processors,
    and adds them to X
    so that X will sum the results of distant rows dist_row
   */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AssembleRowValues(const Vector<T2>& Xrow, Vector<T2>& X) const
  {
    AssembleValues(Xrow, global_row_to_recv, ptr_global_row_to_recv,
		   proc_row_to_recv,
                   local_row_to_send, proc_row_to_send, X);
  }


  //! Sends/receives values of Xcol on distant columns and adds them to X
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AssembleColValues(const Vector<T2>& Xcol, Vector<T2>& X) const
  {
    AssembleValues(Xcol, global_col_to_recv, ptr_global_col_to_recv,
		   proc_col_to_recv,
                   local_col_to_send, proc_col_to_send, X);
  }
  
  
  //! assembles the vector (adds values of shared rows)
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AssembleVec(Vector<T2>& X) const
  {
    AssembleVector(X, MPI::SUM, *ProcSharingRows, *SharingRowNumbers,
                   *comm_, nodl_scalar_, nb_unknowns_scal_, 14);
  }


  //! assembles the vector (by taking the minimum instead 
  //! of summing for AssembleVec
  /*!
    The minimal Xproc is search, if equality the minimal X is searched
   */
  template<class T, class Prop, class Storage, class Allocator>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::AssembleVecMin(Vector<int>& X, Vector<int>& Xproc) const
  {
    AssembleVectorMin(X, Xproc, *ProcSharingRows, *SharingRowNumbers,
                      *comm_, nodl_scalar_, nb_unknowns_scal_, 13);
  }
  
  
  //! Y = Y + alpha A X with only distant columns of A
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Storage2, class Allocator2,
           class T4, class Storage4, class Allocator4>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddCol(const class_SeldonNoTrans& Trans,
              const Vector<T2, Storage2, Allocator2>& X,
              Vector<T4, Storage4, Allocator4>& Y) const
  {
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          int jloc = dist_col(i).Index(j);
          Y(i) += dist_col(i).Value(j)*X(jloc);
        }
  }


  //! Y = Y + alpha A^T X with only distant rows of A
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Storage2, class Allocator2,
           class T4, class Storage4, class Allocator4>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddRow(const class_SeldonTrans& Trans,
              const Vector<T2, Storage2, Allocator2>& X,
              Vector<T4, Storage4, Allocator4>& Y) const
  {
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          int jloc = dist_row(i).Index(j);
          Y(i) += dist_row(i).Value(j)*X(jloc);
        }
  }


  //! Y = Y + alpha A^H X with only distant rows of A
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Storage2, class Allocator2,
           class T4, class Storage4, class Allocator4>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddRow(const class_SeldonConjTrans& Trans,
              const Vector<T2, Storage2, Allocator2>& X,
              Vector<T4, Storage4, Allocator4>& Y) const
  {
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          int jloc = dist_row(i).Index(j);
          Y(i) += conjugate(dist_row(i).Value(j))*X(jloc);
        }
  }


  //! Y = Y + alpha A X with only distant rows of A
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Storage2, class Allocator2,
           class T4, class Storage4, class Allocator4>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddRow(const class_SeldonNoTrans& Trans,
              const Vector<T2, Storage2, Allocator2>& X,
              Vector<T4, Storage4, Allocator4>& Y) const
  {
    T4 zero; SetComplexZero(zero);
    Y.Reallocate(global_row_to_recv.GetM());
    Y.Fill(zero);
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          int jrow = dist_row(i).Index(j);
          Y(jrow) += dist_row(i).Value(j)*X(i);
        }
  }


  //! Y = Y + alpha A^T X with only distant columns of A
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Storage2, class Allocator2,
           class T4, class Storage4, class Allocator4>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddCol(const class_SeldonTrans& Trans,
              const Vector<T2, Storage2, Allocator2>& X,
              Vector<T4, Storage4, Allocator4>& Y) const
  {
    T4 zero; SetComplexZero(zero);
    Y.Reallocate(global_col_to_recv.GetM());
    Y.Fill(zero);
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          int jrow = dist_col(i).Index(j);
          Y(jrow) += dist_col(i).Value(j)*X(i);
        }
  }


  //! Y = Y + alpha A^H X with only distant columns of A
  template<class T, class Prop, class Storage, class Allocator>
  template<class T2, class Storage2, class Allocator2,
           class T4, class Storage4, class Allocator4>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::MltAddCol(const class_SeldonConjTrans& Trans,
              const Vector<T2, Storage2, Allocator2>& X,
              Vector<T4, Storage4, Allocator4>& Y) const
  {
    T4 zero; SetComplexZero(zero);
    Y.Reallocate(global_col_to_recv.GetM());
    Y.Fill(zero);
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          int jrow = dist_col(i).Index(j);
          Y(jrow) += conjugate(dist_col(i).Value(j))*X(i);
        }
  }
  

  //! grouping all the local rows of the matrix into CSR form
  /*!
    Creation of a sparse matrix B, that will contain the local rows
    of the current matrix. The column numbers are global.
    Values placed on non-local rows are ignored
  */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Allocator0>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetDistributedRows(Matrix<T0, General, ArrayRowSparse,
		       Allocator0>& B) const
  {    
    int m = this->GetM();
    int n = this->GetGlobalM();
    
    Copy(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), B);
    
    // now, we are using global numbers
    // and removing lower part of the matrix if symmetric
    B.Resize(m, n);
    const IVect& RowNumber = this->GetGlobalRowNumber();
    for (int i = 0; i < m; i++)
      {
	int size_row = B.GetRowSize(i);
        size_row += dist_col(i).GetM();
	IVect index(size_row);
	Vector<T0> value(size_row);
	int nb = 0;
	int num_row = RowNumber(i);
	if (IsSymmetricMatrix(*this))
          {
            // local values
            for (int j = 0; j < B.GetRowSize(i); j++)
              {
                int num_col = RowNumber(B.Index(i, j));
                if (num_row <= num_col)
                  {
                    index(nb) = num_col;
                    value(nb) = B.Value(i, j);
                    nb++;
                  }
              }
            
            // distant values
            for (int j = 0; j < dist_col(i).GetM(); j++)
              {
                int num_col = dist_col(i).Index(j);
                if (this->local_number_distant_values)
                  num_col = global_col_to_recv(num_col);
                
                if (num_row <= num_col)
                  {
                    index(nb) = num_col;
                    value(nb) = dist_col(i).Value(j);
                    nb++;
                  }
              }
          }
        else
          {
            // local values
            for (int j = 0; j < B.GetRowSize(i); j++)
              {
                int num_col = RowNumber(B.Index(i, j));
                index(nb) = num_col;
                value(nb) = B.Value(i, j);
                nb++;
              }
            
            // distant values
            for (int j = 0; j < dist_col(i).GetM(); j++)
              {
                int num_col = dist_col(i).Index(j);
                if (this->local_number_distant_values)
                  num_col = global_col_to_recv(num_col);
                
                index(nb) = num_col;
                value(nb) = dist_col(i).Value(j);
                nb++;
              }
          }
        
        Seldon::Assemble(nb, index, value);
        
	B.ReallocateRow(i, nb);
        for (int j = 0; j < nb; j++)
          {
            B.Index(i, j) = index(j);
            B.Value(i, j) = value(j);
          }
      }
  }
  
  
  //! grouping all the local columns of the matrix into CSC form
  /*!
    Creation of a sparse matrix B, that will contain the local columns
    of the current matrix. The row numbers are global.
    Values placed on non-local columns are ignored
  */
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0, class Allocator0>
  void DistributedMatrix<T, Prop, Storage, Allocator>::
  GetDistributedColumns(Matrix<T0, General, ArrayColSparse, Allocator0>& B,
                        bool sym_pattern) const
  {
    int m = this->GetGlobalM();
    int n = this->GetN();
    
    // conversion to CSC format of local part and symmetrisation of pattern
    IVect Ptr, Ind; Vector<T0> Val;
    Prop sym;
    ConvertToCSC(*this, sym, Ptr, Ind, Val, sym_pattern);
    
    // for row numbers, we put global numbers and we add some distant entries
    // (i.e entries with local columns,
    // and null values by symmetry of local rows )
    B.Clear(); B.Reallocate(m, n);
    const IVect& RowNumber = this->GetGlobalRowNumber();
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
        size_col += dist_row(i).GetM();
        if (sym_pattern)
          size_col += dist_col(i).GetM();
        
	IVect index(size_col);
	Vector<T0> value(size_col);
	int nb = 0;
        // local values
        for (int j = Ptr(i); j < Ptr(i+1); j++)
          {
            index(nb) = RowNumber(Ind(j));
            value(nb) = Val(j);
            nb++;
          }
        
        // distant values
        for (int j = 0; j < dist_row(i).GetM(); j++)
          {
            index(nb) = dist_row(i).Index(j);
            if (this->local_number_distant_values)    
              index(nb) = global_row_to_recv(index(nb));
            
            value(nb) = dist_row(i).Value(j);
            nb++;
          }
        
        // values due to symmetrisation of pattern
        if (sym_pattern)
          for (int j = 0; j < dist_col(i).GetM(); j++)
            {
              index(nb) = dist_col(i).Index(j);
              if (this->local_number_distant_values)    
                index(nb) = global_col_to_recv(index(nb));
            
              value(nb) = 0;
              nb++;
            }
        
        Seldon::Assemble(nb, index, value);
        
	B.ReallocateColumn(i, nb);
        for (int j = 0; j < nb; j++)
          {
            B.Index(i, j) = index(j);
            B.Value(i, j) = value(j);
          }
      }
  }
   
  
  //! matrix vector product with a distributed matrix
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3& beta,
              Vector<T4, Storage4, Allocator4>& Yres, bool assemble)
  {
    const MPI::Comm& comm = M.GetCommunicator();
    bool proceed_distant_row = true, proceed_distant_col = true;
    if (comm.Get_size() == 1)
      {
        proceed_distant_row = false;
        proceed_distant_col = false;
      }
    else
      {
        if (!M.IsReadyForMltAdd())
          {
            // preparing the matrix vector product
            // this method will be called once for
	    // the first matrix-vector product
            const_cast<DistributedMatrix<T1, Prop1,
					 Storage1, Allocator1>& >(M)
              .PrepareMltAdd();
          }
        
        if (M.GetMaxDataSizeDistantCol() == 0)
          proceed_distant_col = false;
        
        if (M.GetMaxDataSizeDistantRow() == 0)
          proceed_distant_row = false;
      }
    
    T3 zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    Vector<T4, Storage4, Allocator4> Y;
    if (beta == zero)
      Y.SetData(Yres.GetM(), Yres.GetData());
    else
      Y.Reallocate(Yres.GetM());
    
    Y.Fill(zero);
        
    // scattering column values
    Vector<T2> Xcol;
    if (proceed_distant_col)
      M.ScatterColValues(X, Xcol);
    
    // local matrix
    Mlt(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(M),
	X, Y);
    
    // adding contributions of distant columns
    if (proceed_distant_col)
      M.MltAddCol(SeldonNoTrans, Xcol, Y);
    
    // contributions of distant rows
    Vector<T4> Yrow;
    if (proceed_distant_row)
      M.MltAddRow(SeldonNoTrans, X, Yrow);

    // assembling row values
    if (proceed_distant_row)
      M.AssembleRowValues(Yrow, Y);
    
    // assembling rows shared between processors
    if (assemble)
      M.AssembleVec(Y);

    if (beta == zero)
      {
        Mlt(alpha, Y);
        Y.Nullify();
      }
    else
      {
        Mlt(beta, Yres);
        Add(alpha, Y, Yres);
      }
  }
  

  //! matrix vector product with a distributed matrix
  template<class T0, class T1, class Prop1, class Storage1,
	   class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const class_SeldonNoTrans& Trans,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3& beta,
              Vector<T4, Storage4, Allocator4>& Y, bool assemble)
  {
    MltAdd(alpha, M, X, beta, Y, assemble);
  }


  //! matrix vector product with a distributed matrix (transpose)
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const class_SeldonTrans& Trans,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3& beta,
              Vector<T4, Storage4, Allocator4>& Yres, bool assemble)
  {
    const MPI::Comm& comm = M.GetCommunicator();
    bool proceed_distant_row = true, proceed_distant_col = true;
    if (comm.Get_size() == 1)
      {
        proceed_distant_row = false;
        proceed_distant_col = false;
      }
    else
      {
        if (!M.IsReadyForMltAdd())
          {
            // preparing the matrix vector product
            // this method will be called once for 
	    // the first matrix-vector product
            const_cast<DistributedMatrix<T1, Prop1,
					 Storage1, Allocator1>& >(M)
              .PrepareMltAdd();
          }
        
        if (M.GetMaxDataSizeDistantCol() == 0)
          proceed_distant_col = false;
        
        if (M.GetMaxDataSizeDistantRow() == 0)
          proceed_distant_row = false;
      }
    
    T3 zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    Vector<T4, Storage4, Allocator4> Y;
    if (beta == zero)
      Y.SetData(Yres.GetM(), Yres.GetData());
    else
      Y.Reallocate(Yres.GetM());
    
    Y.Fill(zero);
    
    // scattering row values
    Vector<T2> Xrow;
    if (proceed_distant_row)
      M.ScatterRowValues(X, Xrow);
    
    // local matrix
    Mlt(Trans, static_cast<const Matrix<T1, Prop1,
	Storage1, Allocator1>& >(M), X, Y);

    // adding contributions of distant rows
    if (proceed_distant_row)
      M.MltAddRow(SeldonTrans, Xrow, Y);
    
    // contributions of distant columns
    Vector<T4> Ycol;
    if (proceed_distant_col)
      M.MltAddCol(SeldonTrans, X, Ycol);

    // assembling row values
    if (proceed_distant_col)
      M.AssembleColValues(Ycol, Y);

    // assembling rows shared between processors
    if (assemble)
      M.AssembleVec(Y);

    if (beta == zero)
      {
        Mlt(alpha, Y);
        Y.Nullify();
      }
    else
      {
        Mlt(beta, Yres);
        Add(alpha, Y, Yres);
      }
  }


  //! matrix vector product with a distributed matrix (conjugate transpose)
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const class_SeldonConjTrans& Trans,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              const Vector<T2, Storage2, Allocator2>& X,
              const T3& beta,
              Vector<T4, Storage4, Allocator4>& Yres, bool assemble)
  {
    const MPI::Comm& comm = M.GetCommunicator();
    bool proceed_distant_row = true, proceed_distant_col = true;
    if (comm.Get_size() == 1)
      {
        proceed_distant_row = false;
        proceed_distant_col = false;
      }
    else
      {
        if (!M.IsReadyForMltAdd())
          {
            // preparing the matrix vector product
            // this method will be called once
	    // for the first matrix-vector product
            const_cast<DistributedMatrix<T1, Prop1,
					 Storage1, Allocator1>& >(M)
              .PrepareMltAdd();
          }
        
        if (M.GetMaxDataSizeDistantCol() == 0)
          proceed_distant_col = false;
        
        if (M.GetMaxDataSizeDistantRow() == 0)
          proceed_distant_row = false;
      }
    
    // scattering row values
    Vector<T2> Xrow;
    if (proceed_distant_row)
      M.ScatterRowValues(X, Xrow);

    T3 zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    Vector<T4, Storage4, Allocator4> Y;
    if (beta == zero)
      Y.SetData(Yres.GetM(), Yres.GetData());
    else
      Y.Reallocate(Yres.GetM());
    
    Y.Fill(zero);
            
    // local matrix
    Mlt(Trans, static_cast<const Matrix<T1, Prop1,
	Storage1, Allocator1>& >(M), X, Y);

    // adding contributions of distant rows
    if (proceed_distant_row)
      M.MltAddRow(SeldonConjTrans, Xrow, Y);
    
    // contributions of distant columns
    Vector<T4> Ycol;
    if (proceed_distant_col)
      M.MltAddCol(SeldonConjTrans, X, Ycol);

    // assembling row values
    if (proceed_distant_col)
      M.AssembleColValues(Ycol, Y);

    // assembling rows shared between processors
    if (assemble)
      M.AssembleVec(Y);

    if (beta == zero)
      {
        Mlt(alpha, Y);
        Y.Nullify();
      }
    else
      {
        Mlt(beta, Yres);
        Add(alpha, Y, Yres);
      }
  }
  

  //! Performs the multiplication of a matrix with a vector.
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& M,
	   const Vector<T1, Storage1, Allocator1>& X,
	   Vector<T2, Storage2, Allocator2>& Y, bool assemble)
  {
    T2 one, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAdd(one, M, X, zero, Y, assemble);
  }

  
  //! Performs the multiplication of a matrix with a vector.
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T3& alpha,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble)
  {
    T3 zero;
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAdd(alpha, M, X, zero, Y, assemble);
  }

  
  // case where alpha = real and Y is a complex vector
  // if this method is not present, Mlt with SeldonTranspose is called
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T3& alpha,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<complex<T3>, Storage3, Allocator3>& Y, bool assemble)
  {
    complex<T3> zero;
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAdd(complex<T3>(alpha, 0), M, X, zero, Y, assemble);
  }

  
  //! computes Y = alpha M X
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(int alpha,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble)
  {
    T3 zero;
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAdd(double(alpha), M, X, zero, Y, assemble);
  }

  
  //! Performs the multiplication of a matrix with a vector.
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonNoTrans& Trans,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble)
  {
    T3 one, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAdd(one, Trans, M, X, zero, Y, assemble);
  }

  
  //! computes Y = M^T X
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonTrans& Trans,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble)
  {
    T3 one, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAdd(one, Trans, M, X, zero, Y, assemble);
  }

  
  //! computes Y = M^H X
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const class_SeldonConjTrans& Trans,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y, bool assemble)
  {
    T3 one, zero;
    SetComplexOne(one);
    SetComplexZero(zero);
    Y.Fill(zero);
    MltAdd(one, Trans, M, X, zero, Y, assemble);
  }


  //! multiplication by a scalar
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1>
  void Mlt(const T0& alpha,
           DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    A *= alpha;
  }
  

  //! matrix-vector product selecting the minimal processor number
  //! and column number for each row
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          int col = A.Index(i, j);
          if (Yproc(col) < Yproc(i))
            {
              Yproc(i) = Yproc(col);
              Y(i) = Y(col);
            }
          else if (Yproc(col) == Yproc(i))
            {
              if (Y(col) < Y(i))
                Y(i) = Y(col);
              else
                Y(col) = Y(i);
            }
          else
            {
              Yproc(col) = Yproc(i);
              Y(col) = Y(i);
            }
        }
  }


  //! matrix-vector product selecting the minimal processor number
  //! and column number for each row
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          int col = A.Index(i, j);
          if (Yproc(col) < Yproc(i))
            {
              Yproc(i) = Yproc(col);
              Y(i) = Y(col);
            }
          else if (Yproc(col) == Yproc(i))
            {
              if (Y(col) < Y(i))
                Y(i) = Y(col);
              else
                Y(col) = Y(i);
            }
          else
            {
              Yproc(col) = Yproc(i);
              Y(col) = Y(i);
            }
        }
  }


#ifdef SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX
  //! matrix-vector product selecting the minimal processor number
  //! and column number for each row
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int col = A.IndexReal(i, j);
            if (Yproc(col) < Yproc(i))
              {
                Yproc(i) = Yproc(col);
                Y(i) = Y(col);
              }
            else if (Yproc(col) == Yproc(i))
              {
                if (Y(col) < Y(i))
                  Y(i) = Y(col);
                else
                  Y(col) = Y(i);
              }
            else
              {
                Yproc(col) = Yproc(i);
                Y(col) = Y(i);
              }
          }

        for (int j = 0; j < A.GetImagRowSize(i); j++)
          {
            int col = A.IndexImag(i, j);
            if (Yproc(col) < Yproc(i))
              {
                Yproc(i) = Yproc(col);
                Y(i) = Y(col);
              }
            else if (Yproc(col) == Yproc(i))
              {
                if (Y(col) < Y(i))
                  Y(i) = Y(col);
                else
                  Y(col) = Y(i);
              }
            else
              {
                Yproc(col) = Yproc(i);
                Y(col) = Y(i);
              }
          }
      }
  }


  //! matrix-vector product selecting the minimal processor number
  //! and column number for each row
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1,
	      ArrayRowSymComplexSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            int col = A.IndexReal(i, j);
            if (Yproc(col) < Yproc(i))
              {
                Yproc(i) = Yproc(col);
                Y(i) = Y(col);
              }
            else if (Yproc(col) == Yproc(i))
              {
                if (Y(col) < Y(i))
                  Y(i) = Y(col);
                else
                  Y(col) = Y(i);
              }
            else
              {
                Yproc(col) = Yproc(i);
                Y(col) = Y(i);
              }
          }

        for (int j = 0; j < A.GetImagRowSize(i); j++)
          {
            int col = A.IndexImag(i, j);
            if (Yproc(col) < Yproc(i))
              {
                Yproc(i) = Yproc(col);
                Y(i) = Y(col);
              }
            else if (Yproc(col) == Yproc(i))
              {
                if (Y(col) < Y(i))
                  Y(i) = Y(col);
                else
                  Y(col) = Y(i);
              }
            else
              {
                Yproc(col) = Yproc(i);
                Y(col) = Y(i);
              }
          }
      }
  }
#endif
  
              
  //! matrix vector product with a distributed matrix (minimal column numbers)
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void MltMin(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              IVect& Y, IVect& Yproc)
  {
    const MPI::Comm& comm = M.GetCommunicator();
    if (!M.IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix<T1, Prop1, Storage1, Allocator1>& >(M)
          .PrepareMltAdd();
      }
    
    // scattering column values
    int N = M.global_col_to_recv.GetM();
    Vector<int> Yrow(N), Yrow_proc(N);
    Vector<int> Xcol, Xcol_proc;
    if (comm.Get_size() > 1)
      {
        M.ScatterColValues(Y, Xcol);
        M.ScatterColValues(Yproc, Xcol_proc);
      }
    
    // local matrix
    const IVect& global = M.GetGlobalRowNumber();
    MltMin(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(M),
           global, Y, Yproc);
    
    // contributions of distant columns
    Yrow.Fill(0); Yrow_proc.Fill(comm.Get_size());
    for (int i = 0; i < M.dist_col.GetM(); i++)
      for (int j = 0; j < M.dist_col(i).GetM(); j++)
        {
          int k = M.dist_col(i).Index(j);
          if (Xcol_proc(k) < Yproc(i))
            {
              Yproc(i) = Xcol_proc(k);
              Y(i) = Xcol(k);
            }
          else if (Xcol_proc(k) == Yproc(i))
            {
              if (Xcol(k) < Y(i))
                Y(i) = Xcol(k);
            }

          if (Yproc(i) < Yrow_proc(k))
            {
              Yrow_proc(k) = Yproc(i);
              Yrow(k) = Y(i);
            }
          else if (Yproc(i) == Yrow_proc(k))
            {
              if (Y(i) < Yrow(k))
                Yrow(k) = Y(i);
            }
        }
    
    M.AssembleValuesMin(Yrow, Yrow_proc,
                        M.global_col_to_recv, M.ptr_global_col_to_recv, 
                        M.proc_col_to_recv, M.local_col_to_send,
			M.proc_col_to_send, Y, Yproc);
    
    // contributions of distant rows
    if (comm.Get_size() > 1)
      {
        M.ScatterRowValues(Y, Xcol);
        M.ScatterRowValues(Yproc, Xcol_proc);
      }

    N = M.global_row_to_recv.GetM();
    Yrow.Reallocate(N), Yrow_proc.Reallocate(N);
    Yrow.Fill(0); Yrow_proc.Fill(comm.Get_size());
    for (int i = 0; i < M.dist_row.GetM(); i++)
      for (int j = 0; j < M.dist_row(i).GetM(); j++)
        {
          int jrow = M.dist_row(i).Index(j);
          if (Xcol_proc(jrow) < Yproc(i))
            {
              Yproc(i) = Xcol_proc(jrow);
              Y(i) = Xcol(jrow);
            }
          else if (Xcol_proc(jrow) == Yproc(i))
            {
              if (Xcol(jrow) < Y(i))
                Y(i) = Xcol(jrow);
            }

          if (Yproc(i) < Yrow_proc(jrow))
            {
              Yrow_proc(jrow) = Yproc(i);
              Yrow(jrow) = Y(i);
            }
          else if (Yproc(i) == Yrow_proc(jrow))
            {
              if (Y(i) < Yrow(jrow))
                Yrow(jrow) = Y(i);
            }
        }

    // assembling row values
    if (comm.Get_size() > 1)
      {
        M.AssembleValuesMin(Yrow, Yrow_proc,
                            M.global_row_to_recv, M.ptr_global_row_to_recv,
                            M.proc_row_to_recv, M.local_row_to_send,
			    M.proc_row_to_send, Y, Yproc);
        
        // assembling rows shared between processors
        M.AssembleVecMin(Y, Yproc);
      }
  }

  
  //! Adds two distributed matrices (B = B + alpha A)
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void Add(const T0& alpha,
	   const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A, 
           DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    // adding local part
    Add(alpha,
	static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(A),
        static_cast<Matrix<T2, Prop2, Storage2, Allocator2>& >(B));
    
    const_cast<DistributedMatrix<T1, Prop1, Storage1, Allocator1>& >(A)
      .SwitchToGlobalNumbers();
    
    B.SwitchToGlobalNumbers();
    
    // then distant part
    for (int i = 0; i < A.dist_row.GetM(); i++)
      for (int j = 0; j < A.dist_row(i).GetM(); j++)
        B.AddRowDistantInteraction(A.dist_row(i).Index(j), i,
                                   A.proc_row(i)(j),
				   alpha*A.dist_row(i).Value(j));
    
    for (int i = 0; i < A.dist_col.GetM(); i++)
      for (int j = 0; j < A.dist_col(i).GetM(); j++)
        B.AddDistantInteraction(i, A.dist_col(i).Index(j),
                                A.proc_col(i)(j),
				alpha*A.dist_col(i).Value(j));
  }
  
  
  //! returns the maximal element (in absolute value) of the matrix
  /*!
    The result will be different from a sequential execution
    since A(i, j) is obtained as a sum of values stored
    in different processors
   */
  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  MaxAbs(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    typename ClassComplexType<T1>::Treal res, amax;
    res = MaxAbs(static_cast<const Matrix<T1, Prop1,
		 Storage1, Allocator1>& >(A));
    
    for (int i = 0; i < A.dist_row.GetM(); i++)
      for (int j = 0; j < A.dist_row(i).GetM(); j++)
        res = max(res, abs(A.dist_row(i).Value(j)));

    for (int i = 0; i < A.dist_col.GetM(); i++)
      for (int j = 0; j < A.dist_col(i).GetM(); j++)
        res = max(res, abs(A.dist_col(i).Value(j)));
    
    // selecting the maximum between processors
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return res;
    
    Vector<int64_t> xtmp;
    MpiAllreduce(comm, &res, xtmp, &amax, 1, MPI::MAX);
    
    return amax;
  }

  
  //! adds contribution of dist_col for GetRowSum
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetRowSumDistantCol(Vector<T0>& vec_sum) const
  {
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        vec_sum(i) += abs(dist_col(i).Value(j));
  }
  
  
  //! adds contribution of dist_row for GetRowSum
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetRowSumDistantRow(Vector<T0>& vec_sum) const
  {
    T0 zero; SetComplexZero(zero);
    Vector<T0> Y(global_row_to_recv.GetM());
    Y.Fill(zero);
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        {
          int jrow = dist_row(i).Index(j);
          Y(jrow) += abs(dist_row(i).Value(j));
        }
    
    AssembleRowValues(Y, vec_sum);
  }
  
  
  //! For each row of the matrix, computation of the sum of absolute values
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetRowSum(Vector<T0>& vec_sum,
                 const DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    GetRowSum(vec_sum,
	      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(A));
    
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    if (!A.IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix<T, Prop, Storage, Allocator>& >(A)
          .PrepareMltAdd();
      }
    
    A.GetRowSumDistantCol(vec_sum);    
    A.GetRowSumDistantRow(vec_sum);
    
    A.AssembleVec(vec_sum);
  }
  


  //! adds contribution of dist_col for GetColSum
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetColSumDistantCol(Vector<T0>& vec_sum) const
  {
    T0 zero; SetComplexZero(zero);
    Vector<T0> Y(global_col_to_recv.GetM());
    Y.Fill(zero);
    for (int i = 0; i < dist_col.GetM(); i++)
      for (int j = 0; j < dist_col(i).GetM(); j++)
        {
          int jrow = dist_col(i).Index(j);
          Y(jrow) += abs(dist_col(i).Value(j));
        }
    
    AssembleColValues(Y, vec_sum);
  }
  
  
  //! adds contribution of dist_row for GetColSum
  template<class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void DistributedMatrix<T, Prop, Storage, Allocator>
  ::GetColSumDistantRow(Vector<T0>& vec_sum) const
  {
    for (int i = 0; i < dist_row.GetM(); i++)
      for (int j = 0; j < dist_row(i).GetM(); j++)
        vec_sum(i) += abs(dist_row(i).Value(j));
  }
  
  
  //! For each column of the matrix, computation of 
  //! the sum of absolute values
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetColSum(Vector<T0>& vec_sum,
                 const DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    GetColSum(vec_sum,
	      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(A));
    
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    if (!A.IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix<T, Prop, Storage, Allocator>& >(A)
          .PrepareMltAdd();
      }
    
    A.GetColSumDistantCol(vec_sum);    
    A.GetColSumDistantRow(vec_sum);
    
    A.AssembleVec(vec_sum);
  }

  
  //! gets sum of absolute value for each row and each colum
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetRowColSum(Vector<T0>& sum_row,
                    Vector<T0>& sum_col,
                    const DistributedMatrix<T, Prop, Storage, Allocator> & A)
  {    
    GetRowColSum(sum_row, sum_col,
                 static_cast<const Matrix<T, Prop, Storage, Allocator>& >(A));
    
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    if (!A.IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix<T, Prop, Storage, Allocator>& >(A)
          .PrepareMltAdd();
      }

    A.GetRowSumDistantCol(sum_row);
    A.GetRowSumDistantRow(sum_row);
    
    A.GetColSumDistantCol(sum_col);    
    A.GetColSumDistantRow(sum_col);
    
    A.AssembleVec(sum_row);
    A.AssembleVec(sum_col);
  }
  
  
  //! returns the 1-norm of A
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  Norm1(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return Norm1(static_cast<const Matrix<T1, Prop1,
		   Storage1, Allocator1>& >(A));
    
    Vector<typename ClassComplexType<T1>::Treal> sum_col;
    GetColSum(sum_col, A);

    typename ClassComplexType<T1>::Treal res, amax;
    amax = 0;
    for (int i = 0; i < sum_col.GetM(); i++)
      amax = max(amax, abs(sum_col(i)));
    
    Vector<int64_t> xtmp;
    MpiAllreduce(comm, &amax, xtmp, &res, 1, MPI::MAX);
    return res;
  }


  //! returns the infinity norm of A
  /*!
    The result is different from a sequential execution (because 
    some non-zero entries may be shared by different processors)
   */
  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  NormInf(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return NormInf(static_cast<const Matrix<T1, Prop1,
		     Storage1, Allocator1>& >(A));
    
    Vector<typename ClassComplexType<T1>::Treal> sum_row;
    GetRowSum(sum_row, A);

    typename ClassComplexType<T1>::Treal res, amax;
    amax = 0;
    for (int i = 0; i < sum_row.GetM(); i++)
      amax = max(amax, abs(sum_row(i)));
    
    Vector<int64_t> xtmp;
    MpiAllreduce(comm, &amax, xtmp, &res, 1, MPI::MAX);
    return res;
  }


  //! replaces A by its transpose
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Transpose(DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    Transpose(static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(A));
    
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    SwapPointer(A.dist_row, A.dist_col);
    SwapPointer(A.proc_row, A.proc_col);
    SwapPointer(A.global_row_to_recv, A.global_col_to_recv);
    SwapPointer(A.ptr_global_row_to_recv, A.ptr_global_col_to_recv);
    SwapPointer(A.local_row_to_send, A.local_col_to_send);
    SwapPointer(A.proc_row_to_recv, A.proc_col_to_recv);
    SwapPointer(A.proc_row_to_send, A.proc_col_to_send);
    
    int itmp = A.size_max_distant_row;
    A.size_max_distant_row = A.size_max_distant_col;
    A.size_max_distant_col = itmp;
  }
  

  //! replaces A by its conjugate
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Conjugate(DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    Conjugate(static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(A));
    
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    for (int i = 0; i < A.dist_row.GetM(); i++)
      for (int j = 0; j < A.dist_row(i).GetM(); j++)
        A.dist_row(i).Value(j) = conjugate(A.dist_row(i).Value(j));
    
    for (int i = 0; i < A.dist_col.GetM(); i++)
      for (int j = 0; j < A.dist_col(i).GetM(); j++)
        A.dist_col(i).Value(j) = conjugate(A.dist_col(i).Value(j));
  }
  
  
  //! computes B = transpose(A)
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Transpose(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
                 DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    Transpose(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(A),
              static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(B));
    
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    B.Init(A);
    B.dist_col = A.dist_row;
    B.proc_col = A.proc_row;
    B.global_col_to_recv = A.global_row_to_recv;
    B.ptr_global_col_to_recv = A.ptr_global_row_to_recv;
    B.local_col_to_send = A.local_row_to_send;
    B.proc_col_to_recv = A.proc_row_to_recv;
    B.proc_col_to_send = A.proc_row_to_send;

    B.dist_row = A.dist_col;
    B.proc_row = A.proc_col;
    B.global_row_to_recv = A.global_col_to_recv;
    B.ptr_global_row_to_recv = A.ptr_global_col_to_recv;
    B.local_row_to_send = A.local_col_to_send;
    B.proc_row_to_recv = A.proc_col_to_recv;
    B.proc_row_to_send = A.proc_col_to_send;
    B.local_number_distant_values = B.local_number_distant_values;
    
    B.size_max_distant_row = A.size_max_distant_col;
    B.size_max_distant_col = A.size_max_distant_row;
  }
  
  
  //! Matrix transposition and conjugation.
  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                     DistributedMatrix<T, Prop, Storage, Allocator>& B)
  {
    Transpose(A, B);
    Conjugate(B);
  }


  //! Matrix transposition and conjugation.
  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    Transpose(A);
    Conjugate(A);
  }
  

  //! computes C = A B
  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T4, class Prop4, class Storage4, class Allocator4>
  void Mlt(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
           const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
           DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return
	Mlt(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(A),
	    static_cast<const Matrix<T2, Prop2, Storage2, Allocator2>& >(B),
	    static_cast<Matrix<T4, Prop4, Storage4, Allocator4>& >(C));
    
    cout << "Mlt not implemented for distributed matrices" << endl;
    abort();
  }

  
  //! computes C = beta C + alpha A B
  template<class T0,
           class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T3,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
              const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
              const T3& beta,
              DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return
	MltAdd(alpha,
	       static_cast<const Matrix<T1, Prop1,Storage1,Allocator1>& >(A),
	       static_cast<const Matrix<T2, Prop2,Storage2,Allocator2>& >(B),
	       beta,
	       static_cast<Matrix<T4, Prop4, Storage4, Allocator4>& >(C));
    
    cout << "MltAdd not implemented for distributed matrices" << endl;
    abort();
  }

  
  //! computes C = beta C + alpha A B (with transposes)
  template<class T0, class TransA,
           class T1, class Prop1, class Storage1, class Allocator1,
           class TransB,
	   class T2, class Prop2, class Storage2, class Allocator2,
           class T3,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T0& alpha, const TransA& transA,
              const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
              const TransB& transB,
              const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
              const T3& beta,
              DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return MltAdd(alpha, transA,
                    static_cast<const Matrix<T1, Prop1,
		    Storage1, Allocator1>& >(A),
                    transB,
                    static_cast<const Matrix<T2, Prop2,
		    Storage2, Allocator2>& >(B),
                    beta,
                    static_cast<Matrix<T4, Prop4,
		    Storage4, Allocator4>& >(C));
    
    cout << "MltAdd not implemented for distributed matrices" << endl;
    abort();
  }
  

  //! X = M(i, :)
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void GetRow(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
              int i, Vector<T1, VectSparse, Allocator1>& X)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return GetRow(static_cast<const Matrix<T0, Prop0,
		    Storage0, Allocator0>& >(A), i, X);
    
    cout << "GetRow not implemented for distributed matrices" << endl;
    abort();
  }


  //! X = M(:, i)
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void GetCol(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
              int i, Vector<T1, VectSparse, Allocator1>& X)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return GetCol(static_cast<const Matrix<T0, Prop0,
		    Storage0, Allocator0>& >(A), i, X);
    
    cout << "GetCol not implemented for distributed matrices" << endl;
    abort();
  }
  
  
  //! M(i, :) = X
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
              int i, DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return SetRow(X, i, static_cast<Matrix<T0, Prop0,
		    Storage0, Allocator0>& >(A));
    
    cout << "SetRow not implemented for distributed matrices" << endl;
    abort();
  }
  

  //! M(:, i) = X
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
              int i, DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return SetCol(X, i, static_cast<Matrix<T0, Prop0,
		    Storage0, Allocator0>& >(A));
    
    cout << "SetCol not implemented for distributed matrices" << endl;
    abort();
  }

  
  //! applies permutation on rows and columns of A
  template<class T, class Prop, class Storage, class Allocator>
  void ApplyPermutation(DistributedMatrix<T, Prop, Storage, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm)
  {
    MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return ApplyPermutation(static_cast<Matrix<T, Prop,
			      Storage, Allocator>& >(A),
                              row_perm, col_perm);
    
    cout << "ApplyPermutation not implemented for distributed matrices"
	 << endl;
    abort();
  }
  

  //! applies permutation on rows and columns of A
  template<class T, class Prop, class Storage, class Allocator>
  void ApplyInversePermutation(DistributedMatrix<T, Prop,
			       Storage, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm)
  {
    MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return ApplyInversePermutation(static_cast<Matrix<T, Prop, Storage,
				     Allocator>& >(A),
                                     row_perm, col_perm);
    
    cout << "ApplyInversePermutation not implemented for distributed matrices"
	 << endl;
    abort();
  }

  
  //! applies successive over-relaxation steps to solve A X = B
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SOR(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   Vector<T2, Storage2, Allocator2>& X,
	   const Vector<T1, Storage1, Allocator1>& B,
	   const T3& omega, int iter, int type_ssor)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return SOR(static_cast<const Matrix<T0, Prop0,
		 Storage0, Allocator0>& >(A),
                 X, B, omega, iter, type_ssor);
    
    cout << "SOR not implemented for distributed matrices" << endl;
    abort();    
  }
  
  
  //! applies successive over-relaxation steps to solve A^T X = B
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SOR(const class_SeldonTrans& transM,
	   const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   Vector<T2, Storage2, Allocator2>& X,
	   const Vector<T1, Storage1, Allocator1>& B,
	   const T3& omega, int iter, int type_ssor)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return SOR(transM,
                 static_cast<const Matrix<T0, Prop0,
		 Storage0, Allocator0>& >(A),
                 X, B, omega, iter, type_ssor);
    
    cout << "SOR not implemented for distributed matrices" << endl;
    abort();
  }
  
  
  //! extracts a set of sparse columns from matrix A
  template<class T1, class Prop, class Storage, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const DistributedMatrix<T1, Prop, Storage, Allocator>& A,
              const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>,
	      VectSparse, Allocator3>& V)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return GetCol(static_cast<const Matrix<T1, Prop,
		    Storage, Allocator>& >(A),
                    col_number, V);
    
    cout << "GetCol not implemented for distributed matrices" << endl;
    abort();
  }
  
  
  //! converts a distributed matrix to another one
  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void Copy(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
            DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    Copy(static_cast<const Matrix<T1, Prop1, Storage1, Allocator1>& >(A),
         static_cast<Matrix<T2, Prop2, Storage2, Allocator2>& >(B));
    
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    B.OverlapRowNumbers = A.OverlapRowNumbers;
    B.OverlapProcNumbers = A.OverlapProcNumbers;
    B.GlobalRowNumbers = A.GlobalRowNumbers;
    B.ProcSharingRows = A.ProcSharingRows;
    B.SharingRowNumbers = A.SharingRowNumbers;
    B.nodl_scalar_ = A.nodl_scalar_;
    B.nb_unknowns_scal_ = A.nb_unknowns_scal_;
    B.nglob_ = A.nglob_;
    B.comm_ = A.comm_;
    B.dist_col = A.dist_col;
    B.dist_row = A.dist_row;
    B.proc_col = A.proc_col;
    B.proc_row = A.proc_row;
    B.global_row_to_recv = A.global_row_to_recv;
    B.global_col_to_recv = A.global_col_to_recv;
    B.ptr_global_row_to_recv = A.ptr_global_row_to_recv;
    B.ptr_global_col_to_recv = A.ptr_global_col_to_recv;
    B.local_row_to_send = A.local_row_to_send;
    B.local_col_to_send = A.local_col_to_send;
    B.proc_col_to_recv = A.proc_col_to_recv;
    B.proc_col_to_send = A.proc_col_to_send;
    B.proc_row_to_recv = A.proc_row_to_recv;
    B.proc_row_to_send = A.proc_row_to_send;
    B.local_number_distant_values = A.local_number_distant_values;
    B.size_max_distant_row = A.size_max_distant_row;
    B.size_max_distant_col = A.size_max_distant_col;
  }
  
  
  //! extracts real part of a matrix
  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void CopyReal(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
		DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    MPI::Comm& comm = B.GetCommunicator();
    if (comm.Get_size() == 1)
      return CopyReal(static_cast<const Matrix<T1, Prop1,
		      Storage1, Allocator1>& >(A),
                      static_cast<Matrix<T2, Prop2,
		      Storage2, Allocator2>& >(B));

    cout << "CopyReal not implemented for distributed matrices" << endl;
    abort();
  }
  
  
  //! extracts a sub-matrix B from matrix A
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                    int m, int n,
		    DistributedMatrix<T, Prop, Storage, Allocator>& B)
  {
    MPI::Comm& comm = B.GetCommunicator();
    if (comm.Get_size() == 1)
      return GetSubMatrix(static_cast<const Matrix<T, Prop,
			  Storage, Allocator>& >(A),
                          m, n, static_cast<Matrix<T, Prop,
			  Storage, Allocator>& >(B));

    cout << "GetSubMatrix not implemented for distributed matrices" << endl;
    abort();    
  }
  
  
  //! extracts a sub-matrix B from matrix A
  template<class T, class Prop, class Storage, class Allocator>
  void GetSubMatrix(const DistributedMatrix<T, Prop, Storage, Allocator>& A,
                    int m1, int m2, int n1, int n2,
                    DistributedMatrix<T, Prop, Storage, Allocator>& B)
  {
    MPI::Comm& comm = B.GetCommunicator();
    if (comm.Get_size() == 1)
      return GetSubMatrix(static_cast<const Matrix<T, Prop,
			  Storage, Allocator>& >(A),
                          m1, m2, n1, n2,
                          static_cast<Matrix<T, Prop,
			  Storage, Allocator>& >(B));

    cout << "GetSubMatrix not implemented for distributed matrices" << endl;
    abort();    
  }
  
  
  //! returns Froebenius norm of matrix A
  template<class T, class Prop, class Storage, class Allocator>
  typename ClassComplexType<T>::Treal
  NormFro(const DistributedMatrix<T, Prop, Storage, Allocator>& A)
  {
    const MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return NormFro(static_cast<const Matrix<T, Prop,
		     Storage, Allocator>& >(A));

    cout << "NormFro not implemented for distributed matrices" << endl;
    abort();    
  }
  

  //! multiplies rows by coefficients
  template<class T, class Storage, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(DistributedMatrix<T, General, Storage, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow)
  {
    ScaleLeftMatrix(static_cast<Matrix<T, General,
		    Storage, Allocator>& >(A), Drow);
    
    MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    if (!A.IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix<T, General, Storage, Allocator>& >(A)
          .PrepareMltAdd();
      }
    
    // scaling distant columns
    for (int i = 0; i < A.dist_col.GetM(); i++)
      A.dist_col(i) *= Drow(i);
    
    // scaling distant rows
    Vector<T1, VectFull, Allocator1> Drow_glob;
    A.ScatterRowValues(Drow, Drow_glob);
    
    for (int i = 0; i < A.dist_row.GetM(); i++)
      for (int j = 0; j < A.dist_row(i).GetM(); j++)
        A.dist_row(i).Value(j) *= Drow_glob(A.dist_row(i).Index(j));
  }


  //! multiplies columns by coefficients
  template<class T, class Storage, class Allocator,
           class T1, class Allocator1>
  void ScaleRightMatrix(DistributedMatrix<T, General, Storage, Allocator>& A,
                        const Vector<T1, VectFull, Allocator1>& Dcol)
  {
    ScaleRightMatrix(static_cast<Matrix<T, General,
		     Storage, Allocator>& >(A), Dcol);
    
    MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;

    if (!A.IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix<T, General, Storage, Allocator>& >(A)
          .PrepareMltAdd();
      }
    
    // scaling distant rows
    for (int i = 0; i < A.dist_row.GetM(); i++)
      A.dist_row(i) *= Dcol(i);
    
    // scaling distant columns
    Vector<T1, VectFull, Allocator1> Dcol_glob;
    A.ScatterColValues(Dcol, Dcol_glob);
    
    for (int i = 0; i < A.dist_col.GetM(); i++)
      for (int j = 0; j < A.dist_col(i).GetM(); j++)
        A.dist_col(i).Value(j) *= Dcol_glob(A.dist_col(i).Index(j));    
  }


  //! multiplies rows and columns by coefficients
  template<class T, class Prop, class Storage, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(DistributedMatrix<T, Prop, Storage, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    ScaleMatrix(static_cast<Matrix<T, Prop, Storage, Allocator>& >(A),
		Drow, Dcol);
    
    MPI::Comm& comm = A.GetCommunicator();
    if (comm.Get_size() == 1)
      return;

    if (!A.IsReadyForMltAdd())
      {
        // preparing the matrix vector product
        // this method will be called once for the first matrix-vector product
        const_cast<DistributedMatrix<T, Prop, Storage, Allocator>& >(A)
          .PrepareMltAdd();
      }
    
    // scaling distant columns
    for (int i = 0; i < A.dist_col.GetM(); i++)
      A.dist_col(i) *= Drow(i);
    
    // scaling distant rows
    Vector<T1, VectFull, Allocator1> Drow_glob;
    A.ScatterRowValues(Drow, Drow_glob);
    
    for (int i = 0; i < A.dist_row.GetM(); i++)
      for (int j = 0; j < A.dist_row(i).GetM(); j++)
        A.dist_row(i).Value(j) *= Drow_glob(A.dist_row(i).Index(j));

    // scaling distant rows
    for (int i = 0; i < A.dist_row.GetM(); i++)
      A.dist_row(i) *= Dcol(i);
    
    // scaling distant columns
    Vector<T1, VectFull, Allocator1> Dcol_glob;
    A.ScatterColValues(Dcol, Dcol_glob);
    
    for (int i = 0; i < A.dist_col.GetM(); i++)
      for (int j = 0; j < A.dist_col(i).GetM(); j++)
        A.dist_col(i).Value(j) *= Dcol_glob(A.dist_col(i).Index(j));    
  }
  
  
  //! assembling a distributed matrix, global rows are given in CSR form
  /*!
    \param[in] A distributed matrix
    \param[in] sym the matrix is assumed to be symmetric,
    only upper part is assembled
    \param[in] comm MPI communicator
    \param[out] row_numbers global row numbers
    \param[out] local_row_numbers local row numbers
    \param[out] PtrA Index of first element stored for each row
    \param[out] IndA global column numbers
    \param[out] ValA values
    \param[in] sym_pattern not used
   */
  template<class MatrixSparse, class Tint, class T>
  void AssembleDistributed(MatrixSparse& A,
			   Symmetric& sym, const MPI::Comm& comm,
                           IVect& row_numbers, IVect& local_row_numbers,
			   Vector<Tint>& PtrA, Vector<Tint>& IndA,
			   Vector<T>& ValA, bool sym_pattern)
  {
    int m = A.GetM();
    int n = A.GetGlobalM();
    
    PtrA.Clear();
    IndA.Clear();
    ValA.Clear();

    // we convert A in ArrayRowSparse format
    Matrix<T, General, ArrayRowSparse> B;
    A.GetDistributedRows(B);
    
    // global to local conversion
    IVect Glob_to_local(n);
    const IVect& RowNumber = A.GetGlobalRowNumber();
    Glob_to_local.Fill(-1);
    for (int i = 0; i < m; i++)
      Glob_to_local(RowNumber(i)) = i;
    
    const IVect& OverlapRowNumber = A.GetOverlapRowNumber();
    const IVect& OverlapProcNumber = A.GetOverlapProcNumber();
    
    // now sending/receiving overlapped rows, and distant rows
    int nb_proc = comm.Get_size();
    int rank = comm.Get_rank();
    Vector<MPI::Request> request(3*nb_proc);
    Vector<bool> request_active(nb_proc);
    request_active.Fill(false);
    IVect nsend_int(nb_proc);
    Vector<IVect> EntierToSend(nb_proc);
    Vector<Vector<T> > FloatToSend(nb_proc);
    Vector<Vector<int64_t> > FloatToSend_tmp(nb_proc);
    for (int i = 0; i < nb_proc; i++)
      {
	int nrow = 0, nsend_float = 0, irow;
	nsend_int(i) = 2;
	for (int j = 0; j < OverlapRowNumber.GetM(); j++)
	  if (OverlapProcNumber(j) == i)
	    {
	      nrow++; nsend_int(i) += 2;
	      irow = OverlapRowNumber(j);
	      nsend_int(i) += B.GetRowSize(irow);
	      nsend_float += B.GetRowSize(irow);
	    }
        
        for (int j = 0; j < A.dist_row.GetM(); j++)
          for (int k = 0; k < A.dist_row(j).GetM(); k++)
            if (A.proc_row(j)(k) == i)
              {
                irow = A.dist_row(j).Index(k);
                if (A.local_number_distant_values)
                  irow = A.global_row_to_recv(irow);
                
                if (irow <= RowNumber(j))
                  {
                    nrow++; nsend_int(i) += 3;
                    nsend_float++;
                  }
              }
        
        if (nrow <= 0)
          nsend_int(i) = 0;
        
        if (i != rank)
          request(i) = comm.Isend(&nsend_int(i), 1, MPI::INTEGER, i, 4);
        
	if (nrow > 0)
	  {
	    request_active(i) = true;
	    EntierToSend(i).Reallocate(nsend_int(i));
	    FloatToSend(i).Reallocate(nsend_float);
	    EntierToSend(i)(0) = nsend_float;
            EntierToSend(i)(1) = nrow;	    
	    nsend_int(i) = 2; nsend_float = 0;
            // storing values and indices of a row shared with processor i
            // processor i is the owner of this row
	    for (int j = 0; j < OverlapRowNumber.GetM(); j++)
	      if (OverlapProcNumber(j) == i)
		{
		  irow = OverlapRowNumber(j);
		  EntierToSend(i)(nsend_int(i)++) = RowNumber(irow);
		  EntierToSend(i)(nsend_int(i)++) = B.GetRowSize(irow);
		  for (int k = 0; k < B.GetRowSize(irow); k++)
		    {
		      EntierToSend(i)(nsend_int(i)++) = B.Index(irow, k);
		      FloatToSend(i)(nsend_float++) = B.Value(irow, k);
		    }
		}

            // storing values of row associated with processor i
            for (int j = 0; j < A.dist_row.GetM(); j++)
              for (int k = 0; k < A.dist_row(j).GetM(); k++)
                if (A.proc_row(j)(k) == i)
                  {
                    irow = A.dist_row(j).Index(k);
                    if (A.local_number_distant_values)
                      irow = A.global_row_to_recv(irow);
                    
                    if (irow <= RowNumber(j))
                      {
                        EntierToSend(i)(nsend_int(i)++) = irow;
                        EntierToSend(i)(nsend_int(i)++) = 1;
                        EntierToSend(i)(nsend_int(i)++) = RowNumber(j);
                        FloatToSend(i)(nsend_float++)
			  = A.dist_row(j).Value(k);
                      }
                  }
        
            // sending all the values and indices stored to the processor i
            // only rows are sent
	    request(i+nb_proc) = 
              comm.Isend(EntierToSend(i).GetData(), nsend_int(i),
			 MPI::INTEGER, i, 5);
	    
            request(i+2*nb_proc) = 
              MpiIsend(comm, FloatToSend(i), FloatToSend_tmp(i),
		       nsend_float, i, 6);
	  }
	else
	  nsend_int(i) = 0;
       
      }

    A.Clear();
    
    // receiving the number of entries
    MPI::Status status; IVect nrecv_int(nb_proc); nrecv_int.Fill(0);
    Vector<IVect> EntierToRecv(nb_proc);
    Vector<Vector<T> > FloatToRecv(nb_proc);
    Vector<int64_t> FloatToRecv_tmp;
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        comm.Recv(&nrecv_int(i), 1, MPI::INTEGER, i, 4, status);
    
    // waiting for sending of nsend_int effective
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        request(i).Wait(status);
    
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            EntierToRecv(i).Reallocate(nrecv_int(i));
            comm.Recv(EntierToRecv(i).GetData(), nrecv_int(i),
                      MPI::INTEGER, i, 5, status);
          }
      }
    
    // waiting for sending of EntierToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (request_active(i))
	request(i+nb_proc).Wait(status);
    
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            int nsend_float = EntierToRecv(i)(0);
            
            FloatToRecv(i).Reallocate(nsend_float);
            MpiRecv(comm, FloatToRecv(i), FloatToRecv_tmp,
		    nsend_float, i, 6, status);
          }
      }
    
    // waiting for sending of FloatToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (request_active(i))
	request(i + 2*nb_proc).Wait(status);
    
    // deleting arrays
    request_active.Fill(false);
    for (int i = 0; i < nb_proc; i++)
      {
        nsend_int(i) = 0;
        EntierToSend(i).Clear();
        FloatToSend(i).Clear();
      }

    // constructing array to know if a column is overlapped 
    // (treated by another processor)
    IVect OverlappedCol(n); OverlappedCol.Fill(-1);
    for (int i = 0; i < OverlapRowNumber.GetM(); i++)
      OverlappedCol(OverlapRowNumber(i)) = i;

    // assembling matrix B with interactions coming from other processors
    for (int i = 0; i < nb_proc; i++)
      {       
        if (nrecv_int(i) > 0)
          {	
            int nrow = EntierToRecv(i)(1);    
            // loop over rows
            nrecv_int(i) = 2; int nsend_float = 0;
            for (int j = 0; j < nrow; j++)
              {
                int iglob = EntierToRecv(i)(nrecv_int(i)++);
                int irow = Glob_to_local(iglob);
                int size_row = EntierToRecv(i)(nrecv_int(i)++);
                IVect index(size_row); Vector<T> values(size_row);
                for (int k = 0; k < size_row; k++)
                  {
                    index(k) = EntierToRecv(i)(nrecv_int(i)++);
                    values(k) = FloatToRecv(i)(nsend_float++);
                  }
                
                // adding to matrix B if the row is not shared
                // otherwise we send the row to the original processor
                if (OverlappedCol(irow) == -1)
                  {
                    if (size_row == 1)
                      B.AddInteraction(irow, index(0), values(0));
                    else
                      B.AddInteractionRow(irow, size_row, index, values);
                  }
                else
                  {
                    int irow_ = OverlappedCol(irow);
                    int proc = OverlapProcNumber(irow_);
                    
                    int offset_int(2), offset_float(0);
                    if (nsend_int(proc) == 0)
                      {
                        nsend_int(proc) = 2;
                        EntierToSend(proc).Reallocate(size_row+4);
                        FloatToSend(proc).Reallocate(size_row);
                        EntierToSend(proc)(0) = 0;
                        EntierToSend(proc)(1) = 0;
                      }
                    else
                      {
                        offset_int = EntierToSend(proc).GetM();
                        offset_float = FloatToSend(proc).GetM();
                        EntierToSend(proc).Resize(size_row+2+offset_int);
                        FloatToSend(proc).Resize(size_row+offset_float);
                      }
                    
                    nsend_int(proc) += size_row+2;
                    EntierToSend(proc)(0) += size_row;
                    EntierToSend(proc)(1)++;
                    EntierToSend(proc)(offset_int++) = iglob;
                    EntierToSend(proc)(offset_int++) = size_row;
                    for (int k = 0; k < size_row; k++)
                      {
                        EntierToSend(proc)(offset_int++) = index(k);
                        FloatToSend(proc)(offset_float++) = values(k);
                      }                    
                  }
              }
          }
      }

    // sending rows overlapped to original processor
    // this case occurs only if distant rows contain overlapped rows
    for (int i = 0; i < nb_proc; i++)
      {
        if (i != rank)
          {
            request(i) = comm.Isend(&nsend_int(i), 1, MPI::INTEGER, i, 4);
            if (nsend_int(i) > 0)
              {
                request_active(i) = true;
                
                request(i+nb_proc) = 
                  comm.Isend(EntierToSend(i).GetData(), nsend_int(i),
			     MPI::INTEGER, i, 5);
                
                request(i+2*nb_proc) = 
                  MpiIsend(comm, FloatToSend(i), FloatToSend_tmp(i),
                           FloatToSend(i).GetM(), i, 6);            
              }
          }
      }
    
    // receiving rows
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        comm.Recv(&nrecv_int(i), 1, MPI::INTEGER, i, 4, status);
    
    // trying to receive EntierToRecv
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            EntierToRecv(i).Reallocate(nrecv_int(i));
            comm.Recv(EntierToRecv(i).GetData(), nrecv_int(i),
                      MPI::INTEGER, i, 5, status);
          }
      }

    // waiting for sending of EntierToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (request_active(i))
	request(i+nb_proc).Wait(status);
        
    // trying to receive FloatToRecv
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            int nsend_float = EntierToRecv(i)(0);
            
            FloatToRecv(i).Reallocate(nsend_float);
            MpiRecv(comm, FloatToRecv(i), FloatToRecv_tmp,
		    nsend_float, i, 6, status);
          }
      }
    
    // waiting for sending of FloatToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (request_active(i))
	request(i + 2*nb_proc).Wait(status);
    
    // deleting arrays
    request_active.Fill(false);
    for (int i = 0; i < nb_proc; i++)
      {
        nsend_int(i) = 0;
        EntierToSend(i).Clear();
        FloatToSend(i).Clear();
      }

    // assembling matrix B with last interactions coming from other processors
    for (int i = 0; i < nb_proc; i++)
      {       
        if (nrecv_int(i) > 0)
          {	
            int ncol = EntierToRecv(i)(1);    
            // loop over rows
            nrecv_int(i) = 2; int nsend_float = 0;
            for (int j = 0; j < ncol; j++)
              {
                int iglob = EntierToRecv(i)(nrecv_int(i)++);
                int irow = Glob_to_local(iglob);
                int size_row = EntierToRecv(i)(nrecv_int(i)++);
                IVect index(size_row); Vector<T> values(size_row);
                for (int k = 0; k < size_row; k++)
                  {
                    index(k) = EntierToRecv(i)(nrecv_int(i)++);
                    values(k) = FloatToRecv(i)(nsend_float++);
                  }
                
                // adding to matrix B if the column is not shared
                // otherwise an error is generated (impossible case)
                if (OverlappedCol(irow) == -1)
                  {
                    if (size_row == 1)
                      B.AddInteraction(irow, index(0), values(0));
                    else
                      B.AddInteractionRow(irow, size_row, index, values);
                  }
                else
                  {
                    cout << "Case impossible, this row is not overlapped"
			 << endl;
                    abort();
                  }
              }
          }
      }    
    
    // now we convert Bh to RowSparse while removing overlapped rows
    Glob_to_local.Clear();
    
    int nloc = m - OverlapRowNumber.GetM();
    PtrA.Reallocate(nloc+1);
    row_numbers.Reallocate(nloc);
    int nrow = 0, nnz = 0;
    PtrA(nrow) = 0;
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
	{
	  row_numbers(nrow) = RowNumber(i);
	  PtrA(nrow+1) = PtrA(nrow) + B.GetRowSize(i);
	  nrow++;
          nnz += B.GetRowSize(i);
	}
    
    IndA.Reallocate(nnz);
    ValA.Reallocate(nnz); nrow = 0; nnz = 0;
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
	{
	  for (int j = 0; j < B.GetRowSize(i); j++)
	    {
              IndA(nnz) = B.Index(i, j);
	      ValA(nnz) = B.Value(i, j);
	      nnz++;
	    }
	  nrow++;
	}
    
    local_row_numbers.Reallocate(nloc);
    nrow = 0;
    for (int i = 0; i < m; i++)
      if (OverlappedCol(i) == -1)
	local_row_numbers(nrow++) = i;
    
  }
  
  
  //! assembling a distributed matrix, global columns are given in CSC form
  /*!
    \param[in] A distributed matrix
    \param[in] prop the matrix is assumed to be unsymmetric
    \param[in] comm MPI communicator
    \param[out] col_numbers global column numbers
    \param[out] local_col_numbers local column numbers
    \param[out] PtrA Index of first element stored for each column
    \param[out] IndA global row numbers
    \param[out] ValA values
    \param[in] sym_pattern if true,
    the pattern of the unsymmetric matrix is symmetrized
   */
  template<class MatrixSparse, class Tint, class T>
  void AssembleDistributed(MatrixSparse& A,
			   General& prop, const MPI::Comm& comm,
                           IVect& col_numbers, IVect& local_col_numbers,
                           Vector<Tint>& PtrA, Vector<Tint>& IndA,
			   Vector<T>& ValA, bool sym_pattern)
  {
    int n = A.GetN();
    int m = A.GetGlobalM();
    
    PtrA.Clear();
    IndA.Clear();
    ValA.Clear();
    
    // we convert A in CSC format
    Matrix<T, General, ArrayColSparse> B;
    A.GetDistributedColumns(B, sym_pattern);
    
    // global to local conversion
    IVect Glob_to_local(m);
    const IVect& RowNumber = A.GetGlobalRowNumber();
    Glob_to_local.Fill(-1);
    for (int i = 0; i < n; i++)
      Glob_to_local(RowNumber(i)) = i;
    
    const IVect& OverlapRowNumber = A.GetOverlapRowNumber();
    const IVect& OverlapProcNumber = A.GetOverlapProcNumber();
    
    // we send to each processor additional entries due to overlapping
    // distant columns or distant rows (because of symmetrisation of patterns
    int nb_proc = comm.Get_size();
    int rank = comm.Get_rank();
    Vector<MPI::Request> request(3*nb_proc);
    Vector<bool> request_active(nb_proc);
    request_active.Fill(false);
    IVect nsend_int(nb_proc);
    Vector<IVect> EntierToSend(nb_proc);
    Vector<Vector<T> > FloatToSend(nb_proc);
    Vector<Vector<int64_t> > FloatToSend_tmp(nb_proc);
    // loop over processors
    for (int i = 0; i < nb_proc; i++)
      {
	int nrow = 0, nsend_float = 0, irow;
	nsend_int(i) = 2;
        for (int j = 0; j < OverlapRowNumber.GetM(); j++)
	  if (OverlapProcNumber(j) == i)
	    {
	      nrow++; nsend_int(i) += 2;
	      irow = OverlapRowNumber(j);
	      nsend_int(i) += B.GetColumnSize(irow);
	      nsend_float += B.GetColumnSize(irow);
	    }
        
        if (sym_pattern)
          for (int j = 0; j < A.dist_row.GetM(); j++)
            for (int k = 0; k < A.dist_row(j).GetM(); k++)
              if (A.proc_row(j)(k) == i)
                {
                  nrow++; nsend_int(i) += 3;
                  nsend_float++;
                }
        
        for (int j = 0; j < A.dist_col.GetM(); j++)
          for (int k = 0; k < A.dist_col(j).GetM(); k++)
            if (A.proc_col(j)(k) == i)
              {
                nrow++; nsend_int(i) += 3;
                nsend_float++;
              }
        
        if (nrow <= 0)
          nsend_int(i) = 0;
        
        if (i != rank)
          request(i) = comm.Isend(&nsend_int(i), 1, MPI::INTEGER, i, 4);
        
        if (nrow > 0)
	  {
	    request_active(i) = true;
	    EntierToSend(i).Reallocate(nsend_int(i));
	    FloatToSend(i).Reallocate(nsend_float);
	    EntierToSend(i)(0) = nsend_float;
            EntierToSend(i)(1) = nrow;	    
	    nsend_int(i) = 2; nsend_float = 0;
            // storing values and indices of a column shared with processor i
            // processor i is the owner of this column
	    for (int j = 0; j < OverlapRowNumber.GetM(); j++)
	      if (OverlapProcNumber(j) == i)
		{
		  irow = OverlapRowNumber(j);
		  EntierToSend(i)(nsend_int(i)++) = RowNumber(irow);
		  EntierToSend(i)(nsend_int(i)++) = B.GetColumnSize(irow);
		  for (int k = 0; k < B.GetColumnSize(irow); k++)
		    {
		      EntierToSend(i)(nsend_int(i)++) = B.Index(irow, k);
		      FloatToSend(i)(nsend_float++) = B.Value(irow, k);
		    }
		}
            
            // storing values to enforce a symmetric pattern
            if (sym_pattern)
              for (int j = 0; j < A.dist_row.GetM(); j++)
                for (int k = 0; k < A.dist_row(j).GetM(); k++)
                  if (A.proc_row(j)(k) == i)
                    {
                      irow = A.dist_row(j).Index(k);
                      if (A.local_number_distant_values)
                        irow = A.global_row_to_recv(irow);

                      EntierToSend(i)(nsend_int(i)++) = irow;
                      EntierToSend(i)(nsend_int(i)++) = 1;
                      EntierToSend(i)(nsend_int(i)++) = RowNumber(j);
                      FloatToSend(i)(nsend_float++) = 0;
                    }
            
            // storing values of columns associated with processor i
            for (int j = 0; j < A.dist_col.GetM(); j++)
              for (int k = 0; k < A.dist_col(j).GetM(); k++)
                if (A.proc_col(j)(k) == i)
                  {
                    irow = A.dist_col(j).Index(k);
                    if (A.local_number_distant_values)
                      irow = A.global_col_to_recv(irow);

                    EntierToSend(i)(nsend_int(i)++) = irow;
                    EntierToSend(i)(nsend_int(i)++) = 1;
                    EntierToSend(i)(nsend_int(i)++) = RowNumber(j);
                    FloatToSend(i)(nsend_float++) = A.dist_col(j).Value(k);
                  }
            
            // sending all the values and indices stored to the processor i
            // only columns are sent
            request(i+nb_proc) = 
              comm.Isend(EntierToSend(i).GetData(), nsend_int(i),
			 MPI::INTEGER, i, 5);
	    
            request(i+2*nb_proc) = 
              MpiIsend(comm, FloatToSend(i), FloatToSend_tmp(i),
		       nsend_float, i, 6);
	  }
	else
	  nsend_int(i) = 0;
       
      }
    
    A.Clear();
    
    // receiving the number of entries
    MPI::Status status; IVect nrecv_int(nb_proc); nrecv_int.Fill(0);
    Vector<IVect> EntierToRecv(nb_proc);
    Vector<Vector<T> > FloatToRecv(nb_proc);
    Vector<int64_t> FloatToRecv_tmp;
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        comm.Recv(&nrecv_int(i), 1, MPI::INTEGER, i, 4, status);
    
    // waiting for sending of nsend_int effective
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        request(i).Wait(status);
    
    // trying to receive EntierToRecv
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            EntierToRecv(i).Reallocate(nrecv_int(i));
            comm.Recv(EntierToRecv(i).GetData(), nrecv_int(i),
                      MPI::INTEGER, i, 5, status);
          }
      }

    // waiting for sending of EntierToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (request_active(i))
	request(i+nb_proc).Wait(status);
        
    // trying to receive FloatToRecv
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            int nsend_float = EntierToRecv(i)(0);
            
            FloatToRecv(i).Reallocate(nsend_float);
            MpiRecv(comm, FloatToRecv(i), FloatToRecv_tmp,
		    nsend_float, i, 6, status);
          }
      }
    
    // waiting for sending of FloatToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (request_active(i))
	request(i + 2*nb_proc).Wait(status);
    
    // deleting arrays
    request_active.Fill(false);
    for (int i = 0; i < nb_proc; i++)
      {
        nsend_int(i) = 0;
        EntierToSend(i).Clear();
        FloatToSend(i).Clear();
        FloatToSend_tmp(i).Clear();
      }

    // constructing array to know if a column is overlapped 
    // (treated by another processor)
    IVect OverlappedCol(n); OverlappedCol.Fill(-1);
    for (int i = 0; i < OverlapRowNumber.GetM(); i++)
      OverlappedCol(OverlapRowNumber(i)) = i;
    
    // assembling matrix B with interactions coming from other processors
    for (int i = 0; i < nb_proc; i++)
      {       
        if (nrecv_int(i) > 0)
          {	
            int ncol = EntierToRecv(i)(1);    
            // loop over columns
            nrecv_int(i) = 2; int nsend_float = 0;
            for (int j = 0; j < ncol; j++)
              {
                int iglob = EntierToRecv(i)(nrecv_int(i)++);
                int icol = Glob_to_local(iglob);
                int size_col = EntierToRecv(i)(nrecv_int(i)++);
                IVect index(size_col); Vector<T> values(size_col);
                for (int k = 0; k < size_col; k++)
                  {
                    index(k) = EntierToRecv(i)(nrecv_int(i)++);
                    values(k) = FloatToRecv(i)(nsend_float++);
                  }
                
                // adding to matrix B if the column is not shared
                // otherwise we send the column to the original processor
                if (OverlappedCol(icol) == -1)
                  {
                    if (size_col == 1)
                      B.AddInteraction(index(0), icol, values(0));
                    else
                      B.AddInteractionColumn(icol, size_col, index, values);
                  }
                else
                  {
                    int icol_ = OverlappedCol(icol);
                    int proc = OverlapProcNumber(icol_);
                    
                    int offset_int(2), offset_float(0);
                    if (nsend_int(proc) == 0)
                      {
                        nsend_int(proc) = 2;
                        EntierToSend(proc).Reallocate(size_col+4);
                        FloatToSend(proc).Reallocate(size_col);
                        EntierToSend(proc)(0) = 0;
                        EntierToSend(proc)(1) = 0;
                      }
                    else
                      {
                        offset_int = EntierToSend(proc).GetM();
                        offset_float = FloatToSend(proc).GetM();
                        EntierToSend(proc).Resize(size_col+2+offset_int);
                        FloatToSend(proc).Resize(size_col+offset_float);
                      }
                    
                    nsend_int(proc) += size_col+2;
                    EntierToSend(proc)(0) += size_col;
                    EntierToSend(proc)(1)++;
                    EntierToSend(proc)(offset_int++) = iglob;
                    EntierToSend(proc)(offset_int++) = size_col;
                    for (int k = 0; k < size_col; k++)
                      {
                        EntierToSend(proc)(offset_int++) = index(k);
                        FloatToSend(proc)(offset_float++) = values(k);
                      }                    
                  }
              }
          }
      }

    // sending columns overlapped to original processor
    // this case occurs only if distant columns contain overlapped columns
    for (int i = 0; i < nb_proc; i++)
      {
        if (i != rank)
          {
            request(i) = comm.Isend(&nsend_int(i), 1, MPI::INTEGER, i, 4);
            if (nsend_int(i) > 0)
              {
                request_active(i) = true;
                
                request(i+nb_proc) = 
                  comm.Isend(EntierToSend(i).GetData(),
			     nsend_int(i), MPI::INTEGER, i, 5);
                
                request(i+2*nb_proc) = 
                  MpiIsend(comm, FloatToSend(i), FloatToSend_tmp(i),
                           FloatToSend(i).GetM(), i, 6);            
              }
          }
      }
    
    // receiving columns
    for (int i = 0; i < nb_proc; i++)
      if (i != rank)
        comm.Recv(&nrecv_int(i), 1, MPI::INTEGER, i, 4, status);
    
    // trying to receive EntierToRecv
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            EntierToRecv(i).Reallocate(nrecv_int(i));
            comm.Recv(EntierToRecv(i).GetData(), nrecv_int(i),
                      MPI::INTEGER, i, 5, status);
          }
      }

    // waiting for sending of EntierToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (request_active(i))
	request(i+nb_proc).Wait(status);
        
    // trying to receive FloatToRecv
    for (int i = 0; i < nb_proc; i++)
      {
        if (nrecv_int(i) > 0)
          {
            int nsend_float = EntierToRecv(i)(0);
            
            FloatToRecv(i).Reallocate(nsend_float);
            MpiRecv(comm, FloatToRecv(i), FloatToRecv_tmp,
		    nsend_float, i, 6, status);
          }
      }
    
    // waiting for sending of FloatToSend effective
    for (int i = 0; i < nb_proc; i++)
      if (request_active(i))
	request(i + 2*nb_proc).Wait(status);
    
    // deleting arrays
    request_active.Fill(false);
    for (int i = 0; i < nb_proc; i++)
      {
        nsend_int(i) = 0;
        EntierToSend(i).Clear();
        FloatToSend(i).Clear();
      }

    // assembling matrix B with last interactions coming from other processors
    for (int i = 0; i < nb_proc; i++)
      {       
        if (nrecv_int(i) > 0)
          {	
            int ncol = EntierToRecv(i)(1);    
            // loop over columns
            nrecv_int(i) = 2; int nsend_float = 0;
            for (int j = 0; j < ncol; j++)
              {
                int iglob = EntierToRecv(i)(nrecv_int(i)++);
                int icol = Glob_to_local(iglob);
                int size_col = EntierToRecv(i)(nrecv_int(i)++);
                IVect index(size_col); Vector<T> values(size_col);
                for (int k = 0; k < size_col; k++)
                  {
                    index(k) = EntierToRecv(i)(nrecv_int(i)++);
                    values(k) = FloatToRecv(i)(nsend_float++);
                  }
                
                // adding to matrix B if the column is not shared
                // otherwise an error is generated (impossible case)
                if (OverlappedCol(icol) == -1)
                  {
                    if (size_col == 1)
                      B.AddInteraction(index(0), icol, values(0));
                    else
                      B.AddInteractionColumn(icol, size_col, index, values);
                  }
                else
                  {
                    cout << "Case impossible, this column is not overlapped"
			 << endl;
                    abort();
                  }
              }
          }
      }    
    
    // now we convert B to ColSparse while removing overlapped columns
    Glob_to_local.Clear();
    
    int nloc = n - OverlapRowNumber.GetM();
    PtrA.Reallocate(nloc+1);
    col_numbers.Reallocate(nloc);
    int ncol = 0, nnz = 0;
    PtrA(ncol) = 0;
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
        {
          col_numbers(ncol) = RowNumber(i);
          PtrA(ncol+1) = PtrA(ncol) + B.GetColumnSize(i);
          ncol++;
          nnz += B.GetColumnSize(i);
        }
    
    IndA.Reallocate(nnz);
    ValA.Reallocate(nnz);
    ncol = 0; nnz = 0;
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
        {
          for (int j = 0; j < B.GetColumnSize(i); j++)
            {
              IndA(nnz) = B.Index(i, j);
              ValA(nnz) = B.Value(i, j);
              nnz++;
            }
          
          ncol++;
        }
    
    local_col_numbers.Reallocate(nloc);
    ncol = 0;
    for (int i = 0; i < n; i++)
      if (OverlappedCol(i) == -1)
	local_col_numbers(ncol++) = i;
    
  }

  
  //! clears distant numbers in a sparse matrix
  template<class TypeDist>
  void EraseDistantEntries(MPI::Comm& comm, const Vector<bool>& IsRowDropped,
                           const Vector<bool>& IsRowDroppedDistant,
                           TypeDist& dist_row, Vector<IVect>& proc_row,
                           TypeDist& dist_col, Vector<IVect>& proc_col)
  {
    typedef typename TypeDist::value_type Vect1;
    typedef typename Vect1::value_type T1;
    for (int i = 0; i < dist_col.GetM(); i++)
      if (IsRowDropped(i))
        {
          dist_col(i).Clear();
          proc_col(i).Clear();
        }

    for (int j = 0; j < dist_row.GetM(); j++)
      {
        int nb = 0;
        for (int iloc = 0; iloc < dist_row(j).GetM(); iloc++)
          if (IsRowDroppedDistant(dist_row(j).Index(iloc)))
            nb++;
        
        if (nb > 0)
          {
            int size = dist_row(j).GetM();
            IVect row_num(size), proc(size); Vector<T1> val(size);
            for (int iloc = 0; iloc < dist_row(j).GetM(); iloc++)
              {
                row_num(iloc) = dist_row(j).Index(iloc);
                val(iloc) = dist_row(j).Value(iloc);
                proc(iloc) = proc_row(j)(iloc);
              }
            
            dist_row(j).Reallocate(size-nb);
            proc_row(j).Reallocate(size-nb);
            nb = 0;
            for (int iloc = 0; iloc < size; iloc++)
              if (!IsRowDroppedDistant(row_num(iloc)))
                {
                  dist_row(j).Index(nb) = row_num(iloc);
                  dist_row(j).Value(nb) = val(iloc);
                  proc_row(j)(nb) = proc(iloc);
                  nb++;
                }
          }
      }   
  }
  
  
  //! clears several columns of a sparse matrix
  /*!
    \param[in] num numbers of the columns to be cleared
    \param[in] comm MPI communicator
    \param[inout] A sparse matrix where columns are erased
    In this functions the column numbers are assumed to be increasing
   */  
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void EraseCol(const IVect& num,
		DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    // erasing columns of the local entries
    EraseCol(num, static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(A));

    MPI::Comm& comm = A.GetCommunicator();

    if (comm.Get_size() == 1)
      return;
    
    if (!A.IsReadyForMltAdd())
      A.PrepareMltAdd();
    
    Vector<bool> IsColDropped(A.GetM());
    IsColDropped.Fill(false);
    for (int i = 0; i < num.GetM(); i++)
      IsColDropped(num(i)) = true;
    
    Vector<bool> IsColDroppedDistant;
    A.ScatterColValues(IsColDropped, IsColDroppedDistant);
    
    EraseDistantEntries(comm, IsColDropped, IsColDroppedDistant,
                        A.dist_col, A.proc_col, A.dist_row, A.proc_row); 

    if (IsSymmetricMatrix(A))
      {
        A.ScatterRowValues(IsColDropped, IsColDroppedDistant);
        
        EraseDistantEntries(comm, IsColDropped, IsColDroppedDistant,
                            A.dist_row, A.proc_row, A.dist_col, A.proc_col); 
      }
  }
  

  //! clears several rows of a sparse matrix
  /*!
    \param[in] num numbers of the rows to be cleared
    \param[in] comm MPI communicator
    \param[inout] A sparse matrix where rows are erased
    In this functions the row numbers are assumed to be increasing
   */  
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void EraseRow(const IVect& num,
		DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    // erasing rows of the local entries
    EraseRow(num, static_cast<Matrix<T1, Prop1, Storage1, Allocator1>& >(A));
    
    MPI::Comm& comm = A.GetCommunicator();    
    if (comm.Get_size() == 1)
      return;
    
    if (!A.IsReadyForMltAdd())
      A.PrepareMltAdd();
    
    Vector<bool> IsRowDropped(A.GetM());
    IsRowDropped.Fill(false);
    for (int i = 0; i < num.GetM(); i++)
      IsRowDropped(num(i)) = true;
    
    Vector<bool> IsRowDroppedDistant;
    A.ScatterRowValues(IsRowDropped, IsRowDroppedDistant);
    
    EraseDistantEntries(comm, IsRowDropped, IsRowDroppedDistant,
                        A.dist_row, A.proc_row, A.dist_col, A.proc_col); 

    if (IsSymmetricMatrix(A))
      {
        A.ScatterColValues(IsRowDropped, IsRowDroppedDistant);
        
        EraseDistantEntries(comm, IsRowDropped, IsRowDroppedDistant,
                            A.dist_col, A.proc_col, A.dist_row, A.proc_row); 
      }

  }


  //! extracts a sub-matrix B from matrix A
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void
  CopySubMatrix(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		const IVect& row, const IVect& col,
		DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B)
  {
    CopySubMatrix(static_cast<const Matrix<T0, Prop0,
		  Storage0, Allocator0>& >(A),
                  row, col, static_cast<Matrix<T1, Prop1,
		  Storage1, Allocator1>& >(B));
    
    MPI::Comm& comm = B.GetCommunicator();
    if (comm.Get_size() == 1)
      return;
    
    if (!A.IsReadyForMltAdd())
      const_cast<DistributedMatrix<T0, Prop0, Storage0, Allocator0>& >(A)
        .PrepareMltAdd();
    
    B.Init(A);
    int m = A.GetM(), n = A.GetN();    
    Vector<bool> RowKept(m), ColKept(n);
    RowKept.Fill(false); ColKept.Fill(false);
    for (int i = 0; i < row.GetM(); i++)
      RowKept(row(i)) = true;
    
    for (int i = 0; i < col.GetM(); i++)
      ColKept(col(i)) = true;
    
    // checking consistency of row/col with symmetry of matrix B
    if (IsSymmetricMatrix(B))
      {
        if (m != n)
          {
            cout << "A is non-symmetric while B is symmetric" << endl;
            abort();
          }
        
        for (int i = 0; i < m; i++)
          if (RowKept(i) ^ ColKept(i))
            {
              cout << "row and col must be identic to obtain "
		   << "a symmetric matrix" << endl;
              abort();
            }
      }

    Vector<bool> RowKeptDistant, ColKeptDistant;
    A.ScatterRowValues(RowKept, RowKeptDistant);
    A.ScatterColValues(ColKept, ColKeptDistant);
    
    // using global numbers for B
    B.SwitchToGlobalNumbers();
    
    // extracting values of dist_col
    B.dist_col.Reallocate(m);
    B.proc_col.Reallocate(m);
    for (int i = 0; i < A.dist_col.GetM(); i++)
      {
        if (RowKept(i))
          {
            int size_row = 0;
            for (int j = 0; j < A.dist_col(i).GetM(); j++)
              if (ColKeptDistant(A.dist_col(i).Index(j)))
                size_row++;
            
            B.dist_col(i).Reallocate(size_row);
            B.proc_col(i).Reallocate(size_row);
            size_row = 0;
            for (int j = 0; j < A.dist_col(i).GetM(); j++)
              if (ColKeptDistant(A.dist_col(i).Index(j)))
                {
                  B.dist_col(i).Value(size_row) = A.dist_col(i).Value(j);
                  B.dist_col(i).Index(size_row)
                    = A.global_col_to_recv(A.dist_col(i).Index(j));
                  
                  B.proc_col(i)(size_row) = A.proc_col(i)(j);
                  size_row++;
                }
          }
        else
          {
            B.dist_col(i).Clear();
            B.proc_col(i).Clear();
          }
      }
    
    // same stuff for distant rows
    B.dist_row.Reallocate(n);
    B.proc_row.Reallocate(n);
    for (int i = 0; i < A.dist_row.GetM(); i++)
      {
        if (ColKept(i))
          {
            int size_col = 0;
            for (int j = 0; j < A.dist_row(i).GetM(); j++)
              if (RowKeptDistant(A.dist_row(i).Index(j)))
                size_col++;
            
            B.dist_row(i).Reallocate(size_col);
            B.proc_row(i).Reallocate(size_col);
            size_col = 0;
            for (int j = 0; j < A.dist_row(i).GetM(); j++)
              if (RowKeptDistant(A.dist_row(i).Index(j)))
                {
                  B.dist_row(i).Value(size_col) = A.dist_row(i).Value(j);
                  B.dist_row(i).Index(size_col)
                    = A.global_row_to_recv(A.dist_row(i).Index(j));
                  
                  B.proc_row(i)(size_col) = A.proc_row(i)(j);
                  size_col++;
                }
          }
        else
          {
            B.dist_row(i).Clear();
            B.proc_row(i).Clear();
          }

      }
  }
  
  
}

#define SELDON_FILE_DISTRIBUTED_MATRIX_CXX
#endif


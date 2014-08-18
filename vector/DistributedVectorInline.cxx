#ifndef SELDON_FILE_DISTRIBUTED_VECTOR_INLINE_CXX

#include "DistributedVector.hxx"

namespace Seldon
{
  
  //! Constructor taking overlapped rows
  /*!
    In the array rows, you store all the row numbers
    which are already counted in another processor.
    For example :
    if processor 0 contains rows (0, 3, 5, 7, 8)
    if processor 1 contains rows (0, 1, 2, 4, 5, 6)
    Then you can set OverlapRowNumbers empty on processor 0,
    and equal to (0, 5) on processor 1.
    If there is no row shared by processors, the array rows
    will be empty for each processor
  */
  template<class T, class Allocator>
  inline DistributedVector<T, Allocator>
  ::DistributedVector(const IVect& rows, const MPI::Comm& comm)
    : OverlapRowNumbers(rows), comm_(comm)
  {
  }
  
  
  //! Copy constructor
  template<class T, class Allocator>
  inline DistributedVector<T, Allocator>::
  DistributedVector(const DistributedVector<T, Allocator>& V)
    : OverlapRowNumbers(V.OverlapRowNumbers), comm_(V.comm_)
  {
    Vector<T, Vect_Full, Allocator>::
      Copy(static_cast<const Vector<T, Vect_Full, Allocator>& >(V));
  }
  
  
  //! setting pointer
  template<class T, class Allocator>
  inline void DistributedVector<T, Allocator>::
  SetData(Vector<T, Vect_Full, Allocator>& x)
  {
    Vector<T, Vect_Full, Allocator>::SetData(x.GetM(), x.GetData());
  }
  
  
  //! setting pointer
  template<class T, class Allocator>
  inline void DistributedVector<T, Allocator>::SetData(int n, T* data)
  {
    Vector<T, Vect_Full, Allocator>::SetData(n, data);
  }
  
  
  //! operator =
  template<class T, class Allocator>
  inline DistributedVector<T, Allocator>& DistributedVector<T, Allocator>::
  operator=(const DistributedVector<T, Allocator>& X)
  {
    Vector<T>::Copy(static_cast<const Vector<T, Vect_Full, Allocator>& >(X));
    return *this;
  }
  
  
  //! returns the number of rows already counted 
  template<class T, class Allocator>
  inline int DistributedVector<T, Allocator>::GetNbOverlap() const
  {
    return OverlapRowNumbers.GetM();
  }
  
  
  //! returns an overlapped row number
  template<class T, class Allocator>
  inline int DistributedVector<T, Allocator>::GetOverlapRow(int i) const
  {
    return OverlapRowNumbers(i);
  }
  
  
  //! returns communicator
  template<class T, class Allocator>
  inline const MPI::Comm& DistributedVector<T, Allocator>
  ::GetCommunicator() const
  {
    return comm_;
  }
  
}

#define SELDON_FILE_DISTRIBUTED_VECTOR_INLINE_CXX
#endif


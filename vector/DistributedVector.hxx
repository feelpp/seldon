#ifndef SELDON_FILE_DISTRIBUTED_VECTOR_HXX

namespace Seldon
{
  
  //! class storing a vector distributed over all the processors
  /*!
    This class is useful when some rows of the vector
    are shared by several processors. Functions DotProd, DotProdConj
    and Norm2 are overloaded in order to take into account this overlap.
   */
  template<class T, class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class DistributedVector : public Vector<T, Vect_Full, Allocator>
  {
  protected :
    //! row numbers shared with other processors
    /*!
      In this array, you store all the row numbers
      which are already counted in another processor.
      For example :
      if processor 0 contains rows (0, 3, 5, 7, 8)
      if processor 1 contains rows (0, 1, 2, 4, 5, 6)
      Then you can set OverlapRowNumbers empty on processor 0,
      and equal to (0, 5) on processor 1
     */
    const IVect& OverlapRowNumbers;    
    //! MPI communicator grouping processors involved in the computation
    const MPI::Comm& comm_;
        
  public :
    // constructors
    DistributedVector(const IVect& rows, const MPI::Comm& comm);
    DistributedVector(const DistributedVector<T, Allocator>& V);
    
    // memory management
    void SetData(Vector<T, Vect_Full, Allocator>& x);
    void SetData(int n, T* data);
    
    // basic functions
    DistributedVector<T, Allocator>&
    operator=(const DistributedVector<T, Allocator>& X);
    
    int GetNbOverlap() const;
    int GetOverlapRow(int i) const;
    const MPI::Comm& GetCommunicator() const;
    
  };
  
  
  // returns X.Y
  template<class T1, class Allocator1>
  T1 DotProd(const DistributedVector<T1, Allocator1>& X,
	     const DistributedVector<T1, Allocator1>& Y);
  
  // returns X' . Y
  template<class T1, class Allocator1>
  T1 DotProdConj(const DistributedVector<T1, Allocator1>& X,
		 const DistributedVector<T1, Allocator1>& Y);  
  
  // returns euclidian norm of x
  template<class T, class Allocator>
  T Norm2(const DistributedVector<complex<T>, Allocator>& x);
  
  // returns euclidian norm of x
  template<class T, class Allocator>
  T Norm2(const DistributedVector<T, Allocator>& x);

}

#define SELDON_FILE_DISTRIBUTED_VECTOR_HXX
#endif


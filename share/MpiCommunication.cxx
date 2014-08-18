#ifndef SELDON_FILE_MPI_COMMUNICATION_CXX

namespace Seldon
{

  inline const MPI::Datatype& GetMpiDataType(const Vector<bool>&)
  {
    return MPI::BOOL;
  }

  inline const MPI::Datatype& GetMpiDataType(const Vector<int>&)
  {
    return MPI::INTEGER;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const Vector<float>&)
  {
    return MPI::FLOAT;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const Vector<complex<float> >&)
  {
    return MPI::FLOAT;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const Vector<double>&)
  {
    return MPI::DOUBLE;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const Vector<complex<double> >&)
  {
    return MPI::DOUBLE;
  }

  inline const MPI::Datatype& GetMpiDataType(const int&)
  {
    return MPI::INTEGER;
  }

  inline const MPI::Datatype& GetMpiDataType(const bool&)
  {
    return MPI::BOOL;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const float&)
  {
    return MPI::FLOAT;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const complex<float>&)
  {
    return MPI::FLOAT;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const double&)
  {
    return MPI::DOUBLE;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const complex<double>&)
  {
    return MPI::DOUBLE;
  }
  
  template<class T>
  int GetRatioMpiDataType(const Vector<T>&)
  {
    return 1;
  }
  
  template<class T>
  int GetRatioMpiDataType(const Vector<complex<T> >&)
  {
    return 2;
  }

  template<class T>
  MPI::Request MpiIsend(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    Vector<T> xvec;
    xvec.SetData(n, x);
    
    MPI::Request rq = MpiIsend(comm, xvec, xtmp, n, proc, tag);
    
    xvec.Nullify();
    return rq;
  }
  
  template<class T>
  MPI::Request MpiIsend(const MPI::Comm& comm, Vector<T>& x,
			Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    return comm.Isend(x.GetData(), n*GetRatioMpiDataType(x),
		      GetMpiDataType(x), proc, tag); 
  }

  template<class T>
  MPI::Request MpiIrecv(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    Vector<T> xvec;
    xvec.SetData(n, x);
    
    MPI::Request rq = MpiIrecv(comm, xvec, xtmp, n, proc, tag);
    
    xvec.Nullify();
    return rq;
  }

  template<class T>
  MPI::Request MpiIrecv(const MPI::Comm& comm, Vector<T>& x,
			Vector<int64_t>& xtmp,
                        int n, int proc, int tag)
  {
    return comm.Irecv(x.GetData(), n*GetRatioMpiDataType(x),
		      GetMpiDataType(x), proc, tag); 
  }
  
  template<class T>
  void MpiCompleteIrecv(T* x, Vector<int64_t>& xtmp, int n)
  {
    Vector<T> xvec;
    xvec.SetData(n, x);
    
    MpiCompleteIrecv(xvec, xtmp, n);
    
    xvec.Nullify();
  }
  
  template<class T>
  void MpiCompleteIrecv(Vector<T>& x, Vector<int64_t>& xtmp, int n)
  {    
  }

  template<class T>
  void MpiSsend(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    Vector<T> xvec;
    xvec.SetData(n, x);
    
    MpiSsend(comm, xvec, xtmp, n, proc, tag);
    
    xvec.Nullify();
  }
  
  template<class T>
  void MpiSsend(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    comm.Ssend(x.GetData(), n*GetRatioMpiDataType(x),
	       GetMpiDataType(x), proc, tag);
  }
                       
  template<class T>
  void MpiSend(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                int n, int proc, int tag)
  {
    Vector<T> xvec;
    xvec.SetData(n, x);
    
    MpiSend(comm, xvec, xtmp, n, proc, tag);
    
    xvec.Nullify();
  }
  
  template<class T>
  void MpiSend(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag)
  {
    comm.Send(x.GetData(), n*GetRatioMpiDataType(x),
	      GetMpiDataType(x), proc, tag);
  }
                     
  template<class T>
  void MpiGather(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                 T* y, int n, int proc)
  {
    Vector<T> xvec, yvec;
    xvec.SetData(n, x);
    if (proc == comm.Get_rank())
      yvec.SetData(n*comm.Get_size(), y);
    
    MpiGather(comm, xvec, xtmp, yvec, n, proc);
    
    xvec.Nullify();
    if (proc == comm.Get_rank())
      yvec.Nullify();
  }
  
  template<class T>
  void MpiGather(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                 Vector<T>& y, int n, int proc)
  {
    comm.Gather(x.GetData(), n*GetRatioMpiDataType(x), GetMpiDataType(x),
                y.GetData(), n*GetRatioMpiDataType(x), GetMpiDataType(x),
		proc);
  }
  
  template<class T>
  void MpiAllreduce(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                    T* y, int n, const MPI::Op& op)
  {
    Vector<T> xvec, yvec;
    xvec.SetData(n, x);
    yvec.SetData(n, y);
    
    MpiAllreduce(comm, xvec, xtmp, yvec, n, op);
    
    xvec.Nullify();
    yvec.Nullify();
  }
  
  template<class T>
  void MpiAllreduce(const MPI::Comm& comm, Vector<T>& x,
		    Vector<int64_t>& xtmp,
                    Vector<T>& y, int n, const MPI::Op& op)
  {
    comm.Allreduce(x.GetData(), y.GetData(), n*GetRatioMpiDataType(x),
                   GetMpiDataType(x), op);
  }
  
  template<class T>
  void MpiReduce(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
                 T* y, int n, const MPI::Op& op, int proc)
  {
    Vector<T> xvec, yvec;
    xvec.SetData(n, x);
    if (proc == comm.Get_rank())
      yvec.SetData(n, y);
    
    MpiReduce(comm, xvec, xtmp, yvec, n, op, proc);
    
    xvec.Nullify();
    if (proc == comm.Get_rank())
      yvec.Nullify();
  }
  
  template<class T>
  void MpiReduce(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
                 Vector<T>& y, int n, const MPI::Op& op, int proc)
  {
    comm.Reduce(x.GetData(), y.GetData(), n*GetRatioMpiDataType(x),
                GetMpiDataType(x), op, proc);
  }
  
  template<class T>
  void MpiRecv(const MPI::Comm& comm, T* x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI::Status& status)
  {
    Vector<T> xvec;
    xvec.SetData(n, x);
    
    MpiRecv(comm, xvec, xtmp, n, proc, tag, status);
    
    xvec.Nullify();
  }
  
  template<class T>
  void MpiRecv(const MPI::Comm& comm, Vector<T>& x, Vector<int64_t>& xtmp,
               int n, int proc, int tag, MPI::Status& status)
  {
    comm.Recv(x.GetData(), n*GetRatioMpiDataType(x), GetMpiDataType(x),
              proc, tag, status);
  }

  template<class T>
  void MpiBcast(const MPI::Comm& comm, T* x,
		Vector<int64_t>& xtmp, int n, int proc)
  {
    Vector<T> xvec;
    xvec.SetData(n, x);
    
    MpiBcast(comm, xvec, xtmp, n, proc);
    
    xvec.Nullify();
  }
  
  template<class T>
  void MpiBcast(const MPI::Comm& comm, Vector<T>& x,
		Vector<int64_t>& xtmp, int n, int proc)
  {
    comm.Bcast(x.GetData(), n, GetMpiDataType(x), proc);
  }
  
}

#define SELDON_FILE_MPI_COMMUNICATION_CXX
#endif

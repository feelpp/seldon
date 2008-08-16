// Copyright (C) 2001-2003 Vivien Mallet
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

// To be included by Seldon.hxx

#ifndef FILE_MATRIX_ARRAY_SPARSE_HXX

namespace Seldon
{
  
  class ArrayRowSparse : public RowSparse
  {
  };
  
  class ArrayColSparse : public ColSparse
  {
  };
  
  class ArrayRowSymSparse : public RowSymSparse
  {
  };
  
  class ArrayColSymSparse : public ColSymSparse
  {
  };
  
  //! Sparse Array-matrix class.
  /*!
    Sparse matrices are defined by: (1) the number of rows and columns;
    (2) the number of non-zero entries; (3) an array of vectors ind
    ind(i) is a vector, which contains indices of columns of the row i
    (4) an array of vectors val : val(i) is a vector, which contains values of
    the row i
  */
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_ArraySparse
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef Vector<T,Vect_Full,Allocator> vect_value;
    typedef Vector<T,Vect_Full,Allocator>* vect_value_ptr;
    
    // Static attributes.
  protected:
    static Allocator allocator_;

    // Attributes.
  protected:
    // Number of rows.
    int m_;
    // Number of columns.
    int n_;
    // Number of non-zero elements.
    int nz_;
    Vector<IVect,Vect_Full,NewAlloc<IVect> > ind;
    Vector<vect_value,Vect_Full,NewAlloc<vect_value> > val;
    
    // Methods.
  protected :
    // Memory management.
    void ClearVector(int i);
    void ReallocateVector(int i,int j);
    void ResizeVector(int i,int j);
    void SwapVector(int i,int i_);
    void ReplaceIndexVector(int i,IVect& new_index);
    
    int GetVectorSize(int i) const;
    void PrintVector(int i) const;
    void AssembleVector(int i);
    
    void AddInteraction(int i, int j, const T& a);
    
    template<class Storage1,class Allocator1>
    void AddInteractionVector(int i, int nb_interac, IVect col_interac,
			      Vector<T, Storage1, Allocator1> val_interac);
    
  public:
    // Constructors.
    Matrix_ArraySparse();
    Matrix_ArraySparse(int i, int j);
    
    // Destructor.
    ~Matrix_ArraySparse();
    void Clear();
    
    // Memory management.
    int GetMemorySize() const;
    void Reallocate(int i,int j);
    void Resize(int i,int j);
    
    // Basic methods.
    int GetM() const;
    int GetN() const;
    int GetNonZeros() const;
    int GetDataSize() const;
    int* GetInd(int i) const;
    pointer GetData(int i) const;
    
    IVect* GetInd() const;
    vect_value* GetData() const;
    
    // Element acess and affectation.
    value_type operator() (int i, int j) const;
    reference operator() (int i, int j);
    
    const_reference Value(int num_row, int i) const;
    reference Value(int num_row, int i);
    int Index(int num_row, int i) const;
    int& Index(int num_row, int i);
    
    void SetData(int, int, int, IVect*, vect_value*);
    void Nullify();
    
    // Convenient functions.
    void Print() const;
    void Assemble();
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);
    
    void SetIdentity();
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    template <class T0>
    Matrix_ArraySparse<T, Prop, Storage, Allocator>& operator= (const T0& x);
    void FillRand();
    
    template <class T_, class Prop_, class Storage_, class Allocator_>
    friend ostream&
    operator <<(ostream& out,
		const Matrix_ArraySparse<T_, Prop_, Storage_, Allocator_>& A);
  };
  
  //! Matrix allocator.
  template <class T, class Prop, class Storage, class Allocator>
  Allocator Matrix_ArraySparse<T, Prop, Storage, Allocator>::allocator_;
  
  
  //! Column-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSparse, Allocator> :
    public Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>
  {
  public:
    Matrix()  throw();
    Matrix(int i, int j);
    
    // Memory management.
    void ClearColumn(int i);
    void ReallocateColumn(int i, int j);
    void ResizeColumn(int i, int j);
    void SwapColumn(int i, int i_);
    void ReplaceIndexColumn(int i, IVect& new_index);
    
    int GetColumnSize(int i) const;
    void PrintColumn(int i) const;
    void AssembleColumn(int i);
    
    void AddInteraction(int i, int j, const T& val);
    
    void AddInteractionRow(int, int, int*, T*);
    void AddInteractionColumn(int, int, int*, T*);
    
    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<T, Vect_Full, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<T, Vect_Full, Alloc1>& val);
  };
  
  
  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSparse, Allocator> :
    public Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>
  {
  public:
    Matrix()  throw();
    Matrix(int i, int j);
    
    // Memory management.
    void ClearRow(int i);
    void ReallocateRow(int i, int j);
    void ResizeRow(int i, int j);
    void SwapRow(int i, int i_);
    void ReplaceIndexRow(int i, IVect& new_index);
    
    int GetRowSize(int i) const;
    void PrintRow(int i) const;
    void AssembleRow(int i);
    
    void AddInteraction(int i, int j, const T& val);
    
    void AddInteractionRow(int, int, int*, T*);
    void AddInteractionColumn(int, int, int*, T*);
    
    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<T, Vect_Full, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<T, Vect_Full, Alloc1>& val);
  };
  
  //! Column-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSymSparse, Allocator>:
    public Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>
  {
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::reference reference;
    
    Matrix()  throw();
    Matrix(int i, int j);
    
    value_type operator() (int i, int j) const;
    reference operator() (int i, int j);
    
    // Memory management.
    void ClearColumn(int i);
    void ReallocateColumn(int i, int j);
    void ResizeColumn(int i, int j);
    void SwapColumn(int i, int i_);
    void ReplaceIndexColumn(int i, IVect& new_index);
    
    int GetColumnSize(int i) const;
    void PrintColumn(int i) const;
    void AssembleColumn(int i);
    
    void AddInteraction(int i, int j, const T& val);
    
    void AddInteractionRow(int, int, int*, T*);
    void AddInteractionColumn(int, int, int*, T*);
    
    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<T, Vect_Full, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<T, Vect_Full, Alloc1>& val);
  };
  
  
  //! Row-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymSparse, Allocator>:
    public Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>
  {
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::reference reference;
    
    Matrix()  throw();
    Matrix(int i, int j);
    
    value_type operator() (int i, int j) const;
    reference operator() (int i, int j);
    
    // Memory management.
    void ClearRow(int i);
    void ReallocateRow(int i, int j);
    void ResizeRow(int i, int j);
    void SwapRow(int i, int i_);
    void ReplaceIndexRow(int i, IVect& new_index);
    
    int GetRowSize(int i) const;
    void PrintRow(int i) const;
    void AssembleRow(int i);
    
    void AddInteraction(int i, int j, const T& val);
    
    void AddInteractionRow(int, int, int*, T*);
    void AddInteractionColumn(int, int, int*, T*);
    
    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<T, Vect_Full, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<T, Vect_Full, Alloc1>& val);
  };
  
} // namespace Seldon

#define FILE_MATRIX_ARRAY_SPARSE_HXX
#endif

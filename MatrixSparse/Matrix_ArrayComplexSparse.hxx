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

// To be included by Seldon.hxx

#ifndef FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX

namespace Seldon
{
  
  class ArrayRowComplexSparse : public RowComplexSparse
  {
  };
  
  
  class ArrayRowSymComplexSparse : public RowSymComplexSparse
  {
  };
  
  
  //! Sparse Array-matrix class.
  /*!
    Sparse matrices are defined by: (1) the number of rows and columns;
    (2) the number of non-zero entries; (3) an array of vectors ind
    ind(i) is a vector, which contains indices of columns of the row i
    (4) an array of vectors val : val(i) is a vector, which contains values of the row i
  */
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_ArrayComplexSparse
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef Vector<T,Vect_Full,Allocator> vect_value;
    
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
    int nz_real_, nz_imag_;
    Vector<IVect,Vect_Full,NewAlloc<IVect> > ind_real;
    Vector<vect_value,Vect_Full,NewAlloc<vect_value> > val_real;
    Vector<IVect,Vect_Full,NewAlloc<IVect> > ind_imag;
    Vector<vect_value,Vect_Full,NewAlloc<vect_value> > val_imag;
    
    // Methods.
  public:
    // Constructors.
    Matrix_ArrayComplexSparse();
    Matrix_ArrayComplexSparse(int i, int j);
    
    // Destructor.
    ~Matrix_ArrayComplexSparse();
    void Clear();
    void ClearRealRow(int i);
    void ClearImagRow(int i);
    
    // Memory management.
    int GetMemorySize() const;
    void Reallocate(int i,int j);
    void Resize(int i,int j);
    void ReallocateRealRow(int i,int j);
    void ReallocateImagRow(int i,int j);
    void ResizeRealRow(int i,int j);
    void ResizeImagRow(int i,int j);
    void SwapRealRow(int i,int i_);
    void SwapImagRow(int i,int i_);
    void ReplaceIndexRealRow(int i,IVect& new_index);
    void ReplaceIndexImagRow(int i,IVect& new_index);
    
    // Basic methods.
    int GetM() const;
    int GetN() const;
    int GetRealNonZeros() const;
    int GetImagNonZeros() const;
    int GetRealDataSize() const;
    int GetImagDataSize() const;
    int GetDataSize() const;
    int GetRealRowSize(int i) const;
    int GetImagRowSize(int i) const;
    int* GetRealInd(int i) const;
    int* GetImagInd(int i) const;
    pointer GetRealData(int i) const;
    pointer GetImagData(int i) const;
    
    // Element acess and affectation.
    complex<T> operator() (int i, int j) const;
    const_reference ValueReal(int num_row,int i) const;
    reference ValueReal(int num_row,int i);
    int IndexReal(int num_row,int i) const;
    int& IndexReal(int num_row,int i);
    const_reference ValueImag(int num_row,int i) const;
    reference ValueImag(int num_row,int i);
    int IndexImag(int num_row,int i) const;
    int& IndexImag(int num_row,int i);
    const IVect& GetRealIndexRow(int i) const;
    const IVect& GetImagIndexRow(int i) const;
    const vect_value& GetRealValueRow(int i) const;
    const vect_value& GetImagValueRow(int i) const;
    
    void AddRealInteraction(int i,int j,const T& a);
    void AddImagInteraction(int i,int j,const T& a);
    
    template<class Storage1,class Allocator1>
    void AddRealInteraction(int i,int nb_interac, IVect col_interac,
			    Vector<T,Storage1,Allocator1> val_interac);
    
    template<class Storage1,class Allocator1>
    void AddImagInteraction(int i,int nb_interac, IVect col_interac,
			    Vector<T,Storage1,Allocator1> val_interac);
    
    // Convenient functions.
    void Print() const;
    void PrintRealRow(int i) const;
    void PrintImagRow(int i) const;
    void Assemble();
    void AssembleRealRow(int i);
    void AssembleImagRow(int i);
    
    void SetIdentity();
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const complex<T0>& x);
    template <class T0>
    Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>& operator=
    (const complex<T0>& x);
    void FillRand();
    
    template <class T_, class Prop_, class Storage_, class Allocator_>
    friend ostream&
    operator <<(ostream& out, const Matrix_ArrayComplexSparse<T_, Prop_,
		Storage_, Allocator_>& A);
  };
  
  
  //! Matrix allocator.
  template <class T, class Prop, class Storage, class Allocator>
  Allocator Matrix_ArrayComplexSparse<T, Prop,
				      Storage, Allocator>::allocator_;
  
  
  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowComplexSparse, Allocator>
    : public Matrix_ArrayComplexSparse<T, Prop, ArrayRowComplexSparse,
				       Allocator>
  {
  public:
    Matrix()  throw();
    Matrix(int i, int j);
    
    void AddInteraction(int i, int j, const complex<T>& val);
    
    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<complex<T>, Vect_Full, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<complex<T>, Vect_Full, Alloc1>& val);
  };
  
  
  //! Row-major symmetric sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
    : public Matrix_ArrayComplexSparse<T, Prop, ArrayRowSymComplexSparse,
				       Allocator>
  {
  public:
    Matrix()  throw();
    Matrix(int i, int j);
    
    complex<T> operator() (int i, int j) const;
    
    void AddInteraction(int i, int j, const complex<T>& val);
    
    template<class Alloc1>
    void AddInteractionRow(int i, int nb, const IVect& col,
			   const Vector<complex<T>, Vect_Full, Alloc1>& val);
    template<class Alloc1>
    void AddInteractionColumn(int i, int nb, const IVect& row,
			      const Vector<complex<T>, Vect_Full, Alloc1>& val);
  };
  
} // namespace Seldon

#define FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX
#endif

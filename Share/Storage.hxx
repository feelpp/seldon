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

#ifndef SELDON_FILE_STORAGE_HXX

namespace Seldon
{


  //////////////////////
  // GENERAL MATRICES //
  //////////////////////


  class ColMajor
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowMajor
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };



  /////////////
  // VECTORS //
  /////////////


  class Vect_Full;
  class Vect_Sparse;


  ////////////
  // SPARSE //
  ////////////


  class ColSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class ColComplexSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowComplexSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class ColSymSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowSymSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class ColSymComplexSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowSymComplexSparse
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };

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
  
  class ArrayRowComplexSparse : public RowComplexSparse
  {
  };
  
  class ArrayRowSymComplexSparse : public RowSymComplexSparse
  {
  };
  
  class ArrayColComplexSparse : public ColComplexSparse
  {
  };
  
  class ArrayColSymComplexSparse : public ColSymComplexSparse
  {
  };
  
  
  ///////////////
  // SYMMETRIC //
  ///////////////


  class ColSymPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowSymPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class ColSym
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowSym
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };



  ///////////////
  // HERMITIAN //
  ///////////////


  class ColHerm
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowHerm
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class ColHermPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };


  class RowHermPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
  };



  ////////////////
  // TRIANGULAR //
  ////////////////


  class ColUpTriang
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static bool UpLo();
  };


  class ColLoTriang
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static bool UpLo();
  };


  class RowUpTriang
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static bool UpLo();
  };


  class RowLoTriang
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static bool UpLo();
  };


  class ColUpTriangPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static bool UpLo();
  };


  class ColLoTriangPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static bool UpLo();
  };


  class RowUpTriangPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static bool UpLo();
  };


  class RowLoTriangPacked
  {
  public:
    static int GetFirst(int i, int j);
    static int GetSecond(int i, int j);
    static bool UpLo();
  };


} // namespace Seldon.

#define SELDON_FILE_STORAGE_HXX
#endif

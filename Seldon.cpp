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


#ifndef SELDON_FILE_SELDON_CPP


#include "Seldon.hxx"

namespace Seldon
{
  template class MallocAlloc<int>;
  template class Vector_Base<int, MallocAlloc<int> >;
  template class Vector<int, Vect_Full, MallocAlloc<int> >;
  template class MallocAlloc<double>;
  template class Vector_Base<double, MallocAlloc<double> >;
  template class Vector<double, Vect_Full, MallocAlloc<double> >;

  template class Matrix_Base<int, MallocAlloc<int> >;
  template class Matrix_Pointers<int, General, RowMajor, MallocAlloc<int> >;
  template class Matrix<int, General, RowMajor, MallocAlloc<int> >;
  template class Matrix_Base<double, MallocAlloc<double> >;
  template class Matrix_Pointers<double, General, RowMajor, MallocAlloc<double> >;
  template class Matrix<double, General, RowMajor, MallocAlloc<double> >;
}


#define SELDON_FILE_SELDON_CPP
#endif

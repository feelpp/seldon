// Copyright (C) 2010 Vivien Mallet
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


#ifndef SELDON_FILE_COMPUTATION_SPARSESOLVER_HXX


namespace Seldon
{


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SparseSolveLU(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                     Vector<T1, Storage1, Allocator1>& Y);


  template <class T, class Prop0, class Allocator0, class Allocator1>
  void SolveLU(Matrix<T, Prop0, ColSparse, Allocator0>& M,
               Vector<T, VectFull, Allocator1>& Y);


  template <class T, class Prop0, class Allocator0, class Allocator1>
  void SolveLU(Matrix<T, Prop0, RowSparse, Allocator0>& M,
               Vector<T, VectFull, Allocator1>& Y);


}  // namespace Seldon.


#define SELDON_FILE_COMPUTATION_SPARSESOLVER_HXX
#endif

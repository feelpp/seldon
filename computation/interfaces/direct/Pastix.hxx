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

#ifndef SELDON_FILE_PASTIX_HXX

// including Pastix headers
extern "C"
{
#define _COMPLEX_H

#include "pastix.h"
}

namespace Seldon
{

  template<class T>
  class MatrixPastix
  {
  protected :
    //! pastix structure
    pastix_data_t* pastix_data;
    //! options (integers)
    int    iparm[64];
    //! options (floats)
    double dparm[64];
    //! number of columns
    int n;
    IVect perm, invp, col_num;
    bool distributed;
    int print_level;

  public :

    MatrixPastix();
    ~MatrixPastix();

    void Clear();

    void CallPastix(const MPI_Comm&, int* colptr, int* row, T* val,
                    T* b, int nrhs);
    void CheckMatrix(const MPI_Comm&, int**, int**, T**);

    void HideMessages();
    void ShowMessages();


    template<class Prop, class Storage, class Allocator>
    void FindOrdering(Matrix<T, Prop, Storage, Allocator> & mat,
		      IVect& numbers, bool keep_matrix = false);

    template<class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T, General, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T, Symmetric, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Allocator2>
    void Solve(Vector<T, VectFull, Allocator2>& x);

    template<class Allocator2, class Transpose_status>
    void Solve(const Transpose_status& TransA,
	       Vector<T, VectFull, Allocator2>& x);

#ifdef SELDON_WITH_MPI
    template<class Prop, class Allocator>
    void
    FactorizeDistributedMatrix(Matrix<T, General, ColSparse, Allocator>& A,
                               const Prop& sym, const IVect& glob_number,
                               bool keep_matrix = false);

    template<class Allocator2>
    void SolveDistributed(Vector<T, Vect_Full, Allocator2>& x,
                          const IVect& glob_num);

    template<class Allocator2, class Transpose_status>
    void SolveDistributed(const Transpose_status& TransA,
			  Vector<T, Vect_Full, Allocator2>& x,
                          const IVect& glob_num);
#endif

  };

}

#define SELDON_FILE_PASTIX_HXX
#endif


// Copyright (C) 2001-2010 Vivien Mallet
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


#include <cppunit/extensions/HelperMacros.h>

#include "Seldon.hxx"
#include "SeldonSolver.hxx"
using namespace Seldon;


class SparseMatrixTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(SparseMatrixTest);
  CPPUNIT_TEST(test_conversion);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;
  int Nelement_;
  int Nloop_;

public:
  void setUp()
  {
  }


  void tearDown()
  {
  }


  void test_conversion()
  {
    m_ = 25;
    n_ = 10;
    Nelement_ = 30;
    Nloop_ = 10;
    conversion();

    m_ = 50;
    n_ = 60;
    Nelement_ = 5;
    Nloop_ = 10;
    conversion();

    m_ = 50;
    n_ = 50;
    Nelement_ = 1;
    Nloop_ = 10;
    conversion();

    m_ = 10;
    n_ = 25;
    Nelement_ = 30;
    Nloop_ = 10;
    conversion();

    m_ = 10;
    n_ = 5;
    Nelement_ = 40;
    Nloop_ = 10;
    conversion();

    m_ = 5;
    n_ = 10;
    Nelement_ = 100;
    Nloop_ = 2;
    conversion();
  }


  void conversion()
  {
    srand(time(NULL));

    int i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        A_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;

        Copy(A_array, A);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < n_; j++)
            {
              CPPUNIT_ASSERT(A_full(i, j) == A(i, j));
              CPPUNIT_ASSERT(A_full(i, j) == A_array(i, j));
            }

        Vector<int> row_index, col_index;
        Vector<double> value;
        ConvertMatrix_to_Coordinates(A_array, row_index, col_index, value);

        ConvertMatrix_from_Coordinates(row_index, col_index, value, A);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < n_; j++)
            CPPUNIT_ASSERT(A_full(i, j) == A(i, j));
      }
  }
};

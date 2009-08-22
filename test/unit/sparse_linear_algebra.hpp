// Copyright (C) 2001-2009 Vivien Mallet
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

#include <cstdlib>
#include <ctime>

#include "Seldon.hxx"
#include "SeldonSolver.hxx"
using namespace Seldon;


class SparseLinearAlgebraTest: public CppUnit::TestFixture
{


  CPPUNIT_TEST_SUITE(SparseLinearAlgebraTest);
  CPPUNIT_TEST(test_mlt);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;
  int p_;
  int Nelement_;
  int Nloop_;

public:
  void setUp()
  {
  }


  void tearDown()
  {
  }


  void test_mlt()
  {
    Nloop_ = 20;

    m_ = 10;
    n_ = 25;
    p_ = 9;
    Nelement_ = 20;
    mlt();

    m_ = 100;
    n_ = 250;
    p_ = 9;
    Nelement_ = 20;
    mlt();

    m_ = 10;
    n_ = 10;
    p_ = 2000;
    Nelement_ = 200;
    mlt();
  }


  void mlt()
  {
    srand(time(NULL));

    int i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        Matrix<double> B_full(n_, p_);
        Matrix<double> C_full(m_, p_);
        A_full.Zero();
        B_full.Zero();
        C_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        Matrix<double, General, ArrayRowSparse> B_array(n_, p_);
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % n_;
            j = rand() % p_;
            value = double(rand());
            B_array.AddInteraction(i, j, value);
            B_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;
        Matrix<double, General, RowSparse> B;
        Matrix<double, General, RowSparse> C;

        Copy(A_array, A);
        Copy(B_array, B);

        Mlt(A, B, C);
        Mlt(A_full, B_full, C_full);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < p_; j++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(C_full(i, j), C(i, j),
                                         1.e-14 * C_full(i, j));
      }
  }


};

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


#include <cppunit/extensions/HelperMacros.h>

#include "Seldon.hxx"
using namespace Seldon;


class MatrixTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(MatrixTest);
  CPPUNIT_TEST(test_reallocate);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;

public:
  void setUp()
  {
    m_ = 10;
    n_ = 25;
  }

  
  void tearDown()
  {
  }


  void test_reallocate()
  {
    Matrix<double> M(m_, n_);
    M.Reallocate(0, 0);
    CPPUNIT_ASSERT(M.GetM() == 0);
    CPPUNIT_ASSERT(M.GetN() == 0);
    M.Reallocate(m_, n_);
    CPPUNIT_ASSERT(M.GetM() == m_);
    CPPUNIT_ASSERT(M.GetN() == n_);
    M.Reallocate(2 * m_, n_);
    CPPUNIT_ASSERT(M.GetM() == 2 * m_);
    CPPUNIT_ASSERT(M.GetN() == n_);
    M.Reallocate(2 * m_, 2 * n_);
    CPPUNIT_ASSERT(M.GetM() == 2 * m_);
    CPPUNIT_ASSERT(M.GetN() == 2 * n_);
    M.Reallocate(m_, 0);
    CPPUNIT_ASSERT(M.GetM() == m_);
    CPPUNIT_ASSERT(M.GetN() == 0);
    M.Reallocate(0, n_);
    CPPUNIT_ASSERT(M.GetM() == 0);
    CPPUNIT_ASSERT(M.GetN() == n_);
    M.Reallocate(0, 0);
    CPPUNIT_ASSERT(M.GetM() == 0);
    CPPUNIT_ASSERT(M.GetN() == 0);
  }
};

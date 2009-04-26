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
  CPPUNIT_TEST(test_constructor);
  CPPUNIT_TEST(test_reallocate);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;

public:
  void setUp()
  {
    m_ = 25;
    n_ = 10;
  }

  
  void tearDown()
  {
  }


  void test_constructor()
  {
    {
      Matrix<@real, General, @storage_rectangular_full> M@storage_rectangular_full_@real(m_, n_);
      Matrix<complex<@real>, General, @storage_rectangular_full> M@storage_rectangular_full_complex_@real(m_, n_);
      CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == m_);
      CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
      CPPUNIT_ASSERT(M@storage_rectangular_full_complex_@real.GetM() == m_);
      CPPUNIT_ASSERT(M@storage_rectangular_full_complex_@real.GetN() == n_);
    }
    {
      Matrix<@real, General, @storage_rectangular_full> M@storage_rectangular_full_@real(0, n_);
      Matrix<complex<@real>, General, @storage_rectangular_full> M@storage_rectangular_full_complex_@real(0, n_);
      CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == 0);
      CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
      CPPUNIT_ASSERT(M@storage_rectangular_full_complex_@real.GetM() == 0);
      CPPUNIT_ASSERT(M@storage_rectangular_full_complex_@real.GetN() == n_);
    }
    {
      Matrix<@real, General, @storage_rectangular_full> M@storage_rectangular_full_@real(m_, 0);
      Matrix<complex<@real>, General, @storage_rectangular_full> M@storage_rectangular_full_complex_@real(m_, 0);
      CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == m_);
      CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == 0);
      CPPUNIT_ASSERT(M@storage_rectangular_full_complex_@real.GetM() == m_);
      CPPUNIT_ASSERT(M@storage_rectangular_full_complex_@real.GetN() == 0);
    }
    {
      Matrix<@real, General, @storage_full> M@storage_full_@real(0, 0);
      Matrix<complex<@real>, General, @storage_full> M@storage_full_complex_@real(0, 0);
      CPPUNIT_ASSERT(M@storage_full_@real.GetM() == 0);
      CPPUNIT_ASSERT(M@storage_full_@real.GetN() == 0);
      CPPUNIT_ASSERT(M@storage_full_complex_@real.GetM() == 0);
      CPPUNIT_ASSERT(M@storage_full_complex_@real.GetN() == 0);
    }
    {
      Matrix<@real, General, @storage_full> M@storage_full_@real(m_, m_);
      Matrix<complex<@real>, General, @storage_full> M@storage_full_complex_@real(m_, m_);
      CPPUNIT_ASSERT(M@storage_full_@real.GetM() == m_);
      CPPUNIT_ASSERT(M@storage_full_@real.GetN() == m_);
      CPPUNIT_ASSERT(M@storage_full_complex_@real.GetM() == m_);
      CPPUNIT_ASSERT(M@storage_full_complex_@real.GetN() == m_);
    }
  }


  void test_reallocate()
  {
    Matrix<@real, General, @storage_rectangular_full> M@storage_rectangular_full_@real(m_, n_);
    Matrix<complex<@real>, General, @storage_rectangular_full> M@storage_rectangular_full_complex_@real(m_, n_);

    M@storage_rectangular_full_@real.Reallocate(0, 0);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == 0);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == 0);
    M@storage_rectangular_full_@real.Reallocate(m_, n_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == m_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
    M@storage_rectangular_full_@real.Reallocate(m_, n_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == m_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
    M@storage_rectangular_full_@real.Reallocate(2 * m_, n_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == 2 * m_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
    M@storage_rectangular_full_@real.Reallocate(2 * m_, 2 * n_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == 2 * m_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == 2 * n_);
    M@storage_rectangular_full_@real.Reallocate(m_, 0);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == m_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == 0);
    M@storage_rectangular_full_@real.Reallocate(0, n_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == 0);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
    M@storage_rectangular_full_@real.Reallocate(0, 0);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == 0);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == 0);
    M@storage_rectangular_full_@real.Reallocate(0, n_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == 0);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
    M@storage_rectangular_full_@real.Reallocate(m_, n_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == m_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
    M@storage_rectangular_full_@real.Reallocate(0, n_);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetM() == 0);
    CPPUNIT_ASSERT(M@storage_rectangular_full_@real.GetN() == n_);
  }
};

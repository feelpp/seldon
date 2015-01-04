// Copyright (C) 2001-2009 Vivien Mallet
// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_ARRAY3D_CXX

#include "Array3D.hxx"

namespace Seldon
{

  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the 3D array stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::Zero()
  {
    array3D_allocator_.memoryset(data_, char(0),
				 GetDataSize()*sizeof(value_type));
  }


  //! Fills the array.
  /*!
    On exit, the 3D array is filled with 1, 2, 3, 4, ... The order of
    those numbers depends on the storage.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::Fill()
  {
    for (int i = 0; i < GetDataSize(); i++)
      SetComplexReal(i, data_[i]);
  }


  //! Fills the 3D array with a given value.
  /*!
    On exit, the 3D array is filled with 'x'.
    \param x the value to fill the 3D array with.
  */
  template <class T, class Allocator>
  template <class T0>
  void Array3D<T, Allocator>::Fill(const T0& x)
  {
    T x_;
    SetComplexReal(x, x_);
    for (int i = 0; i < GetDataSize(); i++)
      data_[i] = x_;
  }


  //! Fills the 3D array randomly.
  /*!
    On exit, the 3D array is filled with random values.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::FillRand()
  {
#ifndef SELDON_WITHOUT_REINIT_RANDOM
    srand(time(NULL));
#endif
    for (int i = 0; i < GetDataSize(); i++)
      SetComplexReal(rand(), this->data_[i]);
  }


  //! Displays the array on the standard output.
  /*!
    Displays elements on the standard output, in text format.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::Print() const
  {
    int i, j, k;

    for (i = 0; i < GetLength1(); i++)
      {
	for (j = 0; j < GetLength2(); j++)
	  {
	    for (k = 0; k < GetLength3(); k++)
	      cout << (*this)(i, j, k) << '\t';
	    cout << endl;
	  }
	cout << endl;
      }
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the 3D array in a file.
  /*!
    Stores the 3D array in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory
    \param FileName output file name.
  */
  template <class T, class Allocator> void Array3D<T, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Array3D::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the 3D array to an output stream.
  /*!
    Writes the 3D array to an output stream in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and array elements are then written in the same order
    as in memory
    \param FileStream output stream.
  */
  template <class T, class Allocator> void Array3D<T, Allocator>
  ::Write(ofstream& FileStream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Array3D::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    if (with_size)
      {
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&length1_)),
			 sizeof(int));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&length2_)),
			 sizeof(int));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&length3_)),
			 sizeof(int));
      }

    FileStream.write(reinterpret_cast<char*>(data_),
		     length23_ * length1_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Array3D::Write(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string("  The output file may have been removed")
		    + " or there is no space left on device.");
#endif

  }


  //! Reads the 3D array from a file.
  /*!
    Reads a 3D array stored in binary format in a file.
    The dimensions of the array are read (i,j, k three integers),
    and array elements are then read in the same order
    as it should be in memory
    \param FileName input file name.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Array3D::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    Read(FileStream, with_size);

    FileStream.close();
  }


  //! Reads the 3D array from an input stream.
  /*!
    Reads a 3D array in binary format from an input stream.
    The dimensions of the array are read (i,j, k three integers),
    and array elements are then read in the same order
    as it should be in memory
    \param FileStream input stream.
  */
  template <class T, class Allocator>
  void Array3D<T, Allocator>
  ::Read(ifstream& FileStream, bool with_size)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Read(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    if (with_size)
      {
	int new_l1, new_l2, new_l3;
	FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	Reallocate(new_l1, new_l2, new_l3);
      }

    FileStream.read(reinterpret_cast<char*>(data_),
		    length23_ * length1_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Array3D::Read(ifstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The intput file may have been removed")
		    + " or may not contain enough data.");
#endif

  }


  //! operator<< overloaded for a 3D array.
  /*!
    \param out output stream.
    \param A the 3D array.
    \return The updated stream.
  */
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
			const Array3D<T, Allocator>& A)
  {
    int i, j, k;

    for (i = 0; i < A.GetLength1(); i++)
      {
	for (j = 0; j < A.GetLength2(); j++)
	  {
	    for (k = 0; k < A.GetLength3(); k++)
	      out << A(i, j, k) << '\t';
	    out << endl;
	  }
	out << endl;
      }

    return out;
  }


  //! Multiplication of all elements of a 3D array by a scalar.
  /*!
    \param[in] alpha scalar by which \a A should be multiplied.
    \param[in,out] A the array to be multiplied by \a alpha.
  */
  template <class T0, class T, class Allocator>
  void Mlt(const T0& alpha, Array3D<T, Allocator>& A)
  {
    T* data = A.GetData();
    for (int i = 0; i < A.GetDataSize(); i++)
      data[i] *= alpha;
  }
  
} // namespace Seldon.

#define SELDON_FILE_ARRAY3D_CXX
#endif

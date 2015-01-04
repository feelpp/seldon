// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu
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


#ifndef SELDON_FILE_MATRIX_COLLECTION_CXX

#include "MatrixCollection.hxx"


namespace Seldon
{


  //////////////////////
  // MATRIXCOLLECTION //
  //////////////////////
  
  
  /************************
   * CONVENIENT FUNCTIONS *
   ************************/
  
  
  //! Displays the matrix collection on the standard output.
  /*!
    Displays elements on the standard output, in text format.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixCollection<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < Mmatrix_; i++)
      {
	for (int j = 0; j < Nmatrix_; j++)
	  cout << GetMatrix(i, j) << endl;
	cout << endl;
      }
  }


  //! Displays an underlying matrix on the standard output.
  /*!
    \param[in] m row index of the underlying matrix.
    \param[in] n column index of the underlying matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixCollection<T, Prop, Storage, Allocator>
  ::Print(int i, int j) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Mmatrix_)
      throw WrongRow("MatrixCollection::Print(int, int)",
                     string("Row index should be in [0, ")
                     + to_str(Mmatrix_ - 1) + "], but is equal to "
                     + to_str(i) + ".");
    if (j < 0 || j >= Nmatrix_)
      throw WrongCol("MatrixCollection::Print(int, int)",
                     string("Column index should be in [0, ")
                     + to_str(Nmatrix_ - 1) + "], but is equal to "
                     + to_str(j) + ".");
#endif

    cout << matrix_(i, j) << endl;
  }


  //! Writes the matrix collection in a file.
  /*! Stores the matrix collection in a file in binary format. The number of
    rows of matrices (integer) and the number of columns of matrices (integer)
    are written, and the underlying matrices are then written in the same
    order as in memory (e.g. row-major storage).
    \param[in] FileName output file name.
    \param[in] with_size if set to 'false', the dimensions of the matrix are
    not saved.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixCollection<T, Prop, Storage, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("MatrixCollection::Write(string FileName, bool)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the matrix collection to an output stream.
  /*! Writes the matrix collection to an output stream in binary format.  The
    number of rows of matrices (integer) and the number of columns of matrices
    (integer) are written, and the underlying matrices are then written in the
    same order as in memory (e.g. row-major storage).
    \param[in,out] FileStream output stream.
    \param[in] with_size if set to 'false', the dimensions of the matrix are
    not saved.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixCollection<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream, bool with_size = true) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("MatrixCollection::Write(ostream& FileStream, bool)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&Mmatrix_)),
                         sizeof(int));
        FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&Nmatrix_)),
                         sizeof(int));
      }

    int i, j;
    for (i = 0; i < this->GetMmatrix(); i++)
      {
        for (j = 0; j < this->GetNmatrix(); j++)
          GetMatrix(i, j).Write(FileStream, with_size);
      }

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("MatrixCollection::Write(ostream& FileStream, bool)",
                    "Output operation failed.");
#endif

  }


  //! Writes the matrix collection in a file.
  /*! Stores the matrix in a file in text format. Only the underlying matrices
    are written, without the dimensions. Each row is written on a single line
    and elements of a row are delimited by tabulations.
    \param[in] FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixCollection<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("MatrixCollection::WriteText(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the matrix collection to an output stream.
  /*! Stores the matrix to an output stream in text format. Only the
    underlying matrices are written, without the dimensions. Each row is
    written on a single line and elements of a row are delimited by
    tabulations.
    \param[in,out] FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixCollection<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("MatrixCollection::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    int i, j;
    for (i = 0; i < this->GetMmatrix(); i++)
      {
        for (j = 0; j < this->GetNmatrix(); j++)
          GetMatrix(i, j).WriteText(FileStream);
        FileStream << endl;
      }

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("MatrixCollection::WriteText(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Reads the matrix collection from a file.
  /*! Reads a matrix collection stored in binary format in a file.  The number
    of rows of matrices (integer) and the number of columns of matrices
    (integer) are read, and the underlying matrices are then read in the same
    order as it should be in memory (e.g. row-major storage).
    \param[in] FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixCollection<T, Prop, Storage, Allocator>::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("MatrixCollection::Read(string FileName)",
                    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix collection from an input stream.
  /*! Reads a matrix collection stored in binary format from a stream.  The
    number of rows of matrices (integer) and the number of columns of matrices
    (integer) are read, and the underlying matrices are then read in the same
    order as it should be in memory (e.g. row-major storage).
    \param[in,out] FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void MatrixCollection<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("MatrixCollection::Read(istream& FileStream)",
                    "The stream is not ready.");
#endif

    int *new_m, *new_n;
    new_m = new int;
    new_n = new int;

    FileStream.read(reinterpret_cast<char*>(new_m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(new_n), sizeof(int));

    this->Reallocate(*new_m, *new_n);

    T working_matrix;
    int i, j;
    for (i = 0; i < *new_m; i++)
      for (j = 0; j < *new_n; j++)
	{
	  working_matrix.Read(FileStream);
	  SetMatrix(i, j, working_matrix);
	  working_matrix.Nullify();
	}


    delete new_n;
    delete new_m;

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("MatrixCollection::Read(istream& FileStream)",
                    "Input operation failed.");
#endif

  }
  
  
} // namespace Seldon.


#define SELDON_FILE_MATRIX_COLLECTION_CXX
#endif

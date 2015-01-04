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


#ifndef SELDON_FILE_HETEROGENEOUS_MATRIX_COLLECTION_CXX

#include "HeterogeneousMatrixCollection.hxx"

namespace Seldon
{


  ///////////////////////////////////
  // HETEROGENEOUSMATRIXCOLLECTION //
  ///////////////////////////////////


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Displays the matrix collection on the standard output.
  /*!
    Displays elements on the standard output, in text format.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Print() const
  {
    for (int i = 0; i < Mlocal_sum_(Mmatrix_); i++)
      {
	for (int j = 0; j < Nlocal_sum_(Nmatrix_); j++)
	  cout << (*this)(i, j) << endl;
	cout << endl;
      }
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
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::Write(string FileName)",
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
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Write(ostream& FileStream, bool with_size = true) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::Write(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&Mmatrix_)),
                         sizeof(int));
        FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&Nmatrix_)),
                         sizeof(int));
      }

    collection_.Write(FileStream, with_size);

    float_dense_m m0a;
    float_sparse_m m1a;
    double_dense_m m2a;
    double_sparse_m m3a;

    int i, j;
    for (i = 0; i < Mmatrix_; i++)
      for (j = 0; j < Nmatrix_; j++)
	{
	  switch (GetType(i, j))
	    {
	    case 0:
	      GetMatrix(i, j, m0a);
	      m0a.Write(FileStream, with_size);
	      m0a.Nullify();
	      break;
	    case 1:
	      throw Undefined("Matrix<FloatDouble, DenseSparseCollection>"
			      "Storage0, Prop1, Storage1, Allocator>"
			      "::Write(ostream& FileStream, bool "
			      "with_size = true) ");
	    case 2:
	      GetMatrix(i, j, m2a);
	      m2a.Write(FileStream, with_size);
	      m2a.Nullify();
	      break;
	    case 3:
	      throw Undefined("Matrix<FloatDouble, DenseSparseCollection>"
			      "Storage0, Prop1, Storage1, Allocator>"
			      "::Write(ostream& FileStream, bool "
			      "with_size = true) ");
	    default:
	      throw WrongArgument("Matrix<FloatDouble, DenseSparseCollection>"
				  "::Write(ostream& FileStream, "
                                  "bool with_size = true) ",
				  "Underlying matrix  A (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	}


#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::Write(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Writes the matrix collection in a file.
  /*! Stores the matrix in a file in text format. Only the underlying matrices
    are written, without the dimensions. Each row is written on a single line
    and elements of a row are delimited by tabulations.
    \param[in] FileName output file name.
  */
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("HeterogeneousMatrixCollection"
                    "::WriteText(string FileName)",
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
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    float_dense_m m0a;
    float_sparse_m m1a;
    double_dense_m m2a;
    double_sparse_m m3a;

    int i, j;
    for (i = 0; i < Mmatrix_; i++)
      for (j = 0; j < Nmatrix_; j++)
	{
	  switch (GetType(i, j))
	    {
	    case 0:
	      GetMatrix(i, j, m0a);
	      m0a.WriteText(FileStream);
	      m0a.Nullify();
	      break;
	    case 1:
	      GetMatrix(i, j, m1a);
	      m1a.WriteText(FileStream);
	      m1a.Nullify();
	      break;
	    case 2:
	      GetMatrix(i, j, m2a);
	      m2a.WriteText(FileStream);
	      m2a.Nullify();
	      break;
	    case 3:
	      GetMatrix(i, j, m3a);
	      m3a.WriteText(FileStream);
	      m3a.Nullify();
	      break;
	    default:
	      throw WrongArgument("Matrix<FloatDouble, DenseSparseCollection>"
				  "::Write(ostream& FileStream, "
                                  "bool with_size = true) ",
				  "Underlying matrix  A (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	  FileStream << endl;
	}

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::WriteText(ostream& FileStream)",
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
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("HeterogeneousMatrixCollection<Prop0, Storage0, Prop1,"
		    " Storage1, Allocator>::Read(string FileName)",
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
  template <class Prop0, class Storage0,
	    class Prop1, class Storage1,
	    template <class U> class Allocator> void
  HeterogeneousMatrixCollection<Prop0, Storage0, Prop1, Storage1, Allocator>
  ::Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection<Prop0, Storage0, Prop1,"
		    " Storage1, Allocator>::Read(istream& FileStream)",
                    "The stream is not ready.");
#endif

    int *new_m, *new_n;
    new_m = new int;
    new_n = new int;

    FileStream.read(reinterpret_cast<char*>(new_m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(new_n), sizeof(int));

    this->Reallocate(*new_m, *new_n);

    collection_.Read(FileStream);

    float_dense_m m0a;
    float_sparse_m m1a;
    double_dense_m m2a;
    double_sparse_m m3a;
    int i, j;
    for (i = 0; i < Mmatrix_; i++)
      for (j = 0; j < Nmatrix_; j++)
	{
          switch (GetType(i, j))
	    {
	    case 0:
	      m0a.Read(FileStream);
	      SetMatrix(i, j, m0a);
	      m0a.Nullify();
	      break;
	    case 1:
	      throw Undefined("Matrix<FloatDouble, DenseSparseCollection>"
			      "Storage0, Prop1, Storage1, Allocator>"
			      "::Read(istream& FileStream)");
	    case 2:
	      m2a.Read(FileStream);
	      SetMatrix(i, j, m2a);
	      m2a.Nullify();
	      break;
	    case 3:
	      throw Undefined("Matrix<FloatDouble, DenseSparseCollection>"
			      "Storage0, Prop1, Storage1, Allocator>"
			      "::Read(istream& FileStream)");
	      break;
	    default:
	      throw WrongArgument("HeterogeneousMatrixCollection<Prop0, "
				  "Storage0, Prop1, Storage1, Allocator>"
				  "::Read(istream& FileStream) ",
				  "Underlying matrix  A (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	}


    delete new_n;
    delete new_m;

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("HeterogeneousMatrixCollection"
                    "::Read(istream& FileStream)",
                    "Input operation failed.");
#endif

  }
  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_HETEROGENEOUS_COLLECTION_CXX
#endif

// Copyright (C) 2010 Lin Wu
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


#ifndef SELDON_FILE_ARRAY_CXX


namespace Seldon
{

  //! Duplicates an array.
  /*!
    \param A array to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Copy(const Array<T, N, Allocator>& A)
  {
    if (N == 3)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2));
    if (N == 4)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3));
    if (N == 5)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4));
    if (N == 6)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5));
    if (N == 7)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6));
    if (N == 8)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6), A.GetLength(7));
    if (N == 9)
      Reallocate(A.GetLength(0), A.GetLength(1), A.GetLength(2),
		 A.GetLength(3), A.GetLength(4), A.GetLength(5),
		 A.GetLength(6), A.GetLength(7), A.GetLength(8));

    array_allocator_.memorycpy(data_, A.GetData(), GetDataSize());
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the array stores complex
    structures, use 'Fill' instead.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Zero()
  {
    array_allocator_.memoryset(data_, char(0),
			       GetDataSize()*sizeof(value_type));
  }


  //! Fills the array.
  /*!
    On exit, the array is filled with 1, 2, 3, 4, ... The order of
    those numbers depends on the storage.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Fill()
  {
    for (int i = 0; i < GetDataSize(); i++)
      SetComplexReal(i, data_[i]);
  }


  //! Fills the array with a given value.
  /*!
    On exit, the array is filled with 'x'.
    \param x the value to fill the array with.
  */
  template <class T, int N, class Allocator>
  template <class T0>
  void Array<T, N, Allocator>::Fill(const T0& x)
  {
    T x_;
    SetComplexReal(x, x_);
    for (int i = 0; i < GetDataSize(); i++)
      data_[i] = x_;
  }


  //! Fills the array randomly.
  /*!
    On exit, the array is filled with random values.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::FillRand()
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
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Print() const
  {
    if (N == 3)
      {
	int i, j, k;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  cout << (*this)(i, j, k) << '\t';
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 4)
      {
	int i, j, k, l;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      cout << (*this)(i, j, k, l) << '\t';
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 5)
      {
	int i, j, k, l, m;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  cout << (*this)(i, j, k, l, m) << '\t';
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 6)
      {
	int i, j, k, l, m, n;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      cout << (*this)(i, j, k, l, m, n) << '\t';
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 7)
      {
	int i, j, k, l, m, n, o;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  cout << (*this)(i, j, k, l, m, n, o)
				       << '\t';
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 8)
      {
	int i, j, k, l, m, n, o, p;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  {
				    for (p = 0; p < GetLength(7); p++)
				      cout << (*this)(i, j, k, l, m, n, o, p)
					   << '\t';
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 9)
      {
	int i, j, k, l, m, n, o, p, q;
	for (i = 0; i < GetLength(0); i++)
	  {
	    for (j = 0; j < GetLength(1); j++)
	      {
		for (k = 0; k < GetLength(2); k++)
		  {
		    for (l = 0; l < GetLength(3); l++)
		      {
			for (m = 0; m < GetLength(4); m++)
			  {
			    for (n = 0; n < GetLength(5); n++)
			      {
				for (o = 0; o < GetLength(6); o++)
				  {
				    for (p = 0; p < GetLength(7); p++)
				      {
					for (q = 0; q < GetLength(8); q++)
					  cout << (*this)(i, j, k, l, m, n,
							  o, p, q)
					       << '\t';
					cout << endl;
				      }
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the array in a file.
  /*!
    Stores the array in a file in binary format. The dimensions (integer)
    are written, and array elements are then written in the same order as in
    memory.
    \param FileName output file name.
  */
  template <class T, int N, class Allocator> void Array<T, N, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Array::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the array to an output stream.
  /*!
    Writes the array to an output stream in binary format. The dimensions
    (integerS) are written, and array elements are then written in the same
    order as in memory.
    \param FileStream output stream.
  */
  template <class T, int N, class Allocator> void Array<T, N, Allocator>
  ::Write(ofstream& FileStream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Array::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    if (with_size)
      {
	for (int i = 0; i < N; i++)
	  FileStream.write(reinterpret_cast<char*>
			   (const_cast<int*>(&length_[i])),
			   sizeof(int));
      }

    FileStream.write(reinterpret_cast<char*>(data_),
		     offset_[N - 1] * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Array::Write(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string("  The output file may have been removed")
		    + " or there is no space left on device.");
#endif

  }


  //! Reads the array from a file.
  /*!
    Reads a array stored in binary format in a file.  The dimensions
    (integers) of the array are read, and array elements are then read in the
    same order as it should be in memory
    \param FileName input file name.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Array::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    Read(FileStream, with_size);

    FileStream.close();
  }


  //! Reads the array from an input stream.
  /*!
    Reads a array in binary format from an input stream.  The dimensions
    (integers) of the array are read, and array elements are then read in the
    same order as it should be in memory.
    \param FileStream input stream.
  */
  template <class T, int N, class Allocator>
  void Array<T, N, Allocator>
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
	if (N == 3)
	  {
	    int new_l1, new_l2, new_l3;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3);
	  }
	if (N == 4)
	  {
	    int new_l1, new_l2, new_l3, new_l4;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4);
	  }
	if (N == 5)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5);
	  }
	if (N == 6)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6);
	  }
	if (N == 7)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7);
	  }
	if (N == 8)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7,
	      new_l8;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l8), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7, new_l8);
	  }
	if (N == 9)
	  {
	    int new_l1, new_l2, new_l3, new_l4, new_l5, new_l6, new_l7,
	      new_l8, new_l9;
	    FileStream.read(reinterpret_cast<char*>(&new_l1), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l2), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l3), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l4), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l5), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l6), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l7), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l8), sizeof(int));
	    FileStream.read(reinterpret_cast<char*>(&new_l9), sizeof(int));
	    Reallocate(new_l1, new_l2, new_l3, new_l4, new_l5, new_l6,
		       new_l7, new_l8, new_l9);
	  }
      }

    FileStream.read(reinterpret_cast<char*>(data_), 
		    offset_[N - 1] * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Array::Read(ifstream& FileStream)",
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
  template <class T, int N, class Allocator>
  ostream& operator << (ostream& out,
			const Array<T, N, Allocator>& A)
  {
    if (N == 3)
      {
	int i, j, k;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  out << A(i, j, k) << '\t';
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 4)
      {
	int i, j, k, l;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      out << A(i, j, k, l) << '\t';
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 5)
      {
	int i, j, k, l, m;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  out << A(i, j, k, l, m) << '\t';
			out << endl;
		      }
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 6)
      {
	int i, j, k, l, m, n;

	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      out << A(i, j, k, l, m, n) << '\t';
			    out << endl;
			  }
			out << endl;
		      }
		    out << endl;
		  }
		out << endl;
	      }
	    out << endl;
	  }
      }

    if (N == 7)
      {
	int i, j, k, l, m, n, o;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  cout << A(i, j, k, l, m, n, o)
				       << '\t';
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 8)
      {
	int i, j, k, l, m, n, o, p;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  {
				    for (p = 0; p < A.GetLength(7); p++)
				      cout << A(i, j, k, l, m, n, o, p)
					   << '\t';
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    if (N == 9)
      {
	int i, j, k, l, m, n, o, p, q;
	for (i = 0; i < A.GetLength(0); i++)
	  {
	    for (j = 0; j < A.GetLength(1); j++)
	      {
		for (k = 0; k < A.GetLength(2); k++)
		  {
		    for (l = 0; l < A.GetLength(3); l++)
		      {
			for (m = 0; m < A.GetLength(4); m++)
			  {
			    for (n = 0; n < A.GetLength(5); n++)
			      {
				for (o = 0; o < A.GetLength(6); o++)
				  {
				    for (p = 0; p < A.GetLength(7); p++)
				      {
					for (q = 0; q < A.GetLength(8); q++)
					  cout << A(i, j, k, l, m, n, o, p, q)
					       << '\t';
					cout << endl;
				      }
				    cout << endl;
				  }
				cout << endl;
			      }
			    cout << endl;
			  }
			cout << endl;
		      }
		    cout << endl;
		  }
		cout << endl;
	      }
	    cout << endl;
	  }
      }

    return out;
  }


} // namespace Seldon.

#define SELDON_FILE_ARRAY_CXX
#endif

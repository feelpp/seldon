// Copyright (C) 2001-2012 Vivien Mallet
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


#ifndef SELDON_FILE_COMMON_CXX

#include "Common.hxx"

template <class T>
inline void PrintArray(T* v, int lgth)
{
  for (int k = 0; k < lgth - 1; k++)
    std::cout << v[k] << " | ";
  std::cout << v[lgth - 1] << std::endl;
}

namespace Seldon
{


  //! Default constructor.
  inline Str::Str()
  {
  }


  //! Copy constructor.
  /*!
    \param[in] s 'Str' instance to be copied.
  */
  inline Str::Str(const Str& s)
  {
    output_ << s;
  }


  //! Conversion to string.
  inline Str::operator std::string() const
  {
    return output_.str();
  }


  //! Adds an element to the string.
  /*!
    \param[in] input element added at the end of the string.
  */
  template <class T>
  inline Str& Str::operator << (const T& input)
  {
    output_ << input;
    return *this;
  }


  //! Adds an element to an instance of 'Str'.
  /*!
    \param[in] s 'Str' instance.
    \param[in] input element added at the end of the string.
  */
  template <class T>
  inline Str operator + (const Str& s, const T& input)
  {
    string s_input = s;
    Str output;
    output << s_input << input;
    return output;
  }


  //! Converts a 'str' instance to an 'ostream' instance.
  inline ostream& operator << (ostream& out, Str& in)
  {
    string output = in;
    out << output;
    return out;
  }


  //! Converts a 'str' instance to an 'ostream' instance.
  inline ostream& operator << (ostream& out, Str in)
  {
    string output = in;
    out << output;
    return out;
  }


  //! Converts most types to string.
  /*!
    \param input variable to be converted.
    \return A string containing 'input'.
  */
  template<typename T>
  inline std::string to_str(const T& input)
  {
    std::ostringstream output;
    output << input;
    return output.str();
  }


  //! Converts string to most types, specially numbers.
  /*!
    \param[in] s string to be converted.
    \param[out] num \a s converted to 'T'.
  */
  template <class T>
  inline void to_num(std::string s, T& num)
  {
    std::istringstream str(s);
    str >> num;
  }


  //! Converts string to most types, specially numbers.
  /*!
    \param[in] s string to be converted.
    \return \a s converted to 'T'.
  */
  template <class T>
  inline T to_num(std::string s)
  {
    T num;
    std::istringstream str(s);
    str >> num;
    return num;
  }


  //! Sets a number to zero.
  /*!
    \param[in,out] number number to be set to zero.
  */
  template <class T>
  inline void SetComplexZero(T& number)
  {
    number = T(0);
  }


  //! Sets a complex number to zero.
  /*!
    \param[in,out] number complex number to be set to zero.
  */
  template <class T>
  inline void SetComplexZero(complex<T>& number)
  {
    number = complex<T>(T(0), T(0));
  }


  //! Sets a number to one.
  /*!
    \param[in,out] number number to be set to one.
  */
  template <class T>
  inline void SetComplexOne(T& number)
  {
    number = T(1);
  }


  //! Sets a complex number to (1, 0).
  /*!
    \param[in,out] number complex number to be set to (1, 0).
  */
  template <class T>
  inline void SetComplexOne(complex<T>& number)
  {
    number = complex<T>(T(1), T(0));
  }


  //! Sets a real number to n.
  /*!
    \param[in,out] number real umber to be set to n.
  */
  template <class T>
  inline void SetComplexReal(int n, T& number)
  {
    number = T(n);
  }
  
    
  //! Sets a complex number to (n, 0).
  /*!
    \param[in,out] number complex number to be set to (n, 0).
  */
  template <class T>
  inline void SetComplexReal(int n, complex<T>& number)
  {
    number = complex<T>(n, 0);
  }
  

  //! Sets a complex number to (n, 0).
  /*!
    \param[in,out] number complex number to be set to (n, 0).
  */
  template <class T>
  inline void SetComplexReal(bool n, complex<T>& number)
  {
    number = complex<T>(n, 0);
  }


  //! Sets a complex number to (x, 0).
  /*!
    \param[in,out] number complex number to be set to (x, 0).
  */  
  template <class T>
  inline void SetComplexReal(const T& x, complex<T>& number)
  {
    number = complex<T>(x, 0);
  }
  
  
  //! Sets a complex number to x.
  /*!
    \param[in,out] number complex number to be set to x.
  */
  template <class T0, class T1>
  inline void SetComplexReal(const T0& x, T1& number)
  {
    number = x;
  }
  
  
  //! Returns true for a complex number
  template <class T>
  inline bool IsComplexNumber(const T& number)
  {
    return false;
  }


  //! Returns true for a complex number
  template <class T>
  inline bool IsComplexNumber(const complex<T>& number)
  {
    return true;
  }

  
  //! returns absolute value of val
  template<class T>
  inline T ComplexAbs(const T& val)
  {
    return abs(val);
  }
  
  
  //! returns modulus of val
  template<class T>
  inline T ComplexAbs(const complex<T>& val)
  {
#if defined(SELDON_WITH_BLAS) && !defined(SELDON_WITH_LAPACK)
    // we choose Blas convention
    return abs(real(val)) + abs(imag(val));
#else
    // otherwise usual modulus
    return abs(val);
#endif
  }


  //! returns the square modulus of z
  template<class T>
  inline T absSquare(const complex<T>& z)
  {
    // more optimal than real(z * conj(z))
    return real(z)*real(z) + imag(z)*imag(z);
  }

  
  //! returns the square modulus of z
  template<class T>
  inline T absSquare(const T& z)
  {
    return z*z;
  }
  
  
  //! returns extension of a string 
  inline string GetExtension(const string& nom)
  {
    size_t index = nom.find_last_of('.');
    if (index == string::npos)
      return string("");
    
    string extension = nom.substr(index+1, nom.size()-index);  
    return extension;
  }

  
  //! returns base of a string
  inline string GetBaseString(const string& nom)
  {
    size_t index = nom.find_last_of('.');
    if (index == string::npos)
      return nom;
    
    string base = nom.substr(0, index);
    return base;
  }
  
    
#ifdef SELDON_WITH_HDF5
  //! Gives for most C types the corresponding HDF5 memory datatype.
  /*!
    \param[in] input variable to analyze.
    \return HDF5 memory type of \a input.
  */
  template <class T>
  inline hid_t GetH5Type(T& input)
  {
    double d;
    float f;
    int i;
    long l;
    char c;
    unsigned char uc;
    long long ll;
    unsigned int ui;
    unsigned short us;
    unsigned long ul;
    unsigned long long ull;

    if (typeid(input) == typeid(d))
      return H5T_NATIVE_DOUBLE;
    if (typeid(input) == typeid(f))
      return H5T_NATIVE_FLOAT;
    if (typeid(input) == typeid(i))
      return H5T_NATIVE_INT;
    if (typeid(input) == typeid(l))
      return H5T_NATIVE_LONG;
    if (typeid(input) == typeid(c))
      return H5T_NATIVE_CHAR;
    if (typeid(input) == typeid(uc))
      return H5T_NATIVE_UCHAR;
    if (typeid(input) == typeid(ll))
      return H5T_NATIVE_LLONG;
    if (typeid(input) == typeid(ui))
      return H5T_NATIVE_UINT;
    if (typeid(input) == typeid(us))
      return H5T_NATIVE_USHORT;
    if (typeid(input) == typeid(ul))
      return H5T_NATIVE_ULONG;
    if (typeid(input) == typeid(ull))
      return H5T_NATIVE_ULLONG;
    else
      throw Error("hid_t GetH5Type(T& input)",
                  "Type has no corresponding native HDF5 datatype.");
  }
#endif


}  // namespace Seldon.

#define SELDON_FILE_COMMON_CXX
#endif

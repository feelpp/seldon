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

#ifndef SELDON_FILE_ALLOCATOR_HXX

namespace Seldon
{


  /////////////////
  // MALLOCALLOC //
  /////////////////


  template <class T>
  class MallocAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  public:

    pointer allocate(int num, void* h = 0);
    void deallocate(pointer data, int num, void* h = 0);
    void* reallocate(pointer data, int num, void* h = 0);
    void memoryset(pointer data, char c, size_t num);
    void memorycpy(pointer datat, pointer datas, size_t num);
  };


  /////////////////
  // CALLOCALLOC //
  /////////////////


  template <class T>
  class CallocAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  public:

    pointer allocate(int num, void* h = 0);
    void deallocate(pointer data, int num, void* h = 0);
    void* reallocate(pointer data, int num, void* h = 0);
    void memoryset(pointer data, char c, size_t num);
    void memorycpy(pointer datat, pointer datas, size_t num);
  };


  //////////////
  // NEWALLOC //
  //////////////


  template <class T>
  class NewAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  public:

    pointer allocate(int num, void* h = 0);
    void deallocate(pointer data, int num, void* h = 0);
    void* reallocate(pointer data, int num, void* h = 0);
    void memoryset(pointer data, char c, size_t num);
    void memorycpy(pointer datat, pointer datas, size_t num);
  };


  //////////////
  // NANALLOC //
  //////////////


  template <class T>
  class NaNAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  public:

    pointer allocate(int num, void* h = 0);
    void deallocate(pointer data, int num, void* h = 0);
    void* reallocate(pointer data, int num, void* h = 0);
    void memoryset(pointer data, char c, size_t num);
    void memorycpy(pointer datat, pointer datas, size_t num);
  };


} // namespace Seldon.

#define SELDON_FILE_ALLOCATOR_HXX
#endif

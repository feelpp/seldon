%module seldon
%{
#include "SeldonHeader.hxx"
%}

%include "std_string.i"

using namespace std;

// Include the header file with above prototypes
%include "SeldonHeader.hxx"
%include "Common/Common.hxx"
%include "Common/Storage.hxx"
%include "Common/Properties.hxx"
%include "Vector/Vector.hxx"
%include "Matrix/Matrix_Base.hxx"
%include "Matrix/Matrix_Pointers.hxx"
%include "Common/Allocator.hxx"

namespace Seldon
{
%extend Vector<double, Vect_Full, MallocAlloc<double> >
{
    double __getitem__(int index) {
        if (index < self->GetM())
          return self->GetData()[index];
        else
          return 0;
    }
    void __setitem__(int index, double value) {
        if (index >= 0 && index < self->GetM()) {
            self->GetData()[index] = value;
        }
    }
    unsigned long __len__() {
          return self->GetM();
    }
}
%extend Matrix<double, General, RowMajor, MallocAlloc<double> >
{
    double __getitem__(PyObject* args)
    {
	int i, j;
	int success = PyArg_ParseTuple(args, "ii", &i, &j);
	if (!success)
	   throw std::out_of_range("Failed!");
	return (*self)(i, j);
    }
    double __getitem__(int i)
    {
    	return (*self)[i];
    }
    void __setitem__(PyObject* args, double value)
    {
	int i, j;
	int success = PyArg_ParseTuple(args, "ii", &i, &j);
	if (!success)
	   throw std::out_of_range("Failed!");
	(*self)(i, j) = value;
    }
    unsigned long __len__()
    {
	return self->GetSize();
    }
}
%template(DoubleMalloc) MallocAlloc<double>;
%template(BaseSeldonVector) Vector_Base<double, MallocAlloc<double> >;
%template(VectorDouble) Vector<double, Vect_Full, MallocAlloc<double> >;
%template(MatrixBaseDouble) Matrix_Base<double, MallocAlloc<double> >;
%template(MatrixPointersDouble) Matrix_Pointers<double, General, RowMajor, MallocAlloc<double> >;
%template(MatrixDouble) Matrix<double, General, RowMajor, MallocAlloc<double> >;
}

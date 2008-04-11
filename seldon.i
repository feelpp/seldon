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
%include "Vector/Vector.hxx"
%include "Common/Allocator.hxx"

namespace Seldon
{
%template(DoubleMalloc) MallocAlloc<double>;
%template(BaseSeldonVector) Vector_Base<double, MallocAlloc<double> >;
%template(VectorDouble) Vector<double, Vect_Full, MallocAlloc<double> >;
}

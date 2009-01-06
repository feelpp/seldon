#define SELDON_DEBUG_LEVEL_4
#define SELDON_DEFAULT_ALLOCATOR NewAlloc
#define SELDON_WITH_ABORT

#include "Seldon.hxx"

using namespace Seldon;

// testing class Vector without Blas

int main()
{
  Vector<double> V(0);
  // Vector<double> U(-2);
  V.Reallocate(0);
  V.Reallocate(-2);
  V(0) = 1.0;
  V.Clear();
}

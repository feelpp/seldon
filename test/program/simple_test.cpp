#define SELDON_DEFAULT_ALLOCATOR NewAlloc
#include "Seldon.hxx"

using namespace Seldon;

int main()
{

  cout << "Seldon: compilation test" << endl;

  Vector<double> V(3);
  V.Fill();

  cout << "Vector: " << V << endl;

  V.PushBack(19);

  cout << "Vector: " << V << endl;

  return 0;

}

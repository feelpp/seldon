#ifndef SELDON_FILE_ANASAZI_HXX

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCG.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockDavidson.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
// #include "AnasaziEpetraAdapter.hpp"
#include "AnasaziMVOPTester.hpp"

#ifdef SELDON_WITH_MPI
//#include "Epetra_MpiComm.h"
#else
//#include "Epetra_SerialComm.h"
#endif

// #include "Epetra_Map.h"

// templated multivector
#include "AnasaziMultiVec.hpp"
#include "MyMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "Teuchos_BLAS.hpp"

namespace Anasazi
{
  using namespace std;
  
  template<class EigenPb, class T>
  class OperatorAnasaziEigen : public Operator<T>
  {
  public :
    //! reference to the object containing eigenproblem data    
    EigenPb& var_eig;
    //! operator stored (M, A, Op)
    int type_operator;
    enum {OPERATOR_A, OPERATOR_M, OPERATOR_OP};
    
    OperatorAnasaziEigen(EigenPb& var, int type) : Operator<T>(), var_eig(var)
    {
      type_operator = type;
    }
    
    // computes Y0 = Operator X0
    void Apply(const MultiVec<T>& X0, MultiVec<T>& Y0) const
    {
      // converting to SeldonMultiVec
      const MyMultiVec<T>& X = static_cast<const MyMultiVec<T>& >(X0);
      MyMultiVec<T>& Y = static_cast<MyMultiVec<T>& >(Y0);
      int n = X.GetVecLength();
      int nvecs = X.GetNumberVecs();
      
      Seldon::Vector<T> xvec, yvec, xvec0;
      xvec.Reallocate(n);
      switch (type_operator)
	{
	case OPERATOR_A :
	case OPERATOR_OP :
	  // loop over vectors
	  for (int p = 0; p < nvecs; p++)
	    {
	      T* y = Y[p];
	      T* x = X[p];
	      xvec0.SetData(n, x);
	      Copy(xvec0, xvec);
	      yvec.SetData(n, y);
	      
	      if (var_eig.DiagonalMass() || var_eig.UseCholeskyFactoForMass())
		{
		  // standard eigenvalue problem
		  if (var_eig.GetComputationalMode() == var_eig.REGULAR_MODE)
		    {
		      if (var_eig.DiagonalMass())
			var_eig.MltInvSqrtDiagonalMass(xvec);
		      else
			var_eig.SolveCholeskyMass(Seldon::SeldonTrans, xvec);
		      
		      var_eig.MltStiffness(xvec, yvec);
		      
		      if (var_eig.DiagonalMass())
			var_eig.MltInvSqrtDiagonalMass(yvec);
		      else
			var_eig.SolveCholeskyMass(Seldon::SeldonNoTrans, yvec);
		    }
		  else
		    {
		      if (var_eig.DiagonalMass())
			var_eig.MltSqrtDiagonalMass(xvec);
		      else
			var_eig.MltCholeskyMass(Seldon::SeldonNoTrans, xvec);
		      
		      var_eig.ComputeSolution(xvec, yvec);
		      
		      if (var_eig.DiagonalMass())
			var_eig.MltSqrtDiagonalMass(yvec);
		      else
			var_eig.MltCholeskyMass(Seldon::SeldonTrans, yvec);
		    }
		}
	      else
		{
		  if (var_eig.GetComputationalMode() == var_eig.INVERT_MODE)
		    {
		      if (var_eig.GetTypeSpectrum() != var_eig.CENTERED_EIGENVALUES)
			{
			  var_eig.MltStiffness(xvec0, xvec);
			  var_eig.ComputeSolution(xvec, yvec);
			}
		      else
			{
			  var_eig.MltMass(xvec0, xvec);
			  var_eig.ComputeSolution(xvec, yvec);
			}
		    }
		  else if (var_eig.GetComputationalMode() == var_eig.REGULAR_MODE)
		    {
		      var_eig.MltStiffness(xvec0, yvec);
		    }
		  else
		    {
		      cout << "not implemented " << endl;
		      abort();
		    }
		}
	      
	      var_eig.IncrementProdMatVect();
	      
	      xvec0.Nullify();
	      yvec.Nullify();
	    }
	  break;
	case OPERATOR_M :
	  // loop over vectors
	  for (int p = 0; p < nvecs; p++)
	    {
	      T* y = Y[p];
	      T* x = X[p];
	      xvec.SetData(n, x);
	      yvec.SetData(n, y);
	      
	      if (var_eig.DiagonalMass() || var_eig.UseCholeskyFactoForMass())
		Copy(xvec, yvec);
	      else
		{
		  if (var_eig.GetComputationalMode() == var_eig.INVERT_MODE)
		    Copy(xvec, yvec);
		  else
		    var_eig.MltMass(xvec, yvec);
		}
	      
	      xvec.Nullify();
	      yvec.Nullify();
	    }
	  break;
	}
    }
    
  };
  
}

#define SELDON_FILE_ANASAZI_HXX
#endif

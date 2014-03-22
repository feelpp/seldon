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
      
      Seldon::Vector<T> xvec, yvec;
      switch (type_operator)
	{
	case OPERATOR_A :
	case OPERATOR_OP :
	  // loop over vectors
	  for (int p = 0; p < nvecs; p++)
	    {
	      T* y = Y[p];
	      T* x = X[p];
	      xvec.SetData(n, x);
	      yvec.SetData(n, y);
	      
	      var_eig.MltStiffness(xvec, yvec);
	      
	      xvec.Nullify();
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
	      
	      var_eig.MltMass(xvec, yvec);
	      
	      xvec.Nullify();
	      yvec.Nullify();
	    }
	  break;
	}
    }
    
  };
	
  /* template<class TypeElement, class TypeEquation>
  class OperatorTraits <double, Epetra_MultiVector, Montjoie::MatrixEigen_Montjoie<TypeElement, TypeEquation> > 
  {
  public :
		
    static void Apply(const Montjoie::MatrixEigen_Montjoie<TypeElement,TypeEquation>& A, const Epetra_MultiVector& x, Epetra_MultiVector& y)
    {
      int numvecs = x.NumVectors();
      int m = A.GetM(), n = A.GetN();
      Vector<double> xs(n), ys; // DISP(numvecs);
      for (int k = 0; k < numvecs; k++)
        {
          for (int i = 0; i < n; i++)
            xs(i) = x[k][i];
				
          ys.SetData(m, y[k]);
				
          // DISP(A.eigenvalue_computation_mode); DISP(A.REGULAR_MODE);
          if (A.MassDiagonal())
            {
              if (A.eigenvalue_computation_mode == A.REGULAR_MODE)
                {
                  // DISP(xs); 
                  for (int i = 0; i < n; i++)
                    xs(i) *= A.Dh(i);
						
                  A.MltStiffness(xs, ys);
                  // A.IncrementProdMatVect();
                  for (int i = 0; i < m; i++)
                    ys(i) *= A.Dh(i);
						
                  // DISP(ys);
                }
              else
                {
                  // DISP(xs);
                  for (int i = 0; i < n; i++)
                    xs(i) *= A.Dh(i);
						
                  // A.ComputeSolution(xs, ys);
                  // A.IncrementProdMatVect();
                  for (int i = 0; i < m; i++)
                    ys(i) *= A.Dh(i);
						
                  // DISP(ys);
                }
            }
          else
            {
            }
				
          ys.Nullify();
        }
			
    }
		
    }; */
	
}

#define SELDON_FILE_ANASAZI_HXX
#endif

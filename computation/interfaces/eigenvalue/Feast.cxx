#ifndef SELDON_FILE_FEAST_CXX
#define SELDON_FILE_FEAST_CXX

namespace Seldon
{
  
  // main function to find eigenvalues and eigenvectors with Feast (MKL implementation)
  template<class EigenProblem, class T, class Allocator1,
	   class Allocator2, class Allocator3>
  void FindEigenvaluesFeast(EigenProblem& var,
                            Vector<T, VectFull, Allocator1>& eigen_values,
                            Vector<T, VectFull, Allocator2>& lambda_imag,
                            Matrix<T, General, ColMajor, Allocator3>& eigen_vectors)
  {
    // initialization of feast parameters
    IVect fpm(128);
    fpm.Fill(0);
    feastinit_(fpm.GetData());
    
    fpm(0) = 0;
    if (var.GetPrintLevel() > 0)
      fpm(0) = 1;

    double tol = var.GetStoppingCriterion();    
    fpm(2) = -log10(tol);
    fpm(6) = -log10(tol);
    
    for (int i = 64; i < fpm.GetM(); i++)
      fpm(i) = 0;
    
    fpm(3) = var.GetNbMaximumIterations();
    
    // initialization of computation
    int n = var.GetM();
    var.Init(n);

    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;
        
    int m0 = var.GetNbAskedEigenvalues();

    Treal emin = var.GetLowerBoundInterval();
    Treal emax = var.GetUpperBoundInterval();
    
    T zero; Tcplx one;
    SetComplexZero(zero);
    SetComplexOne(one);

    // work arrays
    Matrix<T, General, ColMajor> work(n, m0);
    Matrix<Tcplx, General, ColMajor> workc(n, m0);
    Vector<Tcplx> xc(n), yc(n);
    Vector<T> aq(m0*m0), sq(m0*m0);
    Vector<T> x(n), y(n);
    
    work.Fill(zero);
    workc.Fill(0);
    xc.Fill(0);
    yc.Fill(0);
    x.Fill(zero);
    y.Fill(zero);
    aq.Fill(zero);
    sq.Fill(zero);

    // output arrays
    Treal epsout(0); Tcplx ze;
    Vector<Treal> res(m0);
    res.Fill(0);
    int loop = 0, m = 0;
    Vector<Treal> lambda(m0);
    lambda.Fill(0);
    eigen_vectors.Reallocate(n, m0);
    eigen_vectors.Fill(zero);
    
    if (var.DiagonalMass())
      {
        // if the mass matrix is diagonal :
        typedef typename EigenProblem::MassValue Tmass;
        Vector<Tmass> Dh_diag;
        // the diagonal is computed
        var.ComputeDiagonalMass(Dh_diag);
        
        // and M^{-1/2} is computed
        var.FactorizeDiagonalMass(Dh_diag);
      }
    else if (var.UseCholeskyFactoForMass())
      {
        // if the user wants to use the Cholesky factorization of the mass matrix:
        
        // computing mass matrix in a convenient form
        var.ComputeMassForCholesky();
        
        // performing Cholesky factorization
        var.FactorizeCholeskyMass();
      }
    else
      var.ComputeMassMatrix();
    
    // evaluating stiffness matrix K
    var.ComputeStiffnessMatrix();

    // main loop (reverse communication interface)
    int ijob = -1, info = 0; DISP(emin); DISP(emax);
    while (ijob != 0)
      { 
        // feast is called
        CallFeast(ijob, n, ze, work, workc, aq, sq, fpm, epsout, loop,
                  emin, emax, m0, lambda, eigen_vectors, m, res, info);
        
        //DISP(ijob); DISP(ze); DISP(loop); DISP(info);
        if (ijob == 10)
          {
            // we have to factorize ze M - K
            var.ComputeAndFactorizeStiffnessMatrix(ze, -one);
          }
        else if (ijob == 11)
          {
            // solves (ze K - M) y = workc(n, m0)
            for (int k = 0; k < m0; k++)
              {
                for (int i = 0; i < n; i++)
                  xc(i) = workc(i, k);
                
                if (var.DiagonalMass())
                  var.MltSqrtDiagonalMass(xc);
                else if (var.UseCholeskyFactoForMass())
                  var.MltCholeskyMass(SeldonNoTrans, xc);
                
                var.ComputeSolution(xc, yc);
                var.IncrementProdMatVect();
                
                if (var.DiagonalMass())
                  var.MltSqrtDiagonalMass(yc);
                else if (var.UseCholeskyFactoForMass())
                  var.MltCholeskyMass(SeldonTrans, yc);
                
                for (int i = 0; i < n; i++)
                  workc(i, k) = yc(i);
              }
          }
        else if (ijob == 20)
          {
            // factorize (ze K - M)^H
            // already done in ijob = 10
          }
        else if (ijob == 21)
          {
            // solves (ze K - M)^H y = workc(n, m0)
            for (int k = 0; k < m0; k++)
              {
                for (int i = 0; i < n; i++)
                  xc(i) = workc(i, k);
                
                if (var.DiagonalMass())
                  var.MltSqrtDiagonalMass(xc);
                else if (var.UseCholeskyFactoForMass())
                  var.MltCholeskyMass(SeldonNoTrans, xc);
                
                Conjugate(xc);
                var.ComputeSolution(SeldonTrans, xc, yc);
                Conjugate(yc);
                
                if (var.DiagonalMass())
                  var.MltSqrtDiagonalMass(yc);
                else if (var.UseCholeskyFactoForMass())
                  var.MltCholeskyMass(SeldonTrans, yc);
                
                for (int i = 0; i < n; i++)
                  workc(i, k) = yc(i);
              }
          }
        else if (ijob == 30)
          {
            // multiplication by matrix A
            int i = fpm(23), j = fpm(23) + fpm(24)-1;
            for (int k = i; k <= j; k++)
              {
                for (int p = 0; p < n; p++)
                  x(p) = eigen_vectors(p, k-1);
                
                if (var.DiagonalMass())
                  var.MltInvSqrtDiagonalMass(x);
                else if (var.UseCholeskyFactoForMass())
                  var.SolveCholeskyMass(SeldonTrans, x);
                
                var.MltStiffness(x, y);
                
                if (var.DiagonalMass())
                  var.MltInvSqrtDiagonalMass(y);
                else if (var.UseCholeskyFactoForMass())
                  var.SolveCholeskyMass(SeldonNoTrans, y);
                
                for (int p = 0; p < n; p++)
                  work(p, k-1) = y(p);
              }
          }
        else if (ijob == 40)
          {
            // multiplication by matrix B
            int i = fpm(23), j = fpm(23) + fpm(24)-1;
            for (int k = i; k <= j; k++)
              {
                for (int p = 0; p < n; p++)
                  x(p) = eigen_vectors(p, k-1);
                
                if (var.DiagonalMass() || var.UseCholeskyFactoForMass())
                  Copy(x, y);
                else
                  var.MltMass(x, y);
                
                for (int p = 0; p < n; p++)
                  work(p, k-1) = y(p);
              }
          }
      }

    if (info != 0)
      {
        cout << "An error occurred during feast eigenvalue solving " << endl;
        cout << "info = " << info << endl;
        abort();
      }
    
    work.Clear();
    workc.Clear();
    
    eigen_values.Reallocate(m);
    lambda_imag.Reallocate(m);
    lambda_imag.Fill(0);
    for (int i = 0; i < m; i++)
      eigen_values(i) = lambda(i);
    
    eigen_vectors.Resize(n, m);
  }
  
}

#endif

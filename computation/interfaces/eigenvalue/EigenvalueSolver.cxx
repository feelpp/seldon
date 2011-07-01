#ifndef SELDON_FILE_EIGENVALUE_SOLVER_CXX

#include "EigenvalueSolver.hxx"

namespace Seldon
{
  
  /******************
   * Initialization *
   ******************/
  
  
  //! default constructor
  template<class T, class MatStiff, class MatMass>
  EigenProblem_Base<T, MatStiff, MatMass>::EigenProblem_Base()
  {
    eigenvalue_computation_mode = 1;
    nb_eigenvalues_wanted = 0;
    // default => we want largest eigenvalues by magnitude
    type_spectrum_wanted = LARGE_EIGENVALUES;
    type_sort_eigenvalues = SORTED_MODULUS;
    
    use_cholesky = false;   
    diagonal_mass = false;
    stopping_criterion = 1e-6;
    nb_maximum_iterations = 1000;
    nb_prod = 0;
    n_ = 0;

    shift = T(0);
    shift_imag = T(0);

    nb_arnoldi_vectors = 0;
    automatic_selection_arnoldi_vectors = true;
    
    print_level = 0;      
    
    complex_system = false;
    Mh = NULL;
    Kh = NULL;
  }
  
  
  //! initialisation of a standard problem
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::Init(int n)
  {
    n_ = n;
    nb_prod = 0;
    if (nb_eigenvalues_wanted >= (n_ - 2))
      {
        cout<<"Too many wanted eigenvalues "<<endl;
        cout << nb_eigenvalues_wanted <<
          " asked eigenvalues, but the rank of the matrix is lower than "
             << n_ << endl;
        
        abort();
      }
      
    if (automatic_selection_arnoldi_vectors)
      nb_arnoldi_vectors = min(n_, 2*nb_eigenvalues_wanted+2);
    
    //cout << "n = " << n << endl;
    //cout << "nb_arnoldi_vectors = " << nb_arnoldi_vectors << endl;
  }
  
  
  //! initialization of a standard eigenvalue problem
  /*!
    Stiffness matrix K is given in argument.
    we will search (lambda, x) such as K x = lambda x
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  InitMatrix(MatStiff& K)
  {
    Kh = &K;
    Mh = NULL;
    this->diagonal_mass = true;
    if ( (!IsSymmetricMatrix(K))
         && (!IsComplexMatrix(K)) && (shift_imag != T(0)) )
      {
        // for real unsymmetric problems, if sigma is complex
        // we have to use mode 3 or 4 in Arpack => generalized problem
        this->diagonal_mass = false;
      }
    
    this->Init(K.GetM());
  }
  
  
  //! initialization of a generalized eigenvalue problem
  /*!
    Mass matrix M and stiffness matrix K are given in argument
    we will search (lambda, x) such as K x = lambda M x
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  InitMatrix(MatStiff& K, MatMass& M)
  {
    Kh = &K;
    Mh = &M;
    this->diagonal_mass = false;
    this->Init(K.GetM());
  }
  
  
  /*******************
   * Basic functions *
   *******************/
  
  
  //! returns the spectral transformation used for evaluation of eigenvalues
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetComputationalMode() const
  {
    return eigenvalue_computation_mode;
  }
  
  
  //! sets the spectral transformation used for evaluation of eigenvalues
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetComputationalMode(int mode)
  {
    eigenvalue_computation_mode = mode;
  }
  
  
  //! returns the number of eigenvalues asked by the user
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetNbAskedEigenvalues() const
  {
    return nb_eigenvalues_wanted;
  }
  
  
  //! sets the number of eigenvalues to compute
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbAskedEigenvalues(int n)
  {
    nb_eigenvalues_wanted = n;
  }
  
  
  //! returns the spectrum desired (large, small eigenvalues, etc)
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetTypeSpectrum() const
  {
    return type_spectrum_wanted;
  }

  
  //! returns how eigenvalues are sorted (real, imaginary part or modulus)
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetTypeSorting() const
  {
    return type_sort_eigenvalues;
  }

  
  //! returns the shift value used
  /*!
    If type_spectrum_wanted is set to CENTERED_EIGENVALUES,
    we search closest eigenvalues to the shift value.
    Matrix (A - (shift + i shift_imag)*I)^{-1} will be used instead of A
  */
  template<class T, class MatStiff, class MatMass>
  T EigenProblem_Base<T, MatStiff, MatMass>::GetShiftValue() const
  {
    return shift;
  }

  
  //! returns the imaginary part of shift value used
  /*!
    If type_spectrum_wanted is set to CENTERED_EIGENVALUES,
    we search closest eigenvalues to the shift value.
    Matrix (A - (shift + i shift_imag)*I)^{-1} will be used instead of A
    shift_imag is accessed only for real unsymmetric problems
  */
  template<class T, class MatStiff, class MatMass>
  T EigenProblem_Base<T, MatStiff, MatMass>::GetImagShiftValue() const
  {
    return shift_imag;
  }
  
  
  //! sets which eigenvalues are searched
  /*!
    You can ask small eigenvalues, large, or eigenvalues
    close to the shift.
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SetTypeSpectrum(int type, const T& val, int type_sort)
  {
    type_spectrum_wanted = type;
    shift = val;
    type_sort_eigenvalues = type_sort;
  }

  
  //! sets which eigenvalues are searched
  /*!
    You can ask small eigenvalues, large, or eigenvalues
    close to the shift.
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SetTypeSpectrum(int type, const complex<T>& val, int type_sort)
  {
    // for real unsymmetric eigenproblems, you can
    // specify a complex shift
    type_spectrum_wanted = type;
    shift = real(val);
    shift_imag = imag(val);

    if (Kh != NULL)
      {
        if ( (!IsSymmetricMatrix(*Kh))
             && (!IsComplexMatrix(*Kh)) && (shift_imag != T(0)) )
          {
            // for real unsymmetric problems, if sigma is complex
            // we have to use mode 3 or 4 in Arpack => generalized problem
            this->diagonal_mass = false;
          }
      }
    
    type_sort_eigenvalues = type_sort;
  }
  
  
  //! indicates the use of Cholesky factorisation in order to 
  //! solve a standard eigenvalue problem instead of a generalized one
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SetCholeskyFactoForMass(bool chol)
  {
    use_cholesky = chol;
  }
  
  
  //! returns true if Cholesky factorisation has to be used for mass matrix
  template<class T, class MatStiff, class MatMass>
  bool EigenProblem_Base<T, MatStiff, MatMass>::
  UseCholeskyFactoForMass() const
  {
    return use_cholesky;
  }
    
  
  //! indicates that the mass matrix is diagonal
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetDiagonalMass(bool diag)
  {
    diagonal_mass = diag;
  }
  
  
  //! returns true if the mass matrix is diagonal
  template<class T, class MatStiff, class MatMass>
  bool EigenProblem_Base<T, MatStiff, MatMass>::DiagonalMass() const
  {
    return diagonal_mass;
  }
  
    
  //! modifies the stopping critertion
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SetStoppingCriterion(double eps)
  {
    stopping_criterion = eps;
  }
  
  
  //! returns the stopping criterion
  template<class T, class MatStiff, class MatMass>
  double EigenProblem_Base<T, MatStiff, MatMass>::
  GetStoppingCriterion() const
  {
    return stopping_criterion;
  }
    
  
  //! sets the maximal number of iterations allowed for the iterative algorithm
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbMaximumIterations(int n)
  {
    nb_maximum_iterations = n;
  }
  
  
  //! returns the maximal number of iterations allowed for the iterative algorithm
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::
  GetNbMaximumIterations() const
  {
    return nb_maximum_iterations;
  }
  
  
  //! returns the number of matrix-vector products performed 
  //! since last call to Init
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::
  GetNbMatrixVectorProducts() const
  {
    return nb_prod;
  }
    
  
  //! returns the number of Arnoldi vectors to use
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetNbArnoldiVectors() const
  {
    return nb_arnoldi_vectors;
  }
  
  
  //! sets the number of Arnoldi vectors to use
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetNbArnoldiVectors(int n)
  {
    automatic_selection_arnoldi_vectors = false;
    nb_arnoldi_vectors = n;
  }
  
  
  //! returns number of rows
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetM() const
  {
    return n_;
  }
    
  
  //! returns number of columns
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetN() const
  {
    return n_;
  }
  
  
  //! returns level of verbosity
  template<class T, class MatStiff, class MatMass>
  int EigenProblem_Base<T, MatStiff, MatMass>::GetPrintLevel() const
  {
    return print_level;
  }
  
  
  //! sets the level of verbosity
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::SetPrintLevel(int lvl)
  {
    print_level = lvl;
  }
  
  
  //! increment of the number of matrix vector products
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::IncrementProdMatVect()
  {
    nb_prod++;
    if (print_level >= 3)
      {
        if (nb_prod%10 == 0)
          cout<<" Iteration number " << nb_prod << endl;
      }
    else if (print_level >= 1)
      {
        if (nb_prod%100 == 0)
          cout<<" Iteration number " << nb_prod << endl;
      }			
  }
  

  //! prints error of initialization and aborts program
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::PrintErrorInit() const
  {
    cout << "InitMatrix has not been called" << endl;
    abort();
  }
  
  
  //! returns true if the matrix is symmetric
  template<class T, class MatStiff, class MatMass>
  bool EigenProblem_Base<T, MatStiff, MatMass>::IsSymmetricProblem() const
  {
    if (Kh != NULL)
      {
        return IsSymmetricMatrix(*Kh);
      }
    else
      PrintErrorInit();
    
    return false;
  }
  
    
  /*********************
   * Mass matrix stuff *
   *********************/
  
  
  //! computation of D^1/2 from D
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  FactorizeDiagonalMass(Vector<MassValue>& D)
  {
    sqrt_diagonal_mass.Reallocate(D.GetM());
    for (int i = 0; i < D.GetM(); i++)
      sqrt_diagonal_mass(i) = sqrt(D(i));
  }
  
  
  //! multiplication of X by D^-1/2
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltInvSqrtDiagonalMass(Vector<T>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) /= sqrt_diagonal_mass(i);
  }
  
  
  //! multiplication of X by D^1/2
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltSqrtDiagonalMass(Vector<T>& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      X(i) *= sqrt_diagonal_mass(i);
  }
    
  
  //! computation of diagonal of mass matrix
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeDiagonalMass(Vector<MassValue>& D)
  {
    if (Mh == NULL)
      {
        // M = identity
        D.Reallocate(this->n_);
        D.Fill(1.0);
      }
    else
      {
        D.Reallocate(this->n_);
        for (int i = 0; i < this->n_; i++)
          D(i) = (*Mh)(i, i);
      }
  }
  
  
  //! computation of mass matrix
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeMassForCholesky()
  {
    // nothing to do, we consider that mass matrix
    // is already computed
  }
  
  
  //! computation of mass matrix M
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::ComputeMassMatrix()
  {
    // mass matrix already computed in Mh
  }
  
  
  //! matrix vector product with mass matrix Y = M X
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltMass(const Vector<T>& X, Vector<T>& Y)
  {
    if (Mh == NULL)
      {
        // default : mass matrix is identity (standard eigenvalue problem)
        Seldon::Copy(X, Y);
      }
    else
      Mlt(*Mh, X, Y);
  }
  
  
  /**************************
   * Stiffness matrix stuff *
   **************************/
  
  
  //! computation of stiffness matrix K
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::ComputeStiffnessMatrix()
  {
    // nothing to do, already computed in Kh
  }
  
  
  //! computation of matrix a M + b*K
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeStiffnessMatrix(const T& a, const T& b)
  {
    // nothing to do, we use Kh and Mh for the matrix vector product
  }
	
  
  //! matrix vector product with stifness matrix Y = K X
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltStiffness(const Vector<T>& X, Vector<T>& Y)
  {
    if (Kh == NULL)
      PrintErrorInit();
    else
      Mlt(*Kh, X, Y);
  }
  
  
  //! matrix vector product with stifness and mass matrix Y = (a M + b K) X
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltStiffness(const T& coef_mass, const T& coef_stiff,
               const Vector<T>& X, Vector<T>& Y)
  {
    if (Kh == NULL)
      PrintErrorInit();
    else
      Mlt(*Kh, X, Y);
    
    if (coef_mass != T(0))
      {
        if (Mh == NULL)
          for (int i = 0; i < Y.GetM(); i++)
            Y(i) += coef_mass*X(i);
        else
          MltAdd(coef_mass, *Mh, X, coef_stiff, Y);
      }
    else
      {
        if (coef_stiff != T(1))
          Mlt(coef_stiff, Y);
      }
  }
    
  
  
  /*************************
   * Functions to overload *
   *************************/

  
  //! computation of matrix a M + b K and factorisation of this matrix
  /*!
    The factorisation process can be also the construction of preconditioning
    if an iterative solver is used to solve linear system (a M + b K) y = x 
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b)
  {
    abort();
  }
  
  
  //! computation of matrix a M + b K and factorisation of this matrix
  /*!
    The factorisation process can be also the construction of preconditioning
    if an iterative solver is used to solve linear system (a M + b K) y = x 
  */
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const complex<T>& a, const complex<T>& b,
                                     bool real_part)
  {
    abort();
  }
  
  
  //! solving the linear system (a M + b K) Y = X
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  ComputeSolution(const Vector<T>& X, Vector<T>& Y)
  {
    abort();
  }
  
  
    //! computation of Cholesky factorisation of M from matrix M
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::FactorizeCholeskyMass()
  {
    abort();
  }
  
  
  //! computation of L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  MltCholeskyMass(const TransStatus& TransA, Vector<T>& X)
  {
    abort();
  }
  
  
  //! computation of L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void EigenProblem_Base<T, MatStiff, MatMass>::
  SolveCholeskyMass(const TransStatus& TransA, Vector<T>& X)
  {
    abort();
  }
  
  
  //! memory release
  template<class T, class MatStiff, class MatMass>
  void EigenProblem_Base<T, MatStiff, MatMass>::Clear()
  {
    sqrt_diagonal_mass.Clear();
  }
  
  
  /********************************************
   * Modification of eigenvalues/eigenvectors *
   ********************************************/
  
  
  //! modification of eigenvectors to take into account scaling by mass matrix
  /*!
    One may desire to use matrix D^-1/2 K D^-1/2 or L^-1 K L^-T instead of K
    in order to solve a standard eigenvalue problem instead of a generalized one.
    => eigenvectors are recovered by multiplying them by matrix D^1/2 or by L^T
      with this function
  */
  template<class EigenPb, class Vector1, class Matrix1, class T0>
  void ApplyScalingEigenvec(EigenPb& var, Vector1& eigen_values, Vector1& lambda_imag,
                            Matrix1& eigen_vectors,
                            const T0& shiftr, const T0& shifti)
  {
    
    if (var.DiagonalMass())
      {
        // scaling to have true eigenvectors
        for (int i = 0; i < var.sqrt_diagonal_mass.GetM(); i++)
          for (int j = 0; j < eigen_vectors.GetN(); j++)
            eigen_vectors(i,j) /= var.sqrt_diagonal_mass(i);
      }      
    else if (var.UseCholeskyFactoForMass())
      {
        Vector<T0> Xcol(eigen_vectors.GetM());
        for (int j = 0; j < eigen_vectors.GetN(); j++)
          {
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              Xcol(i) = eigen_vectors(i,j);
            
            var.SolveCholeskyMass(SeldonTrans, Xcol);
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              eigen_vectors(i,j) = Xcol(i);
          }
      }
    
    if (var.GetComputationalMode() != var.REGULAR_MODE)
      {
        if ((var.DiagonalMass())|| (var.UseCholeskyFactoForMass())
            || (var.eigenvalue_computation_mode == var.INVERT_MODE))
          {
            // shift-invert mode, we have to modify eigenvalues
            for (int i = 0; i < eigen_values.GetM(); i++)
              {
                if ((eigen_values(i) == 0) && (lambda_imag(i) == 0))
                  {
                    eigen_values(i) = 0;
                    lambda_imag(i) = 0;
                  }
                else
                  {
                    complex<T0> val = 1.0/complex<T0>(eigen_values(i), lambda_imag(i))
                      + complex<T0>(shiftr, shifti);
                    
                    eigen_values(i) = real(val);
                    lambda_imag(i) = imag(val);
                  }
              }
        
          }
        else if (var.GetImagShiftValue() != T0(0))
          {
            int n = eigen_vectors.GetM();
            Vector<T0> X(n), Ax(n), Mx(n), Y(n);
            int j = 0;
            Ax.Fill(T0(0));
            Mx.Fill(T0(0));
            while (j < eigen_values.GetM())
              {
                if (lambda_imag(j) == T0(0))
                  {
                    // real eigenvalue
                    // lambda is retrieve by computing Rayleigh quotient
                    for (int i = 0; i < eigen_vectors.GetM(); i++)
                      X(i) = eigen_vectors(i,j);
                    
                    var.MltMass(X, Mx);
                    var.MltStiffness(X, Ax);
                    eigen_values(j) = DotProd(X, Ax)/DotProd(X, Mx);
                    
                    // next eigenvalue
                    j++;
                  }
                else
                  {
                    // conjugate pair of eigenvalues
                    for (int i = 0; i < eigen_vectors.GetM(); i++)
                      {
                        X(i) = eigen_vectors(i, j);
                        Y(i) = eigen_vectors(i, j+1);
                      }
                    
                    // complex Rayleigh quotient
                    var.MltStiffness(X, Ax);
                    T0 numr = DotProd(X, Ax);
                    T0 numi = DotProd(Y, Ax);
                    
                    var.MltStiffness(Y, Ax);
                    numr += DotProd(Y, Ax);
                    numi -= DotProd(X, Ax);
                    
                    var.MltMass(X, Mx);
                    T0 denr = DotProd(X, Mx);
                    T0 deni = DotProd(Y, Mx);
                    
                    var.MltMass(Y, Mx);
                    denr += DotProd(Y, Mx);
                    deni -= DotProd(X, Mx);
                    
                    complex<T0> val = complex<T0>(numr, numi)/complex<T0>(denr, deni);
                    
                    eigen_values(j) = real(val);
                    eigen_values(j+1) = real(val);

                    lambda_imag(j) = -imag(val);
                    lambda_imag(j+1) = imag(val);
                    
                    // next eigenvalue
                    j += 2;
                  }
              }
          }
      }
  }


  //! modification of eigenvectors to take into account scaling by mass matrix
  /*!
    One may desire to use matrix D^-1/2 K D^-1/2 or L^-1 K L^-T instead of K
    in order to solve a standard eigenvalue problem instead of a generalized one.
    => eigenvectors are recovered by multiplying them by matrix D^1/2 or by L^T
      with this function
  */
  template<class EigenPb, class Vector1, class Matrix1, class T0>
  void ApplyScalingEigenvec(EigenPb& var, Vector1& eigen_values, Vector1& lambda_imag,
                            Matrix1& eigen_vectors,
                            const complex<T0>& shiftr, const complex<T0>& shifti)
  {
    
    if (var.DiagonalMass())
      {
        // scaling to have true eigenvectors
        for (int i = 0; i < var.sqrt_diagonal_mass.GetM(); i++)
          for (int j = 0; j < eigen_vectors.GetN(); j++)
            eigen_vectors(i,j) /= var.sqrt_diagonal_mass(i);
      }      
    else if (var.UseCholeskyFactoForMass())
      {
        Vector<complex<T0> > Xcol(eigen_vectors.GetM());
        for (int j = 0; j < eigen_vectors.GetN(); j++)
          {
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              Xcol(i) = eigen_vectors(i,j);
            
            var.SolveCholeskyMass(SeldonTrans, Xcol);
            for (int i = 0; i < eigen_vectors.GetM(); i++)
              eigen_vectors(i,j) = Xcol(i);
          }
      }
    
    if (var.GetComputationalMode() != var.REGULAR_MODE)
      {
        if ((var.DiagonalMass())|| (var.UseCholeskyFactoForMass())
            || (var.eigenvalue_computation_mode == var.INVERT_MODE))
          {
            // shift-invert mode, we have to modify eigenvalues
            for (int i = 0; i < eigen_values.GetM(); i++)
              {
                complex<T0> val = 1.0/eigen_values(i) + shiftr;
                
                eigen_values(i) = val;
              }
            
          }
      }
  }

  
  /***********************
   * Sorting eigenvalues *
   ***********************/
  
  
  //! sorting eigenvalues
  template<class T, class Storage1, class Storage2,
           class Alloc1, class Alloc2>
  void SortEigenvalues(Vector<T>& lambda_r, Vector<T>& lambda_i,
                       Matrix<T, General, Storage1, Alloc1>& eigen_old,
                       Matrix<T, General, Storage2, Alloc2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const T& shift_r, const T& shift_i)
  {
    int n = min(lambda_r.GetM(), eigen_old.GetN());;
				 
    IVect permutation(n);
    permutation.Fill();
    eigen_new.Reallocate(eigen_old.GetM(), n);
    
    // creating a vector that can be sorted
    Vector<T>  L(n);
    if (type_spectrum == EigenProblem_Base<T>::CENTERED_EIGENVALUES)
      {
        // eigenvalues closest to shift are placed at beginning
        switch (type_sort)
          {
          case EigenProblem_Base<T>::SORTED_REAL :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(shift_r - lambda_r(i));
            }
            break;
          case EigenProblem_Base<T>::SORTED_IMAG :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(shift_i - lambda_i(i));
            }
            break;
          case EigenProblem_Base<T>::SORTED_MODULUS :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(complex<T>(shift_r - lambda_r(i), shift_i - lambda_i(i)));
            }
            break;
          }
      }
    else
      {
        // smallest eigenvalues are placed at beginning
        switch (type_sort)
          {
          case EigenProblem_Base<T>::SORTED_REAL :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(lambda_r(i));
            }
            break;
          case EigenProblem_Base<T>::SORTED_IMAG :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(lambda_i(i));
            }
            break;
          case EigenProblem_Base<T>::SORTED_MODULUS :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(complex<T>(lambda_r(i), lambda_i(i)));
            }
            break;            
          }
      }
    
    // sorting L, and retrieving permutation array
    Sort(L, permutation);
				 
    // permuting eigenvalues and eigenvectors
    Vector<T> oldLambda_r = lambda_r, oldLambda_i = lambda_i;
    for (int i = 0; i < n; i++)
      {
        lambda_r(i) = oldLambda_r(permutation(i));
        lambda_i(i) = oldLambda_i(permutation(i));
        for (int j = 0; j < eigen_old.GetM(); j++)
          eigen_new(j, i) = eigen_old(j, permutation(i));
      }
    
  }


  //! sorting eigenvalues
  template<class T, class Storage1, class Storage2,
           class Alloc1, class Alloc2>
  void SortEigenvalues(Vector<complex<T> >& lambda_r, Vector<complex<T> >& lambda_i,
                       Matrix<complex<T>, General, Storage1, Alloc1>& eigen_old,
                       Matrix<complex<T>, General, Storage2, Alloc2>& eigen_new,
                       int type_spectrum, int type_sort,
                       const complex<T>& shift_r, const complex<T>& shift_i)
  {
    // complex case, ignoring lambda_i and shift_i
    int n = min(lambda_r.GetM(), eigen_old.GetN());;
    
    IVect permutation(n);
    permutation.Fill();
    eigen_new.Reallocate(eigen_old.GetM(), n);
    
    // creating a vector that can be sorted
    Vector<T>  L(n);
    if (type_spectrum == EigenProblem_Base<T>::CENTERED_EIGENVALUES)
      {
        // eigenvalues closest to shift are placed at beginning
        switch (type_sort)
          {
          case EigenProblem_Base<T>::SORTED_REAL :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(real(shift_r - lambda_r(i)));
            }
            break;
          case EigenProblem_Base<T>::SORTED_IMAG :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(imag(shift_r - lambda_r(i)));
            }
            break;
          case EigenProblem_Base<T>::SORTED_MODULUS :
            {
              for (int i = 0; i < n; i++)
                L(i) = 1.0/abs(shift_r - lambda_r(i));
            }
            break;
          }
      }
    else
      {
        // smallest eigenvalues are placed at beginning
        switch (type_sort)
          {
          case EigenProblem_Base<T>::SORTED_REAL :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(real(lambda_r(i)));
            }
            break;
          case EigenProblem_Base<T>::SORTED_IMAG :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(imag(lambda_r(i)));
            }
            break;
          case EigenProblem_Base<T>::SORTED_MODULUS :
            {
              for (int i = 0; i < n; i++)
                L(i) = abs(lambda_r(i));
            }
            break;            
          }
      }
    
    // sorting L, and retrieving permutation array
    Sort(L, permutation);
				 
    // permuting eigenvalues and eigenvectors
    Vector<complex<T> > oldLambda_r = lambda_r, oldLambda_i = lambda_i;
    for (int i = 0; i < n; i++)
      {
        lambda_r(i) = oldLambda_r(permutation(i));
        lambda_i(i) = oldLambda_i(permutation(i));
        for (int j = 0; j < eigen_old.GetM(); j++)
          eigen_new(j, i) = eigen_old(j, permutation(i));
      }
    
  }
  
  
  /*********************
   * DenseEigenProblem *
   *********************/
  
  
  //! default constructor
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  DenseEigenProblem() : EigenProblem_Base<T, Matrix<T, Prop, Storage>,
                                          Matrix<Tmass, PropM, StorageM> >()
  {    
  }
  
  
  //! Cholesky factorisation of mass matrix
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  FactorizeCholeskyMass()
  {
    if (this->Mh == NULL)
      {
        mat_chol.Reallocate(this->n_, this->n_);
        mat_chol.SetIdentity();
      }
    else
      {
        mat_chol = *(this->Mh);
        GetCholesky(mat_chol);
        Xchol_real.Reallocate(mat_chol.GetM());
        Xchol_imag.Reallocate(mat_chol.GetM());
      }
  }
  
  
  //! computation of L X or L^T X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransStatus, class T0>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  MltCholeskyMass(const TransStatus& transA, Vector<T0>& X)
  {
    MltCholesky(transA, mat_chol, X);
  }
  
  
  //! computation of L^-1 X or L^-T X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransStatus, class T0>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  SolveCholeskyMass(const TransStatus& transA, Vector<T0>& X)
  {
    SolveCholesky(transA, mat_chol, X);
  }
  
  
  //! computation of L X or L^T X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransStatus>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  MltCholeskyMass(const TransStatus& transA, Vector<complex<double> >& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    MltCholesky(transA, mat_chol, Xchol_real);
    MltCholesky(transA, mat_chol, Xchol_imag);

    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<double>(Xchol_real(i), Xchol_imag(i));    
  }
  
  
  //! computation of L^-1 X or L^-T X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  template<class TransStatus>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  SolveCholeskyMass(const TransStatus& transA, Vector<complex<double> >& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    SolveCholesky(transA, mat_chol, Xchol_real);
    SolveCholesky(transA, mat_chol, Xchol_imag);

    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<double>(Xchol_real(i), Xchol_imag(i));
  }

  
  //! computation and factorisation of a M + b K
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b)
  {
    if (this->Kh == NULL)
      this->PrintErrorInit();
    
    this->complex_system = false;
    // computation of mat_lu = a M + b K
    mat_lu = *(this->Kh);
    Mlt(b, mat_lu);
    if (this->Mh == NULL)
      {
        for (int i = 0; i < this->n_; i++)
          mat_lu(i, i) += a;
      }
    else
      Add(a, *(this->Mh), mat_lu);
    
    // factorisation
    GetLU(mat_lu, pivot);
  }
  
  
  //! computation and factorisation of a M + b K
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeStiffnessMatrix(const complex<T>& a, const complex<T>& b,
                                     bool real_part)
  {
    ComputeAndFactorizeComplexMatrix(a, b);
  }
  
  
  //! computation and factorisation of a M + b K
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeComplexMatrix(const complex<double>& a, const complex<double>& b,
                                   bool real_part)
  {
    if (this->Kh == NULL)
      this->PrintErrorInit();
 
    this->complex_system = true;
    // inverse of (a M + b K), then we take real_part or imaginary part
    Matrix<complex<double>, Prop, Storage> InvMat(this->n_, this->n_);
    for (int i = 0; i < this->n_; i++)
      for (int j = 0; j < this->n_; j++)
        InvMat(i, j) = (*this->Kh)(i, j);
    
    Mlt(b, InvMat);
    if (this->Mh == NULL)
      {
        for (int i = 0; i < this->n_; i++)
          InvMat(i, i) += a;
      }
    else
      {
        for (int i = 0; i < this->n_; i++)
          for (int j = 0; j < this->n_; j++)
            InvMat(i, j) += a * (*this->Mh)(i, j);
      }
    
    // inverse
    GetInverse(InvMat);
  
    // then extracting real or imaginary part
    mat_lu.Reallocate(this->n_, this->n_);
    if (real_part)
      {
        for (int i = 0; i < this->n_; i++)
          for (int j = 0; j < this->n_; j++)
            mat_lu(i, j) = real(InvMat(i, j));
      }
    else
      {
        for (int i = 0; i < this->n_; i++)
          for (int j = 0; j < this->n_; j++)
            mat_lu(i, j) = imag(InvMat(i, j));
      }
  }
  

  //! computation and factorisation of a M + b K
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM> template<class T0>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeAndFactorizeComplexMatrix(const complex<T0>& a,
                                   const complex<T0>& b, bool real_p)
  {
    // this case should not appear
    cout << "Case not handled" << endl;
    abort();
  }
  
  
  //! solution of (a M + b K) Y = X
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  ComputeSolution(const Vector<T>& X, Vector<T>& Y)
  {
    if (this->complex_system)
      Mlt(mat_lu, X, Y);
    else
      {
        Copy(X, Y);
        SolveLU(mat_lu, pivot, Y);
      }
  }
   
  
  //! clearing variables used for eigenvalue resolution
  template<class T, class Prop, class Storage,
           class Tmass, class PropM, class StorageM>
  void DenseEigenProblem<T, Prop, Storage, Tmass, PropM, StorageM>::
  Clear()
  {
    EigenProblem_Base<T, Matrix<T, Prop, Storage>,
      Matrix<Tmass, PropM, StorageM> >::Clear();
    
    mat_lu.Clear();
    mat_chol.Clear();
  }
  
  
  /**********************
   * SparseEigenProblem *
   **********************/
  
  
  //! default constructor
  template<class T, class MatStiff, class MatMass>
  SparseEigenProblem<T, MatStiff, MatMass>::SparseEigenProblem()
    : EigenProblem_Base<T, MatStiff, MatMass>()
  {
    imag_solution = false;
    mat_lu.RefineSolution();
    mat_lu_cplx.RefineSolution();
  }
    
  
  //! computation of Cholesky factorisation of M from matrix M
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  FactorizeCholeskyMass()
  {
    if (this->Mh == NULL)
      this->PrintErrorInit();
    
    if (this->print_level > 0)
      chol_facto_mass_matrix.ShowMessages();
    
    chol_facto_mass_matrix.Factorize(*this->Mh, true);    
    
    if (this->print_level < 2)
      chol_facto_mass_matrix.HideMessages();
    
    Xchol_real.Reallocate(this->n_);
    Xchol_imag.Reallocate(this->n_);
  }
  
  
  //! computation of L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus, class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  MltCholeskyMass(const TransStatus& TransA, Vector<T0>& X)
  {
    chol_facto_mass_matrix.Mlt(TransA, X);
  }
  
  
  //! computation of L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus, class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  SolveCholeskyMass(const TransStatus& TransA, Vector<T0>& X)
  {
    chol_facto_mass_matrix.Solve(TransA, X);
  }


  //! computation of L X or L^T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  MltCholeskyMass(const TransStatus& TransA, Vector<complex<double> >& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    chol_facto_mass_matrix.Mlt(TransA, Xchol_real);
    chol_facto_mass_matrix.Mlt(TransA, Xchol_imag);
    
    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<double>(Xchol_real(i), Xchol_imag(i));
    
  }
  
  
  //! computation of L^-1 X or L^-T x if M = L L^T
  template<class T, class MatStiff, class MatMass>
  template<class TransStatus>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  SolveCholeskyMass(const TransStatus& TransA, Vector<complex<double> >& X)
  {
    for (int i = 0; i < X.GetM(); i++)
      {
        Xchol_real(i) = real(X(i));
        Xchol_imag(i) = imag(X(i));
      }
    
    chol_facto_mass_matrix.Solve(TransA, Xchol_real);
    chol_facto_mass_matrix.Solve(TransA, Xchol_imag);
    
    for (int i = 0; i < X.GetM(); i++)
      X(i) = complex<double>(Xchol_real(i), Xchol_imag(i));

  }
  
  
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const T& a, const T& b)
  {
    this->complex_system = false;
    if (this->Kh == NULL)
      this->PrintErrorInit();
    
    if (this->print_level > 0)
      mat_lu.ShowMessages();
    
    // forming a M + b K
    if (IsSymmetricMatrix(*this->Kh))
      {
        Matrix<T, Symmetric, ArrayRowSymSparse> A;
        Copy(*(this->Kh), A);
        Mlt(b, A);
        if (this->Mh == NULL)
          {
            for (int i = 0; i < this->n_; i++)
              A.AddInteraction(i, i, a);
          }
        else
          Add(a, *(this->Mh), A);
        
        mat_lu.Factorize(A);
      } 
    else
      {
        Matrix<T, General, ArrayRowSparse> A;
        Copy(*(this->Kh), A);
        Mlt(b, A);
        if (this->Mh == NULL)
          {
            for (int i = 0; i < this->n_; i++)
              A.AddInteraction(i, i, a);
          }
        else
          Add(a, *(this->Mh), A);
        
        mat_lu.Factorize(A);
      }      

    if (this->print_level < 2)
      mat_lu.HideMessages();
  }
  
  
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeStiffnessMatrix(const complex<T>& a,
                                     const complex<T>& b,
                                     bool real_p)
  {
    ComputeAndFactorizeComplexMatrix(a, b, real_p);
  }
  
  
  template<class T, class MatStiff, class MatMass> template<class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeComplexMatrix(const complex<T0>& a,
                                   const complex<T0>& b,
                                   bool real_p)
  {
    // this case should not appear
    cout << "Case not handled" << endl;
    abort();
  }

  
  template<class T, class MatStiff, class MatMass> 
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeAndFactorizeComplexMatrix(const complex<double>& a,
                                   const complex<double>& b,
                                   bool real_p)
  {
    this->complex_system = true;
    if (this->Kh == NULL)
      this->PrintErrorInit();
    
    imag_solution = !real_p;
    // forming a M + b K
    Matrix<Complexe, General, ArrayRowSparse> A;
    Copy(*(this->Kh), A);
    Mlt(b, A);
    if (this->Mh == NULL)
      {
        for (int i = 0; i < this->n_; i++)
          A.AddInteraction(i, i, a);
      }
    else
      Add(a, *(this->Mh), A);
    
    if (this->print_level > 0)
      mat_lu_cplx.ShowMessages();

    mat_lu_cplx.Factorize(A);

    if (this->print_level < 2)
      mat_lu_cplx.HideMessages();
  }

  
  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeSolution(const Vector<T>& X, Vector<T>& Y)
  {
    if (this->complex_system)
      {
        ComputeComplexSolution(X, Y);
      }
    else
      {
        Copy(X, Y);
        mat_lu.Solve(Y);
      }
  }
  

  template<class T, class MatStiff, class MatMass> template<class T0>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeComplexSolution(const Vector<T0>& X, Vector<T0>& Y)
  {
    Vector<complex<T0> > Xcplx(this->n_);
    for (int i = 0; i < this->n_; i++)
      Xcplx(i) = X(i);
    
    mat_lu_cplx.Solve(Xcplx);
    
    if (imag_solution)
      for (int i = 0; i < this->n_; i++)
        Y(i) = imag(Xcplx(i));
    else
      for (int i = 0; i < this->n_; i++)
        Y(i) = real(Xcplx(i));
    
  }


  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::
  ComputeComplexSolution(const Vector<complex<double> >& X,
                         Vector<complex<double> >& Y)
  {
    // this case should not appear
    cout << "Case not handled" << endl;
    abort();
  }


  template<class T, class MatStiff, class MatMass>
  void SparseEigenProblem<T, MatStiff, MatMass>::Clear()
  {
    mat_lu.Clear();
    mat_lu_cplx.Clear();
    chol_facto_mass_matrix.Clear();
    Xchol_real.Clear();
    Xchol_imag.Clear();    
  }
  
#ifndef SELDON_WITH_COMPILED_LIBRARY
  int TypeEigenvalueSolver::default_solver(0);  
#endif
  
}

#define SELDON_FILE_EIGENVALUE_SOLVER_CXX
#endif

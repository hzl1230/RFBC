#ifndef POISSON_H
#define POISSON_H

#include <vector>
#include <utility>
#include <iomanip>

#include "espic_type.h"

#ifdef USE_MKL
#include "Eigen/PardisoSupport"
#endif

enum PoissonSolverType { LU, PCG };

class Poisson {
  public:
#ifdef USE_MKL
    typedef SpMatCSR SpMat;   // CSR format
#else
    typedef SpMatCSC SpMat;   // CSC format
#endif

    /* Constructors */
    Poisson(class Mesh*);

    /* Desctructor */
    ~Poisson();

    /* Public methods */

    struct Result {       // returned value by Poisson's solver
      bool linear;        // true if linear solver
      bool iterative;     // true if iterative solver
      bool converged;     // true if solved successfully
      int iterations;     // # of iterations in the end
      int iterations_newton;  // # of iterations of newton-raphson
      Real error;
      Real tolerance;
    };
    Result result;

    const class Mesh& get_mesh() const { return mesh; }
    const Vector& solution() const { return sol; }

    /* Driver to Poisson's solve function */
    void solve(const class ScalarField*, class ScalarField*);

    void showinstantphi(Real& phi);

  private:
    const class Mesh& mesh;
    bool   linear;        // true if linear solver
    int    solver_type;   // LU, PCG
    int    solver_maxit;  // max iterations for iterative solver
    int    newton_maxit;  // max iterations for newton-raphson of nonlinear solvers
    Real   solver_rtol;   // residual tolerance
    Real   newton_rtol;
    bool   isothermal;    // true if isothermal, ignored if linear is true
    Real   gamma;         // polytropic coefficient
    Real   phi_ref;       // potential at ref point for hybrid simulation
    Real   n_ref;         // electron density at ref point for hybrid simulation
    Real   temp_ref;      // electron temperature at ref point for hybrid simulation
    Vector b;             // right-hand-side
    Vector sol;           // solution 
    SpMat  A;             // sparse matrix A
    SpMat  J;             // sparse matrix J, for non-linear solver
    std::vector<bool> is_fixed;
    std::vector< std::pair<Index, Real> > bc_id_val[2];  // Dirichlet & Neumann node ids and values
    Real timestep;
    Real Time;
    std::vector<bool> is_rf;

    /* LU and PCG solver handle */
#ifdef USE_MKL
    Eigen::PardisoLU< SpMat>* ptrLUSolver;
#else
    Eigen::SparseLU< SpMat >* ptrLUSolver;
#endif
    Eigen::BiCGSTAB< SpMat, Eigen::IncompleteLUT<Real> >* ptrBiCGSolver;

    typedef void (Poisson::*PtrSolve)(const class ScalarField*, class ScalarField*);  // ptr to solve function
    PtrSolve ptrsolve;

    /* RHS and Jacobian setting handle for nonlinear solver */
    typedef void (Poisson::*PtrSetRhsNonlin)(const class ScalarField&, const Vector&);  // ptr to set_rhs_nonlin
    PtrSetRhsNonlin ptrsetRhs;

    typedef void (Poisson::*PtrSetJacobianNonlin)(const Vector&);  // ptr to set_jacobian_nonlin
    PtrSetJacobianNonlin ptrsetJac;

    void init();
    void init_lu();           // initialize LU solver
    // void init_pcg();          // initialize PCG solver

    void set_matrix_2d();     // assemble coeff matrix for 2d simulations 
    void set_matrix_axi();    // assemble coeff matrix for axi-symm simulations 

    void check_dirichlet_nodes_2d(std::vector<bool>&); 
    void check_objects_2d(std::vector<bool>&); 
    void check_boundary(const Index, const Index, class std::vector<int>&);

    void lu_solve_lin(const class ScalarField*, class ScalarField*);    // solve Ax=b for x by LU

    void set_rhs_lin(const class ScalarField&);

    // Discharge model
    // void SetSimpleCCP(Real dt);

};

inline std::ostream& operator<< (std::ostream& os, const Poisson::Result& r)
{
  if (r.linear) {
    if (r.iterative) {
      std::string msg = r.converged ? "converged" : "failed to converge";
      os << "Iterative solver: " << msg 
        << " after " << std::setw(3) << r.iterations << " iters, err/tol = ";
      os.precision(1);
      os << std::fixed << std::scientific << r.error << "/" << r.tolerance;
    }
    else {
      std::string msg = r.converged ? "succeeded!" : "failed!";
      os << "Direct solver: " << msg;
    }
  }
  else {
    std::string msg = r.converged ? "converged" : "failed to converge";
    os << "Newton iterations: " << msg 
      << " after " << std::setw(2) << r.iterations_newton << " iters, err/tol = ";
    os.precision(1);
    os << std::fixed << std::scientific << r.error << "/" << r.tolerance << "; ";

    if (r.iterative) {
      os << "Iterative solver: total iters = "
        << std::setw(3) << r.iterations;
    }
    else {
      std::string msg = r.converged ? "succeeded!" : "failed!";
      os << "Direct solver: " << msg;
    }
  }

  return os;
}

#endif

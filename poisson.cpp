#include <iostream>
#include <functional>
#include <cstring>
#include <cmath>
#include <limits>

#include "espic_info.h"
#include "espic_math.h"
#include "poisson.h"
#include "mesh.h"
#include "modeldef.h"
#include "solverdef.h"
#include "scalar_field.h"

enum BoundaryId {XLO, XHI, YLO, YHI, ZLO, ZHI};
Poisson::Poisson(Mesh* m)
  : mesh(*m), 
    linear(true), 
    solver_type(0),
    solver_maxit(1000),
    newton_maxit(100),
    solver_rtol(1e-4),
    newton_rtol(1e-6),
    isothermal(true),
    gamma(1.),
    phi_ref (0.),
    n_ref (1.),
    temp_ref(1.),
    b(m->num_nodes()),
    sol(m->num_nodes()),
    A(m->num_nodes(), m->num_nodes()), 
    J(), 
    is_fixed(m->num_nodes(), false),
    timestep(0.1),
    is_rf(m->num_nodes(), false),
    ptrLUSolver(nullptr), 
    ptrBiCGSolver(nullptr), 
    ptrsolve(nullptr),
    ptrsetRhs(nullptr),
    ptrsetJac(nullptr)
{
  std::cout << "Intializing Poisson's solver: \n";
  solver_type = LU;
  init();
  Time = 0.;
}

/* Destructor */
Poisson::~Poisson()
{
  if (ptrLUSolver   != nullptr) delete ptrLUSolver;
  if (ptrBiCGSolver != nullptr) delete ptrBiCGSolver;
}

void Poisson::solve(const ScalarField* r, ScalarField* x)
{
  // solve Ax = b for x
  return ( (this->*ptrsolve)(r, x) );
}

void Poisson::showinstantphi(Real& phi)
{
  ScalarField* sf = new ScalarField(mesh);
  set_rhs_lin(*sf);
}

/* ----------------- End Public Methods ----------------- */


/* ---------------- Begin Private Methods ---------------- */

void Poisson::init()
{
  result.linear = linear;
  result.iterative = !(solver_type == LU);
  result.iterations_newton = 0;
  result.converged = false;
  if (result.iterative) {
    result.iterations = solver_maxit;
    result.error = -99999;
    result.tolerance = solver_rtol;
  }
  else {
    result.iterations = 0;
    result.error = 0;
    result.tolerance = 0;
  }

  if (2 == mesh.dimension()) { set_matrix_2d(); }
  else if (5 == mesh.dimension()) { set_matrix_axi(); }
  else if (3 == mesh.dimension()) {
    espic_error("3D setup is NOT available yet");
  }

  switch (solver_type) {
    case LU :
      init_lu();
      break;
    // case PCG :
    //   init_pcg();
    //   break;
    default :
      espic_error("Unknown type of Poisson's solver");
  }

  // initialization for nonlinear solver
  if (!linear) {
      espic_error("Non-Linear solver is not supported yet");
  }
}

/* ------------------------------------------------------- */

void Poisson::set_matrix_2d()
{
  /* Begin: Find and set all nodes with fixed value */
  check_dirichlet_nodes_2d(is_fixed);

  check_objects_2d(is_fixed);

  std::cout << "Assembling coefficient matrix for two-dimensional Cartesian mesh.\n";

  /* --------------------------------------------------------- */
  /* ----------- Compute coefficients for matrix A ----------- */
  /* --------------------------------------------------------- */
  Index nd[3] = {mesh.num_nodes(0), mesh.num_nodes(1), mesh.num_nodes(2)};  // # of nodes
  Index nc[3] = {mesh.num_cells(0), mesh.num_cells(1), mesh.num_cells(2)};  // # of cells
  Index nid_pshift[3] = {nc[0], nc[1]*nd[0], nc[2]*nd[1]*nd[0]};   // # of nodes shifted when handling periodic BC
  const int nelem = 5;        // 2D: 5-point discretization scheme
  Index row, col[nelem];      // row, col ids for non-zero entries in a row of matrix A 
  Real val[nelem];            // values for non-zero entries in a row of matrix A
  std::vector<Real> dx[3];    // cell sizes
  std::vector<Tp> Atp;        // all non-zero entries of matrix A in the Triplet form
  const Index NOT_USE = std::numeric_limits<Index>::min();    // a flag for unused nodes

  std::function < Real(Index) > fdx = [&] (Index i) -> Real { return mesh.dx(i,0); };
  std::function < Real(Index) > fdy = [&] (Index i) -> Real { return mesh.dy(0,i); };
  std::function < Real(Index) > fdz = [&] (Index i) -> Real { return mesh.dz(0,0); };
  decltype(fdx) fa[3] = { fdx, fdy, fdz };

  // compute cell sizes
  Real dx_bc[3][2];     // dx at boundaries, used for matrix assembly
  for (int a = 0; a < 3; a++) { 
    dx[a].resize(nc[a]+2);
    for (Index i = -1; i <= nc[a]; ++i) dx[a][i+1] = fa[a](i);

    dx_bc[a][0] = 2*dx[a][0];       // "LO" boundary in X, Y, Z
    dx_bc[a][1] = 2*dx[a][nc[a]];   // "HI" boundary in X, Y, Z
  }

  // 2D case; row/x-major storage of mesh nodes
  std::vector<int> bound;
  for (Index iy = 0; iy < nd[1]; iy++) {
    for (Index ix = 0; ix < nd[0]; ix++) {
      row = iy*nd[0] + ix;    // node id, nid

      // nodes with fixed potential 
      if (is_fixed[row]) { Atp.push_back(Tp(row, row, 1.0)); continue; }

      check_boundary(ix, iy, bound);

      Real c;
      col[0] = row - 1;     col[1] = row + 1;
      col[2] = row - nd[0]; col[3] = row + nd[0];
      col[4] = row;

      c = -2.0/(dx[0][ix] + dx[0][ix+1]);
      val[0] = c/dx[0][ix]; val[1] = c/dx[0][ix+1];

      c = -2.0/(dx[1][iy] + dx[1][iy+1]);
      val[2] = c/dx[1][iy]; val[3] = c/dx[1][iy+1];

      val[4] = -(val[0]+val[1]+val[2]+val[3]);

      for (auto it = bound.cbegin(); it != bound.cend(); ++it) {  // *it = XLO, XHI, YLO or YHI
        int sgn = 1 - 2*((*it)%2);
        Real a;

        switch (mesh.fbc_type(*it)) {
          case Mesh::FBCType::dirichlet:
            break;
          case Mesh::FBCType::neumann :
            col[*it] = NOT_USE;
            val[*it + sgn] += val[*it];
            a = sgn*val[*it]*dx_bc[(*it)/2][(*it)%2]*mesh.fbc_value(*it);
            bc_id_val[1].push_back(std::make_pair(row, a));
            break;
          case Mesh::FBCType::symmetric :
            col[*it] = NOT_USE;
            val[*it + sgn] += val[*it];
            break;
          case Mesh::FBCType::periodic : 
            col[*it] += sgn*nid_pshift[(*it)/2];
        }   // end switch 
      }

      for (int k = 0; k < nelem; ++k) {
        if (col[k] != NOT_USE) 
          Atp.push_back(Tp(row, col[k], val[k]));
      }
    }   // end for (Index ix = 0; ix < nd[0]; ++ix)
  }   // end for (Index iy = 0; iy < nd[1]; ++iy)

  A.setFromTriplets(Atp.begin(), Atp.end());
}

/* ------------------------------------------------------- */

void Poisson::set_matrix_axi()
{ }

/* ------------------------------------------------------- */

void Poisson::check_dirichlet_nodes_2d(std::vector<bool>& is_fixed)
{
  Index nx = mesh.num_nodes(0), ny = mesh.num_nodes(1);

  /* ylo & yhi boundary */
  for (Index k = 0; k < 2; ++k) {
    Index bid = 2 + k;
    Index iy = k*(ny-1);

    if (Mesh::FBCType::dirichlet == mesh.fbc_type(bid)) {
      for (Index ix = 0; ix < nx; ++ix) {
        Index nid = iy*nx + ix;
          
        is_fixed[nid] = true;
        bc_id_val[0].push_back(std::make_pair(nid, mesh.fbc_value(bid)));
      }
    }
  }

  /* xlo & xhi boundary */
  for (Index k = 0; k < 2; ++k) {
    Index bid = k;
    Index ix = k*(nx-1);

    if (Mesh::FBCType::dirichlet == mesh.fbc_type(bid)) {
      for (Index iy = 0; iy < ny; ++iy) {
        Index nid = iy*nx + ix;
          
        is_fixed[nid] = true;
        bc_id_val[0].push_back(std::make_pair(nid, mesh.fbc_value(bid)));
      }
    }
  }
}

/* ------------------------------------------------------- */

void Poisson::check_objects_2d(std::vector<bool>& is_fixed)
{
  std::cout << "Checking conductors in simulation domain: " 
    << mesh.num_conductors() << " conductors in simulation domain.\n";
  if (mesh.num_conductors() == 0) return;
  
  Index nx = mesh.num_nodes(0), ny = mesh.num_nodes(1), nid;
  Real phi;
  bool rf;

  for (Index iy = 0; iy < ny; iy++) {
    for (Index ix = 0; ix < nx; ix++) {
      if (mesh.is_fixed_potential(ix, iy)) {
        nid = iy*nx + ix;
        is_fixed[nid] = true;
        const Conductor* conductor = mesh.get_conductor(ix, iy);
        phi = conductor->get_potential();
        rf = conductor->is_rf();
        bc_id_val[0].push_back(std::make_pair(nid, phi));
        if (rf) is_rf[nid] = true;
      }
    }
  }   // for (Index iy = 0; iy < ny; ++iy)
}

/* ------------------------------------------------------- */

void Poisson::check_boundary(const Index ix,
                             const Index iy,
                             std::vector<int>& bound)
{
  bound.clear();

  if (0 == ix) bound.push_back(XLO);
  else if (mesh.num_cells(0) == ix) bound.push_back(XHI);

  if (0 == iy) bound.push_back(YLO);
  else if (mesh.num_cells(1) == iy) bound.push_back(YHI);

  return;
}

/* ------------------------------------------------------- */

void Poisson::set_rhs_lin(const ScalarField& r)
{
  // initialize rhs with input
  r.copy_to(b);

  // set rf source BC
  Time += timestep;
  Real omegaT = ESPIC::PI2 * 0.0076 * Time;

  // Dirichlet BC
  for (auto it = bc_id_val[0].cbegin(); it != bc_id_val[0].cend(); ++it) {
    if (is_rf[it->first]) 
      b(it->first) = (it->second) * cos(omegaT);

    else 
      b(it->first) = it->second;
  }

  // Neumann BC
  for (auto it = bc_id_val[1].cbegin(); it != bc_id_val[1].cend(); ++it) {
    b(it->first) += it->second;
  }

  for (Index i = 0; i < mesh.num_nodes(); i++) 
      std::cout << b(i) << " ";
  std::cout << std::endl;
}

/* ------------------------------------------------------- */

void Poisson::init_lu()
{
#ifdef USE_MKL
  ptrLUSolver = new Eigen::PardisoLU< SpMat > ();
#else
  ptrLUSolver = new Eigen::SparseLU< SpMat > ();
#endif

  if (linear) {
    ptrLUSolver->compute(A);

    if (ptrLUSolver->info() != Eigen::Success)
      espic_error("Failed in LU decomposition");

    ptrsolve = &Poisson::lu_solve_lin;      // specify solve function

    std::cout << "LU linear solver set.\n";
  }
  else {
    espic_error("Non-linear solver not setup yet.");   // specify solve function
  }

}

/*------------------------------------------------------*/

void Poisson::lu_solve_lin(const ScalarField* r, ScalarField* x)
{
  /* ----------------------------------------------- */
  /* Solve Ax = b for x with LU decomposition method */
  /* ----------------------------------------------- */

  set_rhs_lin(*r);   // copy r to b and set boundary condition

  *x = ptrLUSolver->solve(b);

  result.converged = (ptrLUSolver->info() == Eigen::Success);

  return;
}

/*------------------------------------------------------*/
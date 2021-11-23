#ifndef _SOLVERDEF_H
#define _SOLVERDEF_H

#include <string>

class SolverDef {
  public:
    /* Constructors */
    SolverDef () :
      solver_type("lu"),
      solver_maxit(1000),
      newton_maxit(100),
      solver_tolerance(1e-4),
      newton_tolerance(1e-6) {}

    std::string solver_type;
    int solver_maxit;
    int newton_maxit;
    Real solver_tolerance;
    Real newton_tolerance;
};

#endif

#ifndef ESPIC_TYPE_H
#define ESPIC_TYPE_H

#include <array>
#include <inttypes.h>

#ifdef USE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include "ThirdParty/Eigen3.3.7/Eigen/Sparse"

typedef double Real;
//typedef float Real;

typedef int Smallint;
typedef int Bigint;


//#ifdef BIGNUMBER
//typedef int64_t bigint;
//#endif


typedef Eigen::SparseMatrix<Real, Eigen::RowMajor> SpMatCSR;
typedef Eigen::SparseMatrix<Real, Eigen::ColMajor> SpMatCSC;
typedef Eigen::Triplet<Real> Tp;
typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<Smallint, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXi;
typedef std::array<Real, 3> Vector3;

// -----------------------------------------------
// *** do NOT typedef unsigned for Index ***
// -----------------------------------------------
//typedef typename Vector::Index Index;   // int64_t
typedef int32_t Index;

// class Background {
// public:
//   Background(std::string& name, Real m, Real q, Real n, Real T)
//   : name(name), mass(m),
//   charge(q), ndens(n), temp(T), vth(sqrt(2.*temp/mass))
//   {}

//   const std::string& name;
//   const Real mass;
//   const Real charge;
//   const Real ndens;
//   const Real temp;
//   const Real vth;
// };

// static Background* background;

#endif

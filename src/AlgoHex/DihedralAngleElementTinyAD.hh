/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once


//== INCLUDES =================================================================

#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/StdVector>
#include<cmath>

#include <CoMISo/NSolver/FiniteElementTinyAD.hh>
#include <CoMISo/NSolver/FiniteElementHessianProjection.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================


// penalize deviation from given dihedral angle alpha at edge e0 of a tetrahedron with edge vectors E=(e0|e1|e2)
// f(J) = w*(exp(s*(n0.dot(n1) - cos(alpha))^2)-1.0)
// with n0 = normalize(J*e1 x J*e0) and n1 = normalize(J*e2 x J*e0)
// matrices J and E are col-major

class DihedralAngleFramesElement_F
{
public:

  // define dimensions
  const static int NV = 9; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 12; // [E(3x3),w,s,alpha] = shape matrix E = _c[0..8], weight w = _c[9], exponential constant s=_c[10], cos (target angle alpha)=_c[11]


  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using MatT = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
    using VecT = Eigen::Matrix<ScalarT, 3, 1>;

    // get matrix representation of J and E
    Eigen::Map<MatT> J((ScalarT *) _x.data());
    Eigen::Map<MatD> E((double *) _c.data());

    // get references to constants
    const double &w = _c[9];
    const double &s = _c[10];
    const double &cos_alpha = _c[11];

    // deform tetrahedron (edge vectors)
    MatT V = J * E;

    // calculate normal vectors
    VecT n0 = V.col(1).cross(V.col(0));
    n0 /= n0.norm();
    VecT n1 = V.col(2).cross(V.col(0));
    n1 /= n1.norm();
    // dot product
    ScalarT dp = n0.dot(n1);

    if (!isfinite(dp))
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    else
      return w * (exp(s * pow(dp - cos_alpha, 2)) - 1.0);
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }
};

// create types for optimization with derivatives and projected hessian
using DihedralAngleFramesElement_TAD = COMISO::FiniteElementTinyAD<DihedralAngleFramesElement_F>;
using DihedralAngleFramesElement_TAD_PH = COMISO::FiniteElementHessianProjection<DihedralAngleFramesElement_TAD>;


class DihedralAngleElement_F
{
public:

  // define dimensions
  const static int NV = 12; // P=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T| [x9,x10,x11]^T]] = 3x4 matrix encoding points of tetrahedron
  const static int NC = 4; // [w,s,cos_alpha,sin_alpha] = weight w = _c[0], exponential constant s=_c[1], cos (target angle alpha)=_c[2], sin (target angle alpha)=_c[3]


  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using MatT = Eigen::Matrix<ScalarT, 3, 4, Eigen::ColMajor>;
    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
    using VecT = Eigen::Matrix<ScalarT, 3, 1>;

    // get matrix representation of J and E
    Eigen::Map<MatT> P((ScalarT *) _x.data());

    // get references to constants
    const double &w = _c[0];
    const double &s = _c[1];
    const double &cos_alpha = _c[2];
    const double &sin_alpha = _c[3];

    VecT e0 = P.col(1) - P.col(0);
    VecT e1 = P.col(2) - P.col(0);
    VecT e2 = P.col(3) - P.col(0);

    // calculate normal vectors
    VecT n0 = e0.cross(e1);
    n0 /= n0.norm();
    VecT n1 = e0.cross(e2);
    n1 /= n1.norm();
    // dot product
    ScalarT dp = n0.dot(n1);

    // cross product
    ScalarT cp = (n0.cross(n1)).dot(e0) / e0.norm();

//    std::cerr << "dp-cos_alpha = " << dp-cos_alpha << std::endl;
//    std::cerr << "cp-sin_alpha = " << cp-sin_alpha << std::endl;

//    if(!isfinite(dp))
//      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    if (!isfinite(dp) || !isfinite(cp))
      return w * (exp(2.0 * s) - 1.0);
    else
      return w * (exp(s * (pow(dp - cos_alpha, 2) + pow(cp - sin_alpha, 2))) - 1.0);
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }
};

// create types for optimization with derivatives and projected hessian
using DihedralAngleElement_TAD = COMISO::FiniteElementTinyAD<DihedralAngleElement_F>;
using DihedralAngleElement_TAD_PH = COMISO::FiniteElementHessianProjection<DihedralAngleElement_TAD>;


//=============================================================================
} // namespace AlgoHex
//=============================================================================


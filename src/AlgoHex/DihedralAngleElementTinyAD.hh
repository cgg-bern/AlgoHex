/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once


//== INCLUDES =================================================================

#include <Eigen/Core>

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
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);
};

extern template double DihedralAngleFramesElement_F::eval_f<double>(const Eigen::Matrix<double, DihedralAngleFramesElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, false> DihedralAngleFramesElement_F::eval_f<TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, false>, DihedralAngleFramesElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, true> DihedralAngleFramesElement_F::eval_f<TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, true>, DihedralAngleFramesElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NC, 1> &_c) const;

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
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);
};

extern template double DihedralAngleElement_F::eval_f<double>(const Eigen::Matrix<double, DihedralAngleElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<DihedralAngleElement_F::NV, double, false> DihedralAngleElement_F::eval_f<TinyAD::Scalar<DihedralAngleElement_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngleElement_F::NV, double, false>, DihedralAngleElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<DihedralAngleElement_F::NV, double, true> DihedralAngleElement_F::eval_f<TinyAD::Scalar<DihedralAngleElement_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngleElement_F::NV, double, true>, DihedralAngleElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleElement_F::NC, 1> &_c) const;

using DihedralAngleElement_TAD = COMISO::FiniteElementTinyAD<DihedralAngleElement_F>;
using DihedralAngleElement_TAD_PH = COMISO::FiniteElementHessianProjection<DihedralAngleElement_TAD>;

//=============================================================================
} // namespace AlgoHex
//=============================================================================

namespace COMISO
{
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::DihedralAngleFramesElement_F>>;
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::DihedralAngleElement_F>>;
} // namespace COMISO

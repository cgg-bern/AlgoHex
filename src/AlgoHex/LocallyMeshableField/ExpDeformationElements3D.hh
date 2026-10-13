/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include <Eigen/Core>

#include <CoMISo/NSolver/FiniteElementTinyAD.hh>
#include <CoMISo/NSolver/FiniteElementHessianProjection.hh>


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//f(J(x0, x1, x2)) = exp(w*(tr(J^TA^TAJ)*(det(AJ)^(-2/3))
class ExpDeformationElement
{
public:

  // define dimensions
  const static int NV = 3; // [ax,ay,az] = 3D points A, B, C, D of a tetrahedron
  const static int NC = 19; // [weight_conformal, bx,by,bz,cx,cy,cz,dx,dy,dz,]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

  // calculate constants from reference tetrahedron
  static void compute_constants(const double _w_conformal,
                                const Eigen::Vector3d &_p0, const Eigen::Vector3d &_p1,
                                const Eigen::Vector3d &_p2, const Eigen::Vector3d &_p3,
                                VecC &_c);

  static double
  eval_fx(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c);
};

extern template double ExpDeformationElement::eval_f<double>(const Eigen::Matrix<double, ExpDeformationElement::NV, 1> &_x, const Eigen::Matrix<double, ExpDeformationElement::NC, 1> &_c) const;
extern template TinyAD::Scalar<ExpDeformationElement::NV, double, false> ExpDeformationElement::eval_f<TinyAD::Scalar<ExpDeformationElement::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<ExpDeformationElement::NV, double, false>, ExpDeformationElement::NV, 1> &_x, const Eigen::Matrix<double, ExpDeformationElement::NC, 1> &_c) const;
extern template TinyAD::Scalar<ExpDeformationElement::NV, double, true> ExpDeformationElement::eval_f<TinyAD::Scalar<ExpDeformationElement::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<ExpDeformationElement::NV, double, true>, ExpDeformationElement::NV, 1> &_x, const Eigen::Matrix<double, ExpDeformationElement::NC, 1> &_c) const;

using EDE_AD = COMISO::FiniteElementTinyAD<ExpDeformationElement>;
using EDE_PH_AD = COMISO::FiniteElementHessianProjection<EDE_AD>;
}
//== EXPLICIT INSTANTIATION ===================================================
// The wrapper chains below are explicitly instantiated once in
// ExpDeformationElements3D.cc; they are only declared extern here so that
// consumer translation units do not re-instantiate them.
namespace COMISO
{
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::ExpDeformationElement>>;
} // namespace COMISO

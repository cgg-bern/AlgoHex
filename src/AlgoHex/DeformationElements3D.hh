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


class SymmetricDirichletElement3D_F
{
public:

  // define dimensions
  const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
  const static int NC = 12; // [weight, size_scale, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

  // NC = 12; // [weight, size_scale, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
  // calculate constants from reference tetrahedron
  static void compute_constants(const double _w, const double _size_scale, const double _epsilon,
                                const Eigen::Vector3d &_p0, const Eigen::Vector3d &_p1, const Eigen::Vector3d &_p2,
                                const Eigen::Vector3d &_p3,
                                VecC &_c);

  static double
  determine_suitable_epsilon(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c,
                             const double _target_penalty);

};

extern template double SymmetricDirichletElement3D_F::eval_f<double>(const Eigen::Matrix<double, SymmetricDirichletElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<SymmetricDirichletElement3D_F::NV, double, false> SymmetricDirichletElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletElement3D_F::NV, double, false>, SymmetricDirichletElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<SymmetricDirichletElement3D_F::NV, double, true> SymmetricDirichletElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletElement3D_F::NV, double, true>, SymmetricDirichletElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletElement3D_F::NC, 1> &_c) const;

using SDE3D_AD = COMISO::FiniteElementTinyAD<SymmetricDirichletElement3D_F>;
using SDE3D_PH_AD = COMISO::FiniteElementHessianProjection<SDE3D_AD>;


class FoldoverFree3DV_F
{
public:

  // define dimensions
  const static int NV = 12; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz] = 3D points A, B, C, D of a tetrahedron
  const static int NC = 13; // [weight_conformal, weight_det, target_det, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

  //  const static int NC = 13; // [weight_conformal, weight_det, target_det, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
  // calculate constants from reference tetrahedron
  static void
  compute_constants(const double _w_conformal, const double _w_det, const double _target_det, const double _epsilon,
                    const Eigen::Vector3d &_p0, const Eigen::Vector3d &_p1, const Eigen::Vector3d &_p2,
                    const Eigen::Vector3d &_p3,
                    VecC &_c);

  static double
  determine_suitable_epsilon(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c);

};

extern template double FoldoverFree3DV_F::eval_f<double>(const Eigen::Matrix<double, FoldoverFree3DV_F::NV, 1> &_x, const Eigen::Matrix<double, FoldoverFree3DV_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<FoldoverFree3DV_F::NV, double, false> FoldoverFree3DV_F::eval_f<TinyAD::Scalar<FoldoverFree3DV_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<FoldoverFree3DV_F::NV, double, false>, FoldoverFree3DV_F::NV, 1> &_x, const Eigen::Matrix<double, FoldoverFree3DV_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<FoldoverFree3DV_F::NV, double, true> FoldoverFree3DV_F::eval_f<TinyAD::Scalar<FoldoverFree3DV_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<FoldoverFree3DV_F::NV, double, true>, FoldoverFree3DV_F::NV, 1> &_x, const Eigen::Matrix<double, FoldoverFree3DV_F::NC, 1> &_c) const;

using FF3DV_AD = COMISO::FiniteElementTinyAD<FoldoverFree3DV_F>;
using FF3DV_PH_AD = COMISO::FiniteElementHessianProjection<FF3DV_AD>;


class FoldoverFreeDual3DV_F
{
public:

  // define dimensions
  const static int NV = 24; // [ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz,  ex,ey,ez, fx,fy,fz, gx,gy,gz, hx,hy,hz] = 3D points A, B, C, D of image tetrahedron and 3D points E, F, G, H of domain tetrahedron
  const static int NC = 4; // [weight_conformal, weight_det, target_det, epsilon]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

};

extern template double FoldoverFreeDual3DV_F::eval_f<double>(const Eigen::Matrix<double, FoldoverFreeDual3DV_F::NV, 1> &_x, const Eigen::Matrix<double, FoldoverFreeDual3DV_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<FoldoverFreeDual3DV_F::NV, double, false> FoldoverFreeDual3DV_F::eval_f<TinyAD::Scalar<FoldoverFreeDual3DV_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<FoldoverFreeDual3DV_F::NV, double, false>, FoldoverFreeDual3DV_F::NV, 1> &_x, const Eigen::Matrix<double, FoldoverFreeDual3DV_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<FoldoverFreeDual3DV_F::NV, double, true> FoldoverFreeDual3DV_F::eval_f<TinyAD::Scalar<FoldoverFreeDual3DV_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<FoldoverFreeDual3DV_F::NV, double, true>, FoldoverFreeDual3DV_F::NV, 1> &_x, const Eigen::Matrix<double, FoldoverFreeDual3DV_F::NC, 1> &_c) const;

using FFD3DV_AD = COMISO::FiniteElementTinyAD<FoldoverFreeDual3DV_F>;
using FFD3DV_PH_AD = COMISO::FiniteElementHessianProjection<FFD3DV_AD>;

class AlignmentElement_F
{
public:

  // define dimensions
  const static int NV = 2;
  const static int NC = 1;

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

};

extern template double AlignmentElement_F::eval_f<double>(const Eigen::Matrix<double, AlignmentElement_F::NV, 1> &_x, const Eigen::Matrix<double, AlignmentElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<AlignmentElement_F::NV, double, false> AlignmentElement_F::eval_f<TinyAD::Scalar<AlignmentElement_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<AlignmentElement_F::NV, double, false>, AlignmentElement_F::NV, 1> &_x, const Eigen::Matrix<double, AlignmentElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<AlignmentElement_F::NV, double, true> AlignmentElement_F::eval_f<TinyAD::Scalar<AlignmentElement_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<AlignmentElement_F::NV, double, true>, AlignmentElement_F::NV, 1> &_x, const Eigen::Matrix<double, AlignmentElement_F::NC, 1> &_c) const;

using AE_AD = COMISO::FiniteElementTinyAD<AlignmentElement_F>;

class FittingElement_F
{
public:

  // define dimensions
  const static int NV = 1;
  const static int NC = 2;

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

};

extern template double FittingElement_F::eval_f<double>(const Eigen::Matrix<double, FittingElement_F::NV, 1> &_x, const Eigen::Matrix<double, FittingElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<FittingElement_F::NV, double, false> FittingElement_F::eval_f<TinyAD::Scalar<FittingElement_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<FittingElement_F::NV, double, false>, FittingElement_F::NV, 1> &_x, const Eigen::Matrix<double, FittingElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<FittingElement_F::NV, double, true> FittingElement_F::eval_f<TinyAD::Scalar<FittingElement_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<FittingElement_F::NV, double, true>, FittingElement_F::NV, 1> &_x, const Eigen::Matrix<double, FittingElement_F::NC, 1> &_c) const;

using FE_AD = COMISO::FiniteElementTinyAD<FittingElement_F>;

class DihedralAngle_F
{
public:

  // define dimensions
  const static int NV = 12;
  const static int NC = 3;

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

};

extern template double DihedralAngle_F::eval_f<double>(const Eigen::Matrix<double, DihedralAngle_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngle_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<DihedralAngle_F::NV, double, false> DihedralAngle_F::eval_f<TinyAD::Scalar<DihedralAngle_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngle_F::NV, double, false>, DihedralAngle_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngle_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<DihedralAngle_F::NV, double, true> DihedralAngle_F::eval_f<TinyAD::Scalar<DihedralAngle_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngle_F::NV, double, true>, DihedralAngle_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngle_F::NC, 1> &_c) const;

using DAE_AD = COMISO::FiniteElementTinyAD<DihedralAngle_F>;
using DAE_PH_AD = COMISO::FiniteElementHessianProjection<DAE_AD>;

class BallBarrier3D_F
{
public:

  // define dimensions
  const static int NV = 3; // [x,y,z] = 3D point location
  const static int NC = 5; // [weight, radius, cx, cy, cz]  with center c=[cx,cy,cz]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

};

extern template double BallBarrier3D_F::eval_f<double>(const Eigen::Matrix<double, BallBarrier3D_F::NV, 1> &_x, const Eigen::Matrix<double, BallBarrier3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<BallBarrier3D_F::NV, double, false> BallBarrier3D_F::eval_f<TinyAD::Scalar<BallBarrier3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<BallBarrier3D_F::NV, double, false>, BallBarrier3D_F::NV, 1> &_x, const Eigen::Matrix<double, BallBarrier3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<BallBarrier3D_F::NV, double, true> BallBarrier3D_F::eval_f<TinyAD::Scalar<BallBarrier3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<BallBarrier3D_F::NV, double, true>, BallBarrier3D_F::NV, 1> &_x, const Eigen::Matrix<double, BallBarrier3D_F::NC, 1> &_c) const;

using BB3D_AD = COMISO::FiniteElementTinyAD<BallBarrier3D_F>;
using BB3D_AD_PH = COMISO::FiniteElementHessianProjection<BB3D_AD>;

class JacobianSmoothnessElement3D_F
{
public:

  // define dimensions
  const static int NV = 24; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz, ...] = 3D points A, B, C, D of two tetrahedra
  const static int NC = 29; // [weight, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z, h1x,... h3z, t1x, ..., t3z, sizing_scale]

  using VecV = Eigen::Matrix<double,NV,1>;
  using VecC = Eigen::Matrix<double,NC,1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

  //  const static int NC = 29; // [weight, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z, h1x,... h3z, t1x, ..., t3z, sizing_scale]
  // calculate constants from reference tetrahedron
  static void compute_constants( const double _w,
                                 const double _sizing_scale,
                                 const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                 const Eigen::Vector3d& _q0, const Eigen::Vector3d& _q1, const Eigen::Vector3d& _q2, const Eigen::Vector3d& _q3,
                                 const Eigen::Matrix3d& _T01,
                                 VecC& _c);

};

extern template double JacobianSmoothnessElement3D_F::eval_f<double>(const Eigen::Matrix<double, JacobianSmoothnessElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, JacobianSmoothnessElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<JacobianSmoothnessElement3D_F::NV, double, false> JacobianSmoothnessElement3D_F::eval_f<TinyAD::Scalar<JacobianSmoothnessElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<JacobianSmoothnessElement3D_F::NV, double, false>, JacobianSmoothnessElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, JacobianSmoothnessElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<JacobianSmoothnessElement3D_F::NV, double, true> JacobianSmoothnessElement3D_F::eval_f<TinyAD::Scalar<JacobianSmoothnessElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<JacobianSmoothnessElement3D_F::NV, double, true>, JacobianSmoothnessElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, JacobianSmoothnessElement3D_F::NC, 1> &_c) const;

using JSE3D_AD = COMISO::FiniteElementTinyAD<JacobianSmoothnessElement3D_F>;

//=============================================================================
} // namespace AlgoHex

namespace COMISO
{
extern template class FiniteElementTinyAD<AlgoHex::SymmetricDirichletElement3D_F>;
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::SymmetricDirichletElement3D_F>>;
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::FoldoverFree3DV_F>>;
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::FoldoverFreeDual3DV_F>>;
extern template class FiniteElementTinyAD<AlgoHex::BallBarrier3D_F>;
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::BallBarrier3D_F>>;
extern template class FiniteElementTinyAD<AlgoHex::JacobianSmoothnessElement3D_F>;
} // namespace COMISO

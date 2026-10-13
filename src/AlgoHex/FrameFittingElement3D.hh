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

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================


/** \class FrameFittingElement3D FrameFittingElement3D.hh
 *
 *  @brief Implementation of a 3d-frame energy element for the construction of integer-grid maps
 *   || Jacobi_f * F - Id||_alpha^2 with f being the mapping, F being a frame and alpha from (0,0.5] an anisotropic norm
 *
 */

// 3D-Frame Fitting Element
// element for 1/2*int_tetrahedron(||nabla f *F - Id||^2)dV = 1/2*V ||(u|v|w)^T G^T F - Id||^2
class FrameFittingElement3D
{
public:

  // define dimensions
  const static int NV = 12;
  const static int NC = 14; // C00, C01, ... , C32, alpha, volume

  typedef Eigen::Matrix<size_t, NV, 1> VecI;
  typedef Eigen::Matrix<double, NV, 1> VecV;
  typedef Eigen::Matrix<double, NC, 1> VecC;
  typedef Eigen::Triplet<double> Triplet;

  double eval_f(const VecV &_x, const VecC &_c) const;

  void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const;

  void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

  // generate required constants for a given input triangle
  static void constants_from_tetrahedron_and_frame(const Eigen::Matrix<double, 3, 4> &_P,
                                                   const Eigen::Matrix3d &_F,
                                                   const double _alpha,
                                                   VecC &_c);

};


class FrameFittingExpElement3D_F
{
public:

  // define dimensions
  const static int NV = 15; // [ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz, su, sv, sw] = 3D points A, B, C, D of a parametric tetrahedron and target sizes su, sv, sw
  const static int NC = 21; // [weight, weight_exp, weight_aniso, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z, f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z,]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
  using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
  using Mat3x4dT = Eigen::Matrix<double, 3, 4, Eigen::ColMajor>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

  // NC = 21; // [weight, weight_exp, weight_aniso, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z, f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z,]
  // calculate constants from reference tetrahedron
  static void compute_constants(const double _w, const double _w_exp, const double _w_aniso,
                                const Eigen::Matrix<double, 3, 4> &_P, const Eigen::Matrix3d &_F,
                                VecC &_c);

  static void estimate_sizing(const Mat3x4dT &_P0, const Mat3x4dT &_P, const Mat3x3dT &_F, Eigen::Map<Vec3dT> &_optimal_sizing);

  static double max_exp_arg(const VecV &_x, const VecC &_c);
};

extern template double FrameFittingExpElement3D_F::eval_f<double>(const Eigen::Matrix<double, FrameFittingExpElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, FrameFittingExpElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<FrameFittingExpElement3D_F::NV, double, false> FrameFittingExpElement3D_F::eval_f<TinyAD::Scalar<FrameFittingExpElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<FrameFittingExpElement3D_F::NV, double, false>, FrameFittingExpElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, FrameFittingExpElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<FrameFittingExpElement3D_F::NV, double, true> FrameFittingExpElement3D_F::eval_f<TinyAD::Scalar<FrameFittingExpElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<FrameFittingExpElement3D_F::NV, double, true>, FrameFittingExpElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, FrameFittingExpElement3D_F::NC, 1> &_c) const;

// add derivatives -- no hessian projection required since convex!
using FrameFittingExpElement3D_TAD = COMISO::FiniteElementTinyAD<FrameFittingExpElement3D_F>;

//=============================================================================
} // namespace AlgoHex
//=============================================================================

namespace COMISO
{
extern template class FiniteElementTinyAD<AlgoHex::FrameFittingExpElement3D_F>;
} // namespace COMISO

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
#include <Eigen/SparseCore>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================


// element for f(x) = c_(n+1)*(c.dot(x)+c_n)^2
template<int NVT>
class LinearLeastSquaresElementND
{
public:

  // define dimensions
  const static int NV = NVT; // [x0,x1,x2, ...]^T
  const static int NC = NVT + 2; // [c0,c1,c2, ..., cn,cn+1] = linear coeffs and cn+1 is weight

  typedef Eigen::Matrix<size_t, NV, 1> VecI;
  typedef Eigen::Matrix<double, NV, 1> VecV;
  typedef Eigen::Matrix<double, NC, 1> VecC;
  typedef Eigen::Triplet<double> Triplet;

  double eval_f(const VecV &_x, const VecC &_c) const;

  void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const;

  void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c) const;

  // _alpha is the penalty factor for deviation which is not along the axis
  // requirs NV=3
  static void
  get_constants_and_indices_frame_fitting(const Eigen::Matrix3d &_F, const int _cell_idx, const double _w, const double _alpha,
                                          std::vector<VecI> &_idxs, std::vector<VecC> &_constants);
};

template<> void LinearLeastSquaresElementND<3>::get_constants_and_indices_frame_fitting(const Eigen::Matrix3d &_F, const int _cell_idx, const double _w, const double _alpha,
                                            std::vector<VecI> &_idxs, std::vector<VecC> &_constants);

using LinearLeastSquaresElement1D = LinearLeastSquaresElementND<1>;
using LinearLeastSquaresElement2D = LinearLeastSquaresElementND<2>;
using LinearLeastSquaresElement3D = LinearLeastSquaresElementND<3>;
using LinearLeastSquaresElement4D = LinearLeastSquaresElementND<4>;
using LinearLeastSquaresElement5D = LinearLeastSquaresElementND<5>;
using LinearLeastSquaresElement6D = LinearLeastSquaresElementND<6>;
using LinearLeastSquaresElement9D = LinearLeastSquaresElementND<9>;

extern template class LinearLeastSquaresElementND<1>;
extern template class LinearLeastSquaresElementND<2>;
extern template class LinearLeastSquaresElementND<3>;
extern template class LinearLeastSquaresElementND<4>;
extern template class LinearLeastSquaresElementND<5>;
extern template class LinearLeastSquaresElementND<6>;
extern template class LinearLeastSquaresElementND<9>;

//=============================================================================
} // namespace AlgoHex
//=============================================================================

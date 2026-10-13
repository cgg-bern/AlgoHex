/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include "LinearLeastSquaresElements.hh"

#include <cfloat>
#include <Eigen/Geometry>

namespace AlgoHex
{

template<int NVT>
double LinearLeastSquaresElementND<NVT>::eval_f(const VecV &_x, const VecC &_c) const
{

  double d = _x.dot(_c.head(NV)) + _c[NV];

  return _c[NV + 1] * d * d;

}

template<int NVT>
void LinearLeastSquaresElementND<NVT>::eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
{

  double d = _x.dot(_c.head(NV)) + _c[NV];

  for (unsigned int i = 0; i < NV; ++i)
    _g[i] = _c[NV + 1] * 2.0 * d * _c[i];

}

template<int NVT>
void LinearLeastSquaresElementND<NVT>::eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const
{

  _triplets.clear();

  for (int i = 0; i < NV; ++i)
    for (int j = 0; j < NV; ++j)
      _triplets.push_back(Triplet(i, j, _c[NV + 1] * 2.0 * _c[i] * _c[j]));

}

template<int NVT>
double LinearLeastSquaresElementND<NVT>::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c) const
{

  return DBL_MAX;

}

template<int NVT>
void LinearLeastSquaresElementND<NVT>::get_constants_and_indices_frame_fitting(const Eigen::Matrix3d &_F, const int _cell_idx, const double _w, const double _alpha,
                                            std::vector<VecI> &_idxs, std::vector<VecC> &_constants)
{

}

template<>
void LinearLeastSquaresElementND<3>::get_constants_and_indices_frame_fitting(const Eigen::Matrix3d &_F, const int _cell_idx, const double _w, const double _alpha,
                                            std::vector<VecI> &_idxs, std::vector<VecC> &_constants)
{
  // assume _F has positive determinant (non-degenerate + right-handed)

  _idxs.clear();
  _idxs.reserve(9);
  _constants.clear();
  _constants.reserve(9);

  int bi = 9 * _cell_idx;

  for (int i = 0; i < 3; ++i)
  {
    int bi2 = bi + 3 * i;
    _idxs.push_back(VecI(bi2, bi2 + 1, bi2 + 2));
    _idxs.push_back(VecI(bi2, bi2 + 1, bi2 + 2));
    _idxs.push_back(VecI(bi2, bi2 + 1, bi2 + 2));

    // create local basis for anisotropic norm
    Eigen::Matrix3d B;
    B.col(0) = _F.col(i);
    B.col(1) = _F.col((i + 1) % 3);
    B.col(2) = _F.col((i + 2) % 3);

    // first axis is axis of frame
    B.col(0) /= (B.col(0)).norm();
    // second axis via cross-product
    B.col(2) = B.col(0).cross(B.col(1));
    B.col(2) /= (B.col(2)).norm();
    // third axis via cross-product
    B.col(1) = B.col(2).cross(B.col(0));
    B.col(1) /= (B.col(1)).norm();

    // get current axis in local coordinates
    Eigen::Vector3d fb = B.transpose() * _F.col(i);

    for (int j = 0; j < 3; ++j)
    {
      double w = _w;
      if (i != j) w *= _alpha;
      VecC c;
      c << B(0, j), B(1, j), B(2, j), -fb[j], w;
      _constants.push_back(c);
    }
  }

}

template class LinearLeastSquaresElementND<1>;
template class LinearLeastSquaresElementND<2>;
template class LinearLeastSquaresElementND<3>;
template class LinearLeastSquaresElementND<4>;
template class LinearLeastSquaresElementND<5>;
template class LinearLeastSquaresElementND<6>;
template class LinearLeastSquaresElementND<9>;

} // namespace AlgoHex

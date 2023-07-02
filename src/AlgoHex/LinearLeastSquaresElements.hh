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

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================


/** \class LinearElement3D LinearElement3D.hh
 *
 *  @brief Implementation of a minus log det element for optimization. f(x0,x1,x2) = w*(c0*x0 + c1*x1 + c2*x2 + c3)^2
 *
 */

//  // element for f(x0,x1,x2) = c4*(c0*x0 + c1*x1 + c2*x2 + c3)^2
//  class LinearLeastSquaresElement3D
//  {
//  public:
//
//    // define dimensions
//    const static int NV = 3; // [x0,x1,x2]^T
//    const static int NC = 5; // [c0,c1,c2,c3,c4] = linear coeffs and c4 is weight
//
//    typedef Eigen::Matrix<size_t, NV, 1> VecI;
//    typedef Eigen::Matrix<double, NV, 1> VecV;
//    typedef Eigen::Matrix<double, NC, 1> VecC;
//    typedef Eigen::Triplet<double> Triplet;
//
//    inline double eval_f(const VecV &_x, const VecC &_c) const
//    {
//      double d = _c[0]*_x[0] + _c[1]*_x[1] + _c[2]*_x[2] + _c[3];
//
//      return _c[4] * d * d;
//    }
//
//    inline void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
//    {
//      double d = _c[0]*_x[0] + _c[1]*_x[1] + _c[2]*_x[2] + _c[3];
//
//      _g[0] = 2.0 * d * _c[0];
//      _g[1] = 2.0 * d * _c[1];
//      _g[2] = 2.0 * d * _c[2];
//    }
//
//    inline void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const
//    {
//      _triplets.clear();
//
//      for(int i=0; i<NV; ++i)
//        for(int j=0; j<NV; ++j)
//          _triplets.push_back( Triplet(i,j, 2.0*_c[i]*_c[j]));
//    }
//
//    inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c) const
//    {
//      return DBL_MAX;
//    }
//
//
//    // _alpha is the penalty factor for deviation which is not along the axis
//    static void get_constants_and_indices_frame_fitting(const Mat3d& _F, const int _cell_idx, const double _w, const double _alpha, std::vector<VecI> _idxs, std::vector<VecC> &_constants)
//    {
//      // assume _F has positive determinant (non-degenerate + right-handed)
//
//      _idxs.clear();
//      _idxs.reserve(9);
//      _constants.clear();
//      _constants.reserve(9);
//
//      for(int i=0; i<3; ++i)
//      {
//        _idxs.push_back(VecI(3*i, 3*i+1, 3*i+2));
//        _idxs.push_back(VecI(3*i, 3*i+1, 3*i+2));
//        _idxs.push_back(VecI(3*i, 3*i+1, 3*i+2));
//
//        // create local basis for anisotropic norm
//        Mat3d B;
//        B.col(0) = _F.col(i);
//        B.col(1) = _F.col((i+1)%3);
//        B.col(2) = _F.col((i+2)%3);
//
//        // first axis is axis of frame
//        B.col(0) /= (B.col(0)).norm();
//        // second axis via cross-product
//        B.col(2)  = B.col(0).cross(B.col(1));
//        B.col(2) /= (B.col(2)).norm();
//        // third axis via cross-product
//        B.col(1)  = B.col(2).cross(B.col(0));
//        B.col(1) /= (B.col(1)).norm();
//
//        // get current axis in local coordinates
//        Vec3d fb = B.transpose()*_F.col(i);
//
//        for(int j=0; j<3; ++j)
//        {
//          VecC c;
//          c << B(0,j), B(1,j), B(2,j), -fb[j], _w*( i==j ? 1.0 : _alpha);
//          _constants.push_back(c);
//        }
//      }
//    }
//  };


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

  inline double eval_f(const VecV &_x, const VecC &_c) const
  {
    double d = _x.dot(_c.head(NV)) + _c[NV];

    return _c[NV + 1] * d * d;
  }

  inline void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
  {
    double d = _x.dot(_c.head(NV)) + _c[NV];

    for (unsigned int i = 0; i < NV; ++i)
      _g[i] = _c[NV + 1] * 2.0 * d * _c[i];
  }

  inline void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const
  {
    _triplets.clear();

    for (int i = 0; i < NV; ++i)
      for (int j = 0; j < NV; ++j)
        _triplets.push_back(Triplet(i, j, _c[NV + 1] * 2.0 * _c[i] * _c[j]));
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c) const
  {
    return DBL_MAX;
  }


  // _alpha is the penalty factor for deviation which is not along the axis
  // requirs NV=3
  static void
  get_constants_and_indices_frame_fitting(const Mat3d &_F, const int _cell_idx, const double _w, const double _alpha,
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
      Mat3d B;
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
      Vec3d fb = B.transpose() * _F.col(i);

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
};

using LinearLeastSquaresElement1D = LinearLeastSquaresElementND<1>;
using LinearLeastSquaresElement2D = LinearLeastSquaresElementND<2>;
using LinearLeastSquaresElement3D = LinearLeastSquaresElementND<3>;
using LinearLeastSquaresElement4D = LinearLeastSquaresElementND<4>;
using LinearLeastSquaresElement5D = LinearLeastSquaresElementND<5>;
using LinearLeastSquaresElement6D = LinearLeastSquaresElementND<6>;
using LinearLeastSquaresElement9D = LinearLeastSquaresElementND<9>;

//=============================================================================
} // namespace AlgoHex
//=============================================================================


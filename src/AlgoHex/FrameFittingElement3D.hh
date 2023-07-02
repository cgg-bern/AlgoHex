/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once


//== INCLUDES =================================================================

#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/StdVector>
#include<cmath>

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

  inline double eval_f(const VecV &_x, const VecC &_c) const
  {
    //     std::cerr << "eval f" << std::endl;
    double f(0);

    // map first 12 coefficients into 3x4 matrix
    Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor> > C((double *) _c.data());

    Eigen::Vector3d v = C * _x.segment(0, 4) - Eigen::Vector3d(_c[12], 0, 0);
    f += v.dot(v);
    v = C * _x.segment(4, 4) - Eigen::Vector3d(0, _c[12], 0);
    f += v.dot(v);
    v = C * _x.segment(8, 4) - Eigen::Vector3d(0, 0, _c[12]);
    f += v.dot(v);

    return 0.5 * _c[13] * f;
  }

  inline void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
  {
//      std::cerr << "eval grad" << std::endl;
    // map first 12 coefficients into 3x4 matrix
    Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor> > C((double *) _c.data());

    Eigen::Matrix<double, 4, 4> CtC = C.transpose() * C;

    _g.segment(0, 4) = _c[13] * (CtC * _x.segment(0, 4) - C.transpose() * Eigen::Vector3d(_c[12], 0, 0));
    _g.segment(4, 4) = _c[13] * (CtC * _x.segment(4, 4) - C.transpose() * Eigen::Vector3d(0, _c[12], 0));
    _g.segment(8, 4) = _c[13] * (CtC * _x.segment(8, 4) - C.transpose() * Eigen::Vector3d(0, 0, _c[12]));
  }

  inline void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const
  {
    //     std::cerr << "eval hessian" << std::endl;

    // map first 12 coefficients into 3x4 matrix
    Eigen::Map<Eigen::Matrix<double, 3, 4, Eigen::RowMajor> > C((double *) _c.data());
    Eigen::Matrix<double, 4, 4> CtC = C.transpose() * C;
    CtC *= _c[13];

    _triplets.clear();
    for (unsigned int i = 0; i < 4; ++i)
      for (unsigned int j = 0; j < 4; ++j)
      {
        _triplets.push_back(Triplet(i, j, CtC(i, j)));
        _triplets.push_back(Triplet(i + 4, j + 4, CtC(i, j)));
        _triplets.push_back(Triplet(i + 8, j + 8, CtC(i, j)));
      }
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }

  // generate required constants for a given input triangle
  static void constants_from_tetrahedron_and_frame(const Eigen::Matrix<double, 3, 4> &_P,
                                                   const Eigen::Matrix3d &_F,
                                                   const double _alpha,
                                                   VecC &_c)
  {
    // get square-root of alpha
    double alpha_sqrt = std::sqrt(_alpha);

    // volume of tetrahedron
    Eigen::Matrix3d E;
    E.col(0) = _P.col(1) - _P.col(0);
    E.col(1) = _P.col(2) - _P.col(0);
    E.col(2) = _P.col(3) - _P.col(0);
    double V = 1.0 / 6.0 * E.determinant();

    if (V <= 0.0)
    {
      std::cerr << "Warning: degenerate or inverted tet with volume " << V << std::endl;
      V = std::abs(V);
    }

    // Gradient matrix
    Eigen::Matrix<double, 3, 4> G;
    for (unsigned int i = 0; i < 4; ++i)
    {
      // get points opposite to pi
      Eigen::Vector3d q0 = _P.col((i + 1) % 4);
      Eigen::Vector3d q1 = _P.col((i + 2) % 4);
      Eigen::Vector3d q2 = _P.col((i + 3) % 4);

      if (i % 2 == 0) // correct orientation
        std::swap(q1, q2);

      // store opposite normal vector with length 2*area
      G.col(i) = (q1 - q0).cross(q2 - q0);
    }

    // divide by 6*volume to obtain gradient
    G /= 6.0 * V;

    Eigen::Matrix<double, 3, 4> C_alpha = _F.transpose() * G;
    C_alpha.row(0) *= alpha_sqrt;

    if (!std::isfinite(C_alpha.squaredNorm()))
    {
      std::cerr << "ERROR: failed generating valid constants for FrameFittingElement3D" << std::endl;
      C_alpha.setZero();
    }

    // store constants
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 4; ++j)
        _c[j + 4 * i] = C_alpha(i, j);

    _c[12] = alpha_sqrt;
    _c[13] = V;
  }
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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
    using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using Mat3x4T = Eigen::Matrix<ScalarT, 3, 4, Eigen::ColMajor>;

    // get constants
    const double &w = _c[0];
    const double &w_exp = _c[1];
    const double &w_aniso = _c[2];

    // get gradient computation matrix
    Eigen::Map<Mat3x3dT> G((double *) &(_c[3]));

    // get frame matrix
    Eigen::Map<Mat3x3dT> F((double *) &(_c[12]));

    // get points of tetrahedron
    Eigen::Map<Mat3x4T> P((ScalarT *) &(_x[0]));

    // get target size
    Eigen::Map<Vec3T> S((ScalarT *) &(_x[12]));

    // get edge vectors
    Mat3x3T E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    // calc Jacobi Matrix of map
//    Mat3x3T JF = E*G.transpose()*F;
    Mat3x3T JF = 2.0 * E * G.transpose() * F;

    // subtract target size
    JF(0, 0) -= S[0];
    JF(1, 1) -= S[1];
    JF(2, 2) -= S[2];

    ScalarT f = exp(w_exp * pow(JF(0, 0), 2)) + exp(w_exp * pow(JF(1, 1), 2)) + exp(w_exp * pow(JF(2, 2), 2)) +
                2.0 * w_aniso *
                (exp(w_exp * pow(JF(1, 0), 2)) + exp(w_exp * pow(JF(2, 0), 2)) + exp(w_exp * pow(JF(1, 2), 2)));

//    const double eps = 1e-4;
//
//    ScalarT f = exp(w_exp*sqrt(eps+pow(JF(0,0),2)));
//    f        += exp(w_exp*sqrt(eps+pow(JF(1,1),2)));
//    f        += exp(w_exp*sqrt(eps+pow(JF(2,2),2)));
//    f        += exp(w_exp*sqrt(eps+pow(JF(1,0),2)))*2.0;
//    f        += exp(w_exp*sqrt(eps+pow(JF(2,0),2)))*2.0;
//    f        += exp(w_exp*sqrt(eps+pow(JF(1,2),2)))*2.0;

    return w * (f - 9.0);
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }

  // NC = 21; // [weight, weight_exp, weight_aniso, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z, f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z,]
  // calculate constants from reference tetrahedron
  static void compute_constants(const double _w, const double _w_exp, const double _w_aniso,
                                const Eigen::Matrix<double, 3, 4> &_P, const Eigen::Matrix3d &_F,
                                VecC &_c)
  {
    _c[0] = _w;
    _c[1] = _w_exp;
    _c[2] = _w_aniso;

    // compute gradients of pwl basis functions g1, g2, g3

    // get edge vectors
    Eigen::Matrix3d E;
    E.col(0) = _P.col(1) - _P.col(0);
    E.col(1) = _P.col(2) - _P.col(0);
    E.col(2) = _P.col(3) - _P.col(0);

    double d = E.determinant();
    double di = 1.0 / d;
    if (!std::isfinite(di) || d < 0.0)
    {
      std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;
      di = 0.0;
    }

    _c.segment<3>(3) = di * (E.col(1)).cross(E.col(2));
    _c.segment<3>(6) = di * (E.col(2)).cross(E.col(0));
    _c.segment<3>(9) = di * (E.col(0)).cross(E.col(1));

    // store frame
    Eigen::Map<Eigen::Matrix3d> F(&(_c[12]));
    F = _F;

    // -------------------
    // test and debug code
    // -------------------
    if (0)
    {
      Eigen::Matrix3d J;
      J = _c.segment<3>(3) * E.col(0).transpose() +
          _c.segment<3>(6) * E.col(1).transpose() +
          _c.segment<3>(9) * E.col(2).transpose();

      J -= Eigen::Matrix3d::Identity();

      double Jn = J.norm();
      if (Jn > 1e-6)
      {
        std::cerr << "ERROR: Jacobi matrix was not computed correctly!" << std::endl;

        J = _c.segment<3>(3) * E.col(0).transpose() +
            _c.segment<3>(6) * E.col(1).transpose() +
            _c.segment<3>(9) * E.col(2).transpose();
        std::cerr << J << std::endl;
      }
      else
        std::cerr << "Jacobi matrix within tolerance ---> " << Jn << std::endl;
    }
  }

  template<class VecT>
  static void estimate_sizing(const Mat3x4dT &_P0, const Mat3x4dT &_P, const Mat3x3dT &_F, VecT &_optimal_sizing)
  {
    // compute gradients of pwl basis functions g1, g2, g3

    // get edge vectors
    Eigen::Matrix3d E0;
    E0.col(0) = _P0.col(1) - _P0.col(0);
    E0.col(1) = _P0.col(2) - _P0.col(0);
    E0.col(2) = _P0.col(3) - _P0.col(0);

    double d = E0.determinant();
    double di = 1.0 / d;
    if (!std::isfinite(di) || d < 0.0)
    {
      std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;
      di = 0.0;
    }
    Mat3x3dT G;
    G.col(0) = di * (E0.col(1)).cross(E0.col(2));
    G.col(1) = di * (E0.col(2)).cross(E0.col(0));
    G.col(2) = di * (E0.col(0)).cross(E0.col(1));

    // get edge vectors
    Mat3x3dT E;
    E.col(0) = _P.col(1) - _P.col(0);
    E.col(1) = _P.col(2) - _P.col(0);
    E.col(2) = _P.col(3) - _P.col(0);

    // calc Jacobi Matrix of map
//    Mat3x3dT JF = E * G.transpose() * _F; // hack
    Mat3x3dT JF = 2.0 * E * G.transpose() * _F;

    _optimal_sizing[0] = JF(0, 0);
    _optimal_sizing[1] = JF(1, 1);
    _optimal_sizing[2] = JF(2, 2);
  }

  static double max_exp_arg(const VecV &_x, const VecC &_c)
  {
    // get constants
    const double &w = _c[0];
    const double &w_exp = _c[1];
    const double &w_aniso = _c[2];

    // get gradient computation matrix
    Eigen::Map<Mat3x3dT> G((double *) &(_c[3]));

    // get frame matrix
    Eigen::Map<Mat3x3dT> F((double *) &(_c[12]));

    // get points of tetrahedron
    Eigen::Map<Mat3x4dT> P((double *) &(_x[0]));

    // get target size
    Eigen::Map<Vec3dT> S((double *) &(_x[12]));

    // get edge vectors
    Mat3x3dT E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    // calc Jacobi Matrix of map
//    Mat3x3dT JF = E*G.transpose()*F;
    Mat3x3dT JF = 2.0 * E * G.transpose() * F;

    // subtract target size
    JF(0, 0) -= S[0];
    JF(1, 1) -= S[1];
    JF(2, 2) -= S[2];

    double m = pow(JF(0, 0), 2);
    m = std::max(m, pow(JF(1, 1), 2));
    m = std::max(m, pow(JF(2, 2), 2));
    m = std::max(m, pow(JF(1, 0), 2));
    m = std::max(m, pow(JF(2, 0), 2));
    m = std::max(m, pow(JF(1, 2), 2));

    return m;
  }
};

// add derivatives -- no hessian projection required since convex!
using FrameFittingExpElement3D_TAD = COMISO::FiniteElementTinyAD<FrameFittingExpElement3D_F>;

//=============================================================================
} // namespace AlgoHex
//=============================================================================


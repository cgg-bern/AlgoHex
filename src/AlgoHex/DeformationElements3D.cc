/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include "DeformationElements3D.hh"
#include <cfloat>
#include <CoMISo/Utils/Polynomials.hh>

namespace AlgoHex
{

   double SymmetricDirichletElement3D_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    using Mat3x3T = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
    using Mat3x4T = Eigen::Matrix<double, 3, 4, Eigen::ColMajor>;

    // get constants
    const double &size_scale = _c[1];
    const double &epsilon = _c[2];

    // get gradient computation matrix
    Eigen::Map<Mat3x3T> G((double *) &(_c[3]));

    // get points of tetrahedron
    Eigen::Map<Mat3x4T> P((double *) &(_x[0]));

    // get edge vectors
    Mat3x3T E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    // get points of tetrahedron
    Eigen::Map<Mat3x4T> dP((double *) &(_v[0]));

    // get edge vectors
    Mat3x3T dE;
    dE.col(0) = dP.col(1) - dP.col(0);
    dE.col(1) = dP.col(2) - dP.col(0);
    dE.col(2) = dP.col(3) - dP.col(0);

    // coefficients computed via Mathematica (+index shift since Mathematica is 1 based instead of 0 based)
    double c0 = -(E(0, 2) * E(1, 1) * E(2, 0)) + E(0, 1) * E(1, 2) * E(2, 0) + E(0, 2) * E(1, 0) * E(2, 1) -
                E(0, 0) * E(1, 2) * E(2, 1) - E(0, 1) * E(1, 0) * E(2, 2) + E(0, 0) * E(1, 1) * E(2, 2);

    double c1 = -(dE(2, 2) * E(0, 1) * E(1, 0)) + dE(2, 1) * E(0, 2) * E(1, 0) + dE(2, 2) * E(0, 0) * E(1, 1) -
                dE(2, 0) * E(0, 2) * E(1, 1) - dE(2, 1) * E(0, 0) * E(1, 2) + dE(2, 0) * E(0, 1) * E(1, 2) +
                dE(1, 2) * E(0, 1) * E(2, 0) - dE(1, 1) * E(0, 2) * E(2, 0) - dE(0, 2) * E(1, 1) * E(2, 0) +
                dE(0, 1) * E(1, 2) * E(2, 0) - dE(1, 2) * E(0, 0) * E(2, 1) + dE(1, 0) * E(0, 2) * E(2, 1) +
                dE(0, 2) * E(1, 0) * E(2, 1) - dE(0, 0) * E(1, 2) * E(2, 1) + dE(1, 1) * E(0, 0) * E(2, 2) -
                dE(1, 0) * E(0, 1) * E(2, 2) - dE(0, 1) * E(1, 0) * E(2, 2) + dE(0, 0) * E(1, 1) * E(2, 2);

    double c2 = -(dE(1, 2) * dE(2, 1) * E(0, 0)) + dE(1, 1) * dE(2, 2) * E(0, 0) + dE(1, 2) * dE(2, 0) * E(0, 1) -
                dE(1, 0) * dE(2, 2) * E(0, 1) - dE(1, 1) * dE(2, 0) * E(0, 2) + dE(1, 0) * dE(2, 1) * E(0, 2) +
                dE(0, 2) * dE(2, 1) * E(1, 0) - dE(0, 1) * dE(2, 2) * E(1, 0) - dE(0, 2) * dE(2, 0) * E(1, 1) +
                dE(0, 0) * dE(2, 2) * E(1, 1) + dE(0, 1) * dE(2, 0) * E(1, 2) - dE(0, 0) * dE(2, 1) * E(1, 2) -
                dE(0, 2) * dE(1, 1) * E(2, 0) + dE(0, 1) * dE(1, 2) * E(2, 0) + dE(0, 2) * dE(1, 0) * E(2, 1) -
                dE(0, 0) * dE(1, 2) * E(2, 1) - dE(0, 1) * dE(1, 0) * E(2, 2) + dE(0, 0) * dE(1, 1) * E(2, 2);

    double c3 = -(dE(0, 2) * dE(1, 1) * dE(2, 0)) + dE(0, 1) * dE(1, 2) * dE(2, 0) + dE(0, 2) * dE(1, 0) * dE(2, 1) -
                dE(0, 0) * dE(1, 2) * dE(2, 1) - dE(0, 1) * dE(1, 0) * dE(2, 2) + dE(0, 0) * dE(1, 1) * dE(2, 2);

    // adjust coefficients
    double s = std::pow(size_scale, 3) * G.determinant();
    c0 *= s;
    c1 *= s;
    c2 *= s;
    c3 *= s;
    // shift w.r.t. epsilon
    c0 += epsilon;

    COMISO::Monomial<3> f;
    f.coeffs() << c0, c1, c2, c3;
    double t;
    if (COMISO::Polynomials::first_root_in_interval(f, 0.0, 10.0, t))
    {
//      // check result
//      std::cerr << "determinant at estimated root = " << (size_scale*(E+t*dE)*G.transpose()).determinant() + epsilon << std::endl;

      // reduce a little bit to avoid infeasibility of max_step
      return 0.99 * t;
    }
    else
      return DBL_MAX;
    return DBL_MAX;
  }

   void SymmetricDirichletElement3D_F::compute_constants(const double _w, const double _size_scale, const double _epsilon,
                                const Eigen::Vector3d &_p0, const Eigen::Vector3d &_p1, const Eigen::Vector3d &_p2,
                                const Eigen::Vector3d &_p3,
                                VecC &_c)
  {
    _c[0] = _w;
    _c[1] = _size_scale;
    _c[2] = _epsilon;

    // compute gradients of pwl basis functions g1, g2, g3

    // get edge vectors
    Eigen::Matrix3d E;
    E.col(0) = _p1 - _p0;
    E.col(1) = _p2 - _p0;
    E.col(2) = _p3 - _p0;

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

   double
  SymmetricDirichletElement3D_F::determine_suitable_epsilon(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c,
                             const double _target_penalty)
  {
    using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
    using Mat3x4dT = Eigen::Matrix<double, 3, 4, Eigen::ColMajor>;

    // get constants
    const double &w = _c[0];
    const double &size_scale = _c[1];
    const double &epsilon = _c[2];

    // get gradient computation matrix
    Eigen::Map<Mat3x3dT> G((double *) &(_c[3]));

    // get points of tetrahedron
    Eigen::Map<Mat3x4dT> P((double *) &(_x[0]));

    // get edge vectors
    Mat3x3dT E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    // calc (scaled) Jacobi Matrix of map
    Mat3x3dT J = size_scale * E * G.transpose();

    // calculate determinant
    double d = J.determinant();

    if (d > 0.0)
      return 0.0;
    else
      return 1.0 / _target_penalty - d;
  }

   double FoldoverFree3DV_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo
    return DBL_MAX;
  }

   void
  FoldoverFree3DV_F::compute_constants(const double _w_conformal, const double _w_det, const double _target_det, const double _epsilon,
                    const Eigen::Vector3d &_p0, const Eigen::Vector3d &_p1, const Eigen::Vector3d &_p2,
                    const Eigen::Vector3d &_p3,
                    VecC &_c)
  {
    _c[0] = _w_conformal;
    _c[1] = _w_det;
    _c[2] = _target_det;
    _c[3] = _epsilon;

    // compute gradients of pwl basis functions g1, g2, g3

    // get edge vectors
    Eigen::Matrix3d E;
    E.col(0) = _p1 - _p0;
    E.col(1) = _p2 - _p0;
    E.col(2) = _p3 - _p0;

    double d = E.determinant();
    double di = 1.0 / d;
    if (!std::isfinite(di) || d < 0.0)
    {
      std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << d << std::endl;
      di = 0.0;
    }

    // multiply weights with size of tet
    _c[0] *= d / 6.0;
    _c[1] *= d / 6.0;

    _c.segment<3>(4) = di * (E.col(1)).cross(E.col(2));
    _c.segment<3>(7) = di * (E.col(2)).cross(E.col(0));
    _c.segment<3>(10) = di * (E.col(0)).cross(E.col(1));

    // -------------------
    // test and debug code
    // -------------------
    if (0)
    {
      Eigen::Matrix3d J;
      J = _c.segment<3>(4) * E.col(0).transpose() +
          _c.segment<3>(7) * E.col(1).transpose() +
          _c.segment<3>(10) * E.col(2).transpose();

      J -= Eigen::Matrix3d::Identity();

      double Jn = J.norm();
      if (Jn > 1e-6)
      {
        std::cerr << "ERROR: Jacobi matrix was not computed correctly!" << std::endl;

        J = _c.segment<3>(4) * E.col(0).transpose() +
            _c.segment<3>(7) * E.col(1).transpose() +
            _c.segment<3>(10) * E.col(2).transpose();
        std::cerr << J << std::endl;
      }
      else
        std::cerr << "Jacobi matrix within tolerance ---> " << Jn << std::endl;
    }
  }

   double
  FoldoverFree3DV_F::determine_suitable_epsilon(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c)
  {
    using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
    using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get constants
    const double &w_conformal = _c[0];
    const double &w_det = _c[1];
    const double &target_det = _c[2];
    const double &epsilon = _c[3];

    Eigen::Map<Vec3dT> g1((double *) &(_c[4]));
    Eigen::Map<Vec3dT> g2((double *) &(_c[7]));
    Eigen::Map<Vec3dT> g3((double *) &(_c[10]));

    // get four input points
    Eigen::Map<Vec3dT> p0((double *) &(_x[0]));
    Eigen::Map<Vec3dT> p1((double *) &(_x[3]));
    Eigen::Map<Vec3dT> p2((double *) &(_x[6]));
    Eigen::Map<Vec3dT> p3((double *) &(_x[9]));

    // calc Jacobi Matrix of map
    Mat3x3dT J = g1 * (p1 - p0).transpose() + g2 * (p2 - p0).transpose() + g3 * (p3 - p0).transpose();

    // calculate determinant
    double d = J.determinant();

    double d_s = d / target_det;

    if (d_s > 0.0)
      return 0.0;
    else
    {
      //   ScalarT d_reg = 0.5 * (d_s + sqrt(pow(d_s, 2) + epsilon * epsilon));
      // return w_conformal * ((J.transpose() * J).trace() / pow(d_reg * target_det, 2.0 / 3.0)) +
      //             w_det * ((pow(d_s, 2) + 1.0) / d_reg);
      return std::sqrt(1e-12 + 4.0 * 1e-2 * d_s * d_s);
    }
  }

   double FoldoverFreeDual3DV_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo
    return DBL_MAX;
  }

   double AlignmentElement_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }

   double FittingElement_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }

   double DihedralAngle_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }

   double BallBarrier3D_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo
    return DBL_MAX;
  }

   double JacobianSmoothnessElement3D_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }

   void JacobianSmoothnessElement3D_F::compute_constants( const double _w,
                                 const double _sizing_scale,
                                 const Eigen::Vector3d& _p0, const Eigen::Vector3d& _p1, const Eigen::Vector3d& _p2, const Eigen::Vector3d& _p3,
                                 const Eigen::Vector3d& _q0, const Eigen::Vector3d& _q1, const Eigen::Vector3d& _q2, const Eigen::Vector3d& _q3,
                                 const Eigen::Matrix3d& _T01,
                                 VecC& _c)
  {
    _c[0] = _w;
    _c[28] = _sizing_scale;

    // compute gradients of pwl basis functions g1, g2, g3
    // get edge vectors
    Eigen::Matrix3d E0;
    E0.col(0) = _p1-_p0;
    E0.col(1) = _p2-_p0;
    E0.col(2) = _p3-_p0;

    double d0  = E0.determinant();
    double di0 = 1.0/d0;
    if(!std::isfinite(di0) || d0 < 0.0)
    {
      std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << d0 << std::endl;
      di0 = 0.0;
    }

    _c.segment<3>(1) = di0*(E0.col(1)).cross(E0.col(2));
    _c.segment<3>(4) = di0*(E0.col(2)).cross(E0.col(0));
    _c.segment<3>(7) = di0*(E0.col(0)).cross(E0.col(1));

    Eigen::Matrix3d E1;
    E1.col(0) = _q1-_q0;
    E1.col(1) = _q2-_q0;
    E1.col(2) = _q3-_q0;

    double d1  = E1.determinant();
    double di1 = 1.0/d1;
    if(!std::isfinite(di1) || d1 < 0.0)
    {
      std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << d1 << std::endl;
      di1 = 0.0;
    }

    _c.segment<3>(10) = di1*(E1.col(1)).cross(E1.col(2));
    _c.segment<3>(13) = di1*(E1.col(2)).cross(E1.col(0));
    _c.segment<3>(16) = di1*(E1.col(0)).cross(E1.col(1));

    // store transition function
    _c.segment<3>(19) = _T01.col(0);
    _c.segment<3>(22) = _T01.col(1);
    _c.segment<3>(25) = _T01.col(2);


    // -------------------
    // test and debug code
    // -------------------
//    if(0)
//    {
//      Eigen::Matrix3d J;
//      J =     _c.segment<3>(3)*E.col(0).transpose() +
//              _c.segment<3>(6)*E.col(1).transpose() +
//              _c.segment<3>(9)*E.col(2).transpose();
//
//      J -= Eigen::Matrix3d::Identity();
//
//      double Jn = J.norm();
//      if(Jn > 1e-6)
//      {
//        std::cerr << "ERROR: Jacobi matrix was not computed correctly!" << std::endl;
//
//        J =     _c.segment<3>(3)*E.col(0).transpose() +
//                _c.segment<3>(6)*E.col(1).transpose() +
//                _c.segment<3>(9)*E.col(2).transpose();
//        std::cerr << J << std::endl;
//      }
//      else
//        std::cerr << "Jacobi matrix within tolerance ---> " << Jn << std::endl;
//    }
  }


}

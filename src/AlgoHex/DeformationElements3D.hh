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
#include <CoMISo/Utils/Polynomials.hh>


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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
    using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
    using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using Mat3x4T = Eigen::Matrix<ScalarT, 3, 4, Eigen::ColMajor>;
    using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get constants
    const double &w = _c[0];
    const double &size_scale = _c[1];
    const double &epsilon = _c[2];

    // get gradient computation matrix
    Eigen::Map<Mat3x3dT> G((double *) &(_c[3]));

    // get points of tetrahedron
    Eigen::Map<Mat3x4T> P((ScalarT *) &(_x[0]));

    // get edge vectors
    Mat3x3T E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

//    ScalarT V = 1.0 / 6.0 * E.determinant();

    // calc (scaled) Jacobi Matrix of map
    Mat3x3T J = size_scale * E * G.transpose();

    // calculate determinant
    ScalarT d = J.determinant();

    if (epsilon == 0.0)
    {
      if (d <= 0.0)
        return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
      else
        return w*(J.squaredNorm() + J.inverse().squaredNorm() - 6.0);
//        return w*(1.0+d)*(J.squaredNorm() + J.inverse().squaredNorm() - 6.0); // symmetrized symmetric Dirichlet
    }
    else
    {
      // regularized determinant
      ScalarT d_reg = d + epsilon;

      if (d_reg <= 0.0)
        return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
      else
      {
        ScalarT di = 1.0 / d_reg;

        // compute regularized inverse
        Mat3x3T Ji;
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
          {
            int i1 = (i + 1) % 3;
            int j1 = (j + 1) % 3;
            int i2 = (i + 2) % 3;
            int j2 = (j + 2) % 3;

            Ji(j, i) = di * (J(i1, j1) * J(i2, j2) - J(i1, j2) * J(i2, j1));
          }

        return w * (J.squaredNorm() + Ji.squaredNorm() - 6.0);
      }
    }
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
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

  // NC = 12; // [weight, size_scale, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
  // calculate constants from reference tetrahedron
  static void compute_constants(const double _w, const double _size_scale, const double _epsilon,
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


  static double
  determine_suitable_epsilon(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c,
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
};

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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
    using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
    using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get constants
    const double &w_conformal = _c[0];
    const double &w_det = _c[1];
    const double &target_det = _c[2];
    const double &epsilon = _c[3];

//    Eigen::Map<Vec3dT> g1( (double*) &(_c[4]));
//    Eigen::Map<Vec3dT> g2( (double*) &(_c[7]));
//    Eigen::Map<Vec3dT> g3( (double*) &(_c[10]));

    Eigen::Map<Mat3x3dT> G((double *) &(_c[4]));

    // get four input points
    Eigen::Map<Vec3T> p0((ScalarT *) &(_x[0]));
    Eigen::Map<Vec3T> p1((ScalarT *) &(_x[3]));
    Eigen::Map<Vec3T> p2((ScalarT *) &(_x[6]));
    Eigen::Map<Vec3T> p3((ScalarT *) &(_x[9]));

    // get edge vectors
    Mat3x3T E;
    E.col(0) = p1 - p0;
    E.col(1) = p2 - p0;
    E.col(2) = p3 - p0;

    // calc Jacobi Matrix of map
    Mat3x3T J = E * G.transpose();

    // calculate determinant
    ScalarT d = J.determinant();

    if (epsilon > 0.0)
    {
      ScalarT d_s = d / target_det;
      ScalarT d_reg = 0.5 * (d_s + sqrt(pow(d_s, 2) + epsilon * epsilon));

//      std::cerr << std::endl
//                << "d=" << d << std::endl
//                << "d_s=" << d_s << std::endl
//                << "d_reg=" << d_reg << std::endl
//                << "epsilon=" << epsilon << std::endl
//                << "target_det=" << target_det << std::endl;

//      return w_conformal*((J.transpose()*J).trace()/pow(d_reg*target_det,2.0/3.0)) + w_det*((pow(d_reg,2)+1.0)/d_reg);
      return w_conformal * ((J.transpose() * J).trace() / pow(d_reg * target_det, 2.0 / 3.0)) +
             w_det * ((pow(d_s, 2) + 1.0) / d_reg);
    }
    else
    {
//      std::cerr << std::endl
//                << "d=" << d << std::endl
//                << "epsilon=" << epsilon << std::endl
//                << "target_det=" << target_det << std::endl;


      if (d <= 0.0)
        return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
      else
        return w_conformal * ((J.transpose() * J).trace() / pow(d, 2.0 / 3.0) - 3.0) +
               w_det * (d / target_det + target_det / d - 2.0);
      //      return w_conformal*((J.transpose()*J).trace()/pow(d,2.0/3.0)) + w_det*(d/target_det + target_det/d);
    }
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo
    return DBL_MAX;
  }

  //  const static int NC = 13; // [weight_conformal, weight_det, target_det, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
  // calculate constants from reference tetrahedron
  static void
  compute_constants(const double _w_conformal, const double _w_det, const double _target_det, const double _epsilon,
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

  static double
  determine_suitable_epsilon(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c)
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
};

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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
    using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
    using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get constants
    const double &w_conformal = _c[0];
    const double &w_det = _c[1];
    const double &target_det = _c[2];
    const double &epsilon = _c[3];


    // get four image input points
    Eigen::Map<Vec3T> p0((ScalarT *) &(_x[0]));
    Eigen::Map<Vec3T> p1((ScalarT *) &(_x[3]));
    Eigen::Map<Vec3T> p2((ScalarT *) &(_x[6]));
    Eigen::Map<Vec3T> p3((ScalarT *) &(_x[9]));

    // get four domain input points
    Eigen::Map<Vec3T> q0((ScalarT *) &(_x[12]));
    Eigen::Map<Vec3T> q1((ScalarT *) &(_x[15]));
    Eigen::Map<Vec3T> q2((ScalarT *) &(_x[18]));
    Eigen::Map<Vec3T> q3((ScalarT *) &(_x[21]));

    // compute gradients of pwl basis functions g1, g2, g3

    // get mesh edge vectors
    Mat3x3T E;
    E.col(0) = q1 - q0;
    E.col(1) = q2 - q0;
    E.col(2) = q3 - q0;

    ScalarT dE = E.determinant();
    ScalarT dEi = 1.0 / dE;
    if (!isfinite(dEi) || dE < 0.0)
    {
//      std::cerr << "ERROR: determinant of edge vectors numerically not valid ---> " << dE << std::endl;
//      dEi = 0.0;
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    }
    // volume of tet
    ScalarT vol = dE / 6.0;

//    Vec3T g1 = dEi*(E.col(1)).cross(E.col(2));
//    Vec3T g2 = dEi*(E.col(2)).cross(E.col(0));
//    Vec3T g3 = dEi*(E.col(0)).cross(E.col(1));
//
//    // calc Jacobi Matrix of map
//    Mat3x3T J = g1*(p1-p0).transpose() + g2*(p2-p0).transpose() + g3*(p3-p0).transpose();

    Mat3x3T G;
    G.col(0) = dEi * (E.col(1)).cross(E.col(2));
    G.col(1) = dEi * (E.col(2)).cross(E.col(0));
    G.col(2) = dEi * (E.col(0)).cross(E.col(1));

    // get parametric edge vectors
    Mat3x3T E2;
    E2.col(0) = p1 - p0;
    E2.col(1) = p2 - p0;
    E2.col(2) = p3 - p0;

    // calc Jacobi Matrix of map
    Mat3x3T J = E2 * G.transpose();

    // calculate determinant
    ScalarT d = J.determinant();

    if (epsilon > 0.0)
    {
      ScalarT d_s = d / target_det;
      ScalarT d_reg = 0.5 * (d_s + sqrt(pow(d_s, 2) + epsilon * epsilon));

//      std::cerr << std::endl
//                << "d=" << d << std::endl
//                << "d_s=" << d_s << std::endl
//                << "d_reg=" << d_reg << std::endl
//                << "epsilon=" << epsilon << std::endl
//                << "target_det=" << target_det << std::endl;

      return vol * (w_conformal * ((J.transpose() * J).trace() / pow(d_reg * target_det, 2.0 / 3.0)) +
                    w_det * ((pow(d_s, 2) + 1.0) / d_reg));
    }
    else
    {
//      std::cerr << std::endl
//                << "d=" << d << std::endl
//                << "epsilon=" << epsilon << std::endl
//                << "target_det=" << target_det << std::endl;


      if (d <= 0.0)
        return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
      else
        return vol * (w_conformal * ((J.transpose() * J).trace() / pow(d, 2.0 / 3.0) - 3.0) +
                      w_det * (d / target_det + target_det / d - 2.0));
    }
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo
    return DBL_MAX;
  }
};

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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    return _c[0] * pow(_x[0] - _x[1], 2);
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }
};

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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    return _c[0] * pow(_x[0] - _c[1], 2);
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }
};

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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
    using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
    using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get constants
    const double &w = _c[0];
    const double &ca = _c[1];
    const double &sa = _c[2];

    // get four input points
    Eigen::Map<Vec3T> p0((ScalarT *) &(_x[0]));
    Eigen::Map<Vec3T> p1((ScalarT *) &(_x[3]));
    Eigen::Map<Vec3T> p2((ScalarT *) &(_x[6]));
    Eigen::Map<Vec3T> p3((ScalarT *) &(_x[9]));

    Vec3T n0 = (p1 - p0).cross(p2 - p0);
    n0 /= n0.norm();

    Vec3T n1 = (p3 - p0).cross(p1 - p0);
    n1 /= n1.norm();

    return w * (pow(n0.dot(n1) - ca, 2) + pow((n0.cross(n1)).norm() - sa, 2));
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }
};

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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
    using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;

    // get constants
    const double &w = _c[0];
    const double &r = _c[1];

    Eigen::Map<Vec3dT> c((double *) &(_c[2]));

    // get four input points
    Eigen::Map<Vec3T> p((ScalarT *) &(_x[0]));

    ScalarT f = 1.0 - (p - c).squaredNorm() / (r * r);

    if (f <= 0.0)
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    else
      return -w * log(f);
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo
    return DBL_MAX;
  }
};

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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Vec3T    = Eigen::Matrix<ScalarT,3,1,Eigen::ColMajor>;
    using Vec3dT   = Eigen::Matrix<double ,3,1,Eigen::ColMajor>;
    using Mat3x3T  = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
    using Mat3x4T  = Eigen::Matrix<ScalarT,3,4,Eigen::ColMajor>;
    using Mat3x3dT = Eigen::Matrix<double ,3,3,Eigen::ColMajor>;

    // get constants
    const double& w          = _c[0];
    const double& size_scale = _c[28];


    // get gradient computation matrix
    Eigen::Map<Mat3x3dT> G0( (double*) &(_c[1]));
    Eigen::Map<Mat3x3dT> G1( (double*) &(_c[10]));
    // get transition function
    Eigen::Map<Mat3x3dT> T01( (double*) &(_c[19]));

    // get points of tetrahedron
    Eigen::Map<Mat3x4T> P0((ScalarT*) &(_x[0]));
    Eigen::Map<Mat3x4T> P1((ScalarT*) &(_x[12]));

    // get edge vectors
    Mat3x3T E0, E1;
    E0.col(0) = P0.col(1)-P0.col(0);
    E0.col(1) = P0.col(2)-P0.col(0);
    E0.col(2) = P0.col(3)-P0.col(0);

    E1.col(0) = P1.col(1)-P1.col(0);
    E1.col(1) = P1.col(2)-P1.col(0);
    E1.col(2) = P1.col(3)-P1.col(0);

//    ScalarT V0 = 1.0/6.0*E0.determinant();
//    ScalarT V1 = 1.0/6.0*E1.determinant();

    // calc (scaled) Jacobi Matrix of map
    Mat3x3T J0 = size_scale*E0*G0.transpose();
    Mat3x3T J1 = size_scale*E1*G1.transpose();

//    if(typeid(ScalarT) == typeid(double))
//    {
//      std::cerr << "w = " << w << std::endl;
//      std::cerr << "J0=" << std::endl << J0 << std::endl;
//      std::cerr << "T01=" << std::endl << T01 << std::endl;
//      std::cerr << "T01J0=" << std::endl << T01*J0 << std::endl;
//      std::cerr << "J1=" << std::endl << J1 << std::endl;
//      std::cerr << "T01*J0-J1=" << std::endl << T01*J0-J1 << std::endl;
//    }

    return w*(T01*J0-J1).squaredNorm();
  }


  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }

  //  const static int NC = 29; // [weight, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z, h1x,... h3z, t1x, ..., t3z, sizing_scale]
  // calculate constants from reference tetrahedron
  static void compute_constants( const double _w,
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
};

using JSE3D_AD = COMISO::FiniteElementTinyAD<JacobianSmoothnessElement3D_F>;


//=============================================================================
} // namespace AlgoHex
//=============================================================================


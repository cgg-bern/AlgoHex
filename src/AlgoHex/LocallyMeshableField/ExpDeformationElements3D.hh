/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/StdVector>
#include<cmath>

#include <CoMISo/NSolver/FiniteElementTinyAD.hh>
#include <CoMISo/NSolver/FiniteElementHessianProjection.hh>
#include <CoMISo/Utils/Polynomials.hh>
#include <AlgoHex/Geometry.hh>


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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Vec3T = Eigen::Matrix<ScalarT, 3, 1, Eigen::ColMajor>;
    using Vec3dT = Eigen::Matrix<double, 3, 1, Eigen::ColMajor>;
    using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using Mat3x3dT = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get constants
    const double &w_conformal = _c[0];

//            Eigen::Map<Mat3x3dT> invG( (double*) &(_c[10]));

    // get four input points
    Eigen::Map<Vec3T> pa((ScalarT *) &(_x[0]));
    Eigen::Map<Vec3dT> pb((double *) &(_c[1]));
    Eigen::Map<Vec3dT> pc((double *) &(_c[4]));
    Eigen::Map<Vec3dT> pd((double *) &(_c[7]));


    // edge vectors
    Mat3x3T E;
    E.col(0) = pb - pa;
    E.col(1) = pc - pa;
    E.col(2) = pd - pa;

    // calc Jacobi Matrix of map
//            Mat3x3T J = E * invG;

    Eigen::Map<Mat3x3dT> G((double *) &(_c[10]));
    Mat3x3T J = E * G.transpose();

    // calculate determinant
    ScalarT d = J.determinant();
//            std::cerr<<" d "<<d<<" invG "<<G.determinant()<<" E "<<E.determinant()<<" trace "<<(J.transpose() * J).trace()<<std::endl;

    if (d <= 0.0)
    {
//                double xx = pa[0];
//                Vec3dT paa((double)pa[0], (double)pa[1], (double)pa[2]);
//                Mat3x3dT E2;
//                E2.col(0) = pb - paa;
//                E2.col(1) = pc - paa;
//                E2.col(2) = pd - paa;

      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    }
    else
      return exp(w_conformal * ((J.transpose() * J).trace() / pow(d, 2.0 / 3.0) - 3.0));
  }

   double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

  // calculate constants from reference tetrahedron
  static void compute_constants(const double _w_conformal,
                                const Eigen::Vector3d &_p0, const Eigen::Vector3d &_p1,
                                const Eigen::Vector3d &_p2, const Eigen::Vector3d &_p3,
                                VecC &_c);

  static double
  eval_fx(const Eigen::Matrix<double, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c);

};

using EDE_AD = COMISO::FiniteElementTinyAD<ExpDeformationElement>;
using EDE_PH_AD = COMISO::FiniteElementHessianProjection<EDE_AD>;
}
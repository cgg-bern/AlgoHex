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
#include <Eigen/Eigenvalues>
#include<cmath>

#include <CoMISo/Utils/Polynomials.hh>
#include <CoMISo/NSolver/FiniteElementTinyAD.hh>
#include <CoMISo/NSolver/FiniteElementHessianProjection.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================


/** \class AMIPSFrameElement3DTinyAD AMIPSFrameElement3DTinyAD.hh
 *
 *  @brief Implementation of a minus log det element for optimization. f(J) = w*(exp(s*( alpha*(tr(J^TA^TAJ)*(det(AJ)^(-2/3)-3) + beta*(det(AJ)+det(AJ)^(-1) -2) ))-1.0)
 *
 */

// element function for f(J) = w*(exp(s*( alpha*(tr(J^TA^TAJ)*(det(AJ)^(-2/3)-3) + beta*(det(AJ)+det(AJ)^(-1) -2) ))-1.0)
// matrix is col-major

class AMIPSFrameElement3D_F
{
public:

  // define dimensions
  const static int NV = 9; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 14; // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], det(A) = _c[9], weight w=_c[10], distortion parameters s=_c[11], alpha=_c[12] and beta=_c[13]


  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using MatT = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get matrix representation of J and A
    Eigen::Map<MatT> J((ScalarT *) _x.data());
    Eigen::Map<MatD> A((double *) _c.data());

    // get references to constants
    const double &dA = _c[9];
    const double &w = _c[10];
    const double &s = _c[11];
    const double &alpha = _c[12];
    const double &beta = _c[13];

    // calculate determinant of AJ
    MatT AJ = A * J;
    ScalarT d = (AJ).determinant();

    if (d <= 0.0)
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    else
    {
      // f(J) = w*(exp(s*( alpha*(tr(J^TA^TAJ)*(det(AJ)^(-2/3)-3) + beta*(det(AJ)+det(AJ)^(-1) -2) ))-1.0)
//      ScalarT tr    = AJ.col(0).squaredNorm() + AJ.col(1).squaredNorm() + AJ.col(2).squaredNorm();
//      ScalarT Emips = tr*pow(d,-2.0/3.0) - 3.0;
      ScalarT Emips = (AJ.transpose() * AJ).trace() * pow(d, -2.0 / 3.0) - 3.0;
      ScalarT Edet = d + 1.0 / d - 2.0;

      return w * (exp(s * (alpha * Emips + beta * Edet)) - 1.0);
    }
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get matrix representation of J and dJ
    Eigen::Map<MatD> J((double *) _x.data());
    Eigen::Map<MatD> dJ((double *) _v.data());

    // coefficients computed via Mathematica (+index shift since Mathematica is 1 based instead of 0 based)
    double c0 = -(J(0, 2) * J(1, 1) * J(2, 0)) + J(0, 1) * J(1, 2) * J(2, 0) + J(0, 2) * J(1, 0) * J(2, 1) -
                J(0, 0) * J(1, 2) * J(2, 1) - J(0, 1) * J(1, 0) * J(2, 2) + J(0, 0) * J(1, 1) * J(2, 2);

    double c1 = -(dJ(2, 2) * J(0, 1) * J(1, 0)) + dJ(2, 1) * J(0, 2) * J(1, 0) + dJ(2, 2) * J(0, 0) * J(1, 1) -
                dJ(2, 0) * J(0, 2) * J(1, 1) - dJ(2, 1) * J(0, 0) * J(1, 2) + dJ(2, 0) * J(0, 1) * J(1, 2) +
                dJ(1, 2) * J(0, 1) * J(2, 0) - dJ(1, 1) * J(0, 2) * J(2, 0) - dJ(0, 2) * J(1, 1) * J(2, 0) +
                dJ(0, 1) * J(1, 2) * J(2, 0) - dJ(1, 2) * J(0, 0) * J(2, 1) + dJ(1, 0) * J(0, 2) * J(2, 1) +
                dJ(0, 2) * J(1, 0) * J(2, 1) - dJ(0, 0) * J(1, 2) * J(2, 1) + dJ(1, 1) * J(0, 0) * J(2, 2) -
                dJ(1, 0) * J(0, 1) * J(2, 2) - dJ(0, 1) * J(1, 0) * J(2, 2) + dJ(0, 0) * J(1, 1) * J(2, 2);

    double c2 = -(dJ(1, 2) * dJ(2, 1) * J(0, 0)) + dJ(1, 1) * dJ(2, 2) * J(0, 0) + dJ(1, 2) * dJ(2, 0) * J(0, 1) -
                dJ(1, 0) * dJ(2, 2) * J(0, 1) - dJ(1, 1) * dJ(2, 0) * J(0, 2) + dJ(1, 0) * dJ(2, 1) * J(0, 2) +
                dJ(0, 2) * dJ(2, 1) * J(1, 0) - dJ(0, 1) * dJ(2, 2) * J(1, 0) - dJ(0, 2) * dJ(2, 0) * J(1, 1) +
                dJ(0, 0) * dJ(2, 2) * J(1, 1) + dJ(0, 1) * dJ(2, 0) * J(1, 2) - dJ(0, 0) * dJ(2, 1) * J(1, 2) -
                dJ(0, 2) * dJ(1, 1) * J(2, 0) + dJ(0, 1) * dJ(1, 2) * J(2, 0) + dJ(0, 2) * dJ(1, 0) * J(2, 1) -
                dJ(0, 0) * dJ(1, 2) * J(2, 1) - dJ(0, 1) * dJ(1, 0) * J(2, 2) + dJ(0, 0) * dJ(1, 1) * J(2, 2);

    double c3 = -(dJ(0, 2) * dJ(1, 1) * dJ(2, 0)) + dJ(0, 1) * dJ(1, 2) * dJ(2, 0) + dJ(0, 2) * dJ(1, 0) * dJ(2, 1) -
                dJ(0, 0) * dJ(1, 2) * dJ(2, 1) - dJ(0, 1) * dJ(1, 0) * dJ(2, 2) + dJ(0, 0) * dJ(1, 1) * dJ(2, 2);

    COMISO::Monomial<3> f;
    f.coeffs() << c0, c1, c2, c3;
    double t;
    if (COMISO::Polynomials::first_root_in_interval(f, 0.0, 10.0, t))
    {
      // reduce a little bit to avoid infeasibility of max_step
      return 0.99 * t;
    }
    else
      return DBL_MAX;
  }
};

// create types for optimization with derivatives and projected hessian
using AMIPSFrameElement3D_TAD = COMISO::FiniteElementTinyAD<AMIPSFrameElement3D_F>;
using AMIPSFrameElement3D_TAD_PH = COMISO::FiniteElementHessianProjection<AMIPSFrameElement3D_TAD>;


class SymmetricDirichletFrameElement3D_F
{
public:

  // define dimensions
  const static int NV = 9; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 10; // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], det(A) = _c[9], weight w=_c[10], distortion parameters s=_c[11], alpha=_c[12] and beta=_c[13]


  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using MatT = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get matrix representation of J and A
    Eigen::Map<MatT> J((ScalarT *) _x.data());
    Eigen::Map<MatD> A((double *) _c.data());

    // get references to constants
    const double &w = _c[9];

    // calculate determinant of AJ
    MatT AJ = A * J;
    ScalarT d = (AJ).determinant();

    if (d <= 0.0)
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    else
      return w * (AJ.squaredNorm() + AJ.inverse().squaredNorm() - 6.0);
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get matrix representation of J and dJ
    Eigen::Map<MatD> J((double *) _x.data());
    Eigen::Map<MatD> dJ((double *) _v.data());

    // coefficients computed via Mathematica (+index shift since Mathematica is 1 based instead of 0 based)
    double c0 = -(J(0, 2) * J(1, 1) * J(2, 0)) + J(0, 1) * J(1, 2) * J(2, 0) + J(0, 2) * J(1, 0) * J(2, 1) -
                J(0, 0) * J(1, 2) * J(2, 1) - J(0, 1) * J(1, 0) * J(2, 2) + J(0, 0) * J(1, 1) * J(2, 2);

    double c1 = -(dJ(2, 2) * J(0, 1) * J(1, 0)) + dJ(2, 1) * J(0, 2) * J(1, 0) + dJ(2, 2) * J(0, 0) * J(1, 1) -
                dJ(2, 0) * J(0, 2) * J(1, 1) - dJ(2, 1) * J(0, 0) * J(1, 2) + dJ(2, 0) * J(0, 1) * J(1, 2) +
                dJ(1, 2) * J(0, 1) * J(2, 0) - dJ(1, 1) * J(0, 2) * J(2, 0) - dJ(0, 2) * J(1, 1) * J(2, 0) +
                dJ(0, 1) * J(1, 2) * J(2, 0) - dJ(1, 2) * J(0, 0) * J(2, 1) + dJ(1, 0) * J(0, 2) * J(2, 1) +
                dJ(0, 2) * J(1, 0) * J(2, 1) - dJ(0, 0) * J(1, 2) * J(2, 1) + dJ(1, 1) * J(0, 0) * J(2, 2) -
                dJ(1, 0) * J(0, 1) * J(2, 2) - dJ(0, 1) * J(1, 0) * J(2, 2) + dJ(0, 0) * J(1, 1) * J(2, 2);

    double c2 = -(dJ(1, 2) * dJ(2, 1) * J(0, 0)) + dJ(1, 1) * dJ(2, 2) * J(0, 0) + dJ(1, 2) * dJ(2, 0) * J(0, 1) -
                dJ(1, 0) * dJ(2, 2) * J(0, 1) - dJ(1, 1) * dJ(2, 0) * J(0, 2) + dJ(1, 0) * dJ(2, 1) * J(0, 2) +
                dJ(0, 2) * dJ(2, 1) * J(1, 0) - dJ(0, 1) * dJ(2, 2) * J(1, 0) - dJ(0, 2) * dJ(2, 0) * J(1, 1) +
                dJ(0, 0) * dJ(2, 2) * J(1, 1) + dJ(0, 1) * dJ(2, 0) * J(1, 2) - dJ(0, 0) * dJ(2, 1) * J(1, 2) -
                dJ(0, 2) * dJ(1, 1) * J(2, 0) + dJ(0, 1) * dJ(1, 2) * J(2, 0) + dJ(0, 2) * dJ(1, 0) * J(2, 1) -
                dJ(0, 0) * dJ(1, 2) * J(2, 1) - dJ(0, 1) * dJ(1, 0) * J(2, 2) + dJ(0, 0) * dJ(1, 1) * J(2, 2);

    double c3 = -(dJ(0, 2) * dJ(1, 1) * dJ(2, 0)) + dJ(0, 1) * dJ(1, 2) * dJ(2, 0) + dJ(0, 2) * dJ(1, 0) * dJ(2, 1) -
                dJ(0, 0) * dJ(1, 2) * dJ(2, 1) - dJ(0, 1) * dJ(1, 0) * dJ(2, 2) + dJ(0, 0) * dJ(1, 1) * dJ(2, 2);

    COMISO::Monomial<3> f;
    f.coeffs() << c0, c1, c2, c3;
    double t;
    if (COMISO::Polynomials::first_root_in_interval(f, 0.0, 10.0, t))
    {
      // reduce a little bit to avoid infeasibility of max_step
      return 0.99 * t;
    }
    else
      return DBL_MAX;
  }
};

// create types for optimization with derivatives and projected hessian
using SymmetricDirichletFrameElement3D_TAD = COMISO::FiniteElementTinyAD<SymmetricDirichletFrameElement3D_F>;
using SymmetricDirichletFrameElement3D_TAD_PH = COMISO::FiniteElementHessianProjection<SymmetricDirichletFrameElement3D_TAD>;
//using  SymmetricDirichletFrameElement3D_TAD_PH = COMISO::FiniteElementHessianProjectionIdentity< SymmetricDirichletFrameElement3D_TAD >;
//using  SymmetricDirichletFrameElement3D_TAD_PH = SymmetricDirichletFrameElement3D_TAD;


class SymmetricDirichletDualFrameElement3D_F
{
public:

  // define dimensions
  const static int NV = 21; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix, P=[p0,p1,p2,p3] = 3x4 matrix of tet vertex positions
  const static int NC = 10; // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], weight w=_c[9]


  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using Mat3x3T = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using Mat3x4T = Eigen::Matrix<ScalarT, 3, 4, Eigen::ColMajor>;
    using Mat3x3D = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get matrix representation of J and A
    Eigen::Map<Mat3x3T> J((ScalarT *) _x.data());
    Eigen::Map<Mat3x3D> A((double *) _c.data());

    // get points of tetrahedron
    Eigen::Map<Mat3x4T> P((ScalarT *) &(_x[9]));

    // get edge vectors
    Mat3x3T E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    ScalarT V = 1.0 / 6.0 * E.determinant();

    // get references to constants
    const double &w = _c[9];

    // calculate determinant of AJ
    Mat3x3T AJ = A * J;
    ScalarT d = (AJ).determinant();

    if (d <= 0.0)
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    else
      return w * V * (AJ.squaredNorm() + AJ.inverse().squaredNorm() - 6.0);
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // only check for frame part and rely on different terms to check validity of domain tetrahedron

    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;

    // get matrix representation of J and dJ
    Eigen::Map<MatD> J((double *) _x.data());
    Eigen::Map<MatD> dJ((double *) _v.data());

    // coefficients computed via Mathematica (+index shift since Mathematica is 1 based instead of 0 based)
    double c0 = -(J(0, 2) * J(1, 1) * J(2, 0)) + J(0, 1) * J(1, 2) * J(2, 0) + J(0, 2) * J(1, 0) * J(2, 1) -
                J(0, 0) * J(1, 2) * J(2, 1) - J(0, 1) * J(1, 0) * J(2, 2) + J(0, 0) * J(1, 1) * J(2, 2);

    double c1 = -(dJ(2, 2) * J(0, 1) * J(1, 0)) + dJ(2, 1) * J(0, 2) * J(1, 0) + dJ(2, 2) * J(0, 0) * J(1, 1) -
                dJ(2, 0) * J(0, 2) * J(1, 1) - dJ(2, 1) * J(0, 0) * J(1, 2) + dJ(2, 0) * J(0, 1) * J(1, 2) +
                dJ(1, 2) * J(0, 1) * J(2, 0) - dJ(1, 1) * J(0, 2) * J(2, 0) - dJ(0, 2) * J(1, 1) * J(2, 0) +
                dJ(0, 1) * J(1, 2) * J(2, 0) - dJ(1, 2) * J(0, 0) * J(2, 1) + dJ(1, 0) * J(0, 2) * J(2, 1) +
                dJ(0, 2) * J(1, 0) * J(2, 1) - dJ(0, 0) * J(1, 2) * J(2, 1) + dJ(1, 1) * J(0, 0) * J(2, 2) -
                dJ(1, 0) * J(0, 1) * J(2, 2) - dJ(0, 1) * J(1, 0) * J(2, 2) + dJ(0, 0) * J(1, 1) * J(2, 2);

    double c2 = -(dJ(1, 2) * dJ(2, 1) * J(0, 0)) + dJ(1, 1) * dJ(2, 2) * J(0, 0) + dJ(1, 2) * dJ(2, 0) * J(0, 1) -
                dJ(1, 0) * dJ(2, 2) * J(0, 1) - dJ(1, 1) * dJ(2, 0) * J(0, 2) + dJ(1, 0) * dJ(2, 1) * J(0, 2) +
                dJ(0, 2) * dJ(2, 1) * J(1, 0) - dJ(0, 1) * dJ(2, 2) * J(1, 0) - dJ(0, 2) * dJ(2, 0) * J(1, 1) +
                dJ(0, 0) * dJ(2, 2) * J(1, 1) + dJ(0, 1) * dJ(2, 0) * J(1, 2) - dJ(0, 0) * dJ(2, 1) * J(1, 2) -
                dJ(0, 2) * dJ(1, 1) * J(2, 0) + dJ(0, 1) * dJ(1, 2) * J(2, 0) + dJ(0, 2) * dJ(1, 0) * J(2, 1) -
                dJ(0, 0) * dJ(1, 2) * J(2, 1) - dJ(0, 1) * dJ(1, 0) * J(2, 2) + dJ(0, 0) * dJ(1, 1) * J(2, 2);

    double c3 = -(dJ(0, 2) * dJ(1, 1) * dJ(2, 0)) + dJ(0, 1) * dJ(1, 2) * dJ(2, 0) + dJ(0, 2) * dJ(1, 0) * dJ(2, 1) -
                dJ(0, 0) * dJ(1, 2) * dJ(2, 1) - dJ(0, 1) * dJ(1, 0) * dJ(2, 2) + dJ(0, 0) * dJ(1, 1) * dJ(2, 2);

    COMISO::Monomial<3> f;
    f.coeffs() << c0, c1, c2, c3;
    double t;
    if (COMISO::Polynomials::first_root_in_interval(f, 0.0, 10.0, t))
    {
      // reduce a little bit to avoid infeasibility of max_step
      return 0.99 * t;
    }
    else
      return DBL_MAX;
  }
};

// create types for optimization with derivatives and projected hessian
using SymmetricDirichletDualFrameElement3D_TAD = COMISO::FiniteElementTinyAD<SymmetricDirichletDualFrameElement3D_F>;
using SymmetricDirichletDualFrameElement3D_TAD_PH = COMISO::FiniteElementHessianProjection<SymmetricDirichletDualFrameElement3D_TAD>;


/*
class OptimalRotationFrameElement3D_F
{
public:

  // define dimensions
  const static int NV = 9; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 5; // ToDo


  using VecV = Eigen::Matrix<double,NV,1>;
  using VecC = Eigen::Matrix<double,NC,1>;

  template<class ScalarT>
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
  {
    using MatT = Eigen::Matrix<ScalarT,3,3,Eigen::ColMajor>;
    using MatD = Eigen::Matrix<double, 3,3,Eigen::ColMajor>;

    // get matrix representation of J and A
    Eigen::Map<MatT> J( (ScalarT*) _x.data() );

    ScalarT d = J.determinant();

    // get references to constants
    const double& w      = _c[0];
    const double& l_u    = _c[1];
    const double& l_v    = _c[2];
    const double& l_w    = _c[3];
    const double& aniso  = _c[4];

    if (d < 0.0)
    {
      std::cerr << "ERROR: optrot element at infinity!!!" << std::endl;
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    }
    else
    {
//      MatT JJt = J*J.transpose();
//      Eigen::SelfAdjointEigenSolver<MatT> es(JJt);

      Eigen::SelfAdjointEigenSolver<MatT> es(J*J.transpose());
      auto S = es.operatorSqrt();

//      {
//        MatT JJt = J * J.transpose();
//
//        std::cerr << "JJt for l = " << l_u << std::endl << JJt << std::endl;
//        std::cerr << "S   for l = " << l_u << std::endl << S << std::endl;
//      }

      return w*(pow(S(0,0)-l_u,2) + pow(S(1,1)-l_v,2) + pow(S(2,2)-l_w,2) + 2.0*pow(S(0,1),2) + 2.0*pow(S(0,2),2) + 2.0*pow(S(1,2),2));
    }
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    return DBL_MAX;
  }
};

// create types for optimization with derivatives and projected hessian
using  OptimalRotationFrameElement3D_TAD = COMISO::FiniteElementTinyAD<OptimalRotationFrameElement3D_F>;
using  OptimalRotationFrameElement3D_TAD_PH = COMISO::FiniteElementHessianProjection< OptimalRotationFrameElement3D_TAD >;
*/

//=============================================================================
} // namespace AlgoHex
//=============================================================================


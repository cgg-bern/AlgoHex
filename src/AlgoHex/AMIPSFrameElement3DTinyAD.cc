/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include "AMIPSFrameElement3DTinyAD.hh"

#include <CoMISo/Utils/Polynomials.hh>
#include <cfloat>

namespace AlgoHex
{

  double AMIPSFrameElement3D_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
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

  double SymmetricDirichletFrameElement3D_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
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

  double SymmetricDirichletDualFrameElement3D_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
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

template<class ScalarT>
ScalarT AMIPSFrameElement3D_F::eval_f(const Eigen::Matrix<ScalarT, AMIPSFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, AMIPSFrameElement3D_F::NC, 1> &_c) const
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

template<class ScalarT>
ScalarT SymmetricDirichletFrameElement3D_F::eval_f(const Eigen::Matrix<ScalarT, SymmetricDirichletFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NC, 1> &_c) const
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

template<class ScalarT>
ScalarT SymmetricDirichletDualFrameElement3D_F::eval_f(const Eigen::Matrix<ScalarT, SymmetricDirichletDualFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NC, 1> &_c) const
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

template double AMIPSFrameElement3D_F::eval_f<double>(const Eigen::Matrix<double, AMIPSFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, AMIPSFrameElement3D_F::NC, 1> &_c) const;
template TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, false> AMIPSFrameElement3D_F::eval_f<TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, false>, AMIPSFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, AMIPSFrameElement3D_F::NC, 1> &_c) const;
template TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, true> AMIPSFrameElement3D_F::eval_f<TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, true>, AMIPSFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, AMIPSFrameElement3D_F::NC, 1> &_c) const;

template double SymmetricDirichletFrameElement3D_F::eval_f<double>(const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NC, 1> &_c) const;
template TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, false> SymmetricDirichletFrameElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, false>, SymmetricDirichletFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NC, 1> &_c) const;
template TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, true> SymmetricDirichletFrameElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, true>, SymmetricDirichletFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NC, 1> &_c) const;

template double SymmetricDirichletDualFrameElement3D_F::eval_f<double>(const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NC, 1> &_c) const;
template TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, false> SymmetricDirichletDualFrameElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, false>, SymmetricDirichletDualFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NC, 1> &_c) const;
template TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, true> SymmetricDirichletDualFrameElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, true>, SymmetricDirichletDualFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NC, 1> &_c) const;

} // namespace AlgoHex

namespace COMISO
{
template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::AMIPSFrameElement3D_F>>;
template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::SymmetricDirichletFrameElement3D_F>>;
template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::SymmetricDirichletDualFrameElement3D_F>>;
} // namespace COMISO

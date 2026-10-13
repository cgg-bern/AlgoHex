/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

//== EXPLICIT INSTANTIATIONS ==================================================
//
//  eval_f definitions and all TinyAD wrapper instantiations of the element(s)
//  declared in DihedralAngleElementTinyAD.hh live here, so every wrapper chain is
//  instantiated exactly once. The matching `extern template` declarations in the
//  header suppress implicit instantiation in consumer translation units.

#include "DihedralAngleElementTinyAD.hh"

#include <cfloat>
#include <math.h>

namespace AlgoHex
{

template<class ScalarT>
ScalarT DihedralAngleFramesElement_F::eval_f(const Eigen::Matrix<ScalarT, DihedralAngleFramesElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NC, 1> &_c) const
{

    using MatT = Eigen::Matrix<ScalarT, 3, 3, Eigen::ColMajor>;
    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
    using VecT = Eigen::Matrix<ScalarT, 3, 1>;

    // get matrix representation of J and E
    Eigen::Map<MatT> J((ScalarT *) _x.data());
    Eigen::Map<MatD> E((double *) _c.data());

    // get references to constants
    const double &w = _c[9];
    const double &s = _c[10];
    const double &cos_alpha = _c[11];

    // deform tetrahedron (edge vectors)
    MatT V = J * E;

    // calculate normal vectors
    VecT n0 = V.col(1).cross(V.col(0));
    n0 /= n0.norm();
    VecT n1 = V.col(2).cross(V.col(0));
    n1 /= n1.norm();
    // dot product
    ScalarT dp = n0.dot(n1);

    if (!isfinite(dp))
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    else
      return w * (exp(s * pow(dp - cos_alpha, 2)) - 1.0);

}

template double DihedralAngleFramesElement_F::eval_f<double>(const Eigen::Matrix<double, DihedralAngleFramesElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NC, 1> &_c) const;
template TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, false> DihedralAngleFramesElement_F::eval_f<TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, false>, DihedralAngleFramesElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NC, 1> &_c) const;
template TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, true> DihedralAngleFramesElement_F::eval_f<TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngleFramesElement_F::NV, double, true>, DihedralAngleFramesElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NC, 1> &_c) const;

double DihedralAngleFramesElement_F::max_feasible_step(const Eigen::Matrix<double, DihedralAngleFramesElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NV, 1> &_v, const Eigen::Matrix<double, DihedralAngleFramesElement_F::NC, 1> &_c)
{
  return DBL_MAX;
}

template<class ScalarT>
ScalarT DihedralAngleElement_F::eval_f(const Eigen::Matrix<ScalarT, DihedralAngleElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleElement_F::NC, 1> &_c) const
{

    using MatT = Eigen::Matrix<ScalarT, 3, 4, Eigen::ColMajor>;
    using MatD = Eigen::Matrix<double, 3, 3, Eigen::ColMajor>;
    using VecT = Eigen::Matrix<ScalarT, 3, 1>;

    // get matrix representation of J and E
    Eigen::Map<MatT> P((ScalarT *) _x.data());

    // get references to constants
    const double &w = _c[0];
    const double &s = _c[1];
    const double &cos_alpha = _c[2];
    const double &sin_alpha = _c[3];

    VecT e0 = P.col(1) - P.col(0);
    VecT e1 = P.col(2) - P.col(0);
    VecT e2 = P.col(3) - P.col(0);

    // calculate normal vectors
    VecT n0 = e0.cross(e1);
    n0 /= n0.norm();
    VecT n1 = e0.cross(e2);
    n1 /= n1.norm();
    // dot product
    ScalarT dp = n0.dot(n1);

    // cross product
    ScalarT cp = (n0.cross(n1)).dot(e0) / e0.norm();

//    std::cerr << "dp-cos_alpha = " << dp-cos_alpha << std::endl;
//    std::cerr << "cp-sin_alpha = " << cp-sin_alpha << std::endl;

//    if(!isfinite(dp))
//      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    if (!isfinite(dp) || !isfinite(cp))
      return w * (exp(2.0 * s) - 1.0);
    else
      return w * (exp(s * (pow(dp - cos_alpha, 2) + pow(cp - sin_alpha, 2))) - 1.0);

}

template double DihedralAngleElement_F::eval_f<double>(const Eigen::Matrix<double, DihedralAngleElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleElement_F::NC, 1> &_c) const;
template TinyAD::Scalar<DihedralAngleElement_F::NV, double, false> DihedralAngleElement_F::eval_f<TinyAD::Scalar<DihedralAngleElement_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngleElement_F::NV, double, false>, DihedralAngleElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleElement_F::NC, 1> &_c) const;
template TinyAD::Scalar<DihedralAngleElement_F::NV, double, true> DihedralAngleElement_F::eval_f<TinyAD::Scalar<DihedralAngleElement_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<DihedralAngleElement_F::NV, double, true>, DihedralAngleElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleElement_F::NC, 1> &_c) const;

double DihedralAngleElement_F::max_feasible_step(const Eigen::Matrix<double, DihedralAngleElement_F::NV, 1> &_x, const Eigen::Matrix<double, DihedralAngleElement_F::NV, 1> &_v, const Eigen::Matrix<double, DihedralAngleElement_F::NC, 1> &_c)
{
  return DBL_MAX;
}

} // namespace AlgoHex

namespace COMISO
{
template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::DihedralAngleFramesElement_F>>;
template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::DihedralAngleElement_F>>;
} // namespace COMISO

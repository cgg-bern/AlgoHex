/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include "BarrierElements.hh"
#include <cfloat>

namespace AlgoHex
{

double ReciprocalBarrierElement_F::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
{
  if (_v[0] > 0.0)
    return DBL_MAX;
  else
  {
    double t = 0.9 * (_c[0] - _x[0]) / _v[0];

    if (!std::isfinite(t))
    {
      std::cerr << "Warning: max_feasbile_step of ReciprocalBarrierElement_F generated t = " << t
                << ", x=" << _x << ", v=" << _v << ", c= " << _c.transpose() << std::endl;
      return DBL_MAX;
    }
    else return t;
  }
}

template<class ScalarT>
ScalarT ReciprocalBarrierElement_F::eval_f(const Eigen::Matrix<ScalarT, ReciprocalBarrierElement_F::NV, 1> &_x, const Eigen::Matrix<double, ReciprocalBarrierElement_F::NC, 1> &_c) const
{

    // get constants
    const double &b = _c[0];
    const double &w = _c[1];

    ScalarT v = _x[0] - b;

    if (v <= 0.0)
      return static_cast<ScalarT>(std::numeric_limits<double>::infinity());
    else
      return w * 1.0 / v;

}

template double ReciprocalBarrierElement_F::eval_f<double>(const Eigen::Matrix<double, ReciprocalBarrierElement_F::NV, 1> &_x, const Eigen::Matrix<double, ReciprocalBarrierElement_F::NC, 1> &_c) const;
template TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, false> ReciprocalBarrierElement_F::eval_f<TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, false>, ReciprocalBarrierElement_F::NV, 1> &_x, const Eigen::Matrix<double, ReciprocalBarrierElement_F::NC, 1> &_c) const;
template TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, true> ReciprocalBarrierElement_F::eval_f<TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, true>, ReciprocalBarrierElement_F::NV, 1> &_x, const Eigen::Matrix<double, ReciprocalBarrierElement_F::NC, 1> &_c) const;

} // namespace AlgoHex

namespace COMISO
{
template class FiniteElementTinyAD<AlgoHex::ReciprocalBarrierElement_F>;
} // namespace COMISO

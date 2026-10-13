/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once


//== INCLUDES =================================================================

#include <Eigen/Core>

#include <CoMISo/NSolver/FiniteElementTinyAD.hh>


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================

// implements f(x) = w*(1/(x-b) for (x-b) > 0 and f(x) = inf for (x-b) <= 0
// this is a barrier for x > b and requires a feasible initial x!

class ReciprocalBarrierElement_F
{
public:

  // define dimensions
  const static int NV = 1; // [x]
  const static int NC = 2; // [b=barrier_value, w=weight]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

};

extern template double ReciprocalBarrierElement_F::eval_f<double>(const Eigen::Matrix<double, ReciprocalBarrierElement_F::NV, 1> &_x, const Eigen::Matrix<double, ReciprocalBarrierElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, false> ReciprocalBarrierElement_F::eval_f<TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, false>, ReciprocalBarrierElement_F::NV, 1> &_x, const Eigen::Matrix<double, ReciprocalBarrierElement_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, true> ReciprocalBarrierElement_F::eval_f<TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<ReciprocalBarrierElement_F::NV, double, true>, ReciprocalBarrierElement_F::NV, 1> &_x, const Eigen::Matrix<double, ReciprocalBarrierElement_F::NC, 1> &_c) const;

// AD Element -- Hessian projection not required since convex!
using ReciprocalBarrierElement_TAD = COMISO::FiniteElementTinyAD<ReciprocalBarrierElement_F>;

//=============================================================================
} // namespace AlgoHex
//=============================================================================

namespace COMISO
{
extern template class FiniteElementTinyAD<AlgoHex::ReciprocalBarrierElement_F>;
} // namespace COMISO

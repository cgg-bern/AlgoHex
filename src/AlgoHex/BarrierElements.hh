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
  inline ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const
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

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
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
};

// AD Element -- Hessian projection not required since convex!
using ReciprocalBarrierElement_TAD = COMISO::FiniteElementTinyAD<ReciprocalBarrierElement_F>;

//=============================================================================
} // namespace AlgoHex
//=============================================================================


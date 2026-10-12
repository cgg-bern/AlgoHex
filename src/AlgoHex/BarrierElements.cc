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


}

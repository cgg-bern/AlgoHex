/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include "AMIPSFrameElement3DTinyAD.hh"
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


}

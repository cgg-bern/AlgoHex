/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "SHCoeffRotation.hh"

#include <Eigen/Core>

namespace AlgoHex
{
namespace SHRotation
{

const Mat9d RB_x_pi2 = (Mat9d() <<
                                0, 0, 0, 0, 0, sqrt(14.) / 4., 0, -sqrt(2.) / 4., 0,
        0, -3. / 4., 0, sqrt(7.) / 4., 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, sqrt(2.) / 4., 0, sqrt(14.) / 4., 0,
        0, sqrt(7.) / 4., 0, 3. / 4., 0, 0, 0, 0, 0,
        0, 0, 0, 0, 3. / 8., 0, sqrt(5.) / 4., 0, sqrt(35.) / 8.,
        -sqrt(14.) / 4., 0, -sqrt(2.) / 4., 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, sqrt(5.) / 4., 0, 1. / 2., 0, -sqrt(7.) / 4.,
        sqrt(2.) / 4., 0, -sqrt(14.) / 4., 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, sqrt(35.) / 8., 0, -sqrt(7.) / 4., 0, 1. / 8.).finished();

const Mat9d RB_x_minus_pi2 = RB_x_pi2.transpose();
const Mat9d RB_y_pi2 = RB_x_minus_pi2 * (RB_z(M_PI / 2) * RB_x_pi2);
const Mat9d RB_y_minus_pi2 = RB_y_pi2.transpose();

template<typename Scalar>
/// Efficiently compute sin and cos of (gamma, 2*gamma, 3*gamma, 4*gamma)
/// using trigonometric identities
struct SinCos4Multiples
{
  constexpr SinCos4Multiples(Scalar gamma)
          : s1(sin(gamma)), c1(cos(gamma)), s2(2 * s1 * c1), c2(2 * c1 * c1 - 1), s3(s1 * c2 + c1 * s2),
            c3(c1 * c2 - s1 * s2), s4(2 * s2 * c2), c4(2 * c2 * c2 - 1) {}

  Scalar s1, c1, s2, c2, s3, c3, s4, c4;
};

Mat9d RB_z(double radians)
{
  auto sc = SinCos4Multiples(radians);
  return (Mat9d() <<
                  sc.c4, 0, 0, 0, 0, 0, 0, 0, sc.s4,
          0, sc.c3, 0, 0, 0, 0, 0, sc.s3, 0,
          0, 0, sc.c2, 0, 0, 0, sc.s2, 0, 0,
          0, 0, 0, sc.c1, 0, sc.s1, 0, 0, 0,
          0, 0, 0, 0, 1, 0, 0, 0, 0,
          0, 0, 0, -sc.s1, 0, sc.c1, 0, 0, 0,
          0, 0, -sc.s2, 0, 0, 0, sc.c2, 0, 0,
          0, -sc.s3, 0, 0, 0, 0, 0, sc.c3, 0,
          -sc.s4, 0, 0, 0, 0, 0, 0, 0, sc.c4
  ).finished();
}


Mat9d RB_x(double radians)
{
  return RB_y_pi2 * RB_z(radians) * RB_y_minus_pi2;
}

Mat9d RB_y(double radians)
{
  return RB_x_minus_pi2 * RB_z(radians) * RB_x_pi2;
}

Mat9d RB_from_quat(const Quaternion &q)
{
  auto angles = q.toRotationMatrix().eulerAngles(2, 1, 0);
  return RB_z(angles[0]) * RB_y(angles[1]) * RB_x(angles[2]);
}

Vec9d rot_z(double radians, const Vec9d &x)
{
  auto sc = SinCos4Multiples(radians);

  return (Vec9d() <<
                  sc.c4 * x[0] + sc.s4 * x[8],
          sc.c3 * x[1] + sc.s3 * x[7],
          sc.c2 * x[2] + sc.s2 * x[6],
          sc.c1 * x[3] + sc.s1 * x[5],
          x[4],
          -sc.s1 * x[3] + sc.c1 * x[5],
          -sc.s2 * x[2] + sc.c2 * x[6],
          -sc.s3 * x[1] + sc.c3 * x[7],
          -sc.s4 * x[0] + sc.c4 * x[8]
  ).finished();

}

}
}  // namespace


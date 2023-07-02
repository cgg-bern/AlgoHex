/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "SHCoeffs.hh"
#include "SHCoeffRotation.hh"

#include <Eigen/Core>

namespace AlgoHex
{

using SHRotation::RB_x_pi2;
using SHRotation::RB_x_minus_pi2;
using SHRotation::RB_y_pi2;
using SHRotation::RB_y_minus_pi2;
using SHRotation::RB_z;
using SHRotation::rot_z;


SHCoeffs SHCoeffs::from_quaternion(const Quaternion &q)
{
  auto a = SHCoeffs::axis_aligned_cross();
  a.rotate(q);
  return a;
}

void SHCoeffs::rotate_x(double radians)
{
  *this = RB_y_pi2
          * (RB_z(radians)
             * (RB_y_minus_pi2
                * vec()));
}

void SHCoeffs::rotate_y(double radians)
{
  *this = RB_x_minus_pi2
          * rot_z(radians,
                  (RB_x_pi2 * vec()));
}

void SHCoeffs::rotate_z(double radians)
{
  *this = rot_z(radians, vec());
  // TODO PERF: inspect compiled code to see if unrolling the matrix multiplication is cheaper
}

void SHCoeffs::rotate_xyz(const Vec3d &radians)
{
  rotate_x(radians[0]);
  rotate_y(radians[1]);
  rotate_z(radians[2]);
}

void SHCoeffs::rotate(const Quaternion &q)
{
  auto angles = q.toRotationMatrix().eulerAngles(2, 1, 0);
  rotate_xyz({angles[2], angles[1], angles[0]});
}

SHCoeffs Ex(const SHCoeffs &a)
{
  return {{
                  -sqrt(2.) * a[7],
                  -sqrt(2.) * a[8] - sqrt(3.5) * a[6],
                  -sqrt(3.5) * a[7] - sqrt(4.5) * a[5],
                  -sqrt(4.5) * a[6] - sqrt(10.) * a[4],
                  sqrt(10.) * a[3],
                  sqrt(4.5) * a[2],
                  sqrt(3.5) * a[1] + sqrt(4.5) * a[3],
                  sqrt(2.) * a[0] + sqrt(3.5) * a[2],
                  sqrt(2.) * a[1]}};
}


SHCoeffs Ey(const SHCoeffs &a)
{
  return {{
                  sqrt(2.) * a[1],
                  -sqrt(2.) * a[0] + sqrt(3.5) * a[2],
                  -sqrt(3.5) * a[1] + sqrt(4.5) * a[3],
                  -sqrt(4.5) * a[2],
                  -sqrt(10.) * a[5],
                  -sqrt(4.5) * a[6] + sqrt(10.) * a[4],
                  -sqrt(3.5) * a[7] + sqrt(4.5) * a[5],
                  -sqrt(2.) * a[8] + sqrt(3.5) * a[6],
                  sqrt(2.) * a[7]}};
}


SHCoeffs Ez(const SHCoeffs &a)
{
  return {{4 * a[8], 3 * a[7], 2 * a[6],
           a[5], 0, -a[3],
           -2 * a[2], -3 * a[1], -4 * a[0]}};
}

std::ostream &operator<<(std::ostream &s, const SHCoeffs &shc)
{
  static const Eigen::IOFormat CommaInitFmt(8, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
  s << shc.vec().transpose().format(CommaInitFmt);
  return s;
}

}

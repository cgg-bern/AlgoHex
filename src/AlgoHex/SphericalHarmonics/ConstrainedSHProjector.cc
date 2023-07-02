/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "ConstrainedSHProjector.hh"
#include "SHCoeffRotation.hh"

namespace AlgoHex
{

ConstrainedSHProjector::
ConstrainedSHProjector(const Vec3d &v)
        : v_to_z_(Quaternion::FromTwoVectors(v, Vec3d(0, 0, 1)))
{
  Mat9d rot = SHRotation::RB_from_quat(v_to_z_);
  r1_ = rot.row(0);
  r9_ = rot.row(8);
}

SHProjectionResult
ConstrainedSHProjector::project(const SHCoeffs &sh_coeffs)
{
  // cf. Section 5.1.1 of
  // Palmer, D., Bommes, D., & Solomon, J. (2020). Algebraic Representations for Volumetric Frame Fields. ACM Transactions on Graphics (TOG), 39(2), 1-17.

  const double c1 = r1_.dot(sh_coeffs);
  const double c9 = r9_.dot(sh_coeffs);

  const double proj = sqrt(12. / 5) / (c1 * c1 + c9 * c9 + 7. / 12);

  const double s4t = c1 * proj;
  const double c4t = c9 * proj;

  const double t = .25 * atan2(s4t, c4t);
  // TODO: maybe we can avoid inverse trig functions here
  Quaternion axis_q(Eigen::AngleAxisd(t, Vec3d::UnitZ()));
  SHCoeffs test_shc = SHCoeffs::from_quaternion(axis_q);

  auto q = v_to_z_.conjugate() * axis_q;

  SHCoeffs shc = SHCoeffs::from_quaternion(q); // TODO: compute directly with a precomputed SHC rotation
  return {shc, q, 0};
}


} // namespace AlgoHex

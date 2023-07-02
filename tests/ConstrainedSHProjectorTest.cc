#include "gtest/gtest.h"
#include <AlgoHex/SphericalHarmonics/ConstrainedSHProjector.hh>

#include <AlgoHex/Config/Version.hh>
#include <iostream>


namespace AlgoHex
{

class ConstrainedSHProjectorTest : public ::testing::Test
{
};

SHCoeffs bruteforce(const Vec3d &axis, const SHCoeffs &qshc)
{
  const int steps = 360;
  const auto to_axis = Quaternion::FromTwoVectors(Vec3d(0, 0, 1), axis);

  double max_dot = -1000;
  SHCoeffs best_shc;

  for (int i = 0; i < steps; ++i)
  {
    const double angle = i * (M_PI / (2 * steps));
    auto axis_rotq = Quaternion(Eigen::AngleAxisd(angle, Vec3d::UnitZ()));
    const auto shc0 = SHCoeffs::from_quaternion(axis_rotq);
    auto shc = shc0;
    shc.rotate(to_axis);
    const auto dot = shc.dot(qshc);
    if (dot > max_dot)
    {
      max_dot = dot;
      best_shc = shc;
    }
  }
  return best_shc;
}

auto expect_aligned_to_axis = [](const Vec3d &axis, SHCoeffs shc) {
  const auto axis_to_z = Quaternion::FromTwoVectors(axis, Vec3d(0, 0, 1));
  shc.rotate(axis_to_z);
  // Now we only need to test if shc is z-aligned
  const double tol = 1e-6;
  EXPECT_NEAR(shc[1], 0., tol);
  EXPECT_NEAR(shc[2], 0., tol);
  EXPECT_NEAR(shc[3], 0., tol);
  EXPECT_NEAR(shc[4], sqrt(7. / 12), tol);
  EXPECT_NEAR(shc[5], 0., tol);
  EXPECT_NEAR(shc[6], 0., tol);
  EXPECT_NEAR(shc[7], 0., tol);

};

TEST_F(ConstrainedSHProjectorTest, basic_test
)
{
auto test = [](Vec3d axis, const SHCoeffs &shc) {
  axis.normalize();
  ConstrainedSHProjector cshp(axis);
  auto bf_shc = bruteforce(axis, shc);
  auto pr = cshp.project(shc);
  expect_aligned_to_axis(axis, pr.shc);
  //expect_aligned_to_axis(axis, bf_shc);
  // TODO: test that the results actually are aligned with axis
  auto bf_err = (shc - bf_shc).squaredNorm();
  auto err = (shc - pr.shc).squaredNorm();
  EXPECT_GE(bf_err * 1.01, err); // add a little bit to bruteforce err to have wiggle room for numerics
};
ConstrainedSHProjector P_y{Vec3d{0, 1, 0}};
SHCoeffs ref({0, 0, 0, 0, std::sqrt(7. / 12.), 0, 0, 0, std::sqrt(5. / 12.)});
SHCoeffs a{
        {0.22373389, 0.089979071, 0.14032772, 0.33282704, 0.43818923, 0.48919913, 0.14274704, -0.25976235, 0.5405575}};
for (
const auto &axis
:
{
Vec3d
{
0., 2., 1.
}
,
Vec3d{
1., 2., 1.},
Vec3d{
1., 2., 0.}})
{
test(axis, a
);
test(axis, ref
);
}
}

} // namespace AlgoHex


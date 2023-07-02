#include "gtest/gtest.h"

#include <AlgoHex/SphericalHarmonics/SHCoeffs.hh>
#include <AlgoHex/SphericalHarmonics/SHProjectorRay.hh>
#include <iostream>

namespace AlgoHex
{

class SphericalHarmonicsTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
  }
};

TEST_F(SphericalHarmonicsTest, sh_from_quaternion
)
{
auto test = [](const Quaternion &q, const SHCoeffs &a) {
  const double tol = 1e-6;

  auto r = SHCoeffs::from_quaternion(q);
  EXPECT_NEAR((a - r).norm(), 0, tol);

  for (const auto &rot:
          {Eigen::AngleAxisd(M_PI, Vec3d::UnitX()),
           Eigen::AngleAxisd(M_PI / 2, Vec3d::UnitX()),
           Eigen::AngleAxisd(M_PI / 2, Vec3d::UnitY()),
           Eigen::AngleAxisd(M_PI / 2, Vec3d::UnitZ())})
  {
    r = SHCoeffs::from_quaternion(q * rot);
    EXPECT_NEAR((a - r).norm(), 0, tol);

    SHProjectorRay proj;
    EXPECT_NEAR((proj.project(r).shc - r).norm(), 0, tol);

  }
};

SHCoeffs ref;

ref = {0, 0, 0, 0, std::sqrt(7. / 12.), 0, 0, 0, std::sqrt(5. / 12.)};
test(Quaternion(Eigen::AngleAxisd(0, Vec3d::UnitX())), ref
);

ref = {0, 0, 0, 0, -0.190941, 0, -0.853913, 0, 0.484123};
test(Quaternion(Eigen::AngleAxisd(M_PI / 4., Vec3d::UnitX())), ref
);

ref = {0, 0, 0, 0, -0.190941, 0, 0.853913, 0, 0.484123};
test(Quaternion(Eigen::AngleAxisd(M_PI / 4., Vec3d::UnitY())), ref
);

ref = {0, 0, 0, 0, 0.763763, 0, 0, 0, -0.645497};
test(Quaternion(Eigen::AngleAxisd(M_PI / 4., Vec3d::UnitZ())), ref
);

}

}

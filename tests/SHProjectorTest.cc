#include "gtest/gtest.h"
#include <AlgoHex/SphericalHarmonics/SHProjectorRay.hh>
#include <AlgoHex/SphericalHarmonics/SHProjectorSDP.hh>
#include <AlgoHex/SphericalHarmonics/OctaToQuat.hh>

#include <AlgoHex/Config/Version.hh>
#include <iostream>


namespace AlgoHex
{

#if ALGOHEX_WITH_MOSEK
class SHProjectorSDPTest : public ::testing::Test
{
protected:
    void SetUp() override {
    }
    SHProjectorSDP p;
};

TEST_F(SHProjectorSDPTest, basic_test)
{
    auto test = [this](const Quaternion &refq, const SHCoeffs &a)
    {
#if 0
        const double tol = 1e-6;
        SHCoeffs res = p.project_only_shc(a);
        EXPECT_NEAR((a-res).norm(), 0, tol);
        Quaternion q = from_projected_shc(res);
        auto r = SHCoeffs::from_quaternion(q);
        EXPECT_NEAR((a-r).norm(), 0, tol);
#endif
    };

    SHCoeffs ref;

    ref = {0, 0, 0, 0, std::sqrt(7. / 12.), 0, 0, 0, std::sqrt(5. / 12.)};
    test(Quaternion(Eigen::AngleAxisd(0, Vec3d::UnitX())), ref);

    ref = {0, 0, 0, 0, -0.190941, 0, -0.853913, 0, 0.484123};
    test(Quaternion(Eigen::AngleAxisd(M_PI/4., Vec3d::UnitX())), ref);

    ref = {0, 0, 0, 0, -0.190941, 0, 0.853913, 0, 0.484123};
    test(Quaternion(Eigen::AngleAxisd(M_PI/4., Vec3d::UnitY())), ref);

    ref = {0, 0, 0, 0, 0.763763, 0, 0, 0, -0.645497};
    test(Quaternion(Eigen::AngleAxisd(M_PI/4., Vec3d::UnitZ())), ref);
}
TEST_F(SHProjectorSDPTest, test_from_projected_shc)
{
    SHCoeffs a {{0.22373389, 0.089979071, 0.14032772, 0.33282704, 0.43818923, 0.48919913, 0.14274704, -0.25976235, 0.5405575}};
    a.normalize();

    const double tol = 1e-6;
    auto res = p.project(a);
    EXPECT_NEAR((a - res.shc).norm(), 0, tol);
    auto rec = SHCoeffs::from_quaternion(res.q);
    EXPECT_NEAR((a - rec).norm(), 0, tol);
}
#endif // ALGOHEX_WITH_MOSEK

} // namespace AlgoHex


/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "OctaToQuat.hh"
#include <random>

namespace AlgoHex
{


static Vec3d octa_tensor_gradient(const SHCoeffs &q, Vec3d v)
{
  const double x = v.x();
  const double y = v.y();
  const double z = v.z();

  Vec9d dx;
  Vec9d dy;
  Vec9d dz;

  using std::pow;
  const double pi = M_PI;

  double x0 = 2 * pow(pi, -0.5) * x * (pow(x, 2) + pow(y, 2) + pow(z, 2));
  double y0 = 2 * pow(pi, -0.5) * y * (pow(x, 2) + pow(y, 2) + pow(z, 2));
  double z0 = 2 * pow(pi, -0.5) * z * (pow(x, 2) + pow(y, 2) + pow(z, 2));
  const double fac = sqrt(189.) / 4.;


  dx << -3.0 / 4.0 * sqrt(35) * y * (-3 * pow(x, 2) + pow(y, 2)) * sqrt(1.0 / pi),
          (9.0 / 4.0) * sqrt(70) * x * y * z * sqrt(1.0 / pi),
          -3.0 / 4.0 * sqrt(5) * y * (3 * pow(x, 2) + pow(y, 2) - 6 * pow(z, 2)) * sqrt(1.0 / pi),
          -9.0 / 4.0 * sqrt(10) * x * y * z * sqrt(1.0 / pi),
          (9.0 / 4.0) * x * (pow(x, 2) + pow(y, 2) - 4 * pow(z, 2)) / sqrt(pi),
          (3.0 / 8.0) * sqrt(10) * z * (-9 * pow(x, 2) - 3 * pow(y, 2) + 4 * pow(z, 2)) * sqrt(1.0 / pi),
          -3.0 / 2.0 * sqrt(5) * x * (pow(x, 2) - 3 * pow(z, 2)) * sqrt(1.0 / pi),
          (9.0 / 8.0) * sqrt(70) * z * (pow(x, 2) - pow(y, 2)) * sqrt(1.0 / pi),
          (3.0 / 4.0) * sqrt(35) * (pow(x, 3) - 3 * x * pow(y, 2)) * sqrt(1.0 / pi);
  dy << (3.0 / 4.0) * sqrt(35) * x * (pow(x, 2) - 3 * pow(y, 2)) * sqrt(1.0 / pi),
          (9.0 / 8.0) * sqrt(70) * z * (pow(x, 2) - pow(y, 2)) * sqrt(1.0 / pi),
          -3.0 / 4.0 * sqrt(5) * x * (pow(x, 2) + 3 * pow(y, 2) - 6 * pow(z, 2)) * sqrt(1.0 / pi),
          (3.0 / 8.0) * sqrt(10) * z * (-3 * pow(x, 2) - 9 * pow(y, 2) + 4 * pow(z, 2)) * sqrt(1.0 / pi),
          (9.0 / 4.0) * y * (pow(x, 2) + pow(y, 2) - 4 * pow(z, 2)) / sqrt(pi),
          -9.0 / 4.0 * sqrt(10) * x * y * z * sqrt(1.0 / pi),
          (3.0 / 2.0) * sqrt(5) * y * (pow(y, 2) - 3 * pow(z, 2)) * sqrt(1.0 / pi),
          -9.0 / 4.0 * sqrt(70) * x * y * z * sqrt(1.0 / pi),
          -3.0 / 4.0 * sqrt(35) * (3 * pow(x, 2) * y - pow(y, 3)) * sqrt(1.0 / pi);
  dz << 0,
          -3.0 / 8.0 * sqrt(70) * y * (-3 * pow(x, 2) + pow(y, 2)) * sqrt(1.0 / pi),
          9 * sqrt(5) * x * y * z * sqrt(1.0 / pi),
          -9.0 / 8.0 * sqrt(10) * y * (pow(x, 2) + pow(y, 2) - 4 * pow(z, 2)) * sqrt(1.0 / pi),
          (-9 * pow(x, 2) * z - 9 * pow(y, 2) * z + 6 * pow(z, 3)) / sqrt(pi),
          -9.0 / 8.0 * sqrt(10) * x * (pow(x, 2) + pow(y, 2) - 4 * pow(z, 2)) * sqrt(1.0 / pi),
          (9.0 / 2.0) * sqrt(5) * z * (pow(x, 2) - pow(y, 2)) * sqrt(1.0 / pi),
          (3.0 / 8.0) * sqrt(70) * x * (pow(x, 2) - 3 * pow(y, 2)) * sqrt(1.0 / pi),
          0;

  return {q.dot(dx) + fac * x0,
          q.dot(dy) + fac * y0,
          q.dot(dz) + fac * z0};
}

/// Based on David Palmer's ARFF Octa2Frames
Quaternion ALGOHEX_EXPORT from_projected_shc(const SHCoeffs &shc)
{
  // Retry until we succeed, TODO: fix this properly.
  auto uni = std::uniform_real_distribution(-1., 1.);
  std::mt19937 mt(0);
  for (int i = 0; i < 1000; ++i)
  {
    Vec3d v1{uni(mt), uni(mt), uni(mt)};
    Vec3d v2{uni(mt), uni(mt), uni(mt)};

    double delta = 1;

    while (delta > 1e-14)
    {
      Vec3d old_v1 = v1;
      Vec3d w1 = octa_tensor_gradient(shc, v1).normalized();
      Vec3d w2 = octa_tensor_gradient(shc, v2);

      v1 = w1.normalized();
      v2 = (w2 - w2.dot(v1) * v1).normalized();

      delta = (old_v1 - v1).squaredNorm();
    }

    Mat3d m;
    m.col(0) = v1;
    m.col(1) = v2;
    m.col(2) = v1.cross(v2);
    Quaternion q = Quaternion(m);
    auto rec = SHCoeffs::from_quaternion(q);
    if ((rec - shc).norm() < 1e-3)
      return q;
  }
  std::cerr << "from_projected_shc failed for " << shc << std::endl;
  assert(false);
  return {0, 0, 0, 0};
}

} // namespace AlgoHex

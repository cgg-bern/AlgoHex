/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "SHProjectorRay.hh"

namespace AlgoHex
{

const std::vector<SHProjectorRay::Seed> SHProjectorRay::seeds5{{
                                                                       Quaternion(Eigen::AngleAxisd(0, Vec3d::UnitX())),
                                                                       Quaternion(Eigen::AngleAxisd(M_PI / 4.,
                                                                                                    Vec3d::UnitX())),
                                                                       Quaternion(Eigen::AngleAxisd(M_PI / 4.,
                                                                                                    Vec3d::UnitY())),
                                                                       Quaternion(Eigen::AngleAxisd(M_PI / 4.,
                                                                                                    Vec3d::UnitZ())),
                                                                       Quaternion(Eigen::AngleAxisd(M_PI / 4.,
                                                                                                    Vec3d(1., 1.,
                                                                                                          0.).normalized()))
                                                               }};

SHProjectionResult
SHProjectorRay::project(const SHCoeffs &sh_coeffs)
{
  return project_with_seeds(sh_coeffs, seeds5);
}

SHProjectionResult
SHProjectorRay::project_with_seed(
        const SHCoeffs &sh_coeffs,
        const SHProjectorRay::Seed &seed
) const
{
  const double sq_grad_norm_threshold = grad_norm_threshold_ * grad_norm_threshold_;

  SHCoeffs cur_shc{seed.shc};
  Quaternion cur_qt{seed.q};

  double old_dot_prd = sh_coeffs.dot(seed.shc);

  LocalEulerGradient grad(sh_coeffs);
  /*
  LocalEulerHessian hessian(sh_coeffs);
          auto H = hessian.eval(cur_shc);
          */

  size_t iter = 0;
  for (; iter < max_iter_; ++iter)
  {
    Vec3d g = grad.eval(cur_shc);

    if (g.squaredNorm() < sq_grad_norm_threshold)
    {


      break;
    }

    g *= 0.125;

    cur_qt = Quaternion(Eigen::AngleAxisd(g(2), Vec3d::UnitZ()))
             * (Quaternion(Eigen::AngleAxisd(g(1), Vec3d::UnitY()))
                * (Quaternion(Eigen::AngleAxisd(g(0), Vec3d::UnitX()))
                   * cur_qt));

    cur_shc.rotate_xyz(g);

    double dot_prd = sh_coeffs.dot(cur_shc);
    if (dot_prd - old_dot_prd < improvement_threshold_)
      break;
    old_dot_prd = dot_prd;
  }

#if 0
  if (iter == max_iter_)
      std::cerr << "\nWarning: SPH projection did not converge!" << std::endl;
#endif
  return {cur_shc, cur_qt, iter};
}

SHProjectionResult SHProjectorRay::project_with_seeds(const SHCoeffs &sh_coeffs,
                                                      const std::vector<SHProjectorRay::Seed> &seeds) const
{
  switch (seed_selection_)
  {
    case ClosestSeed: return project_using_closest_seed(sh_coeffs, seeds);
    case BestSeed: return project_using_best_seed(sh_coeffs, seeds);
    default: assert(false);
      return {};
  }
}

std::vector<SHProjectorRay::Seed>::const_iterator
SHProjectorRay::find_closest_seed(const SHCoeffs &sh_coeffs,
                                  const std::vector<SHProjectorRay::Seed> &seeds) const
{
  std::vector<double> dotprods(seeds.size());
  std::transform(seeds.begin(), seeds.end(),
                 dotprods.begin(),
                 [&](const Seed &seed) { return seed.shc.dot(sh_coeffs); });
  int seed_idx = std::distance(dotprods.begin(), std::max_element(dotprods.begin(), dotprods.end()));
  return seeds.begin() + seed_idx;
}

SHProjectionResult SHProjectorRay::project_using_closest_seed(const SHCoeffs &sh_coeffs,
                                                              const std::vector<SHProjectorRay::Seed> &seeds) const
{
  return project_with_seed(sh_coeffs, *find_closest_seed(sh_coeffs, seeds));
}

SHProjectionResult SHProjectorRay::project_using_best_seed(const SHCoeffs &sh_coeffs,
                                                           const std::vector<SHProjectorRay::Seed> &seeds) const
{
  double best_dotprod = -2.;
  SHProjectionResult best_result;
  for (const auto &seed: seeds)
  {
    auto result = project_with_seed(sh_coeffs, seed);
    double dotprod = sh_coeffs.dot(result.shc);
    if (dotprod > best_dotprod)
    {
      best_dotprod = dotprod;
      best_result = result;
    }
  }
  assert(best_dotprod > -2.);

  return best_result;
}

LocalEulerGradient::LocalEulerGradient(const SHCoeffs &q)
{
  M_.row(0) = -Ex(q).vec().transpose();
  M_.row(1) = -Ey(q).vec().transpose();
  M_.row(2) = -Ez(q).vec().transpose();
}

Vec3d LocalEulerGradient::eval(const SHCoeffs &a) const
{
  return M_ * a.vec();
}

LocalEulerHessian::LocalEulerHessian(const SHCoeffs &q)
{
  M_ <<
     -2 * q[0] - sqrt(7) * q[2], -11 * q[1] / 2 - 3 * sqrt(7) * q[3] / 2, -sqrt(7) * q[0] - 8 * q[2],
          -3 * sqrt(7) * q[1] / 2 - 29 * q[3] / 2, -10 * q[4] - 3 * sqrt(5) * q[6], -9 * q[5] / 2 -
                                                                                    3 * sqrt(7) * q[7] / 2,
          -3 * sqrt(5) * q[4] - 8 * q[6] - sqrt(7) * q[8], -3 * sqrt(7) * q[5] / 2 - 11 * q[7] / 2, -sqrt(7) * q[6] -
                                                                                                    2 * q[8],
          -sqrt(7) * q[6] - 2 * q[8], -3 * sqrt(7) * q[5] / 2 - 3 * q[7] / 2, -3 * sqrt(5) * q[4] - q[6] +
                                                                              sqrt(7) * q[8], 9 * q[5] / 2 +
                                                                                              3 * sqrt(7) * q[7] / 2,
          -3 * sqrt(5) * q[2], -3 * sqrt(7) * q[1] / 2 + 11 * q[3] / 2, -sqrt(7) * q[0] + q[2], 3 * q[1] / 2 +
                                                                                                3 * sqrt(7) * q[3] / 2,
          2 * q[0] + sqrt(7) * q[2],
          -4 * sqrt(2) * q[1], -3 * sqrt(2) * q[0] - 3 * sqrt(14) * q[2] / 2, -sqrt(14) * q[1] - 3 * sqrt(2) * q[3],
          -3 * sqrt(2) * q[2] / 2, 0, -sqrt(10) * q[4] - 3 * sqrt(2) * q[6] / 2, -3 * sqrt(2) * q[5] - sqrt(14) * q[7],
          -3 * sqrt(14) * q[6] / 2 - 3 * sqrt(2) * q[8], -4 * sqrt(2) * q[7],
          -2 * q[0] + sqrt(7) * q[2], -11 * q[1] / 2 + 3 * sqrt(7) * q[3] / 2, sqrt(7) * q[0] - 8 * q[2],
          3 * sqrt(7) * q[1] / 2 - 9 * q[3] / 2, -10 * q[4] + 3 * sqrt(5) * q[6], -29 * q[5] / 2 +
                                                                                  3 * sqrt(7) * q[7] / 2,
          3 * sqrt(5) * q[4] - 8 * q[6] + sqrt(7) * q[8], 3 * sqrt(7) * q[5] / 2 - 11 * q[7] / 2, sqrt(7) * q[6] -
                                                                                                  2 * q[8],
          -4 * sqrt(2) * q[7], -3 * sqrt(14) * q[6] / 2 + 3 * sqrt(2) * q[8], -3 * sqrt(2) * q[5] + sqrt(14) * q[7],
          -sqrt(10) * q[4] + 3 * sqrt(2) * q[6] / 2, 0, -3 * sqrt(2) * q[2] / 2, -sqrt(14) * q[1] + 3 * sqrt(2) * q[3],
          -3 * sqrt(2) * q[0] + 3 * sqrt(14) * q[2] / 2, 4 * sqrt(2) * q[1],
          -16 * q[0], -9 * q[1], -4 * q[2], -q[3], 0, -q[5], -4 * q[6], -9 * q[7], -16 * q[8];
}

Mat3d LocalEulerHessian::eval(const SHCoeffs &a) const
{
  auto vals = M_ * a.vec();
  return (Mat3d() << vals[0], vals[1], vals[2],
          vals[1], vals[3], vals[4],
          vals[2], vals[4], vals[5]).finished();
}


}

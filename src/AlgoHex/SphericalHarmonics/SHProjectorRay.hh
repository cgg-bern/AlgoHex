/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "SHProjectorInterface.hh"
#include <AlgoHex/Config/Export.hh>

namespace AlgoHex
{

class ALGOHEX_EXPORT LocalEulerGradient
{
public:
  LocalEulerGradient(const SHCoeffs &q);

  Vec3d eval(const SHCoeffs &a) const;

private:
  Eigen::Matrix<double, 3, 9> M_;
};

class ALGOHEX_EXPORT LocalEulerHessian
{
public:
  LocalEulerHessian(const SHCoeffs &q);

  Mat3d eval(const SHCoeffs &a) const;

private:
  Eigen::Matrix<double, 6, 9> M_;
};


/// Implementation of linearized projection from
/// "On Smooth Frame Field Design", Ray & Sokolov
class ALGOHEX_EXPORT SHProjectorRay : public SHProjectorInterface
{
public:
  struct Seed
  {
    Seed() = default; // uninitialized!
    Seed(const SHCoeffs &_shc, const Quaternion &_q)
            : shc(_shc), q(_q) {}

    Seed(const Quaternion &_q)
            : shc(SHCoeffs::from_quaternion(_q)), q(_q) {}

    SHCoeffs shc;
    Quaternion q;
  };

  // SHProjectorInterface methods:
  SHProjectionResult project(const SHCoeffs &sh_coeffs) override;

  SHProjectionResult project_with_seed(const SHCoeffs &sh_coeffs, const Seed &seed) const;

  SHProjectionResult project_with_seeds(const SHCoeffs &sh_coeffs, const std::vector<Seed> &seeds) const;

  // settings:
  enum Mode
  {
    ClosestSeed,   // choose closest seed (Original Ray) - fast
    BestSeed       // project with every seed, pick best result (Heng Liu) - slow, but fewer failures

  };

  void set_seed_selection(Mode value) { seed_selection_ = value; }

  void set_max_iter(int value) { max_iter_ = value; }

  void set_improvement_threshold(double value) { improvement_threshold_ = value; }

  void set_grad_norm_threshold(double value) { grad_norm_threshold_ = value; }

  std::vector<Seed>::const_iterator
  find_closest_seed(const SHCoeffs &sh_coeffs,
                    const std::vector<Seed> &seeds) const;

private:
  static const std::vector<Seed> seeds5;

  SHProjectionResult project_using_closest_seed(const SHCoeffs &sh_coeffs,
                                                const std::vector<Seed> &seeds) const;

  SHProjectionResult project_using_best_seed(const SHCoeffs &sh_coeffs,
                                             const std::vector<Seed> &seeds) const;

  Mode seed_selection_ = ClosestSeed;
  int max_iter_ = 10000;
  double improvement_threshold_ = 0;
  double grad_norm_threshold_ = 1e-8;

};

} // namespace AlgoHex

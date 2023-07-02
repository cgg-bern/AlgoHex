/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include "TypeDef.hh"
#include <iostream>
#include <iomanip>
#include "assert.h"

namespace AlgoHex
{
template<typename Point>
static Vec3d ovm2eigen(const Point &_p)
{
  return Eigen::Vector3d(_p[0], _p[1], _p[2]);
}

static OpenVolumeMesh::Vec3d eigen2ovm(const Vec3d &_p)
{
  return OpenVolumeMesh::Vec3d(_p[0], _p[1], _p[2]);
}


static double signed_angle(const Vec2d &_a, const Vec2d &_b)
{
  double s = 1.0 / (_a.norm() * _b.norm());
  double cos_alpha = _a.dot(_b) * s;
  double sin_alpha = (_a[0] * _b[1] - _a[1] * _b[0]) * s;

  if (!std::isfinite(cos_alpha))
    std::cerr << "Warning: signed_angle has been called with degenerate vectors, a=" << _a.transpose() << ", b="
              << _b.transpose() << std::endl;

  // clamp
  if (cos_alpha > 1.0) cos_alpha = 1.0;
  if (cos_alpha < -1.0) cos_alpha = -1.0;

  double alpha = std::acos(cos_alpha);

  if (sin_alpha > 0) return alpha;
  else return -alpha;
}


static void complement_to_right_handed_orthonormal_frame(const Vec3d &_base_dir, Vec3d &_x0, Vec3d &_x1, Vec3d &_x2)
{
  _x0 = _base_dir / _base_dir.norm();

  if (!std::isfinite(_x0[0] + _x0[1] + _x0[2]))
  {
    std::cerr << "Warning: complement_to_right_handed_frame received degenerate base direction "
              << _base_dir.transpose() << std::endl;
    _x0 = Vec3d(1, 0, 0);
  }

  if (std::abs(_x0[0]) > std::abs(_x0[1]))
    _x1 = _x0.cross(Vec3d(0, 1, 0));
  else
    _x1 = _x0.cross(Vec3d(1, 0, 0));

  _x1 /= _x1.norm();
  _x2 = _x0.cross(_x1);
}

static double tetrahedron_quality(const Vec3d &_p0, const Vec3d &_p1, const Vec3d &_p2, const Vec3d &_p3)
{
  Mat3d E;
  E.col(0) = _p1 - _p0;
  E.col(1) = _p2 - _p0;
  E.col(2) = _p3 - _p0;

  double V = 1.0 / 6.0 * E.determinant();

  double mean_el =
          (E.col(0).squaredNorm() + E.col(1).squaredNorm() + E.col(2).squaredNorm() + (_p2 - _p1).squaredNorm() +
           (_p3 - _p2).squaredNorm() + (_p1 - _p3).squaredNorm()) / 6.0;

  return (12.0 / std::sqrt(2.0)) * V / std::pow(mean_el, 3.0 / 2.0);
}

static double avg_tetrahedron_edge_length(const Vec3d &_p0, const Vec3d &_p1, const Vec3d &_p2, const Vec3d &_p3)
{
  double avg = 0.0;
  avg += (_p1 - _p0).norm();
  avg += (_p2 - _p0).norm();
  avg += (_p3 - _p0).norm();
  avg += (_p2 - _p1).norm();
  avg += (_p3 - _p2).norm();
  avg += (_p1 - _p3).norm();

  return avg / 6.0;
}

static void construct_uniform_tetrahedron(Vec3d &_p0, Vec3d &_p1, Vec3d &_p2, Vec3d &_p3, const double _edge_length)
{
  // first vertex in origin
  _p0 = Vec3d(0., 0., 0.);
  // second vertex aligned to x axis
  _p1 = Vec3d(_edge_length, 0., 0.);
  // third vertex in xy-plane
  _p2 = Vec3d(0.5 * _edge_length, 0.5 * std::sqrt(3) * _edge_length, 0.);
  // third vertex in barycenter w.r.t. xy-plane and height according to Pythagorean theorem
  _p3 = Vec3d(0.5 * _edge_length, 1.0 / 6.0 * std::sqrt(3) * _edge_length, std::sqrt(6.) / 3.0 * _edge_length);
}

} // namespace AlgoHex

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
template<typename VecT>
static Vec3d anyVec2eigen(const VecT &_p)
{
  return Eigen::Vector3d(_p[0], _p[1], _p[2]);
}

template<typename VecT>
static OpenVolumeMesh::Vec3d anyVec2ovm(const VecT &_p)
{
  return OpenVolumeMesh::Vec3d(_p[0], _p[1], _p[2]);
}

template<typename VecT>
static TetMesh::PointT anyVec2point(const VecT &_p)
{
  return TetMesh::PointT(_p[0], _p[1], _p[2]);
}

double signed_angle(const Vec2d &_a, const Vec2d &_b);

void complement_to_right_handed_orthonormal_frame(const Vec3d &_base_dir, Vec3d &_x0, Vec3d &_x1, Vec3d &_x2);

double tetrahedron_quality(const Vec3d &_p0, const Vec3d &_p1, const Vec3d &_p2, const Vec3d &_p3);

double avg_tetrahedron_edge_length(const Vec3d &_p0, const Vec3d &_p1, const Vec3d &_p2, const Vec3d &_p3);

void construct_uniform_tetrahedron(Vec3d &_p0, Vec3d &_p1, Vec3d &_p2, Vec3d &_p3, const double _edge_length);
} // namespace AlgoHex

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

TetMesh::PointT
ALGOHEX_EXPORT get_vertex_normal(const TetMesh &mesh, const VH _vh);

template<typename PointT>
double dihedral_angle(const PointT &pt0, const PointT &pt1, const PointT &n0, const PointT &n1);

extern template ALGOHEX_EXPORT double dihedral_angle(const OVM::Geometry::Vec3d&, const OVM::Geometry::Vec3d&, const OVM::Geometry::Vec3d&, const OVM::Geometry::Vec3d&);
extern template ALGOHEX_EXPORT double dihedral_angle(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);

double
ALGOHEX_EXPORT dihedral_angle_in_cell(const TetMesh &mesh, const HEH heh, const HFH hfh);

double
ALGOHEX_EXPORT dihedral_angle_from_halfface_to_halfface(const TetMesh &mesh, const HEH heh, const HFH hfh_s, const HFH hfh_e);

//dihedral angle of boundary halffaces
double
ALGOHEX_EXPORT dihedral_angle(const TetMesh &mesh, const HEH _heh);

double ALGOHEX_EXPORT get_cell_volume(const std::vector<TetMesh::PointT> &cell_points);

double ALGOHEX_EXPORT get_cell_volume(const TetMesh &mesh, const CH ch);

double ALGOHEX_EXPORT get_mesh_volume(const TetMesh &_mesh);

double ALGOHEX_EXPORT tetmesh_volume(const TetMesh &mesh);

bool ALGOHEX_EXPORT check_cells_volume(const TetMesh &mesh);

void ALGOHEX_EXPORT check_mesh_dihedral_angle(const TetMesh &mesh);

bool ALGOHEX_EXPORT check_mesh_quality(const TetMesh &mesh);

void
ALGOHEX_EXPORT scale_mesh(TetMesh &_mesh, double _s);

double
ALGOHEX_EXPORT avg_local_edge_length(TetMesh &_mesh, const VH _vh);


} // namespace AlgoHex

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include "TypeDef.hh"
#include <vector>

namespace AlgoHex
{

using VFullConstraints = std::vector<std::pair<VH, Quaternion>>;
VFullConstraints &get_vertex_full_constraints(TetMesh &mesh);

using VNormalConstraints = std::vector<std::pair<VH, TetMesh::PointT>>;
VNormalConstraints &get_vertex_normal_constraints(TetMesh &mesh);

using CPartialConstraints = std::vector<std::pair<CH, TetMesh::PointT>>;
CPartialConstraints &get_cell_partial_constraints(TetMesh &mesh);

}

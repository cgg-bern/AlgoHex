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

template<typename MeshT>
VFullConstraints &get_vertex_full_constraints(MeshT &mesh)
{
  auto prop = mesh.template request_mesh_property<VFullConstraints>("vertex full constraints");
  mesh.set_persistent(prop);
  return prop[OVM::MeshHandle(0)];
}

template<typename MeshT>
using VNormalConstraints = std::vector<std::pair<VH, typename MeshT::PointT>>;

template<typename MeshT>
VNormalConstraints<MeshT> &get_vertex_normal_constraints(MeshT &mesh)
{
  auto prop = mesh.template request_mesh_property<VNormalConstraints<MeshT>>("vertex normal constraints");
  mesh.set_persistent(prop);
  return prop[OVM::MeshHandle(0)];
}


template<typename MeshT>
using CPartialConstraints = std::vector<std::pair<CH, typename MeshT::PointT>>;

template<typename MeshT>
CPartialConstraints<MeshT> &get_cell_partial_constraints(MeshT &mesh)
{
  auto prop = mesh.template request_mesh_property<CPartialConstraints<MeshT>>("cell partial constraints");
  mesh.set_persistent(prop);
  return prop[OVM::MeshHandle(0)];
}

}

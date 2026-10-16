/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include "FieldConstraints.hh"

namespace AlgoHex
{

 VFullConstraints &get_vertex_full_constraints(TetMesh &mesh)
{
  auto prop = mesh.request_mesh_property<VFullConstraints>("vertex full constraints");
  mesh.set_persistent(prop);
  return prop[OVM::MeshHandle(0)];
}

 VNormalConstraints &get_vertex_normal_constraints(TetMesh &mesh)
{
  auto prop = mesh.request_mesh_property<VNormalConstraints>("vertex normal constraints");
  mesh.set_persistent(prop);
  return prop[OVM::MeshHandle(0)];
}

 CPartialConstraints &get_cell_partial_constraints(TetMesh &mesh)
{
  auto prop = mesh.request_mesh_property<CPartialConstraints>("cell partial constraints");
  mesh.set_persistent(prop);
  return prop[OVM::MeshHandle(0)];
}


}

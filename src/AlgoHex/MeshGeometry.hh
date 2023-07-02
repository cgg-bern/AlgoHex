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

template<class MeshT>
typename MeshT::PointT
get_vertex_normal(const MeshT &mesh, const VH _vh)
{
  typename MeshT::PointT vn(0, 0, 0);

  for (auto vohe_it = mesh.voh_iter(_vh); vohe_it.valid(); ++vohe_it)
    for (auto hehf_it = mesh.hehf_iter(*vohe_it); hehf_it.valid(); ++hehf_it)
      if (mesh.is_boundary(*hehf_it)) // TODO: also add faces that are internal interfaces
      {
        vn += mesh.normal(*hehf_it);
      }
  return vn.normalized();
}

template<typename PointT>
double
dihedral_angle(const PointT &pt0, const PointT &pt1, const PointT &n0, const PointT &n1)
{
  auto dprd = n0.dot(n1);
  dprd = std::min(1., std::max(-1., dprd));
  //should be normalized, but we need only the sign
  auto sign = n0.cross(n1).dot(pt1 - pt0);

  return sign >= 0 ? M_PI - std::acos(dprd) : M_PI + std::acos(dprd);
}

template<class MeshT>
double
dihedral_angle_in_cell(const MeshT &mesh, const HEH heh, const HFH hfh)
{
  assert(mesh.incident_cell(hfh).is_valid());

  VH vh0 = mesh.halfedge(heh).from_vertex();
  VH vh1 = mesh.halfedge(heh).to_vertex();
  auto pt0 = mesh.vertex(vh0);
  auto pt1 = mesh.vertex(vh1);

  auto n1 = mesh.normal(hfh);
  auto hf_adj = mesh.adjacent_halfface_in_cell(hfh, heh);
  auto n0 = mesh.normal(hf_adj);

  return dihedral_angle(pt0, pt1, n0, n1);
}

template<class MeshT>
double
dihedral_angle_from_halfface_to_halfface(const MeshT &mesh, const HEH heh, const HFH hfh_s, const HFH hfh_e)
{
  auto hfh_it = hfh_s;
  double angle = 0.;
  do
  {
    if (hfh_it == hfh_e)
      break;

    angle += dihedral_angle_in_cell(mesh, heh, hfh_it);

    hfh_it = mesh.adjacent_halfface_in_cell(hfh_it, heh);
    if (!hfh_it.is_valid())
    {
      std::cerr << "Error: adjacent halfface is invalid!" << std::endl;
      break;
    }
    hfh_it = mesh.opposite_halfface_handle(hfh_it);

  }
  while (hfh_it != hfh_s);

  return angle;
}

//dihedral angle of boundary halffaces
template<class MeshT>
double
dihedral_angle(const MeshT &mesh, const HEH _heh)
{
  if (!mesh.is_boundary(_heh))
    return -M_PI;

  HFH hfh1 = *mesh.hehf_iter(_heh),
          hfh0 = *mesh.hehf_iter(mesh.opposite_halfedge_handle(_heh));

  return dihedral_angle(mesh.vertex(mesh.halfedge(_heh).from_vertex()), mesh.vertex(mesh.halfedge(_heh).to_vertex()),
                        mesh.normal(hfh0), mesh.normal(hfh1));
}

template<typename Vec3>
double get_cell_volume(const std::vector<Vec3> &cell_points)
{
  Vec3 v01, v02, v03;
  v01 = cell_points[1] - cell_points[0];
  v02 = cell_points[2] - cell_points[0];
  v03 = cell_points[3] - cell_points[0];

  return v01 % v02 | v03 / 6.0;
}

template<class MeshT>
double get_cell_volume(const MeshT &mesh, const CH ch)
{
  auto cvhs = mesh.get_cell_vertices(ch);
  std::vector<typename MeshT::PointT> cell_pts;
  cell_pts.reserve(4);
  for (const auto vh: cvhs)
    cell_pts.push_back(mesh.vertex(vh));

  return get_cell_volume(cell_pts);
}

template<class MeshT>
double get_mesh_volume(const MeshT &_mesh)
{
  double V = 0.0;
  for (auto ch: _mesh.cells())
    V += get_cell_volume(_mesh, ch);

  return V;
}

template<class MeshT>
double tetmesh_volume(const MeshT &mesh)
{
  double vol(0.);
  for (const auto ch: mesh.cells())
    vol += get_cell_volume(mesh, ch);

  return vol;
}

template<class MeshT>
bool check_cells_volume(const MeshT &mesh)
{
  double minvol = std::numeric_limits<double>::max(), maxvol = std::numeric_limits<double>::min();
  CH minch(-1);

  for (const auto ch: mesh.cells())
  {
    double vol = get_cell_volume(mesh, ch);
    if (vol < minvol)
    {
      minvol = vol;
      minch = ch;
    }

    if (vol > maxvol)
      maxvol = vol;
  }

  std::cout << "Min/Max cell volume: " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << minvol
            << " / " << maxvol << std::setprecision(6) << ". Min volume cell: " << minch << std::endl;


  if (minvol <= 1e-16)
  {
    std::cout << "Error: (nearly) degenerate cell or flipped cell. Volume "
              << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << minvol
              << std::setprecision(6) << std::endl;

    return false;
  }

  return true;
}

template<class MeshT>
void check_mesh_dihedral_angle(const MeshT &mesh)
{
  double min_angle = std::numeric_limits<double>::max(), max_angle = std::numeric_limits<double>::min();
  CH maxch(-1);

  for (const auto eh: mesh.edges())
  {
    auto heh0 = mesh.halfedge_handle(eh, 0);

    for (auto hehf_it = mesh.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
    {
      auto ch = mesh.incident_cell(*hehf_it);
      if (ch.is_valid())
      {
        auto da = dihedral_angle_in_cell(mesh, heh0, *hehf_it);

        if (da > max_angle)
        {
          max_angle = da;
          maxch = ch;
        }

        if (da < min_angle)
          min_angle = da;
      }
    }
  }


  std::cout << "Min/Max cell dihedral angle: " << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
            << min_angle * 180. / M_PI
            << " / " << max_angle * 180. / M_PI << std::setprecision(6) << " Max dihedral angle cell: " << maxch
            << std::endl;
}


template<class MeshT>
bool check_mesh_quality(const MeshT &mesh)
{
  bool valid = check_cells_volume(mesh);
  if (!valid)
    return false;
  check_mesh_dihedral_angle(mesh);

  return true;
}

template<class MeshT>
void
scale_mesh(MeshT &_mesh, double _s)
{
  VIt v_it = _mesh.vertices_begin();
  VIt v_end = _mesh.vertices_end();
  for (; v_it != v_end; ++v_it)
  {
    _mesh.set_vertex(*v_it, _mesh.vertex(*v_it) * _s);
  }
}

template<class MeshT>
double
avg_local_edge_length(MeshT &_mesh, const VH _vh)
{
  double l = 0;
  int n = 0;
  VOHEIt vhe_it(_vh, &_mesh);
  for (; vhe_it.valid(); ++vhe_it)
  {
    l += _mesh.vector(*vhe_it).norm();
    ++n;
  }

  if (n == 0)
    return 0.0;
  else
    return l / double(n);
}


template<class VecT>
Vec3d
vec2eigen(const VecT &_v)
{
  return Vec3d(_v[0], _v[1], _v[2]);
}

template<class VecT>
OpenVolumeMesh::Vec3d
vec2ovm(const VecT &_v)
{
  return OpenVolumeMesh::Vec3d(_v[0], _v[1], _v[2]);
}

} // namespace AlgoHex

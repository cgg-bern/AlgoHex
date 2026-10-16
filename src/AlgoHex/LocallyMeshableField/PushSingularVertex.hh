/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <set>
#include "DetachAtInvalidSingularNode.hh"
#include <vector>
#include <AlgoHex/TypeDef.hh>
#include "FieldAngleCalculator.hh"

namespace AlgoHex
{

class ALGOHEX_EXPORT PushSingularVertex : public MeshProperties
{
public:
  using MeshProperties::valence_;
  using MeshProperties::trans_prop_;
  using MeshProperties::cell_quaternions_;

  using MeshProperties::feature_fprop_;
  using MeshProperties::feature_face_edge_;
  using MeshProperties::feature_face_vertex_;

  using MeshProperties::feature_edge_;
  using MeshProperties::feature_edge_vertex_;

  using MeshProperties::feature_node_;

  using MeshProperties::sgl_vt_;
  using MeshProperties::target_length_;


  using Point = TetMesh::PointT;


  PushSingularVertex(TetMesh &_mesh) : MeshProperties(_mesh),
                                      mesh_(_mesh),
                                      vo_(mesh_),
                                      sge_(mesh_),
                                      fis_(mesh_),
                                      fac_(mesh_),
                                      es_(mesh_)
  {
  }

public:
  bool push_singular_vertex(const VH _vh);

  bool is_move_inside_candidate(const VH _vh, bool _push_feature_vertex = true, bool _preprocess = false);

private:
  VH add_vertex_boundary_buffer(const VH _vh);

  void determine_matchings(const VH _vh_new, const VH _vh_orig, std::map<HEH, int> &_e_to_trans_idx,
                           std::map<HEH, int> &_hfh_to_trans_idx);

  int check_new_edge_index(const VH _vh, std::map<HEH, int> &_heh_to_trans_idx, std::map<HEH, int> &_hfh_to_trans_idx);

  bool get_original_boundary_halfedge_trans_idx(const VH _vh, std::map<HEH, int> &e_to_idx);

  int is_singular_edge_parallel_to_normal_direction(const HEH _heh, const HFH _hfh, const std::set<CH> _onering_chs);

  bool has_interior_singular_edge_tangent_to_surface(const VH _vh, const std::set<CH> &_or_chs,
                                                     const std::vector<HEH> &_itr_hehs);

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

private:
  TetMesh &mesh_;

  VertexOptimize vo_;

  AlgoHex::SingularGraphExtraction sge_;

  DetachAtInvalidSingularNode fis_;

  FieldAngleCalculator fac_;

  EdgeSplit es_;

  TransitionQuaternion tq_;
};
}

//=============================================================================
//=============================================================================

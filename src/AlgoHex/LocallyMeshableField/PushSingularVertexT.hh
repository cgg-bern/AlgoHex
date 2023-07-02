/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <set>
#include <vector>
#include <AlgoHex/TypeDef.hh>
#include "FieldAngleCalculatorT.hh"

namespace AlgoHex
{

template<class MeshT>
class PushSingularVertexT : public MeshPropertiesT<MeshT>
{
public:
  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;

  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;
  using MeshPropertiesT<MeshT>::feature_face_vertex_;

  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_edge_vertex_;

  using MeshPropertiesT<MeshT>::feature_node_;

  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::target_length_;


  using Point = typename MeshT::PointT;


  PushSingularVertexT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
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
  MeshT &mesh_;

  VertexOptimizeT<MeshT> vo_;

  AlgoHex::SingularGraphExtractionT<MeshT> sge_;

  DetachAtInvalidSingularNodeT<MeshT> fis_;

  FieldAngleCalculatorT<MeshT> fac_;

  EdgeSplitT<MeshT> es_;

  TransitionQuaternion tq_;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(PUSHSINGULARVERTEXT_C)
#define PUSHSINGULARVERTEXT_TEMPLATES

#include "PushSingularVertexT_impl.hh"

#endif
//=============================================================================


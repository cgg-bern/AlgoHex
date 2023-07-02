/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "ArcZippingT.hh"
#include <AlgoHex/SingularGraphExtractionT.hh>
#include "TetRemeshingT.hh"

namespace AlgoHex
{

template<class MeshT>
class EnsureEdgeMeshabilityT : public virtual MeshPropertiesT<MeshT>, public ArcZippingT<MeshT>
{
public:
  using ArcZippingT<MeshT>::valence_;
  using ArcZippingT<MeshT>::trans_prop_;
  using ArcZippingT<MeshT>::cell_quaternions_;

  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;
  using MeshPropertiesT<MeshT>::feature_face_vertex_;

  using ArcZippingT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_edge_vertex_;

  using MeshPropertiesT<MeshT>::feature_node_;

  using MeshPropertiesT<MeshT>::sgl_vt_;

  EnsureEdgeMeshabilityT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh), ArcZippingT<MeshT>(_mesh),
                                         mesh_(_mesh),
                                         sge_(mesh_),
                                         es_(mesh_),
                                         fs_(mesh_) {}


public:
  bool fix_edge_with_non_meshable_footprint(const EH _eh);

  bool fix_interior_compound_singular_edge(const HEH _heh, const bool _split = false);

  bool attempt_matching_adjustment_of_feature_face_sector(const HEH _heh, const HFH _hfh_s, const HFH _hfh_e,
                                                          const int _trans, const bool _allow_high_val_sge = false);

private:
  std::pair<double, int> best_transition_idx(const HEH _heh, const HFH _hfh, const std::vector<int> &_cd_vtrans);

  double frame_smoothness_energy_at_vertex(const VH _vh) const;

  double frame_smoothness_energy_at_edge(const HEH _heh) const;

  double frame_smoothness_energy_at_face(const FH _fh) const;

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

private:
  MeshT &mesh_;
  TransitionQuaternion tq_;

  AlgoHex::SingularGraphExtractionT<MeshT> sge_;
  EdgeSplitT<MeshT> es_;
  FaceSplitT<MeshT> fs_;

};
}



//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ENSUREEDGEMESHABILITYT_C)
#define ENSUREEDGEMESHABILITYT_TEMPLATES

#include "EnsureEdgeMeshabilityT_impl.hh"

#endif
//=============================================================================




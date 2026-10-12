/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once
#include <AlgoHex/TypeDef.hh>


#include "ArcZipping.hh"
#include <AlgoHex/SingularGraphExtraction.hh>
#include "TetRemeshing.hh"

namespace AlgoHex
{

class ALGOHEX_EXPORT EnsureEdgeMeshability : public virtual MeshProperties, public ArcZipping
{
public:
  using ArcZipping::valence_;
  using ArcZipping::trans_prop_;
  using ArcZipping::cell_quaternions_;

  using MeshProperties::feature_fprop_;
  using MeshProperties::feature_face_edge_;
  using MeshProperties::feature_face_vertex_;

  using ArcZipping::feature_edge_;
  using MeshProperties::feature_edge_vertex_;

  using MeshProperties::feature_node_;

  using MeshProperties::sgl_vt_;

  EnsureEdgeMeshability(TetMesh &_mesh);

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

  TetMesh &mesh_;
  TransitionQuaternion tq_;

  AlgoHex::SingularGraphExtraction sge_;
  EdgeSplit es_;
  FaceSplit fs_;

};
}



//=============================================================================
//=============================================================================

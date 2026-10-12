/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once


#include <unordered_set>
#include <AlgoHex/TypeDef.hh>

namespace AlgoHex
{
class ALGOHEX_EXPORT DuplicateOneRingMesh
{
public:
  using VecAxis = std::pair<Vec3d, int>;

  DuplicateOneRingMesh(TetMesh &_fmesh, TetMesh &_tmesh);

  ~DuplicateOneRingMesh() {}

  void copy_one_ring(const std::vector<VH> &_vhs);

  bool check_cell_singularity() const;

  VH onering_mesh_vertex_handle(const VH _fvh) const;

  HEH onering_mesh_halfedge_handle(const HEH _fheh) const;

  HFH original_halfface_handle(const HFH _hfh) const;

  EH original_edge_handle(const EH _eh) const;

  void get_halfedge_axis_in_cell_in_onering_mesh();

private:
  void copy_mesh(const std::vector<VH> &_vhs);

  void copy_properties();

  TetMesh &fmesh_;
  TetMesh &tmesh_;
  CP<Eigen::Quaterniond> cell_qts_1r_;
  HFP<int> trans_prop_1r_;
  EP<int> valence_1r_;
  FP<int> feature_fprop_1r_;
  EP<int> feature_eprop_1r_;

  VP<VH> tvh_fvh_;
  HEP<HEH> theh_fheh_;
  HFP<HFH> thfh_fhfh_;
  CP<CH> tch_fch_;
};
}
//=============================================================================
//=============================================================================

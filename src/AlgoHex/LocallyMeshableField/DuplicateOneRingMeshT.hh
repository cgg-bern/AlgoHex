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
template<class MeshT>
class DuplicateOneRingMeshT
{
public:
  using VecAxis = std::pair<Vec3d, int>;

  DuplicateOneRingMeshT(MeshT &_fmesh, MeshT &_tmesh) : fmesh_(_fmesh), tmesh_(_tmesh),
                                                        cell_qts_1r_(
                                                                tmesh_.template request_cell_property<Eigen::Quaterniond>(
                                                                        "FrameFieldQuaternions")),
                                                        trans_prop_1r_(tmesh_.template request_halfface_property<int>(
                                                                "HalffaceTransiton")),
                                                        valence_1r_(tmesh_.template request_edge_property<int>(
                                                                "edge_valance")),
                                                        feature_fprop_1r_(tmesh_.template request_face_property<int>(
                                                                "AlgoHex::FeatureFaces", false)),
                                                        feature_eprop_1r_(tmesh_.template request_edge_property<int>(
                                                                "AlgoHex::FeatureEdges", false)),
                                                        tvh_fvh_(tmesh_.template request_vertex_property<VH>(
                                                                "vertex handle in original mesh")),
                                                        theh_fheh_(tmesh_.template request_halfedge_property<HEH>(
                                                                "halfedge handle in original mesh")),
                                                        thfh_fhfh_(tmesh_.template request_halfface_property<HFH>(
                                                                "halfface handle in original mesh")),
                                                        tch_fch_(tmesh_.template request_cell_property<CH>(
                                                                "cell handle in original mesh"))
  {
    tmesh_.set_persistent(tch_fch_, true);
    tmesh_.set_persistent(theh_fheh_, true);
  }

  ~DuplicateOneRingMeshT() {}

  void copy_one_ring(const std::vector<VH> &_vhs)
  {
    copy_mesh(_vhs);
    copy_properties();
  }

  bool check_cell_singularity() const;

  VH onering_mesh_vertex_handle(const VH _fvh) const;

  HEH onering_mesh_halfedge_handle(const HEH _fheh) const;

  HFH original_halfface_handle(const HFH _hfh) const;

  EH original_edge_handle(const EH _eh) const;

  void get_halfedge_axis_in_cell_in_onering_mesh();

private:
  void copy_mesh(const std::vector<VH> &_vhs);

  void copy_properties();

private:
  MeshT &fmesh_;
  MeshT &tmesh_;
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
#if defined(INCLUDE_TEMPLATES) && !defined(DUPLICATEONERINGMESHT_C)
#define DUPLICATEONERINGMESHT_TEMPLATES

#include "DuplicateOneRingMeshT_impl.hh"

#endif
//=============================================================================




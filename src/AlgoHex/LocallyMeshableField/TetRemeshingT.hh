/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <AlgoHex/TransitionQuaternionEigen.hh>
#include <AlgoHex/SingularGraphExtractionT.hh>
#include "QuaternionsSmoothing.hh"
#include "MeshProperties.hh"
#include "RemeshingAssist.hh"
#include "CommonFuncs.hh"
#include "VertexOptimizeT.hh"

namespace AlgoHex
{
template<class MeshT>
class CellSplitT
{
public:
  CellSplitT(MeshT &_mesh) :
          mesh_(_mesh),
          target_length_(mesh_.template request_vertex_property<double>("target_length")),
          cell_quaternions_(mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions")) {}

  ~CellSplitT() {}

public:
  VH cell_split(const CH _ch);

  bool is_split_ok(const CH _ch) const;

private:
  void update_properties(const VH _vh_new, const Quaternion &_old_qt);


private:
  MeshT &mesh_;
  VP<double> target_length_;
  CP<Quaternion> cell_quaternions_;
};


template<class MeshT>
class FaceSplitT
{
public:
  FaceSplitT(MeshT &_mesh) :
          mesh_(_mesh),
          target_length_(mesh_.template request_vertex_property<double>("target_length")),
          trans_prop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
          cell_quaternions_(mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions")),
          feature_face_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureFaceVertices")),
          feature_face_edge_(mesh_.template request_edge_property<bool>("AlgoHex::FeatureFaceEdges")),
          feature_fprop_(mesh_.template request_face_property<int>("AlgoHex::FeatureFaces")),
          visited_fprop_(mesh_.template request_face_property<bool>("visited faces")) {}

  ~FaceSplitT() {}

public:
  VH face_split(const FH _fh);

  bool is_split_ok(const FH _fh) const;

private:
  void update_properties(std::vector<VH> &_orig_face_vhs, const VH _vh_new, const int _trans,
                         const int _inv_trans, const std::vector<Quaternion> &_old_qts,
                         const int _ff_id, const bool _is_visited_face);


private:
  MeshT &mesh_;
  VP<double> target_length_;

  HFP<int> trans_prop_;

  CP<Quaternion> cell_quaternions_;

  VP<bool> feature_face_vertex_;
  EP<bool> feature_face_edge_;
  FP<int> feature_fprop_;
  FP<bool> visited_fprop_;
};

template<class MeshT>
class EdgeSplitT
{
public:
  using Point = typename MeshT::PointT;

  EdgeSplitT(MeshT &_mesh) :
          mesh_(_mesh),
          vo_(mesh_),
          valence_(mesh_.template request_edge_property<int>("edge_valance")),
          trans_prop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
          cell_quaternions_(mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions")),
          sgl_vt_(mesh_.template request_vertex_property<int>("singular_vertex")),
          feature_edge_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureEdgeVertices")),
          feature_edge_(mesh_.template request_edge_property<int>("AlgoHex::FeatureEdges")),
          feature_face_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureFaceVertices")),
          feature_face_edge_(mesh_.template request_edge_property<bool>("AlgoHex::FeatureFaceEdges")),
          feature_fprop_(mesh_.template request_face_property<int>("AlgoHex::FeatureFaces")),
          target_length_(mesh_.template request_vertex_property<double>("target_length")),
          sg_edge_pairs_(mesh_.template request_edge_property<int>("separate edge pairs")),
          keep_on_feature_face_(mesh_.template request_edge_property<bool>("sge kept on feature face")),
          visited_fprop_(mesh_.template request_face_property<bool>("visited faces")) {}

  ~EdgeSplitT() {}

public:
  VH edge_split(const EH _eh, const bool _check_energy = false);

  bool is_split_ok(const EH _eh) const;

  void set_split_with_optimization() { opt_split_ = true; };

private:
  void
  update_properties(const VH _vh_f, const VH _vh_t, const VH _vh_new, const int _valence, const int _pair_id,
                    const bool _is_kept_on_ff,
                    const int _fe_id, const std::vector<int> &_feature_face_ids,
                    const std::vector<bool> &_is_visited_face,
                    const std::vector<std::vector<VH> > &_hfs_vhs, const std::vector<int> &_hfs_trans,
                    const std::vector<int> &_opp_hfs_trans, const std::vector<Quaternion> &_old_qts,
                    std::map<std::vector<VH>, int> &_hfvs_axis);

  std::vector<HFH>
  get_new_halffaces(const std::vector<VH> &_orig_vhs, const VH _vh_f, const VH _vh_t, const VH _vh_new) const;

  bool check_energy(const HEH _heh, Point &_split_pt) const;


private:
  MeshT &mesh_;
  AlgoHex::TransitionQuaternion tq_;

  VertexOptimizeT<MeshT> vo_;
  EP<int> valence_;
  HFP<int> trans_prop_;

  CP<Quaternion> cell_quaternions_;

  VP<int> sgl_vt_;
  VP<bool> feature_edge_vertex_;
  EP<int> feature_edge_;
  VP<bool> feature_face_vertex_;
  EP<bool> feature_face_edge_;
  FP<int> feature_fprop_;
  VP<double> target_length_;
  EP<int> sg_edge_pairs_;
  EP<bool> keep_on_feature_face_;
  FP<bool> visited_fprop_;

  bool opt_split_ = false;
};


template<class MeshT>
class EdgeCollapseT : public MeshPropertiesT<MeshT>
{
public:
  using Point = typename MeshT::PointT;

  EdgeCollapseT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                                mesh_(_mesh),
                                vo_(mesh_),
                                valence_(mesh_.template request_edge_property<int>("edge_valance")),
                                trans_prop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
                                cell_quaternions_(mesh_.template request_cell_property<Quaternion>(
                                        "FrameFieldQuaternions")),
                                sgl_vt_(mesh_.template request_vertex_property<int>("singular_vertex")),
                                feature_node_(
                                        mesh_.template request_vertex_property<int>("AlgoHex::FeatureVertices")),
                                feature_edge_vertex_(
                                        mesh_.template request_vertex_property<bool>("AlgoHex::FeatureEdgeVertices")),
                                feature_edge_(
                                        mesh_.template request_edge_property<int>("AlgoHex::FeatureEdges")),
                                feature_face_vertex_(mesh_.template request_vertex_property<bool>(
                                        "AlgoHex::FeatureFaceVertices")),
                                feature_face_edge_(
                                        mesh_.template request_edge_property<bool>("AlgoHex::FeatureFaceEdges")),
                                feature_fprop_(
                                        mesh_.template request_face_property<int>("AlgoHex::FeatureFaces")),
                                keep_on_feature_face_(
                                        mesh_.template request_edge_property<bool>("sge kept on feature face")) {}

  ~EdgeCollapseT() {}

  enum CollapseType
  {
    CollapseOK, DOFError, LinkError
  };

public:
  VH edge_collapse(const HEH _heh, const bool _check_energy = true);

  VH edge_collapse_allow_inverted(const HEH _heh);

  CollapseType is_collapse_ok(const HEH _heh) const;

  //topology check
  CollapseType is_topology_ok(const HEH _heh) const;

  std::vector<EH> get_blocking_edges(const EH _eh) const;

  void set_post_remeshing(const bool _post_remesh) { post_remesh_ = _post_remesh; }

  void set_collapse_with_optimization(const bool _opt_collapse = true) { opt_collapse_ = _opt_collapse; };

  EH opposite_edge_in_cell(const EH _eh, const CH _ch) const;

private:
  //get the new tets after edge collapse(just simulate, not real collapsing)
  void surviving_tets(const HEH _heh, std::vector<std::vector<VH> > &_v_vhs) const;

  void collect_axis_aligned_to_feature_face(const HEH _heh, std::map<std::vector<VH>, int> &_inc_hfv_to_axis,
                                            std::map<std::vector<VH>, int> &_ninc_hfv_to_axis);

  void update_properties(const EH _eh, const VH _vh_new, const std::vector<std::vector<VH> > &_hfs_vhs,
                         const std::vector<int> &_hfs_trans, std::vector<int> &_vals, const std::vector<VH> &_f_vhs,
                         const std::vector<VH> &_k_vhs, const std::vector<VH> &_fe_vhs, const std::vector<int> &_fe_ids,
                         const std::vector<std::vector<VH>> &_ff_vhs, const std::vector<int> &_ff_ids,
                         std::map<std::vector<VH>, int> &_inc_hfvhs_axis,
                         std::map<std::vector<VH>, int> &_ninc_hfvhs_axis);

  HFH opposite_halfface_to_vertex_in_cell(const VH _vh, const CH _ch) const;

  VH opposite_vertex_to_halfface_in_cell(const HFH _hfh) const;

  bool is_singular_node(const VH _vh) const;

  bool is_zipper_node(const VH _vh) const;

  bool is_zipper_node_collapsable(const VH _vh) const;

  bool is_boundary_edge_type_changed(const HEH _heh);

  bool has_boundary_valence_minus_one(const VH _vh) const;

  bool is_complex_singular_edge_after(const HEH _heh) const;

  int halfedge_transition_index(const HEH _heh, const HFH _hfh_s) const;

  bool check_energy(const HEH _heh, const std::vector<std::vector<VH>> &_new_cell_vhs, Point &_collapse_pt,
                    const bool _check_energy) const;

  int n_incident_special_edges(const VH _vh) const;

private:
  MeshT &mesh_;
  VertexOptimizeT<MeshT> vo_;

  AlgoHex::TransitionQuaternion tq_;
  EP<int> valence_;
  HFP<int> trans_prop_;

  CP<Quaternion> cell_quaternions_;


  VP<int> sgl_vt_;
  VP<int> feature_node_;
  VP<bool> feature_edge_vertex_;
  EP<int> feature_edge_;
  VP<bool> feature_face_vertex_;
  EP<bool> feature_face_edge_;
  FP<int> feature_fprop_;

  EP<bool> keep_on_feature_face_;

  //
  bool post_remesh_ = false;

  bool opt_collapse_ = false;
};


template<class MeshT>
class EdgeSwapT
{
public:
  using Point = typename MeshT::PointT;

  EdgeSwapT(MeshT &_mesh) :
          mesh_(_mesh),
          valence_(mesh_.template request_edge_property<int>("edge_valance")),
          trans_prop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
          cell_quaternions_(mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions")),
          sgl_vt_(mesh_.template request_vertex_property<int>("singular_vertex")),
          feature_edge_(mesh_.template request_edge_property<int>("AlgoHex::FeatureEdges")),
          feature_edge_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureEdgeVertices")),
          feature_face_edge_(mesh_.template request_edge_property<bool>("AlgoHex::FeatureFaceEdges")),
          feature_face_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureFaceVertices")),
          feature_fprop_(mesh_.template request_face_property<int>("AlgoHex::FeatureFaces")) {}

  ~EdgeSwapT() {}

public:
  std::set<EH> edge_swap(const EH _eh);

  bool is_swap_ok(const EH _eh) const;

  void set_post_remeshing(const bool _post_remesh) { post_remesh_ = _post_remesh; }

public:
  //make sure the transitions of the deleted halffaces are identity
  bool local_coordinate_transformation(const HEH _heh, const HFH _hf_start, const bool _is_bdy_e);

  std::vector<std::vector<VH>>
  collect_axis_aligned_to_feature_face(const HEH _heh, std::map<std::vector<VH>, int> &hfv_aligned_axis);

  bool has_confict_feature_face_constraints_in_cell(const std::vector<VH> &_cell_vhs,
                                                    std::map<std::vector<VH>, int> &hfv_aligned_axis);

  void collect_axis_aligned_to_feature(const HEH _heh, std::map<EH, int> &feh_aligned_axis);

  std::vector<std::vector<VH>> store_new_feature_faces_and_index(const HEH _heh, int &_ff_id) const;

  bool is_boundary_edge_type_changed(const EH _eh) const;

  void update_properties(const std::vector<FH> &_new_ft_fhs, const int _ff_id, const std::set<EH> &_new_ehs,
                         const std::set<CH> &_new_chs, const Quaternion &_qt,
                         std::map<std::vector<VH>, int> &_bdy_aligned_axis, std::map<EH, int> &feh_aligned_axis);

  //type: three tets cluster turns to two tets cluster after edge swap
  HFH swap_edge_32(const EH _eh, const std::vector<VH> &_equatorial_vertices, std::set<CH> &_new_chs);

  bool check_feature_face_constraints_32(const EH _eh, const std::vector<VH> &_equatorial_vertices,
                                         std::map<std::vector<VH>, int> &hfv_aligned_axis);

  //type: four tets cluster turns to four tets cluster after edge swap(2 different configurations)
  EH swap_edge_44(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh,
                  std::set<CH> &_new_chs);

  bool check_feature_face_constraints_44(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh,
                                         std::map<std::vector<VH>, int> &hfv_aligned_axis);

  //type: five tets cluster turns to six tets cluster after edge swap(5 different configurations)
  std::vector<EH> swap_edge_56(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh,
                               std::set<CH> &_new_chs);

  bool check_feature_face_constraints_56(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh,
                                         std::map<std::vector<VH>, int> &hfv_aligned_axis);

  //ccw
  std::vector<std::vector<VH>> get_equatorial_vertices(const EH _eh) const;

  //
  bool check_new_edges(const std::vector<VH> &_eq_vhs, const int _cnf) const;

  int find_the_best_configuration(const EH _eh_sw, const std::vector<VH> &_eq_vhs,
                                  std::vector<std::vector<Point> > &_old_cells_points) const;

  //get the new tets after edge swap
  std::vector<std::vector<VH> >
  surviving_tets_32(const EH _eh, const std::vector<VH> &_equatorial_vertices) const;

  std::vector<std::vector<VH> >
  surviving_tets_44(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh) const;

  std::vector<std::vector<VH> >
  surviving_tets_56(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh) const;

  HFH find_start_halfface(const HEH _heh, const bool _is_bdy_e) const;

private:
  MeshT &mesh_;
  EP<int> valence_;
  HFP<int> trans_prop_;

  CP<Quaternion> cell_quaternions_;

  VP<int> sgl_vt_;
  EP<int> feature_edge_;
  VP<bool> feature_edge_vertex_;
  EP<bool> feature_face_edge_;
  VP<bool> feature_face_vertex_;
  FP<int> feature_fprop_;

  AlgoHex::TransitionQuaternion tq_;

  //
  bool post_remesh_ = false;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(TETREMESHINGT_C)
#define TETREMESHINGT_TEMPLATES

#include "TetRemeshingT_impl.hh"

#endif

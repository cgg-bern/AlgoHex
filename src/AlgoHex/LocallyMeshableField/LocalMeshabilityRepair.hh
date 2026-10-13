/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once
#include <AlgoHex/TypeDef.hh>


#include <set>
#include "CellInfo.hh"

#include "EnsureEdgeMeshability.hh"
#include "DetachAtInvalidSingularNode.hh"
#include "FixConstrainedZipperNode.hh"
#include "FixZipperNode.hh"
#include "TransversalUnzipping.hh"


namespace AlgoHex
{

class ALGOHEX_EXPORT LocalMeshabilityRepair : virtual public MeshProperties, public ArcZipping
{
public:
  using MeshProperties::valence_;
  using MeshProperties::trans_prop_;
  using MeshProperties::cell_quaternions_;

  using MeshProperties::sgl_vt_;
  using MeshProperties::target_length_;

  using MeshProperties::feature_edge_;
  using MeshProperties::feature_edge_vertex_;
  using MeshProperties::feature_node_;

  using MeshProperties::feature_fprop_;
  using MeshProperties::feature_face_edge_;
  using MeshProperties::feature_face_vertex_;

  using Point = TetMesh::PointT;
  using DI = std::pair<double, int>;
  using DV = std::pair<double, VH>;
  using IHEH = std::pair<int, HEH>;


  using DHFINFO = std::pair<double, HALFFACEINFO>;

  using ICH = std::pair<int, CH>;

  using VecAxis = std::pair<Vec3d, int>;

  LocalMeshabilityRepair(TetMesh &_mesh);

  void set_first_stage_iters(const int _first_iters) { first_iters_ = _first_iters; }

  void set_second_stage_inner_iters(const int _second_inner_iters) { second_inner_iters_ = _second_inner_iters; }

  void set_dual_path_weight_wrt_field_direction(const double _dp_weight);

  void set_connect_other_zipper_node() { ftp_.set_connect_other_zipper_node(); }

  void set_unzip_zipper_nodes_from_longest_arc(
          const bool _longest_first) { fix_tps_on_longest_arc_first_ = _longest_first; };

  void set_max_search_dist_coeff_to_other_zipper_node(
          const double _max_tp_dist) { ftp_.set_max_search_dist_coeff_to_other_zipper_node(_max_tp_dist); }

  void set_peak_coefficient(const double _peak_coeff) { ftp_.set_peak_coefficient(_peak_coeff); }

  void set_segment_length_coefficient(const double _seg_coeff) { ftp_.set_segment_length_coefficient(_seg_coeff); }

  void set_absolute_segment_length_coefficient(const double _absl_seg_coeff);

  void set_merge_zipper_node() { merge_zipper_node_ = true; }

  void set_push_singular_circle() { push_singular_circle_ = true; }

  //_angle is in degree
  void set_angle_threshold(const double _angle) { ftp_.set_angle_threshold(_angle); }

  void preprocess(bool _check_alignment = true);

  int fix_first_stage(bool _check_alignment = true);

  int fix_second_stage(bool _check_alignment = true);

  void check_alignments();

  void split_for_local_meshability_test();

  //
  int unzip_zipper_nodes();

  std::set<VH> get_fixable_zipper_nodes(const bool _include_all) { return ftp_.get_fixable_zipper_nodes(_include_all); }

  int uniform_singular_arc_valence();

  void restore_singular_arc_valence();

  //fix singular vertex touching boundary (two or more arcs)
  int push_singular_vertices(bool _push_feature_vertex = true, bool _push_non_feature_boundary_singular_vertex = false);

  //
  void fix_edges_with_non_meshable_footprint();

  int fix_compound_singular_edges(const bool _first_stage);

  //fix invalid singular nodes
  int fix_interior_invalid_singular_nodes();

  int fix_boundary_invalid_singular_nodes();

  int fix_boundary_invalid_singular_nodes_start_from_bdy_sge();

  int fix_invalid_singular_vertices_on_feature_surface();

  int fix_invalid_singular_nodes_on_feature_surface(const bool _include_bdy = false);

  bool fix_vertex_on_tangential_singular_arc_at_feature_face(const VH _vh);

  int fix_fully_constrained_parabolic_sectors();

  bool fix_fully_constrained_parabolic_sectors_at_vertex(const VH _vh, const bool _print = false);

  int fix_constrained_zipper_nodes();

  int flatten_non_singular_feature_edges();

  bool flatten_non_singular_feature_edge(const HEH _heh);

  void fix_misalignment_at_feature_edges();

  void fix_misalignment_at_feature_edge(const EH _eh);

  void fix_misalignment_at_feature_face_sectors();

  void fix_misalignment_at_feature_face_sector(const EH _eh);

  void fix_misalignment_at_feature_face_sector(const HEH _heh, const HFH _hfh_s, const HFH _hfh_e);

  void remove_singular_triangles();

  int merge_zipper_nodes();

  //if an isolated zipper node is on a singular arc of two singular edges and they touch the boundary, zip it
  int zip_zipper_nodes();

  bool has_invalid_singular_edge() const;

  void non_meshable_vertex_numbers(int &_n_nm_vts, int &_n_tps, int &_n_nm_nodes) const;

  int sign(const int _axis) const { return _axis % 2 == 0 ? 1. : -1.; }

  void visualize_dominant_axis_on_singular_arc(const HEH _heh, std::map<CH, int> &_c_axis, TetMesh &_pm);

  int n_singular_edges_in_face(const FH _fh) const;

  bool is_removable_singular_triangle(const FH _fh) const;

  int n_incident_interior_singular_edges(const VH _vh) const;

  int matching_at_boundary_vertex(const VH _vh) const;

  void export_local_configuration(const VH vh, const std::string &filename) const;

private:
  bool is_alignment_consistent_at_feature_edge(const EH _eh) const;

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

  std::set<CH> get_onering_cells(const VH _vh) const;

private:
  TetMesh &mesh_;
  EP<int> sg_edge_pairs_;
  EP<bool> keep_on_feature_face_;

  AlgoHex::SingularGraphExtraction sge_;
  EdgeSplit es_;
  FaceSplit fs_;

  EdgeCollapse ec_;

  FieldAngleCalculator fac_;

  EnsureEdgeMeshability eem_;
  DetachAtInvalidSingularNode fis_;
  FixConstrainedZipperNode fctp_;
  FixZipperNode ftp_;
  TransversalUnzipping tu_;
  PushSingularVertex psv_;

  TransitionQuaternion tq_;

  double scale_ = 0.;

  bool merge_zipper_node_ = false;
  bool push_singular_circle_ = false;
  bool fix_tps_on_longest_arc_first_ = true;

  int iter_ = 0;
  int first_iters_ = 19;
  int second_inner_iters_ = 4;
};
}

//=============================================================================
//=============================================================================

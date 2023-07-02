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

#include "EnsureEdgeMeshabilityT.hh"
#include "DetachAtInvalidSingularNodeT.hh"
#include "FixConstrainedZipperNodeT.hh"
#include "FixZipperNodeT.hh"
#include "TransversalUnzippingT.hh"


#ifdef ALGOHEX_VERBOSE
#define ALGOHEX_DEBUG_ONLY(x) x
#else
#define ALGOHEX_DEBUG_ONLY(x)
#endif


namespace AlgoHex
{

template<class MeshT>
class LocalMeshabilityRepairT : virtual public MeshPropertiesT<MeshT>, public ArcZippingT<MeshT>
{
public:
  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;

  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::target_length_;

  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_edge_vertex_;
  using MeshPropertiesT<MeshT>::feature_node_;

  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;
  using MeshPropertiesT<MeshT>::feature_face_vertex_;

  using Point = typename MeshT::PointT;
  using DI = std::pair<double, int>;
  using DV = std::pair<double, VH>;
  using IHEH = std::pair<int, HEH>;


  using DHFINFO = std::pair<double, HALFFACEINFO>;

  using ICH = std::pair<int, CH>;

  using VecAxis = std::pair<Vec3d, int>;


  LocalMeshabilityRepairT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh), ArcZippingT<MeshT>(_mesh),
                                          mesh_(_mesh),
                                          sg_edge_pairs_(
                                                  mesh_.template request_edge_property<int>("separate edge pairs")),
                                          keep_on_feature_face_(
                                                  mesh_.template request_edge_property<bool>("sge kept on feature face",
                                                                                             false)),
                                          sge_(mesh_),
                                          es_(mesh_),
                                          fs_(mesh_),
                                          ec_(mesh_),
                                          fac_(mesh_),
                                          eem_(mesh_),
                                          fis_(mesh_),
                                          fctp_(mesh_),
                                          ftp_(mesh_),
                                          tu_(mesh_),
                                          psv_(mesh_)
  {
    //get mesh volume
    double vol = 0.;
    for (const auto &ch: mesh_.cells())
      vol += get_cell_volume(mesh_, ch);

    scale_ = 0.1 * std::pow(vol, 1. / 3.);
  }

  void set_first_stage_iters(const int _first_iters) { first_iters_ = _first_iters; }

  void set_second_stage_inner_iters(const int _second_inner_iters) { second_inner_iters_ = _second_inner_iters; }

  void set_dual_path_weight_wrt_field_direction(const double _dp_weight)
  {
    ftp_.set_dual_path_weight_wrt_field_direction(_dp_weight);
  };

  void set_connect_other_zipper_node() { ftp_.set_connect_other_zipper_node(); }

  void set_unzip_zipper_nodes_from_longest_arc(
          const bool _longest_first) { fix_tps_on_longest_arc_first_ = _longest_first; };

  void set_max_search_dist_coeff_to_other_zipper_node(
          const double _max_tp_dist) { ftp_.set_max_search_dist_coeff_to_other_zipper_node(_max_tp_dist); }

  void set_peak_coefficient(const double _peak_coeff) { ftp_.set_peak_coefficient(_peak_coeff); }

  void set_segment_length_coefficient(const double _seg_coeff) { ftp_.set_segment_length_coefficient(_seg_coeff); }

  void set_absolute_segment_length_coefficient(const double _absl_seg_coeff)
  {
    ftp_.set_absolute_segment_length_coefficient(_absl_seg_coeff);
  }

  void set_merge_zipper_node() { merge_zipper_node_ = true; }

  void set_push_singular_circle() { push_singular_circle_ = true; }

  //_angle is in degree
  void set_angle_threshold(const double _angle) { ftp_.set_angle_threshold(_angle); }

  void preprocess(bool _check_alignment = true);

  int fix_first_stage(bool _check_alignment = true);

  int fix_second_stage(bool _check_alignment = true);

  void check_alignments()
  {
    fac_.check_field_alignments_at_feature_faces();
    fac_.check_field_alignment_at_feature_edges();
  };

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

  void visualize_dominant_axis_on_singular_arc(const HEH _heh, std::map<CH, int> &_c_axis, MeshT &_pm);

  int n_singular_edges_in_face(const FH _fh) const;

  bool is_removable_singular_triangle(const FH _fh) const;

  int n_incident_interior_singular_edges(const VH _vh) const;

  int matching_at_boundary_vertex(const VH _vh) const;

  void export_local_configuration(const VH vh, const std::string &filename) const;

private:
  bool is_alignment_consistent_at_feature_edge(const EH _eh) const;

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

  std::set<CH> get_onering_cells(const VH _vh) const
  {
    std::set<CH> or_chs;
    for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
      or_chs.insert(*vc_it);

    return or_chs;
  }

private:
  MeshT &mesh_;
  EP<int> sg_edge_pairs_;
  EP<bool> keep_on_feature_face_;

  AlgoHex::SingularGraphExtractionT<MeshT> sge_;
  EdgeSplitT<MeshT> es_;
  FaceSplitT<MeshT> fs_;

  EdgeCollapseT<MeshT> ec_;

  FieldAngleCalculatorT<MeshT> fac_;

  EnsureEdgeMeshabilityT<MeshT> eem_;
  DetachAtInvalidSingularNodeT<MeshT> fis_;
  FixConstrainedZipperNodeT<MeshT> fctp_;
  FixZipperNodeT<MeshT> ftp_;
  TransversalUnzippingT<MeshT> tu_;
  PushSingularVertexT<MeshT> psv_;

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
#if defined(INCLUDE_TEMPLATES) && !defined(FIXLOCALLYNONMESHABLEVERTICEST_C)
#define FIXLOCALLYNONMESHABLEVERTICEST_TEMPLATES

#include "LocalMeshabilityRepairT_impl.hh"

#endif
//=============================================================================




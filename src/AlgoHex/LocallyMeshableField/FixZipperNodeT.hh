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
#include "TetRemeshingT.hh"
#include "ArcZippingT.hh"
#include "PushSingularVertexT.hh"
#include "MeshOptimizationT.hh"
#include "LocalMeshabilityCheckerWithFrames.hh"


namespace AlgoHex
{

template<class MeshT>
class FixZipperNodeT : public virtual MeshPropertiesT<MeshT>, public ArcZippingT<MeshT>
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
  using MeshPropertiesT<MeshT>::target_length_;


  using Point = typename MeshT::PointT;
  using DCINFO = std::pair<double, CELLINFO>;
  using DV = std::pair<double, VH>;
  using ICH = std::pair<int, CH>;
  using CHI = std::pair<CH, int>;
  using DCHI = std::pair<double, CHI>;

  using UVW = std::array<int, 3>;


  FixZipperNodeT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh), ArcZippingT<MeshT>(_mesh),
                                 mesh_(_mesh),
                                 sg_edge_pairs_(mesh_.template request_edge_property<int>("separate edge pairs")),
                                 keep_on_feature_face_(
                                         mesh_.template request_edge_property<bool>("sge kept on feature face", false)),
                                 sge_(mesh_),
                                 es_(mesh_),
                                 ec_(mesh_),
                                 fs_(mesh_),
                                 fac_(mesh_),
                                 fis_(mesh_),
                                 psv_(mesh_),
                                 fctp_(mesh_),
                                 lmcwf_(mesh_),
                                 mo_(mesh_),
                                 visited_fprop_(mesh_.template request_face_property<bool>("visited faces"))
  {
    mesh_.set_persistent(sg_edge_pairs_, true);
    mesh_.set_persistent(visited_fprop_, true);
  }

  bool unzip_zipper_node(const VH _vh, const HEH _heh_s, const std::set<VH> &_fixable_vhs, VH &_other_tp);

  bool merge_zipper_node(const VH _vh);

  bool zip_zipper_node(const VH _vh);

  std::set<VH> get_fixable_zipper_nodes(const bool _include_all = true);

  //find zipper node that is not on a circle
  bool is_zipper_node(const VH _vh, const bool _fix_on_circle = true) const;

  bool is_fixable_constrained_zipper_node(const VH _vh, HEH &_sg_heh);

  void set_dual_path_weight_wrt_field_direction(const double _w) { dpath_weight_ = _w; }

  void set_connect_other_zipper_node(const bool _connect = true) { connect_other_zipper_node_ = _connect; }

  void set_max_search_dist_coeff_to_other_zipper_node(const double _max_tp_dist) { max_tp_dist_ = _max_tp_dist; }

  void set_peak_coefficient(const double _peak_coeff) { peak_coeff_ = _peak_coeff; }

  void set_segment_length_coefficient(const double _seg_coeff) { seg_coeff_ = _seg_coeff; }

  void set_absolute_segment_length_coefficient(const double _absl_seg_coeff) { absl_seg_coeff_ = _absl_seg_coeff; }

//  void save_cells(const std::vector<CH> &_v_chs, int id) const;

  void set_angle_threshold(const double _angle) { angle_thr_ = _angle / 180. * M_PI; }

  void set_increase_angle_threshold(const double _angle) { inc_val_thr_ = std::cos(_angle / 180. * M_PI); }

private:
  CH find_start_cell(const VH _vh, HEH &_heh, HFH &_hfh) const;

  std::pair<double, CH> find_start_cell_at_halfedge(const VH _vh, const HEH _heh, HFH &_hfh) const;

  CELLINFO next_cell_info(const CELLINFO &ci_cur, const int _u, const int _v, const int _w, const HFH next_hfh,
                          const double _vw_weight, const double _length_scale, bool &_increase) const;

  std::vector<CH>
  shortest_dual_path_to_boundary(const HEH _heh_s, const HFH _hfh_s, const std::set<VH> &_fixable_vhs, VH &_other_tp,
                                 const double _w);

  std::vector<CH>
  shortest_dual_path_to_boundary(const CH _ch_s, int _axis, const std::set<VH> &_fixable_vhs, VH &_other_tp,
                                 const double _w);

  bool is_cell_found(const CH _ch, const int _trans, int _u) const;

  bool is_cell_found_wrt_zipper_node(const CELLINFO &ci_cur, const std::set<VH> &_fixable_tps, VH &_other_tp, int _u);

  std::vector<CH>
  shortest_dual_path_from_cell_to_target(const HEH _heh_s, const CH _ch_s, const std::set<CH> &_target_chs,
                                         const std::set<CH> &_dp_cells, const std::set<VH> &_excl_vh,
                                         const bool _spv_free) const;

  std::vector<CH>
  get_special_vertex_free_path(const std::vector<CH> &_input_path, const HEH _heh_s, const VH _other_tp);

  std::vector<CH> get_handle_free_dual_path(const std::vector<CH> &_input_path, const HEH _heh_s, const VH _other_tp);

  std::vector<EH>
  get_edge_path_on_cell_path(const std::vector<CH> &_dpath, const std::set<EH> &_excluded_ehs, const VH _vh_s,
                             const std::set<VH> &_target_vhs) const;

  std::vector<EH>
  get_edge_path_on_cell_path_wrt_field(const std::vector<CH> &_dpath, const std::set<EH> &_excluded_ehs, const VH _vh_s,
                                       const HEH _heh_s,
                                       const std::set<VH> &_target_vhs, const bool _opp_dir);

  std::vector<EH>
  shortest_path_v_to_target(const std::set<EH> &_all_eset, const std::set<EH> &_excl_eset, const VH _vh_f,
                            const std::set<VH> &_target_vhs, const VH _other_tp) const;

  std::set<EH> remedy_edge_path(const std::vector<CH> &_dpath, const std::vector<EH> &_raw_e_path) const;

  int split_for_second_edge_path(const std::vector<CH> &_dpath, const std::set<EH> &_raw_e_path, const EH _sg_eh,
                                 std::set<CH> &_dp_cells);

  std::map<CH, UVW> propagate_uvw_in_cells(const HEH _heh_s, const HFH _hfh_s, const std::set<CH> &_cells);

  std::map<CH, UVW> propagate_uvw_in_cells(const CH _ch_s, const UVW &_uvw, const std::set<CH> &_cells);

  //to vertex is turning point
  Vec3d get_difference_vector_at_zipper_node(const HEH _heh_s, const CH _ch_s, const std::vector<CH> &_input_path,
                                             std::map<CH, UVW> &_cell_axis);

  std::set<HFH> bounded_surface(const std::set<CH> &cells, const std::set<EH> &edges, const HEH _sg_heh) const;

  std::set<HFH>
  bounded_surface(const std::set<HFH> &_bound_hfhs, const std::set<EH> &_excl_edges, const HEH _sg_heh) const;

  bool is_zipper_node(const VH _vh, const HEH _heh0, const HEH _heh1);

public:
  int get_zipper_nodes_on_singular_arc(const std::vector<HEH> &_sg_hehs, std::vector<std::pair<VH, int>> &_tp_vhs);

  std::set<VH> filter_zipper_nodes(std::set<VH> &_vhs, const bool _include_all = true);

  std::set<VH> get_fixable_zipper_nodes_on_arc(const std::vector<HEH> &_sg_hehs,
                                               const std::vector<std::pair<VH, int>> &_tp_vi) const;

  std::vector<VH> get_non_meshable_zipper_nodes_on_singular_arc(const std::vector<HEH> &_sg_hehs,
                                                                const std::vector<std::pair<VH, int>> &_tp_vi);

  bool is_singular_arc_meshable(const std::vector<HEH> &_sg_hehs, const std::vector<std::pair<VH, int>> &_tp_vi);

private:
  void smooth_relevant_quaternions(const std::vector<VH> &_sg_vhs, std::map<CH, Quaternion> &_old_qtns);

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

  int n_unvisited_halffaces_at_manifold(const std::set<HFH> &_all_hfhs, const std::set<HFH> &_cut_hfhs, const HEH _heh,
                                        HFH &_hfh) const;

  bool is_valid_interior_singular_node_by_valence(const VH _vh) const
  {
    int nidx = this->node_index(_vh);
    if (nidx > 0 && nidx != 20 && nidx != 10)
      return true;

    return false;
  }

  bool is_other_special_vertex(const VH _vh, const std::set<VH> &_excl_vhs) const
  {
    if ((sgl_vt_[_vh] != 0 || feature_edge_vertex_[_vh]) && _excl_vhs.find(_vh) == _excl_vhs.end())
      return true;

    return false;
  }

  int n_bad_vertices_in_cell(const CH _ch, const std::set<VH> &_excl_vhs) const;

  int n_bad_vertices_in_face(const FH _fh, const std::set<VH> &_excl_vhs) const;

  int n_path_edges_in_face(const FH _fh, const std::set<EH> &_epath, EH &_not_path_eh) const;

  int n_path_vertices_in_face(const FH _fh, const std::set<VH> &_epath_vhs) const;

  bool is_singular_circle(const HEH _heh) const;

  bool is_feature_cell(const CH _ch) const;

  bool check_field_alignment(VH _vh, const std::set<EH> &_sg_ehs, const std::set<HFH> &_bound_hfhs,
                             const bool _save = false) const;


  bool is_mergable(const VH _vh, std::vector<HEH> &_v_hehs) const;

  bool intersect_with_triangle(const ACG::Geometry::PlaneT<double> &_plane, const Point &_p0, const Point &_p1,
                               const Point &_p2, const Point &_pa, const Point &_pb) const;

  Point compute_bary_coord(const Point &_p0, const Point &_p1, const Point &_p2, const Point &_px) const;


  bool is_zipable(const VH _vh, std::set<EH> &_sg_ehs) const;


private:
  MeshT &mesh_;
  EP<int> sg_edge_pairs_;
  EP<bool> keep_on_feature_face_;

  AlgoHex::SingularGraphExtractionT<MeshT> sge_;

  EdgeSplitT<MeshT> es_;
  EdgeCollapseT<MeshT> ec_;
  FaceSplitT<MeshT> fs_;

  FieldAngleCalculatorT<MeshT> fac_;
  DetachAtInvalidSingularNodeT<MeshT> fis_;
  PushSingularVertexT<MeshT> psv_;
  FixConstrainedZipperNodeT<MeshT> fctp_;
  LocalMeshabilityCheckerWithFrames<MeshT> lmcwf_;
  MeshOptimizationT<MeshT> mo_;


  FP<bool> visited_fprop_;

  TransitionQuaternion tq_;

  double inc_val_thr_ = -0.09;
  double dpath_weight_ = 100.;
  bool _edge_path_wrt_field = true;
  bool connect_other_zipper_node_ = false;
  double max_tp_dist_ = 1.1;

  double peak_coeff_ = 0.1;
  double seg_coeff_ = 0.4;
  double absl_seg_coeff_ = 5.;

  double angle_thr_ = 5.;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(FIXTURNINGPOINTT_C)
#define FIXTURNINGPOINTT_TEMPLATES

#include "FixZipperNodeT_impl.hh"

#endif
//=============================================================================


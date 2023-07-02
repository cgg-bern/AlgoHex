/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <AlgoHex/TypeDef.hh>
#include "TetRemeshingT.hh"
#include "FieldAngleCalculatorT.hh"

namespace AlgoHex
{
template<class MeshT>
class SeparableSingularArcFinderT : public MeshPropertiesT<MeshT>
{
public:
  using Point = typename MeshT::PointT;
  using HEPT = std::pair<Point, Point>;

  using Parameter = std::map<CH, std::map<VH, Point>>;
  using CellAxis = std::pair<CH, int>;



  //properties
  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;
  using MeshPropertiesT<MeshT>::feature_node_;
  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_edge_vertex_;
  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;
  using MeshPropertiesT<MeshT>::feature_face_vertex_;


  SeparableSingularArcFinderT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                                              mesh_(_mesh),
                                              keep_on_feature_face_(mesh_.template request_edge_property<bool>(
                                                      "sge kept on feature face")),
                                              fac_(mesh_),
                                              es_(mesh_) {}

  //with matchings
  std::set<HFH> find_anti_parallel_singular_edge_with_disk_path(const VH _vh, const std::set<CH> &_or_chs,
                                                                std::map<CH, int> &_cell_he_axis,
                                                                std::map<CH, int> &_cell_nm_axis,
                                                                std::set<HEH> &_he_pairs);

  std::set<HFH>
  find_anti_parallel_singular_edge_with_disk_path(const VH _vh, const HEH _heh_s, const std::vector<HEH> &_excl_hehs,
                                                  const std::set<CH> &_or_chs,
                                                  std::map<CH, int> &_cell_he_axis, std::map<CH, int> &_cell_nm_axis,
                                                  std::set<HEH> &_he_pairs);

  std::set<HFH>
  find_anti_parallel_singular_edge_with_disk_path(const VH _vh, const std::set<CH> &_or_chs, const HEH _heh_s,
                                                  const std::set<HEH> &_candidate_hehs,
                                                  std::map<CH, int> &_cell_he_axis, std::map<CH, int> &_cell_nm_axis,
                                                  const double _angle_threshold, std::set<HEH> &_he_pair);

  //
  std::set<HFH> find_boundary_with_disk_path(const VH _vh, const std::set<CH> &_or_chs, const std::set<HEH> &_fs_hehs,
                                             std::map<CH, int> &_cell_nm_axis, std::set<HEH> &_he_pairs,
                                             const bool _is_anti_prl, const bool _start_from_bdy);

  std::set<HFH> find_boundary_with_disk_path(const VH _vh, const std::set<CH> &_or_chs,
                                             const HEH _heh_s, const std::set<CH> &_candidate_chs,
                                             std::map<CH, int> &_cell_nm_axis, std::set<HEH> &_he_pair,
                                             const bool _is_anti_prl);

  //
  std::set<HFH> find_zipper_node_edge_pair_with_disk_path(const VH _vh, const std::set<CH> &_or_chs,
                                                          std::map<CH, int> &_cell_he_axis,
                                                          std::map<CH, int> &_cell_nm_axis, std::set<HEH> &_he_pairs);

  std::set<HFH>
  find_zipper_node_edge_pair_with_disk_path(const VH _vh, const HEH _heh_s, const std::vector<HEH> &_excl_hehs,
                                            const std::set<CH> &_or_chs,
                                            std::map<CH, int> &_cell_he_axis, std::map<CH, int> &_cell_nm_axis,
                                            std::set<HEH> &_he_pairs);

  std::set<HFH> find_zipper_node_edge_pair_with_disk_path(const VH _vh, const std::set<CH> &_or_chs, const HEH _heh_s,
                                                          const std::set<HEH> &_candidate_hehs,
                                                          std::map<CH, int> &_cell_he_axis,
                                                          std::map<CH, int> &_cell_nm_axis,
                                                          const double _angle_threshold_tp, std::set<HEH> &_he_pair);

  //
  std::set<HFH> find_disk_to_feature_face(const VH _vh, const std::set<CH> &_or_chs,
                                          std::map<CH, int> &_cell_nm_axis, std::set<HEH> &_he_pairs,
                                          const bool _is_anti_prl);

  std::set<HFH> find_disk_to_feature_face(const VH _vh, const std::set<CH> &_or_chs,
                                          const HEH _heh_s, std::map<CH, int> &_cell_nm_axis, std::set<HEH> &_he_pair,
                                          const bool _is_anti_prl);

  //allow
  HEH find_anti_parallel_singular_edge_with_dual_path(const HEH _heh, const CH _ch_s, std::map<CH, int> &_cell_he_axis,
                                                      std::map<CH, int> &_cell_nm_axis, const int _target_axis,
                                                      const std::set<CH> &_onering_chs,
                                                      const std::set<HEH> &_candidate_hehs,
                                                      const double _angle_threshold, std::vector<CH> &_dpath);

  HEH is_anti_parallel_found(const CellAxis &_ca, const VH _vhf, const HEH _heh0, std::map<CH, int> &_cell_he_axis);

  CellAxis next_cell_axis(const CellAxis &_ca_cur, const HFH _hf_opp, double &_inc_val, bool &_increasing) const;

  //
  CH
  find_boundary_with_dual_path(const HEH _heh, const CH _ch_s, std::map<CH, int> &_cell_nm_axis, const int _target_axis,
                               const std::set<CH> &_onering_chs, const std::set<CH> &_candidate_chs,
                               std::vector<CH> &_dpath, const bool _is_anti_prl);

  bool is_boundary_found(const CellAxis &_ca, std::map<CH, int> &_cell_nm_axis, const bool _is_anti_prl);


  //
  CH find_dual_path_to_feature_face(const HEH _heh, const CH _ch_s, std::map<CH, int> &_cell_nm_axis,
                                    const int _target_axis,
                                    const std::set<CH> &_onering_chs, std::vector<CH> &_dpath, const bool _is_anti_prl);

  bool is_target_feature_face_found(const CellAxis &_ca, const FH _inc_fh, std::map<CH, int> &_cell_nm_axis,
                                    const bool _is_anti_prl);

  //
  HEH find_zipper_node_edge_pair_with_dual_path(const HEH _heh, const CH _ch_s, std::map<CH, int> &_cell_he_axis,
                                                std::map<CH, int> &_cell_nm_axis,
                                                const int _target_axis, const std::set<CH> &_onering_chs,
                                                const std::set<HEH> &_candidate_hehs, const double _angle_threshold_tp,
                                                std::vector<CH> &_dpath);

  HEH
  is_zipper_node_edge_pair_found(const CellAxis &_ca, const VH _vhf, const HEH _heh0, std::map<CH, int> &_cell_he_axis,
                                 std::map<CH, int> &_cell_nm_axis);

  EH opposite_edge_in_cell(const EH _eh, const CH _ch) const;

  std::vector<CH>
  get_singular_vertex_free_path(const VH _vh, const std::vector<CH> &_dpath, const std::vector<HEH> &_sghehs,
                                bool _to_boundary = false);

  std::vector<CH> get_handle_free_dual_path(const std::vector<CH> &_input_path, const std::vector<HEH> &_sghehs);

  std::vector<CH>
  shortest_dual_path_from_cell_to_target(const CH _ch_s, const std::set<CH> &_target_chs, const std::set<CH> &_dp_cells,
                                         const std::set<VH> &_excld_vhs, const bool _no_sgv = true) const;

  std::vector<EH> bad_edges_in_cell2(const CH _ch, const std::set<VH> &_excld_sg_vhs) const;

  std::vector<EH> bad_edges_in_cell(const CH _ch, const std::set<VH> &_excld_sg_vhs) const;

  int n_bad_vertices_in_cell(const CH _ch, const std::set<VH> &_excld_sg_vhs) const;

  std::set<HFH> bounded_surface(const std::vector<CH> &_cells, const std::set<HEH> &_sg_hehs,
                                const bool _excld_feature_face_edge = false) const;

  int n_unvisited_halffaces_at_manifold(const std::set<HFH> &_all_hfhs, const std::set<HFH> &_cut_hfhs, const HEH _heh,
                                        HFH &_hfh) const;

  template<typename PT>
  Vec3d eigen_point(const PT &_pt) const { return Vec3d(_pt[0], _pt[1], _pt[2]); }

  Point ovm_point(const Vec3d &_pt) const { return Point(_pt[0], _pt[1], _pt[2]); }

  int check_parameterization(Parameter &_pm);

  void parametric_tet(Parameter &_pm, const CH _ch, Eigen::Matrix<double, 3, 4> &_P);

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

  HEH singular_halfedge_in_cell(const VH _vhf, const CH _ch, const HEH _he_excl) const;

  std::set<HEH> get_candidate_halfedges(const std::set<HEH> &_hehs, const HEH _heh, int _valence, bool _same_valence);

  double get_angle_threshold(const VH _vh, const HEH _heh, const std::set<CH> &_or_chs) const;

  double get_angle_threshold_zipper_node(const VH _vh, const HEH _heh, const std::set<CH> &_or_chs) const;

  //check feature face sector w.r.t normals
  bool is_valid_sector_after_fix(const VH _vh, const HEH _ff_sgheh, const HEH _nff_sgheh);

  //detach wrt boundary
  std::set<CH>
  get_candidate_cells(const VH _vh, const HEH _sg_heh, const std::set<HEH> &_ft_hehs, const bool _is_anti_prl);

  std::set<CH> get_blocking_cells(const HEH _heh_s, std::set<FH> &_blk_fhs, const bool _is_anti_prl = true);

  std::vector<HEH> sort_starting_halfedges_for_antiparallel(const VH _vh);

  std::pair<HEH, HEH> most_parallel_halfedge_pair_in_set(const std::set<HEH> &_all_hehs, std::map<CH, int>& _cell_region_id);

  bool is_feature_cell(const CH _ch) const;

private:
  MeshT &mesh_;
  EP<bool> keep_on_feature_face_;

  FieldAngleCalculatorT<MeshT> fac_;

  EdgeSplitT<MeshT> es_;

  double triangle_elen_ = 0.;

  AlgoHex::TransitionQuaternion tq_;

public:
  bool find_feature_face_sges_only_ = false;
  bool same_sector_sges_only_ = false;
  double angle_threshold_ = -0.174; // -10 degree
  double angle_threshold_tp_ = 0.707; //45 degree
  bool enable_blocking_ = true;
  bool enable_detach_kept_sge_ = false;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(SEPARABLESINGULARARCFINDERT_C)
#define SEPARABLESINGULARARCFINDERT_TEMPLATES

#include "SeparableSingularArcFinderT_impl.hh"

#endif
//=============================================================================




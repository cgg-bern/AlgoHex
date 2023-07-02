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
#include "ArcZippingT.hh"
#include "TetRemeshingT.hh"
#include "TransversalUnzippingT.hh"


namespace AlgoHex
{

template<class MeshT>
class FixConstrainedZipperNodeT : public virtual MeshPropertiesT<MeshT>, public ArcZippingT<MeshT>
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

  using Point = typename MeshT::PointT;
  using DHFINFO = std::pair<double, HALFFACEINFO>;


  struct GDINFO
  {
  public:
    GDINFO(const HEH _heh, const HFH _hfh, const int _normal_axis, const int _search_axis, const int _rt_axis,
           const bool _ccw = false) :
            heh(_heh), hfh(_hfh), normal_axis(_normal_axis), search_axis(_search_axis), rt_axis(_rt_axis), ccw(_ccw) {}

    friend std::ostream &operator<<(std::ostream &out, const GDINFO &gd)
    {
      return out << "halfedge: " << gd.heh << ", halfface: " << gd.hfh << " normal axis: " << gd.normal_axis
                 << " search axis: " << gd.search_axis << " rotation axis: " << gd.rt_axis << std::endl;
    }

    HEH heh;
    HFH hfh;
    int normal_axis;
    int search_axis;
    int rt_axis;
    bool ccw;
  };


  FixConstrainedZipperNodeT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh), ArcZippingT<MeshT>(_mesh),
                                            mesh_(_mesh),
                                            sg_edge_pairs_(
                                                    mesh_.template request_edge_property<int>("separate edge pairs")),
                                            keep_on_feature_face_(mesh_.template request_edge_property<bool>(
                                                    "sge kept on feature face", false)),
                                            sge_(mesh_),
                                            fac_(mesh_),
                                            es_(mesh_),
                                            fs_(mesh_),
                                            fis_(mesh_),
                                            tu_(mesh_)
  {
    mesh_.set_persistent(keep_on_feature_face_, true);
  }

public:
  void set_second_stage() { second_stage_ = true; };

  int fix_constrained_zipper_node(const VH _vh);

  bool fix_constrained_zipper_node1(const HEH _heh, const VH _vh_target, const int _fix_type,
                                    const bool _keep_on_ff = false);

  bool fix_constrained_zipper_node2(const HEH _heh, const int _fix_type, const bool _keep_on_ff);

  std::vector<GDINFO>
  search_for_guides_along_feature_arc(const HEH _sg_heh, const VH _target_vh, int &_end_type, HEH &_end_ft_he,
                                      HFH &_end_ft_hf, const bool _open_sec, const bool _same_dir, const bool _cross);

  std::vector<GDINFO> search_for_guides_along_singular_arc(const HEH _sg_heh, int &_end_type);

  std::vector<GDINFO>
  search_for_guides_on_feature_surface(const HEH _sg_heh, int &_end_type, HEH &_end_ft_he, HFH &_end_ft_hf,
                                       const bool _open_sec, const bool _same_dir);

  std::vector<GDINFO>
  search_for_guides_on_feature_surface(const HEH _sg_heh, const CH _ch, const VH _vh, const int _search_ax,
                                       int &_end_type, HEH &_end_ft_he, HFH &_end_ft_hf, const bool _open_sec);

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

  int fixable_constrained_zipper_node_type(const VH _vh, HEH &_heh);

  int fixable_constrained_zipper_node_type(const HEH _heh);

  int singular_arc_tangentially_touching_feature_face_type(const VH _vh);

  bool is_singular_arc_crossing_feature_face(const VH _vh) const;

  bool is_singularity_crossing_feature_face(const VH _vh) const;

  bool is_separable_at_vertex(const VH _vh, std::set<HEH> &_ft_hehs);

  bool is_valid_feature_face_sector_after(const HEH _heh, const HFH _hfh_s, const HFH _hfh_e, int _trans);


private:
  void get_special_sectors_at_special_edge_vertex(const VH vh, std::vector<std::vector<HEH>> &v_sec_hehs,
                                                  std::vector<std::set<FH>> &v_sec_fhs) const;

  //end type:
  VH unzip_arc_with_guiding_info(const HEH _heh, const std::vector<GDINFO> &_guides, const int _end_type,
                                 const HEH _end_ft_heh, const HFH _end_ft_hfh, VH &_vh_hint);

  VH extend_with_zipper_node(const HEH _sg_heh);

  std::vector<GDINFO> get_guides_cw(const GDINFO &_gd, const VH _target_vh, const bool _open_sec,
                                    int &_end_type, HEH &_end_ft_he, HFH &_end_ft_hf);

  std::vector<GDINFO> get_guides_ccw(const GDINFO &_gd, const VH _target_vh, const bool _open_sec,
                                     int &_end_type, HEH &_end_ft_he, HFH &_end_ft_hf);

  GDINFO get_next_parallel_special_halfedge_halfface_on_feature_surface_cw(const GDINFO &_gd);

  int find_orthogonal_sector(const GDINFO &_gd, HEH &_heh, HFH &_hfh, int &_sec_angle, const bool _open_sec) const;

  void get_orthogonal_sectors(const HEH _heh, std::vector<std::vector<HEH>> &_v_cvx_orth_sec_hehs,
                              std::vector<std::set<FH>> &_v_cvx_orth_sec_fhs,
                              std::vector<std::vector<HEH>> &_v_ccv_orth_sec_hehs,
                              std::vector<std::set<FH>> &_v_ccv_orth_sec_fhs) const;

  bool has_orthogonal_sector_cw(const GDINFO &_gd, HEH &_heh, HFH &_hfh, int &_angle) const;

  std::tuple<HEH, HFH, int> next_outgoing_special_hehf_with_trans_cw(const HEH _heh, const HFH _hfh) const;

  std::tuple<HEH, HFH, int> next_incoming_special_hehf_with_trans_ccw(const HEH _heh, const HFH _hfh) const;

  std::pair<HEH, HFH> next_flat_sector_with_trans(const HEH _heh, const HFH _hfh, const bool _cw) const;

  GDINFO get_next_parallel_special_halfedge_halfface_on_feature_surface_ccw(const GDINFO &_gd);

  bool has_orthogonal_sector_ccw(const GDINFO &_gd, HEH &_heh, HFH &_hfh, int &_angle) const;

  bool has_special_edge_in_direction_in_sector(const VH _vh, const HFH _hfh_s, const int _ax) const;

  std::tuple<double, HFH, int>
  find_furthest_halfface_in_direction_in_sector(const VH _vh, const HFH _hfh_s, const int _ax) const;

  bool
  has_convexly_parallel_sector_in_direction(const VH _vh, const CH _ch_s, const int _dir, const bool _check_spe) const;

  std::tuple<double, HFH, int>
  furthest_parallel_sector_in_direction(const VH _vh, const CH _ch_s, const int _dir, const bool _check_spe) const;

  bool has_detachable_singular_edge_at_vertex(const CH _ch_s, int _rt_axis, const VH _vh_query);

  bool has_detachable_non_feature_face_singular_edge_at_vertex(const CH _ch_s, int _rt_axis, const VH _vh_query);

private:
  std::vector<FH> shortest_face_path_on_feature_surface(const GDINFO &_gd, VH &_other_vh, const double _weight);

  std::vector<FH>
  shortest_dual_path_from_source_face_to_target_face(const FH _fh_s, const FH _fh_t, const std::set<FH> &_dp_fhs,
                                                     const std::set<VH> &_included_vhs) const;

  std::vector<FH>
  get_special_vertex_free_path_on_feature_surface(const std::vector<FH> &_input_path, const VH _vh_s, const VH _vh_t);

  HALFFACEINFO next_halfface_info(const HALFFACEINFO &hi_cur, const int _u, const int _v, const int _w,
                                  const HEH _next_heh, const HFH _next_hfh, const double _vw_weight,
                                  const double _length_scale, double &_incr_val) const;

  bool is_face_found(const HALFFACEINFO &hi_cur, int _u) const;

  bool is_feature_edge_orthogonal_to_search_direction(const HALFFACEINFO &hi_cur, const EH _eh, int _u) const;

  HFH next_halfface_on_feature_surface(const HEH _heh, const HFH _hfh_s) const;

  EH common_edge(const FH _fh0, const FH _fh1) const;

  VH split_edge_with_face_property_update(const EH _eh, std::set<FH> &_dp_fhs);

  int split_edges_with_face_property_update(std::set<FH> &_dp_fhs, const std::set<EH> &_ehs, std::set<FH> &_target_fhs);

  bool is_other_feature_edge_vertex(const VH _vh, const std::set<VH> &_excluded_vhs) const
  {
    if (feature_edge_vertex_[_vh])
    {
      if (_excluded_vhs.find(_vh) != _excluded_vhs.end())
        return false;
    }
    else
      return false;

    return true;
  }

  bool is_other_singular_vertex(const VH _vh, const std::set<VH> &_excluded_vhs) const
  {
    if (sgl_vt_[_vh] != 0)
    {
      if (_excluded_vhs.find(_vh) != _excluded_vhs.end())
        return false;
    }
    else
      return false;

    return true;
  }

  int n_bad_vertices_in_face(const FH _fh, const std::set<VH> &_excl_vhs) const
  {
    int n = 0;
    for (auto fv_it = mesh_.fv_iter(_fh); fv_it.valid(); ++fv_it)
      if ((sgl_vt_[*fv_it] != 0 || feature_edge_vertex_[*fv_it]) && _excl_vhs.find(*fv_it) == _excl_vhs.end())
        n++;


    return n;
  }

private:
  MeshT &mesh_;
  EP<int> sg_edge_pairs_;
  EP<bool> keep_on_feature_face_;

  SingularGraphExtractionT<MeshT> sge_;
  FieldAngleCalculatorT<MeshT> fac_;

  EdgeSplitT<MeshT> es_;
  FaceSplitT<MeshT> fs_;
  DetachAtInvalidSingularNodeT<MeshT> fis_;
  TransversalUnzippingT<MeshT> tu_;

  TransitionQuaternion tq_;

  bool second_stage_ = false;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(FIXCONSTRAINEDZIPPERNODET_C)
#define FIXCONSTRAINEDZIPPERNODET_TEMPLATES

#include "FixConstrainedZipperNodeT_impl.hh"

#endif
//=============================================================================

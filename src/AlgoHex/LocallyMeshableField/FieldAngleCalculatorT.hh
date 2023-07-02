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
#include <AlgoHex/TransitionQuaternionEigen.hh>
#include "CommonFuncs.hh"
#include "MeshProperties.hh"

namespace AlgoHex
{

template<class MeshT>
class FieldAngleCalculatorT : public MeshPropertiesT<MeshT>
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


  FieldAngleCalculatorT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                                        mesh_(_mesh)
  {
  }

public:
  struct HALFFACEINFO2
  {
  public:
    explicit HALFFACEINFO2(const HFH _hfh, const double _u_dist, const int _trans) : hfh(_hfh), u_dist(_u_dist),
                                                                                     trans(_trans) {}

    friend std::ostream &operator<<(std::ostream &out, const HALFFACEINFO2 &ci)
    {
      return out << "halfface: " << ci.hfh << ", u dist: " << ci.u_dist << " trans prod " << ci.trans << std::endl;
    }

    HFH hfh;
    double u_dist;
    //matching from the start cell to the current cell
    int trans;
  };

public:
//  void feature_face_angles_at_feature_edge(const HEH _heh) const;

  void print_feature_sector_angles(const VH _vh) const;

  bool check_field_alignments_at_feature_faces() const;

  bool check_field_alignment_at_feature_edges() const;

  bool check_field_alignment_at_feature_edge(const EH _eh) const;

  bool check_feature_face_sectors_at_regular_edge(const EH _eh) const;

  static double
  average_field_edge_angle(const MeshT &_mesh, const TransitionQuaternion &_tq, const CP<Quaternion> &_cell_quaternions,
                           const HFP<int> &_trans_prop, const EP<int> &_valence,
                           const EP<int> &_feature_edge, const EH _eh);

//
  static int feature_face_sector_angle_at_halfedge(const MeshT &mesh, const TransitionQuaternion &tq,
                                                   const CP<Quaternion> &cell_quaternions,
                                                   const HFP<int> &trans_prop, const HEH heh, const HFH hfh_s,
                                                   const HFH hfh_e);

  //feature edge sector w.r.t. edge direction
  static int field_angle_status_at_feature_sector(const MeshT &mesh, const TransitionQuaternion &tq,
                                                  const CP<Quaternion> &cell_quaternions,
                                                  const HFP<int> &trans_prop, const EP<int> &valence,
                                                  const std::vector<HEH> &_sector_hehs,
                                                  const std::set<FH> &_sector_fhs);

  static double sector_angle(const MeshT &mesh, const std::vector<HEH> &hehs);

//
  int field_angle_of_feature_edge_sector(const HEH _heh, const HFH _hfh, const bool _cw = true) const;

  std::tuple<HEH, HFH, int>
  next_outgoing_special_hehf_with_trans_cw(const HEH _heh, const HFH _hfh, std::vector<HEH> &_sec_hehs) const;

  std::tuple<HEH, HFH, int>
  next_incoming_special_hehf_with_trans_ccw(const HEH _heh, const HFH _hfh, std::vector<HEH> &_sec_hehs) const;

  bool is_axis_in_cell_parallel_to_any_surface(const VH _vh, const int _axis, const CH _ch);

  bool has_any_convexly_orthogonal_face(const VH _vh, const FH _fh) const;

  bool has_any_convexly_orthogonal_face_cw(const VH _vh, const HFH _hfh) const;

  bool has_any_convexly_orthogonal_face_ccw(const VH _vh, const HFH _hfh) const;

  bool is_axis_in_cell_parallel_to_any_surface_convexly(const VH _vh, const int _axis, const CH _ch);

  bool has_special_edge_in_direction(const VH _vh, const CH _ch_s, const int _ax_s);

  bool has_singular_edge_in_direction(const HEH _heh_query, const CH _ch_query, int _axis_query, const VH _vh_query,
                                      CH &_ch_found, HFH &_hfh_found, HEH &_heh_found, int &_axis_found) const;

  bool
  has_singular_edge_in_direction_longest(const HEH _heh_query, const CH _ch_query, int _axis_query, const VH _vh_query,
                                         CH &_ch_found, HFH &_hfh_found, HEH &_heh_found, int &_axis_found) const;

  bool is_feature_edge_orthogonal_to_any_sector(const HEH _ft_heh);

  bool is_singular_edge_tangent_to_any_surface(const std::set<CH> &_or_chs, const HEH _sg_heh) const;

  bool is_singular_edge_orthogonal_to_any_surface(const std::set<CH> &_or_chs, const HEH _sg_heh);

  bool is_singular_edge_tangent_to_all_surface(const HEH _sg_heh) const;

  int angle_of_axis_in_cell_and_normal_direction(const int _e_axis, const CH _sg_ch, const HFH _ft_hfh,
                                                 const std::set<CH> _onering_chs) const;

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

  std::vector<int> rotation_axes_expressed_in_common_tet(const VH _vh, const std::vector<HEH> &_sg_hehs);

public:
  int axis_in_chs_expressed_in_cht(const CH _ch_s, const CH _ch_t, int _axis,
                                   const std::set<CH> _onering_chs) const;

  int axis_after_transition_along_dual_path(const int _axis, const std::vector<HFH> _dpath,
                                            const bool _reversed = false) const;

  std::vector<HFH> get_dual_path_between_cells_in_onering(const CH _ch_s, const CH _ch_t,
                                                          const std::set<CH> &_onering_chs) const;

  std::map<CH, HFH> spanning_tree_in_onering(const CH _ch_seed, const std::set<CH> &_onering_chs) const;

  std::vector<HFH> get_dual_path_to_cell_in_onering(const CH _ch_t, std::map<CH, HFH> &_pre_hf);

private:
  HALFFACEINFO2 next_halfface_info(const HALFFACEINFO2 &hi_cur, const HEH _he_s, const int _u) const;

private:
  MeshT &mesh_;

  TransitionQuaternion tq_;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(FIELDANGLECALCULATORT_C)
#define FIELDANGLECALCULATORT_TEMPLATES

#include "FieldAngleCalculatorT_impl.hh"

#endif
//=============================================================================





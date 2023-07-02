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
#include <AlgoHex/SingularGraphExtractionT.hh>
#include "ArcZippingT.hh"
#include "SeparableSingularArcFinderT.hh"
//#include "LocalMeshabilityCheckerWithFrames.hh"
#include "OneringParameterizationChecker.hh"

namespace AlgoHex
{

template<class MeshT>
class DetachAtInvalidSingularNodeT : public virtual MeshPropertiesT<MeshT>, public ArcZippingT<MeshT>
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
  using MeshPropertiesT<MeshT>::sgl_vt_;


  DetachAtInvalidSingularNodeT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh), ArcZippingT<MeshT>(_mesh),
                                               mesh_(_mesh),
                                               sge_(mesh_),
                                               es_(mesh_),
                                               ssaf_(mesh_),
                                               opc_(mesh_)
//        ,lmcwf_(mesh_)
  {}

public:
  bool fix_interior_invalid_singular_node(const VH _vh);

  bool fix_invalid_singular_node_on_feature_surface(const VH _vh);

  bool fix_boundary_invalid_singular_node(const VH _vh, bool _start_from_bdy_sge = false);

  bool
  detach_singular_arc_from_singular_node(const VH _vh, const bool _find_featureface_sges, const bool _same_sector_sges);

  bool detach_singular_arc_from_singular_node(const VH _vh, const HEH _heh, const std::vector<HEH> &_excl_hehs,
                                              const bool _find_featureface_sges, const bool _same_sector_sges,
                                              const bool _split = true);
//        bool detach_singular_arc_from_singular_node(const VH _vh, MeshT& _mesh, std::vector<OVM::Vec4f>& _mesh_color,
//                                                    const bool _find_featureface_sge, const bool _same_sector_sges);

  bool detach_singular_arc_from_singular_node_wrt_boundary(const VH _vh, const bool _is_anti_prl,
                                                           const bool _start_from_bdy);

  bool
  detach_singular_arcs_of_zipper_node(const VH _vh, const bool _find_featureface_sge, const bool _same_sector_sges);

  bool detach_singular_arcs_of_zipper_node(const VH _vh, const HEH _heh_s, const std::vector<HEH> &_excl_hehs,
                                           const bool _find_featureface_sges, const bool _same_sector_sges,
                                           const bool _split = true);
//        bool detach_singular_arcs_of_zipper_node(const VH _vh, MeshT& _mesh, std::vector<OVM::Vec4f>& _mesh_color,
//                                                  const bool _find_featureface_sge, const bool _same_sector_sges);

  //similar as boundary, but create a zipper node on the other side
  bool detach_singualr_arc_from_singular_node_wrt_feature_face(const VH _vh, const bool _is_anti_prl);

  bool detaching_singularity(const VH _vh, const std::set<FH> &cut_fhs, const std::set<EH> &_bound_ehs);

  int n_incident_singular_edges(const VH _vh) const;

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

public:
  bool is_locally_meshable_simplified(const VH _vh, int _n_sges);

  bool is_locally_meshable(const VH _vh);

  void set_only_feature_face_singular_edges_as_targets(
          const bool _feature_face_sges) { ssaf_.find_feature_face_sges_only_ = _feature_face_sges; }

  void set_same_sector_singular_edges_as_targets(
          const bool _same_sector_sges) { ssaf_.same_sector_sges_only_ = _same_sector_sges; }

  void set_angle_threhold(const double _angle) { ssaf_.angle_threshold_ = _angle; }

  void set_zipper_node_angle_threhold(const double _angle) { ssaf_.angle_threshold_tp_ = _angle; }

  void set_enable_blocking(const bool _enable_blocking) { ssaf_.enable_blocking_ = _enable_blocking; }

  void
  set_enable_detach_kept_singular_edge(const bool _enable_detach) { ssaf_.enable_detach_kept_sge_ = _enable_detach; }

  bool is_detaching_kept_singular_edge_enabled() { return ssaf_.enable_detach_kept_sge_; }

  bool is_valid_face_sector_after_fix(const VH _vh, const HEH _ff_sgheh, const HEH _nff_sgheh)
  {
    return ssaf_.is_valid_sector_after_fix(_vh, _ff_sgheh, _nff_sgheh);
  }

private:
  bool is_singular_circle_noise(const HEH _heh) const;

  bool is_detachable_arc(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_nff_val_u, int _n_int) const;

  bool is_detachable_arc_to_ffe(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_ffe_val_1, int _n_ffe_val1,
                                int _n_ffe_val_u) const;

  bool is_detachale_arc_wrt_boundary(int _n_int_val_1, int _n_int_val1, int _n_int, int _n_bdy, int _n_ft) const;

  bool
  is_detachale_arc_wrt_boundary_from_boundary_singular_edge(int _n_int_val_1, int _n_int_val1, int _n_int, int _n_bdy,
                                                            int _n_ft) const;

  bool is_detachable_arc_wrt_feature_face(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_int) const;

  bool is_detachable_zipper_node(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_int) const;

  bool
  is_detachable_zipper_node_to_ffe(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_ffe_val_1, int _n_ffe_val1,
                                   int _n_ffe_val_u) const;

private:
  MeshT &mesh_;

  AlgoHex::SingularGraphExtractionT<MeshT> sge_;
  EdgeSplitT<MeshT> es_;

  SeparableSingularArcFinderT<MeshT> ssaf_;

  AlgoHex::OneringParameterizationChecker<MeshT> opc_;
//        AlgoHex::LocalMeshabilityCheckerWithFrames<MeshT> lmcwf_;

  AlgoHex::TransitionQuaternion tq_;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(FIXINVALIDSINGULARNODET_C)
#define FIXINVALIDSINGULARNODET_TEMPLATES

#include "DetachAtInvalidSingularNodeT_impl.hh"

#endif
//=============================================================================

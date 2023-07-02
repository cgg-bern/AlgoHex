/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "MeshProperties.hh"
#include <AlgoHex/TransitionQuaternionEigen.hh>

namespace AlgoHex
{
template<class MeshT>
class ArcZippingT : public virtual MeshPropertiesT<MeshT>
{
public:
  using Point = typename MeshT::PointT;

  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;
  using MeshPropertiesT<MeshT>::feature_edge_;


  ArcZippingT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                              mesh_(_mesh) {}

public:
  int n_cut_faces(const EH _eh, const std::set<FH> &_visited_fhs, FH &_fh) const;

  int n_cut_faces(const HEH _heh, HFH &_hfh) const;

  std::set<FH>
  generate_cut_faces(const VH _vh, const std::set<CH> _onering_chs, const std::vector<HEH> &_bdy_excluded_hehs,
                     const std::vector<HEH> &_interior_excluded_hehs, const std::set<EH> &_bound_ehs) const;

  void coordinate_transformation_in_cell(const HFH hfh);

  void regular_edge_closing(const HEH _heh, const HFH _hfh);

  //change the matchings on onering cut faces (bound edges can be regular)
  bool zip_edges_onering(const VH _vh, const EH start_eh, const std::set<EH> &_bound_ehs, const std::set<FH> &_cut_fhs);

  std::set<EH> zipping_onering(const VH _vh, const std::set<FH> &cut_fhs, const std::set<EH> &_bound_ehs);

  bool zipping_surface(const std::set<HFH> &_surface, const std::set<EH> &_edges, const HEH _sg_heh,
                       const int _new_he_trans);

  std::vector<VH> shortest_path_vertex_to_vertex_on_boundary(const VH &_vh_s, const VH &_vh_t) const;

  std::set<FH>
  generate_cut_surface_with_bounding_edges(const std::set<EH> &_bdy_ehs, const EH &_sg_eh0, const EH &_sg_eh1,
                                           const VH &_vh_c) const;

private:
  MeshT &mesh_;

  TransitionQuaternion tq_;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ARCZIPPINGT_C)
#define ARCZIPPINGT_TEMPLATES

#include "ArcZippingT_impl.hh"

#endif
//=============================================================================



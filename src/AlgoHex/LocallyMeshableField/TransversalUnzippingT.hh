/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "ArcZippingT.hh"
#include <AlgoHex/SingularGraphExtractionT.hh>
#include "TetRemeshingT.hh"

namespace AlgoHex
{

template<class MeshT>
class TransversalUnzippingT : public virtual MeshPropertiesT<MeshT>, public ArcZippingT<MeshT>
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

  TransversalUnzippingT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh), ArcZippingT<MeshT>(_mesh),
                                        mesh_(_mesh),
                                        sge_(mesh_),
                                        es_(mesh_),
                                        fs_(mesh_) {}


public:
  std::vector<HEH> repair_fully_constrained_parabolic_sector(const VH _vh, const std::vector<HEH> &_sector_hehs,
                                                             const std::set<FH> &_sec_fhs, bool _shrink_sector = false);

  HEH transversal_unzipping(const HEH _heh_s, const HFH _ft_hfh, const bool _shrink_sector);

  void get_parabolic_feature_sectors_at_vertex(const VH _vh, std::vector<std::vector<HEH>> &_zero_sectors_hehs,
                                               std::vector<std::set<FH>> &_zero_sectors_fhs);

  bool is_parabolic_sector(const std::vector<HEH> &_zero_sector_hehs, const std::set<FH> &_zero_sectors_fhs) const;

private:
  MeshT &mesh_;
  TransitionQuaternion tq_;

  AlgoHex::SingularGraphExtractionT<MeshT> sge_;
  EdgeSplitT<MeshT> es_;
  FaceSplitT<MeshT> fs_;

};
}



//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(TRANSVERSALUNZIPPINGT_C)
#define TRANSVERSALUNZIPPINGT_TEMPLATES

#include "TransversalUnzippingT_impl.hh"

#endif
//=============================================================================


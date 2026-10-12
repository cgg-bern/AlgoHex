/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once
#include <AlgoHex/TypeDef.hh>


#include "ArcZipping.hh"
#include <AlgoHex/SingularGraphExtraction.hh>
#include "TetRemeshing.hh"

namespace AlgoHex
{

class ALGOHEX_EXPORT TransversalUnzipping : public virtual MeshProperties, public ArcZipping
{
public:
  using ArcZipping::valence_;
  using ArcZipping::trans_prop_;
  using ArcZipping::cell_quaternions_;

  using MeshProperties::feature_fprop_;
  using MeshProperties::feature_face_edge_;
  using MeshProperties::feature_face_vertex_;

  using ArcZipping::feature_edge_;
  using MeshProperties::feature_edge_vertex_;

  using MeshProperties::feature_node_;

  using MeshProperties::sgl_vt_;

  TransversalUnzipping(TetMesh &_mesh);

  std::vector<HEH> repair_fully_constrained_parabolic_sector(const VH _vh, const std::vector<HEH> &_sector_hehs,
                                                             const std::set<FH> &_sec_fhs, bool _shrink_sector = false);

  HEH transversal_unzipping(const HEH _heh_s, const HFH _ft_hfh, const bool _shrink_sector);

  void get_parabolic_feature_sectors_at_vertex(const VH _vh, std::vector<std::vector<HEH>> &_zero_sectors_hehs,
                                               std::vector<std::set<FH>> &_zero_sectors_fhs);

  bool is_parabolic_sector(const std::vector<HEH> &_zero_sector_hehs, const std::set<FH> &_zero_sectors_fhs) const;

private:
  TetMesh &mesh_;
  TransitionQuaternion tq_;

  AlgoHex::SingularGraphExtraction sge_;
  EdgeSplit es_;
  FaceSplit fs_;

};
}



//=============================================================================
//=============================================================================

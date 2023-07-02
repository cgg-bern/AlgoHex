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
#include <queue>
#include <AlgoHex/SingularGraphExtractionT.hh>
#include "../QuaternionsSmoothing.hh"
#include "../SplitHelperT.hh"


namespace AlgoHex
{

template<class MeshT>
class CollapseApproachT : public virtual MeshPropertiesT<MeshT>
{
public:
  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;

  CollapseApproachT(MeshT &_mesh) :
          MeshPropertiesT<MeshT>(_mesh),
          mesh_(_mesh),
          sge_(mesh_),
          ec_(mesh_),
          es_(mesh_),
          vo_(mesh_)
  {
    this->initialize_singular_vertex_property();
  }

  void collapse_approach();

  void matching_adjustment();

  void remove_zigzag();

  void collapse_complex_singular_edges();

  //debug
  void collapse_complex_singular_edge(const HEH &he_cur);

  int n_legal_singular_edges_in_face(const FH &_fh) const;

  int n_singular_edges_in_face(const FH &_fh) const;

  std::vector<std::pair<double, int>> sorted_matchings(const HFH &_hfh) const;

  int n_complex_singular_edges() const;

private:
  MeshT &mesh_;
  AlgoHex::SingularGraphExtractionT<MeshT> sge_;
  EdgeCollapseT<MeshT> ec_;
  EdgeSplitT<MeshT> es_;
  VertexOptimizeT<MeshT> vo_;

  AlgoHex::TransitionQuaternion tq_;

};
}


//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(COLLAPSEPAPERT_C)
#define COLLAPSEPAPERT_TEMPLATES

#include "CollapsePaperT_impl.hh"

#endif
//=============================================================================



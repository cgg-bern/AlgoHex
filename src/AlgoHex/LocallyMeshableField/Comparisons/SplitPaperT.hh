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
#include "../LocalMeshabilityRepairT.hh"


namespace AlgoHex
{

template<class MeshT>
class SplitApproachT : public virtual MeshPropertiesT<MeshT>
{
public:
  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;
  using MeshPropertiesT<MeshT>::feature_face_vertex_;
  using MeshPropertiesT<MeshT>::feature_edge_vertex_;

  SplitApproachT(MeshT &_mesh) :
          MeshPropertiesT<MeshT>(_mesh),
          mesh_(_mesh),
          sge_(mesh_),
          pp_(mesh_),
          flnmv_(mesh_),
          eem_(mesh_),
          fac_(mesh_),
          es_(mesh_)
  {
    this->initialize_singular_vertex_property();
  }

  void split_approach();

  void remove_zigzag();

  void remove_zigzag2();

  void split_complex_singular_edges();

  int n_incident_complex_singular_edges(const VH &_vh) const;

  bool is_removable_zigzag(const FH &_fh) const;

  bool is_removable_zigzag2(const FH &_fh) const;

  int n_legal_singular_edges_in_face(const FH &_fh) const;

  int n_singular_edges_in_face(const FH &_fh) const;

  bool check_zigzag();

  int n_complex_singular_edges();


private:
  MeshT &mesh_;

  AlgoHex::SingularGraphExtractionT<MeshT> sge_;
  AlgoHex::TransitionQuaternion tq_;
  SplitHelperT<MeshT> pp_;
  LocalMeshabilityRepairT<MeshT> flnmv_;
  EnsureEdgeMeshabilityT<MeshT> eem_;
  FieldAngleCalculatorT<MeshT> fac_;
  EdgeSplitT<MeshT> es_;
};

}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(SPLITPAPERT_C)
#define SPLITPAPERT_TEMPLATES

#include "SplitPaperT_impl.hh"

#endif
//=============================================================================



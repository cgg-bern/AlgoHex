/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once
#include <AlgoHex/TypeDef.hh>


#include <set>
#include <vector>
#include <queue>
#include <AlgoHex/SingularGraphExtraction.hh>
#include "../QuaternionsSmoothing.hh"
#include "../SplitHelper.hh"
#include "../LocalMeshabilityRepair.hh"


namespace AlgoHex
{

class ALGOHEX_EXPORT SplitApproach : public virtual MeshProperties
{
public:
  using MeshProperties::valence_;
  using MeshProperties::trans_prop_;
  using MeshProperties::sgl_vt_;
  using MeshProperties::feature_edge_;
  using MeshProperties::feature_fprop_;
  using MeshProperties::feature_face_edge_;
  using MeshProperties::cell_quaternions_;
  using MeshProperties::feature_face_vertex_;
  using MeshProperties::feature_edge_vertex_;

  SplitApproach(TetMesh &_mesh);

  void split_approach();

  void remove_zigzag();

  void remove_zigzag2();

  void split_complex_singular_edges();

  int n_incident_complex_singular_edges(const VH &_vh) const;

  bool is_removable_zigzag(const FH &_fh) const;

  int n_legal_singular_edges_in_face(const FH &_fh) const;

  int n_singular_edges_in_face(const FH &_fh) const;

  bool check_zigzag();

  int n_complex_singular_edges();

private:
  TetMesh &mesh_;

  AlgoHex::SingularGraphExtraction sge_;
  AlgoHex::TransitionQuaternion tq_;
  SplitHelper pp_;
  LocalMeshabilityRepair flnmv_;
  EnsureEdgeMeshability eem_;
  FieldAngleCalculator fac_;
  EdgeSplit es_;
};

}

//=============================================================================
//=============================================================================

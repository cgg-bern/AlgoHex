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


namespace AlgoHex
{

class ALGOHEX_EXPORT CollapseApproach : public virtual MeshProperties
{
public:
  using MeshProperties::valence_;
  using MeshProperties::trans_prop_;
  using MeshProperties::sgl_vt_;
  using MeshProperties::feature_edge_;
  using MeshProperties::feature_fprop_;
  using MeshProperties::cell_quaternions_;
  using MeshProperties::feature_face_edge_;

  CollapseApproach(TetMesh &_mesh) :
          MeshProperties(_mesh),
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
  TetMesh &mesh_;
  AlgoHex::SingularGraphExtraction sge_;
  EdgeCollapse ec_;
  EdgeSplit es_;
  VertexOptimize vo_;

  AlgoHex::TransitionQuaternion tq_;

};
}


//=============================================================================
//=============================================================================

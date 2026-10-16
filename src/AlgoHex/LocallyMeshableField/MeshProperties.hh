/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include <AlgoHex/TypeDef.hh>

namespace AlgoHex
{
class ALGOHEX_EXPORT MeshProperties
{
public:
  MeshProperties(TetMesh &_mesh);



  void initialize_singular_vertex_property();


  void import_feature_properties();


  void initialize_feature_vertex_property();


  void initialize_target_edge_length();


  void smooth_target_edge_length(const int _iter = 5, const double _damping = 0.8);


  void localize_properties();


  int n_incident_feature_edges(const VH _vh) const;


  void update_singular_vertex_property(const VH _vh);



  //heuristically get the node index by counting edge valence
  int node_index(const VH _vh) const;


  bool is_singular_node(const VH _vh) const;



protected:
  TetMesh &mesh_;

  //input
  EP<int> valence_;
  HFP<int> trans_prop_;
  CP <Quaternion> cell_quaternions_;

  VP<int> sgl_vt_;
  VP<double> target_length_;

  VP<int> feature_node_;

  VP<bool> feature_edge_vertex_;
  EP<int> feature_edge_;

  VP<bool> feature_face_vertex_;
  EP<bool> feature_face_edge_;
  FP<int> feature_fprop_;
};
}

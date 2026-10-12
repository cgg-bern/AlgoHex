/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once
#include <AlgoHex/TypeDef.hh>


#include "MeshProperties.hh"
#include "TetRemeshing.hh"
#include "CommonFuncs.hh"
#include <AlgoHex/MeshGeometry.hh>


namespace AlgoHex
{
class ALGOHEX_EXPORT SplitHelper : public MeshProperties
{
public:
  using MeshProperties::valence_;
  using MeshProperties::sgl_vt_;
  using MeshProperties::feature_edge_vertex_;
  using MeshProperties::feature_edge_;
  using MeshProperties::feature_fprop_;
  using MeshProperties::feature_face_edge_;
  using MeshProperties::feature_face_vertex_;

  SplitHelper(TetMesh &_mesh);

  void split_for_single_alignment_feature_face();

  bool is_edge_to_split_fe(const EH _eh, const bool _skip_two_sgv = false);

  bool is_edge_to_split_ffe(const EH _eh, const bool _skip_two_sgv = false);

  //it doesn't split regular edge with two singular vertices if true
  void split_for_single_alignment(const bool _skip_two_sgv = false);

  void split_for_good_tet_quality_at_singular_edges();

  void split_special_edges_with_two_special_nodes();

  //don't split singular edges
  static void
  split_one_ring(TetMesh &mesh, EdgeSplit &es, const VH _vh, EP<int> &valence, EP<bool> &feature_face_edge,
                 VP<bool> &feature_face_vertex, EP<int> &feature_edge, VP<int> &sgl_vt, VP<bool> &feature_edge_vertex,
                 const bool collect_garbage = false);

  //singular edges can be split
  static void split_one_ring_all(TetMesh &mesh, EdgeSplit &es, const VH _vh,
                                 EP<int> &valence, EP<bool> &feature_face_edge, VP<bool> &feature_face_vertex,
                                 EP<int> &feature_edge, VP<int> &sgl_vt, VP<bool> &feature_edge_vertex,
                                 const bool collect_garbage = false);

  //singular edges can be split
  static void split_one_ring_all(TetMesh &mesh, EdgeSplit &es, const VH _vh,
                                 EP<int> &valence, EP<bool> &feature_face_edge, VP<bool> &feature_face_vertex,
                                 EP<int> &feature_edge, VP<int> &sgl_vt, VP<bool> &feature_edge_vertex,
                                 std::set<EH> &_excl_ehs, const bool collect_garbage = false);

  static void
  split_singular_edges_with_singular_node(TetMesh &mesh, EdgeSplit &es, EP<int> &valence, const VH _vh,
                                          const bool collect_garbage = false);

  static void
  split_feature_face_singular_edges_with_singular_nodes(TetMesh &mesh, EdgeSplit &es, EP<bool> &feature_face_edge,
                                                        EP<int> &valence, const bool collect_garbage = false);

  static void split_onering_of_zipper_nodes_arc(TetMesh &mesh, EdgeSplit &es, EP<int> &valence, VP<int> &sgl_vt,
                                                std::set<VH> &_tps, const bool collect_garbage = false);


  static int split_faces_with_cell_property_update(TetMesh &mesh, FaceSplit &_fs, std::set<CH> &_dp_cells,
                                                   const std::set<FH> &_fhs, std::set<CH> &_start_cells,
                                                   std::set<CH> &_target_cells, std::set<FH> &_excl_fhs);

  static int split_edges_with_cell_property_update(TetMesh &mesh, EdgeSplit &es, std::set<CH> &_dp_cells,
                                                   const std::set<EH> &_ehs, std::set<CH> &_start_cells,
                                                   std::set<CH> &_target_cells, std::set<EH> &_excl_ehs);

  static int split_edge_cells_at_high_valence(TetMesh &mesh, CellSplit &cs, const HEH heh, FP<int> &feature_fprop,
                                              EP<int> &valence, EP<bool> &feature_face_edge);

  static void split_for_dof(TetMesh &mesh);

  static void sort_singular_arc_label(const TetMesh &mesh, const EP<int> &valence, EP<int> &label);

  //split if the surface mesh of the onering has >= two singular edges
  static void
  split_edges_for_local_meshability_test(TetMesh &mesh, EdgeSplit &es, EP<int> &valence, EP<int> &feature_edge);

  int n_incident_special_edges(const VH &_vh) const;

private:
  TetMesh &mesh_;
  EdgeSplit es_;
  FaceSplit fs_;
  VertexOptimize vo_;
};
}

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once
#include <AlgoHex/TypeDef.hh>


#include <AlgoHex/Util/json.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include "DuplicateOneRingMesh.hh"
#include "../FrameFieldOptimizer3D.hh"
#include <AlgoHex/MeshGeometry.hh>
#include "TetRemeshing.hh"
#include <AlgoHex/SingularGraphExtraction.hh>
#include "QuaternionsSmoothing.hh"
#include "FieldAngleCalculator.hh"

#ifdef ALGOHEX_VERBOSE
#define ALGOHEX_DEBUG_ONLY(x) x
#else
#define ALGOHEX_DEBUG_ONLY(x)
#endif

namespace AlgoHex
{
class ALGOHEX_EXPORT LocalMeshabilityCheckerWithFrames : public MeshProperties
{
public:
  using MeshProperties::valence_;
  using MeshProperties::trans_prop_;
  using MeshProperties::cell_quaternions_;
  using MeshProperties::sgl_vt_;
  using MeshProperties::feature_node_;
  using MeshProperties::feature_edge_;
  using MeshProperties::feature_edge_vertex_;
  using MeshProperties::feature_fprop_;
  using MeshProperties::feature_face_edge_;

  using PairDE = std::pair<double, EH>;

  LocalMeshabilityCheckerWithFrames(TetMesh &_mesh);

  ~LocalMeshabilityCheckerWithFrames() = default;

public:
  nlohmann::json &json_data();

  void optimize_onering_mesh(TetMesh &_or_mesh, const VH _vh, int _opt_iter = 10);

  bool is_locally_meshable_edge(const EH _eh, const bool _align_sg = false);

  bool is_locally_meshable(const VH _vh, const bool _align_sg = false, const int _n_split = 1);

  bool check_parameterization(TetMesh &_mesh, const VH _vh, const int _orig_vh_id, bool _save_config);

  double check_local_meshability(const bool _align_to_sge = false);

  bool is_on_feature_vertex(const VH _vh) const;

  bool is_on_feature_edge(const VH _vh) const;

  bool is_on_feature_face(const VH _vh) const;

  bool is_on_singular_node(const VH _vh) const;

  bool is_zipper_node(const VH _vh) const;

  bool is_on_singular_arc(const VH _vh) const;

  bool is_singular_vertex(const VH _vh) const;

  void enable_save_locally_non_meshable(const std::string _filename_base = "locally_non_meshable_");

  double len_complex_singular_edges() const;

  double len_singular_edges() const;

  int split_edges(TetMesh &_or_mesh, const VH _vh);

  void normalize_one_ring_mesh(TetMesh &_mesh, const VH _vh_center) const;

  void remove_locally_irrelevant_constraints(TetMesh &_mesh, const VH _vh_center) const;

  bool &verbose();

private:
  TetMesh &mesh_;
  AlgoHex::TransitionQuaternion tq_;

  bool save_locally_non_meshable_ = false;
  std::string save_locally_non_meshable_filename_base_;

  int n_opt_ = 0;
  int n_fixed_ = 0;

  bool verbose_ = false;

  // json storage
  nlohmann::json json_data_;
};
}
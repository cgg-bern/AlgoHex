/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include <AlgoHex/Util/json.hh>
#include <AlgoHex/TypeDef.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include "DuplicateOneRingMesh.hh"
#include "../Parametrization3D.hh"
#include <AlgoHex/MeshGeometry.hh>
#include <Base/Debug/DebConfig.hh>
#include "TetRemeshing.hh"
#include "../SingularGraphExtraction.hh"
#include "QuaternionsSmoothing.hh"


namespace AlgoHex
{
class ALGOHEX_EXPORT OneringParameterizationChecker : public MeshProperties
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

  using PairDE = std::pair<double, EH>;

  OneringParameterizationChecker(TetMesh &_mesh);

  ~OneringParameterizationChecker() = default;

public:
  nlohmann::json &json_data();

  void optimize_onering_mesh(TetMesh &_or_mesh, const VH _vh, int _opt_iter = 10);

  int split_edges(TetMesh &_or_mesh, const VH _vh);

  bool is_locally_meshable(const VH _vh, const bool _align_sg = false);

  bool is_locally_meshable2(TetMesh &_onering, const VH _vh, const int _opt_iters, const bool _align_sg = false);

  bool check_parameterization(TetMesh &_mesh, const VH _vh);

  double check_local_meshability(const bool _align_to_sge = false);

  bool is_on_feature_vertex(const VH _vh) const;

  bool is_on_feature_edge(const VH _vh) const;

  bool is_on_feature_face(const VH _vh) const;

  bool is_on_singular_node(const VH _vh) const;

  bool is_zipper_node(const VH _vh) const;

  bool is_on_singular_arc(const VH _vh) const;

  bool is_singular_vertex(const VH _vh) const;

  double len_complex_singular_edges() const;

  double len_singular_edges() const;

  void enable_save_locally_non_meshable(const std::string _filename_base = "locally_non_meshable_");

private:
  TetMesh &mesh_;
  AlgoHex::TransitionQuaternion tq_;

  bool save_locally_non_meshable_ = false;
  std::string save_locally_non_meshable_filename_base_;

  // json storage
  nlohmann::json json_data_;
};
}

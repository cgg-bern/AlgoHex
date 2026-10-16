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
#include "LocalMeshabilityCheckerWithFrames.hh"
#include "OneringParameterizationChecker.hh"


namespace AlgoHex
{
class ALGOHEX_EXPORT LocalMeshabilityChecker : public MeshProperties
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

  LocalMeshabilityChecker(TetMesh &_mesh);


  ~LocalMeshabilityChecker() = default;


public:
  nlohmann::json &json_data();


  double check_edge_local_meshability(const bool _align_to_sge = false);


  double check_local_meshability(const bool _align_to_sge = false);


  bool check_special_vertices_local_meshability(const bool _align_to_sge = false);



  bool is_on_feature_vertex(const VH _vh) const;



  bool is_on_feature_edge(const VH _vh) const;



  bool is_on_feature_face(const VH _vh) const;


  bool is_on_singular_node(const VH _vh) const;


  bool is_zipper_node(const VH _vh) const;



  bool is_on_singular_arc(const VH _vh) const;


  bool is_singular_vertex(const VH _vh) const;


  void enable_save_locally_non_meshable(const std::string _filename_base = "locally_non_meshable_");


  int n_complex_singular_edges() const;


  int n_singular_edges() const;


  double len_complex_singular_edges() const;


  double len_singular_edges() const;


  bool &verbose();


private:
  TetMesh &mesh_;
  OneringParameterizationChecker opc_;
  LocalMeshabilityCheckerWithFrames lmcwf_;

  bool save_locally_non_meshable_ = false;
  std::string save_locally_non_meshable_filename_base_;

  // json storage
  nlohmann::json json_data_;
};
}

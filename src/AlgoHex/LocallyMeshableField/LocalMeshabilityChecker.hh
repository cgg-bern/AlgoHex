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
template<class MeshT>
class LocalMeshabilityChecker : public MeshPropertiesT<MeshT>
{
public:
  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;
  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::feature_node_;
  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_edge_vertex_;
  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;

  LocalMeshabilityChecker(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                                          mesh_(_mesh),
                                          opc_(mesh_),
                                          lmcwf_(mesh_) {}

  ~LocalMeshabilityChecker() = default;


public:
  nlohmann::json &json_data() { return json_data_; }

  double check_edge_local_meshability(const bool _align_to_sge = false)
  {
    Debug::ScopedOutputLevel sol(0);

    std::cout << "##### Check local meshability of edges ..." << std::endl;
    int n(0.0);
    int n_meshable(0.0);

    int n_on_feature_edge(0);
    int n_on_feature_edge_meshable(0);
    int n_on_feature_face(0);
    int n_on_feature_face_meshable(0);

    int n_on_singular_arc(0);
    int n_on_singular_arc_meshable(0);

    // check all vertices
    for (auto ehi: mesh_.edges())
    {
      // global counter
      ++n;

      // only use frame-field based checker
      bool ilm = lmcwf_.is_locally_meshable_edge(ehi, _align_to_sge);

      if (ilm)
      {
        ++n_meshable;
      }
      else
        std::cerr << "edge " << ehi << " is not locally meshable." << std::endl;


      // count on feature edge
      if (feature_edge_[ehi] > 0)
      {
        ++n_on_feature_edge;
        if (ilm)
          ++n_on_feature_edge_meshable;
      }

      // count on feature face
      if (feature_face_edge_[ehi])
      {
        ++n_on_feature_face;
        if (ilm)
          ++n_on_feature_face_meshable;
      }

      // count on singular arc
      if (valence_[ehi] != 0)
      {
        ++n_on_singular_arc;
        if (ilm)
          ++n_on_singular_arc_meshable;
      }
    }

    // avoid NAN
    std::cerr << "#Edges          = " << std::setw(5) << n << ", #meshable = " << std::setw(5) << n_meshable << " ("
              << double(n_meshable) / double(n) * 100.0 << "%)\n";
    std::cerr << "#on feature edge   = " << std::setw(5) << n_on_feature_edge << ", #meshable = " << std::setw(5)
              << n_on_feature_edge_meshable << " ("
              << double(n_on_feature_edge_meshable) / double(n_on_feature_edge) * 100.0 << "%)\n";
    std::cerr << "#on feature face   = " << std::setw(5) << n_on_feature_face << ", #meshable = " << std::setw(5)
              << n_on_feature_face_meshable << " ("
              << double(n_on_feature_face_meshable) / double(n_on_feature_face) * 100.0 << "%)\n";
    std::cerr << "#on singular arc   = " << std::setw(5) << n_on_singular_arc << ", #meshable = " << std::setw(5)
              << n_on_singular_arc_meshable << " ("
              << double(n_on_singular_arc_meshable) / double(n_on_singular_arc) * 100.0 << "%)\n";

    json_data_["percentage_meshable_edges"] = double(n_meshable) / double(n);
    json_data_["non-meshable edges"] = n - n_meshable;
    json_data_["all edges"] = n;
    json_data_["feature edges"] = n_on_feature_edge - n_on_feature_edge_meshable;
    json_data_["complex singular edges"] = n_complex_singular_edges();
    json_data_["singular edges"] = n_singular_edges();
    // return percentage of meshable vertices
    return double(n_meshable) / double(n);
  }

  double check_local_meshability(const bool _align_to_sge = false)
  {
    Debug::ScopedOutputLevel sol(0);

    std::cout << "##### Check local meshability ..." << std::endl;
    int n(0.0);
    int n_meshable(0.0);

    int n_on_feature_vertex(0);
    int n_on_feature_vertex_meshable(0);
    int n_on_feature_edge(0);
    int n_on_feature_edge_meshable(0);
    int n_on_feature_face(0);
    int n_on_feature_face_meshable(0);

    int n_on_singular_vertex(0);
    int n_on_singular_vertex_meshable(0);
    int n_on_singular_node(0);
    int n_on_singular_node_meshable(0);
    int n_zipper_node(0);
    int n_zipper_node_meshable(0);
    int n_on_singular_arc(0);
    int n_on_singular_arc_meshable(0);

    //save path
    if (save_locally_non_meshable_)
      lmcwf_.enable_save_locally_non_meshable(save_locally_non_meshable_filename_base_);
//            if(save_locally_non_meshable_)
//                opc_.enable_save_locally_non_meshable(save_locally_non_meshable_filename_base_);

    // check all vertices
    for (VIt v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
    {
      // global counter
      ++n;
//                bool ilm = opc_.is_locally_meshable(*v_it, _align_to_sge);
//                if(!ilm)
//                    ilm = lmcwf_.is_locally_meshable(*v_it, _align_to_sge);

      // only use frame-field based checker
      std::cerr << "check vh " << *v_it << " ";
      bool ilm = lmcwf_.is_locally_meshable(*v_it, _align_to_sge);

      if (ilm)
      {
        ++n_meshable;
      }
      else
        std::cerr << "vertex " << *v_it << " is not locally meshable." << std::endl;

      // count on feature vertex
      if (is_on_feature_vertex(*v_it))
      {
        ++n_on_feature_vertex;
        if (ilm)
          ++n_on_feature_vertex_meshable;
      }

      // count on feature edge
      if (is_on_feature_edge(*v_it))
      {
        ++n_on_feature_edge;
        if (ilm)
          ++n_on_feature_edge_meshable;
      }

      // count on feature face
      if (is_on_feature_face(*v_it))
      {
        ++n_on_feature_face;
        if (ilm)
          ++n_on_feature_face_meshable;
      }

      // count on singular node
      if (is_singular_vertex(*v_it))
      {
        ++n_on_singular_vertex;
        if (ilm)
          ++n_on_singular_vertex_meshable;
      }

      // count on singular node
      if (is_on_singular_node(*v_it))
      {
        ++n_on_singular_node;
        if (ilm)
          ++n_on_singular_node_meshable;
      }

      // count turning points
      if (is_zipper_node(*v_it))
      {
        ++n_zipper_node;
        if (ilm)
          ++n_zipper_node_meshable;
      }

      // count on singular arc
      if (is_on_singular_arc(*v_it))
      {
        ++n_on_singular_arc;
        if (ilm)
          ++n_on_singular_arc_meshable;
      }
    }
    std::cerr << std::endl;
    // avoid NAN
    double feature_vertex_ratio = double(n_on_feature_vertex_meshable) / double(n_on_feature_vertex);
    if (!std::isfinite(feature_vertex_ratio))
      feature_vertex_ratio = 1.0;

    std::cerr << "#vertices          = " << std::setw(5) << n << ", #meshable = " << std::setw(5) << n_meshable << " ("
              << double(n_meshable) / double(n) * 100.0 << "%)\n";
    std::cerr << "#on feature node = " << std::setw(5) << n_on_feature_vertex << ", #meshable = " << std::setw(5)
              << n_on_feature_vertex_meshable << " (" << feature_vertex_ratio * 100.0 << "%)\n";
    std::cerr << "#on feature edge   = " << std::setw(5) << n_on_feature_edge << ", #meshable = " << std::setw(5)
              << n_on_feature_edge_meshable << " ("
              << double(n_on_feature_edge_meshable) / double(n_on_feature_edge) * 100.0 << "%)\n";
    std::cerr << "#on feature face   = " << std::setw(5) << n_on_feature_face << ", #meshable = " << std::setw(5)
              << n_on_feature_face_meshable << " ("
              << double(n_on_feature_face_meshable) / double(n_on_feature_face) * 100.0 << "%)\n";
    std::cerr << "#on singular node  = " << std::setw(5) << n_on_singular_node << ", #meshable = " << std::setw(5)
              << n_on_singular_node_meshable << " ("
              << double(n_on_singular_node_meshable) / double(n_on_singular_node) * 100.0 << "%)\n";
    std::cerr << "#turning point  = " << std::setw(5) << n_zipper_node << ", #meshable = " << std::setw(5)
              << n_zipper_node_meshable << " (" << double(n_zipper_node_meshable) / double(n_zipper_node) * 100.0
              << "%)\n";
    std::cerr << "#on singular arc   = " << std::setw(5) << n_on_singular_arc << ", #meshable = " << std::setw(5)
              << n_on_singular_arc_meshable << " ("
              << double(n_on_singular_arc_meshable) / double(n_on_singular_arc) * 100.0 << "%)\n";
    std::cerr << "#singular vertices   = " << std::setw(5) << n_on_singular_vertex << ", #meshable = " << std::setw(5)
              << n_on_singular_vertex_meshable << " ("
              << double(n_on_singular_vertex_meshable) / double(n_on_singular_vertex) * 100.0 << "%)\n";
    std::cerr << "#feature vertices   = " << std::setw(5) << n_on_feature_edge << ", #meshable = " << std::setw(5)
              << n_on_feature_edge_meshable << " ("
              << double(n_on_feature_edge_meshable) / double(n_on_feature_edge) * 100.0 << "%)\n";

    json_data_["percentage_meshable_vertices"] = double(n_meshable) / double(n);
    json_data_["non-meshable"] = n - n_meshable;
    json_data_["all vertices"] = n;
    json_data_["singular nodes"] = n_on_singular_node - n_on_singular_node_meshable;
    json_data_["turning points"] = n_zipper_node;
    json_data_["singular vertices"] = n_on_singular_vertex - n_on_singular_vertex_meshable;
    json_data_["feature vertices"] = n_on_feature_edge - n_on_feature_edge_meshable;
    json_data_["complex singular edges"] = len_complex_singular_edges();
    json_data_["singular edges"] = len_singular_edges();
    // return percentage of meshable vertices
    return double(n_meshable) / double(n);
  }

  bool check_special_vertices_local_meshability(const bool _align_to_sge = false)
  {
    Debug::ScopedOutputLevel sol(0);

    std::cout << "##### Check local meshability of special vertice..." << std::endl;

    // check all special vertices
    for (VIt v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
    {
      if (sgl_vt_[*v_it] || feature_edge_vertex_[*v_it])
      {
        std::cerr << "check vertex " << *v_it;
//                    bool ilm = opc_.is_locally_meshable(*v_it, _align_to_sge);
//                    if (!ilm)
        bool ilm = lmcwf_.is_locally_meshable(*v_it, _align_to_sge);

        if (!ilm)
        {
          std::cerr << "vertex " << *v_it << " is not locally meshable." << std::endl;
          return false;
        }
      }
    }

    return true;
  }


  bool is_on_feature_vertex(const VH _vh) const
  {
    // is feature vertex?
    return feature_node_[_vh];
  }


  bool is_on_feature_edge(const VH _vh) const
  {
    // is on feature edge?
    VOHEIt vhe_it(_vh, &mesh_);
    for (; vhe_it.valid(); ++vhe_it)
      if (feature_edge_[mesh_.edge_handle(*vhe_it)])
        return true;

    return false;
  }


  bool is_on_feature_face(const VH _vh) const
  {
    // is on feature triangle?
    VFIt vf_it(_vh, &mesh_);
    for (; vf_it.valid(); ++vf_it)
      if (feature_fprop_[*vf_it])
        return true;

    return false;
  }

  bool is_on_singular_node(const VH _vh) const
  {
    std::vector<int> sv;

    VOHEIt vhe_it(_vh, &mesh_);
    for (; vhe_it.valid(); ++vhe_it)
    {
      int val = valence_[mesh_.edge_handle(*vhe_it)];
      if (val != 0)
        sv.push_back(val);
    }

    if (sv.empty()) return false;
    else if (sv.size() != 2) return true;
    else if (sv.size() == 2)
      return (sv[0] != sv[1]);

    return false;
  }

  bool is_zipper_node(const VH _vh) const
  {
    std::vector<int> sv;

    VOHEIt vhe_it(_vh, &mesh_);
    for (; vhe_it.valid(); ++vhe_it)
    {
      int val = valence_[mesh_.edge_handle(*vhe_it)];
      if (val != 0)
        sv.push_back(val);
    }

    if (sv.size() == 2)
    {
      if ((sv[0] == 1 && sv[1] == -1) || (sv[0] == -1 && sv[1] == 1))
        return true;
    }

    return false;
  }


  bool is_on_singular_arc(const VH _vh) const
  {
    std::vector<int> sv;

    VOHEIt vhe_it(_vh, &mesh_);
    for (; vhe_it.valid(); ++vhe_it)
    {
      int val = valence_[mesh_.edge_handle(*vhe_it)];
      if (val != 0)
        sv.push_back(val);
    }

    if (sv.size() != 2) return false;
    else return (sv[0] == sv[1]);
  }

  bool is_singular_vertex(const VH _vh) const
  {
    return sgl_vt_[_vh] != 0;
  }

  void enable_save_locally_non_meshable(const std::string _filename_base = "locally_non_meshable_")
  {
    save_locally_non_meshable_ = true;
    save_locally_non_meshable_filename_base_ = _filename_base;
  }

  int n_complex_singular_edges() const
  {
    int n = 0;
    for (const auto ehi: mesh_.edges())
    {
      if (valence_[ehi] == std::numeric_limits<int>::max())
        n++;
    }

    return n;
  }

  int n_singular_edges() const
  {
    int n = 0;
    for (const auto ehi: mesh_.edges())
    {
      if (valence_[ehi] != 0)
        n++;
    }

    return n;
  }

  double len_complex_singular_edges() const
  {
    double len_cse = 0.;
    for (const auto ehi: mesh_.edges())
    {
      if (valence_[ehi] == std::numeric_limits<int>::max())
        len_cse += mesh_.length(ehi);
    }

    return len_cse;
  }

  double len_singular_edges() const
  {
    int len_se = 0.;
    for (const auto ehi: mesh_.edges())
    {
      if (valence_[ehi] != 0)
        len_se++;
    }

    return len_se;
  }

  bool &verbose() { return lmcwf_.verbose(); }

private:
  MeshT &mesh_;
  OneringParameterizationChecker <MeshT> opc_;
  LocalMeshabilityCheckerWithFrames <MeshT> lmcwf_;

  bool save_locally_non_meshable_ = false;
  std::string save_locally_non_meshable_filename_base_;

  // json storage
  nlohmann::json json_data_;
};
}

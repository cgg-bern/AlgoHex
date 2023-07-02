/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include <AlgoHex/Util/json.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include "DuplicateOneRingMeshT.hh"
#include "../FrameFieldOptimizer3DT.hh"
#include <AlgoHex/MeshGeometry.hh>
#include "TetRemeshingT.hh"
#include <AlgoHex/SingularGraphExtractionT.hh>
#include "QuaternionsSmoothing.hh"
#include "FieldAngleCalculatorT.hh"

#ifdef ALGOHEX_VERBOSE
#define ALGOHEX_DEBUG_ONLY(x) x
#else
#define ALGOHEX_DEBUG_ONLY(x)
#endif

namespace AlgoHex
{
template<class MeshT>
class LocalMeshabilityCheckerWithFrames : public MeshPropertiesT<MeshT>
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

  using PairDE = std::pair<double, EH>;

  LocalMeshabilityCheckerWithFrames(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                                                    mesh_(_mesh) {}

  ~LocalMeshabilityCheckerWithFrames() = default;


public:
  nlohmann::json &json_data() { return json_data_; }

  void optimize_onering_mesh(MeshT &_or_mesh, const VH _vh, int _opt_iter = 10)
  {
    auto chax = _or_mesh.template request_cell_property<std::pair<HEH, int> >("EdgeAxis");

    ALGOHEX_DEBUG_ONLY(check_mesh_quality(_or_mesh);)

    //fix feature edge vertices or singular vertices
    auto fixed_vt_prop = _or_mesh.template request_vertex_property<int>("fixed vertex", 0);
    auto ic_sge_prop = _or_mesh.template request_vertex_property<EH>("incident singular edge");
//            auto orig_bdy_e = _or_mesh.template request_edge_property<bool>("original boundary edge", false);

    _or_mesh.set_persistent(fixed_vt_prop, true);
    auto or_valence = _or_mesh.template request_edge_property<int>("edge_valance");
    auto or_trans_prop = _or_mesh.template request_halfface_property<int>("HalffaceTransiton");
    auto or_cell_quaternions = _or_mesh.template request_cell_property<Eigen::Quaterniond>(
            "FrameFieldQuaternions");
    auto fe_prop = _or_mesh.template request_edge_property<int>("AlgoHex::FeatureEdges");
    auto ff_prop = _or_mesh.template request_face_property<int>("AlgoHex::FeatureFaces");

    auto be_prop = _or_mesh.template request_edge_property<bool>("feature face edges", false);
    _or_mesh.set_persistent(be_prop, true);
    auto bv_prop = _or_mesh.template request_vertex_property<bool>("feature face vertices", false);
    _or_mesh.set_persistent(bv_prop, true);

    //
    for (const auto &fhi: _or_mesh.faces())
    {
      if (ff_prop[fhi])
      {
        for (auto fe_it = _or_mesh.fe_iter(fhi); fe_it.valid(); ++fe_it)
          be_prop[*fe_it] = true;
      }
    }
    for (const auto &fhi: _or_mesh.faces())
    {
      if (ff_prop[fhi])
      {
        for (auto fv_it = _or_mesh.fv_iter(fhi); fv_it.valid(); ++fv_it)
          bv_prop[*fv_it] = true;
      }
    }



//            for(int j=0; j<3; ++j) {
//                split_edges(_or_mesh, _vh);

    //vertex relocation
    //store incident special edge number
    for (const auto vhi: _or_mesh.vertices())
      fixed_vt_prop[vhi] = 0;

    for (const auto &ehi: _or_mesh.edges())
    {
      if ((or_valence[ehi] >= -2 && or_valence[ehi] <= 4 && or_valence[ehi] != 0) || fe_prop[ehi])
      {
        auto vhf = _or_mesh.edge(ehi).from_vertex();
        auto vht = _or_mesh.edge(ehi).to_vertex();
        fixed_vt_prop[vhf]++;
        fixed_vt_prop[vht]++;

        ic_sge_prop[vhf] = ehi;
        ic_sge_prop[vht] = ehi;
      }
    }

    //fix boundary vertex
    for (const auto &fhi: _or_mesh.faces())
    {
      if (ff_prop[fhi])
      {
        for (auto fv_it = _or_mesh.fv_iter(fhi); fv_it.valid(); ++fv_it)
        {
          if (fixed_vt_prop[*fv_it] == 0)
            fixed_vt_prop[*fv_it] = -1;
        }
      }
    }

    //collect alignment
    std::map<HFH, int> aligned_axis;
    for (const auto &chi: _or_mesh.cells())
    {
      for (auto chf_it = _or_mesh.chf_iter(chi); chf_it.valid(); ++chf_it)
      {
        auto hfh_opp = _or_mesh.opposite_halfface_handle(*chf_it);
        if (ff_prop[_or_mesh.face_handle(hfh_opp)])
        {
          auto nm = _or_mesh.normal(hfh_opp);
          int axis = AxisAlignmentHelpers::closest_axis(or_cell_quaternions[chi],
                                                        nm).second;

          aligned_axis[hfh_opp] = axis;
        }
      }
    }

    VertexOptimizeT vo(_or_mesh);
    vo.set_energy_type(VertexOptimizeT<MeshT>::WORST_QUALITY);

    for (int i = 0; i < _opt_iter; ++i)
    {
      for (const auto &vhi: _or_mesh.vertices())
      {
        if (!bv_prop[vhi])
        {
          if (fixed_vt_prop[vhi] == 0)
            vo.interior_vertex_optimize(vhi);
          else if (fixed_vt_prop[vhi] == 1)
          {
            EH ehi = ic_sge_prop[vhi];
            auto pt_s = _or_mesh.vertex(_or_mesh.edge(ehi).from_vertex());
            auto pt_e = _or_mesh.vertex(_or_mesh.edge(ehi).to_vertex());

            vo.optimize_vertex_on_edge(vhi, pt_s, pt_e, 20);
          }
        }
        else
        {
          if (fixed_vt_prop[vhi] == 1)
          {
            EH ehi = ic_sge_prop[vhi];
            if (be_prop[ehi])
            {
              auto pt_s = _or_mesh.vertex(_or_mesh.edge(ehi).from_vertex());
              auto pt_e = _or_mesh.vertex(_or_mesh.edge(ehi).to_vertex());

              vo.optimize_vertex_on_edge(vhi, pt_s, pt_e, 20);
            }
          }
        }
      }
    }

    std::map<CH, std::vector<std::pair<Vec3d, int>>> alignments;
    std::set<FH> bdy_fhs;
    for (const auto &fhi: _or_mesh.faces())
      if (_or_mesh.is_boundary(fhi))
        bdy_fhs.insert(fhi);

    //align to face
    for (auto&[hfi, ax]: aligned_axis)
    {
      CH ch = _or_mesh.incident_cell(_or_mesh.opposite_halfface_handle(hfi));

      auto nm = _or_mesh.normal(hfi);

      alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), ax);
    }

    //align to feature and singular edges
    for (const auto chi: _or_mesh.cells())
    {
      if (chax[chi].first.is_valid())
      {
        auto dir = _or_mesh.vertex(_or_mesh.halfedge(chax[chi].first).to_vertex()) -
                   _or_mesh.vertex(_or_mesh.halfedge(chax[chi].first).from_vertex());
        dir.normalize();

        alignments[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), chax[chi].second);
      }
    }

    for (int i = 0; i < 10; ++i)
      QuaternionSmoothing::smooth_field_quaternion(_or_mesh, or_trans_prop, _or_mesh.cells(), bdy_fhs, alignments, tq_,
                                                   or_cell_quaternions);

    ALGOHEX_DEBUG_ONLY(check_mesh_quality(_or_mesh);
                               std::cerr << std::endl;)
  }

  bool is_locally_meshable_edge(const EH _eh, const bool _align_sg = false)
  {
    bool valid = true;

    if (!feature_face_edge_[_eh])
    {
      if (valence_[_eh] != std::numeric_limits<int>::max())
        valid = true;
      else
        valid = false;
    }
    else
    {
      //store old quaterions
      std::map<CH, Quaternion> old_cell_quaternions;
      if (_align_sg)
      {
        std::set<CH> chs;
        for (auto ec_it = mesh_.ec_iter(_eh); ec_it.valid(); ++ec_it)
          old_cell_quaternions[*ec_it] = cell_quaternions_[*ec_it];

        QuaternionSmoothing::optimize_quaternions_wrt_prescribed_valence(mesh_, trans_prop_, tq_, cell_quaternions_,
                                                                         valence_, feature_edge_, feature_fprop_, chs,
                                                                         true);
      }

      std::vector<HFH> ft_hfhs;
      HEH heh0 = mesh_.halfedge_handle(_eh, 0);
      for (auto hehf_it = mesh_.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
      {
        if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
        {
          ft_hfhs.push_back(*hehf_it);
        }
      }

      int fsec_agl = 0;
      if (!mesh_.is_boundary(ft_hfhs[0]))
      {
        fsec_agl = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                       trans_prop_, heh0, ft_hfhs[0],
                                                                                       ft_hfhs[1]);

        if (fsec_agl == 0)
        {
          valid = false;
        }
//                    da1 = dihedral_angle_from_halfface_to_halfface(mesh_, heh0, ft_hfhs[0], ft_hfhs[1]);

      }
      int fsec_agl2 = 0;
      if (!mesh_.is_boundary(ft_hfhs[1]))
      {
        fsec_agl2 = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                        trans_prop_, heh0, ft_hfhs[1],
                                                                                        ft_hfhs[0]);

        if (fsec_agl2 == 0)
        {
          valid = false;
        }
//                    da2 = dihedral_angle_from_halfface_to_halfface(mesh_, heh0, ft_hfhs[1], ft_hfhs[0]);
      }

      int sum = fsec_agl + fsec_agl2;

      int target_agl = valence_[_eh] + 4;
      if (mesh_.is_boundary(_eh))
        target_agl -= 2;

      if (sum != target_agl)
      {
        valid = false;
      }

      if (_align_sg)
      {
        //recover quaternions
        for (auto ec_it = mesh_.ec_iter(_eh); ec_it.valid(); ++ec_it)
          cell_quaternions_[*ec_it] = old_cell_quaternions[*ec_it];
      }
    }

    return valid;
  }

  bool is_locally_meshable(const VH _vh, const bool _align_sg = false, const int _n_split = 1)
  {
    //store old quaterions
    std::map<CH, Quaternion> old_cell_quaternions;
    if (_align_sg)
    {
      for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
        old_cell_quaternions[*vc_it] = cell_quaternions_[*vc_it];

      QuaternionSmoothing::optimize_quaternions_wrt_prescribed_valence(mesh_, trans_prop_, tq_, cell_quaternions_,
                                                                       valence_, feature_edge_, feature_fprop_, _vh,
                                                                       true);
    }

    //check one ring parameterization
    bool valid = false;
    MeshT onering;
    DuplicateOneRingMeshT<MeshT> dp(mesh_, onering);
    dp.copy_one_ring(std::vector<VH>{_vh});

    if (dp.check_cell_singularity())
    {
      if (verbose_)
      {
        std::cerr << "checking vertex: " << _vh << " one ring vh " << dp.onering_mesh_vertex_handle(_vh)
                  << std::endl;
        check_mesh_quality(onering);
      }

      VH vh_or = dp.onering_mesh_vertex_handle(_vh);
      // cutout ball where all edges are of unit length
      normalize_one_ring_mesh(onering, vh_or);
      // remove all constraints that are not incident to vh_or
      remove_locally_irrelevant_constraints(onering, vh_or);

      valid = check_parameterization(onering, vh_or, _vh.idx(), save_locally_non_meshable_);

      // alwasy verbose for non-meshable!
      if (!valid && !verbose_)
      {
        verbose_ = true;
        check_parameterization(onering, vh_or, _vh.idx(), false);
        verbose_ = false;
      }

      if (!valid)
      {
        dp.get_halfedge_axis_in_cell_in_onering_mesh();
        optimize_onering_mesh(onering, vh_or, 5);
        check_parameterization(onering, vh_or, _vh.idx(), false);
        verbose_ = true;
      }

      int n_iter = 0;
//              {
//                // write file
//                std::stringstream ss;
//                ss << "mesh_split_" << n_iter << ".ovm";
//
//                OpenVolumeMesh::IO::FileManager fm;
//                fm.writeFile(ss.str(), onering);
//              }

      while (!valid && n_iter < _n_split)
      {
        ++n_iter;
        std::cerr << "################ Refinement Iter " << n_iter << " with #vertices before = "
                  << onering.n_vertices();
        split_edges(onering, vh_or);
        normalize_one_ring_mesh(onering, vh_or);

//                // write file
//                  std::stringstream ss;
//                  ss << "mesh_split_" << n_iter << ".ovm";
//
//                  OpenVolumeMesh::IO::FileManager fm;
//                  fm.writeFile(ss.str(), onering);

        std::cerr << " and #vertices after = " << onering.n_vertices() << std::endl;
        valid = check_parameterization(onering, vh_or, _vh.idx(), false);
      }
      verbose_ = false;


//                if (!valid) {
//                    dp.get_halfedge_axis_in_cell_in_onering_mesh();
//
//                    n_opt_++;
//                    optimize_onering_mesh(onering, vh_or, 10);
//                    valid = check_parameterization(onering, vh_or, _vh.idx(), true);
//
//                    std::cerr << " is valid? " << valid << std::endl;
//                    if (!valid)
//                        check_mesh_dihedral_angle(onering);
//                    else
//                        n_fixed_++;
//                }
    }

    if (_align_sg)
    {
      //recover quaternions
      for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
        cell_quaternions_[*vc_it] = old_cell_quaternions[*vc_it];
    }

    return valid;
  }

  bool check_parameterization(MeshT &_mesh, const VH _vh, const int _orig_vh_id, bool _save_config)
  {
    //check parameterization
    AlgoHex::FrameFieldOptimizer3DT<MeshT> ffopt(_mesh);
    ffopt.import_frames_from_quaternion_property_onering(_vh);
    if (verbose_)
    {
      ffopt.check_frame_rotation_angles();
      ffopt.check_valence_consistency();
    }

    MeshT export_mesh;
    bool is_valid = ffopt.is_locally_meshable(_vh, export_mesh, verbose_);

    if (_save_config && !is_valid && save_locally_non_meshable_)
    {
      // export local configuration
      std::stringstream ss;
      ss << save_locally_non_meshable_filename_base_ << "_nlm_vertex_" << _orig_vh_id << ".ovm";

      // write file
      OpenVolumeMesh::IO::FileManager fm;
      fm.writeFile(ss.str(), export_mesh);
    }

    return is_valid;
  }


  double check_local_meshability(const bool _align_to_sge = false)
  {
//            Debug::ScopedOutputLevel sol(0);

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

    // check all vertices
    for (VIt v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
    {
      // global counter
      ++n;
      bool ilm = is_locally_meshable(*v_it, _align_to_sge);
      if (ilm)
        ++n_meshable;

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

      // count zipper nodes
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
    std::cerr << "Number of optimized local mesh: " << n_opt_ << " number of fixed: " << n_fixed_ << std::endl;

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
    std::cerr << "#zipper node  = " << std::setw(5) << n_zipper_node << ", #meshable = " << std::setw(5)
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
    json_data_["zipper nodes"] = n_zipper_node;
    json_data_["singular vertices"] = n_on_singular_vertex - n_on_singular_vertex_meshable;
    json_data_["feature vertices"] = n_on_feature_edge - n_on_feature_edge_meshable;
    json_data_["complex singular edges"] = len_complex_singular_edges();
    json_data_["singular edges"] = len_singular_edges();
    // return percentage of meshable vertices
    return double(n_meshable) / double(n);
  }

  bool is_on_feature_vertex(const VH _vh) const
  {
    // is feature vertex?
    return feature_node_[_vh] > 0;
  }


  bool is_on_feature_edge(const VH _vh) const
  {
    // is on feature edge?
    VOHEIt vhe_it(_vh, &mesh_);
    for (; vhe_it.valid(); ++vhe_it)
      if (feature_edge_[mesh_.edge_handle(*vhe_it)] > 0)
        return true;

    return false;
  }


  bool is_on_feature_face(const VH _vh) const
  {
    // is on feature triangle?
    VFIt vf_it(_vh, &mesh_);
    for (; vf_it.valid(); ++vf_it)
      if (feature_fprop_[*vf_it] > 0)
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

  int split_edges(MeshT &_or_mesh, const VH _vh)
  {
    std::cout << "Splitting edges ..." << std::endl;
    EdgeSplitT<MeshT> es(_or_mesh);

    auto ve_prop = _or_mesh.template request_edge_property<int>("vertex incident edges", false);
    for (auto ve_it = _or_mesh.ve_iter(_vh); ve_it.valid(); ++ve_it)
      ve_prop[*ve_it] = true;

    std::set<PairDE, std::greater<PairDE>> ehs_split;
    double avr_elen = 0.;
    for (const auto ehi: _or_mesh.edges())
      if (!ve_prop[ehi])
      {
        double len = _or_mesh.length(ehi);
        ehs_split.emplace(len, ehi);
        avr_elen += len;
      }

    avr_elen /= (double) ehs_split.size();

    int n = 0;
    while (!ehs_split.empty() && n < 1000)
    {
      auto de_cur = *ehs_split.begin();
      ehs_split.erase(ehs_split.begin());

      if (es.is_split_ok(de_cur.second))
      {
        if (de_cur.first != _or_mesh.length(de_cur.second))
          continue;

        auto vh_f = _or_mesh.edge(de_cur.second).from_vertex();
        auto vh_t = _or_mesh.edge(de_cur.second).to_vertex();

//            if(de_cur.first < 4./3.*avr_elen)
//              continue;

        auto vh_mid = es.edge_split(de_cur.second);
        if (vh_mid.is_valid())
        {
          auto eh_mid = _or_mesh.edge_handle(_or_mesh.find_halfedge(_vh, vh_mid));
          ve_prop[eh_mid] = true;

          n++;

//              //push
//              for(auto ve_it = _or_mesh.ve_iter(vh_mid); ve_it.valid(); ++ve_it)
//                if(!ve_prop[*ve_it])
//                  ehs_split.emplace(_or_mesh.length(*ve_it), *ve_it);
        }
      }
    }

    std::cout << "Split " << n << " edges." << std::endl;
    _or_mesh.collect_garbage();

    return n;
  }

  void normalize_one_ring_mesh(MeshT &_mesh, const VH _vh_center) const
  {
    // translate center vertex to origin
    auto p = _mesh.vertex(_vh_center);

    for (auto vh: _mesh.vertices())
      _mesh.set_vertex(vh, _mesh.vertex(vh) - p);

    // scale all edges to unit length
    for (auto vh: _mesh.vertices())
      if (vh != _vh_center)
        _mesh.set_vertex(vh, _mesh.vertex(vh) / _mesh.vertex(vh).norm());
  }


  void remove_locally_irrelevant_constraints(MeshT &_mesh, const VH _vh_center) const
  {
    auto fvp(_mesh.template request_vertex_property<int>("AlgoHex::FeatureVertices"));
    auto fep(_mesh.template request_edge_property<int>("AlgoHex::FeatureEdges"));
    auto ffp(_mesh.template request_face_property<int>("AlgoHex::FeatureFaces"));
    auto evp(_mesh.template request_edge_property<int>("edge_valance"));


    for (auto vh: _mesh.vertices())
      if (vh != _vh_center)
        fvp[vh] = false;

    for (auto eh: _mesh.edges())
    {
      bool is_incident = false;
      auto evhs = _mesh.edge_vertices(eh);
      for (auto vh: evhs)
        if (vh == _vh_center)
          is_incident = true;

      if (!is_incident)
      {
        fep[eh] = false;
        evp[eh] = 0;
      }
    }

    for (auto fh: _mesh.faces())
    {
      bool is_incident = false;
      auto fvhs = _mesh.get_halfface_vertices(_mesh.halfface_handle(fh, 0));
      for (auto vh: fvhs)
        if (vh == _vh_center)
          is_incident = true;

      if (!is_incident)
        ffp[fh] = false;
    }
  }

  bool &verbose() { return verbose_; }

private:
  MeshT &mesh_;
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
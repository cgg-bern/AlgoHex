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
#include "DuplicateOneRingMeshT.hh"
#include "../Parametrization3DT.hh"
#include <AlgoHex/MeshGeometry.hh>
#include <Base/Debug/DebConfig.hh>
#include "TetRemeshingT.hh"
#include "../SingularGraphExtractionT.hh"
#include "QuaternionsSmoothing.hh"


namespace AlgoHex
{
template<class MeshT>
class OneringParameterizationChecker : public MeshPropertiesT<MeshT>
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

  using PairDE = std::pair<double, EH>;

  OneringParameterizationChecker(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                                                 mesh_(_mesh) {}

  ~OneringParameterizationChecker() = default;


public:
  nlohmann::json &json_data() { return json_data_; }

  int num_invalid_singular_node() const
  {
    int n_invalid = 0, n_checked = 0;
    for (const auto vhi: mesh_.vertices())
    {
      int node_id = MeshPropertiesT<MeshT>::node_index(vhi);
      if (node_id != 0)
      {
        if (!is_node_onering_parameterization_valid(vhi))
        {
          n_invalid++;
        }

        n_checked++;
      }
    }
    std::cout << "Number of invalid singular nodes: " << n_invalid << " of " << n_checked << std::endl;

    return n_invalid;
  }

  void optimize_onering_mesh(MeshT &_or_mesh, const VH _vh, int _opt_iter = 10)
  {
    auto chax = _or_mesh.template request_cell_property<std::pair<HEH, int> >("EdgeAxis");
//            std::cerr<<"cax before split ";
//            for(auto chi : _or_mesh.cells())
//                if(chax[chi].first.is_valid())
//                    std::cerr<<"chi "<<chi<<" he "<<_or_mesh.edge_handle(chax[chi].first)<<" ax "<<chax[chi].second<<std::endl;

//            RemedyFrameMatchingInconsistencyT<MeshT> rfmi(_or_mesh);
//            rfmi.remedy_faces();
//            std::cerr<<"cax after split ";
//            for(auto chi : _or_mesh.cells())
//                if(chax[chi].first.is_valid())
//                    std::cerr<<"chi "<<chi<<" he "<<_or_mesh.edge_handle(chax[chi].first)<<" ax "<<chax[chi].second<<std::endl;
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

//            //movable boundary feature vertex is -2
//            for (const auto &ehi : _or_mesh.edges()) {
//                if (fe_prop[ehi]) {
//                    auto vhf = _or_mesh.edge(ehi).from_vertex();
//                    auto vht = _or_mesh.edge(ehi).to_vertex();
//                    fixed_vt_prop[vhf]--;
//                    fixed_vt_prop[vht]--;
//
//                    ic_sge_prop[vhf] = ehi;
//                    ic_sge_prop[vht] = ehi;
//                }
//            }

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
//            std::set<CH> visited_chs;
//            for (const auto &chi : _or_mesh.cells()) {
//                if(visited_chs.find(chi) != visited_chs.end())
//                    continue;
//
//                for (auto ce_it = _or_mesh.ce_iter(chi); ce_it.valid(); ++ce_it) {
//                    if (std::abs(or_valence[*ce_it]) == 1 ||
//                        std::abs(or_valence[*ce_it]) == 2) {
//                        HEH heh0 = _or_mesh.halfedge_handle(*ce_it, 0);
//                        auto dir = _or_mesh.vertex(_or_mesh.halfedge(heh0).to_vertex()) -
//                                   _or_mesh.vertex(_or_mesh.halfedge(heh0).from_vertex());
//                        dir.normalize();
//
//                        auto hehf_it = _or_mesh.hehf_iter(heh0);
//                        CH ch_s(-1);
//                        bool e_bdy = _or_mesh.is_boundary(*ce_it);
//                        if(!e_bdy)
//                            ch_s = _or_mesh.incident_cell(_or_mesh.opposite_halfface_handle(*hehf_it));
//                        else
//                            ch_s = _or_mesh.incident_cell(*hehf_it);
//
//                        int axis = halfedge_rotational_axis_idx(_or_mesh, tq_, or_cell_quaternions, or_trans_prop, or_valence, heh0, *hehf_it);
//                        if (or_valence[*ce_it] == 1)
//                            axis = axis % 2 == 0 ? axis + 1 : axis - 1;
//
//                        alignments[ch_s].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
//                        visited_chs.insert(ch_s);
//
//                        if(e_bdy)
//                            ++hehf_it;
//
//                        for (; hehf_it.valid(); ++hehf_it) {
//                            auto chi = _or_mesh.incident_cell(*hehf_it);
//                            if (chi.is_valid()) {
//                                axis = tq_.axis_after_transition(axis, or_trans_prop[*hehf_it]);
//
//                                alignments[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
//                                visited_chs.insert(chi);
//                            }
//                        }
//                    } else if (fe_prop[*ce_it]) {
//                        auto dir = _or_mesh.vertex(_or_mesh.edge(*ce_it).to_vertex()) -
//                                   _or_mesh.vertex(_or_mesh.edge(*ce_it).from_vertex());
//                        dir.normalize();
//
//                        int axis = AxisAlignmentHelpers::closest_axis(or_cell_quaternions[chi],
//                                                                                 dir).second;
//                        alignments[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
//                        visited_chs.insert(chi);
//                    }
//                }
//            }

    for (int i = 0; i < 10; ++i)
      QuaternionSmoothing::smooth_field_quaternion(_or_mesh, or_trans_prop, _or_mesh.cells(), bdy_fhs, alignments, tq_,
                                                   or_cell_quaternions);
//            }

    ALGOHEX_DEBUG_ONLY(check_mesh_quality(_or_mesh);
                               std::cerr << std::endl;)
  }

  int split_edges(MeshT &_or_mesh, const VH _vh)
  {
    ALGOHEX_DEBUG_ONLY(std::cout << "Splitting edges ..." << std::endl;)
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

        if (de_cur.first < 4. / 3. * avr_elen)
          continue;

        auto vh_mid = es.edge_split(de_cur.second);
        if (vh_mid.is_valid())
        {
          auto eh_mid = _or_mesh.edge_handle(_or_mesh.halfedge(_vh, vh_mid));
          ve_prop[eh_mid] = true;

          n++;

          //push
          for (auto ve_it = _or_mesh.ve_iter(vh_mid); ve_it.valid(); ++ve_it)
            if (!ve_prop[*ve_it])
              ehs_split.emplace(_or_mesh.length(*ve_it), *ve_it);
        }
      }
    }

    std::cout << "Split " << n << " edges." << std::endl;
    _or_mesh.collect_garbage();


    return n;
  }

  bool is_locally_meshable(const VH _vh, const bool _align_sg = false)
  {
    //store old quaterions
    std::map<CH, Quaternion> old_cell_quaternions;
    if (_align_sg)
    {
      for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
        old_cell_quaternions[*vc_it] = cell_quaternions_[*vc_it];

      //align quaternions to singular or feature edges or boundary faces
//                QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_, feature_fprop_, feature_edge_, _vh, true);
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
      ALGOHEX_DEBUG_ONLY(
              std::cerr << "checking vertex: " << _vh << " one ring vh " << dp.onering_mesh_vertex_handle(_vh)
                        << std::endl;)
      ALGOHEX_DEBUG_ONLY(check_mesh_quality(onering);)

//            std::cerr<<"before: "<<std::endl;
      auto or_trans_prop = onering.template request_halfface_property<int>("HalffaceTransiton");
      auto or_cell_quaternions = onering.template request_cell_property<Eigen::Quaterniond>(
              "FrameFieldQuaternions");

//            for(auto fhi : onering.faces())
//                std::cerr<<" hfh "<<onering.halfface_handle(fhi, 0)<<" trans "<<or_trans_prop[onering.halfface_handle(fhi, 0)]<<std::endl;
//
//            for(auto chi : onering.cells())
//                std::cerr<<" ch "<<chi<<" qtn "<<or_cell_quaternions[chi].coeffs()<<std::endl;

      VH vh_or = dp.onering_mesh_vertex_handle(_vh);
      valid = check_parameterization(onering, vh_or);


      if (!valid)
      {
        dp.get_halfedge_axis_in_cell_in_onering_mesh();


        optimize_onering_mesh(onering, vh_or, 10);
        valid = check_parameterization(onering, vh_or);

        ALGOHEX_DEBUG_ONLY(std::cerr << " is valid? " << valid << std::endl;)
        ALGOHEX_DEBUG_ONLY(if (!valid)
                             check_mesh_dihedral_angle(onering);)
      }
    }


    if (!valid && save_locally_non_meshable_)
    {
      // export local configuration
      std::stringstream ss;
      ss << save_locally_non_meshable_filename_base_ << "_nlm_vertex_" << _vh.idx() << ".ovm";

      // write file
      OpenVolumeMesh::IO::FileManager fm;
      fm.writeFile(ss.str(), onering);
    }

    if (_align_sg)
    {
      //recover quaternions
      for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
        cell_quaternions_[*vc_it] = old_cell_quaternions[*vc_it];
    }

    return valid;
  }

  bool is_locally_meshable2(MeshT &_onering, const VH _vh, const int _opt_iters, const bool _align_sg = false)
  {
    //store old quaterions
    std::map<CH, Quaternion> old_cell_quaternions;
    if (_align_sg)
    {
      for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
        old_cell_quaternions[*vc_it] = cell_quaternions_[*vc_it];

      //align quaternions to singular or feature edges or boundary faces
//                QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_, feature_fprop_, feature_edge_, _vh, true);
      QuaternionSmoothing::optimize_quaternions_wrt_prescribed_valence(mesh_, trans_prop_, tq_, cell_quaternions_,
                                                                       valence_, feature_edge_, feature_fprop_, _vh,
                                                                       true);
    }

    //check one ring parameterization
    DuplicateOneRingMeshT<MeshT> dp(mesh_, _onering);
    dp.copy_one_ring(std::vector<VH>{_vh});
    if (!dp.check_cell_singularity())
      return false;

    ALGOHEX_DEBUG_ONLY(std::cerr << "checking vertex: " << _vh << " one ring vh " << dp.onering_mesh_vertex_handle(_vh)
                                 << std::endl;)
    ALGOHEX_DEBUG_ONLY(check_mesh_quality(_onering);)

//            std::cerr<<"before: "<<std::endl;
    auto or_trans_prop = _onering.template request_halfface_property<int>("HalffaceTransiton");
    auto or_cell_quaternions = _onering.template request_cell_property<Eigen::Quaterniond>("FrameFieldQuaternions");

//            for(auto fhi : onering.faces())
//                std::cerr<<" hfh "<<onering.halfface_handle(fhi, 0)<<" trans "<<or_trans_prop[onering.halfface_handle(fhi, 0)]<<std::endl;
//
//            for(auto chi : onering.cells())
//                std::cerr<<" ch "<<chi<<" qtn "<<or_cell_quaternions[chi].coeffs()<<std::endl;

    VH vh_or = dp.onering_mesh_vertex_handle(_vh);
    bool valid = check_parameterization(_onering, vh_or);

    if (!valid)
    {
      dp.get_halfedge_axis_in_cell_in_onering_mesh();
      optimize_onering_mesh(_onering, vh_or, _opt_iters);
      valid = check_parameterization(_onering, vh_or);

      ALGOHEX_DEBUG_ONLY(std::cerr << " is valid? " << valid << std::endl;)
      ALGOHEX_DEBUG_ONLY(if (!valid)
                           check_mesh_dihedral_angle(_onering);)
    }


    if (!valid && save_locally_non_meshable_)
    {
      // export local configuration
      std::stringstream ss;
      ss << save_locally_non_meshable_filename_base_ << "_nlm_vertex_" << _vh.idx() << ".ovm";

      // write file
      OpenVolumeMesh::IO::FileManager fm;
      fm.writeFile(ss.str(), _onering);
    }

    if (_align_sg)
    {
      //recover quaternions
      for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
        cell_quaternions_[*vc_it] = old_cell_quaternions[*vc_it];
    }

    return valid;
  }

  bool check_parameterization(MeshT &_mesh, const VH _vh)
  {
    auto or_valence = _mesh.template request_edge_property<int>("edge_valance");
    auto or_trans_prop = _mesh.template request_halfface_property<int>("HalffaceTransiton");
    //store old matchings
    std::map<HFH, int> old_trans;
    for (const auto hfhi: _mesh.halffaces())
      old_trans[hfhi] = or_trans_prop[hfhi];

    //check parameterization
    OpenVolumeMesh::StatusAttrib status(_mesh);

    Parametrization3DT<MeshT> pm(_mesh, status);
    pm.import_frames_from_quaternion_field();

    pm.set_max_stiffening_iters(10);

    pm.parametrize(1000, 1.0, false);

    int n_invalid = pm.number_of_invalid_parametric_tets();
    if (n_invalid > 0)
      std::cerr << "Vertex: " << _vh << " has " << n_invalid << " invalid cells!" << std::endl;
    //check the singular edge angle in parameterization domain
    bool angle_ok = true;
    for (auto voh_it = _mesh.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      EH eh = _mesh.edge_handle(*voh_it);
      int val = or_valence[eh];
      double ref = 2 * M_PI;
      if (val >= -2 && val <= 4 && val != 0)
        ref += M_PI_2 * val;
      else if (val != 0)
        ALGOHEX_DEBUG_ONLY(
                std::cerr << "Warning: complex singular edge of valence " << val << " in one ring!" << std::endl;);

      double ea = pm.edge_angle(_mesh.edge_handle(*voh_it));
      if (_mesh.is_boundary(*voh_it))
        ref -= M_PI;

      bool angle_i_ok = std::abs(ref - ea) < 0.035;

      if (!angle_i_ok)
      {
        ALGOHEX_DEBUG_ONLY(
                std::cerr << " edge: " << _mesh.edge_handle(*voh_it) << " of valence: " << val << " edge angle: "
                          << ea << " target angle: " << ref << std::endl;);

        angle_ok = false;
      }
    }

    bool valid = n_invalid == 0 && angle_ok;
    //restore matchings because pm combed the field
    if (!valid)
    {
      for (const auto hfhi: _mesh.halffaces())
        or_trans_prop[hfhi] = old_trans[hfhi];
    }

    return valid;
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

  void enable_save_locally_non_meshable(const std::string _filename_base = "locally_non_meshable_")
  {
    save_locally_non_meshable_ = true;
    save_locally_non_meshable_filename_base_ = _filename_base;
  }


private:
  MeshT &mesh_;
  AlgoHex::TransitionQuaternion tq_;

  bool save_locally_non_meshable_ = false;
  std::string save_locally_non_meshable_filename_base_;

  // json storage
  nlohmann::json json_data_;
};
}

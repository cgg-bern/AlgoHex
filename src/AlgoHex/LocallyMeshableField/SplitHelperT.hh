/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "MeshProperties.hh"
#include "TetRemeshingT.hh"


namespace AlgoHex
{
template<class MeshT>
class SplitHelperT : public MeshPropertiesT<MeshT>
{
public:
  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::feature_edge_vertex_;
  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;
  using MeshPropertiesT<MeshT>::feature_face_vertex_;

  SplitHelperT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                               mesh_(_mesh),
                               es_(mesh_),
                               fs_(mesh_),
                               vo_(mesh_) {}

  //
  void split_for_single_alignment_feature_face()
  {
    std::set<FH> fhs;
    for (const auto fh: mesh_.faces())
    {
      if (feature_fprop_[fh] == 0)
      {//split non-feature face which has three feature face vertices
        int count = 0;
        for (auto fv_it = mesh_.fv_iter(fh); fv_it.valid(); ++fv_it)
        {
          if (feature_face_vertex_[*fv_it])
            count++;
          else
            break;
        }

        if (count == 3)
          fhs.insert(fh);
      }
      else
      {//split feature face which has >= 2 feature edges
        int count = 0;
        for (auto fe_it = mesh_.fe_iter(fh); fe_it.valid(); ++fe_it)
        {
          if (feature_edge_[*fe_it] > 0)
            count++;
        }

        if (count >= 2)
          fhs.insert(fh);
      }
    }

    for (const auto fhi: fhs)
      if (!mesh_.is_deleted(fhi))
        fs_.face_split(fhi);

    //split non-feature face edge that has two feature face vertices
    std::vector<EH> ehs3;
    for (const auto eh: mesh_.edges())
    {
      if (!feature_face_edge_[eh])
      {
        VH vh_f = mesh_.edge(eh).from_vertex();
        VH vh_t = mesh_.edge(eh).to_vertex();
        if (feature_face_vertex_[vh_f] && feature_face_vertex_[vh_t])
          ehs3.push_back(eh);
      }
    }

    for (const auto ehi: ehs3)
      if (!mesh_.is_deleted(ehi))
      {
        es_.edge_split(ehi);
      }


    mesh_.collect_garbage();
  }


  bool is_edge_to_split_fe(const EH _eh, const bool _skip_two_sgv = false)
  {
    VH vh0 = mesh_.edge(_eh).from_vertex();
    VH vh1 = mesh_.edge(_eh).to_vertex();

    if (valence_[_eh] == 0 && feature_edge_[_eh] == 0)
    {
      if ((feature_edge_vertex_[vh1] && feature_edge_vertex_[vh0]) ||
          (sgl_vt_[vh1] && feature_edge_vertex_[vh0]) ||
          (feature_edge_vertex_[vh1] && sgl_vt_[vh0]))
        return true;

      if (!_skip_two_sgv && sgl_vt_[vh1] && sgl_vt_[vh0])
        return true;
    }
    else if (valence_[_eh] == 0 && feature_edge_[_eh] > 0)
    {
      if (!_skip_two_sgv && sgl_vt_[vh1] && sgl_vt_[vh0])
        return true;
    }
    else if (valence_[_eh] != 0 && feature_edge_[_eh] == 0)
    {
      if (feature_edge_vertex_[vh1] && feature_edge_vertex_[vh0])
        return true;
    }
    else
    {//singular or feature edges which has two special nodes
      if (n_incident_special_edges(vh1) > 2 && n_incident_special_edges(vh0) > 2)
        return true;
    }

    return false;
  }

  bool is_edge_to_split_ffe(const EH _eh, const bool _skip_two_sgv = false)
  {
    VH vh0 = mesh_.edge(_eh).from_vertex();
    VH vh1 = mesh_.edge(_eh).to_vertex();

    if (valence_[_eh] == 0 && !feature_face_edge_[_eh])
    {
      if ((feature_face_vertex_[vh1] && feature_face_vertex_[vh0]) ||
          (sgl_vt_[vh1] && feature_face_vertex_[vh0]) ||
          (feature_face_vertex_[vh1] && sgl_vt_[vh0]))
        return true;

      if (!_skip_two_sgv && sgl_vt_[vh1] && sgl_vt_[vh0])
        return true;
    }
    else if (valence_[_eh] == 0 && feature_face_edge_[_eh])
    {
      if (!_skip_two_sgv && sgl_vt_[vh1] && sgl_vt_[vh0])
        return true;
    }
    else if (valence_[_eh] != 0 && !feature_face_edge_[_eh])
    {
      if (feature_face_vertex_[vh1] && feature_face_vertex_[vh0])
        return true;
    }
    else
    {//singular or feature edges which has two special nodes
      if (n_incident_special_edges(vh1) > 2 && n_incident_special_edges(vh0) > 2)
        return true;
    }

    return false;
  }

  //it doesn't split regular edge with two singular vertices if true
  void split_for_single_alignment(const bool _skip_two_sgv = false)
  {
    //split faces with three edges being special
    std::set<FH> fhs_split;
    for (const auto fhi: mesh_.faces())
    {
      int n_se = 0;
      for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
        if (valence_[*fe_it] != 0 || feature_edge_[*fe_it] > 0)
          n_se++;

      if (n_se == 3)
        fhs_split.insert(fhi);
    }

    for (const auto fhi: fhs_split)
      if (fs_.is_split_ok(fhi))
        fs_.face_split(fhi);

    split_for_single_alignment_feature_face();

    std::vector<EH> ehs1;
    //find feature face edges to split
    for (const auto ehi: mesh_.edges())
      if (is_edge_to_split_fe(ehi, _skip_two_sgv))
        ehs1.push_back(ehi);

    for (const auto ehi: ehs1)
      if (!mesh_.is_deleted(ehi))
      {
        es_.edge_split(ehi);
      }

    //find interior edges to split
    ehs1.clear();
    for (const auto ehi: mesh_.edges())
      if (is_edge_to_split_ffe(ehi, _skip_two_sgv))
        ehs1.push_back(ehi);

    for (const auto ehi: ehs1)
      if (!mesh_.is_deleted(ehi))
      {
        es_.edge_split(ehi);
      }

    mesh_.collect_garbage();
  }

  void split_for_good_tet_quality_at_singular_edges()
  {
    double max_da = -10;
    std::vector<EH> ehs;

    //split cells incident to singular edges with large dihedral angle
    for (const auto ehi: mesh_.edges())
    {
      if (valence_[ehi] != 0)
      {
        auto heh0 = mesh_.halfedge_handle(ehi, 0);
        VH vh0 = mesh_.edge(ehi).from_vertex();
        VH vh1 = mesh_.edge(ehi).to_vertex();
        auto pt0 = mesh_.vertex(vh0);
        auto pt1 = mesh_.vertex(vh1);

        for (auto hehf_it = mesh_.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
        {
          CH chi = mesh_.incident_cell(*hehf_it);
          if (chi.is_valid())
          {
            auto n1 = mesh_.normal(*hehf_it);
            auto hf_adj = mesh_.adjacent_halfface_in_cell(*hehf_it, heh0);
            auto n0 = mesh_.normal(hf_adj);

            auto da = dihedral_angle(pt0, pt1, n0, n1);
            if (da > M_PI_2)
            {
              if (da > max_da)
                max_da = da;
              auto prev_heh = mesh_.prev_halfedge_in_halfface(heh0, *hehf_it);
              auto bt_hf = mesh_.adjacent_halfface_in_cell(*hehf_it, prev_heh);
              auto opposite_edge = mesh_.edge_handle(
                      mesh_.next_halfedge_in_halfface(mesh_.opposite_halfedge_handle(prev_heh), bt_hf));

              ehs.push_back(opposite_edge);
            }
          }
        }
      }
    }

    std::cerr << "Max dihedral angle in cells adjacent to singular edges: " << max_da / M_PI * 180. << " Will split: "
              << ehs.size() << " edges." << std::endl;

    //split and relocate new vertex
    for (const auto ehi: ehs)
      if (!mesh_.is_deleted(ehi))
      {
        auto new_vh = es_.edge_split(ehi);

        if (sgl_vt_[new_vh] == 0)
          vo_.vertex_optimize(new_vh);
      }

    mesh_.collect_garbage();
  }


  void split_special_edges_with_two_special_nodes()
  {
    std::vector<EH> ehs1;
    //find singular edges which have two singular or feature vertices
    for (const auto ehi: mesh_.edges())
      if (abs(valence_[ehi]) == 1 || feature_edge_[ehi] > 0)
      {
        VH vh0 = mesh_.edge(ehi).from_vertex();
        VH vh1 = mesh_.edge(ehi).to_vertex();
        if (n_incident_special_edges(vh0) > 2 && n_incident_special_edges(vh1) > 2)
          ehs1.push_back(ehi);
      }

    for (const auto ehi: ehs1)
      if (!mesh_.is_deleted(ehi))
      {
        es_.edge_split(ehi);
      }

    mesh_.collect_garbage();
  }

  //don't split singular edges
  static void
  split_one_ring(MeshT &mesh, EdgeSplitT<MeshT> &es, const VH _vh, EP<int> &valence, EP<bool> &feature_face_edge,
                 VP<bool> &feature_face_vertex, EP<int> &feature_edge, VP<int> &sgl_vt, VP<bool> &feature_edge_vertex,
                 const bool collect_garbage = false)
  {
    std::set<EH> ehs;
    for (auto vc_it = mesh.vc_iter(_vh); vc_it.valid(); ++vc_it)
    {
      for (auto ce_it = mesh.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
      {
        VH vh0 = mesh.edge(*ce_it).from_vertex();
        VH vh1 = mesh.edge(*ce_it).to_vertex();
        if (valence[*ce_it] == 0)
        {
          if (!feature_face_edge[*ce_it])
          {//split regular edges with two singular/ff vhs
            if ((sgl_vt[vh0] || feature_face_vertex[vh0]) &&
                (sgl_vt[vh1] || feature_face_vertex[vh1]))
              ehs.insert(*ce_it);
          }
          else
          {
            if (sgl_vt[vh1] && sgl_vt[vh0])
            {//split regular edges with two singular vhs
              ehs.insert(*ce_it);
            }

            if (feature_edge[*ce_it] == 0)
            {
              if (feature_edge_vertex[vh1] && feature_edge_vertex[vh0])
              {//split regular edges with two feature edge vhs
                ehs.insert(*ce_it);
              }
            }
          }
        }
      }
    }

    for (const auto &ehi: ehs)
      if (!mesh.is_deleted(ehi))
      {
        es.edge_split(ehi);
      }

    if (collect_garbage && !ehs.empty())
      mesh.collect_garbage();
  }

  //singular edges can be split
  static void split_one_ring_all(MeshT &mesh, EdgeSplitT<MeshT> &es, const VH _vh,
                                 EP<int> &valence, EP<bool> &feature_face_edge, VP<bool> &feature_face_vertex,
                                 EP<int> &feature_edge, VP<int> &sgl_vt, VP<bool> &feature_edge_vertex,
                                 const bool collect_garbage = false)
  {
    FaceSplitT<MeshT> fs(mesh);
    std::set<FH> fhs_split, fhs_all;
    for (auto vc_it = mesh.vc_iter(_vh); vc_it.valid(); ++vc_it)
    {
      for (auto cf_it = mesh.cf_iter(*vc_it); cf_it.valid(); ++cf_it)
        fhs_all.insert(*cf_it);
    }
    //split faces with three edges being special
    for (const auto fhi: fhs_all)
    {
      int n_se = 0;
      for (auto fe_it = mesh.fe_iter(fhi); fe_it.valid(); ++fe_it)
        if (valence[*fe_it] != 0 || feature_edge[*fe_it] > 0)
          n_se++;

      if (n_se == 3)
        fhs_split.insert(fhi);
    }

    for (const auto fhi: fhs_split)
      if (fs.is_split_ok(fhi))
        fs.face_split(fhi);

    //split edges
    std::set<EH> ehs;
    for (auto vc_it = mesh.vc_iter(_vh); vc_it.valid(); ++vc_it)
    {
      for (auto ce_it = mesh.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
      {
        VH vh0 = mesh.edge(*ce_it).from_vertex();
        VH vh1 = mesh.edge(*ce_it).to_vertex();
        if (valence[*ce_it] == 0)
        { //split regular edges with two singular vhs
          if (sgl_vt[vh1] && sgl_vt[vh0])
          {
            ehs.insert(*ce_it);
          }

          if (!feature_face_edge[*ce_it])
          {//split regular edges with two singular/ff vhs
            if ((sgl_vt[vh0] || feature_face_vertex[vh0]) &&
                (sgl_vt[vh1] || feature_face_vertex[vh1]))
              ehs.insert(*ce_it);
          }
        }
        else
        {
          if (!feature_face_edge[*ce_it])
          {//split non-ffe singular edges with two feature face vertex
            if (feature_face_vertex[vh0] && feature_face_vertex[vh1])
              ehs.insert(*ce_it);
          }
        }

        if (feature_edge[*ce_it] == 0)
        {//split non fe with two fe vhs
          if (feature_edge_vertex[vh1] && feature_edge_vertex[vh0])
          {
            ehs.insert(*ce_it);
          }
        }
      }
    }


    for (const auto &ehi: ehs)
    {
      if (!mesh.is_deleted(ehi))
      {
        es.edge_split(ehi);
      }
    }

    if (collect_garbage && !ehs.empty() && !fhs_split.empty())
      mesh.collect_garbage();
  }

  //singular edges can be split
  static void split_one_ring_all(MeshT &mesh, EdgeSplitT<MeshT> &es, const VH _vh,
                                 EP<int> &valence, EP<bool> &feature_face_edge, VP<bool> &feature_face_vertex,
                                 EP<int> &feature_edge, VP<int> &sgl_vt, VP<bool> &feature_edge_vertex,
                                 std::set<EH> &_excl_ehs, const bool collect_garbage = false)
  {
    FaceSplitT<MeshT> fs(mesh);
    std::set<FH> fhs_split, fhs_all;
    for (auto vc_it = mesh.vc_iter(_vh); vc_it.valid(); ++vc_it)
    {
      for (auto cf_it = mesh.cf_iter(*vc_it); cf_it.valid(); ++cf_it)
        fhs_all.insert(*cf_it);
    }
    //split faces with three edges being special
    for (const auto fhi: fhs_all)
    {
      int n_se = 0;
      for (auto fe_it = mesh.fe_iter(fhi); fe_it.valid(); ++fe_it)
        if (valence[*fe_it] != 0 || feature_edge[*fe_it] > 0)
          n_se++;

      if (n_se == 3)
        fhs_split.insert(fhi);
    }

    for (const auto fhi: fhs_split)
      if (fs.is_split_ok(fhi))
        fs.face_split(fhi);

    //split edges
    std::set<EH> ehs;
    for (auto vc_it = mesh.vc_iter(_vh); vc_it.valid(); ++vc_it)
    {
      for (auto ce_it = mesh.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
      {
        VH vh0 = mesh.edge(*ce_it).from_vertex();
        VH vh1 = mesh.edge(*ce_it).to_vertex();
        if (valence[*ce_it] == 0)
        { //split regular edges with two singular vhs
          if (sgl_vt[vh1] && sgl_vt[vh0])
          {
            ehs.insert(*ce_it);
          }

          if (!feature_face_edge[*ce_it])
          {//split regular edges with two singular/ff vhs
            if ((sgl_vt[vh0] || feature_face_vertex[vh0]) &&
                (sgl_vt[vh1] || feature_face_vertex[vh1]))
              ehs.insert(*ce_it);
          }
        }
        else
        {
          if (!feature_face_edge[*ce_it])
          {//split non-ffe singular edges with two feature face vertex
            if (feature_face_vertex[vh0] && feature_face_vertex[vh1])
              ehs.insert(*ce_it);
          }
        }

        if (feature_edge[*ce_it] == 0)
        {//split non fe with two fe vhs
          if (feature_edge_vertex[vh1] && feature_edge_vertex[vh0])
          {
            ehs.insert(*ce_it);
          }
        }
      }
    }


    for (const auto &ehi: ehs)
    {
      if (!mesh.is_deleted(ehi) && _excl_ehs.find(ehi) == _excl_ehs.end())
      {
        es.edge_split(ehi);
      }
    }

    if (collect_garbage && !ehs.empty() && !fhs_split.empty())
      mesh.collect_garbage();
  }

  static void
  split_singular_edges_with_singular_node(MeshT &mesh, EdgeSplitT<MeshT> &es, EP<int> &valence, const VH _vh,
                                          const bool collect_garbage = false)
  {
    if (n_incident_singular_edges(mesh, valence, _vh) > 2)
    {
      std::set<EH> ehs;
      //find halfedge whose to_vertex is singular node
      for (auto voh_it = mesh.voh_iter(_vh); voh_it.valid(); ++voh_it)
      {
        EH ve = mesh.edge_handle(*voh_it);
        if (valence[ve] != 0)
        {
          VH vht = mesh.halfedge(*voh_it).to_vertex();
          if (n_incident_singular_edges(mesh, valence, vht) > 2)
            ehs.insert(ve);
        }
      }


      for (const auto &ehi: ehs)
        if (!mesh.is_deleted(ehi))
        {
          es.edge_split(ehi);
        }

      if (collect_garbage && !ehs.empty())
        mesh.collect_garbage();
    }
  }

  static void
  split_feature_face_singular_edges_with_singular_nodes(MeshT &mesh, EdgeSplitT<MeshT> &es, EP<bool> &feature_face_edge,
                                                        EP<int> &valence, const bool collect_garbage = false)
  {
    std::set<EH> ehs;

    for (const auto ehi: mesh.edges())
    {
      if (valence[ehi] != 0 && feature_face_edge[ehi])
      {
        VH vht = mesh.edge(ehi).to_vertex();
        VH vhf = mesh.edge(ehi).from_vertex();

        if (n_incident_singular_edges(mesh, valence, vht) > 2 &&
            n_incident_singular_edges(mesh, valence, vhf) > 2)
          ehs.insert(ehi);
      }
    }

    for (const auto &ehi: ehs)
      if (!mesh.is_deleted(ehi))
      {
        es.edge_split(ehi);
      }

    if (collect_garbage && !ehs.empty())
      mesh.collect_garbage();
  }

  static void split_onering_of_zipper_nodes_arc(MeshT &mesh, EdgeSplitT<MeshT> &es, EP<int> &valence, VP<int> &sgl_vt,
                                                std::set<VH> &_tps, const bool collect_garbage = false)
  {
    std::set<EH> ehs;

    while (!_tps.empty())
    {
      auto vh_cur = *_tps.begin();

      HEH heh(-1);
      for (auto voh_it = mesh.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
        if (valence[mesh.edge_handle(*voh_it)] != 0)
        {
          heh = *voh_it;
          break;
        }

      auto sg_vhs = get_vertices_on_singular_arc(mesh, valence, heh);

      //find edges to split
      std::set<CH> or_cells;
      for (auto i = 0u; i < sg_vhs.size(); ++i)
      {
        for (auto vc_it = mesh.vc_iter(sg_vhs[i]); vc_it.valid(); ++vc_it)
          or_cells.insert(*vc_it);
      }

      for (const auto chi: or_cells)
      {
        for (auto ce_it = mesh.ce_iter(chi); ce_it.valid(); ++ce_it)
        {
          VH vh0 = mesh.edge(*ce_it).from_vertex();
          VH vh1 = mesh.edge(*ce_it).to_vertex();
          if (valence[*ce_it] == 0)
          {
            if (sgl_vt[vh1] && sgl_vt[vh0])
            {
              ehs.insert(*ce_it);
            }
          }
        }
      }

      for (const auto vhi: sg_vhs)
        _tps.erase(vhi);
    }

    for (const auto &ehi: ehs)
      if (!mesh.is_deleted(ehi))
      {
        es.edge_split(ehi);
      }

    if (collect_garbage && !ehs.empty())
      mesh.collect_garbage();
  }


  static int split_faces_with_cell_property_update(MeshT &mesh, FaceSplitT<MeshT> &_fs, std::set<CH> &_dp_cells,
                                                   const std::set<FH> &_fhs, std::set<CH> &_start_cells,
                                                   std::set<CH> &_target_cells, std::set<FH> &_excl_fhs)
  {
    int n = 0;
    for (const auto &fh: _fhs)
    {
      if (_fs.is_split_ok(fh))
      {
        auto hfh0 = mesh.halfface_handle(fh, 0);
        auto hfh1 = mesh.halfface_handle(fh, 1);

        auto ch0 = mesh.incident_cell(hfh0);
        auto ch1 = mesh.incident_cell(hfh1);

        //dp cells
        bool isdp0 = false, isdp1 = false;
        bool istg0 = false, istg1 = false;
        bool isst0 = false, isst1 = false;
        bool isex = false;
        if (ch0.is_valid() && _dp_cells.find(ch0) != _dp_cells.end())
        {
          isdp0 = true;
          _dp_cells.erase(ch0);

          //target cells
          if (_target_cells.find(ch0) != _target_cells.end())
          {
            istg0 = true;
            _target_cells.erase(ch0);
          }
          //start cells
          if (_start_cells.find(ch0) != _start_cells.end())
          {
            isst0 = true;
            _start_cells.erase(ch0);
          }
        }
        if (ch1.is_valid() && _dp_cells.find(ch1) != _dp_cells.end())
        {
          isdp1 = true;
          _dp_cells.erase(ch1);

          //target cells
          if (_target_cells.find(ch1) != _target_cells.end())
          {
            istg1 = true;
            _target_cells.erase(ch1);
          }
          //start cells
          if (_start_cells.find(ch1) != _start_cells.end())
          {
            isst1 = true;
            _start_cells.erase(ch1);
          }
        }

        if (_excl_fhs.find(fh) != _excl_fhs.end())
        {
          isex = true;
          _excl_fhs.erase(fh);
        }

        auto hfv0 = mesh.get_halfface_vertices(hfh0);
        auto hfv1 = mesh.get_halfface_vertices(hfh1);

        auto vh_new = _fs.face_split(fh);
        n++;

        if (isdp0)
        {
          for (int i = 0; i < 3; ++i)
          {
            VH vh_tmp = hfv0[i];
            hfv0[i] = vh_new;
            auto ch0 = mesh.incident_cell(mesh.find_halfface(hfv0));
            _dp_cells.insert(ch0);
            //update target cells
            if (istg0)
              _target_cells.insert(ch0);
            //update start cells
            if (isst0)
              _start_cells.insert(ch0);

            hfv0[i] = vh_tmp;
          }
        }

        if (isdp1)
        {
          for (int i = 0; i < 3; ++i)
          {
            VH vh_tmp = hfv1[i];
            hfv1[i] = vh_new;
            auto ch1 = mesh.incident_cell(mesh.find_halfface(hfv1));

            _dp_cells.insert(ch1);
            //update target cells
            if (istg1)
              _target_cells.insert(ch1);
            //update start cells
            if (isst1)
              _start_cells.insert(ch1);

            hfv1[i] = vh_tmp;
          }
        }

        if (isex)
        {
          for (int i = 0; i < 3; ++i)
          {
            VH vh_tmp = hfv1[i];
            hfv1[i] = vh_new;
            auto fh_i = mesh.face_handle(mesh.find_halfface(hfv1));
            _excl_fhs.insert(fh_i);

            hfv1[i] = vh_tmp;
          }
        }
      }
    }

    return n;
  }


  static int split_edges_with_cell_property_update(MeshT &mesh, EdgeSplitT<MeshT> &es, std::set<CH> &_dp_cells,
                                                   const std::set<EH> &_ehs, std::set<CH> &_start_cells,
                                                   std::set<CH> &_target_cells, std::set<EH> &_excl_ehs)
  {
    int n = 0;
    for (const auto eh: _ehs)
    {
      if (es.is_split_ok(eh))
      {
        std::vector<std::vector<VH>> dp_cvhs, tg_cvhs, st_cvhs;
        auto heh = mesh.halfedge_handle(eh, 0);
        auto vh_f = mesh.halfedge(heh).from_vertex();
        auto vh_t = mesh.halfedge(heh).to_vertex();

        for (auto hehf_it = mesh.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
        {
          auto ch = mesh.incident_cell(*hehf_it);

          if (ch.is_valid() && _dp_cells.find(ch) != _dp_cells.end())
          {
            auto cvhs = mesh.get_cell_vertices(*hehf_it, heh);
            dp_cvhs.push_back(cvhs);
            _dp_cells.erase(ch);

            //target cells
            if (_target_cells.find(ch) != _target_cells.end())
            {
              tg_cvhs.push_back(cvhs);
              _target_cells.erase(ch);
            }
            //start cells
            if (_start_cells.find(ch) != _start_cells.end())
            {
              st_cvhs.push_back(cvhs);
              _start_cells.erase(ch);
            }
          }
        }


        auto vh_new = es.edge_split(eh);
        n++;

        for (auto &cvhs: dp_cvhs)
        {
          auto chfh = mesh.find_halfface(std::vector<VH>{vh_new, cvhs[1], cvhs[2]});
          _dp_cells.insert(mesh.incident_cell(chfh));

          chfh = mesh.find_halfface(std::vector<VH>{cvhs[0], vh_new, cvhs[2]});
          _dp_cells.insert(mesh.incident_cell(chfh));
        }

        //update target cells
        for (auto &cvhs: tg_cvhs)
        {
          auto chfh = mesh.find_halfface(std::vector<VH>{vh_new, cvhs[1], cvhs[2]});
          auto ch0 = mesh.incident_cell(chfh);
          _target_cells.insert(ch0);

          chfh = mesh.find_halfface(std::vector<VH>{cvhs[0], vh_new, cvhs[2]});
          auto ch1 = mesh.incident_cell(chfh);
          _target_cells.insert(ch1);
        }

        //update start cells
        for (auto &cvhs: st_cvhs)
        {
          auto chfh = mesh.find_halfface(std::vector<VH>{vh_new, cvhs[1], cvhs[2]});
          auto ch0 = mesh.incident_cell(chfh);
          _start_cells.insert(ch0);

          chfh = mesh.find_halfface(std::vector<VH>{cvhs[0], vh_new, cvhs[2]});
          auto ch1 = mesh.incident_cell(chfh);
          _start_cells.insert(ch1);
        }

        //update edge set
        if (_excl_ehs.find(eh) != _excl_ehs.end())
        {
          _excl_ehs.erase(eh);
          _excl_ehs.insert(mesh.edge_handle(mesh.find_halfedge(vh_f, vh_new)));
          _excl_ehs.insert(mesh.edge_handle(mesh.find_halfedge(vh_t, vh_new)));
        }
      }
    }

    return n;
  }


  static int split_edge_cells_at_high_valence(MeshT &mesh, CellSplitT<MeshT> &cs, const HEH heh, FP<int> &feature_fprop,
                                              EP<int> &valence, EP<bool> &feature_face_edge)
  {
    int n_split = 0;

    bool do_split = false;
    EH ehi = mesh.edge_handle(heh);
    if (feature_face_edge[ehi] && valence[ehi] > 0)
    {
      //find feature faces
      std::vector<HFH> ft_hfhs;
      for (auto hehf_it = mesh.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
      {
        if (feature_fprop[mesh.face_handle(*hehf_it)] > 0)
          ft_hfhs.push_back(*hehf_it);
      }
      int n_incf = 0;
      HFH hf_itt = ft_hfhs[0];
      do
      {
        //next incident halfface
        hf_itt = mesh.adjacent_halfface_in_cell(hf_itt, heh);

        if (!hf_itt.is_valid())
        {
          std::cerr << "Error: adjacent halfface is invalid!" << std::endl;
          break;
        }

        hf_itt = mesh.opposite_halfface_handle(hf_itt);

        if (feature_fprop[mesh.face_handle(hf_itt)] > 0)
          break;

        n_incf++;
      }
      while (hf_itt != ft_hfhs[1]);

      if (valence[ehi] + 1 >= 2 * n_incf)
        do_split = true;


      //other side
      n_incf = 0;
      hf_itt = ft_hfhs[1];
      do
      {
        //next incident halfface
        hf_itt = mesh.adjacent_halfface_in_cell(hf_itt, heh);

        if (!hf_itt.is_valid())
        {
          std::cerr << "Error: adjacent halfface is invalid!" << std::endl;
          break;
        }

        hf_itt = mesh.opposite_halfface_handle(hf_itt);

        if (feature_fprop[mesh.face_handle(hf_itt)] > 0)
          break;

        n_incf++;
      }
      while (hf_itt != ft_hfhs[0]);

      if (valence[ehi] + 1 >= 2 * n_incf)
        do_split = true;
    }
    else
    {
      int n_ifs = 0;
      for (auto hehf_it = mesh.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
      {
        n_ifs++;
      }

      if (valence[ehi] + 4 >= n_ifs)
        do_split = true;
    }

    if (do_split)
    {
      std::set<CH> chs_split;
      for (auto ec_it = mesh.ec_iter(ehi); ec_it.valid(); ++ec_it)
      {
        chs_split.insert(*ec_it);
      }

      //split
      for (const auto &chi: chs_split)
        cs.cell_split(chi);

      n_split = chs_split.size();
    }

    return n_split;
  }


  static void split_for_dof(MeshT &mesh)
  {
    std::cerr << "###Splitting w.r.t. feature..." << std::endl;
    auto feature_fprop = mesh.template request_face_property<int>("AlgoHex::FeatureFaces");
    auto feature_edge_prop = mesh.template request_edge_property<int>("AlgoHex::FeatureEdges");
    auto feature_face_vertex = mesh.template request_vertex_property<bool>("AlgoHex::FeatureFaceVertices");
    auto feature_face_edge = mesh.template request_edge_property<bool>("AlgoHex::FeatureFaceEdges");
    auto feature_edge_vertex = mesh.template request_vertex_property<bool>("AlgoHex::FeatureEdgeVertices");

    FaceSplitT<MeshT> fs(mesh);
    std::set<FH> fhs_split;
    for (const auto fh: mesh.faces())
    {
      if (feature_fprop[fh] == 0)
      {//split non-feature face which has three feature face vertices
        int count = 0;
        for (auto fv_it = mesh.fv_iter(fh); fv_it.valid(); ++fv_it)
        {
          if (feature_face_vertex[*fv_it])
            count++;
          else
            break;
        }

        if (count == 3)
          fhs_split.insert(fh);
      }
      else
      {//split feature face which has >= 2 feature edges
        int count = 0;
        for (auto fe_it = mesh.fe_iter(fh); fe_it.valid(); ++fe_it)
        {
          if (feature_edge_prop[*fe_it] > 0)
            count++;
        }

        if (count >= 2)
          fhs_split.insert(fh);
      }
    }

    for (const auto fhi: fhs_split)
      if (fs.is_split_ok(fhi))
        fs.face_split(fhi);

//        mesh.collect_garbage();

    //split non-feature face edge that has two feature face vertices
    EdgeSplitT<MeshT> es(mesh);
    std::set<EH> ehs_split;
    for (const auto eh: mesh.edges())
    {
      if (!feature_face_edge[eh])
      {
        VH vh_f = mesh.edge(eh).from_vertex();
        VH vh_t = mesh.edge(eh).to_vertex();
        if (feature_face_vertex[vh_f] && feature_face_vertex[vh_t])
          ehs_split.insert(eh);
      }
      else
      {
        if (feature_edge_prop[eh] == 0)
        {
          VH vh0 = mesh.edge(eh).from_vertex();
          VH vh1 = mesh.edge(eh).to_vertex();
          if (feature_edge_vertex[vh1] && feature_edge_vertex[vh0])
          {
            ehs_split.insert(eh);
          }
        }
      }
    }

    for (const auto ehi: ehs_split)
      if (es.is_split_ok(ehi))
      {
        es.edge_split(ehi, false);
      }

    mesh.collect_garbage();
  }


  static void sort_singular_arc_label(const MeshT &mesh, const EP<int> &valence, EP<int> &label)
  {
    for (auto e_it = mesh.e_iter(); e_it.valid(); ++e_it)
      label[*e_it] = 0;

    int i_label = 1;
    std::vector<bool> edge_conquered(mesh.n_edges(), false);
    for (auto e_it = mesh.e_iter(); e_it.valid(); ++e_it)
    {
      if (valence[*e_it] != 0 && !edge_conquered[(*e_it).idx()])
      {
        auto he0 = mesh.halfedge_handle(*e_it, 0);
        edge_conquered[(*e_it).idx()] = true;
        label[*e_it] = i_label;

        auto vh_s = mesh.halfedge(he0).from_vertex();
        auto vh_t = mesh.halfedge(he0).to_vertex();

        std::stack<VH> st_vhs;
        if (!is_singular_node(mesh, valence, vh_t))
          st_vhs.push(vh_t);
        if (!is_singular_node(mesh, valence, vh_s))
          st_vhs.push(vh_s);
        //trace back until the boundary or an edge which is conquered, then trace front
        while (!st_vhs.empty())
        {
          VH vh_cur = st_vhs.top();
          st_vhs.pop();

          for (auto voh_it = mesh.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
          {
            auto eh_og = mesh.edge_handle(*voh_it);

            if (valence[eh_og] != 0 && !edge_conquered[eh_og.idx()])
            {
              edge_conquered[eh_og.idx()] = true;
              label[eh_og] = i_label;

              auto vh_next = mesh.halfedge(*voh_it).to_vertex();
              if (!is_singular_node(mesh, valence, vh_next))
                st_vhs.push(vh_next);
            }
          }
        }

        i_label++;
      }
    }
  }


  //split if the surface mesh of the onering has >= two singular edges
  static void
  split_edges_for_local_meshability_test(MeshT &mesh, EdgeSplitT<MeshT> &es, EP<int> &valence, EP<int> &feature_edge)
  {
    auto label = mesh.template request_edge_property<int>("singular arc label", 0);
    mesh.set_persistent(label, true);
    sort_singular_arc_label(mesh, valence, label);

    std::queue<VH> vh_que;
//        for(const auto vhi : mesh.vertices())
//            vh_que.push(vhi);

    int n_split = 0;

    auto sgl_vt = mesh.template request_vertex_property<int>("singular_vertex");
    //store singular nodes
    std::set<VH> sg_nodes;
    for (const auto vhi: mesh.vertices())
    {
      if (sgl_vt[vhi] > 0)
      {
        if (n_incident_singular_edges(mesh, valence, vhi) > 2)
          sg_nodes.insert(vhi);
      }
    }

    //split any singular edge with two sg nodes
    std::set<EH> ehs_split;
    for (const auto ehi: mesh.edges())
    {
      if (valence[ehi] != 0)
      {
        auto vhf = mesh.edge(ehi).from_vertex();
        auto vht = mesh.edge(ehi).to_vertex();

        if (sg_nodes.find(vhf) != sg_nodes.end() &&
            sg_nodes.find(vht) != sg_nodes.end())
        {
          ehs_split.insert(ehi);
        }
      }
    }
    for (const auto &ehi: ehs_split)
    {
      if (!mesh.is_deleted(ehi))
      {
        es.edge_split(ehi);
      }
    }

    //split if onering boundary has more than one sg arc labels
    for (const auto vhi: mesh.vertices())
      vh_que.push(vhi);
    while (!vh_que.empty() && n_split < 1000000)
    {
      VH vh_cur = vh_que.front();
      vh_que.pop();

      std::set<EH> all_ehs, inc_ehs;
      for (auto ve_it = mesh.ve_iter(vh_cur); ve_it.valid(); ++ve_it)
        inc_ehs.insert(*ve_it);

      std::set<EH> sg_ehs;
      std::set<int> labels;
      for (auto vc_it = mesh.vc_iter(vh_cur); vc_it.valid(); ++vc_it)
        for (auto ce_it = mesh.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
          if (valence[*ce_it] != 0 && inc_ehs.find(*ce_it) == inc_ehs.end())
          {
            sg_ehs.insert(*ce_it);
            labels.insert(label[*ce_it]);
          }

      if (labels.size() >= 2)
      {
        std::set<EH> ehs_split;
        std::map<VH, int> vh_int;
        for (auto ehi: sg_ehs)
        {
          auto vhf = mesh.edge(ehi).from_vertex();
          auto vht = mesh.edge(ehi).to_vertex();

          if (vh_int.find(vhf) == vh_int.end())
            vh_int[vhf] = 1;
          else
            vh_int[vhf] += 1;

          if (vh_int.find(vht) == vh_int.end())
            vh_int[vht] = 1;
          else
            vh_int[vht] += 1;
        }


        if (ehs_split.empty())
        {
          //split if incident to more than one sge
          for (auto&[vhi, num]: vh_int)
          {
            if (num > 1 && sg_nodes.find(vhi) == sg_nodes.end())
            {
              auto eh0 = mesh.edge_handle(mesh.find_halfedge(vh_cur, vhi));
              if (eh0.is_valid() && valence[eh0] == 0)
                ehs_split.insert(eh0);
            }
          }
        }

        if (!ehs_split.empty())
          vh_que.push(vh_cur);

        for (auto ehi: ehs_split)
        {
          if (es.is_split_ok(ehi))
          {
            vh_que.push(es.edge_split(ehi));
            n_split++;
          }
        }

        ALGOHEX_DEBUG_ONLY(if (n_split % 10000 == 0) std::cerr << "performed " << n_split << " split - a" << std::endl;)
      }
    }

    std::cerr << "split " << n_split << " edges for local meshability test - a " << std::endl;

    n_split = 0;
    //split
    for (const auto vhi: mesh.vertices())
      vh_que.push(vhi);

    while (!vh_que.empty() && n_split < 1000000)
    {
      VH vh_cur = vh_que.front();
      vh_que.pop();

      std::set<EH> all_ehs, inc_ehs;
      for (auto ve_it = mesh.ve_iter(vh_cur); ve_it.valid(); ++ve_it)
        inc_ehs.insert(*ve_it);

      std::set<EH> sg_ehs;
      std::set<int> labels;
      for (auto vc_it = mesh.vc_iter(vh_cur); vc_it.valid(); ++vc_it)
        for (auto ce_it = mesh.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
          if (valence[*ce_it] != 0 && inc_ehs.find(*ce_it) == inc_ehs.end())
          {
            sg_ehs.insert(*ce_it);
            labels.insert(label[*ce_it]);
          }

      //it on the one ring boundary, there are two arcs
      if (labels.size() >= 2)
      {
        std::set<EH> ehs_split;
        std::map<VH, int> vh_int;
        for (auto ehi: sg_ehs)
        {
          auto vhf = mesh.edge(ehi).from_vertex();
          auto vht = mesh.edge(ehi).to_vertex();

          if (vh_int.find(vhf) == vh_int.end())
            vh_int[vhf] = 1;
          else
            vh_int[vhf] += 1;

          if (vh_int.find(vht) == vh_int.end())
            vh_int[vht] = 1;
          else
            vh_int[vht] += 1;
        }

        //split if open end
        for (auto&[vhi, num]: vh_int)
        {
          if (num == 1 && sg_nodes.find(vhi) == sg_nodes.end())
          {
            auto eh0 = mesh.edge_handle(mesh.find_halfedge(vh_cur, vhi));
            if (eh0.is_valid() && valence[eh0] == 0)
              ehs_split.insert(eh0);
          }
        }

        if (!ehs_split.empty())
          vh_que.push(vh_cur);

        for (auto ehi: ehs_split)
        {
          if (es.is_split_ok(ehi))
          {
            vh_que.push(es.edge_split(ehi));
            n_split++;
          }
        }

        ALGOHEX_DEBUG_ONLY(if (n_split % 10000 == 0) std::cerr << "performed " << n_split << " split - b" << std::endl;)
      }
    }

    std::cerr << "split " << n_split << " edges for local meshability test - b" << std::endl;


    mesh.collect_garbage();

    std::cerr << "\n#Tet mesh data after split:";
    std::cerr << "\n#Tet mesh vertex number: " << mesh.n_vertices();
    std::cerr << "\n#Tet mesh edge number: " << mesh.n_edges();
    std::cerr << "\n#Tet mesh face number: " << mesh.n_faces();
    std::cerr << "\n#Tet mesh cell number: " << mesh.n_cells() << std::endl;
  }

  int n_incident_special_edges(const VH &_vh) const
  {
    int n_spc = 0;
    for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
      if (valence_[*ve_it] != 0 || feature_edge_[*ve_it] > 0)
        n_spc++;

    return n_spc;
  }

private:
  MeshT &mesh_;
  EdgeSplitT<MeshT> es_;
  FaceSplitT<MeshT> fs_;
  VertexOptimizeT<MeshT> vo_;
};
}

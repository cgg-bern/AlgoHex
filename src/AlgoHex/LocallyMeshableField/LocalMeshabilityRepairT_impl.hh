/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define FIXLOCALLYNONMESHABLEVERTICEST_C

#include "LocalMeshabilityRepairT.hh"

namespace AlgoHex
{

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::preprocess(bool _check_alignment)
{
  std::cerr << "####Preprocessing..." << std::endl;

  std::cerr << "####Fix misalignment at feature edges ..." << std::endl;

  fix_misalignment_at_feature_edges();

  std::cerr << "####Fix misalignment at feature face sectors ..." << std::endl;

  fix_misalignment_at_feature_face_sectors();
  if (_check_alignment)
    check_alignments();

  std::cerr << "####Pushing boundary singular edges which are not on feature arcs to interior ..." << std::endl;
  push_singular_vertices(false, true);
  if (_check_alignment)
    check_alignments();

  std::cerr << "####Removing singular triangles ..." << std::endl;
  remove_singular_triangles();
  if (_check_alignment)
    check_alignments();
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::fix_first_stage(bool _check_alignment)
{
  int n_fixed = 0;
  if (iter_ < first_iters_)
  {
    std::cerr << "####Fixing local invalid at first stage..." << std::endl;
    if (_check_alignment)
      check_alignments();

    //Ensure edge meshability: fixing non-meshable footprint
    if (iter_ == 0)
    {
      std::cerr << "####Fix non-meshable footprint on ff ..." << std::endl;
      fix_edges_with_non_meshable_footprint();
    }

    //Ensure edge meshability: fixing compound singular edges
    std::cerr << "####Fix complex singular edges non ff ..." << std::endl;
    fix_compound_singular_edges(true);
    if (_check_alignment)
      check_alignments();

    //Fix invalid singular nodes
    if (iter_ >= 4)
    { //let singular circle shrink before fixing
      if (merge_zipper_node_)
        n_fixed += merge_zipper_nodes();

      std::cerr << "########Fix interior(non ffv) invalid nodes... " << std::endl;
      n_fixed += fix_interior_invalid_singular_nodes();
    }

    std::cerr << "########Fix invalid nodes on boundary... " << std::endl;
    n_fixed += fix_boundary_invalid_singular_nodes();
    if (_check_alignment)
      check_alignments();

    std::cerr << "########Fix invalid singular vertices on feature surface... " << std::endl;
    n_fixed += fix_invalid_singular_vertices_on_feature_surface();
    if (_check_alignment)
      check_alignments();

    std::cerr << "########Fix fully constrained parabolic sectors... " << std::endl;
    n_fixed += fix_fully_constrained_parabolic_sectors();
    if (_check_alignment)
      check_alignments();

    //Fix constrained zipper nodes
    std::cerr << "########Fix constrained zero sectors1... " << std::endl;
    n_fixed += fix_constrained_zipper_nodes();
    if (_check_alignment)
      check_alignments();

    {
      std::cerr << "########Push sgv to interior... " << std::endl;
      n_fixed += push_singular_vertices(true, false);

      //last priority: start from a boundary singular edge if pushing to interior is not possible, e.g. i28_tire
      n_fixed += fix_boundary_invalid_singular_nodes_start_from_bdy_sge();

      if (_check_alignment)
        check_alignments();
    }

    //remove zipper nodes by zipping in onering.
    if (iter_ >= 4)
    {
      std::cerr << "########Remove zipper nodes... " << std::endl;
      n_fixed += zip_zipper_nodes();
      if (_check_alignment)
        check_alignments();
    }
  }

  return n_fixed;
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::fix_second_stage(bool _check_alignment)
{
  int n_fixed = 0;
  int second_iters = iter_ - first_iters_;
  fis_.set_enable_detach_kept_singular_edge(true);
  fctp_.set_second_stage();


  if (second_iters % second_inner_iters_ == 0 && second_iters >= 0)
  {
    std::cerr << "####Fixing local invalid singular nodes at second stage..." << std::endl;
    if (_check_alignment)
      check_alignments();

    //Optional: remove zipper nodes by matching adjustment.
    std::cerr << "########Remove zipper nodes... " << std::endl;
    n_fixed += zip_zipper_nodes();
    if (_check_alignment)
      check_alignments();


    //Fix complex singular edges
    std::cerr << "####Fix complex singular edges ..." << std::endl;
    n_fixed += fix_compound_singular_edges(false);
    if (_check_alignment)
      check_alignments();

    //Fix invalid singular nodes
    std::cerr << "########Fix interior(non ffv) invalid nodes... " << std::endl;
    n_fixed += fix_interior_invalid_singular_nodes();
    if (_check_alignment)
      check_alignments();
    //
    std::cerr << "########Fix invalid nodes on boundary... " << std::endl;
    n_fixed += fix_boundary_invalid_singular_nodes();
    if (_check_alignment)
      check_alignments();

    std::cerr << "########Fix invalid nodes on feature surface... " << std::endl;
    n_fixed += fix_invalid_singular_nodes_on_feature_surface();

    //Fix zero feature sectors
    std::cerr << "########Fix fully constrained parabolic sectors... " << std::endl;
    n_fixed += fix_fully_constrained_parabolic_sectors();
    if (_check_alignment)
      check_alignments();

    //Fix constrained zipper nodes
    std::cerr << "########Fix constrained tps... " << std::endl;
    n_fixed += fix_constrained_zipper_nodes();
    if (_check_alignment)
      check_alignments();


    std::cerr << "########Fix zipper nodes... " << std::endl;
    n_fixed += unzip_zipper_nodes();
    n_fixed += fix_invalid_singular_nodes_on_feature_surface();

    if (_check_alignment)
      check_alignments();
  }

  iter_++;

  return n_fixed;
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::split_for_local_meshability_test()
{
  SplitHelperT<MeshT>::split_edges_for_local_meshability_test(mesh_, es_, valence_, feature_edge_);
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::fix_interior_invalid_singular_nodes()
{
  std::vector<VH> cd_vhs;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] != 0 && !feature_face_vertex_[vhi])
    {
      int n_inc_sge = n_incident_interior_singular_edges(vhi);
      if (n_inc_sge >= 3)
        cd_vhs.push_back(vhi);
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "fixing invalid interior singular nodes..." << std::endl;)
  int n_fixed = 0;
  for (const auto vhi: cd_vhs)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "\nfixing vh: " << vhi;)
    bool fixed = fis_.fix_interior_invalid_singular_node(vhi);
    if (fixed)
      n_fixed++;
  }

  mesh_.collect_garbage();

  return n_fixed;
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::fix_boundary_invalid_singular_nodes_start_from_bdy_sge()
{
  std::vector<VH> cd_vhs;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] != 0 && mesh_.is_boundary(vhi))
    {
      if (feature_node_[vhi])
        cd_vhs.push_back(vhi);
    }
  }

  ALGOHEX_DEBUG_ONLY(
          std::cerr << "###Fixing invalid boundary singular nodes (start from boundary sge)..." << std::endl;)

  int n_fixed = 0;
  for (const auto vhi: cd_vhs)
  {
//            std::cerr<<"\nfixing vh: "<<vhi<<std::endl;
    bool fixed = fis_.fix_boundary_invalid_singular_node(vhi, false);
    if (!fixed)
      fixed = fis_.fix_boundary_invalid_singular_node(vhi, true);

    if (fixed)
      n_fixed++;
  }

  mesh_.collect_garbage();

  ALGOHEX_DEBUG_ONLY(
          std::cerr << "###Fixed  singular nodes " << n_fixed << " (start from boundary sge)..." << std::endl;)


  return n_fixed;
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::fix_boundary_invalid_singular_nodes()
{
  std::set<std::pair<int, VH>> que;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] != 0 && mesh_.is_boundary(vhi))
    {
      int n_sge = n_incident_singular_edges(mesh_, valence_, vhi);
      que.emplace(n_sge, vhi);
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "###Fixing invalid boundary singular nodes..." << std::endl;)

  int n_fixed = 0, n_pushed = 0;
  while (!que.empty())
  {
    auto ivh_cur = *que.begin();
    que.erase(que.begin());

    int n_sge_new = n_incident_singular_edges(mesh_, valence_, ivh_cur.second);
    if (n_sge_new != ivh_cur.first)
      continue;

    bool fixed = fis_.fix_boundary_invalid_singular_node(ivh_cur.second, false);
    if (!fixed)
    {
      if (!feature_node_[ivh_cur.second] && psv_.is_move_inside_candidate(ivh_cur.second, true, false))
      {

        //move singular vertex to the interior
        fixed = psv_.push_singular_vertex(ivh_cur.second);

        n_pushed++;
      }
    }

    if (fixed)
    {
      n_fixed++;

      //push
      for (auto voh_it = mesh_.voh_iter(ivh_cur.second); voh_it.valid(); ++voh_it)
      {
        auto vht = mesh_.halfedge(*voh_it).to_vertex();
        if (sgl_vt_[vht] != 0 && feature_face_vertex_[vht])
        {
          int n_sge = n_incident_singular_edges(mesh_, valence_, vht);
          que.emplace(n_sge, vht);
        }
      }
    }
  }

  mesh_.collect_garbage();

  ALGOHEX_DEBUG_ONLY(
          std::cerr << "###Fixed " << n_fixed << " boundary singular vertices, pushed " << n_pushed << std::endl;)


  return n_fixed;
}


template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::
fix_invalid_singular_vertices_on_feature_surface()
{
  std::set<std::pair<int, VH>> que;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] != 0 && feature_face_vertex_[vhi] && !mesh_.is_boundary(vhi))
    {
      int n_sge = n_incident_singular_edges(mesh_, valence_, vhi);
      que.emplace(n_sge, vhi);
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "###Fixing invalid singular vertices on feature surface..." << std::endl;)
  int n_fixed = 0;
  while (!que.empty())
  {
    auto ivh_cur = *que.begin();
    que.erase(que.begin());

    int n_sge_new = n_incident_singular_edges(mesh_, valence_, ivh_cur.second);
    if (n_sge_new != ivh_cur.first)
      continue;

    bool fixed = false;
    //try fix vts on sg arc
    if (fctp_.singular_arc_tangentially_touching_feature_face_type(ivh_cur.second) == 1)
    {
      fixed = fix_vertex_on_tangential_singular_arc_at_feature_face(ivh_cur.second);
    }

    //if not fixed on sg arc, fix sg nodes
    if (!fixed)
      fixed = fis_.fix_invalid_singular_node_on_feature_surface(ivh_cur.second);


    if (fixed)
    {
      n_fixed++;

      //push
      for (auto voh_it = mesh_.voh_iter(ivh_cur.second); voh_it.valid(); ++voh_it)
      {
        auto vht = mesh_.halfedge(*voh_it).to_vertex();
        if (sgl_vt_[vht] != 0 && feature_face_vertex_[vht])
        {
          int n_sge = n_incident_singular_edges(mesh_, valence_, vht);
          que.emplace(n_sge, vht);
        }
      }
    }
  }


  ALGOHEX_DEBUG_ONLY(
          std::cerr << "\nFixed " << n_fixed << " singular vertices on feature surface(not bdy)!" << std::endl;)


  return n_fixed;
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::fix_invalid_singular_nodes_on_feature_surface(const bool _include_bdy)
{
  std::vector<VH> cd_vhs;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] != 0 && feature_face_vertex_[vhi] && (_include_bdy || !mesh_.is_boundary(vhi)))
      cd_vhs.push_back(vhi);
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "###Fixing invalid singular nodes on feature surface..." << std::endl;)
  int n_fixed = 0;
  for (const auto vhi: cd_vhs)
  {
    bool fixed = fis_.fix_invalid_singular_node_on_feature_surface(vhi);

    if (fixed)
      n_fixed++;
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "\nFixed " << n_fixed << " singular nodes on feature surface(not bdy)!" << std::endl;)

  return n_fixed;
}

template<class MeshT>
bool
LocalMeshabilityRepairT<MeshT>::fix_vertex_on_tangential_singular_arc_at_feature_face(const VH _vh)
{
  auto sg_hehs = get_singular_halfedges(mesh_, valence_, _vh);

  if (sg_hehs.size() != 2)
    return false;

  std::vector<EH> sg_ehs{mesh_.edge_handle(sg_hehs[0]), mesh_.edge_handle(sg_hehs[1])};
  for (const auto ehi: sg_ehs)
  {
    int old_val = valence_[ehi];
    sge_.compute_edge_valence(ehi);
    if (valence_[ehi] != old_val)
      std::cerr << "\nWarning: edge valence changed!" << std::endl;
  }

  if (!feature_face_edge_[sg_ehs[0]] && !feature_face_edge_[sg_ehs[1]])
  {
    if (fis_.detach_singular_arc_from_singular_node(_vh, false, true))
      return true;

    if (fis_.detach_singular_arcs_of_zipper_node(_vh, false, true))
      return true;
  }
  else if ((feature_face_edge_[sg_ehs[0]] && !feature_face_edge_[sg_ehs[1]]) ||
           (!feature_face_edge_[sg_ehs[0]] && feature_face_edge_[sg_ehs[1]]))
  {
    HEH ff_heh = feature_face_edge_[sg_ehs[0]] ? sg_hehs[0] : sg_hehs[1];
    HEH nff_heh = feature_face_edge_[sg_ehs[0]] ? sg_hehs[1] : sg_hehs[0];

    if (fis_.is_valid_face_sector_after_fix(_vh, ff_heh, nff_heh))
    {
      if (fis_.detach_singular_arc_from_singular_node(_vh, true, true))
        return true;

      return fis_.detach_singular_arcs_of_zipper_node(_vh, true, true);
    }
    else
    {//change matching in the sector on the other side (make sure the singular arc is not on feature surface)
      auto or_chs = get_onering_cells(_vh);
      std::set<HFH> bound_ft_hfhs;
      auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, *mesh_.hec_iter(nff_heh),
                                               bound_ft_hfhs);

      HFH hf_s, hf_e;
      for (auto hehf_it = mesh_.hehf_iter(ff_heh); hehf_it.valid(); ++hehf_it)
      {
        if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
        {
          if (bound_ft_hfhs.find(*hehf_it) != bound_ft_hfhs.end())
          {
            hf_e = *hehf_it;
          }
          else if (bound_ft_hfhs.find(mesh_.opposite_halfface_handle(*hehf_it)) != bound_ft_hfhs.end())
          {
            hf_s = *hehf_it;
          }
        }
      }

      if (!hf_s.is_valid() || !hf_e.is_valid())
      {
        std::cerr << "Error: couldn't find start halfface and end halfface at feature face sector!" << std::endl;
        return false;
      }

      HFH hf_change = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_s, ff_heh));
      int old_trans = trans_prop_[hf_change];

      int e_trans = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, ff_heh, hf_change);

      int new_trans = tq_.mult_transitions_idx(old_trans, tq_.inverse_transition_idx(e_trans));

      return eem_.attempt_matching_adjustment_of_feature_face_sector(ff_heh, hf_s, hf_e, new_trans);
    }
  }

  return false;
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::
push_singular_vertices(bool _push_feature_vertex, bool _push_non_feature_boundary_singular_vertex)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "###Pushing boundary singular vertices to interior..." << std::endl;)

  int n = 0;

  std::queue<VH> que;
  for (const auto &vh: mesh_.vertices())
  {
    if (sgl_vt_[vh] != 0 && mesh_.is_boundary(vh))
    {
      que.push(vh);
    }
  }

  while (!que.empty())
  {
    auto vh_cur = que.front();
    que.pop();

    if (!psv_.is_move_inside_candidate(vh_cur, _push_feature_vertex, _push_non_feature_boundary_singular_vertex))
      continue;

    //
    if (n_incident_feature_singular_edges(mesh_, feature_edge_, valence_, vh_cur) >= 2)
    {
      int n_int_sge = n_incident_interior_singular_edges(vh_cur);
      int n_all_sge = n_incident_singular_edges(mesh_, valence_, vh_cur);
      if (n_int_sge >= 2)
      {
        if (fis_.detach_singular_arc_from_singular_node(vh_cur, false, true))
          continue;
        else if (fis_.detach_singular_arcs_of_zipper_node(vh_cur, false, true))
          continue;
      }
      else if (n_int_sge == 1 && n_all_sge == 2)
      {
        if (fis_.detach_singular_arc_from_singular_node(vh_cur, true, true))
          continue;
        else if (fis_.detach_singular_arcs_of_zipper_node(vh_cur, true, true))
          continue;
      }

      if (fis_.detach_singular_arc_from_singular_node_wrt_boundary(vh_cur, true, false))
        continue;
    }

    //move singular vertex to the interior
    if (psv_.push_singular_vertex(vh_cur))
    {
      n++;

      //push
      for (auto voh_it = mesh_.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
      {
        auto veh = mesh_.edge_handle(*voh_it);
        if (valence_[veh] != 0)
        {
          que.push(mesh_.halfedge(*voh_it).to_vertex());
        }
      }
    }
  }

  if (push_singular_circle_)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "\nPushing boundary sg circle vertices!" << std::endl;)

    int n_circle_vhs = 0;
    std::set<EH> bdy_sgehs;
    for (const auto ehi: mesh_.edges())
    {
      if (valence_[ehi] == 1 && mesh_.is_boundary(ehi))
      {
        bdy_sgehs.insert(ehi);
      }
    }

    std::set<VH> circle_vhs;
    while (!bdy_sgehs.empty())
    {
      auto eh_cur = *bdy_sgehs.begin();
      bdy_sgehs.erase(bdy_sgehs.begin());

      auto sg_hehs = get_halfedges_of_singular_arc_on_feature_face(mesh_, valence_, feature_face_edge_,
                                                                   mesh_.halfedge_handle(eh_cur, 0), true);

      //it's a circle
      VH vh0 = mesh_.halfedge(sg_hehs[0]).from_vertex();
      VH vh1 = mesh_.halfedge(sg_hehs.back()).to_vertex();
      if (vh0 == vh1)
      {
        int n_sge = n_incident_singular_edges(mesh_, valence_, vh0);
        if (n_sge == 2)
        {
          for (const auto &hei: sg_hehs)
          {
            if (mesh_.is_boundary(hei))
            {
              circle_vhs.insert(mesh_.halfedge(hei).from_vertex());
              circle_vhs.insert(mesh_.halfedge(hei).to_vertex());

              bdy_sgehs.erase(mesh_.edge_handle(hei));
            }
          }
        }
      }
    }

    for (const auto vhi: circle_vhs)
    {
      psv_.push_singular_vertex(vhi);
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "\nPushed " << circle_vhs.size() << " boundary sg circle vertices!" << std::endl;)
  }

  mesh_.collect_garbage();

  ALGOHEX_DEBUG_ONLY(std::cerr << "\nPushed " << n << " boundary singular vertices!" << std::endl;)

  return n;
}

template<class MeshT>
int LocalMeshabilityRepairT<MeshT>::flatten_non_singular_feature_edges()
{
  int n_spd = 0;

  std::queue<VH> v_que;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] && feature_edge_vertex_[vhi] && !mesh_.is_boundary(vhi))
      v_que.push(vhi);
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "###Flattening... " << std::endl;)

  while (!v_que.empty())
  {
    auto vh_cur = v_que.front();
    v_que.pop();

    std::set<HEH> ft_hehs;
    if (fctp_.is_separable_at_vertex(vh_cur, ft_hehs))
    {
      //separate
      for (const auto hehi: ft_hehs)
      {
        //the feature edge will be split, store the to_vertex in ahead
        VH vht = mesh_.halfedge(hehi).to_vertex();
//        std::cerr << " \nfe vhs " << vh_cur << " " << vht << std::endl;
        bool suc = flatten_non_singular_feature_edge(hehi);
        if (suc)
        {
          v_que.push(vht);
          n_spd++;
        }
        else
        {
          std::cerr << " failed to separate.";
        }
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "###Flattened " << n_spd << " edges " << std::endl;)

  fix_invalid_singular_nodes_on_feature_surface();

  return n_spd;
}


template<class MeshT>
bool LocalMeshabilityRepairT<MeshT>::flatten_non_singular_feature_edge(const HEH _heh)
{
  std::vector<HFH> ft_hfhs;
  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
      ft_hfhs.push_back(*hehf_it);
  }
  if (ft_hfhs.size() != 2)
  {
    std::cerr << "Error: non-manifold feature surface at edge!" << std::endl;
    return false;
  }

  int fsec_agl = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                     trans_prop_, _heh, ft_hfhs[0],
                                                                                     ft_hfhs[1]);

  VH vhf = mesh_.halfedge(_heh).from_vertex();
  VH vht = mesh_.halfedge(_heh).to_vertex();

  {
    //singular free start face
    HFH hf_change = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(ft_hfhs[0], _heh));
    HEH he_open = mesh_.prev_halfedge_in_halfface(_heh, hf_change);

    VH vh_opp = mesh_.halfedge(he_open).from_vertex();

//            if (sgl_vt_[vh_opp] != 0)
//split to better prevent valence flipping
    {
      VH vhm = fs_.face_split(mesh_.face_handle(hf_change));

      std::vector<VH> hfvhs;
      hfvhs.push_back(mesh_.halfedge(_heh).from_vertex());
      hfvhs.push_back(mesh_.halfedge(_heh).to_vertex());
      hfvhs.push_back(vhm);

      hf_change = mesh_.find_halfface(hfvhs);
      he_open = mesh_.prev_halfedge_in_halfface(_heh, hf_change);
    }
    //get rotation axis
    CH ch0 = mesh_.incident_cell(ft_hfhs[0]);
    Point dir = mesh_.vertex(mesh_.halfedge(_heh).to_vertex()) -
                mesh_.vertex(mesh_.halfedge(_heh).from_vertex());
    int rt_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch0], dir).second;
    if (fsec_agl == 1 || fsec_agl == 0)
      rt_axis = negate(rt_axis);

    int rt_trans = rt_axis + 4;
    int new_trans = tq_.mult_transitions_idx(trans_prop_[hf_change],
                                             rt_trans);

    bool fixed = eem_.attempt_matching_adjustment_of_feature_face_sector(_heh, ft_hfhs[0], ft_hfhs[1], new_trans);

    //the other side
    int fsec_agl1 = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                        trans_prop_, _heh, ft_hfhs[1],
                                                                                        ft_hfhs[0]);

    HFH hf_change1 = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(ft_hfhs[1], _heh));
    HEH he_open1 = mesh_.prev_halfedge_in_halfface(_heh, hf_change1);

    VH vhf_opp1 = mesh_.halfedge(he_open1).from_vertex();
//            if (sgl_vt_[vhf_opp1] != 0)
    {
      VH vhm = fs_.face_split(mesh_.face_handle(hf_change1));

      std::vector<VH> hfvhs;
      hfvhs.push_back(mesh_.halfedge(_heh).from_vertex());
      hfvhs.push_back(mesh_.halfedge(_heh).to_vertex());
      hfvhs.push_back(vhm);

      hf_change1 = mesh_.find_halfface(hfvhs);
      he_open1 = mesh_.prev_halfedge_in_halfface(_heh, hf_change1);
    }

    CH ch1 = mesh_.incident_cell(ft_hfhs[1]);
    int rt_axis1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch1], dir).second;
    if (fsec_agl1 == 1 || fsec_agl1 == 0)
      rt_axis1 = negate(rt_axis1);

    int rt_trans1 = rt_axis1 + 4;
    int new_trans1 = tq_.mult_transitions_idx(trans_prop_[hf_change1],
                                              rt_trans1);

    bool fixed1 = eem_.attempt_matching_adjustment_of_feature_face_sector(_heh, ft_hfhs[1], ft_hfhs[0], new_trans1);


    if (fixed)
    {
      if (sgl_vt_[vhf])
      {
        if (std::abs(valence_[mesh_.edge_handle(he_open)]) == 1)
        {
          bool suc = fis_.detach_singular_arc_from_singular_node(vhf,
                                                                 mesh_.opposite_halfedge_handle(he_open),
                                                                 std::vector<HEH>{}, false, true);
          if (!suc)
            fis_.detach_singular_arcs_of_zipper_node(vhf, mesh_.opposite_halfedge_handle(he_open),
                                                     std::vector<HEH>{}, false, true);
        }
        else
        {
          bool suc = fis_.detach_singular_arc_from_singular_node(vhf, false, true);
          if (!suc)
            fis_.detach_singular_arcs_of_zipper_node(vhf, false, true);
        }
      }
    }

    if (fixed1)
    {
      if (sgl_vt_[vhf])
      {
        if (std::abs(valence_[mesh_.edge_handle(he_open1)]) == 1)
        {
          bool suc = fis_.detach_singular_arc_from_singular_node(vhf,
                                                                 mesh_.opposite_halfedge_handle(he_open1),
                                                                 std::vector<HEH>{}, false, true);
          if (!suc)
            fis_.detach_singular_arcs_of_zipper_node(vhf, mesh_.opposite_halfedge_handle(he_open1),
                                                     std::vector<HEH>{}, false, true);
        }
        else
        {
          bool suc = fis_.detach_singular_arc_from_singular_node(vhf, false, true);
          if (!suc)
            fis_.detach_singular_arcs_of_zipper_node(vhf, false, true);
        }
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "Separated edge " << mesh_.edge_handle(_heh) << std::endl;)

    return true;
  }

  return false;
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::unzip_zipper_nodes()
{
  //split
  int n_fixed = 0, n_fixed_ctp = 0;
  std::set<VH> cd_vhs;
  for (const auto &vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi])
    {
      int n1 = 0, n_1 = 0, n_all = 0;
      for (auto ve_it = mesh_.ve_iter(vhi); ve_it.valid(); ++ve_it)
      {
        {
          if (valence_[*ve_it] == 1)
            n1++;

          if (valence_[*ve_it] == -1)
            n_1++;

          if (valence_[*ve_it] != 0)
            n_all++;
        }
      }

      if (!mesh_.is_boundary(vhi))
      {
        if (n1 == 1 && n_1 == 1 && n_all == 2)
        {
          cd_vhs.insert(vhi);
        }
      }
      else
      {
        if ((n1 == 1 || n_1 == 1) && n_all == 1)
        {
          if (ftp_.is_zipper_node(vhi))
            cd_vhs.insert(vhi);
        }
      }
    }
  }

  SplitHelperT<MeshT>::split_onering_of_zipper_nodes_arc(mesh_, es_, valence_, sgl_vt_, cd_vhs);


  std::set<VH> fixable_tps = ftp_.get_fixable_zipper_nodes(true);

  ALGOHEX_DEBUG_ONLY(std::cerr << "\nthere are: " << fixable_tps.size() << " zipper nodes" << std::endl;)
  ALGOHEX_DEBUG_ONLY(for (const auto &vh: fixable_tps)
                       std::cerr << " " << vh;
                             std::cerr << std::endl;)


  //in case of garbage collection
  auto is_tp_prop = mesh_.template request_vertex_property<bool>("zipper_node", false);
  auto length_prop = mesh_.template request_vertex_property<double>("tp_arc_length", 0.0);

  mesh_.set_persistent(length_prop, true);
  mesh_.set_persistent(is_tp_prop, true);
  for (auto vhi: fixable_tps)
    is_tp_prop[vhi] = true;


  //store zipper nodes' arc length
  for (const auto &vhi: fixable_tps)
  {
    HEH heh_s(-1);
    for (auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
    {
      if (valence_[mesh_.edge_handle(*voh_it)] != 0)
      {
        heh_s = *voh_it;
        break;
      }
    }
    if (!heh_s.is_valid())
      continue;

    auto sghes = get_halfedges_on_singular_arc(mesh_, valence_, heh_s);
    double len = 0.;
    for (const auto &hehi: sghes)
      len += mesh_.length(hehi);

    length_prop[vhi] = len;
  }


  //reset visited face property after each fix round
  auto visited_fprop = mesh_.template request_face_property<bool>("visited faces");
  for (const auto fhi: mesh_.faces())
    visited_fprop[fhi] = false;

  for (auto i = 0u; i < fixable_tps.size(); ++i)
  {
    std::set<VH> fixable_tps_new;

    //get tp to process
    VH vh_cur(-1);
    double max_len = -1., min_len = std::numeric_limits<double>::max();
    for (auto vhi: mesh_.vertices())
    {
      if (is_tp_prop[vhi])
      {
        if (fix_tps_on_longest_arc_first_)
        {
          if (length_prop[vhi] > max_len)
          {
            vh_cur = vhi;
            max_len = length_prop[vhi];
          }
        }
        else
        {
          if (length_prop[vhi] < min_len)
          {
            vh_cur = vhi;
            min_len = length_prop[vhi];
          }
        }

        fixable_tps_new.insert(vhi);
      }
    }

    if (!vh_cur.is_valid())
      break;

    ALGOHEX_DEBUG_ONLY(std::cerr << "\nfixing " << n_fixed << " tp/ctp: " << vh_cur << std::endl;)
    bool is_tp = ftp_.is_zipper_node(vh_cur);
    if (is_tp)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "\nfixing zipper node " << vh_cur << std::endl;)

      VH other_tp(-1);
      bool fixed = ftp_.unzip_zipper_node(vh_cur, HEH(-1), fixable_tps_new, other_tp);

      if (fixed)
      {
        if (other_tp.is_valid())
        {
          is_tp_prop[other_tp] = false;
        }

        n_fixed++;
      }
    }
    else
    {
      HEH sg_heh;
      if (ftp_.is_fixable_constrained_zipper_node(vh_cur, sg_heh))
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "\nfixing constrained zipper node " << vh_cur << std::endl;)

        VH other_tp(-1);
        bool fixed = ftp_.unzip_zipper_node(vh_cur, mesh_.opposite_halfedge_handle(sg_heh), fixable_tps_new, other_tp);

        if (fixed)
        {
          if (other_tp.is_valid())
          {
            is_tp_prop[other_tp] = false;
          }

          n_fixed++;
        }
      }
    }
    is_tp_prop[vh_cur] = false;

    if (n_fixed >= 500)
      break;
  }

  fac_.check_field_alignments_at_feature_faces();
  fac_.check_field_alignment_at_feature_edges();
  ALGOHEX_DEBUG_ONLY(std::cerr << "Fixed " << n_fixed << " zipper nodes" << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "Fixed " << n_fixed_ctp << " constrained zipper nodes" << std::endl;)

  mesh_.set_persistent(length_prop, false);
  mesh_.set_persistent(is_tp_prop, false);

  return n_fixed;
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::uniform_singular_arc_valence()
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Uniforming valence on singular arcs..." << std::endl;)
  auto flipped_valence = mesh_.template request_edge_property<bool>("flipped valence");
  mesh_.set_persistent(flipped_valence, true);

  std::set<VH> cd_vhs;
  for (const auto &vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi])
    {
      int nid = MeshPropertiesT<MeshT>::node_index(vhi);
      if (nid == 10 || nid == 12)
        cd_vhs.insert(vhi);
    }
  }

  int n_fixed = 0;
  while (!cd_vhs.empty())
  {
    auto vh_cur = *cd_vhs.begin();
    cd_vhs.erase(cd_vhs.begin());

    HEH heh(-1);
    for (auto voh_it = mesh_.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
      if (valence_[mesh_.edge_handle(*voh_it)] != 0)
      {
        heh = *voh_it;
        break;
      }

    if (!heh.is_valid())
      continue;

    auto sg_hehs = get_halfedges_on_singular_arc(mesh_, valence_, heh);
    if (sg_hehs.empty())
      continue;

    for (const auto hehi: sg_hehs)
    {
      auto ehi = mesh_.edge_handle(hehi);
      int old_val = valence_[ehi];
      int real_val = sge_.calculate_edge_valence(ehi);
      if (valence_[ehi] != old_val)
        std::cerr << "\nWarning: edge valence wrong! It is " << old_val << " should be " << real_val << " vhs "
                  << mesh_.edge(ehi).from_vertex() << " " << mesh_.edge(ehi).to_vertex() << std::endl;
    }

    int n_tp = 0;
    for (auto i = 0u; i < sg_hehs.size() - 1; ++i)
    {
      VH vht = mesh_.halfedge(sg_hehs[i]).to_vertex();
      int n_id = MeshPropertiesT<MeshT>::node_index(vht);
      if (n_id == 10 || n_id == 12)
      {
        cd_vhs.erase(vht);
        n_tp++;
      }
    }

    if (n_tp != 2 && n_tp > 0)
      std::cerr << "Warning: " << n_tp << " zipper nodes on singular arc!" << std::endl;

    int n_val_n = 0, n_val_p = 0;
    int abs_val = std::abs(valence_[mesh_.edge_handle(sg_hehs[0])]);
    double len_n = 0., len_p = 0.;
    for (const auto hehi: sg_hehs)
    {
      auto ehi = mesh_.edge_handle(hehi);
      if (valence_[ehi] == -abs_val)
        len_n += mesh_.length(ehi);
      else if (valence_[ehi] == abs_val)
        len_p += mesh_.length(ehi);
    }

    int val = len_n > len_p ? -abs_val : abs_val;

    //change valence
    for (const auto hehi: sg_hehs)
    {
      auto ehi = mesh_.edge_handle(hehi);
      if (valence_[ehi] == -val)
        flipped_valence[ehi] = true;
      valence_[ehi] = val;
    }

    n_fixed++;
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "Fixed " << n_fixed << std::endl;)

  return n_fixed;
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::restore_singular_arc_valence()
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Restore flipped valence on singular arcs." << std::endl;)
  auto flipped_valence = mesh_.template request_edge_property<bool>("flipped valence");

  for (const auto ehi: mesh_.edges())
  {
    if (flipped_valence[ehi])
    {
      valence_[ehi] = -valence_[ehi];
      flipped_valence[ehi] = false;
    }
  }

  mesh_.set_persistent(flipped_valence, false);
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::non_meshable_vertex_numbers(int &_n_nm_vts, int &_n_tps, int &_n_nm_nodes) const
{
  _n_nm_vts = 0, _n_tps = 0, _n_nm_nodes = 0;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] != 0)
    {
      int node_idx = MeshPropertiesT<MeshT>::node_index(vhi);

      if (node_idx == 20 && n_incident_singular_edges(mesh_, valence_, vhi) <= 4)
        _n_nm_nodes++;
      else if (node_idx == 10)
        _n_tps++;
    }
  }

  _n_nm_vts = _n_nm_nodes + _n_tps;
}

template<class MeshT>
bool
LocalMeshabilityRepairT<MeshT>::has_invalid_singular_edge() const
{
  for (const auto ehi: mesh_.edges())
  {
    if (valence_[ehi] == std::numeric_limits<int>::max() || valence_[ehi] == std::numeric_limits<int>::lowest())
    {
      return true;
    }
  }

  return false;
}


//TODO: support non-manifold
template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::fix_edges_with_non_meshable_footprint()
{
  int n_ce_removed = 0;
  std::set<VH> split_vhs;
  for (const auto eh: mesh_.edges())
    if ((valence_[eh] == -2 || valence_[eh] == 2) && feature_face_edge_[eh] && !keep_on_feature_face_[eh])
    {
      split_vhs.insert(mesh_.edge(eh).from_vertex());
      split_vhs.insert(mesh_.edge(eh).to_vertex());
    }

  for (const auto vhi: split_vhs)
    SplitHelperT<MeshT>::split_one_ring(mesh_, es_, vhi, valence_, feature_face_edge_, feature_face_vertex_,
                                        feature_edge_, sgl_vt_, feature_edge_vertex_, false);


  std::vector<EH> cpl_ehs;
  for (const auto eh: mesh_.edges())
    if ((valence_[eh] == -2 || valence_[eh] == 2) && feature_face_edge_[eh] && !keep_on_feature_face_[eh])
    {
      cpl_ehs.push_back(eh);
    }

  std::set<VH> vhs;
  for (const auto ehi: cpl_ehs)
  {
    auto fixed = eem_.fix_edge_with_non_meshable_footprint(ehi);
    if (fixed)
    {
      vhs.insert(mesh_.edge(ehi).from_vertex());
      vhs.insert(mesh_.edge(ehi).to_vertex());

      n_ce_removed++;
    }
  }

  for (const auto vhi: vhs)
  {
    bool suc = fis_.detach_singular_arc_from_singular_node(vhi, false, true);
    if (!suc)
      fis_.detach_singular_arcs_of_zipper_node(vhi, false, true);
  }

  mesh_.collect_garbage();

  ALGOHEX_DEBUG_ONLY(std::cerr << "Fixed " << n_ce_removed << " edges with non-meshable footprint." << std::endl;)
}


template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::
fix_compound_singular_edges(const bool _first_stage)
{
  std::vector<HEH> cp_hehs;
  for (const auto ehi: mesh_.edges())
  {
    if (valence_[ehi] == std::numeric_limits<int>::max() || std::abs(valence_[ehi]) >= 2)
    {
      if (mesh_.is_boundary(ehi))
        continue;

      if (_first_stage)
      {
        auto vh0 = mesh_.edge(ehi).from_vertex();
        auto vh1 = mesh_.edge(ehi).to_vertex();
        if (valence_[ehi] == std::numeric_limits<int>::max())
        {//try to shrink compound singular edges first
          if ((feature_face_vertex_[vh0] &&
               n_interior_complex_singular_edge(mesh_, valence_, vh0) >= 1) ||
              (feature_face_vertex_[vh1] && n_interior_complex_singular_edge(mesh_, valence_, vh1) >= 1))
            cp_hehs.push_back(mesh_.halfedge_handle(ehi, 0));
        }
        else if (std::abs(valence_[ehi]) >= 2 && valence_[ehi] != std::numeric_limits<int>::lowest() &&
                 !keep_on_feature_face_[ehi]) //high valence sges as a result of repair
          cp_hehs.push_back(mesh_.halfedge_handle(ehi, 0));
      }
      else if (valence_[ehi] ==
               std::numeric_limits<int>::max()) //in the second stage, we allow high valence as a result of fix
        cp_hehs.push_back(mesh_.halfedge_handle(ehi, 0));
    }
  }

  int n_fixed = 0;
  for (const auto hehi: cp_hehs)
    if (eem_.fix_interior_compound_singular_edge(hehi))
      n_fixed++;

  ALGOHEX_DEBUG_ONLY(std::cerr << "fixed " << n_fixed << " interior complex singular edges" << std::endl;)
  return n_fixed;
}


template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::fix_fully_constrained_parabolic_sectors()
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Fixing zero feature sectors...\n";)

  //preprocess: split singular edges that have two singular nodes
  SplitHelperT<MeshT>::split_feature_face_singular_edges_with_singular_nodes(mesh_, es_, feature_face_edge_, valence_);

  std::queue<VH> que;
  for (const auto &vh: mesh_.vertices())
  {
    if (feature_edge_vertex_[vh])
      que.push(vh);
  }

  //fix zero sectors
  int n_fixed = 0;
  while (!que.empty())
  {
    auto vh_cur = que.front();
    que.pop();

    bool suc = fix_fully_constrained_parabolic_sectors_at_vertex(vh_cur);
    if (suc)
    {
      n_fixed++;

      que.push(vh_cur);
    }
  }

  mesh_.collect_garbage();

  ALGOHEX_DEBUG_ONLY(std::cerr << "fixed " << n_fixed << " zero sectors" << std::endl;)

  return n_fixed;
}

template<class MeshT>
bool
LocalMeshabilityRepairT<MeshT>::fix_fully_constrained_parabolic_sectors_at_vertex(const VH _vh, const bool _print)
{
  if (_print)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Defect feature edges: " << std::endl;)
    for (auto &ehi: mesh_.edges())
    {
      if (feature_edge_[ehi] > 0 && !mesh_.is_boundary(ehi))
      {
        std::vector<HFH> hfs;
        HEH heh0 = mesh_.halfedge_handle(ehi, 0);
        for (auto hehf_it = mesh_.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
        {
          if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
          {
            hfs.push_back(*hehf_it);
          }
        }

        if (hfs.size() < 2)
          continue;

        int fsec_agl0 = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                            cell_quaternions_,
                                                                                            trans_prop_, heh0, hfs[0],
                                                                                            hfs[1]);
        ALGOHEX_DEBUG_ONLY(if (fsec_agl0 != 2)
                             std::cerr << " " << ehi;)
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "\nZero sectors fhs: " << std::endl;)
    for (const auto &vhi: mesh_.vertices())
    {
      if (feature_node_[vhi] && !mesh_.is_boundary(vhi))
      {
        std::vector<std::vector<HEH>> v_sec_hehs;
        std::vector<std::set<FH>> v_sec_fhs;
        get_feature_sectors_at_feature_edge_vertex(mesh_, feature_fprop_, feature_edge_,
                                                   feature_edge_vertex_, valence_, vhi, v_sec_hehs,
                                                   v_sec_fhs);

        for (auto j = 0u; j < v_sec_hehs.size(); ++j)
        {
          if (v_sec_hehs[j].size() < 2 || feature_edge_[mesh_.edge_handle(v_sec_hehs[j].front())] == 0 ||
              feature_edge_[mesh_.edge_handle(v_sec_hehs[j].back())] == 0)
            continue;

          int ea_st = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_,
                                                                                         trans_prop_,
                                                                                         valence_,
                                                                                         v_sec_hehs[j], v_sec_fhs[j]);

          ALGOHEX_DEBUG_ONLY(if (ea_st == 0)
                             {
                               for (auto fhi: v_sec_fhs[j])
                                 std::cerr << " " << fhi;
                             })
        }
      }
    }

    //debug
    {
      for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
      {
        auto ve = mesh_.edge_handle(*voh_it);
        if (feature_edge_[ve] > 0 || (valence_[ve] != 0 && feature_face_edge_[ve]))
        {
          std::vector<HFH> ft_hfhs;
          for (auto hehf_it = mesh_.hehf_iter(*voh_it); hehf_it.valid(); ++hehf_it)
          {
            if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
              ft_hfhs.push_back(*hehf_it);
          }

          double da0 = -1., da1 = -1.;
          int fsec_agl = -1000;
          if (!mesh_.is_boundary(ft_hfhs[0]))
          {
            fsec_agl = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                           cell_quaternions_,
                                                                                           trans_prop_, *voh_it,
                                                                                           ft_hfhs[0],
                                                                                           ft_hfhs[1]);
            da0 = dihedral_angle_from_halfface_to_halfface(mesh_, *voh_it, ft_hfhs[0], ft_hfhs[1]);
          }
          int fsec_agl2 = -1000;
          if (!mesh_.is_boundary(ft_hfhs[1]))
          {
            fsec_agl2 = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                            cell_quaternions_,
                                                                                            trans_prop_, *voh_it,
                                                                                            ft_hfhs[1],
                                                                                            ft_hfhs[0]);
            da1 = dihedral_angle_from_halfface_to_halfface(mesh_, *voh_it, ft_hfhs[1], ft_hfhs[0]);
          }


          ALGOHEX_DEBUG_ONLY(std::cerr << "Feature face sector angle0 " << fsec_agl
                                       << " ff sghe " << mesh_.halfedge(*voh_it).from_vertex() << " "
                                       << mesh_.halfedge(*voh_it).to_vertex()
                                       << " fs " << mesh_.face_handle(ft_hfhs[0]) << " fe "
                                       << mesh_.face_handle(ft_hfhs[1])
                                       << " val " << valence_[mesh_.edge_handle(*voh_it)]
                                       << " sec agl opp " << fsec_agl2 << " dihedral angle0 "
                                       << da0 << " da1" << da1
                                       << std::endl;)

        }
      }
    }
    ALGOHEX_DEBUG_ONLY(std::cerr << "\ndone" << std::endl;)
  }

  std::vector<std::vector<HEH>> zero_sectors_hehs;
  std::vector<std::set<FH>> zero_sectors_fhs;
  tu_.get_parabolic_feature_sectors_at_vertex(_vh, zero_sectors_hehs, zero_sectors_fhs);


  std::set<VH> zr_sec_vhs;
  if (!zero_sectors_hehs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "\nvertex: " << _vh << std::endl;)
    ALGOHEX_DEBUG_ONLY(std::cerr << "fixing sec: " << mesh_.edge_handle(zero_sectors_hehs[0][0]) << " "
                                 << mesh_.edge_handle(zero_sectors_hehs[0].back())
                                 << std::endl;)


    VH vhf = mesh_.halfedge(zero_sectors_hehs[0][0]).from_vertex();
    if (vhf != _vh)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "zero sector is not valid, continue..." << std::endl;)
      return false;
    }

    //if it doesn't create a constrained zipper node, leave it to second stage
    if (!sgl_vt_[_vh])
    {
      bool has_cv_orth = fac_.has_any_convexly_orthogonal_face(_vh, *zero_sectors_fhs[0].begin());
      if (!has_cv_orth && iter_ < first_iters_)
        return false;
    }


    ALGOHEX_DEBUG_ONLY(std::cerr << "\nfixing zero sector by creating a zipper node " << std::endl;)
    auto new_hehs = tu_.repair_fully_constrained_parabolic_sector(vhf, zero_sectors_hehs[0], zero_sectors_fhs[0]);

    if (!new_hehs.empty())
    {
      //detach
      for (const auto hehi: new_hehs)
      {
        bool dsa = fis_.detach_singular_arc_from_singular_node(_vh, hehi, std::vector<HEH>{},
                                                               false, true);
        if (!dsa)
          dsa = fis_.detach_singular_arcs_of_zipper_node(_vh, hehi, std::vector<HEH>{},
                                                         false, true);
      }

      if (mesh_.is_boundary(_vh))
        fis_.fix_boundary_invalid_singular_node(_vh, false);
    }
    else
      return false;

    return true;
  }


  return false;
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::fix_constrained_zipper_nodes()
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "########Flatten non-singular feature edges... " << std::endl;)
  //flaten edges such that constrained zipper nodes are at feature nodes
  flatten_non_singular_feature_edges();

  //
  std::vector<VH> cd_vhs;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] != 0 && feature_edge_vertex_[vhi])
    {
      if (feature_node_[vhi] == 0)
      {
        int n_sge = n_incident_singular_edges(mesh_, valence_, vhi);
        if (n_sge > 2)
          continue;
      }
      cd_vhs.push_back(vhi);
    }
  }

  int n_fixed = 0;
  for (const auto vhi: cd_vhs)
  {
    HEH heh(-1);
    int type = fctp_.fixable_constrained_zipper_node_type(vhi, heh);

    if (type > 0)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "\nFixing constrained zipper node at vh " << vhi << std::endl;)

      int nf = fctp_.fix_constrained_zipper_node(vhi);
      n_fixed += nf;
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "\nFixed " << n_fixed << " constrained zipper nodes!" << std::endl;)

  return n_fixed;
}

//remove singular triangle by changing the matching
template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::remove_singular_triangles()
{
  std::set<FH> sg_fhs;
  for (const auto eh: mesh_.edges())
  {
    if (valence_[eh] != 0)
    {
      for (auto ef_it = mesh_.ef_iter(eh); ef_it.valid(); ++ef_it)
      {
        if (is_removable_singular_triangle(*ef_it))
        {
          sg_fhs.insert(*ef_it);
        }
      }
    }
  }

  int n_removed = 0;
  while (!sg_fhs.empty() && n_removed < 10000)
  {
    auto fh_cur = *sg_fhs.begin();
    sg_fhs.erase(sg_fhs.begin());

    if (mesh_.is_boundary(fh_cur))
      continue;

    if (n_singular_edges_in_face(fh_cur) == 3)
    {
      HFH hf0 = mesh_.halfface_handle(fh_cur, 0);

      auto old_trans = trans_prop_[hf0];

      auto hes = mesh_.halfface(hf0).halfedges();

      std::vector<int> e_trans;

      for (const auto &hehi: hes)
        if (valence_[mesh_.edge_handle(hehi)] != 0)
          e_trans.push_back(EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, hehi, hf0));

      bool accept = false;
      for (const auto trans: e_trans)
      {
        trans_prop_[hf0] = tq_.mult_transitions_idx(old_trans, tq_.inverse_transition_idx(trans));
        trans_prop_[mesh_.opposite_halfface_handle(hf0)] = tq_.inverse_transition_idx(trans_prop_[hf0]);

        //check
        int n_regular = 0, n_complex = 0;
        for (const auto &heh: hes)
        {
          int val = sge_.calculate_edge_valence(mesh_.edge_handle(heh));
          if (val == 0)
            n_regular++;
          if (abs(val) > 1)
            n_complex++;
        }
        if (n_complex == 0 && n_regular == 3)
        {
          //check alignment of feature edge direction
          bool aligned = true;
          for (const auto hehi: hes)
          {
            if (feature_edge_[mesh_.edge_handle(hehi)] > 0)
            {
              if (!is_alignment_consistent_at_feature_edge(mesh_.edge_handle(hehi)))
              {
                aligned = false;
                break;
              }
            }
          }

          if (aligned)
          {
            accept = true;
            break;
          }
        }
      }

      if (!accept)
      {
        trans_prop_[hf0] = old_trans;
        trans_prop_[mesh_.opposite_halfface_handle(hf0)] = tq_.inverse_transition_idx(trans_prop_[hf0]);
      }
      else
      {
        for (const auto &heh: hes)
          sge_.compute_edge_valence(mesh_.edge_handle(heh));

        //update singular vertex
        for (auto hfv_it = mesh_.hfv_iter(hf0); hfv_it.valid(); ++hfv_it)
          MeshPropertiesT<MeshT>::update_singular_vertex_property(*hfv_it);

        //update quaternions
        std::set<EH> ehs;
        for (auto hfe_it = mesh_.hfe_iter(hf0); hfe_it.valid(); ++hfe_it)
          ehs.insert(*hfe_it);

        QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                           feature_fprop_, feature_edge_, ehs);

        n_removed++;

        //push neigbour
        for (const auto hehi: hes)
          for (auto hef_it = mesh_.hef_iter(hehi); hef_it.valid(); ++hef_it)
            if (is_removable_singular_triangle(*hef_it))
              sg_fhs.insert(*hef_it);
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "Fixed " << n_removed << " singular triangles" << std::endl;)


  std::vector<HEH> cd_hehs;
  for (auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] > 0 && feature_edge_vertex_[vhi] && !mesh_.is_boundary(vhi))
    {
      for (auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
      {
        auto ve = mesh_.edge_handle(*voh_it);
        if (feature_face_edge_[ve] && valence_[ve] == 0)
        {
          auto vh_to = mesh_.halfedge(*voh_it).to_vertex();
          if (sgl_vt_[vh_to] && !feature_edge_vertex_[vh_to])
          {
            cd_hehs.push_back(*voh_it);
          }
        }
      }
    }
  }

  std::vector<std::pair<HEH, std::vector<HFH>>> v_hehfs;
  for (auto hei: cd_hehs)
  {
    std::vector<HFH> cd_hfs;
    for (auto hehf_it = mesh_.hehf_iter(hei); hehf_it.valid(); ++hehf_it)
    {
      if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
        continue;

      auto he_n = mesh_.next_halfedge_in_halfface(hei, *hehf_it);
      auto he_nn = mesh_.next_halfedge_in_halfface(he_n, *hehf_it);
      if (std::abs(valence_[mesh_.edge_handle(he_n)]) == 1 && std::abs(valence_[mesh_.edge_handle(he_nn)]) == 1)
      {
        cd_hfs.push_back(*hehf_it);
      }
    }

    if (cd_hfs.size() == 2)
    {
      v_hehfs.emplace_back(hei, std::vector<HFH>{cd_hfs[0], cd_hfs[1]});
    }
  }


  int n_removed_quads = 0;
  for (const auto &hehfs: v_hehfs)
  {
    int old_trans0 = trans_prop_[hehfs.second[0]];
    auto he_n0 = mesh_.next_halfedge_in_halfface(hehfs.first, hehfs.second[0]);
    auto he_nn0 = mesh_.next_halfedge_in_halfface(he_n0, hehfs.second[0]);
    int e_trans0 = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, he_n0, hehfs.second[0]);
    trans_prop_[hehfs.second[0]] = tq_.mult_transitions_idx(old_trans0, tq_.inverse_transition_idx(e_trans0));
    trans_prop_[mesh_.opposite_halfface_handle(hehfs.second[0])] = tq_.inverse_transition_idx(
            trans_prop_[hehfs.second[0]]);

    int tt0 = sge_.calculate_edge_valence(mesh_.edge_handle(he_n0));
    if (tt0 != 0)
    {
      std::cerr << "Error: matching adjustment is wrong" << std::endl;
    }

    int old_trans1 = trans_prop_[hehfs.second[1]];
    auto he_n1 = mesh_.next_halfedge_in_halfface(hehfs.first, hehfs.second[1]);
    auto he_nn1 = mesh_.next_halfedge_in_halfface(he_n1, hehfs.second[1]);
    int e_trans1 = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, he_n1, hehfs.second[1]);
    trans_prop_[hehfs.second[1]] = tq_.mult_transitions_idx(old_trans1, tq_.inverse_transition_idx(e_trans1));
    trans_prop_[mesh_.opposite_halfface_handle(hehfs.second[1])] = tq_.inverse_transition_idx(
            trans_prop_[hehfs.second[1]]);

    int tt1 = sge_.calculate_edge_valence(mesh_.edge_handle(he_n1));
    if (tt1 != 0)
    {
      std::cerr << "Error: matching adjustment is wrong" << std::endl;
    }

    int t0 = sge_.calculate_edge_valence(mesh_.edge_handle(he_nn0));
    int t1 = sge_.calculate_edge_valence(mesh_.edge_handle(he_nn1));
    int t2 = sge_.calculate_edge_valence(mesh_.edge_handle(hehfs.first));

    if (t0 == 0 && t1 == 0 && t2 == 0)
    {
      //check sector angle at feature face edge
      if (fac_.check_feature_face_sectors_at_regular_edge(mesh_.edge_handle(hehfs.first)))
      {

        valence_[mesh_.edge_handle(he_n0)] = 0;
        valence_[mesh_.edge_handle(he_nn0)] = 0;
        valence_[mesh_.edge_handle(he_n1)] = 0;
        valence_[mesh_.edge_handle(he_nn1)] = 0;
        valence_[mesh_.edge_handle(hehfs.first)] = 0;

        //update singular vertex
        for (auto hfv_it = mesh_.hfv_iter(hehfs.second[0]); hfv_it.valid(); ++hfv_it)
          MeshPropertiesT<MeshT>::update_singular_vertex_property(*hfv_it);
        for (auto hfv_it = mesh_.hfv_iter(hehfs.second[1]); hfv_it.valid(); ++hfv_it)
          MeshPropertiesT<MeshT>::update_singular_vertex_property(*hfv_it);

        //update quaternions
        std::set<EH> ehs;
        ehs.insert(mesh_.edge_handle(he_n0));
        ehs.insert(mesh_.edge_handle(he_nn0));
        ehs.insert(mesh_.edge_handle(he_n1));
        ehs.insert(mesh_.edge_handle(he_nn1));
        ehs.insert(mesh_.edge_handle(hehfs.first));


        QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                           feature_fprop_,
                                                           feature_edge_, ehs);

        n_removed_quads++;
      }
      else
      {
        trans_prop_[hehfs.second[0]] = old_trans0;
        trans_prop_[mesh_.opposite_halfface_handle(hehfs.second[0])] = tq_.inverse_transition_idx(
                trans_prop_[hehfs.second[0]]);

        trans_prop_[hehfs.second[1]] = old_trans1;
        trans_prop_[mesh_.opposite_halfface_handle(hehfs.second[1])] = tq_.inverse_transition_idx(
                trans_prop_[hehfs.second[1]]);
      }
    }
    else
    {
      trans_prop_[hehfs.second[0]] = old_trans0;
      trans_prop_[mesh_.opposite_halfface_handle(hehfs.second[0])] = tq_.inverse_transition_idx(
              trans_prop_[hehfs.second[0]]);

      trans_prop_[hehfs.second[1]] = old_trans1;
      trans_prop_[mesh_.opposite_halfface_handle(hehfs.second[1])] = tq_.inverse_transition_idx(
              trans_prop_[hehfs.second[1]]);
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "Fixed " << n_removed_quads << " singular quads" << std::endl;)

}

template<class MeshT>
bool
LocalMeshabilityRepairT<MeshT>::is_removable_singular_triangle(const FH _fh) const
{
  if (n_singular_edges_in_face(_fh) == 3)
  {
//            std::vector<EH> ffsge, nffsge;
    for (auto fe_it = mesh_.fe_iter(_fh); fe_it.valid(); ++fe_it)
      if (feature_face_edge_[*fe_it] && !mesh_.is_boundary(*fe_it))
        return false;

    for (auto fv_it = mesh_.fv_iter(_fh); fv_it.valid(); ++fv_it)
    {
      if (feature_node_[*fv_it] && !mesh_.is_boundary(*fv_it))
        return false;
    }

    return true;
  }

  return false;
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::fix_misalignment_at_feature_edges()
{
  for (const auto ehi: mesh_.edges())
  {
    if (feature_edge_[ehi] > 0)
    {
      fix_misalignment_at_feature_edge(ehi);
    }
  }
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::fix_misalignment_at_feature_edge(const EH _eh)
{
  HEH heh0 = mesh_.halfedge_handle(_eh, 0);
  auto he_dir = mesh_.vertex(mesh_.halfedge(heh0).to_vertex()) - mesh_.vertex(mesh_.halfedge(heh0).from_vertex());
  he_dir.normalize();

  int fe_axis = -1, nm_axis = -1;
  HEHFIt hehf_it = mesh_.hehf_iter(heh0);
  for (auto hehf_it2 = mesh_.hehf_iter(heh0); hehf_it2.valid(); ++hehf_it2)
  {
    if (feature_fprop_[mesh_.face_handle(*hehf_it2)] > 0)
    {
      hehf_it = hehf_it2;
      break;
    }
  }
  if (feature_fprop_[mesh_.face_handle(*hehf_it)] == 0)
  {
    std::cerr << "Warning: no incident feature faces at the feature edge " << _eh << std::endl;
    return;
  }

  CH ch_s = mesh_.incident_cell(*hehf_it);
  if (!ch_s.is_valid())
    return;


  fe_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], he_dir).second;
  nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], mesh_.normal(*hehf_it)).second;

  ++hehf_it;
  std::vector<HFH> hf_changed;
  int trans_all = 0;
  for (; hehf_it.valid(); ++hehf_it)
  {
    CH chi = mesh_.incident_cell(*hehf_it);
    if (!chi.is_valid())
      continue;

    trans_all = tq_.mult_transitions_idx(trans_prop_[*hehf_it], trans_all);
    fe_axis = tq_.axis_after_transition(fe_axis, trans_prop_[*hehf_it]);
    int cur_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[chi], he_dir).second;

    if (cur_ax != fe_axis)
    {
      std::cerr << "Warning: feature aligned axis is not correct at edge " << _eh << " feax: " << fe_axis
                << " cur ax: " << cur_ax << std::endl;

      int temp_nm_axis = tq_.axis_after_transition(nm_axis, trans_all);
      int trans_i = -1;
      for (int i = 0; i < 24; ++i)
      {
        if (tq_.axis_after_transition(temp_nm_axis, i) == temp_nm_axis
            && tq_.axis_after_transition(fe_axis, i) == cur_ax)
        {
          trans_i = i;
          break;
        }
      }

      if (trans_i != -1)
      {
        trans_all = tq_.mult_transitions_idx(trans_i, trans_all);
        trans_prop_[*hehf_it] = tq_.mult_transitions_idx(trans_i, trans_prop_[*hehf_it]);
        trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = tq_.inverse_transition_idx(trans_prop_[*hehf_it]);
        fe_axis = cur_ax;
        hf_changed.push_back(*hehf_it);
        ALGOHEX_DEBUG_ONLY(std::cerr << "Fixed, transi " << trans_i << std::endl;)
      }
    }
  }

  std::set<EH> ehs;
  for (const auto hfhi: hf_changed)
    for (auto hfe_it = mesh_.hfe_iter(hfhi); hfe_it.valid(); ++hfe_it)
      ehs.insert(*hfe_it);

  for (const auto ehi: ehs)
    sge_.compute_edge_valence(ehi);

  std::set<VH> vhs;
  for (const auto ehi: ehs)
  {
    auto vhf = mesh_.edge(ehi).from_vertex();
    auto vht = mesh_.edge(ehi).to_vertex();
    vhs.insert(vhf);
    vhs.insert(vht);
  }

  for (const auto vhi: vhs)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);
}


template<class MeshT>
bool
LocalMeshabilityRepairT<MeshT>::is_alignment_consistent_at_feature_edge(const EH _eh) const
{
  HEH heh0 = mesh_.halfedge_handle(_eh, 0);
  auto he_dir = mesh_.vertex(mesh_.halfedge(heh0).to_vertex()) - mesh_.vertex(mesh_.halfedge(heh0).from_vertex());
  he_dir.normalize();

  int fe_axis = -1, nm_axis = -1;
  HEHFIt hehf_it = mesh_.hehf_iter(heh0);
  for (auto hehf_it2 = mesh_.hehf_iter(heh0); hehf_it2.valid(); ++hehf_it2)
  {
    if (feature_fprop_[mesh_.face_handle(*hehf_it2)] > 0)
    {
      hehf_it = hehf_it2;
      break;
    }
  }
  CH ch_s = mesh_.incident_cell(*hehf_it);
  fe_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], he_dir).second;

  ++hehf_it;
  std::vector<HFH> hf_changed;
  for (; hehf_it.valid(); ++hehf_it)
  {
    CH chi = mesh_.incident_cell(*hehf_it);
    if (!chi.is_valid())
      continue;

    fe_axis = tq_.axis_after_transition(fe_axis, trans_prop_[*hehf_it]);
    int cur_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[chi], he_dir).second;

    if (cur_ax != fe_axis)
      return false;
  }

  return true;
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::
fix_misalignment_at_feature_face_sectors()
{
  std::vector<EH> ft_ehs;
  for (const auto ehi: mesh_.edges())
    if (feature_face_edge_[ehi] && !mesh_.is_boundary(ehi))
      ft_ehs.push_back(ehi);

  for (const auto ehi: ft_ehs)
    fix_misalignment_at_feature_face_sector(ehi);
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::
fix_misalignment_at_feature_face_sector(const EH _eh)
{
  HEH heh0 = mesh_.halfedge_handle(_eh, 0);
  std::vector<HFH> ft_hfhs;
  for (auto hehf_it = mesh_.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
  {
    if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
      ft_hfhs.push_back(*hehf_it);
  }
  if (ft_hfhs.size() != 2)
  {
    std::cerr << "Error: non-manifold feature surface at edge!" << std::endl;
    return;
  }

  fix_misalignment_at_feature_face_sector(heh0, ft_hfhs[0], ft_hfhs[1]);
  fix_misalignment_at_feature_face_sector(heh0, ft_hfhs[1], ft_hfhs[0]);
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::
fix_misalignment_at_feature_face_sector(const HEH _heh, const HFH _hfh_s, const HFH _hfh_e)
{
  if (mesh_.is_boundary(_hfh_s))
    return;

  int fsec_agl = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                     trans_prop_, _heh, _hfh_s, _hfh_e);
  EH eh = mesh_.edge_handle(_heh);
  if ((fsec_agl == 1 || fsec_agl == 3))
  {
    if (feature_edge_[eh] > 0)
    {
      auto da = dihedral_angle_from_halfface_to_halfface(mesh_, _heh, _hfh_s, _hfh_e);
      if ((da < M_PI && fsec_agl == 1) || (da > M_PI && fsec_agl == 3)) //valid case
        return;
    }
    //get rotation axis
    HFH hf_start = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(_hfh_s, _heh));

    CH ch0 = mesh_.incident_cell(_hfh_s);
    Point n0 = mesh_.normal(mesh_.opposite_halfface_handle(_hfh_s));
    int closest_n0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch0], n0).second;

    CH ch1 = mesh_.incident_cell(mesh_.opposite_halfface_handle(_hfh_e));
    Point n1 = mesh_.normal(_hfh_e);
    int closest_n1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch1], n1).second;

    int half_etrans = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, _heh, hf_start, _hfh_e);
    int T_inv_id = tq_.inverse_transition_idx(half_etrans);

    int n1_in_q0 = tq_.axis_after_transition(closest_n1, T_inv_id);

    int rt_axis = the_third_axis((AxisAlignment) closest_n0, (AxisAlignment) n1_in_q0);

    int rt_trans = rt_axis + 4;

    int new_trans = tq_.mult_transitions_idx(trans_prop_[hf_start], rt_trans);

    bool allow_cpe = !mesh_.is_boundary(eh) && feature_edge_[eh] == 0;
    eem_.attempt_matching_adjustment_of_feature_face_sector(_heh, _hfh_s, _hfh_e, new_trans, allow_cpe);
  }
  else if (fsec_agl == 0 || fsec_agl == 4)
  {
    if (feature_edge_[eh] > 0)
    {
      CH ch_i = mesh_.incident_cell(_hfh_s);
      auto dir_i = mesh_.vertex(mesh_.halfedge(_heh).to_vertex()) -
                   mesh_.vertex(mesh_.halfedge(_heh).from_vertex());
      int ax_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_i],
                                                    dir_i.normalize()).second;
      if (fsec_agl == 0)
        ax_i = negate(ax_i);

      int rt_trans = ax_i + 4;

      HFH hf_start = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(_hfh_s, _heh));

      int new_trans = tq_.mult_transitions_idx(trans_prop_[hf_start], rt_trans);
      eem_.attempt_matching_adjustment_of_feature_face_sector(_heh, _hfh_s, _hfh_e, new_trans);
    }
  }
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::matching_at_boundary_vertex(const VH _vh) const
{
  HEH heh_s(-1);
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    if (feature_edge_[mesh_.edge_handle(*voh_it)] > 0)
    {
      heh_s = *voh_it;
      break;
    }

  if (!heh_s.is_valid())
    heh_s = *mesh_.voh_iter(_vh);

  int trans = 0;
  HEH he_it = heh_s;
  while (1)
  {
    auto secs = get_boundary_halfedges_in_feature_sector_ccw(mesh_, feature_edge_, valence_, he_it);

    for (auto i = 0u; i < secs.size() - 1; ++i)
    {
      for (auto hehf_it = mesh_.hehf_iter(secs[i]); hehf_it.valid(); ++hehf_it)
        if (!mesh_.is_boundary(mesh_.face_handle(*hehf_it)))
          trans = tq_.mult_transitions_idx(trans_prop_[*hehf_it], trans);
    }

    he_it = secs.back();
    if (he_it == heh_s)
      break;
  }


  return trans;
}


template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::
merge_zipper_nodes()
{
  ec_.set_collapse_with_optimization(true);

  std::set<VH> fixable_tps = ftp_.get_fixable_zipper_nodes(false);
  ALGOHEX_DEBUG_ONLY(std::cerr << "Mergin tps, fixable tp size " << fixable_tps.size() << std::endl;)

  int n_mg = 0;
  for (const auto vhi: fixable_tps)
  {
    if (!feature_face_vertex_[vhi])
    {
      bool mg = ftp_.merge_zipper_node(vhi);
      if (mg)
        n_mg++;
    }
  }

  ec_.set_collapse_with_optimization(false);

  mesh_.collect_garbage();

  ALGOHEX_DEBUG_ONLY(std::cerr << "Merged " << n_mg << " zipper nodes" << std::endl;)

  return n_mg;
}

template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::visualize_dominant_axis_on_singular_arc(const HEH _heh, std::map<CH, int> &_c_axis,
                                                                        MeshT &_pm)
{
  VH vh0 = mesh_.halfedge(_heh).from_vertex();
  double scale = target_length_[vh0];

  auto sg_hehs = get_halfedges_on_singular_arc(mesh_, valence_, _heh);

  std::set<CH> or_cells;
  for (auto i = 0u; i < sg_hehs.size(); ++i)
  {
    for (auto vc_it = mesh_.vc_iter(mesh_.halfedge(sg_hehs[i]).from_vertex()); vc_it.valid(); ++vc_it)
      or_cells.insert(*vc_it);
  }
  for (auto vc_it = mesh_.vc_iter(mesh_.halfedge(sg_hehs.back()).to_vertex()); vc_it.valid(); ++vc_it)
    or_cells.insert(*vc_it);

  std::set<EH> ehs;
  for (auto chi: or_cells)
  {
    for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
    {
      VH vh0 = mesh_.edge(*ce_it).from_vertex();
      VH vh1 = mesh_.edge(*ce_it).to_vertex();
      if (valence_[*ce_it] == 0)
      {
        if (sgl_vt_[vh1] && sgl_vt_[vh0])
        {
          ehs.insert(*ce_it);
        }
      }
    }
  }

  for (const auto &ehi: ehs)
    if (!mesh_.is_deleted(ehi))
    {
      es_.edge_split(ehi);
    }
  mesh_.collect_garbage();

  HEH heh_on(-1);
  for (auto voh_it = mesh_.voh_iter(vh0); voh_it.valid(); ++voh_it)
  {
    if (valence_[mesh_.edge_handle(*voh_it)] != 0)
    {
      heh_on = *voh_it;
      break;
    }
  }
  sg_hehs = get_halfedges_on_singular_arc(mesh_, valence_, heh_on);

  //check
  std::map<CH, int> c_axis;
  auto hehf_it = mesh_.hehf_iter(sg_hehs[0]);
  CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
  int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   sg_hehs[0], *hehf_it);
  int val = valence_[mesh_.edge_handle(sg_hehs[0])];
  if (val == 1)
    rt_axis = negate(rt_axis);
  c_axis[ch_s] = rt_axis;

  for (; hehf_it.valid(); ++hehf_it)
  {
    auto chi = mesh_.incident_cell(*hehf_it);
    if (chi == ch_s)
      break;

    CH chopp = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
    if (c_axis.find(chopp) == c_axis.end())
      std::cerr << "Error: couldn't find start cell " << chopp << std::endl;

    c_axis[chi] = tq_.axis_after_transition(c_axis[chopp], trans_prop_[*hehf_it]);
  }

  for (auto i = 0u; i < sg_hehs.size() - 1; ++i)
  {
    CH ch_prev = *mesh_.hec_iter(sg_hehs[i]);
    int ax_prev = c_axis[ch_prev];

    VH vh_t = mesh_.halfedge(sg_hehs[i]).to_vertex();

    std::set<CH> v_or_cells;
    for (auto vc_it = mesh_.vc_iter(vh_t); vc_it.valid(); ++vc_it)
      v_or_cells.insert(*vc_it);

    auto hfhi_it = mesh_.hehf_iter(sg_hehs[i + 1]);
    CH chopp_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hfhi_it));

    int axis_i = fac_.axis_in_chs_expressed_in_cht(ch_prev, chopp_s, ax_prev, v_or_cells);
    c_axis[chopp_s] = axis_i;

    for (; hfhi_it.valid(); ++hfhi_it)
    {
      auto chi = mesh_.incident_cell(*hfhi_it);
      if (chi == chopp_s)
        break;

      CH chopp = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hfhi_it));
      if (c_axis.find(chopp) == c_axis.end())
        std::cerr << "Error: couldn't find cell " << chopp << std::endl;

      c_axis[chi] = tq_.axis_after_transition(c_axis[chopp], trans_prop_[*hfhi_it]);
    }
  }
  _c_axis = c_axis;


  for (auto i = 0u; i < sg_hehs.size(); ++i)
  {
    VH vh_f = mesh_.halfedge(sg_hehs[i]).from_vertex();
    VH vh_t = mesh_.halfedge(sg_hehs[i]).to_vertex();
    Point mid_pt = (mesh_.vertex(vh_f) + mesh_.vertex(vh_t)) / 2.;

    CH ch_inc = *mesh_.hec_iter(sg_hehs[i]);
    auto dm_axis = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_inc], (AxisAlignment) c_axis[ch_inc]);


    auto vm = _pm.add_vertex(mid_pt);
    auto v1 = _pm.add_vertex(mid_pt + Point(dm_axis[0], dm_axis[1], dm_axis[2]) * scale);
    _pm.add_edge(vm, v1);
  }
}

template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::zip_zipper_nodes()
{
  std::vector<VH> cd_vhs;
  for (const auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] && !mesh_.is_boundary(vhi) && MeshPropertiesT<MeshT>::node_index(vhi) == 10)
      cd_vhs.push_back(vhi);
  }

  int nrm = 0;
  for (const auto vhi: cd_vhs)
  {
    bool rm = ftp_.zip_zipper_node(vhi);
    if (rm)
      nrm++;
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "Removed " << nrm << " zipper nodes" << std::endl;)

  return nrm;
}


template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::n_singular_edges_in_face(const FH _fh) const
{
  int n_sg = 0;
  for (auto fe_it = mesh_.fe_iter(_fh); fe_it.valid(); ++fe_it)
    if (valence_[*fe_it] != 0)
      n_sg++;

  return n_sg;
}


template<class MeshT>
int
LocalMeshabilityRepairT<MeshT>::n_incident_interior_singular_edges(const VH _vh) const
{
  int n_sg = 0;
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (valence_[*ve_it] != 0 && !feature_face_edge_[*ve_it])
      n_sg++;

  return n_sg;
}


template<class MeshT>
void
LocalMeshabilityRepairT<MeshT>::export_local_configuration(const VH vh, const std::string &filename) const
{
  MeshT tetmesh;

  auto eval = tetmesh.template request_edge_property<int>("edge_valance");
  tetmesh.set_persistent(eval, true);

  auto vfeat = tetmesh.template request_vertex_property<int>("AlgoHex::FeatureVertices");
  tetmesh.set_persistent(vfeat, true);

  auto efeat = tetmesh.template request_edge_property<int>("AlgoHex::FeatureEdges");
  tetmesh.set_persistent(efeat, true);

  auto ffeat = tetmesh.template request_face_property<int>("AlgoHex::FeatureFaces");
  tetmesh.set_persistent(ffeat, true);

  auto hftrans = tetmesh.template request_halfface_property<int>("HalffaceTransiton");
  tetmesh.set_persistent(hftrans, true);

  auto qtns = tetmesh.template request_cell_property<Quaternion>("FrameFieldQuaternions");
  tetmesh.set_persistent(qtns, true);


  std::map<VH, VH> vm;
  std::map<VH, VH> vm_inv;

  // generate mesh and set cell properties
  for (auto vc_it = mesh_.vc_iter(vh); vc_it.valid(); ++vc_it)
  {
    CH ch = *vc_it;
    std::vector<VH> vhs;
    TVIt tv_it(ch, &mesh_);
    for (; tv_it.valid(); ++tv_it)
    {
      auto m_it = vm.find(*tv_it);
      if (m_it == vm.end())
      {
        VH vh_new = tetmesh.add_vertex(mesh_.vertex(*tv_it));
        vm[*tv_it] = vh_new;
        vm_inv[vh_new] = *tv_it;
        vhs.push_back(vh_new);
      }
      else vhs.push_back(m_it->second);
    }

    CH ch_new = tetmesh.add_cell(vhs);

    // store initial frame
    qtns[ch_new] = cell_quaternions_[*vc_it];
  }

  // set vertex properties
  for (VIt v_it = tetmesh.vertices_begin(); v_it != tetmesh.vertices_end(); ++v_it)
  {
    VH vh_orig = vm_inv[*v_it];
    vfeat[*v_it] = feature_node_[vh_orig];
  }

  // set edge properties
  for (EIt e_it = tetmesh.edges_begin(); e_it != tetmesh.edges_end(); ++e_it)
  {
    // get original edge handle
    EH eh = *e_it;
    VH vh0 = tetmesh.halfedge(mesh_.halfedge_handle(eh, 0)).to_vertex();
    VH vh1 = tetmesh.halfedge(mesh_.halfedge_handle(eh, 1)).to_vertex();
    VH vho0 = vm_inv[vh0];
    VH vho1 = vm_inv[vh1];
    HEH heh = mesh_.halfedge(vho0, vho1);
    if (!heh.is_valid())
      std::cerr << "ERROR: could not obtain halfedge of original mesh in export_local_configuration!!!" << std::endl;
    EH eho = mesh_.edge_handle(heh);

    // copy valence property
    eval[*e_it] = valence_[eho];
    efeat[*e_it] = feature_edge_[eho];
  }

  // set face properties
  for (FIt f_it = tetmesh.faces_begin(); f_it != tetmesh.faces_end(); ++f_it)
  {
    FH fh = *f_it;
    HFH hfh0 = tetmesh.halfface_handle(fh, 0);
    HFH hfh1 = tetmesh.halfface_handle(fh, 1);

    auto f0 = tetmesh.halfface(hfh0);
    std::vector<VH> vhs0;
    vhs0.push_back(tetmesh.halfedge(f0.halfedges()[0]).to_vertex());
    vhs0.push_back(tetmesh.halfedge(f0.halfedges()[1]).to_vertex());
    vhs0.push_back(tetmesh.halfedge(f0.halfedges()[2]).to_vertex());

    // map vertex indices
    std::vector<VH> vhs0_orig;
    vhs0_orig.push_back(vm_inv[vhs0[0]]);
    vhs0_orig.push_back(vm_inv[vhs0[1]]);
    vhs0_orig.push_back(vm_inv[vhs0[2]]);

    // get corresponding halfface in original mesh
    HFH hfh0_orig = mesh_.halfface(vhs0_orig);
    if (!hfh0_orig.is_valid())
      std::cerr << "ERROR: could not obtain original halfface in export_local_configuration!!!" << std::endl;
    HFH hfh1_orig = mesh_.opposite_halfface_handle(hfh0_orig);
    FH fh_orig = mesh_.face_handle(hfh0_orig);

    // copy properties
    ffeat[fh] = feature_fprop_[fh_orig];
    hftrans[hfh0] = trans_prop_[hfh0_orig];
    hftrans[hfh1] = trans_prop_[hfh1_orig];
  }

  // write file
  OpenVolumeMesh::IO::FileManager fm;
  fm.writeFile(filename + ".ovm", tetmesh);

  //
  std::ofstream f_write(filename + ".qtn");
  for (const auto ch: tetmesh.cells())
  {
    // write to file
    f_write << qtns[ch].w() << " ";
    f_write << qtns[ch].x() << " ";
    f_write << qtns[ch].y() << " ";
    f_write << qtns[ch].z() << " ";
  }
  f_write.close();
}

}

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "PushSingularVertexT.hh"

namespace AlgoHex
{
template<class MeshT>
bool
PushSingularVertexT<MeshT>::push_singular_vertex(const VH _vh)
{
  SplitHelperT<MeshT>::split_one_ring(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                      feature_edge_, sgl_vt_, feature_edge_vertex_);

  //map singular edge handle to transition index
  std::map<HEH, int> e_trans_idxes, hfh_to_trans_idx;
  bool success = get_original_boundary_halfedge_trans_idx(_vh, e_trans_idxes);

  //not suitable for pushing inside
  if (!success)
  {
    ALGOHEX_DEBUG_ONLY(std::cout << "skip vertex: " << _vh << std::endl;)
    return false;
  }

  int sge_stautus = check_new_edge_index(_vh, e_trans_idxes, hfh_to_trans_idx);

//        if(sge_stautus == 1 && feature_node_[_vh]) {
//            std::cerr<<"new valence -2 or +2 singular edge at singular node, skip vertex!"<<std::endl;
//            return false;
//        }
  if (sge_stautus == 2)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "new valence 10 singular edge will be created, skip vertex!" << std::endl;)
    return false;
  }

  auto vh_new = add_vertex_boundary_buffer(_vh);

  if (!vh_new.is_valid())
    return false;

//        determine_matchings(vh_new, _vh, e_trans_idxes);
  determine_matchings(vh_new, _vh, e_trans_idxes, hfh_to_trans_idx);

  //update edge valence
  std::set<EH> rel_ehs;
  for (auto ve_it = mesh_.ve_iter(vh_new); ve_it.valid(); ++ve_it)
    rel_ehs.insert(*ve_it);
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
    rel_ehs.insert(*ve_it);

  for (const auto &eh: rel_ehs)
    sge_.compute_edge_valence(eh);

  //update singular vertex property
  for (auto vv_it = mesh_.vv_iter(vh_new); vv_it.valid(); ++vv_it)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(*vv_it);

  for (auto vv_it = mesh_.vv_iter(_vh); vv_it.valid(); ++vv_it)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(*vv_it);


  //update feature edge property
  auto heh_new = mesh_.find_halfedge(_vh, vh_new);
  for (auto hehf_it = mesh_.hehf_iter(heh_new); hehf_it.valid(); ++hehf_it)
  {
    auto prev_heh = mesh_.prev_halfedge_in_halfface(heh_new, *hehf_it);
    auto next_heh = mesh_.next_halfedge_in_halfface(heh_new, *hehf_it);

    feature_edge_[mesh_.edge_handle(next_heh)] = feature_edge_[mesh_.edge_handle(prev_heh)];
    feature_edge_[mesh_.edge_handle(prev_heh)] = 0;

    if (mesh_.is_boundary(mesh_.edge_handle(prev_heh)))
      std::cerr << "Error: edge " << mesh_.edge_handle(prev_heh) << " should be in the interior" << std::endl;
  }

  //update feature edge vertex property
  if (feature_edge_vertex_[_vh])
  {
    feature_edge_vertex_[vh_new] = true;
    feature_edge_vertex_[_vh] = false;
  }

  //update feature node property
  if (feature_node_[_vh] > 0)
  {
    feature_node_[vh_new] = feature_node_[_vh];
    feature_node_[_vh] = false;
  }

  QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                      feature_fprop_, feature_edge_, vh_new);

  return true;
}

template<class MeshT>
bool
PushSingularVertexT<MeshT>::is_move_inside_candidate(const VH _vh, bool _push_feature_vertex, bool _preprocess)
{
  if (mesh_.is_deleted(_vh))
    return false;

  if (sgl_vt_[_vh] && mesh_.is_boundary(_vh))
  {
    std::vector<HEH> ft_hehs, bdy_hehs, itr_hehs;
    get_special_halfedges_at_vertex(mesh_, feature_edge_, valence_, _vh, ft_hehs, bdy_hehs, itr_hehs);


    int n_ft = ft_hehs.size();

    if (itr_hehs.size() == 1 && bdy_hehs.size() == 0)
      return false;

    //push singular vts on non-feature boundary or invalid sg arc vts on feature arc to the interior
    if (_preprocess)
    {
      if (!_push_feature_vertex && feature_edge_vertex_[_vh])
        return false;

      if (n_ft == 0)
        return true;
      else if (n_ft == 2)
      {
        if ((bdy_hehs.size() == 0 && itr_hehs.size() == 2))
          return true;
      }

      return false;
    }
    else
    {
      if (!_push_feature_vertex && feature_edge_vertex_[_vh])
        return false;

      //in case of zipper nodes after fixing zero sectors
      if (feature_node_[_vh] > 0 && bdy_hehs.size() == 0 && itr_hehs.size() >= 2)
      {
        if (n_interior_complex_singular_edge(mesh_, valence_, _vh) > 0)
          return false;

        std::set<CH> or_chs;
        for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
          or_chs.insert(*vc_it);
        bool has_prl = has_interior_singular_edge_tangent_to_surface(_vh, or_chs, itr_hehs);
        if (has_prl)
          return true;
        else
          return false;
      }

      //if invalid
      if (itr_hehs.size() > 0)
      {
        //including sg arc orthogonally touch ff
        if (feature_node_[_vh] == 0 && itr_hehs.size() == 2 && bdy_hehs.size() == 0)
        {
          return true;
        }

        std::set<CH> or_chs;
        for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
          or_chs.insert(*vc_it);
        bool is_prl = has_interior_singular_edge_tangent_to_surface(_vh, or_chs, itr_hehs);
        if (is_prl)
        {
          if (!fis_.is_locally_meshable(_vh))
            return true;
        }
      }
    }
  }

  return false;
}

template<class MeshT>
bool
PushSingularVertexT<MeshT>::has_interior_singular_edge_tangent_to_surface(const VH _vh, const std::set<CH> &_or_chs,
                                                                          const std::vector<HEH> &_itr_hehs)
{
  //check if any interior singular edge is parallel to surface
  for (const auto &hehi: _itr_hehs)
  {
    for (auto vhf_it = mesh_.vhf_iter(_vh); vhf_it.valid(); ++vhf_it)
    {
      if (mesh_.is_boundary(*vhf_it))
      {
        int is_prl = is_singular_edge_parallel_to_normal_direction(
                mesh_.opposite_halfedge_handle(hehi), *vhf_it, _or_chs);
        if (is_prl == 0)
        {
          return true;
        }
      }
    }
  }

  return false;
}

template<class MeshT>
int
PushSingularVertexT<MeshT>::is_singular_edge_parallel_to_normal_direction(const HEH _sg_heh, const HFH _hfh,
                                                                          const std::set<CH> _onering_chs)
{
  HFH hfh_s = *mesh_.hehf_iter(_sg_heh);
  int e_rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                     valence_, _sg_heh, hfh_s);

  if (valence_[mesh_.edge_handle(_sg_heh)] == 1)
    e_rt_axis = negate(e_rt_axis);

  CH ch_sg = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh_s));

  CH ch_bdy = mesh_.incident_cell(mesh_.opposite_halfface_handle(_hfh));
  auto nm = mesh_.normal(_hfh);
  int nm_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_bdy], nm).second;

  VH vh_f = mesh_.halfedge(_sg_heh).from_vertex();
  int nm_in_sgch = fac_.axis_in_chs_expressed_in_cht(ch_bdy, ch_sg, nm_ax, _onering_chs);

  if (nm_in_sgch / 2 == e_rt_axis / 2)
  {
    if (nm_in_sgch != e_rt_axis) //opposite direction
      return 1;
    else //same direction
      return 2;
  }

  return 0;
}

template<class MeshT>
bool
PushSingularVertexT<MeshT>::
get_original_boundary_halfedge_trans_idx(const VH _vh, std::map<HEH, int> &e_to_idx)
{
  for (auto vih_it = mesh_.vih_iter(_vh); vih_it.valid(); ++vih_it)
  {
    auto eh = mesh_.edge_handle(*vih_it);
    if (mesh_.is_boundary(*vih_it) && valence_[eh] != 0)
    {
      if (valence_[eh] == -2 || valence_[eh] == 2)
      {
        std::cout << "Warning: complex boundary singular edge in the one ring neigbourhood." << std::endl;
        return false;
      }

      //find start halfface
      HFH hf_start = *mesh_.hehf_iter(*vih_it),
              hf_start_opp = mesh_.opposite_halfface_handle(hf_start),
              hf_end = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(mesh_.opposite_halfedge_handle(*vih_it)));

      int idx = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, *vih_it, hf_start);

      //compute transition product
      CH ch0 = mesh_.incident_cell(hf_start);
      Point n0 = mesh_.normal(hf_start_opp);
      int closest_n0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch0], n0).second;

      CH ch1 = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_end));
      Point n1 = mesh_.normal(hf_end);
      int closest_n1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch1], n1).second;
      int n1_in_q0 = tq_.axis_after_transition(closest_n1, tq_.inverse_transition_idx(idx));

      if (n1_in_q0 / 2 == closest_n0 / 2)
      {
        VH vh_f = mesh_.edge(eh).from_vertex();
        VH vh_t = mesh_.edge(eh).to_vertex();

        if (n_incident_singular_edges(mesh_, valence_, vh_f) == 1 ||
            n_incident_singular_edges(mesh_, valence_, vh_t) == 1)
        {
          sge_.compute_edge_valence(eh);
          MeshPropertiesT<MeshT>::update_singular_vertex_property(vh_f);
          MeshPropertiesT<MeshT>::update_singular_vertex_property(vh_t);
        }
        else
          std::cerr << "Warning: should be a complex boundary singular edge or regular edge! Current valence of edge "
                    << eh << "is not consistent!" << std::endl;

        return false;
      }

      //boundary edge rotation axis to transition index
      int rt_axis = the_third_axis((AxisAlignment) n1_in_q0, (AxisAlignment) closest_n0);
      int edge_idx = rt_axis + 4;
      e_to_idx[*vih_it] = edge_idx;
    }
  }


  return true;
}

template<class MeshT>
int
PushSingularVertexT<MeshT>::check_new_edge_index(const VH _vh, std::map<HEH, int> &_heh_to_trans_idx,
                                                 std::map<HEH, int> &_hfh_to_trans_idx)
{
  //calculate halfface transition index
  for (auto vih_it = mesh_.vih_iter(_vh); vih_it.valid(); ++vih_it)
  {
    if (mesh_.is_boundary(*vih_it))
    {
      int interior_trans = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, *vih_it,
                                                                           *mesh_.hehf_iter(*vih_it));
      //overall transition of the halfedge if starting from current cell
      int edge_idx = 0;
      if (valence_[mesh_.edge_handle(*vih_it)] != 0)
        edge_idx = _heh_to_trans_idx[*vih_it];

      //calculate the transition
      _hfh_to_trans_idx[*vih_it] = tq_.mult_transitions_idx(edge_idx, tq_.inverse_transition_idx(interior_trans));
    }
  }

  //circulate incoming halfedges in ccw
  HEH he_s;
  for (auto vih_it = mesh_.vih_iter(_vh); vih_it.valid(); ++vih_it)
    if (mesh_.is_boundary(*vih_it))
    {
      he_s = *vih_it;
      break;
    }

  assert(he_s.is_valid());

  HEH he_it = he_s;
  int trans = 0;
  do
  {
    trans = tq_.mult_transitions_idx(_hfh_to_trans_idx[he_it], trans);

    HFH hf_bdy = *mesh_.hehf_iter(he_it);
    he_it = mesh_.next_halfedge_in_halfface(he_it, hf_bdy);
    he_it = mesh_.opposite_halfedge_handle(he_it);
  }
  while (he_it != he_s);

  if ((trans < 10 && trans >= 4) || trans == 0)
    return 0;
  else if (trans < 4 && trans >= 1)
  {
    return 1;
  }


  return 2;
}

template<class MeshT>
VH
PushSingularVertexT<MeshT>::add_vertex_boundary_buffer(const VH _vh)
{
  if (!mesh_.is_boundary(_vh))
    return VH(-1);

  auto pt_old = mesh_.vertex(_vh);
  //TODO: make sure the initial start vertex generates positive tet volume
  auto normal = vertex_normal_uniform_weighted(mesh_, _vh);

  auto pt_in = pt_old - 0.001 * target_length_[_vh] * normal;

  //check if original incident cells have positive volume with the new position
  bool has_invalid = false;
  auto cells_pts1 = vo_.get_vertex_cells_points(_vh);


  auto cells_pts = cells_pts1;
  //check if new created cells have positive volume with the new position
  //get boundary halffaces
  std::set<HFH> hfhs;
  for (auto vhf_it = mesh_.vhf_iter(_vh); vhf_it.valid(); ++vhf_it)
  {
    if (mesh_.is_boundary(*vhf_it))
    {
      hfhs.insert(*vhf_it);

      auto hfh_vhs = mesh_.get_halfface_vertices(*vhf_it);
      std::vector<Point> cell_points;
      cell_points.reserve(4);

      cell_points.push_back(pt_in);
      for (const auto &vh: hfh_vhs)
        cell_points.push_back(mesh_.vertex(vh));

      cells_pts.push_back(cell_points);
    }
  }

  for (auto &cpts: cells_pts)
  {
    cpts[0] = pt_in;

    if (!RemeshingAssist<MeshT>::is_tet_valid(cpts))
    {
      ALGOHEX_DEBUG_ONLY(
              std::cerr << "Couldn't find a valid start with normal estimation. New tets are degenerate!" << _vh
                        << " is sgl? " << sgl_vt_[_vh] << std::endl;)
      int n_bdy = 0, n_val1 = 0, n_int = 0, n_feature = 0;
      for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
      {
        if (valence_[*ve_it] != 0)
        {
          if (feature_edge_[*ve_it] > 0)
            n_feature++;

          if (mesh_.is_boundary(*ve_it))
            n_bdy++;
          else
            n_int++;
        }
      }
      ALGOHEX_DEBUG_ONLY(
              std::cerr << "n_bdy sg: " << n_bdy << " n_int sg: " << n_int << " n_feature " << n_feature << std::endl;)

      has_invalid = true;
      break;
    }
  }

  //sample points on opposite faces
  if (has_invalid)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Sample points on the opposite triangles in original tets..." << std::endl;)
    std::vector<Point> s_points;
    for (auto &c_pts: cells_pts1)
    {
      RemeshingAssist<MeshT>::sample_points(c_pts[1], c_pts[2], c_pts[3], 0, 2, s_points);
    }

    for (const auto &pt_s_i: s_points)
    {
      auto dir = (pt_s_i - pt_old).normalize();
      auto pt_in2 = pt_old + 0.001 * target_length_[_vh] * dir;

      bool valid_pt = true;
      for (auto &cpts: cells_pts)
      {
        cpts[0] = pt_in2;

        if (!RemeshingAssist<MeshT>::is_tet_valid(cpts))
        {
          valid_pt = false;
          break;
        }
      }

      if (valid_pt)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "Found a valid start point " << pt_old << " " << pt_in2 << std::endl;)

        has_invalid = false;
        pt_in = pt_in2;
        break;
      }
    }
  }

  if (has_invalid)
  {
    std::cerr << "Warning: sample points failed, still couldn't find a valid start point..." << std::endl;
    return VH(-1);
  }

  //add new vertex
  auto vh_new = mesh_.add_vertex(pt_old);

  target_length_[vh_new] = target_length_[_vh];

  //add new cells
  for (const auto hfh: hfhs)
  {
    auto hfvhs = mesh_.get_halfface_vertices(hfh);
    hfvhs.push_back(vh_new);
    auto ch_new = mesh_.add_cell(hfvhs);

    //update quaternions as initial starting point in optimization later
    cell_quaternions_[ch_new] = cell_quaternions_[mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh))];
  }

  //set pt_in to _vh
  mesh_.set_vertex(_vh, pt_in);

  //mark it as non ffv
  feature_face_vertex_[_vh] = false;

  //optimize vertex position;
  vo_.vertex_optimize(_vh);


  return vh_new;
}


template<class MeshT>
void
PushSingularVertexT<MeshT>::determine_matchings(const VH _vh_new, const VH _vh_orig,
                                                std::map<HEH, int> &_e_to_trans_idx,
                                                std::map<HEH, int> &_hfh_to_trans_idx)
{
  //find boundary faces
  for (auto vhf_it = mesh_.vhf_iter(_vh_new); vhf_it.valid(); ++vhf_it)
    if (mesh_.is_boundary(*vhf_it))
    {
      trans_prop_[*vhf_it] = -1;
      trans_prop_[mesh_.opposite_halfface_handle(*vhf_it)] = -1;
    }


  auto heh = mesh_.find_halfedge(_vh_orig, _vh_new);
  for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
  {
    auto heh_prev = mesh_.prev_halfedge_in_halfface(heh, *hehf_it);

    //set matching of new faces
    int trans = _hfh_to_trans_idx[heh_prev];
    if (_hfh_to_trans_idx.find(heh_prev) == _hfh_to_trans_idx.end())
      std::cerr << "Error: could not find halfedge trans map" << std::endl;
    trans_prop_[*hehf_it] = trans;
    trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = tq_.inverse_transition_idx(trans);

    //find bottom face
    HFH hfh_btm = mesh_.adjacent_halfface_in_cell(*hehf_it, heh_prev);
    FH fh_btm = mesh_.face_handle(hfh_btm);

    //set matching
    trans_prop_[hfh_btm] = 0;
    trans_prop_[mesh_.opposite_halfface_handle(hfh_btm)] = 0;


    //update feature face properties
    auto heh_nxt = mesh_.next_halfedge_in_halfface(heh, *hehf_it);
    auto fh_top = mesh_.face_handle(mesh_.adjacent_halfface_in_cell(*hehf_it, heh_nxt));
    feature_fprop_[fh_top] = feature_fprop_[fh_btm];
    feature_fprop_[fh_btm] = 0;

    for (auto fe_it = mesh_.fe_iter(fh_btm); fe_it.valid(); ++fe_it)
      feature_face_edge_[*fe_it] = false;
    for (auto fv_it = mesh_.fv_iter(fh_btm); fv_it.valid(); ++fv_it)
      feature_face_vertex_[*fv_it] = false;

    if (feature_fprop_[fh_top] > 0)
    {
      for (auto fe_it = mesh_.fe_iter(fh_top); fe_it.valid(); ++fe_it)
        feature_face_edge_[*fe_it] = true;

      for (auto fv_it = mesh_.fv_iter(fh_top); fv_it.valid(); ++fv_it)
        feature_face_vertex_[*fv_it] = true;
    }
  }
}


}
/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define ARCZIPPINGT_C

#include "ArcZippingT.hh"

namespace AlgoHex
{
template<class MeshT>
std::set<FH> ArcZippingT<MeshT>::generate_cut_faces(const VH _vh, const std::set<CH> _onering_chs,
                                                    const std::vector<HEH> &_bdy_excluded_hehs,
                                                    const std::vector<HEH> &_interior_excluded_hehs,
                                                    const std::set<EH> &_bound_ehs) const
{
  //
  std::map<CH, int> visited_chs;
  std::vector<CH> visited_chs_other;
  std::set<FH> visited_fhs;


  //visit boundary excluded cells. same component
  int id = 0;
  for (const auto &heh: _bdy_excluded_hehs)
  {
    for (auto hec_it = mesh_.hec_iter(heh); hec_it.valid(); ++hec_it)
    {
      visited_chs[*hec_it] = id;
    }
  }

  for (const auto &heh: _bdy_excluded_hehs)
  {
    for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
      visited_fhs.insert(mesh_.face_handle(*hehf_it));
  }

  //visit interior excluded cells. different component
  for (const auto &heh: _interior_excluded_hehs)
  {
    id++;
    for (auto hec_it = mesh_.hec_iter(heh); hec_it.valid(); ++hec_it)
    {
      visited_chs[*hec_it] = id;
    }
  }

  //faces
  for (const auto &heh: _interior_excluded_hehs)
  {
    for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
      visited_fhs.insert(mesh_.face_handle(*hehf_it));
  }


  std::queue<HFH> que;

  if (visited_chs.empty())
  {
    auto ch0 = *mesh_.vc_iter(_vh);
    visited_chs[ch0] = 0;
    for (auto i = 0u; i < mesh_.cell(ch0).halffaces().size(); ++i)
      que.push(mesh_.cell(ch0).halffaces()[i]);
  }
  else
  {
    for (auto&[chi, idx]: visited_chs)
      for (auto i = 0u; i < mesh_.cell(chi).halffaces().size(); ++i)
        que.push(mesh_.cell(chi).halffaces()[i]);
  }


  while (!que.empty())
  {
    auto hf_cur = que.front();
    que.pop();

    auto hf_opp = mesh_.opposite_halfface_handle(hf_cur);
    CH ch_next = mesh_.incident_cell(hf_opp);
    // only process if not on boundary
    if (ch_next.is_valid())
    {
      //check if cell is not visited, or is in the one ring cells
      if (_onering_chs.find(ch_next) != _onering_chs.end() && visited_chs.find(ch_next) == visited_chs.end())
      {
        // remove face from cut graph
        visited_fhs.insert(mesh_.face_handle(hf_cur));

        //process cell
        CH ch_cur = mesh_.incident_cell(hf_cur);
        visited_chs[ch_next] = visited_chs[ch_cur];

        //push next
        for (unsigned int i = 0; i < mesh_.cell(ch_next).halffaces().size(); ++i)
          que.push(mesh_.cell(ch_next).halffaces()[i]);
      }
    }
  }

  //connect components if there are
  if (!_interior_excluded_hehs.empty())
  {
    std::map<int, int> neu_idx;
    for (int i = 0; i <= id; ++i)
      neu_idx[i] = i;

    for (auto vf_it = mesh_.vf_iter(_vh); vf_it.valid(); ++vf_it)
    {
      CH ch0 = mesh_.incident_cell(mesh_.halfface_handle(*vf_it, 0));
      CH ch1 = mesh_.incident_cell(mesh_.halfface_handle(*vf_it, 1));

      int id0 = neu_idx[visited_chs[ch0]];
      int id1 = neu_idx[visited_chs[ch1]];

      if (id0 < id1)
      {
        neu_idx[visited_chs[ch1]] = id0;
        visited_fhs.insert(*vf_it);
      }
      else if (id0 > id1)
      {
        neu_idx[visited_chs[ch0]] = id1;
        visited_fhs.insert(*vf_it);
      }
    }

    for (auto&[i, id]: neu_idx)
      if (id != 0)
        std::cerr << "components not connected" << std::endl;
  }


  //shrink cut surface
  std::set<EH> or_ehs;
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
    or_ehs.insert(*ve_it);

  std::queue<EH> qeh;
  FH cut_fh;
  //initialize queue with edges (1) not on boundary, (2) not singularity and (3) adjacent to one cut face
  for (const auto eh: or_ehs)
    if (!mesh_.is_boundary(eh) && _bound_ehs.find(eh) == _bound_ehs.end() &&
        n_cut_faces(mesh_, eh, visited_fhs, cut_fh) == 1)
      qeh.push(eh);

  // shrink while possible
  while (!qeh.empty())
  {
    // get next element
    EH eh = qeh.front();
    qeh.pop();

    // still a valid candidate?
    if (n_cut_faces(mesh_, eh, visited_fhs, cut_fh) == 1)
    {
      // remove from cut
      visited_fhs.insert(cut_fh);
      // push new candidates
      auto hehs = mesh_.face(cut_fh).halfedges();
      for (auto i = 0u; i < hehs.size(); ++i)
      {
        eh = mesh_.edge_handle(hehs[i]);
        if (!mesh_.is_boundary(eh) && or_ehs.find(eh) != or_ehs.end() && _bound_ehs.find(eh) == _bound_ehs.end())
          qeh.push(eh);
      }
    }
  }

  //get cut faces
  std::set<FH> cut_fhs;
  for (auto vf_it = mesh_.vf_iter(_vh); vf_it.valid(); ++vf_it)
    if (!mesh_.is_boundary(*vf_it) && visited_fhs.find(*vf_it) == visited_fhs.end())
      cut_fhs.insert(*vf_it);

  return cut_fhs;
}

template<class MeshT>
void ArcZippingT<MeshT>::coordinate_transformation_in_cell(const HFH hfh)
{
  CH ch = mesh_.incident_cell(hfh);

  if (!ch.is_valid())
    return;

  // get matching
  int trans = trans_prop_[hfh];
  int trans_opp = trans_prop_[mesh_.opposite_halfface_handle(hfh)];


  for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
  {
    if (mesh_.is_boundary(mesh_.face_handle(*chf_it)))
      continue;

    trans_prop_[*chf_it] = tq_.mult_transitions_idx(trans_opp, trans_prop_[*chf_it]);
    trans_prop_[mesh_.opposite_halfface_handle(*chf_it)] = tq_.inverse_transition_idx(trans_prop_[*chf_it]);
  }

  // transform quaternion
  // q_t = q_ch^-1 * q_chopp
  // q_ch^-1 = q_t*q_chopp^-1
  // q_t^-1*q_t = Id
  // multiply q_t^-1 on both sides
  // q_ch -> q_ch*q_t
  cell_quaternions_[ch] = cell_quaternions_[ch] * tq_.transition(trans);
}


template<class MeshT>
void ArcZippingT<MeshT>::regular_edge_closing(const HEH _he, const HFH _hf)
{
  auto hf_it = mesh_.adjacent_halfface_in_cell(_hf, _he);
  hf_it = mesh_.opposite_halfface_handle(hf_it);
  int idx = 0;
  while (hf_it != _hf)
  {
    idx = tq_.mult_transitions_idx(trans_prop_[hf_it], idx);
    //ccw
    hf_it = mesh_.adjacent_halfface_in_cell(hf_it, _he);
    hf_it = mesh_.opposite_halfface_handle(hf_it);
  }

  trans_prop_[_hf] = tq_.inverse_transition_idx(idx);
  trans_prop_[mesh_.opposite_halfface_handle(_hf)] = idx;
}


template<class MeshT>
bool ArcZippingT<MeshT>::zip_edges_onering(const VH _vh, const EH start_eh, const std::set<EH> &_bound_ehs,
                                           const std::set<FH> &_cut_fhs)
{
  HFH hf_s;
  auto sg_heh = mesh_.halfedge_handle(start_eh, 0);
  for (auto hehf_it = mesh_.hehf_iter(sg_heh); hehf_it.valid(); ++hehf_it)
    if (_cut_fhs.find(mesh_.face_handle(*hehf_it)) != _cut_fhs.end())
      hf_s = *hehf_it;

  if (!hf_s.is_valid())
  {
    std::cout << "Warning: could not find cut face adjacent to halfedge. Edge: " << start_eh << std::endl;
    return false;
  }
  int trans_orig = trans_prop_[hf_s];


  //get original transition of the singular edge start from hf_s
  int e_trans = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, sg_heh, hf_s);


  if (!mesh_.is_boundary(start_eh))
  {
    //reset matchings on the cutfaces
    for (const auto fh: _cut_fhs)
    {
      trans_prop_[mesh_.halfface_handle(fh, 0)] = -1;
      trans_prop_[mesh_.halfface_handle(fh, 1)] = -1;
    }

    trans_prop_[hf_s] = tq_.mult_transitions_idx(trans_orig, tq_.inverse_transition_idx(e_trans));
    trans_prop_[mesh_.opposite_halfface_handle(hf_s)] = tq_.inverse_transition_idx(trans_prop_[hf_s]);
  }
  else
  {
    //tb ta = tf   tn tb ta = treal  tn tx ta = 0  tx = tn-1 ta-1  tn = treal ta-1 tb-1  tx = tb ta treal-1 ta-1
    // tr tu = tx  tr to = tb  tr = tb to-1  tu = tr-1 tx = to tb-1 tx = to ta treal-1 ta-1
    auto hf_bdy = *mesh_.hehf_iter(sg_heh);
    int full_e_trans = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, sg_heh, hf_bdy);
    int full_trans =
            EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                               sg_heh, hf_bdy) + 4;

    int t_a = tq_.mult_transitions_idx(tq_.inverse_transition_idx(e_trans), full_e_trans);
    int t_u = tq_.mult_transitions_idx(tq_.inverse_transition_idx(full_trans), tq_.inverse_transition_idx(t_a));
    t_u = tq_.mult_transitions_idx(t_a, t_u);
    t_u = tq_.mult_transitions_idx(trans_orig, t_u);

    //reset matchings on the cutfaces
    for (const auto fh: _cut_fhs)
    {
      trans_prop_[mesh_.halfface_handle(fh, 0)] = -1;
      trans_prop_[mesh_.halfface_handle(fh, 1)] = -1;
    }

    trans_prop_[hf_s] = t_u;
    trans_prop_[mesh_.opposite_halfface_handle(hf_s)] = tq_.inverse_transition_idx(t_u);
  }

  auto fh_s = mesh_.face_handle(hf_s);

  std::set<EH> patch_ehs, interior_patch_ehs;
  for (const auto fhi: _cut_fhs)
  {
    for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
      patch_ehs.insert(*fe_it);
  }
  for (const auto ehi: patch_ehs)
  {
    if (!mesh_.is_boundary(ehi) && _bound_ehs.find(ehi) == _bound_ehs.end())
      interior_patch_ehs.insert(ehi);
  }

  std::queue<EH> que;
  for (auto fe_it = mesh_.fe_iter(fh_s); fe_it.valid(); ++fe_it)
    if (interior_patch_ehs.find(*fe_it) != interior_patch_ehs.end() && valence_[*fe_it] == 0)
      que.push(*fe_it);

  while (!que.empty())
  {
    auto eh = que.front();
    que.pop();

    if (mesh_.is_boundary(eh))
      continue;

    auto he = mesh_.halfedge_handle(eh, 0);
    HFH hf;
    if (n_cut_faces(he, hf) == 1)
    {
      regular_edge_closing(he, hf);

      for (auto hfe_it = mesh_.hfe_iter(hf); hfe_it.valid(); ++hfe_it)
        if (interior_patch_ehs.find(*hfe_it) != interior_patch_ehs.end() && valence_[*hfe_it] == 0)
          que.push(*hfe_it);
    }
  }

  return true;
}


template<class MeshT>
std::set<EH>
ArcZippingT<MeshT>::zipping_onering(const VH _vh, const std::set<FH> &cut_fhs, const std::set<EH> &_bound_ehs)
{
  if (cut_fhs.empty())
    return std::set<EH>{};

  //try to find an admissable singular edge
  EH eh_open(-1);
  for (const auto ehi: _bound_ehs)
    if ((valence_[ehi] == -1 || valence_[ehi] == 1) && !mesh_.is_boundary(ehi))
      eh_open = ehi;

  if (!eh_open.is_valid() && !mesh_.is_boundary(*_bound_ehs.begin()))
    eh_open = *_bound_ehs.begin();

  //boundary singular edge
  if (!eh_open.is_valid())
  {
    for (const auto ehi: _bound_ehs)
      if ((valence_[ehi] == -1 || valence_[ehi] == 1))
        eh_open = ehi;
  }

  if (!eh_open.is_valid())
    return std::set<EH>{};

  bool suc = zip_edges_onering(_vh, eh_open, _bound_ehs, cut_fhs);

  if (!suc)
    return std::set<EH>{};

  //check
  for (const auto fhi: cut_fhs)
  {
    if (mesh_.is_boundary(fhi))
      continue;
    auto hfh0 = mesh_.halfface_handle(fhi, 0);
    if (trans_prop_[hfh0] == -1)
    {
      std::cerr << "Error: matching of face " << fhi << " is not solved!" << std::endl;

      return std::set<EH>{};
    }
  }

  //for updating edge valence
  std::set<EH> cf_ehs;
  for (const auto &cf: cut_fhs)
  {
    for (auto fe_it = mesh_.fe_iter(cf); fe_it.valid(); ++fe_it)
      cf_ehs.insert(*fe_it);
  }

  return cf_ehs;
}

template<class MeshT>
bool
ArcZippingT<MeshT>::zipping_surface(const std::set<HFH> &_surface, const std::set<EH> &_edges, const HEH _sg_heh,
                                    const int _new_he_trans)
{
  HFH hf_s;
  for (auto hehf_it = mesh_.hehf_iter(_sg_heh); hehf_it.valid(); ++hehf_it)
    if (_surface.find(*hehf_it) != _surface.end() ||
        _surface.find(mesh_.opposite_halfface_handle(*hehf_it)) != _surface.end())
    {
      hf_s = *hehf_it;
      break;
    }

  if (!hf_s.is_valid())
  {
    std::cout << "Error: couldn't find incident halfface of edge " << mesh_.edge_handle(_sg_heh) << " on surface!"
              << std::endl;
    return false;
  }

  int orig_hf_s_trans = trans_prop_[hf_s];
  int e_trans = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, _sg_heh, hf_s);


  //store old
  std::map<HFH, int> old_hf_trans;
  for (const auto hfh: _surface)
  {
    if (trans_prop_[hfh] != -1)
    {
      old_hf_trans.insert(std::make_pair(hfh, trans_prop_[hfh]));
      HFH hf_opp = mesh_.opposite_halfface_handle(hfh);
      old_hf_trans.insert(std::make_pair(hf_opp, trans_prop_[hf_opp]));
    }
  }

  //t_e t_x = t_new_e
  //t_x = t_e^-1 t_new_e
  if (!mesh_.is_boundary(_sg_heh))
  {
    //reset transition on curfaces
    for (const auto hfh: _surface)
    {
      trans_prop_[hfh] = -1;
      trans_prop_[mesh_.opposite_halfface_handle(hfh)] = -1;
    }

    int t_x = tq_.mult_transitions_idx(tq_.inverse_transition_idx(e_trans), _new_he_trans);
    trans_prop_[hf_s] = tq_.mult_transitions_idx(orig_hf_s_trans, t_x);
    trans_prop_[mesh_.opposite_halfface_handle(hf_s)] = tq_.inverse_transition_idx(trans_prop_[hf_s]);
  }
  else
  {
    auto hf_bdy = *mesh_.hehf_iter(_sg_heh);
    int full_e_trans = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, _sg_heh, hf_bdy);

    int t_r = tq_.mult_transitions_idx(tq_.inverse_transition_idx(full_e_trans), e_trans);
    t_r = tq_.mult_transitions_idx(_new_he_trans, t_r);
    t_r = tq_.mult_transitions_idx(tq_.inverse_transition_idx(e_trans), t_r);

    //reset transition on curfaces
    for (const auto hfh: _surface)
    {
      trans_prop_[hfh] = -1;
      trans_prop_[mesh_.opposite_halfface_handle(hfh)] = -1;
    }

    trans_prop_[hf_s] = tq_.mult_transitions_idx(orig_hf_s_trans, t_r);
    trans_prop_[mesh_.opposite_halfface_handle(hf_s)] = tq_.inverse_transition_idx(trans_prop_[hf_s]);
  }


  std::queue<EH> que;
  for (auto hfe_it = mesh_.hfe_iter(hf_s); hfe_it.valid(); ++hfe_it)
    que.push(*hfe_it);

  while (!que.empty())
  {
    auto e_cur = que.front();
    que.pop();

    if (mesh_.is_boundary(e_cur) || _edges.find(e_cur) != _edges.end())
      continue;

    HFH hf_uk;
    auto he_cur = mesh_.halfedge_handle(e_cur, 0);
    int n_cf = n_cut_faces(he_cur, hf_uk);

    if (n_cf == 1)
    {
      regular_edge_closing(he_cur, hf_uk);

      for (auto hfe_it = mesh_.hfe_iter(hf_uk); hfe_it.valid(); ++hfe_it)
        que.push(*hfe_it);
    }
  }


  //check transitions on surface
  bool suc = true;
  for (auto&[hfh, trans]: old_hf_trans)
  {
    if (trans_prop_[hfh] == -1)
    {
      suc = false;
      std::cerr << "Warning: unsolved face." << std::endl;
      break;
    }
  }

  //restore
  if (!suc)
  {
    for (auto&[hfh, trans]: old_hf_trans)
      trans_prop_[hfh] = trans;
  }

  return suc;
}


template<class MeshT>
int ArcZippingT<MeshT>::n_cut_faces(const EH _eh, const std::set<FH> &_visited_fhs, FH &_fh) const
{
  int n = 0;
  auto heh = mesh_.halfedge_handle(_eh, 0);
  for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
  {
    auto fh = mesh_.face_handle(*hehf_it);
    if (_visited_fhs.find(fh) == _visited_fhs.end())
    {
      _fh = mesh_.face_handle(*hehf_it);
      n++;
    }
  }

  return n;
}

template<class MeshT>
int ArcZippingT<MeshT>::n_cut_faces(const HEH _heh, HFH &_hfh) const
{
  int n = 0;
  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    if (trans_prop_[*hehf_it] == -1)
    {
      _hfh = *hehf_it;
      n++;
    }
  }

  return n;
}

template<class MeshT>
std::vector<VH>
ArcZippingT<MeshT>::
shortest_path_vertex_to_vertex_on_boundary(const VH &_vh_s, const VH &_vh_t) const
{
  std::vector<VH> path;

  //distance to _vh0
  std::vector<double> dist(mesh_.n_vertices(), std::numeric_limits<double>::max());
  //previous vertex in the path
  std::vector<VH> pre_vh(mesh_.n_vertices(), VH(-1));

  //visited vertices
  std::set<VH> visited_vhs;
  visited_vhs.insert(_vh_s);

  dist[_vh_s.idx()] = 0.;
  pre_vh[_vh_s.idx()] = VH(-1);

  std::queue<VH> que;
  que.push(_vh_s);

  while (!que.empty())
  {
    VH vh_cur = que.front();
    que.pop();

    if (vh_cur == _vh_t)
      break;

    for (auto voh_it = mesh_.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
    {
      VH tvh = mesh_.halfedge(*voh_it).to_vertex();
      if (!mesh_.is_boundary(tvh))
        continue;

      double distance = dist[vh_cur.idx()] + mesh_.length(*voh_it);

      //if new distance is shorter
      if (distance < dist[tvh.idx()])
      {
        dist[tvh.idx()] = distance;
        pre_vh[tvh.idx()] = vh_cur;

        visited_vhs.insert(tvh);

        que.push(tvh);
      }
    }
  }

  if (dist[_vh_t.idx()] < std::numeric_limits<int>::max())
  {
    path.push_back(_vh_t);

    //find the path
    std::queue<VH> qp;
    qp.push(_vh_t);
    while (!qp.empty())
    {
      VH vh = qp.front();
      qp.pop();

      VH vh_pre = pre_vh[vh.idx()];
      if (vh_pre != VH(-1))
      {
        path.push_back(vh_pre);

        qp.push(vh_pre);
      }
    }
  }

  std::reverse(path.begin(), path.end());

  return path;
}


template<class MeshT>
std::set<FH>
ArcZippingT<MeshT>::
generate_cut_surface_with_bounding_edges(const std::set<EH> &_bdy_ehs, const EH &_sg_eh0, const EH &_sg_eh1,
                                         const VH &_vh_c) const
{
  //vertices on the boundary path
  std::set<VH> bdy_vhs;
  for (const auto ehi: _bdy_ehs)
  {
    bdy_vhs.insert(mesh_.edge(ehi).from_vertex());
    bdy_vhs.insert(mesh_.edge(ehi).to_vertex());
  }

  //grow cell cluster from the center vertex until all incident cells of bdy_ehs are in the cluster
  //incident cells
  std::set<CH> target_chs;
  for (const auto ehi: _bdy_ehs)
    for (auto ec_it = mesh_.ec_iter(ehi); ec_it.valid(); ++ec_it)
      target_chs.insert(*ec_it);


  std::queue<CH> que;

  for (auto vc_it = mesh_.vc_iter(_vh_c); vc_it.valid(); ++vc_it)
    que.push(*vc_it);

//        std::set<VH> visited_vhs;
  std::set<CH> visited_chs;

  while (!que.empty())
  {
    auto ch_cur = que.front();
    que.pop();

    if (visited_chs.find(ch_cur) != visited_chs.end())
      continue;

    visited_chs.insert(ch_cur);

    bool all_found = std::includes(visited_chs.begin(), visited_chs.end(), target_chs.begin(), target_chs.end());

    if (all_found)
      break;

    auto hfhs = mesh_.cell(ch_cur).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hfh_opp = mesh_.opposite_halfface_handle(hfhs[i]);
      if (mesh_.is_boundary(hfh_opp))
        continue;

      auto ch_opp = mesh_.incident_cell(hfh_opp);
      que.push(ch_opp);
    }
  }

  //find all boundary edges except for the path
  std::set<EH> other_ball_bdy_ehs;
  std::set<FH> ball_bdy_fhs, sge0_fhs, sge1_fhs;
  for (const auto chi: visited_chs)
  {
    for (auto chf_it = mesh_.chf_iter(chi); chf_it.valid(); ++chf_it)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(*chf_it);
      auto ch_opp = mesh_.incident_cell(hf_opp);

      if (!ch_opp.is_valid() || visited_chs.find(ch_opp) == visited_chs.end())
      {
        ball_bdy_fhs.insert(mesh_.face_handle(*chf_it));
      }
    }
  }

  for (auto ef_it = mesh_.ef_iter(_sg_eh0); ef_it.valid(); ++ef_it)
  {
    sge0_fhs.insert(*ef_it);
  }
  for (auto ef_it = mesh_.ef_iter(_sg_eh1); ef_it.valid(); ++ef_it)
  {
    sge1_fhs.insert(*ef_it);
  }

  for (const auto fhi: ball_bdy_fhs)
  {
    for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
    {
      if (_bdy_ehs.find(*fe_it) == _bdy_ehs.end())
        other_ball_bdy_ehs.insert(*fe_it);
    }
  }

  //connect all incident cells
  std::set<CH> visited_chs2;
  for (const auto &ehi: other_ball_bdy_ehs)
  {
    for (auto ec_it = mesh_.ec_iter(ehi); ec_it.valid(); ++ec_it)
      if (visited_chs.find(*ec_it) != visited_chs.end())
        visited_chs2.insert(*ec_it);
  }

  //faces
  std::set<FH> visited_fhs;
  for (const auto &ehi: other_ball_bdy_ehs)
  {
    for (auto ef_it = mesh_.ef_iter(ehi); ef_it.valid(); ++ef_it)
      visited_fhs.insert(*ef_it);
  }

  //construct spanning tree
  std::queue<HFH> hf_que;

  //push to queue
  for (const auto &chi: visited_chs2)
    for (auto i = 0u; i < mesh_.cell(chi).halffaces().size(); ++i)
      hf_que.push(mesh_.cell(chi).halffaces()[i]);

  while (!hf_que.empty())
  {
    auto hf_cur = mesh_.opposite_halfface_handle(hf_que.front());
    hf_que.pop();

    // only process if not on boundary
    if (!mesh_.is_boundary(hf_cur))
    {
      //check if cell is not visited, or is in the one ring cells
      CH ch = mesh_.incident_cell(hf_cur);
      if (visited_chs.find(ch) != visited_chs.end() && visited_chs2.find(ch) == visited_chs2.end())
      {
        // remove face from cut graph
        visited_fhs.insert(mesh_.face_handle(hf_cur));

        //process cell
        visited_chs2.insert(ch);
        for (unsigned int i = 0; i < mesh_.cell(ch).halffaces().size(); ++i)
          hf_que.push(mesh_.cell(ch).halffaces()[i]);
      }
    }
  }

  //generate cut surface that crosses the path
  //cluster boundary edges
  other_ball_bdy_ehs.insert(_bdy_ehs.begin(), _bdy_ehs.end());

  //shrink cut surface
  std::queue<EH> qeh;
  FH cut_fh;
  //initialize queue with edges (1) not on boundary, (2) not singularity and (3) adjacent to one cut face
  for (const auto chi: visited_chs)
  {
    for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
      if (!mesh_.is_boundary(*ce_it) && valence_[*ce_it] == 0 &&
          n_cut_faces(*ce_it, visited_fhs, cut_fh) == 1)
        qeh.push(*ce_it);
  }

  // shrink while possible
  while (!qeh.empty())
  {
    // get next element
    EH eh = qeh.front();
    qeh.pop();

    if (valence_[eh] != 0 || mesh_.is_boundary(eh) || other_ball_bdy_ehs.find(eh) != other_ball_bdy_ehs.end())
      continue;

    // still a valid candidate?
    if (n_cut_faces(eh, visited_fhs, cut_fh) == 1)
    {
      // remove from cut
      visited_fhs.insert(cut_fh);
      // push new candidates
      auto hehs = mesh_.face(cut_fh).halfedges();
      for (auto i = 0u; i < hehs.size(); ++i)
      {
        eh = mesh_.edge_handle(hehs[i]);
        if (!mesh_.is_boundary(eh) && valence_[eh] == 0)
          qeh.push(eh);
      }
    }
  }

  //get cut faces
  std::set<FH> all_fhs;
  for (const auto chi: visited_chs)
    for (auto cf_it = mesh_.cf_iter(chi); cf_it.valid(); ++cf_it)
      all_fhs.insert(*cf_it);


  //not on ball boundary and not visited
  std::set<FH> cut_fhs;
  for (const auto fhi: all_fhs)
    if (ball_bdy_fhs.find(fhi) == ball_bdy_fhs.end() && visited_fhs.find(fhi) == visited_fhs.end())
      cut_fhs.insert(fhi);

  return cut_fhs;
}
}

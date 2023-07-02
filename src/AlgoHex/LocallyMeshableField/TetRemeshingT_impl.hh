/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define TETREMESHINGT_C

#include "TetRemeshingT.hh"

namespace AlgoHex
{
template<class MeshT>
VH
CellSplitT<MeshT>::cell_split(const CH _ch)
{
  if (!is_split_ok(_ch))
    return VH(-1);

  Quaternion orig_qt = cell_quaternions_[_ch];

  auto cvhs = mesh_.get_cell_vertices(_ch);

  auto vh_new = mesh_.add_vertex(mesh_.barycenter(_ch));

  mesh_.delete_cell(_ch);

  VH vh_tmp = cvhs[0];
  cvhs[0] = vh_new;
  mesh_.add_cell(cvhs);
  cvhs[0] = vh_tmp;

  vh_tmp = cvhs[1];
  cvhs[1] = vh_new;
  mesh_.add_cell(cvhs);
  cvhs[1] = vh_tmp;

  vh_tmp = cvhs[2];
  cvhs[2] = vh_new;
  mesh_.add_cell(cvhs);
  cvhs[2] = vh_tmp;

  vh_tmp = cvhs[3];
  cvhs[3] = vh_new;
  mesh_.add_cell(cvhs);
  cvhs[3] = vh_tmp;

  update_properties(vh_new, orig_qt);

  return vh_new;
}

template<class MeshT>
bool
CellSplitT<MeshT>::is_split_ok(const CH _ch) const
{
  if (!_ch.is_valid() || mesh_.is_deleted(_ch))
    return false;

  return true;
}


template<class MeshT>
void
CellSplitT<MeshT>::update_properties(const VH _vh_new, const Quaternion &_old_qt)
{
  double tl = 0.;

  for (auto vc_it = mesh_.vc_iter(_vh_new); vc_it.valid(); ++vc_it)
    cell_quaternions_[*vc_it] = _old_qt;

  for (auto vv_it = mesh_.vv_iter(_vh_new); vv_it.valid(); ++vv_it)
    tl += target_length_[*vv_it];
  //
  target_length_[_vh_new] = tl / 4.;
}

//===============================================================================================================//

template<class MeshT>
VH
FaceSplitT<MeshT>::face_split(const FH _fh)
{
  //invalid face handle
  if (!is_split_ok(_fh))
    return VH(-1);


  auto hfh0 = mesh_.halfface_handle(_fh, 0);
  auto hfh_vhs = mesh_.get_halfface_vertices(hfh0);
  auto hfh1 = mesh_.halfface_handle(_fh, 1);

  //transition stored in halfface property
  int trans = trans_prop_[hfh0];
  int inv_trans = trans_prop_[hfh1];

  std::vector<Quaternion> orig_cell_qts;
  orig_cell_qts.reserve(2);

  //quaternions stored in cells property
  CH ch0 = mesh_.incident_cell(hfh0);
  if (ch0.is_valid())
    orig_cell_qts.push_back(cell_quaternions_[ch0]);
  else
    orig_cell_qts.push_back(Quaternion(1, 0, 0, 0));

  CH ch1 = mesh_.incident_cell(hfh1);
  if (ch1.is_valid())
    orig_cell_qts.push_back(cell_quaternions_[ch1]);
  else
    orig_cell_qts.push_back(Quaternion(1, 0, 0, 0));


  int ff_id = feature_fprop_[_fh];

  bool is_visited_face = visited_fprop_[_fh];

  auto vh_new = mesh_.split_face(_fh, mesh_.barycenter(_fh));


  update_properties(hfh_vhs, vh_new, trans, inv_trans, orig_cell_qts, ff_id, is_visited_face);

  return vh_new;
}

template<class MeshT>
bool
FaceSplitT<MeshT>::is_split_ok(const FH _fh) const
{
  if (!_fh.is_valid() || mesh_.is_deleted(_fh))
    return false;

  return true;
}


template<class MeshT>
void
FaceSplitT<MeshT>::update_properties(std::vector<VH> &_orig_face_vhs, const VH _vh_new, const int _trans,
                                     const int _inv_trans, const std::vector<Quaternion> &_old_qts, const int _ff_id,
                                     const bool _is_visited_face)
{
  double tl = 0.;

  for (int i = 0; i < 3; ++i)
  {
    //replace the vertex vha with the new vertex
    auto vha = _orig_face_vhs[i];
    _orig_face_vhs[i] = _vh_new;

    auto hfhi = mesh_.find_halfface(_orig_face_vhs);
    if (!hfhi.is_valid())
    {
      std::cerr << "Error: no halfface exist after face split!" << std::endl;
    }
    auto hfhi_opp = mesh_.opposite_halfface_handle(hfhi);

    //update transition property
    trans_prop_[hfhi] = _trans;
    trans_prop_[hfhi_opp] = _inv_trans;

    //update cell quaternion property
    CH ch0 = mesh_.incident_cell(hfhi);
    if (ch0.is_valid())
      cell_quaternions_[ch0] = _old_qts[0];

    CH ch1 = mesh_.incident_cell(hfhi_opp);
    if (ch1.is_valid())
      cell_quaternions_[ch1] = _old_qts[1];

    //restore original face vertex
    _orig_face_vhs[i] = vha;

    //target length
    tl += target_length_[vha];

    //feature face
    if (_ff_id > 0)
    {
      auto fhi = mesh_.face_handle(hfhi);
      feature_fprop_[fhi] = _ff_id;

      for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
        feature_face_edge_[*fe_it] = true;
      for (auto fv_it = mesh_.fv_iter(fhi); fv_it.valid(); ++fv_it)
        feature_face_vertex_[*fv_it] = true;
    }

    //visited face after fixing zipper node
    if (_is_visited_face)
    {
      auto fhi = mesh_.face_handle(hfhi);
      visited_fprop_[fhi] = _is_visited_face;
    }
  }

  //
  target_length_[_vh_new] = tl / 3.;
}

//===============================================================================================================//

template<class MeshT>
VH
EdgeSplitT<MeshT>::edge_split(const EH _eh, const bool _check_energy)
{
  std::vector<HFH> ops_hfhs;

  auto heh = mesh_.halfedge_handle(_eh, 0);
  auto vh_f = mesh_.halfedge(heh).from_vertex();
  auto vh_t = mesh_.halfedge(heh).to_vertex();

  //check energy
  Point split_pt = (mesh_.vertex(vh_f) + mesh_.vertex(vh_t)) / 2.;

  if (_check_energy)
  {
    bool better_energy = check_energy(heh, split_pt);
    if (!better_energy)
      return VH(-1);
  }

  //get reference objects for updating properties
  std::vector<std::vector<VH> > hfs_vhs;
  std::vector<int> hfs_trans, opp_hfs_trans;
  std::vector<int> feature_face_ids;
  std::vector<bool> is_visited_face;
  std::vector<Quaternion> orig_cell_qts;

  std::map<std::vector<VH>, int> alignments;
  for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
  {
    auto hfvhs = mesh_.get_halfface_vertices(*hehf_it);
    hfs_vhs.push_back(hfvhs);

    FH fhi = mesh_.face_handle(*hehf_it);
    //feature face property
    feature_face_ids.push_back(feature_fprop_[fhi]);

    //visited face property
    is_visited_face.push_back(visited_fprop_[fhi]);

    hfs_trans.push_back(trans_prop_[*hehf_it]);
    opp_hfs_trans.push_back(trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)]);

    if (!mesh_.is_boundary(*hehf_it))
    {
      orig_cell_qts.push_back(cell_quaternions_[mesh_.incident_cell(*hehf_it)]);
    }

    //for field alignment
    if (opt_split_)
    {
      if (feature_fprop_[fhi] > 0)
      {
        auto chi = mesh_.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          auto nm = mesh_.normal(*hehf_it);
          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[chi], nm).second;
          alignments.insert(std::make_pair(hfvhs, axis));
        }

        //opposite halfface
        auto hfh_opp = mesh_.opposite_halfface_handle(*hehf_it);
        auto chi_opp = mesh_.incident_cell(hfh_opp);
        if (chi_opp.is_valid())
        {
          auto hfvhs_opp = mesh_.get_halfface_vertices(hfh_opp);
          auto nm_opp = mesh_.normal(hfh_opp);
          int axis_opp = AxisAlignmentHelpers::closest_axis(cell_quaternions_[chi_opp], nm_opp).second;
          alignments.insert(std::make_pair(hfvhs_opp, axis_opp));
        }
      }
    }
  }


  int val = valence_[_eh];

  int fe_id = feature_edge_[_eh];

  int pid = sg_edge_pairs_[_eh];

  bool kept_on_ff = keep_on_feature_face_[_eh];

  VH vh = mesh_.split_edge(_eh);

  if (opt_split_)
    mesh_.set_vertex(vh, split_pt);

  update_properties(vh_f, vh_t, vh, val, pid, kept_on_ff, fe_id, feature_face_ids, is_visited_face, hfs_vhs, hfs_trans,
                    opp_hfs_trans, orig_cell_qts,
                    alignments);

  return vh;
}

template<class MeshT>
bool
EdgeSplitT<MeshT>::is_split_ok(const EH _eh) const
{
  if (!_eh.is_valid() || mesh_.is_deleted(_eh))
    return false;

  return true;
}


template<class MeshT>
void
EdgeSplitT<MeshT>::update_properties(const VH _vh_f, const VH _vh_t, const VH _vh_new, const int _valence,
                                     const int _pair_id,
                                     const bool _is_kept_on_ff, const int _fe_id,
                                     const std::vector<int> &_feature_face_ids,
                                     const std::vector<bool> &_is_visited_face,
                                     const std::vector<std::vector<VH> > &_hfs_vhs, const std::vector<int> &_hfs_trans,
                                     const std::vector<int> &_opp_hfs_trans, const std::vector<Quaternion> &_old_qts,
                                     std::map<std::vector<VH>, int> &_hfvs_axis)
{
  auto new_eh0 = mesh_.edge_handle(mesh_.find_halfedge(_vh_f, _vh_new));
  auto new_eh1 = mesh_.edge_handle(mesh_.find_halfedge(_vh_new, _vh_t));

  //update valence
  valence_[new_eh0] = _valence;
  valence_[new_eh1] = _valence;

  //update pair id
  sg_edge_pairs_[new_eh0] = _pair_id;
  sg_edge_pairs_[new_eh1] = _pair_id;

  //singular edge kept on feature face
  keep_on_feature_face_[new_eh0] = _is_kept_on_ff;
  keep_on_feature_face_[new_eh1] = _is_kept_on_ff;

  //update halfface transitions and cell matrices
  int j = 0;
  for (auto i = 0u; i < _hfs_vhs.size(); ++i)
  {
    auto hfs = get_new_halffaces(_hfs_vhs[i], _vh_f, _vh_t, _vh_new);
    bool is_bdy = mesh_.is_boundary(hfs[0]);
    for (const auto hfh: hfs)
    {
      trans_prop_[hfh] = _hfs_trans[i];
      trans_prop_[mesh_.opposite_halfface_handle(hfh)] = _opp_hfs_trans[i];

      if (!is_bdy)
      {
        auto chi = mesh_.incident_cell(hfh);
        cell_quaternions_[chi] = _old_qts[j];
      }

      //feature face property
      if (_feature_face_ids[i] > 0)
      {
        auto fhi = mesh_.face_handle(hfh);
        feature_fprop_[fhi] = _feature_face_ids[i];

        for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
          feature_face_edge_[*fe_it] = true;
        for (auto fv_it = mesh_.fv_iter(fhi); fv_it.valid(); ++fv_it)
          feature_face_vertex_[*fv_it] = true;
      }

      if (_is_visited_face[i])
      {
        auto fhi = mesh_.face_handle(hfh);
        visited_fprop_[fhi] = _is_visited_face[i];
      }
    }

    if (!is_bdy)
      j++;
  }

  //update singular vertex property
  if (_valence != 0)
  {
    if (mesh_.is_boundary(new_eh0))
      sgl_vt_[_vh_new] = 2;
    else
      sgl_vt_[_vh_new] = 1;
  }

  //update target length property
  target_length_[_vh_new] = (target_length_[_vh_f] + target_length_[_vh_t]) / 2.;

  //update feature edge property
  feature_edge_[new_eh0] = _fe_id;
  feature_edge_[new_eh1] = _fe_id;

  if (_fe_id > 0)
    feature_edge_vertex_[_vh_new] = true;

  //alignment
  if (opt_split_)
  {
    std::set<CH> cells;
    for (auto vc_it = mesh_.vc_iter(_vh_new); vc_it.valid(); ++vc_it)
      cells.insert(*vc_it);
    std::set<FH> bdy_fhs;
    for (const auto ch: cells)
      for (auto cf_it = mesh_.cf_iter(ch); cf_it.valid(); ++cf_it)
        if (mesh_.is_boundary(*cf_it))
          bdy_fhs.insert(*cf_it);

    std::map<CH, std::vector<std::pair<Vec3d, int>>> alignment;
    //feature face alignment
    std::map<HFH, int> hfh_to_axis;
    for (auto&[hfvhs, axis]: _hfvs_axis)
    {
      auto hfs = get_new_halffaces(hfvhs, _vh_f, _vh_t, _vh_new);

      for (const auto hfhi: hfs)
      {
        hfh_to_axis.insert(std::make_pair(hfhi, axis));
      }
    }

    for (auto&[hfh, axis]: hfh_to_axis)
    {
      auto ch_new = mesh_.incident_cell(hfh);
      if (!ch_new.is_valid())
        continue;

      auto nm = mesh_.normal(hfh);
      alignment[ch_new].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
    }

    //side faces
    std::set<FH> inc_fhs;
    for (auto vf_it = mesh_.vf_iter(_vh_new); vf_it.valid(); ++vf_it)
      inc_fhs.insert(*vf_it);

    for (const auto ch: cells)
      for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
      {
        auto cfh = mesh_.face_handle(*chf_it);
        if (feature_fprop_[cfh] > 0 && inc_fhs.find(cfh) == inc_fhs.end())
        {
          auto nm = mesh_.normal(*chf_it);
          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], nm).second;
          alignment[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
        }
      }

    //feature edge alignment
    for (const auto ch: cells)
    {
      for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
      {
        if (feature_edge_[*ce_it] > 0)
        {
          auto dir = mesh_.vertex(mesh_.edge(*ce_it).to_vertex()) - mesh_.vertex(mesh_.edge(*ce_it).from_vertex());
          dir.normalize();

          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], dir).second;
          alignment[ch].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
        }
      }
    }
//
//            //check
//            for(const auto ch : cells) {
//                if(alignment[ch].size() == 2) {
//                    if(alignment[ch][0].second == alignment[ch][1].second) {
//                        std::cerr<<"Error: cell "<<ch<<" has more than one constraints of the same axis in splittig edge at vh "<<_vh_new<<std::endl;
//                        for(auto [nm, ax] : alignment[ch]) {
//                            std::cerr<<" dir "<<nm[0]<<" "<<nm[1]<<" "<<nm[2]<<" ax "<<ax;
//                        }
//                        std::cerr<<std::endl;
//                    }
//                } else if(alignment[ch].size() > 2) {
//                    std::cerr<<"Error: cell "<<ch<<" has "<<alignment[ch].size()<<" constraints in splittig edge at vh "<<_vh_new<<std::endl;
//                    for(auto [nm, ax] : alignment[ch]) {
//                        std::cerr<<" dir "<<nm[0]<<" "<<nm[1]<<" "<<nm[2]<<" ax "<<ax;
//                    }
//                    std::cerr<<std::endl;
//                }
//            }

    for (int i = 0; i < 10; ++i)
      QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, cells, bdy_fhs, alignment, tq_,
                                                   cell_quaternions_);
  }
}

template<class MeshT>
std::vector<HFH>
EdgeSplitT<MeshT>::get_new_halffaces(const std::vector<VH> &_orig_vhs, const VH _vh_f, const VH _vh_t,
                                     const VH _vh_new) const
{
  std::vector<VH> new_hf0_vhs, new_hf1_vhs;
  new_hf0_vhs.reserve(3);
  new_hf1_vhs.reserve(3);
  for (auto vh: _orig_vhs)
  {
    if (vh == _vh_t)
      new_hf0_vhs.push_back(_vh_new);
    else
      new_hf0_vhs.push_back(vh);

    if (vh == _vh_f)
      new_hf1_vhs.push_back(_vh_new);
    else
      new_hf1_vhs.push_back(vh);
  }

  std::vector<HFH> hfs;
  hfs.reserve(2);
  hfs.push_back(mesh_.find_halfface(new_hf0_vhs));
  hfs.push_back(mesh_.find_halfface(new_hf1_vhs));

  return hfs;
}

template<class MeshT>
bool
EdgeSplitT<MeshT>::check_energy(const HEH _heh, Point &_split_pt) const
{
  std::vector<std::vector<Point> > new_cells_points, old_cells_points;
  for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
  {
    auto cvhs = mesh_.get_cell_vertices(*hec_it);
    std::vector<Point> cell_points;
    cell_points.reserve(4);
    for (const auto &vhi: cvhs)
      cell_points.push_back(mesh_.vertex(vhi));
    old_cells_points.push_back(cell_points);
  }

  new_cells_points.reserve(2 * old_cells_points.size());

  std::vector<HFH> new_cell_hfhs;
  new_cell_hfhs.reserve(new_cells_points.size());

  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    if (!mesh_.is_boundary(*hehf_it))
    {
      auto he_next = mesh_.next_halfedge_in_halfface(_heh, *hehf_it);
      new_cell_hfhs.push_back(mesh_.adjacent_halfface_in_cell(*hehf_it, he_next));

      auto he_prev = mesh_.prev_halfedge_in_halfface(_heh, *hehf_it);
      new_cell_hfhs.push_back(mesh_.adjacent_halfface_in_cell(*hehf_it, he_prev));
    }
  }

  for (const auto &hfhi: new_cell_hfhs)
  {
    auto hf_vhs = mesh_.get_halfface_vertices(mesh_.opposite_halfface_handle(hfhi));
    std::vector<Point> cell_points;
    cell_points.reserve(4);
    cell_points.push_back(_split_pt);
    for (const auto &vhi: hf_vhs)
      cell_points.push_back(mesh_.vertex(vhi));

    new_cells_points.push_back(cell_points);
  }

  if (opt_split_)
  {
    auto eh = mesh_.edge_handle(_heh);

    if (feature_edge_[eh] > 0 || valence_[eh] != 0)
    {
      Point pt_s = mesh_.vertex(mesh_.edge(eh).from_vertex());
      Point pt_e = mesh_.vertex(mesh_.edge(eh).to_vertex());

      std::vector<Point> s_points;
      for (int i = 0; i <= 10; ++i)
      {
        s_points.push_back((double) i / (double) 10 * pt_s + (double) (10 - i) / (double) 10 * pt_e);
      }

      _split_pt = vo_.optimize_vertex_on_edge_worst(s_points, new_cells_points);
    }
    else if (feature_face_edge_[eh])
    {
      std::vector<Point> s_points;
      s_points.push_back(_split_pt);
      for (auto ef_it = mesh_.ef_iter(eh); ef_it.valid(); ++ef_it)
      {
        if (feature_fprop_[*ef_it] > 0)
        {
          auto hf_vhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(*ef_it, 0));
          RemeshingAssist<MeshT>::sample_points(mesh_.vertex(hf_vhs[0]), mesh_.vertex(hf_vhs[1]),
                                                mesh_.vertex(hf_vhs[2]), 0, 2, s_points);
        }
      }

      _split_pt = vo_.optimize_feature_face_vertex_location_worst(s_points, new_cells_points);

    }
    else
    {
      _split_pt = vo_.optimize_interior_vertex_location_worst(new_cells_points);
    }

    for (auto &c_pts: new_cells_points)
      c_pts[0] = _split_pt;
  }


  if (!RemeshingAssist<MeshT>::is_quality_improved(new_cells_points, old_cells_points, true))
    return false;

  return true;
}

//===============================================================================================================//
template<class MeshT>
VH
EdgeCollapseT<MeshT>::edge_collapse(const HEH _heh, const bool _check_energy)
{
  VH vh_f = mesh_.halfedge(_heh).from_vertex();
  VH vh_t = mesh_.halfedge(_heh).to_vertex();

  std::vector<std::vector<VH> > v_vhs;
  surviving_tets(_heh, v_vhs);

  Point collapse_pt = mesh_.vertex(vh_t);
  bool better_energy = check_energy(_heh, v_vhs, collapse_pt, _check_energy);
  if (!better_energy)
    return VH(-1);


  //3. check if the collapse changes the boundary singular edge status
  if (mesh_.is_boundary(_heh))
  {
    if (is_boundary_edge_type_changed(_heh))
      return VH(-1);
  }

  //4. don't do it if one of the new regular edges connects to two feature edges or feature face
  bool is_vt_bdy = mesh_.is_boundary(vh_t);
  for (auto voh_it = mesh_.voh_iter(vh_f); voh_it.valid(); ++voh_it)
  {
    VH vhff = mesh_.halfedge(*voh_it).to_vertex();
    if (vhff == vh_t)
      continue;
    EH veh = mesh_.edge_handle(*voh_it);

    if (post_remesh_)
    {
      if (valence_[veh] == 0 && sgl_vt_[vh_t] != 0 && sgl_vt_[vhff] != 0)
        return VH(-1);
    }

    if (feature_edge_[veh] == 0 && feature_edge_vertex_[vh_t] && feature_edge_vertex_[vhff])
      return VH(-1);

    if (feature_face_vertex_[vh_t] && !feature_face_edge_[veh] && feature_face_vertex_[vhff])
      if (!mesh_.find_halfedge(vh_t, vhff).is_valid())
        return VH(-1);
  }

  //don't collapse if the operation creates a complex singular edge
  int val = valence_[mesh_.edge_handle(_heh)];
  if (!is_vt_bdy && (val == 1 || val == -1))
  {
    if (is_complex_singular_edge_after(_heh))
      return VH(-1);
  }

  //get reference for updating properties
  //for transitions
  std::vector<std::vector<VH> > hfs_vhs;
  std::vector<int> hfs_trans;
  std::vector<Quaternion> hfs_trans_q, opp_hfs_trans_q;

  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    if (is_vt_bdy && mesh_.is_boundary(*hehf_it))
      continue;
    CH cur_ch = mesh_.incident_cell(*hehf_it);
    auto hfh_top = opposite_halfface_to_vertex_in_cell(vh_f, cur_ch);
    auto hfh_bottom = opposite_halfface_to_vertex_in_cell(vh_t, cur_ch);

    auto vhs = mesh_.get_halfface_vertices(hfh_top);
    hfs_vhs.push_back(vhs);

    int trans = tq_.mult_transitions_idx(trans_prop_[mesh_.opposite_halfface_handle(hfh_top)],
                                         trans_prop_[hfh_bottom]);
    hfs_trans.push_back(tq_.inverse_transition_idx(trans));
  }

  //for updating edge valence
  std::vector<int> vals;
  std::vector<VH> f_vhs, k_vhs;
  auto eh = mesh_.edge_handle(_heh);

  //
  for (auto ec_it = mesh_.ec_iter(eh); ec_it.valid(); ++ec_it)
  {
    auto eh_opp = opposite_edge_in_cell(eh, *ec_it);
    if (valence_[eh_opp] > 1)
      return VH(-1);
  }

  if (valence_[eh] != 0)
  {
    for (auto voh_it = mesh_.voh_iter(vh_f); voh_it.valid(); ++voh_it)
    {
      auto ve = mesh_.edge_handle(*voh_it);
      int val = valence_[ve];
      if (val != 0)
      {
        auto vh_tt = mesh_.halfedge(*voh_it).to_vertex();
        f_vhs.push_back(vh_tt);
        vals.push_back(val);

        if (keep_on_feature_face_[ve])
          k_vhs.push_back(vh_tt);
      }
    }
  }

  //for updating feature edge property
  std::vector<VH> feature_vhs;
  std::vector<int> fe_ids;
  if (feature_edge_[eh] > 0)
  {
    for (auto voh_it = mesh_.voh_iter(vh_f); voh_it.valid(); ++voh_it)
    {
      if (feature_edge_[mesh_.edge_handle(*voh_it)] > 0)
      {
        feature_vhs.push_back(mesh_.halfedge(*voh_it).to_vertex());
        fe_ids.push_back(feature_edge_[mesh_.edge_handle(*voh_it)]);
      }
    }
  }

  //for updating feature face property
  std::vector<std::vector<VH>> ff_vhs;
  std::vector<int> ff_ids;
  if (feature_face_vertex_[vh_f] || feature_face_vertex_[vh_t])
  {
    std::set<FH> e_inc_ff;
    for (auto hef_it = mesh_.hef_iter(_heh); hef_it.valid(); ++hef_it)
      if (feature_fprop_[*hef_it] > 0)
        e_inc_ff.insert(*hef_it);
    std::set<FH> inc_ff;
    for (auto vf_it = mesh_.vf_iter(vh_f); vf_it.valid(); ++vf_it)
      if (feature_fprop_[*vf_it] > 0 && e_inc_ff.find(*vf_it) == e_inc_ff.end())
      {
        auto hfvhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it, 0), vh_f);
        ff_vhs.push_back(hfvhs);
        ff_ids.push_back(feature_fprop_[*vf_it]);
      }

    for (auto vf_it = mesh_.vf_iter(vh_t); vf_it.valid(); ++vf_it)
      if (feature_fprop_[*vf_it] > 0 && e_inc_ff.find(*vf_it) == e_inc_ff.end())
      {
        auto hfvhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it, 0), vh_t);
        ff_vhs.push_back(hfvhs);
        ff_ids.push_back(feature_fprop_[*vf_it]);
      }
  }

  //store the aligned normal axis of neighbouring feature faces
  std::map<std::vector<VH>, int> inc_hfv_to_axis, ninc_hfv_to_axis;
  collect_axis_aligned_to_feature_face(_heh, inc_hfv_to_axis, ninc_hfv_to_axis);


  VH vh_survive = mesh_.collapse_edge(_heh);

  if (opt_collapse_)
    mesh_.set_vertex(vh_survive, collapse_pt);

  update_properties(eh, vh_survive, hfs_vhs, hfs_trans, vals, f_vhs, k_vhs, feature_vhs, fe_ids, ff_vhs, ff_ids,
                    inc_hfv_to_axis, ninc_hfv_to_axis);

  return vh_survive;
}

template<class MeshT>
VH
EdgeCollapseT<MeshT>::edge_collapse_allow_inverted(const HEH _heh)
{
  VH vh_f = mesh_.halfedge(_heh).from_vertex();
  VH vh_t = mesh_.halfedge(_heh).to_vertex();
  auto eh = mesh_.edge_handle(_heh);

  if (feature_face_vertex_[vh_f] && !feature_face_vertex_[vh_t])
    return VH(-1);

  if (feature_edge_vertex_[vh_f] && !feature_edge_vertex_[vh_t])
    return VH(-1);

  if (feature_node_[vh_f])
    return VH(-1);

  std::vector<std::vector<VH> > v_vhs;
  surviving_tets(_heh, v_vhs);

  Point collapse_pt = mesh_.vertex(vh_t);

  //3. check if the collapse changes the boundary singular edge status
  if (mesh_.is_boundary(_heh))
  {
    if (is_boundary_edge_type_changed(_heh))
    {
      std::cerr << "bdy " << std::endl;
      return VH(-1);
    }
  }

  //4. don't do it if one of the new regular edges connects to two feature edges or feature face
  bool is_vt_bdy = mesh_.is_boundary(vh_t);

  //don't collapse if the operation creates a complex singular edge
  int val = valence_[mesh_.edge_handle(_heh)];
  if (!is_vt_bdy && (val == 1 || val == -1))
  {
    if (is_complex_singular_edge_after(_heh))
    {
      std::cerr << "d " << std::endl;

      return VH(-1);
    }
  }

  //get reference for updating properties
  //for transitions
  std::vector<std::vector<VH> > hfs_vhs;
  std::vector<int> hfs_trans;
  std::vector<Quaternion> hfs_trans_q, opp_hfs_trans_q;

  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    if (is_vt_bdy && mesh_.is_boundary(*hehf_it))
      continue;
    CH cur_ch = mesh_.incident_cell(*hehf_it);
    auto hfh_top = opposite_halfface_to_vertex_in_cell(vh_f, cur_ch);
    auto hfh_bottom = opposite_halfface_to_vertex_in_cell(vh_t, cur_ch);

    auto vhs = mesh_.get_halfface_vertices(hfh_top);
    hfs_vhs.push_back(vhs);

    int trans = tq_.mult_transitions_idx(trans_prop_[mesh_.opposite_halfface_handle(hfh_top)],
                                         trans_prop_[hfh_bottom]);
    hfs_trans.push_back(tq_.inverse_transition_idx(trans));
  }

  //for updating edge valence
  std::vector<int> vals;
  std::vector<VH> f_vhs, k_vhs;

  if (valence_[eh] != 0)
  {
    for (auto voh_it = mesh_.voh_iter(vh_f); voh_it.valid(); ++voh_it)
    {
      auto ve = mesh_.edge_handle(*voh_it);
      int val = valence_[ve];
      if (val != 0)
      {
        auto vh_tt = mesh_.halfedge(*voh_it).to_vertex();
        f_vhs.push_back(vh_tt);
        vals.push_back(val);

        if (keep_on_feature_face_[ve])
          k_vhs.push_back(vh_tt);
      }
    }
  }

  //for updating feature edge property
  std::vector<VH> feature_vhs;
  std::vector<int> fe_ids;
  if (feature_edge_[eh] > 0)
  {
    for (auto voh_it = mesh_.voh_iter(vh_f); voh_it.valid(); ++voh_it)
    {
      if (feature_edge_[mesh_.edge_handle(*voh_it)] > 0)
      {
        feature_vhs.push_back(mesh_.halfedge(*voh_it).to_vertex());
        fe_ids.push_back(feature_edge_[mesh_.edge_handle(*voh_it)]);
      }
    }
  }

  //for updating feature face property
  std::vector<std::vector<VH>> ff_vhs;
  std::vector<int> ff_ids;
  if (feature_face_vertex_[vh_f] || feature_face_vertex_[vh_t])
  {
    std::set<FH> e_inc_ff;
    for (auto hef_it = mesh_.hef_iter(_heh); hef_it.valid(); ++hef_it)
      if (feature_fprop_[*hef_it] > 0)
        e_inc_ff.insert(*hef_it);
    std::set<FH> inc_ff;
    for (auto vf_it = mesh_.vf_iter(vh_f); vf_it.valid(); ++vf_it)
      if (feature_fprop_[*vf_it] > 0 && e_inc_ff.find(*vf_it) == e_inc_ff.end())
      {
        auto hfvhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it, 0), vh_f);
        ff_vhs.push_back(hfvhs);
        ff_ids.push_back(feature_fprop_[*vf_it]);
      }

    for (auto vf_it = mesh_.vf_iter(vh_t); vf_it.valid(); ++vf_it)
      if (feature_fprop_[*vf_it] > 0 && e_inc_ff.find(*vf_it) == e_inc_ff.end())
      {
        auto hfvhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it, 0), vh_t);
        ff_vhs.push_back(hfvhs);
        ff_ids.push_back(feature_fprop_[*vf_it]);
      }
  }

  //store the aligned normal axis of neighbouring feature faces
  std::map<std::vector<VH>, int> inc_hfv_to_axis, ninc_hfv_to_axis;
  collect_axis_aligned_to_feature_face(_heh, inc_hfv_to_axis, ninc_hfv_to_axis);


  VH vh_survive = mesh_.collapse_edge(_heh);

  update_properties(eh, vh_survive, hfs_vhs, hfs_trans, vals, f_vhs, k_vhs, feature_vhs, fe_ids, ff_vhs, ff_ids,
                    inc_hfv_to_axis, ninc_hfv_to_axis);

  return vh_survive;
}

template<class MeshT>
bool
EdgeCollapseT<MeshT>::check_energy(const HEH _heh, const std::vector<std::vector<VH>> &_new_cell_vhs,
                                   Point &_collapse_pt, const bool _check_energy) const
{
  VH vh_f = mesh_.halfedge(_heh).from_vertex();
  VH vh_t = mesh_.halfedge(_heh).to_vertex();

  std::vector<std::vector<Point> > new_cells_points, old_cells_points;
  new_cells_points.reserve(_new_cell_vhs.size());

  //1. check if the new tets are degenerated
  for (const auto &vhs: _new_cell_vhs)
  {
    std::vector<Point> cell_points;
    cell_points.reserve(4);
    for (const auto &vhi: vhs)
      cell_points.push_back(mesh_.vertex(vhi));

    if (!RemeshingAssist<MeshT>::is_tet_valid(cell_points))
      return false;

    new_cells_points.push_back(cell_points);
  }

  //2. check the energy
  if (!opt_collapse_)
  {
    if (_check_energy)
    {
      for (auto vc_it = mesh_.vc_iter(vh_f); vc_it.valid(); ++vc_it)
      {
        auto cvhs = mesh_.get_cell_vertices(*vc_it);
        std::vector<Point> cell_points;
        cell_points.reserve(4);
        for (const auto &vhi: cvhs)
          cell_points.push_back(mesh_.vertex(vhi));
        old_cells_points.push_back(cell_points);
      }

      if (!RemeshingAssist<MeshT>::is_quality_improved(new_cells_points, old_cells_points, true))
        return false;
    }
  }
  else
  {
    //get new cells points and old cells points
    std::set<CH> collapsing_cells;
    for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
      collapsing_cells.insert(*hec_it);

    std::set<CH> incident_cells;
    for (auto vc_it = mesh_.vc_iter(vh_t); vc_it.valid(); ++vc_it)
      incident_cells.insert(*vc_it);

    std::vector<CH> surviving_cells2;
    std::set_difference(incident_cells.begin(), incident_cells.end(), collapsing_cells.begin(),
                        collapsing_cells.end(),
                        std::inserter(surviving_cells2, surviving_cells2.end()));


    //old cells points
    for (auto vc_it = mesh_.vc_iter(vh_f); vc_it.valid(); ++vc_it)
    {
      auto cvhs = mesh_.get_cell_vertices(*vc_it);
      std::vector<Point> cell_points;
      cell_points.reserve(4);
      for (const auto &vhi: cvhs)
        cell_points.push_back(mesh_.vertex(vhi));
      old_cells_points.push_back(cell_points);
    }

    for (const auto &chi: surviving_cells2)
    {
      auto cvhs = mesh_.get_cell_vertices(chi);
      std::vector<Point> cell_points;
      cell_points.reserve(4);
      for (const auto &vhi: cvhs)
        cell_points.push_back(mesh_.vertex(vhi));
      old_cells_points.push_back(cell_points);
    }

    //new cell points (start from vh_t)
    new_cells_points.reserve(_new_cell_vhs.size() + surviving_cells2.size());
    for (const auto &chi: surviving_cells2)
    {
      auto cvhs = mesh_.get_cell_vertices(chi, vh_t);
      std::vector<Point> cell_points;
      cell_points.reserve(4);
      for (const auto &vhi: cvhs)
        cell_points.push_back(mesh_.vertex(vhi));
      new_cells_points.push_back(cell_points);
    }

    //optimize new position of vh_t
    auto eh = mesh_.edge_handle(_heh);
    if (feature_edge_vertex_[vh_t] > 0 || sgl_vt_[vh_t])
    {//special vertex
      if (feature_node_[vh_t] > 0 || n_incident_special_edges(vh_t) != 2)
      {
        _collapse_pt = mesh_.vertex(vh_t);
      }
      else
      {//
        std::set<EH> spc_ehs;
        for (auto ve_it = mesh_.ve_iter(vh_t); ve_it.valid(); ++ve_it)
        {
          if (feature_edge_[*ve_it] > 0 || valence_[*ve_it] != 0)
            spc_ehs.insert(*ve_it);
        }
        if ((feature_edge_[eh] > 0 || valence_[eh] != 0) &&
            (feature_node_[vh_f] == 0 && n_incident_special_edges(vh_f) == 2))
        {//more space
          for (auto ve_it = mesh_.ve_iter(vh_f); ve_it.valid(); ++ve_it)
          {
            if (feature_edge_[*ve_it] > 0 || valence_[*ve_it] != 0)
              spc_ehs.insert(*ve_it);
          }
        }

        std::vector<Point> s_points;
        for (const auto &ehi: spc_ehs)
        {
          Point pt_s = mesh_.vertex(mesh_.edge(ehi).from_vertex());
          Point pt_e = mesh_.vertex(mesh_.edge(ehi).to_vertex());

          for (int i = 0; i <= 10; ++i)
          {
            s_points.push_back(
                    (double) i / (double) 10 * pt_s + (double) (10 - i) / (double) 10 * pt_e);
          }
        }

        _collapse_pt = vo_.optimize_vertex_on_edge_worst(s_points, new_cells_points);
      }
    }
    else if (feature_face_vertex_[vh_t])
    {//feature face vertex
      std::vector<Point> s_points;
      s_points.push_back(mesh_.vertex(vh_t));
      for (auto vf_it = mesh_.vf_iter(vh_t); vf_it.valid(); ++vf_it)
      {
        if (feature_fprop_[*vf_it] > 0)
        {
          auto hf_vhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it, 0));
          RemeshingAssist<MeshT>::sample_points(mesh_.vertex(hf_vhs[0]), mesh_.vertex(hf_vhs[1]),
                                                mesh_.vertex(hf_vhs[2]), 0, 2, s_points);
        }
      }
      if (feature_face_edge_[eh])
      {//safe. if ffe, vh_f must be a regular vertex if allowed to collapse
        for (auto vf_it = mesh_.vf_iter(vh_f); vf_it.valid(); ++vf_it)
        {
          if (feature_fprop_[*vf_it] > 0)
          {
            auto hf_vhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it, 0));
            RemeshingAssist<MeshT>::sample_points(mesh_.vertex(hf_vhs[0]), mesh_.vertex(hf_vhs[1]),
                                                  mesh_.vertex(hf_vhs[2]), 0, 2, s_points);
          }
        }
      }

      _collapse_pt = vo_.optimize_feature_face_vertex_location_worst(s_points, new_cells_points);
    }
    else
    {//regular vertex
      _collapse_pt = vo_.optimize_interior_vertex_location_worst(new_cells_points);
    }

    for (auto &c_pts: new_cells_points)
      c_pts[0] = _collapse_pt;

    if (!RemeshingAssist<MeshT>::is_quality_improved(new_cells_points, old_cells_points, true))
      return false;
  }


  return true;
}


template<class MeshT>
void
EdgeCollapseT<MeshT>::collect_axis_aligned_to_feature_face(const HEH _heh,
                                                           std::map<std::vector<VH>, int> &_inc_hfv_to_axis,
                                                           std::map<std::vector<VH>, int> &_ninc_hfv_to_axis)
{
  std::set<HFH> visited_hfhs;
  VH vh_f = mesh_.halfedge(_heh).from_vertex();
  VH vh_t = mesh_.halfedge(_heh).to_vertex();
  //in case the deleted bottom halfface has non-identity transition
  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    CH cur_ch = mesh_.incident_cell(*hehf_it);
    if (!cur_ch.is_valid())
      continue;

    auto hfh_top = opposite_halfface_to_vertex_in_cell(vh_f, cur_ch);
    auto hfh_bottom = opposite_halfface_to_vertex_in_cell(vh_t, cur_ch);

    if (feature_face_vertex_[vh_t])
    {
      auto top_fh = mesh_.face_handle(hfh_top);
      if (feature_fprop_[top_fh] > 0)
      {
        //top hfh
        visited_hfhs.insert(hfh_top);
        auto hfvhs = mesh_.get_halfface_vertices(hfh_top, vh_t);
        auto nm = mesh_.normal(hfh_top);

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[cur_ch], nm).second;

//                    if(trans_prop_[hfh_bottom] != 0) {
        axis = tq_.axis_after_transition(axis, trans_prop_[mesh_.opposite_halfface_handle(hfh_bottom)]);
//                    }

        _inc_hfv_to_axis.insert(std::make_pair(hfvhs, axis));

        //opposite halfface
        auto hfh_top_opp = mesh_.opposite_halfface_handle(hfh_top);
        visited_hfhs.insert(hfh_top_opp);

        auto ch_top_opp = mesh_.incident_cell(hfh_top_opp);
        if (!ch_top_opp.is_valid())
          continue;
        auto hfvhs_opp = mesh_.get_halfface_vertices(hfh_top_opp, vh_t);
        auto nm_opp = mesh_.normal(hfh_top_opp);

        int axis_opp = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_top_opp], nm_opp).second;

        _inc_hfv_to_axis.insert(std::make_pair(hfvhs_opp, axis_opp));
      }
    }
  }

  //feature faces incident to edge vertices
  if (feature_face_vertex_[vh_f])
  {
    for (auto vhf_it = mesh_.vhf_iter(vh_f); vhf_it.valid(); ++vhf_it)
    {
      auto fhi = mesh_.face_handle(*vhf_it);
      if (feature_fprop_[fhi] > 0 && visited_hfhs.find(*vhf_it) == visited_hfhs.end())
      {
        visited_hfhs.insert(*vhf_it);

        auto ch = mesh_.incident_cell(*vhf_it);
        if (!ch.is_valid())
          continue;

        auto hfvhs = mesh_.get_halfface_vertices(*vhf_it, vh_f);
        auto nm = mesh_.normal(*vhf_it);

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], nm).second;

        _inc_hfv_to_axis.insert(std::make_pair(hfvhs, axis));
      }
    }
  }

  if (feature_face_vertex_[vh_t])
  {
    for (auto vhf_it = mesh_.vhf_iter(vh_t); vhf_it.valid(); ++vhf_it)
    {
      auto fhi = mesh_.face_handle(*vhf_it);
      if (feature_fprop_[fhi] > 0 && visited_hfhs.find(*vhf_it) == visited_hfhs.end())
      {
        visited_hfhs.insert(*vhf_it);

        auto ch = mesh_.incident_cell(*vhf_it);
        if (!ch.is_valid())
          continue;

        auto hfvhs = mesh_.get_halfface_vertices(*vhf_it, vh_t);
        auto nm = mesh_.normal(*vhf_it);
        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], nm).second;

        _inc_hfv_to_axis.insert(std::make_pair(hfvhs, axis));
      }
    }
  }

  //feature faces not incident to edge vertices
  for (auto vc_it = mesh_.vc_iter(vh_f); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0 && visited_hfhs.find(*chf_it) == visited_hfhs.end())
      {
        visited_hfhs.insert(*chf_it);

        auto nm = mesh_.normal(*chf_it);
        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it],
                                                      nm).second;

        auto hfvhs = mesh_.get_halfface_vertices(*chf_it);

        _ninc_hfv_to_axis.insert(std::make_pair(hfvhs, axis));
      }
    }
  }

  for (auto vc_it = mesh_.vc_iter(vh_t); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0 && visited_hfhs.find(*chf_it) == visited_hfhs.end())
      {
        visited_hfhs.insert(*chf_it);

        auto nm = mesh_.normal(*chf_it);
        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it],
                                                      nm).second;

        auto hfvhs = mesh_.get_halfface_vertices(*chf_it);

        _ninc_hfv_to_axis.insert(std::make_pair(hfvhs, axis));
      }
    }
  }
}


template<class MeshT>
bool
EdgeCollapseT<MeshT>::is_boundary_edge_type_changed(const HEH _heh)
{
  VH vh_f = mesh_.halfedge(_heh).from_vertex();
  VH vh_t = mesh_.halfedge(_heh).to_vertex();

  std::set<HFH> gone_hfhs;
  std::vector<HFH> survive_hfhs;
  auto hfh_l = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(mesh_.opposite_halfedge_handle(_heh)));
  auto hfh_r = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(_heh));

  gone_hfhs.insert(hfh_l);
  gone_hfhs.insert(hfh_r);

  for (auto vhf_it = mesh_.vhf_iter(vh_f); vhf_it.valid(); ++vhf_it)
  {
    if (mesh_.is_boundary(*vhf_it))
    {
      if (gone_hfhs.count(*vhf_it) == 0)
        survive_hfhs.push_back(*vhf_it);
    }
  }

  //get new face normal
  std::map<HFH, Point> normal_map;
  for (const auto &hfh: survive_hfhs)
  {
    std::vector<VH> hfvhs = mesh_.get_halfface_vertices(hfh, vh_f);
    for (int i = 0; i < 3; ++i)
      if (hfvhs[i] == vh_f)
        hfvhs[i] = vh_t;
    auto normal = face_normal(mesh_, hfvhs);

    normal_map[hfh] = normal;
  }

  auto heh_l = mesh_.prev_halfedge_in_halfface(_heh, hfh_l);
  auto hfh_ll = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh_l));
  normal_map[hfh_l] = normal_map[hfh_ll];

  auto heh_r = mesh_.next_halfedge_in_halfface(mesh_.opposite_halfedge_handle(_heh), hfh_r);
  auto hfh_rr = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh_r));
  normal_map[hfh_r] = normal_map[hfh_rr];


  //get edges of boundary faces
  std::set<EH> rel_ehs;
  for (const auto &hfh: gone_hfhs)
  {
    for (auto hfe_it = mesh_.hfe_iter(hfh); hfe_it.valid(); ++hfe_it)
      rel_ehs.insert(*hfe_it);
  }
  for (const auto &hfh: survive_hfhs)
  {
    for (auto hfe_it = mesh_.hfe_iter(hfh); hfe_it.valid(); ++hfe_it)
      rel_ehs.insert(*hfe_it);
  }

  //remove three vanishing edges
  rel_ehs.erase(mesh_.edge_handle(heh_l));
  rel_ehs.erase(mesh_.edge_handle(_heh));
  rel_ehs.erase(mesh_.edge_handle(heh_r));


  for (const auto ehi: rel_ehs)
  {
    auto heh0 = mesh_.halfedge_handle(ehi, 0);
    auto heh1 = mesh_.opposite_halfedge_handle(heh0);
    auto vh0 = mesh_.halfedge(heh0).from_vertex();
    auto vh1 = mesh_.halfedge(heh1).from_vertex();

    auto hfhi_l = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh1));
    auto hfhi_r = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh0));

    Point newnl, newnr;
    if (normal_map.find(hfhi_r) == normal_map.end())
      newnr = mesh_.normal(hfhi_r);
    else
      newnr = normal_map[hfhi_r];

    if (normal_map.find(hfhi_l) == normal_map.end())
      newnl = mesh_.normal(hfhi_l);
    else
      newnl = normal_map[hfhi_l];

    if (!mesh_.is_boundary(hfhi_r) || !mesh_.is_boundary(hfhi_l))
    {
      std::cerr << "Error: not boundary! " << mesh_.face_handle(hfhi_r) << " " << mesh_.face_handle(hfhi_l)
                << std::endl;

      std::cerr << "edge: " << ehi << " incident cells: ";
      for (auto ec_it = mesh_.ec_iter(ehi); ec_it.valid(); ++ec_it)
        std::cerr << " " << *ec_it << " is bdy? " << mesh_.is_boundary(*ec_it);
      std::cerr << " incident faces: ";
      for (auto ef_it = mesh_.ef_iter(ehi); ef_it.valid(); ++ef_it)
        std::cerr << " " << *ef_it << " is bdy? " << mesh_.is_boundary(*ef_it);
      std::cerr << std::endl;
    }
    //check normal angle changes
    auto angle = dihedral_angle(mesh_.vertex(vh0), mesh_.vertex(vh1), newnl, newnr);
    auto old_angle = dihedral_angle(mesh_.vertex(vh0), mesh_.vertex(vh1), mesh_.normal(hfhi_l), mesh_.normal(hfhi_r));

    if (!std::isfinite(angle) || !std::isfinite(old_angle))
    {
      for (auto nm: normal_map)
        std::cout << "nm: " << nm.first << "->" << nm.second << std::endl;
      std::cout << " angle wrong: " << angle << " " << old_angle << " nl " << newnl << " nr: " << newnr << " fh: "
                << hfhi_l
                << " " << hfhi_r << std::endl;
    }


    int old_type = RemeshingAssist<MeshT>::boundary_edge_type(old_angle);
    int new_type = RemeshingAssist<MeshT>::boundary_edge_type(angle);

    //the edge dihedral angle changed
    //allow flatten sharp edge but not marked as sg
    if (old_type != new_type && (new_type != 0 || valence_[ehi] != 0))
    {
      return true;
    }
  }

  return false;
}

//topology check
template<class MeshT>
typename EdgeCollapseT<MeshT>::CollapseType
EdgeCollapseT<MeshT>::is_collapse_ok(const HEH _heh) const
{
  auto eh = mesh_.edge_handle(_heh);
  //if edge is deleted
  if (mesh_.is_deleted(eh))
    return DOFError;

  auto from_vh = mesh_.halfedge(_heh).from_vertex();
  auto to_vh = mesh_.halfedge(_heh).to_vertex();

  //if from vertex is feature node
  if (feature_node_[from_vh] > 0)
    return DOFError;

  //if from vertex is feature edge vertex but the edge is not feature
  if (feature_edge_vertex_[from_vh] && feature_edge_[eh] == 0)
    return DOFError;

  //if from vertex is feature face vertex but the edge is not feature face edge
  if (feature_face_vertex_[from_vh] && !feature_face_edge_[eh])
    return DOFError;

  bool is_vhf_bdy = mesh_.is_boundary(from_vh);
  bool is_vht_bdy = mesh_.is_boundary(to_vh);
  bool is_e_bdy = mesh_.is_boundary(eh);

  if (sgl_vt_[from_vh] && valence_[eh] == 0)
    return DOFError;


  //if from vertex is boundary but the edge is not
  if (is_vhf_bdy && !is_e_bdy)
    return DOFError;

  //1. if to vertex is on feature face (not boundary) but from vertex is a zipper node
  //2. if to vertex is on boundary face, and from vertex is a zipper node, the other to vertex is feature node
  //3. if to vertex is on feature face and is feature node, but from vertec is a zipper node

  if (valence_[eh] != 0 && !feature_face_edge_[eh] && feature_face_vertex_[to_vh])
  {
    int tid = MeshPropertiesT<MeshT>::node_index(from_vh);
    if (tid == 10)
    {
      //1
      if (!is_vht_bdy)
        return DOFError;
      else
      {
        if (feature_node_[to_vh] > 0)
          return DOFError;
        else
        {
          for (auto ve_it = mesh_.ve_iter(from_vh); ve_it.valid(); ++ve_it)
          {
            if (valence_[*ve_it] != 0 && *ve_it != eh)
            {
              VH other_to_vh = mesh_.edge(*ve_it).from_vertex() == from_vh ? mesh_.edge(*ve_it).to_vertex() :
                               mesh_.edge(*ve_it).from_vertex();
              if (feature_node_[other_to_vh] > 0)
                return DOFError;
            }
          }
        }
      }

    }
  }

  bool is_vhf_node = is_singular_node(from_vh);
  if (post_remesh_ && is_vhf_node)
    return DOFError;

  if (is_vhf_node && feature_node_[to_vh] > 0)
    return DOFError;

  bool is_vht_node = is_singular_node(to_vh);
  //don't collapse the admissable singular edges if the singular node will touch the boundary
  bool ad_val = std::abs(valence_[eh]) == 1 || std::abs(valence_[eh]) == 2;
  if (ad_val && is_vhf_node && is_vht_bdy)
    return DOFError;

  //don't collapse the admissable singular edges if both vertices are singular node
  if (ad_val && is_vht_node && is_vhf_node)
  {
    int fid = MeshPropertiesT<MeshT>::node_index(from_vh);
    int tid = MeshPropertiesT<MeshT>::node_index(to_vh);
    if (fid < 10 && fid > 0 && tid < 10 && tid > 0)
      return DOFError;
  }

//        //don't collapse singular edge if vhf is a zipper node
  // and one of the vht is feature node and no hope to fix with guiding features
  if (valence_[eh] != 0)
  {
    if (is_vht_bdy && is_zipper_node(from_vh))
    {
      for (auto voh_it = mesh_.voh_iter(from_vh); voh_it.valid(); ++voh_it)
      {
        if (valence_[mesh_.edge_handle(*voh_it)] != 0)
        {
          auto vht = mesh_.halfedge(*voh_it).to_vertex();
          if (feature_node_[vht] > 0
//                        && !has_boundary_valence_minus_one(vht)
                  )
          {
            return DOFError;
          }
        }
      }
    }
  }

  //prevent an interior edge having two ffvs after edge collapse
  if (feature_face_vertex_[to_vh] && !feature_face_vertex_[from_vh])
  {
    for (const auto vvh: mesh_.vertex_vertices(from_vh))
    {
      if (feature_face_vertex_[vvh] && vvh != to_vh)
      {
        HEH heh = mesh_.find_halfedge(vvh, to_vh);
        if (!heh.is_valid())
          return DOFError;
        else
        {
          if (!feature_face_edge_[mesh_.edge_handle(heh)])
            return DOFError;
          else
          {
            for (auto voh_it = mesh_.voh_iter(vvh); voh_it.valid(); ++voh_it)
            {
              auto ve = mesh_.edge_handle(*voh_it);
              if (feature_face_edge_[ve])
              {
                auto vvhi = mesh_.halfedge(*voh_it).to_vertex();
                if (vvhi != to_vh)
                {
                  HEH heh_other = mesh_.find_halfedge(vvhi, from_vh);
                  if (heh_other.is_valid())
                    return DOFError;
                }
              }
            }
          }
        }
      }
    }
  }


  return is_topology_ok(_heh);
}

template<class MeshT>
typename EdgeCollapseT<MeshT>::CollapseType
EdgeCollapseT<MeshT>::
is_topology_ok(const HEH _heh) const
{
  auto from_vh = mesh_.halfedge(_heh).from_vertex();
  auto to_vh = mesh_.halfedge(_heh).to_vertex();
  //Link condition check. This is maybe not necessary since the positive tet volume check implicitly take care of this.
  //1.vertex in face
  //                / \
        //               / | \
        //              /  .  \
        //             / /   \ \
        //            /_________\
        //check the intersection of the one-rings of from_vh and to_vh,
  std::set<VH> vvh_from, vvh_to;
  for (const auto vh: mesh_.vertex_vertices(from_vh))
    vvh_from.insert(vh);
  for (const auto vh: mesh_.vertex_vertices(to_vh))
    vvh_to.insert(vh);
  std::vector<VH> vvh_ints;
  std::set_intersection(vvh_from.begin(), vvh_from.end(), vvh_to.begin(), vvh_to.end(),
                        std::inserter(vvh_ints, vvh_ints.end()));
////
  std::vector<VH> hf_vhs(3);
  hf_vhs[0] = from_vh;
  hf_vhs[1] = to_vh;
  for (const auto vh: vvh_ints)
  {
    hf_vhs[2] = vh;
    auto hfh = mesh_.find_halfface(hf_vhs);

    if (!hfh.is_valid() || mesh_.is_deleted(hfh))
    {
      return LinkError;
    }
  }

  //2. vertex in tetrahedron
  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    VH r = opposite_vertex_to_halfface_in_cell(*hehf_it);
    VH v = opposite_vertex_to_halfface_in_cell(mesh_.opposite_halfface_handle(*hehf_it));
    VH s = mesh_.halfedge(mesh_.next_halfedge_in_halfface(_heh, *hehf_it)).to_vertex();

    if (r.is_valid() && v.is_valid())
    {
      if (mesh_.find_halfedge(r, v).is_valid())
      {
        auto hf_qrv = mesh_.find_halfface(std::vector<VH>{to_vh, r, v});
        if (hf_qrv.is_valid())
        {
          auto inc_ch = mesh_.incident_cell(hf_qrv);
          if (inc_ch.is_valid())
          {
            auto is_s = opposite_vertex_to_halfface_in_cell(hf_qrv);
            if (is_s == s)
            {
              //find the second tet
              auto hf_srv = mesh_.find_halfface(std::vector<VH>{s, r, v});
              if (hf_srv.is_valid())
              {
                auto inc_ch2 = mesh_.incident_cell(hf_srv);
                if (inc_ch2.is_valid())
                {
                  auto is_p = opposite_vertex_to_halfface_in_cell(hf_srv);
                  if (is_p == from_vh)
                  {
                    return LinkError;
                  }
                }
              }
            }
          }
        }
      }
    }

  }
////

  //valance check: no valance-2 vertex after collapsing is allowed
  if (vvh_from.size() + vvh_to.size() - vvh_ints.size() - 2 < 3)
    return LinkError;

  return CollapseOK;
}

template<class MeshT>
std::vector<EH>
EdgeCollapseT<MeshT>::
get_blocking_edges(const EH _eh) const
{
  std::vector<EH> blk_ehs;

  auto from_vh = mesh_.edge(_eh).from_vertex();
  auto to_vh = mesh_.edge(_eh).to_vertex();

  //check the link condition
  //1.
  std::set<VH> vvh_from, vvh_to;
  for (const auto vh: mesh_.vertex_vertices(from_vh))
    vvh_from.insert(vh);
  for (const auto vh: mesh_.vertex_vertices(to_vh))
    vvh_to.insert(vh);
  std::vector<VH> vvh_ints;
  std::set_intersection(vvh_from.begin(), vvh_from.end(), vvh_to.begin(), vvh_to.end(),
                        std::inserter(vvh_ints, vvh_ints.end()));

  std::vector<VH> hf_vhs(3);
  hf_vhs[0] = from_vh;
  hf_vhs[1] = to_vh;
  for (const auto vh: vvh_ints)
  {
    hf_vhs[2] = vh;
    auto hfh = mesh_.find_halfface(hf_vhs);

    if (!hfh.is_valid() || mesh_.is_deleted(hfh))
    {
      blk_ehs.push_back(mesh_.edge_handle(mesh_.find_halfedge(from_vh, vh)));
//                blk_ehs.push_back(mesh_.edge_handle(mesh_.halfedge(to_vh, vh)));
    }
  }


  //2. vertex in tetrahedron
  auto heh = mesh_.halfedge_handle(_eh, 0);
  for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
  {
    VH r = opposite_vertex_to_halfface_in_cell(*hehf_it);
    VH v = opposite_vertex_to_halfface_in_cell(mesh_.opposite_halfface_handle(*hehf_it));
    VH s = mesh_.halfedge(mesh_.next_halfedge_in_halfface(heh, *hehf_it)).to_vertex();

    if (r.is_valid() && v.is_valid())
    {
      if (mesh_.find_halfedge(r, v).is_valid())
      {
        auto he_adj = mesh_.find_halfedge(to_vh, v);
        if (he_adj.is_valid() && !mesh_.is_deleted(he_adj))
        {
          auto hf_qrv = mesh_.find_halfface(std::vector<VH>{to_vh, r, v});
          if (hf_qrv.is_valid())
          {
            auto inc_ch = mesh_.incident_cell(hf_qrv);
            if (inc_ch.is_valid())
            {
              auto is_s = opposite_vertex_to_halfface_in_cell(hf_qrv);
              if (is_s == s)
              {
                //find the second tet
                auto hf_srv = mesh_.find_halfface(std::vector<VH>{s, r, v});
                if (hf_srv.is_valid())
                {
                  auto inc_ch2 = mesh_.incident_cell(hf_srv);
                  if (inc_ch2.is_valid())
                  {
                    auto is_p = opposite_vertex_to_halfface_in_cell(hf_srv);
                    if (is_p == from_vh)
                    {
                      blk_ehs.push_back(mesh_.edge_handle(he_adj));
                      auto he_adj2 = mesh_.find_halfedge(from_vh, v);
                      if (he_adj2.is_valid() && !mesh_.is_deleted(he_adj2))
                        blk_ehs.push_back(mesh_.edge_handle(he_adj2));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return blk_ehs;
}

//get the new tets after edge collapse(just simulate, not real collapsing)
template<class MeshT>
void
EdgeCollapseT<MeshT>::surviving_tets(const HEH _heh, std::vector<std::vector<VH> > &_v_vhs) const
{
  auto from_vh = mesh_.halfedge(_heh).from_vertex();
  auto to_vh = mesh_.halfedge(_heh).to_vertex();

  // find cells that will collapse, i.e. are incident to the collapsing halfedge
  std::set<CH> collapsing_cells;
  for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
    collapsing_cells.insert(*hec_it);

  std::set<CH> incident_cells;
  for (auto vc_it = mesh_.vc_iter(from_vh); vc_it.valid(); ++vc_it)
    incident_cells.insert(*vc_it);

  std::vector<CH> surviving_cells;
  std::set_difference(incident_cells.begin(), incident_cells.end(), collapsing_cells.begin(), collapsing_cells.end(),
                      std::inserter(surviving_cells, surviving_cells.end()));

  for (const auto chi: surviving_cells)
  {
    auto cell_vertices = mesh_.get_cell_vertices(chi, from_vh);
    cell_vertices[0] = to_vh;
    _v_vhs.push_back(cell_vertices);
  }
}


template<class MeshT>
void
EdgeCollapseT<MeshT>::update_properties(const EH _eh, const VH _vh_new, const std::vector<std::vector<VH> > &_hfs_vhs,
                                        const std::vector<int> &_hfs_trans, std::vector<int> &_vals,
                                        const std::vector<VH> &_f_vhs, const std::vector<VH> &_k_vhs,
                                        const std::vector<VH> &_fe_vhs, const std::vector<int> &_fe_ids,
                                        const std::vector<std::vector<VH>> &_ff_vhs, const std::vector<int> &_ff_ids,
                                        std::map<std::vector<VH>, int> &_inc_hfvhs_axis,
                                        std::map<std::vector<VH>, int> &_ninc_hfvhs_axis)
{
  //update halfface transition
  int i = 0;
  for (const auto &vhs: _hfs_vhs)
  {
    auto hfh = mesh_.find_halfface(vhs);
    if (!hfh.is_valid())
      std::cout << "\nError: halfface is invalid after edge collapse!";

    if (!mesh_.is_boundary(mesh_.face_handle(hfh)))
    {
      trans_prop_[hfh] = _hfs_trans[i];
      trans_prop_[mesh_.opposite_halfface_handle(hfh)] = tq_.inverse_transition_idx(_hfs_trans[i]);
    }
    else
    {
      trans_prop_[hfh] = -1;
      trans_prop_[mesh_.opposite_halfface_handle(hfh)] = -1;
    }

    ++i;
  }

  //update valence
  for (auto i = 0u; i < _f_vhs.size(); ++i)
  {
    auto heh = mesh_.find_halfedge(_f_vhs[i], _vh_new);
    if (heh.is_valid())
    {
      valence_[mesh_.edge_handle(heh)] = _vals[i];
    }
  }

  //update kept on feature face
  for (auto i = 0u; i < _k_vhs.size(); ++i)
  {
    auto heh = mesh_.find_halfedge(_k_vhs[i], _vh_new);
    if (heh.is_valid())
    {
      keep_on_feature_face_[mesh_.edge_handle(heh)] = true;
    }
  }

  //update feature edges
  for (auto i = 0u; i < _fe_vhs.size(); ++i)
  {
    auto heh = mesh_.find_halfedge(_fe_vhs[i], _vh_new);
    if (heh.is_valid())
    {
      feature_edge_[mesh_.edge_handle(heh)] = _fe_ids[i];
      feature_edge_vertex_[mesh_.halfedge(heh).from_vertex()] = true;
      feature_edge_vertex_[mesh_.halfedge(heh).to_vertex()] = true;
    }
  }

  //update feature faces
  for (auto i = 0u; i < _ff_vhs.size(); ++i)
  {
    auto hfvnew = _ff_vhs[i];
    hfvnew[0] = _vh_new;

    auto hf_new = mesh_.find_halfface(hfvnew);
    if (hf_new.is_valid())
    {
      feature_fprop_[mesh_.face_handle(hf_new)] = _ff_ids[i];

      for (auto hfe_it = mesh_.hfe_iter(hf_new); hfe_it.valid(); ++hfe_it)
        feature_face_edge_[*hfe_it] = true;

      for (auto vhi: hfvnew)
        feature_face_vertex_[vhi] = true;
    }
  }

  //smooth quaternions
  std::vector<CH> cells;
  for (auto vc_it = mesh_.vc_iter(_vh_new); vc_it.valid(); ++vc_it)
    cells.push_back(*vc_it);

  std::set<FH> bdy_fhs;
  for (const auto ch: cells)
    for (auto cf_it = mesh_.cf_iter(ch); cf_it.valid(); ++cf_it)
      if (mesh_.is_boundary(*cf_it))
        bdy_fhs.insert(*cf_it);

  std::map<CH, std::vector<std::pair<Vec3d, int>>> alignment, alignment2;
  //feature face alignment
  std::map<HFH, int> hfh_to_axis;
  for (auto&[hfvhs, axis]: _inc_hfvhs_axis)
  {
    auto hfvhs_temp = hfvhs;
    hfvhs_temp[0] = _vh_new;
    auto hf_new = mesh_.find_halfface(hfvhs_temp);
    if (hf_new.is_valid())
    {
//                std::cerr<<"hfvhs: "<<hfvhs[0]<<" "<<hfvhs[1]<<" "<<hfvhs[2]<<std::endl;

      if (feature_fprop_[mesh_.face_handle(hf_new)] == 0)
        std::cerr << "Error: new halfface is not feature!" << std::endl;
      hfh_to_axis.insert(std::make_pair(hf_new, axis));
    }
  }
//        std::cerr<<"ninc ";
  for (auto&[hfvhs, axis]: _ninc_hfvhs_axis)
  {
    auto hf_new = mesh_.find_halfface(hfvhs);
    if (hf_new.is_valid())
    {
//                std::cerr<<"hfvhs: "<<hfvhs[0]<<" "<<hfvhs[1]<<" "<<hfvhs[2]<<std::endl;
      if (feature_fprop_[mesh_.face_handle(hf_new)] == 0)
        std::cerr << "Error: new halfface is not feature!" << std::endl;
      hfh_to_axis.insert(std::make_pair(hf_new, axis));
    }
  }

  for (auto&[hfh, axis]: hfh_to_axis)
  {
    auto ch_new = mesh_.incident_cell(hfh);
    if (!ch_new.is_valid())
      continue;

    auto nm = mesh_.normal(hfh);
    alignment[ch_new].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
  }

  alignment2 = alignment;
//        //check
//        for(const auto ch : cells) {
//            for(auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it){
//                auto fhi = mesh_.face_handle(*chf_it);
//                if (feature_fprop_[fhi] > 0) {
//                    if(alignment.find(ch) == alignment.end()) {
//                        std::cerr<<"Error: feature face cell has no alignment!"<<std::endl;
//                    } else {
//                        auto normal = mesh_.normal(*chf_it);
//                        Vec3d nm(normal[0], normal[1],normal[2]);
//                        Vec3d aligned_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch], (AxisAlignment)alignment[ch].front().second);
//                        Vec3d aligned_vec_real = alignment[ch].front().first;
//
//                        double prod = nm.dot(aligned_vec);
//                        if(prod < 0.5) {
//                            std::cerr<<"Warning: in collapsing edge "<<_eh<<" before clpsing ch "<<ch<<" fhi "<<fhi <<" dpd: "<<prod<<" vhnew "<<_vh_new<<" #alignments "<<alignment[ch].size()<<" normal "<<normal<<" aligned axis "<<alignment[ch].front().second<<" vec qtn "<<aligned_vec.transpose()<<" vec real "<<aligned_vec_real.transpose()<<std::endl;
//                        }
//                    }
//                }
//            }
//            if(alignment[ch].size() > 1) {
//                std::cerr<<"Error: cell "<<ch<<" has more than one feeature face constraints at vh "<<_vh_new;
//                for(auto [nm, ax] : alignment[ch]) {
//                    std::cerr<<" dir "<<nm[0]<<" "<<nm[1]<<" "<<nm[2]<<" ax "<<ax;
//                }
//                std::cerr<<std::endl;
//            }
//        }


  //feature edge alignment
  for (const auto ch: cells)
  {
    for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
    {
      if (feature_edge_[*ce_it] > 0)
      {
        auto dir = mesh_.vertex(mesh_.edge(*ce_it).to_vertex()) - mesh_.vertex(mesh_.edge(*ce_it).from_vertex());
        dir.normalize();

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], dir).second;
        alignment[ch].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
      }
    }
  }

  //check
  for (const auto ch: cells)
  {
    if (alignment[ch].size() == 2)
    {
      if (alignment[ch][0].second == alignment[ch][1].second)
      {
        std::cerr << "Error: cell " << ch << " has more than one constraints of the same axis in collapsing edge "
                  << _eh << std::endl;
        for (auto[nm, ax]: alignment[ch])
        {
          std::cerr << " dir " << nm[0] << " " << nm[1] << " " << nm[2] << " ax " << ax;
        }
        std::cerr << std::endl;
      }
    }
    else if (alignment[ch].size() > 2)
    {
      std::cerr << "Error: cell " << ch << " has " << alignment[ch].size() << " constraints in collapsing edge " << _eh
                << std::endl;
      for (auto[nm, ax]: alignment[ch])
      {
        std::cerr << " dir " << nm[0] << " " << nm[1] << " " << nm[2] << " ax " << ax;
      }
      std::cerr << std::endl;
    }
  }

  for (int i = 0; i < 10; ++i)
    QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, cells, bdy_fhs, alignment, tq_, cell_quaternions_);


  //check
  for (const auto ch: cells)
  {
    for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0)
      {
        if (alignment2.find(ch) == alignment2.end())
        {
          std::cerr << "Error: feature face cell has no alignment2!" << std::endl;
        }
        else
        {
          auto normal = mesh_.normal(*chf_it);
          Vec3d nm(normal[0], normal[1], normal[2]);
          Vec3d aligned_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch],
                                                                      (AxisAlignment) alignment2[ch].front().second);
          Vec3d aligned_vec_real = alignment2[ch].front().first;

          double prod = nm.dot(aligned_vec);
          if (prod < 0.5)
          {
            std::cerr << "Warning22: in collapsing edge " << _eh << " before clpsing ch " << ch << " fhi " << fhi
                      << " dpd: " << prod
                      << " vhnew " << _vh_new << " #alignments " << alignment2[ch].size() << " normal " << normal
                      << " aligned axis " << alignment2[ch].front().second << " vec qtn " << aligned_vec.transpose()
                      << " vec real " << aligned_vec_real.transpose() << std::endl;
          }
        }
      }
    }
    if (alignment2[ch].size() > 1)
    {
      std::cerr << "Error: cell " << ch << " has more than one feeature face constraints at vh " << _vh_new;
      for (auto[nm, ax]: alignment2[ch])
      {
        std::cerr << " dir " << nm[0] << " " << nm[1] << " " << nm[2] << " ax " << ax;
      }
      std::cerr << std::endl;
    }
  }
}

template<class MeshT>
HFH
EdgeCollapseT<MeshT>::opposite_halfface_to_vertex_in_cell(const VH _vh, const CH _ch) const
{
  if (!_ch.is_valid())
  {
    std::cout << "\nError: Invalid cell!";
    return HFH(-1);
  }
  auto cvhs = mesh_.get_cell_vertices(_ch, _vh);
  std::vector<VH> fvhs = {cvhs[2], cvhs[1], cvhs[3]};

  return mesh_.find_halfface(fvhs);
}

template<class MeshT>
VH
EdgeCollapseT<MeshT>::opposite_vertex_to_halfface_in_cell(const HFH _hfh) const
{
  if (mesh_.is_boundary(_hfh))
    return VH(-1);

  CH ch = mesh_.incident_cell(_hfh);
  std::set<VH> fvhs, cvhs;
  for (auto hfv_it = mesh_.hfv_iter(_hfh); hfv_it.valid(); ++hfv_it)
    fvhs.insert(*hfv_it);

  for (auto cv_it = mesh_.cv_iter(ch); cv_it.valid(); ++cv_it)
    if (fvhs.find(*cv_it) == fvhs.end())
      return *cv_it;

  return VH(-1);
}

template<class MeshT>
bool
EdgeCollapseT<MeshT>::is_singular_node(const VH _vh) const
{
  int n = 0;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    if (valence_[mesh_.edge_handle(*voh_it)] != 0)
      n++;

  if (n > 2)
    return true;
  return false;
}

template<class MeshT>
bool
EdgeCollapseT<MeshT>::is_zipper_node(const VH _vh) const
{
  int n_val1 = 0, n_val_1 = 0;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto veh = mesh_.edge_handle(*voh_it);
    if (valence_[veh] == -1)
      n_val_1++;
    else if (valence_[veh] == 1)
      n_val1++;
  }

  if (n_val_1 == 1 || n_val1 == 1)
    return true;

  return false;
}

template<class MeshT>
int
EdgeCollapseT<MeshT>::
n_incident_special_edges(const VH _vh) const
{
  int n = 0;
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (valence_[*ve_it] != 0 || feature_edge_[*ve_it] > 0)
      n++;

  return n;
}

template<class MeshT>
bool
EdgeCollapseT<MeshT>::is_zipper_node_collapsable(const VH _vh) const
{
  int n_val1 = 0, n_val_1 = 0;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto veh = mesh_.edge_handle(*voh_it);
    if (valence_[veh] == -1)
      n_val_1++;
    else if (valence_[veh] == 1)
      n_val1++;
  }

  if (n_val_1 < 2 && n_val_1 + n_val1 == 2)
    return true;

  return false;
}


template<class MeshT>
bool
EdgeCollapseT<MeshT>::has_boundary_valence_minus_one(const VH _vh) const
{
  int n_val_1 = 0;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto veh = mesh_.edge_handle(*voh_it);
    if (valence_[veh] == -1)
      return true;
  }

  return false;
}


template<class MeshT>
bool
EdgeCollapseT<MeshT>::is_complex_singular_edge_after(const HEH _heh) const
{
  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    HEH heh_next = mesh_.next_halfedge_in_halfface(_heh, *hehf_it);
    //should not be a boundary edge
    HEH heh_next_opp = mesh_.opposite_halfedge_handle(heh_next);

    HEH heh_prev = mesh_.prev_halfedge_in_halfface(_heh, *hehf_it);

    int trans_next = halfedge_transition_index(heh_next_opp, mesh_.opposite_halfface_handle(*hehf_it));
    trans_next = tq_.mult_transitions_idx(trans_next, trans_prop_[*hehf_it]);

    HFH hfh_bottom = mesh_.adjacent_halfface_in_cell(*hehf_it, heh_prev);
    int trans_prev = halfedge_transition_index(heh_prev, mesh_.opposite_halfface_handle(hfh_bottom));
    trans_prev = tq_.mult_transitions_idx(tq_.inverse_transition_idx(trans_prop_[*hehf_it]), trans_prev);

    int trans_new = tq_.mult_transitions_idx(trans_next, trans_prev);

    if (
//                    trans_new > 9 ||
            (trans_new < 4 && trans_new > 0))
      return true;
  }

  return false;
}


template<class MeshT>
int
EdgeCollapseT<MeshT>::halfedge_transition_index(const HEH _heh, const HFH _hfh_s) const
{
  auto hfh_it = _hfh_s;
  //if starting from boundary face, go to the neighbour
  if (mesh_.is_boundary(mesh_.face_handle(hfh_it)))
    hfh_it = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hfh_it, _heh));

  int idx = 0;
  do
  {
    if (mesh_.is_boundary(hfh_it))
      break;

    idx = tq_.mult_transitions_idx(trans_prop_[hfh_it], idx);

    auto hfh_adj = mesh_.adjacent_halfface_in_cell(hfh_it, _heh);
    if (!hfh_adj.is_valid())
    {
      std::cerr << "Error: adjacent halfface is invalid!" << std::endl;
      break;
    }
    hfh_it = mesh_.opposite_halfface_handle(hfh_adj);
  }
  while (hfh_it != _hfh_s);


  return idx;
}

template<class MeshT>
EH
EdgeCollapseT<MeshT>::opposite_edge_in_cell(const EH _eh, const CH _ch) const
{
  VH vh0 = mesh_.edge(_eh).from_vertex();
  VH vh1 = mesh_.edge(_eh).to_vertex();

  std::vector<VH> vhs;
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
  {
    if (*cv_it != vh0 && *cv_it != vh1)
    {
      vhs.push_back(*cv_it);
    }
  }

  return mesh_.edge_handle(mesh_.find_halfedge(vhs[0], vhs[1]));
}



//==================================================================================================================//
template<class MeshT>
std::set<EH>
EdgeSwapT<MeshT>::edge_swap(const EH _eh)
{
  bool is_bdy = mesh_.is_boundary(_eh);
  bool is_ft = feature_face_edge_[_eh];
  //don't do swap if the dihedral angle changes a lot
  if (is_bdy)
  {
    if (is_boundary_edge_type_changed(_eh))
      return std::set<EH>();
  }

  //0. get equatorial vertices
  auto v_equatorial_vertices = get_equatorial_vertices(_eh);
  if (v_equatorial_vertices[0].size() < 3 || v_equatorial_vertices[0].size() > 5)
  {
    if (v_equatorial_vertices[0].empty())
      std::cerr << "Error: no incident faces at edge " << _eh << std::endl;
    return std::set<EH>();
  }

  if (!v_equatorial_vertices[1].empty())
    if (v_equatorial_vertices[1].size() < 3 || v_equatorial_vertices[1].size() > 5)
      return std::set<EH>();


  if (feature_face_edge_[_eh])
  {
    //if a regular edge after swap connect two feature edge vertices
    if (feature_edge_vertex_[v_equatorial_vertices[0][0]] && feature_edge_vertex_[v_equatorial_vertices[0].back()])
      return std::set<EH>();

    //if a regular edge after swap connect two singular vertices
    if (post_remesh_)
    {
      if (sgl_vt_[v_equatorial_vertices[0][0]] && sgl_vt_[v_equatorial_vertices[0].back()])
        return std::set<EH>();
    }
  }

  auto cur_heh = mesh_.halfedge_handle(_eh, 0);
  auto vh_from = mesh_.halfedge(cur_heh).from_vertex();
  auto vh_to = mesh_.halfedge(cur_heh).to_vertex();

  //
  for (auto hehf_it = mesh_.hehf_iter(cur_heh); hehf_it.valid(); ++hehf_it)
  {
    HEH he_next = mesh_.next_halfedge_in_halfface(cur_heh, *hehf_it);
    EH eh_n = mesh_.edge_handle(he_next);
    HEH he_prev = mesh_.prev_halfedge_in_halfface(cur_heh, *hehf_it);
    EH eh_p = mesh_.edge_handle(he_prev);


    if ((valence_[eh_n] > 1 && valence_[eh_n] < std::numeric_limits<int>::max())
        || (valence_[eh_p] > 1 && valence_[eh_p] < std::numeric_limits<int>::max()))
      return std::set<EH>();
  }


  //1. check if energy is improved after and find the best configuration
  std::vector<int> v_ref_pos;
  for (const auto &equatorial_vertices: v_equatorial_vertices)
  {
    if (equatorial_vertices.empty())
      continue;

    int ref_pos = -1;
    std::vector<std::vector<Point> > old_cells_points;
    old_cells_points.reserve(equatorial_vertices.size());

    if (v_equatorial_vertices[1].empty())
    {
      for (const auto &hec: mesh_.halfedge_cells(cur_heh))
      {
        auto cell_vertices = mesh_.get_cell_vertices(hec);
        std::vector<Point> cell_points;
        cell_points.reserve(4);
        for (const auto &vh: cell_vertices)
          cell_points.push_back(mesh_.vertex(vh));

        old_cells_points.push_back(cell_points);
      }
    }
    else
    {
      for (int ii = 0; ii < (int) equatorial_vertices.size() - 1; ++ii)
      {
        std::vector<Point> cell_points;
        cell_points.reserve(4);

        cell_points.push_back(mesh_.vertex(vh_from));
        cell_points.push_back(mesh_.vertex(vh_to));
        cell_points.push_back(mesh_.vertex(equatorial_vertices[ii]));
        cell_points.push_back(mesh_.vertex(equatorial_vertices[ii + 1]));

        old_cells_points.push_back(cell_points);
      }
    }

    ref_pos = find_the_best_configuration(_eh, equatorial_vertices, old_cells_points);
    v_ref_pos.push_back(ref_pos);
  }

  //both sides should have better quality
  for (auto pos: v_ref_pos)
    if (pos == -1)
      return std::set<EH>();

  //2. swap
  std::set<EH> new_ehs;
  if (!v_ref_pos.empty())
  {
    HFH hf_start = find_start_halfface(cur_heh, is_bdy);
    bool is_transformed = local_coordinate_transformation(cur_heh, hf_start, is_bdy);
    if (!is_transformed)
      return new_ehs;

    //initialize new cell quaternions with this quaternion
    Quaternion qt = cell_quaternions_[*mesh_.hec_iter(cur_heh)];

    int ff_id = 0;
    auto new_ft_hfs_vhs = store_new_feature_faces_and_index(cur_heh, ff_id);

    //store axis aligned to boundary normal in case the mesh geometry changes
    std::map<std::vector<VH>, int> hfv_aligned_axis;
    auto ft_hfs_vhs = collect_axis_aligned_to_feature_face(cur_heh, hfv_aligned_axis);

    //store axis aligned to feature edge
    std::map<EH, int> he_aligned_axis;
    collect_axis_aligned_to_feature(cur_heh, he_aligned_axis);

    std::vector<FH> new_ft_fhs;
    std::set<CH> new_chs;
    for (auto ii = 0u; ii < v_equatorial_vertices.size(); ++ii)
    {
      //do the swap
      if (v_equatorial_vertices[ii].size() == 3)
      {
        bool suc = check_feature_face_constraints_32(_eh, v_equatorial_vertices[ii], hfv_aligned_axis);
        if (suc)
        {
          swap_edge_32(_eh, v_equatorial_vertices[ii], new_chs);
        }
        else
        {
          std::cerr << "face alignment conflict 32: " << _eh << std::endl;
          return std::set<EH>();
        }
      }
      else if (v_equatorial_vertices[ii].size() == 4)
      {
        if (check_feature_face_constraints_44(_eh, v_equatorial_vertices[ii],
                                              v_equatorial_vertices[ii][v_ref_pos[ii]], hfv_aligned_axis))
        {
          new_ehs.insert(
                  swap_edge_44(_eh, v_equatorial_vertices[ii], v_equatorial_vertices[ii][v_ref_pos[ii]],
                               new_chs));
        }
        else
        {
          std::cerr << "face alignment conflict 44: " << _eh << std::endl;

          return std::set<EH>();
        }
      }
      else if (v_equatorial_vertices[ii].size() == 5)
      {
        if (check_feature_face_constraints_56(_eh, v_equatorial_vertices[ii], v_equatorial_vertices[ii][v_ref_pos[ii]],
                                              hfv_aligned_axis))
        {
          auto new_ehs_part = swap_edge_56(_eh, v_equatorial_vertices[ii], v_equatorial_vertices[ii][v_ref_pos[ii]],
                                           new_chs);
          new_ehs.insert(new_ehs_part.begin(), new_ehs_part.end());
        }
        else
        {
          std::cerr << "face alignment conflict 56: " << _eh << std::endl;

          return std::set<EH>();
        }
      }
    }

    //push new feature face edge
    if (feature_face_edge_[_eh])
    {
      auto ff_he = mesh_.find_halfedge(v_equatorial_vertices[0][0], v_equatorial_vertices[0].back());
      new_ehs.insert(mesh_.edge_handle(ff_he));
    }

    for (const auto &vhs: new_ft_hfs_vhs)
    {
      auto hfhi = mesh_.find_halfface(vhs);
      if (hfhi.is_valid())
        new_ft_fhs.push_back(mesh_.face_handle(hfhi));
      else
      {
        std::cerr << "Error: cannot find feature face with vertices! EH " << _eh << " is bdy " << is_bdy << " is ffe "
                  << is_ft << std::endl;
      }
    }
    //update properties
    update_properties(new_ft_fhs, ff_id, new_ehs, new_chs, qt, hfv_aligned_axis, he_aligned_axis);

    //check
//            for(const auto chi : new_chs) {
//                int nn = 0;
//                for(auto cf_it = mesh_.cf_iter(chi); cf_it.valid(); ++cf_it) {
//                    if(mesh_.is_boundary(*cf_it))
//                        nn++;
//                }
//
//                if(nn > 1) {
//                    std::cerr << "Error:check b: swap edge " << _eh << " results in " << nn << " bdy faces; case "
//                              << equatorial_vertices.size() << "   ";
//                    for(auto vhi : equatorial_vertices) {
//                        std::cerr<<" vh "<<vhi<<" bdy: "<<mesh_.is_boundary(vhi);
//                    }
//                    std::cerr<<std::endl;
//                }
//            }

//            for(const auto chi : new_chs) {
//                auto hfhs = mesh_.cell(chi).halffaces();
//                for(int i=0; i<4; ++i) {
//                    auto hfh_opp = mesh_.opposite_halfface_handle(hfhs[i]);
//                    if(mesh_.is_boundary(hfh_opp)) {
//                        auto nm = mesh_.normal(hfh_opp);
//                        Vec3d nmeg(nm[0], nm[1], nm[2]);
//
//                        // determine alignment axis
//                        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[chi], nmeg).second;
//                        Vec3d axis_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[chi], axis);
//
//                        if (nmeg.dot(axis_vec) < nmeg.norm() * 0.99985)
//                            std::cerr << "Warning FF alignment: face " << mesh_.face_handle(hfh_opp) << " normal deviates by "
//                                      << std::acos(nmeg.dot(axis_vec) / nmeg.norm()) * 180.0 / M_PI
//                                      << " degree... Happened in swapping edge "<<mesh_.edge(_eh).from_vertex()<<" "<<mesh_.edge(_eh).to_vertex()<<" is bdy? "<<mesh_.is_boundary(_eh)<<std::endl;
//                    }
//                }
//
//            }
  }

  return new_ehs;
}

template<class MeshT>
bool
EdgeSwapT<MeshT>::is_swap_ok(const EH _eh) const
{
  if (!_eh.is_valid() || mesh_.is_deleted(_eh))
    return false;

  if (feature_edge_[_eh] > 0 || valence_[_eh] != 0)
    return false;

  return true;
}

template<class MeshT>
int
EdgeSwapT<MeshT>::find_the_best_configuration(const EH _eh_sw, const std::vector<VH> &_eq_vhs,
                                              std::vector<std::vector<Point> > &_old_cells_points) const
{
  //index of the first vertex in eq vhs
  int ref_pos = -1;

  //number of different configurations
  int num_cnfgs = 0;
  if (_eq_vhs.size() == 3)
    num_cnfgs = 1;
  else if (_eq_vhs.size() == 4)
    num_cnfgs = 2;
  else if (_eq_vhs.size() == 5)
    num_cnfgs = 5;

  //compare energy of different configurations
  for (int i = 0; i < num_cnfgs; ++i)
  {
    std::vector<std::vector<VH> > new_cells_vertices;
    if (num_cnfgs == 1)
    {
      if (check_new_edges(_eq_vhs, i))
        new_cells_vertices = surviving_tets_32(_eh_sw, _eq_vhs);
    }
    else if (num_cnfgs == 2)
    {
      if (check_new_edges(_eq_vhs, i))
        new_cells_vertices = surviving_tets_44(_eh_sw, _eq_vhs, _eq_vhs[i]);
    }
    else if (num_cnfgs == 5)
    {
      if (check_new_edges(_eq_vhs, i))
        new_cells_vertices = surviving_tets_56(_eh_sw, _eq_vhs, _eq_vhs[i]);
    }

    if (new_cells_vertices.empty())
      continue;

    bool has_degenerated_cell = false;
    std::vector<std::vector<Point> > new_cells_points;
    new_cells_points.reserve(new_cells_vertices.size());

    for (const auto &cvhs: new_cells_vertices)
    {
      std::vector<Point> cell_points;
      cell_points.reserve(4);

      for (const auto &vh: cvhs)
        cell_points.push_back(mesh_.vertex(vh));

      if (!RemeshingAssist<MeshT>::is_tet_valid(cell_points))
      {
        has_degenerated_cell = true;
        break;
      }
      new_cells_points.push_back(cell_points);
    }

    if (has_degenerated_cell)
      continue;

    if (RemeshingAssist<MeshT>::is_quality_improved(new_cells_points, _old_cells_points, true))
    {
      _old_cells_points.swap(new_cells_points);
      ref_pos = i;
    }
  }

  return ref_pos;
}

template<class MeshT>
bool
EdgeSwapT<MeshT>::
check_new_edges(const std::vector<VH> &_eq_vhs, const int _cnf) const
{
  int eqsize = _eq_vhs.size();
  if (eqsize == 3)
  {
    for (int j = 0; j < eqsize; ++j)
    {
      auto hehj = mesh_.find_halfedge(_eq_vhs[j], _eq_vhs[(j + 1) % eqsize]);
      if (!hehj.is_valid())
        continue;
      //case 1: check if any regular edge connects to two feature face vertices
      if (!feature_face_edge_[mesh_.edge_handle(hehj)] &&
          (feature_face_vertex_[_eq_vhs[j]] && feature_face_vertex_[_eq_vhs[(j + 1) % eqsize]]))
        return false;
    }
  }
  else if (eqsize == 4)
  {
    //1.5 check if the new regular edge connects to two singular vertex
    if (post_remesh_)
    {
      if (sgl_vt_[_eq_vhs[_cnf]] && sgl_vt_[_eq_vhs[(_cnf + 2) % 4]])
        return false;
    }
    //1.5 check if the new regular edge connects to two feature face vertices
    if (feature_face_vertex_[_eq_vhs[_cnf]] && feature_face_vertex_[_eq_vhs[(_cnf + 2) % 4]])
      return false;
  }
  else if (eqsize == 5)
  {
    //check if the new regular edge connects to two singular vertex
    if (post_remesh_)
    {
      if ((sgl_vt_[_eq_vhs[_cnf]] && sgl_vt_[_eq_vhs[(_cnf + 2) % 5]]) ||
          (sgl_vt_[_eq_vhs[_cnf]] && sgl_vt_[_eq_vhs[(_cnf + 3) % 5]]))
        return false;
    }
    //check if the new regular edge connects to two feature face vertices
    if ((feature_face_vertex_[_eq_vhs[_cnf]] && feature_face_vertex_[_eq_vhs[(_cnf + 2) % 5]]) ||
        (feature_face_vertex_[_eq_vhs[_cnf]] && feature_face_vertex_[_eq_vhs[(_cnf + 3) % 5]]))
      return false;;
  }

  return true;
}

//make sure the transitions of the deleted halffaces are identity
template<class MeshT>
bool
EdgeSwapT<MeshT>::local_coordinate_transformation(const HEH _heh, const HFH _hf_start, const bool _is_bdy_e)
{
  auto hf_it = _hf_start;
  if (_is_bdy_e)
  {
    if (mesh_.is_boundary(hf_it))
      return false;
    hf_it = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_it, _heh));

    while (hf_it != _hf_start)
    {
      if (mesh_.is_boundary(hf_it))
        break;

      auto trans = trans_prop_[hf_it];
      auto trans_opp = trans_prop_[mesh_.opposite_halfface_handle(hf_it)];
      auto ch = mesh_.incident_cell(hf_it);
      auto hfs = mesh_.cell(ch).halffaces();

      for (int i = 0; i < 4; ++i)
      {
        if (mesh_.is_boundary(mesh_.face_handle(hfs[i])))
          continue;
        auto hfh_opp = mesh_.opposite_halfface_handle(hfs[i]);
        trans_prop_[hfs[i]] = tq_.mult_transitions_idx(trans_opp, trans_prop_[hfs[i]]);
        trans_prop_[hfh_opp] = tq_.inverse_transition_idx(trans_prop_[hfs[i]]);
      }

      // transform quaternion and matrices
      // q_t = q_ch^-1 * q_chopp
      // q_ch^-1 = q_t*q_chopp^-1
      // q_t^-1*q_t = Id
      // multiply q_t^-1 on both sides
      // q_ch -> q_ch*q_t
      cell_quaternions_[ch] = cell_quaternions_[ch] * tq_.transition(trans);

      //move to next
      hf_it = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_it, _heh));
    }
  }
  else
  {
    do
    {
      auto trans = trans_prop_[hf_it];
      auto trans_opp = trans_prop_[mesh_.opposite_halfface_handle(hf_it)];
      auto ch = mesh_.incident_cell(hf_it);
      auto hfs = mesh_.cell(ch).halffaces();

      for (int i = 0; i < 4; ++i)
      {
        if (mesh_.is_boundary(mesh_.face_handle(hfs[i])))
          continue;
        auto hfh_opp = mesh_.opposite_halfface_handle(hfs[i]);
        trans_prop_[hfs[i]] = tq_.mult_transitions_idx(trans_opp, trans_prop_[hfs[i]]);
        trans_prop_[hfh_opp] = tq_.inverse_transition_idx(trans_prop_[hfs[i]]);
      }

      // transform quaternion and matrices
      // q_t = q_ch^-1 * q_chopp
      // q_ch^-1 = q_t*q_chopp^-1
      // q_t^-1*q_t = Id
      // multiply q_t^-1 on both sides
      // q_ch -> q_ch*q_t
      cell_quaternions_[ch] = cell_quaternions_[ch] * tq_.transition(trans);

      //move to next
      hf_it = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_it, _heh));
    }
    while (hf_it != _hf_start);
  }

  return true;
}

template<class MeshT>
std::vector<std::vector<VH>>
EdgeSwapT<MeshT>::collect_axis_aligned_to_feature_face(const HEH _heh, std::map<std::vector<VH>, int> &hfv_aligned_axis)
{
  std::vector<std::vector<VH>> fthfs_vhs;
  if (!feature_face_edge_[mesh_.edge_handle(_heh)])
  {//feature face can only be cluster boundary
    for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
    {
      std::vector<HFH> fhfhs;
      for (auto chf_it = mesh_.chf_iter(*hec_it); chf_it.valid(); ++chf_it)
      {
        auto fhi = mesh_.face_handle(*chf_it);
        if (feature_fprop_[fhi] > 0)
        {
          fhfhs.push_back(*chf_it);
        }
      }

      if (fhfhs.size() == 1)
      {
        auto hf_vhs = mesh_.get_halfface_vertices(fhfhs[0]);

        auto nm = mesh_.normal(fhfhs[0]);
        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*hec_it], nm).second;

        hfv_aligned_axis.insert(std::make_pair(hf_vhs, axis));
        fthfs_vhs.push_back(hf_vhs);
      }
      else if (fhfhs.size() > 1)
      {
        std::cerr << "Warning: " << fhfhs.size() << " feature faces in a cell!" << std::endl;
      }
    }
  }
  else
  {//new feature faces will be created after swap
    std::vector<HFH> fhfhs;
    for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
      if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
      {
        fhfhs.push_back(*hehf_it);
      }
    }

    if (fhfhs.size() != 2)
    {
      std::cerr << "Error: non-manifold feature surface at edge " << mesh_.edge_handle(_heh)
                << " in edge swap!" << std::endl;
      return std::vector<std::vector<VH>>{};
    }

    //since there's no non-identity matching, we just pick one
    auto nm = mesh_.normal(fhfhs[0]);
    int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[mesh_.incident_cell(fhfhs[0])], nm).second;

    VH vh_f = mesh_.halfedge(_heh).from_vertex();
    VH vh_t = mesh_.halfedge(_heh).to_vertex();

    auto heh_fr = mesh_.next_halfedge_in_halfface(_heh, fhfhs[0]);
    auto heh_tl = mesh_.next_halfedge_in_halfface(_heh, fhfhs[1]);

    auto vh_r = mesh_.halfedge(heh_fr).to_vertex();
    auto vh_l = mesh_.halfedge(heh_tl).to_vertex();

    std::vector<VH> vhs0{vh_r, vh_f, vh_l}, vhs1{vh_r, vh_l, vh_t};
    hfv_aligned_axis.insert(std::make_pair(vhs0, axis));
    hfv_aligned_axis.insert(std::make_pair(vhs1, axis));

    fthfs_vhs.push_back(vhs0);
    fthfs_vhs.push_back(vhs1);

    if (!mesh_.is_boundary(_heh))
    {
      auto hfh_opp = mesh_.opposite_halfface_handle(fhfhs[0]);
      auto nm_opp = mesh_.normal(hfh_opp);
      int axis_opp = AxisAlignmentHelpers::closest_axis(
              cell_quaternions_[mesh_.incident_cell(hfh_opp)], nm_opp).second;

      int axis_ngt = axis % 2 == 0 ? axis + 1 : axis - 1;
      if (axis_opp != axis_ngt)
      {
        std::cerr << "Error: different aligned axis at feature face after coordinate transformation!" << std::endl;
      }
      std::vector<VH> vhs2{vh_f, vh_r, vh_l}, vhs3{vh_r, vh_t, vh_l};
      hfv_aligned_axis.insert(std::make_pair(vhs2, axis_opp));
      hfv_aligned_axis.insert(std::make_pair(vhs3, axis_opp));
      fthfs_vhs.push_back(vhs2);
      fthfs_vhs.push_back(vhs3);
    }
  }

  return fthfs_vhs;
}

template<class MeshT>
std::vector<std::vector<VH>>
EdgeSwapT<MeshT>::store_new_feature_faces_and_index(const HEH _heh, int &_ff_id) const
{
  _ff_id = 0;

  std::vector<std::vector<VH>> fthfs_vhs;
  if (feature_face_edge_[mesh_.edge_handle(_heh)])
  {//new feature faces will be created after swap
    std::vector<HFH> fhfhs;
    for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
      if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
      {
        fhfhs.push_back(*hehf_it);
      }
    }

    if (fhfhs.size() != 2)
    {
      std::cerr << "Error: non-manifold feature surface at edge " << mesh_.edge_handle(_heh)
                << " in edge swap!" << std::endl;
      return std::vector<std::vector<VH>>{};
    }

    VH vh_f = mesh_.halfedge(_heh).from_vertex();
    VH vh_t = mesh_.halfedge(_heh).to_vertex();

    auto heh_fr = mesh_.next_halfedge_in_halfface(_heh, fhfhs[0]);
    auto heh_tl = mesh_.next_halfedge_in_halfface(_heh, fhfhs[1]);

    auto vh_r = mesh_.halfedge(heh_fr).to_vertex();
    auto vh_l = mesh_.halfedge(heh_tl).to_vertex();

    std::vector<VH> vhs0{vh_r, vh_f, vh_l}, vhs1{vh_r, vh_l, vh_t};

    fthfs_vhs.push_back(vhs0);
    fthfs_vhs.push_back(vhs1);

    _ff_id = feature_fprop_[mesh_.face_handle(fhfhs[0])];
  }

  return fthfs_vhs;
}

template<class MeshT>
void
EdgeSwapT<MeshT>::collect_axis_aligned_to_feature(const HEH _heh, std::map<EH, int> &feh_aligned_axis)
{
  for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
  {
    for (auto ce_it = mesh_.ce_iter(*hec_it); ce_it.valid(); ++ce_it)
    {
      if (feature_edge_[*ce_it] > 0)
      {
        VH vhf = mesh_.edge(*ce_it).from_vertex();
        VH vht = mesh_.edge(*ce_it).to_vertex();
        auto dir = mesh_.vertex(vht) - mesh_.vertex(vhf);
        dir.normalize();

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*hec_it], dir).second;
        feh_aligned_axis[*ce_it] = axis;
      }
    }
  }
}

template<class MeshT>
bool EdgeSwapT<MeshT>::
has_confict_feature_face_constraints_in_cell(const std::vector<VH> &_cell_vhs,
                                             std::map<std::vector<VH>, int> &hfv_aligned_axis)
{
  //copy
  std::map<std::set<VH>, int> hfvs_axes;
  for (auto&[vhs, ax]: hfv_aligned_axis)
  {
    hfvs_axes.insert(std::make_pair(std::set<VH>{vhs[0], vhs[1], vhs[2]}, ax));
  }

  std::vector<int> num_cons(6, 0);

  //cell halffaces with alignment constraints
  auto it = hfvs_axes.find(std::set<VH>{_cell_vhs[0], _cell_vhs[1], _cell_vhs[2]});
  if (it != hfvs_axes.end())
  {
    int ax = it->second;
    if (ax < 6 && ax >= 0)
      num_cons[ax]++;
    else
      return true;
  }

  it = hfvs_axes.find(std::set<VH>{_cell_vhs[1], _cell_vhs[2], _cell_vhs[3]});
  if (it != hfvs_axes.end())
  {
    int ax = it->second;
    if (ax < 6 && ax >= 0)
      num_cons[ax]++;
    else
      return true;
  }

  it = hfvs_axes.find(std::set<VH>{_cell_vhs[2], _cell_vhs[3], _cell_vhs[0]});
  if (it != hfvs_axes.end())
  {
    int ax = it->second;
    if (ax < 6 && ax >= 0)
      num_cons[ax]++;
    else
      return true;
  }

  it = hfvs_axes.find(std::set<VH>{_cell_vhs[0], _cell_vhs[1], _cell_vhs[3]});
  if (it != hfvs_axes.end())
  {
    int ax = it->second;
    if (ax < 6 && ax >= 0)
      num_cons[ax]++;
    else
      return true;
  }

  for (int i = 0; i < 6; ++i)
    if (num_cons[i] > 1)
    {
      return true;
    }

  return false;
}

//
//              vh_t
//             / |  \
    //       vh_l  \ |  /  vh_r
//              vh_f
//
template<class MeshT>
bool
EdgeSwapT<MeshT>::is_boundary_edge_type_changed(const EH _eh) const
{
  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  auto heh1 = mesh_.halfedge_handle(_eh, 1);

  auto vh_f = mesh_.halfedge(heh0).from_vertex();
  auto vh_t = mesh_.halfedge(heh0).to_vertex();

  auto hfh_r = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh0));
  auto hfh_l = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh1));

  assert(mesh_.is_boundary(hfh_l) && mesh_.is_boundary(hfh_r));

  auto heh_fr = mesh_.next_halfedge_in_halfface(heh1, hfh_r);
  auto heh_tl = mesh_.next_halfedge_in_halfface(heh0, hfh_l);

  auto vh_r = mesh_.halfedge(heh_fr).to_vertex();
  auto vh_l = mesh_.halfedge(heh_tl).to_vertex();

  auto n_up = face_normal(mesh_, vh_l, vh_r, vh_t);
  auto n_down = face_normal(mesh_, vh_r, vh_l, vh_f);


  auto angle = dihedral_angle(mesh_.vertex(vh_l), mesh_.vertex(vh_r), n_up, n_down);
  auto old_angle = dihedral_angle(mesh_.vertex(vh_f), mesh_.vertex(vh_t), mesh_.normal(hfh_l), mesh_.normal(hfh_r));


  int old_type = RemeshingAssist<MeshT>::boundary_edge_type(old_angle);
  int new_type = RemeshingAssist<MeshT>::boundary_edge_type(angle);

  //the edge type changed
  if (old_type != new_type && new_type != 0)
  {
    return true;
  }


  std::map<HEH, Point> normal_map;
  normal_map[heh_fr] = n_down;
  normal_map[heh_tl] = n_up;


  std::vector<HEH> check_hehs;
  check_hehs.reserve(4);

  //check the other four edges (maybe not necessary)
  auto heh_rt = mesh_.next_halfedge_in_halfface(heh_fr, hfh_r);
  auto heh_lf = mesh_.next_halfedge_in_halfface(heh_tl, hfh_l);
  normal_map[heh_lf] = n_down;
  normal_map[heh_rt] = n_up;

  check_hehs.push_back(heh_fr);
  check_hehs.push_back(heh_tl);
  check_hehs.push_back(heh_lf);
  check_hehs.push_back(heh_rt);

  std::map<HEH, Point> normal_map_old;
  normal_map_old[heh_fr] = mesh_.normal(hfh_r);
  normal_map_old[heh_rt] = normal_map_old[heh_fr];

  normal_map_old[heh_tl] = mesh_.normal(hfh_l);
  normal_map_old[heh_lf] = normal_map_old[heh_tl];

  normal_map_old[mesh_.opposite_halfedge_handle(heh_fr)] = mesh_.normal(
          mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh_fr)));
  normal_map_old[mesh_.opposite_halfedge_handle(heh_rt)] = mesh_.normal(
          mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh_rt)));
  normal_map_old[mesh_.opposite_halfedge_handle(heh_tl)] = mesh_.normal(
          mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh_tl)));
  normal_map_old[mesh_.opposite_halfedge_handle(heh_lf)] = mesh_.normal(
          mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heh_lf)));

  for (const auto heh: check_hehs)
  {
    auto vh_f = mesh_.halfedge(heh).from_vertex();
    auto vh_t = mesh_.halfedge(heh).to_vertex();

    auto angle = dihedral_angle(mesh_.vertex(vh_f), mesh_.vertex(vh_t), normal_map[heh],
                                normal_map_old[mesh_.opposite_halfedge_handle(heh)]);
    auto old_angle = dihedral_angle(mesh_.vertex(vh_f), mesh_.vertex(vh_t), normal_map_old[heh],
                                    normal_map_old[mesh_.opposite_halfedge_handle(heh)]);


    int old_type = RemeshingAssist<MeshT>::boundary_edge_type(old_angle);
    int new_type = RemeshingAssist<MeshT>::boundary_edge_type(angle);

    //the edge type changed
    if (old_type != new_type && (new_type != 0 || valence_[mesh_.edge_handle(heh)] != 0))
    {
      return true;
    }
  }

  return false;
}


template<class MeshT>
void
EdgeSwapT<MeshT>::update_properties(const std::vector<FH> &_new_ft_fhs, const int _ff_id, const std::set<EH> &_new_ehs,
                                    const std::set<CH> &_new_chs, const Quaternion &_qt,
                                    std::map<std::vector<VH>, int> &_ff_aligned_axis,
                                    std::map<EH, int> &feh_aligned_axis)
{
  //update valence
  for (const auto ehi: _new_ehs)
    valence_[ehi] = 0;

  //update feature faces
  for (const auto fhi: _new_ft_fhs)
  {
    feature_fprop_[fhi] = _ff_id;

    for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
      feature_face_edge_[*fe_it] = true;

    for (auto fv_it = mesh_.fv_iter(fhi); fv_it.valid(); ++fv_it)
      feature_face_vertex_[*fv_it] = true;
  }

  //update transitions for new faces (swap_32)
//        trans_prop_[_new_hfh] = 0;
//        trans_prop_[mesh_.opposite_halfface_handle(_new_hfh)] = 0;

  //other cases
  for (const auto ehi: _new_ehs)
  {
    for (auto hehf_it = mesh_.hehf_iter(mesh_.halfedge_handle(ehi, 0)); hehf_it.valid(); ++hehf_it)
      if (!mesh_.is_boundary(mesh_.face_handle(*hehf_it)))
      {
        trans_prop_[*hehf_it] = 0;
        trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = 0;
      }
      else
      {
        trans_prop_[*hehf_it] = -1;
        trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = -1;
      }
  }

  //initialize
  for (const auto chi: _new_chs)
  {
    cell_quaternions_[chi] = _qt;
  }

  //smooth cell quaternions
  std::set<FH> bdy_fhs;
  for (const auto ch: _new_chs)
    for (auto cf_it = mesh_.cf_iter(ch); cf_it.valid(); ++cf_it)
      if (mesh_.is_boundary(*cf_it))
        bdy_fhs.insert(*cf_it);

  //feature face aligned axis
  std::map<CH, std::vector<std::pair<Vec3d, int>>> alignment;

  for (auto&[vhs, aid]: _ff_aligned_axis)
  {
    auto hfh = mesh_.find_halfface(vhs);
    if (!hfh.is_valid())
    {
      std::cerr << "vhs: " << vhs[0] << " " << vhs[1] << " " << vhs[2] << " " << aid << std::endl;
      std::cerr << "Error: no halfface corresponds to vertices" << std::endl;
//                continue;
    }

    auto nm = mesh_.normal(hfh);
    auto ch = mesh_.incident_cell(hfh);

    alignment[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), aid);
  }

  //check
  for (const auto ch: _new_chs)
  {
    if (alignment[ch].size() >= 2)
    {
      std::cerr << "Warning: cell " << ch << " has " << alignment[ch].size() << " feature face constraints"
                << std::endl;
      for (auto[nm, ax]: alignment[ch])
      {
        std::cerr << " dir " << nm[0] << " " << nm[1] << " " << nm[2] << " ax " << ax;
      }
      std::cerr << std::endl;
    }
  }

  //feature aligned axis
  for (const auto ch: _new_chs)
  {
    for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
    {
      if (feature_edge_[*ce_it] > 0)
      {
        int ax = feh_aligned_axis[*ce_it];
        auto dir = mesh_.vertex(mesh_.edge(*ce_it).to_vertex()) - mesh_.vertex(mesh_.edge(*ce_it).from_vertex());
        dir.normalize();

        alignment[ch].emplace_back(Vec3d(dir[0], dir[1], dir[2]), ax);
      }
    }
  }


  for (const auto ch: _new_chs)
  {
    if (alignment[ch].size() > 2)
    {
      std::cerr << "Warning: " << alignment[ch].size() << " constraints in cell in edge swap " << ch << std::endl;
      for (auto[nm, ax]: alignment[ch])
      {
        std::cerr << " dir " << nm[0] << " " << nm[1] << " " << nm[2] << " ax " << ax;
      }
      std::cerr << std::endl;
    }
  }

  for (int i = 0; i < 10; ++i)
    QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, _new_chs, bdy_fhs, alignment, tq_,
                                                 cell_quaternions_);
}


//type: three tets cluster turns to two tets cluster after edge swap
template<class MeshT>
HFH
EdgeSwapT<MeshT>::swap_edge_32(const EH _eh, const std::vector<VH> &_equatorial_vertices, std::set<CH> &_new_chs)
{
  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  auto from_vh = mesh_.halfedge(heh0).from_vertex();
  auto to_vh = mesh_.halfedge(heh0).to_vertex();

  //delete edge (at interior feature face, if not deleted)
  if (!mesh_.is_deleted(_eh))
    mesh_.delete_edge(_eh);

  //add new cells
  auto ch0 = mesh_.add_cell(_equatorial_vertices[0], _equatorial_vertices[1], _equatorial_vertices[2], to_vh);
  auto ch1 = mesh_.add_cell(_equatorial_vertices[0], _equatorial_vertices[2], _equatorial_vertices[1], from_vh);

  _new_chs.insert(ch0);
  _new_chs.insert(ch1);


  return mesh_.find_halfface(_equatorial_vertices);
}

template<class MeshT>
bool
EdgeSwapT<MeshT>::check_feature_face_constraints_32(const EH _eh, const std::vector<VH> &_equatorial_vertices,
                                                    std::map<std::vector<VH>, int> &hfv_aligned_axis)
{
  if (!hfv_aligned_axis.empty())
  {
    auto heh0 = mesh_.halfedge_handle(_eh, 0);
    auto from_vh = mesh_.halfedge(heh0).from_vertex();
    auto to_vh = mesh_.halfedge(heh0).to_vertex();

    bool has_conf = has_confict_feature_face_constraints_in_cell(
            std::vector<VH>{_equatorial_vertices[0], _equatorial_vertices[1], _equatorial_vertices[2], to_vh},
            hfv_aligned_axis);

    if (has_conf)
      return false;
    else
    {
      has_conf = has_confict_feature_face_constraints_in_cell(
              std::vector<VH>{_equatorial_vertices[0], _equatorial_vertices[1], _equatorial_vertices[2], from_vh},
              hfv_aligned_axis);
      return !has_conf;
    }
  }

  return true;
}


//type: four tets cluster turns to four tets cluster after edge swap(2 different configurations)
template<class MeshT>
EH
EdgeSwapT<MeshT>::swap_edge_44(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh,
                               std::set<CH> &_new_chs)
{
  auto it = std::find(_equatorial_vertices.begin(), _equatorial_vertices.end(), _ref_vh);
  if (it == _equatorial_vertices.end())
    return EH(-1);

  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  auto from_vh = mesh_.halfedge(heh0).from_vertex();
  auto to_vh = mesh_.halfedge(heh0).to_vertex();

//        //delete old cells
//        for(auto hec_it = mesh_.hec_iter(heh0); hec_it.valid(); ++hec_it)
//            mesh_.delete_cell(*hec_it);

  //delete edge (at interior feature face, if not deleted)
  if (!mesh_.is_deleted(_eh))
    mesh_.delete_edge(_eh);

  //add new cells
  EH new_eh(-1);
  CH new_ch;
  if (_ref_vh == _equatorial_vertices[0] || _ref_vh == _equatorial_vertices[2])
  {
    new_ch = mesh_.add_cell(_equatorial_vertices[0], _equatorial_vertices[1],
                            _equatorial_vertices[2], to_vh);
    _new_chs.insert(new_ch);

    new_ch = mesh_.add_cell(_equatorial_vertices[0], _equatorial_vertices[2],
                            _equatorial_vertices[1], from_vh);
    _new_chs.insert(new_ch);

    new_ch = mesh_.add_cell(_equatorial_vertices[2], _equatorial_vertices[3],
                            _equatorial_vertices[0], to_vh);
    _new_chs.insert(new_ch);

    new_ch = mesh_.add_cell(_equatorial_vertices[2], _equatorial_vertices[0],
                            _equatorial_vertices[3], from_vh);
    _new_chs.insert(new_ch);

    new_eh = mesh_.edge_handle(mesh_.find_halfedge(_equatorial_vertices[0], _equatorial_vertices[2]));
  }
  else if (_ref_vh == _equatorial_vertices[1] || _ref_vh == _equatorial_vertices[3])
  {
    new_ch = mesh_.add_cell(_equatorial_vertices[1], _equatorial_vertices[2],
                            _equatorial_vertices[3], to_vh);
    _new_chs.insert(new_ch);

    new_ch = mesh_.add_cell(_equatorial_vertices[1], _equatorial_vertices[3],
                            _equatorial_vertices[2], from_vh);
    _new_chs.insert(new_ch);

    new_ch = mesh_.add_cell(_equatorial_vertices[3], _equatorial_vertices[0],
                            _equatorial_vertices[1], to_vh);
    _new_chs.insert(new_ch);

    new_ch = mesh_.add_cell(_equatorial_vertices[3], _equatorial_vertices[1],
                            _equatorial_vertices[0], from_vh);
    _new_chs.insert(new_ch);

//            mesh_.halfedge(_equatorial_vertices[1], _equatorial_vertices[3]);
    new_eh = mesh_.edge_handle(mesh_.find_halfedge(_equatorial_vertices[1], _equatorial_vertices[3]));
  }


  return new_eh;
}

template<class MeshT>
bool
EdgeSwapT<MeshT>::check_feature_face_constraints_44(const EH _eh, const std::vector<VH> &_equatorial_vertices,
                                                    const VH _ref_vh, std::map<std::vector<VH>, int> &hfv_aligned_axis)
{
  if (!hfv_aligned_axis.empty())
  {
    auto it = std::find(_equatorial_vertices.begin(), _equatorial_vertices.end(), _ref_vh);
    if (it == _equatorial_vertices.end())
      return false;

    auto heh0 = mesh_.halfedge_handle(_eh, 0);
    auto from_vh = mesh_.halfedge(heh0).from_vertex();
    auto to_vh = mesh_.halfedge(heh0).to_vertex();

    if (_ref_vh == _equatorial_vertices[0] || _ref_vh == _equatorial_vertices[2])
    {
      bool has_conf = has_confict_feature_face_constraints_in_cell(
              std::vector<VH>{_equatorial_vertices[0], _equatorial_vertices[1],
                              _equatorial_vertices[2], to_vh}, hfv_aligned_axis);

      if (has_conf)
        return false;
      else
      {
        has_conf = has_confict_feature_face_constraints_in_cell(
                std::vector<VH>{_equatorial_vertices[0], _equatorial_vertices[2], _equatorial_vertices[1], from_vh},
                hfv_aligned_axis);

        if (has_conf)
          return false;
        else
        {
          has_conf = has_confict_feature_face_constraints_in_cell(
                  std::vector<VH>{_equatorial_vertices[2], _equatorial_vertices[3], _equatorial_vertices[0], to_vh},
                  hfv_aligned_axis);

          if (has_conf)
            return false;
          else
          {
            has_conf = has_confict_feature_face_constraints_in_cell(
                    std::vector<VH>{_equatorial_vertices[2], _equatorial_vertices[0], _equatorial_vertices[3], from_vh},
                    hfv_aligned_axis);

            return !has_conf;
          }
        }
      }
    }
    else if (_ref_vh == _equatorial_vertices[1] || _ref_vh == _equatorial_vertices[3])
    {
      bool has_conf = has_confict_feature_face_constraints_in_cell(
              std::vector<VH>{_equatorial_vertices[1], _equatorial_vertices[2],
                              _equatorial_vertices[3], to_vh}, hfv_aligned_axis);

      if (has_conf)
        return false;
      else
      {
        has_conf = has_confict_feature_face_constraints_in_cell(
                std::vector<VH>{_equatorial_vertices[1], _equatorial_vertices[3], _equatorial_vertices[2], from_vh},
                hfv_aligned_axis);

        if (has_conf)
          return false;
        else
        {
          has_conf = has_confict_feature_face_constraints_in_cell(
                  std::vector<VH>{_equatorial_vertices[3], _equatorial_vertices[0], _equatorial_vertices[1], to_vh},
                  hfv_aligned_axis);

          if (has_conf)
            return false;
          else
          {
            has_conf = has_confict_feature_face_constraints_in_cell(
                    std::vector<VH>{_equatorial_vertices[3], _equatorial_vertices[1], _equatorial_vertices[0], from_vh},
                    hfv_aligned_axis);

            return !has_conf;
          }
        }
      }
    }
  }

  return true;
}

//type: five tets cluster turns to six tets cluster after edge swap(5 different configurations)
template<class MeshT>
std::vector<EH>
EdgeSwapT<MeshT>::swap_edge_56(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh,
                               std::set<CH> &_new_chs)
{
  std::vector<EH> new_ehs;
  int pos = -1;
  for (int i = 0; i < 5; ++i)
    if (_ref_vh == _equatorial_vertices[i])
      pos = i;
  if (pos == -1)
    return new_ehs;

  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  auto from_vh = mesh_.halfedge(heh0).from_vertex();
  auto to_vh = mesh_.halfedge(heh0).to_vertex();

  //delete edge (at interior feature face, if not deleted)
  if (!mesh_.is_deleted(_eh))
    mesh_.delete_edge(_eh);

  //add new cells
  CH new_ch;
  new_ch = mesh_.add_cell(_equatorial_vertices[pos], _equatorial_vertices[(pos + 1) % 5],
                          _equatorial_vertices[(pos + 2) % 5], to_vh);
  _new_chs.insert(new_ch);

  new_ch = mesh_.add_cell(_equatorial_vertices[pos], _equatorial_vertices[(pos + 2) % 5],
                          _equatorial_vertices[(pos + 1) % 5], from_vh);
  _new_chs.insert(new_ch);

  new_ch = mesh_.add_cell(_equatorial_vertices[pos], _equatorial_vertices[(pos + 2) % 5],
                          _equatorial_vertices[(pos + 3) % 5], to_vh);
  _new_chs.insert(new_ch);

  new_ch = mesh_.add_cell(_equatorial_vertices[pos], _equatorial_vertices[(pos + 3) % 5],
                          _equatorial_vertices[(pos + 2) % 5], from_vh);
  _new_chs.insert(new_ch);

  new_ch = mesh_.add_cell(_equatorial_vertices[pos], _equatorial_vertices[(pos + 3) % 5],
                          _equatorial_vertices[(pos + 4) % 5], to_vh);
  _new_chs.insert(new_ch);

  new_ch = mesh_.add_cell(_equatorial_vertices[pos], _equatorial_vertices[(pos + 4) % 5],
                          _equatorial_vertices[(pos + 3) % 5], from_vh);
  _new_chs.insert(new_ch);

  new_ehs.push_back(mesh_.edge_handle(mesh_.find_halfedge(_equatorial_vertices[pos],
                                                          _equatorial_vertices[(pos + 2) % 5])));
  new_ehs.push_back(mesh_.edge_handle(mesh_.find_halfedge(_equatorial_vertices[pos],
                                                          _equatorial_vertices[(pos + 3) % 5])));


  return new_ehs;
}

template<class MeshT>
bool
EdgeSwapT<MeshT>::check_feature_face_constraints_56(const EH _eh, const std::vector<VH> &_equatorial_vertices,
                                                    const VH _ref_vh, std::map<std::vector<VH>, int> &hfv_aligned_axis)
{
  if (!hfv_aligned_axis.empty())
  {
    int pos = -1;
    for (int i = 0; i < 5; ++i)
      if (_ref_vh == _equatorial_vertices[i])
        pos = i;
    if (pos == -1)
      return false;

    auto heh0 = mesh_.halfedge_handle(_eh, 0);
    auto from_vh = mesh_.halfedge(heh0).from_vertex();
    auto to_vh = mesh_.halfedge(heh0).to_vertex();

    bool has_conf = has_confict_feature_face_constraints_in_cell(
            std::vector<VH>{_equatorial_vertices[pos], _equatorial_vertices[(pos + 1) % 5],
                            _equatorial_vertices[(pos + 2) % 5], to_vh}, hfv_aligned_axis);

    if (has_conf)
      return false;
    else
    {
      has_conf = has_confict_feature_face_constraints_in_cell(
              std::vector<VH>{_equatorial_vertices[pos], _equatorial_vertices[(pos + 2) % 5],
                              _equatorial_vertices[(pos + 1) % 5], from_vh}, hfv_aligned_axis);

      if (has_conf)
        return false;
      else
      {
        has_conf = has_confict_feature_face_constraints_in_cell(
                std::vector<VH>{_equatorial_vertices[pos], _equatorial_vertices[(pos + 2) % 5],
                                _equatorial_vertices[(pos + 3) % 5], to_vh}, hfv_aligned_axis);

        if (has_conf)
          return false;
        else
        {
          has_conf = has_confict_feature_face_constraints_in_cell(
                  std::vector<VH>{_equatorial_vertices[pos], _equatorial_vertices[(pos + 3) % 5],
                                  _equatorial_vertices[(pos + 2) % 5], from_vh}, hfv_aligned_axis);

          if (has_conf)
            return false;
          else
          {
            has_conf = has_confict_feature_face_constraints_in_cell(
                    std::vector<VH>{_equatorial_vertices[pos], _equatorial_vertices[(pos + 3) % 5],
                                    _equatorial_vertices[(pos + 4) % 5], to_vh}, hfv_aligned_axis);

            if (has_conf)
              return false;
            else
            {
              has_conf = has_confict_feature_face_constraints_in_cell(
                      std::vector<VH>{_equatorial_vertices[pos], _equatorial_vertices[(pos + 4) % 5],
                                      _equatorial_vertices[(pos + 3) % 5], from_vh},
                      hfv_aligned_axis);

              return !has_conf;
            }
          }
        }
      }
    }
  }

  return true;
}

//ccw
template<class MeshT>
std::vector<std::vector<VH>>
EdgeSwapT<MeshT>::get_equatorial_vertices(const EH _eh) const
{
  auto heh0 = mesh_.halfedge_handle(_eh, 0);

  std::vector<std::vector<VH>> v_equatorial_vhs(2);
  //try to get feature faces if it is a feature face edge
  HFH hfh0(-1);
  if (feature_face_edge_[_eh])
  {
    for (auto hehf_it = mesh_.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
      if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
      {
        hfh0 = *hehf_it;
        break;
      }
  }
  else
    hfh0 = *mesh_.hehf_iter(heh0);

  auto hf_it = hfh0;

  do
  {
    v_equatorial_vhs[0].push_back(mesh_.halfedge(mesh_.next_halfedge_in_halfface(heh0, hf_it)).to_vertex());
    if (v_equatorial_vhs[0].size() == 1000)
    {
      std::cout << "Error: broken mesh1! Something bad happened in the remeshing!" << std::endl;
      break;
    }

    hf_it = mesh_.adjacent_halfface_in_cell(hf_it, heh0);
    if (!hf_it.is_valid())
    {
      std::cout << "Error: adjacent halfface in cell is invalid! Error in swapping edge " << _eh
                << "! feature face edge " << feature_face_edge_[_eh] << " bdy " << mesh_.is_boundary(_eh) << std::endl;

      for (auto ec_it = mesh_.ec_iter(_eh); ec_it.valid(); ++ec_it)
        if (!mesh_.is_valid(*ec_it) || mesh_.is_deleted(*ec_it))
          std::cout << "Error: adjacent cell of edge " << _eh << " is invalid!" << std::endl;

      for (auto ef_it = mesh_.ef_iter(_eh); ef_it.valid(); ++ef_it)
        if (!mesh_.is_valid(*ef_it) || mesh_.is_deleted(*ef_it))
          std::cout << "Error: adjacent face of edge " << _eh << " is invalid!" << std::endl;

      int n_ff = 0;
      for (auto ef_it = mesh_.ef_iter(_eh); ef_it.valid(); ++ef_it)
        if (mesh_.is_valid(*ef_it) && !mesh_.is_deleted(*ef_it))
          n_ff++;

      std::cout << n_ff << " adjacent feature face of edge " << _eh << "!" << std::endl;

      break;
    }
    hf_it = mesh_.opposite_halfface_handle(hf_it);

    if (feature_fprop_[mesh_.face_handle(hf_it)] > 0)
    {
      v_equatorial_vhs[0].push_back(mesh_.halfedge(mesh_.next_halfedge_in_halfface(heh0, hf_it)).to_vertex());
      break;
    }
  }
  while (hf_it != hfh0);


  if (feature_face_edge_[_eh] && !mesh_.is_boundary(hf_it))
  {
    do
    {
      v_equatorial_vhs[1].push_back(mesh_.halfedge(mesh_.next_halfedge_in_halfface(heh0, hf_it)).to_vertex());
      if (v_equatorial_vhs[1].size() == 1000)
      {
        std::cout << "Error: broken mesh2! Something bad happened in the remeshing!" << std::endl;
        break;
      }

      hf_it = mesh_.adjacent_halfface_in_cell(hf_it, heh0);
      if (!hf_it.is_valid())
      {
        std::cout << "Error: adjacent halfface in cell is invalid2! Something bad happened in the remeshing!"
                  << std::endl;
        break;
      }
      hf_it = mesh_.opposite_halfface_handle(hf_it);

      if (feature_fprop_[mesh_.face_handle(hf_it)] > 0)
      {
        v_equatorial_vhs[1].push_back(mesh_.halfedge(mesh_.next_halfedge_in_halfface(heh0, hf_it)).to_vertex());
        break;
      }
    }
    while (hf_it != hfh0);
  }

  return v_equatorial_vhs;
}

//get the new tets after edge swap
template<class MeshT>
std::vector<std::vector<VH> >
EdgeSwapT<MeshT>::surviving_tets_32(const EH _eh, const std::vector<VH> &_equatorial_vertices) const
{
  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  auto from_vh = mesh_.halfedge(heh0).from_vertex();
  auto to_vh = mesh_.halfedge(heh0).to_vertex();

  std::vector<std::vector<VH> > new_cells_vertices;
  new_cells_vertices.reserve(2);

  new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[0], _equatorial_vertices[1],
                                               _equatorial_vertices[2], to_vh});

  new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[0], _equatorial_vertices[2],
                                               _equatorial_vertices[1], from_vh});

  return new_cells_vertices;
}

template<class MeshT>
std::vector<std::vector<VH> >
EdgeSwapT<MeshT>::surviving_tets_44(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh) const
{
  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  auto from_vh = mesh_.halfedge(heh0).from_vertex();
  auto to_vh = mesh_.halfedge(heh0).to_vertex();

  std::vector<std::vector<VH> > new_cells_vertices;
  new_cells_vertices.reserve(4);

  if (_ref_vh == _equatorial_vertices[0] || _ref_vh == _equatorial_vertices[2])
  {
    new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[0],
                                                 _equatorial_vertices[1], _equatorial_vertices[2], to_vh});
    new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[0], _equatorial_vertices[2],
                                                 _equatorial_vertices[1], from_vh});
    new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[2], _equatorial_vertices[3],
                                                 _equatorial_vertices[0], to_vh});
    new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[2], _equatorial_vertices[0],
                                                 _equatorial_vertices[3], from_vh});
  }
  else if (_ref_vh == _equatorial_vertices[1] || _ref_vh == _equatorial_vertices[3])
  {
    new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[1], _equatorial_vertices[2],
                                                 _equatorial_vertices[3], to_vh});
    new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[1], _equatorial_vertices[3],
                                                 _equatorial_vertices[2], from_vh});
    new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[3], _equatorial_vertices[0],
                                                 _equatorial_vertices[1], to_vh});
    new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[3], _equatorial_vertices[1],
                                                 _equatorial_vertices[0], from_vh});
  }

  return new_cells_vertices;
}

template<class MeshT>
std::vector<std::vector<VH> >
EdgeSwapT<MeshT>::surviving_tets_56(const EH _eh, const std::vector<VH> &_equatorial_vertices, const VH _ref_vh) const
{
  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  auto from_vh = mesh_.halfedge(heh0).from_vertex();
  auto to_vh = mesh_.halfedge(heh0).to_vertex();

  std::vector<std::vector<VH> > new_cells_vertices;
  int pos = -1;
  for (int i = 0; i < 5; ++i)
    if (_ref_vh == _equatorial_vertices[i])
      pos = i;
  if (pos == -1)
    return new_cells_vertices;

  new_cells_vertices.reserve(6);

  new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[(pos)],
                                               _equatorial_vertices[(pos + 1) % 5], _equatorial_vertices[(pos + 2) % 5],
                                               to_vh});
  new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[(pos)],
                                               _equatorial_vertices[(pos + 2) % 5], _equatorial_vertices[(pos + 1) % 5],
                                               from_vh});
  new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[(pos)],
                                               _equatorial_vertices[(pos + 2) % 5], _equatorial_vertices[(pos + 3) % 5],
                                               to_vh});
  new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[(pos)],
                                               _equatorial_vertices[(pos + 3) % 5], _equatorial_vertices[(pos + 2) % 5],
                                               from_vh});
  new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[(pos)],
                                               _equatorial_vertices[(pos + 3) % 5], _equatorial_vertices[(pos + 4) % 5],
                                               to_vh});
  new_cells_vertices.push_back(std::vector<VH>{_equatorial_vertices[(pos)],
                                               _equatorial_vertices[(pos + 4) % 5], _equatorial_vertices[(pos + 3) % 5],
                                               from_vh});

  return new_cells_vertices;
}

template<class MeshT>
HFH
EdgeSwapT<MeshT>::
find_start_halfface(const HEH _heh, const bool _is_bdy_e) const
{
  HFH hfh = *mesh_.hehf_iter(_heh);
  if (_is_bdy_e)
  {
    for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
      auto ch = mesh_.incident_cell(*hehf_it);
      if (ch.is_valid())
      {
        if (mesh_.is_boundary(ch))
        {
          hfh = *hehf_it;
          break;
        }
      }
    }
  }

  return hfh;
}

}

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define ENSUREEDGEMESHABILITYT_C

#include "EnsureEdgeMeshabilityT.hh"
#include "EdgeMonodromyHelperT.hh"
#include "FieldAngleCalculatorT.hh"
#include "../AxisAlignment.hh"

namespace AlgoHex
{
template<class MeshT>
bool
EnsureEdgeMeshabilityT<MeshT>::
fix_interior_compound_singular_edge(const HEH _heh, const bool _split)
{
  VH vht = mesh_.halfedge(_heh).to_vertex();

  SplitHelperT<MeshT>::split_one_ring(mesh_, es_, vht, valence_, feature_face_edge_, feature_face_vertex_,
                                      feature_edge_, sgl_vt_,
                                      feature_edge_vertex_);


  EH eh = mesh_.edge_handle(_heh);

  HFH hf_best(-1);
  int trans_best = -1;
  if (valence_[eh] == std::numeric_limits<int>::max())
  {
    double smoothest = 100000.;
    for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
      HEH heh_next = mesh_.next_halfedge_in_halfface(_heh, *hehf_it);
      VH vh_next = mesh_.halfedge(heh_next).to_vertex();
      if (!sgl_vt_[vh_next])
      {
        //change the matching of it
        //t = tn ... t1 t0

        //singular edge transition index
        int e_trans = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, _heh, *hehf_it);

        //store new candidate transition indices of the halfface
        std::vector<int> cd_vtrans;
        for (int id = 4; id < 10; ++id)
        {
          int e_trans_new = tq_.mult_transitions_idx(e_trans, id);
          if (e_trans_new < 10 && e_trans_new > 3)
            cd_vtrans.push_back(tq_.mult_transitions_idx(trans_prop_[*hehf_it], id));
        }

        //if complex singular edge cannot be split into two singular edges
        if (cd_vtrans.empty())
        {
          for (int id = 4; id < 10; ++id)
          {
            cd_vtrans.push_back(tq_.mult_transitions_idx(trans_prop_[*hehf_it], id));
          }
        }

        auto bt = best_transition_idx(_heh, *hehf_it, cd_vtrans);

        if (bt.first < smoothest && bt.first >= 0.)
        {
          smoothest = bt.first;
          hf_best = *hehf_it;
          trans_best = bt.second;
        }
      }
    }
  }
  else if (std::abs(valence_[eh]) >= 2)
  {
    for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
      HEH heh_next = mesh_.next_halfedge_in_halfface(_heh, *hehf_it);
      VH vh_next = mesh_.halfedge(heh_next).to_vertex();
      if (!sgl_vt_[vh_next])
      {
        hf_best = *hehf_it;
        break;
      }
    }

    std::map<CH, Quaternion> old_cell_quaternions;
    for (auto vc_it = mesh_.vc_iter(vht); vc_it.valid(); ++vc_it)
      old_cell_quaternions[*vc_it] = cell_quaternions_[*vc_it];

    //align quaternions to singular or feature edges or feature faces
    QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                        feature_fprop_, feature_edge_, vht, true);

    int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                     valence_, _heh, hf_best);

    //calculate real valence with aligned field
    int val = sge_.calculate_edge_valence(eh);

    //safe check
    bool is_ft_node5 = false;
    if (feature_node_[vht])
    {
      if (n_incident_feature_edges(mesh_, feature_edge_, vht) >= 5)
        is_ft_node5 = true;
    }
    if (is_ft_node5 && val == -2)
    {
      std::cerr << "Warning: at feature node the valence of singular edge is -2" << std::endl;
      val = 2;
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "old valence " << val << "   ddddd ";)
    //store new candidate transition indices of the halfface
    int old_trans = trans_prop_[hf_best];
    std::vector<int> vtrans{rt_axis + 4, negate(rt_axis) + 4};
    for (int j = 0; j < 2; ++j)
    {
      trans_prop_[hf_best] = tq_.mult_transitions_idx(trans_prop_[hf_best], vtrans[j]);
      trans_prop_[mesh_.opposite_halfface_handle(hf_best)] = tq_.inverse_transition_idx(trans_prop_[hf_best]);

      int val_new = sge_.calculate_edge_valence(eh);
      ALGOHEX_DEBUG_ONLY(std::cerr << " new valence " << val_new << "   ddddd ";)

      if (std::abs(val_new) == std::abs(val) - 1)
      {
        trans_best = trans_prop_[hf_best];
        ALGOHEX_DEBUG_ONLY(std::cerr << " found " << trans_best << "   ddddd ";)

        for (auto vc_it = mesh_.vc_iter(vht); vc_it.valid(); ++vc_it)
          cell_quaternions_[*vc_it] = old_cell_quaternions[*vc_it];

        break;
      }
      //restore matching
      trans_prop_[hf_best] = old_trans;
      trans_prop_[mesh_.opposite_halfface_handle(hf_best)] = tq_.inverse_transition_idx(trans_prop_[hf_best]);
    }

    if (trans_best == -1)
    {
      std::cerr << "Error: couldn't find a matching to untangle singular edge of valence +/-2" << std::endl;
      return false;
    }
  }


  if (hf_best.is_valid())
  {
    if (_split)
    {
      VH vh_f = mesh_.halfedge(_heh).from_vertex();
      VH vh_m = fs_.face_split(mesh_.face_handle(hf_best));
      std::vector<VH> hf_vhs{vh_f, vht, vh_m};
      hf_best = mesh_.find_halfface(hf_vhs);
    }

    trans_prop_[hf_best] = trans_best;
    trans_prop_[mesh_.opposite_halfface_handle(hf_best)] = tq_.inverse_transition_idx(trans_prop_[hf_best]);

    //update edge valence
    for (auto hfe_it = mesh_.hfe_iter(hf_best); hfe_it.valid(); ++hfe_it)
      sge_.compute_edge_valence(*hfe_it);

    //update singular vertex property
    for (auto hfv_it = mesh_.hfv_iter(hf_best); hfv_it.valid(); ++hfv_it)
      MeshPropertiesT<MeshT>::update_singular_vertex_property(*hfv_it);

    QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                        feature_fprop_, feature_edge_, vht);

    return true;
  }

  return false;
}

template<class MeshT>
bool
EnsureEdgeMeshabilityT<MeshT>::
fix_edge_with_non_meshable_footprint(const EH _eh)
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
    return false;
  }

  HFH hf_start = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(ft_hfhs[0], heh0));
  if (feature_fprop_[mesh_.face_handle(hf_start)] > 0)
  {
    std::cerr << "Error: no DOF for fixing boundary complex edge!" << std::endl;
    return false;
  }

  HFH hf_s = ft_hfhs[0], hf_e = ft_hfhs[1];
  auto ch_s = mesh_.incident_cell(ft_hfhs[0]);

  int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   heh0, hf_start);
  int e_axis = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, heh0, ch_s, rt_axis);
  auto e_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_s], (AxisAlignment) e_axis);
  e_dir.normalize();
  auto nm = mesh_.normal(ft_hfhs[0]);

  bool fixed = false;
  if (abs(e_dir.dot(ovm2eigen(nm))) > 0.98)
  {//singular edge is orthogonal to surface
    ALGOHEX_DEBUG_ONLY(
            std::cerr << "High/low valence " << valence_[_eh] << " singular edge is orthogonal to feature surface"
                      << std::endl;)
    //change the matching
    int mult_trans = rt_axis / 2 + 1;

    int newtrans0 = tq_.mult_transitions_idx(trans_prop_[hf_start], tq_.inverse_transition_idx(mult_trans));
    fixed = attempt_matching_adjustment_of_feature_face_sector(heh0, hf_s, hf_e, newtrans0);
  }
  else //singular edge is tangent to surface
  {
    bool increase = false;
    if (!mesh_.is_boundary(heh0))
    {
      //coordinate transformation, making matchings on feature face identity
      ArcZippingT<MeshT>::coordinate_transformation_in_cell(ft_hfhs[0]);
      ArcZippingT<MeshT>::coordinate_transformation_in_cell(ft_hfhs[1]);

      int fsa0 = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                     trans_prop_, heh0, ft_hfhs[0],
                                                                                     ft_hfhs[1]);
      auto da0 = dihedral_angle_from_halfface_to_halfface(mesh_, heh0, ft_hfhs[0], ft_hfhs[1]);
      int fsa1 = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                     trans_prop_, heh0, ft_hfhs[1],
                                                                                     ft_hfhs[0]);
      ALGOHEX_DEBUG_ONLY(std::cerr << " fsa0 " << fsa0 << " fsa1 " << fsa1 << " da0 " << da0 << " ";)


      if (fsa0 + fsa1 == 2)
      {//increase sector
        increase = true;
        if (fsa0 == 0 || (fsa0 == 1 && da0 > M_PI))
        {
          hf_s = ft_hfhs[0];
          hf_e = ft_hfhs[1];
        }
        else
        {
          hf_s = ft_hfhs[1];
          hf_e = ft_hfhs[0];
        }
      }
      else if (fsa0 + fsa1 == 6)
      {//shrink
        if (fsa0 == 4 || (fsa0 == 3 && da0 < M_PI))
        {
          hf_s = ft_hfhs[0];
          hf_e = ft_hfhs[1];
        }
        else
        {
          hf_s = ft_hfhs[1];
          hf_e = ft_hfhs[0];
        }
      }
      else
      {
        ALGOHEX_DEBUG_ONLY(
                std::cerr << "Feature edge: " << _eh << " agl sec0 " << fsa0 << " sec1 " << fsa1 << std::endl;)

        return false;
      }
      ALGOHEX_DEBUG_ONLY(std::cerr << " hfs " << hf_s << " hfe " << hf_e << " ";)

      //update reference axis
      ch_s = mesh_.incident_cell(hf_s);
      hf_start = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_s, heh0));
      rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   heh0,
                                                                   hf_start);
      e_axis = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, heh0, ch_s, rt_axis);
    }
    else
    {
      if (valence_[_eh] == -2)
        increase = true;
    }

    //change the matching
    int mult_trans = 0;
    if (!increase)
      rt_axis = negate(e_axis);
    else
      rt_axis = e_axis;

    mult_trans = rt_axis + 4;

    int newtrans0 = tq_.mult_transitions_idx(trans_prop_[hf_start], tq_.inverse_transition_idx(mult_trans));
    fixed = attempt_matching_adjustment_of_feature_face_sector(heh0, hf_s, hf_e, newtrans0);
  }

  if (fixed)
    return true;

  return false;
}


template<class MeshT>
bool
EnsureEdgeMeshabilityT<MeshT>::
attempt_matching_adjustment_of_feature_face_sector(const HEH _heh, const HFH _hfh_s, const HFH _hfh_e, const int _trans,
                                                   const bool _allow_high_val_sge)
{
//        if(trans_prop_[_hfh_s] > 0) {
//            //coordinate transformation, making matchings on feature face identity
//            coordinate_transformation_in_cell(mesh_, _hfh_s, tq_, cell_quaternions_, trans_prop_);
//        }
//        if(trans_prop_[_hfh_e] > 0) {
//            coordinate_transformation_in_cell(mesh_, _hfh_e, tq_, cell_quaternions_, trans_prop_);
//        }

  int fsa_before = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                       trans_prop_, _heh, _hfh_s,
                                                                                       _hfh_e);

  HFH hf_change = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(_hfh_s, _heh));

  VH vh_to = mesh_.halfedge(mesh_.next_halfedge_in_halfface(_heh, hf_change)).to_vertex();
  if (sgl_vt_[vh_to] != 0)
  {
    VH vhf = mesh_.halfedge(_heh).from_vertex();
    VH vht = mesh_.halfedge(_heh).to_vertex();
    auto hfvhs = mesh_.get_halfface_vertices(hf_change);
    VH vh_new = fs_.face_split(mesh_.face_handle(hf_change));

    //get new face
    if (vh_new.is_valid())
    {
      for (auto i = 0u; i < hfvhs.size(); ++i)
      {
        if (hfvhs[i] == vh_to)
        {
          hfvhs[i] = vh_new;
          break;
        }
      }

      hf_change = mesh_.find_halfface(hfvhs);
    }
  }


  int old_trans = trans_prop_[hf_change];

  trans_prop_[hf_change] = _trans;
  trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(_trans);


  //check
  {
    auto ehi = mesh_.edge_handle(_heh);
    int val = sge_.calculate_edge_valence(ehi);

    if (!_allow_high_val_sge && abs(val) > 1)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "Complex singular edge appears, eh: " << ehi << " valence: " << val
                                   << " restore..." << std::endl;)

      //restore matching
      trans_prop_[hf_change] = old_trans;
      trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(old_trans);

      return false;
    }

    int fsa_after = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                        trans_prop_, _heh, _hfh_s,
                                                                                        _hfh_e);

    if (feature_edge_[mesh_.edge_handle(_heh)] == 0)
    {
      if (val == 0 && fsa_after != 2)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "Feature face sector fix failed! Sector angle before " << fsa_before
                                     << " ff sghe " << mesh_.halfedge(_heh).from_vertex() << " "
                                     << mesh_.halfedge(_heh).to_vertex()
                                     << " fs " << mesh_.face_handle(_hfh_s) << " fe " << mesh_.face_handle(_hfh_e)
                                     << " val " << valence_[mesh_.edge_handle(_heh)]
                                     << " sec agl after " << fsa_after << " dihedral angle "
                                     << dihedral_angle_from_halfface_to_halfface(mesh_, _heh, _hfh_s, _hfh_e)
                                     << std::endl;)

        //restore matching
        trans_prop_[hf_change] = old_trans;
        trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(old_trans);

        return false;
      }
    }
    else
    {
      auto da = dihedral_angle_from_halfface_to_halfface(mesh_, _heh, _hfh_s, _hfh_e);
      if ((da >= M_PI && fsa_after < 2) || fsa_after < 1)
      {
        if (fsa_after < 2)
        {
          ALGOHEX_DEBUG_ONLY(
                  std::cerr << "Feature face sector fix failed at feature edge! Sector angle before " << fsa_before
                            << "ff sghe " << mesh_.halfedge(_heh).from_vertex() << " "
                            << mesh_.halfedge(_heh).to_vertex()
                            << " fs " << mesh_.face_handle(_hfh_s) << " fe " << mesh_.face_handle(_hfh_e) << " val "
                            << valence_[mesh_.edge_handle(_heh)]
                            << " sec agl after " << fsa_after << " dihedral angle "
                            << dihedral_angle_from_halfface_to_halfface(mesh_, _heh, _hfh_s, _hfh_e) << std::endl;)

          //restore matching
          trans_prop_[hf_change] = old_trans;
          trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(old_trans);

          return false;
        }
      }
    }
  }

  std::vector<HFH> ft_hfhs{_hfh_s, _hfh_e};
  for (auto hfhi: ft_hfhs)
  {
    if (!mesh_.is_boundary(mesh_.face_handle(hfhi)))
    {
      CH ch0 = mesh_.incident_cell(hfhi);
      int axis0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch0], mesh_.normal(hfhi)).second;
      CH ch1 = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfhi));
      int axis1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch1], mesh_.normal(hfhi)).second;
      int ax1in0 = tq_.axis_after_transition(axis1, trans_prop_[hfhi]);
      if (axis0 != ax1in0)
      {
        std::cerr << "Error: feature face alignment wrong after fixing sector! ch0 " << ch0 << " ch1 " << ch1 << " fh "
                  << mesh_.face_handle(hfhi) << " ax0 " << axis0 << " ax1 " << axis1 << " ax1in0 " << ax1in0
                  << std::endl;
        //restore matching
        trans_prop_[hf_change] = old_trans;
        trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(old_trans);

        return false;
      }
    }
  }


  //update valence
  for (auto hfe_it = mesh_.hfe_iter(hf_change); hfe_it.valid(); ++hfe_it)
    sge_.compute_edge_valence(*hfe_it);

  for (auto hfv_it = mesh_.hfv_iter(hf_change); hfv_it.valid(); ++hfv_it)
  {
    MeshPropertiesT<MeshT>::update_singular_vertex_property(*hfv_it);
  }

  return true;
}

template<class MeshT>
std::pair<double, int>
EnsureEdgeMeshabilityT<MeshT>::best_transition_idx(const HEH _heh, const HFH _hfh, const std::vector<int> &_cd_vtrans)
{
  VH vht = mesh_.halfedge(_heh).to_vertex();
  VH vhf = mesh_.halfedge(_heh).from_vertex();

  std::set<CH> chs;
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    chs.insert(*vc_it);
  for (auto vc_it = mesh_.vc_iter(vht); vc_it.valid(); ++vc_it)
    chs.insert(*vc_it);

  bool is_ft_node5 = false;
  if (feature_node_[vht])
  {
    if (n_incident_feature_edges(mesh_, feature_edge_, vht) >= 5)
      is_ft_node5 = true;
  }

  HFH hf_opp = mesh_.opposite_halfface_handle(_hfh);

  //store old matching
  std::map<HFH, int> old_hf_trans;
  old_hf_trans[_hfh] = trans_prop_[_hfh];
  old_hf_trans[hf_opp] = trans_prop_[hf_opp];

  //store old quaternions
  std::map<CH, Quaternion> old_qtns;
  for (const auto &chi: chs)
    old_qtns[chi] = cell_quaternions_[chi];

  double best_field_smoothness = std::numeric_limits<double>::max();
  int best_trans = 0;

  EH eh = mesh_.edge_handle(_heh);
  for (auto idx: _cd_vtrans)
  {
//                std::cerr << "candidate trans: " << idx << std::endl;

    trans_prop_[_hfh] = idx;
    trans_prop_[hf_opp] = tq_.inverse_transition_idx(idx);

    QuaternionSmoothing::optimize_quaternions_in_cells(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                       feature_fprop_,
                                                       feature_edge_, chs);


    double fs = frame_smoothness_energy_at_edge(_heh);

    //restore matching and quaternions
    trans_prop_[_hfh] = old_hf_trans[_hfh];
    trans_prop_[hf_opp] = old_hf_trans[hf_opp];

    for (const auto &chi: chs)
      cell_quaternions_[chi] = old_qtns[chi];

    if (fs < best_field_smoothness)
    {
      best_field_smoothness = fs;
      best_trans = idx;
    }
  }
  return {best_field_smoothness, best_trans};
}

template<class MeshT>
double EnsureEdgeMeshabilityT<MeshT>::
frame_smoothness_energy_at_vertex(const VH _vh) const
{
  double energy = 0.;
  for (auto vf_it = mesh_.vf_iter(_vh); vf_it.valid(); ++vf_it)
  {
    energy += frame_smoothness_energy_at_face(*vf_it);
  }

  return energy;
}

template<class MeshT>
double EnsureEdgeMeshabilityT<MeshT>::
frame_smoothness_energy_at_edge(const HEH _heh) const
{
  double energy = 0.;
  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    auto fhi = mesh_.face_handle(*hehf_it);
    energy += frame_smoothness_energy_at_face(fhi);
  }

  return energy;
}

template<class MeshT>
double EnsureEdgeMeshabilityT<MeshT>::
frame_smoothness_energy_at_face(const FH _fh) const
{
  double energy = 0.;
  if (mesh_.is_boundary(_fh))
    return energy;

  auto hfh0 = mesh_.halfface_handle(_fh, 0);
  auto hfh1 = mesh_.opposite_halfface_handle(hfh0);

  auto ch0 = mesh_.incident_cell(hfh0);
  auto ch1 = mesh_.incident_cell(hfh1);

  Quaternion tranq = tq_.transition(trans_prop_[hfh0]);

  Quaternion q0in1 = cell_quaternions_[ch0] * tranq;

  Quaternion rt = cell_quaternions_[ch1] * q0in1.conjugate();

  energy += 2 * std::acos(std::min(std::abs(rt.w()), 1.));

  return energy;
}
}

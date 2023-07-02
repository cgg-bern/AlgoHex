/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define FIELDANGLECALCULATORT_C

#include "FieldAngleCalculatorT.hh"
#include "EdgeMonodromyHelperT.hh"

namespace AlgoHex
{
//template<class MeshT>
//void FieldAngleCalculatorT<MeshT>::
//feature_face_angles_at_feature_edge(const HEH _heh) const
//{
//  std::vector<HFH> hfs;
//  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
//  {
//    if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
//    {
//      hfs.push_back(*hehf_it);
//    }
//  }
//
//  int fsec_agl0 = -1, fsec_agl1 = -1;
//  if (!mesh_.is_boundary(hfs[0]))
//  {
//    fsec_agl0 = feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
//                                                      trans_prop_, _heh, hfs[0], hfs[1]);
//  }
//
//  if (!mesh_.is_boundary(hfs[1]))
//  {
//    fsec_agl1 = feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
//                                                      trans_prop_, _heh, hfs[1], hfs[0]);
//  }
//
//  ALGOHEX_DEBUG_ONLY(std::cerr << "ff sghe " << mesh_.halfedge(_heh).from_vertex() << " "
//                               << mesh_.halfedge(_heh).to_vertex()
//                               << " fs " << mesh_.face_handle(hfs[0]) << " fe " << mesh_.face_handle(hfs[1])
//                               << " " << valence_[mesh_.edge_handle(_heh)] << " sec agl0 " << fsec_agl0
//                               << " sec agl1 " << fsec_agl1 << " dihedral angle0 "
//                               << dihedral_angle_from_halfface_to_halfface(mesh_, _heh, hfs[0], hfs[1]) * 180. / M_PI
//                               << std::endl;)
//}

template<class MeshT>
void FieldAngleCalculatorT<MeshT>::print_feature_sector_angles(const VH _vh) const
{
  //store
  std::vector<std::vector<HEH>> v_sec_hehs;
  std::vector<std::set<FH>> v_sec_fhs;
  get_feature_sectors_at_feature_edge_vertex(mesh_, feature_fprop_, feature_edge_, feature_edge_vertex_, valence_, _vh,
                                             v_sec_hehs, v_sec_fhs);

  for (auto j = 0u; j < v_sec_hehs.size(); ++j)
  {
    if (v_sec_hehs[j].size() < 2 || feature_edge_[mesh_.edge_handle(v_sec_hehs[j].front())] == 0 ||
        feature_edge_[mesh_.edge_handle(v_sec_hehs[j].back())] == 0)
      continue;

    int ea_st = field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                     v_sec_hehs[j], v_sec_fhs[j]);

    ALGOHEX_DEBUG_ONLY(std::cerr << " sector: " << mesh_.edge_handle(v_sec_hehs[j].front()) << " "
                                 << mesh_.edge_handle(v_sec_hehs[j].back())
                                 << " ag " << ea_st << std::endl;)
  }
}

template<class MeshT>
bool FieldAngleCalculatorT<MeshT>::check_field_alignments_at_feature_faces() const
{

  MeshT tmp0;

  for (const auto fhi: mesh_.faces())
  {
    if (feature_fprop_[fhi] > 0)
    {
      if (mesh_.is_boundary(fhi))
        continue;

      HFH hf0 = mesh_.halfface_handle(fhi, 0);
      CH ch0 = mesh_.incident_cell(hf0);
      int axis0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch0], mesh_.normal(hf0)).second;
      CH ch1 = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf0));
      int axis1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch1], mesh_.normal(hf0)).second;
      int ax1in0 = tq_.axis_after_transition(axis1, trans_prop_[hf0]);
      if (axis0 != ax1in0)
      {
        std::cerr << "Error: feature face alignment wrong! ch0 " << ch0 << " ch1 " << ch1 << " fh " << fhi << " ax0 "
                  << axis0 << " ax1 " << axis1 << " ax1in0 " << ax1in0 << std::endl;

        //debug
        {
          std::vector<VH> vhs;
          for (auto fv_it = mesh_.fv_iter(fhi); fv_it.valid(); ++fv_it)
          {
            vhs.push_back(tmp0.add_vertex(mesh_.vertex(*fv_it)));
          }
          tmp0.add_face(vhs);
        }
      }
    }
  }

  for (const auto ehi: mesh_.edges())
  {
    if (feature_face_edge_[ehi] && valence_[ehi] == 0 && feature_edge_[ehi] == 0 && !mesh_.is_boundary(ehi))
    {
      HEH he0 = mesh_.halfedge_handle(ehi, 0);
      HFH hf_s;
      for (auto hehf_it = mesh_.hehf_iter(he0); hehf_it.valid(); ++hehf_it)
      {
        if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
        {
          hf_s = *hehf_it;
          break;
        }
      }

      HFH hf_it = hf_s;
      HFH hf_e(-1);
      CH ch_s = mesh_.incident_cell(hf_it);
      int ax0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], mesh_.normal(hf_s)).second;
      int trans = 0;
      do
      {
        hf_it = mesh_.adjacent_halfface_in_cell(hf_it, he0);
        hf_it = mesh_.opposite_halfface_handle(hf_it);

        if (feature_fprop_[mesh_.face_handle(hf_it)] > 0)
        {
          hf_e = mesh_.opposite_halfface_handle(hf_it);
          break;
        }

        trans = tq_.mult_transitions_idx(trans_prop_[hf_it], trans);
      }
      while (hf_it != hf_s);

      CH ch_e = mesh_.incident_cell(hf_e);
      int ax1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_e], mesh_.normal(hf_e)).second;

      int ax0in1 = tq_.axis_after_transition(ax0, trans);
      if (ax0in1 != ax1)
      {
        std::cerr << "Error: at ffe " << ehi << " axes aligned to ff are not parallel ch0 " << ch_s << " ax0 " << ax0
                  << " ch1 " << ch_e << " ax1 " << ax1 << " ax0in1 " << ax0in1 << std::endl;

        //debug
        {
          auto vh0 = tmp0.add_vertex(mesh_.vertex(mesh_.edge(ehi).from_vertex()));
          auto vh1 = tmp0.add_vertex(mesh_.vertex(mesh_.edge(ehi).to_vertex()));

          tmp0.add_edge(vh0, vh1);
        }
      }

      //other side
      HEH he1 = mesh_.halfedge_handle(ehi, 1);
      hf_s = mesh_.opposite_halfface_handle(hf_s);
      hf_it = hf_s;
      ch_s = mesh_.incident_cell(hf_it);
      ax0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], mesh_.normal(hf_s)).second;
      trans = 0;
      do
      {
        hf_it = mesh_.adjacent_halfface_in_cell(hf_it, he1);
        hf_it = mesh_.opposite_halfface_handle(hf_it);

        if (feature_fprop_[mesh_.face_handle(hf_it)] > 0)
        {
          hf_e = mesh_.opposite_halfface_handle(hf_it);
          break;
        }

        trans = tq_.mult_transitions_idx(trans_prop_[hf_it], trans);
      }
      while (hf_it != hf_s);

      ch_e = mesh_.incident_cell(hf_e);
      ax1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_e], mesh_.normal(hf_e)).second;

      ax0in1 = tq_.axis_after_transition(ax0, trans);
      if (ax0in1 != ax1)
      {
        std::cerr << "Error: at ffe " << ehi << " axes aligned to ff are not parallel ch0 " << ch_s << " ax0 " << ax0
                  << " ch1 " << ch_e << " ax1 " << ax1 << " ax0in1 " << ax0in1 << std::endl;

        //debug
        {
          auto vh0 = tmp0.add_vertex(mesh_.vertex(mesh_.edge(ehi).from_vertex()));
          auto vh1 = tmp0.add_vertex(mesh_.vertex(mesh_.edge(ehi).to_vertex()));

          tmp0.add_edge(vh0, vh1);
        }
      }
    }
  }

  std::time_t result = std::time(nullptr);

  return true;
}


template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::
check_field_alignment_at_feature_edges() const
{
  int n_invalid = 0;
  for (const auto ehi: mesh_.edges())
  {
    if (feature_edge_[ehi] > 0)
    {
      if (!check_field_alignment_at_feature_edge(ehi))
      {
        std::cerr << "Error: feature aligned axis is not correct at edge " << ehi << " valence " << valence_[ehi]
                  << std::endl;
        n_invalid++;
      }
    }
  }

  return n_invalid > 0;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::check_field_alignment_at_feature_edge(const EH _eh) const
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
  if (feature_fprop_[mesh_.face_handle(*hehf_it)] == 0 || !ch_s.is_valid())
  {
    return true;
  }

  fe_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], he_dir).second;

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
      return false;
  }

  return true;
}

template<class MeshT>
bool FieldAngleCalculatorT<MeshT>::check_feature_face_sectors_at_regular_edge(const EH _eh) const
{
  HEH he0 = mesh_.halfedge_handle(_eh, 0);
  HFH hf_s;
  for (auto hehf_it = mesh_.hehf_iter(he0); hehf_it.valid(); ++hehf_it)
  {
    if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
    {
      hf_s = *hehf_it;
      break;
    }
  }

  HFH hf_it = hf_s;
  HFH hf_e(-1);
  CH ch_s = mesh_.incident_cell(hf_it);
  int ax0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], mesh_.normal(hf_s)).second;
  int trans = 0;
  do
  {
    hf_it = mesh_.adjacent_halfface_in_cell(hf_it, he0);
    hf_it = mesh_.opposite_halfface_handle(hf_it);

    if (feature_fprop_[mesh_.face_handle(hf_it)] > 0)
    {
      hf_e = mesh_.opposite_halfface_handle(hf_it);
      break;
    }

    trans = tq_.mult_transitions_idx(trans_prop_[hf_it], trans);
  }
  while (hf_it != hf_s);

  CH ch_e = mesh_.incident_cell(hf_e);
  int ax1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_e], mesh_.normal(hf_e)).second;

  int ax0in1 = tq_.axis_after_transition(ax0, trans);
  if (ax0in1 != ax1)
  {
    return false;
  }

  //other side
  HEH he1 = mesh_.halfedge_handle(_eh, 1);
  hf_s = mesh_.opposite_halfface_handle(hf_s);
  hf_it = hf_s;
  ch_s = mesh_.incident_cell(hf_it);
  ax0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], mesh_.normal(hf_s)).second;
  trans = 0;
  do
  {
    hf_it = mesh_.adjacent_halfface_in_cell(hf_it, he1);
    hf_it = mesh_.opposite_halfface_handle(hf_it);

    if (feature_fprop_[mesh_.face_handle(hf_it)] > 0)
    {
      hf_e = mesh_.opposite_halfface_handle(hf_it);
      break;
    }

    trans = tq_.mult_transitions_idx(trans_prop_[hf_it], trans);
  }
  while (hf_it != hf_s);

  ch_e = mesh_.incident_cell(hf_e);
  ax1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_e], mesh_.normal(hf_e)).second;

  ax0in1 = tq_.axis_after_transition(ax0, trans);
  if (ax0in1 != ax1)
  {
    return false;
  }

  return true;
}

template<class MeshT>
double
FieldAngleCalculatorT<MeshT>::average_field_edge_angle(const MeshT &_mesh, const TransitionQuaternion &_tq,
                                                       const CP <Quaternion> &_cell_quaternions,
                                                       const HFP<int> &_trans_prop, const EP<int> &_valence,
                                                       const EP<int> &_feature_edge, const EH _eh)
{
  auto vhe = _mesh.halfedge_handle(_eh, 0);
  auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells(_mesh, _tq, _cell_quaternions,
                                                                                 _trans_prop, _valence,
                                                                                 _feature_edge, vhe);
  VH vh0 = _mesh.halfedge(vhe).from_vertex();
  VH vh1 = _mesh.halfedge(vhe).to_vertex();
  auto vn = _mesh.vertex(vh1) - _mesh.vertex(vh0);
  vn.normalize();

  double avr_agl = 0.;
  for (auto&[chi, axi]: cell_edge_axis)
  {
    auto ax_dir = AxisAlignmentHelpers::quaternion_vector(_cell_quaternions[chi], AxisAlignment(axi));
//    std::cerr << "  in cell " << chi << " angle: " << std::acos(ax_dir.dot(ovm2eigen(vn))) * 180. / M_PI << " \n";
    avr_agl += std::acos(ax_dir.dot(ovm2eigen(vn))) * 180. / M_PI;
  }
  avr_agl /= (double) cell_edge_axis.size();

  return avr_agl;
}

template<class MeshT>
int
FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(const MeshT &mesh, const TransitionQuaternion &tq,
                                                                    const CP <Quaternion> &cell_quaternions,
                                                                    const HFP<int> &trans_prop, const HEH heh,
                                                                    const HFH hfh_s, const HFH hfh_e)
{
  using Point = typename MeshT::PointT;

  HFH hf_start = mesh.opposite_halfface_handle(mesh.adjacent_halfface_in_cell(hfh_s, heh));

  CH ch0 = mesh.incident_cell(hfh_s);
  Point n0 = mesh.normal(mesh.opposite_halfface_handle(hfh_s));
  int closest_n0 = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch0], n0).second;

  CH ch1 = mesh.incident_cell(mesh.opposite_halfface_handle(hfh_e));
  Point n1 = mesh.normal(hfh_e);
  int closest_n1 = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch1], n1).second;

  int T_inv_id = tq.inverse_transition_idx(
          EdgeMonodromyHelperT::halfedge_transition_index(mesh, tq, trans_prop, heh, hf_start, hfh_e));

  int n1_in_q0 = tq.axis_after_transition(closest_n1, T_inv_id);
  if (closest_n0 != n1_in_q0)
  {
    if (closest_n0 / 2 == n1_in_q0 / 2)
    {
      auto da = dihedral_angle_from_halfface_to_halfface(mesh, heh, hfh_s, hfh_e);
      if (da < 0.8 * M_PI) //TODO: use sector angle
        return 0;
      else if (da >= 0.8 * M_PI)
        return 4;
    }
    else
    {
      int third_axis_id = the_third_axis((AxisAlignment) closest_n0, (AxisAlignment) n1_in_q0);
      auto third_axis = AxisAlignmentHelpers::quaternion_vector(cell_quaternions[ch0],
                                                                (AxisAlignment) third_axis_id);

      VH vh0 = mesh.halfedge(heh).from_vertex();
      VH vh1 = mesh.halfedge(heh).to_vertex();
      Point vn = mesh.vertex(vh0) - mesh.vertex(vh1);

      double sign = Point(third_axis(0), third_axis(1), third_axis(2)) | vn;
      if (sign < -0.)
        return 3;
      else if (sign > 0.)
        return 1;
      else
      {
        std::cout << "Warning: feature face sector is ambiguous! Maybe it's a complex singular edge." << std::endl;
        return -100;
      }
    }
  }
  else
    return 2;


  //never reach here
  return -100;
}

//view from inside, halfedges in sector is in ccw
template<class MeshT>
int
FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(const MeshT &mesh, const TransitionQuaternion &tq,
                                                                   const CP <Quaternion> &cell_quaternions,
                                                                   const HFP<int> &trans_prop, const EP<int> &valence,
                                                                   const std::vector<HEH> &_sector_hehs,
                                                                   const std::set<FH> &_sector_fhs)
{
  using Point = typename MeshT::PointT;

  if (_sector_hehs.empty())
    return -1;

  //only one feature edge
  if (_sector_hehs[0] == _sector_hehs.back())
    return -1;

  //check rotational axes
  int e_dir0 = -1, e_dir1 = -1;

  std::set<FH> visited_fhs;
  //find start and end halfface
  HFH hf_s0(-1), hf_s1(-1);
  for (auto hehf_it = mesh.hehf_iter(_sector_hehs[0]); hehf_it.valid(); ++hehf_it)
  {
    auto fhi = mesh.face_handle(*hehf_it);
    if (_sector_fhs.find(fhi) != _sector_fhs.end())
    {
      hf_s0 = *hehf_it;
      break;
    }
  }

  for (auto hehf_it = mesh.hehf_iter(_sector_hehs.back()); hehf_it.valid(); ++hehf_it)
  {
    auto fhi = mesh.face_handle(*hehf_it);
    if (_sector_fhs.find(fhi) != _sector_fhs.end())
    {
      hf_s1 = *hehf_it;
      visited_fhs.insert(fhi);

      break;
    }
  }
  if (!hf_s0.is_valid() || !hf_s1.is_valid())
    return -1;

  //edge direction of the first halfedge
  int val0 = valence[mesh.edge_handle(_sector_hehs[0])];
  CH ch_s0 = mesh.incident_cell(hf_s0);
  bool is_bdy0 = mesh.is_boundary(_sector_hehs[0]);
  if (!is_bdy0)
    ch_s0 = mesh.incident_cell(mesh.opposite_halfface_handle(hf_s0));;

  if (val0 != 0)
  {
    e_dir0 = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh, tq, cell_quaternions, trans_prop, valence,
                                                                _sector_hehs[0], hf_s0);
    e_dir0 = EdgeMonodromyHelperT::halfedge_axis_idx(mesh, cell_quaternions, _sector_hehs[0], ch_s0, e_dir0);

    //express the dir in ch_s0
    if (!is_bdy0)
    {
      e_dir0 = tq.axis_after_transition(e_dir0, trans_prop[hf_s0]);
    }
  }
  else
  {
    ch_s0 = mesh.incident_cell(hf_s0);
    Point v0 = mesh.vertex(mesh.halfedge(_sector_hehs[0]).to_vertex()) -
               mesh.vertex(mesh.halfedge(_sector_hehs[0]).from_vertex());
    e_dir0 = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch_s0], v0.normalize()).second;
  }

  //rt axis of the second halfedge
  int val1 = valence[mesh.edge_handle(_sector_hehs.back())];
  bool is_bdy1 = mesh.is_boundary(_sector_hehs.back());
  CH ch_s1 = mesh.incident_cell(mesh.opposite_halfface_handle(hf_s1));

  if (valence[mesh.edge_handle(_sector_hehs.back())] != 0)
  {
    //if bdy, ignores the start halfface, but start from bdy halfface
    e_dir1 = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh, tq, cell_quaternions, trans_prop, valence,
                                                                _sector_hehs.back(), hf_s1);;

    if (is_bdy1)
    {
      int trans1 = EdgeMonodromyHelperT::halfedge_transition_index(mesh, tq, trans_prop, _sector_hehs.back(),
                                                                   *mesh.hehf_iter(_sector_hehs.back()));
      e_dir1 = tq.axis_after_transition(e_dir1, trans1);
    }

    e_dir1 = EdgeMonodromyHelperT::halfedge_axis_idx(mesh, cell_quaternions, _sector_hehs.back(), ch_s1, e_dir1);
  }
  else
  {
    Point v1 = mesh.vertex(mesh.halfedge(_sector_hehs.back()).to_vertex()) -
               mesh.vertex(mesh.halfedge(_sector_hehs.back()).from_vertex());
    e_dir1 = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch_s1], v1.normalize()).second;
  }

  //express e_dir0 in ch1
  int trans = 0;
  for (int i = (int) _sector_hehs.size() - 2; i > 0; --i)
  {
    //halfedge halfface should be oriented
    HFH hf_ss(-1);
    for (auto hehf_it = mesh.hehf_iter(_sector_hehs[i]); hehf_it.valid(); ++hehf_it)
    {
      auto fhi = mesh.face_handle(*hehf_it);
      if (_sector_fhs.find(fhi) != _sector_fhs.end() && visited_fhs.find(fhi) != visited_fhs.end())
      {
        hf_ss = *hehf_it;
        break;
      }
    }

    if (!hf_ss.is_valid())
      return -1;

    int n = 0;
    HFH hf_it = hf_ss;
    do
    {
      hf_it = mesh.adjacent_halfface_in_cell(hf_it, _sector_hehs[i]);
      if (!hf_it.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! he_s: " << _sector_hehs[i] << " hfh_s " << hf_ss
                  << std::endl;
        return -1;
      }

      hf_it = mesh.opposite_halfface_handle(hf_it);

      auto fh_it = mesh.face_handle(hf_it);
      if (_sector_fhs.find(fh_it) != _sector_fhs.end() && visited_fhs.find(fh_it) == visited_fhs.end())
      {
        visited_fhs.insert(fh_it);
        break;
      }

      trans = tq.mult_transitions_idx(trans_prop[hf_it], trans);

      n++;

      if (n > 100)
      {
        std::cerr << "Error: mesh is likely changed in measuring zero sector angle" << std::endl;
        return -1;
      }

    }
    while (hf_it != hf_ss);
  }

  int eax1_in0 = tq.axis_after_transition(e_dir1, trans);

  if (eax1_in0 / 2 == e_dir0 / 2)
  {
    if (eax1_in0 == e_dir0)
    {
      auto sc_agl = sector_angle(mesh, _sector_hehs);
      if (sc_agl < M_PI)
        return 0;//zero sector
      else
        return 4;//2 PI
    }
    else //parallel
      return 2;
  }
  else
  { //corner
    auto nm = mesh.normal(hf_s0);
    int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions[mesh.incident_cell(hf_s0)],
                                                     nm.normalize()).second;

    int ax3 = the_third_axis((AxisAlignment) nm_axis, (AxisAlignment) e_dir0);

    int n_eax1_in0 = eax1_in0 % 2 == 0 ? eax1_in0 + 1 : eax1_in0 - 1;
    if (ax3 == eax1_in0)
      return 1;
    else if (ax3 == n_eax1_in0)
      return 3;
    else
    {//e.g. singular edge is snapped to feature edge
      std::cerr << "Warning: sector classification is wrong. Probably there is a singular edge in the sector. "
                << std::endl;
      return -1;
    }
  }

  return -1;
}

template<class MeshT>
double FieldAngleCalculatorT<MeshT>::sector_angle(const MeshT &mesh, const std::vector<HEH> &hehs)
{
  double angle_sum = 0.;
  for (auto i = 0u; i < hehs.size() - 1; ++i)
  {
    auto dir0 = mesh.vertex(mesh.halfedge(hehs[i]).to_vertex()) -
                mesh.vertex(mesh.halfedge(hehs[i]).from_vertex());
    dir0.normalize();
    auto dir1 = mesh.vertex(mesh.halfedge(hehs[i + 1]).to_vertex()) -
                mesh.vertex(mesh.halfedge(hehs[i + 1]).from_vertex());
    dir1.normalize();

    auto dprd = dir0 | dir1;
    dprd = std::min(1., std::max(-1., dprd));

    angle_sum += std::acos(dprd);
  }

  return angle_sum;
}

template<class MeshT>
int
FieldAngleCalculatorT<MeshT>::field_angle_of_feature_edge_sector(const HEH _heh, const HFH _hfh, const bool _cw) const
{
  using Point = typename MeshT::PointT;

  EH eh = mesh_.edge_handle(_heh);
  if (feature_edge_[eh] == 0)
    return -1;

  //rotational axes
  int e_dir0 = -1, e_dir1 = -1;

  FH fh = mesh_.face_handle(_hfh);
  if (feature_fprop_[fh] == 0)
    return -1;

  //edge direction of the first halfedge
  int val0 = valence_[eh];
  if (abs(val0) == 1)
  {
    e_dir0 = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                _heh, _hfh);
    if (val0 == 1)
      e_dir0 = e_dir0 % 2 == 0 ? e_dir0 + 1 : e_dir0 - 1;
    //express the dir in the cell incident to the feature face
    if (!mesh_.is_boundary(_heh))
    {
      e_dir0 = tq_.axis_after_transition(e_dir0, trans_prop_[_hfh]);
    }
  }
  else
  {
    Point v0 = mesh_.vertex(mesh_.halfedge(_heh).to_vertex()) - mesh_.vertex(mesh_.halfedge(_heh).from_vertex());
    e_dir0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[mesh_.incident_cell(_hfh)],
                                                v0.normalize()).second;
  }
  //outgoing
  if (_cw)
    e_dir0 = negate(e_dir0);

  std::vector<HEH> sec_hehs;

  HEH he_e;
  HFH hf_e;
  int trans = 0;
  if (_cw)
    std::tie(he_e, hf_e, trans) = next_outgoing_special_hehf_with_trans_cw(_heh, _hfh, sec_hehs);
  else
    std::tie(he_e, hf_e, trans) = next_incoming_special_hehf_with_trans_ccw(_heh, _hfh, sec_hehs);

  if (!he_e.is_valid())
    return -1;


  //rt axis of the second halfedge
  int val1 = valence_[mesh_.edge_handle(he_e)];
  if (abs(val1) == 1)
  {
    e_dir1 = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                he_e, hf_e);
    if (val1 == 1)
      e_dir1 = e_dir1 % 2 == 0 ? e_dir1 + 1 : e_dir1 - 1;
    if (!mesh_.is_boundary(he_e))
    {
      e_dir1 = tq_.axis_after_transition(e_dir1, trans_prop_[hf_e]);
    }
  }
  else
  {
    Point v1 = mesh_.vertex(mesh_.halfedge(he_e).to_vertex()) - mesh_.vertex(mesh_.halfedge(he_e).from_vertex());
    e_dir1 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[mesh_.incident_cell(hf_e)], v1.normalize()).second;
  }

  if (!_cw)
    e_dir1 = negate(e_dir1);

  //express e_dir0 in ch1
  int eax0_in1 = tq_.axis_after_transition(e_dir0, trans);

  if (eax0_in1 / 2 == e_dir1 / 2)
  {
    if (eax0_in1 == e_dir1)
    {
      auto sc_agl = sector_angle(mesh_, sec_hehs);
      if (sc_agl < M_PI)
        return 0;//zero sector
      else
        return 4;//2 PI
    }
    else //parallel
      return 2;
  }
  else
  { //corner
    HFH hf_bdy1 = mesh_.opposite_halfface_handle(hf_e);
    auto nm = mesh_.normal(hf_bdy1);
    int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[mesh_.incident_cell(hf_e)],
                                                     nm.normalize()).second;

    int ax3 = the_third_axis((AxisAlignment) e_dir1, (AxisAlignment) nm_axis);
    if (!_cw)
      ax3 = negate(ax3);

    int n_eax0_in1 = eax0_in1 % 2 == 0 ? eax0_in1 + 1 : eax0_in1 - 1;
    if (ax3 == eax0_in1)
      return 1;
    else if (ax3 == n_eax0_in1)
      return 3;
    else
    {//e.g. singular edge is snapped to feature edge
      std::cerr << "Warning: sector classification is wrong. Probably there is a singular edge in the sector. "
                << std::endl;
      return -1;
    }
  }

  return -1;
}

//input is incoming halfedge and incident halfface
template<class MeshT>
std::tuple<HEH, HFH, int>
FieldAngleCalculatorT<MeshT>::next_outgoing_special_hehf_with_trans_cw(const HEH _heh, const HFH _hfh,
                                                                       std::vector<HEH> &_sec_hehs) const
{
  _sec_hehs.push_back(mesh_.opposite_halfedge_handle(_heh));

  HEH he_itt = _heh;
  HFH hf_itt = _hfh;
  int trans_prd = 0;
  do
  {
    he_itt = mesh_.next_halfedge_in_halfface(he_itt, hf_itt);

    //store sector halfedges
    _sec_hehs.push_back(he_itt);

    //if special edge, return
    EH eh_itt = mesh_.edge_handle(he_itt);
    if (feature_edge_[eh_itt] > 0 || valence_[eh_itt] != 0)
    {
      return std::make_tuple(he_itt, hf_itt, trans_prd);
    }

    //move to next halfface on feature surface
    HFH hf_it_in = hf_itt;
    do
    {
      //next incident halfface
      hf_it_in = mesh_.adjacent_halfface_in_cell(hf_it_in, he_itt);
      hf_it_in = mesh_.opposite_halfface_handle(hf_it_in);

      if (!hf_it_in.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! he_s: " << he_itt << " hfh_s " << hf_itt << std::endl;
        return std::make_tuple(HEH(-1), HFH(-1), -1);
      }
      if (feature_fprop_[mesh_.face_handle(hf_it_in)] > 0 || hf_it_in == hf_itt)
        break;

      trans_prd = tq_.mult_transitions_idx(trans_prop_[hf_it_in], trans_prd);
    }
    while (hf_it_in != hf_itt);

    he_itt = mesh_.opposite_halfedge_handle(he_itt);
    hf_itt = mesh_.opposite_halfface_handle(hf_it_in);
  }
  while (he_itt != _heh);

  return std::make_tuple(HEH(-1), HFH(-1), -1);
}

//input is outgoing
template<class MeshT>
std::tuple<HEH, HFH, int>
FieldAngleCalculatorT<MeshT>::next_incoming_special_hehf_with_trans_ccw(const HEH _heh, const HFH _hfh,
                                                                        std::vector<HEH> &_sec_hehs) const
{
  HEH he_itt = _heh;
  HFH hf_itt = _hfh;
  int trans_prd = 0;
  do
  {
    he_itt = mesh_.prev_halfedge_in_halfface(he_itt, hf_itt);

    //store sector halfedges
    _sec_hehs.push_back(he_itt);

    //if special edge, return
    EH eh_itt = mesh_.edge_handle(he_itt);
    if (feature_edge_[eh_itt] > 0 || valence_[eh_itt] != 0)
    {
      return std::make_tuple(he_itt, hf_itt, trans_prd);
    }

    //move to next halfface on feature surface
    HFH hf_it_in = hf_itt;
    do
    {
      //next incident halfface
      hf_it_in = mesh_.adjacent_halfface_in_cell(hf_it_in, he_itt);
      hf_it_in = mesh_.opposite_halfface_handle(hf_it_in);

      if (!hf_it_in.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! he_s: " << he_itt << " hfh_s " << hf_itt << std::endl;
        trans_prd = -1;
        return std::make_tuple(HEH(-1), HFH(-1), -1);
      }
      if (feature_fprop_[mesh_.face_handle(hf_it_in)] > 0 || hf_it_in == hf_itt)
        break;

      trans_prd = tq_.mult_transitions_idx(trans_prop_[hf_it_in], trans_prd);
    }
    while (hf_it_in != hf_itt);

    he_itt = mesh_.opposite_halfedge_handle(he_itt);
    hf_itt = mesh_.opposite_halfface_handle(hf_it_in);
  }
  while (he_itt != _heh);

  return std::make_tuple(HEH(-1), HFH(-1), -1);
}


template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::is_axis_in_cell_parallel_to_any_surface(const VH _vh, const int _axis, const CH _ch)
{
  auto or_chs = get_onering_cells(mesh_, _vh);
  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, _ch, ft_hfhs);

  for (const auto hfhi: ft_hfhs)
  {
    int is_prl = angle_of_axis_in_cell_and_normal_direction(_axis, _ch, hfhi, por_chs);
    if (is_prl == 0)
    {
      return true;
    }
  }

  return false;
}

//search till an orthogonal sector is hit
template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::has_any_convexly_orthogonal_face(const VH _vh, const FH _fh) const
{
  std::vector<HFH> hfhs{mesh_.halfface_handle(_fh, 0), mesh_.halfface_handle(_fh, 1)};

  for (const auto hfhi: hfhs)
  {
    if (!mesh_.is_boundary(hfhi))
    {
      if (has_any_convexly_orthogonal_face_cw(_vh, hfhi))
        return true;

      if (has_any_convexly_orthogonal_face_ccw(_vh, hfhi))
        return true;
    }
  }


  return false;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::has_any_convexly_orthogonal_face_cw(const VH _vh, const HFH _hfh) const
{
  CH ch_ft = mesh_.incident_cell(_hfh);
  auto nm = mesh_.normal(_hfh);
  int nm_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_ft], nm).second;


  std::set<FH> visited_fhs;
  visited_fhs.insert(mesh_.face_handle(_hfh));

  std::queue<HALFFACEINFO2> que;
  que.push(HALFFACEINFO2(_hfh, 0., 0));

  while (!que.empty())
  {
    auto hfi_cur = que.front();
    que.pop();

    //find incident halfedge
    HEH he_s;
    for (auto hfhe_it = mesh_.hfhe_iter(hfi_cur.hfh); hfhe_it.valid(); ++hfhe_it)
      if (mesh_.halfedge(*hfhe_it).from_vertex() == _vh)
      {
        he_s = *hfhe_it;
        break;
      }


    if (he_s.is_valid())
    {
      auto next_hfi = next_halfface_info(hfi_cur, he_s, nm_ax);
      auto next_fh = mesh_.face_handle(next_hfi.hfh);
      if (visited_fhs.find(next_fh) == visited_fhs.end())
      {
        visited_fhs.insert(next_fh);
        int nm_in_chi = tq_.axis_after_transition(nm_ax, next_hfi.trans);

        CH chi = mesh_.incident_cell(next_hfi.hfh);
        auto nm_cur = mesh_.normal(next_hfi.hfh);
        int nm_ax_cur = AxisAlignmentHelpers::closest_axis(cell_quaternions_[chi], nm_cur).second;

        if (nm_in_chi / 2 != nm_ax_cur / 2)
        {
          if (next_hfi.u_dist > 0.)
            return true;
          else if (next_hfi.u_dist < 0.)
            return false;
        }

        que.push(next_hfi);
      }
    }
  }

  return false;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::has_any_convexly_orthogonal_face_ccw(const VH _vh, const HFH _hfh) const
{
  CH ch_ft = mesh_.incident_cell(_hfh);
  auto nm = mesh_.normal(_hfh);
  int nm_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_ft], nm).second;


  std::set<FH> visited_fhs;
  visited_fhs.insert(mesh_.face_handle(_hfh));

  std::queue<HALFFACEINFO2> que;
  que.push(HALFFACEINFO2(_hfh, 0., 0));

  while (!que.empty())
  {
    auto hfi_cur = que.front();
    que.pop();

    //find incident halfedge
    HEH he_s;
    for (auto hfhe_it = mesh_.hfhe_iter(hfi_cur.hfh); hfhe_it.valid(); ++hfhe_it)
      if (mesh_.halfedge(*hfhe_it).to_vertex() == _vh)
      {
        he_s = *hfhe_it;
        break;
      }


    if (he_s.is_valid())
    {
      auto next_hfi = next_halfface_info(hfi_cur, he_s, nm_ax);
      auto next_fh = mesh_.face_handle(next_hfi.hfh);
      if (visited_fhs.find(next_fh) == visited_fhs.end())
      {
        visited_fhs.insert(next_fh);
        int nm_in_chi = tq_.axis_after_transition(nm_ax, next_hfi.trans);

        CH chi = mesh_.incident_cell(next_hfi.hfh);
        auto nm_cur = mesh_.normal(next_hfi.hfh);
        int nm_ax_cur = AxisAlignmentHelpers::closest_axis(cell_quaternions_[chi], nm_cur).second;

        if (nm_in_chi / 2 != nm_ax_cur / 2)
        {
          if (next_hfi.u_dist > 0.)
            return true;
          else if (next_hfi.u_dist < 0.)
            return false;
        }

        que.push(next_hfi);
      }
    }
  }

  return false;
}

template<class MeshT>
typename FieldAngleCalculatorT<MeshT>::HALFFACEINFO2
FieldAngleCalculatorT<MeshT>::next_halfface_info(const HALFFACEINFO2 &hi_cur, const HEH _he_s, const int _u) const
{
  //the first segment
  auto f_bct = mesh_.barycenter(mesh_.face_handle(hi_cur.hfh));
  auto e_bct = mesh_.barycenter(mesh_.edge_handle(_he_s));
  Vec3d dir = ovm2eigen(e_bct - f_bct);

  auto ch_cur = mesh_.incident_cell(hi_cur.hfh);
  //u expressed in the current cell
  int u_cur = tq_.axis_after_transition(_u, hi_cur.trans);
  Vec3d u_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_cur], (AxisAlignment) u_cur);

  //first segment
  double u_dist = hi_cur.u_dist;

  u_dist += u_dir.dot(dir);

  //find next feature face
  FH fh_cur = mesh_.face_handle(hi_cur.hfh);
  HFH hfh_next;
  for (auto hehf_it = mesh_.hehf_iter(_he_s); hehf_it.valid(); ++hehf_it)
  {
    auto fh_next = mesh_.face_handle(*hehf_it);
    if (feature_fprop_[fh_next] > 0 && fh_next != fh_cur)
    {
      hfh_next = *hehf_it;
      break;
    }
  }

  HFH hf_n_opp = mesh_.opposite_halfface_handle(hfh_next);
  //the second segment
  auto ch_next = mesh_.incident_cell(hf_n_opp);

  auto f_next_bct = mesh_.barycenter(mesh_.face_handle(hf_n_opp));
  dir = ovm2eigen(f_next_bct - e_bct);

  HFH hfh_ns = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hi_cur.hfh, _he_s));
  int trans_next = 0;
  if (feature_fprop_[mesh_.face_handle(hfh_ns)] == 0)
    trans_next = tq_.mult_transitions_idx(
            EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, _he_s, hfh_ns, hfh_next),
            hi_cur.trans);
  //u expressed in next cell
  int u_n = tq_.axis_after_transition(_u, trans_next);
  Vec3d un_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_next], (AxisAlignment) u_n);

  u_dist += un_dir.dot(dir);

  HALFFACEINFO2 next_hi(hf_n_opp, u_dist, trans_next);

  return next_hi;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::has_special_edge_in_direction(const VH _vh, const CH _ch_s, const int _ax_s)
{
  //one ring cells
  std::set<CH> or_chs = get_onering_cells(mesh_, _vh);

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, _ch_s, ft_hfhs);

  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    EH ve = mesh_.edge_handle(*voh_it);
    if (feature_edge_[ve] > 0 || (feature_face_edge_[ve] && valence_[ve] != 0))
    {
      CH chi;
      for (auto hec_it = mesh_.hec_iter(*voh_it); hec_it.valid(); ++hec_it)
      {
        if (por_chs.find(*hec_it) != por_chs.end())
        {
          chi = *hec_it;
          break;
        }
      }

      if (chi.is_valid())
      {
        int ax_in_i = axis_in_chs_expressed_in_cht(_ch_s, chi, _ax_s, por_chs);

        auto dir0 = mesh_.vertex(mesh_.halfedge(*voh_it).to_vertex()) -
                    mesh_.vertex(mesh_.halfedge(*voh_it).from_vertex());
        int ax_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[chi],
                                                      dir0.normalize()).second;

        if (ax_i == ax_in_i)
        {
          return true;
        }
      }
    }
  }

  return false;
}


template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::has_singular_edge_in_direction(const HEH _heh_query, const CH _ch_query, int _axis_query,
                                                             const VH _vh_query,
                                                             CH &_ch_found, HFH &_hfh_found, HEH &_heh_found,
                                                             int &_axis_found) const
{
  //one ring cells
  std::set<CH> or_chs = get_onering_cells(mesh_, _vh_query);

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh_query, or_chs, _ch_query, ft_hfhs);

  //find valence > 0
  for (auto voh_it = mesh_.voh_iter(_vh_query); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh_.edge_handle(*voh_it);
    if (valence_[ve] > 0 && *voh_it != _heh_query)
    {
      CH ch_i(-1);
      HFH hf_i(-1);
      int ax_i = -1;
      if (feature_face_edge_[ve])
      {
        //same sector
        HFH hf_start, hf_end;
        for (auto hehf_it = mesh_.hehf_iter(*voh_it); hehf_it.valid(); ++hehf_it)
        {
          if (ft_hfhs.find(*hehf_it) != ft_hfhs.end())
          {
            hf_start = *hehf_it;
          }
          else if (ft_hfhs.find(mesh_.opposite_halfface_handle(*hehf_it)) != ft_hfhs.end())
          {
            hf_end = *hehf_it;
          }
        }

        ch_i = mesh_.incident_cell(hf_start);
        hf_i = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_start, *voh_it));
        ax_i = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  *voh_it,
                                                                  hf_i);
      }
      else
      {
        hf_i = *mesh_.hehf_iter(*voh_it);
        ch_i = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_i));
        ax_i = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  *voh_it,
                                                                  hf_i);
      }
      //e dir
      int e_ax_i = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, *voh_it, ch_i, ax_i);

      int e_ax_i_in_ss = axis_in_chs_expressed_in_cht(ch_i, _ch_query, e_ax_i, por_chs);

      if (e_ax_i_in_ss == _axis_query)
      {
        _ch_found = ch_i;
        _hfh_found = hf_i;
        _heh_found = *voh_it;
        _axis_found = e_ax_i;

        ALGOHEX_DEBUG_ONLY(std::cerr << "found sg eh: " << mesh_.edge_handle(_heh_found) << " val "
                                     << valence_[mesh_.edge_handle(_heh_found)] << " ch_found" << _ch_found
                                     << " e_ax_i: "
                                     << _axis_found << " rt axi " << ax_i << std::endl;)

        return true;
      }
    }
  }

  return false;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::has_singular_edge_in_direction_longest(const HEH _heh_query, const CH _ch_query,
                                                                     int _axis_query, const VH _vh_query,
                                                                     CH &_ch_found, HFH &_hfh_found, HEH &_heh_found,
                                                                     int &_axis_found) const
{
  //one ring cells
  std::set<CH> or_chs = get_onering_cells(mesh_, _vh_query);

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh_query, or_chs, _ch_query, ft_hfhs);

  double max_len = -1.;
  //find valence > 0
  for (auto voh_it = mesh_.voh_iter(_vh_query); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh_.edge_handle(*voh_it);
    if (valence_[ve] > 0 && *voh_it != _heh_query)
    {
      CH ch_i(-1);
      HFH hf_i(-1);
      int ax_i = -1;
      if (feature_face_edge_[ve])
      {
        //same sector
        HFH hf_start, hf_end;
        for (auto hehf_it = mesh_.hehf_iter(*voh_it); hehf_it.valid(); ++hehf_it)
        {
          if (ft_hfhs.find(*hehf_it) != ft_hfhs.end())
          {
            hf_start = *hehf_it;
          }
          else if (ft_hfhs.find(mesh_.opposite_halfface_handle(*hehf_it)) != ft_hfhs.end())
          {
            hf_end = *hehf_it;
          }
        }

        ch_i = mesh_.incident_cell(hf_start);
        hf_i = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_start, *voh_it));
        ax_i = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  *voh_it,
                                                                  hf_i);
      }
      else
      {
        hf_i = *mesh_.hehf_iter(*voh_it);
        ch_i = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_i));
        ax_i = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  *voh_it,
                                                                  hf_i);
      }
      //e dir
      int e_ax_i = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, *voh_it, ch_i, ax_i);

      int e_ax_i_in_ss = axis_in_chs_expressed_in_cht(ch_i, _ch_query, e_ax_i, por_chs);

      if (e_ax_i_in_ss == _axis_query)
      {

        auto sg_hehs_i = get_halfedges_on_singular_arc_directional(mesh_, valence_, *voh_it);
        bool has_diff_val = false;
        for (auto hehi: sg_hehs_i)
          if (valence_[mesh_.edge_handle(hehi)] != valence_[mesh_.edge_handle(sg_hehs_i[0])])
          {
            has_diff_val = true;
            break;
          }
        if (has_diff_val)
        {
          std::cerr << "Different valence on sg arc, continue..." << std::endl;
          continue;
        }


        double len = 0.;
        for (auto heh_j: sg_hehs_i)
          len += mesh_.length(heh_j);

        std::cerr << "\nsg_hehs size: " << sg_hehs_i.size() << " maxlen " << max_len << std::endl;
        if (len > max_len)
        {
          max_len = len;
          _ch_found = ch_i;
          _hfh_found = hf_i;
          _heh_found = *voh_it;
          _axis_found = e_ax_i;

          ALGOHEX_DEBUG_ONLY(std::cerr << "found sg eh: " << mesh_.edge_handle(_heh_found) << " val "
                                       << valence_[mesh_.edge_handle(_heh_found)] << " ch_found"
                                       << _ch_found << " e_ax_i: "
                                       << _axis_found << " rt axi " << ax_i << std::endl;)
        }
      }
    }
  }

  if (_heh_found.is_valid())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "longest found sg eh: " << mesh_.edge_handle(_heh_found) << " val "
                                 << valence_[mesh_.edge_handle(_heh_found)] << " ch_found"
                                 << _ch_found << " e_ax_i: "
                                 << _axis_found << " rt axi " << std::endl;)
    return true;
  }

  return false;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::is_feature_edge_orthogonal_to_any_sector(const HEH _ft_heh)
{
  VH vhf = mesh_.halfedge(_ft_heh).from_vertex();

  auto or_chs = get_onering_cells(mesh_, vhf);
  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vhf, or_chs, *mesh_.hec_iter(_ft_heh), ft_hfhs);
  CH ch_s(-1);
  for (auto hec_it = mesh_.hec_iter(_ft_heh); hec_it.valid(); ++hec_it)
  {
    if (por_chs.find(*hec_it) != por_chs.end())
    {
      ch_s = *hec_it;
      break;
    }
  }

  auto e_dir = mesh_.vertex(mesh_.halfedge(_ft_heh).to_vertex()) - mesh_.vertex(mesh_.halfedge(_ft_heh).from_vertex());
  int fe_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], e_dir.normalize()).second;

  for (const auto hfhi: ft_hfhs)
  {
    int is_prl = angle_of_axis_in_cell_and_normal_direction(fe_axis, ch_s, hfhi, por_chs);
    if (is_prl == 1 || is_prl == 2)
    {
      return true;
    }
  }

  return false;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::is_singular_edge_tangent_to_any_surface(const std::set<CH> &_or_chs,
                                                                      const HEH _sg_heh) const
{
  VH vh_f = mesh_.halfedge(_sg_heh).from_vertex();
  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vh_f, _or_chs, *mesh_.hec_iter(_sg_heh), ft_hfhs);
  CH ch_s(-1);
  HFH hf_s(-1);
  for (auto hehf_it = mesh_.hehf_iter(_sg_heh); hehf_it.valid(); ++hehf_it)
  {
    auto chi = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
    if (por_chs.find(chi) != por_chs.end())
    {
      ch_s = chi;
      hf_s = *hehf_it;
      break;
    }
  }

  int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   _sg_heh, hf_s);
  for (const auto hfhi: ft_hfhs)
  {
    int is_prl = angle_of_axis_in_cell_and_normal_direction(rt_axis, ch_s, hfhi, por_chs);
    if (is_prl == 0)
    {
      return true;
    }
  }

  return false;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::is_singular_edge_orthogonal_to_any_surface(const std::set<CH> &_or_chs, const HEH _sg_heh)
{
  VH vh_f = mesh_.halfedge(_sg_heh).from_vertex();

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vh_f, _or_chs, *mesh_.hec_iter(_sg_heh), ft_hfhs);

  HFH hf_s = *mesh_.hehf_iter(_sg_heh);
  CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_s));

  int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   _sg_heh, hf_s);

  for (const auto hfhi: ft_hfhs)
  {
    int is_prl = angle_of_axis_in_cell_and_normal_direction(rt_axis, ch_s, hfhi, por_chs);
    if (is_prl == 1 || is_prl == 2)
    {
      return true;
    }
  }

  return false;
}

template<class MeshT>
bool
FieldAngleCalculatorT<MeshT>::is_singular_edge_tangent_to_all_surface(const HEH _sg_heh) const
{
  VH vh_f = mesh_.halfedge(_sg_heh).from_vertex();

  auto or_chs = get_onering_cells(mesh_, vh_f);

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vh_f, or_chs, *mesh_.hec_iter(_sg_heh), ft_hfhs);
  CH ch_s(-1);
  HFH hf_s(-1);
  HEH he_s = mesh_.opposite_halfedge_handle(_sg_heh);
  for (auto hehf_it = mesh_.hehf_iter(he_s); hehf_it.valid(); ++hehf_it)
  {
    auto chi = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
    if (por_chs.find(chi) != por_chs.end())
    {
      ch_s = chi;
      hf_s = *hehf_it;
      break;
    }
  }

  int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   he_s, hf_s);
  for (const auto hfhi: ft_hfhs)
  {
    int is_prl = angle_of_axis_in_cell_and_normal_direction(rt_axis, ch_s, hfhi, por_chs);
    if (is_prl == 1 || is_prl == 2)
    {
      return false;
    }
  }

  return true;
}

template<class MeshT>
std::vector<int>
FieldAngleCalculatorT<MeshT>::rotation_axes_expressed_in_common_tet(const VH _vh, const std::vector<HEH> &_sg_hehs)
{
  std::vector<int> axes;

  auto v_size = _sg_hehs.size();

  if (v_size < 2)
    return axes;

  //one ring cells
  std::set<CH> or_chs;
  for (auto cv_it = mesh_.vc_iter(_vh); cv_it.valid(); ++cv_it)
    or_chs.insert(*cv_it);

  //rotational axis of the first halfedge
  HFH hfh_t = *mesh_.hehf_iter(_sg_hehs[0]);
  int axis_s0 = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   _sg_hehs[0], hfh_t);

  axes.push_back(axis_s0);

  CH ch_t(-1);
  if (mesh_.is_boundary(_sg_hehs[0]))
    ch_t = mesh_.incident_cell(hfh_t);
  else
    ch_t = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh_t));

  auto pre_hf = spanning_tree_in_onering(ch_t, or_chs);

  for (auto i = 1u; i < v_size; ++i)
  {
    HFH hfh_s = *mesh_.hehf_iter(_sg_hehs[i]);
    int axis_si = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                     valence_, _sg_hehs[i], hfh_s);

    CH ch_s(-1);
    if (mesh_.is_boundary(_sg_hehs[i]))
      ch_s = mesh_.incident_cell(hfh_s);
    else
      ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh_s));

    auto dpath = get_dual_path_to_cell_in_onering(ch_s, pre_hf);

    int axis_i = axis_after_transition_along_dual_path(axis_si, dpath, true);

    axes.push_back(axis_i);
  }

  return axes;
}

template<class MeshT>
int
FieldAngleCalculatorT<MeshT>::angle_of_axis_in_cell_and_normal_direction(const int _axis, const CH _ch,
                                                                         const HFH _ft_hfh,
                                                                         const std::set<CH> _onering_chs) const
{
  CH ch_ft = mesh_.incident_cell(_ft_hfh);
  auto nm = mesh_.normal(_ft_hfh);
  int nm_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_ft], nm).second;

  int nm_in_sgch = axis_in_chs_expressed_in_cht(ch_ft, _ch, nm_ax, _onering_chs);

  if (nm_in_sgch / 2 == _axis / 2)
  {
    if (nm_in_sgch != _axis) //opposite direction
      return 1;
    else //same direction
      return 2;
  }

  return 0;
}

//one ring
template<class MeshT>
int
FieldAngleCalculatorT<MeshT>::axis_in_chs_expressed_in_cht(const CH _ch_s, const CH _ch_t, int _axis,
                                                           const std::set<CH> _onering_chs) const
{
  auto dpath = get_dual_path_between_cells_in_onering(_ch_s, _ch_t, _onering_chs);

  return axis_after_transition_along_dual_path(_axis, dpath);
}


//one ring
template<class MeshT>
int
FieldAngleCalculatorT<MeshT>::axis_after_transition_along_dual_path(const int _axis, const std::vector<HFH> _dpath,
                                                                    const bool _reversed) const
{
  //transition of dual path
  int dp_trans_idx = 0;
  for (const auto hfh: _dpath)
    dp_trans_idx = tq_.mult_transitions_idx(trans_prop_[hfh], dp_trans_idx);

  if (_reversed)
    dp_trans_idx = tq_.inverse_transition_idx(dp_trans_idx);
//std::cerr<<" trans: "<<dp_trans_idx<<" axis: "<<_axis;
  return tq_.axis_after_transition(_axis, dp_trans_idx);
}

template<class MeshT>
std::vector<HFH>
FieldAngleCalculatorT<MeshT>::
get_dual_path_between_cells_in_onering(const CH _ch_s, const CH _ch_t, const std::set<CH> &_onering_chs) const
{
  auto pre_hf = spanning_tree_in_onering(_ch_s, _onering_chs);

  std::vector<HFH> dpath;

  //get the path
  std::queue<HFH> que_dpath;
  que_dpath.push(pre_hf[_ch_t]);
  while (!que_dpath.empty())
  {
    auto hf = que_dpath.front();
    que_dpath.pop();
    if (hf == HFH(-1))
      break;

    dpath.push_back(mesh_.opposite_halfface_handle(hf));
    que_dpath.push(pre_hf[mesh_.incident_cell(hf)]);
  }

  std::reverse(dpath.begin(), dpath.end());


  return dpath;
}


template<class MeshT>
std::map<CH, HFH>
FieldAngleCalculatorT<MeshT>::
spanning_tree_in_onering(const CH _ch_seed, const std::set<CH> &_onering_chs) const
{
  std::map<CH, HFH> pre_hf;
  pre_hf[_ch_seed] = HFH(-1);

  std::queue<CH> que;
  que.push(_ch_seed);

  while (!que.empty())
  {
    auto ch_cur = que.front();
    que.pop();

    auto hfs = mesh_.cell(ch_cur).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(hfs[i]);
      auto ch_opp = mesh_.incident_cell(hf_opp);

      if (!ch_opp.is_valid())
        continue;

      if (_onering_chs.find(ch_opp) != _onering_chs.end() && pre_hf.find(ch_opp) == pre_hf.end())
      {
        pre_hf[ch_opp] = hfs[i];
        que.push(ch_opp);
      }
    }
  }

  return pre_hf;
}

template<class MeshT>
std::vector<HFH>
FieldAngleCalculatorT<MeshT>::
get_dual_path_to_cell_in_onering(const CH _ch_t, std::map<CH, HFH> &_pre_hf)
{
  std::vector<HFH> dpath;

  //get the path
  std::queue<HFH> que_dpath;
  que_dpath.push(_pre_hf[_ch_t]);
  while (!que_dpath.empty())
  {
    auto hf = que_dpath.front();
    que_dpath.pop();
    if (hf == HFH(-1))
      break;

    dpath.push_back(mesh_.opposite_halfface_handle(hf));
    que_dpath.push(_pre_hf[mesh_.incident_cell(hf)]);
  }

  std::reverse(dpath.begin(), dpath.end());


  return dpath;
}
}

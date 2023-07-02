/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define TRANSVERSALUNZIPPINGT_C

#include "TransversalUnzippingT.hh"

namespace AlgoHex
{
template<class MeshT>
std::vector<HEH>
TransversalUnzippingT<MeshT>::repair_fully_constrained_parabolic_sector(const VH _vh,
                                                                        const std::vector<HEH> &_sector_hehs,
                                                                        const std::set<FH> &_sec_fhs,
                                                                        bool _shrink_sector)
{
  auto eh0 = mesh_.edge_handle(_sector_hehs[0]);
  auto eh1 = mesh_.edge_handle(_sector_hehs.back());

  //if it's a zero sector between two singular edges, it'll be pushed inside
  if (feature_edge_[eh0] == 0 && feature_edge_[eh1] == 0)
    return std::vector<HEH>{};

  ALGOHEX_DEBUG_ONLY(std::cerr << "sector eh: " << eh0 << " " << eh1 << std::endl;)

  if (n_interior_complex_singular_edge(mesh_, valence_, _vh))
    return std::vector<HEH>{};

  if (_sector_hehs.size() <= 2)
  {
    std::cerr << "Error: no edge between two feature edges in the boundary zero-sector!" << std::endl;
    return std::vector<HEH>{};
  }

  //find a regular edge in the sector
  HEH he_s(-1);
  for (int i = 1; i < (int) _sector_hehs.size() - 1; ++i)
    if (valence_[mesh_.edge_handle(_sector_hehs[i])] == 0)
    {
      he_s = _sector_hehs[i];
      break;
    }

  if (!he_s.is_valid())
    return std::vector<HEH>{};

  int n_split = 0;

  VH vhf = mesh_.halfedge(he_s).from_vertex();
  VH vht = mesh_.halfedge(he_s).to_vertex();
  if (sgl_vt_[vht] != 0)
  {
    VH vh_new = es_.edge_split(mesh_.edge_handle(he_s));
    he_s = mesh_.find_halfedge(vhf, vh_new);

    n_split++;
  }

  //find a feature face
  std::vector<HFH> ft_hfhs;
  for (auto hehf_it = mesh_.hehf_iter(he_s); hehf_it.valid(); ++hehf_it)
  {
    if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
      ft_hfhs.push_back(*hehf_it);
  }

  if (ft_hfhs.size() != 2)
  {
    std::cerr << "Error: non-manifold feature surface!" << std::endl;
    return std::vector<HEH>{};
  }


  std::vector<HEH> outgoing_sg_hehs;
  for (const auto hfhi: ft_hfhs)
  {
    HEH sg_he = transversal_unzipping(he_s, hfhi, _shrink_sector);
    if (sg_he.is_valid())
      outgoing_sg_hehs.push_back(sg_he);
  }

  return outgoing_sg_hehs;
}


template<class MeshT>
HEH
TransversalUnzippingT<MeshT>::transversal_unzipping(const HEH _heh_s, const HFH _ft_hfh, const bool _shrink_sector)
{
  if (!_heh_s.is_valid())
    return HEH(-1);
  if (mesh_.is_boundary(_ft_hfh))
    return HEH(-1);

  int n_split = 0;
  //singular free start edge
  HEH he_s = _heh_s;
  HFH hf_s = _ft_hfh;

  //singular free start face
  HFH hf_adj = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_s, he_s));

  HEH he_open = mesh_.prev_halfedge_in_halfface(he_s, hf_adj);
  VH vhf_open = mesh_.halfedge(he_open).from_vertex();

  if (sgl_vt_[vhf_open] != 0)
  {
    VH vhm = fs_.face_split(mesh_.face_handle(hf_adj));
    n_split++;

    std::vector<VH> hfvhs;
    hfvhs.push_back(mesh_.halfedge(he_s).from_vertex());
    hfvhs.push_back(mesh_.halfedge(he_s).to_vertex());
    hfvhs.push_back(vhm);

    hf_adj = mesh_.find_halfface(hfvhs);
    he_open = mesh_.prev_halfedge_in_halfface(he_s, hf_adj);
  }


  //get rotational axis
  CH ch_s = mesh_.incident_cell(hf_s);
  auto nm = mesh_.normal(hf_s);
  int nm_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s], nm).second;


  //transition
  int trans_e_new = nm_ax + 4;

  if (_shrink_sector)
    trans_e_new = tq_.inverse_transition_idx(trans_e_new);

  std::set<HFH> surface{hf_adj};
  std::set<EH> bound_ehs;
  for (auto hfe_it = mesh_.hfe_iter(hf_adj); hfe_it.valid(); ++hfe_it)
    bound_ehs.insert(*hfe_it);

  ArcZippingT<MeshT>::zipping_surface(surface, std::set<EH>{}, he_open, trans_e_new);

  //update edge valence
  for (auto hfe_it = mesh_.hfe_iter(hf_adj); hfe_it.valid(); ++hfe_it)
    sge_.compute_edge_valence(*hfe_it);

  //update singular vertex
  std::set<VH> sgvhs;
  for (const auto &eh: bound_ehs)
  {
    sgvhs.insert(mesh_.edge(eh).from_vertex());
    sgvhs.insert(mesh_.edge(eh).to_vertex());
  }
  for (const auto &vh: sgvhs)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vh);

  QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                     feature_fprop_, feature_edge_, bound_ehs);

  return mesh_.opposite_halfedge_handle(he_open);
}

template<class MeshT>
void
TransversalUnzippingT<MeshT>::get_parabolic_feature_sectors_at_vertex(const VH _vh,
                                                                      std::vector<std::vector<HEH>> &_zero_sectors_hehs,
                                                                      std::vector<std::set<FH>> &_zero_sectors_fhs)
{
  if (!feature_edge_vertex_[_vh])
    return;

  if (n_interior_complex_singular_edge(mesh_, valence_, _vh) > 0)
    return;

  SplitHelperT<MeshT>::split_one_ring(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                      feature_edge_, sgl_vt_,
                                      feature_edge_vertex_);

  std::vector<HEH> ft_hehs, bdy_hehs, itr_hehs;
  get_special_halfedges_at_vertex(mesh_, feature_edge_, valence_, _vh, ft_hehs, bdy_hehs, itr_hehs);

  if (ft_hehs.size() <= 1)
    return;

  //get feature(singular) sectors on feature surface
  std::vector<std::vector<HEH>> v_sec_hehs;
  std::vector<std::set<FH>> v_sec_fhs;
  get_feature_sectors_at_feature_edge_vertex(mesh_, feature_fprop_, feature_edge_, feature_edge_vertex_, valence_,
                                             _vh, v_sec_hehs, v_sec_fhs);

  //store zero sectors
  for (auto j = 0u; j < v_sec_hehs.size(); ++j)
  {
    auto eh_start = mesh_.edge_handle(v_sec_hehs[j].front());
    auto eh_end = mesh_.edge_handle(v_sec_hehs[j].back());
    if (v_sec_hehs[j].size() < 2 || (feature_edge_[eh_start] == 0 && valence_[eh_start] == 0) ||
        (feature_edge_[eh_end] == 0 && valence_[eh_end] == 0))
      continue;

    if (is_parabolic_sector(v_sec_hehs[j], v_sec_fhs[j]))
    {
      _zero_sectors_hehs.push_back(v_sec_hehs[j]);
      _zero_sectors_fhs.push_back(v_sec_fhs[j]);
    }
  }

  //output info
  if (!_zero_sectors_hehs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Feature vertex " << _vh << ": " << std::endl;)
    for (auto j = 0u; j < v_sec_hehs.size(); ++j)
    {
      auto eh_start = mesh_.edge_handle(v_sec_hehs[j].front());
      auto eh_end = mesh_.edge_handle(v_sec_hehs[j].back());
      if (v_sec_hehs[j].size() < 2 || (feature_edge_[eh_start] == 0 && valence_[eh_start] == 0) ||
          (feature_edge_[eh_end] == 0 && valence_[eh_end] == 0))
        continue;

      int ea_st = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_,
                                                                                     cell_quaternions_,
                                                                                     trans_prop_, valence_,
                                                                                     v_sec_hehs[j],
                                                                                     v_sec_fhs[j]);

      ALGOHEX_DEBUG_ONLY(std::cerr << " sector: " << mesh_.edge_handle(v_sec_hehs[j].front()) << " "
                                   << mesh_.edge_handle(v_sec_hehs[j].back())
                                   << " ag " << ea_st << std::endl;)
    }
  }
}

template<class MeshT>
bool
TransversalUnzippingT<MeshT>::is_parabolic_sector(const std::vector<HEH> &_zero_sector_hehs,
                                                  const std::set<FH> &_zero_sectors_fhs) const
{
  if (_zero_sector_hehs.empty())
  {
    return false;
  }

  for (const auto hehi: _zero_sector_hehs)
  {
    if (mesh_.is_deleted(hehi))
      return false;
  }

  int ea_st = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_,
                                                                                 trans_prop_, valence_,
                                                                                 _zero_sector_hehs, _zero_sectors_fhs);
  if (ea_st == 0)
    return true;

  return false;
}
}
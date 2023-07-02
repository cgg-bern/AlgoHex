/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define DUPLICATEONERINGMESHT_C

#include "DuplicateOneRingMeshT.hh"
#include "EdgeMonodromyHelperT.hh"

namespace AlgoHex
{
template<class MeshT>
void DuplicateOneRingMeshT<MeshT>::copy_mesh(const std::vector<VH> &_vhs)
{
  std::unordered_set<int> onering_vhs;
  for (const auto vh: _vhs)
  {
    if (onering_vhs.count(vh.idx()) == 0)
      onering_vhs.insert(vh.idx());
    for (auto vv_it = fmesh_.vv_iter(vh); vv_it.valid(); ++vv_it)
      if (onering_vhs.count((*vv_it).idx()) == 0)
        onering_vhs.insert((*vv_it).idx());
  }

  std::unordered_set<int> onering_chs;
  for (const auto vh: _vhs)
  {
    for (auto vc_it = fmesh_.vc_iter(vh); vc_it.valid(); ++vc_it)
      if (onering_chs.count((*vc_it).idx()) == 0)
        onering_chs.insert((*vc_it).idx());
  }


  std::map<VH, VH> fvh_tvh;
  //add vertices
  for (const auto vh: onering_vhs)
  {
    fvh_tvh[VH(vh)] = tmesh_.add_vertex(fmesh_.vertex(VH(vh)));
    tvh_fvh_[fvh_tvh[VH(vh)]] = VH(vh);
  }


  //add cells
  for (const auto ch: onering_chs)
  {
    std::vector<VH> t_cvhs;
    auto f_cvhs = fmesh_.get_cell_vertices(CH(ch));
    for (auto cvh: f_cvhs)
      t_cvhs.push_back(fvh_tvh[cvh]);

    tch_fch_[tmesh_.add_cell(t_cvhs)] = CH(ch);
  }

  //set properties
  //halfedge
  for (const auto heh: tmesh_.halfedges())
  {
    auto fheh = fmesh_.find_halfedge(tvh_fvh_[tmesh_.halfedge(heh).from_vertex()],
                                     tvh_fvh_[tmesh_.halfedge(heh).to_vertex()]);
    theh_fheh_[heh] = fheh;
  }
  tmesh_.set_persistent(theh_fheh_, true);

  //halfface
  for (const auto hfh: tmesh_.halffaces())
  {
    std::vector<VH> hfvs;
    for (auto hfv_it = tmesh_.hfv_iter(hfh); hfv_it.valid(); ++hfv_it)
      hfvs.push_back(tvh_fvh_[*hfv_it]);
    thfh_fhfh_[hfh] = fmesh_.find_halfface(hfvs);
  }
}


template<class MeshT>
void DuplicateOneRingMeshT<MeshT>::copy_properties()
{
//        for(auto cprop_it = fmesh_.cell_props_begin(); cprop_it != fmesh_.cell_props_end(); ++cprop_it) {
//            auto name =  (*cprop_it)->name();
//            tmesh_. template request_cell_property<decltype(*(*cprop_it)->begin())>(name);
//        }

  tmesh_.set_persistent(cell_qts_1r_, true);
  tmesh_.set_persistent(trans_prop_1r_, true);
  tmesh_.set_persistent(valence_1r_, true);

  tmesh_.set_persistent(feature_fprop_1r_, true);
  tmesh_.set_persistent(feature_eprop_1r_, true);


  auto cell_quaternionsf = fmesh_.template request_cell_property<Eigen::Quaterniond>("FrameFieldQuaternions");
  auto trans_propf = fmesh_.template request_halfface_property<int>("HalffaceTransiton");
  auto valencef = fmesh_.template request_edge_property<int>("edge_valance");
  auto feature_edgef = fmesh_.template request_edge_property<int>("AlgoHex::FeatureEdges");
  auto feature_facef = fmesh_.template request_face_property<int>("AlgoHex::FeatureFaces");


  for (auto ch: tmesh_.cells())
  {
    auto fch = tch_fch_[ch];
    cell_qts_1r_[ch] = cell_quaternionsf[fch];
  }

  for (auto hfh: tmesh_.halffaces())
  {
    if (tmesh_.is_boundary(tmesh_.face_handle(hfh)))
      trans_prop_1r_[hfh] = 0;
    else
    {
      auto fhfh = thfh_fhfh_[hfh];
      trans_prop_1r_[hfh] = trans_propf[fhfh];
    }
  }

  for (auto eh: tmesh_.edges())
  {
    auto fheh = theh_fheh_[tmesh_.halfedge_handle(eh, 0)];
    valence_1r_[eh] = valencef[fmesh_.edge_handle(fheh)];
    feature_eprop_1r_[eh] = feature_edgef[fmesh_.edge_handle(fheh)];
  }

  for (auto fh: tmesh_.faces())
  {
    auto ffh = fmesh_.face_handle(thfh_fhfh_[tmesh_.halfface_handle(fh, 0)]);
    feature_fprop_1r_[fh] = feature_facef[ffh];
  }
}

template<class MeshT>
bool DuplicateOneRingMeshT<MeshT>::check_cell_singularity() const
{
  //check if a cell contains > 1 singular edges
  for (auto ch: tmesh_.cells())
  {
    int count = 0;
    for (auto eh: tmesh_.cell_edges(ch))
      if (valence_1r_[eh] != 0)
        count++;

    if (count > 1)
    {
      std::cerr << "\nError: cell " << ch << " contains " << count << " singular edges" << std::endl;
      return false;
    }
  }

  return true;
}

template<class MeshT>
VH DuplicateOneRingMeshT<MeshT>::onering_mesh_vertex_handle(const VH _fvh) const
{
  for (const auto vh: tmesh_.vertices())
  {
    if (tvh_fvh_[vh] == _fvh)
      return vh;
  }

  return VH(-1);
}

template<class MeshT>
HEH DuplicateOneRingMeshT<MeshT>::onering_mesh_halfedge_handle(const HEH _fheh) const
{
  for (const auto heh: tmesh_.halfedges())
  {
    if (theh_fheh_[heh] == _fheh)
      return heh;
  }

  return HEH(-1);
}

template<class MeshT>
HFH DuplicateOneRingMeshT<MeshT>::original_halfface_handle(const HFH _hfh) const
{
  return thfh_fhfh_[_hfh];
}

template<class MeshT>
EH DuplicateOneRingMeshT<MeshT>::original_edge_handle(const EH _eh) const
{
  return fmesh_.edge_handle(theh_fheh_[tmesh_.halfedge_handle(_eh, 0)]);
}


//there should be only one singular (feature) edge in one cell
template<class MeshT>
void DuplicateOneRingMeshT<MeshT>::get_halfedge_axis_in_cell_in_onering_mesh()
{
  TransitionQuaternion tq;
  auto cell_quaternionsf = fmesh_.template request_cell_property<Eigen::Quaterniond>("FrameFieldQuaternions");
  auto trans_propf = fmesh_.template request_halfface_property<int>("HalffaceTransiton");
  auto valencef = fmesh_.template request_edge_property<int>("edge_valance");
  auto feature_edgef = fmesh_.template request_edge_property<int>("AlgoHex::FeatureEdges");

  auto che_axis_1r = tmesh_.template request_cell_property<std::pair<HEH, int> >("EdgeAxis",
                                                                                 std::make_pair(HEH(-1), -1));
  tmesh_.set_persistent(che_axis_1r, true);

  //cell halfedg axis in original mesh
  std::set<EH> eset;
  for (const auto chi: tmesh_.cells())
    for (auto ce_it = tmesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
      if ((valence_1r_[*ce_it] >= -2 && valence_1r_[*ce_it] <= 4 && valence_1r_[*ce_it] != 0)
          || feature_eprop_1r_[*ce_it] > 0)
        eset.insert(*ce_it);

  //store cell halfedge axis in original mesh
  std::map<CH, std::pair<HEH, int>> c_hax;
  for (const auto ehi: eset)
  {
    HEH heh0 = tmesh_.halfedge_handle(ehi, 0);
    HEH heh0_org = theh_fheh_[heh0];

    std::map<CH, std::vector<VecAxis >> cell_alignments;
    if ((valence_1r_[ehi] >= -2 && valence_1r_[ehi] <= 4 && valence_1r_[ehi] != 0) || feature_eprop_1r_[ehi] > 0)
    {
      //align to edge
      auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells(fmesh_, tq, cell_quaternionsf,
                                                                                     trans_propf, valencef,
                                                                                     feature_edgef, heh0_org);
//      auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells_wrt_valence(fmesh_, tq, cell_quaternionsf, trans_propf, valencef, feature_edgef, heh0_org);

      for (auto&[chi, axi]: cell_edge_axis)
      {
//                    cell_alignments[chi].emplace_back(eigen_dir, axi);
        c_hax[chi] = std::make_pair(heh0_org, axi);;
      }
    }
  }

  for (const auto chi: tmesh_.cells())
  {
    for (auto ce_it = tmesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
    {
      if ((valence_1r_[*ce_it] >= -2 && valence_1r_[*ce_it] <= 4 && valence_1r_[*ce_it] != 0) ||
          feature_eprop_1r_[*ce_it] > 0)
      {
        HEH heh0 = tmesh_.halfedge_handle(*ce_it, 0);

        HEH heh0_org = theh_fheh_[heh0];
        CH ch_org = tch_fch_[chi];

        if (c_hax[ch_org].first == heh0_org)
          che_axis_1r[chi] = std::make_pair(heh0, c_hax[ch_org].second);
        else
          che_axis_1r[chi] = std::make_pair(tmesh_.opposite_halfedge_handle(heh0), c_hax[ch_org].second);

        break;
      }
    }
  }
}

}

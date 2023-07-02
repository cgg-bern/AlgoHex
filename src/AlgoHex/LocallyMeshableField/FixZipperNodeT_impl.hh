/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define FIXTURNINGPOINTT_C

#include "FixZipperNodeT.hh"
#include "CellInfo.hh"
//#include "TetmeshOperationsT.hh"
#include "AlgoHex/BinarySpacePartitionTrees/PlaneT.hh"
#include "CommonFuncs.hh"
#include "QuaternionsSmoothing.hh"


namespace AlgoHex
{
template<class MeshT>
bool
FixZipperNodeT<MeshT>::
unzip_zipper_node(const VH _vh, const HEH _heh_s, const std::set<VH> &_fixable_vhs, VH &_other_tp)
{
  SplitHelperT<MeshT>::split_one_ring(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                      feature_edge_, sgl_vt_, feature_edge_vertex_);

  HEH he_s;
  HFH hf_s;
  CH ch_s;
  if (!_heh_s.is_valid())
    ch_s = find_start_cell(_vh, he_s, hf_s);
  else
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "heh start is provided " << mesh_.edge_handle(_heh_s) << std::endl;)
    he_s = _heh_s;
    ch_s = find_start_cell_at_halfedge(_vh, he_s, hf_s).second;
  }

  if (!ch_s.is_valid() || mesh_.is_deleted(ch_s))
  {
    return false;
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "dpath weight: " << dpath_weight_ << " eh_s: " << mesh_.edge_handle(he_s) << " hf_s: "
                               << mesh_.face_handle(hf_s) << " ch_s: " << ch_s << std::endl;)

  auto path0 = shortest_dual_path_to_boundary(he_s, hf_s, _fixable_vhs, _other_tp, dpath_weight_);

  auto path1 = get_handle_free_dual_path(path0, he_s, _other_tp);

  auto path2 = get_special_vertex_free_path(path1, he_s, _other_tp);

  if (path2.empty())
    return false;

  //find edge paths
  std::set<EH> eset;
  //exclude edges from edge path
  std::set<EH> excluded_ehs;
  auto eh_s = mesh_.edge_handle(he_s);
  excluded_ehs.insert(eh_s);

  for (auto ce_it = mesh_.ce_iter(path2.back()); ce_it.valid(); ++ce_it)
    if (mesh_.is_boundary(*ce_it))
      excluded_ehs.insert(*ce_it);

  //target vertices
  std::set<VH> target_vhs;
  if (_other_tp.is_valid())
    target_vhs.insert(_other_tp);
  else
  {
    for (auto cv_it = mesh_.cv_iter(path2.back()); cv_it.valid(); ++cv_it)
      if (mesh_.is_boundary(*cv_it))
        target_vhs.insert(*cv_it);
  }


  auto vh_f = mesh_.halfedge(he_s).from_vertex();
  auto vh_t = mesh_.halfedge(he_s).to_vertex();

  //first edge path
  std::vector<EH> ep1_r;
  if (!feature_face_vertex_[_vh])
    ep1_r = get_edge_path_on_cell_path_wrt_field(path2, excluded_ehs, vh_f, he_s, target_vhs, false);
  else //only one incident singular edge
    ep1_r = get_edge_path_on_cell_path(path2, excluded_ehs, vh_f, target_vhs);

  ALGOHEX_DEBUG_ONLY(std::cerr << "raw epath: ";
                             for (const auto &eh: ep1_r)
                               std::cerr << " " << eh;
                             std::cerr << std::endl;)

  auto ep1 = remedy_edge_path(path2, ep1_r);

  ALGOHEX_DEBUG_ONLY(std::cerr << "after remedy epath: ";
                             for (const auto &eh: ep1)
                               std::cerr << " " << eh;
                             std::cerr << std::endl;)

  if (ep1.empty())
  {
    std::cout << "Error: couldn't find the first edge path!" << std::endl;
    return false;
  }

  //exclude the first edge path
  excluded_ehs.insert(ep1.begin(), ep1.end());

  //exclude the last vertex on ep1 if ends on boundary
  if (!_other_tp.is_valid())
  {
    target_vhs.erase(mesh_.edge(ep1_r.back()).from_vertex());
    target_vhs.erase(mesh_.edge(ep1_r.back()).to_vertex());
  }

  std::set<CH> dp_cells;
  int n_split = split_for_second_edge_path(path2, ep1, mesh_.edge_handle(he_s), dp_cells);


  //get bounding surface of the stripe
  std::set<HFH> all_bound_hfs, all_hfs;
  for (const auto ch: dp_cells)
    for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
    {
      all_hfs.insert(mesh_.opposite_halfface_handle(*chf_it));
    }
  //bounding halffaces (no boundary halfface)
  for (const auto hf: all_hfs)
  {
    if (dp_cells.find(mesh_.incident_cell(hf)) == dp_cells.end() && !mesh_.is_boundary(hf))
      all_bound_hfs.insert(hf);
  }

  std::set<EH> all_bound_ehs;
  for (const auto hf: all_bound_hfs)
    for (auto hfe_it = mesh_.hfe_iter(hf); hfe_it.valid(); ++hfe_it)
      all_bound_ehs.insert(*hfe_it);

  std::set<EH> excl_ehs(ep1.begin(), ep1.end());

  auto ep2 = shortest_path_v_to_target(all_bound_ehs, excl_ehs, vh_t, target_vhs, _other_tp);


  if (ep2.empty())
  {
    std::cout << "Error: couldn't find the second edge path!" << std::endl;
    return false;
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << " e path2: ";
                             for (const auto eh: ep2)
                               std::cerr << " " << eh;
                             std::cerr << std::endl;)

  excluded_ehs.insert(ep2.begin(), ep2.end());

  HEH he_surf = he_s;
  auto surface = bounded_surface(all_bound_hfs, excluded_ehs, he_surf);

  if (surface.empty())
  {
    //try the other side
    he_surf = mesh_.opposite_halfedge_handle(he_s);
    surface = bounded_surface(all_bound_hfs, excluded_ehs, he_surf);

    if (surface.empty())
    {
      std::cout << "Error: couldn't find the surface!" << std::endl;
      return false;
    }
  }

//        std::cerr << "Check alignment before matching adj" << std::endl;
//  check_field_alignment(_vh, excluded_ehs, surface, false);

  bool suc = ArcZippingT<MeshT>::zipping_surface(surface, excluded_ehs, he_surf, 0);
  if (suc)
  {
    //block surface
    for (const auto &hf: surface)
    {
      visited_fprop_[mesh_.face_handle(hf)] = true;
    }

    //update edge valence
    for (const auto &eh: excluded_ehs)
    {
      if (!mesh_.is_deleted(eh))
        sge_.compute_edge_valence(eh);
      else
        std::cerr << "Warning: edge " << eh << " on the path is deleted" << std::endl;
    }

    QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                       feature_fprop_, feature_edge_,
                                                       excluded_ehs);

    //update singular vertex
    std::set<VH> sgvhs;
    for (const auto &eh: excluded_ehs)
    {
      sgvhs.insert(mesh_.edge(eh).from_vertex());
      sgvhs.insert(mesh_.edge(eh).to_vertex());
    }
    for (const auto &vh: sgvhs)
      MeshPropertiesT<MeshT>::update_singular_vertex_property(vh);

    //post process (push boundary vertices inside the volume if there are)
    VertexOptimizeT<MeshT> vo(mesh_);
    for (const auto &vh: sgvhs)
    {
      if (sgl_vt_[vh] && feature_face_vertex_[vh] && vh != vh_f &&
          n_incident_singular_edges(mesh_, valence_, vh) > 1)
      {
        if (mesh_.is_boundary(vh))
        {
          psv_.push_singular_vertex(vh);
        }
      }
    }

    if (_other_tp.is_valid())
    {
      std::vector<HEH> end_hehs;
      for (auto voh_it = mesh_.voh_iter(_other_tp); voh_it.valid(); ++voh_it)
      {
        auto ve = mesh_.edge_handle(*voh_it);
        if (excluded_ehs.find(ve) != excluded_ehs.end())
          end_hehs.push_back(*voh_it);
      }
      if (end_hehs.size() != 2)
        std::cerr << "Error: incident path edges number at the other tp is not 2" << std::endl;

      bool suc = fis_.detach_singular_arc_from_singular_node(_other_tp, end_hehs[0], end_hehs, false, true);
      if (!suc)
        fis_.detach_singular_arcs_of_zipper_node(_other_tp, end_hehs[0], end_hehs, false, true);
    }

    //update singular edge pairs
    //in case a singular edge becomes regular
    for (const auto &eh: excluded_ehs)
      if (valence_[eh] == 0)
        sg_edge_pairs_[eh] = 0;

    //update pair id
    auto iter1 = std::find_if(ep1.begin(), ep1.end(),
                              [&](const EH &edge) { return valence_[edge] == -1 || valence_[edge] == 1; });
    auto iter2 = std::find_if(ep2.begin(), ep2.end(),
                              [&](const EH &edge) { return valence_[edge] == -1 || valence_[edge] == 1; });

    if (iter1 != ep1.end() && iter2 != ep2.end())
    {
      EH eh1 = *iter1;
      auto sgehs1 = get_edges_on_singular_arc(mesh_, valence_, mesh_.halfedge_handle(eh1, 0));
      ALGOHEX_DEBUG_ONLY(std::cerr << " eh1 " << eh1 << "sgehs1 size " << sgehs1.size();)

      EH eh2 = *iter2;
      auto sgehs2 = get_edges_on_singular_arc(mesh_, valence_, mesh_.halfedge_handle(eh2, 0));
      ALGOHEX_DEBUG_ONLY(std::cerr << " eh2 " << eh2 << " sgehs2 size " << sgehs2.size() << std::endl;)

      ALGOHEX_DEBUG_ONLY(std::cerr << "update pairs" << std::endl;)
      int cur_id = 0;
      for (const auto &ehi: sgehs1)
      {
        if (sg_edge_pairs_[ehi] > 0)
        {
          cur_id = sg_edge_pairs_[ehi];
          break;
        }
      }


      int max_pair_id = *std::max_element(sg_edge_pairs_.begin(), sg_edge_pairs_.end()) + 1;

      if (cur_id == 0)
        cur_id = max_pair_id;

      for (const auto &ehi: sgehs1)
        if (valence_[ehi] != 0)
        {
          sg_edge_pairs_[ehi] = cur_id;
        }

      cur_id++;
      for (const auto &ehi: sgehs2)
        if (valence_[ehi] != 0)
        {
          sg_edge_pairs_[ehi] = cur_id;
        }

      //relocate sg vhs
//                double scale = mo_.get_scale();
//                mo_.remesh(cur_id-1, 1. / scale, 1. / scale, 1. / scale, 1. / scale, 1. / scale, 2);
    }
    mesh_.collect_garbage();
    return true;
  }

  mesh_.collect_garbage();
  return false;
}

template<class MeshT>
bool
FixZipperNodeT<MeshT>::zip_zipper_node(const VH _vh)
{
  std::set<EH> sg_ehs;
  if (!is_zipable(_vh, sg_ehs))
    return false;

  EH eh0 = *sg_ehs.begin();
  EH eh1 = *sg_ehs.rbegin();

  VH vh_s = mesh_.edge(eh0).from_vertex() == _vh ? mesh_.edge(eh0).to_vertex() : mesh_.edge(eh0).from_vertex();
  VH vh_t = mesh_.edge(eh1).from_vertex() == _vh ? mesh_.edge(eh1).to_vertex() : mesh_.edge(eh1).from_vertex();

  auto pvhs = ArcZippingT<MeshT>::shortest_path_vertex_to_vertex_on_boundary(vh_s, vh_t);

  //if there is feature edge vertex on the path, skip
  int n_fev = 0;
  for (auto i = 1u; i < pvhs.size() - 1; ++i)
    if (feature_edge_vertex_[pvhs[i]])
      n_fev++;
  if (n_fev >= 1)
    return false;

  std::set<EH> bdy_ehs;
  for (auto i = 0u; i < pvhs.size() - 1; ++i)
  {
    bdy_ehs.insert(mesh_.edge_handle(mesh_.find_halfedge(pvhs[i], pvhs[i + 1])));
  }
  auto cutfhs = ArcZippingT<MeshT>::generate_cut_surface_with_bounding_edges(bdy_ehs, *sg_ehs.begin(), *sg_ehs.rbegin(),
                                                                             _vh);

  return fis_.detaching_singularity(_vh, cutfhs, sg_ehs);
}

template<class MeshT>
bool
FixZipperNodeT<MeshT>::
merge_zipper_node(const VH _vh)
{
  std::vector<HEH> hehs_clps;
  auto mcase = is_mergable(_vh, hehs_clps);
  if (!mcase)
    return false;

  std::set<CH> onering_chs;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    onering_chs.insert(*vc_it);

  std::set<EH> split_ehs;
  std::vector<HEH> tp_sghehs;
  HEH heha, hehb;
  double min_len = 10000000000.;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto val = valence_[mesh_.edge_handle(*voh_it)];
    if (val == 1 || val == -1)
    {
      double elen = mesh_.length(*voh_it);
      if (elen < min_len)
        min_len = elen;

      tp_sghehs.push_back(*voh_it);

      if (val == 1)
        heha = *voh_it;
    }
  }

  Point pttp0 = mesh_.vertex(mesh_.halfedge(tp_sghehs[0]).to_vertex());
  Point pttp1 = mesh_.vertex(mesh_.halfedge(tp_sghehs[1]).to_vertex());
  Point pt0 = mesh_.vertex(_vh);
  Point pt1 = pt0 + 100 * (pttp0 - pt0);
  Point pt2 = pt0 + 100 * (pttp1 - pt0);

  ACG::Geometry::PlaneT<double> plane(pt0, pt1, pt2);

  for (const auto hehi: hehs_clps)
  {
    VH vhb = mesh_.halfedge(hehi).to_vertex();


    for (auto voh_it = mesh_.voh_iter(vhb); voh_it.valid(); ++voh_it)
    {
      auto val = valence_[mesh_.edge_handle(*voh_it)];
      if (val == 1)
      {
        hehb = *voh_it;
        break;
      }
    }

    if (!heha.is_valid() || !hehb.is_valid())
      return false;

    //if they are on the same arc, return false
    auto sg_vhs1 = get_vertices_on_singular_arc(mesh_, valence_, heha);
    auto result1 = std::find(sg_vhs1.begin(), sg_vhs1.end(), vhb);
    if (result1 != std::end(sg_vhs1))
      return false;

    //if two arcs share common vertex, return false
    auto sg_vhs2 = get_vertices_on_singular_arc(mesh_, valence_, hehb);
    if (sg_vhs1[0] == sg_vhs2[0] || sg_vhs1[0] == sg_vhs2.back()
        || sg_vhs1.back() == sg_vhs2[0] || sg_vhs1.back() == sg_vhs2.back())
      return false;


    //check rt axis
    int rtax_a = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, heha,
                                                                    *mesh_.hehf_iter(heha));
    int rtax_b = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, hehb,
                                                                    *mesh_.hehf_iter(hehb));
//std::cerr<<"rta: "<<rtax_a<<" rtb: "<<rtax_b<<std::endl;
    CH ch_common = *mesh_.hec_iter(hehi);

    int rt_ax_a_in_com = fac_.axis_in_chs_expressed_in_cht(
            mesh_.incident_cell(mesh_.opposite_halfface_handle(*mesh_.hehf_iter(heha))),
            ch_common, rtax_a, onering_chs);


    std::set<CH> onering_b_chs;
    for (auto vc_it = mesh_.vc_iter(vhb); vc_it.valid(); ++vc_it)
      onering_b_chs.insert(*vc_it);

    int rt_ax_b_in_com = fac_.axis_in_chs_expressed_in_cht(
            mesh_.incident_cell(mesh_.opposite_halfface_handle(*mesh_.hehf_iter(hehb))),
            ch_common, rtax_b, onering_b_chs);

//        std::cerr<<"    rta common: "<<rt_ax_a_in_com<<" rtb common: "<<rt_ax_b_in_com<<" case: "<<mcase;

    //sg arc parallel to zipper node edges
    if (rt_ax_a_in_com / 2 == rt_ax_b_in_com / 2)
      return false;


    bool ints_arc = false;
    auto sg_hehs = get_halfedges_on_singular_arc(mesh_, valence_, hehb);
    for (const auto &sg_heh: sg_hehs)
    {
      if (intersect_with_triangle(plane, pt0, pt1, pt2,
                                  mesh_.vertex(mesh_.halfedge(sg_heh).from_vertex()),
                                  mesh_.vertex(mesh_.halfedge(sg_heh).to_vertex())))
      {
        ints_arc = true;
//                    std::cerr<<"Intersected plane "<<pt0<<" "<<pt1<<" "<<pt2<<std::endl;
//                    std::cerr<<"Intersected seg "<<mesh_.vertex(mesh_.halfedge(sg_heh).from_vertex())<<" "<<mesh_.vertex(mesh_.halfedge(sg_heh).to_vertex())<<std::endl;

        break;
      }
    }

    if (!ints_arc)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "No intersection with zipper node sector" << std::endl;)
      return false;
    }

//        std::cerr<<" collapse ok: "<<(ec_.is_topology_ok(heh_clps)== ec_.CollapseOK)<<" "<<std::endl;

    VH vh_new(-1);
    if (ec_.is_topology_ok(hehi) == ec_.CollapseOK)
      vh_new = ec_.edge_collapse(hehi, true);

    if (!vh_new.is_valid())
    {
      HEH heh_opp = mesh_.opposite_halfedge_handle(hehi);
      if (ec_.is_topology_ok(heh_opp) == ec_.CollapseOK)
        vh_new = ec_.edge_collapse(heh_opp, true);
    }

    if (vh_new.is_valid())
    {
      //update edge length and push
      for (auto ve_it = mesh_.ve_iter(vh_new); ve_it.valid(); ++ve_it)
      {
        if (!mesh_.is_deleted(*ve_it))
        {
          //update valence
          if (!mesh_.is_boundary(*ve_it) || valence_[*ve_it] != 0)
            sge_.compute_edge_valence(*ve_it);
        }
      }

      for (auto vv_it = mesh_.vv_iter(vh_new); vv_it.valid(); ++vv_it)
        MeshPropertiesT<MeshT>::update_singular_vertex_property(*vv_it);
      MeshPropertiesT<MeshT>::update_singular_vertex_property(vh_new);

      return true;
    }
    else
    {
      for (auto voh_it = mesh_.voh_iter(vhb); voh_it.valid(); ++voh_it)
      {
        auto val = valence_[mesh_.edge_handle(*voh_it)];
        if (val == 1 && mesh_.length(*voh_it) > min_len)
        {
          split_ehs.insert(mesh_.edge_handle(*voh_it));
        }
      }
    }
  }

  for (const auto ehi: split_ehs)
  {
    if (es_.is_split_ok(ehi))
      es_.edge_split(ehi);
  }

  return false;
}


template<class MeshT>
bool
FixZipperNodeT<MeshT>::is_zipable(const VH _vh, std::set<EH> &_sg_ehs) const
{
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh_.edge_handle(*voh_it);
    if (valence_[ve] == -1 || valence_[ve] == 1)
    {
      VH vh_t = mesh_.halfedge(*voh_it).to_vertex();
      if (!mesh_.is_boundary(vh_t))
        return false;

      if (feature_node_[vh_t])
        return false;

      _sg_ehs.insert(ve);
    }
  }

  if (_sg_ehs.size() == 2)
  {
    if ((valence_[*_sg_ehs.begin()] == 1 && valence_[*(++_sg_ehs.begin())] == -1) ||
        (valence_[*_sg_ehs.begin()] == -1 && valence_[*(++_sg_ehs.begin())] == 1))
      return true;
  }

  return false;
}


template<class MeshT>
bool
FixZipperNodeT<MeshT>::
is_mergable(const VH _vh, std::vector<HEH> &_v_hehs) const
{
  if (mesh_.is_deleted(_vh))
    return false;

  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto veh = mesh_.edge_handle(*voh_it);
    if (valence_[veh] == 0)
    {
      auto vht = mesh_.halfedge(*voh_it).to_vertex();
      if (sgl_vt_[vht] == 1)
      {
        int n_val_1 = 0, n_val1 = 0, n_all = 0;
        for (auto ve_it = mesh_.ve_iter(vht); ve_it.valid(); ++ve_it)
        {
          if (valence_[*ve_it] == -1)
            n_val_1++;
          if (valence_[*ve_it] == 1)
            n_val1++;
          if (valence_[*ve_it] != 0)
            n_all++;
        }

        if (n_all == 2)
        {
          if (n_val1 == 2)
          {
            _v_hehs.push_back(*voh_it);
          }
        }
      }
    }
  }

  if (!_v_hehs.empty())
    return true;

  return false;
}

template<class MeshT>
bool
FixZipperNodeT<MeshT>::
intersect_with_triangle(const ACG::Geometry::PlaneT<double> &_plane, const Point &_p0, const Point &_p1,
                        const Point &_p2, const Point &_pa, const Point &_pb) const
{
  Point pt;
  double ratio = 0.;
  bool ints = _plane.intersect_linesegment(_pa, _pb, pt, ratio);
  if (ints)
  {
    auto bcs = compute_bary_coord(_p0, _p1, _p2, pt);
    for (auto bc: bcs)
      if (bc > 1 || bc < 0)
        return false;

//            std::cerr<<"bc "<<bcs[0]<<" "<<bcs[1]<<" "<<bcs[2]<<std::endl;
//            std::cerr<<"Intersected point "<<pt<<std::endl;

    return true;
  }

  return false;
}

template<class MeshT>
typename MeshT::PointT
FixZipperNodeT<MeshT>::compute_bary_coord(const Point &_p0, const Point &_p1, const Point &_p2, const Point &_px) const
{
  Point v0 = _p1 - _p0, v1 = _p2 - _p0, v2 = _px - _p0;
  auto d00 = v0.dot(v0);
  auto d01 = v0.dot(v1);
  auto d11 = v1.dot(v1);
  auto d20 = v2.dot(v0);
  auto d21 = v2.dot(v1);
  auto denom = d00 * d11 - d01 * d01;
  auto alpha = (d11 * d20 - d01 * d21) / denom;
  auto beta = (d00 * d21 - d01 * d20) / denom;
  auto gama = 1.0 - alpha - beta;

  return Point(alpha, beta, gama);
}


template<class MeshT>
bool
FixZipperNodeT<MeshT>::
check_field_alignment(VH _vh, const std::set<EH> &_sg_ehs, const std::set<HFH> &_bound_hfhs, const bool _save) const
{
  std::set<FH> ffhs;
  std::set<VH> sgvhs;

  for (auto ehi: _sg_ehs)
  {
    sgvhs.insert(mesh_.edge(ehi).from_vertex());
    sgvhs.insert(mesh_.edge(ehi).to_vertex());
  }

  for (auto vhi: sgvhs)
  {
    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
    {
      for (auto cf_it = mesh_.cf_iter(*vc_it); cf_it.valid(); ++cf_it)
        if (feature_fprop_[*cf_it] > 0)
          ffhs.insert(*cf_it);
    }
  }

  MeshT tmp0;

  for (auto fhi: ffhs)
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
      std::cerr << "Error: feature face alignment wrong! ch0 " << ch0 << " ch1 " << ch1 << " fh " << fhi
                << " ax0 " << axis0 << " ax1 " << axis1 << " ax1in0 " << ax1in0 << std::endl;

      //debug
      if (_save)
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

  return true;
}


template<class MeshT>
bool
FixZipperNodeT<MeshT>::is_fixable_constrained_zipper_node(const VH _vh, HEH &_sg_heh)
{
  if (sgl_vt_[_vh] != 0 && feature_face_vertex_[_vh])
  {
    int n_ff_spe = n_incident_special_edges_on_feature_face(mesh_, feature_face_edge_, feature_edge_, valence_, _vh);
    if (n_ff_spe >= 3 || feature_node_[_vh])
    {
      std::vector<HEH> sg_hehs;
      for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
      {
        EH ve = mesh_.edge_handle(*voh_it);
        if (valence_[ve] == -1 && !feature_face_edge_[ve])
        {
          sg_hehs.push_back(*voh_it);
        }
      }

      for (const auto hehi: sg_hehs)
      {
        int ctp_type = fctp_.fixable_constrained_zipper_node_type(hehi);
        if (ctp_type == 5 || ctp_type == 2)
        {
          _sg_heh = hehi;

          return true;
        }
      }
    }
  }

  return false;
}

template<class MeshT>
bool
FixZipperNodeT<MeshT>::is_zipper_node(const VH _vh, const bool _fix_on_circle) const
{
  std::vector<HEH> hehs;
  int n_sge = 0;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    EH eh = mesh_.edge_handle(*voh_it);
    if (valence_[eh] == 1 || valence_[eh] == -1)
      hehs.push_back(*voh_it);

    if (valence_[eh] != 0)
      n_sge++;
  }
  if (n_sge != (int) hehs.size())
    return false;

  if (hehs.size() == 1)
  {
    if (!mesh_.is_boundary(_vh))
      return false;

    HFH hfh_s = *mesh_.hehf_iter(hehs[0]);
    CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh_s));

    //transition axis of sg edge
    int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                     valence_, hehs[0], hfh_s);

    //halfedge direction
    int e_aixs = rt_axis;
    if (valence_[mesh_.edge_handle(hehs[0])] == 1)
      e_aixs = negate(rt_axis);


    //find boundary hf
    HFH hfh_t(-1);
    for (auto vhf_it = mesh_.vhf_iter(_vh); vhf_it.valid(); ++vhf_it)
    {
      if (mesh_.is_boundary(*vhf_it))
      {
        hfh_t = *vhf_it;
        break;
      }
    }

    if (hfh_t.is_valid())
    {
      CH ch_t = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh_t));
      Point nm = mesh_.normal(hfh_t);
      int nm_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_t],
                                                     nm.normalize()).second;

      std::set<CH> onering_chs;
      for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
        onering_chs.insert(*vc_it);
      int e_ax_in_t = fac_.axis_in_chs_expressed_in_cht(ch_s, ch_t, e_aixs, onering_chs);

      if (e_ax_in_t == nm_ax)
        return true;
    }
  }
  else if (hehs.size() == 2)
  {
    //should be of different valence
    if (valence_[mesh_.edge_handle(hehs[0])] == valence_[mesh_.edge_handle(hehs[1])])
      return false;

    for (const auto &heh: hehs)
      if (mesh_.is_boundary(heh))
        return false;

    //it's a circle
    if (!_fix_on_circle && is_singular_circle(hehs[0]))
      return false;

    return true;
  }

  return false;
}

template<class MeshT>
bool
FixZipperNodeT<MeshT>::is_zipper_node(const VH _vh, const HEH _heh0, const HEH _heh1)
{
  //
  int val0 = valence_[mesh_.edge_handle(_heh0)];
  int val1 = valence_[mesh_.edge_handle(_heh1)];

  if (!(val0 == -1 && val1 == 1) && !(val0 == 1 && val1 == -1))
    return false;

  //check halfedge direction
  int rtax_a = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  _heh0, *mesh_.hehf_iter(_heh0));
  if (val0 == 1)
    rtax_a = negate(rtax_a);

  int rtax_b = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  _heh1, *mesh_.hehf_iter(_heh1));
  if (val1 == 1)
    rtax_b = negate(rtax_b);

  //start cell
  CH cha, chb;
  if (mesh_.is_boundary(_heh0))
    cha = mesh_.incident_cell(*mesh_.hehf_iter(_heh0));
  else
    cha = mesh_.incident_cell(mesh_.opposite_halfface_handle(*mesh_.hehf_iter(_heh0)));

  if (mesh_.is_boundary(_heh1))
    chb = mesh_.incident_cell(*mesh_.hehf_iter(_heh1));
  else
    chb = mesh_.incident_cell(mesh_.opposite_halfface_handle(*mesh_.hehf_iter(_heh1)));


  std::set<CH> onering_chs;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    onering_chs.insert(*vc_it);

  //express in chb
  int rt_ax_a_in_b = fac_.axis_in_chs_expressed_in_cht(cha, chb, rtax_a, onering_chs);


  //sg arc parallel to zipper node edges
  if (rt_ax_a_in_b == rtax_b)
    return true;

  return false;
}

template<class MeshT>
std::set<VH>
FixZipperNodeT<MeshT>::get_fixable_zipper_nodes(const bool _include_all)
{
  //update valence
  sge_.get_edges_valence();

  std::set<VH> cd_vhs;
  std::vector<VH> hv_nodes, hv_nodes2;
  std::vector<VH> rel_vhs;
  for (const auto &vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi])
    {
      MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);
      rel_vhs.push_back(vhi);
    }
  }

  for (const auto vhi: rel_vhs)
  {
    if (is_zipper_node(vhi))
      cd_vhs.insert(vhi);
    else
    {
      HEH sg_heh;
      if (is_fixable_constrained_zipper_node(vhi, sg_heh))
        cd_vhs.insert(vhi);
    }

    //fix one zipper node of the ones connected to high valence feature nodes at a time
    if (feature_node_[vhi])
    {
      if (n_incident_feature_edges(mesh_, feature_edge_, vhi) > 5)
      {
        hv_nodes.push_back(vhi);
      }
      else if (n_incident_special_edges_kept_on_feature_face(mesh_, feature_face_edge_, feature_edge_, valence_,
                                                             keep_on_feature_face_, vhi) >= 6)
      {
        hv_nodes2.push_back(vhi);
      }
    }
  }

//	std::cerr<<"angle threshold: "<<angle_thr_<<std::endl;
  auto flt_tps = filter_zipper_nodes(cd_vhs, _include_all);

  //TODO: the rt axes are not safe if has orthogonal sge
  for (const auto vhi: hv_nodes)
  {
    std::vector<HEH> sg_hehs;
    for (auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
    {
      if (std::abs(valence_[mesh_.edge_handle(*voh_it)]) == 1)
      {
        sg_hehs.push_back(*voh_it);
      }
    }

    auto raw_axes = fac_.rotation_axes_expressed_in_common_tet(vhi, sg_hehs);
    if (raw_axes.empty())
      continue;

    std::vector<int> v_e_axes;
    for (auto i = 0u; i < sg_hehs.size(); ++i)
    {
      if (valence_[mesh_.edge_handle(sg_hehs[i])] == 1)
      {
        v_e_axes.push_back(negate(raw_axes[i]));
      }
      else
        v_e_axes.push_back(raw_axes[i]);
    }

    std::vector<std::vector<HEH>> ax_hehs(6);
    for (auto j = 0u; j < v_e_axes.size(); ++j)
      for (int i = 0; i < 6; ++i)
      {
        if (v_e_axes[j] == i)
          ax_hehs[i].push_back(sg_hehs[j]);
      }

    //find the fixable tp on the longest sg arc and erase the rest
    for (int i = 0; i < 6; ++i)
    {
      VH best_tp(-1);
      double max_len = -1.;
      bool found_first = false;
      for (const auto hehi: ax_hehs[i])
      {
        auto sg_hehs_i = get_halfedges_on_singular_arc_with_intertior_surface(mesh_, valence_, feature_face_vertex_,
                                                                              hehi);
        double len = 0.;
        for (auto sg_hehi: sg_hehs_i)
          len += mesh_.length(sg_hehi);

        std::vector<std::pair<VH, int>> arc_tp_vi;
        get_zipper_nodes_on_singular_arc(sg_hehs_i, arc_tp_vi);

//                    std::cerr<<"tps on arc: ";
//                    for(const auto& [vhj, id] : arc_tp_vi)
//                        std::cerr<<" "<<vhj;
//                    std::cerr<<std::endl;

        VH fxb_tp(-1);
        for (const auto&[vhj, id]: arc_tp_vi)
        {
          if (flt_tps.find(vhj) != flt_tps.end())
          {
            fxb_tp = vhj;
            break;
          }
        }

        if (fxb_tp.is_valid())
        {
          if (len > max_len)
          {
            best_tp = fxb_tp;
            max_len = len;
          }

          //erase
          for (const auto&[vhj, id]: arc_tp_vi)
            flt_tps.erase(vhj);
        }
      }

      if (best_tp.is_valid())
        flt_tps.insert(best_tp);
      ALGOHEX_DEBUG_ONLY(std::cerr << "best tp " << best_tp << std::endl;)
    }
  }

  for (const auto vhi: hv_nodes2)
  {
    std::vector<HEH> sg_hehs;
    for (auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
    {
      EH ehi = mesh_.edge_handle(*voh_it);
      if (std::abs(valence_[ehi]) == 1 && !feature_face_edge_[ehi])
      {
        sg_hehs.push_back(*voh_it);
      }
    }

    if (sg_hehs.empty())
      continue;

    for (const auto hehi: sg_hehs)
    {
      //find the start guide
      EH eh_query = mesh_.edge_handle(hehi);
      HFH hf_query = *mesh_.hehf_iter(hehi);
      CH ch_query = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_query));
      int sg_rt_ax = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                        valence_, sg_hehs[0], hf_query);
      //search in singular halfedge direction
      int search_ax = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, hehi, ch_query, sg_rt_ax);


      //find a start feature halfedge
      HEH he_s(-1);
      HFH hf_s(-1), hf_e(-1);
      CH ch_s(-1);
      int search_ax_s = -1;
      int rt_ax_s = -1;

      bool has_in_dir = fac_.has_singular_edge_in_direction(hehi, ch_query, search_ax, vhi, ch_s, hf_s, he_s,
                                                            search_ax_s);

      if (!has_in_dir)
      {
        continue;
      }

      if (!keep_on_feature_face_[mesh_.edge_handle(he_s)])
        continue;

      //one ring cells
      std::set<CH> or_chs = get_onering_cells(mesh_, vhi);

      std::set<HFH> ft_hfhs;
      auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vhi, or_chs, ch_query, ft_hfhs);

      int n_nffsge = n_incident_non_ffe_singular_edges_in_same_region(mesh_, feature_face_edge_, valence_, vhi,
                                                                      por_chs);

      if (n_nffsge > 1)
        continue;
      auto sg_hehs_i = get_halfedges_on_singular_arc_with_intertior_surface(mesh_, valence_, feature_face_vertex_,
                                                                            hehi);

      std::vector<std::pair<VH, int>> arc_tp_vi;
      get_zipper_nodes_on_singular_arc(sg_hehs_i, arc_tp_vi);

      for (const auto&[vhj, id]: arc_tp_vi)
        flt_tps.erase(vhj);
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "Filter done" << std::endl;)

  return flt_tps;
}


template<class MeshT>
int
FixZipperNodeT<MeshT>::
get_zipper_nodes_on_singular_arc(const std::vector<HEH> &_sg_hehs, std::vector<std::pair<VH, int>> &_tp_vi)
{
  if (_sg_hehs.empty())
    return -1;

  std::vector<VH> sg_vhs;
  for (const auto hehi: _sg_hehs)
  {
    sg_vhs.push_back(mesh_.halfedge(hehi).from_vertex());
  }
  if (mesh_.halfedge(_sg_hehs.back()).to_vertex() != sg_vhs.front())
    sg_vhs.push_back(mesh_.halfedge(_sg_hehs.back()).to_vertex());

  int last_tp_pos = 0;

  HEH he_tmp(-1);
  if (is_zipper_node(sg_vhs[0]))
    _tp_vi.emplace_back(sg_vhs[0], 0);
  else if (is_fixable_constrained_zipper_node(sg_vhs[0], he_tmp))
  {
    _tp_vi.emplace_back(sg_vhs[0], 0);
  }

  for (auto i = 1u; i < sg_vhs.size() - 1; ++i)
    if (MeshPropertiesT<MeshT>::node_index(sg_vhs[i]) == 10)
    {
      _tp_vi.emplace_back(sg_vhs[i], i);
      last_tp_pos = i;
    }

  if (is_zipper_node(sg_vhs.back()))
  {
    last_tp_pos = sg_vhs.size() - 1;
    _tp_vi.emplace_back(sg_vhs.back(), last_tp_pos);
  }
  else if (is_fixable_constrained_zipper_node(sg_vhs.back(), he_tmp))
  {
    last_tp_pos = sg_vhs.size() - 1;
    _tp_vi.emplace_back(sg_vhs.back(), last_tp_pos);
  }

  return last_tp_pos;
}


template<class MeshT>
std::set<VH>
FixZipperNodeT<MeshT>::filter_zipper_nodes(std::set<VH> &_vhs, const bool _include_all)
{
  std::set<VH> fix_tps;
  while (!_vhs.empty())
  {
    auto vh_cur = *_vhs.begin();

    //find the start halfedge
    HEH heh(-1);
    if (is_zipper_node(vh_cur))
    {
      for (auto voh_it = mesh_.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
      {
        auto ve = mesh_.edge_handle(*voh_it);
        if (fabs(valence_[ve]) == 1 && !feature_face_edge_[ve])
        {
          heh = *voh_it;
          break;
        }
      }
    }
    else
    {
      HEH sg_heh;
      if (is_fixable_constrained_zipper_node(vh_cur, sg_heh))
        heh = sg_heh;
    }

    if (!heh.is_valid())
    {
      _vhs.erase(_vhs.begin());
      continue;
    }

    //process tps on the sg arc
    auto sg_hehs = get_halfedges_on_singular_arc_with_intertior_surface(mesh_, valence_, feature_face_vertex_, heh);

    std::vector<std::pair<VH, int>> tp_vi;
    int last_pos = get_zipper_nodes_on_singular_arc(sg_hehs, tp_vi);

    auto fix_tps_i = get_fixable_zipper_nodes_on_arc(sg_hehs, tp_vi);

    fix_tps.insert(fix_tps_i.begin(), fix_tps_i.end());

    //erase
    for (const auto&[vhi, id]: tp_vi)
      _vhs.erase(vhi);

    //if no zp node found by filter. Check if there are tps which cannot be removed by uniforming valence of the arc. e.g. elliptic sector
    if (_include_all)
    {
      if (fix_tps_i.empty() && !tp_vi.empty())
      {
        auto fix_tps_j = get_non_meshable_zipper_nodes_on_singular_arc(sg_hehs, tp_vi);
        fix_tps.insert(fix_tps_j.begin(), fix_tps_j.end());
      }
    }
  }

  return fix_tps;
}

template<class MeshT>
std::set<VH> FixZipperNodeT<MeshT>::get_fixable_zipper_nodes_on_arc(const std::vector<HEH> &_sg_hehs,
                                                                    const std::vector<std::pair<VH, int>> &_tp_vi) const
{
  std::set<VH> fix_tps;

  auto tp_vi = _tp_vi;
  VH vh0 = mesh_.halfedge(_sg_hehs[0]).from_vertex();
  VH vh1 = mesh_.halfedge(_sg_hehs.back()).to_vertex();


  //compute dominant cell axis in one ring
  std::map<CH, int> cell_axis;
  auto hehf_it = mesh_.hehf_iter(_sg_hehs[0]);
  CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
  int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   _sg_hehs[0], *hehf_it);
  int val = valence_[mesh_.edge_handle(_sg_hehs[0])];
  if (val == 1)
    rt_axis = negate(rt_axis);
  cell_axis[ch_s] = rt_axis;

  for (; hehf_it.valid(); ++hehf_it)
  {
    auto chi = mesh_.incident_cell(*hehf_it);
    if (chi == ch_s)
      break;

    CH chopp = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
    ALGOHEX_DEBUG_ONLY(if (cell_axis.find(chopp) == cell_axis.end())
      std::cerr << "error: couldn't find start cell " << chopp << std::endl;)

    cell_axis[chi] = tq_.axis_after_transition(cell_axis[chopp], trans_prop_[*hehf_it]);
  }

  for (auto i = 0u; i < _sg_hehs.size() - 1; ++i)
  {
    CH ch_prev = *mesh_.hec_iter(_sg_hehs[i]);
    int ax_prev = cell_axis[ch_prev];

    VH vh_t = mesh_.halfedge(_sg_hehs[i]).to_vertex();

    std::set<CH> v_or_cells;
    for (auto vc_it = mesh_.vc_iter(vh_t); vc_it.valid(); ++vc_it)
      v_or_cells.insert(*vc_it);

    auto hfhi_it = mesh_.hehf_iter(_sg_hehs[i + 1]);
    CH chopp_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hfhi_it));

    int axis_i = fac_.axis_in_chs_expressed_in_cht(ch_prev, chopp_s, ax_prev, v_or_cells);
    cell_axis[chopp_s] = axis_i;

    for (; hfhi_it.valid(); ++hfhi_it)
    {
      auto chi = mesh_.incident_cell(*hfhi_it);
      if (chi == chopp_s)
        break;

      CH chopp = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hfhi_it));
      ALGOHEX_DEBUG_ONLY(if (cell_axis.find(chopp) == cell_axis.end())
        std::cerr << "error: couldn't find cell " << chopp << std::endl;)

      cell_axis[chi] = tq_.axis_after_transition(cell_axis[chopp], trans_prop_[*hfhi_it]);
    }
  }


  //compute coordinates in frame space
  std::map<VH, double> v_coords;

  v_coords[vh0] = 0.;
  for (auto i = 0u; i < _sg_hehs.size(); ++i)
  {
    VH vh_f = mesh_.halfedge(_sg_hehs[i]).from_vertex();
    VH vh_t = mesh_.halfedge(_sg_hehs[i]).to_vertex();

    CH ch_inc = *mesh_.hec_iter(_sg_hehs[i]);
    auto dm_axis = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_inc],
                                                           (AxisAlignment) cell_axis[ch_inc]);
    auto he_dir = ovm2eigen(mesh_.vertex(vh_t) - mesh_.vertex(vh_f));
    double disti = dm_axis.dot(he_dir);

    v_coords[vh_t] = v_coords[vh_f] + disti;
  }

  //in case of circle, v_coords[vh0] may be overwritten
  double vc0 = 0.;
  double vc1 = v_coords[vh1];


  //average segment length
  double avrg_seg_lgt = 0.;
  avrg_seg_lgt += std::fabs(vc0 - v_coords[tp_vi[0].first]);
  for (auto i = 0u; i < tp_vi.size() - 1; ++i)
  {
    avrg_seg_lgt += std::fabs(v_coords[tp_vi[i + 1].first] - v_coords[tp_vi[i].first]);
  }
  avrg_seg_lgt += std::fabs(vc1 - v_coords[tp_vi.back().first]);
  avrg_seg_lgt /= (double) (tp_vi.size() + 1);
  ALGOHEX_DEBUG_ONLY(std::cerr << "avr seg length: " << avrg_seg_lgt << std::endl;)

  double peak_threshold = peak_coeff_ * avrg_seg_lgt;

  ALGOHEX_DEBUG_ONLY(std::cerr << " raw tp vhs on arc: ";
                             for (auto[vhi, id]: tp_vi)
                               std::cerr << " " << vhi;
                             std::cerr << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << " vh0: " << vh0 << " coord: " << vc0;
                             for (auto[vhi, id]: tp_vi)
                               std::cerr << " vh: " << vhi << " coord: " << v_coords[vhi];
                             std::cerr << " vh1: " << vh1 << " coord: " << vc1;
                             std::cerr << std::endl;)

  bool is_v0_btp = tp_vi[0].first == vh0 && feature_face_vertex_[vh0];
  bool is_v1_btp = tp_vi.back().first == vh1 && feature_face_vertex_[vh1];

  //erase noise if at feature. make sure the peak tp is correct
  //vh0 is feature zipper node
  if (is_v0_btp && tp_vi.size() > 1)
  {
    if (!feature_face_vertex_[tp_vi[1].first])
    {
      //erase the first two if they are too close
      if (std::fabs(v_coords[tp_vi[0].first] - v_coords[tp_vi[1].first]) <= peak_threshold)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "vh0 is tp on ff";)
        tp_vi.erase(tp_vi.begin());
        tp_vi.erase(tp_vi.begin());
        ALGOHEX_DEBUG_ONLY(std::cerr << "erased first two vhs" << std::endl;)
      }
    }
  }

  if (is_v1_btp && tp_vi.size() > 1)
  {
    if (!feature_face_vertex_[tp_vi[tp_vi.size() - 2].first])
    {
      //erase the last two if they are too close
      if (std::fabs(v_coords[tp_vi.back().first] - v_coords[tp_vi[tp_vi.size() - 2].first]) <= peak_threshold)
      {
        tp_vi.pop_back();
        tp_vi.pop_back();
        ALGOHEX_DEBUG_ONLY(std::cerr << "erased last two vhs" << std::endl;)
      }
    }
  }

  if (!tp_vi.empty())
  {
    //take the peak if possible
    std::vector<DV> dvs;
    for (const auto &[vhi, id]: tp_vi)
      dvs.emplace_back(v_coords[vhi], vhi);
    std::sort(dvs.begin(), dvs.end());


    //min len of max positive pt to end pts
    double peak_pos = std::min(v_coords[dvs.back().second] - vc0,
                               v_coords[dvs.back().second] - vc1);
    double peak_neg = std::min(vc0 - v_coords[dvs[0].second],
                               vc1 - v_coords[dvs[0].second]);
    ALGOHEX_DEBUG_ONLY(std::cerr << "peak_pos: " << peak_pos << " peak_neg: " << peak_neg << " peak_threshold "
                                 << peak_threshold;)
    double max_peak_dist = std::max(peak_neg, peak_pos);
    VH vh_peak(-1);
    if (max_peak_dist == peak_pos)
      vh_peak = dvs.back().second;
    else if (max_peak_dist == peak_neg)
      vh_peak = dvs[0].second;

    if (max_peak_dist >= 0.)
    {
      if (tp_vi.size() % 2 == 1)
        fix_tps.insert(vh_peak);
      else if (tp_vi.size() > 1)
      {//tp is the max on the arc
        if (peak_pos >= peak_threshold)
          fix_tps.insert(vh_peak);
      }
    }

    //if incident segments at zipper node are large enough
    double seg_threshold = seg_coeff_ * avrg_seg_lgt;
    double av_len = 0.;
    for (const auto&[vhi, id]: tp_vi)
      av_len += target_length_[vhi];
    av_len /= (double) tp_vi.size();
    double absl_seg_threshold = absl_seg_coeff_ * av_len;
    ALGOHEX_DEBUG_ONLY(std::cerr << "\naverage edge length: " << av_len << " seg_threshold " << seg_threshold
                                 << " thrhold wrt target length " << absl_seg_threshold;
                               std::cerr << "\ntps size: " << tp_vi.size() << " sgvhs size " << _sg_hehs.size() + 1;
                               std::cerr << " remaining tps: ";)
    if (tp_vi.size() == 1)
    {
      double dist0 = std::fabs(v_coords[tp_vi[0].first] - vc0);
      double dist1 = std::fabs(vc1 - v_coords[tp_vi[0].first]);
      ALGOHEX_DEBUG_ONLY(std::cerr << "Tp size=1. seg " << seg_threshold << " dist0: " << dist0 << " dist1: " << dist1
                                   << std::endl;)

      if (feature_face_vertex_[vh0])//if one segment ends at the boundary and the other is large enough, take it
      {
        if (feature_face_vertex_[vh1])
        {
          fix_tps.insert(tp_vi[0].first);
          ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi[0].first << " ";)
        }
        else if (dist1 > seg_threshold || dist1 > absl_seg_threshold)
        {
          fix_tps.insert(tp_vi[0].first);
          ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi[0].first << " ";)
        }
      }
      else
      {
        if (feature_face_vertex_[vh1])
        {
          if (dist0 > seg_threshold || dist0 > absl_seg_threshold)
          {
            fix_tps.insert(tp_vi[0].first);
            ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi[0].first << " ";)
          }
        }
        else
        {
          if ((dist0 > seg_threshold && dist1 > seg_threshold) ||
              (dist0 > absl_seg_threshold && dist1 > absl_seg_threshold))
          {
            fix_tps.insert(tp_vi[0].first);
            ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi[0].first << " ";)
          }
        }

      }
    }
    else if (tp_vi.size() > 1)
    {
      if (tp_vi.size() * 5 < (_sg_hehs.size() + 1) || tp_vi.size() <= 2)
      {//in case too many tps on arc
        //beginning
        double dist0 = std::fabs(v_coords[tp_vi[0].first] - vc0);
        double dist1 = std::fabs(v_coords[tp_vi[1].first] - v_coords[tp_vi[0].first]);
        if (feature_face_vertex_[vh0])//if one segment ends at the boundary and the other is large enough, take it
        {
          if (dist1 > seg_threshold || dist1 > absl_seg_threshold)
          {
            ALGOHEX_DEBUG_ONLY(std::cerr << "Tp start to bdy: seg_threshold: " << seg_threshold << " dist1: " << dist1
                                         << " considered: " << dist1 << std::endl;)
            fix_tps.insert(tp_vi[0].first);
            ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi[0].first << " ";)
          }
        }
        else
        {
          if ((dist0 > seg_threshold && dist1 > seg_threshold) ||
              (dist0 > absl_seg_threshold && dist1 > absl_seg_threshold))
          {
            fix_tps.insert(tp_vi[0].first);
            ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi[0].first << " ";)
          }
        }

        //middle
        for (auto i = 1u; i < tp_vi.size() - 1; ++i)
        {
          dist0 = std::fabs(v_coords[tp_vi[i].first] - v_coords[tp_vi[i - 1].first]);
          dist1 = std::fabs(v_coords[tp_vi[i + 1].first] - v_coords[tp_vi[i].first]);
          if ((dist0 > seg_threshold && dist1 > seg_threshold) ||
              (dist0 > absl_seg_threshold && dist1 > absl_seg_threshold))
          {
            fix_tps.insert(tp_vi[i].first);
            ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi[i].first << " ";)

          }
        }

        //end
        dist0 = std::fabs(v_coords[tp_vi.back().first] - vc1);
        dist1 = std::fabs(v_coords[tp_vi[tp_vi.size() - 2].first] - v_coords[tp_vi.back().first]);
        if (feature_face_vertex_[vh1])
        {//if one segment ends at the boundary and the other is large enough, take it
          if (dist1 > seg_threshold || dist1 > absl_seg_threshold)
          {
            ALGOHEX_DEBUG_ONLY(std::cerr << "Tp end to bdy: seg_threshold: " << seg_threshold << " dist0: " << dist0
                                         << " dist considered: " << dist1 << std::endl;)
            fix_tps.insert(tp_vi.back().first);

            ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi.back().first << " ";)
          }
        }
        else
        {
          if ((dist0 > seg_threshold && dist1 > seg_threshold) ||
              (dist0 > absl_seg_threshold && dist1 > absl_seg_threshold))
          {
            fix_tps.insert(tp_vi.back().first);
            ALGOHEX_DEBUG_ONLY(std::cerr << tp_vi.back().first << " ";)
          }
        }
      }
    }
//                std::cerr<<std::endl;

    ALGOHEX_DEBUG_ONLY(std::cerr << "filtered tps on arc ";
                               for (auto vhi: fix_tps)
                                 std::cerr << " " << vhi;
                               std::cerr << std::endl;)
  }

  return fix_tps;
}

template<class MeshT>
std::vector<VH> FixZipperNodeT<MeshT>::get_non_meshable_zipper_nodes_on_singular_arc(const std::vector<HEH> &_sg_hehs,
                                                                                     const std::vector<std::pair<VH, int>> &_tp_vi)
{
  std::vector<VH> sg_vhs;
  for (auto j = 0u; j < _sg_hehs.size(); ++j)
  {
    sg_vhs.push_back(mesh_.halfedge(_sg_hehs[j]).from_vertex());
  }
  sg_vhs.push_back(mesh_.halfedge(_sg_hehs.back()).to_vertex());


  //split
  for (const auto vhi: sg_vhs)
    SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, vhi, valence_, feature_face_edge_, feature_face_vertex_,
                                            feature_edge_,
                                            sgl_vt_, feature_edge_vertex_);

  std::map<CH, Quaternion> old_qtns;

  //uniform valence
  std::map<EH, int> old_eh_val;
  int n_val_n = 0, n_val_p = 0;
  int abs_val = std::abs(valence_[mesh_.edge_handle(_sg_hehs[0])]);
  double len_n = 0., len_p = 0.;
  for (const auto hehi: _sg_hehs)
  {
    auto ehi = mesh_.edge_handle(hehi);
    old_eh_val[ehi] = valence_[ehi];

    if (valence_[ehi] == -abs_val)
      len_n += mesh_.length(ehi);
    else if (valence_[ehi] == abs_val)
      len_p += mesh_.length(ehi);
  }

  int val = len_n > len_p ? -abs_val : abs_val;

  for (const auto hehi: _sg_hehs)
  {
    auto ehi = mesh_.edge_handle(hehi);
    valence_[ehi] = val;
  }

  //smooth quaternions
  smooth_relevant_quaternions(sg_vhs, old_qtns);

  //check local meshability
  std::vector<std::pair<VH, int>> all_nm_vi;
  for (auto j = 1u; j < sg_vhs.size() - 1; ++j)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "check vh: " << sg_vhs[j] << std::endl;)
    bool suc = lmcwf_.is_locally_meshable(sg_vhs[j], true);
    if (!suc)
      all_nm_vi.emplace_back(sg_vhs[j], j);
  }

  std::set<VH> tp_set;
  for (const auto&[vhi, id]: _tp_vi)
    tp_set.insert(vhi);

  std::vector<VH> nm_vhs;
  for (const auto&[vhi, id]: all_nm_vi)
    if (tp_set.find(vhi) != tp_set.end())
      nm_vhs.push_back(vhi);

  //return the closest tp
  if (!all_nm_vi.empty() && nm_vhs.empty())
  {
    VH best_tp(-1);
    int best_id = -1;
    int min_dist = std::numeric_limits<int>::max();
    for (const auto&[vhi, id]: _tp_vi)
    {
      int dist = fabs(id - all_nm_vi[0].second);
      if (dist < min_dist)
      {
        min_dist = dist;
        best_tp = vhi;
        best_id = id;
      }
    }

    if (best_tp.is_valid())
    {
      bool found_other_tp = false;
      if (best_id > all_nm_vi[0].second)
      {
        for (const auto&[vhi, id]: _tp_vi)
          if (id < all_nm_vi[0].second)
          {
            found_other_tp = true;
            break;
          }
      }
      else if (best_id < all_nm_vi[0].second)
      {
        for (const auto&[vhi, id]: _tp_vi)
          if (id > all_nm_vi[0].second)
          {
            found_other_tp = true;
            break;
          }
      }

      if (found_other_tp)
      {
        nm_vhs.push_back(best_tp);
        ALGOHEX_DEBUG_ONLY(std::cerr << "found tps bounding non meshable vt" << std::endl;)
      }
      else
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "not found tps bounding non meshable vt" << std::endl;)
      }
    }
  }


  ALGOHEX_DEBUG_ONLY(std::cerr << "Nm tps: ";
  for (auto vhi: nm_vhs)
    std::cerr << " " << vhi;
  std::cerr << std::endl;)

  //recover valence
  for (const auto hehi: _sg_hehs)
  {
    auto ehi = mesh_.edge_handle(hehi);
    valence_[ehi] = old_eh_val[ehi];
  }

  //recover quaternions
  for (const auto&[chi, qtn]: old_qtns)
    cell_quaternions_[chi] = qtn;

  return nm_vhs;
}

template<class MeshT>
bool FixZipperNodeT<MeshT>::is_singular_arc_meshable(const std::vector<HEH> &_sg_hehs,
                                                     const std::vector<std::pair<VH, int>> &_tp_vi)
{
  if (_tp_vi.size() % 2 == 1)
    return false;

  if (_tp_vi.empty())
    return true;

  if (_tp_vi.size() == 2)
  {
    if (is_singular_circle(_sg_hehs.front()))
      return false;
  }

  std::set<VH> sg_vhs;
  for (auto j = 0u; j < _sg_hehs.size(); ++j)
  {
    sg_vhs.insert(mesh_.halfedge(_sg_hehs[j]).from_vertex());
  }
  sg_vhs.insert(mesh_.halfedge(_sg_hehs.back()).to_vertex());


  //split
  for (const auto vhi: sg_vhs)
    SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, vhi, valence_, feature_face_edge_, feature_face_vertex_,
                                            feature_edge_,
                                            sgl_vt_, feature_edge_vertex_);

  std::map<CH, Quaternion> old_qtns;
  //uniform valence
  std::map<EH, int> old_eh_val;
  int n_val_n = 0, n_val_p = 0;
  int abs_val = std::abs(valence_[mesh_.edge_handle(_sg_hehs[0])]);
  double len_n = 0., len_p = 0.;
  for (const auto hehi: _sg_hehs)
  {
    auto ehi = mesh_.edge_handle(hehi);
    old_eh_val[ehi] = valence_[ehi];

    if (valence_[ehi] == -abs_val)
      len_n += mesh_.length(ehi);
    else if (valence_[ehi] == abs_val)
      len_p += mesh_.length(ehi);
  }

  int val = len_n > len_p ? -abs_val : abs_val;

  for (const auto hehi: _sg_hehs)
  {
    auto ehi = mesh_.edge_handle(hehi);
    valence_[ehi] = val;
  }

  //smooth quaternions
  std::vector<VH> vec_sg_vhs(sg_vhs.begin(), sg_vhs.end());
  smooth_relevant_quaternions(vec_sg_vhs, old_qtns);
//        QuaternionSmoothing::optimize_quaternions_wrt_prescribed_valence(mesh_, trans_prop_, tq_, cell_quaternions_, valence_, feature_edge_, feature_fprop_, sg_cells, true);


  //check local meshability
  for (auto vhi: sg_vhs)
  {
    std::cerr << "check vh: " << vhi << std::endl;
    bool suc = lmcwf_.is_locally_meshable(vhi, true, 1);
    if (!suc)
    {
      //recover valence
      for (const auto hehi: _sg_hehs)
      {
        auto ehi = mesh_.edge_handle(hehi);
        valence_[ehi] = old_eh_val[ehi];
      }

      //recover quaternions
      for (const auto&[chi, qtn]: old_qtns)
        cell_quaternions_[chi] = qtn;

      return false;
    }
  }


  //recover valence
  for (const auto hehi: _sg_hehs)
  {
    auto ehi = mesh_.edge_handle(hehi);
    valence_[ehi] = old_eh_val[ehi];
  }

  //recover quaternions
  for (const auto&[chi, qtn]: old_qtns)
    cell_quaternions_[chi] = qtn;

  return true;
}

template<class MeshT>
void
FixZipperNodeT<MeshT>::smooth_relevant_quaternions(const std::vector<VH> &_sg_vhs, std::map<CH, Quaternion> &_old_qtns)
{

  std::set<CH> sg_cells;
  for (const auto vhi: _sg_vhs)
    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
      sg_cells.insert(*vc_it);

  for (const auto chi: sg_cells)
    _old_qtns.insert(std::make_pair(chi, cell_quaternions_[chi]));
  std::queue<CH> que;
  for (const auto chi: sg_cells)
    que.push(chi);

  int i = 0;
  while (!que.empty() && i < 1000000)
  {
    auto ch_cur = que.front();
    que.pop();
//            std::cerr<<" "<<ch_cur;
    Quaternion qt_old = cell_quaternions_[ch_cur];

    //store qtn for recovery
    if (_old_qtns.find(ch_cur) == _old_qtns.end())
      _old_qtns[ch_cur] = qt_old;

    QuaternionSmoothing::optimize_quaternions_wrt_prescribed_valence(mesh_, trans_prop_, tq_, cell_quaternions_,
                                                                     valence_, feature_edge_, feature_fprop_,
                                                                     std::set<CH>{ch_cur}, true);

    Quaternion rt = cell_quaternions_[ch_cur] * qt_old.conjugate();
    Eigen::AngleAxis<double> aa(rt);

    double angle = aa.angle();

    if (angle > angle_thr_)
    {
      for (auto cc_it = mesh_.cc_iter(ch_cur); cc_it.valid(); ++cc_it)
        que.push(*cc_it);
    }

    i++;
  }
}

template<class MeshT>
typename AlgoHex::CELLINFO
FixZipperNodeT<MeshT>::next_cell_info(const CELLINFO &ci_cur, const int _u, const int _v, const int _w,
                                      const HFH next_hfh, const double _vw_weight, const double _length_scale,
                                      bool &_increase) const
{
  //the first segment
  auto c_bct = mesh_.barycenter(ci_cur.ch);
  auto f_bct = mesh_.barycenter(mesh_.face_handle(next_hfh));
  Vec3d dir = ovm2eigen(f_bct - c_bct);

  //u expressed in the current cell
  int u_cur = tq_.axis_after_transition(_u, ci_cur.trans);
  Vec3d u_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ci_cur.ch], (AxisAlignment) u_cur);
  //v expressed in the current cell
  int v_cur = tq_.axis_after_transition(_v, ci_cur.trans);
  Vec3d v_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ci_cur.ch], (AxisAlignment) v_cur);
  //w expressed in the current cell
  int w_cur = tq_.axis_after_transition(_w, ci_cur.trans);
  Vec3d w_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ci_cur.ch], (AxisAlignment) w_cur);

  //first segment
  Vec3d param_next = ci_cur.param;
  param_next[_u / 2] += u_dir.dot(dir);
  param_next[_v / 2] += v_dir.dot(dir);
  param_next[_w / 2] += w_dir.dot(dir);


  //check if it's increasing in u dir
  Vec3d face_dir = ovm2eigen(mesh_.normal(next_hfh));
  double inc_val = face_dir.dot(u_dir);
  if (inc_val < inc_val_thr_)
    _increase = false;
  else
    _increase = true;

  //the second segment
  auto ch_next = mesh_.incident_cell(next_hfh);

  auto c_next_bct = mesh_.barycenter(ch_next);
  dir = ovm2eigen(c_next_bct - f_bct);

  int trans_next = tq_.mult_transitions_idx(trans_prop_[next_hfh], ci_cur.trans);
  //u expressed in next cell
  int u_n = tq_.axis_after_transition(_u, trans_next);
  Vec3d un_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_next], (AxisAlignment) u_n);
  //v expressed in next cell
  int v_n = tq_.axis_after_transition(_v, trans_next);
  Vec3d vn_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_next], (AxisAlignment) v_n);
  //w expressed in next cell
  int w_n = tq_.axis_after_transition(_w, trans_next);
  Vec3d wn_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_next], (AxisAlignment) w_n);


  param_next[_u / 2] += un_dir.dot(dir);
  param_next[_v / 2] += vn_dir.dot(dir);
  param_next[_w / 2] += wn_dir.dot(dir);

  //length of the segment
  double scale = _length_scale * _length_scale;
  double dist = sqrt(scale * param_next[_u / 2] * param_next[_u / 2]
                     + scale * _vw_weight *
                       (param_next[_v / 2] * param_next[_v / 2] + param_next[_w / 2] * param_next[_w / 2]));

  CELLINFO next_ci(ch_next, param_next);
  next_ci.trans = trans_next;
  next_ci.length = dist;


//        if(_increase)
//        std::cerr<<"    inc val "<<inc_val<<" "<<next_ci;

  return next_ci;
}

template<class MeshT>
CH
FixZipperNodeT<MeshT>::find_start_cell(const VH _vh, HEH &_heh, HFH &_hfh) const
{
  std::vector<HEH> sghehs;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    if (valence_[mesh_.edge_handle(*voh_it)] == 1 || valence_[mesh_.edge_handle(*voh_it)] == -1)
    {
      if (!feature_face_edge_[mesh_.edge_handle(*voh_it)])
      {
        sghehs.push_back(mesh_.opposite_halfedge_handle(*voh_it));
      }
    }


  //take the cell whose barycenter is the furthest
  double max_dist = std::numeric_limits<double>::lowest();
  CH ch_max(-1);
  for (const auto hehi: sghehs)
  {
    HFH hfhi;
    auto dch = find_start_cell_at_halfedge(_vh, hehi, hfhi);
    if (dch.first > max_dist && dch.second.is_valid())
    {
      max_dist = dch.first;
      ch_max = dch.second;
      _hfh = hfhi;
      _heh = hehi;
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << " ch max " << ch_max << " dist " << max_dist << std::endl;)

  return ch_max;
}

template<class MeshT>
std::pair<double, CH>
FixZipperNodeT<MeshT>::find_start_cell_at_halfedge(const VH _vh, const HEH _heh, HFH &_hfh) const
{
  //take the cell whose barycenter is the furthest
  double max_dist = std::numeric_limits<double>::lowest();
  CH ch_max(-1);
  //get halfedge direction
  auto hf_it = mesh_.hehf_iter(_heh);
  int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   _heh, *hf_it);
  int e_dir_ax = valence_[mesh_.edge_handle(_heh)] == 1 ? negate(rt_axis) : rt_axis;

  CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hf_it));
  Vec3d e_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_s], (AxisAlignment) e_dir_ax);

  //barycenter in frame space
  Point pt0 = mesh_.vertex(_vh);
  Point pt = mesh_.barycenter(ch_s);
  Vec3d vdir = ovm2eigen(pt - pt0);

  max_dist = vdir.dot(e_dir);
  _hfh = *hf_it;
  ch_max = ch_s;

  for (; hf_it.valid(); ++hf_it)
  {
    CH ch_it = mesh_.incident_cell(*hf_it);
    if (ch_it == ch_s)
      break;

    //halfedge direction in current cell
    e_dir_ax = tq_.axis_after_transition(e_dir_ax, trans_prop_[*hf_it]);
    e_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_it], (AxisAlignment) e_dir_ax);

    pt = mesh_.barycenter(ch_it);
    vdir = ovm2eigen(pt - pt0);

    double dist = vdir.dot(e_dir);

    if (dist > max_dist)
    {
      ch_max = ch_it;
      max_dist = dist;
      _hfh = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(*hf_it, _heh));
    }
  }

  return {max_dist, ch_max};
}

template<class MeshT>
std::vector<CH>
FixZipperNodeT<MeshT>::shortest_dual_path_to_boundary(const HEH _heh_s, const HFH _hfh_s,
                                                      const std::set<VH> &_fixable_vhs, VH &_other_tp, const double _w)
{
  //transition axis of sg edge
  int trans_idx = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, _heh_s, _hfh_s);
  if (trans_idx == -1)
    return std::vector<CH>{};

  //rotational axis direction
  int axis = tq_.rotation_axis_2d(trans_idx);
  //edge direction in parametrization coordinates
  int valence = valence_[mesh_.edge_handle(_heh_s)];
  if (valence == 1)
    axis = negate(axis);

  CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(_hfh_s));

  return shortest_dual_path_to_boundary(ch_s, axis, _fixable_vhs, _other_tp, _w);
}

template<class MeshT>
std::vector<CH>
FixZipperNodeT<MeshT>::shortest_dual_path_to_boundary(const CH _ch_s, int _axis, const std::set<VH> &_fixable_vhs,
                                                      VH &_other_tp, const double _w)
{
  //edge direction in parametrization coordinates
  int u = _axis;
  int v = (u + 2) % 6;
  int w = the_third_axis((AxisAlignment) u, (AxisAlignment) v);

  CELLINFO chsinfo(_ch_s, Vec3d(0, 0, 0));


  ALGOHEX_DEBUG_ONLY(std::cerr << "init " << chsinfo << " u: " << u;
                             std::cerr << " u dir: " << mesh_.barycenter(_ch_s) << " "
                                       << (ovm2eigen(mesh_.barycenter(_ch_s)) +
                                           AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[_ch_s],
                                                                                   (AxisAlignment) u)).transpose()
                                       << std::endl;)

  std::vector<double> min_dist(mesh_.n_cells(), std::numeric_limits<double>::infinity());
  std::vector<CH> prv_ch(mesh_.n_cells(), CH(-1));

//        std::map<CH, CELLINFO> ch_chinfo;
  min_dist[_ch_s.idx()] = 0.;
//        ch_chinfo.insert(std::make_pair(_ch_s, chsinfo));

  ALGOHEX_DEBUG_ONLY(std::cerr << "Floodfill: " << _w << std::endl;)

  CH ch_target(-1);
  std::set<DCINFO> que;
  que.insert(DCINFO(0, chsinfo));
  double search_dist = std::numeric_limits<double>::infinity();
  while (!que.empty())
  {
    auto dc_cur = *que.begin();
    que.erase(que.begin());

//            std::cout<<" "<<dc_cur.second.ch;
    if (is_cell_found(dc_cur.second.ch, dc_cur.second.trans, u))
    {
      ch_target = dc_cur.second.ch;
      if (ch_target != _ch_s)
      {
        ALGOHEX_DEBUG_ONLY(
                std::cerr << " u " << u << " nm ax " << tq_.axis_after_transition(u, dc_cur.second.trans)
                          << " ";
                std::cerr << "found first boundary cell " << dc_cur.second.ch << " dist "
                          << dc_cur.second.length << std::endl;)

        //set the maximum search distance
        search_dist = max_tp_dist_ * dc_cur.second.length;
        break;
      }
    }

    if (connect_other_zipper_node_)
    {
      if (is_cell_found_wrt_zipper_node(dc_cur.second, _fixable_vhs, _other_tp, u))
      {
        ch_target = dc_cur.second.ch;
        if (ch_target != _ch_s)
        {
          ALGOHEX_DEBUG_ONLY(
                  std::cerr << "found other tp " << _other_tp << " in cell " << ch_target << std::endl;)

          break;
        }
      }
    }

    auto chfs = mesh_.cell(dc_cur.second.ch).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);

      auto ch_next = mesh_.incident_cell(hf_opp);
//                std::cout<<"     hf_opp "<<hf_opp<<" ch_next: "<<ch_next<<std::endl;

      if (!ch_next.is_valid())
        continue;

      FH fh_next = mesh_.face_handle(hf_opp);
      //feature face normal should agree with the search direction
      if (feature_fprop_[fh_next] > 0)
      {
        auto normal = mesh_.normal(hf_opp);
        auto ca = AxisAlignmentHelpers::closest_axis(cell_quaternions_[dc_cur.second.ch], normal);

        if (ca.second != tq_.axis_after_transition(u, dc_cur.second.trans))
          continue;
      }

      if (visited_fprop_[fh_next])
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "blocked " << fh_next;)
        continue;
      }

      //floodfill
      if (min_dist[ch_next.idx()] != std::numeric_limits<double>::infinity() || prv_ch[ch_next.idx()] != CH(-1))
        continue;

      bool increase = false;
      auto ci_next = next_cell_info(dc_cur.second, u, v, w, hf_opp, _w, 1., increase);

      if (!increase)
        continue;

      //absolute distance to origin
      min_dist[ch_next.idx()] = ci_next.length;

      prv_ch[ch_next.idx()] = dc_cur.second.ch;

      que.insert(std::make_pair(min_dist[ch_next.idx()], ci_next));
    }
  }

  if (connect_other_zipper_node_ && !_other_tp.is_valid())
  {
    CH ch_target2(-1);
    std::vector<double> min_dist2(mesh_.n_cells(), std::numeric_limits<double>::infinity());
    std::vector<CH> prv_ch2(mesh_.n_cells(), CH(-1));
    min_dist2[_ch_s.idx()] = 0.;
    prv_ch2[_ch_s.idx()] = CH(-1);
    que.insert(DCINFO(0, chsinfo));
    while (!que.empty())
    {
      auto dc_cur = *que.begin();

      que.erase(que.begin());

      if (dc_cur.second.length > search_dist)
        continue;

//                std::cout << "search dist b " << search_dist << " " << dc_cur.second;

      if (is_cell_found_wrt_zipper_node(dc_cur.second, _fixable_vhs, _other_tp, u))
      {
        if (ch_target2 != _ch_s)
        {
          ch_target2 = dc_cur.second.ch;
          ALGOHEX_DEBUG_ONLY(
                  std::cerr << "found other tp" << _other_tp << " in cell " << ch_target2 << " dist "
                            << dc_cur.second.length << std::endl;)

          break;
        }
      }

      auto chfs = mesh_.cell(dc_cur.second.ch).halffaces();
      for (int i = 0; i < 4; ++i)
      {
        auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);

        auto ch_next = mesh_.incident_cell(hf_opp);

        if (!ch_next.is_valid())
          continue;

        FH fh_next = mesh_.face_handle(hf_opp);

        //feature face normal should agree with the search direction
        if (feature_fprop_[fh_next] > 0)
        {
          auto normal = mesh_.normal(hf_opp);
          auto ca = AxisAlignmentHelpers::closest_axis(cell_quaternions_[dc_cur.second.ch],
                                                       normal);

          if (ca.second != tq_.axis_after_transition(u, dc_cur.second.trans))
            continue;
        }

        if (visited_fprop_[fh_next])
          continue;

        //floodfill
        if (min_dist2[ch_next.idx()] != std::numeric_limits<double>::infinity() || prv_ch2[ch_next.idx()] != CH(-1))
          continue;

        bool increase = false;
        auto ci_next = next_cell_info(dc_cur.second, u, v, w, hf_opp, _w, 0.5, increase);


        if (!increase)
          continue;

        //absolute distance to origin
        min_dist2[ch_next.idx()] = ci_next.length;

        prv_ch2[ch_next.idx()] = dc_cur.second.ch;

        que.insert(std::make_pair(min_dist2[ch_next.idx()], ci_next));
      }
    }

    //find cells
    //get the path
    if (_other_tp.is_valid())
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "Other tp " << _other_tp;)
      std::vector<CH> rpath, path;
      CH ch_it = ch_target2;
      while (ch_it.idx() >= 0 && rpath.size() < 1000000)
      {
        rpath.push_back(ch_it);
        ALGOHEX_DEBUG_ONLY(std::cerr << " " << ch_it;)

        ch_it = prv_ch2[ch_it.idx()];
      }

      if (!rpath.empty() && rpath.back() == _ch_s)
      {
        path.reserve(rpath.size());
        for (auto r_it = rpath.rbegin(); r_it != rpath.rend(); ++r_it)
          path.push_back(*r_it);
      }

      ALGOHEX_DEBUG_ONLY(std::cout << " Path0: ";
                                 for (const auto &ch: path)
                                   std::cout << " " << ch;
                                 std::cerr << std::endl;)

      return path;
    }
  }


  //find cells
  //get the path
  std::vector<CH> rpath, path;
  CH ch_it = ch_target;

  while (ch_it.idx() >= 0 && rpath.size() < 1000000)
  {
    rpath.push_back(ch_it);
    ch_it = prv_ch[ch_it.idx()];
  }

  if (!rpath.empty() && rpath.back() == _ch_s)
  {
    path.reserve(rpath.size());
    for (auto r_it = rpath.rbegin(); r_it != rpath.rend(); ++r_it)
      path.push_back(*r_it);
  }


  ALGOHEX_DEBUG_ONLY(std::cerr << " Path0: ";
                             for (const auto &ch: path)
                               std::cerr << " " << ch;
                             std::cerr << std::endl;)

  return path;
}

template<class MeshT>
std::vector<CH>
FixZipperNodeT<MeshT>::
shortest_dual_path_from_cell_to_target(const HEH _heh_s, const CH _ch_s, const std::set<CH> &_target_chs,
                                       const std::set<CH> &_dp_cells, const std::set<VH> &_excl_vh,
                                       const bool _spv_free) const
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "finding dual path from: " << _ch_s << " to targets: ";
                             for (auto chi: _target_chs)
                               std::cerr << " " << chi;
                             std::cerr << std::endl;)

  HFH hf_s(-1);
  for (auto hehf_it = mesh_.hehf_iter(_heh_s); hehf_it.valid(); ++hehf_it)
  {
    if (mesh_.incident_cell(*hehf_it) == _ch_s)
    {
      hf_s = *hehf_it;
      break;
    }
  }

  if (!hf_s.is_valid())
    std::cerr << "Error: start cell " << _ch_s << "is not incident to halfedge " << _heh_s << std::endl;

  //search axis in _ch_s

  //transition axis of sg edge
  HFH hf_s0 = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_s, _heh_s));

  int trans_idx = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, _heh_s, hf_s0);

  //rotational axis direction
  int e_ax = tq_.rotation_axis_2d(trans_idx);
  //edge direction in parametrization coordinates
  int valence = valence_[mesh_.edge_handle(_heh_s)];
  if (valence == 1)
    e_ax = negate(e_ax);

  //
  bool to_tp = _excl_vh.size() == 3;
  ALGOHEX_DEBUG_ONLY(std::cerr << "to other tp? " << to_tp << std::endl;)

  std::map<CH, int> min_dist;
  min_dist.insert(std::make_pair(_ch_s, 0));

  std::map<CH, CH> previous;
  previous.insert(std::make_pair(_ch_s, CH(-1)));

//std::vector<CH> tp_chs;
  auto vh_t = mesh_.halfedge(_heh_s).to_vertex();


  std::set<DCHI> que;
  que.insert(std::make_pair(0., std::make_pair(_ch_s, 0)));
  CH ch_t(-1);
  while (!que.empty() && !ch_t.is_valid())
  {
    auto dct_cur = *que.begin();
    que.erase(que.begin());

    auto ch_cur = dct_cur.second.first;
    int trans_cur = dct_cur.second.second;
    ALGOHEX_DEBUG_ONLY(std::cerr << " " << ch_cur;)

    auto chfs = mesh_.cell(ch_cur).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);

      auto c_next = mesh_.incident_cell(hf_opp);
      if (!c_next.is_valid())
        continue;

      //if not in the candidate cells
      if (_dp_cells.find(c_next) == _dp_cells.end())
        continue;

      //feature face normal should agree with the search direction
      auto fh_next = mesh_.face_handle(hf_opp);
      if (feature_fprop_[fh_next] > 0)
      {
        auto normal = mesh_.normal(hf_opp);
        auto ca = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_cur], normal);

        if (ca.second != tq_.axis_after_transition(e_ax, trans_cur))
        {
          ALGOHEX_DEBUG_ONLY(
                  std::cerr << "found feature face's normal tangential to search direction, continue..." << std::endl;)
          continue;
        }
      }


      //check if found
      int trans_next = tq_.mult_transitions_idx(trans_prop_[hf_opp], trans_cur);
      if (_target_chs.find(c_next) != _target_chs.end())
      {
        //TODO: verify axis if it reaches other tp
        if ((!to_tp && is_cell_found(c_next, trans_next, e_ax)) || to_tp)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "found " << c_next << std::endl;)
          ch_t = c_next;
          min_dist[c_next] = min_dist[ch_cur] + 1;
          previous[c_next] = ch_cur;
          break;
        }
        else
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "found in target set, but axis doesnot match " << std::endl;)
        }
      }

      //if contains special vertices
      if (_spv_free && (n_bad_vertices_in_cell(c_next, _excl_vh) > 0 || mesh_.is_boundary(c_next)))
        continue;

      int ax_next = tq_.axis_after_transition(e_ax, trans_next);

      //if not visited
      double dist = min_dist[ch_cur] + 1;
      if (min_dist.find(c_next) == min_dist.end() || dist < min_dist[c_next])
      {
        que.erase(std::make_pair(min_dist[c_next], std::make_pair(c_next, trans_next)));

        min_dist[c_next] = dist;
        previous[c_next] = ch_cur;

        que.insert(std::make_pair(min_dist[c_next], std::make_pair(c_next, trans_next)));
      }
    }
  }

  //find cells
  //get the path
  if (ch_t.is_valid())
  {
    std::vector<CH> rpath;
    CH ch_it = ch_t;
    while (ch_it.is_valid())
    {
      rpath.push_back(ch_it);
      ch_it = previous[ch_it];
    }

    std::vector<CH> path(rpath.rbegin(), rpath.rend());

    ALGOHEX_DEBUG_ONLY(std::cerr << "dpath: ";
                               for (const auto ch: path)
                                 std::cout << " " << ch;
                               std::cerr << std::endl;)

    return path;
  }

  return std::vector<CH>{};
}

//template<class MeshT>
//void
//FixZipperNodeT<MeshT>::
//save_cells(const std::vector<CH> &_v_chs, int id) const
//{
//  MeshT tm;
//  for (const auto chi: _v_chs)
//  {
//    auto cvhs = mesh_.get_cell_vertices(chi);
//    std::vector<VH> cvhs_new;
//    for (auto vhi: cvhs)
//      cvhs_new.push_back(tm.add_vertex(mesh_.vertex(vhi)));
//
//    tm.add_cell(cvhs_new);
//  }
//
//  OpenVolumeMesh::IO::FileManager fm;
//  fm.template writeFile("/path/cells" + std::to_string(id) + ".ovm",
//                        tm);
//}


template<class MeshT>
bool
FixZipperNodeT<MeshT>::
is_cell_found(const CH _ch, const int _trans, int _u) const
{
  if (!mesh_.is_boundary(_ch))
    return false;

  auto chfs = mesh_.cell(_ch).halffaces();
  for (int i = 0; i < 4; ++i)
  {
    auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);
    if (mesh_.is_boundary(hf_opp))
    {
      auto normal = mesh_.normal(hf_opp);
      auto ca = AxisAlignmentHelpers::closest_axis(cell_quaternions_[_ch], normal);

      if (ca.second == tq_.axis_after_transition(_u, _trans))
        return true;
    }
  }

  return false;
}

template<class MeshT>
bool
FixZipperNodeT<MeshT>::
is_cell_found_wrt_zipper_node(const CELLINFO &ci_cur, const std::set<VH> &_fixable_tps, VH &_other_tp, int _u)
{
  for (auto cv_it = mesh_.cv_iter(ci_cur.ch); cv_it.valid(); ++cv_it)
  {
    if (sgl_vt_[*cv_it] != 0 && !mesh_.is_boundary(*cv_it) &&
        //            !feature_face_vertex_[*cv_it] &&
        _fixable_tps.find(*cv_it) != _fixable_tps.end())
    {

      std::vector<HEH> cd_hehs;
      for (auto voh_it = mesh_.voh_iter(*cv_it); voh_it.valid(); ++voh_it)
      {
        EH vei = mesh_.edge_handle(*voh_it);
        if (!feature_face_edge_[vei] && std::fabs(valence_[vei]) == 1)
        {
          cd_hehs.push_back(*voh_it);
        }
      }

      if (cd_hehs.size() > 2)
        return false;

      auto or_chs = get_onering_cells(mesh_, *cv_it);
      std::set<HFH> ft_hfhs;
      auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, *cv_it, or_chs, ci_cur.ch, ft_hfhs);

      for (auto vohe: cd_hehs)
      {
        HFH hf_t = *mesh_.hehf_iter(vohe);
        CH ch_t = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_t));
        if (por_chs.find(ch_t) == por_chs.end())
          continue;

        int rt_ax = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                       valence_, vohe, hf_t);
        if (valence_[mesh_.edge_handle(vohe)] == 1)
          rt_ax = negate(rt_ax);

        int u_axis = tq_.axis_after_transition(_u, ci_cur.trans);

        int u_axis_in_t = fac_.axis_in_chs_expressed_in_cht(ci_cur.ch, ch_t, u_axis, por_chs);

        if (rt_ax == u_axis_in_t)
        {
          _other_tp = *cv_it;
          return true;
        }

      }
    }
  }

  return false;
}

template<class MeshT>
std::vector<CH>
FixZipperNodeT<MeshT>::get_special_vertex_free_path(const std::vector<CH> &_input_path, const HEH _heh_s,
                                                    const VH _other_tp)
{
  if (_input_path.empty())
  {
    std::cerr << "Warning: input dual path is empty!" << std::endl;
    return _input_path;
  }

  //split number
  int n_split_all = 0, n_split = 0;

  std::set<CH> dp_cells;
  for (const auto &ch: _input_path)
  {
    dp_cells.insert(ch);
  }

  //recover dual path
  std::vector<FH> dp_ffhs;
  for (auto i = 0u; i < _input_path.size() - 1; ++i)
  {
    auto fhi = common_face(mesh_, _input_path[i], _input_path[i + 1]);
    if (fhi.is_valid() && feature_fprop_[fhi] > 0)
    {
      dp_ffhs.push_back(fhi);
    }
  }


  CH ch_s = *_input_path.begin();
  //store target cells
  CH ch_t = _input_path.back();
  if (!_other_tp.is_valid() && !mesh_.is_boundary(ch_t))
    std::cerr << "Error: the last cell is not boundary" << std::endl;

  CellSplitT<MeshT> cs(mesh_);

  //split ch_t if it's at feature
  if (_other_tp.is_valid())
  {
    for (int j = 0; j < 2; ++j)
    {
      if (is_feature_cell(ch_t))
      {
        if (cs.is_split_ok(ch_t))
        {
          dp_cells.erase(ch_t);

          n_split++;
          auto vh_new = cs.cell_split(ch_t);

          std::set<CH> inc_cells;
          for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
            inc_cells.insert(*vc_it);

          ch_t = CH(-1);
          for (auto vc_it = mesh_.vc_iter(_other_tp); vc_it.valid(); ++vc_it)
            if (!is_feature_cell(*vc_it) && inc_cells.find(*vc_it) != inc_cells.end())
            {
              ch_t = *vc_it;
              break;
            }

          for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
            dp_cells.insert(*vc_it);

          //found a non feature cell
          if (ch_t.is_valid())
            break;
          else
          { //keep splitting
            for (auto vc_it = mesh_.vc_iter(_other_tp); vc_it.valid(); ++vc_it)
              if (inc_cells.find(*vc_it) != inc_cells.end())
              {
                ch_t = *vc_it;
                break;
              }
          }
        }
      }
      else
        break;
    }
  }

  //split ch_s if it's at feature
  FH other_ff(-1);
  FH first_dpf(-1);
  if (_input_path.size() > 1)
    first_dpf = common_face(mesh_, _input_path[0], _input_path[1]);

  if (first_dpf.is_valid())
  {
    for (auto cf_it = mesh_.cf_iter(ch_s); cf_it.valid(); ++cf_it)
    {
      if (*cf_it != first_dpf && feature_fprop_[*cf_it] > 0)
      {
        other_ff = *cf_it;
        break;
      }
    }
  }

  if (other_ff.is_valid())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << " split chs " << ch_s << std::endl;)

    if (cs.is_split_ok(ch_s))
    {
      dp_cells.erase(ch_s);

      n_split++;
      auto vh_new = cs.cell_split(ch_s);

      ch_s = CH(-1);
      for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
      {
        if (!is_feature_cell(*hec_it))
        {
          ch_s = *hec_it;
          break;
        }
      }

      for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
        dp_cells.insert(*vc_it);
    }
  }

  if (!ch_s.is_valid() || !ch_t.is_valid())
  {
    std::cerr << "Error: chs or cht is invalid " << ch_s << " " << ch_t << std::endl;
    return std::vector<CH>{};
  }

  //excluded cells from splitting
  std::set<CH> excl_chs;
  excl_chs.insert(ch_s);
  excl_chs.insert(ch_t);

  std::set<FH> excl_fhs;
  std::set<EH> excl_ehs;

  for (auto cf_it = mesh_.cf_iter(ch_t); cf_it.valid(); ++cf_it)
    if (feature_fprop_[*cf_it] > 0)
    {
      excl_fhs.insert(*cf_it);
    }
  for (auto fhi: dp_ffhs)
  {
    excl_fhs.insert(fhi);

    auto hfhi = mesh_.halfface_handle(fhi, 0);
    auto chi = mesh_.incident_cell(hfhi);
    if (chi.is_valid())
      excl_chs.insert(chi);

    auto hfhi_opp = mesh_.halfface_handle(fhi, 1);
    auto chi_opp = mesh_.incident_cell(hfhi_opp);
    if (chi_opp.is_valid())
      excl_chs.insert(chi_opp);
  }


  VH vh_s0 = mesh_.halfedge(_heh_s).from_vertex();
  VH vh_s1 = mesh_.halfedge(_heh_s).to_vertex();


  std::set<CH> target_chs{ch_t};
  std::set<CH> start_chs{ch_s};

  //w.r.t feature faces
  //split cells which contain feature face(s)
  std::vector<CH> split_chs;
  for (const auto &ch: dp_cells)
  {
    if (is_feature_cell(ch) && excl_chs.find(ch) == excl_chs.end())
    {
      split_chs.push_back(ch);
    }
  }
//        std::cerr<<std::endl;
  for (const auto &ch: split_chs)
  {
    if (cs.is_split_ok(ch))
    {
      dp_cells.erase(ch);

      n_split++;
      auto vh_new = cs.cell_split(ch);

      for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
        dp_cells.insert(*vc_it);
    }
  }
  n_split_all += n_split;
  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " feature cells " << std::endl;)
  n_split = 0;

  //split cells that have feature face edges
  std::set<EH> all_ehs;
  for (const auto &ch: dp_cells)
  {
    for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
      if (feature_face_edge_[*ce_it])
        all_ehs.insert(*ce_it);
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "ffes in dpcells ";
                             for (auto ehi: all_ehs)
                               std::cerr << ehi << " ";
                             std::cerr << std::endl;)
  //remove edges of target cell and start cell if it's boundary
  for (const auto chi: excl_chs)
    for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
      if (feature_face_edge_[*ce_it])
        all_ehs.erase(*ce_it);

  //split faces to get enough DOF
  std::vector<EH> fes_split;
  std::set<FH> ff_eh_fhs;
  for (const auto eh: all_ehs)
  {
    //if the incident feature face is on the dual path
    std::vector<FH> v_ffhs;
    for (auto ef_it = mesh_.ef_iter(eh); ef_it.valid(); ++ef_it)
    {
      auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(*ef_it, 0));
      auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(*ef_it, 1));

      //incident to cells on path
      if (dp_cells.find(ch0) != dp_cells.end() || dp_cells.find(ch1) != dp_cells.end())
      {
        //not dual path face
        if (feature_fprop_[*ef_it] == 0 || (feature_fprop_[*ef_it] > 0 && excl_fhs.find(*ef_it) == excl_fhs.end()))
          v_ffhs.push_back(*ef_it);
      }
    }

    if (!v_ffhs.empty())
    {
      fes_split.push_back(eh);
      for (auto ef_it = mesh_.ef_iter(eh); ef_it.valid(); ++ef_it)
      {
        auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(*ef_it, 0));
        auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(*ef_it, 1));

        if (dp_cells.find(ch0) != dp_cells.end() || dp_cells.find(ch1) != dp_cells.end())
          ff_eh_fhs.insert(*ef_it);
      }
    }
  }

  n_split += SplitHelperT<MeshT>::split_faces_with_cell_property_update(mesh_, fs_, dp_cells, ff_eh_fhs, start_chs,
                                                                        target_chs, excl_fhs);
  //remove the cells that have feature face edges
  for (const auto eh: fes_split)
  {
    for (auto ec_it = mesh_.ec_iter(eh); ec_it.valid(); ++ec_it)
      dp_cells.erase(*ec_it);
  }
  n_split_all += n_split;
  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " faces with feature face edges" << std::endl;)
  n_split = 0;

  //w.r.t singularity and feature
  std::set<VH> excl_vhs;
  excl_vhs.insert(vh_s0);
  excl_vhs.insert(vh_s1);
  if (_other_tp.is_valid())
    excl_vhs.insert(_other_tp);
  //split cells that have four singular/feature edge vertices
  split_chs.clear();
  for (const auto &ch: dp_cells)
    if (n_bad_vertices_in_cell(ch, excl_vhs) == 4)
      split_chs.push_back(ch);
  for (const auto &ch: split_chs)
    if (cs.is_split_ok(ch))
    {
      dp_cells.erase(ch);

      bool target_fd = target_chs.find(ch) != target_chs.end();

      n_split++;
      auto vh_new = cs.cell_split(ch);
      for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
        dp_cells.insert(*vc_it);

      //update target cells
      if (target_fd)
      {
        target_chs.erase(ch);

        for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
          if (!_other_tp.is_valid() && mesh_.is_boundary(*vc_it))
            target_chs.insert(*vc_it);
      }
    }

  n_split_all += n_split;
  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " cells " << std::endl;)
  n_split = 0;

  //split face which has three bad singular vertices
  std::set<FH> fhs3;
  for (const auto &ch: dp_cells)
  {
    if (n_bad_vertices_in_cell(ch, excl_vhs) == 3)
    {
      for (auto cf_it = mesh_.cf_iter(ch); cf_it.valid(); ++cf_it)
      {
        if (n_bad_vertices_in_face(*cf_it, excl_vhs) == 3)
          fhs3.insert(*cf_it);
      }
    }
  }

  n_split += SplitHelperT<MeshT>::split_faces_with_cell_property_update(mesh_, fs_, dp_cells, fhs3, start_chs,
                                                                        target_chs, excl_fhs);
  n_split_all += n_split;

  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " faces " << std::endl;)
  n_split = 0;

  //split regular edges that has two singular vertex
  std::set<EH> ehs2;
  for (const auto &ch: dp_cells)
  {
    if (n_bad_vertices_in_cell(ch, excl_vhs) == 2)
    {
      for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
      {
        if (valence_[*ce_it] == 0 && feature_edge_[*ce_it] == 0)
        {
          auto from_vh = mesh_.edge(*ce_it).from_vertex();
          auto to_vh = mesh_.edge(*ce_it).to_vertex();

          bool fv_bad = is_other_special_vertex(from_vh, excl_vhs);
          bool tv_bad = is_other_special_vertex(to_vh, excl_vhs);
          if (fv_bad && tv_bad)
            ehs2.insert(*ce_it);
        }
      }
    }
  }

  n_split += SplitHelperT<MeshT>::split_edges_with_cell_property_update(mesh_, es_, dp_cells, ehs2, start_chs,
                                                                        target_chs, excl_ehs);
  n_split_all += n_split;

  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " s2 edges " << std::endl;)
  n_split = 0;

  //split regular edges that has one singular or feature vertex
  std::set<EH> ehs1;
  for (const auto &ch: dp_cells)
  {
    if (n_bad_vertices_in_cell(ch, excl_vhs) >= 1)
    {
      for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
      {
        if (valence_[*ce_it] == 0 && feature_edge_[*ce_it] == 0)
        {
          auto from_vh = mesh_.edge(*ce_it).from_vertex();
          auto to_vh = mesh_.edge(*ce_it).to_vertex();

          bool fv_bad = is_other_special_vertex(from_vh, excl_vhs);
          bool tv_bad = is_other_special_vertex(to_vh, excl_vhs);
          if (fv_bad || tv_bad)
            ehs1.insert(*ce_it);
        }
      }
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "split s1 ehs: ";
                             for (auto ehi: ehs1)
                               std::cerr << " " << ehi;
                             std::cerr << std::endl;)
  if (ehs1.size() > 10000)
  {
    std::cerr << "Waring: too many s1 edges to split on this dual path. Try later!" << std::endl;
    return std::vector<CH>{};
  }
  n_split += SplitHelperT<MeshT>::split_edges_with_cell_property_update(mesh_, es_, dp_cells, ehs1, start_chs,
                                                                        target_chs, excl_ehs);
  n_split_all += n_split;

  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " s1 edges " << std::endl;
                             std::cerr << "split all " << n_split_all << " times " << std::endl;)

  if (n_split_all > 0)
  {
    CH ch_s;
    for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
      if (dp_cells.find(*hec_it) != dp_cells.end())
      {
        ch_s = *hec_it;
        break;
      }

    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " target chs: ";)
    //get target cells
    std::set<CH> real_target_chs;
    if (_other_tp.is_valid())
    {
      for (auto vc_it = mesh_.vc_iter(_other_tp); vc_it.valid(); ++vc_it)
      {
        if (target_chs.find(*vc_it) != target_chs.end() && n_bad_vertices_in_cell(*vc_it, excl_vhs) == 0)
        {
          real_target_chs.insert(*vc_it);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << *vc_it;)
        }
      }
      ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
    }
    else
    {
      for (const auto chi: target_chs)
      {
        if (mesh_.is_boundary(chi) && n_bad_vertices_in_cell(chi, excl_vhs) == 0)
        {
          real_target_chs.insert(chi);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << chi;)
        }
      }
      ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
    }

    if (ch_s.is_valid() && !real_target_chs.empty())
      return shortest_dual_path_from_cell_to_target(_heh_s, ch_s, real_target_chs, dp_cells, excl_vhs, true);
    else
      return std::vector<CH>{};
  }

  return _input_path;
}

template<class MeshT>
std::vector<CH>
FixZipperNodeT<MeshT>::get_handle_free_dual_path(const std::vector<CH> &_input_path, const HEH _heh_s,
                                                 const VH _other_tp)
{
  if (_input_path.empty())
  {
    return _input_path;
  }

  VH vh_s0 = mesh_.halfedge(_heh_s).from_vertex();
  VH vh_s1 = mesh_.halfedge(_heh_s).to_vertex();

  std::set<VH> excl_vhs;
  excl_vhs.insert(vh_s0);
  excl_vhs.insert(vh_s1);
  if (_other_tp.is_valid())
    excl_vhs.insert(_other_tp);

  std::set<FH> blocking_fhs;
  std::map<FH, int> step_f;

  for (auto i = 0u; i < _input_path.size(); ++i)
    for (auto cf_it = mesh_.cf_iter(_input_path[i]); cf_it.valid(); ++cf_it)
      step_f[*cf_it] = -1;
  for (auto i = 0u; i < _input_path.size(); ++i)
  {
    for (auto cf_it = mesh_.cf_iter(_input_path[i]); cf_it.valid(); ++cf_it)
    {
      if (step_f[*cf_it] == -1 || step_f[*cf_it] == (int) i - 1)
      {
        step_f[*cf_it] = i;
      }
      else
      {
        blocking_fhs.insert(*cf_it);
      }
    }
  }

  std::vector<CH> output_path = _input_path;
  if (!blocking_fhs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "blocking fhs: ";
                               for (auto fhi: blocking_fhs)
                                 std::cerr << " " << fhi;
                               std::cerr << std::endl;)

    //store target cells
    CH ch_s = output_path[0];
    //cannot be split
    CH ch_t = output_path.back();
    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " ch_t: " << ch_t << " ehs: " << mesh_.edge_handle(_heh_s);)

    std::set<CH> start_chs{ch_s};

    CellSplitT<MeshT> cs(mesh_);

    //split cells which contain blocking face(s)
    ALGOHEX_DEBUG_ONLY(std::cerr << "split cell: ";)
    std::set<CH> split_chs;
    for (const auto &fhi: blocking_fhs)
    {
      auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 0));
      auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 1));

      split_chs.insert(ch0);
      split_chs.insert(ch1);
    }

    std::set<CH> dp_cells;
    for (const auto &chi: output_path)
      dp_cells.insert(chi);

    ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
    int n_split = 0;
    for (const auto &ch: split_chs)
    {
      if (cs.is_split_ok(ch))
      {
        dp_cells.erase(ch);

        n_split++;
        auto vh_new = cs.cell_split(ch);

        for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
        {
          dp_cells.insert(*vc_it);
        }

        if (ch == ch_s)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "Split start cell " << ch << std::endl;)
          start_chs.erase(ch);
          for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
          {
            start_chs.insert(*vc_it);
          }
        }

        if (ch == ch_t)
        {
          std::cerr << "Warning: handle found at end cell " << ch << std::endl;
        }
      }
    }

    for (const auto &fhi: blocking_fhs)
    {
      auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 0));
      auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 1));

      dp_cells.erase(ch0);
      dp_cells.erase(ch1);

      start_chs.erase(ch0);
      start_chs.erase(ch1);
    }
    ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " blocking cells. " << std::endl;)

    //find start cell
    ch_s.reset();
    bool bdy_se = mesh_.is_boundary(_heh_s);
    for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
    {
      if (start_chs.find(*hec_it) != start_chs.end())
      {
        if (bdy_se)
        {
          if (mesh_.is_boundary(*hec_it))
          {
            ch_s = *hec_it;
            break;
          }
        }
        else
        {
          ch_s = *hec_it;
          break;
        }
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " target chs: " << ch_t;)
    //get target cells
    std::set<CH> real_target_chs{ch_t};

    if (ch_s.is_valid() && !real_target_chs.empty())
      output_path = shortest_dual_path_from_cell_to_target(_heh_s, ch_s, real_target_chs, dp_cells, excl_vhs, false);
    else
      return std::vector<CH>{};
  }

  std::set<EH> blocking_ehs;
  std::map<EH, int> step_e;

  for (auto i = 0u; i < output_path.size(); ++i)
    for (auto ce_it = mesh_.ce_iter(output_path[i]); ce_it.valid(); ++ce_it)
      step_e[*ce_it] = -1;
  for (auto i = 0u; i < output_path.size(); ++i)
  {
    for (auto ce_it = mesh_.ce_iter(output_path[i]); ce_it.valid(); ++ce_it)
    {
      if (step_e[*ce_it] == -1 || step_e[*ce_it] == (int) i - 1)
      {
        step_e[*ce_it] = i;
      }
      else
      {
        blocking_ehs.insert(*ce_it);
      }
    }
  }

  if (!blocking_ehs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "blocking ehs: ";
                               for (auto ehi: blocking_ehs)
                                 std::cerr << " " << ehi;
                               std::cerr << std::endl;)

    //store target cells
    CH ch_s = output_path[0];
    CH ch_t = output_path.back();
    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " ch_t: " << ch_t << " ehs: " << mesh_.edge_handle(_heh_s);)

    std::set<CH> target_chs{ch_t};
    std::set<CH> start_chs{ch_s};

    //find faces incident to blocking ehs and the faces are not in ch_s
    std::set<FH> cell_fhs;
    for (const auto &chi: output_path)
      for (auto cf_it = mesh_.cf_iter(chi); cf_it.valid(); ++cf_it)
        cell_fhs.insert(*cf_it);

    //exclude edges of ch_s and ch_t
    for (auto cf_it = mesh_.cf_iter(ch_s); cf_it.valid(); ++cf_it)
      cell_fhs.erase(*cf_it);

    //faces to split
    std::set<FH> fhs_split;
    for (const auto &ehi: blocking_ehs)
    {
      for (auto ef_it = mesh_.ef_iter(ehi); ef_it.valid(); ++ef_it)
        if (cell_fhs.find(*ef_it) != cell_fhs.end())
          fhs_split.insert(*ef_it);
    }

    std::set<CH> dp_cells;
    for (const auto &chi: output_path)
      dp_cells.insert(chi);


    //split faces
    std::set<FH> excl_fhs;
    int n_split = SplitHelperT<MeshT>::split_faces_with_cell_property_update(mesh_, fs_, dp_cells, fhs_split, start_chs,
                                                                             target_chs, excl_fhs);
    ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " faces " << std::endl;)

    //remove cells incident to blocking edges
    std::set<CH> protected_chs;
    protected_chs.insert(ch_s);
    for (auto chf_it = mesh_.chf_iter(ch_s); chf_it.valid(); ++chf_it)
    {
      auto hfh_opp = mesh_.opposite_halfface_handle(*chf_it);
      auto ch_opp = mesh_.incident_cell(hfh_opp);
      if (ch_opp.is_valid() && dp_cells.find(ch_opp) != dp_cells.end())
        protected_chs.insert(ch_opp);
    }

    for (const auto ehi: blocking_ehs)
      for (auto ec_it = mesh_.ec_iter(ehi); ec_it.valid(); ++ec_it)
      {
        if (protected_chs.find(*ec_it) == protected_chs.end())
          dp_cells.erase(*ec_it);
      }

    ALGOHEX_DEBUG_ONLY(std::cerr << "dp cells after split: ";
                               for (const auto &ch: dp_cells)
                               {
                                 std::cerr << " " << ch;
                               }
                               std::cerr << std::endl;)

    //find start cell
    ch_s.reset();
    bool bdy_se = mesh_.is_boundary(_heh_s);
    for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
    {
      if (start_chs.find(*hec_it) != start_chs.end())
      {
        if (bdy_se)
        {
          if (mesh_.is_boundary(*hec_it))
          {
            ch_s = *hec_it;
            break;
          }
        }
        else
        {
          ch_s = *hec_it;
          break;
        }
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " target chs: ";)
    //get target cells
    std::set<CH> real_target_chs;
    if (_other_tp.is_valid())
    {
      for (auto vc_it = mesh_.vc_iter(_other_tp); vc_it.valid(); ++vc_it)
      {
        if (target_chs.find(*vc_it) != target_chs.end())
        {
          //contains no blocking vertex
          bool has_bl_ehs = false;
          for (auto ce_it = mesh_.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
          {
            if (blocking_ehs.find(*ce_it) != blocking_ehs.end())
            {
              has_bl_ehs = true;
              break;
            }
          }

          if (!has_bl_ehs)
          {
            real_target_chs.insert(*vc_it);
            ALGOHEX_DEBUG_ONLY(std::cerr << " " << *vc_it;)
          }
        }
      }
    }
    else
    {
      for (const auto chi: target_chs)
      {
        if (mesh_.is_boundary(chi))
        {
          //contains no blocking vertex
          bool has_bl_ehs = false;
          for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
          {
            if (blocking_ehs.find(*ce_it) != blocking_ehs.end())
            {
              has_bl_ehs = true;
              break;
            }
          }

          if (!has_bl_ehs)
          {
            real_target_chs.insert(chi);
            ALGOHEX_DEBUG_ONLY(std::cerr << " " << chi;)
          }
        }
      }
      ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
    }

    if (ch_s.is_valid() && !real_target_chs.empty())
      output_path = shortest_dual_path_from_cell_to_target(_heh_s, ch_s, real_target_chs, dp_cells, excl_vhs, false);
    else
      return std::vector<CH>{};
  }


  //try to find blocking vertices
  std::set<VH> blocking_vhs;
  std::map<VH, int> step_v;

  for (auto i = 0u; i < output_path.size(); ++i)
    for (auto cv_it = mesh_.cv_iter(output_path[i]); cv_it.valid(); ++cv_it)
      step_v[*cv_it] = -1;


  for (auto i = 0u; i < output_path.size(); ++i)
  {
    for (auto cv_it = mesh_.cv_iter(output_path[i]); cv_it.valid(); ++cv_it)
    {
      if (step_v[*cv_it] == -1 || step_v[*cv_it] == (int) i - 1)
      {
        step_v[*cv_it] = i;
      }
      else
      {
        blocking_vhs.insert(*cv_it);
      }
    }
  }


  if (!blocking_vhs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "blocking vhs: ";
                               for (auto vhi: blocking_vhs)
                                 std::cerr << " " << vhi;
                               std::cerr << std::endl;)



    //store target cells
    CH ch_s = output_path[0];
    CH ch_t = output_path.back();
    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " ch_t: " << ch_t << " ehs: " << mesh_.edge_handle(_heh_s);)

    std::set<CH> target_chs{ch_t};
    std::set<CH> start_chs{ch_s};

    //find edges incident to blocking vhs and the edges are not in ch_s
    std::set<EH> cell_ehs;
    for (const auto &chi: output_path)
      for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
        cell_ehs.insert(*ce_it);

    //exclude edges of ch_s and ch_t
    for (auto ce_it = mesh_.ce_iter(ch_s); ce_it.valid(); ++ce_it)
      cell_ehs.erase(*ce_it);
//            for (auto ce_it = mesh_.ce_iter(ch_t); ce_it.valid(); ++ce_it)
//                cell_ehs.erase(*ce_it);

    //edges to split
    std::set<EH> ehs1;
    for (const auto &vhi: blocking_vhs)
    {
      for (auto ve_it = mesh_.ve_iter(vhi); ve_it.valid(); ++ve_it)
        if (cell_ehs.find(*ve_it) != cell_ehs.end())
          ehs1.insert(*ve_it);
    }

    std::set<CH> dp_cells;
    for (const auto &chi: output_path)
      dp_cells.insert(chi);

    //split
    std::set<EH> excl_ehs;
    int n_split = SplitHelperT<MeshT>::split_edges_with_cell_property_update(mesh_, es_, dp_cells, ehs1, start_chs,
                                                                             target_chs, excl_ehs);
    ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " edges " << std::endl;)

    //remove cells incident to blocking vhs if not
    std::set<CH> protected_chs;
    protected_chs.insert(ch_s);
    for (auto chf_it = mesh_.chf_iter(ch_s); chf_it.valid(); ++chf_it)
    {
      auto hfh_opp = mesh_.opposite_halfface_handle(*chf_it);
      auto ch_opp = mesh_.incident_cell(hfh_opp);
      if (ch_opp.is_valid() && dp_cells.find(ch_opp) != dp_cells.end())
        protected_chs.insert(ch_opp);
    }
    for (const auto vhi: blocking_vhs)
      for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
      {
        if (protected_chs.find(*vc_it) == protected_chs.end())
          dp_cells.erase(*vc_it);
      }


    //find start cell
    ch_s.reset();

    for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
    {
      if (start_chs.find(*hec_it) != start_chs.end())
      {
        if (mesh_.is_boundary(_heh_s) && mesh_.is_boundary(*hec_it))
          ch_s = *hec_it;
        else
          ch_s = *hec_it;
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " target chs: ";)
    //get target cells
    std::set<CH> real_target_chs;
    if (_other_tp.is_valid())
    {
      for (auto vc_it = mesh_.vc_iter(_other_tp); vc_it.valid(); ++vc_it)
      {
        if (target_chs.find(*vc_it) != target_chs.end())
        {
          //contains no blocking vertex
          bool has_bl_vts = false;
          for (auto cv_it = mesh_.cv_iter(*vc_it); cv_it.valid(); ++cv_it)
          {
            if (blocking_vhs.find(*cv_it) != blocking_vhs.end())
            {
              has_bl_vts = true;
              break;
            }
          }

          if (!has_bl_vts)
          {
            real_target_chs.insert(*vc_it);
            ALGOHEX_DEBUG_ONLY(std::cerr << " " << *vc_it;)
          }
        }
      }
      ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
    }
    else
    {
      for (const auto chi: target_chs)
      {
        if (mesh_.is_boundary(chi))
        {
          //contains no blocking vertex
          bool has_bl_vts = false;
          for (auto cv_it = mesh_.cv_iter(chi); cv_it.valid(); ++cv_it)
          {
            if (blocking_vhs.find(*cv_it) != blocking_vhs.end())
            {
              has_bl_vts = true;
              break;
            }
          }

          if (!has_bl_vts)
          {
            real_target_chs.insert(chi);
            ALGOHEX_DEBUG_ONLY(std::cerr << " " << chi;)
          }
        }
      }
      ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
    }


    if (ch_s.is_valid() && !real_target_chs.empty())
      return shortest_dual_path_from_cell_to_target(_heh_s, ch_s, real_target_chs, dp_cells, excl_vhs, false);
    else
      return std::vector<CH>{};

  }


  return output_path;
}


//eset: path should be subset
//vset: path not touch the vertex of the set
template<class MeshT>
std::vector<EH>
FixZipperNodeT<MeshT>::get_edge_path_on_cell_path_wrt_field(const std::vector<CH> &_dpath,
                                                            const std::set<EH> &_excluded_ehs,
                                                            const VH _vh_s, const HEH _heh_s,
                                                            const std::set<VH> &_target_vhs, const bool _opp_dir)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "finding epath wrt difference vector " << std::endl;)
  std::set<CH> all_cells;
  all_cells.insert(_dpath.begin(), _dpath.end());

  HFH hfh_s;
  for (auto hehf_it = mesh_.hehf_iter(_heh_s); hehf_it.valid(); ++hehf_it)
  {
    if (all_cells.find(mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it))) != all_cells.end())
    {
      hfh_s = *hehf_it;
      break;
    }
  }
  CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh_s));
  ALGOHEX_DEBUG_ONLY(std::cerr << "heh_s: " << _heh_s << " hf_s: " << hfh_s << " ch_s: " << ch_s << std::endl;)

  auto cell_axes = propagate_uvw_in_cells(_heh_s, hfh_s, all_cells);

  auto uvw_orig = cell_axes[ch_s];

  //get difference vector expressed in pm coordinates of ch_s
  Vec3d diff_vec = get_difference_vector_at_zipper_node(_heh_s, ch_s, _dpath, cell_axes);

  if (_opp_dir)
    diff_vec = -diff_vec;

  std::map<VH, double> min_dist;
  //avoid handle
  std::map<VH, VH> previous;
  min_dist.insert(std::make_pair(_vh_s, 0.));
  previous.insert(std::make_pair(_vh_s, VH(-1)));

  VH vh_ee = mesh_.halfedge(_heh_s).to_vertex();
  Vec3d u_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_s], (AxisAlignment) uvw_orig[0]);
  Vec3d v_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_s], (AxisAlignment) uvw_orig[1]);
  Vec3d w_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_s], (AxisAlignment) uvw_orig[2]);
  Vec3d dv_in_g =
          diff_vec[uvw_orig[0] / 2] * u_dir + diff_vec[uvw_orig[1] / 2] * v_dir + diff_vec[uvw_orig[2] / 2] * v_dir;


  //edge path shouldn't touch these vertices
  std::set<VH> excluded_vhs;
  for (const auto ehi: _excluded_ehs)
  {
    if (!mesh_.is_boundary(ehi))
    {
      excluded_vhs.insert(mesh_.edge(ehi).from_vertex());
      excluded_vhs.insert(mesh_.edge(ehi).to_vertex());
    }
  }


  VH vh_e(-1);
  VH vh_it = _vh_s;
  std::set<VH> visited_vhs;
  visited_vhs.insert(vh_it);
  for (auto i = 0u; i < _dpath.size(); ++i)
  {
    auto uvw = cell_axes[_dpath[i]];

    double min_dist = 10000000000.;
    VH vh_min(-1);

    std::set<EH> cehs;
    for (auto ce_it = mesh_.ce_iter(_dpath[i]); ce_it.valid(); ++ce_it)
    {
      cehs.insert(*ce_it);
    }

    std::set<VH> cd_vhs;

    if (i < _dpath.size() - 1)
    {
      FH fh_common = common_face(mesh_, _dpath[i], _dpath[i + 1]);
      for (auto fv_it = mesh_.fv_iter(fh_common); fv_it.valid(); ++fv_it)
        cd_vhs.insert(*fv_it);
    }
    else if (i == _dpath.size() - 1)
      cd_vhs = _target_vhs;

    if (cd_vhs.find(vh_it) != cd_vhs.end())
      continue;

    for (auto voh_it = mesh_.voh_iter(vh_it); voh_it.valid(); ++voh_it)
    {
      auto ceh = mesh_.edge_handle(*voh_it);
      if (cehs.find(ceh) != cehs.end())
      {
        if (_excluded_ehs.find(ceh) != _excluded_ehs.end())
          continue;

        auto vh_t = mesh_.halfedge(*voh_it).to_vertex();
        if (cd_vhs.find(vh_t) == cd_vhs.end())
          continue;

        if (visited_vhs.find(vh_t) != visited_vhs.end())
          continue;

        if (excluded_vhs.find(vh_t) != excluded_vhs.end())
          continue;

        Vec3d u_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[_dpath[i]], (AxisAlignment) uvw[0]);
        Vec3d v_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[_dpath[i]], (AxisAlignment) uvw[1]);
        Vec3d w_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[_dpath[i]], (AxisAlignment) uvw[2]);


        auto e_dir = ovm2eigen(mesh_.vertex(vh_t) - mesh_.vertex(vh_it));

        //express the edge in the pm coordinates of ch_s
        Vec3d e_dir_param;
        e_dir_param[uvw_orig[0] / 2] = e_dir.dot(u_dir);
        e_dir_param[uvw_orig[1] / 2] = e_dir.dot(v_dir);
        e_dir_param[uvw_orig[2] / 2] = e_dir.dot(w_dir);

        double dprod = e_dir_param.dot(diff_vec);
        if (dprod < min_dist)
        {
          min_dist = dprod;
          vh_min = vh_t;
          if (_target_vhs.find(vh_t) != _target_vhs.end())
          {
            vh_e = vh_t;
          }
        }
      }
    }

    previous[vh_min] = vh_it;
    vh_it = vh_min;
    visited_vhs.insert(vh_min);
  }

  std::vector<EH> path;
  if (vh_e.is_valid())
  {
    //get the path
    std::vector<VH> rvpath;
    VH vh_it = vh_e;
    while (vh_it.is_valid())
    {
      rvpath.push_back(vh_it);

      vh_it = previous[vh_it];
    }

    if (!rvpath.empty())
    {
      path.reserve(rvpath.size() - 1);
      for (int i = rvpath.size() - 1; i > 0; --i)
        path.push_back(mesh_.edge_handle(mesh_.find_halfedge(rvpath[i], rvpath[i - 1])));

      ALGOHEX_DEBUG_ONLY(std::cerr << " edge path: ";
                                 for (const auto eh: path)
                                   std::cout << " " << eh;
                                 std::cerr << std::endl;)
    }
    else
      std::cout << "Error: could not find edge path!" << std::endl;
  }


  return path;
}

template<class MeshT>
std::set<EH>
FixZipperNodeT<MeshT>::remedy_edge_path(const std::vector<CH> &_dpath, const std::vector<EH> &_raw_e_path) const
{
  std::set<FH> path_fhs;
  for (auto i = 0u; i < _dpath.size() - 1; ++i)
  {
    FH cfh = common_face(mesh_, _dpath[i], _dpath[i + 1]);
    assert(cfh.is_valid());
    path_fhs.insert(cfh);
  }

  std::set<EH> path_ehs;
  for (const auto ehi: _raw_e_path)
    path_ehs.insert(ehi);

  std::queue<FH> que;
  for (const auto fhi: path_fhs)
    que.push(fhi);

  int n = 0;
  while (!que.empty() && n < 10000)
  {
    auto fh_cur = que.front();
    que.pop();


    EH np_eh(-1);
    int n_pe = n_path_edges_in_face(fh_cur, path_ehs, np_eh);

    if (n_pe == 2)
    {
      for (auto fe_it = mesh_.fe_iter(fh_cur); fe_it.valid(); ++fe_it)
      {
        if (*fe_it == np_eh)
          path_ehs.insert(*fe_it);
        else
          path_ehs.erase(*fe_it);
      }

      for (auto ef_it = mesh_.ef_iter(np_eh); ef_it.valid(); ++ef_it)
        if (path_fhs.find(*ef_it) != path_fhs.end())
          que.push(*ef_it);

      n++;
    }
  }

  //check
  for (auto fhi: path_fhs)
  {
    EH np_eh;
    int n_pe = n_path_edges_in_face(fhi, path_ehs, np_eh);
    if (n_pe > 1)
    {
      std::cerr << "Error: face " << fhi << " has " << n_pe << " path edges" << std::endl;
    }
  }

  return path_ehs;
}

//_all_eset: path should be subset
//_excl_eset: path not touch the vertices of the eset
template<class MeshT>
std::vector<EH>
FixZipperNodeT<MeshT>::shortest_path_v_to_target(const std::set<EH> &_all_eset, const std::set<EH> &_excl_eset,
                                                 const VH _vh_f, const std::set<VH> &_target_vhs,
                                                 const VH _other_tp) const
{
  std::set<VH> excl_vhs;
  for (const auto ehi: _excl_eset)
  {
    excl_vhs.insert(mesh_.edge(ehi).from_vertex());
    excl_vhs.insert(mesh_.edge(ehi).to_vertex());
  }
  excl_vhs.erase(_other_tp);

  std::map<VH, double> min_dist;
  std::map<VH, VH> previous;
  min_dist.insert(std::make_pair(_vh_f, 0.));
  previous.insert(std::make_pair(_vh_f, VH(-1)));
//std::cerr<<"start vh: "<<_vh_f<<" go: ";
  std::set<std::pair<double, VH>> que;
  que.insert(std::make_pair(0., _vh_f));
  VH vh_t(-1);
  while (!que.empty() && !vh_t.is_valid())
  {
    auto dv_cur = *que.begin();
    que.erase(que.begin());

    //find vertex in target vertex set
    if (_target_vhs.find(dv_cur.second) != _target_vhs.end())
    {
      vh_t = dv_cur.second;
      break;
    }

    for (auto voh_it = mesh_.voh_iter(dv_cur.second); voh_it.valid(); ++voh_it)
    {
      auto vvh = mesh_.halfedge(*voh_it).to_vertex();

      if (_all_eset.find(mesh_.edge_handle(*voh_it)) != _all_eset.end() && excl_vhs.find(vvh) == excl_vhs.end())
      {
        double dist = min_dist[dv_cur.second] + mesh_.length(*voh_it);

        if (min_dist.find(vvh) == min_dist.end() || dist < min_dist[vvh])
        {
          que.erase(std::make_pair(min_dist[vvh], vvh));

          min_dist[vvh] = dist;
          previous[vvh] = dv_cur.second;

          que.insert(std::make_pair(min_dist[vvh], vvh));
        }
      }
    }
  }
  std::vector<EH> path;

  if (vh_t.is_valid())
  {
    //get the path
    std::vector<VH> rvpath;

    VH vh_it = vh_t;
    while (vh_it.is_valid())
    {
      rvpath.push_back(vh_it);
      vh_it = previous[vh_it];
    }

    if (!rvpath.empty())
    {
      path.reserve(rvpath.size() - 1);
      for (int i = rvpath.size() - 1; i > 0; --i)
        path.push_back(mesh_.edge_handle(mesh_.find_halfedge(rvpath[i], rvpath[i - 1])));
    }
    else
      std::cout << "Error: could not find edge path!" << std::endl;
  }

  return path;
}


template<class MeshT>
int
FixZipperNodeT<MeshT>::split_for_second_edge_path(const std::vector<CH> &_dpath, const std::set<EH> &_raw_e_path,
                                                  const EH _sg_eh, std::set<CH> &_dp_cells)
{
  std::set<FH> path_fhs, all_fhs;
  for (const auto chi: _dpath)
  {
    for (auto cf_it = mesh_.cf_iter(chi); cf_it.valid(); ++cf_it)
      path_fhs.insert(*cf_it);
  }

  std::set<EH> path_ehs;
  for (const auto ehi: _raw_e_path)
    path_ehs.insert(ehi);

  std::set<VH> path_vhs;
  for (const auto ehi: _raw_e_path)
  {
    path_vhs.insert(mesh_.edge(ehi).from_vertex());
    path_vhs.insert(mesh_.edge(ehi).to_vertex());
  }

  std::set<EH> split_ehs;
  for (const auto fhi: path_fhs)
  {
    int nfv = n_path_vertices_in_face(fhi, path_vhs);
    if (nfv == 3)
    {
      for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
      {
        if (path_ehs.find(*fe_it) == path_ehs.end() && *fe_it != _sg_eh)
          split_ehs.insert(*fe_it);
      }
    }
  }

  _dp_cells.insert(_dpath.begin(), _dpath.end());
  std::set<CH> target_chs, start_cells;
  std::set<EH> excl_ehs;
  SplitHelperT<MeshT>::split_edges_with_cell_property_update(mesh_, es_, _dp_cells, split_ehs, start_cells, target_chs,
                                                             excl_ehs);


  return (int) split_ehs.size();
}


template<class MeshT>
std::map<CH, typename FixZipperNodeT<MeshT>::UVW>
FixZipperNodeT<MeshT>::propagate_uvw_in_cells(const HEH _heh_s, const HFH _hfh_s, const std::set<CH> &_cells)
{
  std::map<CH, UVW> cell_axes;

  if (!_hfh_s.is_valid())
    return cell_axes;

  CH ch_start = mesh_.incident_cell(mesh_.opposite_halfface_handle(_hfh_s));

  //rotational axis direction
  int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                   _heh_s, _hfh_s);
  //if it's not a singular edge
  if (rt_axis == -1)
  {
    VH vhf = mesh_.halfedge(_heh_s).from_vertex();
    VH vht = mesh_.halfedge(_heh_s).to_vertex();
    Point dir = mesh_.vertex(vht) - mesh_.vertex(vhf);
    dir.normalize();

    rt_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_start], dir).second;
  }

  //edge direction in parametrization coordinates
  int u = rt_axis;
  int valence = valence_[mesh_.edge_handle(_heh_s)];
  if (valence == 1)
    u = negate(u);
  int v = (u + 2) % 6;
  int w = the_third_axis((AxisAlignment) u, (AxisAlignment) v);

  return propagate_uvw_in_cells(ch_start, UVW{u, v, w}, _cells);
}

template<class MeshT>
std::map<CH, typename FixZipperNodeT<MeshT>::UVW>
FixZipperNodeT<MeshT>::propagate_uvw_in_cells(const CH _ch_s, const UVW &_uvw, const std::set<CH> &_cells)
{
  std::map<CH, UVW> cell_axes;

  cell_axes.insert(std::make_pair(_ch_s, _uvw));

  // process dual spanning tree
  std::queue<HFH> qhfh;
  // add seed
  for (unsigned int i = 0; i < mesh_.cell(_ch_s).halffaces().size(); ++i)
    qhfh.push(mesh_.cell(_ch_s).halffaces()[i]);

  // grow tree until all cells in one ring visited
  while (!qhfh.empty())
  {
    // get opposite of next halfface
    auto hfh0 = qhfh.front();
    auto hfh1 = mesh_.opposite_halfface_handle(hfh0);

    auto fh = mesh_.face_handle(hfh0);
    auto ch_cur = mesh_.incident_cell(hfh0);
    qhfh.pop();

    // only process if not on boundary and in shell chs
    if (!mesh_.is_boundary(fh))
    {
      // check if cell already visited
      CH ch_new = mesh_.incident_cell(hfh1);
      if (_cells.find(ch_new) != _cells.end() && cell_axes.count(ch_new) == 0)
      {
        // get matching
        int trans = trans_prop_[hfh1];

        // transform axis
        auto u_new = tq_.axis_after_transition(cell_axes[ch_cur][0], trans);
        auto v_new = tq_.axis_after_transition(cell_axes[ch_cur][1], trans);
        auto w_new = tq_.axis_after_transition(cell_axes[ch_cur][2], trans);

        cell_axes[ch_new] = UVW{u_new, v_new, w_new};

        // process halffaces to neighbors
        for (auto i = 0u; i < mesh_.cell(ch_new).halffaces().size(); ++i)
        {
          // get new half-face
          auto hfh_new = mesh_.cell(ch_new).halffaces()[i];
          // push as candidate
          qhfh.push(hfh_new);
        }
      }
    }
  }

  return cell_axes;
}


//1.
template<class MeshT>
Vec3d FixZipperNodeT<MeshT>::
get_difference_vector_at_zipper_node(const HEH _heh_s, const CH _ch_s, const std::vector<CH> &_input_path,
                                     std::map<CH, UVW> &_cell_axis)
{
  VH vht = mesh_.halfedge(_heh_s).to_vertex();
  VH vhf = mesh_.halfedge(_heh_s).from_vertex();
  EH eh_s = mesh_.edge_handle(_heh_s);
  HEH sgheh_other(-1);
  for (auto voh_it = mesh_.voh_iter(vht); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh_.edge_handle(*voh_it);
    if (valence_[ve] != 0 && ve != eh_s)
    {
      sgheh_other = *voh_it;
      break;
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "chs " << _ch_s << " vht " << vht << " node id " <<
                               MeshPropertiesT<MeshT>::node_index(vht) << " val eh_s "
                               << valence_[mesh_.edge_handle(_heh_s)] << " heh other " << sgheh_other << std::endl;)
//
  //find in the dual path the furthest cell that is in the onering
  std::set<CH> or_chs;
  for (auto vc_it = mesh_.vc_iter(vht); vc_it.valid(); ++vc_it)
    or_chs.insert(*vc_it);

  //first vec
  auto uvw = _cell_axis[_ch_s];
  Vec3d u_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[_ch_s], (AxisAlignment) uvw[0]);
  Vec3d v_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[_ch_s], (AxisAlignment) uvw[1]);
  Vec3d w_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[_ch_s], (AxisAlignment) uvw[2]);

  auto e_dir0 = ovm2eigen(mesh_.vertex(vhf) - mesh_.vertex(vht));

  Vec3d e_dir_param;
  e_dir_param[uvw[0] / 2] = e_dir0.dot(u_dir);
  e_dir_param[uvw[1] / 2] = e_dir0.dot(v_dir);
  e_dir_param[uvw[2] / 2] = e_dir0.dot(w_dir);

  //seecond vec
  auto or_cell_axis = propagate_uvw_in_cells(_ch_s, _cell_axis[_ch_s], or_chs);

  CH ch2 = *mesh_.hec_iter(sgheh_other);
  auto uvw2 = or_cell_axis[ch2];
  Vec3d u_dir2 = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch2], (AxisAlignment) uvw2[0]);
  Vec3d v_dir2 = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch2], (AxisAlignment) uvw2[1]);
  Vec3d w_dir2 = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch2], (AxisAlignment) uvw2[2]);

  auto e_dir2 = ovm2eigen(mesh_.vertex(mesh_.halfedge(sgheh_other).to_vertex()) - mesh_.vertex(vht));
  Vec3d e_dir_param2;
  //expressed in ch_s
  e_dir_param2[uvw[0] / 2] = e_dir2.dot(u_dir2);
  e_dir_param2[uvw[1] / 2] = e_dir2.dot(v_dir2);
  e_dir_param2[uvw[2] / 2] = e_dir2.dot(w_dir2);

  //difference vector in param coordinates
  Vec3d v_diff = e_dir_param2 - e_dir_param;
  v_diff[uvw[0] / 2] = 0.;

  //project
  return v_diff;
}

//eset: path should be subset
//vset: path not touch the vertex of the set
template<class MeshT>
std::vector<EH>
FixZipperNodeT<MeshT>::get_edge_path_on_cell_path(const std::vector<CH> &_dpath, const std::set<EH> &_excluded_ehs,
                                                  const VH _vh_s, const std::set<VH> &_target_vhs) const
{
  std::vector<double> min_dist(mesh_.n_vertices(), std::numeric_limits<double>::infinity());
  //avoid handle
  std::vector<VH> previous(mesh_.n_vertices(), VH(-1));
  min_dist[_vh_s.idx()] = 0.;

  //edge path shouldn't touch these vertices
  std::set<VH> excluded_vhs;
  for (const auto ehi: _excluded_ehs)
  {
    if (!mesh_.is_boundary(ehi))
    {
      excluded_vhs.insert(mesh_.edge(ehi).from_vertex());
      excluded_vhs.insert(mesh_.edge(ehi).to_vertex());
    }
  }


  VH vh_e(-1);
  for (auto i = 0u; i < _dpath.size(); ++i)
  {
    for (auto che_it = mesh_.che_iter(_dpath[i]); che_it.valid(); ++che_it)
    {
      auto ceh = mesh_.edge_handle(*che_it);
      if (_excluded_ehs.find(ceh) != _excluded_ehs.end())
        continue;

      auto vh_f = mesh_.halfedge(*che_it).from_vertex();
      auto vh_t = mesh_.halfedge(*che_it).to_vertex();

      if ((excluded_vhs.find(vh_f) != excluded_vhs.end() && vh_f != _vh_s) ||
          excluded_vhs.find(vh_t) != excluded_vhs.end())
        continue;

      if (min_dist[vh_f.idx()] < std::numeric_limits<double>::infinity())
      {
        double dist = min_dist[vh_f.idx()] + 1.;
        if (dist < min_dist[vh_t.idx()])
        {
          min_dist[vh_t.idx()] = dist;
          previous[vh_t.idx()] = vh_f;

          if (_target_vhs.find(vh_t) != _target_vhs.end())
            vh_e = vh_t;
        }
      }
    }
  }

  std::vector<EH> path;
  if (vh_e.is_valid())
  {
    //get the path
    std::vector<VH> rvpath;
    VH vh_it = vh_e;
    while (vh_it.is_valid())
    {
      rvpath.push_back(vh_it);
      vh_it = previous[vh_it.idx()];
    }

    if (!rvpath.empty())
    {
      path.reserve(rvpath.size() - 1);
      for (int i = rvpath.size() - 1; i > 0; --i)
        path.push_back(mesh_.edge_handle(mesh_.find_halfedge(rvpath[i], rvpath[i - 1])));
    }
    else
      std::cout << "Error: could not find edge path!" << std::endl;
  }

  return path;
}


template<class MeshT>
std::set<HFH>
FixZipperNodeT<MeshT>::
bounded_surface(const std::set<CH> &_cells, const std::set<EH> &_edges, const HEH _sg_heh) const
{
  std::set<HFH> surface;

  std::set<HFH> all_bound_hfs;
  std::set<HFH> all_hfs;

  for (const auto ch: _cells)
    for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
    {
      all_hfs.insert(mesh_.opposite_halfface_handle(*chf_it));
    }

  //bound halffaces (no boundary halfface)
  for (const auto hf: all_hfs)
  {
    if (_cells.find(mesh_.incident_cell(hf)) == _cells.end() && !mesh_.is_boundary(hf))
      all_bound_hfs.insert(hf);
  }

  return bounded_surface(all_bound_hfs, _edges, _sg_heh);
}


template<class MeshT>
std::set<HFH>
FixZipperNodeT<MeshT>::
bounded_surface(const std::set<HFH> &_bound_hfhs, const std::set<EH> &_excl_edges, const HEH _sg_heh) const
{
  std::set<HFH> surface;


  if (_excl_edges.empty())
  {
    std::set<EH> sf_ehs;
    for (const auto hfh: _bound_hfhs)
    {
      for (auto hfe_it = mesh_.hfe_iter(hfh); hfe_it.valid(); ++hfe_it)
      {
        sf_ehs.insert(*hfe_it);
      }
    }
    //check if bd hfhs has interior singular edges
    for (const auto ehi: sf_ehs)
    {
      if (!mesh_.is_boundary(ehi) && valence_[ehi] != 0)
      {
        std::cout << "Warning: interior singular edges exist on surface. Skip..." << std::endl;
        return std::set<HFH>{};
      }
    }

    return _bound_hfhs;
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "bounded ehs : ";
                             for (const auto e: _excl_edges)
                               std::cout << " " << e;
                             std::cout << std::endl;)

  //start halfface
  HFH hfh;
  for (auto hehf_it = mesh_.hehf_iter(_sg_heh); hehf_it.valid(); ++hehf_it)
    if (_bound_hfhs.find(*hehf_it) != _bound_hfhs.end())
    {
      hfh = *hehf_it;
      break;
    }
  surface.insert(hfh);

  std::queue<HEH> que;
  for (auto hfhe_it = mesh_.hfhe_iter(hfh); hfhe_it.valid(); ++hfhe_it)
    que.push(mesh_.opposite_halfedge_handle(*hfhe_it));

  while (!que.empty())
  {
    auto he_cur = que.front();
    que.pop();

    //touch the boundary edges of the surface
    if (mesh_.is_boundary(he_cur) || (_excl_edges.find(mesh_.edge_handle(he_cur)) != _excl_edges.end()))
      continue;

    HFH hf_uk;
    int n_cut_fhs = n_unvisited_halffaces_at_manifold(_bound_hfhs, surface, he_cur, hf_uk);

    if (n_cut_fhs == 1)
    {
      surface.insert(hf_uk);

      for (auto hfhe_it = mesh_.hfhe_iter(hf_uk); hfhe_it.valid(); ++hfhe_it)
      {
        que.push(mesh_.opposite_halfedge_handle(*hfhe_it));
      }
    }
  }

  //check if it is manifold
  std::set<EH> sf_ehs;
  for (const auto hfh: surface)
  {
    for (auto hfe_it = mesh_.hfe_iter(hfh); hfe_it.valid(); ++hfe_it)
    {
      sf_ehs.insert(*hfe_it);
    }
  }

  for (const auto e: sf_ehs)
  {
    if (_excl_edges.find(e) == _excl_edges.end())
    {
      int n = 0;
      for (auto hehf_it = mesh_.hehf_iter(mesh_.halfedge_handle(e, 0)); hehf_it.valid(); ++hehf_it)
      {
        if (surface.find(*hehf_it) != surface.end())
          n++;
      }
      for (auto hehf_it = mesh_.hehf_iter(mesh_.halfedge_handle(e, 1)); hehf_it.valid(); ++hehf_it)
      {
        if (surface.find(*hehf_it) != surface.end())
          n++;
      }

      if (n != 2)
      {
        std::cout << "Warning: interior edge " << e << " has " << n << " faces on surface" << std::endl;
        return std::set<HFH>{};
      }
    }
    else if (_excl_edges.find(e) != _excl_edges.end())
    {
      int n = 0;
      for (auto hehf_it = mesh_.hehf_iter(mesh_.halfedge_handle(e, 0)); hehf_it.valid(); ++hehf_it)
      {
        if (surface.find(*hehf_it) != surface.end())
          n++;
      }
      for (auto hehf_it = mesh_.hehf_iter(mesh_.halfedge_handle(e, 1)); hehf_it.valid(); ++hehf_it)
      {
        if (surface.find(*hehf_it) != surface.end())
          n++;
      }

      if (n != 1)
      {
        std::cout << "Warning: boundary edge " << e << " has " << n << " faces on surface" << std::endl;
        return std::set<HFH>{};

      }
    }
  }


  return surface;
}

template<class MeshT>
int
FixZipperNodeT<MeshT>::
n_unvisited_halffaces_at_manifold(const std::set<HFH> &_all_hfhs, const std::set<HFH> &_cut_hfhs, const HEH _heh,
                                  HFH &_hfh) const
{
  int n_all = 0, n_unvisited = 0;
  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    auto fh = mesh_.face_handle(*hehf_it);
    if (_all_hfhs.find(*hehf_it) != _all_hfhs.end())
    {
      n_all++;

      //unvisited
      if (_cut_hfhs.count(*hehf_it) == 0)
      {
        _hfh = *hehf_it;
        n_unvisited++;
      }
    }
  }

  if (n_all != 1)
    n_unvisited = 0;

  return n_unvisited;
}

template<class MeshT>
int
FixZipperNodeT<MeshT>::n_bad_vertices_in_cell(const CH _ch, const std::set<VH> &_excl_vhs) const
{
  int n = 0;
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
  {
    if ((sgl_vt_[*cv_it] != 0 || feature_edge_vertex_[*cv_it]) && _excl_vhs.find(*cv_it) == _excl_vhs.end())
      n++;
  }
  return n;
}

template<class MeshT>
int
FixZipperNodeT<MeshT>::n_bad_vertices_in_face(const FH _fh, const std::set<VH> &_excl_vhs) const
{
  int n = 0;
  for (auto fv_it = mesh_.fv_iter(_fh); fv_it.valid(); ++fv_it)
    if ((sgl_vt_[*fv_it] != 0 || feature_edge_vertex_[*fv_it]) && _excl_vhs.find(*fv_it) == _excl_vhs.end())
      n++;


  return n;
}

template<class MeshT>
bool
FixZipperNodeT<MeshT>::is_singular_circle(const HEH _heh) const
{
  if (valence_[mesh_.edge_handle(_heh)] != 0)
  {
    auto sg_hehs = get_halfedges_on_singular_arc(mesh_, valence_, _heh);

    //it's a circle
    if (mesh_.halfedge(sg_hehs[0]).from_vertex() == mesh_.halfedge(sg_hehs.back()).to_vertex())
      return true;
  }

  return false;
}

template<class MeshT>
bool FixZipperNodeT<MeshT>::
is_feature_cell(const CH _ch) const
{
  for (auto cf_it = mesh_.cf_iter(_ch); cf_it.valid(); ++cf_it)
  {
    if (feature_fprop_[*cf_it] > 0)
    {
      return true;
    }
  }

  return false;
}

template<class MeshT>
int
FixZipperNodeT<MeshT>::n_path_edges_in_face(const FH _fh, const std::set<EH> &_epath, EH &_not_path_eh) const
{
  int n_pe = 0;
  for (auto fe_it = mesh_.fe_iter(_fh); fe_it.valid(); ++fe_it)
    if (_epath.find(*fe_it) != _epath.end())
      n_pe++;
    else
      _not_path_eh = *fe_it;

  return n_pe;
}

template<class MeshT>
int
FixZipperNodeT<MeshT>::n_path_vertices_in_face(const FH _fh, const std::set<VH> &_epath_vhs) const
{
  int n_pv = 0;
  for (auto fv_it = mesh_.fv_iter(_fh); fv_it.valid(); ++fv_it)
    if (_epath_vhs.find(*fv_it) != _epath_vhs.end())
      n_pv++;

  return n_pv;
}
}

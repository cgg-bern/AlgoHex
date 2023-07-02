/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define FIXCONSTRAINEDZIPPERNODET_C

#include "FixConstrainedZipperNodeT.hh"
#include <tuple>
#include "CommonFuncs.hh"


namespace AlgoHex
{
//First, detach whatever possible
//Then, fix constrained zipper node with feature arc as guidance
//(At interior feature face vertex, detach sg arc on feature edge before search a surface strip)
//Search for a surface strip if necessary
template<class MeshT>
int
FixConstrainedZipperNodeT<MeshT>::fix_constrained_zipper_node(const VH _vh)
{
  //update valence
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
  {
    if (valence_[*ve_it] != 0)
    {
      sge_.compute_edge_valence(*ve_it);
    }
  }
  auto non_ff_sges = get_non_feature_singular_halfedges(mesh_, feature_face_edge_, valence_, _vh);
  int n_non_ffsge = non_ff_sges.size();
  while (n_non_ffsge >= 2)
  {
    bool dsa = fis_.detach_singular_arc_from_singular_node(_vh, false, true);
    if (!dsa)
      fis_.detach_singular_arcs_of_zipper_node(_vh, false, true);

    auto non_ff_sges_after = get_non_feature_singular_halfedges(mesh_, feature_face_edge_, valence_, _vh);
    int n_non_ffsge_after = non_ff_sges_after.size();

    std::swap(non_ff_sges, non_ff_sges_after);

    if (n_non_ffsge_after == n_non_ffsge)
      break;

    n_non_ffsge = n_non_ffsge_after;
  }

  int n_fixed = 0;

  std::set<HEH> heh_set;
  for (const auto hei: non_ff_sges)
    heh_set.insert(hei);
//        int i=0;
  bool is_bdy_v = mesh_.is_boundary(_vh);
  while (!heh_set.empty())
  {
    auto he_cur = *heh_set.begin();
    heh_set.erase(heh_set.begin());

//            ALGOHEX_DEBUG_ONLY(std::cerr<<"###Fixing ctp at sg edge "<<mesh_.edge_handle(he_cur)<<" at vh "<<mesh_.halfedge(he_cur).from_vertex()<<std::endl;)
//            std::cerr<<"###Fixing ctp at sg edge (with guiding arcs) "<<mesh_.edge_handle(he_cur)<<" at vh "<<mesh_.halfedge(he_cur).from_vertex()<<std::endl;

    int type = fixable_constrained_zipper_node_type(he_cur);
//            ALGOHEX_DEBUG_ONLY(std::cerr <<" type "<<type<<std::endl;)
//            std::cerr <<" type "<<type<<std::endl;

    if (type > 0)
    {
//                ALGOHEX_DEBUG_ONLY(std::cerr<<" fixable ctp type "<<type<<std::endl;)
//                std::cerr<<" fixable ctp type "<<type<<std::endl;

      bool push_to_interior = false;//for singular arc that is long, push to the interior
      if (is_bdy_v)
      {
        auto sg_hehs = get_halfedges_on_singular_arc(mesh_, valence_, he_cur);
        if (type > 2 && type < 6 && sg_hehs.size() > 20)
        {
          push_to_interior = true;
          break;
        }
      }

      if (!push_to_interior)
      {
        bool fixed = fix_constrained_zipper_node1(he_cur, VH(-1), type, true);
        if (fixed)
        {
          auto non_ff_sges2 = get_non_feature_singular_halfedges(mesh_, feature_face_edge_, valence_, _vh);
          for (const auto hei: non_ff_sges2)
            heh_set.insert(hei);

          n_fixed++;
//                        i++;
        }
      }
      else
      {
        ALGOHEX_DEBUG_ONLY(
                std::cerr << "Long singular arc at boundary singular vertex, push to interior later." << std::endl;)
      }
    }
    else
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << " not fixable ctp " << std::endl;)
    }
  }


  //search face strip
  auto non_ff_sges_after = get_non_feature_singular_halfedges(mesh_, feature_face_edge_, valence_, _vh);
  for (const auto hei: non_ff_sges_after)
    heh_set.insert(hei);

  while (!heh_set.empty())
  {
    auto he_cur = *heh_set.begin();
    heh_set.erase(heh_set.begin());

//            ALGOHEX_DEBUG_ONLY(std::cerr << " fixing ctp at sg edge " << mesh_.edge_handle(he_cur);)
//            std::cerr<<"###Fixing ctp at sg edge (face strip) "<<mesh_.edge_handle(he_cur)<<" at vh "<<mesh_.halfedge(he_cur).from_vertex()<<std::endl;
    int type = fixable_constrained_zipper_node_type(he_cur);
//            ALGOHEX_DEBUG_ONLY(std::cerr << " type " << type << std::endl;)
//            std::cerr<<" fixable ctp type "<<type<<std::endl;

    if (type > 0)
    {
//                ALGOHEX_DEBUG_ONLY(std::cerr << " fixable ctp type (no guiding special edge) " << type << std::endl;)
//                std::cerr << " fixable ctp type (no guiding special edge) " << type << std::endl;
      bool fixed = fix_constrained_zipper_node2(he_cur, type, true);
      if (fixed)
      {
        auto non_ff_sges2 = get_non_feature_singular_halfedges(mesh_, feature_face_edge_, valence_, _vh);
        for (const auto &hei: non_ff_sges2)
          heh_set.insert(hei);

        n_fixed++;
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "Fixed " << n_fixed << " ctp at vertex " << _vh << std::endl;)

  return n_fixed;
}

template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::fix_constrained_zipper_node1(const HEH _heh, const VH _vh_target, const int _fix_type,
                                                               const bool _keep_on_ff)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "\nFix constrained zipper node at edge (searching special edge as guide)"
                               << mesh_.edge_handle(_heh) << " at vh " << mesh_.halfedge(_heh).from_vertex() << " val "
                               << valence_[mesh_.edge_handle(_heh)] << " fix type " << _fix_type << std::endl;)

  if (!_heh.is_valid())
  {
    std::cerr << "Error: sge not valid " << std::endl;
    return false;
  }

  ALGOHEX_DEBUG_ONLY(if (valence_[mesh_.edge_handle(_heh)] == 0)
                     {
                       std::cerr << "Error: regular edge" << std::endl;
                       return false;
                     })

  if (_fix_type == 1 || _fix_type == 2)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Extend constrained tp to the other side. " << std::endl;)
    VH vh_next = extend_with_zipper_node(_heh);
    if (!vh_next.is_valid())
      return false;

    return true;
  }
  else if (_fix_type == 3 || _fix_type == 4 || _fix_type == 6)
  {
    std::vector<GDINFO> gds;
    int end_type = 0;
    HEH end_ft_heh(-1);
    HFH end_ft_hfh(-1);
    auto eh_s = mesh_.edge_handle(_heh);

    if (valence_[eh_s] == 1)
    {
      if (_fix_type == 3 || _fix_type == 4)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "Unzip in singular edge direction"
                                     << std::endl;)
        gds = search_for_guides_along_feature_arc(_heh, VH(-1), end_type, end_ft_heh, end_ft_hfh, true,
                                                  true,
                                                  false);

        ALGOHEX_DEBUG_ONLY(std::cerr << "Final gd: "
                                     << std::endl;
                                   for (const auto &gd: gds)
                                   {
                                     std::cerr << mesh_.edge_handle(gd.heh) << " ";
                                   })

        //project to singular arc
        if (second_stage_ && gds.empty())
        {
          VH vhf = mesh_.halfedge(_heh).from_vertex();
          if (feature_node_[vhf] &&
              n_incident_special_edges_kept_on_feature_face(mesh_, feature_face_edge_, feature_edge_, valence_,
                                                            keep_on_feature_face_, vhf) >= 6)
          {
            ALGOHEX_DEBUG_ONLY(std::cerr << "Try unzipping to singular arc on feature face... " << std::endl;)

            gds = search_for_guides_along_singular_arc(_heh, end_type);
            ALGOHEX_DEBUG_ONLY(std::cerr << "Final gd: "
                                         << std::endl;
                                       for (const auto &gd: gds)
                                       {
                                         std::cerr << mesh_.edge_handle(gd.heh) << " ";
                                       })
          }
        }
      }
    }

    if (second_stage_ && _fix_type == 6)
    {
      VH vhf = mesh_.halfedge(_heh).from_vertex();
      if (feature_node_[vhf] &&
          n_incident_special_edges_kept_on_feature_face(mesh_, feature_face_edge_, feature_edge_, valence_,
                                                        keep_on_feature_face_, vhf) >= 6)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "Try unzipping interior singular arc... " << std::endl;)

        gds = search_for_guides_along_singular_arc(_heh, end_type);
        ALGOHEX_DEBUG_ONLY(std::cerr << "Final gd: " << std::endl;
                                   for (const auto &gd: gds)
                                   {
                                     std::cerr << mesh_.edge_handle(gd.heh) << " ";
                                   })
      }
    }

    //split the regular arc to one val+1 and one val-1
    if (!gds.empty())
    {
      if (_keep_on_ff)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "keep sge on face ";)
        for (const auto gdi: gds)
        {
          auto ehi = mesh_.edge_handle(gdi.heh);
          if (feature_face_edge_[ehi])
            keep_on_feature_face_[ehi] = true;
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << ehi;)
        }
      }
      VH vh_hint(-1);
      VH vh_next = unzip_arc_with_guiding_info(_heh, gds, end_type, end_ft_heh, end_ft_hfh, vh_hint);

      if (_fix_type == 6)
      {
        if (vh_next.is_valid() && vh_hint.is_valid())
        {
          HEH heh0, heh1;
          for (auto voh_it = mesh_.voh_iter(vh_next); voh_it.valid(); ++voh_it)
          {
            if (valence_[mesh_.edge_handle(*voh_it)] != 0)
            {
              heh0 = *voh_it;
              break;
            }
          }

          for (auto voh_it = mesh_.voh_iter(vh_hint); voh_it.valid(); ++voh_it)
          {
            if (valence_[mesh_.edge_handle(*voh_it)] != 0)
            {
              heh1 = *voh_it;
              break;
            }
          }
          ALGOHEX_DEBUG_ONLY(std::cerr << "update sg pair id, heh0 " << heh0 << " heh1 " << heh1 << std::endl;)
          if (heh0.is_valid() && heh1.is_valid())
          {
            int max_pair_id = *std::max_element(sg_edge_pairs_.begin(), sg_edge_pairs_.end()) + 1;
            auto sgehs0 = get_edges_on_singular_arc(mesh_, valence_, heh0);
            for (const auto &ehi: sgehs0)
              sg_edge_pairs_[ehi] = max_pair_id;

            max_pair_id++;
            auto sgehs1 = get_edges_on_singular_arc(mesh_, valence_, heh1);
            for (const auto &ehi: sgehs1)
              sg_edge_pairs_[ehi] = max_pair_id;

            ALGOHEX_DEBUG_ONLY(std::cerr << "successfully updated sg pair id." << std::endl;)
          }
        }
      }


      if (vh_next.is_valid())
        return true;
    }
  }

  return false;
}

template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::fix_constrained_zipper_node2(const HEH _heh, const int _fix_type,
                                                               const bool _keep_on_ff)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "\nFix constrained zipper node at edge (searching face strip on feature surface)"
                               << mesh_.edge_handle(_heh) << " at vh " << mesh_.halfedge(_heh).from_vertex() << " val "
                               << valence_[mesh_.edge_handle(_heh)] << " fix type " << _fix_type << std::endl;)

  if (!_heh.is_valid())
  {
    std::cerr << "Error: sge not valid " << std::endl;
    return false;
  }
  if (mesh_.is_deleted(_heh))
  {
    std::cerr << "Error: sge is deleted " << std::endl;
    return false;
  }
  if (valence_[mesh_.edge_handle(_heh)] == 0)
  {
    std::cerr << "Error: regular edge" << std::endl;
    return false;
  }

  if (_fix_type == 3 || _fix_type == 4)
  {
    std::vector<GDINFO> gds;
    int end_type = 0;
    HEH end_ft_heh(-1);
    HFH end_ft_hfh(-1);
    auto eh_s = mesh_.edge_handle(_heh);

    if (valence_[eh_s] == 1)
    {
      //face strip
      if (gds.empty() && _fix_type != 6)
      {
        gds = search_for_guides_on_feature_surface(_heh, end_type, end_ft_heh, end_ft_hfh, true, true);
      }
    }

    //split the regular arc to 1 val+1 and 1 val-1
    if (!gds.empty())
    {
      if (_keep_on_ff)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "keep sge on face ";)
        for (const auto gdi: gds)
        {
          auto ehi = mesh_.edge_handle(gdi.heh);
          if (feature_face_edge_[ehi])
            keep_on_feature_face_[ehi] = true;
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << ehi;)
        }
      }
      VH vh_hint(-1);
      VH vh_next = unzip_arc_with_guiding_info(_heh, gds, end_type, end_ft_heh, end_ft_hfh, vh_hint);

      if (vh_next.is_valid())
        return true;
    }
  }

  return false;
}

//0 - not fixable
//1 - singular arc of val+1, convexly orthogonal but sector angle >= 2(__\|); val-1, convexly orthogonal
//2 - val+1, concavely orthogonal(__\ )
//                               (   |)
//3 - val+1, convexly orthogonal, sector angle < 2, reverse the unzipping and project to feature face or feature(singular) edge
//4 - parallel to all sectors
//5 - val-1, concavely_orthogonal, unzip to boundary;
//6 - val+1, not parallel to any face, but has val>0 singular edge in the same direction. Reverse the unzipping and project to the singular edge.
template<class MeshT>
int
FixConstrainedZipperNodeT<MeshT>::fixable_constrained_zipper_node_type(const VH _vh, HEH &_heh)
{
  _heh = HEH(-1);

  if (sgl_vt_[_vh] == 0 || !feature_face_vertex_[_vh])
    return 0;

  //update valence and sgv prop
  std::set<VH> vhs_set{_vh};
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
  {
    sge_.compute_edge_valence(*ve_it);
    vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
    vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
  }
  for (const auto &vh: vhs_set)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vh);

  std::vector<HEH> sg_hehs;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto ehi = mesh_.edge_handle(*voh_it);
    if (valence_[ehi] == 1 && !feature_face_edge_[ehi]) //
      sg_hehs.push_back(*voh_it);
  }
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto ehi = mesh_.edge_handle(*voh_it);
    if (valence_[ehi] == -1 && !feature_face_edge_[ehi]) //
      sg_hehs.push_back(*voh_it);
  }

  auto or_chs = get_onering_cells(mesh_, _vh);
  for (const auto hehi: sg_hehs)
  {
    auto ehi = mesh_.edge_handle(hehi);

    //can try both search directions.
    if (fac_.is_singular_edge_tangent_to_all_surface(hehi))
    {
      _heh = hehi;
      return 4;
    }
    else
    {
      bool is_tan_af = fac_.is_singular_edge_tangent_to_any_surface(or_chs, hehi);
      ALGOHEX_DEBUG_ONLY(std::cerr << " is tangent to any face at vh?" << is_tan_af << std::endl;)
      if (is_tan_af)
      {
        //direction may not be correct if other sge exist
        std::set<HFH> ft_hfhs;
        auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, *mesh_.hec_iter(hehi),
                                                 ft_hfhs);
        if (n_incident_non_ffe_singular_edges_in_same_region(mesh_, feature_face_edge_, valence_, _vh,
                                                             por_chs) > 1)
        {
          if (fis_.is_locally_meshable(_vh))
          {
            return 0;
          }
          else
            or_chs = get_onering_cells(mesh_, _vh);
        }

        //
        HFH hf_s = *mesh_.hehf_iter(hehi);
        CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_s));

        int e_rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                           valence_, hehi, hf_s);
        if (valence_[ehi] == 1)
          e_rt_axis = negate(e_rt_axis);
        bool has_cvx_prl = has_convexly_parallel_sector_in_direction(_vh, ch_s, e_rt_axis, false);
        if (!has_cvx_prl)
        {
          ALGOHEX_DEBUG_ONLY(
                  std::cerr << "sector parallel in the opposite dir of edge " << ehi << ". val " << valence_[ehi]
                            << "valid case" << std::endl;)
          continue;
        }
        ALGOHEX_DEBUG_ONLY(std::cerr << " edge " << ehi << " val " << valence_[ehi] << std::endl;)

        //find sectors that orthogonal to edge dir
        std::vector<std::vector<HEH>> v_cvx_orth_sec_hehs, v_ccv_orth_sec_hehs;
        std::vector<std::set<FH>> v_cvx_orth_sec_fhs, v_ccv_orth_sec_fhs;
        get_orthogonal_sectors(hehi, v_cvx_orth_sec_hehs, v_cvx_orth_sec_fhs, v_ccv_orth_sec_hehs,
                               v_ccv_orth_sec_fhs);

        if (valence_[ehi] == 1)
        {
          if (!v_ccv_orth_sec_hehs.empty())
          {//transversal unzipping, creating an isolated zipper node
            _heh = hehi;
            return 2;
          }
          else
          {
            if (v_cvx_orth_sec_hehs.empty())
            {
              std::cerr << "Warning: has orthogonal sector at vertex " << _vh << ", but not found."
                        << n_incident_singular_edges(mesh_, valence_, _vh)
                        << " incident singular edges." << std::endl;
              continue;
            }

            if (!mesh_.is_boundary(_vh))
            {
              for (auto j = 0u; j < v_cvx_orth_sec_hehs.size(); ++j)
              {
                int angle = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_,
                                                                                               cell_quaternions_,
                                                                                               trans_prop_, valence_,
                                                                                               v_cvx_orth_sec_hehs[j],
                                                                                               v_cvx_orth_sec_fhs[j]);
                ALGOHEX_DEBUG_ONLY(std::cerr << " cvx sec agl: " << angle << " ";)

                if (angle >= 2)
                {
                  if (angle == 2)
                  {
                    auto sa = FieldAngleCalculatorT<MeshT>::sector_angle(mesh_, v_cvx_orth_sec_hehs[j]);
                    if (sa < M_PI)
                    {//better hex
                      ALGOHEX_DEBUG_ONLY(std::cerr
                                                 << " convexly orthogonal sec found and there is sec angle ==2, can extend to other side"
                                                 << std::endl;)

                      _heh = hehi;

                      return 1;
                    }
                  }
                  if (angle > 2)
                  {
                    ALGOHEX_DEBUG_ONLY(std::cerr
                                               << " convexly orthogonal sec found and there is sec angle " << angle
                                               << ", can extend to other side"
                                               << std::endl;)

                    _heh = hehi;
                    return 1;
                  }
                }
              }
            }

            _heh = hehi;

            return 3;
          }
        }
        else if (valence_[ehi] == -1)
        {
          if (!v_cvx_orth_sec_hehs.empty())
          {//val -1 arc can be extended to other side
            _heh = hehi;

            return 1;
          }
          else
          { //fix tp
            if (v_ccv_orth_sec_hehs.empty())
            {
              ALGOHEX_DEBUG_ONLY(std::cerr << " Warning: has orthogonal sector, but not found."
                                           << n_incident_singular_edges(mesh_, valence_, _vh)
                                           << " incident singular edges." << std::endl;)
              continue;
            }

            _heh = hehi;
            return 5;
          }
        }
      }
      else if (feature_node_[_vh] && n_incident_feature_edges(mesh_, feature_edge_, _vh) >= 5)
      {
        HFH hf_s = *mesh_.hehf_iter(hehi);
        CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_s));

        int e_rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                           valence_, hehi, hf_s);

        e_rt_axis = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, hehi, ch_s, e_rt_axis);

        HEH he_found(-1);
        HFH hf_found(-1), hf_e(-1);
        CH ch_found(-1);
        int search_ax_found = -1;
        ALGOHEX_DEBUG_ONLY(std::cerr << "eh " << ehi << "valence_ " << valence_[ehi] << std::endl;)

        if(valence_[ehi] > 0) {
          bool has_in_dir = fac_.has_singular_edge_in_direction(_heh, ch_s, e_rt_axis, _vh, ch_found,
                                                                hf_found, he_found, search_ax_found);
          ALGOHEX_DEBUG_ONLY(std::cerr << "found singular edge in the same direction? " << has_in_dir << std::endl;)
          if (has_in_dir)
          {
            _heh = hehi;

            return 6;
          }
        } else if(valence_[ehi] < 0) {//reverse direction
          bool has_in_dir = fac_.has_singular_edge_in_direction(_heh, ch_s, negate(e_rt_axis), _vh, ch_found,
                                                                hf_found, he_found, search_ax_found);
          ALGOHEX_DEBUG_ONLY(std::cerr << "negative valence! found singular edge in the reverse direction? " << has_in_dir << std::endl;)
          if (has_in_dir)
          {
            _heh = hehi;

            return 6;
          }
        }
      }
    }
  }

  return 0;
}


//0 - not fixable
//1 - val+1, convexly orthogonal but sector angle >= 2; val-1, convexly orthogonal
//2 - val+1, concavely orthogonal; val-1, concavely but sector angle >= 2
//3 - val+1, convexly orthogonal, sector angle < 2, project to feature face or feature(singular) edge
//4 - parallel to all sectors
//5 - val-1, concavely_orthogonal, fix tp; val+1, concavely orthogonal, push inside
//6 - val+1, not parallel to any face, but has val>0 singular edge in the same direction. Project to the singular edge.
template<class MeshT>
int
FixConstrainedZipperNodeT<MeshT>::fixable_constrained_zipper_node_type(const HEH _heh)
{
  auto eh = mesh_.edge_handle(_heh);
  if (mesh_.is_deleted(eh))
  {
    std::cerr << "Error: received a deleted edge " << eh << std::endl;
    return 0;
  }

  VH vh_f = mesh_.halfedge(_heh).from_vertex();
  if (sgl_vt_[vh_f] == 0 || !feature_face_vertex_[vh_f])
  {
    std::cerr << "Error: received edge whose from vertex is non-singular" << vh_f << std::endl;

    return 0;
  }

  if (feature_face_edge_[eh])
  {
    std::cerr << "Received feature face edge " << eh << std::endl;
    return 0;
  }

  if (fabs(valence_[eh]) != 1)
  {
    std::cerr << "Regular or complex singular edge " << eh << std::endl;

    return 0;
  }

  //
  sge_.compute_edge_valence(eh);
  MeshPropertiesT<MeshT>::update_singular_vertex_property(mesh_.edge(eh).from_vertex());
  MeshPropertiesT<MeshT>::update_singular_vertex_property(mesh_.edge(eh).to_vertex());

  //can try both search directions. In the end, either no tangentially crossing or end at at feature node, turned to the second casee
  if (fac_.is_singular_edge_tangent_to_all_surface(_heh))
    return 4;
  else
  {
    auto or_chs = get_onering_cells(mesh_, vh_f);
    bool is_tan_af = fac_.is_singular_edge_tangent_to_any_surface(or_chs, _heh);
    ALGOHEX_DEBUG_ONLY(std::cerr << " is tangent to any face " << is_tan_af << std::endl;)
    if (is_tan_af)
    {
      //direction may not be correct if other sge exist
      std::set<HFH> ft_hfhs;
      auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vh_f, or_chs, *mesh_.hec_iter(_heh),
                                               ft_hfhs);
      int n_inc_nff_sge = n_incident_non_ffe_singular_edges_in_same_region(mesh_, feature_face_edge_, valence_, vh_f,
                                                                           por_chs);
      if (n_inc_nff_sge > 1)
      {
        if (fis_.is_locally_meshable(vh_f))
        {
          return 0;
        }
      }

      //
      HFH hf_s = *mesh_.hehf_iter(_heh);
      CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_s));

      int e_rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                         valence_, _heh, hf_s);
      if (valence_[mesh_.edge_handle(_heh)] == 1)
        e_rt_axis = negate(e_rt_axis);
      bool has_cvx_prl = has_convexly_parallel_sector_in_direction(vh_f, ch_s, e_rt_axis, false);
      if (!has_cvx_prl)
      {
        ALGOHEX_DEBUG_ONLY(
                std::cerr << "sector parallel in the opposite dir of edge " << eh << ". valid case" << std::endl;)
        return 0;
      }

      //find sectors that orthogonal to edge dir
      std::vector<std::vector<HEH>> v_cvx_orth_sec_hehs, v_ccv_orth_sec_hehs;
      std::vector<std::set<FH>> v_cvx_orth_sec_fhs, v_ccv_orth_sec_fhs;
      get_orthogonal_sectors(_heh, v_cvx_orth_sec_hehs, v_cvx_orth_sec_fhs, v_ccv_orth_sec_hehs, v_ccv_orth_sec_fhs);

      if (valence_[eh] == 1)
      {
        if (!v_ccv_orth_sec_hehs.empty())
        {
          if (mesh_.is_boundary(vh_f) && n_inc_nff_sge > 1)
          {//prefer pushing to interior than creating another zipper node
            return 5;
          }

          return 2;
        }
        else
        {
          if (v_cvx_orth_sec_hehs.empty())
          {
            std::cerr << " Warning: has orthogonal sector, but not found" << std::endl;
            return 0;
          }

          if (!mesh_.is_boundary(vh_f))
          {
            for (auto j = 0u; j < v_cvx_orth_sec_hehs.size(); ++j)
            {
              int angle = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_,
                                                                                             cell_quaternions_,
                                                                                             trans_prop_, valence_,
                                                                                             v_cvx_orth_sec_hehs[j],
                                                                                             v_cvx_orth_sec_fhs[j]);
              ALGOHEX_DEBUG_ONLY(std::cerr << " cvx sec agl: " << angle << " ";)

              if (angle == 2)
              {
                auto sa = FieldAngleCalculatorT<MeshT>::sector_angle(mesh_, v_cvx_orth_sec_hehs[j]);
                if (sa < M_PI)
                {//better hex
                  ALGOHEX_DEBUG_ONLY(std::cerr
                                             << " convexly ortho sec found and there is sec angle ==2, can extend to other side"
                                             << std::endl;)

                  return 1;
                }
              }
              if (angle > 2)
              {
                ALGOHEX_DEBUG_ONLY(std::cerr
                                           << " convexly ortho sec found and there is sec angle " << angle
                                           << ", can extend to other side"
                                           << std::endl;)

                return 1;
              }
            }
          }

          return 3;
        }
      }
      else if (valence_[eh] == -1)
      {//val -1 arc can be detached from the node if convexly orthogonal
        if (!v_cvx_orth_sec_hehs.empty())
          return 1;
        else
        {
          if (v_ccv_orth_sec_hehs.empty())
          {
            std::cerr << " Warning: has orthogonal sector, but not found" << std::endl;
            return 0;
          }

          return 5;
        }
      }
    }
    else if (feature_node_[vh_f] && n_incident_feature_edges(mesh_, feature_edge_, vh_f) >= 5)
    {
      HFH hf_s = *mesh_.hehf_iter(_heh);
      CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_s));

      int e_rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                         valence_, _heh, hf_s);
      e_rt_axis = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, _heh, ch_s, e_rt_axis);

      HEH he_found(-1);
      HFH hf_found(-1), hf_e(-1);
      CH ch_found(-1);
      int search_ax_found = -1;
      ALGOHEX_DEBUG_ONLY(std::cerr << "eh " << eh << "valence_ " << valence_[eh] << std::endl;)

      if(valence_[eh] > 0) {
        bool has_in_dir = fac_.has_singular_edge_in_direction(_heh, ch_s, e_rt_axis, vh_f, ch_found, hf_found, he_found,
                                                              search_ax_found);
        ALGOHEX_DEBUG_ONLY(std::cerr << "found singular edge in the same direction? " << has_in_dir << std::endl;)
        if (has_in_dir)
          return 6;
      } else if(valence_[eh] < 0) {//reverse direction
        bool has_in_dir = fac_.has_singular_edge_in_direction(_heh, ch_s, negate(e_rt_axis), vh_f, ch_found, hf_found, he_found,
                                                              search_ax_found);
        ALGOHEX_DEBUG_ONLY(std::cerr << "negative valence! found singular edge in the reverse direction? " << has_in_dir << std::endl;)
        if (has_in_dir)
          return 6;
      }
    }
  }

  return 0;
}

template<class MeshT>
int
FixConstrainedZipperNodeT<MeshT>::singular_arc_tangentially_touching_feature_face_type(const VH _vh)
{
  if (sgl_vt_[_vh] == 0 || !feature_face_vertex_[_vh] || mesh_.is_boundary(_vh))
    return 0;

  //get singular halfedges
  auto sg_hehs = get_singular_halfedges(mesh_, valence_, _vh);
  if (sg_hehs.size() != 2)
    return 0;

  //if both are feature face edges, skip
  auto eh0 = mesh_.edge_handle(sg_hehs[0]);
  auto eh1 = mesh_.edge_handle(sg_hehs[1]);
  auto is_eh0_ffe = feature_face_edge_[eh0];
  auto is_eh1_ffe = feature_face_edge_[eh1];
  if (is_eh0_ffe && is_eh1_ffe)
    return 0;

  //if kept on feature face
  if (feature_edge_vertex_[_vh])
  {
    if ((is_eh0_ffe || is_eh1_ffe) && (keep_on_feature_face_[eh0] || keep_on_feature_face_[eh1]))
    {
      bool detach = false;
      for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
      {
        if (feature_face_edge_[*ve_it] && keep_on_feature_face_[*ve_it] && valence_[*ve_it] == 0)
        {
          detach = true;
          break;
        }
      }

      if (!detach)
        return 0;
    }
//            if ((is_eh0_ffe && keep_on_feature_face_[eh0]) ||
//                (feature_face_edge_[eh1] && keep_on_feature_face_[eh1]))
//                return 0;
  }

  //if cross
  auto or_chs = get_onering_cells(mesh_, _vh);
  bool is_crs = is_singular_arc_crossing_feature_face(_vh);
  if (is_crs)
    return 0;

  ALGOHEX_DEBUG_ONLY(std::cerr << "vh " << _vh << " eho ff " << is_eh0_ffe << " eh1 ff " << is_eh1_ffe;)
  //if at feature edge vertex, one singular edge is on feature face and there is a sector orthogonal to the sge direction
  if (feature_edge_vertex_[_vh])
  {
    if ((is_eh0_ffe && !is_eh1_ffe) || (!is_eh0_ffe && is_eh1_ffe))
    {
      if (valence_[eh0] != 4 && valence_[eh1] != 4)
      {
        HEH nffheh(-1);
        if (!is_eh0_ffe)
          nffheh = sg_hehs[0];
        else if (!is_eh1_ffe)
          nffheh = sg_hehs[1];

        bool is_oth = fac_.is_singular_edge_orthogonal_to_any_surface(or_chs, nffheh);
        ALGOHEX_DEBUG_ONLY(std::cerr << " is orth " << is_oth;)
        if (is_oth)
        {
          int idx = MeshPropertiesT<MeshT>::node_index(_vh);
          ALGOHEX_DEBUG_ONLY(std::cerr << " nid " << idx << "   ";)
          if (idx == 10)
            return 2;

          if (fis_.is_locally_meshable(_vh))
            return 2;
        }
      }
    }
  }


  return 1;
}

template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::is_singular_arc_crossing_feature_face(const VH _vh) const
{
  //get singular halfedges
  auto sg_hehs = get_singular_halfedges(mesh_, valence_, _vh);

  if (sg_hehs.size() != 2)
    return false;

  EH other_sge = mesh_.edge_handle(sg_hehs[1]);
  if (feature_face_edge_[mesh_.edge_handle(sg_hehs[0])] || feature_face_edge_[other_sge])
    return false;

  auto or_chs = get_onering_cells(mesh_, _vh);
  std::set<HFH> bound_ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, *mesh_.hec_iter(sg_hehs[0]),
                                           bound_ft_hfhs);

  for (const auto chi: por_chs)
  {
    for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
    {
      if (*ce_it == other_sge)
        return false;
    }
  }

  return true;
}

template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::is_singularity_crossing_feature_face(const VH _vh) const
{
  //get singular halfedges
  auto sg_hehs = get_singular_halfedges(mesh_, valence_, _vh);

  if (sg_hehs.size() < 2)
    return false;

  auto or_chs = get_onering_cells(mesh_, _vh);
  std::set<HFH> bound_ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, *mesh_.hec_iter(sg_hehs[0]),
                                           bound_ft_hfhs);

  std::set<EH> sgehs;
  for (const auto chi: por_chs)
  {
    for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
    {
      if (valence_[*ce_it] != 0)
        sgehs.insert(*ce_it);
    }
  }

  for (auto i = 1u; i < sg_hehs.size(); ++i)
    if (sgehs.find(mesh_.edge_handle(sg_hehs[i])) == sgehs.end())
      return true;

  return false;
}

template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::
is_separable_at_vertex(const VH _vh, std::set<HEH> &_ft_hehs)
{
  if (sgl_vt_[_vh] == 0 || !feature_edge_vertex_[_vh] || mesh_.is_boundary(_vh))
    return false;

  //get singular halfedges
  auto sg_hehs = get_singular_halfedges(mesh_, valence_, _vh);

  if (sg_hehs.size() < 2)
    return false;
  else if (sg_hehs.size() == 2)
  {
    bool is_crs = is_singular_arc_crossing_feature_face(_vh);
    if (!is_crs)
      return false;
  }
  else if (sg_hehs.size() > 2)
  {
    bool crs = is_singularity_crossing_feature_face(_vh);
    if (!crs)
      return false;
  }

  //feature face sector angle
  std::vector<HEH> v_hes;
//        std::vector<HFH> v_hfs;
  std::vector<int> v_scagls;
//        std::vector<int> v_das;

  for (auto voh_it2 = mesh_.voh_iter(_vh); voh_it2.valid(); ++voh_it2)
  {
    auto ve = mesh_.edge_handle(*voh_it2);
    if (feature_edge_[ve] && valence_[ve] == 0)
    {
      std::vector<HFH> ft_hfhs;
      for (auto hehf_it = mesh_.hehf_iter(*voh_it2); hehf_it.valid(); ++hehf_it)
      {
        if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
        {
          ft_hfhs.push_back(*hehf_it);
        }
      }
      int fsec_agl = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_, cell_quaternions_,
                                                                                         trans_prop_, *voh_it2,
                                                                                         ft_hfhs[0],
                                                                                         ft_hfhs[1]);


      if (fsec_agl == 2 || fsec_agl == 1)
      {
        v_hes.push_back(*voh_it2);
        v_scagls.push_back(fsec_agl);
      }
      else if (fsec_agl == 3)
      {
        v_hes.push_back(*voh_it2);
        v_scagls.push_back(4 - fsec_agl);
      }
    }
  }

  //if equal
  if (!feature_node_[_vh] && v_hes.size() == 2)
  {
    if (v_scagls[0] == v_scagls[1])
      return false;
  }

  auto or_chs = get_onering_cells(mesh_, _vh);

  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    EH ehi = mesh_.edge_handle(*voh_it);
    if (fabs(valence_[ehi]) == 1 && !feature_face_edge_[ehi])
    {
      //find the start guide
      EH eh_ss = mesh_.edge_handle(*voh_it);
      HFH hf_ss = *mesh_.hehf_iter(*voh_it);
      CH ch_ss = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_ss));
      ALGOHEX_DEBUG_ONLY(
              std::cerr << "At vh " << _vh << " Searching on sg feature halfedge with defect angle wrt sge " << ehi
                        << std::endl;)
      int sg_rt_ax = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                        valence_, *voh_it, hf_ss);


      std::set<HFH> ft_hfhs;
      auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, ch_ss, ft_hfhs);

      //find from the set of feature edges candidate feature edges
      for (auto j = 0u; j < v_hes.size(); ++j)
      {
        HFH ft_hf_start, ft_hf_end;
        for (auto hehf_it = mesh_.hehf_iter(v_hes[j]); hehf_it.valid(); ++hehf_it)
        {
          if (ft_hfhs.find(*hehf_it) != ft_hfhs.end())
          {
            ft_hf_start = *hehf_it;
          }
          else if (ft_hfhs.find(mesh_.opposite_halfface_handle(*hehf_it)) != ft_hfhs.end())
          {
            ft_hf_end = *hehf_it;
          }
        }

        CH ch_i = mesh_.incident_cell(ft_hf_start);
        auto dir_i = mesh_.vertex(mesh_.halfedge(v_hes[j]).to_vertex()) -
                     mesh_.vertex(mesh_.halfedge(v_hes[j]).from_vertex());
        int ax_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_i],
                                                      dir_i.normalize()).second;
//                    int ax_i_in_ss = fac_.axis_in_chs_expressed_in_cht(ch_i, ch_ss, ax_i, or_chs);
        int rt_ax_in_i = fac_.axis_in_chs_expressed_in_cht(ch_ss, ch_i, sg_rt_ax, por_chs);

        ALGOHEX_DEBUG_ONLY(
                std::cerr << "test eh: " << mesh_.edge_handle(v_hes[j]) << " ch_i " << ch_i << " ax_i " << ax_i
                          << " ax_i_in_ss: "
                          << rt_ax_in_i << " ch_ss " << ch_ss << " rt ax " << sg_rt_ax << " v_scagls " << v_scagls[j]
                          << std::endl;)

        if (rt_ax_in_i / 2 == ax_i / 2)
        {
          if (v_scagls[j] == 1)
          {
            //check if cancelling the sge will result in a valid feature sector
            int trans_i = rt_ax_in_i + 4;
            HFH hf_change = mesh_.opposite_halfface_handle(
                    mesh_.adjacent_halfface_in_cell(ft_hf_start, v_hes[j]));
            int new_trans = tq_.mult_transitions_idx(trans_prop_[hf_change], trans_i);
            int old_trans = trans_prop_[hf_change];

            //apply change
            trans_prop_[hf_change] = new_trans;
            trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(
                    new_trans);

            int fsec_after = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                                 cell_quaternions_,
                                                                                                 trans_prop_, v_hes[j],
                                                                                                 ft_hf_start,
                                                                                                 ft_hf_end);

            ALGOHEX_DEBUG_ONLY(std::cerr << "Feature face sector! Sector agl after " << fsec_after << " fhs "
                                         << mesh_.face_handle(ft_hf_start) << " " << mesh_.face_handle(ft_hf_end)
                                         << std::endl;)

            if (fsec_after == 2)
            {
              _ft_hehs.insert(v_hes[j]);
            }

            //restore matching
            trans_prop_[hf_change] = old_trans;
            trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(
                    old_trans);
          }
        }
      }

//                same direction, no need to check the other sge
      if (sg_hehs.size() == 2 && !_ft_hehs.empty())
        return true;
    }
  }

  if (!_ft_hehs.empty())
    return true;

  return false;
}


template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::is_valid_feature_face_sector_after(const HEH _heh, const HFH _hfh_s, const HFH _hfh_e,
                                                                     int _trans)
{
  int fsec_agl_before = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                            cell_quaternions_,
                                                                                            trans_prop_, _heh, _hfh_s,
                                                                                            _hfh_e);


  HFH hf_change = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(_hfh_s, _heh));
  int old_trans = trans_prop_[hf_change];

  int new_trans = tq_.mult_transitions_idx(old_trans, _trans);
  trans_prop_[hf_change] = new_trans;
  trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(new_trans);

  int fsec_agl_after = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                           cell_quaternions_,
                                                                                           trans_prop_, _heh, _hfh_s,
                                                                                           _hfh_e);


  ALGOHEX_DEBUG_ONLY(std::cerr << "ff sghe " << mesh_.halfedge(_heh).from_vertex() << " "
                               << mesh_.halfedge(_heh).to_vertex()
                               << " fs " << mesh_.face_handle(_hfh_s) << " fe " << mesh_.face_handle(_hfh_e) << " fc "
                               << mesh_.face_handle(hf_change) << " trans " << _trans << " old tran "
                               << old_trans << " new tran " << new_trans
                               << " " << valence_[mesh_.edge_handle(_heh)] << " sec agl before " << fsec_agl_before
                               << " sec agl after " << fsec_agl_after << " dihedral angle "
                               << dihedral_angle_from_halfface_to_halfface(mesh_, _heh, _hfh_s, _hfh_e) * 180. / M_PI
                               << std::endl;)

  //restore matching
  trans_prop_[hf_change] = old_trans;
  trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(old_trans);

  auto da = dihedral_angle_from_halfface_to_halfface(mesh_, _heh, _hfh_s, _hfh_e);
  if (fsec_agl_after > 4)//TODO: if higher valence, it's valid
    return false;

  if (da >= M_PI)
  {
    if (fsec_agl_after < 2)
      return false;
  }
  else if (fsec_agl_after < 1)
    return false;

  return true;
}

template<class MeshT>
std::vector<typename FixConstrainedZipperNodeT<MeshT>::GDINFO>
FixConstrainedZipperNodeT<MeshT>::search_for_guides_along_feature_arc(const HEH _sg_heh, const VH _target_vh,
                                                                      int &_end_type, HEH &_end_ft_he, HFH &_end_ft_hf,
                                                                      const bool _open_sec, const bool _same_dir,
                                                                      const bool _cross)
{
  VH vh_f = mesh_.halfedge(_sg_heh).from_vertex();

  _end_type = -1;
  if (!feature_edge_vertex_[vh_f])
    return std::vector<GDINFO>{};

  //find the start guide
  EH eh_ss = mesh_.edge_handle(_sg_heh);
  HFH hf_ss = *mesh_.hehf_iter(_sg_heh);
  CH ch_ss = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_ss));
  ALGOHEX_DEBUG_ONLY(std::cerr << "Searching guiding halfedges..." << std::endl;
                             std::cerr << "hes " << _sg_heh << " ch_ss " << ch_ss << " hf_ss " << hf_ss << std::endl;)
  int sg_rt_ax = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, _sg_heh, hf_ss);

  //measure with sg_rt_ax

  //search in singular halfedge direction or opposite
  int search_ax = sg_rt_ax;
  if (valence_[mesh_.edge_handle(_sg_heh)] == 1)
    search_ax = negate(sg_rt_ax);
  if (!_same_dir)
    search_ax = negate(search_ax);

  //in the direction of rt axis, this will increase sector angle by +90
  int rt_ax = search_ax;
  if (_open_sec)
    rt_ax = negate(rt_ax);

  bool negated = search_ax == negate(rt_ax);

  //one ring cells
  std::set<CH> or_chs = get_onering_cells(mesh_, vh_f);

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vh_f, or_chs, ch_ss, ft_hfhs);

  //find a start feature halfedge
  HEH he_s(-1);
  HFH hf_s(-1), hf_e(-1);
  CH ch_s(-1);
  int search_ax_s = -1;
  int rt_ax_s = -1;
  for (auto voh_it = mesh_.voh_iter(vh_f); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh_.edge_handle(*voh_it);
    if (feature_edge_[ve])
    {
      //don't create +2 or -2
      if ((valence_[ve] == -1 && !_open_sec) || (valence_[ve] == 1 && _open_sec))
        continue;

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

      CH ch_i = mesh_.incident_cell(hf_start);
      auto dir0 = mesh_.vertex(mesh_.halfedge(*voh_it).to_vertex()) -
                  mesh_.vertex(mesh_.halfedge(*voh_it).from_vertex());
      int ax_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_i],
                                                    dir0.normalize()).second;
      int ax_i_in_ss = fac_.axis_in_chs_expressed_in_cht(ch_i, ch_ss, ax_i, por_chs);

      ALGOHEX_DEBUG_ONLY(std::cerr << "test eh: " << mesh_.edge_handle(*voh_it) << " ch_i" << ch_i << " search_ax_i: "
                                   << ax_i << " axi in s: " << ax_i_in_ss << " ch_s " << ch_s << " search ax s "
                                   << search_ax << std::endl;)

      if (ax_i_in_ss == search_ax)
      {
        int fsec_agl_before = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                                  cell_quaternions_,
                                                                                                  trans_prop_, *voh_it,
                                                                                                  hf_start, hf_end);

        int rt_ax_i = ax_i;
        if (negated)
          rt_ax_i = negate(rt_ax_i);

        int trans = rt_ax_i + 4;

        bool is_valid_sec = is_valid_feature_face_sector_after(*voh_it, hf_start, hf_end, trans);
        ALGOHEX_DEBUG_ONLY(std::cerr << " valid " << is_valid_sec << std::endl;)

        //check the other side of the feature face sector
        if (is_valid_sec && _cross)
        {
          HFH hf_i_opp = mesh_.opposite_halfface_handle(hf_start);
          int ax_i_other = tq_.axis_after_transition(ax_i, trans_prop_[hf_i_opp]);
          if (negated)
            ax_i_other = negate(ax_i_other);

          trans = ax_i_other + 4;
          is_valid_sec = is_valid_feature_face_sector_after(mesh_.opposite_halfedge_handle(*voh_it), hf_i_opp,
                                                            mesh_.opposite_halfface_handle(hf_end), trans);
        }


        if (is_valid_sec)
        {
          ch_s = ch_i;
          hf_s = hf_start;
          hf_e = hf_end;
          he_s = *voh_it;
          search_ax_s = ax_i;
          rt_ax_s = rt_ax_i;

          ALGOHEX_DEBUG_ONLY(
                  std::cerr << "found fe ch: " << mesh_.edge_handle(he_s) << " ch_s" << ch_s << " search_ax_s: "
                            << search_ax << " rt ax: " << rt_ax_s << std::endl;)

          break;
        }
      }
    }
  }

  if (he_s.is_valid())
  {
    std::vector<std::vector<GDINFO>> v_guides;
    int end_cw = -1, end_ccw = -1;

    //normal axis
    int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s],
                                                     mesh_.normal(hf_s)).second;

    ALGOHEX_DEBUG_ONLY(
            std::cerr << "bdy ch: " << ch_s << " he ax: " << search_ax << " nm ax: " << nm_axis << std::endl;)

    std::vector<HFH> inc_hfhs;
    HEH end_ft_heh_cw(-1);
    HFH end_ft_hfh_cw(-1);
    auto guides = get_guides_cw(GDINFO(he_s, hf_s, nm_axis, search_ax_s, rt_ax_s, false), _target_vh,
                                _open_sec, end_cw, end_ft_heh_cw, end_ft_hfh_cw);
    ALGOHEX_DEBUG_ONLY(std::cerr << "cw end type " << end_cw << " guiding edges: ";
                               for (auto gd: guides)
                                 std::cerr << " " << mesh_.edge_handle(gd.heh);
                               std::cerr << std::endl;)

    if (!guides.empty())
    {
      v_guides.push_back(guides);
    }


    //test ccw
    HEH end_ft_heh_ccw(-1);
    HFH end_ft_hfh_ccw(-1);
    HFH hf_e_opp = mesh_.opposite_halfface_handle(hf_e);
    HEH he_s_opp = mesh_.opposite_halfedge_handle(he_s);

    CH ch_i = mesh_.incident_cell(hf_e_opp);
    int ax_in_i = fac_.axis_in_chs_expressed_in_cht(ch_s, ch_i, search_ax_s, por_chs);

    auto dir0 = mesh_.vertex(mesh_.halfedge(he_s_opp).from_vertex()) -
                mesh_.vertex(mesh_.halfedge(he_s_opp).to_vertex());
    int ax_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_i],
                                                  dir0.normalize()).second;
    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_i: " << ch_i << " he ax: " << ax_in_i << " fe ax: " << ax_i << std::endl;
                               std::cerr << " search ax dir " << mesh_.vertex(mesh_.halfedge(he_s_opp).to_vertex())
                                         << " "
                                         << (ovm2eigen(mesh_.vertex(mesh_.halfedge(he_s_opp).to_vertex())) +
                                             AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_i],
                                                                                     (AxisAlignment) ax_i)).transpose();)

    if (ax_i == ax_in_i)
    {
      int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_i],
                                                       mesh_.normal(hf_e_opp)).second;
      int rt_ax_in_i = ax_in_i;
      if (negated)
        rt_ax_in_i = negate(rt_ax_in_i);


      guides = get_guides_ccw(GDINFO(he_s_opp, hf_e_opp, nm_axis, ax_in_i, rt_ax_in_i, true), _target_vh,
                              _open_sec, end_ccw, end_ft_heh_ccw, end_ft_hfh_ccw);
      ALGOHEX_DEBUG_ONLY(std::cerr << "ccw end type " << end_ccw << "guiding edges: ";
                                 for (auto gd: guides)
                                   std::cerr << " " << mesh_.edge_handle(gd.heh);
                                 std::cerr << std::endl;)

      if (!guides.empty())
      {
        if (v_guides.empty() || (end_cw + 1) / 2 == (end_ccw + 1) / 2)
          v_guides.push_back(guides);
        else if ((end_cw + 1) / 2 > (end_ccw + 1) / 2)
          v_guides[0] = guides;
      }
    }


    //take the better
    if (v_guides.size() == 2)
    {
      if (v_guides[0].size() > v_guides[1].size())
      {//take shorter
        v_guides[0].swap(v_guides[1]);
        ALGOHEX_DEBUG_ONLY(std::cerr << " ccw is shorter" << std::endl;)
      }
      else if (v_guides[0].size() == v_guides[1].size())
      {
        VH vh_t = mesh_.halfedge(_sg_heh).to_vertex();
        VH vh_other;
        for (auto voh_it = mesh_.voh_iter(vh_t); voh_it.valid(); ++voh_it)
        {
          auto ve = mesh_.edge_handle(*voh_it);
          if (valence_[ve] != 0 && ve != mesh_.edge_handle(_sg_heh))
          {
            vh_other = mesh_.halfedge(*voh_it).to_vertex();
            break;
          }
        }

        auto v0 = (mesh_.vertex(vh_other) - mesh_.vertex(vh_f)).normalize();

        double max_dist = -100000.;
        int pos = 0;
        for (auto i = 0u; i < v_guides.size(); ++i)
        {
          auto pt_m = (mesh_.vertex(mesh_.halfedge(v_guides[i].front().heh).to_vertex()) +
                       mesh_.vertex(mesh_.halfedge(v_guides[i].front().heh).from_vertex())) / 2.;
          auto pt_bc = mesh_.barycenter(mesh_.face_handle(v_guides[i].front().hfh));
          auto vi = pt_bc - pt_m;
          vi.normalize();

          double dist = vi.dot(v0);
          if (dist > max_dist)
          {
            pos = i;
            max_dist = dist;
          }
        }

        ALGOHEX_DEBUG_ONLY(std::cerr << " cw and ccw same size. closer guides " << pos << std::endl;)
        if (pos == 1)
        {
          v_guides[0].swap(v_guides[1]);
        }
      }
    }

    if (!v_guides.empty())
    {
      if (v_guides[0].back().ccw)
      {
        _end_ft_he = end_ft_heh_ccw;
        _end_ft_hf = end_ft_hfh_ccw;
        _end_type = end_ccw;
      }
      else
      {
        _end_ft_he = end_ft_heh_cw;
        _end_ft_hf = end_ft_hfh_cw;
        _end_type = end_cw;
      }

      if (_end_type == 6)
      {
        VH vh_last = mesh_.halfedge(v_guides[0].back().heh).to_vertex();
        if (v_guides[0].back().ccw)
          vh_last = mesh_.halfedge(v_guides[0].back().heh).from_vertex();
        auto guides2 = search_for_guides_on_feature_surface(v_guides[0].back().heh,
                                                            mesh_.incident_cell(v_guides[0].back().hfh), vh_last,
                                                            v_guides[0].back().search_axis,
                                                            _end_type, _end_ft_he, _end_ft_hf, _open_sec);

        v_guides[0].insert(v_guides[0].end(), guides2.begin(), guides2.end());
      }

      return v_guides[0];
    }
  }

  return std::vector<GDINFO>{};
}

template<class MeshT>
std::vector<typename FixConstrainedZipperNodeT<MeshT>::GDINFO>
FixConstrainedZipperNodeT<MeshT>::search_for_guides_along_singular_arc(const HEH _sg_heh, int &_end_type)
{
  VH vh_f = mesh_.halfedge(_sg_heh).from_vertex();

  _end_type = -1;
  if (!feature_node_[vh_f])
    return std::vector<GDINFO>{};

  //find the start guide
  EH eh_query = mesh_.edge_handle(_sg_heh);
  HFH hf_query = *mesh_.hehf_iter(_sg_heh);
  CH ch_query = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_query));
  ALGOHEX_DEBUG_ONLY(std::cerr << "Searching guiding singular arc..." << std::endl;)
  int sg_rt_ax = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, _sg_heh, hf_query);
  ALGOHEX_DEBUG_ONLY(std::cerr << "hes " << _sg_heh << " ch_query " << ch_query << " hf_query " << hf_query << " rt ax "
                               << sg_rt_ax;)

  //search in singular halfedge direction
  int search_ax = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, _sg_heh, ch_query, sg_rt_ax);

  if(valence_[eh_query] < 0)
    search_ax = negate(search_ax);

  int rt_ax = sg_rt_ax;

  bool negated = rt_ax == negate(search_ax);

  //one ring cells
  std::set<CH> or_chs = get_onering_cells(mesh_, vh_f);

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vh_f, or_chs, ch_query, ft_hfhs);

  //find a start feature halfedge
  HEH he_s(-1);
  HFH hf_s(-1), hf_e(-1);
  CH ch_s(-1);
  int search_ax_s = -1;
  int rt_ax_s = -1;

  bool has_in_dir = fac_.has_singular_edge_in_direction_longest(_sg_heh, ch_query, search_ax, vh_f, ch_s, hf_s, he_s,
                                                                search_ax_s);

  if (!has_in_dir)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "no sge in direction!" << std::endl;)
    return std::vector<GDINFO>{};
  }
  //when at feature node, it's not flat, e.g. i21b
  if (n_incident_feature_edges(mesh_, feature_edge_, vh_f) <= 5)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << " he kept ? " << keep_on_feature_face_[mesh_.edge_handle(he_s)];)
    if (!keep_on_feature_face_[mesh_.edge_handle(he_s)])
      return std::vector<GDINFO>{};

    int n_nffsge = n_incident_non_ffe_singular_edges_in_same_region(mesh_, feature_face_edge_, valence_, vh_f, por_chs);

    if (n_nffsge > 1)
      return std::vector<GDINFO>{};
  }


  if (has_in_dir)
  {
    auto sg_hehs = get_halfedges_on_singular_arc_directional(mesh_, valence_, he_s);
    VH vh_end = mesh_.halfedge(sg_hehs.back()).to_vertex();
    bool legal_end = false;
    if (mesh_.is_boundary(vh_end))
    {
      legal_end = true;
    }

    //TODO: end at singular node which has detachable singular edge and has high valence edge of the same type
    if (n_incident_singular_edges(mesh_, valence_, vh_end) > 2 && feature_node_[vh_end])
      legal_end = true;

    if (legal_end)
    {
      //split
      for (const auto hehi: sg_hehs)
      {
        auto vhfi = mesh_.halfedge(hehi).from_vertex();
        SplitHelperT<MeshT>::split_one_ring(mesh_, es_, vhfi, valence_, feature_face_edge_, feature_face_vertex_,
                                            feature_edge_, sgl_vt_, feature_edge_vertex_);
      }
      SplitHelperT<MeshT>::split_one_ring(mesh_, es_, vh_end, valence_, feature_face_edge_, feature_face_vertex_,
                                          feature_edge_, sgl_vt_, feature_edge_vertex_);

      //get guides
      _end_type = 0;

      std::vector<GDINFO> guides;

      CH ch_prev = ch_s;

      for (const auto hehi: sg_hehs)
      {
        auto ehi = mesh_.edge_handle(hehi);

        CH ch_i(-1);
        HFH hf_i(-1);
        int ax_i = -1;
        if (feature_face_edge_[ehi])
        {
          //same sector
          VH vh_i = mesh_.halfedge(hehi).from_vertex();
          //one ring cells
          std::set<CH> or_chs_i = get_onering_cells(mesh_, vh_i);

          std::set<HFH> ft_hfhs_i;
          auto por_chs_i = get_part_of_onering_cells(mesh_, feature_fprop_, vh_i, or_chs_i, ch_prev, ft_hfhs_i);


          HFH ft_hf;
          for (auto hehf_it = mesh_.hehf_iter(hehi); hehf_it.valid(); ++hehf_it)
          {
            if (ft_hfhs_i.find(*hehf_it) != ft_hfhs_i.end())
            {
              ft_hf = *hehf_it;
            }
          }

          ch_i = mesh_.incident_cell(ft_hf);
          hf_i = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(ft_hf, hehi));
        }
        else
        {
          hf_i = *mesh_.hehf_iter(hehi);
          ch_i = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_i));


          if (valence_[mesh_.edge_handle(hehi)] >= 2)
          {
            //find face with min rt angle
            double min_ra = std::numeric_limits<double>::max();
            for (auto hehf_it = mesh_.hehf_iter(hehi); hehf_it.valid(); ++hehf_it)
            {
              double ra = sge_.frame_smoothness_energy_at_face(mesh_.face_handle(*hehf_it));
              if (ra < min_ra)
              {
                hf_i = *hehf_it;
                ch_i = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
                min_ra = ra;

              }
//                                std::cerr<<"rt agl "<<ra*180./M_PI<<std::endl;
            }
//                            std::cerr<<"min rt agl "<<min_ra*180./M_PI<<std::endl;
          }
        }

        ax_i = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  hehi, hf_i);
//                    std::cerr<<"negated? "<<negated<<std::endl;
        int e_ax_i = negated ? negate(ax_i) : ax_i;
        int trans_ax_i = ax_i;
        if (valence_[mesh_.edge_handle(hehi)] % 4 == 3)
        {
          trans_ax_i = negate(trans_ax_i);
        }
        else if (valence_[mesh_.edge_handle(hehi)] % 4 == 2)
        {
          e_ax_i = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, hehi, ch_i, ax_i);
          trans_ax_i = negate(e_ax_i);
        }

        ch_prev = ch_i;

        hf_i = mesh_.adjacent_halfface_in_cell(mesh_.opposite_halfface_handle(hf_i), hehi);

        guides.emplace_back(hehi, hf_i, -1, e_ax_i, trans_ax_i, false);
      }

      return guides;
    }
    else
    {
      std::cerr << "Warning: projected singular arc doesn't end at boundary." << std::endl;
    }
  }

  return std::vector<GDINFO>{};
}

template<class MeshT>
std::vector<typename FixConstrainedZipperNodeT<MeshT>::GDINFO>
FixConstrainedZipperNodeT<MeshT>::search_for_guides_on_feature_surface(const HEH _sg_heh, int &_end_type,
                                                                       HEH &_end_ft_he, HFH &_end_ft_hf,
                                                                       const bool _open_sec, const bool _same_dir)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Searching a face strip on feature surface..." << std::endl;)
  _end_type = -1;

  //find the start guide
  HFH hf_ss = *mesh_.hehf_iter(_sg_heh);
  CH ch_ss = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_ss));
  int sg_rt_ax = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, _sg_heh, hf_ss);

  //search in singular halfedge direction or opposite
  int search_ax = sg_rt_ax;
  if (valence_[mesh_.edge_handle(_sg_heh)] == 1)
    search_ax = negate(sg_rt_ax);
  if (!_same_dir)
    search_ax = negate(search_ax);

  VH vh_f = mesh_.halfedge(_sg_heh).from_vertex();

  return search_for_guides_on_feature_surface(_sg_heh, ch_ss, vh_f, search_ax, _end_type, _end_ft_he, _end_ft_hf,
                                              _open_sec);
}


template<class MeshT>
std::vector<typename FixConstrainedZipperNodeT<MeshT>::GDINFO>
FixConstrainedZipperNodeT<MeshT>::search_for_guides_on_feature_surface(const HEH _sg_heh, const CH _ch, const VH _vh,
                                                                       const int _search_ax, int &_end_type,
                                                                       HEH &_end_ft_he, HFH &_end_ft_hf,
                                                                       const bool _open_sec)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Searching a face strip on feature surface..." << std::endl;)
  std::vector<GDINFO> fn_gds;
  _end_type = -1;
  VH other_vh(-1);

  //in the direction of rt axis, this will increase sector angle by +90
  int rt_ax = _search_ax;
  if (_open_sec)
    rt_ax = negate(rt_ax);

  bool negated = _search_ax == negate(rt_ax);

  if (fac_.has_special_edge_in_direction(_vh, _ch, _search_ax))
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Found feature edge in the same search direction, return..." << std::endl;)

    return fn_gds;
  }

  auto dhi_max = furthest_parallel_sector_in_direction(_vh, _ch, _search_ax, true);

  ALGOHEX_DEBUG_ONLY(std::cerr << " vh " << _vh << " max_dist  " << std::get<0>(dhi_max) << std::endl;)


  if (std::get<0>(dhi_max) <= 0.174)
  {//10 degree
    ALGOHEX_DEBUG_ONLY(
            std::cerr << " Start facee not found in the edge direction!" << std::get<0>(dhi_max) << std::endl;)
    return fn_gds;
  }


  if (std::get<1>(dhi_max).is_valid())
  {
    HFH hf_s = std::get<1>(dhi_max);
    HEH he_s(-1);
    for (auto hfhe_it = mesh_.hfhe_iter(hf_s); hfhe_it.valid(); ++hfhe_it)
    {
      if (mesh_.halfedge(*hfhe_it).from_vertex() == _vh)
      {
        he_s = *hfhe_it;
        break;
      }
    }

    //split start edge
    auto hf_vhs = mesh_.get_halfface_vertices(hf_s, _vh);
    if (feature_edge_vertex_[hf_vhs[1]] || sgl_vt_[hf_vhs[1]])
    {
      VH vh_m = fs_.face_split(mesh_.face_handle(hf_s));
      hf_vhs[1] = vh_m;
      hf_s = mesh_.find_halfface(hf_vhs);
      he_s = mesh_.find_halfedge(_vh, vh_m);
    }
    if (feature_edge_vertex_[hf_vhs[2]] || sgl_vt_[hf_vhs[2]])
    {
      VH vh_m = fs_.face_split(mesh_.face_handle(hf_s));
      hf_vhs[2] = vh_m;
      hf_s = mesh_.find_halfface(hf_vhs);
    }

    ALGOHEX_DEBUG_ONLY(
            std::cerr << " he_s " << he_s << "hf_s " << hf_s << " max dist " << std::get<0>(dhi_max) << " ax "
                      << std::get<2>(dhi_max) << std::endl;)

    //normal axis
    CH ch_s = mesh_.incident_cell(hf_s);
    int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_s],
                                                     mesh_.normal(hf_s)).second;

    //find guides on feature surface (no halfedge information yet)
    int search_ax_s = std::get<2>(dhi_max);
    int rt_ax_s = search_ax_s;
    if (negated)
      rt_ax_s = negate(search_ax_s);

    auto dpath = shortest_face_path_on_feature_surface(GDINFO(he_s, hf_s, nm_axis, search_ax_s, rt_ax_s), other_vh,
                                                       100.);

    if (dpath.empty())
      return std::vector<GDINFO>{};

    //get singular edge free path on the feature surface
    auto dpath2 = get_special_vertex_free_path_on_feature_surface(dpath, _vh, other_vh);

    if (dpath2.empty())
    {
      std::cerr << "Error: dpath face after split is empty" << std::endl;

      return std::vector<GDINFO>{};
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "Searching for guides after split" << std::endl;)
    std::vector<GDINFO> guides;
    guides.emplace_back(HEH(-1), hf_s, nm_axis, search_ax_s, rt_ax_s);

    HFH hf_cur = hf_s;
    HEH he_cur(-1);
    int trans = 0;
    int rt_ax_i = -1;
    for (auto i = 0u; i < dpath2.size() - 1; ++i)
    {
      auto cm_eh = common_edge(dpath2[i], dpath2[i + 1]);
      if (!cm_eh.is_valid())
      {
        std::cerr << "Warning: no common edge found on the dual path at position " << i << " of " << dpath2.size()
                  << std::endl;
        return fn_gds;
      }

      for (auto hfhe_it = mesh_.hfhe_iter(hf_cur); hfhe_it.valid(); ++hfhe_it)
      {
        if (mesh_.edge_handle(*hfhe_it) == cm_eh)
        {
          he_cur = *hfhe_it;
          break;
        }
      }

      if (!he_cur.is_valid())
      {
        std::cerr << "Warning: no common edge found in the first face" << std::endl;
        return fn_gds;
      }

      HFH hf_end(-1);
      for (auto hehf_it = mesh_.hehf_iter(he_cur); hehf_it.valid(); ++hehf_it)
      {
        if (mesh_.face_handle(*hehf_it) == dpath2[i + 1])
        {
          hf_end = *hehf_it;
          break;
        }
      }

      auto hf_start = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_cur, he_cur));
      trans = tq_.mult_transitions_idx(
              EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, he_cur, hf_start, hf_end),
              trans);

      int nm_ax_cur = tq_.axis_after_transition(nm_axis, trans);
      int he_ax_cur = tq_.axis_after_transition(search_ax_s, trans);
      rt_ax_i = he_ax_cur;
      if (negated)
        rt_ax_i = negate(rt_ax_i);

      hf_cur = mesh_.opposite_halfface_handle(hf_end);

      guides.emplace_back(HEH(-1), hf_cur, nm_ax_cur, he_ax_cur, rt_ax_i);
    }

    //map
    std::map<HEH, GDINFO> he_to_gd;
    for (const auto &gd: guides)
    {
      for (auto hfhe_it = mesh_.hfhe_iter(gd.hfh); hfhe_it.valid(); ++hfhe_it)
      {
        he_to_gd.emplace(*hfhe_it, gd);
      }
    }

    ALGOHEX_DEBUG_ONLY(
            for (auto&[hehi, gdi]: he_to_gd)
              std::cerr << "hehi " << hehi << " " << gdi;
    )

    //find a bounding arc
    std::set<EH> end_ehs;
    HEH he_end(-1);
    for (auto hfhe_it = mesh_.hfhe_iter(guides.back().hfh); hfhe_it.valid(); ++hfhe_it)
    {
      auto hfe = mesh_.edge_handle(*hfhe_it);
      if (feature_edge_[hfe] || valence_[hfe] != 0)
      {
        he_end = *hfhe_it;
        end_ehs.insert(hfe);
      }
    }

    if (he_end.is_valid())
    {
      HFH hf_end(-1);
      for (auto hehf_it = mesh_.hehf_iter(he_end); hehf_it.valid(); ++hehf_it)
      {
        if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0 && *hehf_it != guides.back().hfh)
        {
          hf_end = *hehf_it;
          break;
        }
      }

      auto hf_start = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(guides.back().hfh, he_end));
      int trans0 = EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, he_end, hf_start, hf_end);

      int nm_ax_trans = tq_.axis_after_transition(guides.back().normal_axis, trans0);
      int he_ax_trans = tq_.axis_after_transition(guides.back().search_axis, trans0);

      int nm_ax_cur = AxisAlignmentHelpers::closest_axis(
              cell_quaternions_[mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_end))],
              mesh_.normal(hf_end)).second;
    }


    if (end_ehs.empty() && !other_vh.is_valid())
    {
      std::cerr << "Error: the face strip is not ending anywhere." << std::endl;
      return fn_gds;
    }
    //mark bounding edges
    std::map<EH, int> eh_count;
    for (const auto &gd: guides)
    {
      for (auto hfe_it = mesh_.hfe_iter(gd.hfh); hfe_it.valid(); ++hfe_it)
      {
        if (eh_count.find(*hfe_it) != eh_count.end())
          eh_count[*hfe_it] += 1;
        else
          eh_count[*hfe_it] = 1;
      }
    }

    ALGOHEX_DEBUG_ONLY(
            std::cerr << "end ehs; ";
            for (auto ec: end_ehs)
              std::cerr << " " << ec;
            std::cerr << std::endl;

            std::cerr << "bd ehs; ";
            for (auto ec: eh_count)
              if (ec.second == 1)
                std::cerr << " " << ec.first;
            std::cerr << std::endl;
    )

    //check if target vertex is on the bounding edges
    if (other_vh.is_valid())
    {
      bool found_tg = false;
      for (auto ec: eh_count)
        if (ec.second == 1)
        {
          if (mesh_.edge(ec.first).from_vertex() == other_vh ||
              mesh_.edge(ec.first).to_vertex() == other_vh)
          {
            found_tg = true;
            break;
          }
        }

      if (!found_tg)
      {
        std::cerr << "Error: target vertex not found on bounding edges of face strip" << std::endl;
        return fn_gds;
      }
    }

    //find a start edge
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto ve = mesh_.edge_handle(*voh_it);
      if (eh_count.find(ve) != eh_count.end() && eh_count[ve] == 1 && he_to_gd.find(*voh_it) != he_to_gd.end())
      {
        he_s = *voh_it;
        break;
      }
    }

    ALGOHEX_DEBUG_ONLY(
            std::cerr << "hes " << mesh_.halfedge(he_s).from_vertex() << " " << mesh_.halfedge(he_s).to_vertex()
                      << std::endl;)

    std::queue<HEH> he_que;
    he_que.push(he_s);

    while (!he_que.empty())
    {
      auto he_cur = he_que.front();
      he_que.pop();

      auto eh_cur = mesh_.edge_handle(he_cur);
//                std::cerr<<" "<<eh_cur;

      //find incident halfface
      if (he_to_gd.find(he_cur) == he_to_gd.end())
      {
        std::cerr << "Error: couldn't find he in gd map, eh " << he_cur << std::endl;
        return std::vector<GDINFO>{};
      }
      GDINFO gdd = he_to_gd.find(he_cur)->second;

      fn_gds.emplace_back(he_cur, gdd.hfh, gdd.normal_axis, gdd.search_axis, gdd.rt_axis);

      VH vh_to = mesh_.halfedge(he_cur).to_vertex();
      if (vh_to == other_vh)
        break;

      for (auto voh_it = mesh_.voh_iter(vh_to); voh_it.valid(); ++voh_it)
      {
        auto ve = mesh_.edge_handle(*voh_it);
        if (eh_count[ve] == 1 && ve != eh_cur)
        {
          if (!other_vh.is_valid())
          {
            if (end_ehs.find(ve) == end_ehs.end())
              he_que.push(*voh_it);
          }
          else
            he_que.push(*voh_it);
        }
      }
    }
  }

  //distinguish end type
  if (!fn_gds.empty())
  {
    if (other_vh.is_valid())
      _end_type = 0;
    else
    {
      ALGOHEX_DEBUG_ONLY(
              std::cerr << "End at fev: " << mesh_.halfedge(fn_gds.back().heh).to_vertex() << " orth patch type";)

      int sec_angle = -1;
      int type = find_orthogonal_sector(fn_gds.back(), _end_ft_he, _end_ft_hf, sec_angle, _open_sec);

      ALGOHEX_DEBUG_ONLY(std::cerr << " end type " << type << std::endl;)

      if (type > 0)
      {
        _end_type = type;
      }
      else
      {
        std::cerr << "Warning: no othogonal sector found at vertex " << mesh_.halfedge(fn_gds.back().heh).to_vertex()
                  << std::endl;
        fn_gds.clear();
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(
          std::cout << " fn_gds cells: ";
          for (const auto &gd: fn_gds)
            std::cout << " " << mesh_.incident_cell(gd.hfh);
          std::cerr << std::endl;
  )

  return fn_gds;
}

template<class MeshT>
VH
FixConstrainedZipperNodeT<MeshT>::unzip_arc_with_guiding_info(const HEH _heh, const std::vector<GDINFO> &_guides,
                                                              const int _end_type, const HEH _end_ft_heh,
                                                              const HFH _end_ft_hfh, VH &_vh_hint)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "changing matching, end type " << _end_type << std::endl;)

  //split cells to prevent 180 degree rotation
  //split incident cells
  int n_cell_split = 0;
  std::set<CH> chs_split;
  CellSplitT<MeshT> cs(mesh_);
  std::set<EH> ehs_set;

  for (const auto gdi: _guides)
  {
    EH ehi = mesh_.edge_handle(gdi.heh);
    ehs_set.insert(ehi);
    sge_.compute_edge_valence(ehi);
    MeshPropertiesT<MeshT>::update_singular_vertex_property(mesh_.edge(ehi).from_vertex());
    MeshPropertiesT<MeshT>::update_singular_vertex_property(mesh_.edge(ehi).to_vertex());

    n_cell_split += SplitHelperT<MeshT>::split_edge_cells_at_high_valence(mesh_, cs, gdi.heh, feature_fprop_, valence_,
                                                                          feature_face_edge_);
  }

  if (n_cell_split > 0)
  {
    //smooth quaternion
    QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                       feature_fprop_, feature_edge_,
                                                       ehs_set);
  }

  std::set<VH> vhs_cd;
  std::set<EH> ehs_cd;
  std::vector<VH> vhs_dt;

  //TODO: exclude existing sges from the detaching
  int n_split = 0;
  int i = 0;
  CH ch_s(-1);
  for (const auto &gd: _guides)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << gd;)
    int rt_axis = gd.rt_axis;
    if (gd.ccw)
      rt_axis = negate(rt_axis);

    int trans = rt_axis + 4;
    HFH hf_mdf = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(gd.hfh, gd.heh));

    if (!hf_mdf.is_valid())
    {
      std::cerr << "Error: halfface for matching change is not found gd" << std::endl;
      return VH(-1);
    }

    //split face/edge if the third vertex in face is singular
    VH vh_thrd = mesh_.halfedge(mesh_.next_halfedge_in_halfface(gd.heh, hf_mdf)).to_vertex();
    if (sgl_vt_[vh_thrd])
    {
      VH vhf = mesh_.halfedge(gd.heh).from_vertex();
      VH vht = mesh_.halfedge(gd.heh).to_vertex();
      auto hfvhs = mesh_.get_halfface_vertices(hf_mdf);
      VH vh_new = fs_.face_split(mesh_.face_handle(hf_mdf));
      n_split++;

      //get new face
      if (vh_new.is_valid())
      {
        for (auto i = 0u; i < hfvhs.size(); ++i)
        {
          if (hfvhs[i] != vhf && hfvhs[i] != vht)
          {
            hfvhs[i] = vh_new;
          }
        }

        hf_mdf = mesh_.find_halfface(hfvhs);
      }
    }

    if (!hf_mdf.is_valid())
    {
      std::cerr << "Error: halfface for matching change is not found" << std::endl;
      return VH(-1);
    }

    //store first cell
    if (i == 0)
      ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_mdf));

    trans_prop_[hf_mdf] = tq_.mult_transitions_idx(trans_prop_[hf_mdf], trans);
    trans_prop_[mesh_.opposite_halfface_handle(hf_mdf)] = tq_.inverse_transition_idx(trans_prop_[hf_mdf]);


    for (auto hfe_it = mesh_.hfe_iter(hf_mdf); hfe_it.valid(); ++hfe_it)
      ehs_cd.insert(*hfe_it);

    for (auto hfv_it = mesh_.hfv_iter(hf_mdf); hfv_it.valid(); ++hfv_it)
      vhs_cd.insert(*hfv_it);

    VH vh0 = mesh_.halfedge(gd.heh).from_vertex();
    vhs_dt.push_back(vh0);
    i++;
  }


  //
  if (!_guides.empty())
  {
    //
    VH vh_end = mesh_.halfedge(_guides.back().heh).to_vertex();
    if (_guides.back().ccw)
    {
      vhs_dt.push_back(mesh_.halfedge(_guides.front().heh).to_vertex());
      vh_end = mesh_.halfedge(_guides.back().heh).from_vertex();
    }
    vhs_dt.push_back(vh_end);

    //update valence
    for (const auto ehi: ehs_cd)
      sge_.compute_edge_valence(ehi);

    for (const auto &vh: vhs_cd)
      MeshPropertiesT<MeshT>::update_singular_vertex_property(vh);

    //smooth quaternion
    QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                       feature_fprop_, feature_edge_,
                                                       ehs_cd);

    //detach at the constrained tp
    VH vh_start = mesh_.halfedge(_heh).from_vertex();
    bool dsa = fis_.detach_singular_arcs_of_zipper_node(vh_start, _heh, std::vector<HEH>{},
                                                        false, true);
    if (!dsa)
      dsa = fis_.detach_singular_arc_from_singular_node(vh_start, _heh, std::vector<HEH>{},
                                                        false, true);
    bool extend = false;
    {
      int val = valence_[mesh_.edge_handle(_guides.back().heh)];
      //convexly orthogonal
      if (((_end_type == 1 || _end_type == 3) && val == 1))
        extend = true;
    }
    ALGOHEX_DEBUG_ONLY(std::cerr << "###end edge " << mesh_.edge_handle(_guides.back().heh) << " val "
                                 << valence_[mesh_.edge_handle(_guides.back().heh)] << " end type " << _end_type
                                 << std::endl;)

    HFH hf_mdf = mesh_.opposite_halfface_handle(
            mesh_.adjacent_halfface_in_cell(_guides.back().hfh, _guides.back().heh));
    ALGOHEX_DEBUG_ONLY(std::cerr << " hf_mdf " << hf_mdf;)

    HEH end_non_ff_sge(-1);
    if (!_guides.back().ccw)
      end_non_ff_sge = mesh_.next_halfedge_in_halfface(_guides.back().heh, hf_mdf);
    else
      end_non_ff_sge = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_in_halfface(_guides.back().heh, hf_mdf));
    ALGOHEX_DEBUG_ONLY(std::cerr << " end_non_ff_sge " << end_non_ff_sge;)

    if (mesh_.is_boundary(vh_end))
    {
      if (end_non_ff_sge.is_valid())
      {
        //detach the non ff sge introduced by changing matching
        bool dsa = fis_.detach_singular_arc_from_singular_node(vh_end, end_non_ff_sge, std::vector<HEH>{},
                                                               false, true);
        if (!dsa)
          dsa = fis_.detach_singular_arcs_of_zipper_node(vh_end, end_non_ff_sge, std::vector<HEH>{},
                                                         false, true);
      }

      //TODO: specified start halfedge
      bool suc = fis_.detach_singular_arc_from_singular_node_wrt_boundary(vh_end, true, false);
      if (!suc)
        suc = fis_.detach_singular_arc_from_singular_node_wrt_boundary(vh_end, false, false);

      for (auto voh_it = mesh_.voh_iter(vh_end); voh_it.valid(); ++voh_it)
      {
        VH vht = mesh_.halfedge(*voh_it).to_vertex();
        if (sgl_vt_[vht] && valence_[mesh_.edge_handle(*voh_it)] == 0 && mesh_.is_boundary(vht))
        {
          _vh_hint = vht;
          break;
        }
      }
    }
    else if (extend)
    {
      //extend the zipper node to the other region if possible (reduce singularity complexity at feature nodes)
      HEH he_end = _guides.back().heh;
      int edir_end = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                        valence_, he_end, hf_mdf);
      if (valence_[mesh_.edge_handle(he_end)] == 1)
      {
        edir_end = negate(edir_end);
      }
      if (_guides.back().ccw)
        edir_end = negate(edir_end);

      if (!_end_ft_heh.is_valid() || mesh_.is_deleted(_end_ft_heh))
      {
        std::cerr << "Error: end feature edge of orthogonal sector is not valid, maybe split!" << std::endl;
        return VH(-1);
      }
      HEH he_ft = _end_ft_heh;
      HFH hf_orth = _end_ft_hfh;

      ALGOHEX_DEBUG_ONLY(std::cerr << "end type " << _end_type << " fh orth " << mesh_.face_handle(hf_orth) << " fe "
                                   << mesh_.edge_handle(he_ft) << " ccw " << _guides.back().ccw
                                   << std::endl;)

      if (he_ft.is_valid() && hf_orth.is_valid())
      {
        HEH he_orth(-1);
        if (_guides.back().ccw)
        {
          he_orth = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_in_halfface(he_ft, hf_orth));
        }
        else
          he_orth = mesh_.next_halfedge_in_halfface(he_ft, hf_orth);

        if (!he_orth.is_valid())
        {
          std::cerr << "Error: non halfedge found in feature sector!" << std::endl;
          return VH(-1);
        }

        //find feature faces
        std::vector<HFH> ft_hfhs;
        for (auto hehf_it = mesh_.hehf_iter(he_orth); hehf_it.valid(); ++hehf_it)
        {
          if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
            ft_hfhs.push_back(*hehf_it);
        }

        if (ft_hfhs.size() != 2)
        {
          std::cerr << "Error: non-manifold feature surface!" << std::endl;
          return VH(-1);
        }

        bool negate = false;
        if (_end_type == 1 || _end_type == 3)
          negate = false;
        else if (_end_type == 2 || _end_type == 4)
          negate = true;

        ALGOHEX_DEBUG_ONLY(std::cerr << "Extend zipper node to other side, shrink sector " << negate << std::endl;)

        for (const auto hfhi: ft_hfhs)
          tu_.transversal_unzipping(he_orth, hfhi, negate);
      }

      if (end_non_ff_sge.is_valid())
      {
        //detach the non ff sge introduced by changing matching
        bool dsa = fis_.detach_singular_arc_from_singular_node(vh_end, end_non_ff_sge, std::vector<HEH>{},
                                                               false, true);
        if (!dsa)
          dsa = fis_.detach_singular_arcs_of_zipper_node(vh_end, end_non_ff_sge, std::vector<HEH>{},
                                                         false, true);
      }
    }


    ALGOHEX_DEBUG_ONLY(std::cerr << " split " << n_split << " edges" << std::endl;)

    for (const auto &vhi: vhs_dt)
    {
      bool suc = fis_.detach_singular_arc_from_singular_node(vhi, false, true);
      if (!suc && feature_face_vertex_[vhi])
        suc = fis_.detach_singular_arc_from_singular_node(vhi, false, false);

      if (!suc)
        suc = fis_.detach_singular_arcs_of_zipper_node(vhi, false, true);
    }

    return vh_end;
  }

  return VH(-1);
}

template<class MeshT>
void
FixConstrainedZipperNodeT<MeshT>::
get_special_sectors_at_special_edge_vertex(const VH vh, std::vector<std::vector<HEH>> &v_sec_hehs,
                                           std::vector<std::set<FH>> &v_sec_fhs) const
{
  if (!feature_edge_vertex_[vh] && sgl_vt_[vh] == 0)
    return;

  HEH sp_heh_s(-1);
  for (auto voh_it = mesh_.voh_iter(vh); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh_.edge_handle(*voh_it);
    if (feature_edge_[ve] || (std::abs(valence_[ve]) == 1
                              && keep_on_feature_face_[ve]
    ))
    {
      sp_heh_s = *voh_it;
      break;
    }
  }
//std::cerr<<"sp_heh_s "<<sp_heh_s;
  if (sp_heh_s.is_valid())
  {
    std::set<FH> visited_fhs;

    HEH he_it = sp_heh_s;
    while (he_it.is_valid())
    {
      std::set<FH> sec_fhs;
      auto sec_hehs = get_halfedges_in_feature_sector(mesh_, feature_fprop_, feature_edge_,
                                                      valence_, he_it, visited_fhs, sec_fhs);


//                std::cerr << " checking ehs " << mesh_.edge_handle(sec_hehs.front()) << " "
//                          << mesh_.edge_handle(sec_hehs.back()) << std::endl;
      v_sec_hehs.push_back(sec_hehs);
      v_sec_fhs.push_back(sec_fhs);

      he_it = sec_hehs.back();

      if (he_it == sp_heh_s)
        break;
    }
  }
}

//1: concavely
//2: convexly
template<class MeshT>
void
FixConstrainedZipperNodeT<MeshT>::
get_orthogonal_sectors(const HEH _heh, std::vector<std::vector<HEH>> &_v_cvx_orth_sec_hehs,
                       std::vector<std::set<FH>> &_v_cvx_orth_sec_fhs,
                       std::vector<std::vector<HEH>> &_v_ccv_orth_sec_hehs,
                       std::vector<std::set<FH>> &_v_ccv_orth_sec_fhs) const
{
  VH vh_f = mesh_.halfedge(_heh).from_vertex();
  auto or_chs = get_onering_cells(mesh_, vh_f);
  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vh_f, or_chs, *mesh_.hec_iter(_heh),
                                           ft_hfhs);

  //find a sector normal is the same as edge dir
  std::vector<std::vector<HEH>> v_sec_hehs;
  std::vector<std::set<FH>> v_sec_fhs;
  get_special_sectors_at_special_edge_vertex(vh_f, v_sec_hehs, v_sec_fhs);

  HFH hf_s = *mesh_.hehf_iter(_heh);
  CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_s));

  int e_rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                     valence_, _heh, hf_s);

  if (valence_[mesh_.edge_handle(_heh)] == 1)
    e_rt_axis = negate(e_rt_axis);

  for (auto j = 0u; j < v_sec_fhs.size(); ++j)
  {
    HFH hf_in_sec;
    auto fh0 = *v_sec_fhs[j].begin();
    auto hfh0 = mesh_.halfface_handle(fh0, 0);
    auto hfh1 = mesh_.halfface_handle(fh0, 1);

    if (ft_hfhs.find(hfh0) != ft_hfhs.end())
    {
      hf_in_sec = hfh0;
    }
    else if (ft_hfhs.find(hfh1) != ft_hfhs.end())
    {
      hf_in_sec = hfh1;
    }

    if (hf_in_sec.is_valid())
    {
      int is_prl = fac_.angle_of_axis_in_cell_and_normal_direction(e_rt_axis, ch_s, hf_in_sec,
                                                                   por_chs);

      //convexly orthogonal: fixing constrained tp will result in a zero sector. So project to surface
      if (is_prl == 2)
      {
        _v_cvx_orth_sec_hehs.push_back(v_sec_hehs[j]);
        _v_cvx_orth_sec_fhs.push_back(v_sec_fhs[j]);
      }
      else if (is_prl == 1)
      {
        _v_ccv_orth_sec_hehs.push_back(v_sec_hehs[j]);
        _v_ccv_orth_sec_fhs.push_back(v_sec_fhs[j]);
      }
    }
  }
}

template<class MeshT>
VH
FixConstrainedZipperNodeT<MeshT>::
extend_with_zipper_node(const HEH _sg_heh)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "\nExtending at edge " << mesh_.edge_handle(_sg_heh) << " val "
                               << valence_[mesh_.edge_handle(_sg_heh)] << std::endl;)

  VH vh_f = mesh_.halfedge(_sg_heh).from_vertex();
  auto or_chs = get_onering_cells(mesh_, vh_f);
  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vh_f, or_chs, *mesh_.hec_iter(_sg_heh),
                                           ft_hfhs);
  auto sg_hehs = incident_non_ffe_singular_edges_in_same_region(mesh_, feature_face_edge_, valence_, vh_f, por_chs);
  if (sg_hehs.size() > 1)
  {
    ALGOHEX_DEBUG_ONLY(
            std::cerr << "Num incident sge in the same region " << sg_hehs.size() << ", continue... " << std::endl;)
    return VH(-1);
  }

  //
  std::vector<std::vector<HEH>> v_cvx_orth_sec_hehs, v_ccv_orth_sec_hehs;
  std::vector<std::set<FH>> v_cvx_orth_sec_fhs, v_ccv_orth_sec_fhs;
  get_orthogonal_sectors(_sg_heh, v_cvx_orth_sec_hehs, v_cvx_orth_sec_fhs, v_ccv_orth_sec_hehs, v_ccv_orth_sec_fhs);


  ALGOHEX_DEBUG_ONLY(
          for (auto j = 0u; j < v_cvx_orth_sec_hehs.size(); ++j)
          {
            std::cerr << " cvx sec : " << mesh_.edge_handle(v_cvx_orth_sec_hehs[j].front()) << " "
                      << mesh_.edge_handle(v_cvx_orth_sec_hehs[j].back()) << std::endl;
          }
          for (auto j = 0u; j < v_ccv_orth_sec_hehs.size(); ++j)
          {
            std::cerr << " ccv sec : " << mesh_.edge_handle(v_ccv_orth_sec_hehs[j].front()) << " "
                      << mesh_.edge_handle(v_ccv_orth_sec_hehs[j].back()) << std::endl;
          }
  )

  int val = valence_[mesh_.edge_handle(_sg_heh)];

  HEH he_orth(-1);
  bool shrink = false;
  if (val == -1)
  {
    if (!v_cvx_orth_sec_hehs.empty())
    {
      he_orth = v_cvx_orth_sec_hehs[0][1];//initialize
      for (auto j = 0u; j < v_cvx_orth_sec_hehs.size(); ++j)
      {//prefer zero sector
        int fd_angle = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_,
                                                                                          trans_prop_,
                                                                                          valence_,
                                                                                          v_cvx_orth_sec_hehs[j],
                                                                                          v_cvx_orth_sec_fhs[j]);

        if (fd_angle == 0)
        {
          he_orth = v_cvx_orth_sec_hehs[j][1];
          break;
        }
      }
    }
    else if (!v_ccv_orth_sec_hehs.empty())
    {
      HEH he_best(-1);//initialize
      int fd_max = 1;
      for (auto j = 0u; j < v_ccv_orth_sec_hehs.size(); ++j)
      {
        int fd_angle = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_,
                                                                                          trans_prop_, valence_,
                                                                                          v_ccv_orth_sec_hehs[j],
                                                                                          v_ccv_orth_sec_fhs[j]);

        if (fd_max < fd_angle)
        {
          he_best = v_ccv_orth_sec_hehs[j][1];
          fd_max = fd_angle;
        }
      }

      if (he_best.is_valid())
      {//shrink
        he_orth = he_best;
        shrink = true;
      }
    }
  }
  else if (val == 1)
  {
    if (!v_ccv_orth_sec_hehs.empty())
    {
      he_orth = v_ccv_orth_sec_hehs[0][1];
      for (auto j = 0u; j < v_ccv_orth_sec_hehs.size(); ++j)
      {
        int fd_angle = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_,
                                                                                          trans_prop_,
                                                                                          valence_,
                                                                                          v_ccv_orth_sec_hehs[j],
                                                                                          v_ccv_orth_sec_fhs[j]);
        if (fd_angle == 0)
        {
          he_orth = v_ccv_orth_sec_hehs[j][1];
          break;
        }
      }
    }
    else if (!v_cvx_orth_sec_hehs.empty())
    {//shrink
      HEH he_best(-1);//initialize
      int fd_max = 1;
      for (auto j = 0u; j < v_cvx_orth_sec_hehs.size(); ++j)
      {
        int fd_angle = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_,
                                                                                          trans_prop_, valence_,
                                                                                          v_cvx_orth_sec_hehs[j],
                                                                                          v_cvx_orth_sec_fhs[j]);

        if (fd_max < fd_angle)
        {
          he_best = v_cvx_orth_sec_hehs[j][1];
          fd_max = fd_angle;
        }
      }

      if (he_best.is_valid())
      {//shrink
        he_orth = he_best;
        shrink = true;
      }
    }
  }

  if (he_orth.is_valid())
  {
    if (feature_edge_[mesh_.edge_handle(he_orth)])
    {
      std::cerr << "Error: no non-feature edge found in feature sector!" << std::endl;
      return VH(-1);
    }

    //find feature faces
    std::vector<HFH> ft_hfhs;
    for (auto hehf_it = mesh_.hehf_iter(he_orth); hehf_it.valid(); ++hehf_it)
    {
      if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
        ft_hfhs.push_back(*hehf_it);
    }

    if (ft_hfhs.size() != 2)
    {
      std::cerr << "Error: non-manifold feature surface!" << std::endl;
      return VH(-1);
    }

    std::vector<HEH> sghehs_new;
    for (const auto hfhi: ft_hfhs)
    {
      HEH sghe_new = tu_.transversal_unzipping(he_orth, hfhi, shrink);
      if (sghe_new.is_valid())
        sghehs_new.push_back(sghe_new);
    }

    //detach the non ff sge introduced by changing matching
    VH vh_end = mesh_.halfedge(_sg_heh).from_vertex();
    for (const auto hei: sghehs_new)
    {
      bool suc = fis_.detach_singular_arc_from_singular_node(vh_end, hei, std::vector<HEH>{},
                                                             false, true);
      if (!suc)
        suc = fis_.detach_singular_arcs_of_zipper_node(vh_end, hei, std::vector<HEH>{},
                                                       false, true);

      if (!suc)
        std::cerr << "Warning: after creating zipper node, the singular arc failed to detach!" << std::endl;
    }

    return vh_end;
  }

  return VH(-1);
}


template<class MeshT>
std::vector<FH>
FixConstrainedZipperNodeT<MeshT>::shortest_face_path_on_feature_surface(const GDINFO &_gd, VH &_vh_target,
                                                                        const double _weight)
{
  bool connect_other_zipper_node = true;

  std::map<HFH, double> min_dist;
  std::map<HFH, HALFFACEINFO> previous_hfinfo;

  //directions in parametrization coordinates
  int u = _gd.search_axis;
  int v = _gd.normal_axis;
  int w = the_third_axis((AxisAlignment) u, (AxisAlignment) v);

  HFH hfh_s = _gd.hfh;
  CH ch_s = mesh_.incident_cell(hfh_s);

  bool is_bdy = mesh_.is_boundary(ch_s);

  ALGOHEX_DEBUG_ONLY(
          std::cout << "\n\nfh_s: " << mesh_.face_handle(hfh_s) << " ch " << ch_s << " main dir: " << u << std::endl;
          std::cout << " u dir: "
                    << AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_s], (AxisAlignment) u).transpose()
                    << std::endl;
  )

  //initialize
  min_dist[hfh_s] = 0.;
  previous_hfinfo.emplace(hfh_s, HALFFACEINFO(HFH(-1), Vec3d(0., 0., 0.)));

  ALGOHEX_DEBUG_ONLY(std::cout << "Floodfill faces: axis " << u << " \n";)
  HALFFACEINFO hfi_target(HFH(-1), Vec3d(0., 0., 0.)),
          hfi_target2(HFH(-1), Vec3d(0., 0., 0.));
  std::set<DHFINFO> que;
  HALFFACEINFO hfhinfo(hfh_s, Vec3d(0., 0., 0.));
  que.insert(DHFINFO(0., hfhinfo));

  double search_dist = std::numeric_limits<double>::infinity();
  bool found_tv = false;

  while (!que.empty())
  {
    auto dc_cur = *que.begin();

    que.erase(que.begin());

    auto fh_cur = mesh_.face_handle(dc_cur.second.hfh);
//            std::cerr<<" "<<fh_cur;
    //for each edge on the current face
    for (auto hfhe_it = mesh_.hfhe_iter(dc_cur.second.hfh); hfhe_it.valid(); ++hfhe_it)
    {
      auto he_opp = mesh_.opposite_halfedge_handle(*hfhe_it);

      //get the next face
      HFH hfh_next;
      for (auto hehf_it = mesh_.hehf_iter(he_opp); hehf_it.valid(); ++hehf_it)
      {
        auto fh_i = mesh_.face_handle(*hehf_it);
        if (feature_fprop_[fh_i] > 0 && fh_i != fh_cur)
        {
          hfh_next = *hehf_it;
          break;
        }
      }

      //if visited
      if (min_dist.find(hfh_next) != min_dist.end())
        continue;

      double inc_value = 0.;
      auto hi_next = next_halfface_info(dc_cur.second, u, v, w, he_opp, hfh_next, _weight, 1., inc_value);

      //angle threshold for matching
      if (inc_value < -0.5)
        continue;

      if (connect_other_zipper_node)
      {
        //check if it ends at singular vertex which has a singular edge available for detaching
        //get singular vertices on face
        std::vector<VH> sg_vhs;
        for (auto hfv_it = mesh_.hfv_iter(hi_next.hfh); hfv_it.valid(); ++hfv_it)
        {
          if (sgl_vt_[*hfv_it] && feature_edge_vertex_[*hfv_it])
          {
            sg_vhs.push_back(*hfv_it);
          }
        }

        if (!sg_vhs.empty())
        {
          int axis_trans = tq_.axis_after_transition(u, hi_next.trans);
          for (auto vhi: sg_vhs)
          {
            bool has_dsgv = has_detachable_non_feature_face_singular_edge_at_vertex(mesh_.incident_cell(hi_next.hfh),
                                                                                    axis_trans,
                                                                                    vhi);
            if (has_dsgv)
            {
              //cannot have incoming feature edge in the same direction as the search direction
              bool has_prl_se = fac_.has_special_edge_in_direction(vhi, mesh_.incident_cell(hi_next.hfh),
                                                                   negate(axis_trans));
              if (!has_prl_se)
              {
                _vh_target = vhi;
                std::cout << " found other vh " << _vh_target << std::endl;
                hfi_target = hi_next;

                min_dist[hfh_next] = hi_next.length;

                previous_hfinfo.emplace(hi_next.hfh, dc_cur.second);
                break;
              }
            }
          }

          if (_vh_target.is_valid())
            break;
        }
      }

      EH eh_next = mesh_.edge_handle(he_opp);

      if (valence_[eh_next] != 0 || feature_edge_[eh_next])
      {//when touches singular or feature edge
        if (is_feature_edge_orthogonal_to_search_direction(dc_cur.second, eh_next, u))
        {
          if (is_face_found(hi_next, u))
          {
            hfi_target = hi_next;
            ALGOHEX_DEBUG_ONLY(
                    std::cout << " found face " << mesh_.face_handle(hfi_target.hfh) << std::endl;)
            min_dist[hfh_next] = hi_next.length;

            previous_hfinfo.emplace(hi_next.hfh, dc_cur.second);
            search_dist = 1.1 * dc_cur.second.length;

            break;
          }
          else if (is_bdy)
            continue;
        }
        else
          continue;
      }

      //flood fill
      //absolute distance to origin

      min_dist[hfh_next] = hi_next.length;

      previous_hfinfo.emplace(hi_next.hfh, dc_cur.second);


      que.insert(std::make_pair(min_dist[hfh_next], hi_next));
    }

    if (hfi_target.hfh.is_valid())
      break;
  }

  if (!_vh_target.is_valid() && hfi_target.hfh.is_valid() && connect_other_zipper_node)
  {
//            std::cout << " continue searching other vh " << _vh_target << std::endl;

    std::map<HFH, double> min_dist2;
    std::map<HFH, HALFFACEINFO> previous_hfinfo2;

    min_dist2[hfh_s] = 0.;
    previous_hfinfo2.emplace(hfh_s, HALFFACEINFO(HFH(-1), Vec3d(0., 0., 0.)));


    que.insert(DHFINFO(0., hfhinfo));

    while (!que.empty())
    {
      auto dc_cur = *que.begin();

      que.erase(que.begin());

      if (dc_cur.second.length > search_dist)
        continue;

      auto fh_cur = mesh_.face_handle(dc_cur.second.hfh);

      for (auto hfhe_it = mesh_.hfhe_iter(dc_cur.second.hfh); hfhe_it.valid(); ++hfhe_it)
      {
        auto he_opp = mesh_.opposite_halfedge_handle(*hfhe_it);

        HFH hfh_next;
        for (auto hehf_it = mesh_.hehf_iter(he_opp); hehf_it.valid(); ++hehf_it)
        {
          auto fh_i = mesh_.face_handle(*hehf_it);
          if (feature_fprop_[fh_i] > 0 && fh_i != fh_cur)
          {
            hfh_next = *hehf_it;
            break;
          }
        }


        //visited?
        if (min_dist2.find(hfh_next) != min_dist2.end())
          continue;

        //next halfface data
        double inc_value = 0.;
        auto hi_next = next_halfface_info(dc_cur.second, u, v, w, he_opp, hfh_next, _weight, 0.5, inc_value);

        if (inc_value < -0.5)
          continue;

        //end at singular vertex which has a singular edge in the same direction
        std::vector<VH> sg_vhs;
        for (auto hfv_it = mesh_.hfv_iter(hi_next.hfh); hfv_it.valid(); ++hfv_it)
        {
          if (sgl_vt_[*hfv_it] && feature_edge_vertex_[*hfv_it])
          {
            sg_vhs.push_back(*hfv_it);
          }
        }

        if (!sg_vhs.empty())
        {
          int axis_trans = tq_.axis_after_transition(u, hi_next.trans);
          for (auto vhi: sg_vhs)
          {
            bool has_dsgv = has_detachable_non_feature_face_singular_edge_at_vertex(mesh_.incident_cell(hi_next.hfh),
                                                                                    axis_trans,
                                                                                    vhi);
            //found
            if (has_dsgv)
            {
              //cannot have incoming feature edge in the same direction as the search direction
              bool has_prl_se = fac_.has_special_edge_in_direction(vhi, mesh_.incident_cell(hi_next.hfh),
                                                                   negate(axis_trans));
              if (!has_prl_se)
              {
                _vh_target = vhi;
                std::cout << " found other vh in the second try " << _vh_target << std::endl;

                hfi_target2 = hi_next;

                min_dist2[hfh_next] = hi_next.length;

                previous_hfinfo2.emplace(hi_next.hfh, dc_cur.second);
                break;
              }
            }
          }
        }

        if (_vh_target.is_valid())
          break;


        //cannot tangentially cross feature patch
        EH eh_next = mesh_.edge_handle(he_opp);
        if (valence_[eh_next] != 0 || feature_edge_[eh_next])
        {//when touches singular or feature edge
          if (!is_feature_edge_orthogonal_to_search_direction(dc_cur.second, eh_next, u))
            continue;
          else
          {
            //face normals should agree
            if (is_face_found(hi_next, u)) //face normal (anti-)parallel to search direction, continue...
              continue;
          }
        }

        //update
        min_dist2[hfh_next] = hi_next.length;

        previous_hfinfo2.emplace(hi_next.hfh, dc_cur.second);


        que.insert(std::make_pair(min_dist2[hfh_next], hi_next));
      }

      if (_vh_target.is_valid())
        break;
    }

    //collect guides
    if (_vh_target.is_valid())
    {
      //if the end face is at the start vertex, return empty
      VH vh_f = mesh_.halfedge(_gd.heh).from_vertex();
      for (auto hfv_it = mesh_.hfv_iter(hfi_target2.hfh); hfv_it.valid(); ++hfv_it)
      {
        if (*hfv_it == vh_f)
        {
          std::cerr << "Warning: end at the start vertex, return empty" << std::endl;
          return std::vector<FH>{};
        }
      }

      //get the path
      std::vector<FH> rpath;
      auto hfi_it = hfi_target2;
      while (hfi_it.hfh.is_valid())
      {
        rpath.push_back(mesh_.face_handle(hfi_it.hfh));
        auto hfinfo = previous_hfinfo2.find(hfi_it.hfh);
        if (hfinfo != previous_hfinfo2.end())
          hfi_it = hfinfo->second;
        else
          break;
      }

      std::vector<FH> path(rpath.rbegin(), rpath.rend());
      return path;
    }
  }


  //collect guides
  if (hfi_target.hfh.is_valid())
  {
    //if the end face is at the start vertex, return empty
    VH vh_f = mesh_.halfedge(_gd.heh).from_vertex();
    for (auto hfv_it = mesh_.hfv_iter(hfi_target.hfh); hfv_it.valid(); ++hfv_it)
    {
      if (*hfv_it == vh_f)
      {
        std::cerr << "Warning: end at the start vertex, return empty" << std::endl;
        return std::vector<FH>{};
      }
    }

    //get the path
    std::vector<FH> rpath;
    auto hfi_it = hfi_target;
    while (hfi_it.hfh.is_valid())
    {
      rpath.push_back(mesh_.face_handle(hfi_it.hfh));
      auto hfinfo = previous_hfinfo.find(hfi_it.hfh);
      if (hfinfo != previous_hfinfo.end())
        hfi_it = hfinfo->second;
      else
        break;
    }

    std::vector<FH> path(rpath.rbegin(), rpath.rend());
    return path;
  }
  else
    std::cerr << "No face strip found" << std::endl;


  return std::vector<FH>();
}

template<class MeshT>
std::vector<FH>
FixConstrainedZipperNodeT<MeshT>::
get_special_vertex_free_path_on_feature_surface(const std::vector<FH> &_input_path, const VH _vh_s, const VH _vh_t)
{
  if (_input_path.empty())
  {
    return _input_path;
  }
  std::set<FH> dp_fhs(_input_path.begin(), _input_path.end());

  if (!_vh_t.is_valid())
    dp_fhs.erase(_input_path.back());

  int n_split_all = 0;
  //exclude edges from splitting
  std::vector<EH> dp_ehs;
  for (auto i = 0u; i < _input_path.size() - 1; ++i)
  {
    EH ce = common_edge(_input_path[i], _input_path[i + 1]);
    if (ce.is_valid() && (feature_edge_[ce] || valence_[ce] != 0))
      dp_ehs.push_back(ce);
  }

  for (auto j = 0u; j < dp_ehs.size(); ++j)
  {
    //split the feature or singular edges on the dual path
    VH vh_t0 = mesh_.edge(dp_ehs[j]).from_vertex();
    VH vh_t1 = mesh_.edge(dp_ehs[j]).to_vertex();
    int n_sge0 = n_incident_singular_edges(mesh_, valence_, vh_t0);
    int n_sge1 = n_incident_singular_edges(mesh_, valence_, vh_t1);
    int val_le = valence_[dp_ehs[j]];
    //either vertex is bad
    if ((val_le == 0 && (sgl_vt_[vh_t0] || sgl_vt_[vh_t1])) ||
        (val_le != 0 && (n_sge0 > 2 || n_sge1 > 2)) ||
        (feature_edge_[dp_ehs[j]] && (feature_node_[vh_t0] || feature_node_[vh_t1])))
    {
      VH vhm = split_edge_with_face_property_update(dp_ehs[j], dp_fhs);
      n_split_all++;

      EH eh_after = dp_ehs[j];
      if ((val_le == 0 && sgl_vt_[vh_t0] && sgl_vt_[vh_t1]) ||
          (val_le != 0 && n_sge0 > 2 && n_sge1 > 2) ||
          (feature_node_[vh_t0] && feature_node_[vh_t1]))
      {//both vertices are bad, split once more
        eh_after = mesh_.edge_handle(mesh_.find_halfedge(vhm, vh_t1));
        VH vhm2 = split_edge_with_face_property_update(eh_after, dp_fhs);
        eh_after = mesh_.edge_handle(mesh_.find_halfedge(vhm, vhm2));

        dp_ehs[j] = eh_after;

        n_split_all++;
      }
      else
      {//either vertex is bad, get the good edge
        VH vhm2 = vh_t0;
        if ((val_le == 0 && sgl_vt_[vh_t0]) || (val_le != 0 && n_sge0 > 2) || (feature_node_[vh_t0]))
          vhm2 = vh_t1;

        eh_after = mesh_.edge_handle(mesh_.find_halfedge(vhm, vhm2));

        dp_ehs[j] = eh_after;
      }
    }
  }

  //exclude these edges from splitting
  std::set<EH> excl_ehs(dp_ehs.begin(), dp_ehs.end());

  //exclude singular or feature edge vertices
  std::set<VH> excl_vhs;
  excl_vhs.insert(_vh_s);

  for (const auto ehi: excl_ehs)
  {
    excl_vhs.insert(mesh_.edge(ehi).from_vertex());
    excl_vhs.insert(mesh_.edge(ehi).to_vertex());
  }

  EH eh_last(-1);
  if (_vh_t.is_valid())
    excl_vhs.insert(_vh_t);
  else         //if end at feature or singular edge
    eh_last = dp_ehs.back();

  ALGOHEX_DEBUG_ONLY(
          std::cerr << "excl vhs ";
          for (auto vhi: excl_vhs)
            std::cerr << " " << vhi;
  )

  //store target face
  std::set<FH> target_fhs;
  for (auto ef_it = mesh_.ef_iter(eh_last); ef_it.valid(); ++ef_it)
  {
    if (feature_fprop_[*ef_it] > 0 && dp_fhs.find(*ef_it) != dp_fhs.end())
    {
      target_fhs.insert(*ef_it);
      break;
    }
  }

  //store edges to split
  std::set<EH> ehs_split;


  //split regular edges that has two feature vertices
  for (const auto &fhi: dp_fhs)
  {
    for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
    {
      if (feature_edge_[*fe_it] || valence_[*fe_it] != 0)
        continue;

      auto from_vh = mesh_.edge(*fe_it).from_vertex();
      auto to_vh = mesh_.edge(*fe_it).to_vertex();

      bool fv_bad = is_other_feature_edge_vertex(from_vh, excl_vhs)
                    || is_other_singular_vertex(from_vh, excl_vhs);
      bool tv_bad = is_other_feature_edge_vertex(to_vh, excl_vhs)
                    || is_other_singular_vertex(to_vh, excl_vhs);
      if (fv_bad && tv_bad)
        ehs_split.insert(*fe_it);
    }
  }

  int n_split = 0;
  n_split = split_edges_with_face_property_update(dp_fhs, ehs_split, target_fhs);
  n_split_all += n_split;
  ALGOHEX_DEBUG_ONLY(
          std::cerr << "split " << n_split << " edges with both vertices being singular or feature edge vertex"
                    << std::endl;)

  ehs_split.clear();
  for (const auto &fhi: dp_fhs)
  {
    for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
    {
      if (feature_edge_[*fe_it] || valence_[*fe_it] != 0)
        continue;

      auto from_vh = mesh_.edge(*fe_it).from_vertex();
      auto to_vh = mesh_.edge(*fe_it).to_vertex();

      bool fv_bad = is_other_feature_edge_vertex(from_vh, excl_vhs)
                    || is_other_singular_vertex(from_vh, excl_vhs);
      bool tv_bad = is_other_feature_edge_vertex(to_vh, excl_vhs)
                    || is_other_singular_vertex(to_vh, excl_vhs);

      if (fv_bad || tv_bad)
        ehs_split.insert(*fe_it);
    }
  }

  n_split = split_edges_with_face_property_update(dp_fhs, ehs_split, target_fhs);
  n_split_all += n_split;

  ALGOHEX_DEBUG_ONLY(
          std::cerr << "split " << n_split << " edges with one vertex being singular or feature " << std::endl;)

  if (n_split_all > 0)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "dp faces after split: ";
                               for (auto fhi: dp_fhs)
                                 std::cerr << " " << fhi;
                               std::cerr << std::endl;
    )

    //get the first face
    FH source_fh(-1);
    for (auto vf_it = mesh_.vf_iter(_vh_s); vf_it.valid(); ++vf_it)
    {
      if (dp_fhs.find(*vf_it) != dp_fhs.end())
      {
        if (n_bad_vertices_in_face(*vf_it, excl_vhs) == 0)
        {
          source_fh = *vf_it;
          break;
        }
      }
    }
    //check if it's in dpath
    if (!source_fh.is_valid())
    {
      std::cerr << "Error: source face " << source_fh << " is not found!" << std::endl;
      return std::vector<FH>{};
    }

    //get target face
    FH target_fh(-1);
    if (_vh_t.is_valid())
    {
      for (auto vf_it = mesh_.vf_iter(_vh_t); vf_it.valid(); ++vf_it)
      {
        if (dp_fhs.find(*vf_it) != dp_fhs.end())
        {
          if (n_bad_vertices_in_face(*vf_it, excl_vhs) == 0)
          {
            target_fh = *vf_it;
            break;
          }
        }
      }
    }
    else if (eh_last.is_valid())
    {
      for (auto ef_it = mesh_.ef_iter(eh_last); ef_it.valid(); ++ef_it)
      {
        if (dp_fhs.find(*ef_it) != dp_fhs.end())
        {
          target_fh = *ef_it;
          break;
        }
      }
    }
    //check if it's in dpath
    if (!target_fh.is_valid())
    {
      std::cerr << "Error: target face is not found!" << std::endl;
      return std::vector<FH>{};
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "\nsource fh: " << source_fh << " target fh: " << target_fh << std::endl;)
    if (source_fh.is_valid() && target_fh.is_valid())
      return shortest_dual_path_from_source_face_to_target_face(source_fh, target_fh, dp_fhs, excl_vhs);
    else
      return std::vector<FH>{};
  }

  std::vector<FH> input_path(_input_path.begin(), _input_path.end());
  if (!_vh_t.is_valid())
    input_path.pop_back();

  return input_path;
}

template<class MeshT>
std::vector<FH>
FixConstrainedZipperNodeT<MeshT>::shortest_dual_path_from_source_face_to_target_face(const FH _fh_s, const FH _fh_t,
                                                                                     const std::set<FH> &_dp_fhs,
                                                                                     const std::set<VH> &_included_vhs) const
{
  //set weight of every face
  //1.has singular vertex: infinity 2.has feature vertex: 10000 3.other: 1
  std::map<FH, double> fh_wgt;
  for (const auto fhi: _dp_fhs)
  {
    bool has_other_sgv = false, has_other_fv = false;
    for (auto fv_it = mesh_.fv_iter(fhi); fv_it.valid(); ++fv_it)
    {
      if (is_other_feature_edge_vertex(*fv_it, _included_vhs))
        has_other_fv = true;

      if (is_other_singular_vertex(*fv_it, _included_vhs))
      {
        has_other_sgv = true;
        break;
      }
    }

    if (has_other_sgv)
      fh_wgt[fhi] = std::numeric_limits<double>::infinity();
    else if (has_other_fv)
      fh_wgt[fhi] = 10000.;
    else
      fh_wgt[fhi] = 1.;
  }

  std::vector<double> min_dist(mesh_.n_faces(), std::numeric_limits<double>::infinity());
  std::vector<FH> previous(mesh_.n_faces(), FH(-1));

  //check
  if (feature_fprop_[_fh_s] == 0)
    std::cerr << "Error: input source face is not feature face" << std::endl;
  ALGOHEX_DEBUG_ONLY(std::cout << "\nfh_s: " << _fh_s << " fh_t: " << _fh_t << std::endl;)

  min_dist[_fh_s.idx()] = 0.;
  previous[_fh_s.idx()] = FH(-1);

  HFH hfh_s = mesh_.halfface_handle(_fh_s, 0);
  //choose boundary halfface
  if (mesh_.is_boundary(_fh_s))
  {
    hfh_s = mesh_.is_boundary(mesh_.halfface_handle(_fh_s, 0)) ? mesh_.halfface_handle(_fh_s, 0) :
            mesh_.halfface_handle(_fh_s, 1);
  }

  bool found = false;

  std::queue<HFH> que;
  que.push(hfh_s);
  while (!que.empty())
  {
    auto hfh_cur = que.front();

    que.pop();

    auto fh_cur = mesh_.face_handle(hfh_cur);
    auto f_cur = fh_cur.idx();

    if (fh_cur == _fh_t)
    {
      ALGOHEX_DEBUG_ONLY(std::cout << " found " << _fh_t << std::endl;)
      found = true;

      break;
    }

    auto hehs = mesh_.halfface(hfh_cur).halfedges();
    for (int i = 0; i < 3; ++i)
    {
      auto he_opp = mesh_.opposite_halfedge_handle(hehs[i]);

      auto hfh_next = next_halfface_on_feature_surface(he_opp, mesh_.opposite_halfface_handle(hfh_cur));
      auto fh_next = mesh_.face_handle(hfh_next);

      if (_dp_fhs.find(fh_next) == _dp_fhs.end())
        continue;

      //check
      if (feature_fprop_[fh_next] == 0)
      {
        std::cerr << "Error: next face is not feature face" << std::endl;
        continue;
      }


      auto f_next = fh_next.idx();

      //
      if (min_dist[f_cur] + fh_wgt[fh_next] < min_dist[f_next])
      {
        min_dist[f_next] = min_dist[f_cur] + fh_wgt[fh_next];

        previous[f_next] = fh_cur;

        que.push(hfh_next);
      }
    }
  }

  //find cells
  //get the path
  std::vector<FH> path;
  if (found)
  {
    std::vector<FH> rpath;

    FH fh_it = _fh_t;
    while (fh_it.is_valid())
    {
      rpath.push_back(fh_it);
      fh_it = previous[fh_it.idx()];
    }

    path.reserve(rpath.size());
    for (auto r_it = rpath.rbegin(); r_it != rpath.rend(); ++r_it)
      path.push_back(*r_it);

  }

  ALGOHEX_DEBUG_ONLY(std::cout << " Face Path: ";
                             for (const auto fh: path)
                               std::cout << " " << fh;
                             std::cerr << std::endl;
  )

  //check
  for (const auto fhi: path)
  {
    for (auto fv_it = mesh_.fv_iter(fhi); fv_it.valid(); ++fv_it)
    {
      if ((sgl_vt_[*fv_it] || feature_edge_vertex_[*fv_it]) && _included_vhs.find(*fv_it) == _included_vhs.end())
      {
        std::cerr << "Error: found special vertex on the dual path" << std::endl;
        return std::vector<FH>{};
      }
    }
  }

  return path;
}

template<class MeshT>
typename AlgoHex::HALFFACEINFO
FixConstrainedZipperNodeT<MeshT>::next_halfface_info(const HALFFACEINFO &hi_cur, const int _u, const int _v,
                                                     const int _w,
                                                     const HEH _next_heh, const HFH _next_hfh, const double _vw_weight,
                                                     const double _length_scale, double &_incr_val) const
{
  //the first segment
  auto f_bct = mesh_.barycenter(mesh_.face_handle(hi_cur.hfh));
  auto e_bct = mesh_.barycenter(mesh_.edge_handle(_next_heh));
  Vec3d dir = ovm2eigen(e_bct - f_bct);

  auto ch_cur = mesh_.incident_cell(hi_cur.hfh);
  //u expressed in the current cell
  int u_cur = tq_.axis_after_transition(_u, hi_cur.trans);
  Vec3d u_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_cur], (AxisAlignment) u_cur);
  //v expressed in the current cell
  int v_cur = tq_.axis_after_transition(_v, hi_cur.trans);
  Vec3d v_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_cur], (AxisAlignment) v_cur);
  //w expressed in the current cell
  int w_cur = tq_.axis_after_transition(_w, hi_cur.trans);
  Vec3d w_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_cur], (AxisAlignment) w_cur);

  //first segment
  Vec3d param_next = hi_cur.param;
  param_next[_u / 2] += u_dir.dot(dir);
  param_next[_v / 2] += v_dir.dot(dir);
  param_next[_w / 2] += w_dir.dot(dir);


  //check if it's increasing in u dir
  Point e_dir = (mesh_.vertex(mesh_.halfedge(_next_heh).to_vertex()) -
                 mesh_.vertex(mesh_.halfedge(_next_heh).from_vertex())).normalize();
  Vec3d o_dir = ovm2eigen(mesh_.normal(hi_cur.hfh) % e_dir);
//        std::cerr<<"normal dir: "<<o_dir.transpose()<<" u dir: "<<u_dir.transpose()<<std::endl;

  _incr_val = o_dir.dot(u_dir);

  //the second segment
//        auto hfh_next = *mesh_.hehf_iter(_next_heh);
  auto ch_next = mesh_.incident_cell(_next_hfh);

  auto f_next_bct = mesh_.barycenter(mesh_.face_handle(_next_hfh));
  dir = ovm2eigen(f_next_bct - e_bct);

  HEH he_cur = mesh_.opposite_halfedge_handle(_next_heh);
  HFH hfh_ns = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hi_cur.hfh, he_cur));

  int trans_next = 0;
  if (feature_fprop_[mesh_.face_handle(hfh_ns)] == 0)
    trans_next = tq_.mult_transitions_idx(
            EdgeMonodromyHelperT::halfedge_transition_index(mesh_, tq_, trans_prop_, he_cur, hfh_ns,
                                                            mesh_.opposite_halfface_handle(_next_hfh)), hi_cur.trans);
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

  double scale = _length_scale * _length_scale;

  double dist = sqrt(scale * param_next[_u / 2] * param_next[_u / 2]
                     + scale * _vw_weight *
                       (param_next[_v / 2] * param_next[_v / 2] + param_next[_w / 2] * param_next[_w / 2]));

  HALFFACEINFO next_hi(_next_hfh, param_next);
  next_hi.trans = trans_next;
  next_hi.length = dist;

  return next_hi;
}

template<class MeshT>
EH
FixConstrainedZipperNodeT<MeshT>::common_edge(const FH _fh0, const FH _fh1) const
{
  std::set<EH> ehs0;
  for (auto fe_it = mesh_.fe_iter(_fh0); fe_it.valid(); ++fe_it)
    ehs0.insert(*fe_it);

  for (auto fe_it = mesh_.fe_iter(_fh1); fe_it.valid(); ++fe_it)
    if (ehs0.find(*fe_it) != ehs0.end())
      return *fe_it;

  return EH(-1);
}

template<class MeshT>
VH
FixConstrainedZipperNodeT<MeshT>::
split_edge_with_face_property_update(const EH _eh, std::set<FH> &_dp_fhs)
{
  if (es_.is_split_ok(_eh))
  {
    std::vector<std::vector<VH>> dp_hfvhs, tg_hfvhs;
    auto heh = mesh_.halfedge_handle(_eh, 0);
    auto vh_f = mesh_.halfedge(heh).from_vertex();
    auto vh_t = mesh_.halfedge(heh).to_vertex();

    for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
    {
      auto fh = mesh_.face_handle(*hehf_it);

      if (feature_fprop_[fh] > 0 && _dp_fhs.find(fh) != _dp_fhs.end())
      {
        auto hfvhs = mesh_.get_halfface_vertices(*hehf_it, heh);
        dp_hfvhs.push_back(hfvhs);
        _dp_fhs.erase(fh);
      }
    }


    auto vh_new = es_.edge_split(_eh);

    for (auto &hfvhs: dp_hfvhs)
    {
      auto hfh0 = mesh_.find_halfface(std::vector<VH>{vh_new, hfvhs[1], hfvhs[2]});
      _dp_fhs.insert(mesh_.face_handle(hfh0));

      hfh0 = mesh_.find_halfface(std::vector<VH>{hfvhs[0], vh_new, hfvhs[2]});
      _dp_fhs.insert(mesh_.face_handle(hfh0));
    }

    return vh_new;
  }

  return VH(-1);
}

template<class MeshT>
int
FixConstrainedZipperNodeT<MeshT>::
split_edges_with_face_property_update(std::set<FH> &_dp_fhs, const std::set<EH> &_ehs, std::set<FH> &_target_fhs)
{
  int n = 0;
  for (const auto eh: _ehs)
  {
    if (es_.is_split_ok(eh))
    {
      std::vector<std::vector<VH>> dp_hfvhs, tg_hfvhs;
      auto heh = mesh_.halfedge_handle(eh, 0);
      auto vh_f = mesh_.halfedge(heh).from_vertex();
      auto vh_t = mesh_.halfedge(heh).to_vertex();

      for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
      {
        auto fh = mesh_.face_handle(*hehf_it);

        if (feature_fprop_[fh] > 0 && _dp_fhs.find(fh) != _dp_fhs.end())
        {
          auto hfvhs = mesh_.get_halfface_vertices(*hehf_it, heh);
          dp_hfvhs.push_back(hfvhs);
          _dp_fhs.erase(fh);

          //target faces
          if (_target_fhs.find(fh) != _target_fhs.end())
          {
            tg_hfvhs.push_back(hfvhs);
            _target_fhs.erase(fh);
          }
        }
      }


      auto vh_new = es_.edge_split(eh);
      n++;

      for (auto &hfvhs: dp_hfvhs)
      {
        auto hfh0 = mesh_.find_halfface(std::vector<VH>{vh_new, hfvhs[1], hfvhs[2]});
        _dp_fhs.insert(mesh_.face_handle(hfh0));

        hfh0 = mesh_.find_halfface(std::vector<VH>{hfvhs[0], vh_new, hfvhs[2]});
        _dp_fhs.insert(mesh_.face_handle(hfh0));
      }

      //update target faces
      for (auto &hfvhs: tg_hfvhs)
      {
        auto hfh0 = mesh_.find_halfface(std::vector<VH>{vh_new, hfvhs[1], hfvhs[2]});
        _target_fhs.insert(mesh_.face_handle(hfh0));

        hfh0 = mesh_.find_halfface(std::vector<VH>{hfvhs[0], vh_new, hfvhs[2]});
        _target_fhs.insert(mesh_.face_handle(hfh0));
      }
    }
  }

  return n;
}

//TODO: fac_.axis_in_chs_expressed_in_cht is not optimum if there's other interior singular edges.
template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::
has_detachable_singular_edge_at_vertex(const CH _ch_s, int _rt_axis, const VH _vh_query)
{
  auto or_chs = get_onering_cells(mesh_, _vh_query);
  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh_query, or_chs, _ch_s, ft_hfhs);

  for (auto voh_it = mesh_.voh_iter(_vh_query); voh_it.valid(); ++voh_it)
  {
    auto ehi = mesh_.edge_handle(*voh_it);
    if (valence_[ehi] != 0 && !feature_edge_[ehi])
    {
      HFH hfhi = *mesh_.hehf_iter(*voh_it);
      CH chi;
      if (!mesh_.is_boundary(ehi))
        chi = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfhi));
      else
        chi = mesh_.incident_cell(hfhi);

      if (por_chs.find(chi) != por_chs.end())
      {
        int rt_axi = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                        valence_, *voh_it, hfhi);

        int rt_ax_in_i = fac_.axis_in_chs_expressed_in_cht(_ch_s, chi, _rt_axis, por_chs);

        if (rt_axi == rt_ax_in_i)
        {
          return true;
        }
      }
    }
  }

  return false;
}

//TODO: fac_.axis_in_chs_expressed_in_cht is not optimum if there's other interior singular edges.
template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::
has_detachable_non_feature_face_singular_edge_at_vertex(const CH _ch_s, int _rt_axis, const VH _vh_query)
{
  auto or_chs = get_onering_cells(mesh_, _vh_query);
  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh_query, or_chs, _ch_s, ft_hfhs);

  for (auto voh_it = mesh_.voh_iter(_vh_query); voh_it.valid(); ++voh_it)
  {
    auto ehi = mesh_.edge_handle(*voh_it);
    if (valence_[ehi] != 0 && !feature_face_edge_[ehi])
    {
      HFH hfhi = *mesh_.hehf_iter(*voh_it);
      CH chi = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfhi));

      if (por_chs.find(chi) != por_chs.end())
      {
        int rt_axi = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                        valence_, *voh_it, hfhi);

        int rt_ax_in_i = fac_.axis_in_chs_expressed_in_cht(_ch_s, chi, _rt_axis, por_chs);

        if (rt_axi == rt_ax_in_i)
        {
          return true;
        }
      }
    }
  }

  return false;
}

template<class MeshT>
std::vector<typename FixConstrainedZipperNodeT<MeshT>::GDINFO>
FixConstrainedZipperNodeT<MeshT>::get_guides_cw(const GDINFO &_gd, const VH _target_vh, const bool _open_sec,
                                                int &_end_type, HEH &_end_ft_he, HFH &_end_ft_hf)
{
  _end_type = -1;
  _end_ft_he = HEH(-1);
  _end_ft_hf = HFH(-1);

  if (!_gd.heh.is_valid())
    return std::vector<GDINFO>{};

  ALGOHEX_DEBUG_ONLY(std::cerr << "start gd: " << _gd;)
  std::vector<GDINFO> guides;
  guides.push_back(_gd);

  auto vh_t = mesh_.halfedge(_gd.heh).to_vertex();
  if (sgl_vt_[vh_t])
  {
    if (has_detachable_singular_edge_at_vertex(mesh_.incident_cell(_gd.hfh), negate(_gd.rt_axis), vh_t))
    {
      _end_type = 0;
      return guides;
    }
  }

  if (feature_node_[vh_t])
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "End at fn: " << vh_t << " orth patch type ";)

    int sec_angle = -1;
    int type = find_orthogonal_sector(_gd, _end_ft_he, _end_ft_hf, sec_angle, _open_sec);

    ALGOHEX_DEBUG_ONLY(std::cerr << type << std::endl;)
    if (type > 0)
    {
      _end_type = type;

      guides.push_back(_gd);
      return guides;
    }
    else if (type == 0)
    {
      bool has_e = fac_.has_special_edge_in_direction(vh_t, mesh_.incident_cell(_gd.hfh), _gd.search_axis);
      if (!has_e)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr
                                   << "End at feature node and no orthogonal sector, also no special edge continues, search on feature face instead"
                                   << std::endl;)
        _end_type = 6;

        guides.push_back(_gd);
        return guides;
      }
    }
  }

  //check loop
  std::set<VH> visited_vhs;
  visited_vhs.insert(mesh_.halfedge(_gd.heh).from_vertex());
  visited_vhs.insert(vh_t);

  EH eh_cur = mesh_.edge_handle(_gd.heh);
  CH ch_cur = mesh_.incident_cell(_gd.hfh);

  std::queue<GDINFO> que;
  que.push(_gd);

  VH vh_end(-1);
  bool found = false;
  bool loop = false;
  while (!que.empty())
  {
    auto gd_cur = que.front();
    que.pop();

    auto next_guide = get_next_parallel_special_halfedge_halfface_on_feature_surface_cw(gd_cur);
    ALGOHEX_DEBUG_ONLY(std::cerr << "next gd: " << next_guide;)
    if (next_guide.heh.is_valid())
    {
      //check loop
      vh_end = mesh_.halfedge(next_guide.heh).to_vertex();

      if (visited_vhs.find(vh_end) != visited_vhs.end())
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "A loop exists in the guiding path. Check the loop..." << std::endl;)
        loop = true;
        guides.push_back(next_guide);

        break;
      }
      else
        visited_vhs.insert(vh_end);


      //check singularity
      if (sgl_vt_[vh_end])
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "End at sgv: " << vh_end << " detachable ";)
        bool detachable = has_detachable_singular_edge_at_vertex(mesh_.incident_cell(next_guide.hfh),
                                                                 negate(next_guide.rt_axis), vh_end);
        ALGOHEX_DEBUG_ONLY(std::cerr << detachable << std::endl;)

        if (detachable)
        {
          guides.push_back(next_guide);
          _end_type = 0;

          found = true;

          break;
        }
      }

      if (feature_node_[vh_end])
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "End at fn: " << vh_end << " orth patch type ";)
        int sec_angle = -1;
        int type = find_orthogonal_sector(next_guide, _end_ft_he, _end_ft_hf, sec_angle, _open_sec);


        ALGOHEX_DEBUG_ONLY(std::cerr << " end type " << type << std::endl;)
        if (type > 0)
        {
          _end_type = type;

          found = true;
          guides.push_back(next_guide);
          break;
        }
        else if (type == 0)
        {
          bool has_e = fac_.has_special_edge_in_direction(vh_end, mesh_.incident_cell(next_guide.hfh),
                                                          next_guide.search_axis);
          if (!has_e)
          {
            _end_type = 6;

            found = true;
            guides.push_back(next_guide);
            break;
          }
        }
      }

      guides.push_back(next_guide);

      que.push(next_guide);
    }
  }

  if (found)
  {
    if (_target_vh.is_valid())
    {
      if (mesh_.halfedge(guides.back().heh).to_vertex() != _target_vh)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "Didn't find target vertex " << _target_vh << std::endl;)
        _end_type = -1;
        guides.clear();
      }
    }

    return guides;
  }

  //check loop
  if (loop)
  {
    int rt_axis = guides.front().rt_axis;

    int trans = rt_axis + 4;
    HFH hf_mdf = mesh_.opposite_halfface_handle(
            mesh_.adjacent_halfface_in_cell(guides.front().hfh, guides.front().heh));
    if (hf_mdf.is_valid())
    {
      int old_trans = trans_prop_[hf_mdf];

      trans_prop_[hf_mdf] = tq_.mult_transitions_idx(trans_prop_[hf_mdf], trans);
      trans_prop_[mesh_.opposite_halfface_handle(hf_mdf)] = tq_.inverse_transition_idx(trans_prop_[hf_mdf]);

      int sec_angle = -1;
      int type = find_orthogonal_sector(guides.back(), _end_ft_he, _end_ft_hf, sec_angle, _open_sec);

      //restore matching
      trans_prop_[hf_mdf] = old_trans;
      trans_prop_[mesh_.opposite_halfface_handle(hf_mdf)] = tq_.inverse_transition_idx(trans_prop_[hf_mdf]);

      ALGOHEX_DEBUG_ONLY(std::cerr << " loop. end type " << type << std::endl;)
      if (type > 0)
      {
        _end_type = type;

        return guides;
      }
      else if (type == 0)
      {
        bool has_e = fac_.has_special_edge_in_direction(vh_end, mesh_.incident_cell(guides.back().hfh),
                                                        guides.back().search_axis);
        if (!has_e)
        {
          _end_type = 6;

          return guides;
        }
      }
    }
  }

  return std::vector<GDINFO>{};
}

template<class MeshT>
std::vector<typename FixConstrainedZipperNodeT<MeshT>::GDINFO>
FixConstrainedZipperNodeT<MeshT>::get_guides_ccw(const GDINFO &_gd, const VH _target_vh, const bool _open_sec,
                                                 int &_end_type, HEH &_end_ft_he, HFH &_end_ft_hf)
{
  _end_type = -1;
  _end_ft_he = HEH(-1);
  _end_ft_hf = HFH(-1);

  if (!_gd.heh.is_valid())
    return std::vector<GDINFO>{};

  ALGOHEX_DEBUG_ONLY(std::cerr << "\nstart gd ccw: " << _gd;)
  std::vector<GDINFO> guides;
  guides.push_back(_gd);

  auto vh_f = mesh_.halfedge(_gd.heh).from_vertex();
  if (sgl_vt_[vh_f])
  {
    if (has_detachable_singular_edge_at_vertex(mesh_.incident_cell(_gd.hfh), negate(_gd.rt_axis), vh_f))
    {
      _end_type = 0;
      return guides;
    }
  }

  if (feature_node_[vh_f])
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "End at fn: " << vh_f << " orth patch type ";)

    int sec_angle = -1;
    int type = find_orthogonal_sector(_gd, _end_ft_he, _end_ft_hf, sec_angle, _open_sec);

    ALGOHEX_DEBUG_ONLY(std::cerr << type << std::endl;)

    if (type > 0)
    {
      _end_type = type;

      guides.push_back(_gd);
      return guides;
    }
    else if (type == 0)
    {
      bool has_e = fac_.has_special_edge_in_direction(vh_f, mesh_.incident_cell(_gd.hfh), _gd.search_axis);
      if (!has_e)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr
                                   << "End at feature node and no orthogonal sector, also no special edge continues, search on feature face instead"
                                   << std::endl;)
        _end_type = 6;

        guides.push_back(_gd);
        return guides;
      }
    }
  }

  //check loop
  std::set<VH> visited_vhs;
  visited_vhs.insert(mesh_.halfedge(_gd.heh).to_vertex());
  visited_vhs.insert(vh_f);

  EH eh_cur = mesh_.edge_handle(_gd.heh);
  CH ch_cur = mesh_.incident_cell(_gd.hfh);

  std::queue<GDINFO> que;
  que.push(_gd);

  VH vh_end(-1);
  bool found = false;
  bool loop = false;
  while (!que.empty())
  {
    auto gd_cur = que.front();
    que.pop();

    auto next_guide = get_next_parallel_special_halfedge_halfface_on_feature_surface_ccw(gd_cur);
    ALGOHEX_DEBUG_ONLY(std::cerr << "next gd: " << next_guide;)
    if (next_guide.heh.is_valid())
    {
      vh_end = mesh_.halfedge(next_guide.heh).from_vertex();

      //check loop
      if (visited_vhs.find(vh_end) != visited_vhs.end())
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "A loop exists in the guiding path. Check the loop..." << std::endl;)
        loop = true;
        guides.push_back(next_guide);

        break;
      }
      else
        visited_vhs.insert(vh_end);


      //check singularity
      if (sgl_vt_[vh_end])
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "End at sgv: " << vh_end << " detachable ";)
        bool detachable = has_detachable_singular_edge_at_vertex(mesh_.incident_cell(next_guide.hfh),
                                                                 negate(next_guide.rt_axis), vh_end);
        ALGOHEX_DEBUG_ONLY(std::cerr << detachable << std::endl;)

        if (detachable)
        {
          guides.push_back(next_guide);
          found = true;
          _end_type = 0;

          break;
        }
      }

      if (feature_node_[vh_end])
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "End at fn: " << vh_end << " orth patch type";)

        int sec_angle = -1;
        int type = find_orthogonal_sector(next_guide, _end_ft_he, _end_ft_hf, sec_angle, _open_sec);

        ALGOHEX_DEBUG_ONLY(std::cerr << " end type " << type << std::endl;)

        if (type > 0)
        {
          _end_type = type;

          found = true;
          guides.push_back(next_guide);
          break;
        }
        else if (type == 0)
        {
          bool has_e = fac_.has_special_edge_in_direction(vh_end, mesh_.incident_cell(next_guide.hfh),
                                                          next_guide.search_axis);
          if (!has_e)
          {
            _end_type = 6;

            found = true;
            guides.push_back(next_guide);
            break;
          }
        }
      }

      guides.push_back(next_guide);

      que.push(next_guide);
    }
  }

  if (found)
  {
    if (_target_vh.is_valid())
    {
      if (mesh_.halfedge(guides.back().heh).from_vertex() != _target_vh)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "Didn't find target vertex " << _target_vh << std::endl;)
        guides.clear();
        _end_type = -1;
      }
    }

    return guides;
  }

  //check loop
  if (loop)
  {
    int rt_axis = guides.front().rt_axis;
    rt_axis = negate(rt_axis);

    int trans = rt_axis + 4;
    HFH hf_mdf = mesh_.opposite_halfface_handle(
            mesh_.adjacent_halfface_in_cell(guides.front().hfh, guides.front().heh));
    if (hf_mdf.is_valid())
    {
      int old_trans = trans_prop_[hf_mdf];

      trans_prop_[hf_mdf] = tq_.mult_transitions_idx(trans_prop_[hf_mdf], trans);
      trans_prop_[mesh_.opposite_halfface_handle(hf_mdf)] = tq_.inverse_transition_idx(trans_prop_[hf_mdf]);

      int sec_angle = -1;
      int type = find_orthogonal_sector(guides.back(), _end_ft_he, _end_ft_hf, sec_angle, _open_sec);

      //restore matching
      trans_prop_[hf_mdf] = old_trans;
      trans_prop_[mesh_.opposite_halfface_handle(hf_mdf)] = tq_.inverse_transition_idx(trans_prop_[hf_mdf]);

      ALGOHEX_DEBUG_ONLY(std::cerr << " loop end type " << type << std::endl;)
      if (type > 0)
      {
        _end_type = type;

        return guides;
      }
      else if (type == 0)
      {
        bool has_e = fac_.has_special_edge_in_direction(vh_end, mesh_.incident_cell(guides.back().hfh),
                                                        guides.back().search_axis);
        if (!has_e)
        {
          _end_type = 6;

          return guides;
        }
      }
    }
  }

  return std::vector<GDINFO>{};
}


template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::has_orthogonal_sector_cw(const GDINFO &_gd, HEH &_heh, HFH &_hfh, int &_angle) const
{
  //initialize
  _heh = HEH(-1);
  _hfh = HFH(-1);
  _angle = -1;

  //axes
  int heh_axis_prev = _gd.search_axis;
  int patch_nm_axis_prev = _gd.normal_axis;
  //find next incoming special halfedge
  HEH he_it = _gd.heh;
  EH eh_it = mesh_.edge_handle(he_it);
  HFH hf_it = _gd.hfh;

  int trans_prd = 0, trans_prd_acc = 0;
  std::set<HEH> visited_fhes;
  if (feature_edge_[eh_it] || (std::abs(valence_[eh_it]) == 1 && feature_face_edge_[eh_it]))
    visited_fhes.insert(he_it);
  do
  {
    std::tie(he_it, hf_it, trans_prd) = next_outgoing_special_hehf_with_trans_cw(he_it, hf_it);
    if (!he_it.is_valid())
      return false;
    if (visited_fhes.find(he_it) != visited_fhes.end())
      break;

    visited_fhes.insert(he_it);

    trans_prd_acc = tq_.mult_transitions_idx(trans_prd, trans_prd_acc);

    CH ch_next = mesh_.incident_cell(hf_it);

    int patch_nm_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next],
                                                             mesh_.normal(hf_it)).second;
    int patch_nm_axis_prev_in_i = tq_.axis_after_transition(patch_nm_axis_prev, trans_prd_acc);

    if (patch_nm_axis_prev_in_i != patch_nm_axis_i)
      break;

    //move to next sector
    HFH hf_itt = hf_it;
    do
    {
      //next incident halfface
      hf_itt = mesh_.adjacent_halfface_in_cell(hf_itt, he_it);
      hf_itt = mesh_.opposite_halfface_handle(hf_itt);

      if (!hf_itt.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! he_s: " << he_it << " hfh_s " << hf_it
                  << std::endl;
        return false;
      }

      if (feature_fprop_[mesh_.face_handle(hf_itt)] > 0)
        break;

      trans_prd_acc = tq_.mult_transitions_idx(trans_prop_[hf_itt], trans_prd_acc);
    }
    while (hf_itt != hf_it);

    //incoming
    he_it = mesh_.opposite_halfedge_handle(he_it);
    hf_it = mesh_.opposite_halfface_handle(hf_itt);

    //compare patch normal axis of next sector with edge axis
    ch_next = mesh_.incident_cell(hf_it);

    int heh_axis_prev_in_i = tq_.axis_after_transition(heh_axis_prev, trans_prd_acc);

    patch_nm_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next],
                                                         mesh_.normal(hf_it)).second;

    //edge axis parallel to normal axis in the next sector
    if (patch_nm_axis_i / 2 == heh_axis_prev_in_i / 2)
    {
      _heh = he_it;
      _hfh = hf_it;

      if (patch_nm_axis_i != heh_axis_prev_in_i) //opposite direction
        _angle = 1;
      else //same direction
        _angle = 2;

      return true;
    }
  }
  while (he_it != _gd.heh);

  return false;
}


//find othogonal sector that can be used in extending tp
//return: 0 - no orthogonal
//        1 - has orthogonal zero sector, good for extending tp, open
//        2 - has orthogonal, 3 quater sector, good for extending tp, shrink
//        3 - has orthogonal, and suitable for extending tp, open
//        4 - has orthogonal, and suitable for extending tp, shrink
//        5 - has orthogonal, but no suitable for extending
template<class MeshT>
int
FixConstrainedZipperNodeT<MeshT>::find_orthogonal_sector(const GDINFO &_gd, HEH &_heh, HFH &_hfh, int &_sec_angle,
                                                         const bool _open_sec) const
{
  //initialize
  _heh = HEH(-1);
  _hfh = HFH(-1);
  _sec_angle = -1;

  HEH he_orth0(-1);
  HFH hf_orth0(-1);
  int otho_angle = -1;
  bool cw = !_gd.ccw;
  if (cw)
    has_orthogonal_sector_cw(_gd, he_orth0, hf_orth0, otho_angle);
  else
    has_orthogonal_sector_ccw(_gd, he_orth0, hf_orth0, otho_angle);


  if (!he_orth0.is_valid())
  {
    std::cerr << "Warning: he_orth is not valid " << std::endl;
    return 0;
  }

  int preferred_fd_agl = 0;
  if ((otho_angle == 1 && _open_sec) || (otho_angle == 2 && !_open_sec))
    preferred_fd_agl = 0;
  else if ((otho_angle == 2 && _open_sec) || (otho_angle == 1 && !_open_sec))
    preferred_fd_agl = 3;

  //first found sector angle
  int sec_agl = fac_.field_angle_of_feature_edge_sector(he_orth0, hf_orth0, cw);

  ALGOHEX_DEBUG_ONLY(std::cerr << "fd angle " << sec_agl << " eh " << mesh_.edge_handle(he_orth0) << " fh: "
                               << mesh_.face_handle(hf_orth0) << " ortho agl " << otho_angle
                               << std::endl;)

  int type = 5;

  if ((otho_angle == 1 && _open_sec) || (otho_angle == 2 && !_open_sec))
  {
    _heh = he_orth0;
    _hfh = hf_orth0;
    _sec_angle = sec_agl;
    type = 3;

    if (sec_agl == preferred_fd_agl)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "found preferred sector of angle " << preferred_fd_agl << " : eh "
                                   << mesh_.edge_handle(_heh) << " fh: " << mesh_.face_handle(_hfh) << std::endl;)
      return 1;
    }
  }
  else if ((otho_angle == 2 && _open_sec) || (otho_angle == 1 && !_open_sec))
  {
    if (sec_agl >= 2)
    {
      _heh = he_orth0;
      _hfh = hf_orth0;
      _sec_angle = sec_agl;

      type = 4;
    }

    if (sec_agl == preferred_fd_agl)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "found preferred sector of angle " << preferred_fd_agl << " : eh "
                                   << mesh_.edge_handle(_heh) << " fh: " << mesh_.face_handle(_hfh) << std::endl;)
      return 2;
    }
  }

  HEH he_it = he_orth0;
  HFH hf_it = hf_orth0;

  do
  {
    //keep searching in cw for a flat zero sector
    auto next_flat_sec = next_flat_sector_with_trans(he_it, hf_it, cw);
    if (next_flat_sec.first.is_valid())
    {
      he_it = next_flat_sec.first;
      hf_it = next_flat_sec.second;
      sec_agl = fac_.field_angle_of_feature_edge_sector(he_it, hf_it, cw);
      ALGOHEX_DEBUG_ONLY(std::cerr << "fd angle next " << sec_agl << " eh " << mesh_.edge_handle(he_it) << " fh: "
                                   << mesh_.face_handle(hf_it) << std::endl;)
      if ((otho_angle == 1 && _open_sec) || (otho_angle == 2 && !_open_sec))
      {
        if (sec_agl == preferred_fd_agl)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "found preferred sector of angle " << preferred_fd_agl
                                       << " in further searching: eh " << mesh_.edge_handle(_heh) << " fh: "
                                       << mesh_.face_handle(_hfh) << std::endl;)

          _heh = he_orth0;
          _hfh = hf_orth0;
          _sec_angle = sec_agl;

          return 1;
        }
      }
      else if ((otho_angle == 2 && _open_sec) || (otho_angle == 1 && !_open_sec))
      {
        if (!_heh.is_valid() && sec_agl >= 2)
        {
          _heh = he_orth0;
          _hfh = hf_orth0;
          _sec_angle = sec_agl;
          type = 4;
        }

        if (sec_agl == preferred_fd_agl)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "found preferred sector of angle " << preferred_fd_agl
                                       << " in further searching: eh " << mesh_.edge_handle(_heh) << " fh: "
                                       << mesh_.face_handle(_hfh) << std::endl;)

          _heh = he_orth0;
          _hfh = hf_orth0;
          _sec_angle = sec_agl;

          return 2;
        }
      }
    }
    else
      break;
  }
  while (he_it.is_valid() && he_it != he_orth0);

  return type;
}

//input gd contains the outgoing halfedge and feature halfface
//output the outgoing halfedge in the next sector if orthogonal
template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::has_orthogonal_sector_ccw(const GDINFO &_gd, HEH &_heh, HFH &_hfh, int &_angle) const
{
  //initialize
  _heh = HEH(-1);
  _hfh = HFH(-1);
  _angle = -1;

  //axes
  int sc_axis_prev = _gd.search_axis;
  int patch_nm_axis_prev = _gd.normal_axis;
  //find next incoming special halfedge
  HEH he_it = _gd.heh;
  EH eh_it = mesh_.edge_handle(he_it);
  HFH hf_it = _gd.hfh;

  int trans_prd = 0, trans_prd_acc = 0;
  std::set<HEH> visited_fhes;
  if (feature_edge_[eh_it] || (std::abs(valence_[eh_it]) == 1 && feature_face_edge_[eh_it]))
    visited_fhes.insert(he_it);
  do
  {
    std::tie(he_it, hf_it, trans_prd) = next_incoming_special_hehf_with_trans_ccw(he_it, hf_it);

    if (!he_it.is_valid())
      return false;
    if (visited_fhes.find(he_it) != visited_fhes.end())
      break;
    visited_fhes.insert(he_it);


    trans_prd_acc = tq_.mult_transitions_idx(trans_prd, trans_prd_acc);

    CH ch_next = mesh_.incident_cell(hf_it);

    int patch_nm_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next],
                                                             mesh_.normal(hf_it)).second;
    int patch_nm_axis_prev_in_i = tq_.axis_after_transition(patch_nm_axis_prev, trans_prd_acc);

    if (patch_nm_axis_prev_in_i != patch_nm_axis_i)
      break;

    //move to next sector
    HFH hf_itt = hf_it;
    do
    {
      //next incident halfface
      hf_itt = mesh_.adjacent_halfface_in_cell(hf_itt, he_it);
      hf_itt = mesh_.opposite_halfface_handle(hf_itt);

      if (!hf_itt.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! he_s: " << he_it << " hfh_s " << hf_it
                  << std::endl;
        return false;
      }

      if (feature_fprop_[mesh_.face_handle(hf_itt)] > 0)
        break;

      trans_prd_acc = tq_.mult_transitions_idx(trans_prop_[hf_itt], trans_prd_acc);
    }
    while (hf_itt != hf_it);

    //outgoing
    he_it = mesh_.opposite_halfedge_handle(he_it);
    hf_it = mesh_.opposite_halfface_handle(hf_itt);

    //compare patch normal axis
    ch_next = mesh_.incident_cell(hf_it);
//            std::cerr<< " ch_next "<<ch_next;

    //compare guiding axis
    int sc_axis_prev_in_i = tq_.axis_after_transition(sc_axis_prev, trans_prd_acc);

    patch_nm_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next],
                                                         mesh_.normal(hf_it)).second;

    //edge axis parallel to normal axis in the next sector
    if (patch_nm_axis_i / 2 == sc_axis_prev_in_i / 2)
    {
      _heh = he_it;
      _hfh = hf_it;

      if (patch_nm_axis_i != sc_axis_prev_in_i) //opposite direction
        _angle = 1;
      else //same direction
        _angle = 2;

      return true;
    }
  }
  while (he_it != _gd.heh);

  return false;
}

template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::has_special_edge_in_direction_in_sector(const VH _vh, const HFH _hfh_s,
                                                                          const int _ax) const
{
  HEH he_in, he_out;
  for (auto hfhe_it = mesh_.hfhe_iter(_hfh_s); hfhe_it.valid(); ++hfhe_it)
  {
    if (mesh_.halfedge(*hfhe_it).to_vertex() == _vh)
    {
      he_in = *hfhe_it;
      break;
    }
  }

  if (he_in.is_valid())
  {
    HEH he_right, he_left;
    HFH hf_right, hf_left;
    int trans_prd = 0;
    if (feature_edge_[mesh_.edge_handle(he_in)])
    {
      he_left = mesh_.opposite_halfedge_handle(he_in);
      hf_left = _hfh_s;

      auto dir0 = mesh_.vertex(mesh_.halfedge(he_left).to_vertex()) -
                  mesh_.vertex(mesh_.halfedge(he_left).from_vertex());
      int ax0 = AxisAlignmentHelpers::closest_axis(cell_quaternions_[mesh_.incident_cell(hf_left)],
                                                   dir0.normalize()).second;
      if (ax0 == _ax)
        return true;
    }

    //right most special edge
    std::tie(he_right, hf_right, trans_prd) = next_outgoing_special_hehf_with_trans_cw(he_in, _hfh_s);
    if (he_right.is_valid())
    {
      auto dir1 = mesh_.vertex(mesh_.halfedge(he_right).to_vertex()) -
                  mesh_.vertex(mesh_.halfedge(he_right).from_vertex());
      int ax1 = AxisAlignmentHelpers::closest_axis(
              cell_quaternions_[mesh_.incident_cell(hf_right)],
              dir1.normalize()).second;
      int axin_1 = tq_.axis_after_transition(_ax, trans_prd);
      if (axin_1 == ax1)
        return true;
    }

    if (!he_left.is_valid())
    {
      he_out = mesh_.next_halfedge_in_halfface(he_in, _hfh_s);
      std::tie(he_left, hf_left, trans_prd) = next_incoming_special_hehf_with_trans_ccw(he_out, _hfh_s);
      if (he_left.is_valid())
      {
        he_left = mesh_.opposite_halfedge_handle(he_left);
        auto dir1 = mesh_.vertex(mesh_.halfedge(he_left).to_vertex()) -
                    mesh_.vertex(mesh_.halfedge(he_left).from_vertex());
        int ax1 = AxisAlignmentHelpers::closest_axis(
                cell_quaternions_[mesh_.incident_cell(hf_left)],
                dir1.normalize()).second;
        int axin_1 = tq_.axis_after_transition(_ax, trans_prd);
        if (axin_1 == ax1)
          return true;
      }
    }

  }

  return false;
}

template<class MeshT>
std::tuple<double, HFH, int>
FixConstrainedZipperNodeT<MeshT>::find_furthest_halfface_in_direction_in_sector(const VH _vh, const HFH _hfh_s,
                                                                                const int _ax) const
{
  HEH he_in;
  for (auto hfhe_it = mesh_.hfhe_iter(_hfh_s); hfhe_it.valid(); ++hfhe_it)
  {
    if (mesh_.halfedge(*hfhe_it).to_vertex() == _vh)
    {
      he_in = *hfhe_it;
      break;
    }
  }

  //search in cw
  double max_dist = -1000000000.;
  HFH hf_best;
  int ax_in_best;

  CH ch_s = mesh_.incident_cell(_hfh_s);
  auto ax_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_s], (AxisAlignment) _ax);
  auto he_vec = mesh_.barycenter(mesh_.face_handle(_hfh_s)) - mesh_.vertex(_vh);
  he_vec.normalize();

  double dist = ovm2eigen(he_vec).dot(ax_vec);
  if (dist > max_dist)
  {
    max_dist = dist;
    hf_best = _hfh_s;
    ax_in_best = _ax;
  }


  HEH he_itt = he_in;
  HFH hf_itt = _hfh_s;
  int trans_prd = 0;
  do
  {
    he_itt = mesh_.next_halfedge_in_halfface(he_itt, hf_itt);

    //if special edge, return
    EH eh_itt = mesh_.edge_handle(he_itt);
    if (feature_edge_[eh_itt] || valence_[eh_itt] != 0)
    {
      break;
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
        return {-1., HFH(-1), -1};
      }
      if (feature_fprop_[mesh_.face_handle(hf_it_in)] > 0 || hf_it_in == hf_itt)
        break;

      trans_prd = tq_.mult_transitions_idx(trans_prop_[hf_it_in], trans_prd);
    }
    while (hf_it_in != hf_itt);

    he_itt = mesh_.opposite_halfedge_handle(he_itt);
    hf_itt = mesh_.opposite_halfface_handle(hf_it_in);

    int ax_itt = tq_.axis_after_transition(_ax, trans_prd);
    auto ax_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[mesh_.incident_cell(hf_itt)],
                                                          (AxisAlignment) ax_itt);
    auto he_vec = mesh_.barycenter(mesh_.face_handle(hf_itt)) - mesh_.vertex(_vh);
    he_vec.normalize();

    double dist = ovm2eigen(he_vec).dot(ax_vec);
    if (dist > max_dist)
    {
      max_dist = dist;
      hf_best = hf_itt;
      ax_in_best = ax_itt;
    }
  }
  while (he_itt != he_in);

  //search in ccw
  HEH he_out = mesh_.next_halfedge_in_halfface(he_in, _hfh_s);
  he_itt = he_out;
  hf_itt = _hfh_s;
  trans_prd = 0;
  do
  {
    he_itt = mesh_.prev_halfedge_in_halfface(he_itt, hf_itt);

    //if special edge, return
    EH eh_itt = mesh_.edge_handle(he_itt);
    if (feature_edge_[eh_itt] || valence_[eh_itt] != 0)
    {
      break;
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
        return {-1., HFH(-1), -1};
      }
      if (feature_fprop_[mesh_.face_handle(hf_it_in)] > 0 || hf_it_in == hf_itt)
        break;

      trans_prd = tq_.mult_transitions_idx(trans_prop_[hf_it_in], trans_prd);
    }
    while (hf_it_in != hf_itt);

    he_itt = mesh_.opposite_halfedge_handle(he_itt);
    hf_itt = mesh_.opposite_halfface_handle(hf_it_in);

    int ax_itt = tq_.axis_after_transition(_ax, trans_prd);
    auto ax_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[mesh_.incident_cell(hf_itt)],
                                                          (AxisAlignment) ax_itt);
    auto he_vec = mesh_.barycenter(mesh_.face_handle(hf_itt)) - mesh_.vertex(_vh);
    he_vec.normalize();

    double dist = ovm2eigen(he_vec).dot(ax_vec);
    if (dist > max_dist)
    {
      max_dist = dist;
      hf_best = hf_itt;
      ax_in_best = ax_itt;
    }
  }
  while (he_itt != he_out);

  return {max_dist, hf_best, ax_in_best};
}


template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::
has_convexly_parallel_sector_in_direction(const VH _vh, const CH _ch_s, const int _dir, const bool _check_spe) const
{
  auto dhi_max = furthest_parallel_sector_in_direction(_vh, _ch_s, _dir, _check_spe);
//std::cerr<<" vh "<<_vh<<" max dist "<<std::get<0>(dhi_max)<<" fh "<<mesh_.face_handle(std::get<1>(dhi_max))<<std::endl;
  if (std::get<0>(dhi_max) > 0.174 && std::get<1>(dhi_max).is_valid()) //10 degree
    return true;

  return false;
}


template<class MeshT>
std::tuple<double, HFH, int>
FixConstrainedZipperNodeT<MeshT>::furthest_parallel_sector_in_direction(const VH _vh, const CH _ch_s, const int _dir,
                                                                        const bool _check_spe) const
{
  //one ring cells
  std::set<CH> or_chs = get_onering_cells(mesh_, _vh);

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, _ch_s, ft_hfhs);

  std::tuple<double, HFH, int> dhi_max(-100000., HFH(-1), _dir);

  std::vector<std::vector<HEH>> v_sec_hehs;
  std::vector<std::set<FH>> v_sec_fhs;
  get_special_sectors_at_special_edge_vertex(_vh, v_sec_hehs, v_sec_fhs);

  for (auto j = 0u; j < v_sec_fhs.size(); ++j)
  {
    HFH hfi0 = mesh_.halfface_handle(*v_sec_fhs[j].begin(), 0);
    HFH hfi1 = mesh_.halfface_handle(*v_sec_fhs[j].begin(), 1);
    CH ch0 = mesh_.incident_cell(hfi0);
    CH ch1 = mesh_.incident_cell(hfi1);

    HFH hfi;
    CH ch_i;
    if (por_chs.find(ch0) != por_chs.end())
    {
      hfi = hfi0;
      ch_i = ch0;
    }
    else
    {
      hfi = hfi1;
      ch_i = ch1;
    }

    int sc_ax_in_i = fac_.axis_in_chs_expressed_in_cht(_ch_s, ch_i, _dir, por_chs);

    int nm_ax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_i],
                                                   mesh_.normal(hfi)).second;
    //orthogonal to face normal
    if (nm_ax / 2 != sc_ax_in_i / 2)
    {
      if (_check_spe && has_special_edge_in_direction_in_sector(_vh, hfi, sc_ax_in_i))
        continue;

      auto dhi_tmp = find_furthest_halfface_in_direction_in_sector(_vh, hfi, sc_ax_in_i);

      if (std::get<1>(dhi_tmp).is_valid() && std::get<0>(dhi_tmp) > std::get<0>(dhi_max))
      {
        dhi_max.swap(dhi_tmp);
      }
    }
  }

  return dhi_max;
}

template<class MeshT>
typename FixConstrainedZipperNodeT<MeshT>::GDINFO
FixConstrainedZipperNodeT<MeshT>::get_next_parallel_special_halfedge_halfface_on_feature_surface_cw(const GDINFO &_gd)
{
  //axes
  int patch_nm_axis_prev = _gd.normal_axis;
  int heh_axis_prev = _gd.search_axis;
  bool negated = _gd.search_axis == negate(_gd.rt_axis);
  //find next incoming special halfedge
  HEH he_it = _gd.heh;
  HFH hf_it = _gd.hfh;

  int trans_prd = 0, trans_prd_acc = 0;
  do
  {
    std::tie(he_it, hf_it, trans_prd) = next_outgoing_special_hehf_with_trans_cw(he_it, hf_it);
    trans_prd_acc = tq_.mult_transitions_idx(trans_prd, trans_prd_acc);

    //compare patch normal axis
    CH ch_next = mesh_.incident_cell(hf_it);

    int patch_nm_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next],
                                                             mesh_.normal(hf_it)).second;
    int patch_nm_axis_prev_in_i = tq_.axis_after_transition(patch_nm_axis_prev, trans_prd_acc);

    if (patch_nm_axis_prev_in_i != patch_nm_axis_i)
      break;

    //compare guiding axis
    auto dir = mesh_.vertex(mesh_.halfedge(he_it).to_vertex()) -
               mesh_.vertex(mesh_.halfedge(he_it).from_vertex());
    int heh_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next], dir.normalize()).second;

    int heh_axis_prev_in_i = tq_.axis_after_transition(heh_axis_prev, trans_prd_acc);

    if (heh_axis_i == heh_axis_prev_in_i)
    {
      int rt_ax = negated ? negate(heh_axis_i) : heh_axis_i;
      return {he_it, hf_it, patch_nm_axis_i, heh_axis_i, rt_ax, false};
    }

    //move to next sector
    HFH hf_itt = hf_it;
    do
    {
      //next incident halfface
      hf_itt = mesh_.adjacent_halfface_in_cell(hf_itt, he_it);
      hf_itt = mesh_.opposite_halfface_handle(hf_itt);

      if (!hf_itt.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! he_s: " << he_it << " hfh_s " << hf_it << std::endl;
        return {HEH(-1), HFH(-1), -1, -1, -1, false};
      }

      if (feature_fprop_[mesh_.face_handle(hf_itt)] > 0)
        break;

      trans_prd_acc = tq_.mult_transitions_idx(trans_prop_[hf_itt], trans_prd_acc);
    }
    while (hf_itt != hf_it);

    //incoming
    he_it = mesh_.opposite_halfedge_handle(he_it);
    hf_it = mesh_.opposite_halfface_handle(hf_itt);
  }
  while (he_it != _gd.heh);

  return {HEH(-1), HFH(-1), -1, -1, -1, false};
}


template<class MeshT>
typename FixConstrainedZipperNodeT<MeshT>::GDINFO
FixConstrainedZipperNodeT<MeshT>::get_next_parallel_special_halfedge_halfface_on_feature_surface_ccw(const GDINFO &_gd)
{
  //axes
  int patch_nm_axis_prev = _gd.normal_axis;
  int sc_axis_prev = _gd.search_axis;
  bool negated = _gd.search_axis == negate(_gd.rt_axis);

  //find next incoming special halfedge
  HEH he_it = _gd.heh;
  HFH hf_it = _gd.hfh;
  int trans_prd = 0, trans_prd_acc = 0;

  do
  {
    std::tie(he_it, hf_it, trans_prd) = next_incoming_special_hehf_with_trans_ccw(he_it, hf_it);

    trans_prd_acc = tq_.mult_transitions_idx(trans_prd, trans_prd_acc);

    CH ch_next = mesh_.incident_cell(hf_it);

    //compare normal axes
    int patch_nm_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next],
                                                             mesh_.normal(hf_it)).second;
    int patch_nm_axis_prev_in_i = tq_.axis_after_transition(patch_nm_axis_prev, trans_prd_acc);

    if (patch_nm_axis_prev_in_i != patch_nm_axis_i)
      break;

    //compare guiding axis
    auto dir = mesh_.vertex(mesh_.halfedge(he_it).from_vertex()) -
               mesh_.vertex(mesh_.halfedge(he_it).to_vertex());
    int search_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next], dir.normalize()).second;

    int sc_axis_prev_in_i = tq_.axis_after_transition(sc_axis_prev, trans_prd_acc);

    if (search_axis_i == sc_axis_prev_in_i)
    {
      int rt_ax = negated ? negate(search_axis_i) : search_axis_i;
      return {he_it, hf_it, patch_nm_axis_i, search_axis_i, rt_ax, true};
    }

    //move to next sector
    HFH hf_itt = hf_it;
    do
    {
      //next incident halfface
      hf_itt = mesh_.adjacent_halfface_in_cell(hf_itt, he_it);
      hf_itt = mesh_.opposite_halfface_handle(hf_itt);

      if (!hf_itt.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! he_s: " << he_it << " hfh_s " << hf_it << std::endl;
        return {HEH(-1), HFH(-1), -1, -1, -1, true};
      }

      if (feature_fprop_[mesh_.face_handle(hf_itt)] > 0)
        break;

      trans_prd_acc = tq_.mult_transitions_idx(trans_prop_[hf_itt], trans_prd_acc);
    }
    while (hf_itt != hf_it);

    //outgoing
    he_it = mesh_.opposite_halfedge_handle(he_it);
    hf_it = mesh_.opposite_halfface_handle(hf_itt);

  }
  while (he_it != _gd.heh);

  return {HEH(-1), HFH(-1), -1, -1, -1, true};
}

//input is incoming halfedge and incident halfface
template<class MeshT>
std::tuple<HEH, HFH, int>
FixConstrainedZipperNodeT<MeshT>::next_outgoing_special_hehf_with_trans_cw(const HEH _heh, const HFH _hfh) const
{
  HEH he_itt = _heh;
  HFH hf_itt = _hfh;
  int trans_prd = 0;
  do
  {
    he_itt = mesh_.next_halfedge_in_halfface(he_itt, hf_itt);

    //if special edge, return
    EH eh_itt = mesh_.edge_handle(he_itt);
    if (feature_edge_[eh_itt] || valence_[eh_itt] != 0)
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

//cw: input is incoming halfedge and incident halfface
//ccw: input is outgoing halfedge and incident halfface
template<class MeshT>
std::pair<HEH, HFH>
FixConstrainedZipperNodeT<MeshT>::next_flat_sector_with_trans(const HEH _heh, const HFH _hfh, const bool _cw) const
{
  //axes
  CH ch_cur = mesh_.incident_cell(_hfh);
  int patch_nm_axis_cur = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_cur], mesh_.normal(_hfh)).second;

  //find next incoming special halfedge
  HEH he_it = _heh;
  HFH hf_it = _hfh;

  int trans_prd = 0;
  if (_cw)
    std::tie(he_it, hf_it, trans_prd) = next_outgoing_special_hehf_with_trans_cw(he_it, hf_it);
  else
    std::tie(he_it, hf_it, trans_prd) = next_incoming_special_hehf_with_trans_ccw(he_it, hf_it);


  //move to next sector
  HFH hf_itt = hf_it;
  do
  {
    //next incident halfface
    hf_itt = mesh_.adjacent_halfface_in_cell(hf_itt, he_it);
    hf_itt = mesh_.opposite_halfface_handle(hf_itt);

    if (!hf_itt.is_valid())
    {
      std::cerr << "Error: adjacent halfface is invalid! he_s: " << he_it << " hfh_s " << hf_it << std::endl;
      return {HEH(-1), HFH(-1)};
    }

    if (feature_fprop_[mesh_.face_handle(hf_itt)] > 0)
      break;

    trans_prd = tq_.mult_transitions_idx(trans_prop_[hf_itt], trans_prd);
  }
  while (hf_itt != hf_it);

  //incoming
  he_it = mesh_.opposite_halfedge_handle(he_it);
  hf_it = mesh_.opposite_halfface_handle(hf_itt);

  //compare patch normal axis
  CH ch_next = mesh_.incident_cell(hf_it);

  int patch_nm_axis_i = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch_next],
                                                           mesh_.normal(hf_it)).second;
  int patch_nm_axis_cur_in_i = tq_.axis_after_transition(patch_nm_axis_cur, trans_prd);

  if (patch_nm_axis_cur_in_i != patch_nm_axis_i)
    return {HEH(-1), HFH(-1)};

  return {he_it, hf_it};
}

//input is outgoing
template<class MeshT>
std::tuple<HEH, HFH, int>
FixConstrainedZipperNodeT<MeshT>::next_incoming_special_hehf_with_trans_ccw(const HEH _heh, const HFH _hfh) const
{
  HEH he_itt = _heh;
  HFH hf_itt = _hfh;
  int trans_prd = 0;
  do
  {
    he_itt = mesh_.prev_halfedge_in_halfface(he_itt, hf_itt);

    //if special edge, return
    EH eh_itt = mesh_.edge_handle(he_itt);
    if (feature_edge_[eh_itt] || valence_[eh_itt] != 0)
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
HFH
FixConstrainedZipperNodeT<MeshT>::next_halfface_on_feature_surface(const HEH _heh, const HFH _hfh_s) const
{
  auto hfh_it = _hfh_s;
  if (mesh_.is_boundary(hfh_it))
    return HFH(-1);

  int idx = 0;
  do
  {
    auto hfh_adj = mesh_.adjacent_halfface_in_cell(hfh_it, _heh);
    if (!hfh_adj.is_valid())
    {
      std::cerr << "Error: adjacent halfface is invalid!" << std::endl;
      idx = -1;
      break;
    }
    hfh_it = mesh_.opposite_halfface_handle(hfh_adj);

    if (feature_fprop_[mesh_.face_handle(hfh_it)] > 0)
    {
      break;
    }
  }
  while (hfh_it != _hfh_s);

  return hfh_it;
}

template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::
is_face_found(const HALFFACEINFO &hi_cur, int _u) const
{
  auto normal = mesh_.normal(hi_cur.hfh);

  auto ca = AxisAlignmentHelpers::closest_axis(cell_quaternions_[mesh_.incident_cell(hi_cur.hfh)], normal);
  int u2 = tq_.axis_after_transition(_u, hi_cur.trans);

  if (ca.second / 2 == u2 / 2)
    return true;

  return false;
}

template<class MeshT>
bool
FixConstrainedZipperNodeT<MeshT>::
is_feature_edge_orthogonal_to_search_direction(const HALFFACEINFO &hi_cur, const EH _eh, int _u) const
{
  if ((valence_[_eh] >= -2 && valence_[_eh] <= 4 && valence_[_eh] != 0) || feature_edge_[_eh])
  {
    //align to edge
    HEH heh0 = mesh_.halfedge_handle(_eh, 0);
    auto axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cell(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                        valence_,
                                                                        feature_edge_, heh0, hi_cur.hfh);

    int u2 = tq_.axis_after_transition(_u, hi_cur.trans);

    if (axis / 2 != u2 / 2)
      return true;
  }

  return false;
}
}


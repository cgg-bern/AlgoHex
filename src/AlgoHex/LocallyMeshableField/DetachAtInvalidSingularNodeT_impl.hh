/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define FIXINVALIDSINGULARNODET_C

#include "DetachAtInvalidSingularNodeT.hh"
#include "QuaternionsSmoothing.hh"
#include "CommonFuncs.hh"

namespace AlgoHex
{
template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::fix_interior_invalid_singular_node(const VH _vh)
{
  int n = 0;
  bool fixed = false;
  while (n < 10)
  {
    int n_sge = n_incident_singular_edges(_vh);
    if (n_sge <= 2)
      break;

    bool is_meshable = is_locally_meshable_simplified(_vh, n_sge);
    if (is_meshable)
      return fixed;

    if (n_invalid_singular_edge(mesh_, valence_, _vh) > 0)
      return false;

    //check if it's a circle attached to an sg arc, let it shrink instead of fixing it
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
      if (is_singular_circle_noise(*voh_it))
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "singular circle noise, let it shrink." << std::endl;)
        return fixed;
      }

    //fix
    bool is_dtd = detach_singular_arc_from_singular_node(_vh, false, true);

    if (!is_dtd)
    {
      bool is_tp_dtd = detach_singular_arcs_of_zipper_node(_vh, false, true);
      if (is_tp_dtd)
        return true;
    }
    else
    {
      fixed = true;
      break;
    }


    n++;
  }

  return fixed;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::fix_invalid_singular_node_on_feature_surface(const VH _vh)
{
  std::queue<VH> que;
  que.push(_vh);
  int n_fixed = 0;
  while (!que.empty())
  {
    auto vh_cur = que.front();
    que.pop();

    std::set<VH> vhs_set{vh_cur};
    for (auto ve_it = mesh_.ve_iter(vh_cur); ve_it.valid(); ++ve_it)
    {
      sge_.compute_edge_valence(*ve_it);
      vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
      vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
    }
    //update singular vertex property
    for (const auto &vhi: vhs_set)
      MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);

    int n_nff_val_1, n_nff_val1, n_nff_val_u, n_nff, n_bdy_val_1, n_bdy_val1, n_bdy_val_u, n_bdy, n_cmpl, n_ff_val_1, n_ff_val1, n_ff_val_u, n_ff, n_fe;
    count_special_edges_at_vertex(mesh_, feature_face_edge_, feature_edge_, valence_, vh_cur, n_nff_val_1, n_nff_val1,
                                  n_nff_val_u, n_nff,
                                  n_bdy_val_1, n_bdy_val1, n_bdy_val_u, n_bdy, n_cmpl, n_ff_val_1, n_ff_val1,
                                  n_ff_val_u, n_ff, n_fe);

    if (n_invalid_singular_edge(mesh_, valence_, _vh) > 0)
      return false;

    int n_all_sge = n_nff + n_ff;
    bool lcm_checked = false, is_meshable = false;
    if (is_detachable_arc(vh_cur, n_nff_val_1, n_nff_val1, n_nff_val_u, n_nff))
    {
      is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
      lcm_checked = true;
      if (is_meshable)
        return n_fixed > 0;

      ALGOHEX_DEBUG_ONLY(std::cerr << "detach singular arc from sg node..." << vh_cur << std::endl;)
      bool arc_suc = detach_singular_arc_from_singular_node(vh_cur, false, true);
      if (arc_suc)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

        n_fixed++;
        que.push(vh_cur);
        continue;
      }
    }


    if (is_detachable_zipper_node(vh_cur, n_nff_val_1, n_nff_val1, n_nff))
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "detach zipper node from sg node..." << vh_cur << std::endl;)
      if (!lcm_checked)
      {
        is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
        lcm_checked = true;
      }
      if (is_meshable)
        return n_fixed > 0;

      bool tp_suc = detach_singular_arcs_of_zipper_node(vh_cur, false, true);
      if (tp_suc)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

        n_fixed++;
        que.push(vh_cur);
        continue;
      }
    }


    if (!feature_face_vertex_[vh_cur])
      return n_fixed > 0;

    if (is_detachable_arc_to_ffe(vh_cur, n_nff_val_1, n_nff_val1, n_ff_val_1, n_ff_val1, n_ff_val_u))
    {
      if (!lcm_checked)
      {
        is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
        lcm_checked = true;
      }
      if (is_meshable)
        return n_fixed > 0;

      ALGOHEX_DEBUG_ONLY(std::cerr << "detach singular arc from sg node (feature face sge)..." << vh_cur << std::endl;)
      bool arc_suc = detach_singular_arc_from_singular_node(vh_cur, true, true);
      if (arc_suc)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

        n_fixed++;
        que.push(vh_cur);
        continue;
      }
//                }
    }

    if (is_detachable_zipper_node_to_ffe(vh_cur, n_nff_val_1, n_nff_val1, n_ff_val_1, n_ff_val1, n_ff_val_u))
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "detach zipper node from sg node (feature face sge)..." << vh_cur << std::endl;)

      if (!lcm_checked)
      {
        is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
        lcm_checked = true;
      }
      if (is_meshable)
        return n_fixed > 0;

      bool tp_suc = detach_singular_arcs_of_zipper_node(vh_cur, true, true);
      if (tp_suc)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

        n_fixed++;
        que.push(vh_cur);
        continue;
      }
    }

    if (is_detachable_arc(vh_cur, n_nff_val_1, n_nff_val1, n_nff_val_u, n_nff))
    {
      if (!lcm_checked)
      {
        is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
        lcm_checked = true;
      }
      if (is_meshable)
        return n_fixed > 0;

      ALGOHEX_DEBUG_ONLY(
              std::cerr << "detach singular arc from sg node (crossing feature face)..." << vh_cur << std::endl;)
      bool arc_suc = detach_singular_arc_from_singular_node(vh_cur, false, false);
      if (arc_suc)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)
        n_fixed++;
        que.push(vh_cur);

        continue;
      }
    }

    if (is_detachable_arc_wrt_feature_face(vh_cur, n_nff_val_1, n_nff_val1, n_nff))
    {
      if (!lcm_checked)
      {
        is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
        lcm_checked = true;
      }
      if (is_meshable)
        return n_fixed > 0;

      ALGOHEX_DEBUG_ONLY(
              std::cerr << "extend to the other region and then detach singular arc or zipper node..." << vh_cur
                        << std::endl;)
      bool arc_suc = detach_singualr_arc_from_singular_node_wrt_feature_face(vh_cur, true);
      if (arc_suc)
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "fixed crossing..." << std::endl;)
        n_fixed++;
        que.push(vh_cur);

        continue;
      }
    }
  }

  return n_fixed > 0;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::fix_boundary_invalid_singular_node(const VH _vh, bool _start_from_bdy_sge)
{
  std::queue<VH> que;
  que.push(_vh);
  int n_fixed = 0;
  while (!que.empty())
  {
    auto vh_cur = que.front();
    que.pop();

    int n_nff_val_1, n_nff_val1, n_nff_val_u, n_nff, n_bdy_val_1, n_bdy_val1, n_bdy_val_u, n_bdy, n_cmpl, n_ff_val_1, n_ff_val1, n_ff_val_u, n_ff, n_fe;
    count_special_edges_at_vertex(mesh_, feature_face_edge_, feature_edge_, valence_, vh_cur, n_nff_val_1, n_nff_val1,
                                  n_nff_val_u, n_nff,
                                  n_bdy_val_1, n_bdy_val1, n_bdy_val_u, n_bdy, n_cmpl, n_ff_val_1, n_ff_val1,
                                  n_ff_val_u, n_ff, n_fe);

    if (n_invalid_singular_edge(mesh_, valence_, _vh) > 0)
      return false;

    //check if it's a circle attached to an sg arc, let it shrink instead of fixing it
    if (n_bdy < 2)
    {
      for (auto voh_it = mesh_.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
        if (is_singular_circle_noise(*voh_it))
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "singular circle noise, let it shrink." << std::endl;)
          return n_fixed;
        }
    }

    int n_all_sge = n_nff + n_ff;
    bool lcm_checked = false, is_meshable = false;
    if (!_start_from_bdy_sge)
    {
      if (is_detachable_arc(vh_cur, n_nff_val_1, n_nff_val1, n_nff_val_u, n_nff))
      {
        is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
        lcm_checked = true;
        if (is_meshable)
          return n_fixed > 0;

        ALGOHEX_DEBUG_ONLY(std::cerr << "detach singular arc from sg node..." << vh_cur << std::endl;)
        bool arc_suc = detach_singular_arc_from_singular_node(vh_cur, false, true);
        if (arc_suc)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

          n_fixed++;
          que.push(vh_cur);
          continue;
        }
      }


      if (is_detachable_zipper_node(vh_cur, n_nff_val_1, n_nff_val1, n_nff))
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "detach zipper node from sg node..." << vh_cur << std::endl;)
        if (!lcm_checked)
        {
          is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
          lcm_checked = true;
        }
        if (is_meshable)
          return n_fixed > 0;

        bool tp_suc = detach_singular_arcs_of_zipper_node(vh_cur, false, true);
        if (tp_suc)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

          n_fixed++;
          que.push(vh_cur);
          continue;
        }
      }

      if (is_detachale_arc_wrt_boundary(n_nff_val_1, n_nff_val1, n_nff, n_bdy, n_fe))
      {
        ALGOHEX_DEBUG_ONLY(
                std::cerr << "detach singular arc (anti-parallel to boundary normal) from sg node..." << vh_cur
                          << std::endl;)
        if (!(n_nff == 1 && (n_nff_val_1 == 1 || n_nff_val1 == 1) && n_fe == 2))
        {//although locally meshable, for global meshability
          if (!lcm_checked)
          {
            is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
            lcm_checked = true;
          }
          if (is_meshable)
            return n_fixed > 0;
        }

        bool arc_bdy_suc = detach_singular_arc_from_singular_node_wrt_boundary(vh_cur, true, false);
        if (arc_bdy_suc)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

          n_fixed++;
          que.push(vh_cur);
          continue;
        }
      }

      if (is_detachable_arc_to_ffe(vh_cur, n_nff_val_1, n_nff_val1, n_bdy_val_1, n_bdy_val1, n_ff_val_u) &&
          (n_bdy_val_1 + n_bdy_val1 > 0))
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "detach singular arc from sg node (to boundary sge)..." << vh_cur << std::endl;)
        if (!lcm_checked)
        {
          is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
          lcm_checked = true;
        }
        if (is_meshable)
          return n_fixed > 0;
        bool arc_suc = detach_singular_arc_from_singular_node(vh_cur, true, true);
        if (arc_suc)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

          n_fixed++;
          que.push(vh_cur);
          continue;
        }
      }

      if (is_detachable_zipper_node_to_ffe(vh_cur, n_nff_val_1, n_nff_val1, n_bdy_val_1, n_bdy_val1, n_ff_val_u))
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "detach zipper node from sg node (to boundary sge)..." << vh_cur << std::endl;)
        if (!lcm_checked)
        {
          is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
          lcm_checked = true;
        }
        if (is_meshable)
          return n_fixed > 0;
        bool tp_suc = detach_singular_arcs_of_zipper_node(vh_cur, true, true);
        if (tp_suc)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

          n_fixed++;
          que.push(vh_cur);
          continue;
        }
      }

      if (is_detachale_arc_wrt_boundary(n_nff_val_1, n_nff_val1, n_nff, n_bdy, n_fe))
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "detach singular arc (parallel to boundary normal) from sg node..." << vh_cur
                                     << std::endl;)

        if (!lcm_checked)
        {
          is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
          lcm_checked = true;
        }
        if (is_meshable)
          return n_fixed > 0;

        bool arc_bdy_suc = detach_singular_arc_from_singular_node_wrt_boundary(vh_cur, false, false);
        if (arc_bdy_suc)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

          n_fixed++;
          que.push(vh_cur);
          continue;
        }
      }
    }
    else
    {
      // in case push vertex to interior is impossible due to geometric issue
      // lowest priority: detach a boundary arc from invalid node
      if (is_detachale_arc_wrt_boundary_from_boundary_singular_edge(n_nff_val_1, n_nff_val1, n_nff, n_bdy, n_fe))
      {
        ALGOHEX_DEBUG_ONLY(std::cerr
                                   << "detach singular arc (anti-parallel to boundary normal) from sg node (start from bdy sge)..."
                                   << vh_cur << std::endl;)

        if (!lcm_checked)
        {
          is_meshable = is_locally_meshable_simplified(vh_cur, n_all_sge);
          lcm_checked = true;
        }
        if (is_meshable)
          return n_fixed > 0;

        bool arc_bdy_suc = detach_singular_arc_from_singular_node_wrt_boundary(vh_cur, true, true);
        if (arc_bdy_suc)
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "fixed..." << std::endl;)

          n_fixed++;
          que.push(vh_cur);
          continue;
        }
      }
    }
  }

  return n_fixed > 0;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::detach_singular_arc_from_singular_node(const VH _vh,
                                                                            const bool _find_featureface_sges,
                                                                            const bool _same_sector_sges)
{
  SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                          feature_edge_, sgl_vt_, feature_edge_vertex_);
  std::set<VH> vhs_set{_vh};
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
  {
    if (valence_[*ve_it] != 0)
    {
      sge_.compute_edge_valence(*ve_it);
      vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
      vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
    }
  }
  //update singular vertex property
  for (const auto &vhi: vhs_set)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);

  std::set<CH> or_chs;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    or_chs.insert(*vc_it);

  //compute halfedge axis per cell
  std::map<CH, int> cell_he_axis, cell_nm_axis;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    int val = valence_[mesh_.edge_handle(*voh_it)];
    if (val == -1 || val == 1 || val == std::numeric_limits<int>::lowest())
    {
      auto hehf_it = mesh_.hehf_iter(*voh_it);
      CH ch_s(-1);
      if (mesh_.is_boundary(*voh_it))
        ch_s = mesh_.incident_cell(*hehf_it);
      else
        ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));

      if (!ch_s.is_valid())
      {
        std::cerr << "Error: invalid start cell at edge " << mesh_.edge_handle(*voh_it) << std::endl;
        return false;
      }

      int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, *voh_it, *hehf_it);
      if (axis < 0 || axis > 5)
      {
        auto eh = mesh_.edge_handle(*voh_it);
        std::cout << "Error: singular edge: " << eh << "of valence " << valence_[eh] << " has " << axis << " axis"
                  << std::endl;
        return false;
      }

      //halfedge axis
      axis = EdgeMonodromyHelperT::halfedge_axis_idx(mesh_, cell_quaternions_, *voh_it, ch_s, axis);

      cell_he_axis[ch_s] = axis;

      for (; hehf_it.valid(); ++hehf_it)
      {
        auto chi = mesh_.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          if (chi == ch_s)
            continue;

          axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
          cell_he_axis[chi] = axis;
        }
      }
    }
  }

  //feature face aligned axis
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
      if (feature_fprop_[mesh_.face_handle(*chf_it)] > 0)
      {
        int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], mesh_.normal(hfopp)).second;
        cell_nm_axis[*vc_it] = nm_axis;
      }
    }
  }


  ALGOHEX_DEBUG_ONLY(
          std::cerr << "\nDetaching anti-parallel singular edge at vh " << _vh << " ... finding feature face edge "
                    << _find_featureface_sges << " same region " << _same_sector_sges << std::endl;)


  //find anti parallel singular edges
  set_only_feature_face_singular_edges_as_targets(_find_featureface_sges);
  set_same_sector_singular_edges_as_targets(_same_sector_sges);
  std::set<HEH> he_pairs;
  auto hfhs = ssaf_.find_anti_parallel_singular_edge_with_disk_path(_vh, or_chs, cell_he_axis, cell_nm_axis, he_pairs);

  if (hfhs.empty())
    return false;

  std::set<FH> disk_fhs;
  for (const auto hfi: hfhs)
    disk_fhs.insert(mesh_.face_handle(hfi));

  std::set<EH> sg_ehs;
  for (const auto hei: he_pairs)
    sg_ehs.insert(mesh_.edge_handle(hei));

  return detaching_singularity(_vh, disk_fhs, sg_ehs);
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::detach_singular_arc_from_singular_node(const VH _vh, const HEH _heh,
                                                                            const std::vector<HEH> &_excl_hehs,
                                                                            const bool _find_featureface_sges,
                                                                            const bool _same_sector_sges,
                                                                            const bool _split)
{

  if (mesh_.halfedge(_heh).from_vertex() != _vh)
  {
    std::cerr << "Error: from vertex of halfedge " << _heh << " is not " << _vh << std::endl;
    return false;
  }

  if (_split)
    SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                            feature_edge_, sgl_vt_, feature_edge_vertex_);

  std::set<VH> vhs_set{_vh};
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
  {
    if (valence_[*ve_it] != 0)
    {
      sge_.compute_edge_valence(*ve_it);
      vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
      vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
    }
  }
  //update singular vertex property
  for (const auto &vhi: vhs_set)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);


  std::set<CH> or_chs;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    or_chs.insert(*vc_it);

  //compute halfedge axis per cell
  std::map<CH, int> cell_he_axis, cell_nm_axis;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    int val = valence_[mesh_.edge_handle(*voh_it)];
    if (val == -1 || val == 1 || val == std::numeric_limits<int>::lowest())
    {
      auto hehf_it = mesh_.hehf_iter(*voh_it);
      CH ch_s(-1);
      if (mesh_.is_boundary(*voh_it))
        ch_s = mesh_.incident_cell(*hehf_it);
      else
        ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));

      if (!ch_s.is_valid())
      {
        std::cerr << "Error: invalid start cell at edge " << mesh_.edge_handle(*voh_it) << std::endl;
        return false;
      }

      int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, *voh_it, *hehf_it);
      if (axis < 0 || axis > 5)
      {
        auto eh = mesh_.edge_handle(*voh_it);
        std::cout << "Error: singular edge: " << eh << "of valence " << valence_[eh] << " has " << axis << " axis"
                  << std::endl;
        return false;
      }

      //halfedge axis
      if (val == 1)
        axis = negate(axis);

      cell_he_axis[ch_s] = axis;

      for (; hehf_it.valid(); ++hehf_it)
      {
        auto chi = mesh_.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          if (chi == ch_s)
            continue;

          axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
          cell_he_axis[chi] = axis;
        }
      }
    }
  }

  //feature face aligned axis
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
      if (feature_fprop_[mesh_.face_handle(*chf_it)] > 0)
      {
        int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], mesh_.normal(hfopp)).second;
        cell_nm_axis[*vc_it] = nm_axis;
      }
    }
  }


  ALGOHEX_DEBUG_ONLY(
          std::cerr << "\nDetaching anti-parallel singular edge at vh " << _vh << " eh_s " << mesh_.edge_handle(_heh)
                    << " ... finding feature face edge " << _find_featureface_sges << " same region "
                    << _same_sector_sges << std::endl;)

  //find anti parallel singular edges
  set_only_feature_face_singular_edges_as_targets(_find_featureface_sges);
  set_same_sector_singular_edges_as_targets(_same_sector_sges);
  std::set<HEH> he_pairs;
  auto hfhs = ssaf_.find_anti_parallel_singular_edge_with_disk_path(_vh, _heh, _excl_hehs, or_chs, cell_he_axis,
                                                                    cell_nm_axis, he_pairs);

  if (hfhs.empty())
    return false;

  std::set<FH> disk_fhs;
  for (const auto hfi: hfhs)
    disk_fhs.insert(mesh_.face_handle(hfi));

  std::set<EH> sg_ehs;
  for (const auto hei: he_pairs)
    sg_ehs.insert(mesh_.edge_handle(hei));

  return detaching_singularity(_vh, disk_fhs, sg_ehs);
}

//    template <class MeshT>
//    bool
//    DetachAtInvalidSingularNodeT<MeshT>::detach_singular_arc_from_singular_node(const VH _vh, MeshT& _mesh, std::vector<OVM::Vec4f>& _mesh_color, const bool _find_featureface_sges, const bool _same_sector_sges) {
//        SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_, feature_edge_, sgl_vt_, feature_edge_vertex_);
//
//        std::set<VH> vhs_set{_vh};
//        for(auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it) {
//            if(valence_[*ve_it] != 0) {
//                sge_.compute_edge_valence(*ve_it);
//                vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
//                vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
//            }
//        }
//        //update singular vertex property
//        for (const auto &vhi : vhs_set)
//            MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);
//
//        std::set<CH> or_chs;
//        for(auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
//            or_chs.insert(*vc_it);
//
//        //compute halfedge axis per cell
//        std::map<CH, int> cell_he_axis, cell_nm_axis;
//        for(auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it) {
//            int val = valence_[mesh_.edge_handle(*voh_it)];
//            if(val == -1 || val == 1 || val == std::numeric_limits<int>::lowest()) {
//                auto hehf_it = mesh_.hehf_iter(*voh_it);
//                CH ch_s(-1);
//                if(mesh_.is_boundary(*voh_it))
//                    ch_s = mesh_.incident_cell(*hehf_it);
//                else
//                    ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
//
//                if(!ch_s.is_valid()) {
//                    std::cerr<<"Error: invalid start cell at edge "<<mesh_.edge_handle(*voh_it)<<std::endl;
//                    return false;
//                }
//
//                int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_, *voh_it, *hehf_it);
//                if (axis < 0 || axis > 5) {
//                    auto eh = mesh_.edge_handle(*voh_it);
//                    std::cout << "Error: singular edge: "<< eh<< "of valence " << valence_[eh] << " has "<<axis<<" axis"
//                              << std::endl;
//                    return false;
//                }
//
//                //halfedge axis
//                if(val == 1)
//                    axis = negate(axis);
//
//                cell_he_axis[ch_s] = axis;
//
//                for (; hehf_it.valid(); ++hehf_it) {
//                    auto chi = mesh_.incident_cell(*hehf_it);
//                    if(chi.is_valid()) {
//                        if(chi == ch_s)
//                            continue;
//
//                        axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
//                        cell_he_axis[chi] = axis;
//                    }
//                }
//            }
//        }
//
//
//        DuplicateOneRingMeshT<MeshT> dp(mesh_, _mesh);
//        dp.copy_one_ring(std::vector<VH>{_vh});
//        _mesh_color.resize(_mesh.n_cells());
//
//        VH vh_or = dp.onering_mesh_vertex_handle(_vh);
//
//        std::vector<OVM::Vec4f> colors;
//        colors.emplace_back(1.,0,0,0.7);
//        colors.emplace_back(.6,0,0,0.7);
//        colors.emplace_back(0.,1.,0,0.7);
//        colors.emplace_back(0,0.6,0,0.7);
//        colors.emplace_back(0,0,1.,0.7);
//        colors.emplace_back(0,0,.6,0.7);
//
//        for(auto vhf_it = _mesh.vhf_iter(vh_or); vhf_it.valid(); ++vhf_it) {
//            CH ch_or = _mesh.incident_cell(*vhf_it);
//            HFH hf_org = dp.original_halfface_handle(*vhf_it);
//            CH ch_org = mesh_.incident_cell(hf_org);
//
//            int axis = -1;
//            if(cell_he_axis.find(ch_org) != cell_he_axis.end()) {
//                axis = cell_he_axis[ch_org];
//                _mesh_color[ch_or.idx()] = colors[axis];
//            } else
//                _mesh_color[ch_or.idx()] = OVM::Vec4f(1,1,1,0);
//        }
//
//        //nm axis
//        for(auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it) {
//            for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it) {
//                HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
//                if (feature_fprop_[mesh_.face_handle(hfopp)] > 0) {
//                    int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it],
//                                                                                mesh_.normal(hfopp)).second;
//                    cell_nm_axis[*vc_it] = nm_axis;
//                }
//            }
//        }
//
//
//        //find anti parallel singular edges
//        set_only_feature_face_singular_edges_as_targets(_find_featureface_sges);
//        set_same_sector_singular_edges_as_targets(_same_sector_sges);
//        std::set<HEH> he_pairs;
//        auto hfhs = ssaf_.find_anti_parallel_singular_edge_with_disk_path(_vh, or_chs, cell_he_axis, cell_nm_axis, he_pairs);
//
//        if(hfhs.empty())
//            return false;
//
//        std::set<FH> disk_fhs;
//        for(const auto hfi : hfhs)
//            disk_fhs.insert(mesh_.face_handle(hfi));
//
//        std::set<EH> sg_ehs;
//        for(const auto hei : he_pairs)
//            sg_ehs.insert(mesh_.edge_handle(hei));
//
//        return detaching_singularity(_vh, disk_fhs, sg_ehs);
//    }


template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::detach_singular_arc_from_singular_node_wrt_boundary(const VH _vh,
                                                                                         const bool _is_anti_prl,
                                                                                         const bool _start_from_bdy)
{
  SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                          feature_edge_, sgl_vt_, feature_edge_vertex_);

  std::set<VH> vhs_set{_vh};
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
  {
    if (valence_[*ve_it] != 0)
    {
      sge_.compute_edge_valence(*ve_it);
      vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
      vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
    }
  }
  //update singular vertex property
  for (const auto &vhi: vhs_set)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);

  //get boundary halffaces of sector a, b, c ... (separated by feature arcs)
  std::vector<HEH> ft_hehs, bdy_hehs, itr_hehs;
  get_special_halfedges_at_vertex(mesh_, feature_edge_, valence_, _vh, ft_hehs, bdy_hehs, itr_hehs);

  std::set<HEH> fs_hehs;
  fs_hehs.insert(ft_hehs.begin(), ft_hehs.end());
  fs_hehs.insert(bdy_hehs.begin(), bdy_hehs.end());

  //results in boundary zipper node
  if (!_is_anti_prl)
  {
    if (itr_hehs.size() == 1 && bdy_hehs.size() == 2)
    {
      int n_bdy_val1 = 0;
      for (const auto hehi: bdy_hehs)
        if (valence_[mesh_.edge_handle(hehi)] == 1)
          n_bdy_val1++;
      if (n_bdy_val1 < 2)
        return false;
    }
  }

  std::set<CH> or_chs;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    or_chs.insert(*vc_it);

  //get normal axis in boundary cell
  std::map<CH, int> cell_nm_axis;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    int val = valence_[mesh_.edge_handle(*voh_it)];
    if (val == -1 || val == 1)
    {
      auto hehf_it = mesh_.hehf_iter(*voh_it);
      CH ch_s(-1);
      if (mesh_.is_boundary(*voh_it))
        ch_s = mesh_.incident_cell(*hehf_it);
      else
        ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));

      if (!ch_s.is_valid())
      {
        std::cerr << "Error: invalid start cell at edge " << mesh_.edge_handle(*voh_it) << std::endl;
        return false;
      }

      int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, *voh_it, *hehf_it);
      if (axis < 0 || axis > 5)
      {
        auto eh = mesh_.edge_handle(*voh_it);
        std::cout << "Error: singular edge: " << eh << "of valence " << valence_[eh] << " has " << axis << " axis"
                  << std::endl;
        return false;
      }

      //halfedge axis
      if (val == 1)
        axis = negate(axis);

      cell_nm_axis[ch_s] = axis;

      for (; hehf_it.valid(); ++hehf_it)
      {
        auto chi = mesh_.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          if (chi == ch_s)
            continue;

          axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
          cell_nm_axis[chi] = axis;
        }
      }
    }
  }

  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
      if (feature_fprop_[mesh_.face_handle(hfopp)] > 0)
      {
        int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it],
                                                         mesh_.normal(hfopp)).second;
        cell_nm_axis[*vc_it] = nm_axis;
      }
    }
  }


  //find anti parallel singular edges
  std::set<HEH> he_pairs;
  auto hfhs = ssaf_.find_boundary_with_disk_path(_vh, or_chs, fs_hehs, cell_nm_axis, he_pairs, _is_anti_prl,
                                                 _start_from_bdy);

  if (hfhs.empty())
    return false;

  std::set<FH> disk_fhs;
  for (const auto hfi: hfhs)
    disk_fhs.insert(mesh_.face_handle(hfi));

  std::set<EH> sg_ehs;
  for (const auto hei: he_pairs)
    sg_ehs.insert(mesh_.edge_handle(hei));

  return detaching_singularity(_vh, disk_fhs, sg_ehs);
}


template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::detach_singular_arcs_of_zipper_node(const VH _vh,
                                                                         const bool _find_featureface_sges,
                                                                         const bool _same_sector_sges)
{
  SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                          feature_edge_, sgl_vt_, feature_edge_vertex_);

  std::set<VH> vhs_set{_vh};
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
  {
    if (valence_[*ve_it] != 0)
    {
      sge_.compute_edge_valence(*ve_it);
      vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
      vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
    }
  }
  //update singular vertex property
  for (const auto &vhi: vhs_set)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);

  std::set<CH> or_chs;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    or_chs.insert(*vc_it);


  //compute halfedge axis per cell
  std::map<CH, int> cell_he_axis, cell_nm_axis;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    int val = valence_[mesh_.edge_handle(*voh_it)];
    if (val == -1 || val == 1 || val == std::numeric_limits<int>::lowest())
    {
      auto hehf_it = mesh_.hehf_iter(*voh_it);
      CH ch_s(-1);
      if (mesh_.is_boundary(*voh_it))
        ch_s = mesh_.incident_cell(*hehf_it);
      else
        ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));

      if (!ch_s.is_valid())
      {
        std::cerr << "Error: invalid start cell at edge " << mesh_.edge_handle(*voh_it) << std::endl;
        return false;
      }

      int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, *voh_it,
                                                                    *hehf_it);
      if (axis < 0 || axis > 5)
      {
        auto eh = mesh_.edge_handle(*voh_it);
        std::cout << "Error: singular edge: " << eh << "of valence " << valence_[eh] << " has " << axis
                  << " axis"
                  << std::endl;
        return false;
      }

      //halfedge axis
      if (val == 1)
        axis = negate(axis);

      cell_he_axis[ch_s] = axis;

      for (; hehf_it.valid(); ++hehf_it)
      {
        auto chi = mesh_.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          if (chi == ch_s)
            continue;

          axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
          cell_he_axis[chi] = axis;
        }
      }
    }
  }

  //feature face aligned axis
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
      if (feature_fprop_[mesh_.face_handle(*chf_it)] > 0)
      {
        int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], mesh_.normal(hfopp)).second;
        cell_nm_axis[*vc_it] = nm_axis;
      }
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "\nDetaching zipper node at vh " << _vh << " ... finding feature face edge "
                               << _find_featureface_sges << " same region " << _same_sector_sges << std::endl;)


  //find anti parallel singular edges
  set_only_feature_face_singular_edges_as_targets(_find_featureface_sges);
  set_same_sector_singular_edges_as_targets(_same_sector_sges);
  std::set<HEH> he_pairs;
  auto hfhs = ssaf_.find_zipper_node_edge_pair_with_disk_path(_vh, or_chs, cell_he_axis, cell_nm_axis, he_pairs);

  if (hfhs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Empty surface" << std::endl;)
    return false;
  }

  std::set<FH> disk_fhs;
  for (const auto hfi: hfhs)
    disk_fhs.insert(mesh_.face_handle(hfi));

  std::set<EH> sg_ehs;
  for (const auto hei: he_pairs)
    sg_ehs.insert(mesh_.edge_handle(hei));

  return detaching_singularity(_vh, disk_fhs, sg_ehs);
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::detach_singular_arcs_of_zipper_node(const VH _vh, const HEH _heh_s,
                                                                         const std::vector<HEH> &_excl_hehs,
                                                                         const bool _find_featureface_sges,
                                                                         const bool _same_sector_sges,
                                                                         const bool _split)
{
  if (mesh_.halfedge(_heh_s).from_vertex() != _vh)
  {
    std::cerr << "Error: from vertex of halfedge " << _heh_s << " is not " << _vh << std::endl;
    return false;
  }

  if (_split)
    SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                            feature_edge_, sgl_vt_, feature_edge_vertex_);
  std::set<VH> vhs_set{_vh};
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
  {
    if (valence_[*ve_it] != 0)
    {
      sge_.compute_edge_valence(*ve_it);
      vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
      vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
    }
  }
  //update singular vertex property
  for (const auto &vhi: vhs_set)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);

  std::set<CH> or_chs;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    or_chs.insert(*vc_it);

  //compute halfedge axis per cell
  std::map<CH, int> cell_he_axis, cell_nm_axis;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    int val = valence_[mesh_.edge_handle(*voh_it)];
    if (val == -1 || val == 1 || val == std::numeric_limits<int>::lowest())
    {
      auto hehf_it = mesh_.hehf_iter(*voh_it);
      CH ch_s(-1);
      if (mesh_.is_boundary(*voh_it))
        ch_s = mesh_.incident_cell(*hehf_it);
      else
        ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));

      if (!ch_s.is_valid())
      {
        std::cerr << "Error: invalid start cell at edge " << mesh_.edge_handle(*voh_it) << std::endl;
        return false;
      }

      int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, *voh_it,
                                                                    *hehf_it);
      if (axis < 0 || axis > 5)
      {
        auto eh = mesh_.edge_handle(*voh_it);
        std::cout << "Error: singular edge: " << eh << "of valence " << valence_[eh] << " has " << axis
                  << " axis"
                  << std::endl;
        return false;
      }

      //halfedge axis
      if (val == 1)
        axis = negate(axis);

      cell_he_axis[ch_s] = axis;

      for (; hehf_it.valid(); ++hehf_it)
      {
        auto chi = mesh_.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          if (chi == ch_s)
            continue;

          axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
          cell_he_axis[chi] = axis;
        }
      }
    }
  }

  //feature face aligned axis
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
      if (feature_fprop_[mesh_.face_handle(*chf_it)] > 0)
      {
        int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], mesh_.normal(hfopp)).second;
        cell_nm_axis[*vc_it] = nm_axis;
      }
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "\nDetaching zipper node at vh " << _vh << " eh_s " << mesh_.edge_handle(_heh_s)
                               << " ... finding feature face edge " << _find_featureface_sges << " same region "
                               << _same_sector_sges << std::endl;)

  //find anti parallel singular edges
  set_only_feature_face_singular_edges_as_targets(_find_featureface_sges);
  set_same_sector_singular_edges_as_targets(_same_sector_sges);
  std::set<HEH> he_pairs;
  auto hfhs = ssaf_.find_zipper_node_edge_pair_with_disk_path(_vh, _heh_s, _excl_hehs, or_chs, cell_he_axis,
                                                              cell_nm_axis, he_pairs);

  if (hfhs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Empty surface" << std::endl;)
    return false;
  }

  std::set<FH> disk_fhs;
  for (const auto hfi: hfhs)
    disk_fhs.insert(mesh_.face_handle(hfi));

  std::set<EH> sg_ehs;
  for (const auto hei: he_pairs)
    sg_ehs.insert(mesh_.edge_handle(hei));

  return detaching_singularity(_vh, disk_fhs, sg_ehs);
}

//    template <class MeshT>
//    bool
//    DetachAtInvalidSingularNodeT<MeshT>::detach_singular_arcs_of_zipper_node(const VH _vh, MeshT& _mesh, std::vector<OVM::Vec4f>& _mesh_color, const bool _find_featureface_sges, const bool _same_sector_sges) {
//        SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_, feature_edge_, sgl_vt_, feature_edge_vertex_);
//
//        std::set<VH> vhs_set{_vh};
//        for(auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it) {
//            if(valence_[*ve_it] != 0) {
//                sge_.compute_edge_valence(*ve_it);
//                vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
//                vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
//            }
//        }
//        //update singular vertex property
//        for (const auto &vhi : vhs_set)
//            MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);
//
//
//        std::set<CH> or_chs;
//        for(auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
//            or_chs.insert(*vc_it);
//
//        //compute halfedge axis per cell
//        std::map<CH, int> cell_he_axis, cell_nm_axis;
//        for(auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it) {
//            auto ve = mesh_.edge_handle(*voh_it);
//            int val = valence_[ve];
//            if(val == -1 || val == 1 || val == std::numeric_limits<int>::lowest()) {
//                auto hehf_it = mesh_.hehf_iter(*voh_it);
//                CH ch_s(-1);
//                if(mesh_.is_boundary(*voh_it))
//                    ch_s = mesh_.incident_cell(*hehf_it);
//                else
//                    ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
//
//                if(!ch_s.is_valid()) {
//                    std::cerr<<"Error: invalid start cell at edge "<<mesh_.edge_handle(*voh_it)<<std::endl;
//                    return false;
//                }
//
//
//                int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_, *voh_it, *hehf_it);
//                if (axis < 0 || axis > 5) {
//                    auto eh = mesh_.edge_handle(*voh_it);
//                    std::cout << "Error: singular edge: "<< eh<< "of valence " << valence_[eh] << " has "<<axis<<" axis"
//                              << std::endl;
//                    return false;
//                }
//
//                //halfedge axis
//                if(val == 1)
//                    axis = negate(axis);
//
//                cell_he_axis[ch_s] = axis;
//
//                for (; hehf_it.valid(); ++hehf_it) {
//                    auto chi = mesh_.incident_cell(*hehf_it);
//                    if(chi.is_valid()) {
//                        if(chi == ch_s)
//                            continue;
//
//                        axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
//                        cell_he_axis[chi] = axis;
//                    }
//                }
//            }
//        }
//
//        //feature face aligned axis
//        for(auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it) {
//            for(auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it) {
//                HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
//                if(feature_fprop_[mesh_.face_handle(*chf_it)] > 0) {
//                    int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], mesh_.normal(hfopp)).second;
//                    cell_nm_axis[*vc_it] = nm_axis;
//                }
//            }
//        }
//
//        DuplicateOneRingMeshT<MeshT> dp(mesh_, _mesh);
//        dp.copy_one_ring(std::vector<VH>{_vh});
//        _mesh_color.resize(_mesh.n_cells());
//
//        VH vh_or = dp.onering_mesh_vertex_handle(_vh);
//
//        std::vector<OVM::Vec4f> colors;
//        colors.emplace_back(1.,0,0,0.7);
//        colors.emplace_back(.6,0,0,0.7);
//        colors.emplace_back(0.,1.,0,0.7);
//        colors.emplace_back(0,0.6,0,0.7);
//        colors.emplace_back(0,0,1.,0.7);
//        colors.emplace_back(0,0,.6,0.7);
//
//        for(auto vhf_it = _mesh.vhf_iter(vh_or); vhf_it.valid(); ++vhf_it) {
//            CH ch_or = _mesh.incident_cell(*vhf_it);
//            HFH hf_org = dp.original_halfface_handle(*vhf_it);
//            CH ch_org = mesh_.incident_cell(hf_org);
//
//            int axis = -1;
//            if(cell_he_axis.find(ch_org) != cell_he_axis.end()) {
//                axis = cell_he_axis[ch_org];
//                _mesh_color[ch_or.idx()] = colors[axis];
//            } else
//                _mesh_color[ch_or.idx()] = OVM::Vec4f(1,1,1,0);
//        }
//
//
//        //find anti parallel singular edges
//        set_only_feature_face_singular_edges_as_targets(_find_featureface_sges);
//        set_same_sector_singular_edges_as_targets(_same_sector_sges);
//        std::set<HEH> he_pairs;
//        auto hfhs = ssaf_.find_zipper_node_edge_pair_with_disk_path(_vh, or_chs, cell_he_axis, cell_nm_axis, he_pairs);
//
//        if(hfhs.empty())
//            return false;
//
//        std::set<FH> disk_fhs;
//        for(const auto hfi : hfhs)
//            disk_fhs.insert(mesh_.face_handle(hfi));
//
//        std::set<EH> sg_ehs;
//        for(const auto hei : he_pairs)
//            sg_ehs.insert(mesh_.edge_handle(hei));
//
//        return detaching_singularity(_vh, disk_fhs, sg_ehs);
//    }

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::detach_singualr_arc_from_singular_node_wrt_feature_face(const VH _vh,
                                                                                             const bool _is_anti_prl)
{
  SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                          feature_edge_, sgl_vt_, feature_edge_vertex_);

  std::set<CH> or_chs;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    or_chs.insert(*vc_it);

  std::set<VH> vhs_set{_vh};
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
  {
    if (valence_[*ve_it] != 0)
    {
      sge_.compute_edge_valence(*ve_it);
      vhs_set.insert(mesh_.edge(*ve_it).from_vertex());
      vhs_set.insert(mesh_.edge(*ve_it).to_vertex());
    }
  }
  //update singular vertex property
  for (const auto &vhi: vhs_set)
    MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);

  //get normal axis in boundary cell
  std::map<CH, int> cell_nm_axis;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    int val = valence_[mesh_.edge_handle(*voh_it)];
    if (val == -1 || val == 1)
    {
      auto hehf_it = mesh_.hehf_iter(*voh_it);
      CH ch_s(-1);
      if (mesh_.is_boundary(*voh_it))
        ch_s = mesh_.incident_cell(*hehf_it);
      else
        ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));

      if (!ch_s.is_valid())
      {
        std::cerr << "Error: invalid start cell at edge " << mesh_.edge_handle(*voh_it) << std::endl;
        return false;
      }

      int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                    valence_, *voh_it, *hehf_it);
      if (axis < 0 || axis > 5)
      {
        auto eh = mesh_.edge_handle(*voh_it);
        std::cout << "Error: singular edge: " << eh << "of valence " << valence_[eh] << " has " << axis << " axis"
                  << std::endl;
        return false;
      }

      //halfedge axis
      if (val == 1)
        axis = negate(axis);

      cell_nm_axis[ch_s] = axis;

      for (; hehf_it.valid(); ++hehf_it)
      {
        auto chi = mesh_.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          if (chi == ch_s)
            continue;

          axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
          cell_nm_axis[chi] = axis;
        }
      }
    }
  }

  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
      if (feature_fprop_[mesh_.face_handle(hfopp)] > 0)
      {
        int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it],
                                                         mesh_.normal(hfopp)).second;
        cell_nm_axis[*vc_it] = nm_axis;
      }
    }
  }


  //find anti parallel singular edges
  std::set<HEH> he_pairs;
  auto hfhs = ssaf_.find_disk_to_feature_face(_vh, or_chs, cell_nm_axis, he_pairs, _is_anti_prl);

  if (hfhs.empty())
    return false;

  std::set<FH> disk_fhs;
  for (const auto hfi: hfhs)
    disk_fhs.insert(mesh_.face_handle(hfi));

  std::set<EH> sg_ehs;
  for (const auto hei: he_pairs)
    sg_ehs.insert(mesh_.edge_handle(hei));

  //store halfface transitions
  std::map<HFH, int> orig_hf_trans;
  for (const auto fhi: disk_fhs)
  {
    auto hfh0 = mesh_.halfface_handle(fhi, 0);
    orig_hf_trans[hfh0] = trans_prop_[hfh0];
    auto hfh1 = mesh_.halfface_handle(fhi, 1);
    orig_hf_trans[hfh1] = trans_prop_[hfh1];
  }

  bool suc = detaching_singularity(_vh, disk_fhs, sg_ehs);

  if (!suc)
    return false;

  //detach other side
  ALGOHEX_DEBUG_ONLY(std::cerr << "detach other side..." << std::endl;)
  HEH sghe_new(-1);
  for (const auto hei: he_pairs)
    if (std::abs(valence_[mesh_.edge_handle(hei)]) == 1)
      sghe_new = hei;

  //
  std::set<EH> cf_ehs;
  for (const auto fhi: disk_fhs)
  {
    for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
    {
      cf_ehs.insert(*fe_it);
    }
  }

  if (sghe_new.is_valid())
  {
    SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                            feature_edge_, sgl_vt_, feature_edge_vertex_, cf_ehs);

    suc = detach_singular_arc_from_singular_node(_vh, sghe_new, std::vector<HEH>{}, false, true, false);

    if (!suc)
      suc = detach_singular_arcs_of_zipper_node(_vh, sghe_new, std::vector<HEH>{}, false, true, false);

    if (suc)
      return true;
  }
  else
    suc = false;

  if (!suc)
  {
    //restore matching
    for (auto &[hfhi, trans]: orig_hf_trans)
    {
      if (mesh_.is_deleted(hfhi))
      {
        std::cerr << "Error: face was split. Screwed!" << std::endl;
        return false;
      }
      trans_prop_[hfhi] = trans;
    }

    //
    for (const auto &eh: cf_ehs)
      sge_.compute_edge_valence(eh);

    std::set<VH> cf_vhs;
    for (const auto &fhi: disk_fhs)
      for (auto fv_it = mesh_.fv_iter(fhi); fv_it.valid(); ++fv_it)
        cf_vhs.insert(*fv_it);

    //update singular vertex property
    for (const auto &vhi: cf_vhs)
      MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);

    //update quaternions
    std::set<CH> cf_chs;
    for (const auto &vhi: cf_vhs)
      for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
        cf_chs.insert(*vc_it);

    QuaternionSmoothing::optimize_quaternions_in_cells(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                       feature_fprop_, feature_edge_, cf_chs);
  }

  return false;
}


template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::
detaching_singularity(const VH _vh, const std::set<FH> &cut_fhs, const std::set<EH> &_bound_ehs)
{
  //store halfface transitions
  std::map<HFH, int> orig_hf_trans;
  for (const auto fhi: cut_fhs)
  {
    auto hfh0 = mesh_.halfface_handle(fhi, 0);
    orig_hf_trans[hfh0] = trans_prop_[hfh0];
    auto hfh1 = mesh_.halfface_handle(fhi, 1);
    orig_hf_trans[hfh1] = trans_prop_[hfh1];
  }

  auto cf_ehs = ArcZippingT<MeshT>::zipping_onering(_vh, cut_fhs, _bound_ehs);

  if (!cf_ehs.empty())
  {
    //check if complex singular edges appear
    for (const auto &eh: cf_ehs)
    {
      int val = sge_.calculate_edge_valence(eh);
      if ((val > 1 || val < -1) && val != std::numeric_limits<int>::lowest())
      {
        std::cerr << "Warining: complex edge appeared! eh: " << eh << " valence " << val << " restore matchings!"
                  << std::endl;

        for (const auto &fhi: cut_fhs)
        {
          auto hfh0 = mesh_.halfface_handle(fhi, 0);
          trans_prop_[hfh0] = orig_hf_trans[hfh0];
          auto hfh1 = mesh_.halfface_handle(fhi, 1);
          trans_prop_[hfh1] = orig_hf_trans[hfh1];
        }

        return false;
      }
    }


    std::set<VH> cf_vhs;
    for (const auto &fhi: cut_fhs)
      for (auto fv_it = mesh_.fv_iter(fhi); fv_it.valid(); ++fv_it)
        cf_vhs.insert(*fv_it);


    //update quaternions
    std::set<CH> cf_chs;

    for (const auto &vhi: cf_vhs)
      for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
        cf_chs.insert(*vc_it);

    for (const auto &eh: cf_ehs)
      sge_.compute_edge_valence(eh);

    QuaternionSmoothing::optimize_quaternions_in_cells(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                       feature_fprop_, feature_edge_, cf_chs);

    for (const auto &vhi: cf_vhs)
      MeshPropertiesT<MeshT>::update_singular_vertex_property(vhi);

    return true;
  }
  else
  {
    for (const auto &fhi: cut_fhs)
    {
      auto hfh0 = mesh_.halfface_handle(fhi, 0);
      trans_prop_[hfh0] = orig_hf_trans[hfh0];
      auto hfh1 = mesh_.halfface_handle(fhi, 1);
      trans_prop_[hfh1] = orig_hf_trans[hfh1];
    }
  }


  return false;
}

template<class MeshT>
int
DetachAtInvalidSingularNodeT<MeshT>::n_incident_singular_edges(const VH _vh) const
{
  int n_sg = 0;
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (valence_[*ve_it] != 0)
      n_sg++;

  return n_sg;
}

//TODO: the circle cannot bound any other singularities
template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::is_singular_circle_noise(const HEH _heh) const
{
  if (valence_[mesh_.edge_handle(_heh)] != 0)
  {
    VH vhf = mesh_.halfedge(_heh).from_vertex();
    auto or_chs = get_onering_cells(mesh_, vhf);
    std::set<HFH> ft_hfhs;
    auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, vhf, or_chs, *mesh_.hec_iter(_heh), ft_hfhs);
    int n_nffe = n_incident_non_ffe_singular_edges_in_same_region(mesh_, feature_face_edge_, valence_, vhf, por_chs);
    if (n_nffe <= 2)
      return false;

    auto sg_hehs = get_halfedges_on_singular_arc_with_intertior_surface(mesh_, valence_, feature_face_vertex_, _heh);

    if (sg_hehs.size() > 20)
      return false;

    //it's a circle
    if (mesh_.halfedge(sg_hehs[0]).from_vertex() == mesh_.halfedge(sg_hehs.back()).to_vertex())
    {
      //there is zipper node on circle
      bool has_tp = false;
      for (auto i = 1u; i < sg_hehs.size(); ++i)
        if (MeshPropertiesT<MeshT>::node_index(mesh_.halfedge(sg_hehs[i]).from_vertex()) == 10)
        {
          has_tp = true;
          break;
        }

      if (has_tp)
      {
        //if there is a feature face vertex
        for (auto i = 1u; i < sg_hehs.size(); ++i)
          if (feature_face_vertex_[mesh_.halfedge(sg_hehs[i]).from_vertex()])
          {
            return false;
          }
        //otherwise
        return true;
      }
    }
  }

  return false;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::is_locally_meshable(const VH _vh)
{
  SplitHelperT<MeshT>::split_one_ring_all(mesh_, es_, _vh, valence_, feature_face_edge_, feature_face_vertex_,
                                          feature_edge_, sgl_vt_, feature_edge_vertex_);

  bool ilm = opc_.is_locally_meshable(_vh, true);
//if (!ilm)
//            ilm = lmcwf_.is_locally_meshable(_vh, true);

//        std::cerr<<"is meshable? "<<ilm<<std::endl;

  return ilm;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::is_locally_meshable_simplified(const VH _vh, int _n_sges)
{
  bool is_meshable = false;
  int nid = MeshPropertiesT<MeshT>::node_index(_vh);
  if ((_n_sges <= 4 && (nid != 20 && nid != 10)) || _n_sges > 4)
  {
    is_meshable = is_locally_meshable(_vh);
  }

  return is_meshable;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::
is_detachable_arc(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_nff_val_u, int _n_int) const
{
  if (sgl_vt_[_vh] != 0)
  {
    if (feature_face_vertex_[_vh])
    {
      if (_n_int_val_1 + _n_nff_val_u >= 2 || _n_int_val1 + _n_nff_val_u >= 2)
        return true;
    }
    else if (_n_int > 2 && _n_int_val_1 + _n_nff_val_u >= 2 && _n_int_val1 + _n_nff_val_u >= 2)
      return true;
  }

  return false;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::
is_detachable_arc_to_ffe(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_ffe_val_1, int _n_ffe_val1,
                         int _n_ffe_val_u) const
{
  if (sgl_vt_[_vh] != 0 && feature_face_vertex_[_vh])
  {
    if ((_n_int_val_1 >= 1 && _n_ffe_val_1 + _n_ffe_val_u >= 1) ||
        (_n_int_val1 >= 1 && _n_ffe_val1 + _n_ffe_val_u >= 1))
      return true;
  }
  return false;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::
is_detachale_arc_wrt_boundary(int _n_int_val_1, int _n_int_val1, int _n_int, int _n_bdy, int _n_ft) const
{
  if (((_n_int_val_1 + _n_int_val1) >= 1 && _n_ft >= 2) ||
      _n_int_val_1 + _n_int_val1 >= 2)
  {
    //no bdy sg arc can be kept
    if (_n_bdy == 1)
      return false;

    return true;
  }

  return false;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::
is_detachale_arc_wrt_boundary_from_boundary_singular_edge(int _n_int_val_1, int _n_int_val1, int _n_int, int _n_bdy,
                                                          int _n_ft) const
{
  if (_n_ft >= 4 && _n_bdy >= 1)
  {

    return true;
  }

  return false;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::
is_detachable_arc_wrt_feature_face(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_int) const
{
  if (feature_edge_vertex_[_vh] && _n_int >= 2 && (_n_int_val_1 >= 1 || _n_int_val1 >= 1))
  {

    return true;
  }

  return false;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::
is_detachable_zipper_node(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_int) const
{
  if (sgl_vt_[_vh] != 0)
  {
    if (feature_face_vertex_[_vh])
    {
      if (_n_int_val_1 >= 1 && _n_int_val1 >= 1)
        return true;
    }
    else if (_n_int > 2 && _n_int_val_1 >= 1 && _n_int_val1 >= 1)
      return true;
  }
  return false;
}

template<class MeshT>
bool
DetachAtInvalidSingularNodeT<MeshT>::
is_detachable_zipper_node_to_ffe(const VH _vh, int _n_int_val_1, int _n_int_val1, int _n_ffe_val_1, int _n_ffe_val1,
                                 int _n_ffe_val_u) const
{
  if (sgl_vt_[_vh] != 0 && feature_face_vertex_[_vh])
  {
    if ((_n_int_val_1 >= 1 && _n_ffe_val1 + _n_ffe_val_u >= 1) ||
        (_n_int_val1 >= 1 && _n_ffe_val_1 + _n_ffe_val_u >= 1))
      return true;
  }
  return false;
}

}

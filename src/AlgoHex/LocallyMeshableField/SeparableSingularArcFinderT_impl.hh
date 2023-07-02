/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define SEPARABLESINGULARARCFINDERT_C

#include "SeparableSingularArcFinderT.hh"
#include "CellInfo.hh"

namespace AlgoHex
{
template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_anti_parallel_singular_edge_with_disk_path(const VH _vh,
                                                                                    const std::set<CH> &_or_chs,
                                                                                    std::map<CH, int> &_cell_he_axis,
                                                                                    std::map<CH, int> &_cell_nm_axis,
                                                                                    std::set<HEH> &_he_pairs)
{
  std::set<HEH> all_cd_hehs;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto voe = mesh_.edge_handle(*voh_it);
    if (valence_[voe] == -1 || valence_[voe] == 1 || valence_[voe] == std::numeric_limits<int>::lowest())
    {
      all_cd_hehs.insert(*voh_it);
    }
  }

  //get candidate edge set
  auto candidate_hehs = all_cd_hehs;
  if (find_feature_face_sges_only_)
  {
    //has kept edge but not singular edge
    bool force_detach = false;
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto voe = mesh_.edge_handle(*voh_it);
      if (keep_on_feature_face_[voe] && valence_[voe] == 0)
      {
        force_detach = true;
        break;
      }
    }

    for (const auto hehi: all_cd_hehs)
    {
      auto ehi = mesh_.edge_handle(hehi);
      if (valence_[ehi] == -1 || valence_[ehi] == 1 || valence_[ehi] == std::numeric_limits<int>::lowest())
      {
        if (!feature_face_edge_[ehi])
          candidate_hehs.erase(hehi);

        if (!force_detach && keep_on_feature_face_[ehi])
          candidate_hehs.erase(hehi);
      }
    }
  }
  else
  {
    for (const auto hehi: all_cd_hehs)
      if (feature_face_edge_[mesh_.edge_handle(hehi)])
        candidate_hehs.erase(hehi);
  }

  if (!candidate_hehs.empty())
  {
    //find dual path to
    auto s_hehs = sort_starting_halfedges_for_antiparallel(_vh);
    ALGOHEX_DEBUG_ONLY(std::cerr << "sorted ehs: ";
                               for (auto hehi: s_hehs)
                                 std::cerr << " " << mesh_.edge_handle(hehi);
                               std::cerr << std::endl;)

    for (const auto hehi: s_hehs)
    {
      auto ve = mesh_.edge_handle(hehi);
      int val = valence_[ve];
      if (val == -1 || val == 1)
      {
        if (!feature_face_edge_[ve] && all_cd_hehs.find(hehi) != all_cd_hehs.end())
        {
          auto real_cd_hehs = get_candidate_halfedges(candidate_hehs, hehi, val, true);

          ALGOHEX_DEBUG_ONLY(std::cerr << "candidate ehs: ";
                                     for (auto hehi: real_cd_hehs)
                                       std::cerr << " " << mesh_.edge_handle(hehi);
                                     std::cerr << std::endl;)

          if (!real_cd_hehs.empty())
          {
            double angle_threshold = get_angle_threshold(_vh, hehi, _or_chs);

            auto path = find_anti_parallel_singular_edge_with_disk_path(_vh, _or_chs, hehi,
                                                                        real_cd_hehs, _cell_he_axis, _cell_nm_axis,
                                                                        angle_threshold, _he_pairs);
            if (!_he_pairs.empty())
              return path;
          }
        }
      }
    }
  }


  return std::set<HFH>{};
}

template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_anti_parallel_singular_edge_with_disk_path(const VH _vh, const HEH _heh_s,
                                                                                    const std::vector<HEH> &_excl_hehs,
                                                                                    const std::set<CH> &_or_chs,
                                                                                    std::map<CH, int> &_cell_he_axis,
                                                                                    std::map<CH, int> &_cell_nm_axis,
                                                                                    std::set<HEH> &_he_pairs)
{
  std::set<HEH> all_cd_hehs;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto voe = mesh_.edge_handle(*voh_it);
    if (valence_[voe] == -1 || valence_[voe] == 1 || valence_[voe] == std::numeric_limits<int>::lowest())
    {
      all_cd_hehs.insert(*voh_it);
    }
  }

  for (const auto hehi: _excl_hehs)
    all_cd_hehs.erase(hehi);

  //get candidate edge set
  auto candidate_hehs = all_cd_hehs;
  if (find_feature_face_sges_only_)
  {
    //has kept edge but not singular edge
    bool force_detach = false;
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto voe = mesh_.edge_handle(*voh_it);
      if (keep_on_feature_face_[voe] && valence_[voe] == 0)
      {
        force_detach = true;
        break;
      }
    }

    for (const auto hehi: all_cd_hehs)
    {
      auto ehi = mesh_.edge_handle(hehi);
      if (valence_[ehi] == -1 || valence_[ehi] == 1 || valence_[ehi] == std::numeric_limits<int>::lowest())
      {
        if (!feature_face_edge_[ehi])
          candidate_hehs.erase(hehi);

        if (!force_detach && keep_on_feature_face_[ehi])
          candidate_hehs.erase(hehi);
      }
    }
  }
  else
  {
    for (const auto hehi: all_cd_hehs)
      if (feature_face_edge_[mesh_.edge_handle(hehi)])
        candidate_hehs.erase(hehi);
  }

  //get candidate edge set
  auto real_cd_hehs = get_candidate_halfedges(candidate_hehs, _heh_s, valence_[mesh_.edge_handle(_heh_s)], true);

  ALGOHEX_DEBUG_ONLY(std::cerr << "candidate ehs: ";
                             for (auto hehi: real_cd_hehs)
                               std::cerr << " " << mesh_.edge_handle(hehi);
                             std::cerr << std::endl;)

  if (!real_cd_hehs.empty())
  {
    double angle_threshold = get_angle_threshold(_vh, _heh_s, _or_chs);

    auto path = find_anti_parallel_singular_edge_with_disk_path(_vh, _or_chs, _heh_s,
                                                                real_cd_hehs, _cell_he_axis, _cell_nm_axis,
                                                                angle_threshold, _he_pairs);
    if (!_he_pairs.empty())
      return path;
  }


  return std::set<HFH>{};
}


template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_anti_parallel_singular_edge_with_disk_path(const VH _vh,
                                                                                    const std::set<CH> &_or_chs,
                                                                                    const HEH _heh_s,
                                                                                    const std::set<HEH> &_candidate_hehs,
                                                                                    std::map<CH, int> &_cell_he_axis,
                                                                                    std::map<CH, int> &_cell_nm_axis,
                                                                                    const double _angle_threshold,
                                                                                    std::set<HEH> &_he_pair)
{
  for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
  {
    int target_axis = _cell_he_axis[*hec_it];
    target_axis = negate(target_axis);

    std::vector<CH> path;
    ALGOHEX_DEBUG_ONLY(
            std::cerr << "find anti parallel to edge " << mesh_.edge_handle(_heh_s) << " ch: " << *hec_it << std::endl;)
//                    HEH other_sghe = find_anti_parallel_singular_edge(*voh_it, *hec_it, _cell_he_axis, target_axis, _cut_fhs, _or_chs, path);

    HEH other_sghe = find_anti_parallel_singular_edge_with_dual_path(_heh_s, *hec_it,
                                                                     _cell_he_axis, _cell_nm_axis, target_axis, _or_chs,
                                                                     _candidate_hehs, _angle_threshold, path);

    if (other_sghe.is_valid())
    {
      _he_pair.insert(_heh_s);
      _he_pair.insert(other_sghe);

      //process path
      std::vector<HEH> v_sghehs{_heh_s, other_sghe};

      auto path1 = get_handle_free_dual_path(path, v_sghehs);
//                ALGOHEX_DEBUG_ONLY(std::cerr<<"dpath after handle free: ";
//                for (const auto ch : path1)
//                    std::cerr<<" "<<ch;
//                std::cerr<<std::endl;)

      auto path0 = get_singular_vertex_free_path(_vh, path1, v_sghehs);
//                ALGOHEX_DEBUG_ONLY(std::cerr<<"dpath after split: ";
//                for (const auto ch : path0)
//                    std::cerr<<" "<<ch;
//                std::cerr<<std::endl;)


      return bounded_surface(path0, _he_pair);

    }
  }
  return std::set<HFH>{};
}

template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_boundary_with_disk_path(const VH _vh, const std::set<CH> &_or_chs,
                                                                 const std::set<HEH> &_fs_hehs,
                                                                 std::map<CH, int> &_cell_nm_axis,
                                                                 std::set<HEH> &_he_pairs, const bool _is_anti_prl,
                                                                 const bool _start_from_bdy)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Detach " << _is_anti_prl << " parallel wrt bdy at vh " << _vh << std::endl;)

  //find dual path to boundary
  std::vector<HEH> sghehs;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto voe = mesh_.edge_handle(*voh_it);
    if (valence_[voe] == -1 || valence_[voe] == 1)
      sghehs.push_back(*voh_it);
  }
  for (const auto hehi: sghehs)
  {
    if (!mesh_.is_deleted(hehi))
    {
      auto is_bdy = mesh_.is_boundary(hehi);
      if (!is_bdy && !_start_from_bdy)
      {
        //exclude boundary feature sectors if the normal and sge direction are orthogonal
        std::set<CH> candidate_cells = get_candidate_cells(_vh, hehi, _fs_hehs, _is_anti_prl);

        ALGOHEX_DEBUG_ONLY(std::cerr << "candidate chs ";
                                   for (auto chi: candidate_cells)
                                     std::cerr << " " << chi;
                                   std::cerr << std::endl;)
        if (!candidate_cells.empty())
        {
          auto path = find_boundary_with_disk_path(_vh, _or_chs, hehi, candidate_cells, _cell_nm_axis,
                                                   _he_pairs, _is_anti_prl);

          if (!_he_pairs.empty())
            return path;
        }
      }
      else if (is_bdy && _start_from_bdy)
      {
        ALGOHEX_DEBUG_ONLY(
                std::cerr << "Finding path to bdy cell from bdy sge " << mesh_.edge_handle(hehi) << std::endl;)
        //update search direction
        auto dir = mesh_.vertex(mesh_.halfedge(hehi).to_vertex()) - mesh_.vertex(mesh_.halfedge(hehi).from_vertex());
        for (auto hec_it = mesh_.hec_iter(hehi); hec_it.valid(); ++hec_it)
        {
          int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*hec_it], dir.normalize()).second;
          _cell_nm_axis[*hec_it] = nm_axis;
        }

        std::set<CH> candidate_cells = get_candidate_cells(_vh, hehi, _fs_hehs, _is_anti_prl);

        for (auto hec_it = mesh_.hec_iter(hehi); hec_it.valid(); ++hec_it)
          candidate_cells.erase(*hec_it);

        std::set<HFH> path;
        if (!candidate_cells.empty())
        {
          path = find_boundary_with_disk_path(_vh, _or_chs, hehi, candidate_cells, _cell_nm_axis,
                                              _he_pairs, _is_anti_prl);

          if (!_he_pairs.empty())
            return path;
        }

        if (path.empty())
        {//restore
          for (auto hec_it = mesh_.hec_iter(hehi); hec_it.valid(); ++hec_it)
            _cell_nm_axis.erase(*hec_it);

          for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
          {
            if (mesh_.is_boundary(*vc_it))
            {
              for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
              {
                HFH hfopp = mesh_.opposite_halfface_handle(*chf_it);
                if (mesh_.is_boundary(hfopp))
                {
                  int nm_axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it],
                                                                   mesh_.normal(hfopp)).second;
                  _cell_nm_axis[*vc_it] = nm_axis;
                }
              }
            }
          }
        }

      }
    }
  }


  return std::set<HFH>{};
}

template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_boundary_with_disk_path(const VH _vh, const std::set<CH> &_or_chs,
                                                                 const HEH _heh_s, const std::set<CH> &_candidate_chs,
                                                                 std::map<CH, int> &_cell_nm_axis,
                                                                 std::set<HEH> &_he_pair, const bool _is_anti_prl)
{
  for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
  {
    //incoming halfedge direction
    int target_axis = _cell_nm_axis[*hec_it];
    target_axis = negate(target_axis);


    std::vector<CH> path;
    ALGOHEX_DEBUG_ONLY(
            std::cerr << "find " << _is_anti_prl << " parallel boundary to edge " << mesh_.edge_handle(_heh_s)
                      << " valence " << valence_[mesh_.edge_handle(_heh_s)] << " ch: " << *hec_it << std::endl;)
    CH end_ch = find_boundary_with_dual_path(_heh_s, *hec_it,
                                             _cell_nm_axis, target_axis, _or_chs, _candidate_chs, path, _is_anti_prl);

    if (end_ch.is_valid())
    {
      _he_pair.insert(_heh_s);

      //process path
      std::vector<HEH> v_sghehs{_heh_s};

      auto path0 = get_handle_free_dual_path(path, v_sghehs);

      auto path1 = get_singular_vertex_free_path(_vh, path0, v_sghehs, true);


      return bounded_surface(path1, _he_pair);
    }
  }
  return std::set<HFH>{};
}


template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_zipper_node_edge_pair_with_disk_path(const VH _vh, const std::set<CH> &_or_chs,
                                                                              std::map<CH, int> &_cell_he_axis,
                                                                              std::map<CH, int> &_cell_nm_axis,
                                                                              std::set<HEH> &_he_pairs)
{
  std::set<HEH> candidate_hehs;
  //get candidate edge set
  if (find_feature_face_sges_only_)
  {
    //has kept edge but not singular edge
    bool force_detach = false;
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto voe = mesh_.edge_handle(*voh_it);
      if (keep_on_feature_face_[voe] && valence_[voe] == 0)
      {
        force_detach = true;
        break;
      }
    }
    ALGOHEX_DEBUG_ONLY(std::cerr << "force detach " << force_detach << std::endl;)
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto voe = mesh_.edge_handle(*voh_it);
      if ((valence_[voe] == -1 || valence_[voe] == 1 || valence_[voe] == std::numeric_limits<int>::lowest())
          && feature_face_edge_[voe] && !(!force_detach && keep_on_feature_face_[voe]))
        candidate_hehs.insert(*voh_it);
    }
  }
  else
  {
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto voe = mesh_.edge_handle(*voh_it);
      if ((valence_[voe] == -1 || valence_[voe] == 1 || valence_[voe] == std::numeric_limits<int>::lowest()) &&
          !feature_face_edge_[voe])
        candidate_hehs.insert(*voh_it);
    }
  }

  if (!candidate_hehs.empty())
  {
    //find dual path to
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto ve = mesh_.edge_handle(*voh_it);
      int val = valence_[ve];
      if (val == -1 || val == 1)
      {
        if (!feature_face_edge_[ve])
        {
          ALGOHEX_DEBUG_ONLY(std::cerr << "start eh " << ve << " ";)
          auto real_cd_hehs = get_candidate_halfedges(candidate_hehs, *voh_it, val, false);
          if (!real_cd_hehs.empty())
          {
            ALGOHEX_DEBUG_ONLY(std::cerr << "real candidate ehs: ";
                                       for (auto hehi: real_cd_hehs)
                                         std::cerr << " " << mesh_.edge_handle(hehi);
                                       std::cerr << std::endl;)

            double angle_threshold_tp = get_angle_threshold_zipper_node(_vh, *voh_it, _or_chs);

            auto path = find_zipper_node_edge_pair_with_disk_path(_vh, _or_chs, *voh_it,
                                                                  real_cd_hehs, _cell_he_axis, _cell_nm_axis,
                                                                  angle_threshold_tp, _he_pairs);

            if (!_he_pairs.empty())
              return path;
          }
        }
      }
    }
  }


  return std::set<HFH>{};
}

template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_zipper_node_edge_pair_with_disk_path(const VH _vh, const HEH _heh_s,
                                                                              const std::vector<HEH> &_excl_hehs,
                                                                              const std::set<CH> &_or_chs,
                                                                              std::map<CH, int> &_cell_he_axis,
                                                                              std::map<CH, int> &_cell_nm_axis,
                                                                              std::set<HEH> &_he_pairs)
{
  std::set<HEH> candidate_hehs;
  //get candidate edge set
  if (find_feature_face_sges_only_)
  {
    //has kept edge but not singular edge
    bool force_detach = false;
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto voe = mesh_.edge_handle(*voh_it);
      if (keep_on_feature_face_[voe] && valence_[voe] == 0)
      {
        force_detach = true;
        break;
      }
    }


    ALGOHEX_DEBUG_ONLY(std::cerr << "force detach " << force_detach << std::endl;)
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto voe = mesh_.edge_handle(*voh_it);
      if ((valence_[voe] == -1 || valence_[voe] == 1 || valence_[voe] == std::numeric_limits<int>::lowest())
          && feature_face_edge_[voe] && !(!force_detach && keep_on_feature_face_[voe]))
        candidate_hehs.insert(*voh_it);
    }
  }
  else
  {
    for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
    {
      auto voe = mesh_.edge_handle(*voh_it);
      if ((valence_[voe] == -1 || valence_[voe] == 1 || valence_[voe] == std::numeric_limits<int>::lowest()) &&
          !feature_face_edge_[voe])
        candidate_hehs.insert(*voh_it);
    }
  }

  //exclude
  for (const auto hehi: _excl_hehs)
    candidate_hehs.erase(hehi);

  auto real_cd_hehs = get_candidate_halfedges(candidate_hehs, _heh_s, valence_[mesh_.edge_handle(_heh_s)], false);
  if (!real_cd_hehs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "real candidate ehs: ";
                               for (auto hehi: real_cd_hehs)
                                 std::cerr << " " << mesh_.edge_handle(hehi);
                               std::cerr << std::endl;)

    double angle_threshold_tp = get_angle_threshold_zipper_node(_vh, _heh_s, _or_chs);

    auto path = find_zipper_node_edge_pair_with_disk_path(_vh, _or_chs, _heh_s,
                                                          real_cd_hehs, _cell_he_axis, _cell_nm_axis,
                                                          angle_threshold_tp, _he_pairs);

    if (!_he_pairs.empty())
      return path;
  }


  return std::set<HFH>{};
}

template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_zipper_node_edge_pair_with_disk_path(const VH _vh, const std::set<CH> &_or_chs,
                                                                              const HEH _heh_s,
                                                                              const std::set<HEH> &_candidate_hehs,
                                                                              std::map<CH, int> &_cell_he_axis,
                                                                              std::map<CH, int> &_cell_nm_axis,
                                                                              const double _angle_threshold_tp,
                                                                              std::set<HEH> &_he_pair)
{
  for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
  {
    int target_axis = _cell_he_axis[*hec_it];
    std::vector<CH> path;
    ALGOHEX_DEBUG_ONLY(
            std::cerr << "find other zipper node edge to edge " << mesh_.edge_handle(_heh_s) << " ch: " << *hec_it
                      << std::endl;)

    HEH other_sghe = find_zipper_node_edge_pair_with_dual_path(_heh_s, *hec_it,
                                                               _cell_he_axis, _cell_nm_axis, target_axis, _or_chs,
                                                               _candidate_hehs, _angle_threshold_tp, path);

    if (other_sghe.is_valid())
    {
      _he_pair.insert(_heh_s);
      _he_pair.insert(other_sghe);

      //process path
      std::vector<HEH> v_sghehs{_heh_s, other_sghe};

      auto path0 = get_handle_free_dual_path(path, v_sghehs);

      auto path1 = get_singular_vertex_free_path(_vh, path0, v_sghehs);

      ALGOHEX_DEBUG_ONLY(std::cerr << "dpath after split: ";
                                 for (const auto ch: path1)
                                   std::cerr << " " << ch;
                                 std::cerr << std::endl;)


      return bounded_surface(path1, _he_pair);
    }
  }
  return std::set<HFH>{};
}


template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_disk_to_feature_face(const VH _vh, const std::set<CH> &_or_chs,
                                                              std::map<CH, int> &_cell_nm_axis,
                                                              std::set<HEH> &_he_pairs, const bool _is_anti_prl)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Extend " << _is_anti_prl << " parallel wrt feature face at vh " << _vh << std::endl;)

  //find dual path to feature face
  std::vector<HEH> sghehs;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto voe = mesh_.edge_handle(*voh_it);
    if (!feature_face_edge_[voe] && (valence_[voe] == -1 || valence_[voe] == 1))
      sghehs.push_back(*voh_it);
  }
  for (const auto hehi: sghehs)
  {
    auto path = find_disk_to_feature_face(_vh, _or_chs, hehi, _cell_nm_axis,
                                          _he_pairs, _is_anti_prl);

    if (!_he_pairs.empty())
      return path;
  }


  return std::set<HFH>{};
}


template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::find_disk_to_feature_face(const VH _vh, const std::set<CH> &_or_chs,
                                                              const HEH _heh_s, std::map<CH, int> &_cell_nm_axis,
                                                              std::set<HEH> &_he_pair, const bool _is_anti_prl)
{
  for (auto hec_it = mesh_.hec_iter(_heh_s); hec_it.valid(); ++hec_it)
  {
    //incoming halfedge direction
    int target_axis = _cell_nm_axis[*hec_it];
    target_axis = negate(target_axis);


    std::vector<CH> path;
    ALGOHEX_DEBUG_ONLY(
            std::cerr << "find " << _is_anti_prl << " parallel to feature face " << mesh_.edge_handle(_heh_s) << " ch: "
                      << *hec_it << std::endl;)
    CH end_ch = find_dual_path_to_feature_face(_heh_s, *hec_it,
                                               _cell_nm_axis, target_axis, _or_chs, path, _is_anti_prl);

    if (end_ch.is_valid())
    {
      _he_pair.insert(_heh_s);

      //process path
      std::vector<HEH> v_sghehs{_heh_s};

      auto path0 = get_handle_free_dual_path(path, v_sghehs);

      auto path1 = get_singular_vertex_free_path(_vh, path0, v_sghehs, true);

      auto bd_hfhs = bounded_surface(path1, _he_pair, true);
      ALGOHEX_DEBUG_ONLY(std::cerr << "bd fhs ";
                                 for (auto hfh: bd_hfhs)
                                   std::cerr << " " << mesh_.face_handle(hfh);
                                 std::cerr << std::endl;)
      //get feature face edge
      HEH ff_he(-1);
      HFH bd_opp_hfh(-1), ff_hf(-1);
      for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
      {
        auto ve = mesh_.edge_handle(*voh_it);
        if (feature_face_edge_[ve])
        {
          for (auto hehf_it = mesh_.hehf_iter(*voh_it); hehf_it.valid(); ++hehf_it)
          {
            if (bd_hfhs.find(mesh_.opposite_halfface_handle(*hehf_it)) != bd_hfhs.end())
            {
              ff_he = *voh_it;
              bd_opp_hfh = *hehf_it;
              break;
            }
          }
        }

        if (ff_he.is_valid())
          break;
      }

      if (ff_he.is_valid())
      {
        HFH hf_it = bd_opp_hfh;
        do
        {
          auto hfh_adj = mesh_.adjacent_halfface_in_cell(hf_it, ff_he);
          if (!hfh_adj.is_valid())
          {
            std::cerr << "Error: adjacent halfface is invalid! he_s: " << ff_he << " hfh_s " << bd_opp_hfh << std::endl;
            break;
          }
          hf_it = mesh_.opposite_halfface_handle(hfh_adj);
          if (feature_fprop_[mesh_.face_handle(hf_it)] > 0)
          {
            ff_hf = hf_it;
            break;
          }
        }
        while (hf_it != bd_opp_hfh);

        if (ff_hf.is_valid())
        {
          //add halfface to set
          HFH hf_other = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(ff_hf, ff_he));
          bd_hfhs.insert(hf_other);

          //he
          HEH he_other = mesh_.opposite_halfedge_handle(mesh_.prev_halfedge_in_halfface(ff_he, hf_other));
          _he_pair.insert(he_other);

          return bd_hfhs;
        }
      }

    }
  }
  return std::set<HFH>{};
}

//===================================================================================================//
template<class MeshT>
EH
SeparableSingularArcFinderT<MeshT>::opposite_edge_in_cell(const EH _eh, const CH _ch) const
{
  VH vhf = mesh_.edge(_eh).from_vertex();
  VH vht = mesh_.edge(_eh).to_vertex();

  std::vector<VH> other_vhs;
  other_vhs.reserve(2);
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
    if (*cv_it != vhf && *cv_it != vht)
      other_vhs.push_back(*cv_it);

  return mesh_.edge_handle(mesh_.halfedge(other_vhs[0], other_vhs[1]));
}

template<class MeshT>
std::set<HFH>
SeparableSingularArcFinderT<MeshT>::
bounded_surface(const std::vector<CH> &_cells, const std::set<HEH> &_sg_hehs, const bool _excld_feature_face_edge) const
{
  if (_cells.empty())
    return std::set<HFH>{};

  std::set<HFH> surface;

  std::set<HFH> all_bound_hfs;
  std::set<HFH> all_hfs;

  for (const auto ch: _cells)
    for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
    {
      all_hfs.insert(mesh_.opposite_halfface_handle(*chf_it));
    }

  std::set<CH> dpath;
  dpath.insert(_cells.begin(), _cells.end());

  //store one ring cells
  VH vh_f = mesh_.halfedge(*(_sg_hehs.begin())).from_vertex();
  std::set<CH> or_cells;
  for (auto vc_it = mesh_.vc_iter(vh_f); vc_it.valid(); ++vc_it)
    or_cells.insert(*vc_it);

  //bound halffaces (no boundary halfface)
  for (const auto hf: all_hfs)
  {
    auto chi = mesh_.incident_cell(hf);
    if (dpath.find(chi) == dpath.end() && or_cells.find(chi) != or_cells.end())
      all_bound_hfs.insert(hf);
  }

  if (_sg_hehs.empty())
    return all_bound_hfs;

  ALGOHEX_DEBUG_ONLY(std::cout << "one ring cells bdy fhs: ";
                             for (const auto hf: all_bound_hfs)
                               std::cout << " " << mesh_.face_handle(hf);
                             std::cout << std::endl;

                             std::cout << "bounded ehs : ";
                             for (const auto e: _sg_hehs)
                               std::cout << " " << mesh_.edge_handle(e);
                             std::cout << std::endl;)

  //start halfface
  HFH hfh(-1);
  for (auto hehf_it = mesh_.hehf_iter(*(_sg_hehs.begin())); hehf_it.valid(); ++hehf_it)
    if (all_bound_hfs.find(*hehf_it) != all_bound_hfs.end())
    {
      hfh = *hehf_it;
      break;
    }

  //bound edges
  std::set<EH> sgehs;
  std::queue<HEH> que;
  if (hfh.is_valid())
  {
    surface.insert(hfh);

    for (const auto hehi: _sg_hehs)
      sgehs.insert(mesh_.edge_handle(hehi));

    for (auto hfhe_it = mesh_.hfhe_iter(hfh); hfhe_it.valid(); ++hfhe_it)
      que.push(mesh_.opposite_halfedge_handle(*hfhe_it));
  }
//std::cerr<<"surface edges: ";
  while (!que.empty())
  {
    auto he_cur = que.front();
    que.pop();

    //touch the boundary edges of the surface
    EH eh_cur = mesh_.edge_handle(he_cur);
    if ((_excld_feature_face_edge && feature_face_edge_[eh_cur]) || mesh_.is_boundary(eh_cur) ||
        (sgehs.find(eh_cur) != sgehs.end()))
      continue;
//std::cerr<<mesh_.edge_handle(he_cur)<<" ";

    HFH hf_uk;
    int n_cut_fhs = n_unvisited_halffaces_at_manifold(all_bound_hfs, surface, he_cur, hf_uk);

    if (n_cut_fhs == 1)
    {
      surface.insert(hf_uk);

      for (auto hfhe_it = mesh_.hfhe_iter(hf_uk); hfhe_it.valid(); ++hfhe_it)
      {
        que.push(mesh_.opposite_halfedge_handle(*hfhe_it));
      }
    }
  }


  return surface;
}


template<class MeshT>
int
SeparableSingularArcFinderT<MeshT>::
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
std::vector<CH>
SeparableSingularArcFinderT<MeshT>::
get_singular_vertex_free_path(const VH _vh, const std::vector<CH> &_dpath, const std::vector<HEH> &_sghehs,
                              bool _to_boundary)
{
  if (_dpath.empty())
    return std::vector<CH>{};

  int n_split = 0, ns = 0;

  std::set<VH> excld_vhs;
  for (const auto heh: _sghehs)
  {
    excld_vhs.insert(mesh_.halfedge(heh).from_vertex());
    excld_vhs.insert(mesh_.halfedge(heh).to_vertex());
  }

  //excluded cells from splitting
  std::set<CH> excl_chs;

  //exclude feature face edges to be split (can be split if feature edge exist)
  std::set<FH> excl_fhs;

  //exclude feature face edges to be split
  std::set<EH> excl_ehs;

  //initialize excluded entities
  for (auto i = 0u; i < _dpath.size() - 1; ++i)
  {
    auto fhi = common_face(mesh_, _dpath[i], _dpath[i + 1]);
    if (fhi.is_valid())
    {
      if (feature_fprop_[fhi] > 0)
      {
        excl_fhs.insert(fhi);

        excl_chs.insert(_dpath[i]);
        excl_chs.insert(_dpath[i + 1]);
      }
    }
  }
  if (_to_boundary)
  {
    excl_chs.insert(_dpath.back());

    for (auto cf_it = mesh_.cf_iter(_dpath.back()); cf_it.is_valid(); ++cf_it)
    {
      if (feature_fprop_[*cf_it] > 0)
      {
        excl_fhs.insert(*cf_it);
      }
    }
  }

  //store candidate cells
  std::set<CH> dp_cells;
  dp_cells.insert(_dpath.begin(), _dpath.end());

  CH ch_s = _dpath.front();
  CH ch_t = _dpath.back();
  ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " ch_t: " << ch_t << std::endl;)

  //store target cells
  std::set<CH> target_chs{ch_t}, start_chs{ch_s};

  //0. split feature cells
  CellSplitT<MeshT> cs(mesh_);
  ALGOHEX_DEBUG_ONLY(
          std::cerr << "sgeh0 " << mesh_.edge_handle(_sghehs[0]) << " sgeh1 " << mesh_.edge_handle(_sghehs.back())
                    << std::endl;)
  if (is_feature_cell(ch_s) && excl_chs.find(ch_s) == excl_chs.end())
  {
    dp_cells.erase(ch_s);
    start_chs.erase(ch_s);

    auto vh_m = cs.cell_split(ch_s);
    n_split++;
    std::set<CH> inc_chs;
    for (auto vc_it = mesh_.vc_iter(vh_m); vc_it.valid(); ++vc_it)
    {
      if (!is_feature_cell(*vc_it))
      {
        inc_chs.insert(*vc_it);
        dp_cells.insert(*vc_it);
      }
    }

    for (auto hec_it = mesh_.hec_iter(_sghehs.front()); hec_it.valid(); ++hec_it)
    {
      if (inc_chs.find(*hec_it) != inc_chs.end())
      {
        ch_s = *hec_it;
        start_chs.insert(ch_s);
        break;
      }
    }
  }
  if (!_to_boundary)
  {
    if (is_feature_cell(ch_t) && excl_chs.find(ch_t) == excl_chs.end())
    {
      dp_cells.erase(ch_t);
      target_chs.erase(ch_t);
      auto vh_m = cs.cell_split(ch_t);
      n_split++;

      std::set<CH> inc_chs;
      for (auto vc_it = mesh_.vc_iter(vh_m); vc_it.valid(); ++vc_it)
      {
        if (!is_feature_cell(*vc_it))
        {
          inc_chs.insert(*vc_it);
          dp_cells.insert(*vc_it);
        }
      }

      for (auto hec_it = mesh_.hec_iter(_sghehs.back()); hec_it.valid(); ++hec_it)
      {
        if (inc_chs.find(*hec_it) != inc_chs.end())
        {
          ch_t = *hec_it;
          target_chs.insert(ch_t);
          break;
        }
      }
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " start or end feature cells " << std::endl;)
  ns = n_split;

  ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " ch_t: " << ch_t << std::endl;)
  if (!ch_s.is_valid() || !ch_t.is_valid())
    return std::vector<CH>{};

  excl_chs.insert(ch_s);
  excl_chs.insert(ch_t);

  FaceSplitT<MeshT> fs(mesh_);

  //w.r.t feature faces
  //1. split cells which contain feature face(s)
//        std::cerr<<"split cell: ";
  std::vector<CH> split_chs;
  for (const auto chi: dp_cells)
  {
    if (is_feature_cell(chi) && excl_chs.find(chi) == excl_chs.end())
    {
//                std::cerr << " " << chi;
      split_chs.push_back(chi);
    }
  }
//        std::cerr<<std::endl;
  for (const auto chi: split_chs)
  {
    if (cs.is_split_ok(chi))
    {
      dp_cells.erase(chi);

      n_split++;
      auto vh_new = cs.cell_split(chi);

      for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
        if (!is_feature_cell(*vc_it))
          dp_cells.insert(*vc_it);
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split << " feature cells " << std::endl;)
  ns = n_split;

  //2. split faces that have feature/singular edges
  std::set<EH> all_ehs;
  for (const auto ch: dp_cells)
  {
    for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
      if (feature_edge_[*ce_it] > 0 || valence_[*ce_it] != 0)
        all_ehs.insert(*ce_it);
  }

  //remove singular edges
  for (const auto hehi: _sghehs)
    all_ehs.erase(mesh_.edge_handle(hehi));

  std::set<FH> ff_eh_fhs;
  for (const auto eh: all_ehs)
  {
    for (auto ef_it = mesh_.ef_iter(eh); ef_it.valid(); ++ef_it)
    {
      auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(*ef_it, 0));
      auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(*ef_it, 1));

      //incident to cells on path
      if (dp_cells.find(ch0) != dp_cells.end() || dp_cells.find(ch1) != dp_cells.end())
      {
        ff_eh_fhs.insert(*ef_it);
      }
    }
  }

  n_split += SplitHelperT<MeshT>::split_faces_with_cell_property_update(mesh_, fs, dp_cells, ff_eh_fhs, start_chs,
                                                                        target_chs, excl_fhs);

  //remove cells that have feature edges
  for (const auto eh: all_ehs)
  {
    for (auto ec_it = mesh_.ec_iter(eh); ec_it.valid(); ++ec_it)
    {
      dp_cells.erase(*ec_it);
      start_chs.erase(*ec_it);
      target_chs.erase(*ec_it);
    }

    //remove excluded faces
    for (auto ef_it = mesh_.ef_iter(eh); ef_it.valid(); ++ef_it)
    {
      if (feature_fprop_[*ef_it] > 0)
        excl_fhs.erase(*ef_it);
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split - ns << " faces with f edge" << std::endl;)
  ns = n_split;


  //3. split faces that have feature face edges
  //excluded edges
  for (const auto fhi: excl_fhs)
  {
    for (auto fe_it = mesh_.fe_iter(fhi); fe_it.valid(); ++fe_it)
      if (feature_edge_[*fe_it] == 0 && valence_[*fe_it] == 0)
        excl_ehs.insert(*fe_it);
  }
  for (const auto hehi: _sghehs)
    excl_ehs.insert(mesh_.edge_handle(hehi));

  all_ehs.clear();
  for (const auto ch: dp_cells)
  {
    for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
      if (feature_face_edge_[*ce_it])
        all_ehs.insert(*ce_it);
  }

  //remove excluded edges
  for (const auto ehi: excl_ehs)
    all_ehs.erase(ehi);


  //split faces to get enough DOF
  ff_eh_fhs.clear();
  for (const auto eh: all_ehs)
  {
    for (auto ef_it = mesh_.ef_iter(eh); ef_it.valid(); ++ef_it)
    {
      //incident cells
      auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(*ef_it, 0));
      auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(*ef_it, 1));
      //any incident cell is in candidate cell set
      if (dp_cells.find(ch0) != dp_cells.end() || dp_cells.find(ch1) != dp_cells.end())
      {
        //no excluded face
        ff_eh_fhs.insert(*ef_it);
      }
    }
  }

  n_split += SplitHelperT<MeshT>::split_faces_with_cell_property_update(mesh_, fs, dp_cells, ff_eh_fhs, start_chs,
                                                                        target_chs, excl_fhs);


  //remove the cells that have feature face edges
  for (const auto eh: all_ehs)
  {
    for (auto ec_it = mesh_.ec_iter(eh); ec_it.valid(); ++ec_it)
    {
      if (excl_chs.find(*ec_it) == excl_chs.end())
      {
        dp_cells.erase(*ec_it);
        start_chs.erase(*ec_it);
        target_chs.erase(*ec_it);
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "split " << n_split - ns << " faces with ff edge" << std::endl;)
  ns = n_split;


  //4. split edges which has feature edge vertices or singular vertices
  std::set<EH> split_ehs;
  std::set<VH> bad_vhs;
  for (const auto chi: dp_cells)
  {
    auto bad_ehs = bad_edges_in_cell2(chi, excld_vhs);
    split_ehs.insert(bad_ehs.begin(), bad_ehs.end());
  }
  for (const auto ehi: split_ehs)
  {
    auto vh00 = mesh_.edge(ehi).from_vertex();
    auto vh11 = mesh_.edge(ehi).to_vertex();
    if ((sgl_vt_[vh00] != 0 || feature_edge_vertex_[vh00]) && excld_vhs.find(vh00) == excld_vhs.end())
      bad_vhs.insert(vh00);
    if ((sgl_vt_[vh11] != 0 || feature_edge_vertex_[vh11]) && excld_vhs.find(vh11) == excld_vhs.end())
      bad_vhs.insert(vh11);
  }
//        std::cerr<<"split ehs; ";
//        for(auto ehi : split_ehs)
//            std::cerr<<" "<<ehi;
//        std::cerr<<std::endl;
  n_split += SplitHelperT<MeshT>::split_edges_with_cell_property_update(mesh_, es_, dp_cells, split_ehs, start_chs,
                                                                        target_chs, excl_ehs);

  //remove the cells that incident to vertices
  for (const auto vhi: bad_vhs)
  {
    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
    {
      dp_cells.erase(*vc_it);
      start_chs.erase(*vc_it);
      target_chs.erase(*vc_it);
    }

    //remove excluded edges incident to bad vhs
    for (auto ve_it = mesh_.ve_iter(vhi); ve_it.valid(); ++ve_it)
    {
      if (feature_face_edge_[*ve_it] && excl_ehs.find(*ve_it) != excl_ehs.end())
        excl_ehs.erase(*ve_it);
    }
  }
  ALGOHEX_DEBUG_ONLY(
          std::cerr << "Split " << n_split - ns << " edges incident to other special(fe) vertices" << std::endl;)
  ns = n_split;

  //5. split edges which has feature face vertices or singular vertices

  //exclude vertices on excluded edges
  for (const auto ehi: excl_ehs)
  {
    excld_vhs.insert(mesh_.edge(ehi).from_vertex());
    excld_vhs.insert(mesh_.edge(ehi).to_vertex());
  }

  split_ehs.clear();
  for (const auto chi: dp_cells)
  {
    auto bad_ehs = bad_edges_in_cell(chi, excld_vhs);
    split_ehs.insert(bad_ehs.begin(), bad_ehs.end());
  }

  bad_vhs.clear();
  for (const auto ehi: split_ehs)
  {
    auto vh00 = mesh_.edge(ehi).from_vertex();
    auto vh11 = mesh_.edge(ehi).to_vertex();
    if ((sgl_vt_[vh00] != 0 || feature_face_vertex_[vh00]) && excld_vhs.find(vh00) == excld_vhs.end())
      bad_vhs.insert(vh00);
    if ((sgl_vt_[vh11] != 0 || feature_face_vertex_[vh11]) && excld_vhs.find(vh11) == excld_vhs.end())
      bad_vhs.insert(vh11);
  }
//        std::cerr<<"split ehs: ";
//for(auto ehi : split_ehs)
//    std::cerr<<" "<<ehi;
//std::cerr<<std::endl;

  n_split += SplitHelperT<MeshT>::split_edges_with_cell_property_update(mesh_, es_, dp_cells, split_ehs, start_chs,
                                                                        target_chs, excl_ehs);

  //remove the cells incident to vertices
  for (const auto vhi: bad_vhs)
  {
    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
    {
      dp_cells.erase(*vc_it);
      start_chs.erase(*vc_it);
      target_chs.erase(*vc_it);
    }
  }
  ALGOHEX_DEBUG_ONLY(
          std::cerr << "Split " << n_split - ns << " edges incident to other special(ffe) vertices" << std::endl;)

  if (n_split > 0)
  {
    ch_s.reset();
    for (auto hec_it = mesh_.hec_iter(_sghehs[0]); hec_it.valid(); ++hec_it)
      if (start_chs.find(*hec_it) != start_chs.end() && !is_feature_cell(*hec_it))
      {
        ch_s = *hec_it;
        break;
      }
    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " target chs ";)

    //get target cells
    std::set<CH> real_target_chs;
    if (_to_boundary)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "to boundary: ";)
      for (const auto chi: target_chs)
      {
        if (is_feature_cell(chi) && n_bad_vertices_in_cell(chi, excld_vhs) == 0)
        {
          real_target_chs.insert(chi);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << chi;)
        }
      }
      ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
    }
    else
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "to sge: ";)

      for (auto hec_it = mesh_.hec_iter(_sghehs.back()); hec_it.valid(); ++hec_it)
        if (target_chs.find(*hec_it) != target_chs.end() && !is_feature_cell(*hec_it) &&
            n_bad_vertices_in_cell(*hec_it, excld_vhs) == 0)
        {
          real_target_chs.insert(*hec_it);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << *hec_it;)
        }
      ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)

    }


    if (ch_s.is_valid() && !real_target_chs.empty())
      return shortest_dual_path_from_cell_to_target(ch_s, real_target_chs, dp_cells, excld_vhs);
    else
      return std::vector<CH>{};
  }

  return _dpath;
}


template<class MeshT>
std::vector<CH>
SeparableSingularArcFinderT<MeshT>::get_handle_free_dual_path(const std::vector<CH> &_input_path,
                                                              const std::vector<HEH> &_sghehs)
{
  if (_input_path.empty())
  {
    return _input_path;
  }

  std::set<VH> excld_vhs;
  for (const auto heh: _sghehs)
  {
    excld_vhs.insert(mesh_.halfedge(heh).from_vertex());
    excld_vhs.insert(mesh_.halfedge(heh).to_vertex());
  }


  //1. find blocking faces
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
  //detach w.r.t bdy
  bool to_boundary = _sghehs.size() == 1;

  if (!blocking_fhs.empty())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "blocking fhs: ";
                               for (auto fhi: blocking_fhs)
                                 std::cerr << " " << fhi;
                               std::cerr << std::endl;)

    //store target cells
    CH ch_s = _input_path[0];
    CH ch_t = _input_path.back();

    ALGOHEX_DEBUG_ONLY(
            std::cerr << "ch_s: " << ch_s << " ch_t: " << ch_t << " is to bdy: " << to_boundary << std::endl;)

    std::set<CH> target_chs{ch_t};
    std::set<CH> start_chs{ch_s};

    //find cells incident to blocking fhs and cell is not ch_s
    std::set<CH> cells;
    for (const auto chi: _input_path)
      cells.insert(chi);

    auto dp_cells = cells;

    //exclude edges of ch_s
    cells.erase(ch_s);

    //faces to split
    std::set<CH> chs_split;
    for (const auto fhi: blocking_fhs)
    {
      auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 0));
      auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 1));

      if (cells.find(ch0) != cells.end())
        chs_split.insert(ch0);

      if (cells.find(ch1) != cells.end())
        chs_split.insert(ch1);
    }

    //split cells
    int n_split = 0;

    CellSplitT<MeshT> cs(mesh_);
    for (const auto chi: chs_split)
      if (cs.is_split_ok(chi))
      {
        dp_cells.erase(chi);

        bool target_fd = target_chs.find(chi) != target_chs.end();

        n_split++;
        auto vh_new = cs.cell_split(chi);
        for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
          dp_cells.insert(*vc_it);

        //update target cells
        if (target_fd)
        {
          target_chs.erase(chi);

          if (to_boundary)
          {
            for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
              if (is_feature_cell(*vc_it))
                target_chs.insert(*vc_it);
          }
          else
          {
            for (auto vc_it = mesh_.vc_iter(vh_new); vc_it.valid(); ++vc_it)
              target_chs.insert(*vc_it);
          }
        }
      }

    for (const auto fhi: blocking_fhs)
    {
      auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 0));
      auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 1));

      dp_cells.erase(ch0);
      dp_cells.erase(ch1);

      target_chs.erase(ch0);
      target_chs.erase(ch1);
    }

    //find start cell
    for (auto hec_it = mesh_.hec_iter(_sghehs[0]); hec_it.valid(); ++hec_it)
    {
      if (start_chs.find(*hec_it) != start_chs.end() && !is_feature_cell(*hec_it))
      {
        ch_s = *hec_it;
        break;
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "candidate target chs: ";
                               for (const auto chi: target_chs)
                               {
                                 std::cerr << " " << chi;
                               }
                               std::cerr << std::endl;
                               std::cerr << "ch_s: " << ch_s << " target chs: ";)

    //get target cells
    std::set<CH> real_target_chs;
    if (to_boundary)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "bdy, eh: " << mesh_.edge_handle(_sghehs.back()) << std::endl;)
      for (const auto chi: target_chs)
        if (is_feature_cell(chi))
        {
          real_target_chs.insert(chi);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << chi;)
        }


    }
    else
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "interior, eh: " << mesh_.edge_handle(_sghehs.back()) << "   ";)
      for (auto hec_it = mesh_.hec_iter(_sghehs.back()); hec_it.valid(); ++hec_it)
      {
        if (target_chs.find(*hec_it) != target_chs.end())
        {
          real_target_chs.insert(*hec_it);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << *hec_it;)
        }
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)

    if (ch_s.is_valid() && !real_target_chs.empty())
      output_path = shortest_dual_path_from_cell_to_target(ch_s, real_target_chs, dp_cells, excld_vhs, false);
    else
      return std::vector<CH>{};
  }


  //2. blocking edges
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

    ALGOHEX_DEBUG_ONLY(
            std::cerr << "ch_s: " << ch_s << " ch_t: " << ch_t << " is to bdy: " << to_boundary << std::endl;)

    std::set<CH> target_chs{ch_t};
    std::set<CH> start_chs{ch_s};

    //find faces incident to blocking ehs and faces are not in ch_s
    std::set<FH> cell_fhs;
    for (const auto chi: output_path)
      for (auto cf_it = mesh_.cf_iter(chi); cf_it.valid(); ++cf_it)
        cell_fhs.insert(*cf_it);

    //exclude edges of ch_s
    for (auto cf_it = mesh_.cf_iter(ch_s); cf_it.valid(); ++cf_it)
      cell_fhs.erase(*cf_it);

    //faces to split
    std::set<FH> fhs_split;
    for (const auto ehi: blocking_ehs)
    {
      for (auto ef_it = mesh_.ef_iter(ehi); ef_it.valid(); ++ef_it)
        if (cell_fhs.find(*ef_it) != cell_fhs.end())
          fhs_split.insert(*ef_it);
    }
//            std::cerr<<"fhs split : ";
//            for (const auto fhi : fhs_split) {
//                std::cerr<<" "<<fhi;
//            }
//            std::cerr<<std::endl;

    std::set<CH> dp_cells;
    for (const auto chi: output_path)
      dp_cells.insert(chi);
//            std::cerr<<"dp cells : ";
//            for (const auto ch : dp_cells) {
//                    std::cerr<<" "<<ch;
//            }
//            std::cerr<<std::endl;

    //split faces
    FaceSplitT<MeshT> fs(mesh_);
    std::set<FH> excl_fhs;
    int n_split = SplitHelperT<MeshT>::split_faces_with_cell_property_update(mesh_, fs, dp_cells, fhs_split, start_chs,
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
        {
          dp_cells.erase(*ec_it);
          target_chs.erase(*ec_it);
          start_chs.erase(*ec_it);
        }
      }


    //find start cell
    for (auto hec_it = mesh_.hec_iter(_sghehs[0]); hec_it.valid(); ++hec_it)
    {
      if (start_chs.find(*hec_it) != start_chs.end() && !is_feature_cell(*hec_it))
      {
        ch_s = *hec_it;
        break;
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "candidate target chs: ";
                               for (const auto chi: target_chs)
                               {
                                 std::cerr << " " << chi;
                               }
                               std::cerr << std::endl;
                               std::cerr << "ch_s: " << ch_s << " real target chs: ";)

    //get target cells
    std::set<CH> real_target_chs;
    if (to_boundary)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "bdy, eh: " << mesh_.edge_handle(_sghehs.back()) << std::endl;)
      for (const auto chi: target_chs)
        if (is_feature_cell(chi))
        {
          real_target_chs.insert(chi);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << chi;)
        }


    }
    else
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "interior, eh: " << mesh_.edge_handle(_sghehs.back()) << "   ";)
      for (auto hec_it = mesh_.hec_iter(_sghehs.back()); hec_it.valid(); ++hec_it)
      {
        if (target_chs.find(*hec_it) != target_chs.end())
        {
          real_target_chs.insert(*hec_it);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << *hec_it;)
        }
      }
    }


    if (ch_s.is_valid() && !real_target_chs.empty())
      output_path = shortest_dual_path_from_cell_to_target(ch_s, real_target_chs, dp_cells, excld_vhs, false);
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

    ALGOHEX_DEBUG_ONLY(std::cerr << "ch_s: " << ch_s << " ch_t: " << ch_t;)

    std::set<CH> target_chs{ch_t};
    std::set<CH> start_chs{ch_s};

    //find edges incident to blocking vhs and the edges are not in ch_s
    std::set<EH> cell_ehs;
    for (const auto chi: output_path)
      for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
        cell_ehs.insert(*ce_it);

    //exclude edges of ch_s
    for (auto ce_it = mesh_.ce_iter(ch_s); ce_it.valid(); ++ce_it)
      cell_ehs.erase(*ce_it);

    //edges to split
    std::set<EH> ehs1;
    for (const auto vhi: blocking_vhs)
    {
      for (auto ve_it = mesh_.ve_iter(vhi); ve_it.valid(); ++ve_it)
        if (cell_ehs.find(*ve_it) != cell_ehs.end())
          ehs1.insert(*ve_it);
    }

    std::set<CH> dp_cells;
    for (const auto chi: output_path)
      dp_cells.insert(chi);

    //split edges
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
        {
          dp_cells.erase(*vc_it);
          target_chs.erase(*vc_it);
          start_chs.erase(*vc_it);
        }
      }


    //find start cell
    for (auto hec_it = mesh_.hec_iter(_sghehs[0]); hec_it.valid(); ++hec_it)
    {
      if (start_chs.find(*hec_it) != start_chs.end() && !is_feature_cell(*hec_it))
      {
        ch_s = *hec_it;
        break;
      }
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "candidate target chs: ";
                               for (const auto chi: target_chs)
                               {
                                 std::cerr << " " << chi;
                               }
                               std::cerr << std::endl;
                               std::cerr << "ch_s: " << ch_s << " real target chs: ";)

    //get target cells
    std::set<CH> real_target_chs;
    if (to_boundary)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "bdy, eh: " << mesh_.edge_handle(_sghehs.back()) << std::endl;)
      for (const auto chi: target_chs)
        if (is_feature_cell(chi))
        {
          real_target_chs.insert(chi);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << chi;)
        }


    }
    else
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "interior, eh: " << mesh_.edge_handle(_sghehs.back()) << "   ";)
      for (auto hec_it = mesh_.hec_iter(_sghehs.back()); hec_it.valid(); ++hec_it)
      {
        if (target_chs.find(*hec_it) != target_chs.end())
        {
          real_target_chs.insert(*hec_it);
          ALGOHEX_DEBUG_ONLY(std::cerr << " " << *hec_it;)
        }
      }
    }

    if (ch_s.is_valid() && !real_target_chs.empty())
      return shortest_dual_path_from_cell_to_target(ch_s, real_target_chs, dp_cells, excld_vhs, false);
    else
      return std::vector<CH>{};

  }


  return output_path;
}

template<class MeshT>
std::vector<EH>
SeparableSingularArcFinderT<MeshT>::
bad_edges_in_cell(const CH _ch, const std::set<VH> &_excld_sg_vhs) const
{
  std::vector<EH> bad_ehs;
  for (auto ce_it = mesh_.ce_iter(_ch); ce_it.valid(); ++ce_it)
  {
    if (valence_[*ce_it] != 0 || feature_face_edge_[*ce_it])
      continue;

    auto vhf = mesh_.edge(*ce_it).from_vertex();
    auto vht = mesh_.edge(*ce_it).to_vertex();

    if (((sgl_vt_[vhf] != 0 || feature_face_vertex_[vhf]) && _excld_sg_vhs.find(vhf) == _excld_sg_vhs.end()) ||
        ((sgl_vt_[vht] != 0 || feature_face_vertex_[vht]) && _excld_sg_vhs.find(vht) == _excld_sg_vhs.end()))
    {
      bad_ehs.push_back(*ce_it);
    }
  }

  return bad_ehs;
}

template<class MeshT>
std::vector<EH>
SeparableSingularArcFinderT<MeshT>::
bad_edges_in_cell2(const CH _ch, const std::set<VH> &_excld_sg_vhs) const
{
  std::vector<EH> bad_ehs;
  for (auto ce_it = mesh_.ce_iter(_ch); ce_it.valid(); ++ce_it)
  {
    if (valence_[*ce_it] != 0 || feature_edge_[*ce_it] > 0)
      continue;

    auto vhf = mesh_.edge(*ce_it).from_vertex();
    auto vht = mesh_.edge(*ce_it).to_vertex();

    if (((sgl_vt_[vhf] != 0 || feature_edge_vertex_[vhf]) && _excld_sg_vhs.find(vhf) == _excld_sg_vhs.end()) ||
        ((sgl_vt_[vht] != 0 || feature_edge_vertex_[vht]) && _excld_sg_vhs.find(vht) == _excld_sg_vhs.end()))
    {
      bad_ehs.push_back(*ce_it);
    }
  }

  return bad_ehs;
}

template<class MeshT>
int
SeparableSingularArcFinderT<MeshT>::
n_bad_vertices_in_cell(const CH _ch, const std::set<VH> &_excld_sg_vhs) const
{
  int n_bad_vhs = 0;
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
  {
    if (((sgl_vt_[*cv_it] != 0 || feature_edge_vertex_[*cv_it]) && _excld_sg_vhs.find(*cv_it) == _excld_sg_vhs.end()))
    {
      n_bad_vhs++;
    }
  }

  return n_bad_vhs;
}


template<class MeshT>
std::vector<CH>
SeparableSingularArcFinderT<MeshT>::
shortest_dual_path_from_cell_to_target(const CH _ch_s, const std::set<CH> &_target_chs, const std::set<CH> &_dp_cells,
                                       const std::set<VH> &_excld_vhs, const bool _no_sgv) const
{
  using ICH = std::pair<int, CH>;

  std::map<CH, int> min_dist;
  min_dist.insert(std::make_pair(_ch_s, 0));

  std::map<CH, CH> previous;
  previous.insert(std::make_pair(_ch_s, CH(-1)));

  std::vector<CH> chs;
  std::set<ICH> que;
  que.insert(ICH(0, _ch_s));

  CH ch_t(-1);
  while (!que.empty() && !ch_t.is_valid())
  {
    auto dc_cur = *que.begin();
    que.erase(que.begin());

//            std::cerr<<" "<<dc_cur.second;

    auto chfs = mesh_.cell(dc_cur.second).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);
      auto c_next = mesh_.incident_cell(hf_opp);
      if (!c_next.is_valid())
        continue;

      if (_target_chs.find(c_next) != _target_chs.end())
      {
        ALGOHEX_DEBUG_ONLY(std::cout << "found " << c_next << std::endl;)
        min_dist[c_next] = min_dist[dc_cur.second] + 1;
        previous[c_next] = dc_cur.second;
        ch_t = c_next;
        break;
      }

      if (_dp_cells.find(c_next) != _dp_cells.end())
      {
        if (_no_sgv && n_bad_vertices_in_cell(c_next, _excld_vhs) > 0)
          continue;

        double dist = min_dist[dc_cur.second] + 1;
        if (min_dist.find(c_next) == min_dist.end() || dist < min_dist[c_next])
        {
          que.erase(std::make_pair(min_dist[c_next], c_next));

          min_dist[c_next] = dist;
          previous[c_next] = dc_cur.second;

          que.insert(std::make_pair(min_dist[c_next], c_next));
        }
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


    return path;
  }
//        std::cerr<<"empty path"<<std::endl;

  return std::vector<CH>{};
}


template<class MeshT>
HEH
SeparableSingularArcFinderT<MeshT>::find_anti_parallel_singular_edge_with_dual_path(const HEH _heh, const CH _ch_s,
                                                                                    std::map<CH, int> &_cell_he_axis,
                                                                                    std::map<CH, int> &_cell_nm_axis,
                                                                                    const int _target_axis,
                                                                                    const std::set<CH> &_onering_chs,
                                                                                    const std::set<HEH> &_candidate_hehs,
                                                                                    const double _angle_threshold,
                                                                                    std::vector<CH> &_dpath)
{
  //
  VH vhf = mesh_.halfedge(_heh).from_vertex();

  //locally align quaternion field
  //store old quaterions
  std::map<CH, Quaternion> old_cell_quaternions;
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    old_cell_quaternions[*vc_it] = cell_quaternions_[*vc_it];

  //align quaternions to singular or feature edges or boundary faces
  QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                      feature_fprop_, feature_edge_, vhf, true);

  //get blocking cells if at zero feature sectors
  std::set<FH> blk_fhs;
  std::set<CH> blk_cells;
  if (enable_blocking_)
  {
    blk_cells = get_blocking_cells(_heh, blk_fhs);
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "Angle thresold(cos) " << _angle_threshold;
                             std::cerr << " Blocking fhs ";
                             for (auto fhi: blk_fhs)
                               std::cerr << " " << fhi;
                             std::cerr << std::endl;)

  //
  std::map<CellAxis, CellAxis> previous;

  std::map<CellAxis, double> previous_incrs;


  //edge direction in parametrization coordinates
  int u = _target_axis;

  CellAxis ca(_ch_s, u);

  previous[ca] = CellAxis(CH(-1), -1);


  ALGOHEX_DEBUG_ONLY(std::cout << "Floodfilling..." << " _target_axis " << _target_axis << "\n";)

  CellAxis ca_target(CH(-1), -1);
  std::queue<CellAxis> que;
  que.push(ca);
  while (!que.empty() && !ca_target.first.is_valid())
  {
    auto ca_cur = que.front();

    que.pop();

    auto chfs = mesh_.cell(ca_cur.first).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);
      auto fh_next = mesh_.face_handle(hf_opp);
      if (feature_fprop_[fh_next] > 0)
      {
        if (same_sector_sges_only_)
          continue;
        else if (_cell_nm_axis[ca_cur.first] != ca_cur.second)
          continue;
      }

      auto ch_next = mesh_.incident_cell(hf_opp);

      if (!same_sector_sges_only_ && blk_fhs.find(fh_next) != blk_fhs.end())
        continue;

      if (_onering_chs.find(ch_next) == _onering_chs.end())
        continue;

      auto c_next = ch_next.idx();

      bool increase = false;
      double inc_val = 0.;
      auto ca_next = next_cell_axis(ca_cur, hf_opp, inc_val, increase);


      if (inc_val < _angle_threshold)
        continue;


      //visited
      if (previous.find(ca_next) != previous.end())
        continue;

      HEH he_test = is_anti_parallel_found(ca_next, vhf, _heh, _cell_he_axis);
      if (he_test.is_valid())
      {
        ALGOHEX_DEBUG_ONLY(std::cout << "found antiparallel " << mesh_.edge_handle(he_test) << " ";)
        if (_candidate_hehs.find(he_test) != _candidate_hehs.end())
        {
          ca_target = ca_next;
          previous[ca_target] = ca_cur;
          previous_incrs[ca_target] = inc_val;

          ALGOHEX_DEBUG_ONLY(std::cout << "found edge in candidate " << mesh_.edge_handle(he_test) << std::endl;)

          break;
        }
      }

      //push to queue
      //floodfill
      previous[ca_next] = ca_cur;
      previous_incrs[ca_next] = inc_val;

      que.push(ca_next);
    }
  }

  //recover quaternions
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    cell_quaternions_[*vc_it] = old_cell_quaternions[*vc_it];

  //
  HEH he_other(-1);

  if (ca_target.first.is_valid())
  {
    he_other = singular_halfedge_in_cell(vhf, ca_target.first, _heh);

    //find cells
    //get the path
    std::vector<CH> rpath;
    CellAxis ca_it = ca_target;
    while (ca_it.first.is_valid())
    {
      rpath.push_back(ca_it.first);

      ca_it = previous[ca_it];
    }

    //check if there are duplicate elements
    std::set<CH> path(rpath.begin(), rpath.end());
    if (rpath.empty() || path.size() != rpath.size())
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "Empty dpath or has duplicate elements, skip." << std::endl;)
      return HEH(-1);
    }
    else
    {//store
      _dpath.reserve(rpath.size());
      for (auto r_it = rpath.rbegin(); r_it != rpath.rend(); ++r_it)
        _dpath.push_back(*r_it);
    }
  }

  return he_other;
}


template<class MeshT>
CH
SeparableSingularArcFinderT<MeshT>::
find_boundary_with_dual_path(const HEH _heh, const CH _ch_s, std::map<CH, int> &_cell_nm_axis, const int _target_axis,
                             const std::set<CH> &_onering_chs, const std::set<CH> &_candidate_chs,
                             std::vector<CH> &_dpath, const bool _is_anti_prl)
{
  VH vhf = mesh_.halfedge(_heh).from_vertex();

  //locally align quaternion field
  //store old quaterions
  std::map<CH, Quaternion> old_cell_quaternions;
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    old_cell_quaternions[*vc_it] = cell_quaternions_[*vc_it];

  //align quaternions to singular or feature edges or boundary faces
  QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                      feature_fprop_, feature_edge_, vhf, true);

  //
  std::map<CellAxis, CellAxis> previous;
  std::map<CellAxis, double> previous_incrs;

  //edge direction in parametrization coordinates
  int u = _target_axis;

  CellAxis ca(_ch_s, u);

  previous[ca] = CellAxis(CH(-1), -1);

  ALGOHEX_DEBUG_ONLY(std::cout << "Floodfilling..." << " \n";)

  CellAxis ca_target(CH(-1), -1);
  std::queue<CellAxis> que;
  que.push(ca);
  while (!que.empty() && !ca_target.first.is_valid())
  {
    auto ca_cur = que.front();

    que.pop();

//            std::cout<<" ch cur: "<<ca_cur.first<<" ax: "<<ca_cur.second<<std::endl;

    auto chfs = mesh_.cell(ca_cur.first).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);

      auto ch_next = mesh_.incident_cell(hf_opp);

      if (_onering_chs.find(ch_next) == _onering_chs.end())
        continue;

      auto c_next = ch_next.idx();

      bool increase = false;
      double inc_val = 0.;
      auto ca_next = next_cell_axis(ca_cur, hf_opp, inc_val, increase);

      if (_is_anti_prl)
        if (!increase)
          continue;

      //visited
      if (previous.find(ca_next) != previous.end())
        continue;

      bool is_found = is_boundary_found(ca_next, _cell_nm_axis, _is_anti_prl);
      if (is_found)
      {
        ALGOHEX_DEBUG_ONLY(std::cout << "found " << _is_anti_prl << " parallel boundary cell " << ca_next.first << " ";)
        if (_candidate_chs.find(ca_next.first) != _candidate_chs.end())
        {
          ca_target = ca_next;
          previous[ca_target] = ca_cur;
          previous_incrs[ca_target] = inc_val;

          ALGOHEX_DEBUG_ONLY(std::cout << "found cell " << ca_next.first << std::endl;)

          break;
        }
      }

      //push to queue
      //floodfill
      previous[ca_next] = ca_cur;
      previous_incrs[ca_next] = inc_val;

      que.push(ca_next);
    }
  }

  //recover quaternions
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    cell_quaternions_[*vc_it] = old_cell_quaternions[*vc_it];

  if (ca_target.first.is_valid())
  {
    //find cells
    //get the path
    std::vector<CH> rpath;
    CellAxis ca_it = ca_target;
    while (ca_it.first.is_valid())
    {
      rpath.push_back(ca_it.first);
      ca_it = previous[ca_it];
    }

    //check if there are duplicate elements
    std::set<CH> path(rpath.begin(), rpath.end());
    if (rpath.empty() || path.size() != rpath.size())
    {
      std::cerr << "Empty dpath or has duplicate elements, skip." << std::endl;
      return CH(-1);
    }
    else
    {//store
      _dpath.reserve(rpath.size());
      for (auto r_it = rpath.rbegin(); r_it != rpath.rend(); ++r_it)
        _dpath.push_back(*r_it);
    }
  }

  return ca_target.first;
}

template<class MeshT>
CH
SeparableSingularArcFinderT<MeshT>::
find_dual_path_to_feature_face(const HEH _heh, const CH _ch_s, std::map<CH, int> &_cell_nm_axis, const int _target_axis,
                               const std::set<CH> &_onering_chs, std::vector<CH> &_dpath, const bool _is_anti_prl)
{
  VH vhf = mesh_.halfedge(_heh).from_vertex();

  //locally align quaternion field
  //store old quaterions
  std::map<CH, Quaternion> old_cell_quaternions;
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    old_cell_quaternions[*vc_it] = cell_quaternions_[*vc_it];

  //align quaternions to singular or feature edges or boundary faces
  QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                      feature_fprop_, feature_edge_, vhf, true);

  //get blocking cells if at zero feature sectors
  std::set<FH> blk_fhs;
  std::set<CH> blk_cells;
  if (enable_blocking_)
  {
    blk_cells = get_blocking_cells(_heh, blk_fhs);
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << " Blocking fhs ";)
  ALGOHEX_DEBUG_ONLY(for (auto fhi: blk_fhs)
                       std::cerr << " " << fhi;
                             std::cerr << std::endl;)

  //
  std::map<CellAxis, CellAxis> previous;
  std::map<CellAxis, double> previous_incrs;

  //edge direction in parametrization coordinates
  int u = _target_axis;

  CellAxis ca(_ch_s, u);

  previous[ca] = CellAxis(CH(-1), -1);

  ALGOHEX_DEBUG_ONLY(std::cout << "Floodfilling..." << " \n";)

  CellAxis ca_target(CH(-1), -1);
  std::queue<CellAxis> que;
  que.push(ca);
  while (!que.empty() && !ca_target.first.is_valid())
  {
    auto ca_cur = que.front();

    que.pop();

    auto chfs = mesh_.cell(ca_cur.first).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);

      auto fh_next = mesh_.face_handle(hf_opp);

      if (blk_fhs.find(fh_next) != blk_fhs.end())
        continue;

      bool is_found = is_target_feature_face_found(ca_cur, fh_next, _cell_nm_axis, _is_anti_prl);
      if (is_found)
      {
        ALGOHEX_DEBUG_ONLY(
                std::cout << "found " << _is_anti_prl << " parallel dual path to feature face" << ca_cur.first
                          << " ";)
        ca_target = ca_cur;

        break;
      }

      if (feature_fprop_[fh_next] > 0)
      {
        continue;
      }

      auto ch_next = mesh_.incident_cell(hf_opp);
      if (!ch_next.is_valid())
        continue;


      if (_onering_chs.find(ch_next) == _onering_chs.end())
        continue;

      auto c_next = ch_next.idx();

      bool increase = false;
      double inc_val = 0.;
      auto ca_next = next_cell_axis(ca_cur, hf_opp, inc_val, increase);

      if (_is_anti_prl)
        if (!increase)
          continue;

      //visited
      if (previous.find(ca_next) != previous.end())
        continue;

      //push to queue
      //floodfill
      previous[ca_next] = ca_cur;
      previous_incrs[ca_next] = inc_val;

      que.push(ca_next);
    }
  }

  //recover quaternions
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    cell_quaternions_[*vc_it] = old_cell_quaternions[*vc_it];

  if (ca_target.first.is_valid())
  {
    //find cells
    //get the path
    std::vector<CH> rpath;
    CellAxis ca_it = ca_target;
    while (ca_it.first.is_valid())
    {
      rpath.push_back(ca_it.first);
      ca_it = previous[ca_it];
    }

    //check if there are duplicate elements
    std::set<CH> path(rpath.begin(), rpath.end());
    if (rpath.empty() || path.size() != rpath.size())
    {
      std::cerr << "Empty dpath or has duplicate elements, skip." << std::endl;
      return CH(-1);
    }
    else
    {//store
      _dpath.reserve(rpath.size());
      for (auto r_it = rpath.rbegin(); r_it != rpath.rend(); ++r_it)
        _dpath.push_back(*r_it);
    }
  }

  return ca_target.first;
}


template<class MeshT>
HEH
SeparableSingularArcFinderT<MeshT>::find_zipper_node_edge_pair_with_dual_path(const HEH _heh, const CH _ch_s,
                                                                              std::map<CH, int> &_cell_he_axis,
                                                                              std::map<CH, int> &_cell_nm_axis,
                                                                              const int _target_axis,
                                                                              const std::set<CH> &_onering_chs,
                                                                              const std::set<HEH> &_candidate_hehs,
                                                                              const double _angle_threshold_tp,
                                                                              std::vector<CH> &_dpath)
{
  //
  VH vhf = mesh_.halfedge(_heh).from_vertex();

  //locally align quaternion field
  //store old quaterions
  std::map<CH, Quaternion> old_cell_quaternions;
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    old_cell_quaternions[*vc_it] = cell_quaternions_[*vc_it];

  //align quaternions to singular or feature edges or boundary faces
  QuaternionSmoothing::optimize_quaternions_at_vertex(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                      feature_fprop_, feature_edge_, vhf, true);

  //
  std::map<CellAxis, CellAxis> previous;
  std::map<CellAxis, double> previous_incrs;

  //edge direction in parametrization coordinates
  int u = _target_axis;

  CellAxis ca(_ch_s, u);

  previous[ca] = CellAxis(CH(-1), -1);


  ALGOHEX_DEBUG_ONLY(std::cerr << "Angle thresold(cos) " << _angle_threshold_tp << std::endl;
                             std::cout << "Floodfilling..." << " target axis " << u << " \n";)

  CellAxis ca_target(CH(-1), -1);
  std::queue<CellAxis> que;
  que.push(ca);
  while (!que.empty() && !ca_target.first.is_valid())
  {
    auto ca_cur = que.front();

    que.pop();

    auto chfs = mesh_.cell(ca_cur.first).halffaces();
    for (int i = 0; i < 4; ++i)
    {
      auto hf_opp = mesh_.opposite_halfface_handle(chfs[i]);

      if (feature_fprop_[mesh_.face_handle(hf_opp)] > 0)
      {
        if (same_sector_sges_only_)
          continue;
        else if (_cell_nm_axis[ca_cur.first] / 2 == ca_cur.second / 2)
          continue;
      }

      auto ch_next = mesh_.incident_cell(hf_opp);

      if (_onering_chs.find(ch_next) == _onering_chs.end())
        continue;

      auto c_next = ch_next.idx();

      bool increase = false;
      double inc_val = 0.;
      auto ca_next = next_cell_axis(ca_cur, hf_opp, inc_val, increase);


      if (inc_val > _angle_threshold_tp)
        continue;

      //visited
      if (previous.find(ca_next) != previous.end())
        continue;

      HEH he_test = is_zipper_node_edge_pair_found(ca_next, vhf, _heh, _cell_he_axis, _cell_nm_axis);

      if (he_test.is_valid())
      {
        if (_candidate_hehs.find(he_test) != _candidate_hehs.end())
        {
          ca_target = ca_next;
          previous[ca_target] = ca_cur;

          previous_incrs[ca_target] = inc_val;

          ALGOHEX_DEBUG_ONLY(std::cout << "found edge " << mesh_.edge_handle(he_test) << std::endl;)

          break;
        }
      }

      //push to queue
      //floodfill
      previous[ca_next] = ca_cur;
      previous_incrs[ca_next] = inc_val;

      que.push(ca_next);
    }
  }

  //recover quaternions
  for (auto vc_it = mesh_.vc_iter(vhf); vc_it.valid(); ++vc_it)
    cell_quaternions_[*vc_it] = old_cell_quaternions[*vc_it];

  //
  HEH he_other(-1);
  if (ca_target.first.is_valid())
  {
    he_other = singular_halfedge_in_cell(vhf, ca_target.first, _heh);

    //find cells
    //get the path
    std::vector<CH> rpath;
    CellAxis ca_it = ca_target;
    while (ca_it.first.is_valid())
    {
      rpath.push_back(ca_it.first);
      ca_it = previous[ca_it];
    }

    //check if there are duplicate elements
    std::set<CH> path(rpath.begin(), rpath.end());
    if (rpath.empty() || path.size() != rpath.size())
    {
      std::cerr << "Empty dpath or has duplicate elements, skip." << std::endl;
      return HEH(-1);
    }
    else
    {//store
      _dpath.reserve(rpath.size());
      for (auto r_it = rpath.rbegin(); r_it != rpath.rend(); ++r_it)
        _dpath.push_back(*r_it);
    }
  }

  return he_other;
}


template<class MeshT>
typename SeparableSingularArcFinderT<MeshT>::CellAxis
SeparableSingularArcFinderT<MeshT>::next_cell_axis(const CellAxis &_ca_cur, const HFH _hf_opp, double &_inc_val,
                                                   bool &_increasing) const
{
  Vec3d axis_dir = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[_ca_cur.first],
                                                           (AxisAlignment) _ca_cur.second);

  //check if it's increasing in u dir
  Vec3d face_dir = eigen_point(mesh_.normal(_hf_opp));
  double prd = face_dir.dot(axis_dir);
  _inc_val = prd;
  if (prd < 0.)
    _increasing = false;
  else
    _increasing = true;

  int next_axis = tq_.axis_after_transition(_ca_cur.second, trans_prop_[_hf_opp]);
  CH next_ch = mesh_.incident_cell(_hf_opp);

  return std::make_pair(next_ch, next_axis);
}

template<class MeshT>
HEH
SeparableSingularArcFinderT<MeshT>::is_anti_parallel_found(const CellAxis &_ca, const VH _vhf, const HEH _heh0,
                                                           std::map<CH, int> &_cell_he_axis)
{
  auto it = _cell_he_axis.find(_ca.first);
  if (it != _cell_he_axis.end() && _ca.second == it->second)
  {
    auto heh = singular_halfedge_in_cell(_vhf, _ca.first, _heh0);
    if (heh.is_valid())
    {
      int val0 = valence_[mesh_.edge_handle(_heh0)];
      int val = valence_[mesh_.edge_handle(heh)];
      if ((val == val0) || val == std::numeric_limits<int>::lowest())
        return heh;
    }
  }

  return HEH(-1);
}

template<class MeshT>
bool
SeparableSingularArcFinderT<MeshT>::is_boundary_found(const CellAxis &_ca, std::map<CH, int> &_cell_nm_axis,
                                                      const bool _is_anti_prl)
{
  if ((_is_anti_prl && _ca.second == _cell_nm_axis[_ca.first])
      || (!_is_anti_prl && negate(_ca.second) == _cell_nm_axis[_ca.first]))
  {
    if (mesh_.is_boundary(_ca.first))
      return true;
  }

  return false;
}

template<class MeshT>
bool
SeparableSingularArcFinderT<MeshT>::is_target_feature_face_found(const CellAxis &_ca, const FH _inc_fh,
                                                                 std::map<CH, int> &_cell_nm_axis,
                                                                 const bool _is_anti_prl)
{
  if (feature_fprop_[_inc_fh] > 0)
  {
    if ((_is_anti_prl && _ca.second == _cell_nm_axis[_ca.first])
        || (!_is_anti_prl && negate(_ca.second) == _cell_nm_axis[_ca.first]))
    {
      return true;
    }
  }

  return false;
}


template<class MeshT>
HEH
SeparableSingularArcFinderT<MeshT>::is_zipper_node_edge_pair_found(const CellAxis &_ca, const VH _vhf, const HEH _heh0,
                                                                   std::map<CH, int> &_cell_he_axis,
                                                                   std::map<CH, int> &_cell_nm_axis)
{
  auto it = _cell_he_axis.find(_ca.first);
  if (it != _cell_he_axis.end() && _ca.second == it->second)
  {
    auto heh = singular_halfedge_in_cell(_vhf, _ca.first, _heh0);
//            if(!mesh_.is_boundary(heh))
    if (heh.is_valid())
    {
      int val0 = valence_[mesh_.edge_handle(_heh0)];
      int val = valence_[mesh_.edge_handle(heh)];
      if ((val == -val0) || val == std::numeric_limits<int>::lowest())
        return heh;
    }
  }

  return HEH(-1);
}


template<class MeshT>
int
SeparableSingularArcFinderT<MeshT>::check_parameterization(Parameter &_pm)
{
  int n = 0;
  for (auto&[chi, vp]: _pm)
  {
    Eigen::Matrix<double, 3, 4> P;
    parametric_tet(_pm, chi, P);
    Eigen::Matrix3d E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    double det = E.determinant();
    if (det < 0.0)
    {
      std::cerr << "Flipped tet " << chi << " has volume " << det << std::endl;
      ++n;
    }
  }

  std::cerr << "Parameterization has " << n << " flipped tets!" << std::endl;

  return n;
}

template<class MeshT>
void
SeparableSingularArcFinderT<MeshT>::
parametric_tet(Parameter &_pm, const CH _ch, Eigen::Matrix<double, 3, 4> &_P)
{
  int i = 0;
  for (auto tv_it = mesh_.tv_iter(_ch); tv_it.valid(); ++tv_it)
  {
    _P.col(i++) = eigen_point(_pm[_ch][*tv_it]);
  }
}


template<class MeshT>
HEH
SeparableSingularArcFinderT<MeshT>::singular_halfedge_in_cell(const VH _vhf, const CH _ch, const HEH _he_excl) const
{
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
  {
    if (*cv_it != _vhf)
    {
      HEH hei = mesh_.find_halfedge(_vhf, *cv_it);
      if (hei.is_valid())
      {
        int vali = valence_[mesh_.edge_handle(hei)];
        if ((vali == -1 || vali == 1 || vali == std::numeric_limits<int>::lowest()) && hei != _he_excl)
          return hei;
      }
    }
  }

  return HEH(-1);
}

template<class MeshT>
std::set<HEH>
SeparableSingularArcFinderT<MeshT>::
get_candidate_halfedges(const std::set<HEH> &_hehs, const HEH _heh, int _valence, bool _same_valence)
{
  //get candidate edge set
  auto candidate_hehs = _hehs;

  std::set<HEH> cd_hehs;
  VH vh_f = mesh_.halfedge(_heh).from_vertex();
  for (const auto hehi: candidate_hehs)
  {
    auto ehi = mesh_.edge_handle(hehi);
    if (hehi != _heh &&
        ((_same_valence && valence_[ehi] == _valence)
         || (!_same_valence && valence_[ehi] == -_valence)
         || (valence_[ehi] == std::numeric_limits<int>::lowest() ||
             _valence == std::numeric_limits<int>::lowest())))//because of aligning to feature surface
    {
      if (feature_face_edge_[ehi])
      {
        bool valid_sector_after = is_valid_sector_after_fix(vh_f, hehi, _heh);
        if (valid_sector_after)
          cd_hehs.insert(hehi);
      }
      else
        cd_hehs.insert(hehi);
    }
  }

  return cd_hehs;
}


template<class MeshT>
double SeparableSingularArcFinderT<MeshT>::
get_angle_threshold_zipper_node(const VH _vh, const HEH _heh, const std::set<CH> &_or_chs) const
{
  double angle_threshold_tp = angle_threshold_tp_;
  if (n_incident_singular_edges_of_valence_one(mesh_, valence_, _vh) <= 2)
    angle_threshold_tp = 1.;

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, _or_chs, *mesh_.hec_iter(_heh),
                                           ft_hfhs);
  int n_nffe = n_incident_non_ffe_singular_edges_in_same_region(mesh_, feature_face_edge_, valence_, _vh,
                                                                por_chs);
  if ((n_nffe <= 2 && !find_feature_face_sges_only_) || (n_nffe <= 1 && find_feature_face_sges_only_))
    angle_threshold_tp = 1.;

  return angle_threshold_tp;
}

template<class MeshT>
double SeparableSingularArcFinderT<MeshT>::
get_angle_threshold(const VH _vh, const HEH _heh, const std::set<CH> &_or_chs) const
{
  double angle_threshold = angle_threshold_;
  if (n_incident_singular_edges_of_valence_one(mesh_, valence_, _vh) <= 2)
    angle_threshold = -1.;

  std::set<HFH> ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, _or_chs, *mesh_.hec_iter(_heh),
                                           ft_hfhs);
  int n_nffe = n_incident_non_ffe_singular_edges_in_same_region(mesh_, feature_face_edge_, valence_, _vh,
                                                                por_chs);
  if ((n_nffe <= 2 && !find_feature_face_sges_only_) || (n_nffe <= 1 && find_feature_face_sges_only_))
    angle_threshold = -1.;

  return angle_threshold;
}

template<class MeshT>
bool SeparableSingularArcFinderT<MeshT>::
is_valid_sector_after_fix(const VH _vh, const HEH _ff_sgheh, const HEH _nff_sgheh)
{
  auto or_chs = get_onering_cells(mesh_, _vh);
  std::set<HFH> bound_ft_hfhs;
  auto por_chs = get_part_of_onering_cells(mesh_, feature_fprop_, _vh, or_chs, *mesh_.hec_iter(_nff_sgheh),
                                           bound_ft_hfhs);

  HFH hf_s, hf_e;
  for (auto hehf_it = mesh_.hehf_iter(_ff_sgheh); hehf_it.valid(); ++hehf_it)
  {
    if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
    {
      if (bound_ft_hfhs.find(*hehf_it) != bound_ft_hfhs.end())
      {
        hf_s = *hehf_it;
      }
      else if (bound_ft_hfhs.find(mesh_.opposite_halfface_handle(*hehf_it)) != bound_ft_hfhs.end())
      {
        hf_e = *hehf_it;
      }
    }
  }

  if (!hf_s.is_valid() || !hf_e.is_valid())
  {
    std::cerr << "Error: couldn't find start halfface and end halfface at feature face sector!" << std::endl;
    return false;
  }

  HFH hf_change = mesh_.opposite_halfface_handle(mesh_.adjacent_halfface_in_cell(hf_s, _ff_sgheh));
  int old_trans = trans_prop_[hf_change];

  int rt_ax = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                 _ff_sgheh, hf_change);
  if (rt_ax == -1)
  {
    std::cerr << "Warning: invalid axis, return" << std::endl;
    return false;
  }

  //get the matching to turn the edge regular
  int e_trans = rt_ax + 4;

  int fsec_agl_before = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                            cell_quaternions_,
                                                                                            trans_prop_, _ff_sgheh,
                                                                                            hf_s, hf_e);

  int new_trans = tq_.mult_transitions_idx(old_trans, tq_.inverse_transition_idx(e_trans));
  trans_prop_[hf_change] = new_trans;
  trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(new_trans);
  int rt_ax2 = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  _ff_sgheh, hf_change);

  int fsec_agl_after = FieldAngleCalculatorT<MeshT>::feature_face_sector_angle_at_halfedge(mesh_, tq_,
                                                                                           cell_quaternions_,
                                                                                           trans_prop_, _ff_sgheh, hf_s,
                                                                                           hf_e);


  ALGOHEX_DEBUG_ONLY(std::cerr << "ff sghe " << mesh_.edge_handle(_ff_sgheh) << "f-t vhs "
                               << mesh_.halfedge(_ff_sgheh).from_vertex() << " "
                               << mesh_.halfedge(_ff_sgheh).to_vertex()
                               << " fs " << mesh_.face_handle(hf_s) << " trans " << trans_prop_[hf_s] << " fe "
                               << mesh_.face_handle(hf_e) << " trans " << trans_prop_[hf_e] << " changed fh "
                               << mesh_.face_handle(hf_change)
                               << " trans " << old_trans << " -> " << new_trans << " val "
                               << valence_[mesh_.edge_handle(_ff_sgheh)] << " rt axes: " << rt_ax << "->" << rt_ax2
                               << " e_trans " << e_trans << " sec agl before " << fsec_agl_before
                               << " sec agl after " << fsec_agl_after << " dihedral angle "
                               << dihedral_angle_from_halfface_to_halfface(mesh_, _ff_sgheh, hf_s, hf_e) << std::endl;)

  //restore matching
  trans_prop_[hf_change] = old_trans;
  trans_prop_[mesh_.opposite_halfface_handle(hf_change)] = tq_.inverse_transition_idx(old_trans);

  if (feature_edge_[mesh_.edge_handle(_ff_sgheh)] == 0)
  {
    if (fsec_agl_after != 2)
      return false;
  }
  else
  {
    auto da = dihedral_angle_from_halfface_to_halfface(mesh_, _ff_sgheh, hf_s, hf_e);
    if (da >= M_PI)
    {
      if (fsec_agl_after < 2)
        return false;
    }
    else if (fsec_agl_after < 1)
      return false;
  }

  return true;
}

template<class MeshT>
std::set<CH>
SeparableSingularArcFinderT<MeshT>::
get_candidate_cells(const VH _vh, const HEH _sg_heh, const std::set<HEH> &_ft_hehs, const bool _is_anti_prl)
{
  std::set<CH> cd_bdy_chs;

  if (!_ft_hehs.empty())
  {
    EH sg_eh = mesh_.edge_handle(_sg_heh);

    HFH hf_sg = *mesh_.hehf_iter(_sg_heh);
    CH ch_sg = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_sg));
    int rt_axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_,
                                                                     valence_, _sg_heh,
                                                                     hf_sg);
    if (valence_[sg_eh] == 1)
      rt_axis = negate(rt_axis);

    auto or_chs = get_onering_cells(mesh_, _vh);

    //store sectors that are orthogonal to sge
    std::vector<std::vector<HEH>> bdy_sector_hehs_0, bdy_sector_hehs_90, bdy_sector_hehs_180, bdy_sector_hehs_270;
    //get feature(singular) sectors on feature surface
    std::vector<std::vector<HEH>> v_sec_hehs;
    std::vector<std::set<FH>> v_sec_fhs;
    get_feature_sectors_at_feature_edge_vertex(mesh_, feature_fprop_, feature_edge_, feature_edge_vertex_, valence_,
                                               _vh, v_sec_hehs, v_sec_fhs);

    for (auto j = 0u; j < v_sec_hehs.size(); ++j)
    {
      HFH hf_m;
      for (auto hehf_it = mesh_.hehf_iter(v_sec_hehs[j][1]); hehf_it.valid(); ++hehf_it)
      {
        if (mesh_.is_boundary(mesh_.opposite_halfface_handle(*hehf_it)))
        {
          hf_m = *hehf_it;
          break;
        }
      }
      if (!hf_m.is_valid())
        continue;

      int is_prl = fac_.angle_of_axis_in_cell_and_normal_direction(rt_axis, ch_sg, hf_m, or_chs);
      if (is_prl == 0)// 0 is orthogonal, sector and sge is parallel
        continue;

      int ea_st = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_,
                                                                                     trans_prop_, valence_,
                                                                                     v_sec_hehs[j], v_sec_fhs[j]);

      if (ea_st == 1)// 90 degree
        bdy_sector_hehs_90.push_back(v_sec_hehs[j]);
      else if (ea_st == 0)
      {
        double sec_angle = FieldAngleCalculatorT<MeshT>::sector_angle(mesh_, v_sec_hehs[j]);
        if (sec_angle < M_PI_2) // 0 degree
          bdy_sector_hehs_0.push_back(v_sec_hehs[j]);
      }
      else if (ea_st == 2)
        bdy_sector_hehs_180.push_back(v_sec_hehs[j]);
      else if (ea_st == 3)
        bdy_sector_hehs_270.push_back(v_sec_hehs[j]);
    }

    //erase
    if ((valence_[sg_eh] == 1 && _is_anti_prl) || (valence_[sg_eh] == -1 && !_is_anti_prl))
    {
      for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
      {
        if (mesh_.is_boundary(*vc_it))
          cd_bdy_chs.insert(*vc_it);
      }

      //shouldn't go to zero sector or 90 degree sector
      for (const auto &sec_hehs: bdy_sector_hehs_0)
      {
        for (auto i = 1u; i < sec_hehs.size() - 1; ++i)
        {
          for (auto hec_it = mesh_.hec_iter(sec_hehs[i]); hec_it.valid(); ++hec_it)
          {
            if (mesh_.is_boundary(*hec_it))
              cd_bdy_chs.erase(*hec_it);
          }
        }
      }

      for (const auto &sec_hehs: bdy_sector_hehs_90)
      {
        for (auto i = 1u; i < sec_hehs.size() - 1; ++i)
        {
          for (auto hec_it = mesh_.hec_iter(sec_hehs[i]); hec_it.valid(); ++hec_it)
          {
            if (mesh_.is_boundary(*hec_it))
              cd_bdy_chs.erase(*hec_it);
          }
        }
      }

      //erase 180 sector if 270 sector exist for better distorsion
      if (!bdy_sector_hehs_270.empty())
      {
        for (const auto &sec_hehs: bdy_sector_hehs_180)
        {
          for (auto i = 1u; i < sec_hehs.size() - 1; ++i)
          {
            for (auto hec_it = mesh_.hec_iter(sec_hehs[i]); hec_it.valid(); ++hec_it)
            {
              if (mesh_.is_boundary(*hec_it))
                cd_bdy_chs.erase(*hec_it);
            }
          }
        }
      }
    }

    //add
    if ((valence_[sg_eh] == -1 && _is_anti_prl) || (valence_[sg_eh] == 1 && !_is_anti_prl))
    {
      if (!bdy_sector_hehs_0.empty())
      {
        for (const auto &sec_hehs: bdy_sector_hehs_0)
        {
          for (auto i = 1u; i < sec_hehs.size() - 1; ++i)
          {
            for (auto hec_it = mesh_.hec_iter(sec_hehs[i]); hec_it.valid(); ++hec_it)
            {
              if (mesh_.is_boundary(*hec_it))
                cd_bdy_chs.insert(*hec_it);
            }
          }
        }
      }
      else if (!bdy_sector_hehs_90.empty())
      {
        for (const auto &sec_hehs: bdy_sector_hehs_90)
        {
          for (auto i = 1u; i < sec_hehs.size() - 1; ++i)
          {
            for (auto hec_it = mesh_.hec_iter(sec_hehs[i]); hec_it.valid(); ++hec_it)
            {
              if (mesh_.is_boundary(*hec_it))
                cd_bdy_chs.insert(*hec_it);
            }
          }
        }
      }
      else
      {
        for (const auto &sec_hehs: bdy_sector_hehs_180)
        {
          for (auto i = 1u; i < sec_hehs.size() - 1; ++i)
          {
            for (auto hec_it = mesh_.hec_iter(sec_hehs[i]); hec_it.valid(); ++hec_it)
            {
              if (mesh_.is_boundary(*hec_it))
                cd_bdy_chs.insert(*hec_it);
            }
          }
        }
      }
    }
  }
  else
  {
    for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    {
      if (mesh_.is_boundary(*vc_it))
        cd_bdy_chs.insert(*vc_it);
    }
  }

  return cd_bdy_chs;
}

template<class MeshT>
std::set<CH>
SeparableSingularArcFinderT<MeshT>::get_blocking_cells(const HEH _heh_s, std::set<FH> &_blk_fhs,
                                                       const bool _is_anti_prl)
{
  _blk_fhs.clear();

  std::set<CH> blk_chs;
  int val = valence_[mesh_.edge_handle(_heh_s)];
  if ((val <= -1 && _is_anti_prl) || (val >= 1 && !_is_anti_prl))
    return blk_chs;

  VH vhf = mesh_.halfedge(_heh_s).from_vertex();
  if (!feature_edge_vertex_[vhf])
    return blk_chs;

  //get feature(singular) sectors on feature surface
  std::vector<std::vector<HEH>> v_sec_hehs;
  std::vector<std::set<FH>> v_sec_fhs;
  get_feature_sectors_at_feature_edge_vertex(mesh_, feature_fprop_, feature_edge_, feature_edge_vertex_, valence_, vhf,
                                             v_sec_hehs, v_sec_fhs);

  for (auto j = 0u; j < v_sec_hehs.size(); ++j)
  {
    int ea_st = FieldAngleCalculatorT<MeshT>::field_angle_status_at_feature_sector(mesh_, tq_, cell_quaternions_,
                                                                                   trans_prop_, valence_, v_sec_hehs[j],
                                                                                   v_sec_fhs[j]);

    if (ea_st >= 0 && ea_st <= 1)
    {
      for (const auto fhi: v_sec_fhs[j])
      {
        _blk_fhs.insert(fhi);

        auto ch0 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 0));
        auto ch1 = mesh_.incident_cell(mesh_.halfface_handle(fhi, 1));

        if (ch0.is_valid())
          blk_chs.insert(ch0);
        if (ch1.is_valid())
          blk_chs.insert(ch1);
      }
    }
  }

  return blk_chs;
}

template<class MeshT>
std::vector<HEH>
SeparableSingularArcFinderT<MeshT>::
sort_starting_halfedges_for_antiparallel(const VH _vh)
{
  std::set<HEH> all_hehs;
  for (auto voh_it = mesh_.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto oeh = mesh_.edge_handle(*voh_it);
    if (valence_[oeh] != 0 & !feature_face_edge_[oeh])
      all_hehs.insert(*voh_it);
  }

  auto cell_region_id = get_region_id_of_onering_cells(mesh_, feature_fprop_, _vh);

  std::vector<HEH> sorted_hehs;
  while (!all_hehs.empty())
  {
    auto he_pair = most_parallel_halfedge_pair_in_set(all_hehs, cell_region_id);
    if (he_pair.first.is_valid())
    {
      sorted_hehs.push_back(he_pair.first);
      sorted_hehs.push_back(he_pair.second);

      all_hehs.erase(he_pair.first);
      all_hehs.erase(he_pair.second);
    }
    else
      break;
  }

  //
  if (!sorted_hehs.empty())
  {
    auto ag0 = FieldAngleCalculatorT<MeshT>::average_field_edge_angle(mesh_, tq_, cell_quaternions_,
                                                                      trans_prop_, valence_,
                                                                      feature_edge_, mesh_.edge_handle(sorted_hehs[0]));
    auto ag1 = FieldAngleCalculatorT<MeshT>::average_field_edge_angle(mesh_, tq_, cell_quaternions_,
                                                                      trans_prop_, valence_,
                                                                      feature_edge_, mesh_.edge_handle(sorted_hehs[1]));


    if (ag1 < ag0)
      std::swap(sorted_hehs[0], sorted_hehs[1]);
  }

  //in case it's odd number of singular edges
  sorted_hehs.insert(sorted_hehs.end(), all_hehs.begin(), all_hehs.end());

  for(const auto hehi : all_hehs)
    if(feature_face_edge_[mesh_.edge_handle(hehi)])
      sorted_hehs.push_back(hehi);

  return sorted_hehs;
}

template<class MeshT>
std::pair<HEH, HEH>
SeparableSingularArcFinderT<MeshT>::
most_parallel_halfedge_pair_in_set(const std::set<HEH> &_all_hehs, std::map<CH, int>& _cell_region_id)
{
  using HEP = std::pair<HEH, HEH>;
  using DHEP = std::pair<double, HEP>;
  std::vector<DHEP> dheps;
  for (const auto hehi: _all_hehs)
  {
    HEH he_start = hehi;
    auto eh_start = mesh_.edge_handle(he_start);
    int rg_id = _cell_region_id[*mesh_.hec_iter(he_start)];
    int val = valence_[eh_start];

    for (const auto hehi2: _all_hehs)
    {
      if ((hehi2 != he_start) && (valence_[mesh_.edge_handle(hehi2)] == val ||
                                  valence_[mesh_.edge_handle(hehi2)] == std::numeric_limits<int>::lowest()))
      {
        int rg_id2 = _cell_region_id[*mesh_.hec_iter(hehi2)];
        if((same_sector_sges_only_ && rg_id == rg_id2) || (!same_sector_sges_only_ && rg_id != rg_id2))
        {
          auto angle = edge_angle(mesh_, he_start, hehi2);
          dheps.emplace_back(std::fabs(angle - M_PI), std::make_pair(he_start, hehi2));
        }
      }
    }
  }

  if (dheps.empty())
    return std::make_pair(HEH(-1), HEH(-1));

  std::sort(dheps.begin(), dheps.end());

  return dheps[0].second;
}


template<class MeshT>
bool
SeparableSingularArcFinderT<MeshT>::
is_feature_cell(const CH _ch) const
{
  for (auto cf_it = mesh_.cf_iter(_ch); cf_it.valid(); ++cf_it)
    if (feature_fprop_[*cf_it] > 0)
    {
      return true;
    }

  return false;
}

}
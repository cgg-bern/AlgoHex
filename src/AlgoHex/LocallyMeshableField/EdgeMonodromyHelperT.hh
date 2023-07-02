/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#ifdef ALGOHEX_VERBOSE
#define ALGOHEX_DEBUG_ONLY(x) x
#else
#define ALGOHEX_DEBUG_ONLY(x)
#endif

#include <AlgoHex/TypeDef.hh>
#include <AlgoHex/TransitionQuaternionEigen.hh>

namespace AlgoHex
{
class EdgeMonodromyHelperT
{
public:
  template<class MeshT>
  static int halfedge_transition_index(const MeshT &_mesh, const TransitionQuaternion &_tq, const HFP<int> &_trans_prop,
                                       const HEH _heh, const HFH _hfh_start)
  {
    ALGOHEX_DEBUG_ONLY(if (!_heh.is_valid())
                       {
                         std::cerr << "Error: input heh is invalid" << std::endl;
                         return -1;
                       })
    ALGOHEX_DEBUG_ONLY(if (!_hfh_start.is_valid())
                       {
                         std::cerr << "Error: input hfh is invalid" << std::endl;
                         return -1;
                       })

    ALGOHEX_DEBUG_ONLY(if (_mesh.is_boundary(_hfh_start))
                       {
                         std::cerr << "Error: start halfface is boundary" << std::endl;
                         return -1;
                       })

    auto hfh_it = _hfh_start;
    //if starting from boundary face, go to the neighbour
    if (_mesh.is_boundary(_mesh.face_handle(hfh_it)))
      hfh_it = _mesh.opposite_halfface_handle(_mesh.adjacent_halfface_in_cell(hfh_it, _heh));

    int idx = 0;
    do
    {
      if (_mesh.is_boundary(hfh_it))
        break;

      idx = _tq.mult_transitions_idx(_trans_prop[hfh_it], idx);

      auto hfh_adj = _mesh.adjacent_halfface_in_cell(hfh_it, _heh);
      if (!hfh_adj.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! Halfedge handle: " << _heh << " , start halfface: "
                  << _hfh_start << std::endl;
        idx = -1;
        break;
      }
      hfh_it = _mesh.opposite_halfface_handle(hfh_adj);
    }
    while (hfh_it != _hfh_start);


    return idx;
  }

  template<class MeshT>
  static int halfedge_transition_index(const MeshT &_mesh, const TransitionQuaternion &_tq, const HFP<int> &_trans_prop,
                                       const HEH _heh, const HFH _hfh_start, const HFH _hfh_end)
  {
    if (!_heh.is_valid())
    {
      std::cerr << "Error: input heh is invalid" << std::endl;
      return -1;
    }
    if (!_hfh_start.is_valid() || !_hfh_end.is_valid())
    {
      std::cerr << "Error: input hfh is invalid" << std::endl;
      return -1;
    }

    if (_mesh.is_boundary(_hfh_start))
    {
      std::cerr << "Error: start halfface is boundary" << std::endl;
      return -1;
    }

    auto hfh_it = _hfh_start;
    //if starting from boundary face, go to the neighbour
    if (_mesh.is_boundary(_mesh.face_handle(hfh_it)))
      hfh_it = _mesh.opposite_halfface_handle(_mesh.adjacent_halfface_in_cell(hfh_it, _heh));

    int idx = 0;
    do
    {
      if (hfh_it == _hfh_end || _mesh.is_boundary(hfh_it))
        break;

      idx = _tq.mult_transitions_idx(_trans_prop[hfh_it], idx);

      auto hfh_adj = _mesh.adjacent_halfface_in_cell(hfh_it, _heh);
      if (!hfh_adj.is_valid())
      {
        std::cerr << "Error: adjacent halfface is invalid! Halfedge handle: " << _heh << " , start halfface: "
                  << _hfh_start << std::endl;
        idx = -1;
        break;
      }
      hfh_it = _mesh.opposite_halfface_handle(hfh_adj);
    }
    while (hfh_it != _hfh_start);


    return idx;
  }

  template<class MeshT>
  static int halfedge_rotational_axis_idx(const MeshT &_mesh, const TransitionQuaternion &_tq,
                                          const CP <Quaternion> &_cell_quaternions,
                                          const HFP<int> &_trans_prop, const EP<int> &_valence, const HEH _heh,
                                          const HFH _hfh_start)
  {
    using Point = typename MeshT::PointT;

    if (_valence[_mesh.edge_handle(_heh)] == 0)
      return -1;

    if (_mesh.is_boundary(_heh))
    {// axis in the start cell
      HFH hf_start = *_mesh.hehf_iter(_heh);
      //should be boundary halfface if oriented
      HFH hf_end = _mesh.opposite_halfface_handle(*_mesh.hehf_iter(_mesh.opposite_halfedge_handle(_heh)));

      CH ch0 = _mesh.incident_cell(hf_start);
      Point n0 = _mesh.normal(_mesh.opposite_halfface_handle(hf_start));
      int closest_n0 = AxisAlignmentHelpers::closest_axis(_cell_quaternions[ch0], n0).second;

      CH ch1 = _mesh.incident_cell(_mesh.opposite_halfface_handle(hf_end));
      Point n1 = _mesh.normal(hf_end);
      int closest_n1 = AxisAlignmentHelpers::closest_axis(_cell_quaternions[ch1], n1).second;

      int T_inv_id = _tq.inverse_transition_idx(halfedge_transition_index(_mesh, _tq, _trans_prop, _heh, hf_start));

      int n1_in_q0 = _tq.axis_after_transition(closest_n1, T_inv_id);
      if (closest_n0 != n1_in_q0)
      {
        if (closest_n0 / 2 == n1_in_q0 / 2)
        {
          std::cerr << "Warning: non-meshable footprint at boundary singular edge, n0: " << closest_n0 << " n1in0: "
                    << n1_in_q0 << ", return edge direction." << std::endl;
          VH vh0 = _mesh.halfedge(_heh).from_vertex();
          VH vh1 = _mesh.halfedge(_heh).to_vertex();
          Point vn = _mesh.vertex(vh1) - _mesh.vertex(vh0);

          return AxisAlignmentHelpers::closest_axis(_cell_quaternions[ch0], vn).second;
        }
        else
        {
          return the_third_axis((AxisAlignment) n1_in_q0, (AxisAlignment) closest_n0);
        }
      }
      else
        return -1;
    }
    else
    { // axis in the specified cell
      auto trans = halfedge_transition_index(_mesh, _tq, _trans_prop, _heh, _hfh_start);
      return _tq.rotation_axis_2d(trans);
    }

    //never reach here
    return -1;
  }

  template<class MeshT>
  static int halfedge_transition_index_in_cell(const MeshT &_mesh, const TransitionQuaternion &_tq,
                                               const CP <Quaternion> &_cell_quaternions,
                                               const HFP<int> &_trans_prop, const EP<int> &_valence, const HEH _heh,
                                               const HFH _hfh_start)
  {
    int trans = -1;
    if (_mesh.is_boundary(_heh))
    {
      if (abs(_valence[_mesh.edge_handle(_heh)]) == 1)
      {
        auto ch_s = _mesh.incident_cell(_mesh.opposite_halfface_handle(_hfh_start));
        if (!ch_s.is_valid())
          return -1;

        //rt axis in the start boundary cell
        int rt_ax = halfedge_rotational_axis_idx(_mesh, _tq, _cell_quaternions, _trans_prop, _valence, _heh,
                                                 _hfh_start);

        HFH hfh_it = *_mesh.hehf_iter(_heh);

        do
        {
          if (_mesh.is_boundary(hfh_it))
            break;

          CH ch_it = _mesh.incident_cell(hfh_it);

          if (ch_it == ch_s)
            break;


          auto hfh_adj = _mesh.adjacent_halfface_in_cell(hfh_it, _heh);
          if (!hfh_adj.is_valid())
          {
            std::cerr << "Error: adjacent halfface is invalid! Halfedge handle: " << _heh << " , start halfface: "
                      << _hfh_start << std::endl;
            break;
          }
          hfh_it = _mesh.opposite_halfface_handle(hfh_adj);
          rt_ax = _tq.axis_after_transition(rt_ax, _trans_prop[hfh_it]);
        }
        while (1);

        return rt_ax + 4;
      }
    }
    else
    {
      trans = halfedge_transition_index(_mesh, _tq, _trans_prop, _heh, _hfh_start);
    }

    return trans;
  }


  template<class MeshT>
  static int get_dominant_axis_in_cell(const MeshT &_mesh, const TransitionQuaternion &_tq,
                                       const CP <Quaternion> &_cell_quaternions,
                                       const HFP<int> &_trans_prop, const HEH _heh)
  {
    auto dir = _mesh.vertex(_mesh.halfedge(_heh).to_vertex()) -
               _mesh.vertex(_mesh.halfedge(_heh).from_vertex());
    dir.normalize();
    Vec3d vn = ovm2eigen(dir);

    std::vector<double> dprds(3, 0.);
    for (int i = 0; i < 3; ++i)
    {
      int axi = 2 * i;
      for (auto hehf_it = _mesh.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
      {
        CH chi = _mesh.incident_cell(_mesh.opposite_halfface_handle(*hehf_it));
        dprds[i] += vn.dot(AxisAlignmentHelpers::quaternion_vector(_cell_quaternions[chi], AxisAlignment(axi)));
        //axis in the next cell
        axi = _tq.axis_after_transition(axi, _trans_prop[*hehf_it]);
      }
    }

    //choose the max
    double max_drd = -1.;
    int max_ax = -1;
    for (int i = 0; i < 3; ++i)
    {
      if (std::abs(dprds[i]) > max_drd)
      {
        max_drd = std::abs(dprds[i]);
        max_ax = 2 * i;
      }
    }

    if (dprds[max_ax / 2] < 0)
      max_ax = max_ax + 1;


    return max_ax;
  }

  template<class MeshT>
  static int
  halfedge_axis_idx(const MeshT &_mesh, const CP <Quaternion> &_cell_quaternions, const HEH _heh, const CH _ch,
                    const int _axis)
  {
    if (_axis < 0 || _axis > 5)
    {
      std::cerr << "Error: input axis is not valid!" << std::endl;
      return -1;
    }

    auto dir = _mesh.vertex(_mesh.halfedge(_heh).to_vertex()) -
               _mesh.vertex(_mesh.halfedge(_heh).from_vertex());

    Vec3d ax_vec = AxisAlignmentHelpers::quaternion_vector(_cell_quaternions[_ch], (AxisAlignment) _axis);

    if (ax_vec.dot(ovm2eigen(dir)) < 0.)
      return _axis % 2 == 0 ? _axis + 1 : _axis - 1;

    return _axis;
  }

  template<class MeshT>
  static int get_halfedge_aligned_axis_in_cell(const MeshT &_mesh, const TransitionQuaternion &_tq,
                                               const CP <Quaternion> &_cell_quaternions,
                                               const HFP<int> &_trans_prop, const EP<int> &_valence,
                                               const EP<int> &_feature_edge, const HEH _heh, const HFH _hfh)
  {
    EH eh = _mesh.edge_handle(_heh);
    auto dir = _mesh.vertex(_mesh.halfedge(_heh).to_vertex()) - _mesh.vertex(_mesh.halfedge(_heh).from_vertex());
    dir.normalize();

    CH ch_s = _mesh.incident_cell(_hfh);

    if (_valence[eh] % 4 != 0)
    {
      HFH hf_s = _mesh.opposite_halfface_handle(_mesh.adjacent_halfface_in_cell(_hfh, _heh));
      //if boundary edge, the specified halfface is not used
      bool e_bdy = _mesh.is_boundary(_heh);
      if (e_bdy)
        ch_s = _mesh.incident_cell(*_mesh.hehf_iter(_heh));

      int axis = halfedge_rotational_axis_idx(_mesh, _tq, _cell_quaternions, _trans_prop, _valence, _heh, hf_s);
      int e_axis = halfedge_axis_idx(_mesh, _cell_quaternions, _heh, ch_s, axis);

      return e_axis;
    }
    else if (_valence[eh] % 4 == 0)
    {
      if (_feature_edge[eh] > 0)
      {
        auto da = AxisAlignmentHelpers::closest_axis(_cell_quaternions[ch_s], ovm2eigen(dir));

        return da.second;
      }
      else if (_valence[eh] != 0)
      {
        int dm_ax = get_dominant_axis_in_cell(_mesh, _tq, _cell_quaternions, _trans_prop, _heh);

        return dm_ax;
      }
    }

    return -1;
  }

  template<class MeshT>
  static std::map<CH, int> get_halfedge_aligned_axis_in_cells(const MeshT &_mesh, const TransitionQuaternion &_tq,
                                                              const CP <Quaternion> &_cell_quaternions,
                                                              const HFP<int> &_trans_prop, const EP<int> &_valence,
                                                              const EP<int> &_feature_edge, const HEH &_heh)
  {
    std::map<CH, int> cell_axis;
    EH eh = _mesh.edge_handle(_heh);
    auto dir = _mesh.vertex(_mesh.halfedge(_heh).to_vertex()) - _mesh.vertex(_mesh.halfedge(_heh).from_vertex());
    dir.normalize();

    if (_valence[eh] % 4 != 0)
    {
      auto hehf_it = _mesh.hehf_iter(_heh);
      CH ch_s(-1);
      bool e_bdy = _mesh.is_boundary(_heh);
      if (!e_bdy)
        ch_s = _mesh.incident_cell(_mesh.opposite_halfface_handle(*hehf_it));
      else
        ch_s = _mesh.incident_cell(*hehf_it);

      int axis = halfedge_rotational_axis_idx(_mesh, _tq, _cell_quaternions, _trans_prop, _valence, _heh, *hehf_it);
      int e_axis = halfedge_axis_idx(_mesh, _cell_quaternions, _heh, ch_s, axis);

      cell_axis[ch_s] = e_axis;

      if (e_bdy)
        ++hehf_it;

      for (; hehf_it.valid(); ++hehf_it)
      {
        auto chi = _mesh.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          if (chi == ch_s)
            break;

          e_axis = _tq.axis_after_transition(e_axis, _trans_prop[*hehf_it]);
          cell_axis[chi] = e_axis;
        }
      }
    }
    else if (_valence[eh] % 4 == 0)
    {
      if (_feature_edge[eh] > 0)
      {
        for (auto hec_it = _mesh.hec_iter(_heh); hec_it.valid(); ++hec_it)
        {
          auto da = AxisAlignmentHelpers::closest_axis(_cell_quaternions[*hec_it], ovm2eigen(dir));
          int e_axis = da.second;
          cell_axis[*hec_it] = e_axis;
        }
      }
      else if (_valence[eh] != 0)
      {
        auto hehf_it = _mesh.hehf_iter(_heh);

        CH ch_s = _mesh.incident_cell(_mesh.opposite_halfface_handle(*hehf_it));
        int dm_ax = get_dominant_axis_in_cell(_mesh, _tq, _cell_quaternions, _trans_prop, _heh);
        cell_axis[ch_s] = dm_ax;

        for (; hehf_it.valid(); ++hehf_it)
        {
          CH chi = _mesh.incident_cell(*hehf_it);
          dm_ax = _tq.axis_after_transition(dm_ax, _trans_prop[*hehf_it]);
          cell_axis[chi] = dm_ax;
        }
      }
    }

    return cell_axis;
  }

  //TODO: support high valence by checking angle defect
  template<class MeshT>
  static std::map<CH, int>
  get_halfedge_aligned_axis_in_cells_wrt_valence(const MeshT &_mesh, const TransitionQuaternion &_tq,
                                                 const CP <Quaternion> &_cell_quaternions,
                                                 const HFP<int> &_trans_prop, const EP<int> &_valence,
                                                 const EP<int> &_feature_edge, const HEH &_heh)
  {
    std::map<CH, int> cell_axis;
    EH eh = _mesh.edge_handle(_heh);
    auto dir = _mesh.vertex(_mesh.halfedge(_heh).to_vertex()) - _mesh.vertex(_mesh.halfedge(_heh).from_vertex());
    dir.normalize();

    if (_valence[eh] % 4 != 0)
    {
      auto hehf_it = _mesh.hehf_iter(_heh);
      CH ch_s(-1);
      bool e_bdy = _mesh.is_boundary(_heh);
      if (!e_bdy)
        ch_s = _mesh.incident_cell(_mesh.opposite_halfface_handle(*hehf_it));
      else
        ch_s = _mesh.incident_cell(*hehf_it);

      int axis = halfedge_rotational_axis_idx(_mesh, _tq, _cell_quaternions, _trans_prop, _valence, _heh, *hehf_it);
      int e_axis = axis;
      if (_valence[eh] == 1)
        e_axis = axis % 2 == 0 ? axis + 1 : axis - 1;
      else if (std::abs(_valence[eh]) == 2)
      {
        e_axis = halfedge_axis_idx(_mesh, _cell_quaternions, _heh, ch_s, axis);
      }

      cell_axis[ch_s] = e_axis;

      if (e_bdy)
        ++hehf_it;

      for (; hehf_it.valid(); ++hehf_it)
      {
        auto chi = _mesh.incident_cell(*hehf_it);
        if (chi.is_valid())
        {
          if (chi == ch_s)
            break;

          e_axis = _tq.axis_after_transition(e_axis, _trans_prop[*hehf_it]);
          cell_axis[chi] = e_axis;
        }
      }
    }
    else if (_valence[eh] % 4 == 0)
    {
      if (_feature_edge[eh] > 0)
      {
        for (auto hec_it = _mesh.hec_iter(_heh); hec_it.valid(); ++hec_it)
        {
          auto da = AxisAlignmentHelpers::closest_axis(_cell_quaternions[*hec_it], ovm2eigen(dir));
          int e_axis = da.second;
          cell_axis[*hec_it] = e_axis;
        }
      }
      else if (_valence[eh] != 0)
      {
        auto hehf_it = _mesh.hehf_iter(_heh);

        CH ch_s = _mesh.incident_cell(_mesh.opposite_halfface_handle(*hehf_it));
        int dm_ax = get_dominant_axis_in_cell(_mesh, _tq, _cell_quaternions, _trans_prop, _heh);
        cell_axis[ch_s] = dm_ax;

        for (; hehf_it.valid(); ++hehf_it)
        {
          CH chi = _mesh.incident_cell(*hehf_it);
          dm_ax = _tq.axis_after_transition(dm_ax, _trans_prop[*hehf_it]);
          cell_axis[chi] = dm_ax;
        }
      }
    }

    return cell_axis;
  }
};
}



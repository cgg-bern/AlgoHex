/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include "EdgeMonodromyHelperT.hh"

namespace AlgoHex
{
class QuaternionSmoothing
{
public:
  using VecAxis = std::pair<Vec3d, int>;

  template<class MeshT>
  static Quaternion construct_quaternion(const Quaternion &_qt, const std::vector<VecAxis> &_directions)
  {
    int ax0 = _directions[0].second / 2, ax1 = _directions[1].second / 2;
    if (ax0 == ax1)
      return _qt;

    int axis2 = the_third_axis((AxisAlignment) _directions[0].second, (AxisAlignment) _directions[1].second);
    int ax2 = axis2 / 2;

    Eigen::Matrix3d R;
    R.col(ax0) = _directions[0].first;
    R.col(ax0).normalize();

    R.col(ax1) = _directions[1].first;
    R.col(ax1).normalize();

    R.col(ax2) = R.col(ax0).cross(R.col(ax1));

    if (_directions[0].second % 2 == 1)
      R.col(ax0) = -R.col(ax0);

    if (_directions[1].second % 2 == 1)
      R.col(ax1) = -R.col(ax1);

    if (axis2 % 2 == 1)
      R.col(ax2) = -R.col(ax2);

    if (!std::isfinite(R.norm()))
    {
      std::cout << "Warning: get_quaternion frame from two axes failed!!!" << std::endl;
      return _qt;
    }

    Quaternion q = Quaternion(R);

    return q;
  }


  template<class MeshT>
  static Quaternion construct_quaternion(const Quaternion &_qt, const std::vector<Vec3d> &_directions)
  {
    auto di0 = AxisAlignmentHelpers::closest_axis(_qt, _directions[0]);
    auto di1 = AxisAlignmentHelpers::closest_axis(_qt, _directions[1]);

    return construct_quaternion<MeshT>(_qt, {di0.second / 2, _directions[0]}, {di1.second / 2, _directions[1]});
  }


  static Quaternion project_quaternion(const Quaternion &_qt, const int _axis, const Vec3d &_direction)
  {
    // rotate _q to coordinate system, where _axis_vector is the x-axis
    Quaternion R;
    // get corresponding unit vector
    int vnr = _axis / 2;
    bool negative = _axis % 2;

    Eigen::Vector3d a(0 == vnr, 1 == vnr, 2 == vnr);

    if (negative)
      a = -a;
    R.setFromTwoVectors(a, _direction);
    Quaternion Rip = R.inverse() * _qt;

    // project rotated quaternion and normalize
    Quaternion qp;
    switch (_axis)
    {
      case 0 :
      case 1:qp = Quaternion(Rip.w(), Rip.x(), 0, 0);
        break;
      case 2 :
      case 3:qp = Quaternion(Rip.w(), 0, Rip.y(), 0);
        break;
      case 4 :
      case 5:qp = Quaternion(Rip.w(), 0, 0, Rip.z());
        break;
      default:return _qt; // no projection
    }

    qp.normalize();

    if (!std::isfinite(qp.squaredNorm()))
    {
      std::cout << "Warning: projection of quaternion onto direction failed!!!" << std::endl;
      return _qt;
    }

    // express solution in original coordinate system
    qp = R * qp;

    // choose sign according to dot product between quaternions
    if (qp.dot(_qt) < 0.0)
      qp = Quaternion(-qp.w(), -qp.x(), -qp.y(), -qp.z());

    return qp;
  }


  template<class MeshT>
  static Quaternion project_quaternion(const Quaternion &_qt, const Vec3d &_direction)
  {
    auto di = AxisAlignmentHelpers::closest_axis(_qt, _direction);

    return project_quaternion(_qt, di.second, _direction);
  }


  //align to feature faces, feature edges (and singular edges if added)
  template<class MeshT, class CT>
  static void
  smooth_field_quaternion(const MeshT &mesh, const HFP<int> &trans_prop, const CT &cells,
                          const HFP<bool> &_is_hf_bdy,
                          const CP<std::vector<Vec3d>> &_ch_alignments, const AlgoHex::TransitionQuaternion &tq,
                          CP<Quaternion> &cell_quaternions)
  {
    for (const auto &ch: cells)
    {
      auto hfs = mesh.cell(ch).halffaces();
      Quaternion qt = cell_quaternions[ch];

      if (_ch_alignments[ch].size() == 2)
      {
        cell_quaternions[ch] = construct_quaternion<MeshT>(qt, _ch_alignments[ch]);
        continue;
      }
      else if (_ch_alignments[ch].size() > 2)
        continue;

      for (const auto &hfi: hfs)
      {
        auto hfh_opp = mesh.opposite_halfface_handle(hfi);
        if (!_is_hf_bdy[hfh_opp])
        {
          Quaternion tmp =
                  cell_quaternions[mesh.incident_cell(hfh_opp)] * tq.transition(trans_prop[hfh_opp]);

          if (qt.dot(tmp) >= 0)
          {
            qt.w() += tmp.w();
            qt.x() += tmp.x();
            qt.y() += tmp.y();
            qt.z() += tmp.z();
          }
          else
          {
            qt.w() -= tmp.w();
            qt.x() -= tmp.x();
            qt.y() -= tmp.y();
            qt.z() -= tmp.z();
          }
        }
      }
      qt.normalize();

      if (!std::isfinite(qt.squaredNorm()))
      {
        std::cerr << "ERROR: average quaternion degenerated!!!!" << std::endl;
        qt = cell_quaternions[ch]; // fallback to old quaternion
      }

      if (_ch_alignments[ch].size() == 1)
      {
        cell_quaternions[ch] = project_quaternion<MeshT>(qt, _ch_alignments[ch][0]);
      }
      else if (_ch_alignments[ch].size() == 0)
        cell_quaternions[ch] = qt;
    }
  }


  //no alignment
//    template<class MeshT, class CT>
//    void smooth_field_quaternion(const MeshT &mesh, const HFP<int> &trans_prop, const CT &cells,
//                                 const AlgoHex::TransitionQuaternion &tq, CP<Quaternion> &cell_quaternions) {
//        for (const auto &ch: cells) {
//            auto hfs = mesh.cell(ch).halffaces();
//            Quaternion qt = cell_quaternions[ch];
//
//            for (const auto &hfi: hfs) {
//                auto hfh_opp = mesh.opposite_halfface_handle(hfi);
//                if (!mesh.is_boundary(hfh_opp)) {
//                    Quaternion tmp = cell_quaternions[mesh.incident_cell(hfh_opp)] * tq.transition(trans_prop[hfh_opp]);
//
//                    if (qt.dot(tmp) >= 0) {
//                        qt.w() += tmp.w();
//                        qt.x() += tmp.x();
//                        qt.y() += tmp.y();
//                        qt.z() += tmp.z();
//                    } else {
//                        qt.w() -= tmp.w();
//                        qt.x() -= tmp.x();
//                        qt.y() -= tmp.y();
//                        qt.z() -= tmp.z();
//                    }
//                }
//            }
//            qt.normalize();
//
//            if (!std::isfinite(qt.squaredNorm())) {
//                std::cerr << "ERROR: average quaternion degenerated!!!!" << std::endl;
//                qt = cell_quaternions[ch]; // fallback to old quaternion
//            }
//
//            cell_quaternions[ch] = qt;
//        }
//    }


//    //only align to boundary
//    template<class MeshT, class CT>
//    void smooth_field_quaternion(const MeshT &mesh, const HFP<int> &trans_prop, const CT &cells,
//                                 const std::set<FH> &_bdy_fhs,
//                                 const AlgoHex::TransitionQuaternion &tq, CP<Quaternion> &cell_quaternions) {
//        for (const auto &ch: cells) {
//            auto hfs = mesh.cell(ch).halffaces();
//            Quaternion qt = cell_quaternions[ch];
//
//            HFH bdy_hfh(-1);
//
//            for (const auto &hfi: hfs) {
//                auto hfh_opp = mesh.opposite_halfface_handle(hfi);
//                if (_bdy_fhs.find(mesh.face_handle(hfi)) == _bdy_fhs.end()) {
//                    Quaternion tmp = cell_quaternions[mesh.incident_cell(hfh_opp)] * tq.transition(trans_prop[hfh_opp]);
//
//                    if (qt.dot(tmp) >= 0) {
//                        qt.w() += tmp.w();
//                        qt.x() += tmp.x();
//                        qt.y() += tmp.y();
//                        qt.z() += tmp.z();
//                    } else {
//                        qt.w() -= tmp.w();
//                        qt.x() -= tmp.x();
//                        qt.y() -= tmp.y();
//                        qt.z() -= tmp.z();
//                    }
//                } else {
//                    bdy_hfh = hfh_opp;
//                }
//            }
//            qt.normalize();
//
//            if (!std::isfinite(qt.squaredNorm())) {
//                std::cerr << "ERROR: average quaternion degenerated!!!!" << std::endl;
//                qt = cell_quaternions[ch]; // fallback to old quaternion
//            }
//
//            if (bdy_hfh.is_valid()) {
//                auto normal = mesh.normal(bdy_hfh);
//                cell_quaternions[ch] = project_quaternion<MeshT>(qt, Eigen::Vector3d(normal[0], normal[1], normal[2]));
//            } else
//                cell_quaternions[ch] = qt;
//        }
//    }


  //project to specified axis
  //align to feature faces, feature edges (and singular edges if added)
  template<class MeshT, class CT>
  static void
  smooth_field_quaternion(const MeshT &mesh, const HFP<int> &trans_prop, const CT &cells,
                          const HFP<bool> &_is_hf_bdy,
                          const CP<std::vector<VecAxis>> &_ch_alignments, const AlgoHex::TransitionQuaternion &tq,
                          CP<Quaternion> &cell_quaternions)
  {
    for (const auto &ch: cells)
    {
      auto hfs = mesh.cell(ch).halffaces();
      Quaternion qt = cell_quaternions[ch];

      if (_ch_alignments[ch].size() == 2)
      {
        cell_quaternions[ch] = construct_quaternion<MeshT>(qt, _ch_alignments[ch]);
        continue;
      }
      else if (_ch_alignments[ch].size() > 2)
        continue;

      for (const auto &hfi: hfs)
      {
        auto hfh_opp = mesh.opposite_halfface_handle(hfi);
        if (!_is_hf_bdy[hfh_opp])
        {
          Quaternion tmp =
                  cell_quaternions[mesh.incident_cell(hfh_opp)] * tq.transition(trans_prop[hfh_opp]);

          if (qt.dot(tmp) >= 0)
          {
            qt.w() += tmp.w();
            qt.x() += tmp.x();
            qt.y() += tmp.y();
            qt.z() += tmp.z();
          }
          else
          {
            qt.w() -= tmp.w();
            qt.x() -= tmp.x();
            qt.y() -= tmp.y();
            qt.z() -= tmp.z();
          }
        }
      }
      qt.normalize();

      if (!std::isfinite(qt.squaredNorm()))
      {
        std::cerr << "ERROR: average quaternion degenerated!!!!" << std::endl;
        qt = cell_quaternions[ch]; // fallback to old quaternion
      }

      if (_ch_alignments[ch].size() == 1)
      {
        cell_quaternions[ch] = project_quaternion(qt, _ch_alignments[ch][0].second,
                                                  _ch_alignments[ch][0].first);
      }
      else if (_ch_alignments[ch].size() == 0)
        cell_quaternions[ch] = qt;
    }
  }


  template<class MeshT, class CT>
  static void smooth_field_quaternion(const MeshT &mesh, const HFP<int> &trans_prop, const CT &cells,
                                      const std::set<FH> &_bdy_fhs,
                                      std::map<CH, std::vector<VecAxis>> &_ch_alignments,
                                      const AlgoHex::TransitionQuaternion &tq, CP<Quaternion> &cell_quaternions)
  {
    for (const auto &ch: cells)
    {
      auto hfs = mesh.cell(ch).halffaces();
      Quaternion qt = cell_quaternions[ch];

      if (_ch_alignments[ch].size() == 2)
      {
        cell_quaternions[ch] = construct_quaternion<MeshT>(qt, _ch_alignments[ch]);
        continue;
      }
      else if (_ch_alignments[ch].size() > 2)
        continue;

      for (const auto &hfi: hfs)
      {
        auto hfh_opp = mesh.opposite_halfface_handle(hfi);
        if (_bdy_fhs.find(mesh.face_handle(hfi)) == _bdy_fhs.end())
        {
          Quaternion tmp = cell_quaternions[mesh.incident_cell(hfh_opp)] * tq.transition(trans_prop[hfh_opp]);

          if (qt.dot(tmp) >= 0)
          {
            qt.w() += tmp.w();
            qt.x() += tmp.x();
            qt.y() += tmp.y();
            qt.z() += tmp.z();
          }
          else
          {
            qt.w() -= tmp.w();
            qt.x() -= tmp.x();
            qt.y() -= tmp.y();
            qt.z() -= tmp.z();
          }
        }
      }

      qt.normalize();

      if (!std::isfinite(qt.squaredNorm()))
      {
        std::cerr << "ERROR: average quaternion degenerated!!!!" << std::endl;
        qt = cell_quaternions[ch]; // fallback to old quaternion
      }

      if (_ch_alignments[ch].size() == 1)
      {
        cell_quaternions[ch] = project_quaternion(qt, _ch_alignments[ch][0].second,
                                                  _ch_alignments[ch][0].first);
      }
      else if (_ch_alignments[ch].size() == 0)
        cell_quaternions[ch] = qt;
    }
  }

//    //smooth cell qtns at edge
//    template<class MeshT, class CT>
//    void smooth_field_quaternion_at_edge(const MeshT &mesh, const HFP<int> &trans_prop, const CT &cells,
//                                 const std::set<FH> &_bdy_fhs,
//                                 std::map<CH, std::vector<VecAxis>> &_ch_alignments,
//                                 const AlgoHex::TransitionQuaternion &tq, CP<Quaternion> &cell_quaternions) {
//        for (const auto &ch: cells) {
//            auto hfs = mesh.cell(ch).halffaces();
//            Quaternion qt = cell_quaternions[ch];
//
//            if (_ch_alignments[ch].size() == 2) {
//                cell_quaternions[ch] = construct_quaternion<MeshT>(qt, _ch_alignments[ch]);
//                continue;
//            } else if (_ch_alignments[ch].size() > 2)
//                continue;
//
//
//            for (const auto &hfi: hfs) {
//                auto hfh_opp = mesh.opposite_halfface_handle(hfi);
//                if (_bdy_fhs.find(mesh.face_handle(hfi)) == _bdy_fhs.end()) {
//                    CH ch_opp = mesh.incident_cell(hfh_opp);
//                    if(cells.find(ch_opp) == cells.end())
//                        continue;
//
//                    Quaternion tmp = cell_quaternions[mesh.incident_cell(hfh_opp)] * tq.transition(trans_prop[hfh_opp]);
//
//                    if (qt.dot(tmp) >= 0) {
//                        qt.w() += tmp.w();
//                        qt.x() += tmp.x();
//                        qt.y() += tmp.y();
//                        qt.z() += tmp.z();
//                    } else {
//                        qt.w() -= tmp.w();
//                        qt.x() -= tmp.x();
//                        qt.y() -= tmp.y();
//                        qt.z() -= tmp.z();
//                    }
//                }
//            }
//
//            qt.normalize();
//
//            if (!std::isfinite(qt.squaredNorm())) {
//                std::cerr << "ERROR: average quaternion degenerated!!!!" << std::endl;
//                qt = cell_quaternions[ch]; // fallback to old quaternion
//            }
//
//            if (_ch_alignments[ch].size() == 1) {
//                cell_quaternions[ch] = project_quaternion(qt, _ch_alignments[ch][0].second,
//                                                          _ch_alignments[ch][0].first);
//            } else if (_ch_alignments[ch].size() == 0)
//                cell_quaternions[ch] = qt;
//        }
//    }


  //
  template<class MeshT>
  static void
  optimize_quaternions(const MeshT &mesh, const HFP<int> &trans_prop, const AlgoHex::TransitionQuaternion &tq,
                       CP<Quaternion> &cell_quaternions, const EP<int> &valence, const FP<int> &feature_face,
                       const EP<int> &feature_edge, const std::set<CH> &_cells,
                       std::map<CH, int> &_bdy_aligned_axis, const bool _align_sg = false)
  {
    std::set<FH> bdy_fhs_set;
    for (const auto &ch: _cells)
      for (auto cf_it = mesh.cf_iter(ch); cf_it.valid(); ++cf_it)
      {
        if (mesh.is_boundary(*cf_it))
          bdy_fhs_set.insert(*cf_it);
      }

    std::map<CH, std::vector<VecAxis >> cell_alignments;
    //align to feature face
    for (auto&[ch, ax]: _bdy_aligned_axis)
    {
      for (auto chf_it = mesh.chf_iter(ch); chf_it.valid(); ++chf_it)
      {
        auto fhi = mesh.face_handle(*chf_it);
        if (feature_face[fhi] > 0)
        {
          auto nm = mesh.normal(*chf_it);
          cell_alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), ax);
        }
      }
    }

    //align to feature edges and singular edges
    std::set<CH> visited_chs;
    for (const auto &ch: _cells)
    {
      if (visited_chs.find(ch) != visited_chs.end())
        continue;

      for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it)
      {
        if ((_align_sg && (valence[*ce_it] >= -2 && valence[*ce_it] <= 4 && valence[*ce_it] != 0)) ||
            feature_edge[*ce_it] > 0)
        {
          HEH heh0 = mesh.halfedge_handle(*ce_it, 0);
          auto dir = mesh.vertex(mesh.halfedge(heh0).to_vertex()) -
                     mesh.vertex(mesh.halfedge(heh0).from_vertex());
          dir.normalize();
          Vec3d eigen_dir = Vec3d(dir[0], dir[1], dir[2]);

          auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells(mesh, tq,
                                                                                         cell_quaternions,
                                                                                         trans_prop,
                                                                                         valence,
                                                                                         feature_edge,
                                                                                         heh0);

          for (auto&[chi, axi]: cell_edge_axis)
          {
            visited_chs.insert(chi);
            cell_alignments[chi].emplace_back(eigen_dir, axi);
          }
        }
      }
    }

    //check
//        for(const auto ch : _cells) {
//            if (cell_alignments[ch].size() > 2)
//                std::cout << "Warning: cell " << ch << " has " << cell_alignments[ch].size() << " alignment constraints!"
//                          << std::endl;
//            else if (cell_alignments[ch].size() == 2)
//                std::cout << "cell " << ch << " has " << cell_alignments[ch].size() << " alignment constraints!"
//                          << std::endl;
//        }

    for (int i = 0; i < 20; ++i)
    {
      smooth_field_quaternion(mesh, trans_prop, _cells, bdy_fhs_set, cell_alignments, tq, cell_quaternions);
    }
  }

//    template<class MeshT>
//    void optimize_quaternions(const MeshT &mesh, const HFP<int> &trans_prop, const AlgoHex::TransitionQuaternion &tq,
//                              CP<Quaternion> &cell_quaternions,
//                              const EP<int> &valence, const EP<int> &feature_edge, const std::set<CH> &_cells,
//                              const bool _align_sg = false) {
//        std::set<FH> bdy_fhs_set;
//        for (const auto &ch: _cells)
//            for (auto cf_it = mesh.cf_iter(ch); cf_it.valid(); ++cf_it) {
//                if (mesh.is_boundary(*cf_it))
//                    bdy_fhs_set.insert(*cf_it);
//            }
//
//        std::map<CH, std::vector<VecAxis >> cell_alignments;
//        //align to boundary
//        for (const auto &ch: _cells) {
//            if (mesh.is_boundary(ch)) {
//                for (auto chf_it = mesh.chf_iter(ch); chf_it.valid(); ++chf_it) {
//                    auto hfh_opp = mesh.opposite_halfface_handle(*chf_it);
//                    if (mesh.is_boundary(hfh_opp)) {
//                        auto nm = mesh.normal(hfh_opp);
//                        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch], ovm2eigen(nm)).second;
//                        cell_alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
//                    }
//                }
//            }
//        }
//
//        //align to feature edges and singular edges
//        std::set<CH> visited_chs;
//        for (const auto &ch: _cells) {
//            if (visited_chs.find(ch) != visited_chs.end())
//                continue;
//
//            for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it) {
//                if((_align_sg && (valence[*ce_it] >= -2 && valence[*ce_it] <= 4 && valence[*ce_it] != 0)) || feature_edge[*ce_it] > 0) {
//                    HEH heh0 = mesh.halfedge_handle(*ce_it, 0);
//                    auto dir = mesh.vertex(mesh.halfedge(heh0).to_vertex()) - mesh.vertex(mesh.halfedge(heh0).from_vertex());
//                    dir.normalize();
//                    Vec3d eigen_dir = Vec3d(dir[0], dir[1], dir[2]);
//
//                    auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells(mesh, tq, cell_quaternions, trans_prop, valence, feature_edge, heh0);
//
//                    for(auto& [chi, axi] : cell_edge_axis) {
//                        visited_chs.insert(chi);
//                        cell_alignments[chi].emplace_back(eigen_dir, axi);
//                    }
//                }
//            }
//        }
//
////        std::set<CH> visited_chs;
////        for (const auto &ch: _cells) {
////            if (visited_chs.find(ch) != visited_chs.end())
////                continue;
////
////            for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it) {
////                if((_align_sg || feature_edge[*ce_it]) && (valence[*ce_it] == 1 || valence[*ce_it] == -1 || std::abs(valence[*ce_it]) == 2)) {
////                    HEH heh0 = mesh.halfedge_handle(*ce_it, 0);
////                    auto dir = mesh.vertex(mesh.halfedge(heh0).to_vertex()) -
////                               mesh.vertex(mesh.halfedge(heh0).from_vertex());
////                    dir.normalize();
////
////                    auto hehf_it = mesh.hehf_iter(heh0);
////                    CH ch_s(-1);
////                    bool e_bdy = mesh.is_boundary(*ce_it);
////                    if (!e_bdy)
////                        ch_s = mesh.incident_cell(mesh.opposite_halfface_handle(*hehf_it));
////                    else
////                        ch_s = mesh.incident_cell(*hehf_it);
////
////                    int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh, tq, cell_quaternions, trans_prop, valence, heh0,
////                                                            *hehf_it);
////                    if (valence[*ce_it] == 1)
////                        axis = axis % 2 == 0 ? axis + 1 : axis - 1;
////
////                    if(std::abs(valence[*ce_it]) == 2) {
////                        Vec3d ax_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions[ch_s], (AxisAlignment)axis);
////                        if(ax_vec.dot(ovm2eigen(dir)) < 0.)
////                            axis = axis%2 == 0 ? axis+1 : axis-1;
////                    }
////
////                    cell_alignments[ch_s].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
////                    visited_chs.insert(ch_s);
////
////                    if (e_bdy)
////                        ++hehf_it;
////
////                    for (; hehf_it.valid(); ++hehf_it) {
////                        auto chi = mesh.incident_cell(*hehf_it);
////                        if (chi.is_valid()) {
////                            if (chi == ch_s)
////                                break;
////
////                            axis = tq.axis_after_transition(axis, trans_prop[*hehf_it]);
////
////                            cell_alignments[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
////                            visited_chs.insert(chi);
////                        }
////                    }
////                } else if(feature_edge[*ce_it] && valence[*ce_it] == 0) {
////                    auto dir = mesh.vertex(mesh.edge(*ce_it).to_vertex()) - mesh.vertex(mesh.edge(*ce_it).from_vertex());
////                    dir.normalize();
////
////                    int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch], ovm2eigen(dir)).second;
////                    cell_alignments[ch].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
////                    visited_chs.insert(ch);
////                }
////            }
////        }
//
//        //check
////        for(const auto ch : _cells) {
////            if (alignment[ch].second > 2)
////                std::cout << "Warning: cell " << ch << " has " << alignment[ch].second << " alignment constraints!"
////                          << std::endl;
//////            else if (alignment[ch].second == 2)
//////                std::cout << "cell " << ch << " has " << alignment[ch].second << " alignment constraints!"
//////                          << std::endl;
////        }
//
//        for (int i = 0; i < 20; ++i) {
//            smooth_field_quaternion(mesh, trans_prop, _cells, bdy_fhs_set, cell_alignments, tq, cell_quaternions);
//        }
//    }

  template<class MeshT>
  static void optimize_quaternions_in_cells(const MeshT &mesh, const HFP<int> &trans_prop,
                                            const AlgoHex::TransitionQuaternion &tq,
                                            CP<Quaternion> &cell_quaternions,
                                            const EP<int> &valence, const FP<int> &feature_face,
                                            const EP<int> &feature_edge,
                                            const std::set<CH> &_cells, const bool _align_sg = false)
  {
    std::set<FH> bdy_fhs_set;
    for (const auto &ch: _cells)
      for (auto cf_it = mesh.cf_iter(ch); cf_it.valid(); ++cf_it)
      {
        if (mesh.is_boundary(*cf_it))
          bdy_fhs_set.insert(*cf_it);
      }

    std::map<CH, std::vector<VecAxis >> cell_alignments;
    //align to feature faces
    for (const auto &ch: _cells)
    {
      for (auto chf_it = mesh.chf_iter(ch); chf_it.valid(); ++chf_it)
      {
        auto hfh_opp = mesh.opposite_halfface_handle(*chf_it);
        if (feature_face[mesh.face_handle(hfh_opp)] > 0)
        {
          auto nm = mesh.normal(hfh_opp);
          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch], ovm2eigen(nm)).second;
          cell_alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
        }
      }
    }

    //align to feature edges and singular edges
    std::set<CH> visited_chs;
    for (const auto &ch: _cells)
    {
      if (visited_chs.find(ch) != visited_chs.end())
        continue;

      for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it)
      {
        if ((_align_sg && (valence[*ce_it] >= -2 && valence[*ce_it] <= 4 && valence[*ce_it] != 0)) ||
            feature_edge[*ce_it] > 0)
        {
          HEH heh0 = mesh.halfedge_handle(*ce_it, 0);
          auto dir = mesh.vertex(mesh.halfedge(heh0).to_vertex()) -
                     mesh.vertex(mesh.halfedge(heh0).from_vertex());
          dir.normalize();
          Vec3d eigen_dir = Vec3d(dir[0], dir[1], dir[2]);

          auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells(mesh, tq,
                                                                                         cell_quaternions,
                                                                                         trans_prop,
                                                                                         valence,
                                                                                         feature_edge,
                                                                                         heh0);

          for (auto&[chi, axi]: cell_edge_axis)
          {
            visited_chs.insert(chi);
            cell_alignments[chi].emplace_back(eigen_dir, axi);
          }
        }
      }
    }

//        std::set<CH> visited_chs;
//        for (const auto &ch: _cells) {
//            if (visited_chs.find(ch) != visited_chs.end())
//                continue;
//
//            for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it) {
//                if((_align_sg || feature_edge[*ce_it]) && (valence[*ce_it] == 1 || valence[*ce_it] == -1 || std::abs(valence[*ce_it]) == 2)) {
//                    HEH heh0 = mesh.halfedge_handle(*ce_it, 0);
//                    auto dir = mesh.vertex(mesh.halfedge(heh0).to_vertex()) -
//                               mesh.vertex(mesh.halfedge(heh0).from_vertex());
//                    dir.normalize();
//
//                    auto hehf_it = mesh.hehf_iter(heh0);
//                    CH ch_s(-1);
//                    bool e_bdy = mesh.is_boundary(*ce_it);
//                    if (!e_bdy)
//                        ch_s = mesh.incident_cell(mesh.opposite_halfface_handle(*hehf_it));
//                    else
//                        ch_s = mesh.incident_cell(*hehf_it);
//
//                    int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh, tq, cell_quaternions, trans_prop, valence, heh0,
//                                                            *hehf_it);
//                    if (valence[*ce_it] == 1)
//                        axis = axis % 2 == 0 ? axis + 1 : axis - 1;
//
//                    if(std::abs(valence[*ce_it]) == 2) {
//                        Vec3d ax_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions[ch_s], (AxisAlignment)axis);
//                        if(ax_vec.dot(ovm2eigen(dir)) < 0.)
//                            axis = axis%2 == 0 ? axis+1 : axis-1;
//                    }
//
//                    cell_alignments[ch_s].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
//                    visited_chs.insert(ch_s);
//
//                    if (e_bdy)
//                        ++hehf_it;
//
//                    for (; hehf_it.valid(); ++hehf_it) {
//                        auto chi = mesh.incident_cell(*hehf_it);
//                        if (chi.is_valid()) {
//                            if (chi == ch_s)
//                                break;
//
//                            axis = tq.axis_after_transition(axis, trans_prop[*hehf_it]);
//
//                            cell_alignments[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
//                            visited_chs.insert(chi);
//                        }
//                    }
//                } else if(feature_edge[*ce_it] && valence[*ce_it] == 0) {
//                    auto dir =
//                            mesh.vertex(mesh.edge(*ce_it).to_vertex()) - mesh.vertex(mesh.edge(*ce_it).from_vertex());
//                    dir.normalize();
//
//                    auto da = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch], ovm2eigen(dir));
//                    int axis = da.second;
//                    cell_alignments[ch].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
//                    visited_chs.insert(ch);
//                }
//            }
//        }

    //check
//        for(const auto ch : _cells) {
//            if (cell_alignments[ch].size() > 2)
//                std::cout << "Warning: cell " << ch << " has " << cell_alignments[ch].size() << " alignment constraints!"
//                          << std::endl;
//            else if (cell_alignments[ch].size() == 2)
//                std::cout << "cell " << ch << " has " << cell_alignments[ch].size() << " alignment constraints!"
//                          << std::endl;
//        }
//        for(auto& [chi, va] : cell_alignments)
//            for(auto vai : va)
//                std::cerr<<"ch "<<chi<<" ax "<<vai.second<<" dir "<<vai.first<<std::endl;

    for (int i = 0; i < 20; ++i)
    {
      smooth_field_quaternion(mesh, trans_prop, _cells, bdy_fhs_set, cell_alignments, tq, cell_quaternions);
    }
  }


  template<class MeshT>
  static void optimize_quaternions_wrt_prescribed_valence(const MeshT &mesh, const HFP<int> &trans_prop,
                                                          const AlgoHex::TransitionQuaternion &tq,
                                                          CP<Quaternion> &cell_quaternions,
                                                          const EP<int> &valence, const EP<int> &feature_edge,
                                                          const FP<int> &feature_face, const std::set<CH> &_cells,
                                                          const bool _align_sg = false)
  {
    std::set<FH> bdy_fhs_set;
    for (const auto ch: _cells)
      for (auto cf_it = mesh.cf_iter(ch); cf_it.valid(); ++cf_it)
      {
        if (mesh.is_boundary(*cf_it))
//      if (feature_face[*cf_it])
          bdy_fhs_set.insert(*cf_it);
      }

    std::map<CH, std::vector<VecAxis >> cell_alignments;
    //align to boundary
    for (const auto ch: _cells)
    {
      for (auto cf_it = mesh.cf_iter(ch); cf_it.valid(); ++cf_it)
      {
//      if (mesh.is_boundary(*cf_it) || feature_face[*cf_it]) {
        if (feature_face[*cf_it] > 0)
        {
          auto hfh0 = mesh.halfface_handle(*cf_it, 0);
          if (ch != mesh.incident_cell(hfh0))
            hfh0 = mesh.opposite_halfface_handle(hfh0);

          auto nm = mesh.normal(hfh0);
          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch], ovm2eigen(nm)).second;
          cell_alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
        }
      }
    }

    std::set<CH> visited_chs;
    for (const auto ch: _cells)
    {
      if (visited_chs.find(ch) != visited_chs.end())
        continue;

      for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it)
      {
        if ((_align_sg || feature_edge[*ce_it]) &&
            (valence[*ce_it] == 1 || valence[*ce_it] == -1 || std::abs(valence[*ce_it]) == 2))
        {
          HEH heh0 = mesh.halfedge_handle(*ce_it, 0);
          auto dir = mesh.vertex(mesh.halfedge(heh0).to_vertex()) -
                     mesh.vertex(mesh.halfedge(heh0).from_vertex());
          dir.normalize();

          auto hehf_it = mesh.hehf_iter(heh0);
          CH ch_s(-1);
          bool e_bdy = mesh.is_boundary(*ce_it);
          if (!e_bdy)
            ch_s = mesh.incident_cell(mesh.opposite_halfface_handle(*hehf_it));
          else
            ch_s = mesh.incident_cell(*hehf_it);

          int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh, tq, cell_quaternions,
                                                                        trans_prop, valence, heh0,
                                                                        *hehf_it);
          if (valence[*ce_it] == 1)
            axis = axis % 2 == 0 ? axis + 1 : axis - 1;

          if (std::abs(valence[*ce_it]) == 2)
          {
            Vec3d ax_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions[ch_s],
                                                                   (AxisAlignment) axis);
            if (ax_vec.dot(ovm2eigen(dir)) < 0.)
              axis = axis % 2 == 0 ? axis + 1 : axis - 1;
          }

          cell_alignments[ch_s].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
          visited_chs.insert(ch_s);

          if (e_bdy)
            ++hehf_it;

          for (; hehf_it.valid(); ++hehf_it)
          {
            auto chi = mesh.incident_cell(*hehf_it);
            if (chi.is_valid())
            {
              if (chi == ch_s)
                break;

              axis = tq.axis_after_transition(axis, trans_prop[*hehf_it]);

              cell_alignments[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
              visited_chs.insert(chi);
            }
          }
        }
        else if (feature_edge[*ce_it] && valence[*ce_it] == 0)
        {
          auto dir = mesh.vertex(mesh.edge(*ce_it).to_vertex()) -
                     mesh.vertex(mesh.edge(*ce_it).from_vertex());
          dir.normalize();

          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch], ovm2eigen(dir)).second;
          cell_alignments[ch].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
          visited_chs.insert(ch);
        }
      }
    }

    for (int i = 0; i < 20; ++i)
    {
      smooth_field_quaternion(mesh, trans_prop, _cells, bdy_fhs_set, cell_alignments, tq, cell_quaternions);
    }
  }

  template<class MeshT>
  static void optimize_quaternions_wrt_prescribed_valence(const MeshT &mesh, const HFP<int> &trans_prop,
                                                          const AlgoHex::TransitionQuaternion &tq,
                                                          CP<Quaternion> &cell_quaternions,
                                                          const EP<int> &valence, const EP<int> &feature_edge,
                                                          const FP<int> &feature_face, const VH _vh_new,
                                                          const bool _align_sg = false)
  {
    std::set<CH> cells;
    for (auto vc_it = mesh.vc_iter(_vh_new); vc_it.valid(); ++vc_it)
      cells.insert(*vc_it);

    optimize_quaternions_wrt_prescribed_valence(mesh, trans_prop, tq, cell_quaternions, valence, feature_edge,
                                                feature_face, cells, _align_sg);
  }

  //
  template<class MeshT>
  static void optimize_quaternions_at_edge(const MeshT &mesh, const HFP<int> &trans_prop,
                                           const AlgoHex::TransitionQuaternion &tq,
                                           CP<Quaternion> &cell_quaternions, const EP<int> &valence,
                                           const FP<int> &feature_face, const EP<int> &feature_edge, const EH eh)
  {
    std::set<CH> e_chs;
    for (auto ec_it = mesh.ec_iter(eh); ec_it.valid(); ++ec_it)
    {
      e_chs.insert(*ec_it);
    }

    std::set<FH> bdy_fhs_set;
    for (const auto &ch: e_chs)
      for (auto cf_it = mesh.cf_iter(ch); cf_it.valid(); ++cf_it)
      {
        if (mesh.is_boundary(*cf_it))
          bdy_fhs_set.insert(*cf_it);
      }

    std::map<CH, std::vector<VecAxis >> cell_alignments;
    //align to feature faces
//        for (const auto &ch: e_chs) {
//            for (auto chf_it = mesh.chf_iter(ch); chf_it.valid(); ++chf_it) {
//                auto hfh_opp = mesh.opposite_halfface_handle(*chf_it);
//                if (feature_face[mesh.face_handle(hfh_opp)] > 0) {
//                    auto nm = mesh.normal(hfh_opp);
//                    int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions[ch], ovm2eigen(nm)).second;
//                    cell_alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
//                }
//            }
//        }

    if ((valence[eh] >= -2 && valence[eh] <= 4 && valence[eh] != 0) || feature_edge[eh] > 0)
    {
      //align to edge
      HEH heh0 = mesh.halfedge_handle(eh, 0);
      auto dir = mesh.vertex(mesh.halfedge(heh0).to_vertex()) -
                 mesh.vertex(mesh.halfedge(heh0).from_vertex());
      dir.normalize();
      Vec3d eigen_dir = Vec3d(dir[0], dir[1], dir[2]);

      auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells(mesh, tq,
                                                                                     cell_quaternions,
                                                                                     trans_prop, valence,
                                                                                     feature_edge, heh0);
//std::cerr<<"opt qtns at eh "<<eh<<" ";
      for (auto&[chi, axi]: cell_edge_axis)
      {
//std::cerr<<"     ch "<<chi<<" "<<axi<<std::endl;
        cell_alignments[chi].emplace_back(eigen_dir, axi);
      }
//std::cerr<<"opt qtns at eh "<<eh<<std::endl;
//
//            for(auto& [chi, axi] : cell_edge_axis) {
//std::cerr<<"     chi "<<chi<<" line "<<mesh.barycenter(chi)<<" "<< (ovm2eigen(mesh.barycenter(chi)) + AxisAlignmentHelpers::quaternion_vector(cell_quaternions[chi], (AxisAlignment)axi)).transpose()<<" "<<axi<<std::endl;
//            }
    }

    for (int i = 0; i < 20; ++i)
    {
//            smooth_field_quaternion_at_edge(mesh, trans_prop, e_chs, bdy_fhs_set, cell_alignments, tq, cell_quaternions);
      smooth_field_quaternion(mesh, trans_prop, e_chs, bdy_fhs_set, cell_alignments, tq, cell_quaternions);

    }
  }


  template<class MeshT>
  static void optimize_quaternions_at_vertex(const MeshT &mesh, const HFP<int> &trans_prop,
                                             const AlgoHex::TransitionQuaternion &tq,
                                             CP<Quaternion> &cell_quaternions,
                                             const EP<int> &valence, const FP<int> &feature_face,
                                             const EP<int> &feature_edge,
                                             const VH _vh_new, const bool _align_sg = false)
  {
    std::set<CH> cells;
    for (auto vc_it = mesh.vc_iter(_vh_new); vc_it.valid(); ++vc_it)
      cells.insert(*vc_it);

    optimize_quaternions_in_cells(mesh, trans_prop, tq, cell_quaternions, valence, feature_face, feature_edge,
                                  cells,
                                  _align_sg);
  }


//    template<class MeshT>
//    void optimize_quaternions(const MeshT &mesh, const HFP<int> &trans_prop, const AlgoHex::TransitionQuaternion &tq,
//                              CP<Quaternion> &cell_quaternions,
//                              const EP<int> &valence, const EP<int> &feature_edge, const VH _vh_new,
//                              const bool _align_sg = false) {
//        std::set<CH> cells;
//        for (auto vc_it = mesh.vc_iter(_vh_new); vc_it.valid(); ++vc_it)
//            cells.insert(*vc_it);
//
//        optimize_quaternions(mesh, trans_prop, tq, cell_quaternions, valence, feature_edge, cells, _align_sg);
//    }

  template<class MeshT>
  static void optimize_quaternions_at_edges(const MeshT &mesh, const HFP<int> &trans_prop,
                                            const AlgoHex::TransitionQuaternion &tq,
                                            CP<Quaternion> &cell_quaternions,
                                            const EP<int> &valence, const FP<int> &feature_face,
                                            const EP<int> &feature_edge, const std::set<EH> &_edges,
                                            const bool _align_sg = false)
  {
    std::set<VH> vhs;
    for (const auto &e: _edges)
    {
      vhs.insert(mesh.edge(e).from_vertex());
      vhs.insert(mesh.edge(e).to_vertex());
    }

    std::set<CH> chs;
    for (const auto &vh: vhs)
    {
      for (auto vc_it = mesh.vc_iter(vh); vc_it.valid(); ++vc_it)
        chs.insert(*vc_it);
    }

//        std::cerr<<"opt qtn chs: ";
//        for(const auto chi : chs)
//            std::cerr<<" "<<chi;
//        std::cerr<<std::endl;

    optimize_quaternions_in_cells(mesh, trans_prop, tq, cell_quaternions, valence, feature_face, feature_edge,
                                  chs,
                                  _align_sg);
  }

  template<class MeshT>
  static void optimize_quaternions_with_alignment_constraints(const MeshT &mesh, const HFP<int> &trans_prop,
                                                              const AlgoHex::TransitionQuaternion &tq,
                                                              CP<Quaternion> &cell_quaternions,
                                                              const std::set<CH> &_cells,
                                                              std::map<CH, std::vector<VecAxis> > &_alignments,
                                                              const int _iter = 20)
  {
    std::set<FH> bdy_fhs_set;
    for (const auto &ch: _cells)
      for (auto cf_it = mesh.cf_iter(ch); cf_it.valid(); ++cf_it)
      {
        if (mesh.is_boundary(*cf_it))
          bdy_fhs_set.insert(*cf_it);
      }


    for (int i = 0; i < _iter; ++i)
      smooth_field_quaternion(mesh, trans_prop, _cells, bdy_fhs_set, _alignments, tq, cell_quaternions);
  }

  template<class MeshT>
  static void optimize_quaternions_with_alignment_constraints(const MeshT &mesh, const HFP<int> &trans_prop,
                                                              const AlgoHex::TransitionQuaternion &tq,
                                                              CP<Quaternion> &cell_quaternions,
                                                              const VH _vh_new,
                                                              std::map<CH, std::vector<VecAxis> > &_alignments,
                                                              const int _iter = 20)
  {
    std::set<CH> cells;
    for (auto vc_it = mesh.vc_iter(_vh_new); vc_it.valid(); ++vc_it)
      cells.insert(*vc_it);

    optimize_quaternions_with_alignment_constraints(mesh, trans_prop, tq, cell_quaternions, cells, _alignments,
                                                    _iter);
  }
};
}

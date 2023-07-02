/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define FRAMEFIELDOPTIMIZATIONT_C

#include "FrameFieldOptimizationT.hh"
#include "EdgeMonodromyHelperT.hh"
#include "CommonFuncs.hh"


namespace AlgoHex
{

template<class MeshT>
void
FrameFieldOptimizationT<MeshT>::
align_quaternions_to_feature(const int _iterations)
{
  std::cout << "Preprocessing quaternions..." << std::endl;

  auto is_hfh_bdy = mesh_.template request_halfface_property<bool>("is_boundary_halfface", false);
  for (const auto hfh: mesh_.halffaces())
    if (mesh_.is_boundary(hfh))
      is_hfh_bdy[hfh] = true;

  auto alignment = mesh_.template request_cell_property<std::vector<VecAxis> >("cell_alignment");

  bool unclear_alignment = collect_alignments(alignment);

  if (unclear_alignment)
  {
    std::cerr
            << "Initial alignment constraints at some feature are not clear. Smooth field without these constraints to remove ambiguity..."
            << std::endl;
    for (int i = 0; i < 3; ++i)
    {
      QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, mesh_.cells(), is_hfh_bdy, alignment, tq_,
                                                   cell_quaternions_);
    }

    unclear_alignment = collect_alignments(alignment);
    if (unclear_alignment)
      std::cerr << "SERIOUS WARNING: unclear alignment constraints!!!" << std::endl;
  }


  for (int i = 0; i < _iterations; ++i)
  {
    QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, mesh_.cells(), is_hfh_bdy, alignment, tq_,
                                                 cell_quaternions_);
  }

  check_feature_face_alignments_of_field();
  std::cerr << "Done!\n";
}

template<class MeshT>
void
FrameFieldOptimizationT<MeshT>::
optimize_quaternions_local(const int k, const int _iterations)
{
  AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::sgaf_field_opt_local);

  std::cout << "Optimizing quaternions..." << std::endl;

  auto cells = get_k_ring_cells_of_singular_graph(mesh_, k + 1);
//        std::cout<<"Get "<<k+1<<" ring cells ("<<cells.size()<<") takes: "<<sw.stop()/1000.<<"s."<<std::endl;


  auto is_hfh_bdy = mesh_.template request_halfface_property<bool>("is_boundary_halfface", false);
  for (const auto hfh: mesh_.halffaces())
    if (mesh_.is_boundary(hfh))
      is_hfh_bdy[hfh] = true;

  auto alignment = mesh_.template request_cell_property<std::vector<VecAxis> >("cell_alignment");

  //align to normal direction
  for (const auto ch: mesh_.cells())
  {
    for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0)
      {
        auto nm = mesh_.normal(*chf_it);

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], nm).second;
        alignment[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
      }
    }

    if (alignment[ch].size() > 1)
      std::cout << "Warning: cell has " << alignment[ch].size() << " boundary alignment!" << std::endl;
  }


  //align to feature edges
  for (const auto ehi: mesh_.edges())
  {
    if (feature_edge_[ehi] > 0)
    {
      HEH heh0 = mesh_.halfedge_handle(ehi, 0);
      auto dir = mesh_.vertex(mesh_.halfedge(heh0).to_vertex()) -
                 mesh_.vertex(mesh_.halfedge(heh0).from_vertex());
      dir.normalize();

      if (valence_[ehi] >= -2 && valence_[ehi] <= 4 && valence_[ehi] != 0)
      {
        HEH heh0 = mesh_.halfedge_handle(ehi, 0);
        auto dir = mesh_.vertex(mesh_.halfedge(heh0).to_vertex()) - mesh_.vertex(mesh_.halfedge(heh0).from_vertex());
        dir.normalize();

        auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells(mesh_, tq_, cell_quaternions_,
                                                                                       trans_prop_, valence_,
                                                                                       feature_edge_, heh0);

        for (auto&[chi, axi]: cell_edge_axis)
          alignment[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axi);
      }
      else
      {
        for (auto hec_it = mesh_.hec_iter(heh0); hec_it.valid(); ++hec_it)
        {
          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*hec_it],
                                                        dir).second;
          alignment[*hec_it].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
        }
      }
    }
  }

  for (int i = 0; i < _iterations; ++i)
  {
    QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, cells, is_hfh_bdy, alignment, tq_,
                                                 cell_quaternions_);
  }
//        std::cout<<"Quaternion smoothing takes: "<<sw.stop()/1000.<<"s."<<std::endl;
}


template<class MeshT>
void
FrameFieldOptimizationT<MeshT>::
optimize_quaternions_global(const int _iterations)
{
  AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::sgaf_field_opt_global);

  std::cout << "Optimize quaternion field with all alignments ..." << std::endl;
  auto is_hfh_bdy = mesh_.template request_halfface_property<bool>("is_boundary_halfface", false);
  for (const auto hfh: mesh_.halffaces())
    if (mesh_.is_boundary(hfh))
      is_hfh_bdy[hfh] = true;

  auto alignment = mesh_.template request_cell_property<std::vector<VecAxis> >("cell_alignment");

//        auto flipped_valence = mesh_.template request_edge_property<bool>("flipped valence");


  //align frames to feature faces
  for (const auto ch: mesh_.cells())
  {
    for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0)
      {
        auto nm = mesh_.normal(*chf_it);

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], nm).second;
        alignment[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
      }
    }
  }

  //align frames to feature edges and singular edges
  for (const auto ehi: mesh_.edges())
  {
    if (valence_[ehi] >= -2 && valence_[ehi] <= 4 && valence_[ehi] != 0)
    {
      HEH heh0 = mesh_.halfedge_handle(ehi, 0);
      auto dir = mesh_.vertex(mesh_.halfedge(heh0).to_vertex()) - mesh_.vertex(mesh_.halfedge(heh0).from_vertex());
      dir.normalize();

      auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells_wrt_valence(mesh_, tq_,
                                                                                                 cell_quaternions_,
                                                                                                 trans_prop_, valence_,
                                                                                                 feature_edge_, heh0);

      for (auto&[chi, axi]: cell_edge_axis)
        alignment[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axi);
    }
    else if (feature_edge_[ehi] > 0)
    {
      HEH heh0 = mesh_.halfedge_handle(ehi, 0);
      auto dir = mesh_.vertex(mesh_.halfedge(heh0).to_vertex()) - mesh_.vertex(mesh_.halfedge(heh0).from_vertex());
      dir.normalize();

      for (auto hec_it = mesh_.hec_iter(heh0); hec_it.valid(); ++hec_it)
      {
        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*hec_it], dir).second;
        alignment[*hec_it].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
      }
    }
  }

  //check
  for (const auto ch: mesh_.cells())
  {
    if (alignment[ch].size() > 2)
    {
      std::cout << "Warning: cell " << ch << " has " << alignment[ch].size()
                << " alignment constraints! They are dir:";
      for (auto alg: alignment[ch])
      {
        std::cerr << alg.first[0] << " " << alg.first[1] << " " << alg.first[2] << " ax: " << alg.second;
      }
      std::cerr << std::endl;
    }
  }

  for (int i = 0; i < _iterations; ++i)
  {
    QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, mesh_.cells(), is_hfh_bdy, alignment, tq_,
                                                 cell_quaternions_);
  }


  std::cerr << "Checking alignment error: \n";
  check_feature_face_alignments_of_field();
  check_edge_alignments_of_field();
  std::cerr << "Done!" << std::endl;

//        mesh_.set_persistent(flipped_valence, false);
}

template<class MeshT>
void
FrameFieldOptimizationT<MeshT>::
optimize_quaternions_global(const double _angle_thr)
{
  //translate
  auto angle_thr = _angle_thr / 180.0 * M_PI;

  AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::sgaf_field_opt_global);

  std::cout << "Optimize quaternion field with all alignments(adaptive) ..." << std::endl;
  auto is_hfh_bdy = mesh_.template request_halfface_property<bool>("is_boundary_halfface", false);
  for (const auto hfh: mesh_.halffaces())
    if (mesh_.is_boundary(hfh))
      is_hfh_bdy[hfh] = true;

  auto alignment = mesh_.template request_cell_property<std::vector<VecAxis> >("cell_alignment");

//        auto flipped_valence = mesh_.template request_edge_property<bool>("flipped valence");


  //align to ff
  for (const auto ch: mesh_.cells())
  {
    for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0)
      {
        auto nm = mesh_.normal(*chf_it);

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], nm).second;
        alignment[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), axis);
      }
    }
  }

  //align to fe and sge
  for (const auto ehi: mesh_.edges())
  {
    if (valence_[ehi] >= -2 && valence_[ehi] <= 4 && valence_[ehi] != 0)
    {
      HEH heh0 = mesh_.halfedge_handle(ehi, 0);
      auto dir = mesh_.vertex(mesh_.halfedge(heh0).to_vertex()) - mesh_.vertex(mesh_.halfedge(heh0).from_vertex());
      dir.normalize();

      auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells_wrt_valence(mesh_, tq_,
                                                                                                 cell_quaternions_,
                                                                                                 trans_prop_, valence_,
                                                                                                 feature_edge_, heh0);

      for (auto&[chi, axi]: cell_edge_axis)
      {
        alignment[chi].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axi);
      }
    }
    else if (feature_edge_[ehi] > 0)
    {
      HEH heh0 = mesh_.halfedge_handle(ehi, 0);
      auto dir = mesh_.vertex(mesh_.halfedge(heh0).to_vertex()) - mesh_.vertex(mesh_.halfedge(heh0).from_vertex());
      dir.normalize();

      for (auto hec_it = mesh_.hec_iter(heh0); hec_it.valid(); ++hec_it)
      {
        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*hec_it], dir).second;
        alignment[*hec_it].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
      }
    }
  }

  //check
  for (const auto ch: mesh_.cells())
  {
    if (alignment[ch].size() > 2)
    {
      std::cout << "Warning: cell " << ch << " has " << alignment[ch].size()
                << " alignment constraints! They are dir:";
      for (auto alg: alignment[ch])
      {
        std::cerr << alg.first[0] << " " << alg.first[1] << " " << alg.first[2] << " ax: " << alg.second;
      }
      std::cerr << std::endl;
    }
  }


  std::priority_queue<DC> que;

  for (auto chi: mesh_.cells())
    if (!alignment[chi].empty())
      que.push(std::make_pair(angle_thr, chi));

  int n = 0;
  while (!que.empty())
  {
    auto dc_cur = que.top();
    que.pop();

    Quaternion qt_old = cell_quaternions_[dc_cur.second];

    QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, std::array<CH, 1>{dc_cur.second}, is_hfh_bdy,
                                                 alignment, tq_, cell_quaternions_);


    Quaternion rt = cell_quaternions_[dc_cur.second] * qt_old.conjugate();
    Eigen::AngleAxis<double> aa(rt);

    double angle = aa.angle();
//            std::cerr<<"  ch "<<dc_cur.second<<" angle "<<angle<<" angle_thr "<<angle_thr<<std::endl;

    if (angle >= angle_thr)
    {
//                std::cerr<<"  ch "<<dc_cur.second<<" angle "<<angle<<" angle_thr "<<angle_thr<<std::endl;
      for (auto cc_it = mesh_.cc_iter(dc_cur.second); cc_it.valid(); ++cc_it)
        que.push(std::make_pair(angle, *cc_it));
    }

    n++;
  }


  check_feature_face_alignments_of_field();
  check_edge_alignments_of_field();

//        mesh_.set_persistent(flipped_valence, false);
}

template<class MeshT>
void
FrameFieldOptimizationT<MeshT>::check_feature_face_alignments_of_field()
{
  std::cerr << "Checking face alignment error: \n";
  for (auto fh: mesh_.faces())
  {
    if (feature_fprop_[fh] > 0)
    {
      // get non-boundary halffaces
      HFH hfh0 = mesh_.halfface_handle(fh, 0);
      HFH hfh1 = mesh_.halfface_handle(fh, 1);
      std::vector<HFH> hfhs;
      if (!mesh_.is_boundary(hfh0))
        hfhs.push_back(hfh0);
      if (!mesh_.is_boundary(hfh1))
        hfhs.push_back(hfh1);

      for (size_t i = 0; i < hfhs.size(); ++i)
      {
        HFH hfh = hfhs[i];
        // get corresponding cell
        CH ch = mesh_.incident_cell(hfh);

        // get vertices
        VH vh0 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[0]).to_vertex();
        VH vh1 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[1]).to_vertex();
        VH vh2 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[2]).to_vertex();

        // get points
        Point p0 = mesh_.vertex(vh0);
        Point p1 = mesh_.vertex(vh1);
        Point p2 = mesh_.vertex(vh2);

        // get outward normal vector
        Vec3d n = acg2eigen((p1 - p0) % (p2 - p0));
        Vec3d e0 = acg2eigen(p1 - p0);
        Vec3d e1 = acg2eigen(p2 - p0);


        // get frame of cell
        Mat3d F = cell_quaternions_[ch].toRotationMatrix();

        // determine alignment axis
        Mat3d J = F.inverse();
        Vec3d Je0 = J * e0;
        Vec3d Je1 = J * e1;
        Vec3d Jn = Je0.cross(Je1);
        AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Jn);
        int axis_idx = int(e_axis) / 2;
        // test whether field is sufficiently aligned
        Vec3d a = AxisAlignmentHelpers::vector(e_axis);
        if (a.dot(Jn) < Jn.norm() * 0.99985)
          std::cerr << "Warning: feature face alignment - face " << fh << " alignment axis deviates by "
                    << std::acos(a.dot(Jn) / Jn.norm()) * 180.0 / M_PI
                    << " degree\n";
      }
    }
  }
  std::cerr << "Done!\n";
}

template<class MeshT>
void
FrameFieldOptimizationT<MeshT>::check_edge_alignments_of_field()
{
  std::cerr << "Checking edge alignment error: \n";
  for (const auto ehi: mesh_.edges())
  {
    // singular edge?
    if (valence_[ehi] != 0 || feature_edge_[ehi] > 0)
    {
      // get all cells incident to the edge
      std::vector<CH> chs;
      for (HECIt hec_it = mesh_.hec_iter(mesh_.halfedge_handle(ehi, 0)); hec_it.valid(); ++hec_it)
        chs.push_back(*hec_it);

      for (size_t i = 0; i < chs.size(); ++i)
      {
        // get curent cell adjacent to edge
        CH ch = chs[i];

        // get adjacent vertices
        VH vh0 = mesh_.edge(ehi).from_vertex();
        VH vh1 = mesh_.edge(ehi).to_vertex();

        // get edge vector
        Vec3d e = acg2eigen(mesh_.vertex(vh1) - mesh_.vertex(vh0));

        // determine alignment axis
        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], e).second;
        Vec3d axis_vec = AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch], (AxisAlignment) axis);

        if (e.dot(axis_vec) < e.norm() * 0.99985)
          std::cerr << "Warning FF alignment: edge " << ehi << " alignment axis deviates by "
                    << std::acos(e.dot(axis_vec) / e.norm()) * 180.0 / M_PI
                    << " degree\n";
      }
    }
  }
  std::cerr << "Done!\n";
}

template<class MeshT>
bool
FrameFieldOptimizationT<MeshT>::
collect_alignments(CP <std::vector<VecAxis>> &_alignments)
{
  for (const auto chi: mesh_.cells())
    _alignments[chi].clear();

  bool unclear_alignment = false;
  //align to feature faces
  for (const auto ch: mesh_.cells())
  {
    for (auto chf_it = mesh_.chf_iter(ch); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0)
      {
        auto nm = mesh_.normal(*chf_it);

        auto dax = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], nm);
        if (std::fabs(dax.first - std::sqrt(2.) / 2.) < 1e-6)
        {
          std::cerr << "Warning: axis aligned to feature face may be ambiguous. Cell " << ch << std::endl;
          unclear_alignment = true;
          continue;
        }
        _alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), dax.second);
      }
    }

    if (_alignments[ch].size() > 1)
      std::cout << "Warning: cell " << ch << " has " << _alignments[ch].size() << " feature face alignment!"
                << std::endl;
  }


  //align to feature edges
  for (const auto ch: mesh_.cells())
  {
    for (auto ce_it = mesh_.ce_iter(ch); ce_it.valid(); ++ce_it)
    {
      if (feature_edge_[*ce_it] > 0)
      {
        auto dir = mesh_.vertex(mesh_.edge(*ce_it).to_vertex()) - mesh_.vertex(mesh_.edge(*ce_it).from_vertex());
        dir.normalize();

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[ch], dir).second;

        //check conflict with feature face alignment
        if (is_feature_cell(ch) && _alignments[ch][0].second / 2 == axis / 2)
        {
          std::cout << "Warning: cell " << ch << " has " << _alignments[ch].size()
                    << " same alignments. Feature face axis: " << _alignments[ch][0].second << " edge axis: " << axis
                    << std::endl;
          unclear_alignment = true;
          continue;
        }

        _alignments[ch].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
      }
    }

    if (_alignments[ch].size() > 2)
      std::cout << "Warning: cell " << ch << " has " << _alignments[ch].size() << " alignments!" << std::endl;
    else if (_alignments[ch].size() == 2)
    {
      if (_alignments[ch][0].second / 2 == _alignments[ch][1].second / 2)
      {
        std::cout << "Warning: cell " << ch << " has " << _alignments[ch].size()
                  << " same alignments, clear the constraints!" << std::endl;
        _alignments[ch].clear();
      }
    }
  }

  return unclear_alignment;
}


template<class MeshT>
void
FrameFieldOptimizationT<MeshT>::create_closest_field_axis(MeshT &_pm)
{
  //clear the old mesh
  // the mesh.clear() causes crash!!!
  for (auto eh: _pm.edges())
  {
    if (!_pm.is_deleted(eh))
    {
      _pm.delete_edge(eh);
    }
  }
  for (auto vh: _pm.vertices())
  {
    if (!_pm.is_deleted(vh))
    {
      _pm.delete_vertex(vh);
    }
  }
  _pm.collect_garbage();
  for (const auto &eh: mesh_.edges())
  {
    if (valence_[eh] >= -2 && valence_[eh] <= 4 && valence_[eh] != 0 && !mesh_.is_boundary(eh))
    {
      auto vhf = mesh_.edge(eh).from_vertex();
      auto mid_pt = mesh_.vertex(vhf) + mesh_.vertex(mesh_.edge(eh).to_vertex());
      mid_pt /= 2.;

      HEH heh0 = mesh_.halfedge_handle(eh, 0);
      auto cell_edge_axis = EdgeMonodromyHelperT::get_halfedge_aligned_axis_in_cells(mesh_, tq_, cell_quaternions_,
                                                                                     trans_prop_, valence_,
                                                                                     feature_edge_, heh0);

      Vec3d dir_avr;
      dir_avr.setZero();
      for (auto&[chi, axi]: cell_edge_axis)
      {
        dir_avr += AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[chi], (AxisAlignment) axi);
      }

      auto vm = _pm.add_vertex(mid_pt);
      auto v1 = _pm.add_vertex(mid_pt + Point(dir_avr[0], dir_avr[1], dir_avr[2]) * target_length_[vhf] / 2);
      auto v0 = _pm.add_vertex(mid_pt - Point(dir_avr[0], dir_avr[1], dir_avr[2]) * target_length_[vhf] / 2);
      _pm.add_edge(vm, v1);
      _pm.add_edge(vm, v0);
    }
  }
}

template<class MeshT>
bool FrameFieldOptimizationT<MeshT>::
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

}
/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define SINGULARGRAPHEXTRACTIONT_C

//== INCLUDES =================================================================

#include "SingularGraphExtractionT.hh"
#include "MeshGeometry.hh"
#include "AxisAlignment.hh"
#include "Geometry.hh"
#include "./LocallyMeshableField/QuaternionsSmoothing.hh"

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== IMPLEMENTATION ==========================================================

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
get_singular_edges_with_refinement(const double _largest_length_ratio)
{
  ScopedStopWatch sw_total(sw::sg_extraction);

  refine_singular_regions(_largest_length_ratio);

  get_singular_edges();
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
get_singular_edges()
{

  {
    ScopedStopWatch sw(sw::interpolate_cq);
    get_field_per_cell();
  }

  {
    ScopedStopWatch sw(sw::extract_sg);
    get_transitions();

    get_edges_valence();
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
refine_singular_regions(const double _largest_length_ratio)
{
  std::cout << "\n#####Generating singular graph... \nsplitting (edge length ratio: " << _largest_length_ratio
            << ")..." << std::endl;

  ScopedStopWatch sw(sw::refinement);

  refine_interior_regions(_largest_length_ratio * av_edge_length_);
  mesh_.collect_garbage();
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
save_cell_quaternions(const std::string &_filename) const
{
  std::ofstream f_write(_filename);

  for (const auto ch: mesh_.cells())
  {
    // write to file
    f_write << cell_quaternion_[ch].w() << " ";
    f_write << cell_quaternion_[ch].x() << " ";
    f_write << cell_quaternion_[ch].y() << " ";
    f_write << cell_quaternion_[ch].z() << " ";
  }

  f_write.close();
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
get_edges_valence()
{
  //reset the valence
  for (const auto eh: mesh_.edges())
    compute_edge_valence(eh);
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
compute_edge_valence(const EH _eh)
{
  valence_[_eh] = calculate_edge_valence(_eh);
}

//-----------------------------------------------------------------------------

template<class MeshT>
int SingularGraphExtractionT<MeshT>::
calculate_edge_valence(const EH _eh)
{
  bool is_bdy = mesh_.is_boundary(_eh);

  HEH he0 = mesh_.halfedge_handle(_eh, 0);
  int idx = compute_edge_transition(he0, is_bdy);
  int old_val = valence_[_eh];

  if (!is_bdy)
  {
    HFH hf_start = *mesh_.hehf_iter(he0);

    if (idx != 0)
    {
      CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_start));

      Point vn = mesh_.vertex(mesh_.halfedge(he0).to_vertex()) -
                 mesh_.vertex(mesh_.halfedge(he0).from_vertex());

      return compute_valence(he0, cell_quaternion_[ch_s], idx, vn);
    }
    else
    {//consider higher valence
      double us_rt_angle = unsigned_rotation_angle_at_edge(he0);
      if (us_rt_angle > M_PI)
      {
        int old_val = valence_[_eh];
        //suppose it's singular of high valence
        valence_[_eh] = 4;

        //store old quaterions
        std::map<CH, Quaternion> old_cell_quaternions;
        for (auto hec_it = mesh_.hec_iter(he0); hec_it.valid(); ++hec_it)
        {
          old_cell_quaternions[*hec_it] = cell_quaternion_[*hec_it];
        }

        QuaternionSmoothing::optimize_quaternions_at_edge(mesh_, trans_prop_, tq_, cell_quaternion_, valence_,
                                                          feature_fprop_, feature_edge_prop_, _eh);

        double sec_val = calc_sector_angle(he0, *mesh_.hehf_iter(he0));

        int valence = 0;
        if (sec_val != -1000)
          valence = sec_val - 4;
        else
        { // undertermined
          valence = old_val;
        }

        //recover qtn
        for (auto hec_it = mesh_.hec_iter(he0); hec_it.valid(); ++hec_it)
          cell_quaternion_[*hec_it] = old_cell_quaternions[*hec_it];

        valence_[_eh] = old_val;

        return valence;
      }

      return 0;
    }
  }
  else
  {
    HFH hf_start = *mesh_.hehf_iter(he0);
    //should be boundary halfface if oriented
    HFH hf_end = mesh_.opposite_halfface_handle(*mesh_.hehf_iter(mesh_.opposite_halfedge_handle(he0)));

    CH ch0 = mesh_.incident_cell(hf_start);
    Point n0 = mesh_.normal(mesh_.opposite_halfface_handle(hf_start));
    int closest_n0 = AxisAlignmentHelpers::closest_axis(cell_quaternion_[ch0], n0).second;

    CH ch1 = mesh_.incident_cell(mesh_.opposite_halfface_handle(hf_end));
    Point n1 = mesh_.normal(hf_end);
    int closest_n1 = AxisAlignmentHelpers::closest_axis(cell_quaternion_[ch1], n1).second;
    int T_inv_id = tq_.inverse_transition_idx(idx);

    int n1_in_q0 = tq_.axis_after_transition(closest_n1, T_inv_id);

    if (closest_n0 != n1_in_q0)
    {
      if (closest_n0 / 2 == n1_in_q0 / 2)
      {
        auto da = dihedral_angle(mesh_, he0);
        if (da < M_PI)
          return -2;
        else if (da >= M_PI)
          return 2;
      }
      else
      {
        int third_axis_id = the_third_axis((AxisAlignment) closest_n0, (AxisAlignment) n1_in_q0);
        Vec3d third_axis = AxisAlignmentHelpers::quaternion_vector(cell_quaternion_[ch0],
                                                                   (AxisAlignment) third_axis_id);

        VH vh0 = mesh_.halfedge(he0).from_vertex();
        VH vh1 = mesh_.halfedge(he0).to_vertex();
        Point vn = mesh_.vertex(vh0) - mesh_.vertex(vh1);

        double sign = Point(third_axis(0), third_axis(1), third_axis(2)) | vn;
        if (sign < -0.)
          return 1;
        else if (sign > 0.)
          return -1;
        else
        {
          std::cout << "Warning: boundary edge direction is orthogonal to rotational axis!" << std::endl;
          return std::numeric_limits<int>::lowest();
        }
      }
    }
    else
      return 0;
  }

  return 0;
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
load_cell_quaternions(const std::string &_filename)
{
  std::ifstream fin;
  fin.open(_filename);

  for (const auto &ch: mesh_.cells())
  {
    fin >> cell_quaternion_[ch].w();
    fin >> cell_quaternion_[ch].x();
    fin >> cell_quaternion_[ch].y();
    fin >> cell_quaternion_[ch].z();
  }

  fin.close();
}

//=====================================================================================================================//

template<class MeshT>
double SingularGraphExtractionT<MeshT>::
unsigned_rotation_angle_at_edge(const HEH _heh) const
{
  double energy = 0.;
  for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
  {
    auto fhi = mesh_.face_handle(*hehf_it);
    energy += std::abs(frame_smoothness_energy_at_face(fhi));
  }

  return energy;
}

//=====================================================================================================================//

template<class MeshT>
int SingularGraphExtractionT<MeshT>::
get_dominant_axis_in_cell(const HEH _heh) const
{
  Point dir = mesh_.vertex(mesh_.halfedge(_heh).to_vertex()) -
              mesh_.vertex(mesh_.halfedge(_heh).from_vertex());
  dir.normalize();
  Vec3d vn = ovm2eigen(dir);

  std::vector<double> dprds(3, 0.);
  for (int i = 0; i < 3; ++i)
  {
    int axi = 2 * i;
    for (auto hehf_it = mesh_.hehf_iter(_heh); hehf_it.valid(); ++hehf_it)
    {
      CH chi = mesh_.incident_cell(*hehf_it);
      dprds[i] += vn.dot(AxisAlignmentHelpers::quaternion_vector(cell_quaternion_[chi], AxisAlignment(axi)));
      //axis in the next cell
      axi = tq_.axis_after_transition(axi, trans_prop_[*hehf_it]);
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
    max_ax = max_ax % 2 == 0 ? max_ax + 1 : max_ax - 1;

  return max_ax;
}

//=====================================================================================================================//

template<class MeshT>
double SingularGraphExtractionT<MeshT>::
frame_smoothness_energy_at_face(const FH _fh) const
{
  double energy = 0.;
  if (mesh_.is_boundary(_fh))
    return energy;

  auto hfh0 = mesh_.halfface_handle(_fh, 0);
  auto hfh1 = mesh_.opposite_halfface_handle(hfh0);

  auto ch0 = mesh_.incident_cell(hfh0);
  auto ch1 = mesh_.incident_cell(hfh1);

  Quaternion tranq = tq_.transition(trans_prop_[hfh0]);

  Quaternion q0in1 = cell_quaternion_[ch0] * tranq;

  Quaternion rt = cell_quaternion_[ch1] * q0in1.conjugate();

  energy += 2 * std::acos(std::min(std::abs(rt.w()), 1.));

  return energy;
}

//=====================================================================================================================//

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
get_field_per_cell()
{
  std::cout << "interpolating for cell quaternions..." << std::endl;

  int num_cells = mesh_.n_cells();
#pragma omp parallel for schedule(dynamic, 16)
  for (int i = 0; i < num_cells; ++i)
  {
    CH ch(i);
    auto cell_vhs = mesh_.get_cell_vertices(ch);
    SHCoeffs qi;
    qi.setZero();

    for (int i = 0; i < 4; ++i)
      qi += shc_[cell_vhs[i]];

    qi.normalize();
    cell_quaternion_[ch] = projector_.project(qi).q;
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
get_transitions()
{
  std::cout << "compute transitions..." << std::endl;

  for (const auto fh: mesh_.faces())
    compute_face_transition(fh);
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
compute_face_transition(const FH _fh)
{
  if (mesh_.is_boundary(_fh))
  {
    trans_prop_[mesh_.halfface_handle(_fh, 0)] = -1;
    trans_prop_[mesh_.halfface_handle(_fh, 1)] = -1;
  }
  else
  {
    HFH hfh0 = mesh_.halfface_handle(_fh, 0);
    HFH hfh1 = mesh_.halfface_handle(_fh, 1);

    int best_trans(0);
    find_best_transition(hfh0, best_trans);
    trans_prop_[hfh0] = best_trans;
    trans_prop_[hfh1] = tq_.inverse_transition_idx(best_trans);
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
int SingularGraphExtractionT<MeshT>::
compute_edge_transition(const HEH _heh, bool _is_boundary) const
{
  int idx = 0;
  if (!_is_boundary)
  {
    HFH hf_start = *mesh_.hehf_iter(_heh);
    HFH hf_it = hf_start;
    do
    {
      if (hf_it.idx() < 0)
      {
        std::cerr << "Error: invalid halfface " << hf_it << " at edge " << mesh_.edge_handle(_heh) << std::endl;
        return -1;
      }
      idx = tq_.mult_transitions_idx(trans_prop_[hf_it], idx);

      //counterclockwise
      hf_it = mesh_.adjacent_halfface_in_cell(hf_it, _heh);
      hf_it = mesh_.opposite_halfface_handle(hf_it);
    }
    while (hf_it != hf_start);
  }
  else
  {
    //find start halfface
    HFH hf_start = *mesh_.hehf_iter(_heh);

    //compute transition product
    HFH hf_it = hf_start;
    hf_it = mesh_.adjacent_halfface_in_cell(hf_it, _heh);
    hf_it = mesh_.opposite_halfface_handle(hf_it);

    while (!mesh_.is_boundary(hf_it) && hf_it.is_valid())
    {
      if (hf_it.idx() < 0)
      {
        std::cerr << "Error: invalid halfface " << hf_it << " at boundary edge " << mesh_.edge_handle(_heh)
                  << std::endl;
        return -1;
      }

      idx = tq_.mult_transitions_idx(trans_prop_[hf_it], idx);

      //counterclockwise
      hf_it = mesh_.adjacent_halfface_in_cell(hf_it, _heh);
      hf_it = mesh_.opposite_halfface_handle(hf_it);
    }
  }

  return idx;
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SingularGraphExtractionT<MeshT>::
refine_interior_regions(const double _edge_length)
{
  int n_split_e = 0, n_split_f = 0;

  std::queue<FH> que;
  for (const auto fhi: mesh_.faces())
    que.push(fhi);

  int max_split = 20000;
  while (!que.empty() && (n_split_e + n_split_f) < max_split)
  {
    auto fh_cur = que.front();
    que.pop();

    if (mesh_.is_deleted(fh_cur))
      continue;

    auto fvhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(fh_cur, 0));

    int idx = compute_tri_idx(quaternion_[fvhs[0]], quaternion_[fvhs[1]], quaternion_[fvhs[2]]);
    if (idx == 0)
      continue;

    EH longest_edge(-1);
    if (longest_edge_length(fh_cur, longest_edge) > _edge_length)
    {
      auto vh_m = split_edge(longest_edge);

      //push incident faces
      for (auto vf_it = mesh_.vf_iter(vh_m); vf_it.valid(); ++vf_it)
        que.push(*vf_it);

      n_split_e++;
    }
    else
    {
      split_face(fh_cur, fvhs);
      n_split_f++;
    }
  }

  std::cout << "performed " << n_split_e << " edge split" << std::endl;
  std::cout << "Performed " << n_split_f << " face split" << std::endl;
}

//-----------------------------------------------------------------------------

template<class MeshT>
double SingularGraphExtractionT<MeshT>::
find_best_transition(const HFH hfh, int &_best_transition) const
{
  // error checking
  if (mesh_.is_boundary(hfh))
  {
    std::cout << "Error: find best transition received boundary facet!!!\n";
    return -1;
  }

  CH ch1 = mesh_.incident_cell(hfh);
  CH ch0 = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh));

  // find closest transition quaternion
  Quaternion t = cell_quaternion_[ch1].conjugate() * cell_quaternion_[ch0];
  _best_transition = tq_.closest_transition_idx(t);
  double measure = 1.0 - fabs(t.dot(tq_.transition(_best_transition)));

  return measure;
}

//-----------------------------------------------------------------------------

template<class MeshT>
int SingularGraphExtractionT<MeshT>::
compute_valence(const HEH _heh, const Quaternion &_q, int _idx, const Point &_vn)
{
  int valence = 0;
  Mat3d f = _q.toRotationMatrix();
  if (_idx < 10 && _idx > 3)
  {
    Point xaxis(f(0, 0), f(1, 0), f(2, 0));
    Point yaxis(f(0, 1), f(1, 1), f(2, 1));
    Point zaxis(f(0, 2), f(1, 2), f(2, 2));
    Point taxis;
    int rotation_axis = tq_.rotation_axis_2d(_idx);
    if (rotation_axis == 0) taxis = xaxis;
    else if (rotation_axis == 1) taxis = -xaxis;
    else if (rotation_axis == 2) taxis = yaxis;
    else if (rotation_axis == 3) taxis = -yaxis;
    else if (rotation_axis == 4) taxis = zaxis;
    else if (rotation_axis == 5) taxis = -zaxis;

    double sign = _vn | taxis;
    if (sign == 0.0)
    {
      std::cout << "Warning: edge direction is orthogonal to rotational axis!" << std::endl;
      valence = std::numeric_limits<int>::lowest();
    }
    else
    {
      //store old quaterions
      std::map<CH, Quaternion> old_cell_quaternions;
      for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
      {
        old_cell_quaternions[*hec_it] = cell_quaternion_[*hec_it];
      }

      EH eh_cur = mesh_.edge_handle(_heh);
      int old_val = valence_[eh_cur];
      valence_[eh_cur] = 1;
      //align quaternions to singular edge
      QuaternionSmoothing::optimize_quaternions_at_edge(mesh_, trans_prop_, tq_, cell_quaternion_, valence_,
                                                        feature_fprop_, feature_edge_prop_, eh_cur);

      double sec_val = calc_sector_angle(_heh, *mesh_.hehf_iter(_heh));

      if (sec_val != -1000)
      {
        valence = sec_val - 4;
        if (valence == 0)
        {
          std::cerr << "Error: non-identity matching product corresponds to valence 0 at edge "
                    << mesh_.edge_handle(_heh) << std::endl;
          valence = std::numeric_limits<int>::lowest();
        }
      }
      else // undertermined
        valence = std::numeric_limits<int>::lowest();

      //recover qtn
      for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
        cell_quaternion_[*hec_it] = old_cell_quaternions[*hec_it];

      valence_[eh_cur] = old_val;

    }
  }
  else if (_idx <= 3 && _idx > 0)
  {
    //store old quaterions
    std::map<CH, Quaternion> old_cell_quaternions;
    for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
    {
      old_cell_quaternions[*hec_it] = cell_quaternion_[*hec_it];
    }

    EH eh_cur = mesh_.edge_handle(_heh);
    int old_val = valence_[eh_cur];
    valence_[eh_cur] = 2;

    //align quaternions to singular edge
    QuaternionSmoothing::optimize_quaternions_at_edge(mesh_, trans_prop_, tq_, cell_quaternion_, valence_,
                                                      feature_fprop_, feature_edge_prop_, eh_cur);

    double sec_val = calc_sector_angle(_heh, *mesh_.hehf_iter(_heh));

    if (sec_val != -1000)
      valence = sec_val - 4;
    else // undertermined
      valence = 2;

    //recover qtn
    for (auto hec_it = mesh_.hec_iter(_heh); hec_it.valid(); ++hec_it)
      cell_quaternion_[*hec_it] = old_cell_quaternions[*hec_it];

    valence_[eh_cur] = old_val;
  }
  else if (_idx >= 10 && _idx < 24)
    valence = std::numeric_limits<int>::max();

  return valence;
}

template<class TetMeshT>
double
SingularGraphExtractionT<TetMeshT>::
calc_sector_angle(const HEH _heh, const HFH _hfh_sector_start) const
{
  HFH hfh0 = _hfh_sector_start;
  FH fh0 = mesh_.face_handle(hfh0);
  CH ch0 = mesh_.incident_cell(hfh0);

  // get adjacent vertices
  VH vh0 = mesh_.halfedge(_heh).from_vertex();
  VH vh1 = mesh_.halfedge(_heh).to_vertex();

  // global coordinate
  Vec3d e = ovm2eigen(mesh_.vertex(vh1) - mesh_.vertex(vh0));
  Vec3d u, v, w;
  complement_to_right_handed_orthonormal_frame(e, u, v, w);

  // determine alignment axis
  Mat3d F0 = cell_quaternion_[ch0].toRotationMatrix();
  Mat3d J0 = F0.inverse();
  Vec3d e_trans = J0 * e;
  AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(e_trans);

  // get two tangent vectors
  AxisAlignment a0_axis = AxisAlignment((int(e_axis) + 2) % 6);
  AxisAlignment b0_axis = AxisAlignment((int(a0_axis) + 2) % 6);

  //correct (e, a0, b0) to right-handed frame if necessary, if sign(e)==-1 -> swap sign of a0
  if (!AxisAlignmentHelpers::is_positive_axis(e_axis))
    a0_axis = reverse_axis(a0_axis);

  // param domain vectors of axes for rotation estimation
  Vec3d a0_trans = AxisAlignmentHelpers::vector(a0_axis);
  Vec3d b0_trans = AxisAlignmentHelpers::vector(b0_axis);

  Vec3d e0_trans = AxisAlignmentHelpers::vector(e_axis);


  double rotation_angle_a(0.0);
  double rotation_angle_b(0.0);
  double dihedral_angle(0.0);

  // collect data for debugging
  std::vector<double> d_angles;
  std::vector<double> ra_angles;
  std::vector<double> rb_angles;
  std::vector<int> tfs;

  int trans_prd = 0;
  for (int i = 0; i < 10000; ++i) // restrict to 10000 sectors to avoid infinite loops
  {
    // add up dihedral angle
    double da = dihedral_angle_in_cell(mesh_, _heh, hfh0);
    dihedral_angle += da;
    // debug helper
    d_angles.push_back(da);

    // get next halfface
    HFH hfh1 = mesh_.adjacent_halfface_in_cell(hfh0, _heh);
    HFH hfh1_opp = mesh_.opposite_halfface_handle(hfh1);
    FH fh1 = mesh_.face_handle(hfh1_opp);
    if (mesh_.is_boundary(fh1))
      break;
    else
    {
      // debug helper
      tfs.push_back(trans_prop_[hfh1]);

      CH ch1 = mesh_.incident_cell(hfh1_opp);

      Vec3d e00 = F0 * e0_trans;


      // transform F0 to coordinate system of F1
      Mat3d F1 = cell_quaternion_[ch1].toRotationMatrix();
      F0 = (F0 * tq_.transition_matrix_int(trans_prop_[hfh1])).eval();
      // transform axes a and b
      auto TM = tq_.transition_matrix_int(trans_prop_[hfh1_opp]); // inverse transition
      a0_trans = (TM * a0_trans).eval();
      b0_trans = (TM * b0_trans).eval();

      e0_trans = (TM * e0_trans).eval();

      // get corresponding frame vectors
      Vec3d a0 = F0 * a0_trans;
      Vec3d b0 = F0 * b0_trans;
      Vec3d a1 = F1 * a0_trans;
      Vec3d b1 = F1 * b0_trans;

      Vec3d e0 = F0 * e0_trans;
      Vec3d e1 = F1 * e0_trans;


      // accumulate 2d rotation angles around edge
      double ra = signed_angle(Vec2d(a0.dot(v), a0.dot(w)), Vec2d(a1.dot(v), a1.dot(w)));
      double rb = signed_angle(Vec2d(b0.dot(v), b0.dot(w)), Vec2d(b1.dot(v), b1.dot(w)));

      rotation_angle_a += ra;
      rotation_angle_b += rb;

      ra_angles.push_back(ra);
      rb_angles.push_back(rb);

      if (hfh1_opp == _hfh_sector_start)
        break;

      // update variables for next step
      std::swap(F0, F1);
      std::swap(hfh0, hfh1_opp);
    }
  }

  double sector_angle = dihedral_angle - rotation_angle_a;
  double sector_valence = std::round(sector_angle / (0.5 * M_PI));

  double sector_angle_b = dihedral_angle - rotation_angle_b;
  double sector_valence_b = std::round(sector_angle_b / (0.5 * M_PI));

  if (sector_valence != sector_valence_b)
  {
    std::cerr << "Warning: sector angles are not equal at edge " << mesh_.edge_handle(_heh) << " sec val a "
              << sector_valence << " b " << sector_valence_b << std::endl;
    return -1000;
  }

  return sector_valence;
}

//-----------------------------------------------------------------------------

template<class MeshT>
inline int SingularGraphExtractionT<MeshT>::
compute_tri_idx(const Quaternion &q0, const Quaternion &q1, const Quaternion &q2) const
{
  return tq_.mult_transitions_idx(tq_.closest_transition_idx(q0.conjugate() * q2),
                                  tq_.mult_transitions_idx(tq_.closest_transition_idx(q2.conjugate() * q1),
                                                           tq_.closest_transition_idx(q1.conjugate() * q0)));
}

//-----------------------------------------------------------------------------

template<class MeshT>
double SingularGraphExtractionT<MeshT>::
longest_edge_length(const FH _fh, EH &_eh) const
{
  double max_length = -1.;
  for (auto fe_it = mesh_.fe_iter(_fh); fe_it.valid(); ++fe_it)
  {
    auto length = mesh_.length(*fe_it);
    if (max_length < length)
    {
      max_length = length;
      _eh = *fe_it;
    }
  }

  return max_length;
}

//-----------------------------------------------------------------------------

template<class MeshT>
OVM::VertexHandle
SingularGraphExtractionT<MeshT>::
split_edge(const EH _eh)
{
  HEH heh = mesh_.halfedge_handle(_eh, 0);

  VH vh0 = mesh_.halfedge(heh).from_vertex();
  VH vh1 = mesh_.halfedge(heh).to_vertex();

  int feature_edge_id = feature_edge_prop_[_eh];

  //store feature face vertices
  std::vector<std::vector<VH>> hfs_vhs;
  std::vector<int> ffs_id;
  for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
  {
    int ff_id = feature_fprop_[mesh_.face_handle(*hehf_it)];
    if (ff_id > 0)
    {
      hfs_vhs.push_back(mesh_.get_halfface_vertices(*hehf_it));
      ffs_id.push_back(ff_id);
    }
  }

  VH vh_new = mesh_.split_edge(_eh);

  //interpolate sph and get quaternion of the new vertex
  shc_[vh_new] = (shc_[vh0] + shc_[vh1]).normalized();

  quaternion_[vh_new] = projector_.project(shc_[vh_new]).q;

  //update feature edge property
  EH eh0 = mesh_.edge_handle(mesh_.find_halfedge(vh0, vh_new));
  feature_edge_prop_[eh0] = feature_edge_id;

  EH eh1 = mesh_.edge_handle(mesh_.find_halfedge(vh1, vh_new));
  feature_edge_prop_[eh1] = feature_edge_id;

  feature_edge_vertex_[vh_new] = feature_edge_id > 0;

  //update feature face property
  for (auto i = 0u; i < hfs_vhs.size(); ++i)
  {
    for (auto &vhi: hfs_vhs[i])
    {
      if (vhi == vh0)
      {
        vhi = vh_new;

        //update
        auto hfh0 = mesh_.find_halfface(hfs_vhs[i]);
        feature_fprop_[mesh_.face_handle(hfh0)] = ffs_id[i];

        for (auto hfe_it = mesh_.hfe_iter(hfh0); hfe_it.valid(); ++hfe_it)
        {
          feature_face_edge_[*hfe_it] = true;

          feature_face_vertex_[mesh_.edge(*hfe_it).from_vertex()] = true;
          feature_face_vertex_[mesh_.edge(*hfe_it).to_vertex()] = true;
        }

        //reset
        vhi = vh0;
      }

      if (vhi == vh1)
      {
        vhi = vh_new;
        auto hfh1 = mesh_.find_halfface(hfs_vhs[i]);

        //update
        feature_fprop_[mesh_.face_handle(hfh1)] = ffs_id[i];
        for (auto hfe_it = mesh_.hfe_iter(hfh1); hfe_it.valid(); ++hfe_it)
        {
          feature_face_edge_[*hfe_it] = true;

          feature_face_vertex_[mesh_.edge(*hfe_it).from_vertex()] = true;
          feature_face_vertex_[mesh_.edge(*hfe_it).to_vertex()] = true;
        }

        //reset
        vhi = vh1;
      }
    }
  }

  return vh_new;
}

//-----------------------------------------------------------------------------

template<class MeshT>
OVM::VertexHandle
SingularGraphExtractionT<MeshT>::
split_face(const FH _fh, const std::vector<VH> &_fvhs)
{
  //halfface vertices
  auto hf_vhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(_fh, 0));

  int ff_id = feature_fprop_[_fh];


  VH vh_new = mesh_.split_face(_fh);

  shc_[vh_new] = (shc_[_fvhs[0]] + shc_[_fvhs[1]] + shc_[_fvhs[2]]).normalized();

  quaternion_[vh_new] = projector_.project(shc_[vh_new]).q;

  //update feature face property
  if (ff_id > 0)
  {
    for (int i = 0; i < 3; ++i)
    {
      //replace the vertex vha with the new vertex
      auto vha = hf_vhs[i];
      hf_vhs[i] = vh_new;

      auto hfhi = mesh_.find_halfface(hf_vhs);

      //update
      feature_fprop_[mesh_.face_handle(hfhi)] = ff_id;
      for (auto hfe_it = mesh_.hfe_iter(hfhi); hfe_it.valid(); ++hfe_it)
      {
        feature_face_edge_[*hfe_it] = true;

        feature_face_vertex_[mesh_.edge(*hfe_it).from_vertex()] = true;
        feature_face_vertex_[mesh_.edge(*hfe_it).to_vertex()] = true;
      }

      //restore original face vertex
      hf_vhs[i] = vha;
    }
  }

  return vh_new;
}

//=============================================================================
} // namespace AlgoHex
//=============================================================================


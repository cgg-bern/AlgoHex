/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define VERTEX_OPTIMIZET_C

#include "VertexOptimizeT.hh"
#include <CoMISo/NSolver/NPTiming.hh>
#include <CoMISo/NSolver/IPOPTSolver.hh>
#include <CoMISo/NSolver/NewtonSolver.hh>
#include <Base/Debug/DebConfig.hh>


namespace AlgoHex
{
template<class MeshT>
bool
VertexOptimizeT<MeshT>::vertex_optimize(const VH _vh)
{
  //check if there is any degenerated tet
  if (has_invalid_incident_tet(_vh))
  {
    std::cout << "Error: has flipped or degenerate tet at vertex " << _vh << std::endl;
    return false;
  }

  //subdivide surface triangles and take the minimum point
  if (feature_face_vertex_[_vh])
  {
    Point old_point;
    old_point = mesh_.vertex(_vh);

    std::map<CH, int> alignments;
    collect_alignments(_vh, alignments);

    optimize_feature_face_vertex(_vh);
    if (has_invalid_incident_tet(_vh))
    {
      std::cout << "Error: boundary regular vertex optimization leads to invalid tets!!! Set to original point."
                << std::endl;
      mesh_.set_vertex(_vh, old_point);
      return false;
    }

    optimize_quaternions(_vh, alignments);
  }
  else
  {
    Point old_point, new_point;
    old_point = mesh_.vertex(_vh);
    if (!optimize_interior_vertex_location(_vh, new_point))
      return false;
    mesh_.set_vertex(_vh, new_point);

    if (has_invalid_incident_tet(_vh))
    {
      std::cout << "Warning: tet flips after vertex optimization. Set to original point." << std::endl;
      mesh_.set_vertex(_vh, old_point);
      return false;
    }
  }

  return true;
}


template<class MeshT>
bool
VertexOptimizeT<MeshT>::interior_vertex_optimize(const VH _vh)
{
  //check if there is any degenerated tet
  if (has_invalid_incident_tet(_vh))
  {
    std::cout << "Error: has flipped or degenerate tet at vertex " << _vh << std::endl;
    return false;
  }

  Point old_point, new_point;
  old_point = mesh_.vertex(_vh);
  if (!optimize_interior_vertex_location(_vh, new_point))
    return false;
  mesh_.set_vertex(_vh, new_point);

  if (has_invalid_incident_tet(_vh))
  {
    std::cout << "Warning: tet flips after vertex optimization. Set to original point." << std::endl;
    mesh_.set_vertex(_vh, old_point);
    return false;
  }

  return true;
}

template<class MeshT>
bool
VertexOptimizeT<MeshT>::is_vertex_optimize_ok(const VH _vh) const
{
  if (!_vh.is_valid() || mesh_.is_deleted(_vh))
    return false;

  return true;
}

template<class MeshT>
bool
VertexOptimizeT<MeshT>::has_invalid_incident_tet(const VH _vh) const
{
  for (auto ch: mesh_.vertex_cells(_vh))
  {
    auto cell_points = get_cell_points(ch);
    if (!RemeshingAssist<MeshT>::is_tet_valid(cell_points))
      return true;
  }

  return false;
}

template<class MeshT>
bool
VertexOptimizeT<MeshT>::optimize_interior_vertex_location(const VH _vh, Point &_new_point, const int _max_iter) const
{
  //will be updated
  std::vector<std::vector<Point> > cells_points = get_vertex_cells_points(_vh);
  if (!cells_points.empty())
  {
    if (energy_type_ == WORST_QUALITY)
    {
      _new_point = optimize_interior_vertex_location_worst(cells_points);
    }
    else if (energy_type_ == AVRG_QUALITY)
      optimize_interior_vertex_location_average(cells_points, _new_point, _max_iter);
  }
  return true;
}

template<class MeshT>
bool
VertexOptimizeT<MeshT>::optimize_interior_vertex_location_average(std::vector<std::vector<Point> > &cells_points,
                                                                  Point &_new_point, const int _max_iter) const
{
  int iter(0);
  Eigen::Vector3d J;
  Eigen::Matrix3d H;
  Eigen::Vector3d x;
  x << cells_points[0][0][0], cells_points[0][0][1], cells_points[0][0][2];
  Eigen::Vector3d dx;
  double old_f(0.);

  while (iter < _max_iter)
  {
    if (!update_newton(cells_points, old_f, J, H))
      break;

    dx = H.colPivHouseholderQr().solve(-J);

    double lambda2 = -J.transpose() * dx;

    if (lambda2 < 1e-9)
      break;

    double step_size(1.);
    step_size = backtracking_line_search(x, cells_points, old_f, J, dx, step_size);

    //update
    x += step_size * dx;
    update_points(x, cells_points);

    iter++;
  }

  if (iter == 0)
    return false;

  _new_point = cells_points[0][0];

  return true;
}

template<class MeshT>
typename MeshT::PointT
VertexOptimizeT<MeshT>::optimize_interior_vertex_location_worst(
        const std::vector<std::vector<Point> > &_cells_points) const
{
  //When used in edge collapse, virtual new cells which are nealy degenerate may fail with Newton Solver. It's not a problem.
  Point new_point = _cells_points[0][0];
  // fe problem
  COMISO::FiniteElementProblem fe_problem(3);

  COMISO::FiniteElementSet<EDE_PH_AD> fe_tet_deform;
  setup_tet_deformation_objective_function(fe_tet_deform, _cells_points);

  fe_problem.add_set(&fe_tet_deform);

  for (unsigned int i = 0; i < 3; ++i)
    fe_problem.x()[i] = _cells_points[0][0][i];
//        std::cerr << "initial energy  at: "<<new_point<<" energy "<< fe_problem.eval_f(fe_problem.x().data()) << std::endl;

  // solve
//        COMISO::NPTiming fe_problem_timings(&fe_problem); // with timings
//        fe_problem.print_objectives();

//==========
  // Newton
  // output level of COMISO (silent=0)
  Debug::ScopedOutputLevel sol(0);
  COMISO::NewtonSolver ns(1e-4, 1e-9, 200);
  ns.solve(&fe_problem);
  new_point[0] = fe_problem.x()[0];
  new_point[1] = fe_problem.x()[1];
  new_point[2] = fe_problem.x()[2];
//        std::cerr << "final energy  at: "<<new_point<<" energy "<< fe_problem.eval_f(fe_problem.x().data()) << std::endl;

  return new_point;

//==========
//        try {
//            COMISO::IPOPTSolver ipsol;
//            ipsol.set_ipopt_option("print_level", 0);
//            ipsol.set_max_iterations(100);
//            ipsol.set_ipopt_option("tol", 1e-4);
////    ipsol.set_ipopt_option("hessian_approximation", "limited-memory");
////    ipsol.set_ipopt_option("limited_memory_max_history", 20);
//            ipsol.solve(&fe_problem_timings);
//
//            new_point[0] = fe_problem.x()[0];
//            new_point[1] = fe_problem.x()[1];
//            new_point[2] = fe_problem.x()[2];
////            std::cerr << "final frame field energy  at: "<<_new_point<<" energy "<< fe_problem.eval_f(fe_problem.x().data()) << std::endl;
//
//        }
//        catch (...)
//        {
//            std::cerr << "ERROR: IPOPT threw an exception!!!!!" << std::endl;
//        }
}

template<class TetMeshT>
void
VertexOptimizeT<TetMeshT>::
setup_tet_deformation_objective_function(COMISO::FiniteElementSet<EDE_PH_AD> &_fe_ff,
                                         const std::vector<std::vector<Point>> &_cells_pts) const
{
  // clear old data
  _fe_ff.instances().clear_elements();

  // determine parameters
  double w_conformal = determin_suitable_weight(_cells_pts);
//        std::cerr<<"w: "<<w_conformal<<std::endl;
  for (const auto &c_pts: _cells_pts)
  {
    COMISO::FiniteElementSet<EDE_PH_AD>::VecI idx;
    idx[0] = 0;
    idx[1] = 1;
    idx[2] = 2;

    COMISO::FiniteElementSet<EDE_PH_AD>::VecC c;

    // add one element per tetrahedron
    Eigen::Matrix<double, 3, 4> P;

    // collect points
    for (int i = 0; i < 4; ++i)
      P.col(i) = ovm2eigen(c_pts[i]);

    Eigen::Matrix3d m;
    m.col(0) = P.col(1) - P.col(0);
    m.col(1) = P.col(2) - P.col(0);
    m.col(2) = P.col(3) - P.col(0);

    EDE_PH_AD::compute_constants(w_conformal, P.col(0), P.col(1), P.col(2), P.col(3), c);

    // add instance to set
    _fe_ff.instances().add_element(idx, c);
  }
}

template<class TetMeshT>
double
VertexOptimizeT<TetMeshT>::
determin_suitable_weight(const std::vector<std::vector<Point>> &_cells_pts) const
{
  // determine parameters
  double w_conformal = 0.;
  double max_fx = -1.;

  for (const auto &c_pts: _cells_pts)
  {
    COMISO::FiniteElementSet<EDE_PH_AD>::VecC c;

    // add one element per tetrahedron
    Eigen::Matrix<double, 3, 4> P;

    // collect points
    for (int i = 0; i < 4; ++i)
      P.col(i) = ovm2eigen(c_pts[i]);

    Eigen::Matrix3d m;
    m.col(0) = P.col(1) - P.col(0);
    m.col(1) = P.col(2) - P.col(0);
    m.col(2) = P.col(3) - P.col(0);

    EDE_PH_AD::compute_constants(1., P.col(0), P.col(1), P.col(2), P.col(3), c);

    double fx = EDE_PH_AD::eval_fx(P.col(0), c);
    if (max_fx < fx)
      max_fx = fx;
  }

  w_conformal = log(1e5) / max_fx;

  return w_conformal;
}


template<class MeshT>
bool
VertexOptimizeT<MeshT>::optimize_feature_face_vertex(const VH _vh, const int _n_subdivision)
{
  std::vector<Point> s_points;
  s_points.push_back(mesh_.vertex(_vh));
  for (auto vf_it = mesh_.vf_iter(_vh); vf_it.valid(); ++vf_it)
  {
    if (feature_fprop_[*vf_it] > 0)
    {
      auto hf_vhs = mesh_.get_halfface_vertices(mesh_.halfface_handle(*vf_it, 0));
      RemeshingAssist<MeshT>::sample_points(mesh_.vertex(hf_vhs[0]), mesh_.vertex(hf_vhs[1]), mesh_.vertex(hf_vhs[2]),
                                            0, _n_subdivision, s_points);
    }
  }

  std::vector<std::vector<Point> > cells_points = get_vertex_cells_points(_vh);

  Point new_pt;
  if (energy_type_ == WORST_QUALITY)
    new_pt = optimize_feature_face_vertex_location_worst(s_points, cells_points);
  else if (energy_type_ == AVRG_QUALITY)
    new_pt = optimize_feature_face_vertex_location_average(s_points, cells_points);

  mesh_.set_vertex(_vh, new_pt);

  return true;
}

template<class MeshT>
typename MeshT::PointT
VertexOptimizeT<MeshT>::optimize_feature_face_vertex_location_worst(const std::vector<Point> sample_points,
                                                                    std::vector<std::vector<Point> > &cells_points) const
{
  double min_f = std::numeric_limits<double>::max();
  Point new_pt(0.);
  for (const auto &pt: sample_points)
  {
    update_points(Eigen::Vector3d(pt[0], pt[1], pt[2]), cells_points);
    double curf = get_worst_f(cells_points);
    if (curf < min_f)
    {
      min_f = curf;
      new_pt = pt;
    }
  }

  return new_pt;
}

template<class MeshT>
typename MeshT::PointT
VertexOptimizeT<MeshT>::optimize_feature_face_vertex_location_average(const std::vector<Point> sample_points,
                                                                      std::vector<std::vector<Point> > &cells_points) const
{
  double min_f = std::numeric_limits<double>::max();
  Point new_pt(0.);
  for (const auto &pt: sample_points)
  {
    update_points(Eigen::Vector3d(pt[0], pt[1], pt[2]), cells_points);
    double curf = get_f(cells_points);
    if (curf < min_f)
    {
      min_f = curf;
      new_pt = pt;
    }
  }

  return new_pt;
}


template<class MeshT>
bool
VertexOptimizeT<MeshT>::optimize_vertex_on_edge(const VH _vh, const Point &_pt_s, const Point &_pt_e,
                                                const int _n_subdivision)
{
  std::vector<Point> s_points;
  int n_total = _n_subdivision + 1;
  for (int i = 0; i <= n_total; ++i)
  {
    s_points.push_back((double) i / (double) n_total * _pt_s + (double) (n_total - i) / (double) n_total * _pt_e);
  }

  std::vector<std::vector<Point> > cells_points = get_vertex_cells_points(_vh);

  double min_f = std::numeric_limits<double>::max();
  Point new_pt;
  if (energy_type_ == WORST_QUALITY)
    new_pt = optimize_vertex_on_edge_worst(s_points, cells_points);
  else if (energy_type_ == AVRG_QUALITY)
    new_pt = optimize_vertex_on_edge_average(s_points, cells_points);

  mesh_.set_vertex(_vh, new_pt);

  return true;
}


template<class MeshT>
typename MeshT::PointT
VertexOptimizeT<MeshT>::optimize_vertex_on_edge_worst(const std::vector<Point> &_sample_points,
                                                      std::vector<std::vector<Point> > &_cells_points) const
{
  double min_f = std::numeric_limits<double>::max();
  Point new_pt = _cells_points[0][0];
  for (const auto &pt: _sample_points)
  {
    update_points(Eigen::Vector3d(pt[0], pt[1], pt[2]), _cells_points);
    double curf = get_worst_f(_cells_points);
    if (curf < min_f)
    {
      min_f = curf;
      new_pt = pt;
    }
  }

  return new_pt;
}

template<class MeshT>
typename MeshT::PointT
VertexOptimizeT<MeshT>::optimize_vertex_on_edge_average(const std::vector<Point> &_sample_points,
                                                        std::vector<std::vector<Point> > &_cells_points) const
{
  double min_f = std::numeric_limits<double>::max();
  Point new_pt(0);
  for (const auto &pt: _sample_points)
  {
    update_points(Eigen::Vector3d(pt[0], pt[1], pt[2]), _cells_points);
    double curf = get_f(_cells_points);
    if (curf < min_f)
    {
      min_f = curf;
      new_pt = pt;
    }
  }

  return new_pt;
}


template<class MeshT>
bool
VertexOptimizeT<MeshT>::update_newton(const std::vector<std::vector<Point> > &_cells_points,
                                      double &_f, Eigen::Vector3d &_g, Eigen::Matrix3d &_H) const
{
  for (int i = 0; i < 3; ++i)
  {
    _g(i) = 0.;
    for (int j = 0; j < 3; ++j)
      _H(i, j) = 0.;
  }
  //update energy
  _f = get_f(_cells_points);

  //update gradient
  double gi[3] = {0.};
  double hi[9] = {0.};
  for (const auto &cell_points: _cells_points)
  {
    RemeshingAssist<MeshT>::compute_IMRM_gradient(cell_points, gi);
    RemeshingAssist<MeshT>::compute_IMRM_hessian(cell_points, hi);
    for (int j = 0; j < 3; j++)
    {
      _g(j) += gi[j];
      _H(j, 0) += hi[j * 3 + 0];
      _H(j, 1) += hi[j * 3 + 1];
      _H(j, 2) += hi[j * 3 + 2];
    }
  }

  if (_f >= 1e50)
  {
    _f = 1e50;
    return false;
  }
  if (_f <= 0)
  {
    std::cout << "\nError: energy is less than 0";
    return false;
  }
  if (!_g.allFinite())
  {
    std::cout << "\nError: gradient is infinite";
    return false;
  }
  if (!_H.allFinite())
  {
    std::cout << "\nError: hessian is infinite";
    return false;
  }

  return true;
}

template<class MeshT>
double
VertexOptimizeT<MeshT>::backtracking_line_search(const Eigen::Vector3d &_x,
                                                 const std::vector<std::vector<Point> > &_cells_points,
                                                 const double _old_f, const Eigen::Vector3d &_g,
                                                 const Eigen::Vector3d &_dx, double _t,
                                                 const double _alpha, const double _beta) const
{
  auto cells_points_cpy = _cells_points;
  double gtdx = _g.transpose() * _dx;

  //update points
  update_points(_x + _t * _dx, cells_points_cpy);
  //update energy
  double new_f = get_f(cells_points_cpy);

  while (!(new_f <= _old_f + _alpha * _t * gtdx))
  {
    _t *= _beta;

    update_points(_x + _t * _dx, cells_points_cpy);
    new_f = get_f(cells_points_cpy);
  }

  return _t;
}

//
template<class MeshT>
double
VertexOptimizeT<MeshT>::get_f(const std::vector<std::vector<Point> > &_cells_points) const
{
  double new_f = 0.;
  for (const auto &cell_points: _cells_points)
    new_f += RemeshingAssist<MeshT>::compute_IMRM_energy(cell_points);

  return new_f;
}

//
template<class MeshT>
double
VertexOptimizeT<MeshT>::get_worst_f(const std::vector<std::vector<Point> > &_cells_points) const
{
  double max_f = -1.;
  for (const auto &cell_points: _cells_points)
  {
    double f = RemeshingAssist<MeshT>::compute_IMRM_energy(cell_points);
    if (max_f < f)
      max_f = f;
  }

  return max_f;
}

//
template<class MeshT>
void
VertexOptimizeT<MeshT>::update_points(const Eigen::Vector3d &_x, std::vector<std::vector<Point> > &_cells_points) const
{
  for (size_t i = 0; i < _cells_points.size(); ++i)
  {
    _cells_points[i][0][0] = _x(0);
    _cells_points[i][0][1] = _x(1);
    _cells_points[i][0][2] = _x(2);
  }
}

template<class MeshT>
std::vector<std::vector<typename MeshT::PointT> >
VertexOptimizeT<MeshT>::get_incident_boundary_faces_points(const VH _vh) const
{
  std::vector<std::vector<Point> > faces_points;
  for (const auto &cvh: mesh_.vertex_cells(_vh))
  {
    for (const auto hfh: mesh_.cell(cvh).halffaces())
    {
      if (mesh_.is_boundary(mesh_.face_handle(hfh)))
      {
        std::vector<Point> hf_points;
        hf_points.reserve(3);
        for (auto hfv_it = mesh_.hfv_iter(hfh); hfv_it.valid(); ++hfv_it)
          hf_points.push_back(mesh_.vertex(*hfv_it));

        faces_points.push_back(hf_points);
      }
    }
  }

  return faces_points;
}

template<class MeshT>
std::vector<std::vector<typename MeshT::PointT> >
VertexOptimizeT<MeshT>::get_vertex_cells_points(const VH _vh) const
{
  std::vector<std::vector<Point> > cells_points;
  for (const auto vch: mesh_.vertex_cells(_vh))
  {
    auto cell_points = get_cell_points(vch, _vh);
    cells_points.push_back(cell_points);
  }

  return cells_points;
}

template<class MeshT>
std::vector<typename MeshT::PointT>
VertexOptimizeT<MeshT>::get_cell_points(const CH _ch) const
{
  auto cell_vertices = mesh_.get_cell_vertices(_ch);
  std::vector<Point> cell_points;
  cell_points.reserve(4);
  for (const auto cv: cell_vertices)
    cell_points.push_back(mesh_.vertex(cv));

  return cell_points;
}

template<class MeshT>
std::vector<typename MeshT::PointT>
VertexOptimizeT<MeshT>::get_cell_points(const CH _ch, const VH _vh) const
{
  auto cell_vertices = mesh_.get_cell_vertices(_ch, _vh);
  std::vector<Point> cell_points;
  cell_points.reserve(4);
  for (const auto cv: cell_vertices)
    cell_points.push_back(mesh_.vertex(cv));

  return cell_points;
}


template<class MeshT>
void
VertexOptimizeT<MeshT>::optimize_quaternions(const VH _vh, std::map<CH, int> &_aligned_axis)
{
  //smooth quaternions
  std::vector<CH> cells;
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    cells.push_back(*vc_it);

  std::set<FH> bdy_fhs;
  for (const auto ch: cells)
    for (auto cf_it = mesh_.cf_iter(ch); cf_it.valid(); ++cf_it)
      if (mesh_.is_boundary(*cf_it))
        bdy_fhs.insert(*cf_it);

  std::map<CH, std::vector<std::pair<Vec3d, int>>> alignments;

  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0)
      {
        auto nm = mesh_.normal(*chf_it);

        alignments[*vc_it].emplace_back(Vec3d(nm[0], nm[1], nm[2]), _aligned_axis[*vc_it]);
      }
    }
  }

  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto ce_it = mesh_.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
    {
      if (feature_edge_[*ce_it] > 0)
      {
        auto dir = mesh_.vertex(mesh_.edge(*ce_it).to_vertex()) - mesh_.vertex(mesh_.edge(*ce_it).from_vertex());
        dir.normalize();

        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], dir).second;
        alignments[*vc_it].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
      }
    }
  }

  for (int i = 0; i < 10; ++i)
    QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, cells, bdy_fhs, alignments, tq_,
                                                 cell_quaternions_);
}

template<class MeshT>
void
VertexOptimizeT<MeshT>::collect_alignments(const VH _vh, std::map<CH, int> &_aligned_axis)
{
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
  {
    for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
    {
      auto fhi = mesh_.face_handle(*chf_it);
      if (feature_fprop_[fhi] > 0)
      {
        auto nm = mesh_.normal(*chf_it);
        int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], nm).second;

        _aligned_axis[*vc_it] = axis;
      }
    }
  }
}
}

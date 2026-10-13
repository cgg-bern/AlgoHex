/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <CoMISo/NSolver/NProblemInterface.hh>
#include <AlgoHex/TypeDef.hh>

namespace AlgoHex
{
class SingularVertexOptProblemQuaternion : public COMISO::NProblemInterface
{
public:
  using Edge = std::pair<int, int>;
  using Face = std::array<int, 3>;
  using VTuple = std::array<int, 3>;
  using Cell = std::array<int, 4>;
  using DualEdge = std::pair<int, int>;

  using Mat = Mat3d;

  SingularVertexOptProblemQuaternion(const std::vector<double> &_points, std::vector<double> &_partial_points);

  ~SingularVertexOptProblemQuaternion() {}

  //
  virtual int n_unknowns() { return partial_points_.size(); }

  virtual void initial_x(double *_x);

  double initial_f();

  virtual double eval_f(const double *_x);

  double initial_f(const double *_x);

  //project px to the plane defined by p0 and n
  template<typename Vec3>
  static void project_on_tangent_plane(const Vec3 &_p0, const Vec3 &_px, const Vec3 &_n,
                                       Vec3 &_pp);

  virtual void eval_gradient(const double *_x, double *_g);

  virtual void eval_hessian(const double *_x, SMatrixNP &_H);

  virtual void store_result(const double *_x);

  //one unknown vertex
  void add_imrm_energy_element_one(const int _v_newid0, const int _vid1, const int _vid2, const int _vid3,
                                   const double _weight);
  //two
  void add_imrm_energy_element_two(const int _v_newid0, const int _v_newid1, const int _vid2, const int _vid3,
                                   const double _weight);
  //three
  void
  add_imrm_energy_element_three(const int _v_newid0, const int _v_newid1, const int _v_newid2, const int _vid3,
                                const double _weight);

  //four
  void
  add_imrm_energy_element_four(const int _v_newid0, const int _v_newid1, const int _v_newid2, const int _v_newid3,
                               const double _weight);

  void
  add_boundary_simrm_element(const int _v_newid0, const int _vid1, const int _vid2, const Vec3d &_n,
                             const double _weight);

  void add_anti_normal_flip_element(const int _v_newid0, const int _vid1, const int _vid2, const Vec3d &_top_pt,
                                    const double _weight);

  void add_perpendicular_element(const int _vid0, const int _vid1, const Vec3d &_n, const double _w);

  void add_shrink_edge_element(const int _vid0, const int _vid1, const double _w);

  void add_curvature_smooth_element(const int _vid0, const int _vid1, const int _vid2, const double _w);

  void add_repulsion_term(const int _vid, const Vec3d &_target_point, const double _w);

private:
  //complete mesh vertex position
  const std::vector<double> &points_;
  //partial mesh vertex position
  std::vector<double> &partial_points_;

  //inverse mean ratio metric
  std::vector<Cell> cells1_, cells2_, cells3_, cells4_;
  std::vector<double> imrm_weights1_, imrm_weights2_, imrm_weights3_, imrm_weights4_;

  //
  Eigen::VectorXd x_imrme_;

  Vec3d g_imrme1_;
  Mat3d h_imrme1_;

  Eigen::VectorXd g_imrme2_;
  Eigen::MatrixXd h_imrme2_;

  Eigen::VectorXd g_imrme3_;
  Eigen::MatrixXd h_imrme3_;

  Eigen::VectorXd g_imrme4_;
  Eigen::MatrixXd h_imrme4_;


  //surface inverse mean ratio
  std::vector<Face> faces_;
  std::vector<Vec3d> normals_;
  std::vector<double> simrm_weights_;
  //
  Eigen::VectorXd x_simrme_;
  Vec3d g_simrme_;
  Mat3d h_simrme_;


  //anti normal flip energy
  std::vector<Face> base_faces_;
  std::vector<Vec3d> top_pts_;
  std::vector<double> c_weights_;
  //
  Eigen::VectorXd x_nme_;
  Vec3d g_nme_;
  Mat3d h_nme_;


  //axis perpendicular to edge direction (dot or cross product)
  std::vector<Edge> ppd_edges_;
  std::vector<Vec3d> ppd_axes_;
  std::vector<double> ppd_weights_;

  Eigen::VectorXd x_pe_;
  Eigen::VectorXd g_pe_;
  Eigen::MatrixXd h_pe_;

  //edge length term
  std::vector<Edge> comp_edges_;
  std::vector<double> cp_weights_;

  Eigen::VectorXd x_cpe_;
  Eigen::VectorXd g_cpe_;
  Eigen::MatrixXd h_cpe_;


  //curvature smoothing
  std::vector<VTuple> vertex_tuples_;
  std::vector<double> cv_weights_;

  Eigen::VectorXd x_cve_;
  Eigen::VectorXd g_cve_;
  Eigen::MatrixXd h_cve_;


  //repulsion term
  std::vector<int> rp_vertices_;
  std::vector<Vec3d> target_points_;
  std::vector<double> rp_weights_;

  Eigen::VectorXd x_rpe_;
  Vec3d g_rpe_;
  Mat3d h_rpe_;
};

extern template void ALGOHEX_EXPORT SingularVertexOptProblemQuaternion::project_on_tangent_plane<OVM::Geometry::Vec3d>(const OVM::Geometry::Vec3d &, const OVM::Geometry::Vec3d &, const OVM::Geometry::Vec3d &, OVM::Geometry::Vec3d &);
extern template void ALGOHEX_EXPORT SingularVertexOptProblemQuaternion::project_on_tangent_plane<Eigen::Vector3d>(const Eigen::Vector3d &, const Eigen::Vector3d &, const Eigen::Vector3d &, Eigen::Vector3d &);

}

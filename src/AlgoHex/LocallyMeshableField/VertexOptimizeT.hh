/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <AlgoHex/TransitionQuaternionEigen.hh>
#include <AlgoHex/SingularGraphExtractionT.hh>
#include <map>
#include "QuaternionsSmoothing.hh"
#include "MeshProperties.hh"
#include "RemeshingAssist.hh"
#include "CommonFuncs.hh"
#include <AlgoHex/TypeDef.hh>
#include <CoMISo/NSolver/FiniteElementProblem.hh>
#include "ExpDeformationElements3D.hh"
#include <AlgoHex/Geometry.hh>


namespace AlgoHex
{
template<class MeshT>
class VertexOptimizeT
{
public:
  using Point = typename MeshT::PointT;

  VertexOptimizeT(MeshT &_mesh) :
          mesh_(_mesh),
          trans_prop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
          cell_quaternions_(mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions")),
          feature_edge_(mesh_.template request_edge_property<int>("AlgoHex::FeatureEdges")),
          feature_fprop_(mesh_.template request_face_property<int>("AlgoHex::FeatureFaces")),
          feature_face_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureFaceVertices")) {}

  ~VertexOptimizeT() {}

  enum ENERGY_TYPE
  {
    WORST_QUALITY, AVRG_QUALITY
  };

public:
  bool vertex_optimize(const VH _vh);

  bool is_vertex_optimize_ok(const VH _vh) const;

  bool has_invalid_incident_tet(const VH _vh) const;

  bool optimize_interior_vertex_location(const VH _vh, Point &_new_point, const int _max_iter = 15) const;

  bool optimize_interior_vertex_location_average(std::vector<std::vector<Point> > &cells_points, Point &_new_point,
                                                 const int _max_iter) const;

  Point optimize_interior_vertex_location_worst(const std::vector<std::vector<Point> > &_cells_points) const;

  bool optimize_feature_face_vertex(const VH _vh, const int _n_subdivision = 2);

  Point optimize_feature_face_vertex_location_worst(const std::vector<Point> sample_points,
                                                    std::vector<std::vector<Point> > &cells_points) const;

  Point optimize_feature_face_vertex_location_average(const std::vector<Point> sample_points,
                                                      std::vector<std::vector<Point> > &cells_points) const;

  bool interior_vertex_optimize(const VH _vh);

  bool
  optimize_vertex_on_edge(const VH _vh, const Point &_pt_s, const Point &_pt_e, const int _n_subdivision = 10);

  Point optimize_vertex_on_edge_worst(const std::vector<Point> &_sample_points,
                                      std::vector<std::vector<Point> > &_cells_points) const;

  Point optimize_vertex_on_edge_average(const std::vector<Point> &_sample_points,
                                        std::vector<std::vector<Point> > &_cells_points) const;

  void set_energy_type(ENERGY_TYPE energyType) { energy_type_ = energyType; }

  std::vector<std::vector<Point> > get_vertex_cells_points(const VH _vh) const;

  std::vector<Point> get_cell_points(const CH _ch) const;

  std::vector<Point> get_cell_points(const CH _ch, const VH _vh) const;

private:

  bool update_newton(const std::vector<std::vector<Point> > &_cells_points,
                     double &_f, Eigen::Vector3d &_g, Eigen::Matrix3d &_H) const;

  double
  backtracking_line_search(const Eigen::Vector3d &_x, const std::vector<std::vector<Point> > &_cells_points,
                           const double _old_f, const Eigen::Vector3d &_g, const Eigen::Vector3d &_dx, double _t,
                           const double _alpha = 0.2, const double _beta = 0.6) const;
  
  double get_f(const std::vector<std::vector<Point> > &_cells_points) const;

  double get_worst_f(const std::vector<std::vector<Point> > &_cells_points) const;

  void update_points(const Eigen::Vector3d &_x, std::vector<std::vector<Point> > &_cells_points) const;

  std::vector<std::vector<Point> > get_incident_boundary_faces_points(const VH _vh) const;

  void optimize_quaternions(const VH _vh, std::map<CH, int> &_aligned_axis);

  void collect_alignments(const VH _vh, std::map<CH, int> &_aligned_axis);

  void setup_tet_deformation_objective_function(COMISO::FiniteElementSet<EDE_PH_AD> &_fe_ff,
                                                const std::vector<std::vector<Point>> &_cells_pts) const;

  double determin_suitable_weight(const std::vector<std::vector<Point>> &_cells_pts) const;


private:
  MeshT &mesh_;
  HFP<int> trans_prop_;
  CP<Quaternion> cell_quaternions_;
  EP<int> feature_edge_;
  FP<int> feature_fprop_;
  VP<bool> feature_face_vertex_;

  bool energy_type_ = AVRG_QUALITY;

  TransitionQuaternion tq_;
};
}
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(VERTEX_OPTIMIZET_C)
#define TETREMESHINGT_TEMPLATES

#include "VertexOptimizeT_impl.hh"

#endif



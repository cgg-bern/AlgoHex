/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <functional>
#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>
#include <CoMISo/NSolver/ConstraintTools.hh>
#include "SingularVertexOptProblem.hh"
#include <CoMISo/NSolver/NPTiming.hh>
#include <AlgoHex/SingularGraphExtractionT.hh>
#include "AlgoHex/BinarySpacePartitionTrees/TriangleBSPT.hh"

#include "TetRemeshingT.hh"
#include "../Stopwatches.hh"
#include "QuaternionsSmoothing.hh"

namespace OVM = OpenVolumeMesh;

namespace AlgoHex
{
template<class MeshT>
class MeshOptimizationT : public MeshPropertiesT<MeshT>
{
public:
  using Point = typename MeshT::PointT;

  using PairDHE = std::pair<double, HEH>;
  using PairDE = std::pair<double, EH>;
  using PairDI = std::pair<double, int>;

  //properties
  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;

  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::feature_node_;
  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_edge_vertex_;
  using MeshPropertiesT<MeshT>::feature_face_vertex_;
  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::feature_face_edge_;
  using MeshPropertiesT<MeshT>::target_length_;


  using SMatD = Eigen::SparseMatrix<double>;
  using VecD = Eigen::VectorXd;
  using T = Eigen::Triplet<double>;

public:

  MeshOptimizationT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                                    mesh_(_mesh),
                                    face_normal_(
                                            sf_mesh_.template request_face_property<Point>("original face normal")),
                                    vertex_normal_(
                                            sf_mesh_.template request_vertex_property<Point>("original vertex normal")),
                                    sharp_face_(sf_mesh_.template request_face_property<bool>("original sharp face")),
                                    sm_feature_edge_(
                                            sf_mesh_.template request_edge_property<int>("surface mesh: feature edge")),
                                    imrm_cell_(
                                            mesh_.template request_cell_property<bool>("cell with different weight")),
                                    sg_edge_pairs_(mesh_.template request_edge_property<int>("separate edge pairs")),
                                    keep_on_feature_face_(
                                            mesh_.template request_edge_property<bool>("sge kept on feature face",
                                                                                       false)),
                                    sge_(mesh_)
  {
    //reference mesh
    sf_mesh_.set_persistent(face_normal_, true);
    sf_mesh_.set_persistent(vertex_normal_, true);
    sf_mesh_.set_persistent(sm_feature_edge_, true);

    this->initialize_feature_vertex_property();
    build_surface_mesh(mesh_, sf_mesh_);
    bsp_ = build_bsp(sf_mesh_);

    //get mesh volume
    for (const auto &ch: mesh_.cells())
      mesh_volume_ += get_cell_volume(mesh_, ch);

    scale_ = 0.1 * std::pow(mesh_volume_, 1. / 3.);

    //cache feature edges in sf mesh
    for (auto ehi: sf_mesh_.edges())
    {
      if (sm_feature_edge_[ehi] > 0)
        sm_fes_[sm_feature_edge_[ehi]].push_back(ehi);
    }
  }

  ~MeshOptimizationT() {}

public:
  void remesh(const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w, const double _rp_w,
              const int k);

  void remesh(const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w, const double _rp_w,
              const int k, std::vector<MeshT> &meshes);

  void remesh(int _pair_id, const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w,
              const double _rp_w, const int k);

  void remesh_post(const int _remesh_iter);

  void remesh_post(const int _remesh_iter, std::vector<MeshT> &_meshes);

  double get_scale() const { return scale_; }

  std::vector<EH> get_singular_edge_pairs(const int _pair_id) const;

  void check_valence_change();

  void singular_vertices_relocation(const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w,
                                    const double _rp_w, const int _max_iter = 100);

  void singular_vertices_relocation(const std::vector<EH> &_sg_ehs, const double _sge_w, const double _rgl_w,
                                    const double _cs_w, const double _tps_w, const double _rp_w,
                                    const int _max_iter = 100);

  int split_edges(const int k, const bool _split_with_opt = false);

  int split_edges(const std::vector<EH> &_sg_ehs, const int k);

  int collapse_edges(const int k, const bool _post_remesh = false, const bool _check_energy = true);

  int collapse_edges(const std::vector<EH> &_sg_ehs, const int k, const bool _check_energy = true);

  int swap_edges(const int k, const bool _post_remesh = false);

  int swap_edges(const std::vector<EH> &_sg_ehs, const int k);

  //postprocess
//        void optimize_singular_vertices();
//        void optimize_key_vertices();

  void optimize_regular_vertices();

  void save_singular_graph_development(MeshT &_sg_mesh, int _iter);

  void save_singular_graph(const std::string &_filename, int _id, const std::set<VH> &_fxb_tps);


private:
  void add_alignment_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                            double _alignment_weight) const;

  void add_alignment_energy(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                            std::map<VH, int> &_vh_to_idx, double _alignment_weight) const;

  void add_edge_shrink_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                              double _cs_weight, double _tps_w) const;

  void add_edge_shrink_energy(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                              std::map<VH, int> &_vh_to_idx, double _cs_weight, double _tps_w) const;

  void add_curvature_smooth_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                                   double _rgl_weight) const;

  void add_curvature_smooth_energy(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                                   std::map<VH, int> &_vh_to_idx, double _rgl_weight) const;

  void add_imrm_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                       double _imrm_weight) const;

  void add_imrm_energy2(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                        double _imrm_weight) const;

  void add_surface_imrm_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                               double _sf_imrm_weight) const;

  void add_surface_imrm_energy(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                               std::map<VH, int> &_vh_to_idx, double _sf_imrm_weight) const;

  void add_anti_normal_flip_term(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                                 std::set<VH> &_bad_vhs, double _normal_weight) const;

  void add_anti_normal_flip_term(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                                 std::map<VH, int> &_vh_to_idx, std::set<VH> &_bad_vhs, double _normal_weight) const;

  //add repulsion energy
  void add_repulsion_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                            const double _rp_weight);

  void add_regularizer(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                       const double _rp_weight);

  void collect_constrained_vertices(std::map<VH, int> &_new_vh_idx, const std::set<VH> &_bad_vhs,
                                    std::vector<VH> &_fixed_vhs, std::vector<VH> &_movable_vhs1,
                                    std::vector<VH> &_movable_vhs2);

  void collect_fixed_vertices(std::map<VH, int> &_new_vh_idx, const std::set<VH> &_bad_vhs, std::vector<VH> &_fixed_vhs);

  void collect_fixed_vertices(std::map<VH, int> &_new_vh_idx, std::vector<VH> &_fixed_vhs);

  void collect_fixed_vertices_post(std::map<VH, int> &_new_vh_idx, std::vector<VH> &_fixed_vhs);

  void setup_singular_vertex_constraints(int _n_unknowns, std::map<VH, int> &_new_vh_idx,
                                         std::vector<COMISO::LinearConstraint> &_constraints,
                                         std::vector<COMISO::NConstraintInterface *> &_constraint_pointers);

  void setup_singular_vertex_constraints(int _n_unknowns, const std::vector<VH> &_movable_vhs2,
                                         const std::vector<VH> &_movable_vhs1, const std::vector<VH> &_fixed_vhs,
                                         std::map<VH, int> &_new_vh_idx, SMatD &_A, VecD &_b);

  Vec3d best_aligned_rotational_axis(const HEH _heh) const;

  Point get_normal_on_reference_mesh(const Point &_pt) const;

  void project_vertex_to_original_surface(const VH _vh, const Point &_orig_pt);

  void project_vertex_to_original_feature_arc(const VH _vh, const Point &_orig_pt);

  bool has_invalid_incident_tet(std::vector<std::vector<Point>> &_cells_pts, const Point &_query_pt) const;

  bool has_invalid_incident_triangle(const VH _vh, const Point &_query_pt) const;

  std::vector<Point> get_cell_points(const CH _ch, const VH _vh) const;

  std::vector<VH> get_cell_vertices(const CH _ch) const;

  std::vector<VH> get_cell_vertices2(std::map<VH, int> &_vh_to_idx, const CH _ch) const;

  std::vector<VH> get_cell_vertices_post(const CH _ch) const;

  int n_singular_vertex_in_cell(const CH _ch) const;

  int n_singular_vertex_in_cell2(std::map<VH, int> &_vh_to_idx, const CH _ch) const;

  int n_incident_singular_edges(const VH _vh) const;

  bool is_movable(const VH _vh) const;

  std::map<VH, std::map<HFH, int>> collect_normal_aligned_axis(const bool _post = false);

  std::map<VH, std::map<HFH, int>> collect_normal_aligned_axis(const std::vector<EH> &_sg_ehs);

  std::pair<double, Point> find_closest_point_on_other_singular_arcs(const VH _vh, const std::set<EH> &_same_arc_ehs,
                                                                     const std::set<EH> &_all_ehs) const;

  void get_edge_color(const EH _eh, OVM::Vec4f &_color) const;

  void get_vertex_color(const VH _vh, OVM::Vec4f &_color) const;


public:
  void set_min_target_length(const double _factor = 0.05);

  void adapt_target_length(const double _damping = 0.5, const int _norm_order = 4);

  void measure_field_smoothness(VP<double> &_smth_vprop, const int _norm_order);

  void use_truncated_newton() { truncated_newton_ = true; }

private:
  double measure_field_smoothness_at_vertex(const VH _vh, const FP<double> &_smth_fprop, int _norm_order) const;

private:
  MeshT &mesh_;
  MeshT sf_mesh_;

  //reference mesh property
  FP<Point> face_normal_;
  VP<Point> vertex_normal_;
  FP<bool> sharp_face_;
  EP<int> sm_feature_edge_;

  //
  CP<bool> imrm_cell_;
  EP<int> sg_edge_pairs_;
  EP<bool> keep_on_feature_face_;

  AlgoHex::TransitionQuaternion tq_;
  AlgoHex::SingularGraphExtractionT<MeshT> sge_;

  std::unique_ptr<AlgoHex::OpenVolumeMeshTriangleBSPT<MeshT>> bsp_;

  double mesh_volume_ = 0.;
  double scale_ = 0.;

  //used in post-remeshing when adapt the target the edge length
  double min_target_len_ = 0.;

  int iter_ = 0;

  bool truncated_newton_ = false;

  std::map<int, std::vector<EH>> sm_fes_;
};

}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(MESHOPTIMIZATIONT_C)
#define MESHOPTIMIZATIONT_TEMPLATES

#include "MeshOptimizationT_impl.hh"

#endif
//=============================================================================


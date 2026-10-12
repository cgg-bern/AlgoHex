/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once
#include <AlgoHex/TypeDef.hh>


#include <functional>
#include <AlgoHex/MeshGeometry.hh>
#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>
#include <CoMISo/NSolver/ConstraintTools.hh>
#include "SingularVertexOptProblem.hh"
#include <CoMISo/NSolver/NPTiming.hh>
#include <AlgoHex/SingularGraphExtraction.hh>
#include "AlgoHex/BinarySpacePartitionTrees/TriangleBSPT.hh"

#include "TetRemeshing.hh"
#include "../Stopwatches.hh"
#include "QuaternionsSmoothing.hh"

namespace OVM = OpenVolumeMesh;

namespace AlgoHex
{
class ALGOHEX_EXPORT MeshOptimization : public MeshProperties
{
public:
  using PairDHE = std::pair<double, HEH>;
  using PairDE = std::pair<double, EH>;
  using PairDI = std::pair<double, int>;

  //properties
  using MeshProperties::valence_;
  using MeshProperties::trans_prop_;
  using MeshProperties::cell_quaternions_;

  using MeshProperties::sgl_vt_;
  using MeshProperties::feature_node_;
  using MeshProperties::feature_edge_;
  using MeshProperties::feature_edge_vertex_;
  using MeshProperties::feature_face_vertex_;
  using MeshProperties::feature_fprop_;
  using MeshProperties::feature_face_edge_;
  using MeshProperties::target_length_;


  using SMatD = Eigen::SparseMatrix<double>;
  using VecD = Eigen::VectorXd;
  using T = Eigen::Triplet<double>;

  MeshOptimization(TetMesh &_mesh);

  ~MeshOptimization() {}

  void remesh(const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w, const double _rp_w,
              const int k);

  void remesh(const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w, const double _rp_w,
              const int k, std::vector<TetMesh> &meshes);

  void remesh(int _pair_id, const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w,
              const double _rp_w, const int k);

  void remesh_post(const int _remesh_iter);

  void remesh_post(const int _remesh_iter, std::vector<TetMesh> &_meshes);

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

  void optimize_regular_vertices();

  void save_singular_graph_development(TetMesh &_sg_mesh, int _iter);

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

  TetMesh &mesh_;
  TetMesh sf_mesh_;

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
  AlgoHex::SingularGraphExtraction sge_;

  std::unique_ptr<AlgoHex::OpenVolumeMeshTriangleBSPT<TetMesh>> bsp_;

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
//=============================================================================

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

//== INCLUDES =================================================================

#include <map>
#include <vector>

#include <AlgoHex/Util/json.hh>

#include <CoMISo/NSolver/LinearConstraint.hh>
#include <CoMISo/NSolver/FiniteElementProblem.hh>
#include <CoMISo/NSolver/VariableType.hh>

#include "TypeDef.hh"
#include "TransitionQuaternionEigen.hh"
#include "FrameFittingElement3D.hh"
#include "DeformationElements3D.hh"
#include "DihedralAngleElementTinyAD.hh"
#include "BarrierElements.hh"
#include "LinearLeastSquaresElements.hh"

#include <MC3D/Types.hpp>
#include <QGP3D/PathConstraint.hpp>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================




/** \class Parametrization3D Parametrization3D.hh

    Brief Description.

    A more elaborate description follows.
*/

class ALGOHEX_EXPORT Parametrization3D
{
public:

  enum OptimizerType
  {
    OPT_COMISO, OPT_GUROBI
  };

  // typedefs for shorter access

  typedef OpenVolumeMesh::VertexHandle VH;
  typedef OpenVolumeMesh::EdgeHandle EH;
  typedef OpenVolumeMesh::HalfEdgeHandle HEH;
  typedef OpenVolumeMesh::CellHandle CH;
  typedef OpenVolumeMesh::FaceHandle FH;
  typedef OpenVolumeMesh::HalfFaceHandle HFH;
  typedef OpenVolumeMesh::EdgeIter EIt;
  typedef OpenVolumeMesh::HalfEdgeIter HEIt;
  typedef OpenVolumeMesh::VertexIter VIt;
  typedef OpenVolumeMesh::CellIter CIt;
  typedef OpenVolumeMesh::FaceIter FIt;
  typedef OpenVolumeMesh::BoundaryFaceIter BFIt;
  typedef OpenVolumeMesh::HalfFaceIter HFIt;
  typedef OpenVolumeMesh::HalfEdgeCellIter HECIt;
  typedef OpenVolumeMesh::HalfFaceVertexIter HFVIt;
  typedef OpenVolumeMesh::HalfEdgeHalfFaceIter HEHFIt;
  typedef OpenVolumeMesh::VertexOHalfEdgeIter VOHEIt;
  typedef OpenVolumeMesh::VertexCellIter VCIt;
  typedef OpenVolumeMesh::CellCellIter CCIt;
  typedef OpenVolumeMesh::TetVertexIter TVIt;

  typedef TetMesh::PointT Point;
  typedef TetMesh::Cell Cell;
  typedef TetMesh::Face Face;

  // Eigen Types
  typedef Eigen::Triplet<double> Trip;
  typedef Eigen::SparseMatrix<double> SparseMatrix;
  typedef Eigen::Quaternion<double> Quaternion;
  typedef Eigen::VectorXd VectorXd;
  typedef Eigen::Matrix<double, 4, 1> Vec4;
  typedef Eigen::Matrix<double, 3, 1> Vec3;
  typedef Eigen::Matrix<double, 2, 1> Vec2;
  typedef Eigen::Matrix<double, 9, 9> Mat9x9;
  typedef Eigen::Matrix<double, 3, 3> Mat3x3;
  typedef Eigen::Matrix<double, 2, 2> Mat2x2;
  typedef Eigen::Matrix<double, 3, 2> Mat3x2;
  typedef Eigen::Matrix<int, 3, 1> Vec3i;
  typedef Eigen::Matrix<int, 4, 1> Vec4i;

  typedef std::pair<unsigned int, COMISO::VariableType> PairUiV;

  using QPathConstraint = qgp3d::PathConstraint;

  // constructors
  Parametrization3D(TetMesh &_m, OpenVolumeMesh::StatusAttrib &_mesh_status);

//  parametrize frame field to integer-grid map
//  input:
//    frame_cprop_( mesh_.request_cell_property<<Mat3x3>("AlgoHex::FrameField"))
//    transition_hfprop_( mesh_.request_halfface_property<int>("HalffaceTransiton")),
//    valence_eprop_( mesh_.request_edge_property<int>("edge_valance")),
//    feature_eprop_(mesh_.request_edge_property<bool>("feature edges")),
//    input requirements:
//        (1) field is aligned to boundary normals, singular edges and feature curves
//        (2) each tet has at most one such alignment constraint (tet mesh is split properly) except for boundary normal + boundary feature
//  output:
//    igm_cprop_(mesh_.request_cell_property< std::map<OpenVolumeMesh::VertexHandle, ACG::Vec3d> >("Parametrization"))

  std::pair<int, int> parametrize_complete(const int _num_hex_cells = 1000, const double _anisotropy_alpha = 1.0,
                                           const bool _with_integer_constraints = true);

  // robust quantization based on motorcycle complex
  // precondition: valid seamless map is available
  std::pair<int, int> parametrize_robust_quantization(const int _num_hex_cells = 1000);

  // as before but without the nonlinear optimization part
  void parametrize(const int _num_hex_cells = 1000, const double _anisotropy_alpha = 1.0,
                   const bool _with_integer_constraints = true);

  double mesh_volume();

  double frame_field_volume();

  void set_optimizer(const OptimizerType _ot) { opt_type_ = _ot; }

  void set_timelimit(const double _t) { timelimit_ = _t; }

  void set_stiffening_weight(const double _w) { stiffening_weight_ = _w; }

  void set_max_stiffening_iters(const unsigned int _n) { max_stiffening_iters_ = _n; }

  nlohmann::json &json_data() { return json_data_; }

  void import_frames_from_quaternion_field(const std::string &_prop_name = "FrameFieldQuaternions");

  void import_frames_from_parametrization();

  void add_boundary_faces_to_feature_surfaces();

  void reset_feature_face_property();

  void scale_frames_to_target_complexity(const int _num_hex_cells);

  void optimize_integer_grid_map_robust_quantization(const int _num_hex_cells);

  void optimize_integer_grid_map_locally_injective(const int _num_hex_cells, const double _anisotropy_alpha,
                                                   const double _epsilon, const bool _with_integer_constraints);

  void quantize_qgp3d(const int _num_hex_cells);

  void export_quantization_constraint_paths(TetMesh &_export_mesh);

  void initialize_helper_variables_via_seamless_map(std::vector<double> &_x, const int _helper_var_base_idx);

  void convertTetrahedralMeshIndexing(const TetMesh &_mesh, TetMesh &_tetMesh);

private:
  void check_input_mesh() const;

  // generate cut surface by initializing
  //   OpenVolumeMesh::FacePropertyT<bool> is_on_cut_fprop_;
  //   std::vector< std::vector<FH> > cut_surface_;
  void generate_cut_surface();

  void comb_octahedral_field();

  double octahedral_field_energy() const;

  void optimize_integer_grid_map(const int _num_hex_cells, const double _anisotropy_alpha,
                                 const bool _with_integer_constraints);

  void store_igm(const double *_x);

  void setup_vertex_indices();

  void check_vertex_indices() const;

  void setup_frame_objective_function(COMISO::FiniteElementSet<FrameFittingElement3D> &_fe_ff, const double _alpha,
                                      const double _sizing_scale_factor) const;

  void setup_deformation_objective_function(COMISO::FiniteElementSet<FF3DV_PH_AD> &_fe_ff, const double _epsilon,
                                            const double _target_det_scaling);

  void setup_symmetric_dirichlet_objective_function(COMISO::FiniteElementSet<SDE3D_PH_AD> &_fe_sdirichlet, const double _w,
                                               const double _sizing_scale, const double _inversion_penalty);

  void setup_jacobian_smoothness_objective_function( COMISO::FiniteElementSet<JSE3D_AD>& _fe_jacobian_smoothness, const double _w, const double _sizing_scale);

  void setup_domain_deformation_objective_function(COMISO::FiniteElementSet<FF3DV_PH_AD> &_fe_ff,
                                                   const int _domain_idx_base);

  void setup_ball_barrier_objective_function(COMISO::FiniteElementSet<BB3D_AD_PH> &_fe_bb, const int _domain_idx_base,
                                             const double _radius_scale = 1.0);

  void setup_dual_deformation_objective_function(COMISO::FiniteElementSet<FFD3DV_PH_AD> &_fe_ffd, const double _epsilon,
                                                 const double _target_det_scaling, const int _domain_idx_base);

  void setup_dihedral_angle_objective_function(COMISO::FiniteElementSet<DihedralAngleElement_TAD_PH> &_fe_dihedral,
                                               const double _w_dihedral, const double _s_dihedral);

  void setup_frame_exp_objective_function(COMISO::FiniteElementSet<FrameFittingExpElement3D_TAD> &_fe_frame_fitting,
                                          const double _sizing_factor, const double _w_exp, const double _mesh_volume,
                                          const double _limit_size, std::vector<double> &_optimal_sizing);

  void setup_size_fitting_objective_function(COMISO::FiniteElementSet<LinearLeastSquaresElement1D> &_fe_size_fitting,
                                             const double _weight, const double _mesh_volume);

  void setup_size_barrier_objective_function(COMISO::FiniteElementSet<ReciprocalBarrierElement_TAD> &_fe_size_barrier,
                                             const double _weight, const double _limit_size, const double _mesh_volume);

  void test_objective_function() const;

  void setup_constraints(std::vector<COMISO::LinearConstraint> &_constraints,
                         std::vector<COMISO::NConstraintInterface *> &_constraint_pointers,
                         std::vector<PairUiV> &_discrete_constraints,
                         const int _helper_idx_offset = 0,
                         const bool _allocate_helper_variables = true);

  void add_feature_constraints(std::vector<COMISO::LinearConstraint> &_constraints,
                               std::vector<COMISO::NConstraintInterface *> &_constraint_pointers);

  void reset_cell_weights();

  int update_cell_weights(const std::vector<double> &_x);

  int tet_vidx(const CH _ch, const VH _vh) const;

  int halfface_vidx(const HFH _hfh, const VH _vh) const;

  HFH opposite_halfface_handle(const CH _ch, const VH _vh) const;

  VH find_singular_or_feature_node() const;

  VH find_vertex_on_singular_or_feature_curve() const;

  VH find_vertex_on_feature_surface() const;

  int incident_singular_edges(const VH _vh) const;

  int incident_feature_edges(const VH _vh) const;

  int n_cut_faces(const EH _eh) const;

  int n_cut_faces(const EH _eh, std::vector<FH> &_fhs) const;

  void parametric_tet(const CH _ch, Mat3x4d &_P);

  void parametric_tet(const CH _ch, const std::vector<double> &_x, Mat3x4d &_P);

  double parametric_volume();

  // get jacobi matrix of p.w.l. map from _P to _U
  void jacobi_matrix(const Mat3x4d &_P, const Mat3x4d &_U, Mat3d &_J);

  // copy cell property
  void copy_cell_property(OVM::CellPropertyT<std::map<VH, Point>> &_source, OVM::CellPropertyT<std::map<VH, Point>> &_target);

  // T(a) = R*a+t
  void transform_constraint(COMISO::LinearConstraint::SVectorNC &_cx,
                            COMISO::LinearConstraint::SVectorNC &_cy,
                            COMISO::LinearConstraint::SVectorNC &_cz,
                            const int _transition_function,
                            const int _translation_idx) const;

  // T(a) = R*a+t as above but t is determined as t = p_to-R*p_from with indices (_bidx_from, _bidx_from+1, _bidx_from+2), ...
  void transform_constraint(COMISO::LinearConstraint::SVectorNC &_cx,
                            COMISO::LinearConstraint::SVectorNC &_cy,
                            COMISO::LinearConstraint::SVectorNC &_cz,
                            const int _transition_function,
                            const int _bidx_from,
                            const int _bidx_to) const;

  // T(a) = R^-1*(a-t)
  void inverse_transform_constraint(COMISO::LinearConstraint::SVectorNC &_cx,
                                    COMISO::LinearConstraint::SVectorNC &_cy,
                                    COMISO::LinearConstraint::SVectorNC &_cz,
                                    const int _transition_function,
                                    const int _translation_idx) const;

  // T(a) = R*a
  void permute_constraint(COMISO::LinearConstraint::SVectorNC &_cx,
                          COMISO::LinearConstraint::SVectorNC &_cy,
                          COMISO::LinearConstraint::SVectorNC &_cz,
                          const int _transition_function) const;

  void check_quantization_path_constraints();

public:
  int number_of_invalid_parametric_tets();

  int number_of_invalid_edge_valencies();

  void export_invalid_parametric_tets(TetMesh &_tetmesh_ref, TetMesh &_tetmesh_param);

  CH add_tet(TetMesh &_tetmesh, Mat3x4d &_P);

  EH add_edge(TetMesh &_tetmesh, Eigen::Matrix<double, 3, 2> &_P);

  void mesh_tet(const CH _ch, Mat3x4d &_P);

  double degenerate_volume_map();

  // interior dihedral angle in param
  double edge_angle(const EH &_eh);

  // interior dihedral angle in mesh
  double mesh_edge_angle(const EH &_eh);

  // as before but in addition to the sum return vector of tet angles
  double edge_angle(const EH &_eh, std::vector<double> &_tet_angles);

  double mesh_edge_angle(const EH &_eh, std::vector<double> &_tet_angles);

  bool param_invalid_edge_valence_valid_tets(const EH &_eh);

  Point face_normal(const Point &_p0, const Point &_p1, const Point &_p2) const;

  int incident_feature_ohalfedges(const VH _vh, std::vector<HEH> &_fheh) const;

  int incident_feature_faces(const VH _vh, std::vector<FH> &_ffh) const;

  inline double clamp(const double _value, const double _min_val, const double _max_val) const;
private:

  // mesh, status and color attribs
  TetMesh &mesh_;
  OpenVolumeMesh::StatusAttrib &mesh_status_;
//  OpenVolumeMesh::ColorAttrib<ACG::Vec4f>& mesh_color_;

  // ##### input properties
  // frame field
  OpenVolumeMesh::CellPropertyT<Mat3x3> frame_cprop_;
  // halfface transitions (for frames -> F0*T = F1)
  OpenVolumeMesh::HalfFacePropertyT<int> transition_hfprop_;
  // integer index of edges (nonzero for singularities)
  OpenVolumeMesh::EdgePropertyT<int> valence_eprop_;
  // feature vertices/edges/faces
  OpenVolumeMesh::VertexPropertyT<int> feature_vprop_;
  OpenVolumeMesh::EdgePropertyT<int> feature_eprop_;
  OpenVolumeMesh::FacePropertyT<int> feature_fprop_;

  // ##### output properties
  // integer-grid map parametrization
  OpenVolumeMesh::CellPropertyT<std::map<OpenVolumeMesh::VertexHandle, Point> > igm_cprop_;

  // ##### intermediate properties
  // face cut surface
  OpenVolumeMesh::FacePropertyT<bool> is_on_cut_fprop_;
  OpenVolumeMesh::CellPropertyT<Vec4i> global_idx_;
  int max_global_idx_;
  int n_helper_variables_;
  OpenVolumeMesh::CellPropertyT<double> cell_weight_;

  // stiffening parameters
  unsigned int max_stiffening_iters_;
  double stiffening_weight_;

  // cut surface decomposed into manifold patches
  std::vector<std::vector<FH> > cut_surface_;

  TransitionQuaternion tq_;

  OptimizerType opt_type_;

  double timelimit_;

  // data for locally-injective optimization
  std::vector<double> x_prev_;

  // quantization constraints from qgp3d if it has been executed, otherwise empty
  std::vector<QPathConstraint> quantization_path_constraints_;

  // json export
  nlohmann::json json_data_;
};




} // namespace AlgoHex
//=============================================================================

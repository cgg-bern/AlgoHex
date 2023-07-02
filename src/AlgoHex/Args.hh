/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

namespace AlgoHex
{

// Global control arguments
struct Args
{
  // Input tetmesh file path
  std::string inFileName;

  // Output hexmesh file path
  std::string outFileName;

  // Save vertex quaternions file path
  std::string vertexQtOutFileName;

  // Load refined tetmesh (with matchings and singular graph) file path
  std::string rfInFileName;

  // Save refined tetmesh (with matchings and singular graph) file path
  std::string rfOutFileName;

  // Save final tetmesh (with matchings and singular graph) file path
  std::string finalTetmeshOutFileName;

  // Load cell quaternions file path
  std::string cellQtInFileName;

  // Save cell quaternions file path
  std::string cellQtOutFileName;

  // Save the final cell quaternions file path
  std::string finalCellQtOutFileName;

  // Save locally non-meshable one-ring tetrahedra file path
  std::string locallyNonMeshableFileName;

  // Load hexex file path
  std::string hexExInFileName;

  // Save the seamless map in .hexex format
  std::string seamlessOutFileName;

  // Save the Integer-grid map in .hexex format
  std::string igmOutFileName;

  // Flag to save files for debugging
  bool save_debug_files = false;

  // Save json file
  std::string jsonOutFileName;

  // Cell-base field
  bool cell_based = false;

  // Paper: Hex Me If You Can
  bool hexme_pipeline = false;

  // Paper: Frame Field Singularity Correction for Automatic Hexahedralization
  bool split_approach = false;

  // Paper: All-Hex Meshing using Singularity-Restricted Field.
  bool collapse_approach = false;

  // Dihedral angle threshold for feature edges
  double dihedral_angle = 70.;

  // Max non-linear octahedral field optimization iteration number
  int max_field_opt_iters = 500;

  // Penalty of normal alignment
  double penalty = 100.;

  // Edge length ratio threshold of average edge length in refining singular regions
  double refine_edge_length_ratio = .6;

  // Target cell number of hexahedral mesh
  int num_hex_cells = 10000;

  // Anisotroy value of hex elements
  double anisotropy = 1.;

  // Stiffening weight in parametrization
  double stiffening_weight = 4.;

  // Solver option. Default is CoMISo, otherwise Gurobi
  int optimizer = 0;

  // Max stiffening iteration number
  int max_stiffening_iters = 1;

  // Timelimit of parametrization
  double param_timelimit = 3600.;

  // Get singular graph without refining the tetmesh
  bool without_refinement = false;

  // Generate locally hex-meshable frame fields
  bool generate_locally_meshable_field = true;

  // Disable optimization of integrable frame field
  bool without_integrable_field_optimization = false;

  // Force to compute feature edges with dihedral angle even if property exists
  bool force_feature_threshold = false;

  // Add full constraints to frames at features
  bool full_constraints = false;

  // scale input mesh
  double scale_input_mesh = 1.0;

  // Number of remeshing in the postprocess to improve the tetmesh quality
  int n_post_remeshing = 3;

  // check for each vertex whether a locally injective map of the 1-ring can be generated
  bool with_local_meshability_test = false;

  // Parameters of singular graph aligned field optimization
  // Iteration number of the second phase
  int second_iters = 15;

  // Iteration number of the first phase
  int first_iters = 13;

  // Inner iteration number of the second phase
  int second_inner_iters = 6;

  // Singular edge alignment weight to field
  double field_alignment_weight = 10.;

  // Complex singular edge shrinkage weight
  double ces_weight = 1.;

  // zipper node shrinkage weight
  double tps_weight = 1.;

  // Regularization weight
  double rgl_weight = 1.;

  // Repulsion energy weight
  double rp_weight = 0.1;

  // Weight of the dual path to follow the field direction
  double dpath_field_weight = 100.;

  // Coefficient of singular arc segment length for filtering zipper nodes
  double tp_seg_coeff = 0.4;

  // Coefficient of the peak distance for filtering zipper nodes
  double tp_peak_coeff = 0.1;

  // Coefficient of the target edge length for filtering zipper nodes
  double abs_seg_coeff = 5.;

  // Edge path w.r.t. field
  bool edge_path_wrt_field = true;

  // Allow merge zipper node to an valence +1 singular arc
  bool merge_zipper_node = true;

  // Push vertice on boundary singular circle to interior
  bool push_bdy_sg_circle = false;

  // Angle difference threshold in quaternion smoothing
  double qtn_angle_thre = 5.;

  // Allow path to end at other zipper nodes
  bool connect_other_tp = true;

  // Fix zipper nodes on shortest singular arcs first (otherwise longest first)
  bool fix_tps_on_shortest_arc_first = false;

  // Maximum searching distance coefficient in finding another zipper node as to the firstly reached boundary
  double max_tp_dist = 1.1;

  // Use truncated newton
  bool truncated_newton = false;

  // After n second iterations, stops the pipeline once all singular/feature vertices are locally hex-meshable
  int early_stop_iter = 3;

  // Save the singularity graph development during the locally meshable field generation in .ovm format
  std::string animationOutFileName;

  // Save the boundary layer of the output hex mesh
  std::string subHexOutFileName;
};

//=============================================================================
} // namespace AlgoHex
//=============================================================================


/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include <CLI/CLI.hpp>
#include "HexMeshing.hh"

#ifdef _WIN32
#  include <windows.h>
#  include <stdlib.h>
#  include <errhandlingapi.h>
#endif


int main(int argc, const char *argv[])
{

#ifdef _WIN32
  SetErrorMode(0);
#endif

  AlgoHex::Args args;

  CLI::App app{"AlgoHex"};
  app.add_flag("--hexme-pipeline", args.hexme_pipeline,
               "Standard frame field based hexmeshing pipeline used in HexMe paper.");
  app.add_flag("--split-paper", args.split_approach,
               "Paper: Frame Field Singularity Correction for Automatic Hexahedralization.");
  app.add_flag("--collapse-paper", args.collapse_approach,
               "Paper: All-Hex Meshing using Singularity-Restricted Field.");

  app.add_option("-i", args.inFileName, "Input tetrahedral mesh in .ovm/.vtk format (optional).");
  app.add_option("-o", args.outFileName,
                 "Output hexhedral mesh in .ovm format (optional, default: *_inFileName +'_hex'+'.ovm').");
  app.add_option("--vertex-quaternion-out-path", args.vertexQtOutFileName,
                 "Save vertex quaternions in .qtn format (optional).");
  app.add_option("--refined-mesh-in-path", args.rfInFileName, "Load refined mesh in .ovm format (optional).");
  app.add_option("--refined-mesh-out-path", args.rfOutFileName, "Save refined mesh in .ovm format (optional).");
  app.add_option("--final-tetmesh-out-path", args.finalTetmeshOutFileName,
                 "Save the tetmesh after the singular graph optimization in .ovm format (optional).");
  app.add_option("--cell-quaternion-in-path", args.cellQtInFileName,
                 "Load cell quaternions in .qtn format (optional).");
  app.add_option("--cell-quaternion-out-path", args.cellQtOutFileName,
                 "Save cell quaternions after the singular graph extraction in .qtn format (optional).");
  app.add_option("--final-cell-quaternion-out-path", args.finalCellQtOutFileName,
                 "Save cell quaternions after the singular graph optimization in .qtn format (optional).");
  app.add_option("--locally-non-meshable-out-path", args.locallyNonMeshableFileName,
                 "Save locally non-meshable one-ring tetrahedra in .ovm format (optional).");
  app.add_option("--hexex-in-path", args.hexExInFileName, "Load HexEx file in .hexex format (optional).");
  app.add_option("--sm-out-path", args.seamlessOutFileName, "Save the seamless map in .hexex format (optional).");
  app.add_option("--igm-out-path", args.igmOutFileName, "Save the Integer-grid map in .hexex format (optional).");
  app.add_flag("--save-debug-files", args.save_debug_files, "Cell base 3D frame field represented as quaternions");

  app.add_flag("--cell-based", args.cell_based, "Cell base 3D frame field represented as quaternions");
  app.add_flag("--force-feature-threshold", args.force_feature_threshold,
               "Refresh feature edges even if the feature edge property exists.");
  app.add_flag("--full-constraints", args.full_constraints, "Fully constrain the frames at features.");
  app.add_option("-d, --dihedral-angle", args.dihedral_angle, "Dihedral angle that defines feature edge.");
  app.add_option("-f, --max-field-opt-iters", args.max_field_opt_iters,
                 "Maximum iteration number of octahedral field optimization.");
  app.add_option("-p, --penalty", args.penalty, "Penalty of the normal alignment.");
  app.add_option("-l, --edge-length-ratio", args.refine_edge_length_ratio,
                 "Ratio of the target edge length threshold to the average tetmesh edge length for singular regions' refinement.");
  app.add_option("-n", args.num_hex_cells, "Number of output hexmesh cells.");
  app.add_option("-a, --anisotropy", args.anisotropy, "Anisotropy of hexmesh elements.");
  app.add_option("-s, --stiffening-weight", args.stiffening_weight, "Stiffening weight in parametrization.");
  app.add_option("--scale-input-mesh", args.scale_input_mesh, "Scale input mesh");
  app.add_option("-z, --optimizer", args.optimizer, "Use which solver to optimize. (0: CoMISo; 1: Gurobi)");
  app.add_option("-m, --max-stiffening-iters", args.max_stiffening_iters,
                 "Maximum iteration number of stiffening in parametrization.");
  app.add_option("-t, --timelimit", args.param_timelimit, "Time limit in seconds for parametrization.");
  app.add_option("-j, --json-out-file", args.jsonOutFileName, "Save json data file in .json format (optional).");
  app.add_flag("--without-refinement", args.without_refinement, "Get singular graph without refining the tetmesh.");

  app.add_flag("--generate-locally-meshable-field", args.generate_locally_meshable_field,
               "Generate locally meshable field.");
  app.add_flag("--without-integrable-field-optimization", args.without_integrable_field_optimization,
               "DO NOT perform an integrability optimization of the frame field.");

  app.add_option("--n-post-remesh", args.n_post_remeshing,
                 "Number of remeshing in the postprocess to improve the tetmesh quality.");
  app.add_flag("--with-local-meshability-test", args.with_local_meshability_test,
               "Checks for each vertex, whether a locally injective map of the 1-ring tets can be generated.");

  //parameters for locally meshable frame fields
  app.add_option("--first-iters", args.first_iters,
                 "Number of iterations of the first phase in singular graph optimization.");
  app.add_option("--second-iters", args.second_iters,
                 "Number of iterations of the second phase in singular graph optimization.");
  app.add_option("--second-inner-iters", args.second_inner_iters,
                 "Number of inner iterations of the second phase in singular graph optimization.");
  app.add_option("--field-alignment-weight", args.field_alignment_weight,
                 "Weight of aligning singular edges to the field.");
  app.add_option("--ces-weight", args.ces_weight,
                 "Weight of shrinking complex singular edge in singular graph optimization.");
  app.add_option("--tps-weight", args.tps_weight,
                 "Weight of shrinking edges at zipper nodes in singular graph optimization.");
  app.add_option("--rgl-weight", args.rgl_weight, "Weight of regularization in singular graph optimization.");
  app.add_option("--rp-weight", args.rp_weight,
                 "Weight of separating twisted singular arc pairs in singular graph optimization.");
  app.add_option("--dpath-weight", args.dpath_field_weight,
                 "Weight of aligning the dual path to the field in fixing zipper nodes.");
  app.add_option("--tp-peak-coeff", args.tp_peak_coeff,
                 "Coefficient of the peak distance for filtering zipper nodes.");
  app.add_option("--tp-seg-coeff", args.tp_seg_coeff,
                 "Coefficient of the singular arc segment length for filtering zipper nodes.");
  app.add_option("--tp-abs-seg-coeff", args.abs_seg_coeff,
                 "Coefficient of the target edge length for filtering zipper nodes.");
  app.add_flag("--epath-wrt-field, --no-epath-wrt-field{false}", args.edge_path_wrt_field,
               "Search edge path w.r.t. the field information in fixing zipper nodes.");
  app.add_flag("--merge-turning-point, --no-merge-turning-point{false}", args.merge_zipper_node,
               "Allow merging zipper nodes to singular arcs of valence +1.");
  app.add_flag("--push-boundary-singular-circle", args.push_bdy_sg_circle,
               "Push vertice on boundary singular circle to interior.");
  app.add_option("--qtn-angle-threshold", args.qtn_angle_thre, "Angle difference threshold in quaternion smoothing.");
  app.add_flag("--connect-other-turning-point, --no-connect-other-turning-point{false}", args.connect_other_tp,
               "Allow the path to end at other zipper nodes.");
  app.add_flag("--fix-turning-point-from-shortest-arc", args.fix_tps_on_shortest_arc_first,
               "Fix zipper nodes on shortest singular arcs first (otherwise longest first).");
  app.add_option("--max-turning-point-dist", args.max_tp_dist,
                 "Maximum searching distance coefficient in finding another zipper node as to the firstly reached boundary.");
  app.add_flag("--truncated-newton", args.truncated_newton,
               "Use truncated newton method to solve the singularity relocation problem.");
  app.add_option("--early-stop-iter", args.early_stop_iter,
                 "After n second iterations, stops the pipeline once all singular/feature vertices are locally hex-meshable.");
  app.add_option("--animation-out-path", args.animationOutFileName,
                 "Save the singularity graph development during the locally meshable field generation in .ovm format (optional).");
  app.add_option("--sub-hexmesh-out-path", args.subHexOutFileName,
                 "Save the boundary layer of the output hex mesh in .ovm format (optional).");

  try
  {
    app.parse(argc, argv);
  }
  catch (const CLI::ParseError &e)
  {
    return app.exit(e);
  }

  //set default outFileName
  if (args.outFileName.empty())
  {
    if (!args.inFileName.empty())
      args.outFileName = args.inFileName.substr(0, args.inFileName.size() - 4) + "_hex.ovm";
    else if (!args.rfInFileName.empty())
      args.outFileName = args.rfInFileName.substr(0, args.rfInFileName.size() - 4) + "_hex.ovm";
    else if (!args.hexExInFileName.empty())
      args.outFileName = args.hexExInFileName.substr(0, args.hexExInFileName.size() - 6) + "_hex.ovm";
    else
    {
      std::cout << "Error: input file is not specified! \nExample usage: ./Build/bin/HexMeshing -i ../demo/HexMeshing/cylinder.ovm \n"
                   <<"Run with '-h' for usage information." << std::endl;
      return -1;
    }
  }

  AlgoHex::hexMeshing(args);

  return 0;
}

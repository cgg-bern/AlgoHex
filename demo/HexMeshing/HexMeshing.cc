/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

//== INCLUDES =================================================================

#include "HexMeshing.hh"

#include <OpenVolumeMesh/FileManager/FileManager.hh>
#include <OpenVolumeMesh/FileManager/VtkColorReader.hh>
#include <AlgoHex/CommonOVMBDecoders/PropertyCodecsSHCoeffs.hh>
#include <OpenVolumeMesh/IO/WriteOptions.hh>
#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <OpenVolumeMesh/IO/ReadOptions.hh>
#include <OpenVolumeMesh/IO/ovmb_read.hh>

#include <AlgoHex/Util/json.hh>

#include <AlgoHex/Stopwatches.hh>


namespace AlgoHex
{
using Point = OpenVolumeMesh::Geometry::Vec3d;
using TetMesh = OVM::TetrahedralGeometryKernel<Point, OVM::TetrahedralMeshTopologyKernel>;

void hexMeshing(Args &args)
{
  if (!args.hexExInFileName.empty() && !args.inFileName.empty())
  {
    hexMeshing_from_seamless_map(args);
  }
  else if(!args.inFileName.empty() || !args.rfInFileName.empty())
  {
    sw::total.reset();
    sw::total.watch().resume();

    //
    TetMesh tetmesh;

    nlohmann::json json_data;
    json_data["InputTetmeshFile"] = args.inFileName;
    if (!args.jsonOutFileName.empty())
      initialize_json_file(args.jsonOutFileName, json_data);

    load_tetmesh(args, tetmesh);

    get_initial_frame_field(args, tetmesh);

    modify_frame_field(args, tetmesh);

    parameterization(args, tetmesh, json_data);

    extract_hexmesh(args, tetmesh);


    // save timings in json file
    sw::total.watch().stop();
    std::cout << sw::total;

    json_data["time_total"] = sw::total.watch().elapsed_ms();
    json_data["time_frame_field_opt"] = sw::field_opt.watch().elapsed_ms();
    json_data["time_singularity_opt"] = sw::lhfg_total.watch().elapsed_ms();
    json_data["time_frame_field_integrability_opt"] = sw::frame_field_int_opt.watch().elapsed_ms();
    json_data["time_parameterization"] = sw::parameterization.watch().elapsed_ms();
    json_data["time_hexex"] = sw::hex_extraction.watch().elapsed_ms();

    std::cerr << "save json data file into " << args.jsonOutFileName << std::endl;
    store_json_file(args, json_data);
  }
}

template<class MeshT>
void load_tetmesh(const Args &args, MeshT &tetmesh)
{
  if (!args.inFileName.empty())
  {
    //identify file format
    auto found = args.inFileName.find_last_of(".");
    auto file_type = args.inFileName.substr(found + 1);

    if (file_type == "vtk")
    {
      OpenVolumeMesh::Reader::VtkColorReader vtkfm;
      vtkfm.readFile(args.inFileName, tetmesh, true, true);
    }
    else if (file_type == "ovm")
    {
      OpenVolumeMesh::IO::FileManager fm;
      fm.readFile(args.inFileName, tetmesh);
    } else if (file_type == "ovmb")
    {
      read_ovmb_file(args.inFileName, tetmesh);
    }
    else
    {
      std::cerr << "Error: the file type is not supported!" << std::endl;
      return;
    }
  }
}

template<class MeshT>
void get_initial_frame_field(const Args &args, MeshT &tetmesh)
{
  //if the input is a pure tetmesh
  //1. compute initial field and extract singular graph
  if (!args.inFileName.empty())
  {
    initialize_feature_properties(args, tetmesh);

    auto feature_edges = get_feature_edges(tetmesh);

    AlgoHex::FaceNormalCache normal_cache{tetmesh};
    if (!args.cell_based)
    {
      AlgoHex::FieldConstraintGenerator constraint_gen{tetmesh, normal_cache};
      if (args.full_constraints)
      {
        constraint_gen.add_full_constraints_from_feature_edges(feature_edges);
      }
      constraint_gen.add_partial_constraints_from_feature_edges(feature_edges);

      AlgoHex::SmoothOctahedralFieldGeneratorT<TetMesh> sofg(tetmesh);
      sofg.solve_spherical_harmonic_coefficients(args.penalty);
      sofg.project();
      sofg.iterate(args.max_field_opt_iters, args.penalty);
      sofg.print_timings();

      //save vertex quaternion
      if (!args.vertexQtOutFileName.empty())
      {
        std::cout << "Saving vertex quaternions..." << std::endl;
        sofg.save_vertex_quaternion(args.vertexQtOutFileName);
      }
    }
    else
    {
      AlgoHex::FieldConstraintGeneratorCellBased constraint_gen{tetmesh, normal_cache};
      constraint_gen.add_partial_constraints(feature_edges);

      AlgoHex::SmoothOctahedralFieldGeneratorCellBasedT<TetMesh> sofg(tetmesh);
      sofg.solve_spherical_harmonic_coefficients(args.penalty);
      sofg.project();
      sofg.iterate(args.max_field_opt_iters, args.penalty);
      sofg.print_timings();
    }

    AlgoHex::SingularGraphExtractionT<TetMesh> sge(tetmesh);
    if (!args.cell_based)
    {
      if (args.without_refinement)
        sge.get_singular_edges();
      else
        sge.get_singular_edges_with_refinement(args.refine_edge_length_ratio);
    }
    else
    {
      sge.get_transitions();
      sge.get_edges_valence();
    }

    std::cout << "Tetmesh cells: " << tetmesh.n_cells() << std::endl;
    std::cout << "Tetmesh faces: " << tetmesh.n_faces() << std::endl;
    std::cout << "Tetmesh edges: " << tetmesh.n_edges() << std::endl;
    std::cout << "Tetmesh vertices: " << tetmesh.n_vertices() << "\n" << std::endl;


    //save mesh, matching and singular graph
    if (!args.rfOutFileName.empty())
    {
      //do not save properties
      auto vqt = tetmesh.template request_vertex_property<Quaternion>("vertex quaternion");
      tetmesh.set_persistent(vqt, false);
      auto shc = tetmesh.template request_vertex_property<SHCoeffs>(
              "scaled vertex spherical harmonic coefficients");
      tetmesh.set_persistent(shc, false);

      std::cout << "Saving tetmesh, field and singular graph..." << std::endl;
      write_ovmb_file(args.rfOutFileName, tetmesh);

      tetmesh.set_persistent(vqt, true);
      tetmesh.set_persistent(shc, true);

      //save cell quaternion
      if (!args.cellQtOutFileName.empty())
      {
        std::cout << "Saving cell quaternions..." << std::endl;
        sge.save_cell_quaternions(args.cellQtOutFileName);
      }
    }
  }
  else if (!args.rfInFileName.empty())
  { //read singular graph, field and tetmesh from data
    //read ovmb file
    read_ovmb_file(args.rfInFileName, tetmesh);

    //read ovm file and load quaternions
//        fm.readFile(args.rfInFileName, tetmesh);
//        AlgoHex::SingularGraphExtractionT<TetMesh> sge(tetmesh);
//        sge.load_cell_quaternions(args.cellQtInFileName);
  }
  else
  {
    std::cerr << "No input is given!" << std::endl;
    return;
  }
}

template<class MeshT>
void modify_frame_field(Args &args, MeshT &tetmesh)
{
  if (args.split_approach)
  {
    AlgoHex::SingularGraphExtractionT<TetMesh> sge(tetmesh);
    AlgoHex::FrameFieldOptimizationT<TetMesh> ffo(tetmesh);
    ffo.align_quaternions_to_feature(100);
    sge.get_transitions();
    sge.get_edges_valence();

    AlgoHex::SplitApproachT<TetMesh> sa(tetmesh);
    sa.split_approach();
    auto zz_free = sa.check_zigzag();
    auto n_cpe = sa.n_complex_singular_edges();
    if (!zz_free || n_cpe > 0)
      std::cerr << "Zigzags or compound singular edges exist!" << std::endl;
  }
  if (args.collapse_approach)
  {
    AlgoHex::SingularGraphExtractionT<TetMesh> sge(tetmesh);
    AlgoHex::FrameFieldOptimizationT<TetMesh> ffo(tetmesh);
    ffo.align_quaternions_to_feature(100);
    sge.get_transitions();
    sge.get_edges_valence();

    AlgoHex::CollapseApproachT<TetMesh> ca(tetmesh);
    ca.collapse_approach();
    auto n_cpe = ca.n_complex_singular_edges();

    if (n_cpe > 0)
      std::cerr << "Compound singular edges exist!" << std::endl;
  }

  if (args.hexme_pipeline || args.split_approach || args.collapse_approach)
  {
    args.generate_locally_meshable_field = false;
    args.without_integrable_field_optimization = true;
  }

  if (args.generate_locally_meshable_field)
  {
    TetMesh sg_mesh;
    LocallyMeshableFieldGenerationT <TetMesh> lmfg(tetmesh);
    lmfg.set_export_path(args.outFileName.substr(0, args.outFileName.size() - 7));
    set_locally_hexmeshable_field_coeffs(args, lmfg);

    lmfg.generate_locally_meshable_field(sg_mesh);

    //save animation
    if(!args.animationOutFileName.empty())
    {
      OpenVolumeMesh::IO::FileManager fm;
      fm.writeFile(args.animationOutFileName, sg_mesh);
    }
  }
  //save final tetmesh
  if (!args.finalTetmeshOutFileName.empty())
    write_ovmb_file(args.finalTetmeshOutFileName, tetmesh);
  if (!args.finalCellQtOutFileName.empty())
    save_cell_quaternions(args.finalCellQtOutFileName, tetmesh);
}

template<class MeshT>
void parameterization(const Args &args, MeshT &tetmesh, nlohmann::json &json_data)
{
  // hack ---> should be earlier
  if (args.scale_input_mesh != 1.0)
  {
    std::cerr << "scale input mesh with factor " << args.scale_input_mesh << std::endl;
    scale_mesh(tetmesh, args.scale_input_mesh);
  }

  OpenVolumeMesh::IO::FileManager fm;

  //parameterize
  OpenVolumeMesh::StatusAttrib status(tetmesh);
  AlgoHex::Parametrization3DT<TetMesh> param(tetmesh, status);
  param.import_frames_from_quaternion_field();
  param.scale_frames_to_target_complexity(args.num_hex_cells);

  AlgoHex::FrameFieldOptimizer3DT<TetMesh> ffopt(tetmesh);

  // check local meshability
  check_local_meshability(args, tetmesh, json_data);

  // Frame Field Optimization
  if (!args.without_integrable_field_optimization)
  {
    if(args.save_debug_files)
    {
      fm.writeFile(args.finalTetmeshOutFileName + "before_tet_deform.ovm", tetmesh);
    }

    // ToDo: integrable field optimization before dual optimization?
    ffopt.optimize_integrable_regularized(0.1, 1e4, 1e4);
    ffopt.optimize_integrable_regularized(0.1, 1e3, 1e5);
    ffopt.optimize_integrable_regularized(0.1, 1e2, 1e6);
    ffopt.optimize_integrable_regularized(0.1, 1e1, 1e7);
    ffopt.optimize_integrable_regularized(0.1, 0.0, 1e8);
    ffopt.optimize_integrable_regularized(0.1, 0.0, 1e9);

    if(args.save_debug_files)
    {
      fm.writeFile(args.finalTetmeshOutFileName + "after_tet_deform.ovm", tetmesh);
    }

    json_data["FrameFieldOptimization"] = ffopt.json_data();
    store_json_file(args, json_data);
  }

  ffopt.check_valence_consistency();
  ffopt.check_frame_rotation_angles();

  // set parameters
  if (args.optimizer == 0)
    param.set_optimizer(AlgoHex::Parametrization3DT<TetMesh>::OPT_COMISO);
  else if (args.optimizer == 1)
    param.set_optimizer(AlgoHex::Parametrization3DT<TetMesh>::OPT_GUROBI);

  param.set_timelimit(args.param_timelimit);
  // set stiffening params
  param.set_stiffening_weight(args.stiffening_weight);
  param.set_max_stiffening_iters(args.max_stiffening_iters);

  std::pair<int, int> invalid_elements_seamless;

  // parametrize seamless (---> check success of integrable field optimization)
  invalid_elements_seamless = param.parametrize_complete(args.num_hex_cells, args.anisotropy, false);

  if(!args.seamlessOutFileName.empty()) {
    // save seamless map
    auto parameters = tetmesh.template request_cell_property<std::map<OpenVolumeMesh::VertexHandle, TetMesh::PointT>>(
            "Parametrization");
    HexEx::HexExtractor hexExtractor(tetmesh, parameters);
    hexExtractor.writeToFile(args.seamlessOutFileName);
    attach_feature_information(tetmesh, args.seamlessOutFileName);
  }

  //store to json file
  json_data["Parametrization"] = param.json_data();
  store_json_file(args, json_data);

  if(args.save_debug_files)
  {
    // debug output
    if (invalid_elements_seamless.first != 0 || invalid_elements_seamless.second != 0)
    {
      TetMesh invalid_ref, invalid_param;
      param.export_invalid_parametric_tets(invalid_ref, invalid_param);

      std::stringstream ss_ref, ss_param;
      ss_ref << args.outFileName << "_invalid_ref.ovm";
      ss_param << args.outFileName << "_invalid_param.ovm";

      fm.writeFile(ss_ref.str(), invalid_ref);
      fm.writeFile(ss_param.str(), invalid_param);
    }
  }

  // generate integer-grid map
  if (invalid_elements_seamless.first == 0 && invalid_elements_seamless.second == 0)
  {
    std::cout << "Info: Parameterizing using robust quantization-based pipeline." << std::endl;
    // update frames to benefit from valid seamless map when generating IGM
    param.import_frames_from_parametrization();
    param.parametrize_robust_quantization(args.num_hex_cells);
  }
  else
  {
    std::cout << "WARNING: Parameterizing without quantization (non-robust pipeline)." << std::endl;
    // fallback if valid seamless map could not be found, or Gurobi not available
    param.parametrize(args.num_hex_cells, args.anisotropy, true);
  }

  json_data["Parametrization"] = param.json_data();
  store_json_file(args, json_data);

  if(args.save_debug_files)
  {
    //debug
    TetMesh viz_quantization_constraint_paths;
    param.export_quantization_constraint_paths(viz_quantization_constraint_paths);
    auto viz_mesh_name = args.outFileName.substr(0, args.outFileName.size() - 7) + "viz_quant.ovm";
    fm.writeFile(viz_mesh_name, viz_quantization_constraint_paths);
  }
}

template<class MeshT>
void extract_hexmesh(const Args &args, MeshT &tetmesh)
{
  ScopedStopWatch sw(sw::hex_extraction);

  std::cout << "\n#####Extracting hexmesh..." << std::endl;
  auto parameters = tetmesh.template request_cell_property<std::map<OpenVolumeMesh::VertexHandle, TetMesh::PointT>>(
          "Parametrization");
  HexEx::HexExtractor hexExtractor(tetmesh, parameters);

  //save IGM
  if (!args.igmOutFileName.empty())
    hexExtractor.writeToFile(args.igmOutFileName);

  hexExtractor.extract();
  HexEx::HexahedralMesh hexmesh;
  hexExtractor.getHexMesh(hexmesh, true);

  std::cout << "Hexmesh cells: " << hexmesh.n_cells() << std::endl;
  std::cout << "Hexmesh faces: " << hexmesh.n_faces() << std::endl;
  std::cout << "Hexmesh edges: " << hexmesh.n_edges() << std::endl;
  std::cout << "Hexmesh vertices: " << hexmesh.n_vertices() << std::endl;

  auto hex_vol = hexmesh_volume(hexmesh);
  auto tet_vol = tetmesh_volume(tetmesh);
  std::cout << "#####Hexmesh volume / tetmesh voluem: " << hex_vol << "/" << tet_vol << " "
                << hex_vol / tet_vol << std::endl;


  //save hexmesh
  if (!args.outFileName.empty())
  {
    OpenVolumeMesh::IO::FileManager fm;
    fm.writeFile(args.outFileName, hexmesh);
  }
  //save the boundary layer
  if(!args.subHexOutFileName.empty())
  {
    OpenVolumeMesh::IO::FileManager fm;
    HexEx::HexahedralMesh sub_mesh;
    create_sub_hexmesh(hexmesh, sub_mesh);
    auto sub_mesh_name = args.outFileName.substr(0, args.outFileName.size() - 7) + "sub_hex.ovm";
    fm.writeFile(sub_mesh_name, sub_mesh);
  }
}

void hexMeshing_from_seamless_map(const Args &args)
{
  if (args.hexExInFileName.empty())
  {
    std::cerr << "Please specify hexex file!" << std::endl;
    return;
  }
  if (args.inFileName.empty())
  {
    std::cerr << "Please specify tetmesh file!" << std::endl;
    return;
  }

  nlohmann::json json_data;

  // read tetmesh
  OpenVolumeMesh::IO::FileManager fm;

  TetMesh tetmesh;
  fm.readFile(args.inFileName, tetmesh);

  // read hexex file
  auto parameters = tetmesh.request_cell_property<std::map<OpenVolumeMesh::VertexHandle, TetMesh::PointT>>(
          "Parametrization");
  tetmesh.set_persistent(parameters);

  HexEx::HexExtractor hexExtractor(args.hexExInFileName);
  TetMesh tmp_mesh;
  hexExtractor.getInputMesh(tmp_mesh);

  //copy pm
  auto tmp_parameters = tmp_mesh.request_cell_property<std::map<OpenVolumeMesh::VertexHandle, TetMesh::PointT>>(
          "Parametrization");
  for (auto c_it = tmp_mesh.cells_begin(); c_it != tmp_mesh.cells_end(); ++c_it)
    for (auto cv_it = tmp_mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it)
    {
      parameters[*c_it][*cv_it] = tmp_parameters[*c_it][*cv_it];
    }


  //parameterize
  OpenVolumeMesh::StatusAttrib status(tetmesh);
  AlgoHex::Parametrization3DT<TetMesh> param(tetmesh, status);
  // set parameters
  if (args.optimizer == 0)
    param.set_optimizer(AlgoHex::Parametrization3DT<TetMesh>::OPT_COMISO);
  else if (args.optimizer == 1)
    param.set_optimizer(AlgoHex::Parametrization3DT<TetMesh>::OPT_GUROBI);

  param.set_timelimit(args.param_timelimit);
  // set stiffening params
  param.set_stiffening_weight(args.stiffening_weight);
  param.set_max_stiffening_iters(args.max_stiffening_iters);

  // update frames to benefit from valid seamless map when generating IGM
  param.import_frames_from_parametrization();

  //re-run seamless map generation
  auto invalid_elements_seamless = param.parametrize_complete(args.num_hex_cells, args.anisotropy, false);

  if (invalid_elements_seamless.first == 0 && invalid_elements_seamless.second == 0 && COMISO_GUROBI_AVAILABLE)
  {
    param.parametrize_robust_quantization(args.num_hex_cells);
  }
  json_data["Parametrization"] = param.json_data();
  store_json_file(args, json_data);

  std::cout << "\n#####Extracting hexmesh..." << std::endl;
  HexEx::HexahedralMesh hexmesh;
  HexEx::HexExtractor hexExtractor2(tetmesh, parameters);
  hexExtractor2.extract();
  hexExtractor2.getHexMesh(hexmesh, true);
  if (!args.outFileName.empty())
  {
    fm.writeFile(args.outFileName, hexmesh);
  }
  std::cout << "Hexmesh cells: " << hexmesh.n_cells() << std::endl;
  std::cout << "Hexmesh faces: " << hexmesh.n_faces() << std::endl;
  std::cout << "Hexmesh edges: " << hexmesh.n_edges() << std::endl;
  std::cout << "Hexmesh vertices: " << hexmesh.n_vertices() << std::endl;
}

template<class MeshT>
void check_local_meshability(const Args &args, MeshT &tetmesh, nlohmann::json &json_data) {
  if (args.with_local_meshability_test || !args.locallyNonMeshableFileName.empty())
  {
    std::cerr << "******** CHECK LOCAL MESHABILITY ************" << std::endl;
    AlgoHex::LocalMeshabilityChecker lmc(tetmesh);
    if (!args.locallyNonMeshableFileName.empty())
    {
      lmc.enable_save_locally_non_meshable(args.locallyNonMeshableFileName);
      if (args.generate_locally_meshable_field)
        lmc.check_local_meshability();
      else
        lmc.check_local_meshability(true);
    }
    else
    {
      if (args.generate_locally_meshable_field)
        lmc.check_local_meshability();
      else
        lmc.check_local_meshability(true);
    }

    std::cerr << "******** CHECK EDGE LOCAL MESHABILITY ************" << std::endl;

    lmc.check_edge_local_meshability(true);

    json_data["LocalMeshability"] = lmc.json_data();
    store_json_file(args, json_data);


    std::cerr << "#non-meshable = " << json_data["LocalMeshability"]["non-meshable"] << std::endl;
    std::cerr << "******** END LOCAL MESHABILITY ************" << std::endl;
  }
}

template<class MeshT>
void
set_locally_hexmeshable_field_coeffs(Args &args, LocallyMeshableFieldGenerationT <MeshT> &lmfg)
{
  if (args.with_local_meshability_test)
    args.field_alignment_weight = 1.;

  lmfg.first_iters_ = args.first_iters;
  lmfg.second_iters_ = args.second_iters;
  lmfg.second_inner_iteers_ = args.second_inner_iters;
  lmfg.dpath_weight_ = args.dpath_field_weight;
  lmfg.peak_coeff_ = args.tp_peak_coeff;
  lmfg.seg_coeff_ = args.tp_seg_coeff;
  lmfg.abs_seg_coeff_ = args.abs_seg_coeff;
  lmfg.field_alignment_weight_ = args.field_alignment_weight;
  lmfg.ce_shrink_weight_ = args.ces_weight;
  lmfg.tp_shrink_weight_ = args.tps_weight;
  lmfg.regularize_weight_ = args.rgl_weight;
  lmfg.repulsion_weight_ = args.rp_weight;
  lmfg.edge_path_wrt_field_ = args.edge_path_wrt_field;
  lmfg.merge_zipper_node_ = args.merge_zipper_node;
  lmfg.push_singular_circle_ = args.push_bdy_sg_circle;
  lmfg.n_post_remeshing_ = args.n_post_remeshing;
  lmfg.qtn_angle_thre_ = args.qtn_angle_thre;
  lmfg.connect_other_zipper_node_ = args.connect_other_tp;
  lmfg.fix_tps_on_longest_arc_first_ = !args.fix_tps_on_shortest_arc_first;
  lmfg.max_tp_dist_ = args.max_tp_dist;
  lmfg.truncated_newton_ = args.truncated_newton;
  lmfg.early_stop_iter_ = args.early_stop_iter;

}


template<class MeshT>
void initialize_feature_properties(const Args &args, MeshT &tetmesh)
{
  if (!args.force_feature_threshold && tetmesh.template edge_property_exists<int>("AlgoHex::FeatureEdges"))
  {
    std::cout << "Found AlgoHex feature edges." << std::endl;
    if (!tetmesh.template face_property_exists<int>("AlgoHex::FeatureFaces"))
    {
      auto feature_fprop = tetmesh.template request_face_property<int>("AlgoHex::FeatureFaces", 0);
      for (const auto fhi: tetmesh.faces())
        if (tetmesh.is_boundary(fhi))
          feature_fprop[fhi] = 1;
        else
          feature_fprop[fhi] = 0;

      tetmesh.set_persistent(feature_fprop, true);
    }
    MeshPropertiesT <MeshT> mp(tetmesh);
    mp.initialize_feature_vertex_property();
  }
  else if (!args.force_feature_threshold && tetmesh.template edge_property_exists<int>("edge_colors"))
  {
    std::cout << "Found prescribed feature colors (from vtk file)." << std::endl;
    if (!tetmesh.template face_property_exists<int>("face_colors"))
      std::cerr << "Error: no input feature face colors" << std::endl;

    MeshPropertiesT <MeshT> mp(tetmesh);
    mp.import_feature_properties();
    mp.initialize_feature_vertex_property();
  }
  else
  {
    std::cout << "No prescribed feature tags found. Set boundary faces to feature faces and edges w.r.t. dihedral angle threshold "<<args.dihedral_angle<<" to feature edges." << std::endl;
    //Feature faces
    auto feature_fprop = tetmesh.template request_face_property<int>("AlgoHex::FeatureFaces", 0);
    for (const auto fhi: tetmesh.faces())
      if (tetmesh.is_boundary(fhi))
        feature_fprop[fhi] = 1;
      else
        feature_fprop[fhi] = 0;

    tetmesh.set_persistent(feature_fprop, true);


    //Feature edges
    AlgoHex::FaceNormalCache normal_cache{tetmesh};
    AlgoHex::DihedralFeatureDetector feature_detector{tetmesh, normal_cache};
    auto feature_edges = feature_detector.compute(args.dihedral_angle);

    //set feature edge property
    auto feature_edge_prop = tetmesh.template request_edge_property<int>("AlgoHex::FeatureEdges", 0);
    tetmesh.set_persistent(feature_edge_prop, true);
    for (const auto eh: feature_edges)
      feature_edge_prop[eh] = 1;

    MeshPropertiesT <MeshT> mp(tetmesh);
    mp.initialize_feature_vertex_property();
  }
}

template<class MeshT>
std::vector<EH> get_feature_edges(MeshT &tetmesh)
{
  std::vector<EH> feature_edges;
  auto feature_edge_prop = tetmesh.template request_edge_property<int>("AlgoHex::FeatureEdges");

  //split for singular alignment
  AlgoHex::SplitHelperT<MeshT>::split_for_dof(tetmesh);

  feature_edges.clear();
  for (const auto eh: tetmesh.edges())
    if (feature_edge_prop[eh])
      feature_edges.push_back(eh);

  return feature_edges;
}


double hexmesh_volume(const HexEx::HexahedralMesh &_hexmesh)
{
  using Point = HexEx::HexahedralMesh::PointT;
  double vol(0);

  for (const auto ch: _hexmesh.cells())
  {
    //cell barycenter
    Point cb(0);
    for (const auto cvh: _hexmesh.cell_vertices(ch))
      cb += _hexmesh.vertex(cvh);
    cb /= 8.;

    //halffaces of the cell
    auto hfhs = _hexmesh.cell(ch).halffaces();
    for (const auto hfh: hfhs)
    {
      Point fb(0);
      for (const auto hfvh: _hexmesh.halfface_vertices(hfh))
        fb += _hexmesh.vertex(hfvh);
      fb /= 4.;

      auto hehs = _hexmesh.halfface(hfh).halfedges();
      for (const auto heh: hehs)
      {
        auto p2 = _hexmesh.vertex(_hexmesh.halfedge(heh).to_vertex());
        auto p3 = _hexmesh.vertex(_hexmesh.halfedge(heh).from_vertex());

        vol += 1. / 6. * ((fb - cb) % (p2 - cb) | (p3 - cb));
      }
    }
  }

  return vol;
}

void initialize_json_file(const std::string &json_path, nlohmann::json &json_data)
{
  json_data["LocalMeshability"]["percentage_meshable_vertices"] = -1;
  json_data["LocalMeshability"]["non-meshable"] = -1;
  json_data["LocalMeshability"]["all vertices"] = -1;
  json_data["LocalMeshability"]["singular nodes"] = -1;
  json_data["LocalMeshability"]["zipper nodes"] = -1;
  json_data["LocalMeshability"]["singular vertices"] = -1;
  json_data["LocalMeshability"]["feature vertices"] = -1;
  json_data["LocalMeshability"]["complex singular edges"] = -1;
  json_data["LocalMeshability"]["singular edges"] = -1;

  json_data["Parametrization"]["final_energy"] = -1;
  json_data["Parametrization"]["parametric_volume"] = -1;
  json_data["Parametrization"]["n_invalid_param_tets"] = -1;
  json_data["Parametrization"]["n_invalid_valencies"] = -1;
  json_data["Parametrization"]["valid_volume"] = -1;

  json_data["Parametrization"]["final_energy_seamless"] = -1;
  json_data["Parametrization"]["parametric_volume_seamless"] = -1;
  json_data["Parametrization"]["n_invalid_param_tets_seamless"] = -1;
  json_data["Parametrization"]["n_invalid_valencies_seamless"] = -1;
  json_data["Parametrization"]["valid_volume_seamless"] = -1;

  json_data["time_total"] = -1;
  json_data["time_frame_field_opt"] = -1;
  json_data["time_singularity_opt"] = -1;
  json_data["time_frame_field_integrability_opt"] = -1;

  std::cerr << "Initialize json data into " << json_path << std::endl;
  std::ofstream jout(json_path);
  jout << json_data;
  jout.close();
}

void create_sub_hexmesh(const HexEx::HexahedralMesh &_hexmesh, HexEx::HexahedralMesh &_sub_mesh)
{
  std::map<VH, VH> vh_old_to_new;
  for (auto c_it = _hexmesh.c_iter(); c_it.valid(); ++c_it)
  {
    if (_hexmesh.is_boundary(*c_it))
    {
      for (auto cv_it = _hexmesh.cv_iter(*c_it); cv_it.valid(); ++cv_it)
      {
        if (vh_old_to_new.find(*cv_it) == vh_old_to_new.end())
        {
          auto vh_new = _sub_mesh.add_vertex(_hexmesh.vertex(*cv_it));
          vh_old_to_new.insert(std::make_pair(*cv_it, vh_new));
        }
      }
    }
  }

  for (auto c_it = _hexmesh.c_iter(); c_it.valid(); ++c_it)
  {
    if (_hexmesh.is_boundary(*c_it))
    {
      std::vector<HFH> sc_hfs;
      std::vector<HFH> hfs = _hexmesh.cell(*c_it).halffaces();
      for (const auto hfi: hfs)
      {
        std::vector<VH> sc_vhs;

        for (auto hfv_it = _hexmesh.hfv_iter(hfi); hfv_it.valid(); ++hfv_it)
          sc_vhs.push_back(vh_old_to_new[*hfv_it]);

        FH sc_fh = _sub_mesh.add_face(sc_vhs);
        sc_hfs.push_back(_sub_mesh.halfface_handle(sc_fh, 0));
      }

      _sub_mesh.add_cell(sc_hfs);
    }
  }

  std::cout << "Subhexmesh cells: " << _sub_mesh.n_cells() << std::endl;
  std::cout << "Subhexmesh faces: " << _sub_mesh.n_faces() << std::endl;
  std::cout << "Subhexmesh edges: " << _sub_mesh.n_edges() << std::endl;
  std::cout << "Subhexmesh vertices: " << _sub_mesh.n_vertices() << std::endl;
}

void store_json_file(const Args &_args, const nlohmann::json &_json_data)
{
  if (!_args.jsonOutFileName.empty())
  {
    std::ofstream jout(_args.jsonOutFileName, std::ofstream::trunc);
    jout << _json_data;
    jout.close();
  }
}


template<class MeshT>
void attach_feature_information(MeshT &mesh, const std::string &_filename)
{
  std::ofstream os(_filename, std::ofstream::app);

// store feature information if available
  std::vector<OpenVolumeMesh::VertexHandle> ft_vertices;
  std::vector<OpenVolumeMesh::EdgeHandle> ft_edges;
  std::vector<OpenVolumeMesh::FaceHandle> ft_faces;

  // collect
  if (mesh.template vertex_property_exists<int>("AlgoHex::FeatureVertices"))
  {
    auto feature_vprop = mesh.template request_vertex_property<int>("AlgoHex::FeatureVertices");
    for (const auto vh: mesh.vertices())
      if (feature_vprop[vh] > 0)
        ft_vertices.push_back(vh);
  }

  if (mesh.template edge_property_exists<int>("AlgoHex::FeatureEdges"))
  {
    auto feature_eprop = mesh.template request_edge_property<int>("AlgoHex::FeatureEdges");
    for (const auto eh: mesh.edges())
      if (feature_eprop[eh] > 0)
        ft_edges.push_back(eh);
  }

  if (mesh.template face_property_exists<int>("AlgoHex::FeatureFaces"))
  {
    auto feature_fprop = mesh.template request_face_property<int>("AlgoHex::FeatureFaces");
    for (const auto fh: mesh.faces())
      if (feature_fprop[fh] > 0)
        ft_faces.push_back(fh);
  }

  std::cerr << "store #feature_vertices = " << ft_vertices.size() << std::endl
            << "store #feature_edges    = " << ft_edges.size() << std::endl
            << "store #feature_faces    = " << ft_faces.size() << std::endl;

  os << ft_vertices.size() << " " << ft_edges.size() << " " << ft_faces.size() << std::endl;

  for (const auto vh: ft_vertices)
    os << vh.idx() << std::endl;

  for (const auto eh: ft_edges)
    os << mesh.edge(eh).from_vertex().idx() << " " << mesh.edge(eh).to_vertex().idx() << std::endl;

  for (const auto fh: ft_faces)
  {
    for (int i = 0; i < 3; ++i)
      os << mesh.halfedge(mesh.face(fh).halfedges()[i]).to_vertex().idx() << " ";
    os << std::endl;
  }
}

/*
void read_feature_information(...)
{
// number of feature vertices/edges/faces stored in file
  int n_ftv(0), n_fte(0), n_ftf(0);

  is >> n_ftv >> n_fte >> n_ftf;

  std::cerr << "read #feature_vertices = " << n_ftv
            << "read #feature_edges = " << n_fte
            << "read #feature_vertices = " << n_ftf << std::endl;

// request/add feature properties
  auto feature_vprop = mesh.template request_vertex_property<bool>("AlgoHex::FeatureVertices", false);
  auto feature_eprop = mesh.template request_edge_property<int>("AlgoHex::FeatureEdges", false);
  auto feature_fprop = mesh.template request_face_property<int>("AlgoHex::FeatureFaces", false);

  mesh.set_persistent(feature_vprop, true);
  mesh.set_persistent(feature_eprop, true);
  mesh.set_persistent(feature_fprop, true);

  for (int i = 0; i < n_ftv; ++i) {
    int vidx;
    is >> vidx;
    feature_vprop[OpenVolumeMesh::VertexHandle(vidx)] = true;
  }

  for (int i = 0; i < n_fte; ++i) {
    int v0idx, v1idx;
    is >> v0idx >> v1idx;

    OpenVolumeMesh::HalfEdgeHandle heh = mesh.halfedge(OpenVolumeMesh::VertexHandle(v0idx),
                                                       OpenVolumeMesh::VertexHandle(v1idx));
    if (!heh.is_valid())
      std::cerr << "ERROR: could not obtain feature edge stored in .hexex file " << v0idx << " " << v1idx << std::endl;
    else {
      auto eh = mesh.edge_handle(heh);
      feature_eprop[eh] = true;
    }
  }

  for (int i = 0; i < n_ftf; ++i) {
    int v0idx, v1idx, v2idx;
    is >> v0idx >> v1idx >> v2idx;

// map vertex indices
    std::vector<OpenVolumeMesh::VertexHandle> vhs;
    vhs.push_back(OpenVolumeMesh::VertexHandle(v0idx));
    vhs.push_back(OpenVolumeMesh::VertexHandle(v1idx));
    vhs.push_back(OpenVolumeMesh::VertexHandle(v2idx));

// get corresponding halfface in original mesh
    OpenVolumeMesh::HalfFaceHandle hfh = mesh.halfface(vhs);
    if (!hfh.is_valid()) std::cerr << "ERROR: could not obtain feature face stored in .hexex file" << std::endl;
    else {
      OpenVolumeMesh::FaceHandle fh = mesh.face_handle(hfh);
      feature_fprop[fh] = true;
    }
  }
}
 */

template<class MeshT>
void load_cell_quaternions(const Args &_args, MeshT &tetmesh)
{
  std::ifstream fin;
  fin.open(_args.cellQtInFileName);

  OpenVolumeMesh::CellPropertyT <Eigen::Quaterniond> quaternion_prop =
          tetmesh.template request_cell_property<Eigen::Quaterniond>("FrameFieldQuaternions");
  tetmesh.set_persistent(quaternion_prop, true);

  auto num_c = tetmesh.n_cells();
  for (unsigned int i = 0; i < num_c; ++i)
  {
    std::vector<double> temp_q(4);
    for (int j = 0; j < 4; j++)
    {
      fin >> temp_q[j];
    }
    Quaternion quaternion(temp_q[0], temp_q[1], temp_q[2], temp_q[3]);
    quaternion_prop[CH(i)] = quaternion;
  }

  fin.close();
}

template<class MeshT>
void save_cell_quaternions(const std::string &_filename, MeshT &_mesh)
{
  std::ofstream f_write(_filename);
  auto quaternion_prop = _mesh.template request_cell_property<Eigen::Quaterniond>("FrameFieldQuaternions");

  for (const auto ch: _mesh.cells())
  {
    // write to file
    f_write << quaternion_prop[ch].w() << " ";
    f_write << quaternion_prop[ch].x() << " ";
    f_write << quaternion_prop[ch].y() << " ";
    f_write << quaternion_prop[ch].z() << " ";
  }

  f_write.close();
}

template<typename MeshT>
void write_ovmb_file(const std::string &filename, MeshT &mesh)
{
  //ovmb
  OVM::IO::WriteOptions wo;
  OVM::IO::ReadOptions ro;
  OVM::IO::PropertyCodecs propCodecs = OVM::IO::g_default_property_codecs;
  OVM::IO::register_eigen_codecs(propCodecs);
  OVM::IO::register_algohex_codecs(propCodecs);

  std::ofstream off(filename.c_str(), std::ios::binary);
  OVM::IO::ovmb_write(off, mesh, wo, propCodecs);
  off.close();
}

template<typename MeshT>
void read_ovmb_file(const std::string &filename, MeshT &mesh)
{
  OVM::IO::ReadOptions ro;
  OVM::IO::PropertyCodecs propCodecs = OVM::IO::g_default_property_codecs;
  OVM::IO::register_eigen_codecs(propCodecs);
  OVM::IO::register_algohex_codecs(propCodecs);

  std::ifstream iff(filename.c_str(), std::ios::binary);
  OVM::IO::ovmb_read(iff, mesh, ro, propCodecs);
  iff.close();
}

//=============================================================================
} // namespace AlgoHex
//=============================================================================

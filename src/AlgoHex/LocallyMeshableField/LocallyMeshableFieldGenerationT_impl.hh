/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define LOCALLYMESHABLEFIELDGENERATIONT_C

#include "LocallyMeshableFieldGenerationT.hh"
#include <AlgoHex/LocallyMeshableField/MeshOptimizationT.hh>
#include <AlgoHex/LocallyMeshableField/FrameFieldOptimizationT.hh>
#include <AlgoHex/LocallyMeshableField/SplitHelperT.hh>
#include <AlgoHex/LocallyMeshableField/LocalMeshabilityRepairT.hh>
#include <AlgoHex/LocallyMeshableField/LocalMeshabilityChecker.hh>
#include "../Stopwatches.hh"

namespace AlgoHex
{

template<class MeshT>
void
LocallyMeshableFieldGenerationT<MeshT>::generate_locally_meshable_field(MeshT &sg_mesh)
{
  //configure
  AlgoHex::SingularGraphExtractionT<MeshT> sge(mesh_);
  AlgoHex::MeshOptimizationT<MeshT> mo(mesh_);
  AlgoHex::FrameFieldOptimizationT<MeshT> ffo(mesh_);
  AlgoHex::LocalMeshabilityRepairT<MeshT> flnmv(mesh_);
  flnmv.set_first_stage_iters(first_iters_);
  flnmv.set_second_stage_inner_iters(second_inner_iteers_);

  flnmv.set_dual_path_weight_wrt_field_direction(dpath_weight_);
  //for filtering turning points
  flnmv.set_peak_coefficient(peak_coeff_);
  flnmv.set_segment_length_coefficient(seg_coeff_);
  flnmv.set_absolute_segment_length_coefficient(abs_seg_coeff_);
  flnmv.set_angle_threshold(qtn_angle_thre_);
  flnmv.set_unzip_zipper_nodes_from_longest_arc(fix_tps_on_longest_arc_first_);

  if (merge_zipper_node_)
    flnmv.set_merge_zipper_node();
  if (push_singular_circle_)
    flnmv.set_push_singular_circle();
  if (connect_other_zipper_node_)
  {
    flnmv.set_connect_other_zipper_node();
    flnmv.set_max_search_dist_coeff_to_other_zipper_node(max_tp_dist_);
  }

  if (truncated_newton_)
    mo.use_truncated_newton();

  ALGOHEX_DEBUG_ONLY(std::cerr << "###Parameters: first_iters: " << first_iters_ << " second_iters: " << second_iters_ << " second_inner_iters: "
            << second_inner_iteers_ << " field_alignment_weight: " << field_alignment_weight_
            << " tp_shrink_weight: " << tp_shrink_weight_ << " regularize_weight: " << regularize_weight_
            << " repulsion_weight: " << repulsion_weight_ << " ce_shrink_weight: " << ce_shrink_weight_
            << " truncated_newton: " << truncated_newton_ << " dpath_weight: " << dpath_weight_
            << " merge_zipper_node: " << merge_zipper_node_
            << " push_singular_circle " << push_singular_circle_ << " connect_other_zipper_node: "
            << connect_other_zipper_node_ << " max_tp_dist: " << max_tp_dist_ << " peak_coeff " << peak_coeff_
            << " seg_coeff: " << seg_coeff_ << " abs_seg_coeff: " << abs_seg_coeff_ << " qtn_angle_thre: "
            << qtn_angle_thre_ << std::endl;)

  //initialize properties
  AlgoHex::SplitHelperT<MeshT> pp(mesh_);
  pp.initialize_target_edge_length();
  pp.initialize_singular_vertex_property();

  //split
  pp.split_for_single_alignment_feature_face();

  if (!check_mesh_quality(mesh_))
    return;

  //align field
  ffo.align_quaternions_to_feature(100);

  //update matching and edge valence
  sge.get_transitions();
  sge.get_edges_valence();
  //update singular vertex property
  pp.initialize_singular_vertex_property();

  //save for animation
  mo.save_singular_graph_development(sg_mesh, 0);

  //preprocess
  //in coarse meshes, it might happen that many singular triangles connect to each other
  flnmv.preprocess();

  pp.split_special_edges_with_two_special_nodes();

  //save for animation
  mo.save_singular_graph_development(sg_mesh, 1);
  //save for debugging
//  MeshT axis_mesh;
  OpenVolumeMesh::IO::FileManager fm;
//  ffo.create_closest_field_axis(axis_mesh);
//  auto am_name = export_name_ + "axis" + std::to_string(0) + ".ovm";
  //fm.writeFile(am_name, axis_mesh);
  //mo.save_singular_graph(export_name_, 0, std::set<VH>{});

  //
  LocalMeshabilityChecker lmc(mesh_);

  int min_iters = first_iters_ + early_stop_iter_ * second_inner_iteers_;

  int i = 0, sg_id = i + 2;
  int n_fixable_tps = 0;

  int n_iters = first_iters_ + second_iters_ * second_inner_iteers_;
  {
    ScopedStopWatch sw(sw::lhfg_total);
    while (i < n_iters)
    {
      std::cout << "######Pipeline (Quaternion) iter " << i << " ..." << std::endl;

      int n_fixed = 0;
      //fix invalid local configuration
      {
        AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::singular_graph_repair);
        //fix at first stage
        n_fixed += flnmv.fix_first_stage();

        //only executed when enters second stage
        n_fixed += flnmv.fix_second_stage();
//                    mo.save_singular_graph(export_name_, sg_id-1);

        std::cerr << "####Fixed!" << std::endl;
      }

      //singularity relocation
      {
        AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::sg_relocation);
        mo.remesh(field_alignment_weight_, regularize_weight_, ce_shrink_weight_, tp_shrink_weight_,
                  repulsion_weight_, 2);

        {
          std::cerr << "####Splitting for DOF ..." << std::endl;
          AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::split_for_dof);
          pp.split_for_single_alignment(true);
          std::cerr << "####Splitting done!" << std::endl;
        }
      }

      //optimize quaternion field without alignments of singular edges
      ffo.optimize_quaternions_local(2, 50);

      std::cerr << "####ff opt done!" << std::endl;

      {
        AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::sgaf_visualization);
        //save singular edges for visualization
        mo.save_singular_graph_development(sg_mesh, sg_id++);
      }

      //check for termination
      if (i >= min_iters)
      {
        ScopedStopWatch ssw(sw::check_for_termination);

        int n_nm_vts, n_nm_nodes, n_tps, n_fixable_tps = std::numeric_limits<int>::max();
        flnmv.non_meshable_vertex_numbers(n_nm_vts, n_tps, n_nm_nodes);
        bool has_cpe = flnmv.has_invalid_singular_edge();
        std::set<VH> fxb_tps;
        if (n_nm_nodes == 0 && !has_cpe)
        {
          pp.split_for_single_alignment();
          fxb_tps = flnmv.get_fixable_zipper_nodes(true);
          n_fixable_tps = fxb_tps.size();
        }
        std::cerr << "#invalid vertices " << n_nm_vts << " invalid nodes(nsge<=4): " << n_nm_nodes
                  << " has complex edge " << has_cpe << " number of tps " << n_tps << " fixable "
                  << n_fixable_tps << " n_fixed " << n_fixed << std::endl;

        //mo.save_singular_graph(export_name_, sg_id-2, fxb_tps);
        //ffo.create_closest_field_axis(axis_mesh);
        //am_name = export_name_ + "axis" + std::to_string(sg_id-2) + ".ovm";
        //fm.writeFile(am_name, axis_mesh);

        //increase weight of shrinking complex singular edges
        if (i % 5 == 0 && i > 0)
        {
          ce_shrink_weight_ *= 2;

          if (n_fixable_tps == 0 && n_tps > 0)
            tp_shrink_weight_ *= 2;
        }

        //early stop
        if (n_nm_nodes == 0 && !has_cpe && n_fixable_tps == 0)
        {

          flnmv.uniform_singular_arc_valence();

          bool all_sv_meshable = lmc.check_special_vertices_local_meshability(true);

          flnmv.restore_singular_arc_valence();

          if (all_sv_meshable)
          {
            std::cerr << "All special vertices are locally meshable. Stop pipeline." << std::endl;
            i++;
            break;
          }
          else
            std::cerr << "Has non-meshable special vertices, continue..." << std::endl;
        }
      }
      else
      {
        // mo.save_singular_graph(export_name_, sg_id - 2, std::set<VH>{});
        // ffo.create_closest_field_axis(axis_mesh);
        // am_name = export_name_ + "axis" + std::to_string(sg_id-2) + ".ovm";
        // fm.writeFile(am_name, axis_mesh);
      }
      i++;
    }


    //uniform valence on singular arcs in case of pairs of turning points
    flnmv.uniform_singular_arc_valence();

    //avoid global inconsistency in the local meshability test
    {
      pp.split_for_single_alignment();
      std::cerr << "###Split for single alignments done!" << std::endl;
      flnmv.split_for_local_meshability_test();
      std::cerr << "###Split for local meshability test done!" << std::endl;
    }

    //optimize quaternion field with all constraints
    ffo.optimize_quaternions_global(100);

    {
      AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::sgaf_visualization);
      mo.save_singular_graph_development(sg_mesh, sg_id++);
    }
  }


  //post remeshing
  if (n_post_remeshing_ > 0)
  {
    AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::post_remesh);

    std::cerr << "###Mesh Quality before post-process: " << std::endl;
    check_mesh_quality(mesh_);
    mo.remesh_post(n_post_remeshing_);
    std::cerr << "###Mesh Quality after post-process: " << std::endl;
    check_mesh_quality(mesh_);

    //split
    pp.split_for_single_alignment();
    std::cerr << "###Split for single alignments done!" << std::endl;

    //only needed in lm test: avoid global inconsistency in the local meshability test
    flnmv.split_for_local_meshability_test();
    std::cerr << "###Split for local meshability test done!" << std::endl;


    mo.optimize_regular_vertices();
    std::cerr << "###Mesh Quality in the end: " << std::endl;
    check_mesh_quality(mesh_);

    ffo.optimize_quaternions_global(100);

    mo.save_singular_graph_development(sg_mesh, sg_id++);
  }
}

//For debugging
template<class MeshT>
void
LocallyMeshableFieldGenerationT<MeshT>::generate_locally_meshable_field(MeshT &sg_mesh, std::vector<MeshT> &meshes,
                                                                        MeshT &axis_mesh,
                                                                        std::vector<MeshT> &axis_meshes)
{
  int max_iter = first_iters_ + second_iters_ * second_inner_iteers_;
  meshes.clear();
  meshes.reserve(max_iter + 3);

  axis_meshes.clear();
  axis_meshes.reserve(max_iter + 3);

  AlgoHex::SingularGraphExtractionT<MeshT> sge(mesh_);
  AlgoHex::MeshOptimizationT<MeshT> mo(mesh_);

  AlgoHex::FrameFieldOptimizationT<MeshT> ffo(mesh_);
  AlgoHex::LocalMeshabilityRepairT<MeshT> flnmv(mesh_);
  flnmv.set_first_stage_iters(first_iters_);
  flnmv.set_second_stage_inner_iters(second_inner_iteers_);

  flnmv.set_dual_path_weight_wrt_field_direction(dpath_weight_);
  //for filtering turning points
  flnmv.set_peak_coefficient(peak_coeff_);
  flnmv.set_segment_length_coefficient(seg_coeff_);
  flnmv.set_absolute_segment_length_coefficient(abs_seg_coeff_);
  flnmv.set_angle_threshold(5.);
  flnmv.set_unzip_zipper_nodes_from_longest_arc(fix_tps_on_longest_arc_first_);

  AlgoHex::SplitHelperT<MeshT> pp(mesh_);
  pp.initialize_target_edge_length();
  pp.initialize_singular_vertex_property();

  //split
  pp.split_for_single_alignment_feature_face();

  if (!check_mesh_quality(mesh_))
    return;


  ffo.align_quaternions_to_feature(100);
  sge.get_transitions();
  sge.get_edges_valence();
  pp.initialize_singular_vertex_property();

  //save for animation
  meshes.push_back(mesh_);
  ffo.create_closest_field_axis(axis_mesh);
  axis_meshes.push_back(axis_mesh);

  mo.save_singular_graph_development(sg_mesh, 0);

  if (connect_other_zipper_node_)
  {
    flnmv.set_connect_other_zipper_node();
    flnmv.set_max_search_dist_coeff_to_other_zipper_node(max_tp_dist_);
  }
  if (push_singular_circle_)
    flnmv.set_push_singular_circle();
  if (merge_zipper_node_)
    flnmv.set_merge_zipper_node();

  //preprocess
  flnmv.preprocess();

  //save for animation
  meshes.push_back(mesh_);
  ffo.create_closest_field_axis(axis_mesh);
  axis_meshes.push_back(axis_mesh);

  mo.save_singular_graph_development(sg_mesh, 1);

  pp.split_special_edges_with_two_special_nodes();

  ScopedStopWatch sw(sw::lhfg_total);

  int n_fixed = 0;
  for (int i = 0; i < max_iter; ++i)
  {
    std::cout << "######Pipeline (Quaternion) iter " << i << " ..." << std::endl;
    if (i >= 0)
    {
      AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::singular_graph_repair);
      std::cerr << "####Fixing local invalid singular nodes at first stage..." << std::endl;
      //fix at first stage
      n_fixed += flnmv.fix_first_stage();

      n_fixed += flnmv.fix_second_stage();
      //save singular edges for visualization
      mo.save_singular_graph_development(sg_mesh, i + 2);

      //save meshes for visualization
      meshes.push_back(mesh_);

      //save axes
      ffo.create_closest_field_axis(axis_mesh);
      axis_meshes.push_back(axis_mesh);
    }

    {
      AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::sg_relocation);
      mo.remesh(field_alignment_weight_, regularize_weight_, ce_shrink_weight_, tp_shrink_weight_, repulsion_weight_, 2,
                meshes);

      {
        std::cerr << "####Splitting for DOF ..." << std::endl;

        AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::split_for_dof);
        pp.split_for_single_alignment(true);
      }
    }

    ffo.optimize_quaternions_local(2, 50);
  }

  flnmv.uniform_singular_arc_valence();
  meshes.push_back(mesh_);
  ffo.create_closest_field_axis(axis_mesh);
  axis_meshes.push_back(axis_mesh);


  {
    AlgoHex::ScopedStopWatch ssw(AlgoHex::sw::sg_relocation);
    AlgoHex::ScopedStopWatch ssw2(AlgoHex::sw::split_for_dof);

    std::cerr << "###Mesh Quality before post-process: " << std::endl;
    check_mesh_quality(mesh_);

    {
      mo.remesh_post(n_post_remeshing_);
      std::cerr << "###Mesh Quality after post-process: " << std::endl;
      check_mesh_quality(mesh_);

      mo.check_valence_change();
    }

    pp.split_for_single_alignment();
  }

  ffo.optimize_quaternions_global(qtn_angle_thre_);
}

template<class MeshT>
void
LocallyMeshableFieldGenerationT<MeshT>::print_timings() const
{
  std::cout << sw::lhfg_total << std::endl;
}
}

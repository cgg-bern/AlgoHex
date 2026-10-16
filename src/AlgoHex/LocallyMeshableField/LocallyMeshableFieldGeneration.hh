/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once
#include <AlgoHex/TypeDef.hh>


#include <vector>
#include <string>


namespace AlgoHex
{
class ALGOHEX_EXPORT LocallyMeshableFieldGeneration
{
public:
  LocallyMeshableFieldGeneration(TetMesh &_mesh) :
          mesh_(_mesh) {}

  void generate_locally_meshable_field(TetMesh &sg_mesh);

  void generate_locally_meshable_field(TetMesh &sg_mesh, std::vector<TetMesh> &meshes, TetMesh &axis_mesh,
                                       std::vector<TetMesh> &axis_meshes);

  void set_export_path(const std::string &filename) { export_name_ = filename; }

  void print_timings() const;

private:
  TetMesh &mesh_;

public:
  int first_iters_ = 12;
  int second_iters_ = 15;
  int second_inner_iteers_ = 4;
  double dpath_weight_ = 100.;

  //for zipper node filtering
  double peak_coeff_ = 0.1;
  double seg_coeff_ = 0.4;
  double abs_seg_coeff_ = 5.;

  //singularity relocation
  double field_alignment_weight_ = 1.;
  double ce_shrink_weight_ = 1.;
  double tp_shrink_weight_ = 0.1;
  double regularize_weight_ = 1.;
  double repulsion_weight_ = 0.1;

  bool edge_path_wrt_field_ = true;
  bool merge_zipper_node_ = true;
  bool push_singular_circle_ = false;
  int n_post_remeshing_ = 0;

  bool connect_other_zipper_node_ = false;
  bool fix_tps_on_longest_arc_first_ = true;
  double max_tp_dist_ = 1.1;

  double max_qts_angle_ = 180.;

  double qtn_angle_thre_ = 5.;

  bool truncated_newton_ = false;

  int early_stop_iter_ = second_iters_;

  std::string export_name_;
};

}

//=============================================================================
//=============================================================================

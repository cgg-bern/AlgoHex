/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once
#include <AlgoHex/TypeDef.hh>

#include <AlgoHex/TransitionQuaternionEigen.hh>

#include "MeshProperties.hh"


namespace OVM = OpenVolumeMesh;

namespace AlgoHex
{
class ALGOHEX_EXPORT FrameFieldOptimization : public MeshProperties
{
public:
  using Quaternion = Eigen::Quaterniond;
  using Point = TetMesh::PointT;

  using MeshProperties::valence_;
  using MeshProperties::trans_prop_;
  using MeshProperties::cell_quaternions_;

  using MeshProperties::sgl_vt_;
  using MeshProperties::feature_edge_;
  using MeshProperties::feature_fprop_;
  using MeshProperties::target_length_;

  using VecAxis = std::pair<Vec3d, int>;
  using DC = std::pair<double, CH>;

  FrameFieldOptimization(TetMesh &_mesh);

  void align_quaternions_to_feature(const int _iterations);

  void optimize_quaternions_local(const int k, const int _iterations);

  void optimize_quaternions_global(const double _angle_thr = 10.);

  void optimize_quaternions_global(const int _iterations);

  void create_closest_field_axis(TetMesh &_pm);

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

private:
  void check_feature_face_alignments_of_field();

  void check_edge_alignments_of_field();

  bool collect_alignments(CP<std::vector<VecAxis> > &_alignments);

  bool is_feature_cell(const CH _ch) const;

private:
  TetMesh &mesh_;
  AlgoHex::TransitionQuaternion tq_;
};
}

//=============================================================================
//=============================================================================

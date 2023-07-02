/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "QuaternionsSmoothing.hh"
#include "MeshProperties.hh"
#include "../Stopwatches.hh"
#include "../AxisAlignment.hh"


namespace OVM = OpenVolumeMesh;

namespace AlgoHex
{
template<class MeshT>
class FrameFieldOptimizationT : public MeshPropertiesT<MeshT>
{
public:
  using Quaternion = Eigen::Quaterniond;
  using Point = typename MeshT::PointT;

  using MeshPropertiesT<MeshT>::valence_;
  using MeshPropertiesT<MeshT>::trans_prop_;
  using MeshPropertiesT<MeshT>::cell_quaternions_;

  using MeshPropertiesT<MeshT>::sgl_vt_;
  using MeshPropertiesT<MeshT>::feature_edge_;
  using MeshPropertiesT<MeshT>::feature_fprop_;
  using MeshPropertiesT<MeshT>::target_length_;

  using VecAxis = std::pair<Vec3d, int>;
  using DC = std::pair<double, CH>;


  FrameFieldOptimizationT(MeshT &_mesh) : MeshPropertiesT<MeshT>(_mesh),
                                          mesh_(_mesh)
  {
  }


  void align_quaternions_to_feature(const int _iterations);

  void optimize_quaternions_local(const int k, const int _iterations);

  void optimize_quaternions_global(const double _angle_thr = 10.);

  void optimize_quaternions_global(const int _iterations);


  void create_closest_field_axis(MeshT &_pm);

  Eigen::Vector3d acg2eigen(const Point &_p) const { return Eigen::Vector3d(_p[0], _p[1], _p[2]); }

  int negate(const int _axis) const { return _axis % 2 == 0 ? _axis + 1 : _axis - 1; }

private:
  void check_feature_face_alignments_of_field();

  void check_edge_alignments_of_field();

  bool collect_alignments(CP<std::vector<VecAxis> > &_alignments);

  bool is_feature_cell(const CH _ch) const;

private:
  MeshT &mesh_;
  AlgoHex::TransitionQuaternion tq_;
};
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(FRAMEFIELDOPTIMIZATIONT_C)
#define FRAMEFIELDOPTIMIZATIONT_TEMPLATES

#include "FrameFieldOptimizationT_impl.hh"

#endif
//=============================================================================







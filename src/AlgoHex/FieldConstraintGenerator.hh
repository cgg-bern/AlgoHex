/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "TypeDef.hh"
#include "FaceNormalCache.hh"
#include "FieldConstraints.hh"

#include <Eigen/SVD>
#include <optional>

namespace AlgoHex
{

class ALGOHEX_EXPORT DihedralFeatureDetector
{
public:
  DihedralFeatureDetector(TetMesh &mesh, FaceNormalCache &normals_);

  std::vector<EH> compute(double angle_threshold_degrees);

private:
  TetMesh &mesh_;
  FaceNormalCache &normals_;
  //OVM::InterfaceAttrib interface_;
};

/// Try to average unoriented normals if they are not too different from one another.
inline std::optional<Point> normal_average(const std::vector<Point> &normals);

class ALGOHEX_EXPORT FieldConstraintGenerator
{
public:
  using Point = TetMesh::PointT;

  FieldConstraintGenerator(TetMesh &mesh, FaceNormalCache &normals_);

  void add_full_constraints_from_feature_edges(const std::vector<EH> &feature_edges);

  void add_partial_constraints_from_feature_edges(const std::vector<EH> &feature_edges);

  void add_alignment_constraints_from_feature_edges(const std::vector<EH> &feature_edges);

  void add_normal_constraints_from_feature_edges(const std::vector<EH> &feature_edges);

protected:
  std::vector<Point> collect_normals(VH vh, const std::vector<EH> &incident_feature_edges);

  Quaternion computeFullConstraint(const std::vector<Point> &normals);

  Quaternion computeFullConstraint(const std::vector<Point> &normals, const VH &vh, const std::vector<EH> &edges);

private:
  TetMesh &mesh_;
  FaceNormalCache &normals_;
  VFullConstraints &vertex_full_constraints_;
  VNormalConstraints &vertex_normal_constraints_;
};


class ALGOHEX_EXPORT FieldConstraintGeneratorCellBased
{
public:
  using Point = TetMesh::PointT;

  FieldConstraintGeneratorCellBased(TetMesh &mesh, FaceNormalCache &normals_);

  void add_partial_constraints(const std::vector<EH> &feature_edges);

  void add_feature_constraints(const std::vector<EH> &feature_edges);

  void add_normal_constraints();

private:
  TetMesh &mesh_;
  FaceNormalCache &normals_;
  CPartialConstraints &cell_partial_constraints_;
};

}

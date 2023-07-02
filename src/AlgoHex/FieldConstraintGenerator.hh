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

template<typename MeshT>
class DihedralFeatureDetector
{
public:
  DihedralFeatureDetector(MeshT &mesh, FaceNormalCache<MeshT> &normals_)
          : mesh_(mesh),
            normals_(normals_)
  //interface_(mesh),
  {
  }

  std::vector<EH> compute(double angle_threshold_degrees)
  {
    // TODO: handle interior interfaces too; this is a bit tricky:
    //       Each interface face has *two* interface halffaces, we
    //       need to capture the right pair (so the angle threshold makes sense)
    //       We also need to handle non-manifold configurations, possibly
    //       involving the mesh surface.
    std::vector<EH> feature_edges;

    double cos_threshold = std::cos(angle_threshold_degrees / 180. * M_PI);

    /// if there is exactly one boundary halfface adjacent to `heh`, return its handle,
    /// otherwise InvalidHalfFaceHandle.
    auto one_incident_boundary_halfface = [this](const HEH heh) -> HFH {
      HFH res = MeshT::InvalidHalfFaceHandle;
      for (const auto hfh: mesh_.halfedge_halffaces(heh))
      {
        if (!mesh_.is_boundary(hfh))
          continue;
        if (res.is_valid())
        { // more than one
          return MeshT::InvalidHalfFaceHandle;
        }
        res = hfh;
      }
      return res;
    };

    for (const auto &eh: mesh_.edges())
    {
      if (!mesh_.is_boundary(eh))
        continue;
      HFH hfh0 = one_incident_boundary_halfface(mesh_.halfedge_handle(eh, 0));
      if (!hfh0.is_valid()) continue;
      HFH hfh1 = one_incident_boundary_halfface(mesh_.halfedge_handle(eh, 1));
      if (!hfh1.is_valid()) continue;

      auto n0 = normals_[hfh0];
      auto n1 = normals_[hfh1];

      if (n0.dot(n1) < cos_threshold)
      {
        feature_edges.push_back(eh);
      }
    }
    return feature_edges;
  }

private:
  MeshT &mesh_;
  FaceNormalCache<MeshT> &normals_;
  //OVM::InterfaceAttrib interface_;
};

/// Try to average unoriented normals if they are not too different from one another.
template<typename Point>
std::optional<Point> normal_average(const std::vector<Point> &normals)
{
  static const double cos_threshold = std::cos(70 / 180. * M_PI);
  // Arbitrarily take one reference normal, if any other deviates by
  // an angle larger than a threshold, do not fix this vertex normal.
  Point ref = normals[0];
  Point sum{0., 0., 0.};

  for (const auto n: normals)
  {
    double dotprod = n.dot(ref);
    if (std::fabs(dotprod) < cos_threshold)
    {
      return {};
    }
    if (dotprod < 0)
    {
      sum -= n;
    }
    else
    {
      sum += n;
    }
  }
  return sum.normalized();
}

template<typename MeshT>
class FieldConstraintGenerator
{
public:
  using Point = typename MeshT::PointT;

  FieldConstraintGenerator(MeshT &mesh, FaceNormalCache<MeshT> &normals_)
          : mesh_(mesh),
            normals_(normals_),
            vertex_full_constraints_(get_vertex_full_constraints(mesh_)),
            vertex_normal_constraints_(get_vertex_normal_constraints(mesh_))
  {
    vertex_full_constraints_.clear();
    vertex_normal_constraints_.clear();
  }


  void add_full_constraints_from_feature_edges(const std::vector<EH> &feature_edges)
  {

    std::map<VH, std::vector<EH>> incident_feature_edges;
    for (const auto &eh: feature_edges)
    {
      auto edge = mesh_.edge(eh);
      incident_feature_edges[edge.from_vertex()].push_back(eh);
      incident_feature_edges[edge.to_vertex()].push_back(eh);
    }
    for (const auto&[vh, edges]: incident_feature_edges)
    {
      // we can only handle up to 3 normals, do not constrain vertices that are, e.g.,
      // a pyramid apex.

      if (edges.size() > 3 || edges.size() <= 1)
        continue;

      auto normals = collect_normals(vh, edges);
      auto q = computeFullConstraint(normals);

      vertex_full_constraints_.emplace_back(vh, q);
    }

  }

  void add_partial_constraints_from_feature_edges(const std::vector<EH> &feature_edges)
  {
    add_normal_constraints_from_feature_edges(feature_edges);
    add_alignment_constraints_from_feature_edges(feature_edges);
  }

  void add_alignment_constraints_from_feature_edges(const std::vector<EH> &feature_edges)
  {
    auto feature_node = mesh_.template request_vertex_property<int>("AlgoHex::FeatureVertices");

    std::map<VH, std::vector<EH>> incident_feature_edges;
    for (const auto &eh: feature_edges)
    {
      auto edge = mesh_.edge(eh);
      auto vhf = edge.from_vertex();
      if (feature_node[vhf] == 0)
        incident_feature_edges[vhf].push_back(eh);

      auto vht = edge.to_vertex();
      if (feature_node[vht] == 0)
        incident_feature_edges[vht].push_back(eh);
    }

    for (const auto&[vh, edges]: incident_feature_edges)
    {
      if (edges.size() > 2 || edges.empty())
        continue;

      Point feature_dir(0.);
      if (edges.size() == 2)
      {
        VH e0_from_vh =
                mesh_.edge(edges[0]).from_vertex() == vh ? mesh_.edge(edges[0]).to_vertex() : mesh_.edge(
                        edges[0]).from_vertex();
        VH e1_to_vh = mesh_.edge(edges[1]).to_vertex() == vh ? mesh_.edge(edges[1]).from_vertex() : mesh_.edge(
                edges[1]).to_vertex();

        auto e_dir0 = mesh_.vertex(vh) - mesh_.vertex(e0_from_vh);
        e_dir0.normalize();
        auto e_dir1 = mesh_.vertex(e1_to_vh) - mesh_.vertex(vh);
        e_dir1.normalize();


        feature_dir = e_dir0 + e_dir1;
      }
      else
        feature_dir = mesh_.vertex(mesh_.edge(edges[0]).to_vertex()) - mesh_.vertex(mesh_.edge(edges[0]).from_vertex());

      feature_dir.normalize();

      vertex_normal_constraints_.emplace_back(vh, feature_dir);
    }

  }

  void add_normal_constraints_from_feature_edges(const std::vector<EH> &feature_edges)
  {
    auto feature_fprop = mesh_.template request_face_property<int>("AlgoHex::FeatureFaces");

    auto incident_to_feature_edge = mesh_.template request_vertex_property<bool>("", false);

    for (const auto &eh: feature_edges)
    {
      auto edge = mesh_.edge(eh);
      incident_to_feature_edge[edge.from_vertex()] = true;
      incident_to_feature_edge[edge.to_vertex()] = true;
    }

    std::map<VH, std::vector<Point>> vertex_normals;

    for (const auto &fh: mesh_.faces())
    {
      if (feature_fprop[fh] == 0)
        continue;

      HFH hfh = mesh_.halfface_handle(fh, 0);
      if (mesh_.is_boundary(fh) && !mesh_.is_boundary(hfh))
      {
        hfh = mesh_.halfface_handle(fh, 1);
      }

      const auto &normal = normals_[hfh];

      for (const auto &vh: mesh_.face_vertices(fh))
      {
        if (incident_to_feature_edge[vh])
          continue;
        vertex_normals[vh].push_back(normal);
      }
    }

    vertex_normal_constraints_.reserve(vertex_normals.size());
    for (const auto &[vh, normals]: vertex_normals)
    {
#if 0
      if (mesh_.is_boundary(vh)) {
          Point sum = std::accumulate(std::next(normals.begin()), normals.end(), normals[0]);
          vertex_normal_constraints_.emplace_back(vh, sum.normalized());
          return;
      }
#endif

      // as we are not handling interior features yet, this may be a feature vertex,
      // in that case we cannot just average the normals.

      auto na = normal_average(normals);
      if (na.has_value())
      {
        vertex_normal_constraints_.emplace_back(vh, *na);
      }


    }


  }

protected:
  std::vector<Point> collect_normals(VH vh, const std::vector<EH> &incident_feature_edges)
  {
    const auto &edges = incident_feature_edges;
    // walk through incident boundary faces and average
    // normals for each feature-edge partitioned segment
    std::vector<Point> normals;
    std::set<HFH> visited;
    HEH cur_heh = mesh_.halfedge_handle(edges[0], 0);
    if (mesh_.from_vertex_handle(cur_heh) != vh)
    {
      cur_heh = mesh_.halfedge_handle(edges[0], 1);
    }
    HFH cur_hfh = MeshT::InvalidHalfFaceHandle;
    // find first boundary hfh:
    std::vector<Point> segment_normals;
    int n_voh = mesh_.valence(vh);
    // TODO: prevent endless loop by testing against the max possible number of edges visited
    while (true)
    {
      for (const auto hfh: mesh_.halfedge_halffaces(cur_heh))
      {
        if (mesh_.is_boundary(hfh) && visited.find(hfh) == visited.end())
        {
          cur_hfh = hfh;
          visited.insert(cur_hfh);
          break;
        }
      }
      segment_normals.push_back(normals_[cur_hfh]);

      HEH next_heh = MeshT::InvalidHalfEdgeHandle;
      for (const auto heh: mesh_.halfface_halfedges(cur_hfh))
      {
        if (heh == cur_heh) continue;
        auto halfedge = mesh_.halfedge(heh);
        if (halfedge.to_vertex() == vh)
        {
          next_heh = mesh_.opposite_halfedge_handle(heh);
          break;
        }
      }
      assert(next_heh.is_valid());
      cur_heh = next_heh;
      if (std::find(edges.begin(), edges.end(), mesh_.edge_handle(cur_heh)) != edges.end())
      {
        // end of sector.
        // TODO: if normals deviate too much, do not constrain this vertex
        // TODO: angle-weighted normals?
        Point nsum = std::accumulate(segment_normals.begin(),
                                     segment_normals.end(),
                                     Point(0., 0., 0.));
        normals.push_back(nsum.normalized());
        segment_normals.clear();
      }
      if (mesh_.edge_handle(cur_heh) == edges[0])
      {
        // roundtrip complete.
        break;
      }

    }
    return normals;
  }

  Quaternion computeFullConstraint(const std::vector<Point> &normals)
  {
    assert(normals.size() == 2 || normals.size() == 3);

    Point third = (normals.size() == 3) ? normals[2] : normals[0].cross(normals[1]).normalized();

    Eigen::Matrix3d mat;
    mat << normals[0][0], normals[1][0], third[0],
            normals[0][1], normals[1][1], third[1],
            normals[0][2], normals[1][2], third[2];

    // compute nearest frame for normals; SVD -> UV^T
    Eigen::JacobiSVD svd{mat, Eigen::ComputeFullU | Eigen::ComputeFullV};
    Eigen::Matrix3d ortho = svd.matrixU() * svd.matrixV().transpose();

    // flip orientation if needed (might be a reflection instead of a rotation):
    ortho.col(2) = ortho.col(0).cross(ortho.col(1));
    return Quaternion{ortho};
  }


  Quaternion computeFullConstraint(const std::vector<Point> &normals, const VH &vh, const std::vector<EH> &edges)
  {
    assert(normals.size() == 2);
    VH e0_from_vh = mesh_.edge(edges[0]).from_vertex() == vh ? mesh_.edge(edges[0]).to_vertex() : mesh_.edge(
            edges[0]).from_vertex();
    VH e1_to_vh = mesh_.edge(edges[1]).to_vertex() == vh ? mesh_.edge(edges[1]).from_vertex() : mesh_.edge(
            edges[1]).to_vertex();

    auto e_dir0 = mesh_.vertex(vh) - mesh_.vertex(e0_from_vh);
    e_dir0.normalize();
    auto e_dir1 = mesh_.vertex(e1_to_vh) - mesh_.vertex(vh);
    e_dir1.normalize();


    auto feature_dir = e_dir0 + e_dir1;
    feature_dir.normalize();

    auto normal_dir = normals[0] + normals[1];
    normal_dir.normalize();

    auto third = normals[0].cross(normals[1]);

    Eigen::Matrix3d mat;
    mat << feature_dir[0], normal_dir[0], third[0],
            feature_dir[1], normal_dir[1], third[1],
            feature_dir[2], normal_dir[2], third[2];

    // compute nearest frame for normals; SVD -> UV^T
    Eigen::JacobiSVD svd{mat, Eigen::ComputeFullU | Eigen::ComputeFullV};
    Eigen::Matrix3d ortho = svd.matrixU() * svd.matrixV().transpose();

    // flip orientation if needed (might be a reflection instead of a rotation):
    ortho.col(2) = ortho.col(0).cross(ortho.col(1));
    return Quaternion{ortho};
  }


private:
  MeshT &mesh_;
  FaceNormalCache<MeshT> &normals_;
  VFullConstraints &vertex_full_constraints_;
  VNormalConstraints<MeshT> &vertex_normal_constraints_;

};


template<typename MeshT>
class FieldConstraintGeneratorCellBased
{
public:
  using Point = typename MeshT::PointT;

  FieldConstraintGeneratorCellBased(MeshT &mesh, FaceNormalCache<MeshT> &normals_)
          : mesh_(mesh),
            normals_(normals_),
            cell_partial_constraints_(get_cell_partial_constraints(mesh_))
  {
    cell_partial_constraints_.clear();
  }

  void add_partial_constraints(const std::vector<EH> &feature_edges)
  {
    add_normal_constraints();
    add_feature_constraints(feature_edges);
  }

  void add_feature_constraints(const std::vector<EH> &feature_edges)
  {
    std::map<CH, std::vector<EH>> incident_feature_edges;
    for (const auto &eh: feature_edges)
    {
      for (auto ec_it = mesh_.ec_iter(eh); ec_it.valid(); ++ec_it)
        incident_feature_edges[*ec_it].push_back(eh);
    }

    for (const auto&[ch, edges]: incident_feature_edges)
    {
      if (edges.size() > 1 || edges.empty())
        continue;

      Point feature_dir =
              mesh_.vertex(mesh_.edge(edges[0]).to_vertex()) - mesh_.vertex(mesh_.edge(edges[0]).from_vertex());

      feature_dir.normalize();

      cell_partial_constraints_.emplace_back(ch, feature_dir);
    }

  }

  void add_normal_constraints()
  {
//            OVM::InterfaceAttrib interface {mesh_};
    auto feature_fprop = mesh_.template request_face_property<int>("AlgoHex::FeatureFaces");
    std::map<CH, std::vector<Point>> cell_normals;

    for (const auto &fh: mesh_.faces())
    {
      if (feature_fprop[fh] == 0)
        continue;

      HFH hfh = mesh_.halfface_handle(fh, 0);
      CH ch_inc = mesh_.incident_cell(hfh);
      if (ch_inc.is_valid())
      {
        const auto &normal = normals_[hfh];
        cell_normals[mesh_.incident_cell(hfh)].push_back(normal);
      }

      HFH hfh_opp = mesh_.halfface_handle(fh, 1);
      ch_inc = mesh_.incident_cell(hfh_opp);
      if (ch_inc.is_valid())
      {
        const auto &normal = normals_[hfh_opp];
        cell_normals[mesh_.incident_cell(hfh_opp)].push_back(normal);
      }
    }

    cell_partial_constraints_.reserve(cell_normals.size());
    for (const auto &[ch, normals]: cell_normals)
    {
      // One cell should has one boundary face. But if there are more,
      // we average them.

      auto na = normal_average(normals);
      if (na.has_value())
      {
        cell_partial_constraints_.emplace_back(ch, *na);
      }
    }
  }


private:
  MeshT &mesh_;
  FaceNormalCache<MeshT> &normals_;
  CPartialConstraints<MeshT> &cell_partial_constraints_;

};

}

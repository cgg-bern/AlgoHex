/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

//== INCLUDES =================================================================

#include <queue>
#include <vector>
#include "TransitionQuaternionEigen.hh"
#include "TypeDef.hh"
#include "SphericalHarmonics/SHCoeffs.hh"
#include "SphericalHarmonics/SHProjectorRay.hh"


#include <AlgoHex/Stopwatches.hh>


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================
namespace OVM = OpenVolumeMesh;

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================




/** \class SingularGraphExtractionT SingularGraphExtractionT.hh

    Brief Description.

    A more elaborate description follows.
*/

template<class MeshT>
class SingularGraphExtractionT
{
public:
  using Mesh = MeshT;
  using Point = typename Mesh::PointT;

  SingularGraphExtractionT(MeshT &_mesh) :
          mesh_(_mesh),
          shc_(mesh_.template request_vertex_property<SHCoeffs>("scaled vertex spherical harmonic coefficients")),
          quaternion_(mesh_.template request_vertex_property<Quaternion>("vertex quaternion", {0., 0., 0., 0.})),
          cell_quaternion_(mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions", {0., 0., 0., 0.})),
          trans_prop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
          valence_(mesh_.template request_edge_property<int>("edge_valance")),
          feature_edge_prop_(mesh_.template request_edge_property<int>("AlgoHex::FeatureEdges")),
          feature_edge_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureEdgeVertices")),
          feature_fprop_(mesh_.template request_face_property<int>("AlgoHex::FeatureFaces")),
          feature_face_vertex_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureFaceVertices")),
          feature_face_edge_(mesh_.template request_edge_property<bool>("AlgoHex::FeatureFaceEdges"))
  {

    //calculate average edge length
    av_edge_length_ = 0.;
    for (auto ehi: mesh_.edges())
    {
      av_edge_length_ += mesh_.length(ehi);
    }
    av_edge_length_ /= mesh_.n_edges();

    initial_num_vts_ = mesh_.n_vertices();

    mesh_.set_persistent(cell_quaternion_, true);
    mesh_.set_persistent(valence_, true);
    mesh_.set_persistent(trans_prop_, true);
    projector_.set_seed_selection(SHProjectorRay::BestSeed);
  }

  ~SingularGraphExtractionT()
  {
    mesh_.set_persistent(shc_, false);
    mesh_.set_persistent(quaternion_, false);
  }

public:
  void get_singular_edges_with_refinement(const double _largest_length_ratio);

  void save_cell_quaternions(const std::string &_filename) const;

  void get_singular_edges();

  void get_edges_valence();

  void compute_edge_valence(const EH _eh);

  int calculate_edge_valence(const EH _eh);

  void get_transitions();

  void compute_face_transition(const FH _fh);

  int compute_edge_transition(const HEH _heh, bool _is_boundary) const;

  void refine_singular_regions(const double _largest_length_ratio);

  double unsigned_rotation_angle_at_edge(const HEH _heh) const;

  double frame_smoothness_energy_at_face(const FH _fh) const;

  int get_dominant_axis_in_cell(const HEH _heh) const;

  void load_cell_quaternions(const std::string &_filename);

private:
  void refine_interior_regions(const double _edge_length);

  int compute_tri_idx(const Quaternion &q0, const Quaternion &q1, const Quaternion &q2) const;

  double longest_edge_length(const FH _fh, EH &_eh) const;

  VH split_edge(const EH _eh);

  VH split_face(const FH _fh, const std::vector<VH> &_fvhs);

  void get_field_per_cell();

  double find_best_transition(const HFH hfh, int &_best_transition) const;

  int compute_valence(const HEH _heh, const Quaternion &_q, int _idx, const Point &_vn);

  double calc_sector_angle(const HEH _heh, const HFH _hfh_sector_start) const;

private:
  Mesh &mesh_;

  //spherical harmonics coefficients per vertex
  VP<SHCoeffs> shc_;

  //quaternion per vertex
  VP<Quaternion> quaternion_;

  //cell quaternions
  CP<Quaternion> cell_quaternion_;

  //TODO: cache matrix
//        //frame field
//        CP<Mat3d> frame_cprop_;

  //matching property
  HFP<int> trans_prop_;

  //valence property
  EP<int> valence_;

  //feature edge property
  EP<int> feature_edge_prop_;
  VP<bool> feature_edge_vertex_;

  //feature face property
  FP<int> feature_fprop_;
  VP<bool> feature_face_vertex_;
  EP<bool> feature_face_edge_;

  SHProjectorRay projector_;

  AlgoHex::TransitionQuaternion tq_;

  double av_edge_length_;

  int initial_num_vts_;
};


//=============================================================================
} // namespace AlgoHex
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(SINGULARGRAPHEXTRACTIONT_C)
#define SINGULARGRAPHEXTRACTION_TEMPLATES

#include "SingularGraphExtractionT_impl.hh"

#endif
//=============================================================================

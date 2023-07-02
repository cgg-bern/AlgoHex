/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

//== INCLUDES =================================================================

#include <map>
#include <vector>
#include <iostream>

#include <AlgoHex/Util/json.hh>
#include <AlgoHex/Stopwatches.hh>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <CoMISo/Config/config.hh>
#include <CoMISo/NSolver/COMISOSolver.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>
#include <CoMISo/NSolver/FiniteElementProblem.hh>

#include "TypeDef.hh"
#include "TransitionQuaternionEigen.hh"
#include "AxisAlignment.hh"
#include "LogDetElement3D.hh"
#include "AMIPSFrameElement3DTinyAD.hh"
#include "AlignmentElements3D.hh"
#include "DeformationElements3D.hh"
#include "IntegrabilityElements3D.hh"
#include "LinearLeastSquaresElements.hh"

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================

/** \class FrameFieldOptimizer3DT FrameFieldOptimizer3DT.hh

    Brief Description.
  
    A more elaborate description follows.
*/

template<class TetMeshT>
class FrameFieldOptimizer3DT
{
public:
  // typedefs for shorter access
  typedef TetMeshT TetMesh;
  typedef typename TetMesh::PointT Point;
  typedef typename TetMesh::Cell Cell;
  typedef typename TetMesh::Face Face;

  struct CellAlignmentInfo
  {
    HFH hfh = HFH(-1);
    int hfh_axis = -1;
    HEH heh = HEH(-1);
    int heh_axis0 = -1;
    int heh_axis1 = -1;

    inline void set_edge_alignment(const HEH _heh, const int _heh_axis0, const int _heh_axis1)
    {
      heh = _heh;
      heh_axis0 = _heh_axis0;
      heh_axis1 = _heh_axis1;
    }

    inline void set_face_alignment(const HFH _hfh, const int _hfh_axis)
    {
      hfh = _hfh;
      hfh_axis = _hfh_axis;
    }
  };

  // constructors
  FrameFieldOptimizer3DT(TetMesh &_m
//          , OpenVolumeMesh::StatusAttrib &_mesh_status
          /*OpenVolumeMesh::ColorAttrib<ACG::Vec4f>& _mesh_color*/)
          : mesh_(_m),
//            mesh_status_(_mesh_status),
//    mesh_color_(_mesh_color),
            frame_cprop_(mesh_.template request_cell_property<Mat3d>("AlgoHex::FrameField")),
            transition_hfprop_(mesh_.template request_halfface_property<int>("HalffaceTransiton")),
            valence_eprop_(mesh_.template request_edge_property<int>("edge_valance")),
            feature_vprop_(mesh_.template request_vertex_property<int>("AlgoHex::FeatureVertices")),
            feature_eprop_(mesh_.template request_edge_property<int>("AlgoHex::FeatureEdges")),
            feature_fprop_(mesh_.template request_face_property<int>("AlgoHex::FeatureFaces")),
            save_locally_non_meshable_(false),
            save_locally_non_meshable_filename_base_("locally_non_meshable_")
  {
    mesh_.set_persistent(valence_eprop_);
    mesh_.set_persistent(transition_hfprop_);
  }

//  optimize frame field to become more integrable
//  input:
//    frame_cprop_( mesh_.template request_cell_property<<Mat3x3>("AlgoHex::FrameField"))
//    transition_hfprop_( mesh_.template request_halfface_property<int>("HalffaceTransiton")),
//    valence_eprop_( mesh_.template request_edge_property<int>("edge_valance")),
//    feature_vprop_(mesh_.template request_vertex_property<bool>("AlgoHex::FeatureVertices", false)),
//    feature_eprop_(mesh_.template request_edge_property<bool>("feature edges", false)),
//    feature_fprop_(mesh_.template request_face_property<bool>("AlgoHex::FeatureFaces", false))
//    input requirements:
//        (1) field is aligned to boundary normals, singular edges and feature curves
//        (2) each tet has at most one such alignment constraint (tet mesh is split properly) except for boundary normal + boundary feature
//  output:
//    frame_cprop_( mesh_.template request_cell_property<<Mat3x3>("AlgoHex::FrameField"))

  void optimize_integrable(const double _anisotropy_alpha = 1.0);

  void optimize_integrable_regularized(const double _anisotropy_alpha = 1.0, const double _w_frame_smoothness = 1e4,
                                       const double _w_integrability = 1e5);

  double check_local_meshability();

  nlohmann::json &json_data() { return json_data_; }

  void enable_save_locally_non_meshable(const std::string _filename_base = "locally_non_meshable_")
  {
    save_locally_non_meshable_ = true;
    save_locally_non_meshable_filename_base_ = _filename_base;
  }

  void import_frames_from_vec3d_properties();

  void import_frames_from_quaternion_property_onering(const VH &_vh);

  void import_quaternions_from_frames();

  void closest_rotation(const Mat3d &_A, Mat3d &_R) const;

  void check_valence_consistency();

  void check_frame_rotation_angles();

private:

  void mesh_tet(const CH _ch, Mat3x4d &_P) const;

  double volume(const CH _ch) const;

  double mesh_volume() const;

  double volume_frame(const CH _ch) const;

//  double interior_angle(const EH _eh) const;

  bool is_on_feature_vertex(const VH _vh) const;

  bool is_on_feature_edge(const VH _vh) const;

  bool is_on_feature_face(const VH _vh) const;

  bool is_on_singular_node(const VH _vh) const;

  bool is_on_singular_arc(const VH _vh) const;

  bool has_non_feature_boundary_face(const EH _eh) const;

  bool check_sector_valencies(const EH _eh);

  HFH find_feature_halfface(const HEH _heh);

  double calc_sector_angle(const HEH _heh, const HFH _hfh_sector_start, HFH &_hfh_sector_end);

  int incident_feature_ohalfedges(const VH _vh, std::vector<HEH> &_fheh) const;

  int incident_feature_faces(const VH _vh, std::vector<FH> &_ffh) const;

  void restore_field_alignment(const std::map<CH, CellAlignmentInfo> &_cell_alignment_info);

public:
  bool is_locally_meshable(const VH _vh, const bool _verbose = false) const;

  bool is_locally_meshable(const VH _vh, TetMesh &_export_tmesh, const bool _verbose = false) const;

private:
  void export_local_configuration(TetMesh &tetmesh, const std::vector<CH> &_chv, const std::vector<EH> &_invalid_edges,
                                  const std::vector<double> &_x0, const std::vector<double> &_x) const;

private:

  // mesh, status and color attribs
  TetMesh &mesh_;
//  OpenVolumeMesh::StatusAttrib &mesh_status_;
//  OpenVolumeMesh::ColorAttrib<ACG::Vec4f>& mesh_color_;

  // ##### input properties
  // frame field input/output
  OpenVolumeMesh::CellPropertyT<Mat3d> frame_cprop_;
  // halfface transitions (for frames -> F0*T = F1)
  OpenVolumeMesh::HalfFacePropertyT<int> transition_hfprop_;
  // integer index of edges (nonzero for singularities)
  OpenVolumeMesh::EdgePropertyT<int> valence_eprop_;
  // feature vertices/edges/faces
  OpenVolumeMesh::VertexPropertyT<int> feature_vprop_;
  OpenVolumeMesh::EdgePropertyT<int> feature_eprop_;
  OpenVolumeMesh::FacePropertyT<int> feature_fprop_;

  // transiton function handling
  TransitionQuaternion tq_;

  bool save_locally_non_meshable_;
  std::string save_locally_non_meshable_filename_base_;

  // json storage
  nlohmann::json json_data_;

  const bool always_export_local_configuration_ = false;
};



//=============================================================================
} // namespace AlgoHex
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(FRAMEFIELDOPTIMIZER3DT_C)
#define FRAMEFIELDOPTIMIZER3DT_TEMPLATES

#include "FrameFieldOptimizer3DT_impl.hh"

#endif


/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define FRAMEFIELDOPTIMIZER3DT_C


//== INCLUDES =================================================================

//#include <ACG/Utils/HaltonColors.hh>

#include <sstream>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include "FrameFieldOptimizer3DT.hh"

#include <Base/Debug/DebConfig.hh>
#include <CoMISo/NSolver/NPDerivativeChecker.hh>
#include <CoMISo/NSolver/NPTiming.hh>
#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/IPOPTSolver.hh>
#include <CoMISo/NSolver/TruncatedNewtonPCG.hh>
#include <CoMISo/NSolver/TrustregionNewtonPCG.hh>
#include <CoMISo/NSolver/AugmentedLagrangianMethod.hh>
#include <CoMISo/NSolver/ConstraintTools.hh>
#include <CoMISo/Utils/Tools.hh>
#include <AlgoHex/Util/StopWatch.hh>

#include <Eigen/Geometry>

#include "Geometry.hh"
#include "LinearLeastSquaresElements.hh"
#include "AMIPSFrameElement3D.hh"
#include "DihedralAngleElementTinyAD.hh"

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== IMPLEMENTATION ==========================================================

template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::
optimize_integrable(const double _anisotropy_alpha)
{
  ScopedStopWatch sw_total(sw::frame_field_int_opt);

  std::cout << "#####Optimize Frame Field..." << std::endl;

  // check whether input is valid
  check_valence_consistency();
  check_frame_rotation_angles();

  const bool amips_active = false;
  const bool log_det_active = false;
  const bool frame_fitting_active = false;
  const bool frame_smoothness_active = false;
  const bool align_penalty_active = false;
  const bool integrability_penalty_active = true;
  const bool dihedral_active = false;
  const bool symmetric_dirichlet_active = true;

  const double w_amips = 1.0;
  const double w_log_det = 1.0;
  const double w_frame_fitting = 1.0;
  const double w_frame_smoothness = 1e3;
  const double w_align = 1.0;
  const double w_integrability = 1e7;
  const double alpha_fitting = 1.0;
  const double w_dihedral = 1.0;
  const double s_dihedral = 4.0;
  const double w_symmetric_dirichlet = 1.0;
//  const double w_optrot   = 0.0;

  // variables are 3x3 frames per cell
  int nv = 9 * mesh_.n_cells();

  if (nv == 0)
  {
    std::cerr << "skip optimize_integrable since mesh does not have cells\n";
    return;
  }

  double VM = mesh_volume();
  double VM_cbrt = std::cbrt(VM);

  // create problem
  COMISO::FiniteElementProblem fe_problem(nv);
  COMISO::FiniteElementProblem fe_problem_complete(nv);

  // set gradient field as initial solution (gradient field is inverse-transpose of frame field)
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    int idx = 9 * c_it->idx();
    Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[idx]));
    G = (frame_cprop_[*c_it].inverse()).transpose();
  }

  // objective terms
  COMISO::FiniteElementSet <AMIPSFrameElement3D_PH> fe_amips("AMIPS");
  COMISO::FiniteElementSet <AMIPSFrameElement3D_TAD_PH> fe_amips_tad_ph("AMIPS_TinyAD");
  COMISO::FiniteElementSet <LogDetElement3D_PH> fe_logdet("LogDet");
  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_frame_fit("FrameFit");
  COMISO::FiniteElementSet <DihedralAngleFramesElement_TAD_PH> fe_dihedral("DihedralAngle");
  COMISO::FiniteElementSet <SymmetricDirichletFrameElement3D_TAD_PH> fe_sdirichlet("SymmetricDirichlet");

  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_face_align("FaceAlignment");
  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_edge_align("EdgeAlignment");
  COMISO::FiniteElementSet <LinearLeastSquaresElement6D> fe_integrability("Integrability");
  COMISO::FiniteElementSet <LinearLeastSquaresElement2D> fe_frame_smoothness("Frame Smoothness");

//  COMISO::FiniteElementSet<OptimalRotationFrameElement3D_TAD_PH> fe_optrot;

  // constraints
  std::vector<COMISO::LinearConstraint> c_face_align;
  std::vector<COMISO::LinearConstraint> c_edge_align;
  std::vector<COMISO::LinearConstraint> c_integrability;

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // calc volume of tet
    double V = std::abs(volume(*c_it));
    double F_det = frame_cprop_[*c_it].determinant();
    double as = std::cbrt(F_det);

    const double offset = 0.001;
//    const double scale = F_det * (1.0+offset);
    const double scale = F_det;

    // log determinant terms
    LogDetElement3D_PH::VecI idxs0;
    LogDetElement3D_PH::VecC c0;
    int i = 9 * c_it->idx();
    idxs0 << i, i + 1, i + 2, i + 3, i + 4, i + 5, i + 6, i + 7, i + 8;
    c0 << w_log_det, scale, offset; // weight, scale, offset
    fe_logdet.instances().add_element(idxs0, c0);

    // frame fitting terms
    std::vector<LinearLeastSquaresElement3D::VecI> idxs1;
    std::vector<LinearLeastSquaresElement3D::VecC> c1;

    double w1 = (w_frame_fitting * V / VM) * (as * as);
    // get 9 instances of LinearElement
    LinearLeastSquaresElement3D::get_constants_and_indices_frame_fitting(
            Eigen::Map<Mat3d>((double *) &(fe_problem.x()[i])),
            c_it->idx(),
            w1,
            alpha_fitting,
            idxs1,
            c1);
    for (size_t j = 0; j < idxs1.size(); ++j)
    {
      fe_frame_fit.instances().add_element(idxs1[j], c1[j]);
    }

    // amips terms
    AMIPSFrameElement3D_PH::VecI idxs2;
    AMIPSFrameElement3D_PH::VecC c2;
    i = 9 * c_it->idx();
    double w2 = w_amips * V / VM;
    idxs2 << i, i + 3, i + 6, i + 1, i + 4, i + 7, i + 2, i + 5, i + 8;
    // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], det(A) = _c[9], weight w=_c[10], distortion parameters s=_c[11], alpha=_c[12] and beta=_c[13]
    c2 << as, 0, 0, 0, as, 0, 0, 0, as, as * as * as, w_amips * V / VM, 1.0, 1.0, 0.5;
    fe_amips.instances().add_element(idxs2, c2);
    fe_amips_tad_ph.instances().add_element(idxs2, c2);

    // symmetric dirichlet terms
    SymmetricDirichletFrameElement3D_TAD_PH::VecC c3;
    c3 << as, 0, 0, 0, as, 0, 0, 0, as, w_symmetric_dirichlet * V / VM;
    fe_sdirichlet.instances().add_element(idxs2, c3);

//    // optimal rotation terms
//    OptimalRotationFrameElement3D_TAD_PH::VecI idxs3;
//    OptimalRotationFrameElement3D_TAD_PH::VecC c3;
//    c3[0] = w_optrot;
//    c3[1] = 1.0/as;
//    c3[2] = 1.0/as;
//    c3[3] = 1.0/as;
//    c3[4] = 1.0;
//    fe_optrot.instances().add_element(idxs2, c3);
  }


  //setup dihedral angle terms
//  if(w_dihedral > 0.0)
  {
    // iterate over all edges
    for (auto eh: mesh_.edges())
    {
//      EH eh = mesh_.edge_handle(*vohe_it);

//      // skip non-singular edges
//      if(valence_eprop_[eh] == 0)
//        continue;

      // skip boundary edges that are not on feature faces
      if (mesh_.is_boundary(eh) && has_non_feature_boundary_face(eh))
        continue;

      // determine target angle
      double target_angle = double(4 + valence_eprop_[eh]) * 0.5 * M_PI;
      if (mesh_.is_boundary(eh))
      {
        target_angle -= M_PI;
      }

      // get total interior angle
      double int_angle = 2.0 * M_PI;

      if (mesh_.is_boundary(eh))
      {
        int_angle = dihedral_angle(mesh_, mesh_.halfedge_handle(eh, 0));
//        std::cerr << "boundary edge with angle " << int_angle*180.0/M_PI << " and target angle " << target_angle*180.0/M_PI << std::endl;
      }
      // determine uniform scaling of input angles to reach target angle
      double scale_angle = target_angle / int_angle;

      HEH heh = mesh_.halfedge_handle(eh, 0);
      HEH heh_opp = mesh_.opposite_halfedge_handle(heh);
      VH vh0 = mesh_.halfedge(heh).from_vertex();
      VH vh1 = mesh_.halfedge(heh).to_vertex();
      Vec3d pt0 = ovm2eigen(mesh_.vertex(vh0));
      Vec3d pt1 = ovm2eigen(mesh_.vertex(vh1));

      double angle(0.0);
      HEHFIt hehf_it(heh, &mesh_);
      for (; hehf_it.valid(); ++hehf_it)
        if (!mesh_.is_boundary(*hehf_it)) // skip boundary halffaces
        {
          // get halfface
          HFH hfh = *hehf_it;
          // get cell
          CH ch = mesh_.incident_cell(hfh);
          // get other halfface in cell
          HFH hf_adj = mesh_.adjacent_halfface_in_cell(hfh, heh);

          // get vertices and points of both halffaces
          HEH he0n = mesh_.next_halfedge_in_halfface(heh, hfh);
          HEH he1n = mesh_.next_halfedge_in_halfface(heh_opp, hf_adj);
          if (!he0n.is_valid() || !he1n.is_valid())
            std::cerr << "ERROR: next_halfedge_in_halfface is invalid while setting up dihedral terms" << std::endl;

          VH vh00 = mesh_.halfedge(he0n).to_vertex();
          VH vh11 = mesh_.halfedge(he1n).to_vertex();
          Vec3d pt00 = ovm2eigen(mesh_.vertex(vh00));
          Vec3d pt11 = ovm2eigen(mesh_.vertex(vh11));

          // get edge vectors
          Vec3d e0 = pt1 - pt0;
          Vec3d e1 = pt00 - pt0;
          Vec3d e2 = pt11 - pt0;

          // determine target angle
          Vec3d n0 = e1.cross(e0);
          n0 /= n0.norm();
          Vec3d n1 = e2.cross(e0);
          n1 /= n1.norm();
          double dp = n0.dot(n1);
          dp = std::max(-1.0, dp);
          dp = std::min(1.0, dp);
          double alpha = scale_angle * std::acos(dp);

          // get base index
          int idx = 9 * ch.idx();

          DihedralAngleFramesElement_TAD_PH::VecI idxs_dihedral;
          DihedralAngleFramesElement_TAD_PH::VecC c_dihedral;
          idxs_dihedral << idx, idx + 3, idx + 6, idx + 1, idx + 4, idx + 7, idx + 2, idx + 5, idx + 8;
          c_dihedral << e0[0], e0[1], e0[2], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], w_dihedral, s_dihedral, std::cos(
                  alpha);
          fe_dihedral.instances().add_element(idxs_dihedral, c_dihedral);
        }
    }
  }


  // setup alignment terms
  std::map<int, int> cell_face_alignment_constraint;

  // (i) face alignment
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (feature_fprop_[*f_it] > 0)
    {
      // get non-boundary halffaces
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);
      std::vector<HFH> hfhs;
      if (!mesh_.is_boundary(hfh0))
        hfhs.push_back(hfh0);
      if (!mesh_.is_boundary(hfh1))
        hfhs.push_back(hfh1);

      for (size_t i = 0; i < hfhs.size(); ++i)
      {
        HFH hfh = hfhs[i];
        // get corresponding cell
        CH ch = mesh_.incident_cell(hfh);

        // base index of frame
        int bidx = 9 * ch.idx();

        // get vertices
        VH vh0 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[0]).to_vertex();
        VH vh1 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[1]).to_vertex();
        VH vh2 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[2]).to_vertex();

        // get points
        Point p0 = mesh_.vertex(vh0);
        Point p1 = mesh_.vertex(vh1);
        Point p2 = mesh_.vertex(vh2);

        // get outward normal vector
        Vec3d n = ovm2eigen((p1 - p0) % (p2 - p0));
        Vec3d e0 = ovm2eigen(p1 - p0);
        Vec3d e1 = ovm2eigen(p2 - p0);

        double A = 0.5 * n.norm();

        // get local basis of triangle
        Vec3d u = ovm2eigen(p1 - p0);
        u /= u.norm();
        Vec3d v = n.cross(u);
        v /= v.norm();

        // get current gradient vectors
        Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[bidx]));

        // determine alignment axis
        Vec3d Je0 = G.transpose() * e0;
        Vec3d Je1 = G.transpose() * e1;
        Vec3d Jn = Je0.cross(Je1);
        AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Jn);
        int axis_idx = int(e_axis) / 2;
        // test whether field is sufficiently aligned
        Vec3d a = AxisAlignmentHelpers::vector(e_axis);
        if (a.dot(Jn) < Jn.norm() * 0.99985)
          std::cerr << "Warning FF alignment: face " << f_it->idx() << " alignment axis deviates by "
                    << std::acos(a.dot(Jn) / Jn.norm()) * 180.0 / M_PI
                    << " degree\n";

        // memorize axis
        cell_face_alignment_constraint[ch.idx()] = axis_idx;

        double w = w_align * A;

        int bidx2 = bidx + 3 * axis_idx;
        LinearLeastSquaresElement3D::VecI idx(bidx2, bidx2 + 1, bidx2 + 2);
        LinearLeastSquaresElement3D::VecC c0, c1;
        c0 << u[0], u[1], u[2], 0.0, w;
        c1 << v[0], v[1], v[2], 0.0, w;
        fe_face_align.instances().add_element(idx, c0);
        fe_face_align.instances().add_element(idx, c1);

        // generate corresponding constraints
//        double w2 = std::sqrt(A);
        double w2 = 1.0;
        COMISO::LinearConstraint::SVectorNC lc0(nv);
        lc0.coeffRef(idx[0]) = w2 * c0[0];
        lc0.coeffRef(idx[1]) = w2 * c0[1];
        lc0.coeffRef(idx[2]) = w2 * c0[2];
        c_face_align.push_back(COMISO::LinearConstraint(lc0, 0.0));

        COMISO::LinearConstraint::SVectorNC lc1(nv);
        lc1.coeffRef(idx[0]) = w2 * c1[0];
        lc1.coeffRef(idx[1]) = w2 * c1[1];
        lc1.coeffRef(idx[2]) = w2 * c1[2];
        c_face_align.push_back(COMISO::LinearConstraint(lc1, 0.0));
      }
    }

  // (ii) edge alignment
  for (EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
  {
    // singular edge?
    if (valence_eprop_[*e_it] != 0 || feature_eprop_[*e_it] > 0)
    {
      // get all cells incident to the edge
      std::vector<CH> chs;
      for (HECIt hec_it = mesh_.hec_iter(mesh_.halfedge_handle(*e_it, 0)); hec_it.valid(); ++hec_it)
        chs.push_back(*hec_it);

      for (size_t i = 0; i < chs.size(); ++i)
      {
        // get curent cell adjacent to edge
        CH ch = chs[i];

        // base index of frame
        int bidx = 9 * ch.idx();

        // get adjacent vertices
        VH vh0 = mesh_.edge(*e_it).from_vertex();
        VH vh1 = mesh_.edge(*e_it).to_vertex();

        // get edge vector
        Vec3d e = ovm2eigen(mesh_.vertex(vh1) - mesh_.vertex(vh0));

        double l = e.norm();

        // get current gradient vectors
        Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[bidx]));

        // determine alignment axis
        Vec3d Je = G.transpose() * e;
        AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Je);
        int axis_idx = int(e_axis) / 2;
        // test whether field is sufficiently aligned
        Vec3d a = AxisAlignmentHelpers::vector(e_axis);
        if (a.dot(Je) < Je.norm() * 0.99985)
          std::cerr << "Warning FF alignment: edge " << e_it->idx() << " alignment axis deviates by "
                    << std::acos(a.dot(Je) / Je.norm()) * 180.0 / M_PI
                    << " degree\n";

        double w = w_align * VM_cbrt / l;
//        double w2 = std::sqrt(VM_cbrt/l);
        double w2 = 1.0;

        // check if face alignment exists for this cell
        int face_alignment_axis = -1;
        auto cfc = cell_face_alignment_constraint.find(ch.idx());
        if (cfc != cell_face_alignment_constraint.end())
          face_alignment_axis = cfc->second; // ---> skip constraint since linearly dependent

        // the two other axes need to be orthogonal!
        for (int j = 0; j < 3; ++j)
          if (j != axis_idx && j != face_alignment_axis)
          {
            int bidx2 = bidx + 3 * j;
            LinearLeastSquaresElement3D::VecI idx(bidx2, bidx2 + 1, bidx2 + 2);
            LinearLeastSquaresElement3D::VecC c;
            c << e[0], e[1], e[2], 0.0, w;
            fe_edge_align.instances().add_element(idx, c);

            // generate corresponding constraints
            COMISO::LinearConstraint::SVectorNC lc(nv);
            lc.coeffRef(idx[0]) = w2 * c[0];
            lc.coeffRef(idx[1]) = w2 * c[1];
            lc.coeffRef(idx[2]) = w2 * c[2];
            c_edge_align.push_back(COMISO::LinearConstraint(lc, 0.0));
          }
      }
    }
  }

  // setup integrability terms
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (!mesh_.is_boundary(*f_it))
    {
      // get non-boundary halffaces
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);

      // get corresponding cells
      CH ch0 = mesh_.incident_cell(hfh0);
      CH ch1 = mesh_.incident_cell(hfh1);

      // calc volumes of tets
      double V0 = std::abs(volume(ch0));
      double V1 = std::abs(volume(ch1));
      double VF = (V0 + V1) * 0.25;

      double F0_det = frame_cprop_[ch0].determinant();
      double F1_det = frame_cprop_[ch1].determinant();
      double as = std::cbrt(0.5 * (F0_det + F1_det));


      // get local basis of face
      // get vertices
      VH vh0 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[0]).to_vertex();
      VH vh1 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[1]).to_vertex();
      VH vh2 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[2]).to_vertex();

      // get points
      Point p0 = mesh_.vertex(vh0);
      Point p1 = mesh_.vertex(vh1);
      Point p2 = mesh_.vertex(vh2);

      // get outward normal vector
      Vec3d n = ovm2eigen((p1 - p0) % (p2 - p0));
      // get local basis of triangle
      Vec3d u = ovm2eigen(p1 - p0);
      u /= u.norm();
      Vec3d v = n.cross(u);
      v /= v.norm();

      // get transition from ch0 to ch1 such that G0*T01 is w.r.t. coords of G1 (transition for gradients is inverse-transpose of the frame version)
      Mat3d T01 = (tq_.transition_matrix_int(tq_.inverse_transition_idx(transition_hfprop_[hfh0]))).transpose();

      // get permutation and signs of transformed gradients
      Vec3i k;
      Vec3d s;

      if (T01(0, 0) != 0)
      {
        k[0] = 0;
        s[0] = T01(0, 0);
      }
      else if (T01(1, 0) != 0)
      {
        k[0] = 1;
        s[0] = T01(1, 0);
      }
      else if (T01(2, 0) != 0)
      {
        k[0] = 2;
        s[0] = T01(2, 0);
      }

      if (T01(0, 1) != 0)
      {
        k[1] = 0;
        s[1] = T01(0, 1);
      }
      else if (T01(1, 1) != 0)
      {
        k[1] = 1;
        s[1] = T01(1, 1);
      }
      else if (T01(2, 1) != 0)
      {
        k[1] = 2;
        s[1] = T01(2, 1);
      }

      if (T01(0, 2) != 0)
      {
        k[2] = 0;
        s[2] = T01(0, 2);
      }
      else if (T01(1, 2) != 0)
      {
        k[2] = 1;
        s[2] = T01(1, 2);
      }
      else if (T01(2, 2) != 0)
      {
        k[2] = 2;
        s[2] = T01(2, 2);
      }


      // use average cell volume to avoid issues with slivers
      double w = (w_integrability * 0.5 / double(mesh_.n_cells())) * (as * as); // uniform weight
//     double w = (w_integrability * VF/VM)*(as*as);   // geometric weight

//      double w2 = std::sqrt(VF/VM_cbrt);
//      double w2 = 1.0/std::sqrt(double(mesh_.n_cells()));

      // setup integrability constraints for all axes
      for (unsigned int i = 0; i < 3; ++i)
      {
        // get base indices
        int bidx0 = 9 * ch0.idx() + k[i] * 3;
        int bidx1 = 9 * ch1.idx() + i * 3;

        LinearLeastSquaresElement6D::VecI idxs;
        idxs << bidx0, bidx0 + 1, bidx0 + 2, bidx1, bidx1 + 1, bidx1 + 2;

        LinearLeastSquaresElement6D::VecC c0, c1;
        c0 << s[i] * u[0], s[i] * u[1], s[i] * u[2], -u[0], -u[1], -u[2], 0.0, w;
        c1 << s[i] * v[0], s[i] * v[1], s[i] * v[2], -v[0], -v[1], -v[2], 0.0, w;

        fe_integrability.instances().add_element(idxs, c0);
        fe_integrability.instances().add_element(idxs, c1);

        // normalize constraint such that violation is relative to expected length
        Vec3d v1(fe_problem.x()[idxs[0]], fe_problem.x()[idxs[1]], fe_problem.x()[idxs[2]]);
        Vec3d v2(fe_problem.x()[idxs[3]], fe_problem.x()[idxs[4]], fe_problem.x()[idxs[5]]);
        double w2 = 1.0 / (0.5 * (v1.norm() + v2.norm()));

        // generate corresponding constraints
        COMISO::LinearConstraint::SVectorNC lc0(nv);
        lc0.coeffRef(idxs[0]) = w2 * c0[0];
        lc0.coeffRef(idxs[1]) = w2 * c0[1];
        lc0.coeffRef(idxs[2]) = w2 * c0[2];
        lc0.coeffRef(idxs[3]) = w2 * c0[3];
        lc0.coeffRef(idxs[4]) = w2 * c0[4];
        lc0.coeffRef(idxs[5]) = w2 * c0[5];
        c_integrability.push_back(COMISO::LinearConstraint(lc0, 0.0));
        COMISO::LinearConstraint::SVectorNC lc1(nv);
        lc1.coeffRef(idxs[0]) = w2 * c1[0];
        lc1.coeffRef(idxs[1]) = w2 * c1[1];
        lc1.coeffRef(idxs[2]) = w2 * c1[2];
        lc1.coeffRef(idxs[3]) = w2 * c1[3];
        lc1.coeffRef(idxs[4]) = w2 * c1[4];
        lc1.coeffRef(idxs[5]) = w2 * c1[5];
        c_integrability.push_back(COMISO::LinearConstraint(lc1, 0.0));
      }

      double w_smooth = (w_frame_smoothness * VF / VM) * (as * as);   // geometric weight

      // setup smoothness terms for all axes
      for (unsigned int i = 0; i < 3; ++i)
      {
        // get base indices
        int bidx0 = 9 * ch0.idx() + k[i] * 3;
        int bidx1 = 9 * ch1.idx() + i * 3;

        for (unsigned int j = 0; j < 3; ++j)
        {
          LinearLeastSquaresElement2D::VecI idxs;
          idxs << bidx0 + j, bidx1 + j;

          LinearLeastSquaresElement2D::VecC c;
          c << s[i], -1.0, 0.0, w_smooth;

          fe_frame_smoothness.instances().add_element(idxs, c);
        }
      }
    }

  // complete objective function
  {
    fe_problem_complete.add_set(&fe_amips);
    fe_problem_complete.add_set(&fe_amips_tad_ph);
    fe_problem_complete.add_set(&fe_logdet);
    fe_problem_complete.add_set(&fe_frame_fit);
    fe_problem_complete.add_set(&fe_sdirichlet);
    fe_problem_complete.add_set(&fe_dihedral);
    fe_problem_complete.add_set(&fe_edge_align);
    fe_problem_complete.add_set(&fe_face_align);
    fe_problem_complete.add_set(&fe_integrability);
    fe_problem_complete.add_set(&fe_frame_smoothness);
  }

  // add all terms to objective function
  if (amips_active)
    fe_problem.add_set(&fe_amips_tad_ph);

  if (log_det_active)
    fe_problem.add_set(&fe_logdet);

  if (frame_fitting_active)
    fe_problem.add_set(&fe_frame_fit);

  if (symmetric_dirichlet_active)
    fe_problem.add_set(&fe_sdirichlet);

  if (align_penalty_active)
  {
    fe_problem.add_set(&fe_face_align);  // used as constraint ?
    fe_problem.add_set(&fe_edge_align);  // used as constraint ?
  }

  if (integrability_penalty_active)
    fe_problem.add_set(&fe_integrability); // used as constraint ?

  if (frame_smoothness_active)
    fe_problem.add_set(&fe_frame_smoothness);

  // only add if activated
  if (dihedral_active)
    fe_problem.add_set(&fe_dihedral);

//  if(w_optrot > 0.0)
//    fe_problem.add_set(&fe_optrot);

  // collect constraints
  std::vector<COMISO::NConstraintInterface *> constraints;
//  constraints.reserve(c_edge_align.size()+c_face_align.size()+c_integrability.size());
  constraints.reserve(c_integrability.size());
//  for(size_t i=0; i<c_edge_align.size(); ++i)
//    constraints.push_back( &(c_edge_align[i]));
//  for(size_t i=0; i<c_face_align.size(); ++i)
//    constraints.push_back( &(c_face_align[i]));
  for (size_t i = 0; i < c_integrability.size(); ++i)
    constraints.push_back(&(c_integrability[i]));

  std::vector<COMISO::LinearConstraint> linear_constraints;
  linear_constraints.reserve(c_edge_align.size() + c_face_align.size());
  for (size_t i = 0; i < c_edge_align.size(); ++i)
    linear_constraints.push_back(c_edge_align[i]);
  for (size_t i = 0; i < c_face_align.size(); ++i)
    linear_constraints.push_back(c_face_align[i]);

//  for(size_t i=0; i<c_integrability.size(); ++i)
//    linear_constraints.push_back( &(c_integrability[i]));


  std::vector<COMISO::NConstraintInterface *> constraints_all;
  constraints_all.reserve(c_edge_align.size() + c_face_align.size() + c_integrability.size());
  for (size_t i = 0; i < c_edge_align.size(); ++i)
    constraints_all.push_back(&(c_edge_align[i]));
  for (size_t i = 0; i < c_face_align.size(); ++i)
    constraints_all.push_back(&(c_face_align[i]));
//  for(size_t i=0; i<c_integrability.size(); ++i)
//    constraints_all.push_back( &(c_integrability[i]));

//  // make integrability constraints linearly independent
//  COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraints_all);

//  // debug constraints
//  {
//    COMISO::TruncatedNewtonPCG::SMatrixD A;
//    COMISO::TruncatedNewtonPCG::VectorD  b;
//
//    COMISO::LinearConstraintConverter::nsolver_to_eigen(linear_constraints, A, b);
//
//    std::cerr << "||A| = " << A.norm() << std::endl;
//    std::cerr << "Constraint Matrix" << std::endl << A << std::endl;
//    std::cerr << "Constraint rhs" << std::endl    << b << std::endl;
//  }


  if (0) //DEBUG
  {
    // check derivatives
    COMISO::NPDerivativeChecker npdc;
    npdc.config().dx = 1e-6;
    npdc.check_all(&fe_problem);
  }

  // print energies at beginning
  {
    std::cerr << "--- initial energies complete --- \n";
    fe_problem_complete.x() = fe_problem.x();
    fe_problem_complete.print_objectives();

    std::cerr << "--- initial energies selected --- \n";
    fe_problem.print_objectives();
    std::cerr << "-------------------------------\n";
  }

  // optimize ***********
  COMISO::NPTiming npp2(&fe_problem); // with timings
  Eigen::Map <Eigen::VectorXd> x0(fe_problem.x().data(), fe_problem.x().size());
  double tol = 2e2 / x0.norm();
  std::cerr << "||x0|| = " << x0.norm() << std::endl;
  std::cerr << "tolerance = " << tol << std::endl;

  if (1)
  {
    COMISO::TruncatedNewtonPCG tn(tol);
    tn.always_update_preconditioner() = true;
//  tn.adaptive_tolerance_modifier() = 0.01;
    tn.max_pcg_iters() = 300;
    tn.allow_warmstart() = true;
    tn.tolerance_gdx() = 0.1;
    tn.solve_projected_normal_equation(&npp2, linear_constraints);
  }

  // test different truncated Newton settings
  if (0)
  {
    std::vector<double> x_bak = fe_problem.x();
    std::vector<int> max_iters;
    max_iters.push_back(10);
    max_iters.push_back(20);
    max_iters.push_back(50);
    max_iters.push_back(100);
    max_iters.push_back(200);
    max_iters.push_back(500);
    max_iters.push_back(1000);

    for (int i = 0; i < max_iters.size(); ++i)
    {
      // restore initial
      fe_problem.x() = x_bak;
      std::cerr << "************ test max_pcg_iters = " << max_iters[i] << std::endl;
      // Truncated Newton with fixed penalty
      tol = 1e2 / x0.norm();
      COMISO::TruncatedNewtonPCG tn2(tol);
      tn2.always_update_preconditioner() = true;
//  tn.adaptive_tolerance_modifier() = 0.01;
      tn2.max_pcg_iters() = max_iters[i];
      COMISO::NPTiming npp3(&fe_problem); // with timings
//    tn.solve_projected_normal_equation(&npp2, linear_constraints);
      tn2.solve_projected_normal_equation(&npp3, linear_constraints);

      // print energies at end
      {
        std::cerr << "max_pcg_iters = " << max_iters[i] << std::endl;
        std::cerr << "--- final energies complete --- \n";
        fe_problem_complete.x() = fe_problem.x();
        fe_problem_complete.print_objectives();

        std::cerr << "--- final energies selected --- \n";
        fe_problem.print_objectives();
        std::cerr << "-------------------------------\n";
      }

    }
  }

  if (0)
  {
    COMISO::TrustregionNewtonPCG trn(tol);
    trn.solve(&npp2, linear_constraints);
  }

//  // ALM
//  COMISO::AugmentedLagrangianMethod alm(10.0, 1e-1, 1e-6, 7);
//  alm.solve_experimental(&npp2, constraints, linear_constraints);

//  // Newton infeasible start
//  COMISO::NewtonSolver ns;
//  ns.solve_infeasible_start(&npp2, linear_constraints);
////  ns.solve_infeasible_start(&npp2, constraints_all);

  // IPOPT
//  try {
//    COMISO::IPOPTSolver ipsol;
////    ipsol.set_ipopt_option("print_level", 0);
////    ipsol.set_max_iterations(100);
////    ipsol.set_ipopt_option("tol", 1e-4);
////    ipsol.set_ipopt_option("hessian_approximation", "limited-memory");
////    ipsol.set_ipopt_option("limited_memory_max_history", 20);
//    ipsol.solve(&npp2, constraints_all);
//  }
//  catch (...)
//  {
//    std::cerr << "ERROR: IPOPT threw an exception!!!!!" << std::endl;
//  }

  // re-solve with higher weights on integrability
  if (0)
  {
    // boost integrability weight
    for (unsigned int i = 0; i < fe_integrability.instances().size(); ++i)
      fe_integrability.instances().c(i)[7] *= 100.0;

    try
    {
      COMISO::IPOPTSolver ipsol;
//    ipsol.set_ipopt_option("print_level", 0);
      ipsol.set_max_iterations(100);
//    ipsol.set_ipopt_option("tol", 1e-4);
//    ipsol.set_ipopt_option("hessian_approximation", "limited-memory");
//    ipsol.set_ipopt_option("limited_memory_max_history", 20);
      ipsol.solve(&npp2, constraints_all);
    }
    catch (...)
    {
      std::cerr << "ERROR: IPOPT threw an exception!!!!!" << std::endl;
    }
  }

//    COMISO::AugmentedLagrangianMethod alm(1e3, 1e-1, 1e-2, 7);
//    COMISO::NPTiming npp_alm3(&fe_problem_alm3); // with timings
//    fe_problem_alm3.x() = x0;
//    alm.solve_experimental(&npp_alm3, constraints, linear_constraints);
//    std::cerr << "******** finished via TinyAD derivatives ***********" << std::endl;
//    fe_problem.x() = fe_problem_alm3.x();

//    // final high-precision solve
//    if(0)
//    {
//      // make integrability constraints linearly independent
//      COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraints_all);
//
//      COMISO::NewtonSolver ns(1e-6, 1e-9, 2000);
////  ns.set_prescribed_constraint_regularization(-1e-9);
////  ns.set_prescribed_hessian_regularization(1e-9);
////  ns.solve_infeasible_start(&fe_problem, linear_constraints);
//      ns.solve_infeasible_start(&npp_alm3, constraints_all);
//
//      fe_problem.x() = fe_problem_alm3.x();
//    }
//  }

  // print energies at end
  {
    std::cerr << "--- final energies complete --- \n";
    fe_problem_complete.x() = fe_problem.x();
    fe_problem_complete.print_objectives();

    std::cerr << "--- final energies selected --- \n";
    fe_problem.print_objectives();
    std::cerr << "-------------------------------\n";
  }

  // export solution to frame field
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    int idx = 9 * c_it->idx();
    frame_cprop_[*c_it] = (Eigen::Map<Mat3d>((double *) &(fe_problem.x()[idx])).inverse()).transpose();
  }
}

template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::
optimize_integrable_regularized(const double _anisotropy_alpha, const double _w_frame_smoothness,
                                const double _w_integrability)
{
  ScopedStopWatch sw_total(sw::frame_field_int_opt);

  std::cout << "#####Optimize Frame Field regularized..." << std::endl;

  // check whether input is valid
  check_valence_consistency();
  check_frame_rotation_angles();

  const bool amips_active = false;
  const bool log_det_active = false;
  const bool frame_fitting_active = false;
  const bool frame_smoothness_active = true;
  const bool align_penalty_active = false;
  const bool integrability_penalty_active = true;
  const bool dihedral_active = false;
  const bool symmetric_dirichlet_active = true;

  const double w_amips = 1.0;
  const double w_log_det = 1.0;
  const double w_frame_fitting = 1.0;
  const double w_frame_smoothness = _w_frame_smoothness;
  const double w_align = 1.0;
  const double w_integrability = _w_integrability;
  const double alpha_fitting = 1.0;
  const double w_dihedral = 1.0;
  const double s_dihedral = 4.0;
  const double w_symmetric_dirichlet = 1.0;
//  const double w_optrot   = 0.0;

  // variables are 3x3 frames per cell
  int nv = 9 * mesh_.n_cells();

  if (nv == 0)
  {
    std::cerr << "skip optimize_integrable since mesh does not have cells\n";
    return;
  }

  double VM = mesh_volume();
  double VM_cbrt = std::cbrt(VM);

  // create problem
  COMISO::FiniteElementProblem fe_problem(nv);
  COMISO::FiniteElementProblem fe_problem_hc(nv);
  COMISO::FiniteElementProblem fe_problem_complete(nv);

  // set gradient field as initial solution (gradient field is inverse-transpose of frame field)
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    int idx = 9 * c_it->idx();
    Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[idx]));
    G = (frame_cprop_[*c_it].inverse()).transpose();
  }

  // objective terms
  COMISO::FiniteElementSet <AMIPSFrameElement3D_PH> fe_amips("AMIPS");
  COMISO::FiniteElementSet <AMIPSFrameElement3D_TAD_PH> fe_amips_tad_ph("AMIPS_TinyAD");
  COMISO::FiniteElementSet <LogDetElement3D_PH> fe_logdet("LogDet");
  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_frame_fit("FrameFit");
  COMISO::FiniteElementSet <DihedralAngleFramesElement_TAD_PH> fe_dihedral("DihedralAngle");
  COMISO::FiniteElementSet <SymmetricDirichletFrameElement3D_TAD_PH> fe_sdirichlet("SymmetricDirichlet");

  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_face_align("FaceAlignment");
  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_edge_align("EdgeAlignment");
  COMISO::FiniteElementSet <LinearLeastSquaresElement6D> fe_integrability("Integrability");
  COMISO::FiniteElementSet <LinearLeastSquaresElement2D> fe_frame_smoothness("Frame Smoothness");

//  COMISO::FiniteElementSet<OptimalRotationFrameElement3D_TAD_PH> fe_optrot;

  // constraints
  std::vector<COMISO::LinearConstraint> c_face_align;
  std::vector<COMISO::LinearConstraint> c_edge_align;
  std::vector<COMISO::LinearConstraint> c_integrability;

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // calc volume of tet
    double V = std::abs(volume(*c_it));
    double F_det = frame_cprop_[*c_it].determinant();
    double as = std::cbrt(F_det);

    const double offset = 0.001;
//    const double scale = F_det * (1.0+offset);
    const double scale = F_det;

    // log determinant terms
    LogDetElement3D_PH::VecI idxs0;
    LogDetElement3D_PH::VecC c0;
    int i = 9 * c_it->idx();
    idxs0 << i, i + 1, i + 2, i + 3, i + 4, i + 5, i + 6, i + 7, i + 8;
    c0 << w_log_det, scale, offset; // weight, scale, offset
    fe_logdet.instances().add_element(idxs0, c0);

    // frame fitting terms
    std::vector<LinearLeastSquaresElement3D::VecI> idxs1;
    std::vector<LinearLeastSquaresElement3D::VecC> c1;

    double w1 = (w_frame_fitting * V / VM) * (as * as);
    // get 9 instances of LinearElement
    LinearLeastSquaresElement3D::get_constants_and_indices_frame_fitting(
            Eigen::Map<Mat3d>((double *) &(fe_problem.x()[i])),
            c_it->idx(),
            w1,
            alpha_fitting,
            idxs1,
            c1);
    for (size_t j = 0; j < idxs1.size(); ++j)
    {
      fe_frame_fit.instances().add_element(idxs1[j], c1[j]);
    }

    // amips terms
    AMIPSFrameElement3D_PH::VecI idxs2;
    AMIPSFrameElement3D_PH::VecC c2;
    i = 9 * c_it->idx();
    double w2 = w_amips * V / VM;
    idxs2 << i, i + 3, i + 6, i + 1, i + 4, i + 7, i + 2, i + 5, i + 8;
    // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], det(A) = _c[9], weight w=_c[10], distortion parameters s=_c[11], alpha=_c[12] and beta=_c[13]
    c2 << as, 0, 0, 0, as, 0, 0, 0, as, as * as * as, w_amips * V / VM, 1.0, 1.0, 0.5;
    fe_amips.instances().add_element(idxs2, c2);
    fe_amips_tad_ph.instances().add_element(idxs2, c2);

    // symmetric dirichlet terms
    SymmetricDirichletFrameElement3D_TAD_PH::VecC c3;
    c3 << as, 0, 0, 0, as, 0, 0, 0, as, w_symmetric_dirichlet * V / VM;
    fe_sdirichlet.instances().add_element(idxs2, c3);

//    // optimal rotation terms
//    OptimalRotationFrameElement3D_TAD_PH::VecI idxs3;
//    OptimalRotationFrameElement3D_TAD_PH::VecC c3;
//    c3[0] = w_optrot;
//    c3[1] = 1.0/as;
//    c3[2] = 1.0/as;
//    c3[3] = 1.0/as;
//    c3[4] = 1.0;
//    fe_optrot.instances().add_element(idxs2, c3);
  }


  //setup dihedral angle terms
//  if(w_dihedral > 0.0)
  {
    // iterate over all edges
    for (auto eh: mesh_.edges())
    {
//      EH eh = mesh_.edge_handle(*vohe_it);

//      // skip non-singular edges
//      if(valence_eprop_[eh] == 0)
//        continue;

      // skip boundary edges that are not on feature faces
      if (mesh_.is_boundary(eh) && has_non_feature_boundary_face(eh))
        continue;

      // determine target angle
      double target_angle = double(4 + valence_eprop_[eh]) * 0.5 * M_PI;
      if (mesh_.is_boundary(eh))
      {
        target_angle -= M_PI;
      }

      // get total interior angle
      double int_angle = 2.0 * M_PI;

      if (mesh_.is_boundary(eh))
      {
        int_angle = dihedral_angle(mesh_, mesh_.halfedge_handle(eh, 0));
//        std::cerr << "boundary edge with angle " << int_angle*180.0/M_PI << " and target angle " << target_angle*180.0/M_PI << std::endl;
      }
      // determine uniform scaling of input angles to reach target angle
      double scale_angle = target_angle / int_angle;

      HEH heh = mesh_.halfedge_handle(eh, 0);
      HEH heh_opp = mesh_.opposite_halfedge_handle(heh);
      VH vh0 = mesh_.halfedge(heh).from_vertex();
      VH vh1 = mesh_.halfedge(heh).to_vertex();
      Vec3d pt0 = ovm2eigen(mesh_.vertex(vh0));
      Vec3d pt1 = ovm2eigen(mesh_.vertex(vh1));

      double angle(0.0);
      HEHFIt hehf_it(heh, &mesh_);
      for (; hehf_it.valid(); ++hehf_it)
        if (!mesh_.is_boundary(*hehf_it)) // skip boundary halffaces
        {
          // get halfface
          HFH hfh = *hehf_it;
          // get cell
          CH ch = mesh_.incident_cell(hfh);
          // get other halfface in cell
          HFH hf_adj = mesh_.adjacent_halfface_in_cell(hfh, heh);

          // get vertices and points of both halffaces
          HEH he0n = mesh_.next_halfedge_in_halfface(heh, hfh);
          HEH he1n = mesh_.next_halfedge_in_halfface(heh_opp, hf_adj);
          if (!he0n.is_valid() || !he1n.is_valid())
            std::cerr << "ERROR: next_halfedge_in_halfface is invalid while setting up dihedral terms" << std::endl;

          VH vh00 = mesh_.halfedge(he0n).to_vertex();
          VH vh11 = mesh_.halfedge(he1n).to_vertex();
          Vec3d pt00 = ovm2eigen(mesh_.vertex(vh00));
          Vec3d pt11 = ovm2eigen(mesh_.vertex(vh11));

          // get edge vectors
          Vec3d e0 = pt1 - pt0;
          Vec3d e1 = pt00 - pt0;
          Vec3d e2 = pt11 - pt0;

          // determine target angle
          Vec3d n0 = e1.cross(e0);
          n0 /= n0.norm();
          Vec3d n1 = e2.cross(e0);
          n1 /= n1.norm();
          double dp = n0.dot(n1);
          dp = std::max(-1.0, dp);
          dp = std::min(1.0, dp);
          double alpha = scale_angle * std::acos(dp);

          // get base index
          int idx = 9 * ch.idx();

          DihedralAngleFramesElement_TAD_PH::VecI idxs_dihedral;
          DihedralAngleFramesElement_TAD_PH::VecC c_dihedral;
          idxs_dihedral << idx, idx + 3, idx + 6, idx + 1, idx + 4, idx + 7, idx + 2, idx + 5, idx + 8;
          c_dihedral << e0[0], e0[1], e0[2], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], w_dihedral, s_dihedral, std::cos(
                  alpha);
          fe_dihedral.instances().add_element(idxs_dihedral, c_dihedral);
        }
    }
  }


  // setup alignment terms
  std::map<int, int> cell_face_alignment_constraint;

  // (i) face alignment
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (feature_fprop_[*f_it] > 0)
    {
      // get non-boundary halffaces
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);
      std::vector<HFH> hfhs;
      if (!mesh_.is_boundary(hfh0))
        hfhs.push_back(hfh0);
      if (!mesh_.is_boundary(hfh1))
        hfhs.push_back(hfh1);

      for (size_t i = 0; i < hfhs.size(); ++i)
      {
        HFH hfh = hfhs[i];
        // get corresponding cell
        CH ch = mesh_.incident_cell(hfh);

        // base index of frame
        int bidx = 9 * ch.idx();

        // get vertices
        VH vh0 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[0]).to_vertex();
        VH vh1 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[1]).to_vertex();
        VH vh2 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[2]).to_vertex();

        // get points
        Point p0 = mesh_.vertex(vh0);
        Point p1 = mesh_.vertex(vh1);
        Point p2 = mesh_.vertex(vh2);

        // get outward normal vector
        Vec3d n = ovm2eigen((p1 - p0) % (p2 - p0));
        Vec3d e0 = ovm2eigen(p1 - p0);
        Vec3d e1 = ovm2eigen(p2 - p0);

        double A = 0.5 * n.norm();

        // get local basis of triangle
        Vec3d u = ovm2eigen(p1 - p0);
        u /= u.norm();
        Vec3d v = n.cross(u);
        v /= v.norm();

        // get current gradient vectors
        Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[bidx]));

        // determine alignment axis
        Vec3d Je0 = G.transpose() * e0;
        Vec3d Je1 = G.transpose() * e1;
        Vec3d Jn = Je0.cross(Je1);
        AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Jn);
        int axis_idx = int(e_axis) / 2;
        // test whether field is sufficiently aligned
        Vec3d a = AxisAlignmentHelpers::vector(e_axis);
        if (a.dot(Jn) < Jn.norm() * 0.99985)
          std::cerr << "Warning FF alignment: face " << f_it->idx() << " alignment axis deviates by "
                    << std::acos(a.dot(Jn) / Jn.norm()) * 180.0 / M_PI
                    << " degree\n";

        // memorize axis
        cell_face_alignment_constraint[ch.idx()] = axis_idx;

        double w = w_align * A;

        int bidx2 = bidx + 3 * axis_idx;
        LinearLeastSquaresElement3D::VecI idx(bidx2, bidx2 + 1, bidx2 + 2);
        LinearLeastSquaresElement3D::VecC c0, c1;
        c0 << u[0], u[1], u[2], 0.0, w;
        c1 << v[0], v[1], v[2], 0.0, w;
        fe_face_align.instances().add_element(idx, c0);
        fe_face_align.instances().add_element(idx, c1);

        // generate corresponding constraints
//        double w2 = std::sqrt(A);
        double w2 = 1.0;
        COMISO::LinearConstraint::SVectorNC lc0(nv);
        lc0.coeffRef(idx[0]) = w2 * c0[0];
        lc0.coeffRef(idx[1]) = w2 * c0[1];
        lc0.coeffRef(idx[2]) = w2 * c0[2];
        c_face_align.push_back(COMISO::LinearConstraint(lc0, 0.0));

        COMISO::LinearConstraint::SVectorNC lc1(nv);
        lc1.coeffRef(idx[0]) = w2 * c1[0];
        lc1.coeffRef(idx[1]) = w2 * c1[1];
        lc1.coeffRef(idx[2]) = w2 * c1[2];
        c_face_align.push_back(COMISO::LinearConstraint(lc1, 0.0));
      }
    }

  // (ii) edge alignment
  for (EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
  {
    // singular edge?
    if (valence_eprop_[*e_it] != 0 || feature_eprop_[*e_it] > 0)
    {
      // get all cells incident to the edge
      std::vector<CH> chs;
      for (HECIt hec_it = mesh_.hec_iter(mesh_.halfedge_handle(*e_it, 0)); hec_it.valid(); ++hec_it)
        chs.push_back(*hec_it);

      for (size_t i = 0; i < chs.size(); ++i)
      {
        // get curent cell adjacent to edge
        CH ch = chs[i];

        // base index of frame
        int bidx = 9 * ch.idx();

        // get adjacent vertices
        VH vh0 = mesh_.edge(*e_it).from_vertex();
        VH vh1 = mesh_.edge(*e_it).to_vertex();

        // get edge vector
        Vec3d e = ovm2eigen(mesh_.vertex(vh1) - mesh_.vertex(vh0));

        double l = e.norm();

        // get current gradient vectors
        Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[bidx]));

        // determine alignment axis
        Vec3d Je = G.transpose() * e;
        AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Je);
        int axis_idx = int(e_axis) / 2;
        // test whether field is sufficiently aligned
        Vec3d a = AxisAlignmentHelpers::vector(e_axis);
        if (a.dot(Je) < Je.norm() * 0.99985)
          std::cerr << "Warning FF alignment: edge " << e_it->idx() << " alignment axis deviates by "
                    << std::acos(a.dot(Je) / Je.norm()) * 180.0 / M_PI
                    << " degree\n";

        double w = w_align * VM_cbrt / l;
//        double w2 = std::sqrt(VM_cbrt/l);
        double w2 = 1.0;

        // check if face alignment exists for this cell
        int face_alignment_axis = -1;
        auto cfc = cell_face_alignment_constraint.find(ch.idx());
        if (cfc != cell_face_alignment_constraint.end())
          face_alignment_axis = cfc->second; // ---> skip constraint since linearly dependent

        // the two other axes need to be orthogonal!
        for (int j = 0; j < 3; ++j)
          if (j != axis_idx && j != face_alignment_axis)
          {
            int bidx2 = bidx + 3 * j;
            LinearLeastSquaresElement3D::VecI idx(bidx2, bidx2 + 1, bidx2 + 2);
            LinearLeastSquaresElement3D::VecC c;
            c << e[0], e[1], e[2], 0.0, w;
            fe_edge_align.instances().add_element(idx, c);

            // generate corresponding constraints
            COMISO::LinearConstraint::SVectorNC lc(nv);
            lc.coeffRef(idx[0]) = w2 * c[0];
            lc.coeffRef(idx[1]) = w2 * c[1];
            lc.coeffRef(idx[2]) = w2 * c[2];
            c_edge_align.push_back(COMISO::LinearConstraint(lc, 0.0));
          }
      }
    }
  }

  // setup integrability terms
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (!mesh_.is_boundary(*f_it))
    {
      // get non-boundary halffaces
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);

      // get corresponding cells
      CH ch0 = mesh_.incident_cell(hfh0);
      CH ch1 = mesh_.incident_cell(hfh1);

      // calc volumes of tets
      double V0 = std::abs(volume(ch0));
      double V1 = std::abs(volume(ch1));
      double VF = (V0 + V1) * 0.25;

      double F0_det = frame_cprop_[ch0].determinant();
      double F1_det = frame_cprop_[ch1].determinant();
      double as = std::cbrt(0.5 * (F0_det + F1_det));


      // get local basis of face
      // get vertices
      VH vh0 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[0]).to_vertex();
      VH vh1 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[1]).to_vertex();
      VH vh2 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[2]).to_vertex();

      // get points
      Point p0 = mesh_.vertex(vh0);
      Point p1 = mesh_.vertex(vh1);
      Point p2 = mesh_.vertex(vh2);

      // get outward normal vector
      Vec3d n = ovm2eigen((p1 - p0) % (p2 - p0));
      // get local basis of triangle
      Vec3d u = ovm2eigen(p1 - p0);
      u /= u.norm();
      Vec3d v = n.cross(u);
      v /= v.norm();

      // get transition from ch0 to ch1 such that G0*T01 is w.r.t. coords of G1 (transition for gradients is inverse-transpose of the frame version)
      Mat3d T01 = (tq_.transition_matrix_int(tq_.inverse_transition_idx(transition_hfprop_[hfh0]))).transpose();

      // get permutation and signs of transformed gradients
      Vec3i k;
      Vec3d s;

      if (T01(0, 0) != 0)
      {
        k[0] = 0;
        s[0] = T01(0, 0);
      }
      else if (T01(1, 0) != 0)
      {
        k[0] = 1;
        s[0] = T01(1, 0);
      }
      else if (T01(2, 0) != 0)
      {
        k[0] = 2;
        s[0] = T01(2, 0);
      }

      if (T01(0, 1) != 0)
      {
        k[1] = 0;
        s[1] = T01(0, 1);
      }
      else if (T01(1, 1) != 0)
      {
        k[1] = 1;
        s[1] = T01(1, 1);
      }
      else if (T01(2, 1) != 0)
      {
        k[1] = 2;
        s[1] = T01(2, 1);
      }

      if (T01(0, 2) != 0)
      {
        k[2] = 0;
        s[2] = T01(0, 2);
      }
      else if (T01(1, 2) != 0)
      {
        k[2] = 1;
        s[2] = T01(1, 2);
      }
      else if (T01(2, 2) != 0)
      {
        k[2] = 2;
        s[2] = T01(2, 2);
      }


      // use average cell volume to avoid issues with slivers
      double w = (w_integrability * 0.5 / double(mesh_.n_cells())) * (as * as); // uniform weight
//     double w = (w_integrability * VF/VM)*(as*as);   // geometric weight

//      double w2 = std::sqrt(VF/VM_cbrt);
//      double w2 = 1.0/std::sqrt(double(mesh_.n_cells()));

      // setup integrability constraints for all axes
      for (unsigned int i = 0; i < 3; ++i)
      {
        // get base indices
        int bidx0 = 9 * ch0.idx() + k[i] * 3;
        int bidx1 = 9 * ch1.idx() + i * 3;

        LinearLeastSquaresElement6D::VecI idxs;
        idxs << bidx0, bidx0 + 1, bidx0 + 2, bidx1, bidx1 + 1, bidx1 + 2;

        LinearLeastSquaresElement6D::VecC c0, c1;
        c0 << s[i] * u[0], s[i] * u[1], s[i] * u[2], -u[0], -u[1], -u[2], 0.0, w;
        c1 << s[i] * v[0], s[i] * v[1], s[i] * v[2], -v[0], -v[1], -v[2], 0.0, w;

        fe_integrability.instances().add_element(idxs, c0);
        fe_integrability.instances().add_element(idxs, c1);

        // normalize constraint such that violation is relative to expected length
        Vec3d v1(fe_problem.x()[idxs[0]], fe_problem.x()[idxs[1]], fe_problem.x()[idxs[2]]);
        Vec3d v2(fe_problem.x()[idxs[3]], fe_problem.x()[idxs[4]], fe_problem.x()[idxs[5]]);
        double w2 = 1.0 / (0.5 * (v1.norm() + v2.norm()));

        // generate corresponding constraints
        COMISO::LinearConstraint::SVectorNC lc0(nv);
        lc0.coeffRef(idxs[0]) = w2 * c0[0];
        lc0.coeffRef(idxs[1]) = w2 * c0[1];
        lc0.coeffRef(idxs[2]) = w2 * c0[2];
        lc0.coeffRef(idxs[3]) = w2 * c0[3];
        lc0.coeffRef(idxs[4]) = w2 * c0[4];
        lc0.coeffRef(idxs[5]) = w2 * c0[5];
        c_integrability.push_back(COMISO::LinearConstraint(lc0, 0.0));
        COMISO::LinearConstraint::SVectorNC lc1(nv);
        lc1.coeffRef(idxs[0]) = w2 * c1[0];
        lc1.coeffRef(idxs[1]) = w2 * c1[1];
        lc1.coeffRef(idxs[2]) = w2 * c1[2];
        lc1.coeffRef(idxs[3]) = w2 * c1[3];
        lc1.coeffRef(idxs[4]) = w2 * c1[4];
        lc1.coeffRef(idxs[5]) = w2 * c1[5];
        c_integrability.push_back(COMISO::LinearConstraint(lc1, 0.0));
      }

      double w_smooth = (w_frame_smoothness * VF / VM) * (as * as);   // geometric weight

      // setup smoothness terms for all axes
      for (unsigned int i = 0; i < 3; ++i)
      {
        // get base indices
        int bidx0 = 9 * ch0.idx() + k[i] * 3;
        int bidx1 = 9 * ch1.idx() + i * 3;

        for (unsigned int j = 0; j < 3; ++j)
        {
          LinearLeastSquaresElement2D::VecI idxs;
          idxs << bidx0 + j, bidx1 + j;

          LinearLeastSquaresElement2D::VecC c;
          c << s[i], -1.0, 0.0, w_smooth;

          fe_frame_smoothness.instances().add_element(idxs, c);
        }
      }

    }

  // complete objective function
  {
    fe_problem_complete.add_set(&fe_amips);
    fe_problem_complete.add_set(&fe_amips_tad_ph);
    fe_problem_complete.add_set(&fe_logdet);
    fe_problem_complete.add_set(&fe_frame_fit);
    fe_problem_complete.add_set(&fe_sdirichlet);
    fe_problem_complete.add_set(&fe_dihedral);
    fe_problem_complete.add_set(&fe_edge_align);
    fe_problem_complete.add_set(&fe_face_align);
    fe_problem_complete.add_set(&fe_integrability);
    fe_problem_complete.add_set(&fe_frame_smoothness);
  }

  // add all terms to objective function
  if (amips_active)
    fe_problem.add_set(&fe_amips_tad_ph);

  if (log_det_active)
    fe_problem.add_set(&fe_logdet);

  if (frame_fitting_active)
    fe_problem.add_set(&fe_frame_fit);

  if (symmetric_dirichlet_active)
    fe_problem.add_set(&fe_sdirichlet);

  if (align_penalty_active)
  {
    fe_problem.add_set(&fe_face_align);  // used as constraint ?
    fe_problem.add_set(&fe_edge_align);  // used as constraint ?
  }

  if (integrability_penalty_active)
    fe_problem.add_set(&fe_integrability); // used as constraint ?

  if (frame_smoothness_active)
    fe_problem.add_set(&fe_frame_smoothness);

  // only add if activated
  if (dihedral_active)
    fe_problem.add_set(&fe_dihedral);

//  if(w_optrot > 0.0)
//    fe_problem.add_set(&fe_optrot);

  // problem with hard constraints
  fe_problem_hc.add_set(&fe_sdirichlet);
  fe_problem_hc.add_set(&fe_frame_smoothness);

  // collect constraints
  std::vector<COMISO::NConstraintInterface *> constraints;
//  constraints.reserve(c_edge_align.size()+c_face_align.size()+c_integrability.size());
  constraints.reserve(c_integrability.size());
//  for(size_t i=0; i<c_edge_align.size(); ++i)
//    constraints.push_back( &(c_edge_align[i]));
//  for(size_t i=0; i<c_face_align.size(); ++i)
//    constraints.push_back( &(c_face_align[i]));
  for (size_t i = 0; i < c_integrability.size(); ++i)
    constraints.push_back(&(c_integrability[i]));

  std::vector<COMISO::LinearConstraint> linear_constraints;
  linear_constraints.reserve(c_edge_align.size() + c_face_align.size());
  for (size_t i = 0; i < c_edge_align.size(); ++i)
    linear_constraints.push_back(c_edge_align[i]);
  for (size_t i = 0; i < c_face_align.size(); ++i)
    linear_constraints.push_back(c_face_align[i]);

//  for(size_t i=0; i<c_integrability.size(); ++i)
//    linear_constraints.push_back( &(c_integrability[i]));


  std::vector<COMISO::NConstraintInterface *> constraints_all;
  constraints_all.reserve(c_edge_align.size() + c_face_align.size() + c_integrability.size());
  for (size_t i = 0; i < c_edge_align.size(); ++i)
    constraints_all.push_back(&(c_edge_align[i]));
  for (size_t i = 0; i < c_face_align.size(); ++i)
    constraints_all.push_back(&(c_face_align[i]));
  for (size_t i = 0; i < c_integrability.size(); ++i)
    constraints_all.push_back(&(c_integrability[i]));

//  // make integrability constraints linearly independent
//  COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraints_all);

//  // debug constraints
//  {
//    COMISO::TruncatedNewtonPCG::SMatrixD A;
//    COMISO::TruncatedNewtonPCG::VectorD  b;
//
//    COMISO::LinearConstraintConverter::nsolver_to_eigen(linear_constraints, A, b);
//
//    std::cerr << "||A| = " << A.norm() << std::endl;
//    std::cerr << "Constraint Matrix" << std::endl << A << std::endl;
//    std::cerr << "Constraint rhs" << std::endl    << b << std::endl;
//  }


  if (0) //DEBUG
  {
    // check derivatives
    COMISO::NPDerivativeChecker npdc;
    npdc.config().dx = 1e-6;
    npdc.check_all(&fe_problem);
  }

  // print energies at beginning
  {
    std::cerr << "--- initial energies complete --- \n";
    fe_problem_complete.x() = fe_problem.x();
    fe_problem_complete.print_objectives();

    std::cerr << "--- initial energies selected --- \n";
    fe_problem.print_objectives();
    std::cerr << "-------------------------------\n";
  }

  // optimize ***********
  COMISO::NPTiming npp2(&fe_problem); // with timings
  Eigen::Map <Eigen::VectorXd> x0(fe_problem.x().data(), fe_problem.x().size());
  double tol = 2e2 / x0.norm();
  std::cerr << "||x0|| = " << x0.norm() << std::endl;
  std::cerr << "tolerance = " << tol << std::endl;

  if (1)
  {
    COMISO::TruncatedNewtonPCG tn(1e-6);
    tn.always_update_preconditioner() = true;
//  tn.adaptive_tolerance_modifier() = 0.01;
    tn.max_pcg_iters() = 300;
    tn.allow_warmstart() = true;
    tn.tolerance_gdx() = 0.1;
    tn.solve_projected_normal_equation(&npp2, linear_constraints);
  }

  if (0)
  {
    COMISO::NPTiming npp3(&fe_problem_hc); // with timings
    COMISO::NewtonSolver ns(1e-6);
    ns.set_linearsolver(COMISO::NewtonSolver::LS_OSQP_DIRECT);
    ns.solve_infeasible_start(&npp3, constraints_all);
    fe_problem.x() = fe_problem_hc.x();
  }


  // test different truncated Newton settings
  if (0)
  {
    std::vector<double> x_bak = fe_problem.x();
    std::vector<int> max_iters;
    max_iters.push_back(10);
    max_iters.push_back(20);
    max_iters.push_back(50);
    max_iters.push_back(100);
    max_iters.push_back(200);
    max_iters.push_back(500);
    max_iters.push_back(1000);

    for (int i = 0; i < max_iters.size(); ++i)
    {
      // restore initial
      fe_problem.x() = x_bak;
      std::cerr << "************ test max_pcg_iters = " << max_iters[i] << std::endl;
      // Truncated Newton with fixed penalty
      tol = 1e2 / x0.norm();
      COMISO::TruncatedNewtonPCG tn2(tol);
      tn2.always_update_preconditioner() = true;
//  tn.adaptive_tolerance_modifier() = 0.01;
      tn2.max_pcg_iters() = max_iters[i];
      COMISO::NPTiming npp3(&fe_problem); // with timings
//    tn.solve_projected_normal_equation(&npp2, linear_constraints);
      tn2.solve_projected_normal_equation(&npp3, linear_constraints);

      // print energies at end
      {
        std::cerr << "max_pcg_iters = " << max_iters[i] << std::endl;
        std::cerr << "--- final energies complete --- \n";
        fe_problem_complete.x() = fe_problem.x();
        fe_problem_complete.print_objectives();

        std::cerr << "--- final energies selected --- \n";
        fe_problem.print_objectives();
        std::cerr << "-------------------------------\n";
      }

    }
  }

  if (0)
  {
    COMISO::TrustregionNewtonPCG trn(tol);
    trn.solve(&npp2, linear_constraints);
  }

//  // ALM
//  COMISO::AugmentedLagrangianMethod alm(10.0, 1e-1, 1e-6, 7);
//  alm.solve_experimental(&npp2, constraints, linear_constraints);

//  // Newton infeasible start
//  COMISO::NewtonSolver ns;
//  ns.solve_infeasible_start(&npp2, linear_constraints);
////  ns.solve_infeasible_start(&npp2, constraints_all);

  // IPOPT
//  try {
//    COMISO::IPOPTSolver ipsol;
////    ipsol.set_ipopt_option("print_level", 0);
////    ipsol.set_max_iterations(100);
////    ipsol.set_ipopt_option("tol", 1e-4);
////    ipsol.set_ipopt_option("hessian_approximation", "limited-memory");
////    ipsol.set_ipopt_option("limited_memory_max_history", 20);
//    ipsol.solve(&npp2, constraints_all);
//  }
//  catch (...)
//  {
//    std::cerr << "ERROR: IPOPT threw an exception!!!!!" << std::endl;
//  }

  // re-solve with higher weights on integrability
  if (0)
  {
    // boost integrability weight
    for (unsigned int i = 0; i < fe_integrability.instances().size(); ++i)
      fe_integrability.instances().c(i)[7] *= 100.0;

    try
    {
      COMISO::IPOPTSolver ipsol;
//    ipsol.set_ipopt_option("print_level", 0);
      ipsol.set_max_iterations(100);
//    ipsol.set_ipopt_option("tol", 1e-4);
//    ipsol.set_ipopt_option("hessian_approximation", "limited-memory");
//    ipsol.set_ipopt_option("limited_memory_max_history", 20);
      ipsol.solve(&npp2, constraints_all);
    }
    catch (...)
    {
      std::cerr << "ERROR: IPOPT threw an exception!!!!!" << std::endl;
    }
  }

//    COMISO::AugmentedLagrangianMethod alm(1e3, 1e-1, 1e-2, 7);
//    COMISO::NPTiming npp_alm3(&fe_problem_alm3); // with timings
//    fe_problem_alm3.x() = x0;
//    alm.solve_experimental(&npp_alm3, constraints, linear_constraints);
//    std::cerr << "******** finished via TinyAD derivatives ***********" << std::endl;
//    fe_problem.x() = fe_problem_alm3.x();

//    // final high-precision solve
//    if(0)
//    {
//      // make integrability constraints linearly independent
//      COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraints_all);
//
//      COMISO::NewtonSolver ns(1e-6, 1e-9, 2000);
////  ns.set_prescribed_constraint_regularization(-1e-9);
////  ns.set_prescribed_hessian_regularization(1e-9);
////  ns.solve_infeasible_start(&fe_problem, linear_constraints);
//      ns.solve_infeasible_start(&npp_alm3, constraints_all);
//
//      fe_problem.x() = fe_problem_alm3.x();
//    }
//  }

  // print energies at end
  {
    std::cerr << "--- final energies complete --- \n";
    fe_problem_complete.x() = fe_problem.x();
    fe_problem_complete.print_objectives();

    std::cerr << "--- final energies selected --- \n";
    fe_problem.print_objectives();
    std::cerr << "-------------------------------\n";
  }

  // export solution to frame field
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    int idx = 9 * c_it->idx();
    frame_cprop_[*c_it] = (Eigen::Map<Mat3d>((double *) &(fe_problem.x()[idx])).inverse()).transpose();
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
double
FrameFieldOptimizer3DT<TetMeshT>::
check_local_meshability()
{
//  ScopedStopWatch sw_total(sw::frame_field_int_opt);

  std::cout << "##### Check local meshability ..." << std::endl;
  int n(0.0);
  int n_meshable(0.0);

  int n_on_feature_vertex(0);
  int n_on_feature_vertex_meshable(0);
  int n_on_feature_edge(0);
  int n_on_feature_edge_meshable(0);
  int n_on_feature_face(0);
  int n_on_feature_face_meshable(0);

  int n_on_singular_node(0);
  int n_on_singular_node_meshable(0);
  int n_on_singular_arc(0);
  int n_on_singular_arc_meshable(0);

  int n_boundary(0);
  int n_boundary_meshable(0);
  int n_interior(0);
  int n_interior_meshable(0);

  std::vector<int> non_meshable_idxs;

  // check all vertices
  for (VIt v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
  {
    // global counter
    ++n;
    bool ilm = is_locally_meshable(*v_it);
    if (ilm)
      ++n_meshable;
    else
      non_meshable_idxs.push_back(v_it->idx());

    // count on feature vertex
    if (is_on_feature_vertex(*v_it))
    {
      ++n_on_feature_vertex;
      if (ilm)
        ++n_on_feature_vertex_meshable;
    }

    // count on feature edge
    if (is_on_feature_edge(*v_it))
    {
      ++n_on_feature_edge;
      if (ilm)
        ++n_on_feature_edge_meshable;
    }

    // count on feature face
    if (is_on_feature_face(*v_it))
    {
      ++n_on_feature_face;
      if (ilm)
        ++n_on_feature_face_meshable;
    }

    // count on singular node
    if (is_on_singular_node(*v_it))
    {
      ++n_on_singular_node;
      if (ilm)
        ++n_on_singular_node_meshable;
    }

    // count on singular arc
    if (is_on_singular_arc(*v_it))
    {
      ++n_on_singular_arc;
      if (ilm)
        ++n_on_singular_arc_meshable;
    }

    // count on singular arc
    if (mesh_.is_boundary(*v_it))
    {
      ++n_boundary;
      if (ilm)
        ++n_boundary_meshable;
    }
    else
    {
      ++n_interior;
      if (ilm)
        ++n_interior_meshable;
    }
  }

  // avoid NAN
  double feature_vertex_ratio = double(n_on_feature_vertex_meshable) / double(n_on_feature_vertex);
  if (!std::isfinite(feature_vertex_ratio))
    feature_vertex_ratio = 1.0;

  std::cerr << "#vertices          = " << std::setw(5) << n << ", #meshable = " << std::setw(5) << n_meshable << " ("
            << double(n_meshable) / double(n) * 100.0 << "%)\n";
  std::cerr << "#on feature vertex = " << std::setw(5) << n_on_feature_vertex << ", #meshable = " << std::setw(5)
            << n_on_feature_vertex_meshable << " (" << feature_vertex_ratio * 100.0 << "%)\n";
  std::cerr << "#on feature edge   = " << std::setw(5) << n_on_feature_edge << ", #meshable = " << std::setw(5)
            << n_on_feature_edge_meshable << " ("
            << double(n_on_feature_edge_meshable) / double(n_on_feature_edge) * 100.0 << "%)\n";
  std::cerr << "#on feature face   = " << std::setw(5) << n_on_feature_face << ", #meshable = " << std::setw(5)
            << n_on_feature_face_meshable << " ("
            << double(n_on_feature_face_meshable) / double(n_on_feature_face) * 100.0 << "%)\n";
  std::cerr << "#on singular node  = " << std::setw(5) << n_on_singular_node << ", #meshable = " << std::setw(5)
            << n_on_singular_node_meshable << " ("
            << double(n_on_singular_node_meshable) / double(n_on_singular_node) * 100.0 << "%)\n";
  std::cerr << "#on singular arc   = " << std::setw(5) << n_on_singular_arc << ", #meshable = " << std::setw(5)
            << n_on_singular_arc_meshable << " ("
            << double(n_on_singular_arc_meshable) / double(n_on_singular_arc) * 100.0 << "%)\n";
  std::cerr << "#boundary          = " << std::setw(5) << n_boundary << ", #meshable = " << std::setw(5)
            << n_boundary_meshable << " (" << double(n_boundary_meshable) / double(n_boundary) * 100.0 << "%)\n";
  std::cerr << "#interior          = " << std::setw(5) << n_interior << ", #meshable = " << std::setw(5)
            << n_interior_meshable << " (" << double(n_interior_meshable) / double(n_interior) * 100.0 << "%)\n";

  std::cerr << "non-meshable idx: ";
  for (auto it = non_meshable_idxs.begin(); it != non_meshable_idxs.end(); ++it)
    std::cerr << *it << ", ";
  std::cerr << std::endl;

  std::cout << "##### END Check local meshability ..." << std::endl;


  json_data_["percentage_meshable_vertices"] = double(n_meshable) / double(n);

  // return percentage of meshable vertices
  return double(n_meshable) / double(n);
}

//-----------------------------------------------------------------------------

template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
is_on_feature_vertex(const VH _vh) const
{
  // is feature vertex?
  return feature_vprop_[_vh] > 0;
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
is_on_feature_edge(const VH _vh) const
{
  // is on feature edge?
  VOHEIt vhe_it(_vh, &mesh_);
  for (; vhe_it.valid(); ++vhe_it)
    if (feature_eprop_[mesh_.edge_handle(*vhe_it)] > 0)
      return true;

  return false;
}

//-----------------------------------------------------------------------------

template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
is_on_feature_face(const VH _vh) const
{
  // is on feature triangle?
  VFIt vf_it(_vh, &mesh_);
  for (; vf_it.valid(); ++vf_it)
    if (feature_fprop_[*vf_it] > 0)
      return true;

  return false;
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
is_on_singular_node(const VH _vh) const
{
  std::vector<int> sv;

  VOHEIt vhe_it(_vh, &mesh_);
  for (; vhe_it.valid(); ++vhe_it)
  {
    int val = valence_eprop_[mesh_.edge_handle(*vhe_it)];
    if (val != 0)
      sv.push_back(val);
  }

  if (sv.empty()) return false;
  else if (sv.size() != 2) return true;
  else if (sv.size() == 2)
    return (sv[0] != sv[1]);

  return false;
}

//-----------------------------------------------------------------------------

template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
is_on_singular_arc(const VH _vh) const
{
  std::vector<int> sv;

  VOHEIt vhe_it(_vh, &mesh_);
  for (; vhe_it.valid(); ++vhe_it)
  {
    int val = valence_eprop_[mesh_.edge_handle(*vhe_it)];
    if (val != 0)
      sv.push_back(val);
  }

  if (sv.size() != 2) return false;
  else return (sv[0] == sv[1]);
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
has_non_feature_boundary_face(const EH _eh) const
{
  HEHFIt hehf_it(mesh_.halfedge_handle(_eh, 0), &mesh_);
  for (; hehf_it.valid(); ++hehf_it)
  {
    FH fh = mesh_.face_handle(*hehf_it);
    if (mesh_.is_boundary(fh))
      if (feature_fprop_[fh] == 0)
        return true;
  }

  return false;
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
is_locally_meshable(const VH _vh, const bool _verbose) const
{
  TetMesh tetmesh;
  bool valid = is_locally_meshable(_vh, tetmesh, _verbose);

  if ((!valid || always_export_local_configuration_) && save_locally_non_meshable_)
  {
    // export local configuration
    std::stringstream ss;
    ss << save_locally_non_meshable_filename_base_ << "_nlm_vertex_" << _vh.idx() << ".ovm";

    // write file
    OpenVolumeMesh::IO::FileManager fm;
    fm.writeFile(ss.str(), tetmesh);
  }

  return valid;
}

//-----------------------------------------------------------------------------

template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
is_locally_meshable(const VH _vh, TetMesh &_export_tmesh, const bool _verbose) const
{
  if (_verbose) std::cerr << "check vertex " << _vh.idx() << std::endl;

  // Todo: why it's not working in Hexmeshing but working in standalone meshability checker?
//    // output level of COMISO (silent=0)
  Debug::ScopedOutputLevel sol(int(_verbose) * 10);
//    Debug::ScopedOutputLevel sol(0);

  // ############# 1. get local mesh
  std::set<VH> vhs;
  std::set<EH> ehs;
  std::set<FH> fhs;
  std::map<CH, int> chm;
  std::vector<CH> chv;

  // collect cells and their sub-elements
  VCIt vc_it(_vh, &mesh_);
  for (; vc_it.valid(); ++vc_it)
  {
    // collect cells
    chm[*vc_it] = chv.size();
    chv.push_back(*vc_it);

    // collect vertices
    CVIt cv_it(*vc_it, &mesh_);
    for (; cv_it.valid(); ++cv_it)
      vhs.insert(*cv_it);

    // collect edges
    CEIt ce_it(*vc_it, &mesh_);
    for (; ce_it.valid(); ++ce_it)
      ehs.insert(*ce_it);

    // collect faces
    CFIt cf_it(*vc_it, &mesh_);
    for (; cf_it.valid(); ++cf_it)
      fhs.insert(*cf_it);
  }

  if (_verbose)
  {
    std::cerr << "#one-ring vertices = " << vhs.size() << std::endl;
    std::cerr << "#one-ring edges    = " << ehs.size() << std::endl;
    std::cerr << "#one-ring faces    = " << fhs.size() << std::endl;
    std::cerr << "#one-ring cells    = " << chv.size() << std::endl;
  }

  // get volume
  double VM(0.0);
  for (size_t i = 0; i < chv.size(); ++i)
  {
    double vol = volume(chv[i]);

    if (vol < 0.0)
      std::cerr << "Warning: cell has invalid volume of " << vol << std::endl;

    VM += vol;
  }
  // cubical root of volume
  double VM_cbrt = std::cbrt(VM);

  // ############ 2. setup problem

  const double w_amips = 0.0;
  const double w_log_det = 0.0;
  const double w_frame_fitting = 0.0;
  const double alpha_fitting = 1.0;
  // dihedral angle terms
  const double w_dihedral = 0.0;
  const double s_dihedral = 4.0;
  const double w_symmetric_dirichlet = 1.0;
  const double w_frame_smoothness = 1.0;

  // variables are 3x3 frames per cell
  int nv = 9 * chv.size();

  if (nv == 0) return false;

  // create problem
  COMISO::FiniteElementProblem fe_problem_complete(nv);
  COMISO::FiniteElementProblem fe_problem(nv);

  // set gradient field as initial solution (gradient field is inverse-transpose of frame field)
  for (size_t i = 0; i < chv.size(); ++i)
  {
    CH ch = chv[i];
    int idx = 9 * i;
    Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[idx]));
    G = (frame_cprop_[ch].inverse()).transpose();
  }

  // objective terms
  COMISO::FiniteElementSet <AMIPSFrameElement3D_PH> fe_amips("AMIPS");
  COMISO::FiniteElementSet <AMIPSFrameElement3D_TAD_PH> fe_amips_tad_ph("AMIPS_TinyAD");
  COMISO::FiniteElementSet <LogDetElement3D_PH> fe_logdet("LogDet");
  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_frame_fit("FrameFit");
  COMISO::FiniteElementSet <DihedralAngleFramesElement_TAD_PH> fe_dihedral("DihedralAngle");
  COMISO::FiniteElementSet <SymmetricDirichletFrameElement3D_TAD_PH> fe_sdirichlet("SymmetricDirichlet");

  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_face_align("FaceAlignment");
  COMISO::FiniteElementSet <LinearLeastSquaresElement3D> fe_edge_align("EdgeAlignment");
  COMISO::FiniteElementSet <LinearLeastSquaresElement6D> fe_integrability("Integrability");
  COMISO::FiniteElementSet <LinearLeastSquaresElement2D> fe_frame_smoothness("Frame Smoothness");


  // constraints
  std::vector<COMISO::LinearConstraint> c_face_align;
  std::vector<COMISO::LinearConstraint> c_edge_align;
  std::vector<COMISO::LinearConstraint> c_integrability;

  // setup objective terms amips etc.
  for (size_t i = 0; i < chv.size(); ++i)
  {
    CH ch = chv[i];
    int idx = 9 * i;

    // calc volume of tet
    double V = std::abs(volume(ch));
    double F_det = frame_cprop_[ch].determinant();
    // metric scaling for element terms representing A=diag(as,as,as) and the objective is that A*J = R \in SO(3)
    double as = std::cbrt(F_det);

    const double offset = 0.001;
    const double scale = F_det;

    // log determinant terms
    LogDetElement3D_PH::VecI idxs0;
    LogDetElement3D_PH::VecC c0;
    idxs0 << idx, idx + 1, idx + 2, idx + 3, idx + 4, idx + 5, idx + 6, idx + 7, idx + 8;
    c0 << w_log_det, scale, offset; // weight, scale, offset
    fe_logdet.instances().add_element(idxs0, c0);

    // frame fitting terms
    std::vector<LinearLeastSquaresElement3D::VecI> idxs1;
    std::vector<LinearLeastSquaresElement3D::VecC> c1;

    double w1 = (w_frame_fitting * V / VM) / (as * as);
    // get 9 instances of LinearElement
    LinearLeastSquaresElement3D::get_constants_and_indices_frame_fitting(
            Eigen::Map<Mat3d>((double *) &(fe_problem.x()[idx])),
            i,
            w1,
            alpha_fitting,
            idxs1,
            c1);
    for (size_t j = 0; j < idxs1.size(); ++j)
      fe_frame_fit.instances().add_element(idxs1[j], c1[j]);

    // amips terms
    AMIPSFrameElement3D_PH::VecI idxs2;
    AMIPSFrameElement3D_PH::VecC c2;
    // get indices corresponding to Jacobi matrix
    idxs2 << idx, idx + 3, idx + 6, idx + 1, idx + 4, idx + 7, idx + 2, idx + 5, idx + 8;
    // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], det(A) = _c[9], weight w=_c[10], distortion parameters s=_c[11], alpha=_c[12] and beta=_c[13]
    c2 << as, 0, 0, 0, as, 0, 0, 0, as, as * as * as, w_amips * V / VM, 1.0, 1.0, 0.5;
    fe_amips.instances().add_element(idxs2, c2);
    fe_amips_tad_ph.instances().add_element(idxs2, c2);

    // symmetric dirichlet terms
    SymmetricDirichletFrameElement3D_TAD_PH::VecC c3;
    c3 << as, 0, 0, 0, as, 0, 0, 0, as, w_symmetric_dirichlet * V / VM;
    fe_sdirichlet.instances().add_element(idxs2, c3);
  }

  //setup dihedral angle terms
  if (w_dihedral > 0.0)
  {
    // iterate over all edges incident to _vh
    VOHEIt vohe_it(_vh, &mesh_);
    for (; vohe_it.valid(); ++vohe_it)
    {
      EH eh = mesh_.edge_handle(*vohe_it);

//      // skip non-singular edges
//      if(valence_eprop_[eh] == 0)
//        continue;

      // skip boundary edges that are not on feature faces
      if (mesh_.is_boundary(eh) && has_non_feature_boundary_face(eh))
        continue;

      // determine target angle
      double target_angle = double(4 + valence_eprop_[eh]) * 0.5 * M_PI;
      if (mesh_.is_boundary(eh))
      {
        target_angle -= M_PI;
      }

      // get total interior angle
      double int_angle = 2.0 * M_PI;

      if (mesh_.is_boundary(eh))
      {
        int_angle = dihedral_angle(mesh_, mesh_.halfedge_handle(eh, 0));
        std::cerr << "boundary edge with angle " << int_angle * 180.0 / M_PI << " and target angle "
                  << target_angle * 180.0 / M_PI << std::endl;
      }
      // determine uniform scaling of input angles to reach target angle
      double scale_angle = target_angle / int_angle;

      HEH heh = *vohe_it;
      HEH heh_opp = mesh_.opposite_halfedge_handle(heh);
      VH vh0 = mesh_.halfedge(heh).from_vertex();
      VH vh1 = mesh_.halfedge(heh).to_vertex();
      Vec3d pt0 = ovm2eigen(mesh_.vertex(vh0));
      Vec3d pt1 = ovm2eigen(mesh_.vertex(vh1));

      double angle(0.0);
      HEHFIt hehf_it(*vohe_it, &mesh_);
      for (; hehf_it.valid(); ++hehf_it)
        if (!mesh_.is_boundary(*hehf_it)) // skip boundary halffaces
        {
          // get halfface
          HFH hfh = *hehf_it;
          // get cell
          CH ch = mesh_.incident_cell(hfh);
          // get other halfface in cell
          HFH hf_adj = mesh_.adjacent_halfface_in_cell(hfh, heh);

          // get vertices and points of both halffaces
          HEH he0n = mesh_.next_halfedge_in_halfface(heh, hfh);
          HEH he1n = mesh_.next_halfedge_in_halfface(heh_opp, hf_adj);
          if (!he0n.is_valid() || !he1n.is_valid())
            std::cerr << "ERROR: next_halfedge_in_halfface is invalid while setting up dihedral terms" << std::endl;

          VH vh00 = mesh_.halfedge(he0n).to_vertex();
          VH vh11 = mesh_.halfedge(he1n).to_vertex();
          Vec3d pt00 = ovm2eigen(mesh_.vertex(vh00));
          Vec3d pt11 = ovm2eigen(mesh_.vertex(vh11));

          // get edge vectors
          Vec3d e0 = pt1 - pt0;
          Vec3d e1 = pt00 - pt0;
          Vec3d e2 = pt11 - pt0;

          // determine target angle
          Vec3d n0 = e1.cross(e0);
          n0 /= n0.norm();
          Vec3d n1 = e2.cross(e0);
          n1 /= n1.norm();
          double dp = n0.dot(n1);
          dp = std::max(-1.0, dp);
          dp = std::min(1.0, dp);
          double alpha = scale_angle * std::acos(dp);

          // get frame
          if (chm.find(ch) == chm.end())
            std::cerr << "ERROR: required cell is not part of local neighborhood" << std::endl;
          // get base index
          int idx = 9 * chm[ch];

          DihedralAngleFramesElement_TAD_PH::VecI idxs_dihedral;
          DihedralAngleFramesElement_TAD_PH::VecC c_dihedral;
          idxs_dihedral << idx, idx + 3, idx + 6, idx + 1, idx + 4, idx + 7, idx + 2, idx + 5, idx + 8;
          c_dihedral << e0[0], e0[1], e0[2], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], w_dihedral, s_dihedral, std::cos(
                  alpha);
          fe_dihedral.instances().add_element(idxs_dihedral, c_dihedral);

          // debug
          if (0)
          {
            DihedralAngleFramesElement_TAD_PH de;
            DihedralAngleFramesElement_TAD_PH::VecV x_dihedral;
            x_dihedral << fe_problem.x()[idxs_dihedral[0]],
                    fe_problem.x()[idxs_dihedral[1]],
                    fe_problem.x()[idxs_dihedral[2]],
                    fe_problem.x()[idxs_dihedral[3]],
                    fe_problem.x()[idxs_dihedral[4]],
                    fe_problem.x()[idxs_dihedral[5]],
                    fe_problem.x()[idxs_dihedral[6]],
                    fe_problem.x()[idxs_dihedral[7]],
                    fe_problem.x()[idxs_dihedral[8]];

            double fx = de.eval_f(x_dihedral, c_dihedral);

            std::cerr << "fx = " << fx << ", valence=" << valence_eprop_[eh] << " --> target angle = "
                      << alpha / M_PI * 180.0 << " degree, e0^T = " << e0.transpose() << ", e1^T = " << e1.transpose()
                      << ", e2^T = " << e2.transpose() << ", volume = " << volume(ch) << std::endl;

          }
        }
    }
  }

  // setup alignment terms
  std::map<int, int> cell_face_alignment_constraint;

  // (i) face alignment
  for (auto f_it = fhs.begin(); f_it != fhs.end(); ++f_it)
    if (feature_fprop_[*f_it] > 0)
    {
      // get non-boundary halffaces
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);
      std::vector<HFH> hfhs;
      if (!mesh_.is_boundary(hfh0) &&
          chm.find(mesh_.incident_cell(hfh0)) != chm.end())
        hfhs.push_back(hfh0);

      if (!mesh_.is_boundary(hfh1) &&
          chm.find(mesh_.incident_cell(hfh1)) != chm.end())
        hfhs.push_back(hfh1);


      for (size_t i = 0; i < hfhs.size(); ++i)
      {
        HFH hfh = hfhs[i];
        // get corresponding cell
        CH ch = mesh_.incident_cell(hfh);

        // base index of frame
        int chidx = chm[ch];
        int bidx = 9 * chidx;

        // get vertices
        VH vh0 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[0]).to_vertex();
        VH vh1 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[1]).to_vertex();
        VH vh2 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[2]).to_vertex();

        // skip constraints that are not incident to checked vertex
        if (vh0 != _vh && vh1 != _vh && vh2 != _vh)
          continue;

        // get points
        Point p0 = mesh_.vertex(vh0);
        Point p1 = mesh_.vertex(vh1);
        Point p2 = mesh_.vertex(vh2);

        // get outward normal vector
        Vec3d n = ovm2eigen((p1 - p0) % (p2 - p0));
        Vec3d e0 = ovm2eigen(p1 - p0);
        Vec3d e1 = ovm2eigen(p2 - p0);

        double A = 0.5 * n.norm();

        // get local basis of triangle
        Vec3d u = ovm2eigen(p1 - p0);
        u /= u.norm();
        Vec3d v = n.cross(u);
        v /= v.norm();

        // get current gradient vectors
        Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[bidx]));

        // determine alignment axis
        Vec3d Je0 = G.transpose() * e0;
        Vec3d Je1 = G.transpose() * e1;
        Vec3d Jn = Je0.cross(Je1);
        AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Jn);
        int axis_idx = int(e_axis) / 2;
        // test whether field is sufficiently aligned
        Vec3d a = AxisAlignmentHelpers::vector(e_axis);
        if (a.dot(Jn) < Jn.norm() * 0.99985)
          std::cerr << "Warning FF alignment: face " << f_it->idx() << " alignment axis deviates by "
                    << std::acos(a.dot(Jn) / Jn.norm()) * 180.0 / M_PI
                    << " degree\n";

        // memorize axis
        cell_face_alignment_constraint[chidx] = axis_idx;

        double w = A;

        int bidx2 = bidx + 3 * axis_idx;
        LinearLeastSquaresElement3D::VecI idx(bidx2, bidx2 + 1, bidx2 + 2);
        LinearLeastSquaresElement3D::VecC c0, c1;
        c0 << u[0], u[1], u[2], 0.0, w;
        c1 << v[0], v[1], v[2], 0.0, w;
        fe_face_align.instances().add_element(idx, c0);
        fe_face_align.instances().add_element(idx, c1);

        // generate corresponding constraints
        double w2 = std::sqrt(A);
        w2 = 1.0;
        COMISO::LinearConstraint::SVectorNC lc0(nv);
        lc0.coeffRef(idx[0]) = w2 * c0[0];
        lc0.coeffRef(idx[1]) = w2 * c0[1];
        lc0.coeffRef(idx[2]) = w2 * c0[2];
        c_face_align.push_back(COMISO::LinearConstraint(lc0, 0.0));

        COMISO::LinearConstraint::SVectorNC lc1(nv);
        lc1.coeffRef(idx[0]) = w2 * c1[0];
        lc1.coeffRef(idx[1]) = w2 * c1[1];
        lc1.coeffRef(idx[2]) = w2 * c1[2];
        c_face_align.push_back(COMISO::LinearConstraint(lc1, 0.0));
      }
    }

  // (ii) edge alignment
  for (auto e_it = ehs.begin(); e_it != ehs.end(); ++e_it)
  {
    // singular or feature edge?
    if (valence_eprop_[*e_it] != 0 || feature_eprop_[*e_it] > 0)
    {
      // get adjacent vertices
      VH vh0 = mesh_.edge(*e_it).from_vertex();
      VH vh1 = mesh_.edge(*e_it).to_vertex();

      // skip constraints that are not incident to checked vertex
      if (vh0 != _vh && vh1 != _vh)
        continue;

      // get all cells incident to the edge
      std::vector<CH> chs;
      for (HECIt hec_it = mesh_.hec_iter(mesh_.halfedge_handle(*e_it, 0)); hec_it.valid(); ++hec_it)
        if (chm.find(*hec_it) != chm.end()) // only consider cells in local neighborhood
          chs.push_back(*hec_it);


      for (size_t i = 0; i < chs.size(); ++i)
      {
        // get curent cell adjacent to edge
        CH ch = chs[i];

        // base index of frame
        int chidx = chm[ch];
        int bidx = 9 * chidx;

        // get edge vector
        Vec3d e = ovm2eigen(mesh_.vertex(vh1) - mesh_.vertex(vh0));

        double l = e.norm();

        // get current gradient vectors
        Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[bidx]));

        // determine alignment axis
        Vec3d Je = G.transpose() * e;
        AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Je);
        int axis_idx = int(e_axis) / 2;
        // test whether field is sufficiently aligned
        Vec3d a = AxisAlignmentHelpers::vector(e_axis);
        if (a.dot(Je) < Je.norm() * 0.99985)
          std::cerr << "Warning FF alignment: edge " << e_it->idx() << " alignment axis deviates by "
                    << std::acos(a.dot(Je) / Je.norm()) * 180.0 / M_PI
                    << " degree\n";

        Je[(axis_idx + 1) % 3] = 0.0;
        Je[(axis_idx + 2) % 3] = 0.0;

        double w = VM_cbrt / l;
        double w2 = std::sqrt(VM_cbrt / l);
        // hack
        w2 = 1.0 / e.norm();

        // check if face alignment exists for this cell
        int face_alignment_axis = -1;
        auto cfc = cell_face_alignment_constraint.find(chidx);
        if (cfc != cell_face_alignment_constraint.end())
          face_alignment_axis = cfc->second; // ---> skip constraint since linearly dependent

        // the two other axes need to be orthogonal!
        // hack: constrain J*e = axis (not active at the moment)
        for (int j = 0; j < 3; ++j)
          if (j != axis_idx && j != face_alignment_axis)
          {
            int bidx2 = bidx + 3 * j;
            LinearLeastSquaresElement3D::VecI idx(bidx2, bidx2 + 1, bidx2 + 2);
            LinearLeastSquaresElement3D::VecC c;
            c << e[0], e[1], e[2], -Je[j], w;
            fe_edge_align.instances().add_element(idx, c);

            // generate corresponding constraints
            COMISO::LinearConstraint::SVectorNC lc(nv);
            lc.coeffRef(idx[0]) = w2 * c[0];
            lc.coeffRef(idx[1]) = w2 * c[1];
            lc.coeffRef(idx[2]) = w2 * c[2];
            c_edge_align.push_back(COMISO::LinearConstraint(lc, -w2 * Je[j]));
          }
      }
    }
  }

  // setup integrability terms
  for (auto f_it = fhs.begin(); f_it != fhs.end(); ++f_it)
    if (!mesh_.is_boundary(*f_it))
    {
      // get non-boundary halffaces
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);

      // get corresponding cells
      CH ch0 = mesh_.incident_cell(hfh0);
      CH ch1 = mesh_.incident_cell(hfh1);

      // need both in local neighborhood
      if (chm.find(ch0) != chm.end() && chm.find(ch1) != chm.end())
      {
        int ch0idx = chm[ch0];
        int ch1idx = chm[ch1];

        // calc volumes of tets
        double V0 = std::abs(volume(ch0));
        double V1 = std::abs(volume(ch1));
        double VF = (V0 + V1) * 0.25;

        double F0_det = frame_cprop_[ch0].determinant();
        double F1_det = frame_cprop_[ch1].determinant();
        double as = std::cbrt(0.5 * (F0_det + F1_det));

        // get local basis of face
        // get vertices
        VH vh0 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[0]).to_vertex();
        VH vh1 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[1]).to_vertex();
        VH vh2 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[2]).to_vertex();

        // get points
        Point p0 = mesh_.vertex(vh0);
        Point p1 = mesh_.vertex(vh1);
        Point p2 = mesh_.vertex(vh2);

        // get outward normal vector
        Vec3d n = ovm2eigen((p1 - p0) % (p2 - p0));
        // get local basis of triangle
        Vec3d u = ovm2eigen(p1 - p0);
        u /= u.norm();
        Vec3d v = n.cross(u);
        v /= v.norm();

        // get transition from ch0 to ch1 such that G0*T01 is w.r.t. coords of G1 (transition for gradients is inverse-transpose of the frame version)
        Mat3d T01 = (tq_.transition_matrix_int(
                tq_.inverse_transition_idx(transition_hfprop_[hfh0]))).transpose();

        // get permutation and signs of transformed gradients
        Vec3i k;
        Vec3d s;

        if (T01(0, 0) != 0)
        {
          k[0] = 0;
          s[0] = T01(0, 0);
        }
        else if (T01(1, 0) != 0)
        {
          k[0] = 1;
          s[0] = T01(1, 0);
        }
        else if (T01(2, 0) != 0)
        {
          k[0] = 2;
          s[0] = T01(2, 0);
        }

        if (T01(0, 1) != 0)
        {
          k[1] = 0;
          s[1] = T01(0, 1);
        }
        else if (T01(1, 1) != 0)
        {
          k[1] = 1;
          s[1] = T01(1, 1);
        }
        else if (T01(2, 1) != 0)
        {
          k[1] = 2;
          s[1] = T01(2, 1);
        }

        if (T01(0, 2) != 0)
        {
          k[2] = 0;
          s[2] = T01(0, 2);
        }
        else if (T01(1, 2) != 0)
        {
          k[2] = 1;
          s[2] = T01(1, 2);
        }
        else if (T01(2, 2) != 0)
        {
          k[2] = 2;
          s[2] = T01(2, 2);
        }

        double w = VF / VM_cbrt;

        double w2 = std::sqrt(VF / VM_cbrt);

        // hack
        w2 = 1.0;

        // setup constraints for all axes
        for (unsigned int i = 0; i < 3; ++i)
        {
          // get base indices
          int bidx0 = 9 * ch0idx + k[i] * 3;
          int bidx1 = 9 * ch1idx + i * 3;

          LinearLeastSquaresElement6D::VecI idxs;
          idxs << bidx0, bidx0 + 1, bidx0 + 2, bidx1, bidx1 + 1, bidx1 + 2;

          LinearLeastSquaresElement6D::VecC c0, c1;
          c0 << s[i] * u[0], s[i] * u[1], s[i] * u[2], -u[0], -u[1], -u[2], 0.0, w;
          c1 << s[i] * v[0], s[i] * v[1], s[i] * v[2], -v[0], -v[1], -v[2], 0.0, w;

          fe_integrability.instances().add_element(idxs, c0);
          fe_integrability.instances().add_element(idxs, c1);

          // generate corresponding constraints
          COMISO::LinearConstraint::SVectorNC lc0(nv);
          lc0.coeffRef(idxs[0]) = w2 * c0[0];
          lc0.coeffRef(idxs[1]) = w2 * c0[1];
          lc0.coeffRef(idxs[2]) = w2 * c0[2];
          lc0.coeffRef(idxs[3]) = w2 * c0[3];
          lc0.coeffRef(idxs[4]) = w2 * c0[4];
          lc0.coeffRef(idxs[5]) = w2 * c0[5];
          c_integrability.push_back(COMISO::LinearConstraint(lc0, 0.0));
          COMISO::LinearConstraint::SVectorNC lc1(nv);
          lc1.coeffRef(idxs[0]) = w2 * c1[0];
          lc1.coeffRef(idxs[1]) = w2 * c1[1];
          lc1.coeffRef(idxs[2]) = w2 * c1[2];
          lc1.coeffRef(idxs[3]) = w2 * c1[3];
          lc1.coeffRef(idxs[4]) = w2 * c1[4];
          lc1.coeffRef(idxs[5]) = w2 * c1[5];
          c_integrability.push_back(COMISO::LinearConstraint(lc1, 0.0));
        }

        // setup smoothness terms
        double w_smooth = (w_frame_smoothness * VF / VM) * (as * as);   // geometric weight

        // setup smoothness terms for all axes
        for (unsigned int i = 0; i < 3; ++i)
        {
          // get base indices
          int bidx0 = 9 * ch0idx + k[i] * 3;
          int bidx1 = 9 * ch1idx + i * 3;

          for (unsigned int j = 0; j < 3; ++j)
          {
            LinearLeastSquaresElement2D::VecI idxs;
            idxs << bidx0 + j, bidx1 + j;

            LinearLeastSquaresElement2D::VecC c;
            c << s[i], -1.0, 0.0, w_smooth;

            fe_frame_smoothness.instances().add_element(idxs, c);
          }
        }
      }
    }
  // add all terms to objective function
  if (_verbose)
  {
    fe_problem_complete.add_set(&fe_amips);
    fe_problem_complete.add_set(&fe_amips_tad_ph);
    fe_problem_complete.add_set(&fe_logdet);
    fe_problem_complete.add_set(&fe_frame_fit);
    fe_problem_complete.add_set(&fe_sdirichlet);
    fe_problem_complete.add_set(&fe_dihedral);
    fe_problem_complete.add_set(&fe_edge_align);
    fe_problem_complete.add_set(&fe_face_align);
    fe_problem_complete.add_set(&fe_integrability);
    fe_problem_complete.add_set(&fe_frame_smoothness);
  }

  // add all terms to objective function
  if (w_amips > 0.0)
    fe_problem.add_set(&fe_amips_tad_ph);
//    fe_problem.add_set(&fe_amips);
  if (w_log_det > 0.0)
    fe_problem.add_set(&fe_logdet);
  if (w_frame_fitting > 0.0)
    fe_problem.add_set(&fe_frame_fit);
  if (w_symmetric_dirichlet > 0.0)
    fe_problem.add_set(&fe_sdirichlet);
  if (w_dihedral > 0.0)
    fe_problem.add_set(&fe_dihedral);
  if (w_frame_smoothness > 0.0)
    fe_problem.add_set(&fe_frame_smoothness);

  if (_verbose)
  {
    std::cerr << "#edge alignment = " << c_edge_align.size() << std::endl;
    std::cerr << "#face alignment = " << c_face_align.size() << std::endl;
    std::cerr << "#integrability  = " << c_integrability.size() << std::endl;
  }

  std::vector<COMISO::LinearConstraint> linear_constraints;
  linear_constraints.reserve(c_edge_align.size() + c_face_align.size());
  for (size_t i = 0; i < c_edge_align.size(); ++i)
    linear_constraints.push_back(c_edge_align[i]);
  for (size_t i = 0; i < c_face_align.size(); ++i)
    linear_constraints.push_back(c_face_align[i]);

  // only for infeasible start Newton solver!!!
  for (size_t i = 0; i < c_integrability.size(); ++i)
    linear_constraints.push_back(c_integrability[i]);

  // collect constraints
  std::vector<COMISO::NConstraintInterface *> constraints;
  constraints.reserve(linear_constraints.size() + c_integrability.size());

//  for(size_t i=0; i<c_integrability.size(); ++i)
//    constraints.push_back( &(c_integrability[i]));

  for (size_t i = 0; i < linear_constraints.size(); ++i)
    constraints.push_back(&(linear_constraints[i]));

//  {
//    // check rank of constraints
//    SMatXd A;
//    VecXd  b;
//    COMISO::LinearConstraintConverter::nsolver_to_eigen(constraints, A, b, nv);
//
//    Eigen::SparseQR<SMatXd,Eigen::COLAMDOrdering<int> > sqr(A);
//    std::cerr << "rank(A) = " << sqr.rank() << std::endl;
//  }

  // make integrability constraints linearly independent
  COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraints);

  // ############# 3. optimize

  // print energies at beginning
  if (_verbose)
  {
    std::cerr << "--- initial energies complete --- \n";
    fe_problem_complete.x() = fe_problem.x();
    fe_problem_complete.print_objectives();

    std::cerr << "--- initial energies selected --- \n";
    fe_problem.print_objectives();
    std::cerr << "-------------------------------\n";
  }

  // store initial solution
  std::vector<double> x0 = fe_problem.x();

//  COMISO::AugmentedLagrangianMethod alm(1e3, 1e-6, 1e-6, 12);
//  alm.set_silent(true);
////  std::cerr << "#constraints: " << constraints.size() << std::endl;
////  std::cerr << "#linear constraints: " << linear_constraints.size() << std::endl;
//  alm.solve_experimental(&fe_problem, constraints, linear_constraints);

  COMISO::NewtonSolver ns(1e-6, 1e-9, 2000);
//  ns.set_prescribed_constraint_regularization(-1e-9);
//  ns.set_prescribed_hessian_regularization(1e-9);
//  ns.solve_infeasible_start(&fe_problem, linear_constraints);
  auto cv = ns.solve_infeasible_start(&fe_problem, constraints);
//  try {
//    COMISO::IPOPTSolver ipsol;
////    ipsol.set_ipopt_option("print_level", 0);
//    ipsol.set_max_iterations(5000);
//    ipsol.solve(&fe_problem, constraints);
//  }
//  catch (...)
//  {
//    std::cerr << "ERROR: IPOPT threw an exception!!!!!" << std::endl;
//  }

  if (_verbose)
  {
    std::cerr << "--- final energies complete ---\n";
    fe_problem_complete.x() = fe_problem.x();
    fe_problem_complete.print_objectives();

    std::cerr << "--- final energies selected ---\n";
    fe_problem.print_objectives();
    std::cerr << "-------------------------------\n";
  }
  // ############# 4. analyze result

  bool valid = true;

  // check local injectivity
  for (size_t i = 0; i < chv.size(); ++i)
  {
    int bidx = 9 * i;

    Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[bidx]));
    if (G.determinant() < 0.0)
    {
      std::cerr << "Warning: local neighborhood of vertex " << _vh.idx() << " is not locally injective: "
                << G.determinant() << std::endl;
      valid = false;
    }
  }

  // check valence of edges
  VOHEIt vohe_it(_vh, &mesh_);
  std::vector<EH> invalid_edges;
  for (; vohe_it.valid(); ++vohe_it)
  {
    EH eh = mesh_.edge_handle(*vohe_it);

    // skip boundary edges that are not on feature faces
    if (mesh_.is_boundary(eh) && has_non_feature_boundary_face(eh))
      continue;

    double target_angle = double(4 + valence_eprop_[eh]) * 0.5 * M_PI;
    if (mesh_.is_boundary(eh))
      target_angle -= M_PI;

//      std::cerr << "check valence of edge " << eh.idx() << " with target angle " << target_angle/M_PI << " Pi " << std::endl;

    HEH heh = *vohe_it;
    VH vh0 = mesh_.halfedge(heh).from_vertex();
    VH vh1 = mesh_.halfedge(heh).to_vertex();
    Vec3d pt0 = ovm2eigen(mesh_.vertex(vh0));
    Vec3d pt1 = ovm2eigen(mesh_.vertex(vh1));

    double angle(0.0);
    HEHFIt hehf_it(*vohe_it, &mesh_);
    for (; hehf_it.valid(); ++hehf_it)
      if (!mesh_.is_boundary(*hehf_it)) // skip boundary halffaces
      {
        HFH hfh = *hehf_it;
        CH ch = mesh_.incident_cell(hfh);
        // get frame
        if (chm.find(ch) == chm.end())
          std::cerr << "ERROR: required cell is not part of local neighborhood" << std::endl;

        int bidx = 9 * chm[ch];
        Eigen::Map <Mat3d> G((double *) &(fe_problem.x()[bidx])); // G=J^T and has gradients as column vectors
        Mat3d Gi = G.inverse();

        Vec3d n1 = ovm2eigen(mesh_.normal(hfh));
        auto hf_adj = mesh_.adjacent_halfface_in_cell(hfh, heh);
        Vec3d n0 = ovm2eigen(mesh_.normal(hf_adj));

        // transform with Jacobian
        Vec3d pt0t = G.transpose() * pt0;
        Vec3d pt1t = G.transpose() * pt1;
        // transform with inverse-transpose Jacobian
        Vec3d n0t = Gi * n0;
        n0t /= n0t.norm();
        Vec3d n1t = Gi * n1;
        n1t /= n1t.norm();

        angle += dihedral_angle(pt0t, pt1t, n0t, n1t);

//                std::cerr << "add angle: " << dihedral_angle(pt0t, pt1t, n0t, n1t) << std::endl;
      }
//        std::cerr << "total angle: " << angle << std::endl;

    if (std::abs(angle - target_angle) > 1e-4)
    {
      std::cerr << "at vertex " << _vh.idx() << " valence of edge " << eh.idx()
                << " was not reproduced correctly. angle = " << angle * 180 / M_PI << " vs. target angle = "
                << target_angle * 180 / M_PI << std::endl;
      invalid_edges.push_back(eh);
      valid = false;
    }
  }

  if (!valid || always_export_local_configuration_)
    export_local_configuration(_export_tmesh, chv, invalid_edges, x0, fe_problem.x());

  return valid;
}

//-----------------------------------------------------------------------------

template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::
export_local_configuration(TetMesh &tetmesh, const std::vector<CH> &_chv, const std::vector<EH> &_invalid_edges,
                           const std::vector<double> &_x0, const std::vector<double> &_x) const
{
  auto einvalid = tetmesh.template request_edge_property<double>("edge_invalid");
  tetmesh.set_persistent(einvalid, true);

  auto eval = tetmesh.template request_edge_property<int>("edge_valance");
  tetmesh.set_persistent(eval, true);

  auto vfeat = tetmesh.template request_vertex_property<int>("AlgoHex::FeatureVertices");
  tetmesh.set_persistent(vfeat, true);

  auto efeat = tetmesh.template request_edge_property<int>("AlgoHex::FeatureEdges");
  tetmesh.set_persistent(efeat, true);

  auto ffeat = tetmesh.template request_face_property<int>("AlgoHex::FeatureFaces");
  tetmesh.set_persistent(ffeat, true);

  auto hftrans = tetmesh.template request_halfface_property<int>("HalffaceTransiton");
  tetmesh.set_persistent(hftrans, true);

  auto cfu = tetmesh.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_u");
  tetmesh.set_persistent(cfu, true);

  auto cfv = tetmesh.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_v");
  tetmesh.set_persistent(cfv, true);

  auto cfw = tetmesh.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_w");
  tetmesh.set_persistent(cfw, true);

  auto cfun = tetmesh.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_new_u");
  tetmesh.set_persistent(cfun, true);

  auto cfvn = tetmesh.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_new_v");
  tetmesh.set_persistent(cfvn, true);

  auto cfwn = tetmesh.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_new_w");
  tetmesh.set_persistent(cfwn, true);

  std::map<VH, VH> vm;
  std::map<VH, VH> vm_inv;

  // generate mesh and set cell properties
  for (size_t i = 0; i < _chv.size(); ++i)
  {
    CH ch = _chv[i];
    std::vector<VH> vhs;
    TVIt tv_it(ch, &mesh_);
    for (; tv_it.valid(); ++tv_it)
    {
      auto m_it = vm.find(*tv_it);
      if (m_it == vm.end())
      {
        VH vh_new = tetmesh.add_vertex(mesh_.vertex(*tv_it));
        vm[*tv_it] = vh_new;
        vm_inv[vh_new] = *tv_it;
        vhs.push_back(vh_new);
      }
      else vhs.push_back(m_it->second);
    }

    CH ch_new = tetmesh.add_cell(vhs);

    // store initial frame
    cfu[ch_new] = eigen2ovm(frame_cprop_[ch].col(0));
    cfv[ch_new] = eigen2ovm(frame_cprop_[ch].col(1));
    cfw[ch_new] = eigen2ovm(frame_cprop_[ch].col(2));

    // store optimized frame
    int idx = 9 * i;
    Eigen::Map <Mat3d> G((double *) &(_x[idx]));
    Mat3d F = G.inverse().transpose();
    cfun[ch_new] = eigen2ovm(F.col(0));
    cfvn[ch_new] = eigen2ovm(F.col(1));
    cfwn[ch_new] = eigen2ovm(F.col(2));
  }

  // set vertex properties
  for (VIt v_it = tetmesh.vertices_begin(); v_it != tetmesh.vertices_end(); ++v_it)
  {
    VH vh_orig = vm_inv[*v_it];
    vfeat[*v_it] = feature_vprop_[vh_orig];
  }

  // set edge properties
  for (EIt e_it = tetmesh.edges_begin(); e_it != tetmesh.edges_end(); ++e_it)
  {
    // get original edge handle
    EH eh = *e_it;
    VH vh0 = tetmesh.halfedge(mesh_.halfedge_handle(eh, 0)).to_vertex();
    VH vh1 = tetmesh.halfedge(mesh_.halfedge_handle(eh, 1)).to_vertex();
    VH vho0 = vm_inv[vh0];
    VH vho1 = vm_inv[vh1];
    HEH heh = mesh_.find_halfedge(vho0, vho1);
    if (!heh.is_valid())
      std::cerr << "ERROR: could not obtain halfedge of original mesh in export_local_configuration!!!" << std::endl;
    EH eho = mesh_.edge_handle(heh);

    // set default color to gray
    einvalid[*e_it] = 0.0;
    // copy valence property
    eval[*e_it] = valence_eprop_[eho];
    efeat[*e_it] = feature_eprop_[eho];
  }

  // mark problematic edges
  for (size_t i = 0; i < _invalid_edges.size(); ++i)
  {
    EH eh = _invalid_edges[i];
    VH vh0 = mesh_.halfedge(mesh_.halfedge_handle(eh, 0)).to_vertex();
    VH vh1 = mesh_.halfedge(mesh_.halfedge_handle(eh, 1)).to_vertex();

    VH vhn0 = vm[vh0];
    VH vhn1 = vm[vh1];

    einvalid[tetmesh.edge_handle(tetmesh.find_halfedge(vhn0, vhn1))] = 1.0;
  }

  // set face properties
  for (FIt f_it = tetmesh.faces_begin(); f_it != tetmesh.faces_end(); ++f_it)
  {
    FH fh = *f_it;
    HFH hfh0 = tetmesh.halfface_handle(fh, 0);
    HFH hfh1 = tetmesh.halfface_handle(fh, 1);

    auto f0 = tetmesh.halfface(hfh0);
    std::vector<VH> vhs0;
    vhs0.push_back(tetmesh.halfedge(f0.halfedges()[0]).to_vertex());
    vhs0.push_back(tetmesh.halfedge(f0.halfedges()[1]).to_vertex());
    vhs0.push_back(tetmesh.halfedge(f0.halfedges()[2]).to_vertex());

    // map vertex indices
    std::vector<VH> vhs0_orig;
    vhs0_orig.push_back(vm_inv[vhs0[0]]);
    vhs0_orig.push_back(vm_inv[vhs0[1]]);
    vhs0_orig.push_back(vm_inv[vhs0[2]]);

    // get corresponding halfface in original mesh
    HFH hfh0_orig = mesh_.find_halfface(vhs0_orig);
    if (!hfh0_orig.is_valid())
      std::cerr << "ERROR: could not obtain original halfface in export_local_configuration!!!" << std::endl;
    HFH hfh1_orig = mesh_.opposite_halfface_handle(hfh0_orig);
    FH fh_orig = mesh_.face_handle(hfh0_orig);

    // copy properties
    ffeat[fh] = feature_fprop_[fh_orig];
    hftrans[hfh0] = transition_hfprop_[hfh0_orig];
    hftrans[hfh1] = transition_hfprop_[hfh1_orig];
  }
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
double
FrameFieldOptimizer3DT<TetMeshT>::
volume(const CH _ch) const
{
  Mat3x4d P;
  mesh_tet(_ch, P);
  Mat3d E;
  E.col(0) = P.col(1) - P.col(0);
  E.col(1) = P.col(2) - P.col(0);
  E.col(2) = P.col(3) - P.col(0);

  return (1.0 / 6.0 * E.determinant());
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
double
FrameFieldOptimizer3DT<TetMeshT>::
volume_frame(const CH _ch) const
{
  Mat3x4d P;
  mesh_tet(_ch, P);
  Mat3d E;
  E.col(0) = P.col(1) - P.col(0);
  E.col(1) = P.col(2) - P.col(0);
  E.col(2) = P.col(3) - P.col(0);

  Mat3d F = frame_cprop_[_ch];

  // we want det[ 1/6(F^-1 * E)] =det[E]/(6 det[F])
  return (E.determinant() / (6.0 * F.determinant()));

//  // deform tet with Jacobian induced by the frame
//  Mat3x3 J = frame_cprop_[*c_it].inverse();
//  Mat3x3 JE = J*E;
//
//  return( 1.0 / 6.0 * JE.determinant() );
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
double
FrameFieldOptimizer3DT<TetMeshT>::
mesh_volume() const
{
  double V(0.0);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    V += std::abs(volume(*c_it));
  }

  return V;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::
mesh_tet(const CH _ch, Mat3x4d &_P) const
{
  int i = 0;
  for (TVIt tv_it = mesh_.tv_iter(_ch); tv_it.valid(); ++tv_it)
  {
    _P.col(i++) = ovm2eigen(mesh_.vertex(*tv_it));
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::
import_frames_from_vec3d_properties()
{
  if (!mesh_.template cell_property_exists<OpenVolumeMesh::Vec3d>("frame_u"))
    std::cerr << "ERROR: import_frames_from_vec3d_properties was called but property frame_u does not exist"
              << std::endl;
  if (!mesh_.template cell_property_exists<OpenVolumeMesh::Vec3d>("frame_v"))
    std::cerr << "ERROR: import_frames_from_vec3d_properties was called but property frame_v does not exist"
              << std::endl;
  if (!mesh_.template cell_property_exists<OpenVolumeMesh::Vec3d>("frame_w"))
    std::cerr << "ERROR: import_frames_from_vec3d_properties was called but property frame_w does not exist"
              << std::endl;


  auto cfu = mesh_.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_u");
  auto cfv = mesh_.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_v");
  auto cfw = mesh_.template request_cell_property<OpenVolumeMesh::Vec3d>("frame_w");

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    OpenVolumeMesh::Vec3d u = cfu[*c_it];
    OpenVolumeMesh::Vec3d v = cfv[*c_it];
    OpenVolumeMesh::Vec3d w = cfw[*c_it];

    frame_cprop_[*c_it](0, 0) = u[0];
    frame_cprop_[*c_it](1, 0) = u[1];
    frame_cprop_[*c_it](2, 0) = u[2];

    frame_cprop_[*c_it](0, 1) = v[0];
    frame_cprop_[*c_it](1, 1) = v[1];
    frame_cprop_[*c_it](2, 1) = v[2];

    frame_cprop_[*c_it](0, 2) = w[0];
    frame_cprop_[*c_it](1, 2) = w[1];
    frame_cprop_[*c_it](2, 2) = w[2];
  }
}


template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::import_frames_from_quaternion_property_onering(const VH &_vh)
{
  auto cell_quaternions = mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions");
  for (auto vc_it = mesh_.vc_iter(_vh); vc_it.valid(); ++vc_it)
    frame_cprop_[*vc_it] = cell_quaternions[*vc_it].toRotationMatrix();;
}


template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::import_quaternions_from_frames()
{
  auto cell_quaternions = mesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions");
  mesh_.set_persistent(cell_quaternions, true);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    Mat3d R;
    closest_rotation(frame_cprop_[*c_it], R);
    cell_quaternions[*c_it] = Quaternion(R);

//    std::cerr << "Frame     : " << std::endl << frame_cprop_[*c_it] << std::endl;
//    std::cerr << "Quaternion: " << std::endl << cell_quaternions[*c_it].toRotationMatrix() << std::endl << std::endl;
  }
}

template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::closest_rotation(const Mat3d &_A, Mat3d &_R) const
{
  Eigen::JacobiSVD <Mat3d> svd(_A, Eigen::ComputeFullU | Eigen::ComputeFullV);

  Mat3d U = svd.matrixU();
  Mat3d V = svd.matrixV();

  // correct orientations if necessary
  if (U.determinant() < 0.0)
    U.col(0) *= -1.0;
  if (V.determinant() < 0.0)
    V.col(0) *= -1.0;

  _R = U * V.transpose();
}


template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::check_valence_consistency()
{
  std::cerr << "---------- Check Valence/Transition Functions consistency -------------" << std::endl;
  bool is_valid = true;

  int n_invalid_sectors(0);

  for (EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
  {
    if (!check_sector_valencies(*e_it))
    {
      ++n_invalid_sectors;
      is_valid = false;
    }

    bool e_boundary = mesh_.is_boundary(*e_it);

    HEH heh = mesh_.halfedge_handle(*e_it, 0);
    HEHFIt hehf_it(heh, &mesh_);

    // skip first face for boundary edges
    if (mesh_.is_boundary(mesh_.opposite_halfface_handle(*hehf_it)))
      ++hehf_it;

    HFH hfh0 = *hehf_it;
    CH ch0 = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh0));
    int total_transition = 0;

    // target valence of edge
    int valence_e = valence_eprop_[*e_it];

    // get adjacent vertices
    VH vh0 = mesh_.halfedge(heh).from_vertex();
    VH vh1 = mesh_.halfedge(heh).to_vertex();

    // get edge vector
    Vec3d e = ovm2eigen(mesh_.vertex(vh1) - mesh_.vertex(vh0));

    // determine alignment axis
    Mat3d J = frame_cprop_[ch0].inverse();
    Vec3d e_trans = J * e;
    AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(e_trans);

    for (; hehf_it.valid(); ++hehf_it)
    {
      // check alignment in cell if singular or feature
      if (valence_e != 0 || feature_eprop_[*e_it] > 0)
      {
        CH ch_cur = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
        Vec3d e_local = frame_cprop_[ch_cur].inverse() * e;

        if ((e_trans.cross(e_local)).norm() / (e_trans.norm() * e_local.norm()) > 1e-3 || e_trans.dot(e_local) < 0.0)
        {
          std::cerr << "ERROR: alignment axis inconsistent " << e_trans.transpose() / e_trans.norm() << std::endl
                    << "                   should align to " << e_local.transpose() / e_local.norm() << std::endl;
          is_valid = false;
        }
      }

      // transform to next cell (if not on boundary)
      if (!mesh_.is_boundary(mesh_.face_handle(*hehf_it)))
      {
        // transform to next cell
        int cur_trans = transition_hfprop_[mesh_.opposite_halfface_handle(*hehf_it)];
        total_transition = tq_.mult_transitions_idx(total_transition, cur_trans);
        // transform alignment axis
        e_trans = tq_.transition_matrix_int(tq_.inverse_transition_idx(cur_trans)) * e_trans;
      }
    }

    if (!tq_.is_2d_rotation(total_transition) && !e_boundary)
    {
      std::cerr << "ERROR: singularity type is not a 2D rotation. type = " << total_transition << std::endl;
      is_valid = false;
    }

    AxisAlignment t_axis = AxisAlignment(tq_.rotation_axis_2d(total_transition));

    if (total_transition != 0 && int(e_axis) / 2 != int(t_axis) / 2 && !e_boundary)
    {
      std::cerr << "ERROR: transition function axis inconsistent to alignment axis. e_axis = " << int(e_axis)
                << ", t_axis = " << int(t_axis) << std::endl;
      is_valid = false;
    }

    // transition rotation is inverse to field rotation
    int valence_t = tq_.rotation_angle_2d(total_transition) / 90;
    if (e_axis != t_axis)
      valence_t *= -1;

    if (valence_e != valence_t && !e_boundary)
    {
      std::cerr << "ERROR: incorrect valence.  edge-property-valence = " << valence_e << ", transition-frame-valence = "
                << valence_t << std::endl;
      is_valid = false;
    }

//      std::cerr << "ede-property-valence = " << valence_e << "transition-frame-valence = " << valence_t << ", total_transition = " << total_transition << ", axis = " << t_axis << ", edge-axis= " << e_axis << std::endl;
  }

  std::cerr << "#invalid sectors = " << n_invalid_sectors << std::endl;
  std::cerr << "---------- Validity of Valence/Transition Functions = " << int(is_valid) << " -------------"
            << std::endl;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
bool
FrameFieldOptimizer3DT<TetMeshT>::
check_sector_valencies(const EH _eh)
{
  bool is_valid = true;

  // get halfedge
  HEH heh = mesh_.halfedge_handle(_eh, 0);
  // get first halfface
  HEHFIt hehf_it(heh, &mesh_);
  HFH hfh0 = *hehf_it;
  FH fh0 = mesh_.face_handle(hfh0);
  if (mesh_.is_boundary(_eh) && feature_fprop_[fh0] == 0)
  {
//    std::cerr << "Warning: detected boundary face without feature tag such that the valence is undefined. face.idx() =  " << fh0.idx() << std::endl;
    return false;
  }
  else
  {
    HFH hfh_feat = find_feature_halfface(heh);
    if (hfh_feat.is_valid())
    {
      hfh0 = hfh_feat;
      fh0 = mesh_.face_handle(hfh0);
    }
  }

  HFH hfh_sector_start = hfh0;
  HFH hfh_sector_end;
  std::vector<double> sector_valencies;

  for (int i = 0; i <= 10000; ++i) // limit number of sectors to 10000 to prevent infinite loops
  {
    // sector angle
    double sa = calc_sector_angle(heh, hfh_sector_start, hfh_sector_end);
    // sector valence
    double sv = sa / (0.5 * M_PI);
    int sv_int = std::round(sv);
//    if(sv_int < 1 || std::abs(sv-double(sv_int)) > 1e-2) // hack to avoid outputs for non-orthognoal frames which are not handled correctly yet
    if (sv_int < 1 || std::abs(sv - double(sv_int)) > 0.5)
    {
      is_valid = false;
    }

    // store valencies
    sector_valencies.push_back(sv);

    FH fh_sector_end = mesh_.face_handle(hfh_sector_end);

    if (mesh_.is_boundary(fh_sector_end) || (fh_sector_end == fh0))
      break;

    hfh_sector_start = hfh_sector_end;
  }

  if (!is_valid || sector_valencies.size() >= 10000)
  {
//    std::cerr << "ERROR: invalid frame field sector with " << sector_valencies.size() << " sector valencies: " << std::endl;
//    for(auto sv : sector_valencies) std::cerr << sv << " ";
//    std::cerr << std::endl;
    std::cerr << "ERROR: invalid frame field sector with " << sector_valencies.size() << " sector valencies at edge "
              << _eh << " hfhs " << hfh_sector_start << " " << hfh_sector_end << std::endl;

  }

  return is_valid;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
HFH
FrameFieldOptimizer3DT<TetMeshT>::
find_feature_halfface(const HEH _heh)
{
  HEHFIt hehf_it(_heh, &mesh_);
  for (; hehf_it.valid(); ++hehf_it)
  {
    if (feature_fprop_[mesh_.face_handle(*hehf_it)] > 0)
      return *hehf_it;
  }
  // return invalid handle if no feature face found
  return HFH(-1);
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
double
FrameFieldOptimizer3DT<TetMeshT>::
calc_sector_angle(const HEH _heh, const HFH _hfh_sector_start, HFH &_hfh_sector_end)
{
  if (mesh_.is_boundary(_hfh_sector_start))
  {
    std::cerr << "Warning: sector start half-face is on boundary, which should never happen" << std::endl;
    return 0.0;
  }

  HFH hfh0 = _hfh_sector_start;
  FH fh0 = mesh_.face_handle(hfh0);
  CH ch0 = mesh_.incident_cell(hfh0);

  int fh0_feat = feature_fprop_[fh0];

  // get adjacent vertices
  VH vh0 = mesh_.halfedge(_heh).from_vertex();
  VH vh1 = mesh_.halfedge(_heh).to_vertex();

  // get edge vector
  Vec3d e = ovm2eigen(mesh_.vertex(vh1) - mesh_.vertex(vh0));
  Vec3d u, v, w;
  complement_to_right_handed_orthonormal_frame(e, u, v, w);

  // determine alignment axis
  Mat3d F0 = frame_cprop_[ch0];
  Mat3d J0 = frame_cprop_[ch0].inverse();
  Vec3d e_trans = J0 * e;
  AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(e_trans);

  // get two tangent vectors
  AxisAlignment a0_axis = AxisAlignment((int(e_axis) + 2) % 6);
  AxisAlignment b0_axis = AxisAlignment((int(a0_axis) + 2) % 6);

  //correct (e, a0, b0) to right-handed frame if necessary, if sign(e)==-1 -> swap sign of a0
  if (!AxisAlignmentHelpers::is_positive_axis(e_axis))
    a0_axis = reverse_axis(a0_axis);

  Vec3d n0 = ovm2eigen(mesh_.normal(hfh0));
  Vec3d n0_trans = F0.transpose() * n0;  // normal vectors transform with the inverse-transposed
  AxisAlignment n0_axis = AxisAlignmentHelpers::get_dominant_axis(n0_trans);

  if (fh0_feat > 0) // check alignment
  {
    Vec3d n0_axis_vec = AxisAlignmentHelpers::vector(n0_axis);
    if (n0_axis_vec.dot(n0_trans) < 0.0 ||
        (n0_axis_vec.cross(n0_trans)).norm() / (n0_axis_vec.norm() * n0_trans.norm()) > 1e-3)
      std::cerr << "Warning: frame field is not aligned to feature face, normal maps to " << n0_trans.transpose()
                << std::endl;
  }

  // prefer axis closest to face normal since it is more stable for regular edges
  if (unoriented_axis(n0_axis) != unoriented_axis(a0_axis))
  {
    std::swap(a0_axis, b0_axis);
    a0_axis = reverse_axis(a0_axis);
  }

  // param domain vectors of axes for rotation estimation
  Vec3d a0_trans = AxisAlignmentHelpers::vector(a0_axis);
  Vec3d b0_trans = AxisAlignmentHelpers::vector(b0_axis);

  double rotation_angle_a(0.0);
  double rotation_angle_b(0.0);
  double dihedral_angle(0.0);

  // collect data for debugging
  std::vector<double> d_angles;
  std::vector<double> ra_angles;
  std::vector<double> rb_angles;
  std::vector<int> tfs;

  for (int i = 0; i < 10000; ++i) // restrict to 10000 sectors to avoid infinite loops
  {
    // add up dihedral angle
    double da = dihedral_angle_in_cell(mesh_, _heh, hfh0);
    dihedral_angle += da;
    // debug helper
    d_angles.push_back(da);

    // get next halfface
    HFH hfh1 = mesh_.adjacent_halfface_in_cell(hfh0, _heh);
    HFH hfh1_opp = mesh_.opposite_halfface_handle(hfh1);
    FH fh1 = mesh_.face_handle(hfh1_opp);
    if (feature_fprop_[fh1] > 0 || mesh_.is_boundary(fh1))
    {
      _hfh_sector_end = hfh1_opp;
      break;
    }
    else
    {
      // debug helper
      tfs.push_back(transition_hfprop_[hfh1]);

      CH ch1 = mesh_.incident_cell(hfh1_opp);
      // transform F0 to coordinate system of F1
      Mat3d F1 = frame_cprop_[ch1];
      F0 = (F0 * tq_.transition_matrix_int(transition_hfprop_[hfh1])).eval();
      // transform axes a and b
      auto TM = tq_.transition_matrix_int(transition_hfprop_[hfh1_opp]); // inverse transition
      a0_trans = (TM * a0_trans).eval();
      b0_trans = (TM * b0_trans).eval();

      // get corresponding frame vectors
      Vec3d a0 = F0 * a0_trans;
      Vec3d b0 = F0 * b0_trans;
      Vec3d a1 = F1 * a0_trans;
      Vec3d b1 = F1 * b0_trans;

      // accumulate 2d rotation angles around edge
      double ra = signed_angle(Vec2d(a0.dot(v), a0.dot(w)), Vec2d(a1.dot(v), a1.dot(w)));
      double rb = signed_angle(Vec2d(b0.dot(v), b0.dot(w)), Vec2d(b1.dot(v), b1.dot(w)));
      rotation_angle_a += ra;
      rotation_angle_b += rb;

      ra_angles.push_back(ra);
      rb_angles.push_back(rb);

      if (hfh1_opp == _hfh_sector_start)
      {
        _hfh_sector_end = hfh1_opp;
        break;
      }

      // update variables for next step
      std::swap(F0, F1);
      std::swap(hfh0, hfh1_opp);
    }
  }

//  fh0 = mesh_.face_handle(_hfh_sector_end);
//  if(feature_fprop_[fh0]) // check alignment
//  {
//    Vec3d n0 = ovm2eigen(mesh_.normal(fh0));
//    Vec3d n0_trans = F.transpose()*n0;  // normal vectors transform with the inverse-transposed
//    AxisAlignment n0_axis = AxisAlignmentHelpers::get_dominant_axis(n0_trans);
//    Vec3d n0_axis_vec = AxisAlignmentHelpers::vector(n0_axis);
//
//    if( n0_axis_vec.dot(n0_trans) < 0.0 || n0_axis_vec.cross(n0_trans)/(n0_axis_vec.norm()*n0_trans.norm()) > 1e-3)
//      std::cerr << "Warning: frame field is not aligned to feature face, normal maps to " << n0_trans.transpose() << std::endl;
//  }

  // TODO handle non-orthogonal frames, which require a correction angle!!!
  double sector_angle = dihedral_angle - rotation_angle_a;
  double sector_valence = std::round(sector_angle / (0.5 * M_PI));

  // output information for invalid sectors
  if (sector_angle < M_PI * 0.25 || std::abs(sector_valence * 0.5 * M_PI - sector_angle) > M_PI * 0.4)
  {
    EH eh = mesh_.edge_handle(_heh);
    std::cerr << "--------- severe invalid sector ----------" << std::endl;
    std::cerr << "edge valence     = " << valence_eprop_[eh] << std::endl;
    std::cerr << "is boundary      = " << int(mesh_.is_boundary(eh)) << std::endl;
    std::cerr << "dihedral angle   = " << dihedral_angle / M_PI * 180.0 << std::endl;
    std::cerr << "rotation angle a = " << rotation_angle_a / M_PI * 180.0 << std::endl;
    std::cerr << "rotation angle b = " << rotation_angle_b / M_PI * 180.0 << std::endl;
    std::cerr << "sector angle     = " << (dihedral_angle - rotation_angle_a) / M_PI * 180.0 << std::endl;
    std::cerr << "dihedral angles: ";
    for (double da: d_angles) std::cerr << da / M_PI * 180.0 << " ";
    std::cerr << std::endl;
    std::cerr << "rot a angles   : ";
    for (double a: ra_angles) std::cerr << a / M_PI * 180.0 << " ";
    std::cerr << std::endl;
    std::cerr << "rot b angles   : ";
    for (double a: rb_angles) std::cerr << a / M_PI * 180.0 << " ";
    std::cerr << std::endl;
    std::cerr << "transitions    : ";
    for (int t: tfs) std::cerr << t << " ";
    std::cerr << std::endl;
  }

  return sector_angle;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::
check_frame_rotation_angles()
{
  double min_angle = DBL_MAX;
  double max_angle = -DBL_MAX;

  std::vector<double> angles(mesh_.n_faces());

  for (auto fh: mesh_.faces())
    if (!mesh_.is_boundary(fh))
    {
      auto hfh0 = mesh_.halfface_handle(fh, 0);
      auto hfh1 = mesh_.halfface_handle(fh, 1);

      auto ch0 = mesh_.incident_cell(hfh0);
      auto ch1 = mesh_.incident_cell(hfh1);

      auto F0 = frame_cprop_[ch0];
      auto F1 = frame_cprop_[ch1];

      int trans01 = transition_hfprop_[hfh0];

      Mat3d A = F1.inverse() * F0 * tq_.transition_matrix_int(trans01);

      Mat3d R;
      closest_rotation(A, R);

      Eigen::AngleAxis<double> aa(R);

      double angle = aa.angle() / M_PI * 180.0;

      if (angle < min_angle) min_angle = angle;
      if (angle > max_angle) max_angle = angle;

      angles.push_back(angle);
    }

  std::cerr << "Min frame rotation angle: " << min_angle << std::endl;
  std::cerr << "Max frame rotation angle: " << max_angle << std::endl;

  std::sort(angles.begin(), angles.end(), std::greater<double>());

  std::cerr << "largest angles ";
  for (int i = 0; i < std::min(int(angles.size()), int(30)); ++i)
    std::cerr << angles[i] << ", ";
  std::cerr << std::endl;
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
int
FrameFieldOptimizer3DT<TetMeshT>::
incident_feature_ohalfedges(const VH _vh, std::vector<HEH> &_fheh) const
{
  _fheh.clear();
  // is on feature edge?
  VOHEIt vhe_it(_vh, &mesh_);
  for (; vhe_it.valid(); ++vhe_it)
    if (feature_eprop_[mesh_.edge_handle(*vhe_it)] > 0)
      _fheh.push_back(*vhe_it);

  return _fheh.size();
}

//-----------------------------------------------------------------------------

template<class TetMeshT>
int
FrameFieldOptimizer3DT<TetMeshT>::
incident_feature_faces(const VH _vh, std::vector<FH> &_ffh) const
{
  _ffh.clear();
  // is on feature triangle?
  VFIt vf_it(_vh, &mesh_);
  for (; vf_it.valid(); ++vf_it)
    if (feature_fprop_[*vf_it] > 0)
      _ffh.push_back(*vf_it);

  return _ffh.size();
}



//-----------------------------------------------------------------------------

template<class TetMeshT>
void
FrameFieldOptimizer3DT<TetMeshT>::
restore_field_alignment(const std::map<CH, CellAlignmentInfo> &_cell_alignment_info)
{
  for (auto mi: _cell_alignment_info)
  {
    CH ch = mi.first;
    CellAlignmentInfo cai = mi.second;

    // get matrix of gradient vectors
    Mat3d G = (frame_cprop_[ch].inverse()).transpose();

    // has face alignment constraint?
    if (cai.hfh.is_valid())
    {
      if (cai.heh.is_valid()) // has both alignment constraints
      {
        // get normal vector
        Vec3d n = vec2eigen(mesh_.normal(cai.hfh));
        n /= n.norm();

        // get edge vector
        Vec3d e = vec2eigen(mesh_.vector(cai.heh));
        e /= e.norm();

        if (!std::isfinite(n.sum()) || !std::isfinite(e.sum()))
        {
          std::cerr << "Warning: restore_field_alignment called with degenerate normal n=" << n.transpose()
                    << ", or edge=" << e.transpose() << std::endl;
        }
        else
        {
          double dp = n.dot(G.col(cai.hfh_axis));
          if (dp <= 0.0)
            std::cerr
                    << "Warning: restore_field_alignment received face alignment constraints with locally inverted axis..."
                    << std::endl;
          G.col(cai.hfh_axis) = n * dp;

          G.col(cai.heh_axis0) -= e * (e.dot(G.col(cai.heh_axis0)));
          G.col(cai.heh_axis1) -= e * (e.dot(G.col(cai.heh_axis1)));
        }
      }
      else // has only face alignment constraint?
      {
        // get normal vector
        Vec3d n = vec2eigen(mesh_.normal(cai.hfh));
        n /= n.norm();

        if (!std::isfinite(n.sum()))
        {
          std::cerr << "Warning: restore_field_alignment called with degenerate normal n =" << n.transpose()
                    << std::endl;
        }
        else
        {
          double dp = n.dot(G.col(cai.hfh_axis));

          if (dp <= 0.0)
            std::cerr
                    << "Warning: restore_field_alignment received face alignment constraints with locally inverted axis..."
                    << std::endl;

          G.col(cai.hfh_axis) = n * dp;
        }
      }
    }
    else if (cai.heh.is_valid()) // has only edge alignment constraint
    {
      // get edge vector
      Vec3d e = vec2eigen(mesh_.vector(cai.heh));
      e /= e.norm();

      if (!std::isfinite(e.sum()))
      {
        std::cerr << "Warning: restore_field_alignment called with degenerate edge e=" << e.transpose() << std::endl;
      }
      else
      {
        G.col(cai.heh_axis0) -= e * (e.dot(G.col(cai.heh_axis0)));
        G.col(cai.heh_axis1) -= e * (e.dot(G.col(cai.heh_axis1)));
      }
    }

    if (G.determinant() <= 0.0)
      std::cerr << "Warning: G is inverted after restore_field_alignment " << std::endl << G << std::endl;

    // store final frame
    frame_cprop_[ch] = G.inverse().transpose();
  }
}

//=============================================================================
} // namespace AlgoHex
//=============================================================================

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define PARAMETRIZATION3DT_C


//== INCLUDES =================================================================

//#include <ACG/Utils/HaltonColors.hh>

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include "Parametrization3DT.hh"
#include "Geometry.hh"
#include "MeshGeometry.hh"
#include "LinearLeastSquaresElements.hh"

#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/IPOPTSolver.hh>
#include <CoMISo/NSolver/TruncatedNewtonPCG.hh>
#include <CoMISo/NSolver/NPDerivativeChecker.hh>
#include <CoMISo/NSolver/NPTiming.hh>
#include <CoMISo/NSolver/ConstraintTools.hh>
#include <CoMISo/NSolver/LinearConstraintConverter.hh>
#include <CoMISo/Utils/Tools.hh>
#include <AlgoHex/Util/StopWatch.hh>
#include "Stopwatches.hh"

#if COMISO_GUROBI_AVAILABLE
#include <CoMISo/NSolver/GUROBISolver.hh>
#endif

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== IMPLEMENTATION ==========================================================


template<class TetMeshT>
std::pair<int, int>
Parametrization3DT<TetMeshT>::
parametrize_complete(const int _num_hex_cells, const double _anisotropy_alpha, bool _with_integer_constraints)
{
  ScopedStopWatch sw_total(sw::parameterization);
  std::cout << "#####Parametrizing (complete pipeline)..." << std::endl;

  // solve quadratic problem
  parametrize(_num_hex_cells, _anisotropy_alpha, _with_integer_constraints);

  OpenVolumeMesh::CellPropertyT<std::map<OpenVolumeMesh::VertexHandle, Point> > best_solution(
          mesh_.template request_cell_property<std::map<OpenVolumeMesh::VertexHandle, Point> >(
                  "ParametrizationCopy"));
  copy_cell_property(igm_cprop_, best_solution);
  std::pair<int, int> best_invalid;
  best_invalid.first = number_of_invalid_parametric_tets();
  best_invalid.second = number_of_invalid_edge_valencies();

  std::vector<std::pair<int, int> > degeneracy_sequence;
  degeneracy_sequence.push_back(best_invalid);

  double epsilon = 0.1;
  double w_exp = 1.0;
  if (best_invalid.first > 0 || best_invalid.second > 0)
//    for(int j=0; j<5; ++j)
  {
    std::cerr << "---> remaining defects: #invalid_tets = " << best_invalid.first << " and #invalid_edge_valencies = "
              << best_invalid.second << std::endl;
    std::cerr << "---> optimize locally injective with epsilon = " << epsilon << std::endl;

    // solve convex problem with sizing DOFs
    optimize_integer_grid_map_locally_injective(_num_hex_cells, _anisotropy_alpha, epsilon, _with_integer_constraints);
    std::pair<int, int> cur_invalid;
    cur_invalid.first = number_of_invalid_parametric_tets();
    cur_invalid.second = number_of_invalid_edge_valencies();

    degeneracy_sequence.push_back(cur_invalid);

    if (cur_invalid >= best_invalid)
    {
      // recover better solution
      copy_cell_property(best_solution, igm_cprop_);
//        break;
    }
    else
    {
      best_invalid = cur_invalid;
      copy_cell_property(igm_cprop_, best_solution);
    }

//      // valid solution found?
//      if(best_invalid.first == 0 && best_invalid.second == 0)
//        break;

//      w_exp *= 2.0;
  }

  // recover better solution
  copy_cell_property(best_solution, igm_cprop_);

  std::cerr << "-----------> parametrization finished with #invalid_tets = " << best_invalid.first
            << " and #invalid_edge_valencies = " << best_invalid.second << std::endl;
  std::cerr << "degeneracy sequence: ";
  for (auto i: degeneracy_sequence) std::cerr << "(#inv-tets = " << i.first << "#inv-edges = " << i.second << ") ";
  std::cerr << std::endl;

  return best_invalid;
}


template<class TetMeshT>
std::pair<int, int>
Parametrization3DT<TetMeshT>::
parametrize_robust_quantization(const int _num_hex_cells)
{
  ScopedStopWatch sw_total(sw::parameterization);
  std::cout << "#####Parametrizing (robust quantization pipeline)..." << std::endl;

  // quantize
  quantize_qgp3d(_num_hex_cells);

  // optional: check validity of path constraints and output difference within seamless map
  check_quantization_path_constraints();

  // construct IGM subject to quantization constraints
  optimize_integer_grid_map_robust_quantization(_num_hex_cells);
//  parametrize(_num_hex_cells, 0.5, true);

  std::pair<int, int> invalid;
  invalid.first = number_of_invalid_parametric_tets();
  invalid.second = number_of_invalid_edge_valencies();

  std::cerr << "-----------> robust quantization parametrization finished with #invalid_tets = " << invalid.first
            << " and #invalid_edge_valencies = " << invalid.second << std::endl;

  return invalid;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
parametrize(const int _num_hex_cells, const double _anisotropy_alpha, bool _with_integer_constraints)
{
  ScopedStopWatch sw_total(sw::parameterization);

  std::cout << "#####Parametrizing..." << std::endl;

//  test_objective_function();
//  return;

  // 0. check input mesh
  check_input_mesh();

  // 1. generate cut surface
  generate_cut_surface();

  // 2. comb octahedral field
  comb_octahedral_field();

  // 3. optimize integer grid map
  optimize_integer_grid_map(_num_hex_cells, _anisotropy_alpha, _with_integer_constraints);
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
check_input_mesh() const
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "check input mesh..." << std::endl;)

  // check edge lengths
  double min_el(DBL_MAX);
  double max_el(0.0);

  for (EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
  {
    VH vh0 = mesh_.halfedge(mesh_.halfedge_handle(*e_it, 0)).to_vertex();
    VH vh1 = mesh_.halfedge(mesh_.halfedge_handle(*e_it, 1)).to_vertex();

    double el = (mesh_.vertex(vh0) - mesh_.vertex(vh1)).norm();

    if (el > max_el) max_el = el;
    if (el < min_el) min_el = el;
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "min/max edge length: " << min_el << " / " << max_el << std::endl;)

  // check volume
  double min_vol(DBL_MAX);
  double max_vol(0.0);

  Eigen::Matrix<double, 3, 4> P;
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // collect points
    unsigned int i = 0;
    for (TVIt tv_it = mesh_.tv_iter(*c_it); tv_it.valid(); ++tv_it)
      P.col(i++) = ovm2eigen(mesh_.vertex(*tv_it));

    Mat3x3 E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    double vol = 1.0 / 6.0 * E.determinant();

    if (vol > max_vol) max_vol = vol;
    if (vol < min_vol) min_vol = vol;
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "min/max volume: " << min_vol << " / " << max_vol << std::endl;)
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
generate_cut_surface()
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "generate cut surface..." << std::endl;)
  // initialize
  //   OpenVolumeMesh::FacePropertyT<bool> is_on_cut_fprop_;
  //   std::vector< std::vector<FH> > cut_surface_;

  // 1. start with trivial cut (all faces are cut)
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    is_on_cut_fprop_[*f_it] = !mesh_.is_boundary(*f_it);

  // only continue if mesh has cells
  if (mesh_.n_cells() > 0)
  {
    // 2. reduce cut-faces by dual spanning tree
    std::vector<bool> cvisited(mesh_.n_cells(), false);
    std::queue<HFH> qhfh;
    // add seed
    CH ch(0);
    cvisited[ch.idx()] = true;
    for (unsigned int i = 0; i < mesh_.cell(ch).halffaces().size(); ++i)
      qhfh.push(mesh_.cell(ch).halffaces()[i]);

    // grow tree until all cells in connected component visited
    int n_removed_by_spanning_tree(0);
    while (!qhfh.empty())
    {
      // get opposite of next halfface
      HFH hfh = mesh_.opposite_halfface_handle(qhfh.front());
      qhfh.pop();

      // only process if not on boundary
      if (!mesh_.is_boundary(hfh))
      {
        // check if cell already visited
        CH ch = mesh_.incident_cell(hfh);
        if (!cvisited[ch.idx()])
        {
          // remove face from cutgraph
          is_on_cut_fprop_[mesh_.face_handle(hfh)] = false;
          ++n_removed_by_spanning_tree;
          // process cell
          cvisited[ch.idx()] = true;
          for (unsigned int i = 0; i < mesh_.cell(ch).halffaces().size(); ++i)
            qhfh.push(mesh_.cell(ch).halffaces()[i]);
        }
      }
    }

    ALGOHEX_DEBUG_ONLY(
            std::cerr << "#removed by spanning tree: " << n_removed_by_spanning_tree << " (#cells " << mesh_.n_cells()
                      << ")" << std::endl;)

    // 3. further shrink cut surface
    std::queue<EH> qeh;
    int n_shrinked(0);
    // initialize queue with edges (1) not on boundary, (2) not singularity and (3) adjacent to one cut face
    for (EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
      if (!mesh_.is_boundary(*e_it) &&
          valence_eprop_[*e_it] == 0 &&
          n_cut_faces(*e_it) == 1)
        qeh.push(*e_it);

    ALGOHEX_DEBUG_ONLY(std::cerr << "#edge queue: " << qeh.size() << std::endl;)
    std::vector<FH> fhs;

    // shrink while possible
    while (!qeh.empty())
    {
      // get next element
      EH eh = qeh.front();
      qeh.pop();

      // still a valid candidate?
      if (n_cut_faces(eh, fhs) == 1)
      {
        // remove from cut
        FH fh = fhs.front();
        is_on_cut_fprop_[fh] = false;
        ++n_shrinked;
        // push new candidates
        for (unsigned int i = 0; i < mesh_.face(fh).halfedges().size(); ++i)
        {
          eh = mesh_.edge_handle(mesh_.face(fh).halfedges()[i]);
          if (!mesh_.is_boundary(eh) &&
              valence_eprop_[eh] == 0 &&
              n_cut_faces(eh) == 1)
            qeh.push(eh);
        }
      }
    }

    {
      for (EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
        if (!mesh_.is_boundary(*e_it) &&
            valence_eprop_[*e_it] == 0 &&
            n_cut_faces(*e_it) == 1)
          std::cerr << "ERROR: there are shrinking candidates left but queue is empty!" << std::endl;
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "#shrinked faces: " << n_shrinked << std::endl;)

    // 4. decompose cut surface into manifold pieces
    cut_surface_.clear();
    std::set<FH> fvisited;
    for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
      if (is_on_cut_fprop_[*f_it] && fvisited.count(*f_it) == 0)
      {
        // grow a new manifold patch
        cut_surface_.resize(cut_surface_.size() + 1);
        cut_surface_.back().push_back(*f_it);
        fvisited.insert(*f_it);

        // init queue
        std::queue<EH> qeh2;
        for (unsigned int i = 0; i < mesh_.face(*f_it).halfedges().size(); ++i)
          qeh2.push(mesh_.edge_handle(mesh_.face(*f_it).halfedges()[i]));

        std::vector<FH> fhs;
        while (!qeh2.empty())
        {
          EH eh = qeh2.front();
          qeh2.pop();

          if (!mesh_.is_boundary(eh) && valence_eprop_[eh] == 0 && n_cut_faces(eh, fhs) == 2)
          {
            for (unsigned int i = 0; i < fhs.size(); ++i)
            {
              // grow to neighbor if possible
              FH fh_new = fhs[i];
              if (fvisited.count(fh_new) == 0)
              {
                // add to current cluster and mark as visited
                cut_surface_.back().push_back(fh_new);
                fvisited.insert(fh_new);

                // push candidate edges
                for (unsigned int j = 0; j < mesh_.face(fh_new).halfedges().size(); ++j)
                  qeh2.push(mesh_.edge_handle(mesh_.face(fh_new).halfedges()[j]));
              }
            }
          }
        }
      }

    ALGOHEX_DEBUG_ONLY(std::cerr << "#2-manifold patches in cut surface: " << cut_surface_.size() << std::endl;)

    // 5. visualize cut surface with alpha channel
//    for(FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
//      if(is_on_cut_fprop_[*f_it])
//        mesh_color_[*f_it] = ACG::Vec4f(1.0,0.0,0.0,1.0);
//      else
//        mesh_color_[*f_it] = ACG::Vec4f(0.5,0.5,0.5,0.1);
//
//    // 6. visualize 2-manifold patches with different colors
//    ACG::HaltonColors hcol;
//    for(unsigned int i=0; i<cut_surface_.size(); ++i)
//    {
//      ACG::Vec4f cur_col = hcol.get_next_color();
//      for(unsigned int j=0; j<cut_surface_[i].size(); ++j)
//        mesh_color_[cut_surface_[i][j]] = cur_col;
//    }
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
comb_octahedral_field()
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "comb octahedral field... ";)
//  std::cerr << "octahedral field energy before combing: " << octahedral_field_energy() << std::endl;

  // only continue if mesh has cells
  if (mesh_.n_cells() > 0)
  {
    // process dual spanning tree
    std::vector<bool> cvisited(mesh_.n_cells(), false);
    std::queue<HFH> qhfh;
    // add seed
    CH ch(0);
    cvisited[ch.idx()] = true;
    for (unsigned int i = 0; i < mesh_.cell(ch).halffaces().size(); ++i)
      qhfh.push(mesh_.cell(ch).halffaces()[i]);

    // grow tree until all cells in connected component visited
    while (!qhfh.empty())
    {
      // get opposite of next halfface
      HFH hfh0 = qhfh.front();
      HFH hfh1 = mesh_.opposite_halfface_handle(hfh0);
      FH fh = mesh_.face_handle(hfh0);
      qhfh.pop();

      // only process if not on boundary
      if (!mesh_.is_boundary(fh) && !is_on_cut_fprop_[fh])
      {
        // check if cell already visited
        CH ch_new = mesh_.incident_cell(hfh1);
        if (!cvisited[ch_new.idx()])
        {
          // mark as visited
          cvisited[ch_new.idx()] = true;

          // get matching
          int trans = transition_hfprop_[hfh0];

          // transform frame
//                        Eigen::Quaterniond qtrans_inv = tq_.inverse_transition(trans);
//                        quaternion_cprop_[ch_new] = quaternion_cprop_[ch_new] * qtrans_inv;
          Mat3x3 trans_inv = tq_.transition_matrix_int(tq_.inverse_transition_idx(trans));
          frame_cprop_[ch_new] = frame_cprop_[ch_new] * trans_inv;


          // process half-faces to neighbors
          for (unsigned int i = 0; i < mesh_.cell(ch_new).halffaces().size(); ++i)
          {
            // get new half-face
            HFH hfh_new = mesh_.cell(ch_new).halffaces()[i];
            HFH hfh_new_opp = mesh_.opposite_halfface_handle(hfh_new);

            // transform matching and inverse matching
            transition_hfprop_[hfh_new] = tq_.mult_transitions_idx(trans, transition_hfprop_[hfh_new]);
            transition_hfprop_[hfh_new_opp] = tq_.inverse_transition_idx(transition_hfprop_[hfh_new]);

            // push as candidate
            qhfh.push(hfh_new);
          }
        }
      }
    }
  }

  // check result
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (!mesh_.is_boundary(*f_it))
    {
      bool on_cut = is_on_cut_fprop_[*f_it];
      int trans = transition_hfprop_[mesh_.halfface_handle(*f_it, 0)];

      if (!on_cut && trans != 0)
        std::cerr << "ERROR: after combing there is still a non-identity transition away from the cut "
                  << trans << std::endl;
    }
  ALGOHEX_DEBUG_ONLY(std::cerr << "done" << std::endl;)
//  std::cerr << "octahedral field energy after combing: " << octahedral_field_energy() << std::endl;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
optimize_integer_grid_map(const int _num_hex_cells, const double _anisotropy_alpha,
                          const bool _with_integer_constraints)
{
  // compute target edge length
  // ToDo: calculate volume based on frames
  //      double sizing_scale_factor = std::cbrt(mesh_volume() / _num_hex_cells);
  double sizing_scale_factor = std::cbrt(frame_field_volume() / _num_hex_cells);
  ALGOHEX_DEBUG_ONLY(std::cerr << "sizing scale factor: " << sizing_scale_factor << std::endl;)
  // create multiple copies of vertices at cut and store indices in global_idx_
  setup_vertex_indices();
  check_vertex_indices();

  // setup objective function as finite-element problem
  reset_cell_weights();
  COMISO::FiniteElementSet<FrameFittingElement3D> fe_ff;
  setup_frame_objective_function(fe_ff, _anisotropy_alpha, sizing_scale_factor);

  // setup constraints
  std::vector<PairUiV> discrete_constraints;
  std::vector<COMISO::LinearConstraint> constraints;
  std::vector<COMISO::NConstraintInterface *> constraint_pointers;
  setup_constraints(constraints, constraint_pointers, discrete_constraints);

  // hack
  // make constraints linearly independent
  COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraint_pointers);

  if (!_with_integer_constraints)
    discrete_constraints.clear();

  // solve and store result
  unsigned int n_complete = 3 * (max_global_idx_ + 1) + n_helper_variables_;
  COMISO::FiniteElementProblem fe_problem(n_complete);
  fe_problem.add_set(&fe_ff);

  //  std::cout << "---- (optional for debugging) Check derivatives of problem..." << std::endl;
//  COMISO::NPDerivativeChecker npd;
//  npd.check_all(&fe_problem);

  // configure comiso solver
  COMISO::COMISOSolver comiso;
  comiso.solver().misolver().set_local_iters(0);
  comiso.solver().misolver().set_inital_full(true);
  comiso.solver().misolver().set_final_full(true);
  comiso.solver().misolver().set_iter_full(false);
  comiso.solver().misolver().set_cg_iters(20);
  comiso.solver().misolver().set_cg_error(1e-3);
  comiso.solver().misolver().set_multiple_rounding();
  comiso.solver().set_support_constraint_rhs_resolve(false); // this saves a ton of memory

//  comiso.solver().misolver().set_noise(0);
//  comiso.solver().set_noisy(0);

//  comiso.solver().misolver().set_gurobi_rounding(true);
//  comiso.solver().misolver().set_gurobi_max_time(120);
#if COMISO_GUROBI_AVAILABLE
  COMISO::GUROBISolver gb;
//  gb.set_GRB_MIPFocus(1); // >> If you are more interested in finding feasible solutions quickly, you can select MIPFocus=1.
//  gb.set_GRB_Heuristics(0.20);
#endif

  // stiffening
//  std::vector<PairUiV> discrete_constraints2;
  StopWatch sw;
  sw.start();
  std::vector<double> best_found = fe_problem.x();
  int n_best_found = mesh_.n_cells();
  ALGOHEX_DEBUG_ONLY(std::cerr << "run stiffening with weight = " << stiffening_weight_ << " and #max iters = "
                               << max_stiffening_iters_ << std::endl;)
  for (unsigned int i = 0; i < max_stiffening_iters_; ++i)
  {
    if (opt_type_ == OPT_COMISO)
      comiso.solve(&fe_problem, constraint_pointers, discrete_constraints, 0.0, false);
    else if (opt_type_ == OPT_GUROBI)
    {
#if COMISO_GUROBI_AVAILABLE
      gb.solve(&fe_problem, constraint_pointers, discrete_constraints, timelimit_);
#else
      std::cerr << "ERROR: Gurobi solver selected, but no Gurobi support compiled in." << std::endl;
      return;
#endif
    }

    int n_invalid = update_cell_weights(fe_problem.x());
    if (n_invalid < n_best_found)
    {
      best_found = fe_problem.x();
      n_best_found = n_invalid;
    }
    ALGOHEX_DEBUG_ONLY(std::cerr << "stiffening iter " << i << " #invalid: " << n_invalid << std::endl;)
    if (n_invalid == 0)
      break;
    else
      setup_frame_objective_function(fe_ff, _anisotropy_alpha, sizing_scale_factor);
  }
  // restore best solution
  fe_problem.x() = best_found;
  ALGOHEX_DEBUG_ONLY(std::cerr << "best stiffening solution has: " << n_best_found << " degenerate or inverted cells"
                               << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "time stiffening: " << sw.stop() / 1000.0 << "s" << std::endl;)


  // solve
//  if(opt_type_ == OPT_COMISO)
//  {
//    comiso.solve(&fe_problem, constraint_pointers, discrete_constraints, 0.0, false, true);
//  }
//  else
//    if(opt_type_ == OPT_GUROBI)
//  {
//#if COMISO_GUROBI_AVAILABLE
//    COMISO::GUROBISolver gb;
//    gb.solve(&fe_problem, constraint_pointers, discrete_constraints, timelimit_);
//#endif
//  }

  // make integer exact
  for (unsigned int i = 0; i < discrete_constraints.size(); ++i)
  {
    if (discrete_constraints[i].second == COMISO::Integer)
    {
      double val = ROUND_MI(fe_problem.x()[discrete_constraints[i].first]);
      if (std::abs(val - fe_problem.x()[discrete_constraints[i].first]) > 1e-3)
        std::cerr << "Warning: found integer variable with large deviation!!! "
                  << std::abs(val - fe_problem.x()[discrete_constraints[i].first]) << std::endl;
      fe_problem.x()[discrete_constraints[i].first] = val;
    }
  }

  // store result to cell property igm_cprop_
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    int i = 0;
    for (TVIt tv_it = mesh_.tv_iter(*c_it); tv_it.valid(); ++tv_it)
    {
      int bidx = 3 * global_idx_[*c_it][i++];
      igm_cprop_[*c_it][*tv_it] = Point(fe_problem.x()[bidx], fe_problem.x()[bidx + 1],
                                        fe_problem.x()[bidx + 2]);
    }
  }

  // check output
  double final_energy = fe_problem.eval_f(fe_problem.x().data());
  double param_vol = parametric_volume();
  int n_invalid_tets = number_of_invalid_parametric_tets();
  int n_invalid_valencies = number_of_invalid_edge_valencies();
  double dvm = degenerate_volume_map();
  double mv = mesh_volume();

  ALGOHEX_DEBUG_ONLY(std::cerr << "with integers            : " << int(_with_integer_constraints) << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "final energy             : " << final_energy << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "parametric volume        : " << param_vol << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "#invalid parametric tets : " << n_invalid_tets << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "#invalid edge valencies  : " << n_invalid_valencies << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "degenerate map volume    : " << 100.0 * dvm / mv << "%"
                               << std::endl;)

  if (_with_integer_constraints)
  {
    json_data_["final_energy"] = final_energy;
    json_data_["parametric_volume"] = param_vol;
    json_data_["n_invalid_param_tets"] = n_invalid_tets;
    json_data_["n_invalid_valencies"] = n_invalid_valencies;
    json_data_["valid_volume"] = 1.0 - dvm / mv;
  }
  else
  {
    json_data_["final_energy_seamless"] = final_energy;
    json_data_["parametric_volume_seamless"] = param_vol;
    json_data_["n_invalid_param_tets_seamless"] = n_invalid_tets;
    json_data_["n_invalid_valencies_seamless"] = n_invalid_valencies;
    json_data_["valid_volume_seamless"] = 1.0 - dvm / mv;
  }

  // preserve solution vector for potential subsequent optimization
  x_prev_.swap(fe_problem.x());
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
optimize_integer_grid_map_locally_injective(const int _num_hex_cells, const double _anisotropy_alpha,
                                            const double _epsilon, const bool _with_integer_constraints)
{
  // compute target edge length
  double sizing_scale_factor = std::cbrt(frame_field_volume() / _num_hex_cells);
  double sizing_scale_factor_deform = std::cbrt(mesh_volume() / _num_hex_cells);
  ALGOHEX_DEBUG_ONLY(std::cerr << "sizing scale factor: " << sizing_scale_factor << std::endl;)
  // create multiple copies of vertices at cut and store indices in global_idx_
  setup_vertex_indices();
  check_vertex_indices();

  // setup objective function as finite-element problem
  reset_cell_weights();
  COMISO::FiniteElementSet<FrameFittingElement3D> fe_ff("Frame Fitting");
  setup_frame_objective_function(fe_ff, _anisotropy_alpha, sizing_scale_factor);

  COMISO::FiniteElementSet<FF3DV_PH_AD> fe_deform("Foldover Free");
  setup_deformation_objective_function(fe_deform, _epsilon, 1.0 / sizing_scale_factor_deform);

  const double w_sdirichlet = 100.0;
  const double inversion_penalty = 1e3;
  COMISO::FiniteElementSet<SDE3D_PH_AD> fe_sdirichlet("Symmetric Dirichlet");
  setup_symmetric_dirichlet_objective_function(fe_sdirichlet, w_sdirichlet, sizing_scale_factor_deform,
                                               inversion_penalty);

  // dihedral angle support is experimental and deactivated by default
  const double w_dihedral = 0.0;
  const double s_dihedral = 4.0;
  COMISO::FiniteElementSet<DihedralAngleElement_TAD_PH> fe_dihedral("Dihedral");
  if (w_dihedral > 0.0)
    setup_dihedral_angle_objective_function(fe_dihedral, w_dihedral, s_dihedral);

  const double w_regularize_helpers = 1e-6;
  COMISO::FiniteElementSet<LinearLeastSquaresElement1D> fe_regularize_helpers("Helper Regularization");

  // setup constraints
  std::vector<PairUiV> discrete_constraints;
  std::vector<COMISO::LinearConstraint> constraints;
  std::vector<COMISO::NConstraintInterface *> constraint_pointers;
  setup_constraints(constraints, constraint_pointers, discrete_constraints, 0, _with_integer_constraints);


  if (!_with_integer_constraints)
    discrete_constraints.clear();

  // solve and store result
  unsigned int n_param = 3 * (max_global_idx_ + 1);
  unsigned int n_complete = n_param + n_helper_variables_;
  COMISO::FiniteElementProblem fe_problem(n_complete), fe_problem_complete(n_complete);

  // dummy problem to showe energies
  fe_problem_complete.add_set(&fe_ff);
  fe_problem_complete.add_set(&fe_deform);
  fe_problem_complete.add_set(&fe_sdirichlet);
  fe_problem_complete.add_set(&fe_dihedral);
  fe_problem_complete.add_set(&fe_regularize_helpers);

  // optimization problem
//  fe_problem.add_set(&fe_ff);
//  fe_problem.add_set(&fe_deform);
  fe_problem.add_set(&fe_sdirichlet);
  fe_problem.add_set(&fe_regularize_helpers);
//  if(w_dihedral > 0.0)
//    fe_problem.add_set(&fe_dihedral);

  // set initial solution
  if (x_prev_.size() < n_complete)
  {
    std::cerr << "ERROR: dimension of previous solution is " << x_prev_.size() << " but expected " << n_complete
              << " -----> abort" << std::endl;
    return;
  }
  else
  {
    x_prev_.resize(n_complete);
    fe_problem.x() = x_prev_;
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "--- initial energies complete --- \n";)
  fe_problem_complete.x() = fe_problem.x();
  fe_problem_complete.print_objectives();

  ALGOHEX_DEBUG_ONLY(std::cerr << "--- initial energies selected --- \n";)
  fe_problem.print_objectives();
  ALGOHEX_DEBUG_ONLY(std::cerr << "-------------------------------\n";)


  // fix discrete constraints to their current value (assume that the integer problem was solved before)
  if (_with_integer_constraints)
  {
    std::set<int> already_added;
    for (auto dc: discrete_constraints)
    {
      int var_idx = dc.first;
      if (already_added.insert(var_idx).second)
      {
        COMISO::LinearConstraint::SVectorNC coeffs(n_complete);
        coeffs.coeffRef(var_idx) = 1.0;
        constraints.push_back(
                COMISO::LinearConstraint(coeffs, -fe_problem.x()[var_idx], COMISO::NConstraintInterface::NC_EQUAL));
      }
    }

    // fix constraint pointers
    constraint_pointers.clear();
    constraint_pointers.reserve(constraints.size());
    for (unsigned int i = 0; i < constraints.size(); ++i)
      constraint_pointers.push_back(&(constraints[i]));

    // ToDo: regularize helper variables
    for (int k = 0; k < n_helper_variables_; ++k)
    {
      LinearLeastSquaresElement1D::VecI i;
      LinearLeastSquaresElement1D::VecC c;
      i[0] = n_param + k;
      c[0] = 1.0;
      c[1] = 0.0;
      c[2] = 1.0;
      fe_regularize_helpers.instances().add_element(i, c);
    }

  }

  // make constraints linearly independent
  COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraint_pointers);

  // convert constraints
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  COMISO::LinearConstraintConverter::nsolver_to_eigen(constraint_pointers, A, b, n_complete);

  // solve
  COMISO::NPTiming fe_problem_timings(&fe_problem); // with timings

  // solve via TruncatedNewtonPCG
  double tol = 1e-2; // ToDo: remove scale dependency
  COMISO::TruncatedNewtonPCG tn(tol);
  tn.always_update_preconditioner() = true;
//  tn.adaptive_tolerance_modifier() = 0.01;
//  tn.max_iters()     = 100;
//  tn.max_pcg_iters() = 300;
  tn.max_iters() = 30;
  tn.max_pcg_iters() = 300;

  // memorize best solution
  int best_nit = mesh_.n_cells();
  int best_nie = mesh_.n_edges();
  std::vector<double> best_x(n_complete);

  // allow several rounds of penalty updates
  for (int i = 0; i < 3; ++i)
  {
    tn.solve_projected_normal_equation(&fe_problem_timings, A, b);
    store_igm(fe_problem.x().data());

    int nit = number_of_invalid_parametric_tets();
    int nie = number_of_invalid_edge_valencies();

    // better than previous solutions?
    if (nit < best_nit && nie <= best_nie)
    {
      best_nit = nit;
      best_nie = nie;
      best_x = fe_problem.x();
    }

    if (nit == 0)
      break;

    // stop if numerics break due to penalty
    if (!std::isfinite(fe_problem.eval_f(fe_problem.x().data())))
      break;

    double penalty_l = 100.0; // corresponds to a factor of 100^2 in the symmetric dirichlet energy

    int n_nonzero_regularizers = 0;
    for (unsigned int j = 0; j < fe_sdirichlet.instances().size(); ++j)
    {
      // update regularizer
      SDE3D_PH_AD::VecV xl;
      for (int k = 0; k < xl.innerSize(); ++k)
        xl[k] = fe_problem.x()[fe_sdirichlet.instances().index(j, k)];

      SDE3D_PH_AD::VecC cl = fe_sdirichlet.instances().c(j);

      double eps_l = SDE3D_PH_AD::determine_suitable_epsilon(xl, cl, penalty_l);
      fe_sdirichlet.instances().c(j)[2] = eps_l;

      if (eps_l > 0.0)
        ++n_nonzero_regularizers;
    }

    ALGOHEX_DEBUG_ONLY(std::cerr << "----- penalty iter      = " << i <<
                                 "\n penalty               = " << penalty_l <<
                                 "\n #invalid tets         = " << nit <<
                                 "\n #invalid edges        = " << nie <<
                                 "\n #nonzero_regularizers = " << n_nonzero_regularizers <<
                                 "\n -----" << std::endl;)
  }

  // recover best solution found in-between
  fe_problem.x() = best_x;

  // Newton
//  COMISO::NewtonSolver ns(1e-4, 1e-9, 2000);
//  ns.solve_infeasible_start(&fe_problem_timings, constraint_pointers);

//  try {
//    COMISO::IPOPTSolver ipsol;
////    ipsol.set_ipopt_option("print_level", 0);
//    ipsol.set_max_iterations(100);
//    ipsol.set_ipopt_option("tol", 1e-4);
////    ipsol.set_ipopt_option("hessian_approximation", "limited-memory");
////    ipsol.set_ipopt_option("limited_memory_max_history", 20);
//    ipsol.solve(&fe_problem_timings, constraint_pointers);
//  }
//  catch (...)
//  {
//    std::cerr << "ERROR: IPOPT threw an exception!!!!!" << std::endl;
//  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "--- final energies complete --- \n";)
  fe_problem_complete.x() = fe_problem.x();
  fe_problem_complete.print_objectives();

  ALGOHEX_DEBUG_ONLY(std::cerr << "--- final energies selected --- \n";)
  fe_problem.print_objectives();
  ALGOHEX_DEBUG_ONLY(std::cerr << "-------------------------------\n";)


  // make integer exact
  for (unsigned int i = 0; i < discrete_constraints.size(); ++i)
  {
    if (discrete_constraints[i].second == COMISO::Integer)
    {
      double val = ROUND_MI(fe_problem.x()[discrete_constraints[i].first]);
      if (std::abs(val - fe_problem.x()[discrete_constraints[i].first]) > 1e-3)
        std::cerr << "Warning: found integer variable with large deviation!!! "
                  << std::abs(val - fe_problem.x()[discrete_constraints[i].first]) << std::endl;
      fe_problem.x()[discrete_constraints[i].first] = val;
    }
  }

  // store result to cell property igm_cprop_
  store_igm(fe_problem.x().data());

  // check output
  double final_energy = fe_problem.eval_f(fe_problem.x().data());
  double param_vol = parametric_volume();
  int n_invalid_tets = number_of_invalid_parametric_tets();
  int n_invalid_valencies = number_of_invalid_edge_valencies();
  double dvm = degenerate_volume_map();
  double mv = mesh_volume();

  ALGOHEX_DEBUG_ONLY(std::cerr << "with integers            : " << int(_with_integer_constraints) << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "final energy             : " << final_energy << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "parametric volume        : " << param_vol << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "#invalid parametric tets : " << n_invalid_tets << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "#invalid edge valencies  : " << n_invalid_valencies << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "degenerate map volume    : " << 100.0 * dvm / mv << "%"
                               << std::endl;)

  if (_with_integer_constraints)
  {
    json_data_["final_energy"] = final_energy;
    json_data_["parametric_volume"] = param_vol;
    json_data_["n_invalid_param_tets"] = n_invalid_tets;
    json_data_["n_invalid_valencies"] = n_invalid_valencies;
    json_data_["valid_volume"] = 1.0 - dvm / mv;
  }
  else
  {
    json_data_["final_energy_seamless"] = final_energy;
    json_data_["parametric_volume_seamless"] = param_vol;
    json_data_["n_invalid_param_tets_seamless"] = n_invalid_tets;
    json_data_["n_invalid_valencies_seamless"] = n_invalid_valencies;
    json_data_["valid_volume_seamless"] = 1.0 - dvm / mv;
  }

  // preserve solution vector for potential subsequent optimization
  x_prev_.swap(fe_problem.x());
}



//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
optimize_integer_grid_map_robust_quantization(const int _num_hex_cells)
{
  // compute target edge length
  double sizing_scale_factor = std::cbrt(frame_field_volume() / _num_hex_cells);
  double sizing_scale_factor_deform = std::cbrt(mesh_volume() / _num_hex_cells);
  std::cerr << "sizing scale factor: " << sizing_scale_factor << std::endl;
  // create multiple copies of vertices at cut and store indices in global_idx_
  setup_vertex_indices();
  check_vertex_indices();

  const double w_sdirichlet = 100.0;
  const double inversion_penalty = 1e3;
  COMISO::FiniteElementSet<SDE3D_PH_AD> fe_sdirichlet("Symmetric Dirichlet");
  setup_symmetric_dirichlet_objective_function(fe_sdirichlet, w_sdirichlet, sizing_scale_factor_deform,
                                               inversion_penalty);

  // De-activated by default for now
  const double w_jacobian_smoothness = 0.0;
  COMISO::FiniteElementSet<JSE3D_AD> fe_jacobian_smoothness("Jacobian Smoothness");
  if(w_jacobian_smoothness > 0.0)
    setup_jacobian_smoothness_objective_function(fe_jacobian_smoothness, w_jacobian_smoothness, sizing_scale_factor_deform);


  // setup constraints
  std::vector<PairUiV> discrete_constraints;
  std::vector<COMISO::LinearConstraint> constraints;
  std::vector<COMISO::NConstraintInterface *> constraint_pointers;
  setup_constraints(constraints, constraint_pointers, discrete_constraints, 0, false);

  // solve and store result
  unsigned int n_param = 3 * (max_global_idx_ + 1);
  unsigned int n_complete = n_param + n_helper_variables_; // identical to n_param in this version
  COMISO::FiniteElementProblem fe_problem(n_complete);

  fe_problem.add_set(&fe_sdirichlet);
  if(w_jacobian_smoothness > 0.0)
    fe_problem.add_set( &fe_jacobian_smoothness);

  // set initial solution
  if (x_prev_.size() < n_param)
  {
    std::cerr << "ERROR: dimension of previous solution is " << x_prev_.size() << " but expected at least " << n_param
              << " -----> abort" << std::endl;
    return;
  }
  else
  {
    for (unsigned int i = 0; i < n_param; ++i)
      fe_problem.x()[i] = x_prev_[i];
  }

//  // check constraint violation
//  for(unsigned int i=0; i<constraints.size(); ++i)
//  {
//    double v = constraints[i].eval_constraint(fe_problem.x().data());
//    double rhs = -constraints[i].b();
//
//    std::cerr << i << ", violation=" << v << ", rhs=" << rhs << std::endl;
//  }

  std::cerr << "--- initial energies --- \n";
  fe_problem.print_objectives();
  std::cerr << "-------------------------\n";

  // make constraints linearly independent
  std::cerr << "#constraints = " << constraint_pointers.size() << std::endl;
  COMISO::ConstraintTools::remove_dependent_linear_constraints_only_linear_equality(constraint_pointers);
  std::cerr << "#independent constraints = " << constraint_pointers.size() << std::endl;

  // convert constraints
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd b;
  COMISO::LinearConstraintConverter::nsolver_to_eigen(constraint_pointers, A, b, n_complete);

  // solve
  COMISO::NPTiming fe_problem_timings(&fe_problem); // with timings

//  // Newton
//  COMISO::NewtonSolver ns(1e-4, 1e-9, 2000);
//  ns.solve_infeasible_start(&fe_problem_timings, constraint_pointers);

  try
  {
    COMISO::IPOPTSolver ipsol;
//    ipsol.set_ipopt_option("print_level", 0);
    ipsol.set_max_iterations(200);
    ipsol.set_ipopt_option("tol", 1e-4);
//    ipsol.set_ipopt_option("hessian_approximation", "limited-memory");
//    ipsol.set_ipopt_option("limited_memory_max_history", 20);
    ipsol.solve(&fe_problem_timings, constraint_pointers);
  }
  catch (...)
  {
    std::cerr << "ERROR: IPOPT threw an exception!!!!!" << std::endl;
  }


//  // solve via TruncatedNewtonPCG
//  double tol = 1e-2; // ToDo: remove scale dependency
//  COMISO::TruncatedNewtonPCG tn(tol);
//  tn.always_update_preconditioner() = true;
////  tn.adaptive_tolerance_modifier() = 0.01;
////  tn.max_iters()     = 100;
////  tn.max_pcg_iters() = 300;
//  tn.max_iters()     = 30;
//  tn.max_pcg_iters() = 300;
//  tn.solve_projected_normal_equation(&fe_problem_timings, A, b);
//  store_igm(fe_problem.x().data());


  std::cerr << "--- final energies --- \n";
  fe_problem.print_objectives();
  std::cerr << "-------------------------------\n";


  // make integer exact
  for (unsigned int i = 0; i < discrete_constraints.size(); ++i)
  {
    if (discrete_constraints[i].second == COMISO::Integer)
    {
      double val = ROUND_MI(fe_problem.x()[discrete_constraints[i].first]);
      if (std::abs(val - fe_problem.x()[discrete_constraints[i].first]) > 1e-3)
        std::cerr << "Warning: found integer variable with large deviation!!! "
                  << std::abs(val - fe_problem.x()[discrete_constraints[i].first]) << std::endl;
      fe_problem.x()[discrete_constraints[i].first] = val;
    }
  }

  // store result to cell property igm_cprop_
  store_igm(fe_problem.x().data());

  // check output
  double final_energy = fe_problem.eval_f(fe_problem.x().data());
  double param_vol = parametric_volume();
  int n_invalid_tets = number_of_invalid_parametric_tets();
  int n_invalid_valencies = number_of_invalid_edge_valencies();
  double dvm = degenerate_volume_map();
  double mv = mesh_volume();

  std::cerr << "final energy             : " << final_energy << std::endl;
  std::cerr << "parametric volume        : " << param_vol << std::endl;
  std::cerr << "#invalid parametric tets : " << n_invalid_tets << std::endl;
  std::cerr << "#invalid edge valencies  : " << n_invalid_valencies << std::endl;
  std::cerr << "degenerate map volume    : " << 100.0 * dvm / mv << "%"
            << std::endl;

  json_data_["final_energy"] = final_energy;
  json_data_["parametric_volume"] = param_vol;
  json_data_["n_invalid_param_tets"] = n_invalid_tets;
  json_data_["n_invalid_valencies"] = n_invalid_valencies;
  json_data_["valid_volume"] = 1.0 - dvm / mv;

  // preserve solution vector for potential subsequent optimization
  x_prev_.swap(fe_problem.x());
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
store_igm(const double *_x)
{
// store result to cell property igm_cprop_
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    int i = 0;
    for (TVIt tv_it = mesh_.tv_iter(*c_it); tv_it.valid(); ++tv_it)
    {
      int bidx = 3 * global_idx_[*c_it][i++];
      igm_cprop_[*c_it][*tv_it] = Point(_x[bidx], _x[bidx + 1], _x[bidx + 2]);
    }
  }
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_vertex_indices()
{
  // reset indices
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
    global_idx_[*c_it] = Vec4i(-1, -1, -1, -1);

  int cur_idx(0);
  for (VIt v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
  {
    for (VCIt vc_it = mesh_.vc_iter(*v_it); vc_it.valid(); ++vc_it)
    {
      int tvidx = tet_vidx(*vc_it, *v_it);
      // error check
      if (tvidx == -1) std::cerr << "ERROR: invalid tvidx" << std::endl;

      // if not set yet, start a new region
      if (global_idx_[*vc_it][tvidx] == -1)
      {
        global_idx_[*vc_it][tvidx] = cur_idx;

        // try to propagate through incident halffaces
        std::queue<HFH> qhfh;
        // get halfface opposite to vertex in cell
        HFH hfh_vopp = opposite_halfface_handle(*vc_it, *v_it);
        for (unsigned int i = 0; i < mesh_.halfface(hfh_vopp).halfedges().size(); ++i)
        {
          // get halfedge
          HEH heh = mesh_.halfface(hfh_vopp).halfedges()[i];
          // push adjacent halfface (incident to *v_it)
          qhfh.push(mesh_.adjacent_halfface_in_cell(hfh_vopp, heh));
        }

        while (!qhfh.empty())
        {
          HFH hfh = qhfh.front();
          FH fh = mesh_.face_handle(hfh);
          qhfh.pop();

          if (!mesh_.is_boundary(fh) && !is_on_cut_fprop_[fh])
          {
            // get cell handle
            CH ch_new = mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh));
            // get local vertex index
            int tvidx_new = tet_vidx(ch_new, *v_it);

            // error check
            if (tvidx_new == -1) std::cerr << "ERROR: tvidx_new is invalid" << std::endl;

            // already visited?
            if (tvidx_new != -1 && global_idx_[ch_new][tvidx_new] == -1)
            {
              // set current index
              global_idx_[ch_new][tvidx_new] = cur_idx;

              //             std::cerr << "set index of cell " << ch_new.idx() << " vertex " << tvidx_new << std::endl;

              // push halffaces incident to *v_it in ch_new
              hfh_vopp = opposite_halfface_handle(ch_new, *v_it);
              for (unsigned int i = 0; i < mesh_.halfface(hfh_vopp).halfedges().size(); ++i)
              {
                // get halfedge
                HEH heh = mesh_.halfface(hfh_vopp).halfedges()[i];
                // push adjacent halfface (incident to *v_it)
                qhfh.push(mesh_.adjacent_halfface_in_cell(hfh_vopp, heh));
              }
            }
          }
        }
        // move on to next index
        ++cur_idx;
      }
    }
  }

  // store index range
  max_global_idx_ = cur_idx - 1;

  ALGOHEX_DEBUG_ONLY(std::cerr << "max_global_idx_ = " << max_global_idx_ << " (#vertices " << mesh_.n_vertices() << ")"
                               << std::endl;)
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
check_vertex_indices() const
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "check vertex indices...";)

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
    for (unsigned int i = 0; i < 4; ++i)
      if (global_idx_[*c_it][i] == -1)
        std::cerr << "ERROR: uninitialized indices found!" << std::endl;

  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (!mesh_.is_boundary(*f_it))
    {
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);
      CH ch0 = mesh_.incident_cell(hfh0);
      CH ch1 = mesh_.incident_cell(hfh1);

      VH vh0 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[0]).to_vertex();
      VH vh1 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[1]).to_vertex();
      VH vh2 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[2]).to_vertex();

      if (!is_on_cut_fprop_[*f_it])
      {
        if (global_idx_[ch0][tet_vidx(ch0, vh0)] != global_idx_[ch1][tet_vidx(ch1, vh0)])
          std::cerr << "ERROR: different global index despite no cutface!!!" << std::endl;
        if (global_idx_[ch0][tet_vidx(ch0, vh1)] != global_idx_[ch1][tet_vidx(ch1, vh1)])
          std::cerr << "ERROR: different global index despite no cutface!!!" << std::endl;
        if (global_idx_[ch0][tet_vidx(ch0, vh2)] != global_idx_[ch1][tet_vidx(ch1, vh2)])
          std::cerr << "ERROR: different global index despite no cutface!!!" << std::endl;
      }
    }

  ALGOHEX_DEBUG_ONLY(std::cerr << " done" << std::endl;)
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_frame_objective_function(COMISO::FiniteElementSet<FrameFittingElement3D> &_fe_ff, const double _alpha,
                               const double _sizing_scale_factor) const
{
  // clear old data
  _fe_ff.instances().clear_elements();

  COMISO::FiniteElementSet<FrameFittingElement3D>::VecI idx;
  COMISO::FiniteElementSet<FrameFittingElement3D>::VecC c;

  // add one element per tetrahedron
  Eigen::Matrix<double, 3, 4> P;
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // collect points
    unsigned int i = 0;
    for (TVIt tv_it = mesh_.tv_iter(*c_it); tv_it.valid(); ++tv_it)
      P.col(i++) = ovm2eigen(mesh_.vertex(*tv_it));

    // set indices
    for (i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 4; ++j)
        idx[4 * i + j] = 3 * global_idx_[*c_it][j] + i;

    // get frame
//            Eigen::Matrix3d F = _target_edgelength * quaternion_cprop_[*c_it].toRotationMatrix();
    // ToDo: handle target_edgelength correctly
    Mat3x3 F = _sizing_scale_factor * frame_cprop_[*c_it];

//    std::cerr << "points:" << std::endl << P << std::endl;
//    std::cerr << "frame :" << std::endl << F << std::endl;

    // construct constants
    FrameFittingElement3D::constants_from_tetrahedron_and_frame(P, F, _alpha, c);

    // augment with cell weighting
    c[13] *= cell_weight_[*c_it];

    // add instance to set
    _fe_ff.instances().add_element(idx, c);
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_deformation_objective_function(COMISO::FiniteElementSet<FF3DV_PH_AD> &_fe_ff, const double _epsilon,
                                     const double _target_det_scaling)
{
  // clear old data
  _fe_ff.instances().clear_elements();

//  // collect tets incident to invalid edges ---> these should always be allowed to flip teomporarily
//  std::set<CH> ch_next_to_invalid;
//  std::set<CH> ch_invalid;
//  for(auto eh : mesh_.edges())
//  {
//    if(param_invalid_edge_valence_valid_tets(eh))
//    {
//      ECIt ec_it(eh, &mesh_);
//      for (; ec_it.valid(); ++ec_it)
//        ch_next_to_invalid.insert(*ec_it);
//    }
//  }

//  std::cerr << "#cells next to invalid edges = " << ch_next_to_invalid.size() << std::endl;

  COMISO::FiniteElementSet<FF3DV_PH_AD>::VecI idx;
  COMISO::FiniteElementSet<FF3DV_PH_AD>::VecC c;

  // add one element per tetrahedron
  Eigen::Matrix<double, 3, 4> P, U;
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // set indices
    for (unsigned int j = 0; j < 4; ++j)
      for (int i = 0; i < 3; ++i)
        idx[3 * j + i] = 3 * global_idx_[*c_it][j] + i;

    // collect points
    mesh_tet(*c_it, P);

    // check validity of tet
    parametric_tet(*c_it, U);
    // get edge matrix
    Mat3x3 E;
    E.col(0) = U.col(1) - U.col(0);
    E.col(1) = U.col(2) - U.col(0);
    E.col(2) = U.col(3) - U.col(0);

    double det = E.determinant();

    // determine parameters
    double epsilon = _epsilon;
    double w_conformal = 1.0;
    double w_det = 1e-3;
    double target_det_scaling = _target_det_scaling;

    // increase weighting of bad quality tets because they are likely to degenerate
//    if(tetrahedron_quality(P.col(0), P.col(1), P.col(2), P.col(3)) < 0.1 || tetrahedron_quality(U.col(0), U.col(1), U.col(2), U.col(3)) < 0.1)

    if (0)
      if (tetrahedron_quality(P.col(0), P.col(1), P.col(2), P.col(3)) < 0.1)
      {
        Vec3d q0, q1, q2, q3;
        construct_uniform_tetrahedron(q0, q1, q2, q3,
                                      avg_tetrahedron_edge_length(P.col(0), P.col(1), P.col(2), P.col(3)));

        P.col(0) = q0;
        P.col(1) = q1;
        P.col(2) = q2;
        P.col(3) = q3;

//      std::cerr << "quality of uniform tetrahedron = " << tetrahedron_quality(P.col(0), P.col(1), P.col(2), P.col(3)) << std::endl;
//      w_conformal *= 10.0;
//      w_det       *= 10.0;
      }

    // determine target determinant (for deformation energy)
    double target_det = std::pow(target_det_scaling, 3);

    // construct constants
    FF3DV_PH_AD::compute_constants(w_conformal, w_det, target_det, epsilon, P.col(0), P.col(1), P.col(2), P.col(3), c);

    // adapt epsilon
    if (det > 0.0)
    {
      c[3] = 0;
    }
    else
    {
//      // memorize invalid
//      ch_invalid.insert(*c_it);

      // determine parameters
      epsilon = _epsilon;
      w_conformal = 1.0;
      w_det = 1e-3;
      target_det = 1.0 * target_det;
      // construct constants
      FF3DV_PH_AD::compute_constants(w_conformal, w_det, target_det, epsilon, P.col(0), P.col(1), P.col(2), P.col(3),
                                     c);

      FF3DV_PH_AD::VecV x;

      for (unsigned int i = 0; i < x.innerSize(); ++i)
        x[i] = x_prev_[idx[i]];

      double eps_heuristic = FF3DV_PH_AD::determine_suitable_epsilon(x, c);


      double eps_adaptive = 1.0;
      FF3DV_PH_AD ff3dv;
      FF3DV_PH_AD::VecC c2 = c;
      // normalized element without volume weighting
      c2[0] = w_conformal;
      c2[1] = w_det;
      c2[3] = eps_adaptive;

      for (unsigned int i = 0; i < 50; ++i)
      {
        if (ff3dv.eval_f(x, c2) < 1e2 / _epsilon) // become more and more strict with each iteration
        {
          eps_adaptive *= 0.5;
          c2[3] = eps_adaptive;
        }
        else break;
      }

      // choose adaptive epsilon unless it is larger than global epsilon
      c[3] = std::min(eps_adaptive, _epsilon);

      FF3DV_PH_AD::VecC c3(c2), c4(c2);
      c3[3] = eps_heuristic;
      c4[3] = _epsilon;

//      std::cerr << "adaptive epsilon = "    << eps_adaptive   << " (energy=" << ff3dv.eval_f(x, c2) << ")"
//                << ", heuristic epsilon = " << eps_heuristic  << " (energy=" << ff3dv.eval_f(x, c3) << ")"
//                << ", fixed epsilon = "     << _epsilon       << " (energy=" << ff3dv.eval_f(x, c4) << ")"
//                << std::endl;


    }

//    if(ch_next_to_invalid.count(*c_it) && !ch_invalid.count(*c_it))
//    {
//      //  const static int NC = 13; // [weight_conformal, weight_det, target_det, epsilon, g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z]
////      c[0] *= 1000.0;
////      c[1] *= 1000.0;
//      c[3] = 0.1;
//    }

    // add instance to set
    _fe_ff.instances().add_element(idx, c);
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_symmetric_dirichlet_objective_function(COMISO::FiniteElementSet<SDE3D_PH_AD> &_fe_sdirichlet, const double _w,
                                             const double _sizing_scale, const double _inversion_penalty)
{
  // clear old data
  _fe_sdirichlet.instances().clear_elements();

  // calculate mesh volume for normalization of energy
  double VM = mesh_volume();
  double penalty = std::sqrt(_inversion_penalty);

  COMISO::FiniteElementSet<SDE3D_PH_AD>::VecI idx;
  COMISO::FiniteElementSet<SDE3D_PH_AD>::VecC c;

  // add one element per tetrahedron
  Eigen::Matrix<double, 3, 4> P, U;
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // set indices
    for (unsigned int j = 0; j < 4; ++j)
      for (int i = 0; i < 3; ++i)
        idx[3 * j + i] = 3 * global_idx_[*c_it][j] + i;

    // collect points
    mesh_tet(*c_it, P);

    // check validity of tet
    parametric_tet(*c_it, U);

    // normalixed weight
    double w = _w * get_cell_volume(mesh_, *c_it) / VM;

    // construct constants
    SDE3D_PH_AD::compute_constants(w, _sizing_scale, 0.0, P.col(0), P.col(1), P.col(2), P.col(3), c);

    Eigen::Map<SDE3D_PH_AD::VecV> x(U.data());
    c[2] = SDE3D_PH_AD::determine_suitable_epsilon(x, c, _inversion_penalty);

    // add instance to set
    _fe_sdirichlet.instances().add_element(idx, c);
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_jacobian_smoothness_objective_function( COMISO::FiniteElementSet<JSE3D_AD>& _fe_jacobian_smoothness, const double _w, const double _sizing_scale)
{
  // clear old data
  _fe_jacobian_smoothness.instances().clear_elements();

  // calculate mesh volume for normalization of energy
  double  VM      = mesh_volume();
  double cVM      = std::cbrt(VM);

  COMISO::FiniteElementSet<JSE3D_AD>::VecI idx;
  COMISO::FiniteElementSet<JSE3D_AD>::VecC c;

  // add one element per interior face
  Eigen::Matrix<double, 3, 4> P0, P1, U0, U1;
  for(auto fh : mesh_.faces())
    if(!mesh_.is_boundary(fh))
    {
      auto hfh0 = mesh_.halfface_handle(fh,0);
      auto hfh1 = mesh_.halfface_handle(fh,1);

      auto ch0 = mesh_.incident_cell(hfh0);
      auto ch1 = mesh_.incident_cell(hfh1);

      // set indices
      for (unsigned int j = 0; j < 4; ++j)
        for (int i = 0; i < 3; ++i)
        {
          idx[     3 * j + i] = 3 * global_idx_[ch0][j] + i;
          idx[12 + 3 * j + i] = 3 * global_idx_[ch1][j] + i;
        }

      // collect points
      mesh_tet(ch0, P0);
      mesh_tet(ch1, P1);

//      // check validity of tet
//      parametric_tet(ch0, U0);
//      parametric_tet(ch1, U1);

      // transition function ( for Jacobians is inverse to the one of frames)
      Eigen::Matrix3d T01 = tq_.transition_matrix_int(tq_.inverse_transition_idx(transition_hfprop_[hfh0]));

      // compute orthogonal distance of barycenters w.r.t. fh
      double d = std::abs((mesh_.barycenter(ch0)-mesh_.barycenter(ch1)) | mesh_.normal(hfh0));

      // normalized weight
      double w = _w * 0.25*(get_cell_volume(mesh_, ch0)+get_cell_volume(mesh_, ch1)) / (d*d*cVM);
//      double w = _w * 0.25*(get_cell_volume(mesh_, ch0)+get_cell_volume(mesh_, ch1)) / VM;

      // construct constants
      JSE3D_AD::compute_constants(w,
                                  _sizing_scale,
                                  P0.col(0), P0.col(1), P0.col(2), P0.col(3),
                                  P1.col(0), P1.col(1), P1.col(2), P1.col(3),
                                  T01,
                                  c);

      // add instance to set
      _fe_jacobian_smoothness.instances().add_element(idx, c);
    }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_domain_deformation_objective_function(COMISO::FiniteElementSet<FF3DV_PH_AD> &_fe_ff, const int _domain_idx_base)
{
  // clear old data
  _fe_ff.instances().clear_elements();

  // determine parameters
  double epsilon = 0.0;
  double w_conformal = 10.0;
  double w_det = 1e-2;
  double target_det = 1.0;

  COMISO::FiniteElementSet<FF3DV_PH_AD>::VecI idx;
  COMISO::FiniteElementSet<FF3DV_PH_AD>::VecC c;

  // add one element per tetrahedron
  Eigen::Matrix<double, 3, 4> P;
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    // set indices
    int j = 0;
    for (TVIt tv_it = mesh_.tv_iter(*c_it); tv_it.valid(); ++tv_it, ++j)
    {
      int domain_idx = _domain_idx_base + 3 * tv_it->idx();
      for (int i = 0; i < 3; ++i)
        idx[3 * j + i] = domain_idx + i;
    }

    // collect points
    mesh_tet(*c_it, P);

    // construct constants
    FF3DV_PH_AD::compute_constants(w_conformal, w_det, target_det, epsilon, P.col(0), P.col(1), P.col(2), P.col(3), c);

    // add instance to set
    _fe_ff.instances().add_element(idx, c);
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_ball_barrier_objective_function(COMISO::FiniteElementSet<BB3D_AD_PH> &_fe_bb, const int _domain_idx_base,
                                      const double _radius_scale)
{
  // clear old data
  _fe_bb.instances().clear_elements();

  COMISO::FiniteElementSet<BB3D_AD_PH>::VecI idx;
  COMISO::FiniteElementSet<BB3D_AD_PH>::VecC c;

  std::vector<HEH> fheh;
  std::vector<FH> ffh;

  for (auto vh: mesh_.vertices())
    if (incident_feature_ohalfedges(vh, fheh) || incident_feature_faces(vh, ffh))
    {
      Point p = mesh_.vertex(vh);
      double r = _radius_scale * avg_local_edge_length(mesh_, vh);

      int i = 3 * vh.idx() + _domain_idx_base;
      idx[0] = i;
      idx[1] = i + 1;
      idx[2] = i + 2;

      //   const static int NC = 5; // [weight, radius, cx, cy, cz]  with center c=[cx,cy,cz]
      c[0] = 1e-1;
      c[1] = r;
      c[2] = p[0];
      c[3] = p[1];
      c[4] = p[2];

      // add instance to set
      _fe_bb.instances().add_element(idx, c);
    }
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_dihedral_angle_objective_function(COMISO::FiniteElementSet<DihedralAngleElement_TAD_PH> &_fe_dihedral,
                                        const double _w_dihedral, const double _s_dihedral)
{
  // clear old data
  _fe_dihedral.instances().clear_elements();

  COMISO::FiniteElementSet<DihedralAngleElement_TAD_PH>::VecI idx;
  COMISO::FiniteElementSet<DihedralAngleElement_TAD_PH>::VecC c;
  c[0] = _w_dihedral;
  c[1] = _s_dihedral;

  for (auto eh: mesh_.edges())
//  if(mesh_.is_boundary(eh) || valence_eprop_[eh] != 0) // HACK
    if (param_invalid_edge_valence_valid_tets(eh))
    {
      // calculate target angle
      double target_angle = 2.0 * M_PI;
      if (mesh_.is_boundary(eh))
        target_angle = M_PI;
      target_angle += 0.5 * M_PI * valence_eprop_[eh];

//    bool valid_angle = false;
//    edge_angle
//    if(std::abs(edge_angle(eh)-target_angle) < 0.01)
//      valid_angle = true;

      // calculate mesh angle
      std::vector<double> tet_angles;
      double mesh_angle = mesh_edge_angle(eh, tet_angles);

      // clamp tet angles
      double min_valid_angle = 10.0 / 180.0 * M_PI;
      double max_valid_angle = 120.0 / 180.0 * M_PI;
      double mesh_angle_clamped = 0.0;
      for (double &a: tet_angles)
      {
        mesh_angle_clamped += clamp(a, min_valid_angle, max_valid_angle);

        std::cerr << "angle clamping " << a / M_PI * 180.0 << " -> "
                  << clamp(a, min_valid_angle, max_valid_angle) / M_PI * 180.0 << std::endl;
      }

      // scaling factor
      double scale_angle = target_angle / mesh_angle_clamped;

      // process all dihedral angles at edge
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
          Vec3d n0 = e0.cross(e1);
          n0 /= n0.norm();
          Vec3d n1 = e0.cross(e2);
          n1 /= n1.norm();
          double dp = n0.dot(n1);
          dp = std::max(-1.0, dp);
          dp = std::min(1.0, dp);
          double cp = n0.cross(n1).dot(e0) / e0.norm();
          cp = std::max(-1.0, cp);
          cp = std::min(1.0, cp);

          double dh = std::atan2(cp, dp);

          double dh_clamp = clamp(dh, min_valid_angle, max_valid_angle);

          double dh_final = dh_clamp * scale_angle;

          std::cerr << "final angle " << dh_final / M_PI * 180.0 << std::endl;

          c[2] = std::cos(dh_final);
          c[3] = std::sin(dh_final);

          c[0] = _w_dihedral * get_cell_volume(mesh_, ch) * 0.25;

          angle += dh_final;

          // get vertex indices in order vh0, vh1, vh00, vh11
          Vec4i gi = global_idx_[ch];

          int i0 = 3 * gi[tet_vidx(ch, vh0)];
          int i1 = 3 * gi[tet_vidx(ch, vh1)];
          int i2 = 3 * gi[tet_vidx(ch, vh00)];
          int i3 = 3 * gi[tet_vidx(ch, vh11)];

          idx << i0, i0 + 1, i0 + 2,
                  i1, i1 + 1, i1 + 2,
                  i2, i2 + 1, i2 + 2,
                  i3, i3 + 1, i3 + 2;

////        // use current angle if valid
////        valid_angle = false;
////        if(valid_angle)
////        {
//////          c[0] = _w_dihedral;
//////          c[1] = _s_dihedral;
////
////          c[0] = 1.0;
////          c[1] = 1.0;
////
////          Vec3d qt0(x_prev_[i0], x_prev_[i0+1], x_prev_[i0+2]);
////          Vec3d qt1(x_prev_[i1], x_prev_[i1+1], x_prev_[i1+2]);
////
////          Vec3d qt00(x_prev_[i2], x_prev_[i2+1], x_prev_[i2+2]);
////          Vec3d qt11(x_prev_[i3], x_prev_[i3+1], x_prev_[i3+2]);
////
////          // get edge vectors
////          Vec3d eq0 = qt1-qt0;
////          Vec3d eq1 = qt00-qt0;
////          Vec3d eq2 = qt11-qt0;
////
////          // determine target angle
////          Vec3d nq0 = eq0.cross(eq1);
////          nq0 /= nq0.norm();
////          Vec3d nq1 = eq0.cross(eq2);
////          nq1 /= nq1.norm();
////          double dq = nq0.dot(nq1);
////          dq = std::max(-1.0,dq);
////          dq = std::min( 1.0,dq);
////          double cq = nq0.cross(nq1).dot(eq0)/eq0.norm();
////          cp = std::max(-1.0,cp);
////          cp = std::min( 1.0,cp);
////
////          c[2] = dq;
////          c[3] = cq;
////
//////          std::cerr << "param_angle " << alpha_q/M_PI*180.0 << "vs. estimated angle " << alpha/M_PI*180.0 << std::endl;
////        }
////        else
//        {
//          c[0] = _w_dihedral;
//          c[1] = _s_dihedral;
//
//          c[0] = 1.0;
//          c[1] = 1.0;
//        }

          // add element
          _fe_dihedral.instances().add_element(idx, c);
        }

//      std::cerr << "target angle = " << target_angle/M_PI*180. << ", angle = " << angle/M_PI*180.0 << std::endl;
    }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_frame_exp_objective_function(COMISO::FiniteElementSet<FrameFittingExpElement3D_TAD> &_fe_frame_fitting,
                                   const double _sizing_factor, const double _w_exp, const double _mesh_volume,
                                   const double _limit_size, std::vector<double> &_optimal_sizing)
{
  // clear old data
  _fe_frame_fitting.instances().clear_elements();

  // allocate memory for sizing estimates
  _optimal_sizing.resize(3 * mesh_.n_cells());

  FrameFittingExpElement3D_TAD::VecI idx;
  FrameFittingExpElement3D_TAD::VecC c;
  Eigen::Matrix<double, 3, 4> P, U;

  int sizing_var_base = 3 * (max_global_idx_ + 1);

  for (auto ch: mesh_.cells())
  {
    // volume of cell
    double V = get_cell_volume(mesh_, ch);

    // collect points
    mesh_tet(ch, P);
    parametric_tet(ch, U);
    // set vertex indices
    for (unsigned int j = 0; j < 4; ++j)
      for (unsigned int i = 0; i < 3; ++i)
        idx[3 * j + i] = 3 * global_idx_[ch][j] + i;

    // set cell sizing indices
    int sidx_base = sizing_var_base + 3 * ch.idx();
    idx[12] = sidx_base;
    idx[13] = sidx_base + 1;
    idx[14] = sidx_base + 2;

    // target frame
    Mat3x3 F = _sizing_factor * frame_cprop_[ch];

//    double w = V / _mesh_volume;  // geometric weight
    double w = 1.0 / double(mesh_.n_cells()); // uniform cell weight

    // construct constants
    FrameFittingExpElement3D_TAD::compute_constants(w, _w_exp, 1.0, P, F, c);

    // add element
    _fe_frame_fitting.instances().add_element(idx, c);

    // determine initialization of sizing
    Eigen::Map<Vec3d> optimal_sizing(&(_optimal_sizing[3 * ch.idx()]));
    FrameFittingExpElement3D_TAD::estimate_sizing(P, U, F, optimal_sizing);

//    std::cerr << "cell " << ch.idx() << " optimal sizing = " << optimal_sizing.transpose() << std::endl;

    // truncate if necessary
    double acceptable_size = _limit_size + 1.0 / 100.0;
    if (optimal_sizing[0] < acceptable_size) optimal_sizing[0] = acceptable_size;
    if (optimal_sizing[1] < acceptable_size) optimal_sizing[1] = acceptable_size;
    if (optimal_sizing[2] < acceptable_size) optimal_sizing[2] = acceptable_size;

    // debug output
    if (0)
    {
      FrameFittingExpElement3D_TAD ffexp;
      FrameFittingExpElement3D_TAD::VecV x;

      Eigen::Map<Eigen::Matrix<double, 12, 1> > PP(P.data());

      x.segment<12>(0) = PP;
      x.segment<3>(12) = optimal_sizing;

      double f = ffexp.eval_f(x, c);
      std::cerr << "eval f = " << f << std::endl;
    }
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_size_fitting_objective_function(COMISO::FiniteElementSet<LinearLeastSquaresElement1D> &_fe_size_fitting,
                                      const double _weight, const double _mesh_volume)
{
  // clear old data
  _fe_size_fitting.instances().clear_elements();

  LinearLeastSquaresElement1D::VecI idx;
  LinearLeastSquaresElement1D::VecC c;

  int sizing_var_base = 3 * (max_global_idx_ + 1);

  for (auto ch: mesh_.cells())
    for (int j = 0; j < 3; ++j)
    {
      idx[0] = sizing_var_base + 3 * ch.idx() + j;
      c[0] = 1.0; // coefficient of x_i
//      c[1] = -1.0; // negative target size
      c[1] = -2.0; // negative target size
      c[2] = _weight * get_cell_volume(mesh_, ch) / _mesh_volume;
      _fe_size_fitting.instances().add_element(idx, c);
    }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_size_barrier_objective_function(COMISO::FiniteElementSet<ReciprocalBarrierElement_TAD> &_fe_size_barrier,
                                      const double _weight, const double _limit_size, const double _mesh_volume)
{
  // clear old data
  _fe_size_barrier.instances().clear_elements();

  ReciprocalBarrierElement_TAD::VecI idx;
  ReciprocalBarrierElement_TAD::VecC c;
  c[0] = _limit_size;

  int sizing_var_base = 3 * (max_global_idx_ + 1);

  for (auto ch: mesh_.cells())
    for (int j = 0; j < 3; ++j)
    {
      idx[0] = sizing_var_base + 3 * ch.idx() + j;

      c[1] = _weight * get_cell_volume(mesh_, ch) / _mesh_volume;

      _fe_size_barrier.instances().add_element(idx, c);
    }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
setup_constraints(std::vector<COMISO::LinearConstraint> &_constraints,
                  std::vector<COMISO::NConstraintInterface *> &_constraint_pointers,
                  std::vector<PairUiV> &_discrete_constraints,
                  const int _helper_idx_offset,
                  const bool _allocate_helper_variables)
{
  // reset and clear old data
  n_helper_variables_ = 0;
  _constraints.clear();
  _constraint_pointers.clear();
  _discrete_constraints.clear();

  // each manifold piece of the cut surface gets three translational variables
  if (_allocate_helper_variables)
    n_helper_variables_ = 3 * cut_surface_.size();
  else
    n_helper_variables_ = 0;

  // temp variables
  int helper_var_base = 3 * (max_global_idx_ + 1) + _helper_idx_offset;
  int n_complete = helper_var_base + n_helper_variables_;
  COMISO::LinearConstraint::SVectorNC coeffs(n_complete);
  COMISO::LinearConstraint::SVectorNC coeffs2(n_complete);
  COMISO::LinearConstraint::SVectorNC coeffs3(n_complete);

  // 1. transition functions
  for (unsigned int i = 0; i < cut_surface_.size(); ++i)
  {
    for (unsigned int j = 0; j < cut_surface_[i].size(); ++j)
    {
      // get face
      FH fh = cut_surface_[i][j];
      // get both halffaces
      HFH hfh0 = mesh_.halfface_handle(fh, 0);
      HFH hfh1 = mesh_.halfface_handle(fh, 1);

      if (mesh_.is_boundary(fh)) continue;

      // get cells
      CH ch0 = mesh_.incident_cell(hfh0);
      CH ch1 = mesh_.incident_cell(hfh1);

      // get transition from ch0 to ch1 (transition for parametrization is inverse to that of the frames)
      Mat3x3 T01 = tq_.transition_matrix_int(tq_.inverse_transition_idx(transition_hfprop_[hfh0]));

      // get vertices
      VH vh0 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[0]).to_vertex();
      VH vh1 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[1]).to_vertex();
      VH vh2 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[2]).to_vertex();

      // get vertex indices
      int idx00 = 3 * global_idx_[ch0][tet_vidx(ch0, vh0)];
      int idx01 = 3 * global_idx_[ch0][tet_vidx(ch0, vh1)];
      int idx02 = 3 * global_idx_[ch0][tet_vidx(ch0, vh2)];

      int idx10 = 3 * global_idx_[ch1][tet_vidx(ch1, vh0)];
      int idx11 = 3 * global_idx_[ch1][tet_vidx(ch1, vh1)];
      int idx12 = 3 * global_idx_[ch1][tet_vidx(ch1, vh2)];

      if (j != 0 || !_allocate_helper_variables) // j!=0 -> constraint without integer translation
      {
        for (unsigned int k = 0; k < 3; ++k)
        {
          coeffs.setZero();
          coeffs2.setZero();

          for (unsigned int m = 0; m < 3; ++m)
          {
            coeffs.coeffRef(idx01 + m) += T01(k, m);
            coeffs.coeffRef(idx00 + m) += -T01(k, m);

            coeffs2.coeffRef(idx02 + m) += T01(k, m);
            coeffs2.coeffRef(idx00 + m) += -T01(k, m);
          }
          coeffs.coeffRef(idx11 + k) += -1.0;
          coeffs.coeffRef(idx10 + k) += 1.0;

          coeffs2.coeffRef(idx12 + k) += -1.0;
          coeffs2.coeffRef(idx10 + k) += 1.0;

          _constraints.push_back(
                  COMISO::LinearConstraint(coeffs, 0.0, COMISO::NConstraintInterface::NC_EQUAL));
          _constraints.push_back(
                  COMISO::LinearConstraint(coeffs2, 0.0, COMISO::NConstraintInterface::NC_EQUAL));
        }
      }
      else // j=0 -> constraint with integer translation
      {
        for (unsigned int k = 0; k < 3; ++k)
        {
          coeffs.setZero();
          coeffs2.setZero();
          coeffs3.setZero();

          // rotate point0
          for (unsigned int m = 0; m < 3; ++m)
          {
            coeffs.coeffRef(idx00 + m) += T01(k, m);
            coeffs2.coeffRef(idx01 + m) += T01(k, m);
            coeffs3.coeffRef(idx02 + m) += T01(k, m);
          }

          // translate point0
          coeffs.coeffRef(helper_var_base + 3 * i + k) += 1.0;
          coeffs2.coeffRef(helper_var_base + 3 * i + k) += 1.0;
          coeffs3.coeffRef(helper_var_base + 3 * i + k) += 1.0;

          // subtract point1
          coeffs.coeffRef(idx10 + k) += -1.0;
          coeffs2.coeffRef(idx11 + k) += -1.0;
          coeffs3.coeffRef(idx12 + k) += -1.0;

          _constraints.push_back(
                  COMISO::LinearConstraint(coeffs, 0.0, COMISO::NConstraintInterface::NC_EQUAL));
          _constraints.push_back(
                  COMISO::LinearConstraint(coeffs2, 0.0, COMISO::NConstraintInterface::NC_EQUAL));
          _constraints.push_back(
                  COMISO::LinearConstraint(coeffs3, 0.0, COMISO::NConstraintInterface::NC_EQUAL));

          // mark translation as integer
          _discrete_constraints.push_back(PairUiV(helper_var_base + 3 * i + k, COMISO::Integer));
        }
      }
    }
  }

  // 2. feature surface alignment
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (feature_fprop_[*f_it] > 0)
    {
      // get non-boundary halfface
      HFH hfh = mesh_.halfface_handle(*f_it, 0);
      if (mesh_.is_boundary(hfh))
        hfh = mesh_.opposite_halfface_handle(hfh);

      // get corresponding cell
      CH ch = mesh_.incident_cell(hfh);

      // get vertices
      VH vh0 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[0]).to_vertex();
      VH vh1 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[1]).to_vertex();
      VH vh2 = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[2]).to_vertex();

      // get vertex indices
      int idx0 = 3 * global_idx_[ch][tet_vidx(ch, vh0)];
      int idx1 = 3 * global_idx_[ch][tet_vidx(ch, vh1)];
      int idx2 = 3 * global_idx_[ch][tet_vidx(ch, vh2)];

      // get points
      Point p0 = mesh_.vertex(vh0);
      Point p1 = mesh_.vertex(vh1);
      Point p2 = mesh_.vertex(vh2);
      // get outward normal vector
      Vec3 n = ovm2eigen((p1 - p0) % (p2 - p0));
      Vec3 e0 = ovm2eigen(p1 - p0);
      Vec3 e1 = ovm2eigen(p2 - p0);

      // get frame of cell
//                Mat3x3 F = quaternion_cprop_[ch].toRotationMatrix();
      Mat3x3 F = frame_cprop_[ch];

      // determine alignment axis
      Mat3x3 J = F.inverse();
      Vec3 Je0 = J * e0;
      Vec3 Je1 = J * e1;
      Vec3 Jn = Je0.cross(Je1);
      AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Jn);
      int axis_idx = int(e_axis) / 2;
      // test whether field is sufficiently aligned
      Vec3 a = AxisAlignmentHelpers::vector(e_axis);
      if (a.dot(Jn) < Jn.norm() * 0.99985)
        std::cerr << "Warning: face" << *f_it << " alignment axis deviates by "
                  << std::acos(a.dot(Jn) / Jn.norm()) * 180.0 / M_PI << " degree\n";

      // make all equal in normal dimensions
      coeffs.setZero();
      coeffs.coeffRef(idx1 + axis_idx) = 1.0;
      coeffs.coeffRef(idx0 + axis_idx) = -1.0;
      _constraints.push_back(COMISO::LinearConstraint(coeffs, 0.0, COMISO::NConstraintInterface::NC_EQUAL));

      coeffs.setZero();
      coeffs.coeffRef(idx2 + axis_idx) = 1.0;
      coeffs.coeffRef(idx1 + axis_idx) = -1.0;
      _constraints.push_back(COMISO::LinearConstraint(coeffs, 0.0, COMISO::NConstraintInterface::NC_EQUAL));

      // mark as integer
      _discrete_constraints.push_back(PairUiV(idx0 + axis_idx, COMISO::Integer));
      _discrete_constraints.push_back(PairUiV(idx1 + axis_idx, COMISO::Integer));
      _discrete_constraints.push_back(PairUiV(idx2 + axis_idx, COMISO::Integer));
    }

  // 3. singularity alignment
  for (EIt e_it = mesh_.e_iter(); e_it.valid(); ++e_it)
  {
    // singular edge?
    if (valence_eprop_[*e_it] != 0 || feature_eprop_[*e_it] > 0)
    {
      // get a cell adjacent to edge
      CH ch = *mesh_.hec_iter(mesh_.halfedge_handle(*e_it, 0));

      // get adjacent vertices
      VH vh0 = mesh_.edge(*e_it).from_vertex();
      VH vh1 = mesh_.edge(*e_it).to_vertex();

      // get vertex indices
      int idx0 = 3 * global_idx_[ch][tet_vidx(ch, vh0)];
      int idx1 = 3 * global_idx_[ch][tet_vidx(ch, vh1)];

      // get edge vector
      Vec3 e = ovm2eigen(mesh_.vertex(vh1) - mesh_.vertex(vh0));

      // get frame of cell
//                Mat3x3 F = quaternion_cprop_[ch].toRotationMatrix();
      Mat3x3 F = frame_cprop_[ch];

      // determine alignment axis
      Mat3x3 J = F.inverse();
      Vec3 Je = J * e;
      AxisAlignment e_axis = AxisAlignmentHelpers::get_dominant_axis(Je);
      int axis_idx = int(e_axis) / 2;
      // test whether field is sufficiently aligned
      Vec3 a = AxisAlignmentHelpers::vector(e_axis);
      if (a.dot(Je) < Je.norm() * 0.99985)
        std::cerr << "Warning: edge alignment axis deviates by " << std::acos(a.dot(Je) / Je.norm()) * 180.0 / M_PI
                  << " degree\n";

      for (int i = 0; i < 3; ++i)
        if (i != axis_idx)
        {
          // make equal on two dimensions
          coeffs.setZero();
          coeffs.coeffRef(idx1 + i) = 1.0;
          coeffs.coeffRef(idx0 + i) = -1.0;
          _constraints.push_back(
                  COMISO::LinearConstraint(coeffs, 0.0, COMISO::NConstraintInterface::NC_EQUAL));

          // mark as integer
          _discrete_constraints.push_back(PairUiV(idx0 + i, COMISO::Integer));
          _discrete_constraints.push_back(PairUiV(idx1 + i, COMISO::Integer));
        }
    }
  }

  // 4. fix origin with least possible distortion
  VH vh_origin = find_singular_or_feature_node();
  if (vh_origin.is_valid())
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "fix origin on singular/feature node" << std::endl;)
  }
  else
  {
    vh_origin = find_vertex_on_singular_or_feature_curve();
    if (vh_origin.is_valid())
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << "fix origin on singular/feature edge vertex" << std::endl;)
    }
    else
    {
      vh_origin = find_vertex_on_feature_surface();
      if (vh_origin.is_valid())
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "fix origin on feature surface vertex" << std::endl;)
      }
      else
      {
        ALGOHEX_DEBUG_ONLY(std::cerr << "fix origin on vertex with id 0" << std::endl;)
        vh_origin = VH(0);
      }
    }
  }

  if (vh_origin.is_valid())
  {
    // get incident cell
    CH ch_origin = *mesh_.vc_iter(vh_origin);
    int bidx = 3 * global_idx_[ch_origin][tet_vidx(ch_origin, vh_origin)];

    coeffs.setZero();
    coeffs.coeffRef(bidx) = 1.0;
    _constraints.push_back(COMISO::LinearConstraint(coeffs, 0.0, COMISO::NConstraintInterface::NC_EQUAL));

    coeffs.setZero();
    coeffs.coeffRef(bidx + 1) = 1.0;
    _constraints.push_back(COMISO::LinearConstraint(coeffs, 0.0, COMISO::NConstraintInterface::NC_EQUAL));

    coeffs.setZero();
    coeffs.coeffRef(bidx + 2) = 1.0;
    _constraints.push_back(COMISO::LinearConstraint(coeffs, 0.0, COMISO::NConstraintInterface::NC_EQUAL));
  }

  // add quantization constraints if available
  int n_quantization_path_constraints = 0;
  // add quantization path constraints if available
  for (auto &pc: quantization_path_constraints_)
  {
    // express coordinates of end-point in the coordinate system of the start point
    COMISO::LinearConstraint::SVectorNC cx(n_complete);
    COMISO::LinearConstraint::SVectorNC cy(n_complete);
    COMISO::LinearConstraint::SVectorNC cz(n_complete);

    bool has_transitions = false;

    // initialize endpoint
    int tvidx = tet_vidx(pc.tetTo, pc.vTo);
    if(tvidx == -1)
    {
      std::cerr << "ERROR: path constraint has vTo not incident to tetTo" << std::endl;
      continue; // skip this constraint
    }

    int bidx_to = 3 * global_idx_[pc.tetTo][tvidx];
    cx.coeffRef(bidx_to) = 1.0;
    cy.coeffRef(bidx_to + 1) = 1.0;
    cz.coeffRef(bidx_to + 2) = 1.0;

    // transform backwards along path
    for (int i = pc.pathHalffaces.size() - 1; i >= 0; --i)
    {
      // get next backward transformation (transition for parametrization is inverse to that of the frames!!!)
      HFH h = mesh_.opposite_halfface_handle(pc.pathHalffaces[i]);

      // get transition function (param transforms inversely to frame-field)
      int tf = tq_.inverse_transition_idx(transition_hfprop_[h]);

      // cut face? ---> apply transition
      if (is_on_cut_fprop_[mesh_.face_handle(h)])
      {
        // get indices of vertex on cut and corresponding cells to compute translation of transition function
        VH vh = mesh_.halfedge(mesh_.halfface(h).halfedges()[0]).to_vertex();
        HFH h_opp = mesh_.opposite_halfface_handle(h);
        CH ch = mesh_.incident_cell(h);
        CH ch_opp = mesh_.incident_cell(h_opp);

        int bidx_from = 3 * global_idx_[ch][tet_vidx(ch, vh)];
        int bidx_to = 3 * global_idx_[ch_opp][tet_vidx(ch_opp, vh)];

        transform_constraint(cx, cy, cz, tf, bidx_from, bidx_to);
      }
    }

    // subtract starting point
    tvidx = tet_vidx(pc.tetFrom, pc.vFrom);
    if(tvidx == -1)
    {
      std::cerr << "ERROR: path constraint has vFrom not incident to tetFrom" << std::endl;
      continue; // skip this constraint
    }

    int bidx_from = 3 * global_idx_[pc.tetFrom][tvidx];
    cx.coeffRef(bidx_from) -= 1.0;
    cy.coeffRef(bidx_from + 1) -= 1.0;
    cz.coeffRef(bidx_from + 2) -= 1.0;

    if (pc.offset[0] != INT_MAX)
    {
      _constraints.push_back(COMISO::LinearConstraint(cx, -pc.offset[0], COMISO::NConstraintInterface::NC_EQUAL));
      ++n_quantization_path_constraints;
    }
    if (pc.offset[1] != INT_MAX)
    {
      _constraints.push_back(COMISO::LinearConstraint(cy, -pc.offset[1], COMISO::NConstraintInterface::NC_EQUAL));
      ++n_quantization_path_constraints;
    }
    if (pc.offset[2] != INT_MAX)
    {
      _constraints.push_back(COMISO::LinearConstraint(cz, -pc.offset[2], COMISO::NConstraintInterface::NC_EQUAL));
      ++n_quantization_path_constraints;
    }
  }
  std::cerr << "#quantization path constraints = " << n_quantization_path_constraints << std::endl;

  // store pointers to constraints
  _constraint_pointers.reserve(_constraints.size());
  for (unsigned int i = 0; i < _constraints.size(); ++i)
    _constraint_pointers.push_back(&(_constraints[i]));

  ALGOHEX_DEBUG_ONLY(
          std::cerr << "#variables            : " << n_complete << " (#helper: " << n_helper_variables_ << ")"
                    << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "#constraints generated: " << _constraints.size() << std::endl;)
  ALGOHEX_DEBUG_ONLY(std::cerr << "#integer variables    : " << _discrete_constraints.size() << std::endl;)
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
transform_constraint(COMISO::LinearConstraint::SVectorNC &_cx,
                     COMISO::LinearConstraint::SVectorNC &_cy,
                     COMISO::LinearConstraint::SVectorNC &_cz,
                     const int _transition_function,
                     const int _translation_idx) const
{
  // permute
  permute_constraint(_cx, _cy, _cz, _transition_function);

  // add translation
  _cx.coeffRef(_translation_idx) += 1.0;
  _cy.coeffRef(_translation_idx + 1) += 1.0;
  _cz.coeffRef(_translation_idx + 2) += 1.0;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
transform_constraint(COMISO::LinearConstraint::SVectorNC &_cx,
                     COMISO::LinearConstraint::SVectorNC &_cy,
                     COMISO::LinearConstraint::SVectorNC &_cz,
                     const int _transition_function,
                     const int _bidx_from,
                     const int _bidx_to) const
{
  // R*c+t = R*c + p_to-R*p_from = R(c-p_from) + p_to

  // subtract p_from
  _cx.coeffRef(_bidx_from) -= 1.0;
  _cy.coeffRef(_bidx_from + 1) -= 1.0;
  _cz.coeffRef(_bidx_from + 2) -= 1.0;

  // permute
  permute_constraint(_cx, _cy, _cz, _transition_function);

  // add p_to
  _cx.coeffRef(_bidx_to) += 1.0;
  _cy.coeffRef(_bidx_to + 1) += 1.0;
  _cz.coeffRef(_bidx_to + 2) += 1.0;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
inverse_transform_constraint(COMISO::LinearConstraint::SVectorNC &_cx,
                             COMISO::LinearConstraint::SVectorNC &_cy,
                             COMISO::LinearConstraint::SVectorNC &_cz,
                             const int _transition_function,
                             const int _translation_idx) const
{
  // subtract translation
  _cx.coeffRef(_translation_idx) -= 1.0;
  _cy.coeffRef(_translation_idx + 1) -= 1.0;
  _cz.coeffRef(_translation_idx + 2) -= 1.0;


  // permute
  permute_constraint(_cx, _cy, _cz, tq_.inverse_transition_idx(_transition_function));
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
permute_constraint(COMISO::LinearConstraint::SVectorNC &_cx,
                   COMISO::LinearConstraint::SVectorNC &_cy,
                   COMISO::LinearConstraint::SVectorNC &_cz,
                   const int _transition_function) const
{
  // permute
  Mat3x3 T01 = tq_.transition_matrix_int(_transition_function);

  COMISO::LinearConstraint::SVectorNC cxn(_cx.size());
  COMISO::LinearConstraint::SVectorNC cyn(_cy.size());
  COMISO::LinearConstraint::SVectorNC czn(_cz.size());

  // get new cx
  if (T01(0, 0) != 0.0)
    cxn = T01(0, 0) * _cx;
  if (T01(0, 1) != 0.0)
    cxn = T01(0, 1) * _cy;
  if (T01(0, 2) != 0.0)
    cxn = T01(0, 2) * _cz;

  // get new cy
  if (T01(1, 0) != 0.0)
    cyn = T01(1, 0) * _cx;
  if (T01(1, 1) != 0.0)
    cyn = T01(1, 1) * _cy;
  if (T01(1, 2) != 0.0)
    cyn = T01(1, 2) * _cz;

  // get new cz
  if (T01(2, 0) != 0.0)
    czn = T01(2, 0) * _cx;
  if (T01(2, 1) != 0.0)
    czn = T01(2, 1) * _cy;
  if (T01(2, 2) != 0.0)
    czn = T01(2, 2) * _cz;

  // swap new ones with current ones
  _cx.swap(cxn);
  _cy.swap(cyn);
  _cz.swap(czn);
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
add_feature_constraints(std::vector<COMISO::LinearConstraint> &_constraints,
                        std::vector<COMISO::NConstraintInterface *> &_constraint_pointers)
{
  int n_image = 3 * (max_global_idx_ + 1);
  int n_domain = 3 * mesh_.n_vertices();
  int n_helper = n_helper_variables_;
  int n_complete = n_image + n_domain + n_helper;

  COMISO::LinearConstraint::SVectorNC coeffs(n_complete);

  // local feature halfedges
  std::vector<HEH> fheh;
  // local feature halffaces
  std::vector<FH> ffh;

  int n_feature_vertices = 0;
  int n_on_feature_arcs = 0;
  int n_on_feature_patches = 0;
  for (auto vh: mesh_.vertices())
  {
    // get incident feature halfedges
    int n_ifh = incident_feature_ohalfedges(vh, fheh);

    // get current position
    Vec3d p = ovm2eigen(mesh_.vertex(vh));

    if (feature_vprop_[vh] > 0 || n_ifh > 2)
    {
      ++n_feature_vertices;
      int i = n_image + 3 * vh.idx();

      for (int j = 0; j < 3; ++j)
      {
        coeffs.setZero();
        coeffs.coeffRef(i + j) = 1.0;
        _constraints.push_back(COMISO::LinearConstraint(coeffs, -p[j], COMISO::NConstraintInterface::NC_EQUAL));
      }
    }
    else if (n_ifh == 2)
    {
      ++n_on_feature_arcs;
      Vec3d t = ovm2eigen(mesh_.vector(fheh[0]) - mesh_.vector(fheh[1]));
      t /= t.norm();

      if (std::isfinite(t[0] + t[1] + t[2]))
      {
        Vec3d u, v, w;
        complement_to_right_handed_orthonormal_frame(t, u, v, w);

        int i = n_image + 3 * vh.idx();
        coeffs.setZero();
        coeffs.coeffRef(i) = v[0];
        coeffs.coeffRef(i + 1) = v[1];
        coeffs.coeffRef(i + 2) = v[2];
        _constraints.push_back(COMISO::LinearConstraint(coeffs, -p.dot(v), COMISO::NConstraintInterface::NC_EQUAL));

        coeffs.setZero();
        coeffs.coeffRef(i) = w[0];
        coeffs.coeffRef(i + 1) = w[1];
        coeffs.coeffRef(i + 2) = w[2];
        _constraints.push_back(COMISO::LinearConstraint(coeffs, -p.dot(w), COMISO::NConstraintInterface::NC_EQUAL));
      }
      else
        std::cerr << "ERROR: tangent vector is undefined t=" << mesh_.vector(fheh[0]) - mesh_.vector(fheh[1])
                  << std::endl;
    }
    else if (incident_feature_faces(vh, ffh))
    {
      ++n_on_feature_patches;

      // estimate feature patch normal
      Vec3d n(0, 0, 0);
      for (auto fh: ffh)
      {
        Vec3d n_cur = ovm2eigen(mesh_.normal(mesh_.halfface_handle(fh, 0)));
        if (n.dot(n_cur) < 0.0)
          n -= n_cur;
        else
          n += n_cur;
      }

      Vec3d nn = n / n.norm();

      if (std::isfinite(nn[0] + nn[1] + nn[2]))
      {
        int i = n_image + 3 * vh.idx();
        coeffs.setZero();
        coeffs.coeffRef(i) = nn[0];
        coeffs.coeffRef(i + 1) = nn[1];
        coeffs.coeffRef(i + 2) = nn[2];
        _constraints.push_back(COMISO::LinearConstraint(coeffs, -p.dot(nn), COMISO::NConstraintInterface::NC_EQUAL));
      }
      else std::cerr << "ERROR: normal vector is undefined t=" << nn.transpose() << std::endl;
    }
  }

  std::cerr << "#FeatureVertices = " << n_feature_vertices << std::endl;
  std::cerr << "#OnFeatureArc    = " << n_on_feature_arcs << std::endl;
  std::cerr << "#OnFeaturePatch  = " << n_on_feature_patches << std::endl;

  // rebuild constraint pointers
  _constraint_pointers.clear();
  _constraint_pointers.reserve(_constraints.size());
  for (unsigned int i = 0; i < _constraints.size(); ++i)
    _constraint_pointers.push_back(&(_constraints[i]));
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
typename Parametrization3DT<TetMeshT>::VH
Parametrization3DT<TetMeshT>::
find_singular_or_feature_node() const
{
  for (VIt v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
  {
    int n = incident_singular_edges(*v_it);

    if (n == 1 || n >= 3)
      return *v_it;

    if (feature_vprop_[*v_it] > 0)
      return *v_it;
  }

  // none found
  return VH(-1);
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
typename Parametrization3DT<TetMeshT>::VH
Parametrization3DT<TetMeshT>::
find_vertex_on_singular_or_feature_curve() const
{
  for (VIt v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
  {
    int n = incident_singular_edges(*v_it);

    if (n == 2)
      return *v_it;

    int m = incident_feature_edges(*v_it);

    if (m == 2)
      return *v_it;
  }

  // none found
  return VH(-1);
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
typename Parametrization3DT<TetMeshT>::VH
Parametrization3DT<TetMeshT>::
find_vertex_on_feature_surface() const
{
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (feature_fprop_[*f_it] > 0)
    {
      HFH hfh = mesh_.halfface_handle(*f_it, 0);

      VH vh = mesh_.halfedge(mesh_.halfface(hfh).halfedges()[0]).to_vertex();

      return vh;
    }

  // none found
  return VH(-1);
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
incident_singular_edges(const VH _vh) const
{
  int n(0);
  for (VOHEIt vohe_it = mesh_.voh_iter(_vh); vohe_it.valid(); ++vohe_it)
  {
    // count number of singular edges
    if (valence_eprop_[mesh_.edge_handle(*vohe_it)] != 0)
      ++n;
  }
  return n;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
incident_feature_edges(const VH _vh) const
{
  int n(0);
  for (VOHEIt vohe_it = mesh_.voh_iter(_vh); vohe_it.valid(); ++vohe_it)
  {
    // count number of singular edges
    if (feature_eprop_[mesh_.edge_handle(*vohe_it)] > 0)
      ++n;
  }
  return n;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
tet_vidx(const CH _ch, const VH _vh) const
{
  int i = 0;
  for (TVIt tv_it = mesh_.tv_iter(_ch); tv_it.valid(); ++tv_it)
  {
    if (_vh == *tv_it)
      return i;
    ++i;
  }
  // not found?
  return -1;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
halfface_vidx(const HFH _hfh, const VH _vh) const
{
  for (unsigned int i = 0; i < mesh_.halfface(_hfh).halfedges().size(); ++i)
  {
    VH vh = mesh_.halfedge(mesh_.halfface(_hfh).halfedges()[i]).to_vertex();
    if (vh == _vh)
      return int(i);
  }

  // not found?
  return -1;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
typename Parametrization3DT<TetMeshT>::HFH
Parametrization3DT<TetMeshT>::
opposite_halfface_handle(const CH _ch, const VH _vh) const
{
  // check for valid input
  assert(tet_vidx(_ch, _vh) != -1);

  for (unsigned int i = 0; i < mesh_.cell(_ch).halffaces().size(); ++i)
  {
    HFH hfh = mesh_.cell(_ch).halffaces()[i];
    if (halfface_vidx(hfh, _vh) == -1)
      return hfh;
  }

  return HFH(-1);
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
n_cut_faces(const EH _eh) const
{
  std::vector<FH> fhs;
  return n_cut_faces(_eh, fhs);
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
n_cut_faces(const EH _eh, std::vector<FH> &_fhs) const
{
  int n(0);
  _fhs.clear();

  HEH heh = mesh_.halfedge_handle(_eh, 0);
  for (HEHFIt hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
  {
    FH fh = mesh_.face_handle(*hehf_it);
    if (is_on_cut_fprop_[fh])
    {
      // store reference
      _fhs.push_back(fh);
      // increment counter
      ++n;
    }
  }

  return n;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
double
Parametrization3DT<TetMeshT>::
octahedral_field_energy() const
{
  double energy(0.0);

  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (!mesh_.is_boundary(*f_it))
    {
      HFH hfh0 = mesh_.halfface_handle(*f_it, 0);
      HFH hfh1 = mesh_.halfface_handle(*f_it, 1);

      CH ch0 = mesh_.incident_cell(hfh0);
      CH ch1 = mesh_.incident_cell(hfh1);

      Mat3x3 trans = tq_.transition_matrix_int(transition_hfprop_[hfh0]);

//              Eigen::Quaterniond tq = tq_.transition(transition_hfprop_[hfh0]);

//                Mat3x3 F0 = (quaternion_cprop_[ch0] * tq).toRotationMatrix();
//                Mat3x3 F1 = quaternion_cprop_[ch1].toRotationMatrix();

      Mat3x3 F0 = frame_cprop_[ch0] * trans;
      Mat3x3 F1 = frame_cprop_[ch1];

      energy += (F0 - F1).squaredNorm();
    }

  return energy;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
test_objective_function() const
{
  COMISO::FiniteElementSet<FrameFittingElement3D>::VecC c;
  COMISO::FiniteElementSet<FrameFittingElement3D>::VecV x;

  FrameFittingElement3D ffe;

  Mat3x4d P;
  P.col(0) = Vec3(-0.1, 0, 0.1);
  P.col(1) = Vec3(1, 0.1, 0.2);
  P.col(2) = Vec3(0, 1.1, 0.2);
  P.col(3) = Vec3(0, 0.1, 1.2);

  Eigen::Matrix3d F;
  F.col(0) = Vec3(1, 0, 0);
  F.col(1) = Vec3(0, 1, 0);
  F.col(2) = Vec3(0, 0, 1);

  // construct constants
  FrameFittingElement3D::constants_from_tetrahedron_and_frame(P, F, 1.0, c);

  for (unsigned int i = 0; i < 3; ++i)
  {
    x[4 * i + 0] = P.col(0)[i];
    x[4 * i + 1] = P.col(1)[i];
    x[4 * i + 2] = P.col(2)[i];
    x[4 * i + 3] = P.col(3)[i];
  }

  double d = ffe.eval_f(x, c);
  std::cerr << "identity transformation should lead to zero energy: " << d << std::endl;

  Vec3 offset(-1, 3, 5);

  for (unsigned int i = 0; i < 3; ++i)
  {
    x[4 * i + 0] = P.col(0)[i] + offset[i];
    x[4 * i + 1] = P.col(1)[i] + offset[i];
    x[4 * i + 2] = P.col(2)[i] + offset[i];
    x[4 * i + 3] = P.col(3)[i] + offset[i];
  }

  d = ffe.eval_f(x, c);
  std::cerr << "translated identity transformation should lead to zero energy: " << d << std::endl;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
double
Parametrization3DT<TetMeshT>::
degenerate_volume_map()
{
  const double eps = 10e-10;

  double V(0.0);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    Mat3x4d P;
    parametric_tet(*c_it, P);
    Mat3x3 E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    if (E.determinant() <= eps)
    {
      mesh_tet(*c_it, P);
      E.col(0) = P.col(1) - P.col(0);
      E.col(1) = P.col(2) - P.col(0);
      E.col(2) = P.col(3) - P.col(0);

      V += 1.0 / 6.0 * E.determinant();
    }
  }

  return V;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
double
Parametrization3DT<TetMeshT>::
mesh_volume()
{
  double V(0.0);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    Mat3x4d P;
    mesh_tet(*c_it, P);
    Mat3x3 E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    V += 1.0 / 6.0 * E.determinant();
  }

  return V;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
double
Parametrization3DT<TetMeshT>::
parametric_volume()
{
  double V(0.0);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    Mat3x4d P;
    parametric_tet(*c_it, P);
    Mat3x3 E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    V += 1.0 / 6.0 * E.determinant();
  }

  return V;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
double
Parametrization3DT<TetMeshT>::
frame_field_volume()
{
  double V(0.0);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    Mat3x4d P;
    mesh_tet(*c_it, P);
    Mat3x3 E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    Mat3x3 F = frame_cprop_[*c_it];

    // we want det[ 1/6(F^-1 * E)] =det[E]/(6 det[F])
    V += E.determinant() / (6.0 * F.determinant());


//    // deform tet with Jacobian induced by the frame
//    Mat3x3 J = frame_cprop_[*c_it].inverse();
//    Mat3x3 JE = J*E;
//
//    V += std::abs(1.0 / 6.0 * JE.determinant());
  }

  return V;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
number_of_invalid_parametric_tets()
{
  int n(0);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    Mat3x4d P;
    parametric_tet(*c_it, P);
    Mat3x3 E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    double d = E.determinant();

    if (d <= 0.0 || !std::isfinite(d))
      ++n;
  }

  return n;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
number_of_invalid_edge_valencies()
{
  int n(0);

  for (auto eh: mesh_.edges())
  {
    // compute expected valence
    double expected_val = 4;
    if (mesh_.is_boundary(eh))
      expected_val = 2;
    expected_val += valence_eprop_[eh];

    // compute valence
    std::vector<double> av;
    double val = 2.0 * edge_angle(eh, av) / M_PI;

    // check for inverted tets in one-ring
    bool valid_one_ring = true;
    for (auto a: av)
      if (a >= M_PI || !std::isfinite(a))
      {
        valid_one_ring = false;
        break;
      }

    // only count wrong valencies with valid one-ring
    if (valid_one_ring && std::abs(val - expected_val) > 0.1)
    {
//      std::cerr << "Warning: valence mismatch: " << val << " should be " << expected_val << " (on boundary = " << int(mesh_.is_boundary(eh)) << ")" << std::endl;
//      std::cerr << "         tet angles = ";
//      for(auto a : av) std::cerr << a/M_PI*180.0 << " ";
//      std::cerr << std::endl;
      ++n;
    }
  }

  return n;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
parametric_tet(const CH _ch, Mat3x4d &_P)
{
  int i = 0;
  for (TVIt tv_it = mesh_.tv_iter(_ch); tv_it.valid(); ++tv_it)
  {
    _P.col(i++) = ovm2eigen(Point(igm_cprop_[_ch][*tv_it]));
  }
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
parametric_tet(const CH _ch, const std::vector<double> &_x, Mat3x4d &_P)
{
  int i = 0;
  for (TVIt tv_it = mesh_.tv_iter(_ch); tv_it.valid(); ++tv_it)
  {
    int bidx = 3 * global_idx_[_ch][i];
    _P.col(i) = ovm2eigen(Point(_x[bidx], _x[bidx + 1], _x[bidx + 2]));
    ++i;
  }
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
jacobi_matrix(const Mat3x4d &_P, const Mat3x4d &_U, Mat3d &_J)
{
  // get edge vectors
  Mat3d E;
  E.col(0) = _P.col(1) - _P.col(0); //q1-q0;
  E.col(1) = _P.col(2) - _P.col(0); //q2-q0;
  E.col(2) = _P.col(3) - _P.col(0); //q3-q0;

  double dE = E.determinant();
  double dEi = 1.0 / dE;

  // compute gradients of pwl basis functions G = [g0 | g1 | g2]
  Mat3d G;
  G.col(0) = dEi * (E.col(1)).cross(E.col(2));
  G.col(1) = dEi * (E.col(2)).cross(E.col(0));
  G.col(2) = dEi * (E.col(0)).cross(E.col(1));

  // get parametric edge vectors
  Mat3d EU;
  EU.col(0) = _U.col(1) - _U.col(0);
  EU.col(1) = _U.col(2) - _U.col(0);
  EU.col(2) = _U.col(3) - _U.col(0);

  _J = EU * G.transpose();

  // debug
  if (0)
  {
    double diff = (_J * E - EU).norm();

    if (diff > 1e-4)
    {
      std::cerr << "Warning: Jacobian is not correct" << std::endl;
      std::cerr << _J * E << std::endl;
      std::cerr << " = J*E vs E_u = " << std::endl;
      std::cerr << EU << std::endl;

      std::cerr << "check 1 = g0*e0 = " << G.col(0).dot(E.col(0)) << std::endl;
      std::cerr << "check 1 = g1*e1 = " << G.col(1).dot(E.col(1)) << std::endl;
      std::cerr << "check 1 = g2*e2 = " << G.col(2).dot(E.col(2)) << std::endl;
    }
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
double
Parametrization3DT<TetMeshT>::
edge_angle(const EH &_eh, std::vector<double> &_tet_angles)
{
  _tet_angles.clear();

  double angle_sum = 0.;
  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  for (auto hehf_it = mesh_.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
  {
    auto ch = mesh_.incident_cell(*hehf_it);
    if (ch.is_valid())
    {
      auto cvhs = mesh_.get_cell_vertices(*hehf_it, heh0);
      std::vector<Point> pm_pts;
      pm_pts.reserve(4);
      for (const auto &vhi: cvhs)
        pm_pts.push_back(igm_cprop_[ch][vhi]);


      auto n0 = face_normal(pm_pts[0], pm_pts[1], pm_pts[3]);
      auto n1 = face_normal(pm_pts[0], pm_pts[2], pm_pts[1]);

      double dh = dihedral_angle(pm_pts[0], pm_pts[1], n0, n1);
      _tet_angles.push_back(dh);

      angle_sum += dh;
    }
  }

  return angle_sum;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
double
Parametrization3DT<TetMeshT>::
mesh_edge_angle(const EH &_eh, std::vector<double> &_tet_angles)
{
  _tet_angles.clear();

  double angle_sum = 0.;
  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  for (auto hehf_it = mesh_.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
  {
    auto ch = mesh_.incident_cell(*hehf_it);
    if (ch.is_valid())
    {
      auto cvhs = mesh_.get_cell_vertices(*hehf_it, heh0);
      std::vector<Point> pm_pts;
      pm_pts.reserve(4);
      for (const auto &vhi: cvhs)
        pm_pts.push_back(mesh_.vertex(vhi));

      auto n0 = face_normal(pm_pts[0], pm_pts[1], pm_pts[3]);
      auto n1 = face_normal(pm_pts[0], pm_pts[2], pm_pts[1]);

      double dh = dihedral_angle(pm_pts[0], pm_pts[1], n0, n1);
      angle_sum += dh;

      _tet_angles.push_back(dh);
    }
  }

  return angle_sum;
}



//-----------------------------------------------------------------------------


template<class TetMeshT>
bool
Parametrization3DT<TetMeshT>::
param_invalid_edge_valence_valid_tets(const EH &_eh)
{
  bool tets_valid = true;

  double angle_sum = 0.;
  auto heh0 = mesh_.halfedge_handle(_eh, 0);
  for (auto hehf_it = mesh_.hehf_iter(heh0); hehf_it.valid(); ++hehf_it)
  {
    auto ch = mesh_.incident_cell(*hehf_it);
    if (ch.is_valid())
    {
      auto cvhs = mesh_.get_cell_vertices(*hehf_it, heh0);
      std::vector<Point> pm_pts;
      pm_pts.reserve(4);
      for (const auto &vhi: cvhs)
        pm_pts.push_back(igm_cprop_[ch][vhi]);

      auto n0 = face_normal(pm_pts[0], pm_pts[1], pm_pts[3]);
      auto n1 = face_normal(pm_pts[0], pm_pts[2], pm_pts[1]);

      double dh = dihedral_angle(pm_pts[0], pm_pts[1], n0, n1);
      angle_sum += dh;

      if (dh > M_PI)
        tets_valid = false;
    }
  }

  // calct target valence
  int target_valence = 4;
  if (mesh_.is_boundary(_eh))
    target_valence = 2;
  target_valence += valence_eprop_[_eh];

  bool edge_valence_valid = (std::abs(double(target_valence) - angle_sum / (0.5 * M_PI)) < 0.01);

  // return true if edge valence is invalid and tets are valie ---> double covering case
  return (!edge_valence_valid && tets_valid);
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
typename Parametrization3DT<TetMeshT>::Point
Parametrization3DT<TetMeshT>::
face_normal(const Point &_p0, const Point &_p1, const Point &_p2) const
{
  Point p01 = _p1 - _p0;
  Point p02 = _p2 - _p0;

  Point normal = p01.cross(p02);
  normal.normalize();

  return normal;
}



//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
mesh_tet(const CH _ch, Mat3x4d &_P)
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
Parametrization3DT<TetMeshT>::
reset_cell_weights()
{
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    cell_weight_[*c_it] = 1.0;
  }
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
update_cell_weights(const std::vector<double> &_x)
{
  int n(0);

//  std::vector<double> wup (mesh_.n_cells(),0.0);
//  std::vector<double> wup2(mesh_.n_cells(),0.0);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    Mat3x4d P;
    parametric_tet(*c_it, _x, P);
    Mat3x3 E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    if (E.determinant() <= 0.0)
    {
      ++n;
      cell_weight_[*c_it] *= stiffening_weight_;
//      wup[c_it->idx()] = (stiffening_weight_-1.0)*cell_weight_[*c_it];
    }
  }

//  // diffuse
//  for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
//  {
//    double avg(0.0);
//    double n(0.0);
//    for(CCIt cc_it = mesh_.cc_iter(*c_it); cc_it.valid(); ++cc_it)
//    {
//      avg += wup[cc_it->idx()];
//      n += 1.0;
//    }
//
//    avg /= n;
//
//    wup2[c_it->idx()] += wup[c_it->idx()] + 0.5*(avg - wup[c_it->idx()]);
//  }
//  std::swap(wup,wup2);
//
//  // update
//  for(CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
//  {
//    cell_weight_[*c_it] += wup[c_it->idx()];
//  }

  return n;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
import_frames_from_quaternion_field(const std::string &_prop_name)
{
  OpenVolumeMesh::CellPropertyT<Eigen::Quaterniond> quaternion_cprop = mesh_.template request_cell_property<Eigen::Quaterniond>(
          "FrameFieldQuaternions");

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    frame_cprop_[*c_it] = quaternion_cprop[*c_it].toRotationMatrix();
  }
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
import_frames_from_parametrization()
{
  std::cerr << "import_frames_from_parametrization..." << std::endl;
  for (auto ch: mesh_.cells())
  {
    // get tet in mesh and param domain
    Mat3x4d P, U;
    mesh_tet(ch, P);
    parametric_tet(ch, U);

    Mat3d J;
    jacobi_matrix(P, U, J);
    Mat3d Ji = J.inverse();

    if (std::isfinite(Ji.sum()))
      frame_cprop_[ch] = Ji;
    else
    {
      std::cerr << "Warning: valid frame could not be extracted from parametrization.... " << std::endl << Ji
                << std::endl;
      frame_cprop_[ch] = Mat3d::Identity();
    }
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
add_boundary_faces_to_feature_surfaces()
{
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    if (mesh_.is_boundary(*f_it))
      feature_fprop_[*f_it] = 1;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
reset_feature_face_property()
{
  for (FIt f_it = mesh_.f_iter(); f_it.valid(); ++f_it)
    feature_fprop_[*f_it] = 0;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
scale_frames_to_target_complexity(const int _num_hex_cells)
{
  // compute scaling based on target complexity
  double sizing_scale_factor = std::cbrt(frame_field_volume() / _num_hex_cells);

  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    frame_cprop_[*c_it] *= sizing_scale_factor;
  }

}


//-----------------------------------------------------------------------------


template<class TetMeshT>
template<class ExportTetMeshT>
void
Parametrization3DT<TetMeshT>::
export_invalid_parametric_tets(ExportTetMeshT &_tetmesh_ref, TetMeshT &_tetmesh_param)
{
  // properties
  OpenVolumeMesh::CellPropertyT<double> cell_type_ref = _tetmesh_ref.template request_cell_property<double>("CellType");
  OpenVolumeMesh::CellPropertyT<double> cell_type_param = _tetmesh_param.template request_cell_property<double>(
          "CellType");
  _tetmesh_ref.set_persistent(cell_type_ref, true);
  _tetmesh_param.set_persistent(cell_type_param, true);

  OpenVolumeMesh::EdgePropertyT<double> edge_type_ref = _tetmesh_ref.template request_edge_property<double>("CellType");
  OpenVolumeMesh::EdgePropertyT<double> edge_type_param = _tetmesh_param.template request_edge_property<double>(
          "CellType");
  _tetmesh_ref.set_persistent(edge_type_ref, true);
  _tetmesh_param.set_persistent(edge_type_param, true);


  // export invalid tets
  std::cerr << "export invalid tets..." << std::endl;
  int n_invalid_tets = 0;
  for (CIt c_it = mesh_.c_iter(); c_it.valid(); ++c_it)
  {
    Mat3x4d P;
    parametric_tet(*c_it, P);
    Mat3x3 E;
    E.col(0) = P.col(1) - P.col(0);
    E.col(1) = P.col(2) - P.col(0);
    E.col(2) = P.col(3) - P.col(0);

    if (E.determinant() <= 0.0)
    {
      auto chp = add_tet(_tetmesh_param, P);
      cell_type_param[chp] = 1.0;

      mesh_tet(*c_it, P);
      auto chr = add_tet(_tetmesh_ref, P);
      cell_type_ref[chr] = 1.0;

      Mat3x4d Q;
      mesh_tet(*c_it, Q);
      double tq = tetrahedron_quality(Q.col(0), Q.col(1), Q.col(2), Q.col(3));

      // mark tet as inverted by adding a point in its barycenter
      auto b_ref = _tetmesh_ref.barycenter(chr);
      auto vh_ref = _tetmesh_ref.add_vertex();
      _tetmesh_ref.set_vertex(vh_ref, b_ref);

      // mark tet as inverted by adding a point in its barycenter
      auto b_param = _tetmesh_param.barycenter(chp);
      auto vh_param = _tetmesh_param.add_vertex();
      _tetmesh_param.set_vertex(vh_param, b_param);

      std::cerr << "exported invalid tet has quality " << tq << std::endl;
      ++n_invalid_tets;
    }
  }

  std::cerr << "export invalid edges..." << std::endl;
  // export invalid edges
  for (auto eh: mesh_.edges())
  {
    if (param_invalid_edge_valence_valid_tets(eh))
    {
      HEH heh = mesh_.halfedge_handle(eh, 0);
      VH vh0 = mesh_.from_vertex_handle(heh);
      VH vh1 = mesh_.to_vertex_handle(heh);
      Eigen::Matrix<double, 3, 2> P;
      P.col(0) = ovm2eigen(mesh_.vertex(vh0));
      P.col(1) = ovm2eigen(mesh_.vertex(vh1));
      auto eh2 = add_edge(_tetmesh_ref, P);

      edge_type_ref[eh2] = 1.0;

      ECIt ec_it(eh, &mesh_);
      for (; ec_it.valid(); ++ec_it)
      {
        Mat3x4d P;
        parametric_tet(*ec_it, P);
        auto chp = add_tet(_tetmesh_param, P);
        cell_type_param[chp] = -1.0;

        mesh_tet(*ec_it, P);
        auto chr = add_tet(_tetmesh_ref, P);
        cell_type_ref[chr] = -1.0;
      }
    }
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
template<class ExportTetMeshT>
CH
Parametrization3DT<TetMeshT>::
add_tet(ExportTetMeshT &_tetmesh, Mat3x4d &_P)
{
  auto vh0 = _tetmesh.add_vertex();
  auto vh1 = _tetmesh.add_vertex();
  auto vh2 = _tetmesh.add_vertex();
  auto vh3 = _tetmesh.add_vertex();

  _tetmesh.set_vertex(vh0, eigen2ovm(_P.col(0)));
  _tetmesh.set_vertex(vh1, eigen2ovm(_P.col(1)));
  _tetmesh.set_vertex(vh2, eigen2ovm(_P.col(2)));
  _tetmesh.set_vertex(vh3, eigen2ovm(_P.col(3)));

  auto ch = _tetmesh.add_cell(vh0, vh1, vh2, vh3);
  return ch;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
template<class ExportTetMeshT>
EH
Parametrization3DT<TetMeshT>::
add_edge(ExportTetMeshT &_tetmesh, Eigen::Matrix<double, 3, 2> &_P)
{
  auto vh0 = _tetmesh.add_vertex();
  auto vh1 = _tetmesh.add_vertex();

  _tetmesh.set_vertex(vh0, eigen2ovm(_P.col(0)));
  _tetmesh.set_vertex(vh1, eigen2ovm(_P.col(1)));

  auto eh = _tetmesh.add_edge(vh0, vh1);
  return eh;
}


//-----------------------------------------------------------------------------

template<class TetMeshT>
int
Parametrization3DT<TetMeshT>::
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
Parametrization3DT<TetMeshT>::
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
double
Parametrization3DT<TetMeshT>::
clamp(const double _value, const double _min_val, const double _max_val) const
{
  if (_value < _min_val)
    return _min_val;
  if (_value > _max_val)
    return _max_val;

  return _value;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
quantize_qgp3d(const int _num_hex_cells)
{
  // tetMesh: your mc3d::TetMesh
  qgp3d::Quantizer quantizer(mesh_);

  // param: your property defining mc3d::Vec3d parameter values for each tet-vertex-pair
  for (auto tet: mesh_.cells())
    for (auto v: mesh_.cell_vertices(tet))
      quantizer.setParam(tet, v, igm_cprop_[tet][v]);

  // isFeatureF, isFeatureE, isFeatureV: your properties marking features
  for (auto f: mesh_.faces())
    if (feature_fprop_[f] > 0)
      quantizer.setFeature(f, true);
  for (auto e: mesh_.edges())
    if (feature_eprop_[e] > 0)
      quantizer.setFeature(e, true);
  for (auto v: mesh_.vertices())
    if (feature_vprop_[v] > 0)
      quantizer.setFeature(v, true);

  // This scaling factor is applied to the parametrization before quantization
  double scaling = std::cbrt(parametric_volume() / _num_hex_cells);
  std::cerr << "sizing scale factor for quantization = " << scaling <<std::endl;

  // This will hold the number of hexahedra implied by the quantization
  int nHexahedra = 0;

  // Compute a quantization and retrieve the corresponding spacing constraints
  quantization_path_constraints_.clear();
  auto result = quantizer.quantize(scaling, quantization_path_constraints_, nHexahedra);
  if (result != qgp3d::Quantizer::RetCode::SUCCESS) {
      std::cerr << "ERROR: QGP3D quantization failed! Code " << result << std::endl;
      throw std::runtime_error("QGP3D failed");
  }

  std::cerr << "#hexahedra after QGP3D quantization = " << nHexahedra << std::endl;
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
initialize_helper_variables_via_seamless_map(std::vector<double> &_x, const int _helper_var_base_idx)
{
  // calculate translation part of transition functions on cut
  for (unsigned int i = 0; i < cut_surface_.size(); ++i)
  {
    // get face
    FH fh = cut_surface_[i][0];
    // get both halffaces
    HFH hfh0 = mesh_.halfface_handle(fh, 0);
    HFH hfh1 = mesh_.halfface_handle(fh, 1);

    if (mesh_.is_boundary(fh)) continue;

    // get cells
    CH ch0 = mesh_.incident_cell(hfh0);
    CH ch1 = mesh_.incident_cell(hfh1);

    // get transition from ch0 to ch1 (transition for parametrization is inverse to that of the frames)
    Mat3x3 T01 = tq_.transition_matrix_int(tq_.inverse_transition_idx(transition_hfprop_[hfh0]));

    // get vertices
    VH vh0 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[0]).to_vertex();
//    VH vh1 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[1]).to_vertex();
//    VH vh2 = mesh_.halfedge(mesh_.halfface(hfh0).halfedges()[2]).to_vertex();

    // get vertex indices
    int idx00 = 3 * global_idx_[ch0][tet_vidx(ch0, vh0)];
//    int idx01 = 3 * global_idx_[ch0][tet_vidx(ch0, vh1)];
//    int idx02 = 3 * global_idx_[ch0][tet_vidx(ch0, vh2)];

    int idx10 = 3 * global_idx_[ch1][tet_vidx(ch1, vh0)];
//    int idx11 = 3 * global_idx_[ch1][tet_vidx(ch1, vh1)];
//    int idx12 = 3 * global_idx_[ch1][tet_vidx(ch1, vh2)];

    Vec3 p0(_x[idx00], _x[idx00 + 1], _x[idx00 + 2]);
    Vec3 p1(_x[idx10], _x[idx10 + 1], _x[idx10 + 2]);

    Vec3 t = p1 - T01 * p0;

    int hvidx = _helper_var_base_idx + 3 * i;
    _x[hvidx] = t[0];
    _x[hvidx + 1] = t[1];
    _x[hvidx + 2] = t[2];
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
convertTetrahedralMeshIndexing(const TetMeshT &_mesh, TetMeshT &_tetMesh)
{
  // add vertices
  _tetMesh.clear(false);

  for (auto v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    _tetMesh.add_vertex(_mesh.vertex(*v_it));

  // add tets
  for (auto c_it = _mesh.cells_begin(); c_it != _mesh.cells_end(); ++c_it)
  {
    auto vertices = _mesh.get_cell_vertices(*c_it);
    _tetMesh.add_cell(vertices);
  }
}


//-----------------------------------------------------------------------------


template<class TetMeshT>
void
Parametrization3DT<TetMeshT>::
check_quantization_path_constraints()
{
  // debug code: compute differences in seamless map
  int j=0;
  int n_valid=0;
  int n_invalid=0;
  for (auto &pc: quantization_path_constraints_)
  {
    bool is_valid = true;
    // debug check whether dual path is valid
    for (int i = 0; i+1 < pc.pathHalffaces.size(); ++i) {
      HFH hfh0 = pc.pathHalffaces[i];
      HFH hfh1 = pc.pathHalffaces[i + 1];
      if (mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh0)) != mesh_.incident_cell(hfh1))
      {
        std::cerr << "ERROR in dual path!!! path_idx = " << j << std::endl;
        is_valid = false;
      }
    }
    // check start and end vertex
    if (!pc.pathHalffaces.empty()) {
      if (pc.tetFrom != mesh_.incident_cell(pc.pathHalffaces[0]))
      {
        std::cerr << "ERROR in dual path start !!! path_idx = " << j << std::endl;
        is_valid = false;
      }

      if (pc.tetTo != mesh_.incident_cell(
              mesh_.opposite_halfface_handle(pc.pathHalffaces[pc.pathHalffaces.size() - 1])))
      {
        std::cerr << "ERROR in dual path end !!! path_idx = " << j << std::endl;
        is_valid = false;
      }
    }

    // check start vertex
    auto vertices = mesh_.get_cell_vertices(pc.tetFrom);
    if(std::find(vertices.begin(),vertices.end(),pc.vFrom) == vertices.end())
    {
      std::cerr << "ERROR: path constraint " << j << " has a start vertex " << pc.vFrom.idx() << " not being incident to start tet " << pc.tetFrom.idx() << std::endl;
      is_valid = false;
    }

    // check end vertex
    vertices = mesh_.get_cell_vertices(pc.tetTo);
    if(std::find(vertices.begin(),vertices.end(),pc.vTo) == vertices.end())
    {
      std::cerr << "ERROR: path constraint " << j << " has an end vertex " << pc.vTo.idx() << " not being incident to end tet " << pc.tetTo.idx() << std::endl;
      is_valid = false;
    }

    if(is_valid)
      ++n_valid;
    else
      ++n_invalid;
    // move to next constraint
    ++j;
  }

  int count = 0;
  for (auto &pc: quantization_path_constraints_)
  {
    // initialize endpoint
    Vec3 p = ovm2eigen(igm_cprop_[pc.tetTo][pc.vTo]);

    // transform backwards along path
    for (int i = (int)pc.pathHalffaces.size() - 1; i >= 0; --i) {
      // get next backward transformation (transition for parametrization is inverse to that of the frames!!!)
      HFH h = mesh_.opposite_halfface_handle(pc.pathHalffaces[i]);

      HFH h_opp = mesh_.opposite_halfface_handle(h);

      // get transition function (param transforms inversely to frame-field)
      int tf = tq_.inverse_transition_idx(transition_hfprop_[h]);
      auto R = tq_.transition_matrix_int(tf);

      // get vertex of halfface
      VH vh = mesh_.halfedge(mesh_.halfface(h).halfedges()[0]).to_vertex();
      CH ch = mesh_.incident_cell(h);
      CH ch_opp = mesh_.incident_cell(h_opp);

      Vec3 t = ovm2eigen(igm_cprop_[ch_opp][vh]) - R * ovm2eigen(igm_cprop_[ch][vh]);

      p = (R * p).eval() + t;
    }

    Vec3 q = ovm2eigen(igm_cprop_[pc.tetFrom][pc.vFrom]);
    Vec3 d = p - q;
    Vec3 diff = d - ovm2eigen(pc.offset);

    if(std::abs(diff[0]) > 1 || std::abs(diff[1]) > 1 || std::abs(diff[2]) > 1)
        std::cerr << "path constraint " << count << " has large diff to seamless map of " << diff << std::endl;

    ++count;
  }

  std::cerr << "#invalid path constraints = " << n_invalid << ", #valid path constraints = " << n_valid << std::endl;
}

//-----------------------------------------------------------------------------


template<class TetMeshT>
template<class ExportTetMeshT>
void
Parametrization3DT<TetMeshT>::
export_quantization_constraint_paths(ExportTetMeshT &_export_mesh)
{

  _export_mesh.clear();
  auto pci_ep = _export_mesh.template request_edge_property<int>("PathConstraintIndex");
  _export_mesh.set_persistent(pci_ep, true);

  int i=0;
  for(auto pc : quantization_path_constraints_)
  {
    // start segment
    Mat3x2 P;
    P.col(0) = ovm2eigen(mesh_.vertex(pc.vFrom));
    P.col(1) = ovm2eigen(mesh_.barycenter(pc.tetFrom));
    auto eh = add_edge(_export_mesh, P);
    pci_ep[eh] = i;

    // end segment
    P.col(0) = ovm2eigen(mesh_.vertex(pc.vTo));
    P.col(1) = ovm2eigen(mesh_.barycenter(pc.tetTo));
    eh = add_edge(_export_mesh, P);
    pci_ep[eh] = i;

    // intermediate segments
    for(auto hfh : pc.pathHalffaces)
      if(!mesh_.is_boundary(mesh_.face_handle(hfh)))
      {
        P.col(0) = ovm2eigen(mesh_.barycenter(mesh_.incident_cell(hfh)));
        P.col(1) = ovm2eigen(mesh_.barycenter(mesh_.incident_cell(mesh_.opposite_halfface_handle(hfh))));
        eh = add_edge(_export_mesh, P);
        pci_ep[eh] = i;
      }
    // increment index
    ++i;
  }
}


//=============================================================================
} // namespace AlgoHex
//=============================================================================

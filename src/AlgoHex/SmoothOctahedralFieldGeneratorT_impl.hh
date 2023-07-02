/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define SMOOTHOCTAHEDRALFIELDGENERATOR_C

//== INCLUDES =================================================================
#include "OpenNL/OpenNL_psm.h"
#include "SmoothOctahedralFieldGeneratorT.hh"
#include "MeshGeometry.hh"
#include "FieldConstraints.hh"
#include <AlgoHex/Util/StopWatch.hh>
#include <AlgoHex/Stopwatches.hh>

#include <iostream>

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== IMPLEMENTATION ==========================================================

template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::
solve_spherical_harmonic_coefficients(const double _penalty)
{
  ScopedStopWatch sw(sw::field_opt);
  {
    ScopedStopWatch sw(sw::linsys_setup);

    size_t n_normal_constraints = get_vertex_normal_constraints(tmesh_).size();

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, 9 * num_vertices_ + 2 * n_normal_constraints);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlSolverParameteri(NL_SOLVER, 0);

    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    add_smooth_term();

    add_normal_constraints(_penalty);
    add_full_constraints(_penalty);

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
  }


  std::cout << "\n#####Solving for Spherical Harmonic Coefficients..." << std::endl;

  {
    ScopedStopWatch sw(sw::linsys_solve);
    nlSolve();
  }

  for (const auto vh: tmesh_.vertices())
  {
    for (int d = 0; d < 9; d++)
      shcs_[vh][d] = nlGetVariable(9 * vh.idx() + d);
  }

  nlDeleteContext(nlGetCurrent());
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::
iterate(int _N, double _penalty, const bool _original_projection)
{
  ScopedStopWatch sw(sw::field_opt);
  size_t n_normal_constraints = get_vertex_normal_constraints(tmesh_).size();

  double prev_energy = compute_energy();
  std::cout << "#####Start Iterating... \niteration 0 energy: " << prev_energy << std::endl;

  for (int I = 1; I <= _N; I++)
  {
    {
      ScopedStopWatch sw(sw::linsys_setup);
      nlNewContext();
      nlSolverParameteri(NL_NB_VARIABLES, 12 * num_vertices_ + 2 * n_normal_constraints);
      nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
      nlSolverParameteri(NL_SOLVER, 0);

      nlBegin(NL_SYSTEM);
      nlBegin(NL_MATRIX);

      add_smooth_term();

      add_normal_constraints(_penalty);
      add_full_constraints(_penalty);

      add_local_optim_constraints(n_normal_constraints, _penalty);

      nlEnd(NL_MATRIX);
      nlEnd(NL_SYSTEM);
    }


    {
      ScopedStopWatch sw(sw::linsys_solve);
      nlSolve();
    }

    for (const auto vh: tmesh_.vertices())
    {
      for (int d = 0; d < 9; d++)
        shcs_[vh][d] = nlGetVariable(9 * vh.idx() + d);
    }

    nlDeleteContext(nlGetCurrent());

    project(_original_projection);

    double energy = compute_energy();
    if ((prev_energy - energy) / prev_energy < 1e-5 || energy < 1e-7)
    {
      std::cout << "iteration " << I << " energy: " << energy << std::endl;
      std::cout << "Stop iteration due to convergence after " << I << " of " << _N << " iterations." << std::endl;
      break;
    }
    else
    {
      std::cout << "iteration " << I << " energy: " << energy << std::endl;
      prev_energy = energy;
    }
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::
project(const bool _original_projection)
{
  ScopedStopWatch sw_total(sw::field_opt);
  ScopedStopWatch sw(sw::projection);

  if (_original_projection)
  {
    projector_.set_seed_selection(SHProjectorRay::ClosestSeed);
  }
  else
  {
    projector_.set_seed_selection(SHProjectorRay::BestSeed);
  }

#pragma omp parallel for schedule (dynamic, 16)
  for (int i = 0; i < (int) num_vertices_; ++i)
  {
    VH vh(i);
    auto result = projector_.project(shcs_[vh].normalized());
    shcs_p_[vh] = result.shc;
    quaternions_[vh] = result.q;
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::
save_vertex_quaternion(const std::string &_filename) const
{
  std::ofstream f_write(_filename);

  for (const auto vh: tmesh_.vertices())
  {
    // write to file
    f_write << quaternions_[vh].w() << " ";
    f_write << quaternions_[vh].x() << " ";
    f_write << quaternions_[vh].y() << " ";
    f_write << quaternions_[vh].z() << " ";
  }

  f_write.close();
}


//================================================================================================================//

template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::
add_smooth_term() const
{
  ScopedStopWatch sw(sw::linsys_setup_smooth);
  for (const auto eh: tmesh_.edges())
  {
    for (int d = 0; d < 9; d++)
    {
      nlBegin(NL_ROW);
      nlCoefficient(9 * (tmesh_.edge(eh).from_vertex().idx()) + d, 1.0);
      nlCoefficient(9 * (tmesh_.edge(eh).to_vertex().idx()) + d, -1.0);
      nlEnd(NL_ROW);
    }
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::
add_normal_constraints(const double _penalty)
{
  ScopedStopWatch sw(sw::linsys_setup_normal_constr);
  int count = 0;
  for (const auto &ncc: normal_constraint_coeffs_)
  {
    for (int d = 0; d < 9; d++)
    {
      nlRowScaling(_penalty);
      nlBegin(NL_ROW);
      nlCoefficient(9 * ncc.vh.idx() + d, 1.);
      nlCoefficient(9 * num_vertices_ + 2 * count + 0, ncc.h0[d]);
      nlCoefficient(9 * num_vertices_ + 2 * count + 1, ncc.h8[d]);
      nlRightHandSide(ncc.h4[d]);
      nlEnd(NL_ROW);
    }
    ++count;
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::
add_full_constraints(const double _penalty)
{
  ScopedStopWatch sw(sw::linsys_setup_full_constr);
  auto &vertex_full_constraints = get_vertex_full_constraints(tmesh_);
  for (const auto &[vh, q]: vertex_full_constraints)
  {
    // TODO: convert to sph coeffs only once in the beginning
    auto sph = SHCoeffs::from_quaternion(q);
    for (int d = 0; d < 9; d++)
    {
      nlRowScaling(_penalty);
      nlBegin(NL_ROW);
      nlCoefficient(9 * vh.idx() + d, 1);
      nlRightHandSide(sph[d]);
      nlEnd(NL_ROW);
    }
  }

}

//-----------------------------------------------------------------------------

template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::
add_local_optim_constraints(const int _num_cv, const double _penalty) const
{
  ScopedStopWatch sw(sw::linsys_setup_lin);
  int a_idx = 0;
  for (const auto vh: tmesh_.vertices())
  {
    const auto shc = shcs_p_[vh];
    const auto c_x = Ex(shc);
    const auto c_y = Ey(shc);
    const auto c_z = Ez(shc);

    for (int d = 0; d < 9; d++)
    {
      nlBegin(NL_ROW);
      nlCoefficient(9 * vh.idx() + d, _penalty);
      nlCoefficient(9 * num_vertices_ + 2 * _num_cv + 3 * a_idx + 0, -_penalty * c_x[d]);
      nlCoefficient(9 * num_vertices_ + 2 * _num_cv + 3 * a_idx + 1, -_penalty * c_y[d]);
      nlCoefficient(9 * num_vertices_ + 2 * _num_cv + 3 * a_idx + 2, -_penalty * c_z[d]);
      nlRightHandSide(_penalty * shcs_p_[vh][d]);
      nlEnd(NL_ROW);
    }
    a_idx++;
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
double SmoothOctahedralFieldGeneratorT<MeshT>::
compute_energy() const
{
  double energy = 0.0;
  for (const auto eh: tmesh_.edges())
  {
    for (int i = 0; i < 9; i++)
    {
      double diff = (shcs_p_[tmesh_.edge(eh).from_vertex()][i] -
                     shcs_p_[tmesh_.edge(eh).to_vertex()][i]);
      energy += diff * diff;
    }
  }

  return energy;
}

//-----------------------------------------------------------------------------

template<class MeshT>
void
SmoothOctahedralFieldGeneratorT<MeshT>::
cache_constrained_9d_coefficients()
{
  auto vnc = get_vertex_normal_constraints(tmesh_);
  normal_constraint_coeffs_.clear();
  normal_constraint_coeffs_.reserve(vnc.size());

  for (const auto &[vh, normal]: vnc)
  {
    normal_constraint_coeffs_.emplace_back(vh, normal);
  }
}

//-----------------------------------------------------------------------------
template<class MeshT>
void SmoothOctahedralFieldGeneratorT<MeshT>::print_timings() const
{
  std::cout << sw::field_opt << std::endl;
#if 0
  std::cout << "Time spent in SPH projection:      " << sw::projection_.elapsed_ms() << "ms" << std::endl;
  std::cout << "Time spent in Linear System setup: " << sw::linsys_setup_.elapsed_ms() << "ms" << std::endl;
  std::cout << "Time spent in Linear System solve: " << sw::linsys_solve_.elapsed_ms() << "ms" << std::endl;
#endif
}
//-----------------------------------------------------------------------------

//Convention: X(alpha)->Y(beta)->Z(gamma)
inline Vec3d get_euler_angle(const double _x, const double _y, const double _z)
{
  Vec3d angles(0, 0, 0);
  if (std::abs(_z) < .99)
  {
    Eigen::Vector3d dir(_x, _y, _z), nz(0, 0, 1);
    Mat3d m;
    m = Eigen::Quaterniond().setFromTwoVectors(nz, dir);
    //(gamma, beta, alpha)
    angles = m.eulerAngles(2, 1, 0);
  }

  return {angles(2), angles(1), angles(0)};
}

//-----------------------------------------------------------------------------

template<typename Point>
NormalConstraintCoeff::NormalConstraintCoeff(
        const OVM::VertexHandle &_vh,
        const Point &normal)
        : vh(_vh)
{
  auto angles = get_euler_angle(normal[0], normal[1], normal[2]);

  h0 = {std::sqrt(5. / 12.), 0, 0, 0, 0, 0, 0, 0, 0};
  h0.rotate_xyz(angles);

  h4 = {0, 0, 0, 0, std::sqrt(7. / 12.), 0, 0, 0, 0};
  h4.rotate_xyz(angles);

  h8 = {0, 0, 0, 0, 0, 0, 0, 0, std::sqrt(5. / 12.)};
  h8.rotate_xyz(angles);

}


//=============================================================================
} // namespace AlgoHex
//=============================================================================


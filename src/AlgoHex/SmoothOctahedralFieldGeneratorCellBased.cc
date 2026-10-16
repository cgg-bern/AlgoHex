/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/


//== INCLUDES =================================================================
#include "OpenNL/OpenNL_psm.h"
#include "SmoothOctahedralFieldGeneratorCellBased.hh"
#include "SmoothOctahedralFieldGenerator.hh"
#include "TypeDef.hh"
#include "FieldConstraints.hh"
#include <AlgoHex/Stopwatches.hh>
#include <libTimekeeper/StopWatchPrinting.hh>

#include <fstream>

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== IMPLEMENTATION ==========================================================

template<typename Point>
PartialConstraintCoeff::PartialConstraintCoeff(
        const OVM::CellHandle &_ch,
        const Point &normal)
        : ch(_ch)
{
  auto angles = get_euler_angle(normal[0], normal[1], normal[2]);

  h0 = {std::sqrt(5. / 12.), 0, 0, 0, 0, 0, 0, 0, 0};
  h0.rotate_xyz(angles);

  h4 = {0, 0, 0, 0, std::sqrt(7. / 12.), 0, 0, 0, 0};
  h4.rotate_xyz(angles);

  h8 = {0, 0, 0, 0, 0, 0, 0, 0, std::sqrt(5. / 12.)};
  h8.rotate_xyz(angles);

}
template ALGOHEX_EXPORT PartialConstraintCoeff::PartialConstraintCoeff(const OVM::CellHandle &, const OVM::Geometry::Vec3d &);
template ALGOHEX_EXPORT PartialConstraintCoeff::PartialConstraintCoeff(const OVM::CellHandle &, const Eigen::Vector3d &);

SmoothOctahedralFieldGeneratorCellBased::SmoothOctahedralFieldGeneratorCellBased(TetMesh &_tmesh) :
        tmesh_(_tmesh),
        shcs_(tmesh_.request_cell_property<SHCoeffs>("cell spherical harmonic coefficients")),
        shcs_p_(tmesh_.request_cell_property<SHCoeffs>("scaled cell spherical harmonic coefficients")),
        quaternions_(tmesh_.request_cell_property<Quaternion>("FrameFieldQuaternions", {0., 0., 0., 0.}))
{
  num_cells_ = tmesh_.n_cells();

  tmesh_.set_persistent(shcs_p_, true);
  tmesh_.set_persistent(quaternions_, true);

  cache_constrained_9d_coefficients();
}

void SmoothOctahedralFieldGeneratorCellBased::
solve_spherical_harmonic_coefficients(const double _penalty)
{
  ScopedStopWatch sw(sw::field_opt);
  {
    ScopedStopWatch sw(sw::linsys_setup);

    size_t n_partial_constraints = get_cell_partial_constraints(tmesh_).size();

    nlNewContext();
    nlSolverParameteri(NL_NB_VARIABLES, 9 * num_cells_ + 2 * n_partial_constraints);
    nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
    nlSolverParameteri(NL_SOLVER, 0);

    nlBegin(NL_SYSTEM);
    nlBegin(NL_MATRIX);

    add_smooth_term();

    add_partial_constraints(_penalty);

    nlEnd(NL_MATRIX);
    nlEnd(NL_SYSTEM);
  }


  std::cout << "\n#####Solving for Spherical Harmonic Coefficients..." << std::endl;

  {
    ScopedStopWatch sw(sw::linsys_solve);
    nlSolve();
  }

  for (const auto ch: tmesh_.cells())
  {
    for (int d = 0; d < 9; d++)
      shcs_[ch][d] = nlGetVariable(9 * ch.idx() + d);
  }

  nlDeleteContext(nlGetCurrent());
}

//-----------------------------------------------------------------------------

void SmoothOctahedralFieldGeneratorCellBased::
iterate(int _N, double _penalty, const bool _original_projection)
{
  ScopedStopWatch sw(sw::field_opt);
  size_t n_partial_constraints = get_cell_partial_constraints(tmesh_).size();

  double prev_energy = compute_energy();
  std::cout << "#####Start Iterating... \niteration 0 energy: " << prev_energy << std::endl;

  for (int I = 1; I <= _N; I++)
  {
    {
      ScopedStopWatch sw(sw::linsys_setup);
      nlNewContext();
      nlSolverParameteri(NL_NB_VARIABLES, 12 * num_cells_ + 2 * n_partial_constraints);
      nlSolverParameteri(NL_LEAST_SQUARES, NL_TRUE);
      nlSolverParameteri(NL_SOLVER, 0);

      nlBegin(NL_SYSTEM);
      nlBegin(NL_MATRIX);

      add_smooth_term();

      add_partial_constraints(_penalty);

      add_local_optim_constraints(n_partial_constraints, _penalty);

      nlEnd(NL_MATRIX);
      nlEnd(NL_SYSTEM);
    }


    {
      ScopedStopWatch sw(sw::linsys_solve);
      nlSolve();
    }

    for (const auto ch: tmesh_.cells())
    {
      for (int d = 0; d < 9; d++)
        shcs_[ch][d] = nlGetVariable(9 * ch.idx() + d);
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

void SmoothOctahedralFieldGeneratorCellBased::
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
  for (int i = 0; i < (int) num_cells_; ++i)
  {
    CH ch(i);
    auto result = projector_.project(shcs_[ch].normalized());
    shcs_p_[ch] = result.shc;
    quaternions_[ch] = result.q;
  }
}

//-----------------------------------------------------------------------------

void SmoothOctahedralFieldGeneratorCellBased::
save_cell_quaternion(const std::string &_filename) const
{
  std::ofstream f_write(_filename);

  for (const auto ch: tmesh_.cells())
  {
    // write to file
    f_write << quaternions_[ch].w() << " ";
    f_write << quaternions_[ch].x() << " ";
    f_write << quaternions_[ch].y() << " ";
    f_write << quaternions_[ch].z() << " ";
  }

  f_write.close();
}


//================================================================================================================//

void SmoothOctahedralFieldGeneratorCellBased::
add_smooth_term() const
{
  ScopedStopWatch sw(sw::linsys_setup_smooth);
  for (const auto fh: tmesh_.faces())
  {
    if (!tmesh_.is_boundary(fh))
    {
      CH ch0 = tmesh_.incident_cell(tmesh_.halfface_handle(fh, 0));
      CH ch1 = tmesh_.incident_cell(tmesh_.halfface_handle(fh, 1));

      for (int d = 0; d < 9; d++)
      {
        nlBegin(NL_ROW);
        nlCoefficient(9 * ch0.idx() + d, 1.0);
        nlCoefficient(9 * ch1.idx() + d, -1.0);
        nlEnd(NL_ROW);
      }
    }
  }
}

//-----------------------------------------------------------------------------

void SmoothOctahedralFieldGeneratorCellBased::
add_partial_constraints(const double _penalty)
{
  ScopedStopWatch sw(sw::linsys_setup_normal_constr);
  int count = 0;
  for (const auto &ncc: partial_constraint_coeffs_)
  {
    for (int d = 0; d < 9; d++)
    {
      nlRowScaling(_penalty);
      nlBegin(NL_ROW);
      nlCoefficient(9 * ncc.ch.idx() + d, 1.);
      nlCoefficient(9 * num_cells_ + 2 * count + 0, ncc.h0[d]);
      nlCoefficient(9 * num_cells_ + 2 * count + 1, ncc.h8[d]);
      nlRightHandSide(ncc.h4[d]);
      nlEnd(NL_ROW);
    }
    ++count;
  }
}

//-----------------------------------------------------------------------------

void SmoothOctahedralFieldGeneratorCellBased::
add_local_optim_constraints(const int _num_cv, const double _penalty) const
{
  ScopedStopWatch sw(sw::linsys_setup_lin);
  int a_idx = 0;
  for (const auto &ch: tmesh_.cells())
  {
    const auto shc = shcs_p_[ch];
    const auto c_x = Ex(shc);
    const auto c_y = Ey(shc);
    const auto c_z = Ez(shc);

    for (int d = 0; d < 9; d++)
    {
      nlBegin(NL_ROW);
      nlCoefficient(9 * ch.idx() + d, _penalty);
      nlCoefficient(9 * num_cells_ + 2 * _num_cv + 3 * a_idx + 0, -_penalty * c_x[d]);
      nlCoefficient(9 * num_cells_ + 2 * _num_cv + 3 * a_idx + 1, -_penalty * c_y[d]);
      nlCoefficient(9 * num_cells_ + 2 * _num_cv + 3 * a_idx + 2, -_penalty * c_z[d]);
      nlRightHandSide(_penalty * shcs_p_[ch][d]);
      nlEnd(NL_ROW);
    }
    a_idx++;
  }
}

//-----------------------------------------------------------------------------

double SmoothOctahedralFieldGeneratorCellBased::
compute_energy() const
{
  double energy = 0.0;
  for (const auto fh: tmesh_.faces())
  {
    if (!tmesh_.is_boundary(fh))
    {
      CH ch0 = tmesh_.incident_cell(tmesh_.halfface_handle(fh, 0));
      CH ch1 = tmesh_.incident_cell(tmesh_.halfface_handle(fh, 1));

      for (int i = 0; i < 9; i++)
      {
        double diff = (shcs_p_[ch0][i] - shcs_p_[ch1][i]);
        energy += diff * diff;
      }
    }
  }

  return energy;
}

//-----------------------------------------------------------------------------

void
SmoothOctahedralFieldGeneratorCellBased::
cache_constrained_9d_coefficients()
{
  auto cnc = get_cell_partial_constraints(tmesh_);
  partial_constraint_coeffs_.clear();
  partial_constraint_coeffs_.reserve(cnc.size());

  for (const auto &[ch, normal]: cnc)
  {
    partial_constraint_coeffs_.emplace_back(ch, normal);
  }
}

//-----------------------------------------------------------------------------

void SmoothOctahedralFieldGeneratorCellBased::print_timings() const
{
  std::cout << sw::field_opt << std::endl;
#if 0
  std::cout << "Time spent in SPH projection:      " << sw::projection_.elapsed_ms() << "ms" << std::endl;
  std::cout << "Time spent in Linear System setup: " << sw::linsys_setup_.elapsed_ms() << "ms" << std::endl;
  std::cout << "Time spent in Linear System solve: " << sw::linsys_solve_.elapsed_ms() << "ms" << std::endl;
#endif

}

} // namespace AlgoHex
//=============================================================================

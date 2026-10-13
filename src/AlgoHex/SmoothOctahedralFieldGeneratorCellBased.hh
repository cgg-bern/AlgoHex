/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

//== INCLUDES =================================================================

#include "TypeDef.hh"

#include <AlgoHex/SphericalHarmonics/SHCoeffs.hh>
#include <AlgoHex/SphericalHarmonics/SHProjectorRay.hh>


//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================
namespace OVM = OpenVolumeMesh;

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================



struct PartialConstraintCoeff
{
  template<typename Point>
  PartialConstraintCoeff(const OVM::CellHandle &_ch, const Point &normal);

  OVM::CellHandle ch;
  SHCoeffs h0, h4, h8;
};

extern template PartialConstraintCoeff::PartialConstraintCoeff(const OVM::CellHandle &, const OVM::Geometry::Vec3d &);
extern template PartialConstraintCoeff::PartialConstraintCoeff(const OVM::CellHandle &, const Eigen::Vector3d &);

/** \class SmoothOctahedralFieldGenerator SmoothOctahedralFieldGeneratorCellBased.hh

    @brief Implementation of boundary aligned smooth octahedral field generation method by Nicolas Ray and Dmitry Sokolov.

    A more elaborate description follows.
*/

class ALGOHEX_EXPORT SmoothOctahedralFieldGeneratorCellBased
{
public:
  using Mesh = TetMesh;
  using Point = typename Mesh::PointT;

  SmoothOctahedralFieldGeneratorCellBased(TetMesh &_tmesh);

  ~SmoothOctahedralFieldGeneratorCellBased()
  {
  }

public:
  //get the spherical harmonic coefficients by solving a least square problem
  void solve_spherical_harmonic_coefficients(const double _penalty);

  void project(const bool _original_projection = true);

  void iterate(int _N, double _penalty, const bool _original_projection = true);

  void save_cell_quaternion(const std::string &_filename) const;

  void print_timings() const;

private:
  void add_smooth_term() const;

  void add_partial_constraints(const double _penalty);

  void add_local_optim_constraints(const int _num_cv, const double _penalty) const;

  double compute_energy() const;

  void cache_constrained_9d_coefficients();

  Mesh &tmesh_;

  //9D coefficients before projection
  CP<SHCoeffs> shcs_;

  //after projection
  CP<SHCoeffs> shcs_p_;

  //vertex quaternion
  CP<Quaternion> quaternions_;

  //cell number
  size_t num_cells_;

  std::vector<PartialConstraintCoeff> partial_constraint_coeffs_;
  SHProjectorRay projector_;

};

//=============================================================================

} // namespace AlgoHex
//=============================================================================

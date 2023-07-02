/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

//== INCLUDES =================================================================

#include <fstream>
#include "TypeDef.hh"

#include <AlgoHex/SphericalHarmonics/SHCoeffs.hh>
#include <AlgoHex/SphericalHarmonics/SHProjectorRay.hh>
#include <AlgoHex/Util/ScopedStopWatch.hh>


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


/** \class SmoothOctahedralFieldGeneratorT SmoothOctahedralFieldGeneratorCellBasedT.hh

    @brief Implementation of boundary aligned smooth octahedral field generation method by Nicolas Ray and Dmitry Sokolov.

    A more elaborate description follows.
*/

template<class MeshT>
class SmoothOctahedralFieldGeneratorCellBasedT
{
public:
  using Mesh = MeshT;
  using Point = typename Mesh::PointT;


  SmoothOctahedralFieldGeneratorCellBasedT(MeshT &_tmesh) :
          tmesh_(_tmesh),
          shcs_(tmesh_.template request_cell_property<SHCoeffs>("cell spherical harmonic coefficients")),
          shcs_p_(tmesh_.template request_cell_property<SHCoeffs>("scaled cell spherical harmonic coefficients")),
          quaternions_(tmesh_.template request_cell_property<Quaternion>("FrameFieldQuaternions", {0., 0., 0., 0.}))
  {
    num_cells_ = tmesh_.n_cells();

    tmesh_.set_persistent(shcs_p_, true);
    tmesh_.set_persistent(quaternions_, true);

    cache_constrained_9d_coefficients();
  }

  ~SmoothOctahedralFieldGeneratorCellBasedT()
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


private:
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
#if defined(INCLUDE_TEMPLATES) && !defined(SMOOTHOCTAHEDRALFIELDGENERATORCELLBASED_C)
#define SMOOTHOCTAHEDRALFIELDGENERATORCELLBASED_TEMPLATES

#include "SmoothOctahedralFieldGeneratorCellBasedT_impl.hh"

#endif
//=============================================================================


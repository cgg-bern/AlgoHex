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



struct NormalConstraintCoeff
{
  NormalConstraintCoeff(const OVM::VertexHandle &_vh, const Point &normal);

  OVM::VertexHandle vh;
  SHCoeffs h0, h4, h8;
};


/** \class SmoothOctahedralFieldGenerator SmoothOctahedralFieldGenerator.hh

    @brief Implementation of boundary aligned smooth octahedral field generation method by Nicolas Ray and Dmitry Sokolov.

    A more elaborate description follows.
*/

class ALGOHEX_EXPORT SmoothOctahedralFieldGenerator
{
public:
  using Mesh = TetMesh;
  using Point = typename Mesh::PointT;

  SmoothOctahedralFieldGenerator(TetMesh &_tmesh);

  ~SmoothOctahedralFieldGenerator();

public:
  //get the spherical harmonic coefficients by solving a least square problem
  void solve_spherical_harmonic_coefficients(const double _penalty);

  void project(const bool _original_projection = true);

  void iterate(int _N, double _penalty, const bool _original_projection = true);

  void save_vertex_quaternion(const std::string &_filename) const;

  void print_timings() const;

private:
  void add_smooth_term() const;

  void add_normal_constraints(const double _penalty);

  void add_full_constraints(double penalty);

  void add_local_optim_constraints(const int _num_cv, const double _penalty) const;

  double compute_energy() const;

  void cache_constrained_9d_coefficients();

  Mesh &tmesh_;

//        //boolean halfface property for interior constraints
//        HFP<bool> interface_;

  //9D coefficients before projection
  VP<SHCoeffs> shcs_;

  //after projection
  VP<SHCoeffs> shcs_p_;

  //vertex quaternion
  VP<Quaternion> quaternions_;

  //vertex number
  size_t num_vertices_;

  std::vector<NormalConstraintCoeff> normal_constraint_coeffs_;
  SHProjectorRay projector_;

};

//=============================================================================

//Convention: X(alpha)->Y(beta)->Z(gamma)
Vec3d get_euler_angle(const double _x, const double _y, const double _z);

} // namespace AlgoHex
//=============================================================================

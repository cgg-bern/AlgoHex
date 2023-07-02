/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <AlgoHex/TypeDef.hh>

namespace AlgoHex
{
template<class MeshT>
class RemeshingAssist
{
public:
  using Point = typename MeshT::PointT;

  //compare the imrm energies
  static bool is_quality_improved(const std::vector<std::vector<Point> > &_new_cells_points,
                                  const std::vector<std::vector<Point> > &_old_cells_points,
                                  bool compare_max_energy = false);

  //compute the Inverse Mean Ratio Metric energy of a tet
  static double compute_IMRM_energy(const std::vector<Point> &_points);

  //compute the Gradient w.r.t. the first point
  static void compute_IMRM_gradient(const std::vector<Point> &_points, double *_gradient);

  //compute the Gradient w.r.t. the first point
  static void compute_IMRM_hessian(const std::vector<Point> &_points, double *_hessian);

  //check the orientation of the tetrahedron
  static bool is_tet_valid(const std::vector<Point> &_points);

  //boundary edge type w.r.t. dihedral angle:
  static int boundary_edge_type(double _angle);

  //sample points
  static void sample_points(const Point &_p0, const Point &_p1, const Point &_p2, const int _n, const int _N,
                            std::vector<Point> &_points);

private:
  //compute the Inverse Mean Ratio Metric energy of a group of tets that share the common vertex
  static double comformalIMRMEnergy(const double *_x);

  static void comformalIMRMGradient(const double *_x, double *_g);

  static void comformalIMRMHessian(const double *_x, double *_h);

//    public:
//        template <typename Vec3>
//        static Point_3 to_CGAL_point(const Vec3& _pt) {
//            return Point_3(_pt[0], _pt[1], _pt[2]);
//        }
//
//        template <typename Vec3>
//        static Point_3S to_CGAL_point_S(const Vec3& _pt) {
//            return Point_3S(_pt[0], _pt[1], _pt[2]);
//        }

};
}
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(REMESHINGASSIST_C)
#define REMESHINGASSIST_TEMPLATES

#include "RemeshingAssist_impl.hh"

#endif
//=============================================================================

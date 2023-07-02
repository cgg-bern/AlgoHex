/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
\*===========================================================================*/


//=============================================================================
//
//  CLASS PlaneT
//
//=============================================================================


#ifndef ACG_PLANE_HH
#define ACG_PLANE_HH


//== INCLUDES =================================================================

#include <cfloat>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include "Matrix4x4T.hh"


//== NAMESPACES ===============================================================


namespace ACG {
namespace Geometry {


//== CLASS DEFINITION =========================================================


/** \class PlaneT PlaneT.hh <ACG/Geometry/Types/Plane.hh>

    This class stores a plane in normal form and provides useful
    algorithms for left/right test, itersection, and so on.
**/

template <typename Scalar>
class PlaneT
{
public:

  /// typedefs
  typedef VectorT<Scalar, 3>  Vec3;
  typedef VectorT<Scalar, 4>  Vec4;
  typedef Matrix4x4T<Scalar>  Mat4x4;


  /// constructor: coefficients
  PlaneT( Scalar _a=0, Scalar _b=0, Scalar _c=0, Scalar _d=0 )
    : coeffs_(_a, _b, _c, _d)
  { HNF(); }


  /// constructor: origin and normal
  PlaneT( const Vec3& _o, const Vec3& _n )
    : coeffs_(_n[0], _n[1], _n[2], -(_n|_o))
  { HNF(); }


  /// constructor: 3 points
  PlaneT( const Vec3& _v0, const Vec3& _v1, const Vec3& _v2 )
  {
    Vec3 n = (_v1-_v0) % (_v2-_v0);
    coeffs_ = Vec4(n[0], n[1], n[2], -(n|_v0));
    HNF();
  }


  /// normal vector
  Vec3 normal() const { return Vec3(coeffs_[0], coeffs_[1], coeffs_[2]); }


  /// coeffitients
  const Vec4& coeffs() const { return coeffs_; }


  /// signed distance point-plane
  Scalar distance( const Vec3& _v ) const
  {
    return ( _v[0]*coeffs_[0] +
	     _v[1]*coeffs_[1] +
	     _v[2]*coeffs_[2] +
  	           coeffs_[3] );
  }


  /// predicate: above plane
  bool operator() ( const Vec3& _v ) const { return distance(_v) > 0; }


  // INTERSECTION
  enum IntersectionTarget { Line, LineSegment, Ray };

  // intersect with (infinite) line
  bool intersect_line( const Vec3& _v0, const Vec3& _v1,
		       Vec3& _v, Scalar& _t ) const
  { return intersect(_v0, _v1, _v, _t, Line); }

  // intersect with ray
  bool intersect_ray( const Vec3& _v0, const Vec3& _v1,
		      Vec3& _v, Scalar& _t ) const
  { return intersect(_v0, _v1, _v, _t, Ray); }

  // intersect with line segment
  bool intersect_linesegment( const Vec3& _v0, const Vec3& _v1,
			      Vec3& _v, Scalar& _t ) const
  { return intersect(_v0, _v1, _v, _t, LineSegment); }

  /**  \brief General intersection function.
    *
    * General intersection between a line/ray and the plane.
    *
    * @param _v0 start vector of the ray
    * @param _v1 end vector of the ray
    * @param _v returned intersection point
    * @param _t returned relative distance from the interection to _v0 compared to _v1
    * @param _target type of intersection to test
    * @return true if an intersection was found
    */
  bool intersect( const Vec3& _v0,
		  const Vec3& _v1,
		  Vec3& _v,
		  Scalar& _t,
		  IntersectionTarget _target )  const
  {
#define SGN(d) ((d>0.0) ? 1 : -1)

    Scalar  d0(distance(_v0)), d1(distance(_v1));
    Scalar  a0(fabs(d0)), a1(fabs(d1));

    // endpoint on plane
    if (a0 < FLT_MIN)  { _v = _v0; _t = 0.0; return true; }
    if (a1 < FLT_MIN)  { _v = _v1; _t = 1.0; return true; }

    // triv accept
    if (SGN(d0) != SGN(d1))
    {
      _t = (a0/(a0+a1));
      _v = _v0*(1.0-_t) + _v1*_t;
      return true;
    }

    // depends on target
    else
    {
      if (_target == LineSegment)  return false;

      if (fabs(d0-d1) < FLT_MIN) return false; // line parallel to plane
      else _t = d0/(d0-d1);

      if (_target == Ray && _t < 0.0)  return false;

      _v = _v0*(1.0-_t) + _v1*_t;
      return true;
    }
#undef SGN
  }


  /// affine transformation of the plane
  bool affine_transformation( const Mat4x4& _M )
  {
    Mat4x4 M(_M);
    if (!M.invert()) return false;
    M.transpone();
    affineTransformation_precomp(M);
    return true;
  }


  /// affine transformation of the plane (M^{-T} precomputed)
  void affine_transformation_precomp( const Mat4x4& _M_inverseTransposed )
  {
    coeffs_ = _M_inverseTransposed*coeffs_;
    HNF();
  }


private:

  void HNF() {
    Scalar n = normal().norm();
    if (n != 0.0) coeffs_ /= n;
  }

  Vec4 coeffs_;
};



typedef PlaneT<float>   Planef;
typedef PlaneT<double>  Planed;


//=============================================================================
} // namespace Geometry
} // namespace ACG
//=============================================================================
#endif // ACG_PLANE_HH defined
//=============================================================================


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



#ifndef GEO_ALGORITHMS_HH
#define GEO_ALGORITHMS_HH


//== INCLUDES =================================================================

#include <cfloat>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <vector>
#include <iostream>



namespace ACG {
namespace Geometry {


//== 3D STUFF =================================================================

 

/// return circumcenter of tetrahedron (_v0,_v1,_v2,_v3)
template<typename Scalar>
bool
circumCenter( const VectorT<Scalar,3>&  _v0,
	      const VectorT<Scalar,3>&  _v1,
	      const VectorT<Scalar,3>&  _v2,
	      const VectorT<Scalar,3>&  _v3,
	      VectorT<Scalar,3>&        _result );


/// return squared radius of circumcircle of tetrahedron (_v0,_v1,_v2,_v3)
template<typename Scalar>
Scalar
circumRadiusSquared( const VectorT<Scalar,3>&  _v0,
		     const VectorT<Scalar,3>&  _v1,
		     const VectorT<Scalar,3>&  _v2,
		     const VectorT<Scalar,3>&  _v3 )
{
  VectorT<Scalar,3> cc;
  return circumCenter(_v0, _v1, _v2, _v3, cc) ? (cc-_v0).sqrnorm() : FLT_MAX;
}


/// return radius of circumcircle of tetrahedron (_v0,_v1,_v2,_v3)
template<typename Scalar>
Scalar
circumRadius( const VectorT<Scalar,3>&  _v0,
	      const VectorT<Scalar,3>&  _v1,
	      const VectorT<Scalar,3>&  _v2,
	      const VectorT<Scalar,3>&  _v3 )
{
  return sqrt(circumRadiusSquared(_v0, _v1, _v2, _v3));
}



/** \brief Get rotation axis and signed angle of rotation between two vectors
 *
 * Get two vectors, _v0 and _v1, and compute rotation axis _v0 % _v1
 * as well the angle between _v0 and _v1. Note that the angle is always negative.
 * We consider the rotation to be performed from _v0 to _v1.
 * The angle is given in degree if not explicitly demanded in radiant
 * (pass false as fifth parameter).
 *
 * @param _v0 The first vector
 * @param _v1 The second vector
 * @param _axis A reference to a vector in which the rotation axis is stored
 * @param _angle A reference to a scalar type in which the signed angle is stored ( in degree)
 * @param _degree Indicates whether the angle should be given in degree or radiant
 */
template<typename Scalar>
bool
rotationOfTwoVectors( const VectorT<Scalar,3>&  _v0,
                      const VectorT<Scalar,3>&  _v1,
                      VectorT<Scalar,3>&  _axis,
                      Scalar& _angle,
                      bool _degree = true);


/** \brief find a vector that's perpendicular to _v
 *
 * This function takes a vector and  generates a new arbitrary
 * vector that is perpendicular to the input vector.
 *
 * @param _v Input vector
 * @return Perpendicular vector
 */
template <typename Scalar>
VectorT<Scalar,3>
perpendicular( const VectorT<Scalar,3>&  _v );


/**  \brief Intersect a ray and a triangle.
  *
  * Computes the intersection point between a ray and a triangle. The orientation of the triangle
  * does not matter. The distance returned in t will be negative if the triangle is not in the
  * direction given but in the opposite direction.
  *
  * @param _o origin of the ray
  * @param _dir direction vector of the ray
  * @param _v0 first point of the triangle
  * @param _v1 second point of the triangle
  * @param _v2 third point of the triangle
  * @param _t returned distance from the origin to the intersection, in units of _dir ( negative if before origin!)
  * @param _u returned first barycentric coordinate of the intersection point in the triangle
  * @param _v returned second barycentric coordinate of the intersection point in the triangle
  * @return true if an intersection was found
  */
template<typename Vec>
bool
triangleIntersection( const Vec&  _o,
                      const Vec&  _dir,
                      const Vec&  _v0,
                      const Vec&  _v1,
                      const Vec&  _v2,
                      typename Vec::value_type& _t,
                      typename Vec::value_type& _u,
                      typename Vec::value_type& _v );
      

/**  \brief Intersect a ray and an axis aligned bounding box
  *
  * Computes the intersection point between a ray and an axis aligned bounding box
  *
  * @param _o     Origin of the ray
  * @param _dir   direction vector of the ray
  * @param _bbmin lower left front corner of the bounding box
  * @param _bbmax upper right back corner of the bounding box
  * @param _t0 if there was an intersection, this value marks the entry point
  * @param _t1 if there was an intersection, this value marks the exit point
  * @return       true if an intersection was found
  */
template<typename Vec>
bool
axisAlignedBBIntersection( const Vec&  _o,
		                       const Vec&  _dir,
		                       const Vec& _bbmin,
		                       const Vec& _bbmax,
		                       typename Vec::value_type& _t0,
		                       typename Vec::value_type& _t1 );


//== 2D STUFF =================================================================

/// orientation of point _p w.r.t. line through _v0,_v1 in 2D
template<typename Scalar>
Scalar
pointLineOrientation( const VectorT<Scalar,2>&  _p,
		      const VectorT<Scalar,2>&  _v0,
		      const VectorT<Scalar,2>&  _v1 )
{
  VectorT<Scalar,2> d1(_p-_v0), d2(_v1-_v0);
  return (d1[0]*d2[1]-d1[1]*d2[0]);
}


/// are 3 vertices in counterclockwise order? in 2D
template<typename Scalar>
bool
isCCW( const VectorT<Scalar,2>&  _v0,
       const VectorT<Scalar,2>&  _v1,
       const VectorT<Scalar,2>&  _v2 )
{
  return ( pointLineOrientation(_v0, _v1, _v2) < Scalar(0.0) );
}


/// are 3 vertices in clockwise order? in 2D
template<typename Scalar>
bool
isCW( const VectorT<Scalar,2>&  _v0,
      const VectorT<Scalar,2>&  _v1,
      const VectorT<Scalar,2>&  _v2 )
{
  return ( pointLineOrientation(_v0, _v1, _v2) > Scalar(0.0) );
}


/// intersect two line segments (_v0,_v1) and (_v2,_v3)
template<typename Scalar>
bool
lineIntersection( const VectorT<Scalar,2>&  _v0,
		  const VectorT<Scalar,2>&  _v1,
		  const VectorT<Scalar,2>&  _v2,
		  const VectorT<Scalar,2>&  _v3,
		  Scalar& _t1,
		  Scalar& _t2 );


//===========================================================================
/** @name Distance Functions ( N-Dimensional )
* @{ */
//===========================================================================     


/** \brief squared distance from point _p to line segment (_v0,_v1)
 *
 * @param _p     Point to test
 * @param _v0    Start of line segment
 * @param _v1    End of line segment
 * @param _min_v Pointer to vector, to get the closest point or 0 if it's not required
 * @return Distance to line segment
 */
template<class Vec>
typename Vec::value_type
distPointLineSquared( const Vec& _p,
                      const Vec& _v0,
                      const Vec& _v1,
                      Vec*       _min_v=0);



/** \brief Compute distance from point to line segment
 *
 * Compute the distance from a point p to a line segment and possibly return
 * closest point on segment
 *
 * @param _p   Point to test
 * @param _v0  Start of line segment
 * @param _v1  End of line segment
 * @param _min_v Pointer to vector, to get the closest point or 0 if it's not required
 * @return Distance to line segment
 */
template<class Vec>
typename Vec::value_type
distPointLine( const Vec& _p,
               const Vec& _v0,
               const Vec& _v1,
               Vec*       _min_v=0 )
{ return sqrt(distPointLineSquared(_p, _v0, _v1, _min_v)); }


/// squared distance from point _p to triangle (_v0, _v1, _v2)
template <class Vec>
typename Vec::value_type
distPointTriangleSquared( const Vec& _p,
                          const Vec& _v0,
                          const Vec& _v1,
                          const Vec& _v2,
                          Vec& _nearestPoint );

/** \brief squared distance from point _p to triangle (_v0, _v1, _v2)
*
*  In the stable version the distance to the longest edge 
*  is returned if the triangle is degenerate.
*
* @param _p   point to test against triangle
* @param _v0  First point of triangle
* @param _v1  Second point of triangle
* @param _v2  Third point of triangle
* @return     Computed distance
*/
template <class Vec>
typename Vec::value_type
distPointTriangleSquaredStable( const Vec& _p,
                                const Vec& _v0,
                                const Vec& _v1,
                                const Vec& _v2,
                                Vec& _nearestPoint );

/// distance from point _p to triangle (_v0, _v1, _v2)
template <class Vec>
typename Vec::value_type
distPointTriangle( const Vec& _p,
                   const Vec& _v0,
                   const Vec& _v1,
                   const Vec& _v2,
                   Vec& _nearestPoint )
{ return sqrt(distPointTriangleSquared(_p, _v0, _v1, _v2, _nearestPoint)); }

/** \brief distance from point _p to triangle (_v0, _v1, _v2)
*
*   In the stable version the distance to the longest edge 
*   is returned if the triangle is degenerate.
* 
* @param _v0  First point of triangle
* @param _v1  Second point of triangle
* @param _v2  Third point of triangle
* @return     Computed distance
 */
template <class Vec>
typename Vec::value_type
distPointTriangleStable( const Vec& _p,
                         const Vec& _v0,
                         const Vec& _v1,
                         const Vec& _v2,
                         Vec& _nearestPoint )
{ return sqrt(distPointTriangleSquaredStable(_p, _v0, _v1, _v2, _nearestPoint)); }

/** \brief Checks the distance from a point to a plane
*
*
* @param _porigin Planes origin
* @param _pnormal Plane normal ( has to be normalized!)
* @param _point   point to test
* @return         distance
*/
template < typename VectorT , typename ValueT >
inline
ValueT 
distPointPlane(const VectorT& _porigin, 
               const VectorT& _pnormal, 
               const VectorT&  _point);                          

          
/** @} */   
          
//===========================================================================
/** @name Distance Functions ( 3-Dimensional )
* @{ */
//===========================================================================                             
                          
/// squared distance of lines (_v00, _v01) and (_v10, _v11)
template<typename Scalar>
Scalar
distLineLineSquared( const VectorT<Scalar,3>& _v00,
                     const VectorT<Scalar,3>& _v01,
                     const VectorT<Scalar,3>& _v10,
                     const VectorT<Scalar,3>& _v11,
                     VectorT<Scalar,3>*       _min_v0=0,
                     VectorT<Scalar,3>*       _min_v1=0,
                     bool                          _fastApprox=false );


/// distance of lines (_v00, _v01) and (_v10, _v11)
template<typename Scalar>
Scalar
distLineLine( const VectorT<Scalar,3>& _v00,
              const VectorT<Scalar,3>& _v01,
              const VectorT<Scalar,3>& _v10,
              const VectorT<Scalar,3>& _v11,
              VectorT<Scalar,3>*       _min_v0=0,
              VectorT<Scalar,3>*       _min_v1=0 )
{
  return sqrt(distLineLineSquared(_v00, _v01, _v10, _v11,
                                  _min_v0, _min_v1));
}

/** @} */   
   

//===========================================================================
/** @name Projection Functions ( N-Dimensional )
* @{ */
//===========================================================================      


/**
project one point to an edge. If its projection is not on the edge itself, the start or the endpoint is returned
@param _start beginning of edge
@param _end   end of the edge
@param _point point to be projected
*/
template < typename VectorT >
VectorT projectToEdge(const VectorT& _start , 
                      const VectorT& _end , 
                      const VectorT& _point );
                      
                      
/** \brief projects a point to a plane
* @param _porigin Planes origin
* @param _pnormal Plane normal ( has to be normalized! )
* @param _point   point to project
* @return         projected point
*/
template < typename VectorT >
inline
VectorT
projectToPlane(const VectorT& _porigin, 
               const VectorT& _pnormal, 
               const VectorT&  _point);                      
              
/** @} */           

//===========================================================================
/** @name Triangle Functions (2D Only!!)
* @{ */
//===========================================================================   

/** \brief return circumcenter of triangle (_v0,_v1,_v2)
*
*/

/// barycentric coord of _p w.r.t. (_u,_v,_w) in 2D
template<typename Scalar>
bool
baryCoord( const VectorT<Scalar,2>&  _p,
           const VectorT<Scalar,2>&  _u,
           const VectorT<Scalar,2>&  _v,
           const VectorT<Scalar,2>&  _w,
           VectorT<Scalar,3>&        _result );
           

/// is point _p in triangle (_v0,_v1,_v2) in 2D
template<typename Scalar>
bool
isInTriangle( const VectorT<Scalar,2>&  _p,
              const VectorT<Scalar,2>&  _u,
              const VectorT<Scalar,2>&  _v,
              const VectorT<Scalar,2>&  _w )
{
  VectorT<Scalar,3> bary;
  if (baryCoord(_p, _u, _v, _w, bary)) 
    return ((bary[0]>0.0) && (bary[1]>0.0) && (bary[2]>0.0));
  else {
    std::cerr << "point in triangle error\n";
    return false;
  }
}

template<typename Scalar>
bool
circumCenter( const VectorT<Scalar,2>&  _v0,
              const VectorT<Scalar,2>&  _v1,
              const VectorT<Scalar,2>&  _v2,
              VectorT<Scalar,2>&        _result );

/** @} */   

//===========================================================================
/** @name Triangle Functions 3-Dimensional
* @{ */
//===========================================================================  
 
/** barycentric coord of _p w.r.t. (_u,_v,_w) in 3D
    _p has to lie in plane (_u,_v,_w) **/
template<typename Scalar>
bool
baryCoord( const VectorT<Scalar,3>&  _p,
           const VectorT<Scalar,3>&  _u,
           const VectorT<Scalar,3>&  _v,
           const VectorT<Scalar,3>&  _w,
           VectorT<Scalar,3>&        _result );


/// construct min. enclosing sphere
template<typename Scalar>
bool
minSphere(const VectorT<Scalar,3>&  _v0,
          const VectorT<Scalar,3>&  _v1,
          const VectorT<Scalar,3>&  _v2,
          VectorT<Scalar,3>&        _center,
          Scalar&                   _radius);


/// return squared radius of min. enclosing circle of triangle (_v0,_v1,_v2)
template<typename Scalar>
Scalar
minRadiusSquared( const VectorT<Scalar,3>&  _v0,
                  const VectorT<Scalar,3>&  _v1,
                  const VectorT<Scalar,3>&  _v2 );

  
/// return radius of min. enclosing circle of triangle (_v0,_v1,_v2)
template<typename Scalar>
Scalar
minRadius( const VectorT<Scalar,3>&  _v0,
           const VectorT<Scalar,3>&  _v1,
           const VectorT<Scalar,3>&  _v2 )
{
  return sqrt(minRadiusSquared(_v0, _v1, _v2));
}


/// return circumcenter of triangle (_v0,_v1,_v2)
template<typename Scalar>
bool
circumCenter( const VectorT<Scalar,3>&  _v0,
              const VectorT<Scalar,3>&  _v1,
              const VectorT<Scalar,3>&  _v2,
              VectorT<Scalar,3>&        _result );


/// return squared radius of circumcircle of triangle (_v0,_v1,_v2)
template<typename Scalar>
Scalar
circumRadiusSquared( const VectorT<Scalar,3>&  _v0,
                     const VectorT<Scalar,3>&  _v1,
                     const VectorT<Scalar,3>&  _v2 );


/// return radius of circumcircle of triangle (_v0,_v1,_v2)
template<typename Scalar>
Scalar
circumRadius( const VectorT<Scalar,3>&  _v0,
              const VectorT<Scalar,3>&  _v1,
              const VectorT<Scalar,3>&  _v2 )
{
  return sqrt(circumRadiusSquared(_v0, _v1, _v2));
}

/**
* test angles in triangle
* return 0 if all smaller than 90?
* return 1 if angle at _p0 ist greater than 90?
* return 2 if angle at _p1 ist greater than 90?
* return 3 if angle at _p2 ist greater than 90?
*/
template<class VectorT>
int isObtuse(const VectorT& _p0,
             const VectorT& _p1,
             const VectorT& _p2 );

/** @} */   

//===========================================================================
/** @name Triangle Functions N-Dimensional
* @{ */
//===========================================================================   


/** \brief return squared area of triangle (_v0, _v1, _v2)
*
* @param _v0  First point of triangle
* @param _v1  Second point of triangle
* @param _v2  Third point of triangl
*/
template <class Vec>
typename Vec::value_type
triangleAreaSquared( const Vec& _v0,
                     const Vec& _v1,
                     const Vec& _v2 );


/** \brief return area of triangle (_v0, _v1, _v2)
*
* @param _v0  First point of triangle
* @param _v1  Second point of triangle
* @param _v2  Third point of triangl
*/
template <class Vec>
typename Vec::value_type
triangleArea( const Vec& _v0,
              const Vec& _v1,
              const Vec& _v2 )
{
  return sqrt(triangleAreaSquared(_v0,_v1,_v2));
}
  

/** \brief return aspect ratio (length/height) of triangle
*
* @param _v0  First point of triangle
* @param _v1  Second point of triangle
* @param _v2  Third point of triangl
*/
template <typename Scalar, int N>
Scalar
aspectRatio( const VectorT<Scalar, N>& _v0,
             const VectorT<Scalar, N>& _v1,
             const VectorT<Scalar, N>& _v2 );

/** \brief return roundness of triangle: 1=equilateral, 0=colinear
*
* @param _v0  First point of triangle
* @param _v1  Second point of triangle
* @param _v2  Third point of triangl
*/
template <typename Scalar, int N>
Scalar
roundness( const VectorT<Scalar, N>& _v0,
           const VectorT<Scalar, N>& _v1,
           const VectorT<Scalar, N>& _v2 );

/** @} */

template<typename Vector>
Vector closestPointLineSegment(Vector x, Vector p1, Vector p2) {
    const auto delta = ((p2-p1)|(x-p1)) / (p2-p1).sqrnorm();
    //std::cout << "\x1b[32mdelta = " << delta << "\x1b[0m" << std::endl;
    if (delta <= 0) {
        //std::cout << "\x1b[43mdelta <= 0\x1b[0m" << std::endl;
        return p1;
    } else if (delta >= 1) {
        //std::cout << "\x1b[43mdelta >= 1\x1b[0m" << std::endl;
        return p2;
    } else if (delta != delta) { // p1 = p2 or x = p1
        //std::cout << "\x1b[43mdelta != delta\x1b[0m" << std::endl;
        return (x-p1).sqrnorm() < (x-p2).sqrnorm() ? p1 : p2;
    } else {
        //std::cout << "\x1b[43mdelta \\in [0, 1]\x1b[0m" << std::endl;
        return (1 - delta) * p1 + delta * p2;
    }
};


//=============================================================================
} // namespace Geometry
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(GEO_ALGORITHMS_C)
#define GEO_ALGORITHMS_TEMPLATES
#include "Algorithms.cc"
#endif
//=============================================================================
#endif // GEO_ALGORITHMS_HH defined
//=============================================================================


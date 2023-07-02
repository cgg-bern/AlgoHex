/*===========================================================================*\
*                                                                            *
*                              OpenFlipper                                   *
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
*                                                                            *
\*===========================================================================*/


//=============================================================================
//
//  CLASS BSPImplT
//
//=============================================================================

#ifndef BSPIMPLT_HH
#define BSPIMPLT_HH


//== INCLUDES =================================================================


#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include "Algorithms.hh"


//== CLASS DEFINITION =========================================================
#include <vector>
#include <ostream>

template <class BSPCore>
class BSPImplT : public BSPCore
{
public: //---------------------------------------------------------------------

  typedef typename BSPCore::Traits      Traits;
  typedef typename BSPCore::Handle      Handle;
  typedef typename BSPCore::Point       Point;
  typedef typename BSPCore::Scalar      Scalar;
  typedef typename BSPCore::Node        Node;
  typedef typename BSPCore::Handles     Handles;
  typedef typename BSPCore::HandleIter  HandleIter;


public: //---------------------------------------------------------------------

  BSPImplT(const Traits& _traits, const Scalar& _infinity = std::numeric_limits<Scalar>::infinity()) :
      BSPCore(_traits),
      infinity_(_infinity) {}
  ~BSPImplT() {}


  /// Store nearest neighbor information
  struct NearestNeighbor
  {
    NearestNeighbor() {}
    NearestNeighbor(Handle _h, Scalar _d) : handle(_h), dist(_d) {}
    Handle  handle;
    Scalar  dist;
  };

  /// Store nearest neighbor information
  typedef  std::vector< std::pair<Handle,Scalar> > RayCollision;

  /// Return handle of the nearest neighbor face
  NearestNeighbor nearest(const Point& _p) const;

  /** \brief intersect mesh with ray
   *
   * This function shots a ray through the mesh and collects all intersected triangles and
   * the handle of the closest face ( non-directional, so no matter of the ray direction, the
   * closest face handle is returned in either direction)
   *
   * @param _p Start point of the ray
   * @param _r Ray direction
   * @return   Collision information
   */
  RayCollision raycollision (const Point& _p, const Point& _r) const;

  /** \brief intersect mesh with ray
   *
   * This function shots a ray through the mesh and collects all intersected triangles and
   * the handle of the closest face ( directional, so the ray direction is taken into account!).
   *
   * Only hits with a distance > 0.0 to the point p will be collected (_p will be skipped!)
   *
   * @param _p Start point of the ray
   * @param _r Ray direction
   * @return   Collision information
   */
  RayCollision directionalRaycollision (const Point& _p, const Point& _r) const;

  /** \brief intersect mesh with ray
   *
   * This function shots a ray through the mesh and determines the first intersected triangle and
   * the handle of the closest face ( directional, so the ray direction is taken into account!).
   *
   * Only hits with a distance > 0.0 to the point p will be collected (_p will be skipped!).
   * Note that for compatibility reasons the return type is still a vector of collisions.
   *
   * @param _p Start point of the ray
   * @param _r Ray direction
   * @return   Collision information
   */
  RayCollision nearestRaycollision(const Point& _p, const Point& _r) const;

  /** \brief intersect mesh with open ball
   *
   * All triangles that have at least one vertex (!) inside the ball are given to the Callback,
   * triangles which have no vertex inside the ball but intersect it MAY be returned. (TODO)
   * Each triangle can be returned up to three times, make sure to handle this, e.g.
   * by putting the values into an std::(unordered_)set.
   *
   * @param _c Center of the ball
   * @param _r Radius of the ball
   * @param _callback Callable that accepts Handle or const Handle&, e.g. (const Handle &h) -> void
   */
  template<class Callback>
  void intersectBall(const Point & _c, Scalar _r, Callback _callback) const;

private: //---------------------------------------------------------------------


  /// Store nearest neighbor information
  struct NearestNeighborData
  {
    Point   ref;
    Scalar  dist;
    Handle  nearest;

    friend std::ostream &operator<< (std::ostream &stream, const NearestNeighborData &data) {
        stream << "[NearestNeghborData instance. ref: " << data.ref << ", dist: " << data.dist << ", nearest: " << data.nearest.idx() << "]";
        return stream;
    }
  };

  /// Store ray collide information
  struct RayCollisionData
  {
    RayCollisionData() : ref(), ray() {}

    Point   ref;
    Point   ray;
    RayCollision hit_handles;
  };


  // Recursive part of nearest()
  void _nearest(Node* _node, NearestNeighborData& _data) const;

  /**  \brief recursive part of raycollision()
   *
   * @param _node The current node in the tree
   * @param _data Data pointer, used to collect the collision information
   */
  void _raycollision_non_directional(Node* _node, RayCollisionData& _data) const;

  /**  \brief recursive part of directionalRaycollision()
   *
   * @param _node The current node in the tree
   * @param _data Data pointer, used to collect the collision information
   */
  void _raycollision_directional(Node* _node, RayCollisionData& _data) const;

  void _raycollision_nearest_directional(Node* _node, RayCollisionData& _data) const;

  template<class Callback>
  void _intersect_ball(const Node& _node, const Point & _c, Scalar _r, Callback _callback) const;

  template<typename T,typename U>
  struct less_pair_second {
      bool operator()(const std::pair<T,U> &left, const std::pair<T,U> &right) {
          return (fabs(left.second) < fabs(right.second));
      }
  };

  const Scalar infinity_;

};

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(BSPIMPLT_C)
#  define BSPIMPLT_TEMPLATES
#  include "BSPImplT_impl.hh"
#endif
//=============================================================================
#endif // BSPIMPLT_HH defined
//=============================================================================

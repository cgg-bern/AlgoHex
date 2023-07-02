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

#define BSPIMPLT_C

//== INCLUDES =================================================================


#include "BSPImplT.hh"
#include <cfloat>
#include <cmath>

//== CLASS DEFINITION =========================================================
#include <vector>
#include <stdexcept>
#include <limits>

template <class BSPCore>
typename BSPImplT<BSPCore>::NearestNeighbor
BSPImplT<BSPCore>::
nearest(const Point& _p) const
{
  NearestNeighborData  data;
  data.ref  = _p;
  data.dist = infinity_;
  if (this->root_ == 0)
    throw std::runtime_error("It seems like the BSP hasn't been built, yet. Did you call build(...)?");
  _nearest(this->root_, data);
  return NearestNeighbor(data.nearest, sqrt(data.dist));
}


//-----------------------------------------------------------------------------


template <class BSPCore>
void
BSPImplT<BSPCore>::
_nearest(Node* _node, NearestNeighborData& _data) const
{
  // terminal node
  if (!_node->left_child_)
  {
    Scalar dist(0);
    for (HandleIter it=_node->begin(); it!=_node->end(); ++it)
    {
      dist = this->traits_.sqrdist(*it, _data.ref);
      if (dist < _data.dist)
      {
        _data.dist = dist;
        _data.nearest = *it;
      }
    }
  }

  // non-terminal node
  else
  {
    Scalar dist = _node->plane_.distance(_data.ref);
    if (dist > 0.0)
    {
      _nearest(_node->left_child_, _data);
      if (dist*dist < _data.dist)
        _nearest(_node->right_child_, _data);
    }
    else
    {
      _nearest(_node->right_child_, _data);
      if (dist*dist < _data.dist)
        _nearest(_node->left_child_, _data);
    }
  }
}

//-----------------------------------------------------------------------------

template <class BSPCore>
typename BSPImplT<BSPCore>::RayCollision
BSPImplT<BSPCore>::
raycollision(const Point& _p, const Point& _r) const
{
  // Prepare the struct for returning the data
  RayCollisionData  data;
  data.ref  = _p;
  data.ray  = _r;
  data.hit_handles.clear();

  _raycollision_non_directional(this->root_, data);

  std::sort(data.hit_handles.begin(), data.hit_handles.end(), less_pair_second<Handle,Scalar>());
  return RayCollision(data.hit_handles);
}

template <class BSPCore>
typename BSPImplT<BSPCore>::RayCollision
BSPImplT<BSPCore>::
directionalRaycollision(const Point& _p, const Point& _r) const {

  // Prepare the struct for returning the data
  RayCollisionData  data;
  data.ref  = _p;
  data.ray  = _r;
  data.hit_handles.clear();

  _raycollision_directional(this->root_, data);

  std::sort(data.hit_handles.begin(), data.hit_handles.end(), less_pair_second<Handle,Scalar>());
  return RayCollision(data.hit_handles);

}

template <class BSPCore>
typename BSPImplT<BSPCore>::RayCollision
BSPImplT<BSPCore>::
nearestRaycollision(const Point& _p, const Point& _r) const {

  // Prepare the struct for returning the data
  RayCollisionData  data;
  data.ref  = _p;
  data.ray  = _r;
  data.hit_handles.clear();

  _raycollision_nearest_directional(this->root_, data);

  return RayCollision(data.hit_handles);
}

template<class BSPCore>
template<class Callback>
void
BSPImplT<BSPCore>::
intersectBall(const Point &_c, Scalar _r, Callback _callback) const
{
    _intersect_ball(*(this->root_), _c, _r, _callback);
}


//-----------------------------------------------------------------------------


template <class BSPCore>
void
BSPImplT<BSPCore>::
_raycollision_non_directional(Node* _node, RayCollisionData& _data) const
{
  // terminal node
  if (!_node->left_child_)
  {
    Scalar dist;
    Point v0, v1, v2;
    Scalar u, v;

    for (HandleIter it=_node->begin(); it!=_node->end(); ++it)
    {
      this->traits_.points(*it, v0, v1, v2);
      if (ACG::Geometry::triangleIntersection(_data.ref, _data.ray, v0, v1, v2, dist, u, v)) {
        // face intersects with ray.
        _data.hit_handles.push_back(std::pair<Handle,Scalar>(*it,dist));
      }
    }
  }

  // non-terminal node
  else
  {
    Scalar tmin, tmax;
    if ( _node->left_child_ && ACG::Geometry::axisAlignedBBIntersection( _data.ref, _data.ray, _node->left_child_->bb_min, _node->left_child_->bb_max, tmin, tmax)) {
      _raycollision_non_directional(_node->left_child_, _data);
    }
    if ( _node->right_child_ && ACG::Geometry::axisAlignedBBIntersection( _data.ref, _data.ray, _node->right_child_->bb_min, _node->right_child_->bb_max, tmin, tmax)) {
      _raycollision_non_directional(_node->right_child_, _data);
    }
  }
}

//-----------------------------------------------------------------------------


template <class BSPCore>
void
BSPImplT<BSPCore>::
_raycollision_directional(Node* _node, RayCollisionData& _data) const
{
  // terminal node
  if (!_node->left_child_)
  {
    Scalar dist;
    Point v0, v1, v2;
    Scalar u, v;

    for (HandleIter it=_node->begin(); it!=_node->end(); ++it)
    {
      this->traits_.points(*it, v0, v1, v2);
      if (ACG::Geometry::triangleIntersection(_data.ref, _data.ray, v0, v1, v2, dist, u, v)) {
        if (dist > 0.0){
          _data.hit_handles.push_back(std::pair<Handle,Scalar>(*it,dist));
        }
      }
    }
  }

  // non-terminal node
  else
  {
    Scalar tmin, tmax;
    if ( _node->left_child_ && ACG::Geometry::axisAlignedBBIntersection( _data.ref, _data.ray, _node->left_child_->bb_min, _node->left_child_->bb_max, tmin, tmax)) {
      _raycollision_directional(_node->left_child_, _data);
    }
    if ( _node->right_child_ && ACG::Geometry::axisAlignedBBIntersection( _data.ref, _data.ray, _node->right_child_->bb_min, _node->right_child_->bb_max, tmin, tmax)) {
      _raycollision_directional(_node->right_child_, _data);
    }
  }
}

//-----------------------------------------------------------------------------


template <class BSPCore>
void
BSPImplT<BSPCore>::
_raycollision_nearest_directional(Node* _node, RayCollisionData& _data) const
{
  // terminal node
  if (!_node->left_child_)
  {
    Scalar dist;
    Point v0, v1, v2;
    Scalar u, v;

    for (HandleIter it=_node->begin(); it!=_node->end(); ++it)
    {
      this->traits_.points(*it, v0, v1, v2);
      if (ACG::Geometry::triangleIntersection(_data.ref, _data.ray, v0, v1, v2, dist, u, v)) {
        if (dist > 0.0){
          _data.hit_handles.push_back(std::pair<Handle,Scalar>(*it, dist));
        }
      }
    }
		// only return the closest hit
		if(!_data.hit_handles.empty()) {
			std::partial_sort(_data.hit_handles.begin(), _data.hit_handles.begin() + 1, _data.hit_handles.end(), less_pair_second<Handle, Scalar>());
			_data.hit_handles.resize(1);
		}
  }

  // non-terminal node
  else
	{
		// determine order of traversal
		Node* first_node = _node->left_child_, *second_node = _node->right_child_;
		if (!_node->plane_(_data.ref)) {
			std::swap(first_node, second_node);
		}

		Scalar tmin, tmax;
		if ( first_node && ACG::Geometry::axisAlignedBBIntersection( _data.ref, _data.ray, first_node->bb_min, first_node->bb_max, tmin, tmax) ) {
			_raycollision_nearest_directional(first_node, _data);
		}
		// if the second node is further away than the closeset hit skip it
		Scalar dist = std::numeric_limits<Scalar>::max();
		if(!_data.hit_handles.empty()) {
			dist = _data.hit_handles.front().second;
		}
		if ( second_node && ACG::Geometry::axisAlignedBBIntersection( _data.ref, _data.ray, second_node->bb_min, second_node->bb_max, tmin, tmax) && (tmin < dist) ) {
			_raycollision_nearest_directional(second_node, _data);
		}
  }
}
template<class BSPCore>
template<class Callback>
void BSPImplT<BSPCore>::
_intersect_ball(const Node &_node,
                const Point &_c,
                Scalar _r,
                Callback _callback) const
{
    // terminal node
    if (!_node.left_child_)
    {
        const double r_sqr = _r * _r;
        for (const auto &fh: _node) {
            const double dist = this->traits_.sqrdist(fh, _c);
            if (dist < r_sqr) {
                _callback(fh);
            }
        }
    }
    else // non-terminal node
    {
        const Scalar dist = _node.plane_.distance(_c);
        const Node &left = *_node.left_child_;
        const Node &right = *_node.right_child_;

        if (dist > -_r){
            _intersect_ball(left, _c, _r, _callback);
        }
        if (dist < _r) {
            _intersect_ball(right, _c, _r, _callback);
        }
    }
}


//=============================================================================


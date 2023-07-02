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
//  CLASS TriangleBSPCoreT
//
//=============================================================================

#define TRIANGLEBSPCORET_C

//== INCLUDES =================================================================


#include "TriangleBSPCoreT.hh"


//== CLASS DEFINITION =========================================================

template <class BSPTraits>
void
TriangleBSPCoreT<BSPTraits>::
build(unsigned int _max_handles, unsigned int _max_depth)
{
  // init
  delete root_;
  root_ = new Node(handles_, 0);

  // delete own handles (don't store them twice)
  handles_ = Handles();

  nodes=1;
  traits_.calculateBoundingBoxRoot (root_);
  // call recursive helper
  _build(root_, _max_handles, _max_depth);

}


//-----------------------------------------------------------------------------


template <class BSPTraits>
void
TriangleBSPCoreT<BSPTraits>::
_build(Node*         _node,
       unsigned int  _max_handles,
       unsigned int  _depth)
{
  // should we stop at this level ?
  if ((_depth == 0) || ((_node->end()-_node->begin()) <= (int)_max_handles))
    return;

  Point median;
  int axis;
  // compute bounding boxes for children
  traits_.calculateBoundingBox (_node, median, axis);

  // construct splitting plane
  const Point XYZ[3] = { Point(1,0,0), Point(0,1,0), Point(0,0,1) };
  _node->plane_ = Plane(median, XYZ[axis]);

  // partition for left and right child
  Handles lhandles, rhandles;
  lhandles.reserve(_node->handles_.size()/2);
  rhandles.reserve(_node->handles_.size()/2);

  HandleIter it;
  Point p0, p1, p2;
  for (it=_node->begin(); it!=_node->end(); ++it)
  {
    traits_.points(*it, p0, p1, p2);

    /* @TODO Remove this comment block. Replaced by the clause below

    if (_node->plane_(p0))  left  = true;
    else                    right = true;
    if (_node->plane_(p1))  left  = true;
    else                    right = true;
    if (_node->plane_(p2))  left  = true;
    else                    right = true;

    if (left)  lhandles.push_back(*it);
    if (right) rhandles.push_back(*it);
    */

    if ( _node->plane_(p0)  || _node->plane_(p1)  || _node->plane_(p2)  ) lhandles.push_back(*it);
    if ( !_node->plane_(p0) || !_node->plane_(p1) || !_node->plane_(p2) ) rhandles.push_back(*it);

  }

  // check it
  if (lhandles.size() == _node->handles_.size() ||
      rhandles.size() == _node->handles_.size())
    return;
  else
    _node->handles_ = Handles();


  // create children
  _node->left_child_  = new Node(lhandles, _node);  lhandles = Handles();
  _node->right_child_ = new Node(rhandles, _node);  rhandles = Handles();
  nodes+=2;

  //save bounding boxes for children
  /*
  _node->left_child_->bb_min  = _node->bb_min;
  _node->left_child_->bb_max  = _node->bb_max;
  _node->left_child_->bb_max[axis] = median [axis];

  _node->right_child_->bb_min = _node->bb_min;
  _node->right_child_->bb_min[axis] = median [axis];
  _node->right_child_->bb_max = _node->bb_max;
  */
  _node->right_child_->bb_min  = _node->bb_min;
  _node->right_child_->bb_max  = _node->bb_max;
  _node->right_child_->bb_max[axis] = median [axis];

  _node->left_child_->bb_min = _node->bb_min;
  _node->left_child_->bb_min[axis] = median [axis];
  _node->left_child_->bb_max = _node->bb_max;

  // recurse to childen
  _build(_node->left_child_,  _max_handles, _depth-1);
  _build(_node->right_child_, _max_handles, _depth-1);
}

//=============================================================================



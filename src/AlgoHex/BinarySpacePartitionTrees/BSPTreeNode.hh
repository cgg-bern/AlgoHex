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
//  CLASS TreeNode
//
//=============================================================================

#ifndef MB_BSPTREENODE_HH
#define MB_BSPTREENODE_HH

//== INCLUDES =================================================================

#include "PlaneT.hh"
//#include <ACG/Geometry/Algorithms.hh>
#include <ostream>

//== CLASS DEFINITION =========================================================

// Node of the tree: contains parent, children and splitting plane
template <class BSPTraits>
struct TreeNode
{
    typedef typename BSPTraits::Handle       Handle;
    typedef typename BSPTraits::Point        Point;
    typedef typename BSPTraits::VertexHandle VertexHandle;
    typedef std::vector<Handle>              Handles;
    typedef typename Handles::iterator       HandleIter;
    typedef typename Handles::const_iterator HandleConstIter;
    typedef typename Point::value_type       Scalar;
    typedef ACG::Geometry::PlaneT<Scalar>    Plane;

    TreeNode(const Handles& _handles, TreeNode* _parent)
            : handles_(_handles),
            parent_(_parent), left_child_(0), right_child_(0) {}
    ~TreeNode()
    {
        delete left_child_;
        delete right_child_;

        if (parent_)
        {
            if (this == parent_->left_child_)
                parent_->left_child_ = 0;
            else
                parent_->right_child_ = 0;
        }
    }

    HandleIter begin() {
        return handles_.begin();
    }

    HandleIter end()   {
        return handles_.end();
    }

    HandleConstIter begin() const {
        return handles_.begin();
    }

    HandleConstIter end() const {
        return handles_.end();
    }

    size_t size() const {
        return handles_.size();
    }

    Handles     handles_;
    TreeNode    *parent_, *left_child_, *right_child_;
    Plane       plane_;
    Point	bb_min, bb_max;
    
    /// This visualizes the bounding boxes
    template< typename MeshT >
    void visualizeTree(MeshT *_object, int _max_depth)
    {
        if (_max_depth > 0 && (left_child_ || right_child_) )
        {
            if (left_child_)
                left_child_->visualizeTree(_object, _max_depth-1);
            if (right_child_)
                right_child_->visualizeTree(_object, _max_depth-1);
        }
        else
        {
            Point size_ = bb_max - bb_min;

            std::vector<VertexHandle> vhandle(8);
            vhandle[0] = _object->add_vertex(bb_min+Point(0.0,0.0,size_[2]));
            vhandle[1] = _object->add_vertex(bb_min+Point(size_[0],0.0,size_[2]));
            vhandle[2] = _object->add_vertex(bb_min+Point(size_[0],size_[1],size_[2]));
            vhandle[3] = _object->add_vertex(bb_min+Point(0.0,size_[1],size_[2]));
            vhandle[4] = _object->add_vertex(bb_min+Point(0.0,0.0,0.0));
            vhandle[5] = _object->add_vertex(bb_min+Point(size_[0],0.0,0.0));
            vhandle[6] = _object->add_vertex(bb_min+Point(size_[0],size_[1],0.0));
            vhandle[7] = _object->add_vertex(bb_min+Point(0.0,size_[1],0.0));


            // generate (quadrilateral) faces
            std::vector<VertexHandle>  face_vhandles;

            face_vhandles.clear();
            face_vhandles.push_back(vhandle[0]);
            face_vhandles.push_back(vhandle[1]);
            face_vhandles.push_back(vhandle[2]);
            face_vhandles.push_back(vhandle[3]);
            _object->add_face(face_vhandles);

            face_vhandles.clear();
            face_vhandles.push_back(vhandle[7]);
            face_vhandles.push_back(vhandle[6]);
            face_vhandles.push_back(vhandle[5]);
            face_vhandles.push_back(vhandle[4]);
            _object->add_face(face_vhandles);

            face_vhandles.clear();
            face_vhandles.push_back(vhandle[1]);
            face_vhandles.push_back(vhandle[0]);
            face_vhandles.push_back(vhandle[4]);
            face_vhandles.push_back(vhandle[5]);
            _object->add_face(face_vhandles);

            face_vhandles.clear();
            face_vhandles.push_back(vhandle[2]);
            face_vhandles.push_back(vhandle[1]);
            face_vhandles.push_back(vhandle[5]);
            face_vhandles.push_back(vhandle[6]);
            _object->add_face(face_vhandles);

            face_vhandles.clear();
            face_vhandles.push_back(vhandle[3]);
            face_vhandles.push_back(vhandle[2]);
            face_vhandles.push_back(vhandle[6]);
            face_vhandles.push_back(vhandle[7]);
            _object->add_face(face_vhandles);

            face_vhandles.clear();
            face_vhandles.push_back(vhandle[0]);
            face_vhandles.push_back(vhandle[3]);
            face_vhandles.push_back(vhandle[7]);
            face_vhandles.push_back(vhandle[4]);
            _object->add_face(face_vhandles);
        }
    }

    private:
    /*
     * Noncopyable because of root_.
     */
    TreeNode(const TreeNode &rhs);
    TreeNode &operator=(const TreeNode &rhs);

};

template<class BSPTraits>
std::ostream &operator<< (std::ostream &stream, const TreeNode<BSPTraits> &node) {
    stream << "[TreeNode instance. Handles: ";
    for (typename TreeNode<BSPTraits>::Handles::const_iterator it = node.handles_.begin(), it_end = node.handles_.end();
            it != it_end; ++it) {
        stream << it->idx();
        if (it < it_end-1) stream << ", ";
    }
    stream << ", parent: " << node.parent_ << ", left_child_: " << node.left_child_
            << ", right_child_: " << node.right_child_ << ", plane_: <not implemented>, bb_min: "
            << node.bb_min << ", bb_max: " << node.bb_max << ", size(): " << node.size() << "]";
    return stream;
}

//=============================================================================
#endif // MB_BSPTREENODE_HH defined
//=============================================================================

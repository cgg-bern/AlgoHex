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

#ifndef TRIANGLEBSPCORET_HH
#define TRIANGLEBSPCORET_HH


//== INCLUDES =================================================================


#include <vector>
#include "PlaneT.hh"
#include <OpenVolumeMesh/Geometry/VectorT.hh>

#include "TriangleBSPT.hh"


//== CLASS DEFINITION =========================================================


template <class BSPTraits>
class TriangleBSPCoreT
{
public: //---------------------------------------------------------------------

  typedef BSPTraits                      Traits;
  typedef typename BSPTraits::Point      Point;
  typedef typename BSPTraits::Handle     Handle;
  typedef typename BSPTraits::Node	     Node;
  typedef typename Point::value_type     Scalar;
  typedef ACG::Geometry::PlaneT<Scalar>  Plane;
  typedef std::vector<Handle>            Handles;
  typedef typename Handles::iterator     HandleIter;


public: //---------------------------------------------------------------------


  /** Constructor: need traits that define the types and
      give us the points by traits_.point(PointHandle) */
  explicit TriangleBSPCoreT(const BSPTraits& _traits) : traits_(_traits), root_(0), nodes(0), n_triangles(0) {}

  /// Destructor
  ~TriangleBSPCoreT() {
      delete root_;
  }


  /// Reserve memory for _n entries
  void reserve(size_t _n) { handles_.reserve(_n); }
  /// Add a handle to the BSP
  void push_back(Handle _h)     { handles_.push_back(_h); ++n_triangles; }

  /**
   * @return size() == 0
   */
  bool empty()     { return n_triangles == 0; }

  /**
   * @return The number of triangles contained in the tree.
   */
  size_t size()     { return n_triangles; }

  /** Build the tree.
   *
   * @param _max_handles Maximum number of vertices per leaf.
   * @param _max_depth   Maximum depth.
   */
  void build(unsigned int _max_handles, unsigned int _max_depth);

    /** \brief Create a PolyMesh object that visualizes the bounding boxes of the BSP tree
     *
     * @param _object     The output mesh which the tree will be written into
     * @param _max_depth  The maximal depth that will be visualized
     */
  template <typename MeshT>
  void visualizeTree(MeshT *_object, int _max_depth)
  {
    root_->visualizeTree(_object, _max_depth-1);
    _object->update_normals();
  }

private:
  /*
   * Noncopyable because of root_.
   */
  TriangleBSPCoreT(const TriangleBSPCoreT &rhs);
  TriangleBSPCoreT &operator=(const TriangleBSPCoreT &rhs);


private: //---------------------------------------------------------------------


  // Recursive part of build()
  void _build(Node*        _node,
	      unsigned int _max_handles,
	      unsigned int _depth);




protected: //-------------------------------------------------------------------


  BSPTraits  traits_;
  Handles    handles_;
  Node*      root_;
  int	       nodes, n_triangles;

};


//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(TRIANGLEBSPCORET_C)
#  define TRIANGLEBSPCORET_TEMPLATES
#  include "TriangleBSPCoreT_impl.hh"
#endif
//=============================================================================
#endif // TRIANGLEBSPCORET_HH defined
//=============================================================================

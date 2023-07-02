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
//  CLASS TriangleBSPT
//
//=============================================================================

#pragma once


//== INCLUDES =================================================================

#include "BSPTreeNode.hh"
#include "TriangleBSPCoreT.hh"
#include "BSPImplT.hh"

#include <list>
//== CLASS DEFINITION =========================================================

namespace AlgoHex {
    template<class BSPTraits>
    class TriangleBSPT : public BSPImplT<TriangleBSPCoreT<BSPTraits> > {
    public:
        typedef BSPImplT<TriangleBSPCoreT<BSPTraits> > Base;
        typedef typename Base::Scalar Scalar;

        TriangleBSPT(const BSPTraits &_traits,
                     const Scalar &_infinity = std::numeric_limits<Scalar>::infinity()) : Base(_traits, _infinity) {}
    };

//== CLASS DEFINITION =========================================================

    template<class Mesh, class SpecificTraits>
    class OVMOMCommonTriangleBSPTraits : public SpecificTraits {
    public:

        typedef typename SpecificTraits::Point Point;
        typedef typename SpecificTraits::Handle Handle;
        typedef typename Point::value_type Scalar;
        typedef std::vector<Handle> Handles;
        typedef typename Handles::iterator HandleIter;
        typedef TreeNode<SpecificTraits> Node;

        explicit OVMOMCommonTriangleBSPTraits(const Mesh &_mesh) : SpecificTraits(_mesh) {}

        Scalar sqrdist(const Handle _h, const Point &_p) const {
            Point p0, p1, p2, q;
            this->points(_h, p0, p1, p2);
            return ACG::Geometry::distPointTriangleSquaredStable(_p, p0, p1, p2, q);
        }

        void calculateBoundingBox(Node *_node, Point &median, int &axis) {
            //determine splitting axis
            HandleIter it_h;
            Point p0, p1, p2, bb_min, bb_max;
            bb_min.vectorize(std::numeric_limits<typename Point::value_type>::infinity());
            bb_max.vectorize(-std::numeric_limits<typename Point::value_type>::infinity());
            std::list<Point> vertices;

            for (it_h = _node->begin(); it_h != _node->end(); ++it_h) {
                this->points(*it_h, p0, p1, p2);
                /*
                 bb_min.minimize(p0);
                 bb_min.minimize(p1);
                 bb_min.minimize(p2);
                 bb_max.maximize(p0);
                 bb_max.maximize(p1);
                 bb_max.maximize(p2);*/

                vertices.push_back(p0);
                vertices.push_back(p1);
                vertices.push_back(p2);
            }
            bb_min = _node->bb_min;
            bb_max = _node->bb_max;

            // split longest side of bounding box
            Point bb = bb_max - bb_min;
            Scalar length = bb[0];
            axis = 0;
            if (bb[1] > length)
                length = bb[(axis = 1)];
            if (bb[2] > length)
                length = bb[(axis = 2)];

            //calculate the median value in axis-direction
            switch (axis) {
                case 0:
                    vertices.sort(x_sort());
                    break;
                case 1:
                    vertices.sort(y_sort());
                    break;
                case 2:
                    vertices.sort(z_sort());
                    break;
            }
            vertices.unique(); ///todo: does this work with Points?!

            size_t size = vertices.size();
            typename std::list<Point>::iterator it_v;
            it_v = vertices.begin();
            std::advance(it_v, size / 2);
            median = *it_v;

        }

        void calculateBoundingBoxRoot(Node *_node) {
            HandleIter it;
            Point p0, p1, p2, bb_min, bb_max;
            bb_min.vectorize(FLT_MAX);
            bb_max.vectorize(-FLT_MAX);
            for (it = _node->begin(); it != _node->end(); ++it) {
                this->points(*it, p0, p1, p2);
                bb_min.minimize(p0);
                bb_min.minimize(p1);
                bb_min.minimize(p2);
                bb_max.maximize(p0);
                bb_max.maximize(p1);
                bb_max.maximize(p2);
            }
            _node->bb_min = bb_min;
            _node->bb_max = bb_max;
        }

    private:
        //functors for sorting in different directions
        struct x_sort {
            bool operator()(const Point &first, const Point &second) { return (first[0] < second[0]); }
        };

        struct y_sort {
            bool operator()(const Point &first, const Point &second) { return (first[1] < second[1]); }
        };

        struct z_sort {
            bool operator()(const Point &first, const Point &second) { return (first[2] < second[2]); }
        };
    };


//== CLASS DEFINITION =========================================================

    template<class Mesh>
    class OMSpecificTriangleBSPTraits {
    public:
        typedef typename Mesh::Point Point;
        typedef typename Mesh::FaceHandle Handle;
        typedef typename Mesh::VertexHandle VertexHandle;

        explicit OMSpecificTriangleBSPTraits(const Mesh &_mesh) : mesh_(_mesh) {}

        /// Returns the points belonging to the face handle _h
        inline void points(const Handle &_h, Point &_p0, Point &_p1, Point &_p2) const {
            const auto &mesh = this->mesh_;
            typename Mesh::CFVIter fv_it(mesh.cfv_iter(_h));
            _p0 = mesh.point(*fv_it);
            ++fv_it;
            _p1 = mesh.point(*fv_it);
            ++fv_it;
            _p2 = mesh.point(*fv_it);
        }

    protected:
        const Mesh &mesh_;
    };

    template<class Mesh>
    using OpenMeshTriangleBSPTraits = OVMOMCommonTriangleBSPTraits<Mesh, OMSpecificTriangleBSPTraits<Mesh>>;


//== CLASS DEFINITION =========================================================
    template<class Mesh>
    class OpenMeshTriangleBSPT
            : public TriangleBSPT<OpenMeshTriangleBSPTraits<Mesh> > {
    public:
        typedef OpenMeshTriangleBSPTraits<Mesh> Traits;
        typedef TriangleBSPT<Traits> Base;
        typedef typename Traits::Scalar Scalar;

        OpenMeshTriangleBSPT(const Mesh &_mesh,
                             const Scalar &_infinity = std::numeric_limits<Scalar>::infinity())
                : Base(Traits(_mesh), _infinity) {}
    };


//#if (defined ENABLE_POLYHEDRALMESH_SUPPORT) \
//    || (defined ENABLE_HEXAHEDRALMESH_SUPPORT) \
//    || (defined ENABLE_TETRAHEDRALMESH_SUPPORT)

#include <OpenVolumeMesh/Core/Handles.hh>

//== CLASS DEFINITION =========================================================

// Make sure all faces of the mesh have valence 3.
    template<class Mesh>
    class OVMSpecificTriangleBSPTraits {
    public:
        typedef typename Mesh::PointT Point;
        typedef OpenVolumeMesh::FaceHandle Handle;
        typedef OpenVolumeMesh::VertexHandle VertexHandle;

        explicit OVMSpecificTriangleBSPTraits(const Mesh &_mesh) : mesh_(_mesh) {}

        /// Returns the points belonging to the face handle _h
        inline void points(const Handle &_h, Point &_p0, Point &_p1, Point &_p2) const {
            const auto &mesh = this->mesh_;
            assert(mesh.valence(_h) == 3);
            auto hfv_it(mesh.hfv_iter(mesh.halfface_handle(_h, 0)));
            _p0 = mesh.vertex(*hfv_it++);
            _p1 = mesh.vertex(*hfv_it++);
            _p2 = mesh.vertex(*hfv_it++);
        }

    protected:
        const Mesh &mesh_;
    };


    template<class Mesh>
    using OpenVolumeMeshTriangleBSPTraits = OVMOMCommonTriangleBSPTraits<Mesh, OVMSpecificTriangleBSPTraits<Mesh>>;

//== CLASS DEFINITION =========================================================
    template<class Mesh>
    class OpenVolumeMeshTriangleBSPT
            : public TriangleBSPT<OpenVolumeMeshTriangleBSPTraits<Mesh> > {
    public:
        typedef OpenVolumeMeshTriangleBSPTraits<Mesh> Traits;
        typedef TriangleBSPT<Traits> Base;
        typedef typename Traits::Scalar Scalar;

        OpenVolumeMeshTriangleBSPT(const Mesh &_mesh,
                                   const Scalar &_infinity = std::numeric_limits<Scalar>::infinity())
                : Base(Traits(_mesh), _infinity) {}
    };

}
//#endif


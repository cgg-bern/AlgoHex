/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

//== INCLUDES =================================================================

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Core/BaseEntities.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
//#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh>
#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.hh>
#include <OpenVolumeMesh/Attribs/StatusAttrib.hh>
#include "OpenVolumeMesh/Attribs/InterfaceAttrib.hh"

#include <Eigen/Geometry>
#include <Eigen/Core>
#include <Eigen/SparseCore>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace OVM = OpenVolumeMesh;

namespace AlgoHex
{
using VH = OVM::VertexHandle;
using HEH = OVM::HalfEdgeHandle;
using EH = OVM::EdgeHandle;
using HFH = OVM::HalfFaceHandle;
using FH = OVM::FaceHandle;
using CH = OVM::CellHandle;

using EIt = OVM::EdgeIter;
using HEIt = OVM::HalfEdgeIter;
using VIt = OVM::VertexIter;
using CIt = OVM::CellIter;
using FIt = OVM::FaceIter;
using BFIt = OVM::BoundaryFaceIter;
using HFIt = OVM::HalfFaceIter;

using VOHEIt = OVM::VertexOHalfEdgeIter;
using VVIt = OVM::VertexVertexIter;
using VEIt = OVM::VertexEdgeIter;
using VFIt = OVM::VertexFaceIter;
using VCIt = OVM::VertexCellIter;
using ECIt = OVM::EdgeCellIter;
using HECIt = OVM::HalfEdgeCellIter;
using HEHFIt = OVM::HalfEdgeHalfFaceIter;
using CVIt = OVM::CellVertexIter;
using CEIt = OVM::CellEdgeIter;
using CFIt = OVM::CellFaceIter;

using TVIt = OVM::TetVertexIter;


template<typename T>
using VP = OVM::VertexPropertyT<T>;
template<typename T>
using HEP = OVM::HalfEdgePropertyT<T>;
template<typename T>
using EP = OVM::EdgePropertyT<T>;
template<typename T>
using HFP = OVM::HalfFacePropertyT<T>;
template<typename T>
using FP = OVM::FacePropertyT<T>;
template<typename T>
using CP = OVM::CellPropertyT<T>;
template<typename T>
using MP = OVM::MeshPropertyT<T>;

using DI = std::pair<double, int>;


using Vec9d = Eigen::Matrix<double, 9, 1>;
using Quaternion = Eigen::Quaterniond;
using Vec2d = Eigen::Vector2d;
using Vec3d = Eigen::Vector3d;
using Vec3i = Eigen::Vector3i;
using Vec4i = Eigen::Vector4i;
using VecXd = Eigen::VectorXd;
using Mat9d = Eigen::Matrix<double, 9, 9>;
using Mat3d = Eigen::Matrix3d;
using Mat3x4d = Eigen::Matrix<double, 3, 4>;

using Triplet = Eigen::Triplet<double>;
using SMatXd = Eigen::SparseMatrix<double>;


//=============================================================================
} // namespace AlgoHex
//=============================================================================

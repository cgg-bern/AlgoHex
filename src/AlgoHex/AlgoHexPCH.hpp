/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

// Precompiled header for the AlgoHex library and demo targets (see
// CMakeLists.txt, step "A" in SPEEDUPS.md). It pre-parses the heavy
// third-party header stack that nearly every AlgoHex translation unit
// includes: OVM geometry/topology kernels + handles (via TypeDef.hh),
// Eigen, the CoMISo NSolver template headers and the QGP3D/MC3D
// property headers.
//
// The include order below matters:
//   1. CoMISo/Config/config.hh is a build-time generated header that
//      defines the COMISO_*_AVAILABLE macros and MUST precede every
//      other CoMISo header.
//   2. TypeDef.hh pulls in the OVM kernels and Eigen Core/Geometry/
//      SparseCore before anything else.

#include <CoMISo/Config/config.hh>

#include <AlgoHex/TypeDef.hh>

#include <Eigen/Dense>

#include <CoMISo/Base/Debug/DebConfig.hh>
#include <CoMISo/NSolver/NProblemInterface.hh>
#include <CoMISo/NSolver/FiniteElementHessianProjection.hh>
#include <CoMISo/NSolver/FiniteElementTinyAD.hh>
#include <CoMISo/NSolver/FiniteElementProblem.hh>
#include <CoMISo/NSolver/NewtonSolver.hh>
#include <CoMISo/NSolver/TruncatedNewtonPCG.hh>
#include <CoMISo/NSolver/TrustregionNewtonPCG.hh>
#include <CoMISo/NSolver/LinearConstraint.hh>
#include <CoMISo/NSolver/ConstraintTools.hh>
#include <CoMISo/NSolver/NPTiming.hh>
#include <CoMISo/Utils/Polynomials.hh>
#include <CoMISo/Utils/Tools.hh>

#include <QGP3D/Quantizer.hpp>
#include <QGP3D/PathConstraint.hpp>

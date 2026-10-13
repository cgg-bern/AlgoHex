/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once


//== INCLUDES =================================================================

#include <Eigen/Core>

#include <CoMISo/NSolver/FiniteElementTinyAD.hh>
#include <CoMISo/NSolver/FiniteElementHessianProjection.hh>

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================


/** \class AMIPSFrameElement3DTinyAD AMIPSFrameElement3DTinyAD.hh
 *
 *  @brief Implementation of a minus log det element for optimization. f(J) = w*(exp(s*( alpha*(tr(J^TA^TAJ)*(det(AJ)^(-2/3)-3) + beta*(det(AJ)+det(AJ)^(-1) -2) ))-1.0)
 *
 */

// element function for f(J) = w*(exp(s*( alpha*(tr(J^TA^TAJ)*(det(AJ)^(-2/3)-3) + beta*(det(AJ)+det(AJ)^(-1) -2) ))-1.0)
// matrix is col-major

class AMIPSFrameElement3D_F
{
public:

  // define dimensions
  const static int NV = 9; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 14; // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], det(A) = _c[9], weight w=_c[10], distortion parameters s=_c[11], alpha=_c[12] and beta=_c[13]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);
};

extern template double AMIPSFrameElement3D_F::eval_f<double>(const Eigen::Matrix<double, AMIPSFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, AMIPSFrameElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, false> AMIPSFrameElement3D_F::eval_f<TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, false>, AMIPSFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, AMIPSFrameElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, true> AMIPSFrameElement3D_F::eval_f<TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<AMIPSFrameElement3D_F::NV, double, true>, AMIPSFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, AMIPSFrameElement3D_F::NC, 1> &_c) const;

// create types for optimization with derivatives and projected hessian
using AMIPSFrameElement3D_TAD = COMISO::FiniteElementTinyAD<AMIPSFrameElement3D_F>;
using AMIPSFrameElement3D_TAD_PH = COMISO::FiniteElementHessianProjection<AMIPSFrameElement3D_TAD>;

class SymmetricDirichletFrameElement3D_F
{
public:

  // define dimensions
  const static int NV = 9; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 10; // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], det(A) = _c[9], weight w=_c[10], distortion parameters s=_c[11], alpha=_c[12] and beta=_c[13]


  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);
};

extern template double SymmetricDirichletFrameElement3D_F::eval_f<double>(const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, false> SymmetricDirichletFrameElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, false>, SymmetricDirichletFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, true> SymmetricDirichletFrameElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletFrameElement3D_F::NV, double, true>, SymmetricDirichletFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletFrameElement3D_F::NC, 1> &_c) const;

using SymmetricDirichletFrameElement3D_TAD = COMISO::FiniteElementTinyAD<SymmetricDirichletFrameElement3D_F>;
using SymmetricDirichletFrameElement3D_TAD_PH = COMISO::FiniteElementHessianProjection<SymmetricDirichletFrameElement3D_TAD>;

class SymmetricDirichletDualFrameElement3D_F
{
public:

  // define dimensions
  const static int NV = 21; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix, P=[p0,p1,p2,p3] = 3x4 matrix of tet vertex positions
  const static int NC = 10; // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], weight w=_c[9]

  using VecV = Eigen::Matrix<double, NV, 1>;
  using VecC = Eigen::Matrix<double, NC, 1>;

  template<class ScalarT>
  ScalarT eval_f(const Eigen::Matrix<ScalarT, NV, 1> &_x, const Eigen::Matrix<double, NC, 1> &_c) const;

  double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

};

extern template double SymmetricDirichletDualFrameElement3D_F::eval_f<double>(const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, false> SymmetricDirichletDualFrameElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, false>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, false>, SymmetricDirichletDualFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NC, 1> &_c) const;
extern template TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, true> SymmetricDirichletDualFrameElement3D_F::eval_f<TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, true>>(const Eigen::Matrix<TinyAD::Scalar<SymmetricDirichletDualFrameElement3D_F::NV, double, true>, SymmetricDirichletDualFrameElement3D_F::NV, 1> &_x, const Eigen::Matrix<double, SymmetricDirichletDualFrameElement3D_F::NC, 1> &_c) const;

using SymmetricDirichletDualFrameElement3D_TAD = COMISO::FiniteElementTinyAD<SymmetricDirichletDualFrameElement3D_F>;
using SymmetricDirichletDualFrameElement3D_TAD_PH = COMISO::FiniteElementHessianProjection<SymmetricDirichletDualFrameElement3D_TAD>;

//=============================================================================
} // namespace AlgoHex
//=============================================================================

namespace COMISO
{
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::AMIPSFrameElement3D_F>>;
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::SymmetricDirichletFrameElement3D_F>>;
extern template class FiniteElementHessianProjection<FiniteElementTinyAD<AlgoHex::SymmetricDirichletDualFrameElement3D_F>>;
} // namespace COMISO

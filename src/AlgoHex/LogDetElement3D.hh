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
#include <Eigen/SparseCore>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================


/** \class LogDetElement3D LogDetElement3D.hh
 *
 *  @brief Implementation of a minus log det element for optimization. f(x0,x1,x2,x3,x4,x5,x6,x7,x8) = - w*log ( det [[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]])
 *
 */

// element for f(x0,x1,x2,x3,x4,x5,x6,x7,x8) = - w*log ( det [[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]])
// matrix is col-major
class LogDetElement3D
{
public:

  // WARNING: this version is not the latest anymore, please use projected version LogDetElement3D_PH instead

  // define dimensions
  const static int NV = 9; // [[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 2; // [w,o] = weight w and offset o

  typedef Eigen::Matrix<size_t, NV, 1> VecI;
  typedef Eigen::Matrix<double, NV, 1> VecV;
  typedef Eigen::Matrix<double, NC, 1> VecC;
  typedef Eigen::Triplet<double> Triplet;

   double eval_f(const VecV &_x, const VecC &_c) const;

   void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const;

   void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const;

   double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

   double determinant(const VecV &_x) const;

   double grad_det(const VecV &_x, const int _i, const int _j) const;

   double hess_det(const VecV &_x, const int _i, const int _j, const int _k, const int _l) const;
};

// as the one before but with projected Hessian (positive semidefinite)
class LogDetElement3D_PH
{
public:

  // define dimensions
  const static int NV = 9; // [[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 3; // [w,o] = weight w, scale s and offset o

  typedef Eigen::Matrix<size_t, NV, 1> VecI;
  typedef Eigen::Matrix<double, NV, 1> VecV;
  typedef Eigen::Matrix<double, NC, 1> VecC;
  typedef Eigen::Triplet<double> Triplet;

   double eval_f(const VecV &_x, const VecC &_c) const;

   void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const;

   void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const;

   double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

   double determinant(const VecV &_x) const;

   double grad_det(const VecV &_x, const int _i, const int _j) const;

   double hess_det(const VecV &_x, const int _i, const int _j, const int _k, const int _l) const;

   bool diagonally_dominant(const Eigen::Matrix<double, 9, 9> &_A) const;

   void project_to_positive_semidefinite(Eigen::Matrix<double, 9, 9> &_A) const;
};


//=============================================================================
} // namespace AlgoHex
//=============================================================================

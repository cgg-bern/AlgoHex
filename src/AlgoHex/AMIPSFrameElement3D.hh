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


/** \class AMIPSFrameElement3D AMIPSFrameElement3D.hh
 *
 *  @brief Implementation of a minus log det element for optimization. f(J) = w*(exp(s*( alpha*(tr(J^TA^TAJ)*(det(AJ)^(-2/3)-3) + beta*(det(AJ)+det(AJ)^(-1) -2) ))-1.0)
 *
 */

// element for f(J) = w*(exp(s*( alpha*(tr(J^TA^TAJ)*(det(AJ)^(-2/3)-3) + beta*(det(AJ)+det(AJ)^(-1) -2) ))-1.0)
// matrix is col-major
class AMIPSFrameElement3D_PH
{
public:

  // define dimensions
  const static int NV = 9; // J=[[x0,x1,x2]^T | [x3,x4,x5]^T | [x6,x7,x8]^T]] = 3x3 matrix
  const static int NC = 14; // [A(3x3),w,s,alpha,beta] = shape matrix A = _c[0..8], det(A) = _c[9], weight w=_c[10], distortion parameters s=_c[11], alpha=_c[12] and beta=_c[13]

  typedef Eigen::Matrix<size_t, NV, 1> VecI;
  typedef Eigen::Matrix<double, NV, 1> VecV;
  typedef Eigen::Matrix<double, NC, 1> VecC;
  typedef Eigen::Triplet<double> Triplet;

   double eval_f(const VecV &_x, const VecC &_c) const;

   void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const;

   void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets);

   double fh(const VecV &_x, const VecC &_c) const;

   double grad_fh(const VecV &_x, const VecC &_c, const int _i, const int _j) const;

   double hess_fh(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const;

   double fe(const VecV &_x, const VecC &_c) const;

   double grad_fe(const VecV &_x, const VecC &_c, const int _i, const int _j) const;

   double hess_fe(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const;

   double fc(const VecV &_x, const VecC &_c) const;

   double grad_fc(const VecV &_x, const VecC &_c, const int _i, const int _j) const;

   double hess_fc(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const;

   double fb(const VecV &_x, const VecC &_c) const;

   double grad_fb(const VecV &_x, const VecC &_c, const int _i, const int _j) const;

   double hess_fb(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const;

   double fa(const VecV &_x, const VecC &_c) const;

   double grad_fa(const VecV &_x, const VecC &_c, const int _i, const int _j) const;

   double hess_fa(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const;

   double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c);

   double determinant(const VecV &_x) const;

   double grad_det(const VecV &_x, const int _i, const int _j) const;

   double hess_det(const VecV &_x, const int _i, const int _j, const int _k, const int _l) const;

   bool diagonally_dominant(const Eigen::Matrix<double, 9, 9> &_A) const;

   void project_to_positive_semidefinite(Eigen::Matrix<double, 9, 9> &_A) const;

private:
  // temporary
  Eigen::Matrix3d gfh_;
  Eigen::Matrix3d AJ_;
};


//=============================================================================
} // namespace AlgoHex
//=============================================================================

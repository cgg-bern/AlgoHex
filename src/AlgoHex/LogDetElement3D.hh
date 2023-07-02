/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once


//== INCLUDES =================================================================

#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/StdVector>
#include<cmath>

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

  inline double eval_f(const VecV &_x, const VecC &_c) const
  {
    double d = determinant(_x);

    if (d <= 0.0)
      return std::numeric_limits<double>::infinity();
    else
      return -_c[0] * std::log(d);
  }

  inline void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
  {
    double a = -_c[0] / determinant(_x);

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        _g[3 * j + i] = a * grad_det(_x, i, j);
  }

  inline void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const
  {
    _triplets.clear();

    double d = determinant(_x);
    double a = -_c[0] / d;
    double a2 = -a / d;

    // get gradient of determinant
    Eigen::Matrix3d gd;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        gd(i, j) = grad_det(_x, i, j);

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        for (unsigned int k = 0; k < 3; ++k)
          for (unsigned int l = 0; l < 3; ++l)
          {
            const int vi = 3 * j + i;
            const int vj = 3 * l + k;

            // only compute upper triangle of symmetric Hessian
            if (vi <= vj)
            {
              double v = a2 * gd(i, j) * gd(k, l) + a * hess_det(_x, i, j, k, l);

              _triplets.push_back(Triplet(vi, vj, v));

              if (vi != vj)
                _triplets.push_back(Triplet(vj, vi, v));
            }
          }
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo: see https://en.wikipedia.org/wiki/Cubic_equation
    // get coefficients in maple via "collect(Determinant(V*t+G), t)"

    return DBL_MAX;
  }

  inline double determinant(const VecV &_x) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > G((double *) &(_x[0]));

    return G.determinant();
  }

  inline double grad_det(const VecV &_x, const int _i, const int _j) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > G((double *) &(_x[0]));

    const int i0 = (_i + 1) % 3;
    const int i1 = (_i + 2) % 3;
    const int j0 = (_j + 1) % 3;
    const int j1 = (_j + 2) % 3;

    double d = G(i0, j0) * G(i1, j1) - G(i0, j1) * G(i1, j0);

    if ((_i + _j) % 1)
      return -d;
    else
      return d;
  }

  inline double hess_det(const VecV &_x, const int _i, const int _j, const int _k, const int _l) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > G((double *) &(_x[0]));

    const int i0 = (_i + 1) % 3;
    const int i1 = (_i + 2) % 3;
    const int j0 = (_j + 1) % 3;
    const int j1 = (_j + 2) % 3;

//      first derivative w.r.t. (i,j)
//      double d = G(i0,j0)*G(i1,j1)-G(i0,j1)*G(i1,j0);
//
//      if( (_i+_j) % 1)
//        return -d;
//      else
//        return d;

    // get sign
    double s = 1.0;
    if ((_i + _j) % 1)
      s = -1.0;

    if (_k == i0)
    {
      if (_l == j0)
        return s * G(i1, j1);
      else if (_l == j1)
        return -s * G(i1, j0);
      else return 0.0;
    }
    else if (_k == i1)
    {
      if (_l == j0)
        return -s * G(i0, j1);
      else if (_l == j1)
        return s * G(i0, j0);
      else return 0.0;
    }
    else return 0.0;
  }
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

  inline double eval_f(const VecV &_x, const VecC &_c) const
  {
    double d = _c[1] * determinant(_x) - _c[2];

    if (d <= 0.0)
      return std::numeric_limits<double>::infinity();
    else
      return -_c[0] * std::log(d);
  }

  inline void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
  {
    double a = -_c[0] * _c[1] / (_c[1] * determinant(_x) - _c[2]);

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        _g[3 * j + i] = a * grad_det(_x, i, j);
  }

  inline void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const
  {
    _triplets.clear();

    double d = _c[1] * determinant(_x) - _c[2];
    double a = -_c[0] * _c[1] / d;
    double a2 = -_c[1] * a / d;

    // get gradient of determinant
    Eigen::Matrix3d gd;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        gd(i, j) = grad_det(_x, i, j);

    Eigen::Matrix<double, 9, 9> H;

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        for (unsigned int k = 0; k < 3; ++k)
          for (unsigned int l = 0; l < 3; ++l)
          {
            const int vi = 3 * j + i;
            const int vj = 3 * l + k;

            // only compute upper triangle of symmetric Hessian
            if (vi <= vj)
            {
              double v = a2 * gd(i, j) * gd(k, l) + a * hess_det(_x, i, j, k, l);

              H(vi, vj) = v;
//                _triplets.push_back(Triplet(vi, vj, v));

              if (vi != vj)
                H(vj, vi) = v;
//                  _triplets.push_back(Triplet(vj, vi, v));
            }
          }

    if (!diagonally_dominant(H))
      project_to_positive_semidefinite(H);

    for (unsigned int i = 0; i < 9; ++i)
      for (unsigned int j = 0; j < 9; ++j)
        _triplets.push_back(Triplet(i, j, H(i, j)));
  }

  inline double max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo: see https://en.wikipedia.org/wiki/Cubic_equation
    // get coefficients in maple via "collect(Determinant(V*t+G), t)"

    return DBL_MAX;
  }

  inline double determinant(const VecV &_x) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > G((double *) &(_x[0]));

    return G.determinant();
  }

  inline double grad_det(const VecV &_x, const int _i, const int _j) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > G((double *) &(_x[0]));

    const int i0 = (_i + 1) % 3;
    const int i1 = (_i + 2) % 3;
    const int j0 = (_j + 1) % 3;
    const int j1 = (_j + 2) % 3;

    double d = G(i0, j0) * G(i1, j1) - G(i0, j1) * G(i1, j0);

    if ((_i + _j) % 1)
      return -d;
    else
      return d;
  }

  inline double hess_det(const VecV &_x, const int _i, const int _j, const int _k, const int _l) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > G((double *) &(_x[0]));

    const int i0 = (_i + 1) % 3;
    const int i1 = (_i + 2) % 3;
    const int j0 = (_j + 1) % 3;
    const int j1 = (_j + 2) % 3;

//      first derivative w.r.t. (i,j)
//      double d = G(i0,j0)*G(i1,j1)-G(i0,j1)*G(i1,j0);
//
//      if( (_i+_j) % 1)
//        return -d;
//      else
//        return d;

    // get sign
    double s = 1.0;
    if ((_i + _j) % 1)
      s = -1.0;

    if (_k == i0)
    {
      if (_l == j0)
        return s * G(i1, j1);
      else if (_l == j1)
        return -s * G(i1, j0);
      else return 0.0;
    }
    else if (_k == i1)
    {
      if (_l == j0)
        return -s * G(i0, j1);
      else if (_l == j1)
        return s * G(i0, j0);
      else return 0.0;
    }
    else return 0.0;
  }

  template<class MatT>
  inline bool diagonally_dominant(const MatT &_A) const
  {
    // only works for square matrices
    assert(_A.cols() == _A.rows());

    for (unsigned int i = 0; i < _A.rows(); ++i)
    {
      double asum(0.0);
      for (unsigned int j = 0; j < _A.cols(); ++j)
        asum += std::abs(_A(i, j));

      double d = std::abs(_A(i, i));

      if (2 * d < asum)
        return false;
    }

    return true;
  }

  template<class MatT>
  inline void project_to_positive_semidefinite(MatT &_A) const
  {
    typename Eigen::SelfAdjointEigenSolver<MatT> es(_A);

    typename Eigen::SelfAdjointEigenSolver<MatT>::RealVectorType ev = es.eigenvalues();
//      ev.cwiseMax(0.0); // ToDo: why does cwiseMax not work?
    for (unsigned int i = 0; i < ev.size(); ++i)
      ev[i] = std::max(ev[i], 0.0);

    _A = es.eigenvectors() * ev.asDiagonal() * es.eigenvectors().transpose();

//      // DEBUG
//      {
//        Eigen::SelfAdjointEigenSolver<MatT> es2(_A);
//
//        std::cerr << "----------------------------------------------\n";
//        std::cerr << "eigenvalues original : " << es.eigenvalues().transpose() << std::endl;
//        std::cerr << "eigenvalues projected: " << es2.eigenvalues().transpose() << std::endl;
//      }
  }
};


//=============================================================================
} // namespace AlgoHex
//=============================================================================


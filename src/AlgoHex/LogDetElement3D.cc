/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include "LogDetElement3D.hh"
#include <cassert>
#include <cfloat>

namespace AlgoHex
{

   double LogDetElement3D::eval_f(const VecV &_x, const VecC &_c) const
  {
    double d = determinant(_x);

    if (d <= 0.0)
      return std::numeric_limits<double>::infinity();
    else
      return -_c[0] * std::log(d);
  }

   void LogDetElement3D::eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
  {
    double a = -_c[0] / determinant(_x);

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        _g[3 * j + i] = a * grad_det(_x, i, j);
  }

   void LogDetElement3D::eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const
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

   double LogDetElement3D::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo: see https://en.wikipedia.org/wiki/Cubic_equation
    // get coefficients in maple via "collect(Determinant(V*t+G), t)"

    return DBL_MAX;
  }

   double LogDetElement3D::determinant(const VecV &_x) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > G((double *) &(_x[0]));

    return G.determinant();
  }

   double LogDetElement3D::grad_det(const VecV &_x, const int _i, const int _j) const
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

   double LogDetElement3D::hess_det(const VecV &_x, const int _i, const int _j, const int _k, const int _l) const
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

   double LogDetElement3D_PH::eval_f(const VecV &_x, const VecC &_c) const
  {
    double d = _c[1] * determinant(_x) - _c[2];

    if (d <= 0.0)
      return std::numeric_limits<double>::infinity();
    else
      return -_c[0] * std::log(d);
  }

   void LogDetElement3D_PH::eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
  {
    double a = -_c[0] * _c[1] / (_c[1] * determinant(_x) - _c[2]);

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        _g[3 * j + i] = a * grad_det(_x, i, j);
  }

   void LogDetElement3D_PH::eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets) const
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

   double LogDetElement3D_PH::max_feasible_step(const VecV &_x, const VecV &_v, const VecC &_c)
  {
    // ToDo: see https://en.wikipedia.org/wiki/Cubic_equation
    // get coefficients in maple via "collect(Determinant(V*t+G), t)"

    return DBL_MAX;
  }

   double LogDetElement3D_PH::determinant(const VecV &_x) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > G((double *) &(_x[0]));

    return G.determinant();
  }

   double LogDetElement3D_PH::grad_det(const VecV &_x, const int _i, const int _j) const
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

   double LogDetElement3D_PH::hess_det(const VecV &_x, const int _i, const int _j, const int _k, const int _l) const
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

   bool LogDetElement3D_PH::diagonally_dominant(const Eigen::Matrix<double, 9, 9> &_A) const
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

   void LogDetElement3D_PH::project_to_positive_semidefinite(Eigen::Matrix<double, 9, 9> &_A) const
  {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 9, 9>> es(_A);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 9, 9>>::RealVectorType ev = es.eigenvalues();
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


}

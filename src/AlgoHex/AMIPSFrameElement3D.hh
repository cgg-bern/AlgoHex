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

  inline double eval_f(const VecV &_x, const VecC &_c) const
  {
    // split function into
    // f(x) = w*(exp(fh(x))-1.0)
    // with fh(x) = s*(alpha*fe(x)-3.0)+beta*(fb(x)+fc(x)-2.0))
    //      fe(x) = fa(x)*fb(x)^(-2/3)
    //      fc(x) = fb(x)^(-1)
    //      fb(x) = det(AJ) = det(A)*det(J)
    //      fa(x) = tr(J^T A^T A J)

    const double &w = _c[10];

    return w * (std::exp(fh(_x, _c)) - 1.0);
  }

  inline void eval_gradient(const VecV &_x, const VecC &_c, VecV &_g) const
  {
    // get constants
    const double &w = _c[10];
    const double z = w * std::exp(fh(_x, _c));

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        _g[3 * j + i] = z * grad_fh(_x, _c, i, j);
  }

  inline void eval_hessian(const VecV &_x, const VecC &_c, std::vector<Triplet> &_triplets)
  {
    _triplets.clear();

    // get constants
    const double &w = _c[10];
    const double z = w * std::exp(fh(_x, _c));

    // get gradient of fh
    // Eigen::Matrix3d gfh;
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        gfh_(i, j) = grad_fh(_x, _c, i, j);

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
              double v = z * (gfh_(i, j) * gfh_(k, l) + hess_fh(_x, _c, i, j, k, l));

              H(vi, vj) = v;
//                _triplets.push_back(Triplet(vi, vj, v));

              if (vi != vj)
                H(vj, vi) = v;
//                _triplets.push_back(Triplet(vj, vi, v));
            }
          }


    if (!diagonally_dominant(H))
      project_to_positive_semidefinite(H);

    for (unsigned int i = 0; i < 9; ++i)
      for (unsigned int j = 0; j < 9; ++j)
        _triplets.push_back(Triplet(i, j, H(i, j)));

  }

  inline double fh(const VecV &_x, const VecC &_c) const
  {
    // fh(x) = s*(alpha*fe(x)-3.0)+beta*(fb(x)+fc(x)-2.0))

    // get constants
    const double &s = _c[11];
    const double &alpha = _c[12];
    const double &beta = _c[13];

    return s * (alpha * (fe(_x, _c) - 3.0) + beta * (fb(_x, _c) + fc(_x, _c) - 2.0));
  }

  inline double grad_fh(const VecV &_x, const VecC &_c, const int _i, const int _j) const
  {
    // get constants
    const double &s = _c[11];
    const double &alpha = _c[12];
    const double &beta = _c[13];

    double sa = s * alpha;
    double sb = s * beta;

    return sa * grad_fe(_x, _c, _i, _j) + sb * (grad_fb(_x, _c, _i, _j) + grad_fc(_x, _c, _i, _j));
  }

  inline double hess_fh(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const
  {
    // get constants
    const double &s = _c[11];
    const double &alpha = _c[12];
    const double &beta = _c[13];

    double sa = s * alpha;
    double sb = s * beta;

    return sa * hess_fe(_x, _c, _i, _j, _k, _l) +
           sb * (hess_fb(_x, _c, _i, _j, _k, _l) + hess_fc(_x, _c, _i, _j, _k, _l));
  }

  inline double fe(const VecV &_x, const VecC &_c) const
  {
    return fa(_x, _c) * std::pow(fb(_x, _c), -2.0 / 3.0);
  }

  inline double grad_fe(const VecV &_x, const VecC &_c, const int _i, const int _j) const
  {
    const double fav = fa(_x, _c);
    const double fbv = fb(_x, _c);

    return std::pow(fbv, -2.0 / 3.0) * grad_fa(_x, _c, _i, _j) -
           2.0 / 3.0 * fav * std::pow(fbv, -5.0 / 3.0) * grad_fb(_x, _c, _i, _j);
  }

  inline double hess_fe(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const
  {
    const double fav = fa(_x, _c);
    const double fbv = fb(_x, _c);

    const double gav0 = grad_fa(_x, _c, _i, _j);
    const double gav1 = grad_fa(_x, _c, _k, _l);
    const double gbv0 = grad_fb(_x, _c, _i, _j);
    const double gbv1 = grad_fb(_x, _c, _k, _l);

    double v = std::pow(fbv, -2.0 / 3.0) * hess_fa(_x, _c, _i, _j, _k, _l);
    v += -2.0 / 3.0 * std::pow(fbv, -5.0 / 3.0) * (gav0 * gbv1 + gav1 * gbv0 + fav * hess_fb(_x, _c, _i, _j, _k, _l));
    v += 10.0 / 9.0 * fav * std::pow(fbv, -8.0 / 3.0) * gbv0 * gbv1;

    return v;
  }

  inline double fc(const VecV &_x, const VecC &_c) const
  {
    return 1.0 / fb(_x, _c);
  }

  inline double grad_fc(const VecV &_x, const VecC &_c, const int _i, const int _j) const
  {
    double fbv = fb(_x, _c);

    return -std::pow(fbv, -2) * grad_fb(_x, _c, _i, _j);
  }

  inline double hess_fc(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const
  {
    double fbv = fb(_x, _c);

    return 2.0 * std::pow(fbv, -3) * grad_fb(_x, _c, _i, _j) * grad_fb(_x, _c, _k, _l) -
           std::pow(fbv, -2) * hess_fb(_x, _c, _i, _j, _k, _l);
  }

  inline double fb(const VecV &_x, const VecC &_c) const
  {
    // fb(x) = det(AJ) = det(A)*det(J)

    // get constants
    const double &dA = _c[9];

    const double d = determinant(_x);
    if (d <= 0.0)
      return std::numeric_limits<double>::infinity();
    else
      return d * dA;
  }

  inline double grad_fb(const VecV &_x, const VecC &_c, const int _i, const int _j) const
  {
    // get constants
    const double &dA = _c[9];

    return dA * grad_det(_x, _i, _j);
  }

  inline double hess_fb(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const
  {
    // get constants
    const double &dA = _c[9];

    return dA * hess_det(_x, _i, _j, _k, _l);
  }

  inline double fa(const VecV &_x, const VecC &_c) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > J((double *) &(_x[0]));
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > A((double *) &(_c[0]));

    Mat3d AJ;
    AJ = A * J;

    return (AJ.transpose() * AJ).trace();
  }

  inline double grad_fa(const VecV &_x, const VecC &_c, const int _i, const int _j) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > J((double *) &(_x[0]));
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > A((double *) &(_c[0]));

    double v(0.0);
    for (int m = 0; m < 3; ++m)
      for (int n = 0; n < 3; ++n)
        v += A(m, n) * J(n, _j) * A(m, _i);

    return 2.0 * v;
  }

  inline double hess_fa(const VecV &_x, const VecC &_c, const int _i, const int _j, const int _k, const int _l) const
  {
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > J((double *) &(_x[0]));
    Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::ColMajor> > A((double *) &(_c[0]));

    if (_j == _l)
    {
      double v(0.0);

      for (int m = 0; m < 3; ++m)
        v += A(m, _k) * A(m, _i);

      return 2.0 * v;
    }

    return 0.0;
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


private:
  // temporary
  Eigen::Matrix3d gfh_;
  Eigen::Matrix3d AJ_;
};


//=============================================================================
} // namespace AlgoHex
//=============================================================================


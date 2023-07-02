/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include <AlgoHex/TypeDef.hh>


namespace AlgoHex
{
//cross product of the edge vector and orthogonal axis
class PerpendicularEnergyCPRD
{
public:
  static int n_unknowns() { return 6; }

  static double eval_f(const Vec3d &_n, const double *_x)
  {
    double helper_0 = _x[0] - _x[3];
    double helper_1 = _x[1] - _x[4];
    double helper_2 = _x[2] - _x[5];

    return std::pow(helper_0 * _n[1] - helper_1 * _n[0], 2) + std::pow(helper_0 * _n[2] - helper_2 * _n[0], 2) +
           std::pow(helper_1 * _n[2] - helper_2 * _n[1], 2);
  }

  static void eval_gradient(const Vec3d &_n, const double *_x, double *_g)
  {
    double helper_0 = _x[0] - _x[3];
    double helper_1 = _x[1] - _x[4];
    double helper_2 = 2 * helper_0 * _n[1] - 2 * helper_1 * _n[0];
    double helper_3 = helper_2 * _n[1];
    double helper_4 = _x[2] - _x[5];
    double helper_5 = helper_0 * _n[2] - helper_4 * _n[0];
    double helper_6 = 2 * _n[2];
    double helper_7 = helper_5 * helper_6;
    double helper_8 = helper_2 * _n[0];
    double helper_9 = helper_1 * _n[2] - helper_4 * _n[1];
    double helper_10 = helper_6 * helper_9;
    double helper_11 = 2 * helper_5 * _n[0];
    double helper_12 = 2 * helper_9 * _n[1];

    _g[0] = helper_3 + helper_7;
    _g[1] = helper_10 - helper_8;
    _g[2] = -helper_11 - helper_12;
    _g[3] = -helper_3 - helper_7;
    _g[4] = -helper_10 + helper_8;
    _g[5] = helper_11 + helper_12;
  }

  static void eval_hessian(const Vec3d &_n, const double *_x, Eigen::MatrixXd &_H)
  {
    double helper_0 = 2 * std::pow(_n[1], 2);
    double helper_1 = 2 * std::pow(_n[2], 2);
    double helper_2 = helper_0 + helper_1;
    double helper_3 = 2 * _n[0];
    double helper_4 = helper_3 * _n[1];
    double helper_5 = -helper_4;
    double helper_6 = helper_3 * _n[2];
    double helper_7 = -helper_6;
    double helper_8 = -helper_0;
    double helper_9 = -helper_1;
    double helper_10 = helper_8 + helper_9;
    double helper_11 = 2 * std::pow(_n[0], 2);
    double helper_12 = helper_1 + helper_11;
    double helper_13 = 2 * _n[1] * _n[2];
    double helper_14 = -helper_13;
    double helper_15 = -helper_11;
    double helper_16 = helper_15 + helper_9;
    double helper_17 = helper_0 + helper_11;
    double helper_18 = helper_15 + helper_8;

    _H(0, 0) = helper_2;
    _H(0, 1) = helper_5;
    _H(0, 2) = helper_7;
    _H(0, 3) = helper_10;
    _H(0, 4) = helper_4;
    _H(0, 5) = helper_6;
    _H(1, 0) = helper_5;
    _H(1, 1) = helper_12;
    _H(1, 2) = helper_14;
    _H(1, 3) = helper_4;
    _H(1, 4) = helper_16;
    _H(1, 5) = helper_13;
    _H(2, 0) = helper_7;
    _H(2, 1) = helper_14;
    _H(2, 2) = helper_17;
    _H(2, 3) = helper_6;
    _H(2, 4) = helper_13;
    _H(2, 5) = helper_18;
    _H(3, 0) = helper_10;
    _H(3, 1) = helper_4;
    _H(3, 2) = helper_6;
    _H(3, 3) = helper_2;
    _H(3, 4) = helper_5;
    _H(3, 5) = helper_7;
    _H(4, 0) = helper_4;
    _H(4, 1) = helper_16;
    _H(4, 2) = helper_13;
    _H(4, 3) = helper_5;
    _H(4, 4) = helper_12;
    _H(4, 5) = helper_14;
    _H(5, 0) = helper_6;
    _H(5, 1) = helper_13;
    _H(5, 2) = helper_18;
    _H(5, 3) = helper_7;
    _H(5, 4) = helper_14;
    _H(5, 5) = helper_17;


    // Compute Eigen decomposition H = V D V^T
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_H);
    const Eigen::MatrixXd &V = solver.eigenvectors();
    Eigen::VectorXd evals = solver.eigenvalues();

    // check eigenvalues against epsilon, for those less that eps (zero)
    const int n = n_unknowns();
    for (int i = 0; i < n; ++i)
    {
      if (evals[i] < 1e-8)
      {
        evals[i] = 1e-8;
      }
    }

    // compute correction matrix M = V * diag(m) * V^T
    _H.noalias() = V * evals.asDiagonal() * V.transpose();
  }
};


class ComplexEdgeLength
{
public:
  static int n_unknowns() { return 6; }

  static double eval_f(const double *_x)
  {
    return std::pow(_x[0] - _x[3], 2) + std::pow(_x[1] - _x[4], 2) + std::pow(_x[2] - _x[5], 2);
  }

  static void eval_gradient(const double *_x, double *_g)
  {
    double helper_0 = 2 * _x[0];
    double helper_1 = 2 * _x[3];
    double helper_2 = 2 * _x[1];
    double helper_3 = 2 * _x[4];
    double helper_4 = 2 * _x[2];
    double helper_5 = 2 * _x[5];

    _g[0] = helper_0 - helper_1;
    _g[1] = helper_2 - helper_3;
    _g[2] = helper_4 - helper_5;
    _g[3] = -helper_0 + helper_1;
    _g[4] = -helper_2 + helper_3;
    _g[5] = -helper_4 + helper_5;
  }

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H)
  {
    _H.setZero();

    _H(0, 0) = 2;
    _H(0, 3) = -2;
    _H(1, 1) = 2;
    _H(1, 4) = -2;
    _H(2, 2) = 2;
    _H(2, 5) = -2;
    _H(3, 0) = -2;
    _H(3, 3) = 2;
    _H(4, 1) = -2;
    _H(4, 4) = 2;
    _H(5, 2) = -2;
    _H(5, 5) = 2;
  }
};


class CurvatureSmooth
{
public:
  static int n_unknowns() { return 9; }

  static double eval_f(const double *_x)
  {
    return std::pow(_x[0] - 2 * _x[3] + _x[6], 2) + std::pow(_x[1] - 2 * _x[4] + _x[7], 2) +
           std::pow(_x[2] - 2 * _x[5] + _x[8], 2);
  }

  static void eval_gradient(const double *_x, double *_g)
  {
    memset(_g, 0, 9 * sizeof(double));

    double helper_0 = 2 * _x[0] - 4 * _x[3] + 2 * _x[6];
    double helper_1 = 2 * _x[1] - 4 * _x[4] + 2 * _x[7];
    double helper_2 = 2 * _x[2] - 4 * _x[5] + 2 * _x[8];

    _g[0] = helper_0;
    _g[1] = helper_1;
    _g[2] = helper_2;
    _g[3] = -4 * _x[0] + 8 * _x[3] - 4 * _x[6];
    _g[4] = -4 * _x[1] + 8 * _x[4] - 4 * _x[7];
    _g[5] = -4 * _x[2] + 8 * _x[5] - 4 * _x[8];
    _g[6] = helper_0;
    _g[7] = helper_1;
    _g[8] = helper_2;
  }

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H)
  {
    _H.setZero();

    _H(0, 0) = 2;
    _H(0, 3) = -4;
    _H(0, 6) = 2;

    _H(1, 1) = 2;
    _H(1, 4) = -4;
    _H(1, 7) = 2;

    _H(2, 2) = 2;
    _H(2, 5) = -4;
    _H(2, 8) = 2;

    _H(3, 0) = -4;
    _H(3, 3) = 8;
    _H(3, 6) = -4;

    _H(4, 1) = -4;
    _H(4, 4) = 8;
    _H(4, 7) = -4;

    _H(5, 2) = -4;
    _H(5, 5) = 8;
    _H(5, 8) = -4;

    _H(6, 0) = 2;
    _H(6, 3) = -4;
    _H(6, 6) = 2;

    _H(7, 1) = 2;
    _H(7, 4) = -4;
    _H(7, 7) = 2;

    _H(8, 2) = 2;
    _H(8, 5) = -4;
    _H(8, 8) = 2;
  }

};


class RepulsionEnergy
{
public:
  static int n_unknowns() { return 3; }

  static double eval_f(const double *_x, const Vec3d &_target_pt)
  {
    double energy(0);
    double dx0 = _x[0] - _target_pt[0];
    double dx1 = _x[1] - _target_pt[1];
    double dx2 = _x[2] - _target_pt[2];

    energy += dx0 * dx0;
    energy += dx1 * dx1;
    energy += dx2 * dx2;

    return 0.5 * energy;
  }

  static void eval_gradient(const double *_x, const Vec3d &_target_pt, double *_g)
  {
    _g[0] = _x[0] - _target_pt[0];
    _g[1] = _x[1] - _target_pt[1];
    _g[2] = _x[2] - _target_pt[2];
  }

  static void eval_hessian(Mat3d &_H)
  {
    _H << 1, 0, 0,
            0, 1, 0,
            0, 0, 1;
  }
};


//Inverse mean ratio metric energy
class IMRMEnergy
{
public:
  static int n_unknowns() { return 3; }

  inline static double eval_f(const double *_x)
  {
    //trace
    double helper_0 = -_x[0];
    double helper_1 = helper_0 + _x[3];
    double helper_2 = -_x[1];
    double helper_3 = helper_2 + _x[7];
    double helper_4 = -_x[2];
    double helper_5 = helper_4 + _x[11];
    double helper_6 = helper_0 + _x[6];
    double helper_7 = helper_2 + _x[10];
    double helper_8 = helper_4 + _x[5];
    double helper_9 = helper_0 + _x[9];
    double helper_10 = helper_2 + _x[4];
    double helper_11 = helper_4 + _x[8];

    return std::log(std::pow(-helper_1 * helper_11 * helper_7 + helper_1 * helper_3 * helper_5 +
                             helper_10 * helper_11 * helper_9 - helper_10 * helper_5 * helper_6 -
                             helper_3 * helper_8 * helper_9 + helper_6 * helper_7 * helper_8,
                             -0.66666666666666663) *
                    (std::pow(_x[0] - _x[3], 2) + std::pow(_x[1] - _x[4], 2) + std::pow(_x[2] - _x[5], 2) +
                     1.3333333333333333 * std::pow(0.5 * _x[0] + 0.5 * _x[3] - _x[9], 2) +
                     1.3333333333333333 * std::pow(0.5 * _x[1] - _x[10] + 0.5 * _x[4], 2) +
                     1.3333333333333333 * std::pow(-_x[11] + 0.5 * _x[2] + 0.5 * _x[5], 2) +
                     1.5000000000000002 * std::pow(
                             0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] - _x[6] +
                             0.33333333333333326 * _x[9], 2) + 1.5000000000000002 * std::pow(
                            0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] +
                            0.33333333333333326 * _x[4] - _x[7], 2) + 1.5000000000000002 * std::pow(
                            0.33333333333333326 * _x[11] + 0.33333333333333348 * _x[2] +
                            0.33333333333333326 * _x[5] - _x[8], 2))) - 0.2310490601866485;
  }

  inline static void eval_gradient(const double *_x, Vec3d &_g)
  {
    //trace
    double helper_0 = -1.4142135623730951 * _x[0] + 1.4142135623730951 * _x[3];
    double helper_1 = -_x[2];
    double helper_2 = helper_1 + _x[11];
    double helper_3 = -_x[1];
    double helper_4 = helper_3 + _x[7];
    double helper_5 = helper_2 * helper_4;
    double helper_6 = -1.4142135623730951 * _x[2] + 1.4142135623730951 * _x[5];
    double helper_7 = helper_3 + _x[10];
    double helper_8 = -_x[0];
    double helper_9 = helper_8 + _x[6];
    double helper_10 = helper_7 * helper_9;
    double helper_11 = -1.4142135623730951 * _x[1] + 1.4142135623730951 * _x[4];
    double helper_12 = helper_1 + _x[8];
    double helper_13 = helper_8 + _x[9];
    double helper_14 = helper_12 * helper_13;
    double helper_15 = helper_0 * helper_12;
    double helper_16 = helper_11 * helper_2;
    double helper_17 = helper_4 * helper_6;
    double helper_18 =
            helper_0 * helper_5 + helper_10 * helper_6 + helper_11 * helper_14 - helper_13 * helper_17 -
            helper_15 * helper_7 - helper_16 * helper_9;
    double helper_19 = std::pow(helper_18, 0.66666666666666663);
    double helper_20 = 1.0 / helper_19;
    double helper_21 = -_x[10];
    double helper_22 = helper_21 + _x[1];
    double helper_23 = -_x[8];
    double helper_24 = 0.66666666666666663 * helper_11;
    double helper_25 = 0.66666666666666663 * helper_6;
    double helper_26 = -_x[9];
    double helper_27 = -_x[11];
    double helper_28 = -_x[6];
    double helper_29 = -_x[7];
    double helper_30 =
            std::pow(helper_1 + _x[5], 2) + std::pow(helper_3 + _x[4], 2) + std::pow(helper_8 + _x[3], 2) +
            1.3333333333333333 * std::pow(helper_21 + 0.5 * _x[1] + 0.5 * _x[4], 2) +
            1.3333333333333333 * std::pow(helper_26 + 0.5 * _x[0] + 0.5 * _x[3], 2) +
            1.3333333333333333 * std::pow(helper_27 + 0.5 * _x[2] + 0.5 * _x[5], 2) + 1.5000000000000002 *
                                                                                      std::pow(helper_23 +
                                                                                               0.33333333333333326 *
                                                                                               _x[11] +
                                                                                               0.33333333333333348 *
                                                                                               _x[2] +
                                                                                               0.33333333333333326 *
                                                                                               _x[5], 2) +
            1.5000000000000002 * std::pow(
                    helper_28 + 0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] +
                    0.33333333333333326 * _x[9], 2) + 1.5000000000000002 * std::pow(
                    helper_29 + 0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] +
                    0.33333333333333326 * _x[4], 2);
    double helper_31 = std::pow(helper_18, -1.6666666666666665) * helper_30;
    double helper_32 = helper_19 / helper_30;
    double helper_33 = helper_27 + _x[2];
    double helper_34 = 0.66666666666666663 * helper_0;
    double helper_35 = helper_29 + _x[1];

    _g[0] = helper_32 * (helper_20 * (3.0 * _x[0] - 1.0 * _x[3] - 1.0000000000000007 * _x[6] -
                                      0.99999999999999978 * _x[9]) + helper_31 *
                                                                     (0.94280904158206336 * helper_12 *
                                                                      helper_22 -
                                                                      0.66666666666666663 * helper_16 -
                                                                      0.66666666666666663 * helper_17 -
                                                                      helper_22 * helper_25 -
                                                                      helper_24 * (helper_23 + _x[2]) +
                                                                      0.94280904158206336 * helper_5));
    _g[1] = helper_32 * (helper_20 * (3.0 * _x[1] - 0.99999999999999978 * _x[10] - 1.0 * _x[4] -
                                      1.0000000000000007 * _x[7]) + helper_31 * (-helper_13 * helper_25 +
                                                                                 0.94280904158206336 *
                                                                                 helper_14 -
                                                                                 0.66666666666666663 *
                                                                                 helper_15 - helper_25 *
                                                                                             (helper_28 +
                                                                                              _x[0]) -
                                                                                 helper_33 * helper_34 +
                                                                                 0.94280904158206336 *
                                                                                 helper_33 * helper_9));
    _g[2] = helper_32 * (helper_20 * (-0.99999999999999978 * _x[11] + 3.0 * _x[2] - 1.0 * _x[5] -
                                      1.0000000000000007 * _x[8]) + helper_31 *
                                                                    (0.94280904158206336 * helper_10 +
                                                                     0.94280904158206336 * helper_13 *
                                                                     helper_35 - helper_24 * helper_9 -
                                                                     helper_24 * (helper_26 + _x[0]) -
                                                                     helper_34 * helper_35 -
                                                                     helper_34 * helper_7));
  }


  inline static void eval_hessian(const double *_x, Mat3d &_H)
  {
    _H.setZero();

    //trace
    double helper_0 = 3.0 * _x[0];
    double helper_1 = 1.0 * _x[3];
    double helper_2 = 1.0000000000000007 * _x[6];
    double helper_3 = 0.99999999999999978 * _x[9];
    double helper_4 = helper_0 - helper_1 - helper_2 - helper_3;
    double helper_5 = -1.4142135623730951 * _x[0] + 1.4142135623730951 * _x[3];
    double helper_6 = -_x[2];
    double helper_7 = helper_6 + _x[11];
    double helper_8 = -_x[1];
    double helper_9 = helper_8 + _x[7];
    double helper_10 = helper_7 * helper_9;
    double helper_11 = -1.4142135623730951 * _x[2] + 1.4142135623730951 * _x[5];
    double helper_12 = helper_8 + _x[10];
    double helper_13 = -_x[0];
    double helper_14 = helper_13 + _x[6];
    double helper_15 = helper_12 * helper_14;
    double helper_16 = -1.4142135623730951 * _x[1] + 1.4142135623730951 * _x[4];
    double helper_17 = helper_6 + _x[8];
    double helper_18 = helper_13 + _x[9];
    double helper_19 = helper_17 * helper_18;
    double helper_20 = helper_17 * helper_5;
    double helper_21 = helper_16 * helper_7;
    double helper_22 = helper_11 * helper_9;
    double helper_23 =
            helper_10 * helper_5 + helper_11 * helper_15 - helper_12 * helper_20 - helper_14 * helper_21 +
            helper_16 * helper_19 - helper_18 * helper_22;
    double helper_24 = std::pow(helper_23, 0.66666666666666663);
    double helper_25 = 1.0 / helper_24;
    double helper_26 = std::pow(helper_13 + _x[3], 2);
    double helper_27 = std::pow(helper_8 + _x[4], 2);
    double helper_28 = std::pow(helper_6 + _x[5], 2);
    double helper_29 = -_x[9];
    double helper_30 = std::pow(helper_29 + 0.5 * _x[0] + 0.5 * _x[3], 2);
    double helper_31 = -_x[10];
    double helper_32 = std::pow(helper_31 + 0.5 * _x[1] + 0.5 * _x[4], 2);
    double helper_33 = -_x[11];
    double helper_34 = std::pow(helper_33 + 0.5 * _x[2] + 0.5 * _x[5], 2);
    double helper_35 = -_x[6];
    double helper_36 = std::pow(
            helper_35 + 0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] + 0.33333333333333326 * _x[9],
            2);
    double helper_37 = -_x[7];
    double helper_38 = std::pow(helper_37 + 0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] +
                                0.33333333333333326 * _x[4], 2);
    double helper_39 = -_x[8];
    double helper_40 = std::pow(helper_39 + 0.33333333333333326 * _x[11] + 0.33333333333333348 * _x[2] +
                                0.33333333333333326 * _x[5], 2);
    double helper_41 = helper_26 + helper_27 + helper_28 + 1.3333333333333333 * helper_30 +
                       1.3333333333333333 * helper_32 + 1.3333333333333333 * helper_34 +
                       1.5000000000000002 * helper_36 + 1.5000000000000002 * helper_38 +
                       1.5000000000000002 * helper_40;
    double helper_42 = std::pow(helper_23, -1.6666666666666665);
    double helper_43 = 0.94280904158206336 * helper_10;
    double helper_44 = helper_31 + _x[1];
    double helper_45 = helper_17 * helper_44;
    double helper_46 = 0.94280904158206336 * helper_45;
    double helper_47 = 0.66666666666666663 * helper_21;
    double helper_48 = helper_16 * (helper_39 + _x[2]);
    double helper_49 = 0.66666666666666663 * helper_48;
    double helper_50 = 0.66666666666666663 * helper_22;
    double helper_51 = helper_11 * helper_44;
    double helper_52 = 0.66666666666666663 * helper_51;
    double helper_53 = helper_43 + helper_46 - helper_47 - helper_49 - helper_50 - helper_52;
    double helper_54 = helper_42 * helper_53;
    double helper_55 = helper_25 * helper_4 + helper_41 * helper_54;
    double helper_56 = 0.44444444444444431 * helper_24 / std::pow(
            0.66666666666666652 * helper_26 + 0.66666666666666652 * helper_27 +
            0.66666666666666652 * helper_28 + 0.88888888888888873 * helper_30 +
            0.88888888888888873 * helper_32 + 0.88888888888888873 * helper_34 + helper_36 + helper_38 +
            helper_40, 2);
    double helper_57 = helper_55 * helper_56;
    double helper_58 = 1.0 / helper_41;
    double helper_59 = std::pow(helper_23, -0.33333333333333337) * helper_58;
    double helper_60 = helper_55 * helper_59;
    double helper_61 = 3.0 * helper_25;
    double helper_62 = std::pow(helper_23, -2.6666666666666665) * helper_41;
    double helper_63 = helper_53 * helper_62;
    double helper_64 = helper_24 * helper_58;
    double helper_65 = 3.0 * _x[1];
    double helper_66 = 0.99999999999999978 * _x[10];
    double helper_67 = 1.0 * _x[4];
    double helper_68 = 1.0000000000000007 * _x[7];
    double helper_69 = -helper_65 + helper_66 + helper_67 + helper_68;
    double helper_70 = helper_33 + _x[2];
    double helper_71 = helper_14 * helper_70;
    double helper_72 = 0.94280904158206336 * helper_71;
    double helper_73 = 0.94280904158206336 * helper_19;
    double helper_74 = helper_5 * helper_70;
    double helper_75 = 0.66666666666666663 * helper_74;
    double helper_76 = 0.66666666666666663 * helper_20;
    double helper_77 = 0.66666666666666663 * helper_11;
    double helper_78 = helper_18 * helper_77;
    double helper_79 = helper_35 + _x[0];
    double helper_80 = helper_77 * helper_79;
    double helper_81 = -helper_72 - helper_73 + helper_75 + helper_76 + helper_78 + helper_80;
    double helper_82 = helper_72 + helper_73 - helper_75 - helper_76 - helper_78 - helper_80;
    double helper_83 = helper_4 * helper_42;
    double helper_84 = helper_65 - helper_66 - helper_67 - helper_68;
    double helper_85 = 1.6666666666666665 * helper_11;
    double helper_86 =
            -helper_18 * helper_85 + 2.3570226039551585 * helper_19 - 1.6666666666666665 * helper_20 +
            2.3570226039551585 * helper_71 - 1.6666666666666665 * helper_74 - helper_79 * helper_85;
    double helper_87 = 0.99999999999999978 * _x[11];
    double helper_88 = 3.0 * _x[2];
    double helper_89 = 1.0 * _x[5];
    double helper_90 = 1.0000000000000007 * _x[8];
    double helper_91 = helper_87 - helper_88 + helper_89 + helper_90;
    double helper_92 = 0.94280904158206336 * helper_15;
    double helper_93 = helper_37 + _x[1];
    double helper_94 = helper_18 * helper_93;
    double helper_95 = 0.94280904158206336 * helper_94;
    double helper_96 = 0.66666666666666663 * helper_5;
    double helper_97 = helper_12 * helper_96;
    double helper_98 = helper_93 * helper_96;
    double helper_99 = 0.66666666666666663 * helper_16;
    double helper_100 = helper_14 * helper_99;
    double helper_101 = helper_29 + _x[0];
    double helper_102 = helper_101 * helper_99;
    double helper_103 = helper_100 + helper_102 - helper_92 - helper_95 + helper_97 + helper_98;
    double helper_104 = -helper_100 - helper_102 + helper_92 + helper_95 - helper_97 - helper_98;
    double helper_105 = -helper_87 + helper_88 - helper_89 - helper_90;
    double helper_106 = 1.6666666666666665 * helper_5;
    double helper_107 = 1.6666666666666665 * helper_16;
    double helper_108 = -helper_101 * helper_107 - helper_106 * helper_12 - helper_106 * helper_93 -
                        helper_107 * helper_14 + 2.3570226039551585 * helper_15 +
                        2.3570226039551585 * helper_94;
    double helper_109 = helper_42 * helper_82;
    double helper_110 = helper_109 * helper_41 + helper_25 * helper_84;
    double helper_111 = helper_110 * helper_56;
    double helper_112 = helper_110 * helper_59;
    double helper_113 = helper_62 * helper_82;
    double helper_114 = helper_104 * helper_42;
    double helper_115 = helper_105 * helper_25 + helper_114 * helper_41;

    _H(0, 0) = helper_57 * (-helper_0 + helper_1 + helper_2 + helper_3) +
               helper_60 * (-helper_43 - helper_46 + helper_47 + helper_49 + helper_50 + helper_52) +
               helper_64 * (2 * helper_4 * helper_54 + helper_61 + helper_63 * (2.3570226039551585 * helper_10 -
                                                                                1.6666666666666665 * helper_21 -
                                                                                1.6666666666666665 * helper_22 +
                                                                                2.3570226039551585 * helper_45 -
                                                                                1.6666666666666665 * helper_48 -
                                                                                1.6666666666666665 *
                                                                                helper_51));
    _H(0, 1) = helper_57 * helper_69 + helper_60 * helper_81 +
               helper_64 * (helper_54 * helper_84 + helper_63 * helper_86 + helper_82 * helper_83);
    _H(0, 2) = helper_103 * helper_60 + helper_57 * helper_91 +
               helper_64 * (helper_104 * helper_83 + helper_105 * helper_54 + helper_108 * helper_63);
    _H(1, 1) = helper_111 * helper_69 + helper_112 * helper_81 +
               helper_64 * (2 * helper_109 * helper_84 + helper_113 * helper_86 + helper_61);
    _H(1, 2) = helper_103 * helper_112 + helper_111 * helper_91 +
               helper_64 * (helper_105 * helper_109 + helper_108 * helper_113 + helper_114 * helper_84);
    _H(2, 2) = helper_103 * helper_115 * helper_59 + helper_115 * helper_56 * helper_91 +
               helper_64 * (helper_104 * helper_108 * helper_62 + 2 * helper_105 * helper_114 + helper_61);


    for (int i = 1; i < 3; ++i)
      for (int j = 0; j < i; ++j)
        _H(i, j) = _H(j, i);

    // Compute Eigen decomposition H = V D V^T
    Eigen::SelfAdjointEigenSolver<Mat3d> solver(_H);
    const Mat3d &V = solver.eigenvectors();
    Vec3d evals = solver.eigenvalues();

    // check eigenvalues against epsilon, for those less that eps (zero)
    const int n = n_unknowns();
    for (int i = 0; i < n; ++i)
    {
      if (evals[i] < 1e-8)
      {
        evals[i] = 1e-8;
      }
    }

    // compute correction matrix M = V * diag(m) * V^T
    _H.noalias() = V * evals.asDiagonal() * V.transpose();
  }
};


//Two free vertices
class IMRMEnergy2
{
public:
  static int n_unknowns() { return 6; }

  static double eval_f(const double *_x)
  {
    //trace
    double helper_0 = -_x[0];
    double helper_1 = helper_0 + _x[3];
    double helper_2 = -_x[1];
    double helper_3 = helper_2 + _x[7];
    double helper_4 = -_x[2];
    double helper_5 = helper_4 + _x[11];
    double helper_6 = helper_0 + _x[6];
    double helper_7 = helper_2 + _x[10];
    double helper_8 = helper_4 + _x[5];
    double helper_9 = helper_0 + _x[9];
    double helper_10 = helper_2 + _x[4];
    double helper_11 = helper_4 + _x[8];

    return std::log(std::pow(-helper_1 * helper_11 * helper_7 + helper_1 * helper_3 * helper_5 +
                             helper_10 * helper_11 * helper_9 - helper_10 * helper_5 * helper_6 -
                             helper_3 * helper_8 * helper_9 + helper_6 * helper_7 * helper_8,
                             -0.66666666666666663) *
                    (std::pow(_x[0] - _x[3], 2) + std::pow(_x[1] - _x[4], 2) + std::pow(_x[2] - _x[5], 2) +
                     1.3333333333333333 * std::pow(0.5 * _x[0] + 0.5 * _x[3] - _x[9], 2) +
                     1.3333333333333333 * std::pow(0.5 * _x[1] - _x[10] + 0.5 * _x[4], 2) +
                     1.3333333333333333 * std::pow(-_x[11] + 0.5 * _x[2] + 0.5 * _x[5], 2) +
                     1.5000000000000002 * std::pow(
                             0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] - _x[6] +
                             0.33333333333333326 * _x[9], 2) + 1.5000000000000002 * std::pow(
                            0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] +
                            0.33333333333333326 * _x[4] - _x[7], 2) + 1.5000000000000002 * std::pow(
                            0.33333333333333326 * _x[11] + 0.33333333333333348 * _x[2] +
                            0.33333333333333326 * _x[5] - _x[8], 2))) - 0.2310490601866485;
  }

  static void eval_gradient(const double *_x, double *_g)
  {
    //trace
    double helper_0 = -1.4142135623730951 * _x[0] + 1.4142135623730951 * _x[3];
    double helper_1 = -_x[2];
    double helper_2 = helper_1 + _x[11];
    double helper_3 = -_x[1];
    double helper_4 = helper_3 + _x[7];
    double helper_5 = helper_2 * helper_4;
    double helper_6 = -1.4142135623730951 * _x[2] + 1.4142135623730951 * _x[5];
    double helper_7 = helper_3 + _x[10];
    double helper_8 = -_x[0];
    double helper_9 = helper_8 + _x[6];
    double helper_10 = helper_7 * helper_9;
    double helper_11 = -1.4142135623730951 * _x[1] + 1.4142135623730951 * _x[4];
    double helper_12 = helper_1 + _x[8];
    double helper_13 = helper_8 + _x[9];
    double helper_14 = helper_12 * helper_13;
    double helper_15 = helper_0 * helper_12;
    double helper_16 = helper_11 * helper_2;
    double helper_17 = helper_4 * helper_6;
    double helper_18 = helper_0 * helper_5 + helper_10 * helper_6 + helper_11 * helper_14 - helper_13 * helper_17 -
                       helper_15 * helper_7 - helper_16 * helper_9;
    double helper_19 = std::pow(helper_18, 0.66666666666666663);
    double helper_20 = 1.0 / helper_19;
    double helper_21 = 0.94280904158206336 * helper_5;
    double helper_22 = -_x[10];
    double helper_23 = helper_22 + _x[1];
    double helper_24 = 0.94280904158206336 * helper_12 * helper_23;
    double helper_25 = -_x[8];
    double helper_26 = 0.66666666666666663 * helper_11;
    double helper_27 = 0.66666666666666663 * helper_6;
    double helper_28 = -_x[9];
    double helper_29 = -_x[11];
    double helper_30 = -_x[6];
    double helper_31 = -_x[7];
    double helper_32 = std::pow(helper_1 + _x[5], 2) + std::pow(helper_3 + _x[4], 2) + std::pow(helper_8 + _x[3], 2) +
                       1.3333333333333333 * std::pow(helper_22 + 0.5 * _x[1] + 0.5 * _x[4], 2) +
                       1.3333333333333333 * std::pow(helper_28 + 0.5 * _x[0] + 0.5 * _x[3], 2) +
                       1.3333333333333333 * std::pow(helper_29 + 0.5 * _x[2] + 0.5 * _x[5], 2) + 1.5000000000000002 *
                                                                                                 std::pow(helper_25 +
                                                                                                          0.33333333333333326 *
                                                                                                          _x[11] +
                                                                                                          0.33333333333333348 *
                                                                                                          _x[2] +
                                                                                                          0.33333333333333326 *
                                                                                                          _x[5], 2) +
                       1.5000000000000002 * std::pow(
                               helper_30 + 0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] +
                               0.33333333333333326 * _x[9], 2) + 1.5000000000000002 * std::pow(
            helper_31 + 0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] + 0.33333333333333326 * _x[4], 2);
    double helper_33 = std::pow(helper_18, -1.6666666666666665) * helper_32;
    double helper_34 = helper_19 / helper_32;
    double helper_35 = helper_29 + _x[2];
    double helper_36 = 0.94280904158206336 * helper_35 * helper_9;
    double helper_37 = 0.94280904158206336 * helper_14;
    double helper_38 = 0.66666666666666663 * helper_0;
    double helper_39 = 0.94280904158206336 * helper_10;
    double helper_40 = helper_31 + _x[1];
    double helper_41 = 0.94280904158206336 * helper_13 * helper_40;

    _g[0] = helper_34 *
            (helper_20 * (3.0 * _x[0] - 1.0 * _x[3] - 1.0000000000000007 * _x[6] - 0.99999999999999978 * _x[9]) +
             helper_33 *
             (-0.66666666666666663 * helper_16 - 0.66666666666666663 * helper_17 + helper_21 - helper_23 * helper_27 +
              helper_24 - helper_26 * (helper_25 + _x[2])));
    _g[1] = helper_34 *
            (helper_20 * (3.0 * _x[1] - 0.99999999999999978 * _x[10] - 1.0 * _x[4] - 1.0000000000000007 * _x[7]) +
             helper_33 * (-helper_13 * helper_27 - 0.66666666666666663 * helper_15 - helper_27 * (helper_30 + _x[0]) -
                          helper_35 * helper_38 + helper_36 + helper_37));
    _g[2] = helper_34 *
            (helper_20 * (-0.99999999999999978 * _x[11] + 3.0 * _x[2] - 1.0 * _x[5] - 1.0000000000000007 * _x[8]) +
             helper_33 *
             (-helper_26 * helper_9 - helper_26 * (helper_28 + _x[0]) - helper_38 * helper_40 - helper_38 * helper_7 +
              helper_39 + helper_41));
    _g[3] = helper_34 * (helper_20 * (-1.0 * _x[0] + 3.0 * _x[3] - 0.99999999999999989 * _x[6] - 1.0 * _x[9]) +
                         helper_33 * (-helper_21 - helper_24));
    _g[4] = helper_34 * (helper_20 * (-1.0 * _x[1] - 1.0 * _x[10] + 3.0 * _x[4] - 0.99999999999999989 * _x[7]) +
                         helper_33 * (-helper_36 - helper_37));
    _g[5] = helper_34 * (helper_20 * (-1.0 * _x[11] - 1.0 * _x[2] + 3.0 * _x[5] - 0.99999999999999989 * _x[8]) +
                         helper_33 * (-helper_39 - helper_41));
  }

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H)
  {
    _H.setZero();

    //trace
    double helper_0 = 3.0 * _x[0];
    double helper_1 = 1.0 * _x[3];
    double helper_2 = 1.0000000000000007 * _x[6];
    double helper_3 = 0.99999999999999978 * _x[9];
    double helper_4 = helper_0 - helper_1 - helper_2 - helper_3;
    double helper_5 = -1.4142135623730951 * _x[0] + 1.4142135623730951 * _x[3];
    double helper_6 = -_x[2];
    double helper_7 = helper_6 + _x[11];
    double helper_8 = -_x[1];
    double helper_9 = helper_8 + _x[7];
    double helper_10 = helper_7 * helper_9;
    double helper_11 = -1.4142135623730951 * _x[2] + 1.4142135623730951 * _x[5];
    double helper_12 = helper_8 + _x[10];
    double helper_13 = -_x[0];
    double helper_14 = helper_13 + _x[6];
    double helper_15 = helper_12 * helper_14;
    double helper_16 = -1.4142135623730951 * _x[1] + 1.4142135623730951 * _x[4];
    double helper_17 = helper_6 + _x[8];
    double helper_18 = helper_13 + _x[9];
    double helper_19 = helper_17 * helper_18;
    double helper_20 = helper_17 * helper_5;
    double helper_21 = helper_16 * helper_7;
    double helper_22 = helper_11 * helper_9;
    double helper_23 = helper_10 * helper_5 + helper_11 * helper_15 - helper_12 * helper_20 - helper_14 * helper_21 +
                       helper_16 * helper_19 - helper_18 * helper_22;
    double helper_24 = std::pow(helper_23, 0.66666666666666663);
    double helper_25 = 1.0 / helper_24;
    double helper_26 = std::pow(helper_13 + _x[3], 2);
    double helper_27 = std::pow(helper_8 + _x[4], 2);
    double helper_28 = std::pow(helper_6 + _x[5], 2);
    double helper_29 = -_x[9];
    double helper_30 = std::pow(helper_29 + 0.5 * _x[0] + 0.5 * _x[3], 2);
    double helper_31 = -_x[10];
    double helper_32 = std::pow(helper_31 + 0.5 * _x[1] + 0.5 * _x[4], 2);
    double helper_33 = -_x[11];
    double helper_34 = std::pow(helper_33 + 0.5 * _x[2] + 0.5 * _x[5], 2);
    double helper_35 = -_x[6];
    double helper_36 = std::pow(
            helper_35 + 0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] + 0.33333333333333326 * _x[9], 2);
    double helper_37 = -_x[7];
    double helper_38 = std::pow(
            helper_37 + 0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] + 0.33333333333333326 * _x[4], 2);
    double helper_39 = -_x[8];
    double helper_40 = std::pow(
            helper_39 + 0.33333333333333326 * _x[11] + 0.33333333333333348 * _x[2] + 0.33333333333333326 * _x[5], 2);
    double helper_41 =
            helper_26 + helper_27 + helper_28 + 1.3333333333333333 * helper_30 + 1.3333333333333333 * helper_32 +
            1.3333333333333333 * helper_34 + 1.5000000000000002 * helper_36 + 1.5000000000000002 * helper_38 +
            1.5000000000000002 * helper_40;
    double helper_42 = std::pow(helper_23, -1.6666666666666665);
    double helper_43 = helper_31 + _x[1];
    double helper_44 = helper_11 * helper_43;
    double helper_45 = 0.66666666666666663 * helper_44;
    double helper_46 = 0.66666666666666663 * helper_21;
    double helper_47 = helper_16 * (helper_39 + _x[2]);
    double helper_48 = 0.66666666666666663 * helper_47;
    double helper_49 = 0.66666666666666663 * helper_22;
    double helper_50 = 0.94280904158206336 * helper_10;
    double helper_51 = helper_17 * helper_43;
    double helper_52 = 0.94280904158206336 * helper_51;
    double helper_53 = helper_50 + helper_52;
    double helper_54 = -helper_45 - helper_46 - helper_48 - helper_49 + helper_53;
    double helper_55 = helper_42 * helper_54;
    double helper_56 = helper_25 * helper_4 + helper_41 * helper_55;
    double helper_57 = 0.44444444444444431 * helper_24 / std::pow(
            0.66666666666666652 * helper_26 + 0.66666666666666652 * helper_27 + 0.66666666666666652 * helper_28 +
            0.88888888888888873 * helper_30 + 0.88888888888888873 * helper_32 + 0.88888888888888873 * helper_34 +
            helper_36 + helper_38 + helper_40, 2);
    double helper_58 = helper_56 * helper_57;
    double helper_59 = -helper_50 - helper_52;
    double helper_60 = 1.0 / helper_41;
    double helper_61 = std::pow(helper_23, -0.33333333333333337) * helper_60;
    double helper_62 = helper_56 * helper_61;
    double helper_63 = 3.0 * helper_25;
    double helper_64 = 2.3570226039551585 * helper_10;
    double helper_65 = 2.3570226039551585 * helper_51;
    double helper_66 = std::pow(helper_23, -2.6666666666666665) * helper_41;
    double helper_67 = helper_54 * helper_66;
    double helper_68 = helper_24 * helper_60;
    double helper_69 = 3.0 * _x[1];
    double helper_70 = 0.99999999999999978 * _x[10];
    double helper_71 = 1.0 * _x[4];
    double helper_72 = 1.0000000000000007 * _x[7];
    double helper_73 = -helper_69 + helper_70 + helper_71 + helper_72;
    double helper_74 = helper_35 + _x[0];
    double helper_75 = 0.66666666666666663 * helper_11;
    double helper_76 = helper_74 * helper_75;
    double helper_77 = helper_33 + _x[2];
    double helper_78 = helper_5 * helper_77;
    double helper_79 = 0.66666666666666663 * helper_78;
    double helper_80 = 0.66666666666666663 * helper_20;
    double helper_81 = helper_18 * helper_75;
    double helper_82 = helper_14 * helper_77;
    double helper_83 = 0.94280904158206336 * helper_82;
    double helper_84 = 0.94280904158206336 * helper_19;
    double helper_85 = -helper_83 - helper_84;
    double helper_86 = helper_76 + helper_79 + helper_80 + helper_81 + helper_85;
    double helper_87 = helper_83 + helper_84;
    double helper_88 = -helper_76 - helper_79 - helper_80 - helper_81 + helper_87;
    double helper_89 = helper_4 * helper_42;
    double helper_90 = helper_69 - helper_70 - helper_71 - helper_72;
    double helper_91 = 2.3570226039551585 * helper_82;
    double helper_92 = 2.3570226039551585 * helper_19;
    double helper_93 = 1.6666666666666665 * helper_11;
    double helper_94 = -helper_18 * helper_93 - 1.6666666666666665 * helper_20 - helper_74 * helper_93 -
                       1.6666666666666665 * helper_78 + helper_91 + helper_92;
    double helper_95 = 0.99999999999999978 * _x[11];
    double helper_96 = 3.0 * _x[2];
    double helper_97 = 1.0 * _x[5];
    double helper_98 = 1.0000000000000007 * _x[8];
    double helper_99 = helper_95 - helper_96 + helper_97 + helper_98;
    double helper_100 = helper_29 + _x[0];
    double helper_101 = 0.66666666666666663 * helper_16;
    double helper_102 = helper_100 * helper_101;
    double helper_103 = helper_37 + _x[1];
    double helper_104 = 0.66666666666666663 * helper_5;
    double helper_105 = helper_103 * helper_104;
    double helper_106 = helper_104 * helper_12;
    double helper_107 = helper_101 * helper_14;
    double helper_108 = 0.94280904158206336 * helper_15;
    double helper_109 = helper_103 * helper_18;
    double helper_110 = 0.94280904158206336 * helper_109;
    double helper_111 = -helper_108 - helper_110;
    double helper_112 = helper_102 + helper_105 + helper_106 + helper_107 + helper_111;
    double helper_113 = helper_108 + helper_110;
    double helper_114 = -helper_102 - helper_105 - helper_106 - helper_107 + helper_113;
    double helper_115 = -helper_95 + helper_96 - helper_97 - helper_98;
    double helper_116 = 2.3570226039551585 * helper_15;
    double helper_117 = 2.3570226039551585 * helper_109;
    double helper_118 = 1.6666666666666665 * helper_5;
    double helper_119 = 1.6666666666666665 * helper_16;
    double helper_120 =
            -helper_100 * helper_119 - helper_103 * helper_118 + helper_116 + helper_117 - helper_118 * helper_12 -
            helper_119 * helper_14;
    double helper_121 = 1.0 * _x[0];
    double helper_122 = 3.0 * _x[3];
    double helper_123 = 0.99999999999999989 * _x[6];
    double helper_124 = 1.0 * _x[9];
    double helper_125 = helper_121 - helper_122 + helper_123 + helper_124;
    double helper_126 = -1.0 * helper_25;
    double helper_127 = -helper_121 + helper_122 - helper_123 - helper_124;
    double helper_128 = -helper_64 - helper_65;
    double helper_129 = 1.0 * _x[1];
    double helper_130 = 1.0 * _x[10];
    double helper_131 = 3.0 * _x[4];
    double helper_132 = 0.99999999999999989 * _x[7];
    double helper_133 = helper_129 + helper_130 - helper_131 + helper_132;
    double helper_134 = -helper_129 - helper_130 + helper_131 - helper_132;
    double helper_135 = 0.94280904158206336 * _x[11];
    double helper_136 = 0.94280904158206336 * _x[8];
    double helper_137 = helper_41 * helper_42;
    double helper_138 = -helper_91 - helper_92;
    double helper_139 = 1.0 * _x[11];
    double helper_140 = 1.0 * _x[2];
    double helper_141 = 3.0 * _x[5];
    double helper_142 = 0.99999999999999989 * _x[8];
    double helper_143 = helper_139 + helper_140 - helper_141 + helper_142;
    double helper_144 = -helper_139 - helper_140 + helper_141 - helper_142;
    double helper_145 = 0.94280904158206336 * _x[10];
    double helper_146 = 0.94280904158206336 * _x[7];
    double helper_147 = -helper_116 - helper_117;
    double helper_148 = helper_137 * helper_88 + helper_25 * helper_90;
    double helper_149 = helper_148 * helper_57;
    double helper_150 = helper_148 * helper_61;
    double helper_151 = helper_42 * helper_90;
    double helper_152 = helper_66 * helper_88;
    double helper_153 = helper_42 * helper_88;
    double helper_154 = 0.94280904158206336 * _x[6];
    double helper_155 = 0.94280904158206336 * _x[9];
    double helper_156 = helper_114 * helper_137 + helper_115 * helper_25;
    double helper_157 = helper_156 * helper_57;
    double helper_158 = helper_156 * helper_61;
    double helper_159 = helper_115 * helper_42;
    double helper_160 = helper_114 * helper_66;
    double helper_161 = helper_114 * helper_42;
    double helper_162 = helper_127 * helper_25 + helper_137 * helper_59;
    double helper_163 = helper_162 * helper_57;
    double helper_164 = helper_162 * helper_61;
    double helper_165 = helper_42 * helper_59;
    double helper_166 = helper_59 * helper_66;
    double helper_167 = helper_127 * helper_42;
    double helper_168 = helper_134 * helper_25 + helper_137 * helper_85;
    double helper_169 = helper_168 * helper_57;
    double helper_170 = helper_168 * helper_61;
    double helper_171 = helper_134 * helper_42;
    double helper_172 = helper_66 * helper_85;
    double helper_173 = helper_144 * helper_42;
    double helper_174 = helper_111 * helper_137 + helper_144 * helper_25;

    _H(0, 0) = helper_58 * (-helper_0 + helper_1 + helper_2 + helper_3) +
               helper_62 * (helper_45 + helper_46 + helper_48 + helper_49 + helper_59) + helper_68 *
                                                                                         (2 * helper_4 * helper_55 +
                                                                                          helper_63 + helper_67 *
                                                                                                      (-1.6666666666666665 *
                                                                                                       helper_21 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_22 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_44 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_47 +
                                                                                                       helper_64 +
                                                                                                       helper_65));
    _H(0, 1) = helper_58 * helper_73 + helper_62 * helper_86 +
               helper_68 * (helper_55 * helper_90 + helper_67 * helper_94 + helper_88 * helper_89);
    _H(0, 2) = helper_112 * helper_62 + helper_58 * helper_99 +
               helper_68 * (helper_114 * helper_89 + helper_115 * helper_55 + helper_120 * helper_67);
    _H(0, 3) = helper_125 * helper_58 + helper_53 * helper_62 +
               helper_68 * (helper_126 + helper_127 * helper_55 + helper_128 * helper_67 + helper_59 * helper_89);
    _H(0, 4) = helper_133 * helper_58 + helper_62 * helper_87 + helper_68 * (helper_134 * helper_55 +
                                                                             helper_137 * (-helper_135 + helper_136) +
                                                                             helper_138 * helper_67 +
                                                                             helper_85 * helper_89);
    _H(0, 5) = helper_113 * helper_62 + helper_143 * helper_58 + helper_68 * (helper_111 * helper_89 +
                                                                              helper_137 * (helper_145 - helper_146) +
                                                                              helper_144 * helper_55 +
                                                                              helper_147 * helper_67);
    _H(1, 1) = helper_149 * helper_73 + helper_150 * helper_86 +
               helper_68 * (2 * helper_151 * helper_88 + helper_152 * helper_94 + helper_63);
    _H(1, 2) = helper_112 * helper_150 + helper_149 * helper_99 +
               helper_68 * (helper_114 * helper_151 + helper_115 * helper_153 + helper_120 * helper_152);
    _H(1, 3) = helper_125 * helper_149 + helper_150 * helper_53 + helper_68 *
                                                                  (helper_127 * helper_153 + helper_128 * helper_152 +
                                                                   helper_137 * (helper_135 - helper_136) +
                                                                   helper_151 * helper_59);
    _H(1, 4) = helper_133 * helper_149 + helper_150 * helper_87 +
               helper_68 * (helper_126 + helper_134 * helper_153 + helper_138 * helper_152 + helper_151 * helper_85);
    _H(1, 5) = helper_113 * helper_150 + helper_143 * helper_149 + helper_68 * (helper_111 * helper_151 +
                                                                                helper_137 * (helper_154 - helper_155) +
                                                                                helper_144 * helper_153 +
                                                                                helper_147 * helper_152);
    _H(2, 2) = helper_112 * helper_158 + helper_157 * helper_99 +
               helper_68 * (2 * helper_114 * helper_159 + helper_120 * helper_160 + helper_63);
    _H(2, 3) = helper_125 * helper_157 + helper_158 * helper_53 + helper_68 *
                                                                  (helper_127 * helper_161 + helper_128 * helper_160 +
                                                                   helper_137 * (-helper_145 + helper_146) +
                                                                   helper_159 * helper_59);
    _H(2, 4) = helper_133 * helper_157 + helper_158 * helper_87 + helper_68 * (helper_134 * helper_161 +
                                                                               helper_137 * (-helper_154 + helper_155) +
                                                                               helper_138 * helper_160 +
                                                                               helper_159 * helper_85);
    _H(2, 5) = helper_113 * helper_158 + helper_143 * helper_157 +
               helper_68 * (helper_111 * helper_159 + helper_126 + helper_144 * helper_161 + helper_147 * helper_160);
    _H(3, 3) = helper_125 * helper_163 + helper_164 * helper_53 +
               helper_68 * (2 * helper_127 * helper_165 + helper_128 * helper_166 + helper_63);
    _H(3, 4) = helper_133 * helper_163 + helper_164 * helper_87 +
               helper_68 * (helper_134 * helper_165 + helper_138 * helper_166 + helper_167 * helper_85);
    _H(3, 5) = helper_113 * helper_164 + helper_143 * helper_163 +
               helper_68 * (helper_111 * helper_167 + helper_144 * helper_165 + helper_147 * helper_166);
    _H(4, 4) = helper_133 * helper_169 + helper_170 * helper_87 +
               helper_68 * (helper_138 * helper_172 + 2 * helper_171 * helper_85 + helper_63);
    _H(4, 5) = helper_113 * helper_170 + helper_143 * helper_169 +
               helper_68 * (helper_111 * helper_171 + helper_147 * helper_172 + helper_173 * helper_85);
    _H(5, 5) = helper_113 * helper_174 * helper_61 + helper_143 * helper_174 * helper_57 +
               helper_68 * (helper_111 * helper_147 * helper_66 + 2 * helper_111 * helper_173 + helper_63);


    for (int i = 1; i < 6; ++i)
      for (int j = 0; j < i; ++j)
        _H(i, j) = _H(j, i);

    // Compute Eigen decomposition H = V D V^T
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_H);
    const Eigen::MatrixXd &V = solver.eigenvectors();
    Eigen::VectorXd evals = solver.eigenvalues();

    // check eigenvalues against epsilon, for those less that eps (zero)
    const int n = n_unknowns();
    for (int i = 0; i < n; ++i)
    {
      if (evals[i] < 1e-8)
      {
        evals[i] = 1e-8;
      }
    }

    // compute correction matrix M = V * diag(m) * V^T
    _H.noalias() = V * evals.asDiagonal() * V.transpose();
  }
};


//Three free vertices
class IMRMEnergy3
{
public:
  static int n_unknowns() { return 9; }

  static double eval_f(const double *_x)
  {
    //trace
    double helper_0 = -_x[0];
    double helper_1 = helper_0 + _x[3];
    double helper_2 = -_x[1];
    double helper_3 = helper_2 + _x[7];
    double helper_4 = -_x[2];
    double helper_5 = helper_4 + _x[11];
    double helper_6 = helper_0 + _x[6];
    double helper_7 = helper_2 + _x[10];
    double helper_8 = helper_4 + _x[5];
    double helper_9 = helper_0 + _x[9];
    double helper_10 = helper_2 + _x[4];
    double helper_11 = helper_4 + _x[8];

    return std::log(std::pow(-helper_1 * helper_11 * helper_7 + helper_1 * helper_3 * helper_5 +
                             helper_10 * helper_11 * helper_9 - helper_10 * helper_5 * helper_6 -
                             helper_3 * helper_8 * helper_9 + helper_6 * helper_7 * helper_8,
                             -0.66666666666666663) *
                    (std::pow(_x[0] - _x[3], 2) + std::pow(_x[1] - _x[4], 2) + std::pow(_x[2] - _x[5], 2) +
                     1.3333333333333333 * std::pow(0.5 * _x[0] + 0.5 * _x[3] - _x[9], 2) +
                     1.3333333333333333 * std::pow(0.5 * _x[1] - _x[10] + 0.5 * _x[4], 2) +
                     1.3333333333333333 * std::pow(-_x[11] + 0.5 * _x[2] + 0.5 * _x[5], 2) +
                     1.5000000000000002 * std::pow(
                             0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] - _x[6] +
                             0.33333333333333326 * _x[9], 2) + 1.5000000000000002 * std::pow(
                            0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] +
                            0.33333333333333326 * _x[4] - _x[7], 2) + 1.5000000000000002 * std::pow(
                            0.33333333333333326 * _x[11] + 0.33333333333333348 * _x[2] +
                            0.33333333333333326 * _x[5] - _x[8], 2))) - 0.2310490601866485;
  }

  static void eval_gradient(const double *_x, double *_g)
  {
    //trace
    double helper_0 = -1.4142135623730951 * _x[0] + 1.4142135623730951 * _x[3];
    double helper_1 = -_x[2];
    double helper_2 = helper_1 + _x[11];
    double helper_3 = -_x[1];
    double helper_4 = helper_3 + _x[7];
    double helper_5 = helper_2 * helper_4;
    double helper_6 = -1.4142135623730951 * _x[2] + 1.4142135623730951 * _x[5];
    double helper_7 = helper_3 + _x[10];
    double helper_8 = -_x[0];
    double helper_9 = helper_8 + _x[6];
    double helper_10 = helper_7 * helper_9;
    double helper_11 = -1.4142135623730951 * _x[1] + 1.4142135623730951 * _x[4];
    double helper_12 = helper_1 + _x[8];
    double helper_13 = helper_8 + _x[9];
    double helper_14 = helper_12 * helper_13;
    double helper_15 = helper_0 * helper_12;
    double helper_16 = helper_11 * helper_2;
    double helper_17 = helper_4 * helper_6;
    double helper_18 = helper_0 * helper_5 + helper_10 * helper_6 + helper_11 * helper_14 - helper_13 * helper_17 -
                       helper_15 * helper_7 - helper_16 * helper_9;
    double helper_19 = std::pow(helper_18, 0.66666666666666663);
    double helper_20 = 1.0 / helper_19;
    double helper_21 = 0.94280904158206336 * helper_5;
    double helper_22 = -_x[10];
    double helper_23 = helper_22 + _x[1];
    double helper_24 = 0.94280904158206336 * helper_12 * helper_23;
    double helper_25 = -_x[8];
    double helper_26 = 0.66666666666666663 * helper_11;
    double helper_27 = 0.66666666666666663 * helper_6;
    double helper_28 = -_x[9];
    double helper_29 = -_x[11];
    double helper_30 = -_x[6];
    double helper_31 = -_x[7];
    double helper_32 = std::pow(helper_1 + _x[5], 2) + std::pow(helper_3 + _x[4], 2) + std::pow(helper_8 + _x[3], 2) +
                       1.3333333333333333 * std::pow(helper_22 + 0.5 * _x[1] + 0.5 * _x[4], 2) +
                       1.3333333333333333 * std::pow(helper_28 + 0.5 * _x[0] + 0.5 * _x[3], 2) +
                       1.3333333333333333 * std::pow(helper_29 + 0.5 * _x[2] + 0.5 * _x[5], 2) + 1.5000000000000002 *
                                                                                                 std::pow(helper_25 +
                                                                                                          0.33333333333333326 *
                                                                                                          _x[11] +
                                                                                                          0.33333333333333348 *
                                                                                                          _x[2] +
                                                                                                          0.33333333333333326 *
                                                                                                          _x[5], 2) +
                       1.5000000000000002 * std::pow(
                               helper_30 + 0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] +
                               0.33333333333333326 * _x[9], 2) + 1.5000000000000002 * std::pow(
            helper_31 + 0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] + 0.33333333333333326 * _x[4], 2);
    double helper_33 = std::pow(helper_18, -1.6666666666666665) * helper_32;
    double helper_34 = helper_19 / helper_32;
    double helper_35 = helper_29 + _x[2];
    double helper_36 = 0.94280904158206336 * helper_35 * helper_9;
    double helper_37 = 0.94280904158206336 * helper_14;
    double helper_38 = 0.66666666666666663 * helper_0;
    double helper_39 = 0.94280904158206336 * helper_10;
    double helper_40 = helper_31 + _x[1];
    double helper_41 = 0.94280904158206336 * helper_13 * helper_40;
    double helper_42 = helper_28 + _x[0];

    _g[0] = helper_34 *
            (helper_20 * (3.0 * _x[0] - 1.0 * _x[3] - 1.0000000000000007 * _x[6] - 0.99999999999999978 * _x[9]) +
             helper_33 *
             (-0.66666666666666663 * helper_16 - 0.66666666666666663 * helper_17 + helper_21 - helper_23 * helper_27 +
              helper_24 - helper_26 * (helper_25 + _x[2])));
    _g[1] = helper_34 *
            (helper_20 * (3.0 * _x[1] - 0.99999999999999978 * _x[10] - 1.0 * _x[4] - 1.0000000000000007 * _x[7]) +
             helper_33 * (-helper_13 * helper_27 - 0.66666666666666663 * helper_15 - helper_27 * (helper_30 + _x[0]) -
                          helper_35 * helper_38 + helper_36 + helper_37));
    _g[2] = helper_34 *
            (helper_20 * (-0.99999999999999978 * _x[11] + 3.0 * _x[2] - 1.0 * _x[5] - 1.0000000000000007 * _x[8]) +
             helper_33 *
             (-helper_26 * helper_42 - helper_26 * helper_9 - helper_38 * helper_40 - helper_38 * helper_7 + helper_39 +
              helper_41));
    _g[3] = helper_34 * (helper_20 * (-1.0 * _x[0] + 3.0 * _x[3] - 0.99999999999999989 * _x[6] - 1.0 * _x[9]) +
                         helper_33 * (-helper_21 - helper_24));
    _g[4] = helper_34 * (helper_20 * (-1.0 * _x[1] - 1.0 * _x[10] + 3.0 * _x[4] - 0.99999999999999989 * _x[7]) +
                         helper_33 * (-helper_36 - helper_37));
    _g[5] = helper_34 * (helper_20 * (-1.0 * _x[11] - 1.0 * _x[2] + 3.0 * _x[5] - 0.99999999999999989 * _x[8]) +
                         helper_33 * (-helper_39 - helper_41));
    _g[6] = helper_34 * (helper_20 *
                         (-1.0000000000000007 * _x[0] - 0.99999999999999989 * _x[3] + 3.0000000000000004 * _x[6] -
                          0.99999999999999989 * _x[9]) + helper_33 * (-helper_26 * helper_35 - helper_27 * helper_7));
    _g[7] = helper_34 * (helper_20 *
                         (-1.0000000000000007 * _x[1] - 0.99999999999999989 * _x[10] - 0.99999999999999989 * _x[4] +
                          3.0000000000000004 * _x[7]) + helper_33 * (-helper_2 * helper_38 - helper_27 * helper_42));
    _g[8] = helper_34 * (helper_20 *
                         (-0.99999999999999989 * _x[11] - 1.0000000000000007 * _x[2] - 0.99999999999999989 * _x[5] +
                          3.0000000000000004 * _x[8]) + helper_33 * (-helper_13 * helper_26 - helper_23 * helper_38));
  }

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H)
  {
    _H.setZero();


    //trace
    double helper_0 = 3.0 * _x[0];
    double helper_1 = 1.0 * _x[3];
    double helper_2 = 1.0000000000000007 * _x[6];
    double helper_3 = 0.99999999999999978 * _x[9];
    double helper_4 = helper_0 - helper_1 - helper_2 - helper_3;
    double helper_5 = -1.4142135623730951 * _x[0] + 1.4142135623730951 * _x[3];
    double helper_6 = -_x[2];
    double helper_7 = helper_6 + _x[11];
    double helper_8 = -_x[1];
    double helper_9 = helper_8 + _x[7];
    double helper_10 = helper_7 * helper_9;
    double helper_11 = -1.4142135623730951 * _x[2] + 1.4142135623730951 * _x[5];
    double helper_12 = helper_8 + _x[10];
    double helper_13 = -_x[0];
    double helper_14 = helper_13 + _x[6];
    double helper_15 = helper_12 * helper_14;
    double helper_16 = -1.4142135623730951 * _x[1] + 1.4142135623730951 * _x[4];
    double helper_17 = helper_6 + _x[8];
    double helper_18 = helper_13 + _x[9];
    double helper_19 = helper_17 * helper_18;
    double helper_20 = helper_17 * helper_5;
    double helper_21 = helper_16 * helper_7;
    double helper_22 = helper_11 * helper_9;
    double helper_23 = helper_10 * helper_5 + helper_11 * helper_15 - helper_12 * helper_20 - helper_14 * helper_21 +
                       helper_16 * helper_19 - helper_18 * helper_22;
    double helper_24 = std::pow(helper_23, 0.66666666666666663);
    double helper_25 = 1.0 / helper_24;
    double helper_26 = std::pow(helper_13 + _x[3], 2);
    double helper_27 = std::pow(helper_8 + _x[4], 2);
    double helper_28 = std::pow(helper_6 + _x[5], 2);
    double helper_29 = -_x[9];
    double helper_30 = std::pow(helper_29 + 0.5 * _x[0] + 0.5 * _x[3], 2);
    double helper_31 = -_x[10];
    double helper_32 = std::pow(helper_31 + 0.5 * _x[1] + 0.5 * _x[4], 2);
    double helper_33 = -_x[11];
    double helper_34 = std::pow(helper_33 + 0.5 * _x[2] + 0.5 * _x[5], 2);
    double helper_35 = -_x[6];
    double helper_36 = std::pow(
            helper_35 + 0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] + 0.33333333333333326 * _x[9], 2);
    double helper_37 = -_x[7];
    double helper_38 = std::pow(
            helper_37 + 0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] + 0.33333333333333326 * _x[4], 2);
    double helper_39 = -_x[8];
    double helper_40 = std::pow(
            helper_39 + 0.33333333333333326 * _x[11] + 0.33333333333333348 * _x[2] + 0.33333333333333326 * _x[5], 2);
    double helper_41 =
            helper_26 + helper_27 + helper_28 + 1.3333333333333333 * helper_30 + 1.3333333333333333 * helper_32 +
            1.3333333333333333 * helper_34 + 1.5000000000000002 * helper_36 + 1.5000000000000002 * helper_38 +
            1.5000000000000002 * helper_40;
    double helper_42 = std::pow(helper_23, -1.6666666666666665);
    double helper_43 = helper_31 + _x[1];
    double helper_44 = helper_11 * helper_43;
    double helper_45 = 0.66666666666666663 * helper_44;
    double helper_46 = 0.66666666666666663 * helper_21;
    double helper_47 = helper_16 * (helper_39 + _x[2]);
    double helper_48 = 0.66666666666666663 * helper_47;
    double helper_49 = 0.66666666666666663 * helper_22;
    double helper_50 = 0.94280904158206336 * helper_10;
    double helper_51 = helper_17 * helper_43;
    double helper_52 = 0.94280904158206336 * helper_51;
    double helper_53 = helper_50 + helper_52;
    double helper_54 = -helper_45 - helper_46 - helper_48 - helper_49 + helper_53;
    double helper_55 = helper_42 * helper_54;
    double helper_56 = helper_25 * helper_4 + helper_41 * helper_55;
    double helper_57 = 0.44444444444444431 * helper_24 / std::pow(
            0.66666666666666652 * helper_26 + 0.66666666666666652 * helper_27 + 0.66666666666666652 * helper_28 +
            0.88888888888888873 * helper_30 + 0.88888888888888873 * helper_32 + 0.88888888888888873 * helper_34 +
            helper_36 + helper_38 + helper_40, 2);
    double helper_58 = helper_56 * helper_57;
    double helper_59 = -helper_50 - helper_52;
    double helper_60 = 1.0 / helper_41;
    double helper_61 = std::pow(helper_23, -0.33333333333333337) * helper_60;
    double helper_62 = helper_56 * helper_61;
    double helper_63 = 3.0 * helper_25;
    double helper_64 = 2.3570226039551585 * helper_10;
    double helper_65 = 2.3570226039551585 * helper_51;
    double helper_66 = std::pow(helper_23, -2.6666666666666665) * helper_41;
    double helper_67 = helper_54 * helper_66;
    double helper_68 = helper_24 * helper_60;
    double helper_69 = 3.0 * _x[1];
    double helper_70 = 0.99999999999999978 * _x[10];
    double helper_71 = 1.0 * _x[4];
    double helper_72 = 1.0000000000000007 * _x[7];
    double helper_73 = -helper_69 + helper_70 + helper_71 + helper_72;
    double helper_74 = helper_35 + _x[0];
    double helper_75 = 0.66666666666666663 * helper_11;
    double helper_76 = helper_74 * helper_75;
    double helper_77 = helper_33 + _x[2];
    double helper_78 = helper_5 * helper_77;
    double helper_79 = 0.66666666666666663 * helper_78;
    double helper_80 = 0.66666666666666663 * helper_20;
    double helper_81 = helper_18 * helper_75;
    double helper_82 = helper_14 * helper_77;
    double helper_83 = 0.94280904158206336 * helper_82;
    double helper_84 = 0.94280904158206336 * helper_19;
    double helper_85 = -helper_83 - helper_84;
    double helper_86 = helper_76 + helper_79 + helper_80 + helper_81 + helper_85;
    double helper_87 = helper_83 + helper_84;
    double helper_88 = -helper_76 - helper_79 - helper_80 - helper_81 + helper_87;
    double helper_89 = helper_4 * helper_42;
    double helper_90 = helper_69 - helper_70 - helper_71 - helper_72;
    double helper_91 = 2.3570226039551585 * helper_82;
    double helper_92 = 2.3570226039551585 * helper_19;
    double helper_93 = 1.6666666666666665 * helper_11;
    double helper_94 = -helper_18 * helper_93 - 1.6666666666666665 * helper_20 - helper_74 * helper_93 -
                       1.6666666666666665 * helper_78 + helper_91 + helper_92;
    double helper_95 = 0.99999999999999978 * _x[11];
    double helper_96 = 3.0 * _x[2];
    double helper_97 = 1.0 * _x[5];
    double helper_98 = 1.0000000000000007 * _x[8];
    double helper_99 = helper_95 - helper_96 + helper_97 + helper_98;
    double helper_100 = helper_29 + _x[0];
    double helper_101 = 0.66666666666666663 * helper_16;
    double helper_102 = helper_100 * helper_101;
    double helper_103 = helper_37 + _x[1];
    double helper_104 = 0.66666666666666663 * helper_5;
    double helper_105 = helper_103 * helper_104;
    double helper_106 = helper_104 * helper_12;
    double helper_107 = helper_101 * helper_14;
    double helper_108 = 0.94280904158206336 * helper_15;
    double helper_109 = helper_103 * helper_18;
    double helper_110 = 0.94280904158206336 * helper_109;
    double helper_111 = -helper_108 - helper_110;
    double helper_112 = helper_102 + helper_105 + helper_106 + helper_107 + helper_111;
    double helper_113 = helper_108 + helper_110;
    double helper_114 = -helper_102 - helper_105 - helper_106 - helper_107 + helper_113;
    double helper_115 = -helper_95 + helper_96 - helper_97 - helper_98;
    double helper_116 = 2.3570226039551585 * helper_15;
    double helper_117 = 2.3570226039551585 * helper_109;
    double helper_118 = 1.6666666666666665 * helper_5;
    double helper_119 = 1.6666666666666665 * helper_16;
    double helper_120 =
            -helper_100 * helper_119 - helper_103 * helper_118 + helper_116 + helper_117 - helper_118 * helper_12 -
            helper_119 * helper_14;
    double helper_121 = 1.0 * _x[0];
    double helper_122 = 3.0 * _x[3];
    double helper_123 = 0.99999999999999989 * _x[6];
    double helper_124 = 1.0 * _x[9];
    double helper_125 = helper_121 - helper_122 + helper_123 + helper_124;
    double helper_126 = -1.0 * helper_25;
    double helper_127 = -helper_121 + helper_122 - helper_123 - helper_124;
    double helper_128 = -helper_64 - helper_65;
    double helper_129 = 1.0 * _x[1];
    double helper_130 = 1.0 * _x[10];
    double helper_131 = 3.0 * _x[4];
    double helper_132 = 0.99999999999999989 * _x[7];
    double helper_133 = helper_129 + helper_130 - helper_131 + helper_132;
    double helper_134 = -helper_129 - helper_130 + helper_131 - helper_132;
    double helper_135 = 0.94280904158206336 * _x[11];
    double helper_136 = -helper_135;
    double helper_137 = 0.94280904158206336 * _x[8];
    double helper_138 = helper_41 * helper_42;
    double helper_139 = -helper_91 - helper_92;
    double helper_140 = 1.0 * _x[11];
    double helper_141 = 1.0 * _x[2];
    double helper_142 = 3.0 * _x[5];
    double helper_143 = 0.99999999999999989 * _x[8];
    double helper_144 = helper_140 + helper_141 - helper_142 + helper_143;
    double helper_145 = -helper_140 - helper_141 + helper_142 - helper_143;
    double helper_146 = 0.94280904158206336 * _x[10];
    double helper_147 = 0.94280904158206336 * _x[7];
    double helper_148 = -helper_116 - helper_117;
    double helper_149 = 1.0000000000000007 * _x[0];
    double helper_150 = 0.99999999999999989 * _x[3];
    double helper_151 = 3.0000000000000004 * _x[6];
    double helper_152 = 0.99999999999999989 * _x[9];
    double helper_153 = helper_149 + helper_150 - helper_151 + helper_152;
    double helper_154 = helper_101 * helper_77;
    double helper_155 = helper_12 * helper_75;
    double helper_156 = helper_154 + helper_155;
    double helper_157 = -1.0000000000000007 * helper_25;
    double helper_158 = -helper_154 - helper_155;
    double helper_159 = -helper_149 - helper_150 + helper_151 - helper_152;
    double helper_160 = -helper_119 * helper_77 - helper_12 * helper_93;
    double helper_161 = 1.0000000000000007 * _x[1];
    double helper_162 = 0.99999999999999989 * _x[10];
    double helper_163 = 0.99999999999999989 * _x[4];
    double helper_164 = 3.0000000000000004 * _x[7];
    double helper_165 = helper_161 + helper_162 + helper_163 - helper_164;
    double helper_166 = helper_104 * helper_7;
    double helper_167 = helper_100 * helper_75;
    double helper_168 = helper_166 + helper_167;
    double helper_169 = -helper_166 - helper_167;
    double helper_170 = -helper_161 - helper_162 - helper_163 + helper_164;
    double helper_171 = 0.94280904158206336 * _x[5];
    double helper_172 = -helper_100 * helper_93 - helper_118 * helper_7;
    double helper_173 = 0.99999999999999989 * _x[11];
    double helper_174 = 1.0000000000000007 * _x[2];
    double helper_175 = 0.99999999999999989 * _x[5];
    double helper_176 = 3.0000000000000004 * _x[8];
    double helper_177 = helper_173 + helper_174 + helper_175 - helper_176;
    double helper_178 = helper_104 * helper_43;
    double helper_179 = helper_101 * helper_18;
    double helper_180 = helper_178 + helper_179;
    double helper_181 = -helper_178 - helper_179;
    double helper_182 = -helper_173 - helper_174 - helper_175 + helper_176;
    double helper_183 = -helper_146;
    double helper_184 = 0.94280904158206336 * _x[4];
    double helper_185 = -helper_118 * helper_43 - helper_119 * helper_18;
    double helper_186 = helper_138 * helper_88 + helper_25 * helper_90;
    double helper_187 = helper_186 * helper_57;
    double helper_188 = helper_186 * helper_61;
    double helper_189 = helper_42 * helper_90;
    double helper_190 = helper_66 * helper_88;
    double helper_191 = helper_42 * helper_88;
    double helper_192 = 0.94280904158206336 * _x[6];
    double helper_193 = 0.94280904158206336 * _x[9];
    double helper_194 = -helper_193;
    double helper_195 = 0.94280904158206336 * _x[3];
    double helper_196 = helper_114 * helper_138 + helper_115 * helper_25;
    double helper_197 = helper_196 * helper_57;
    double helper_198 = helper_196 * helper_61;
    double helper_199 = helper_115 * helper_42;
    double helper_200 = helper_114 * helper_66;
    double helper_201 = helper_114 * helper_42;
    double helper_202 = helper_127 * helper_25 + helper_138 * helper_59;
    double helper_203 = helper_202 * helper_57;
    double helper_204 = helper_202 * helper_61;
    double helper_205 = helper_42 * helper_59;
    double helper_206 = helper_59 * helper_66;
    double helper_207 = helper_127 * helper_42;
    double helper_208 = -0.99999999999999989 * helper_25;
    double helper_209 = 0.94280904158206336 * _x[2];
    double helper_210 = 0.94280904158206336 * _x[1];
    double helper_211 = helper_134 * helper_25 + helper_138 * helper_85;
    double helper_212 = helper_211 * helper_57;
    double helper_213 = helper_211 * helper_61;
    double helper_214 = helper_134 * helper_42;
    double helper_215 = helper_66 * helper_85;
    double helper_216 = helper_42 * helper_85;
    double helper_217 = 0.94280904158206336 * _x[0];
    double helper_218 = helper_111 * helper_138 + helper_145 * helper_25;
    double helper_219 = helper_218 * helper_57;
    double helper_220 = helper_218 * helper_61;
    double helper_221 = helper_111 * helper_42;
    double helper_222 = helper_111 * helper_66;
    double helper_223 = helper_145 * helper_42;
    double helper_224 = helper_138 * helper_158 + helper_159 * helper_25;
    double helper_225 = helper_224 * helper_57;
    double helper_226 = helper_224 * helper_61;
    double helper_227 = 3.0000000000000004 * helper_25;
    double helper_228 = helper_159 * helper_42;
    double helper_229 = helper_158 * helper_66;
    double helper_230 = helper_158 * helper_42;
    double helper_231 = helper_138 * helper_169 + helper_170 * helper_25;
    double helper_232 = helper_231 * helper_57;
    double helper_233 = helper_231 * helper_61;
    double helper_234 = helper_169 * helper_42;
    double helper_235 = helper_169 * helper_66;
    double helper_236 = helper_181 * helper_42;
    double helper_237 = helper_138 * helper_181 + helper_182 * helper_25;

    _H(0, 0) = helper_58 * (-helper_0 + helper_1 + helper_2 + helper_3) +
               helper_62 * (helper_45 + helper_46 + helper_48 + helper_49 + helper_59) + helper_68 *
                                                                                         (2 * helper_4 * helper_55 +
                                                                                          helper_63 + helper_67 *
                                                                                                      (-1.6666666666666665 *
                                                                                                       helper_21 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_22 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_44 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_47 +
                                                                                                       helper_64 +
                                                                                                       helper_65));
    _H(0, 1) = helper_58 * helper_73 + helper_62 * helper_86 +
               helper_68 * (helper_55 * helper_90 + helper_67 * helper_94 + helper_88 * helper_89);
    _H(0, 2) = helper_112 * helper_62 + helper_58 * helper_99 +
               helper_68 * (helper_114 * helper_89 + helper_115 * helper_55 + helper_120 * helper_67);
    _H(0, 3) = helper_125 * helper_58 + helper_53 * helper_62 +
               helper_68 * (helper_126 + helper_127 * helper_55 + helper_128 * helper_67 + helper_59 * helper_89);
    _H(0, 4) = helper_133 * helper_58 + helper_62 * helper_87 + helper_68 * (helper_134 * helper_55 +
                                                                             helper_138 * (helper_136 + helper_137) +
                                                                             helper_139 * helper_67 +
                                                                             helper_85 * helper_89);
    _H(0, 5) = helper_113 * helper_62 + helper_144 * helper_58 + helper_68 * (helper_111 * helper_89 +
                                                                              helper_138 * (helper_146 - helper_147) +
                                                                              helper_145 * helper_55 +
                                                                              helper_148 * helper_67);
    _H(0, 6) = helper_153 * helper_58 + helper_156 * helper_62 +
               helper_68 * (helper_157 + helper_158 * helper_89 + helper_159 * helper_55 + helper_160 * helper_67);
    _H(0, 7) = helper_165 * helper_58 + helper_168 * helper_62 + helper_68 * (helper_138 * (helper_135 - helper_171) +
                                                                              helper_169 * helper_89 +
                                                                              helper_170 * helper_55 +
                                                                              helper_172 * helper_67);
    _H(0, 8) = helper_177 * helper_58 + helper_180 * helper_62 + helper_68 * (helper_138 * (helper_183 + helper_184) +
                                                                              helper_181 * helper_89 +
                                                                              helper_182 * helper_55 +
                                                                              helper_185 * helper_67);
    _H(1, 1) = helper_187 * helper_73 + helper_188 * helper_86 +
               helper_68 * (2 * helper_189 * helper_88 + helper_190 * helper_94 + helper_63);
    _H(1, 2) = helper_112 * helper_188 + helper_187 * helper_99 +
               helper_68 * (helper_114 * helper_189 + helper_115 * helper_191 + helper_120 * helper_190);
    _H(1, 3) = helper_125 * helper_187 + helper_188 * helper_53 + helper_68 *
                                                                  (helper_127 * helper_191 + helper_128 * helper_190 +
                                                                   helper_138 * (helper_135 - helper_137) +
                                                                   helper_189 * helper_59);
    _H(1, 4) = helper_133 * helper_187 + helper_188 * helper_87 +
               helper_68 * (helper_126 + helper_134 * helper_191 + helper_139 * helper_190 + helper_189 * helper_85);
    _H(1, 5) = helper_113 * helper_188 + helper_144 * helper_187 + helper_68 * (helper_111 * helper_189 +
                                                                                helper_138 * (helper_192 + helper_194) +
                                                                                helper_145 * helper_191 +
                                                                                helper_148 * helper_190);
    _H(1, 6) = helper_153 * helper_187 + helper_156 * helper_188 + helper_68 * (helper_138 * (helper_136 + helper_171) +
                                                                                helper_158 * helper_189 +
                                                                                helper_159 * helper_191 +
                                                                                helper_160 * helper_190);
    _H(1, 7) = helper_165 * helper_187 + helper_168 * helper_188 +
               helper_68 * (helper_157 + helper_169 * helper_189 + helper_170 * helper_191 + helper_172 * helper_190);
    _H(1, 8) = helper_177 * helper_187 + helper_180 * helper_188 + helper_68 * (helper_138 * (helper_193 - helper_195) +
                                                                                helper_181 * helper_189 +
                                                                                helper_182 * helper_191 +
                                                                                helper_185 * helper_190);
    _H(2, 2) = helper_112 * helper_198 + helper_197 * helper_99 +
               helper_68 * (2 * helper_114 * helper_199 + helper_120 * helper_200 + helper_63);
    _H(2, 3) = helper_125 * helper_197 + helper_198 * helper_53 + helper_68 *
                                                                  (helper_127 * helper_201 + helper_128 * helper_200 +
                                                                   helper_138 * (helper_147 + helper_183) +
                                                                   helper_199 * helper_59);
    _H(2, 4) = helper_133 * helper_197 + helper_198 * helper_87 + helper_68 * (helper_134 * helper_201 +
                                                                               helper_138 * (-helper_192 + helper_193) +
                                                                               helper_139 * helper_200 +
                                                                               helper_199 * helper_85);
    _H(2, 5) = helper_113 * helper_198 + helper_144 * helper_197 +
               helper_68 * (helper_111 * helper_199 + helper_126 + helper_145 * helper_201 + helper_148 * helper_200);
    _H(2, 6) = helper_153 * helper_197 + helper_156 * helper_198 + helper_68 * (helper_138 * (helper_146 - helper_184) +
                                                                                helper_158 * helper_199 +
                                                                                helper_159 * helper_201 +
                                                                                helper_160 * helper_200);
    _H(2, 7) = helper_165 * helper_197 + helper_168 * helper_198 + helper_68 * (helper_138 * (helper_194 + helper_195) +
                                                                                helper_169 * helper_199 +
                                                                                helper_170 * helper_201 +
                                                                                helper_172 * helper_200);
    _H(2, 8) = helper_177 * helper_197 + helper_180 * helper_198 +
               helper_68 * (helper_157 + helper_181 * helper_199 + helper_182 * helper_201 + helper_185 * helper_200);
    _H(3, 3) = helper_125 * helper_203 + helper_204 * helper_53 +
               helper_68 * (2 * helper_127 * helper_205 + helper_128 * helper_206 + helper_63);
    _H(3, 4) = helper_133 * helper_203 + helper_204 * helper_87 +
               helper_68 * (helper_134 * helper_205 + helper_139 * helper_206 + helper_207 * helper_85);
    _H(3, 5) = helper_113 * helper_204 + helper_144 * helper_203 +
               helper_68 * (helper_111 * helper_207 + helper_145 * helper_205 + helper_148 * helper_206);
    _H(3, 6) = helper_153 * helper_203 + helper_156 * helper_204 +
               helper_68 * (helper_158 * helper_207 + helper_159 * helper_205 + helper_160 * helper_206 + helper_208);
    _H(3, 7) = helper_165 * helper_203 + helper_168 * helper_204 + helper_68 * (helper_138 * (helper_136 + helper_209) +
                                                                                helper_169 * helper_207 +
                                                                                helper_170 * helper_205 +
                                                                                helper_172 * helper_206);
    _H(3, 8) = helper_177 * helper_203 + helper_180 * helper_204 + helper_68 * (helper_138 * (helper_146 - helper_210) +
                                                                                helper_181 * helper_207 +
                                                                                helper_182 * helper_205 +
                                                                                helper_185 * helper_206);
    _H(4, 4) = helper_133 * helper_212 + helper_213 * helper_87 +
               helper_68 * (helper_139 * helper_215 + 2 * helper_214 * helper_85 + helper_63);
    _H(4, 5) = helper_113 * helper_213 + helper_144 * helper_212 +
               helper_68 * (helper_111 * helper_214 + helper_145 * helper_216 + helper_148 * helper_215);
    _H(4, 6) = helper_153 * helper_212 + helper_156 * helper_213 + helper_68 * (helper_138 * (helper_135 - helper_209) +
                                                                                helper_158 * helper_214 +
                                                                                helper_159 * helper_216 +
                                                                                helper_160 * helper_215);
    _H(4, 7) = helper_165 * helper_212 + helper_168 * helper_213 +
               helper_68 * (helper_169 * helper_214 + helper_170 * helper_216 + helper_172 * helper_215 + helper_208);
    _H(4, 8) = helper_177 * helper_212 + helper_180 * helper_213 + helper_68 * (helper_138 * (helper_194 + helper_217) +
                                                                                helper_181 * helper_214 +
                                                                                helper_182 * helper_216 +
                                                                                helper_185 * helper_215);
    _H(5, 5) = helper_113 * helper_220 + helper_144 * helper_219 +
               helper_68 * (2 * helper_145 * helper_221 + helper_148 * helper_222 + helper_63);
    _H(5, 6) = helper_153 * helper_219 + helper_156 * helper_220 + helper_68 * (helper_138 * (helper_183 + helper_210) +
                                                                                helper_158 * helper_223 +
                                                                                helper_159 * helper_221 +
                                                                                helper_160 * helper_222);
    _H(5, 7) = helper_165 * helper_219 + helper_168 * helper_220 + helper_68 * (helper_138 * (helper_193 - helper_217) +
                                                                                helper_169 * helper_223 +
                                                                                helper_170 * helper_221 +
                                                                                helper_172 * helper_222);
    _H(5, 8) = helper_177 * helper_219 + helper_180 * helper_220 +
               helper_68 * (helper_181 * helper_223 + helper_182 * helper_221 + helper_185 * helper_222 + helper_208);
    _H(6, 6) = helper_153 * helper_225 + helper_156 * helper_226 +
               helper_68 * (2 * helper_158 * helper_228 + helper_160 * helper_229 + helper_227);
    _H(6, 7) = helper_165 * helper_225 + helper_168 * helper_226 +
               helper_68 * (helper_169 * helper_228 + helper_170 * helper_230 + helper_172 * helper_229);
    _H(6, 8) = helper_177 * helper_225 + helper_180 * helper_226 +
               helper_68 * (helper_181 * helper_228 + helper_182 * helper_230 + helper_185 * helper_229);
    _H(7, 7) = helper_165 * helper_232 + helper_168 * helper_233 +
               helper_68 * (2 * helper_170 * helper_234 + helper_172 * helper_235 + helper_227);
    _H(7, 8) = helper_177 * helper_232 + helper_180 * helper_233 +
               helper_68 * (helper_170 * helper_236 + helper_182 * helper_234 + helper_185 * helper_235);
    _H(8, 8) = helper_177 * helper_237 * helper_57 + helper_180 * helper_237 * helper_61 +
               helper_68 * (helper_181 * helper_185 * helper_66 + 2 * helper_182 * helper_236 + helper_227);


    for (int i = 1; i < 9; ++i)
      for (int j = 0; j < i; ++j)
        _H(i, j) = _H(j, i);

    // Compute Eigen decomposition H = V D V^T
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_H);
    const Eigen::MatrixXd &V = solver.eigenvectors();
    Eigen::VectorXd evals = solver.eigenvalues();

    // check eigenvalues against epsilon, for those less that eps (zero)
    const int n = n_unknowns();
    for (int i = 0; i < n; ++i)
    {
      if (evals[i] < 1e-8)
      {
        evals[i] = 1e-8;
      }
    }

    // compute correction matrix M = V * diag(m) * V^T
    _H.noalias() = V * evals.asDiagonal() * V.transpose();
  }
};


//Four vertices
class IMRMEnergy4
{
public:
  static int n_unknowns() { return 12; }

  static double eval_f(const double *_x)
  {
    //trace
    double helper_0 = -_x[0];
    double helper_1 = helper_0 + _x[3];
    double helper_2 = -_x[1];
    double helper_3 = helper_2 + _x[7];
    double helper_4 = -_x[2];
    double helper_5 = helper_4 + _x[11];
    double helper_6 = helper_0 + _x[6];
    double helper_7 = helper_2 + _x[10];
    double helper_8 = helper_4 + _x[5];
    double helper_9 = helper_0 + _x[9];
    double helper_10 = helper_2 + _x[4];
    double helper_11 = helper_4 + _x[8];

    return std::log(std::pow(-helper_1 * helper_11 * helper_7 + helper_1 * helper_3 * helper_5 +
                             helper_10 * helper_11 * helper_9 - helper_10 * helper_5 * helper_6 -
                             helper_3 * helper_8 * helper_9 + helper_6 * helper_7 * helper_8,
                             -0.66666666666666663) *
                    (std::pow(_x[0] - _x[3], 2) + std::pow(_x[1] - _x[4], 2) + std::pow(_x[2] - _x[5], 2) +
                     1.3333333333333333 * std::pow(0.5 * _x[0] + 0.5 * _x[3] - _x[9], 2) +
                     1.3333333333333333 * std::pow(0.5 * _x[1] - _x[10] + 0.5 * _x[4], 2) +
                     1.3333333333333333 * std::pow(-_x[11] + 0.5 * _x[2] + 0.5 * _x[5], 2) +
                     1.5000000000000002 * std::pow(
                             0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] - _x[6] +
                             0.33333333333333326 * _x[9], 2) + 1.5000000000000002 * std::pow(
                            0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] +
                            0.33333333333333326 * _x[4] - _x[7], 2) + 1.5000000000000002 * std::pow(
                            0.33333333333333326 * _x[11] + 0.33333333333333348 * _x[2] +
                            0.33333333333333326 * _x[5] - _x[8], 2))) - 0.2310490601866485;
  }

  static void eval_gradient(const double *_x, double *_g)
  {
    //trace
    double helper_0 = -1.0 * _x[3];
    double helper_1 = -1.4142135623730951 * _x[0] + 1.4142135623730951 * _x[3];
    double helper_2 = -_x[2];
    double helper_3 = helper_2 + _x[11];
    double helper_4 = -_x[1];
    double helper_5 = helper_4 + _x[7];
    double helper_6 = helper_3 * helper_5;
    double helper_7 = -1.4142135623730951 * _x[2] + 1.4142135623730951 * _x[5];
    double helper_8 = helper_4 + _x[10];
    double helper_9 = -_x[0];
    double helper_10 = helper_9 + _x[6];
    double helper_11 = helper_10 * helper_8;
    double helper_12 = -1.4142135623730951 * _x[1] + 1.4142135623730951 * _x[4];
    double helper_13 = helper_2 + _x[8];
    double helper_14 = helper_9 + _x[9];
    double helper_15 = helper_13 * helper_14;
    double helper_16 = helper_1 * helper_13;
    double helper_17 = helper_12 * helper_3;
    double helper_18 = helper_5 * helper_7;
    double helper_19 = helper_1 * helper_6 - helper_10 * helper_17 + helper_11 * helper_7 + helper_12 * helper_15 -
                       helper_14 * helper_18 - helper_16 * helper_8;
    double helper_20 = std::pow(helper_19, 0.66666666666666663);
    double helper_21 = 1.0 / helper_20;
    double helper_22 = 0.94280904158206336 * helper_6;
    double helper_23 = -_x[10];
    double helper_24 = helper_23 + _x[1];
    double helper_25 = 0.94280904158206336 * helper_13 * helper_24;
    double helper_26 = -_x[8];
    double helper_27 = helper_26 + _x[2];
    double helper_28 = 0.66666666666666663 * helper_12;
    double helper_29 = 0.66666666666666663 * helper_7;
    double helper_30 = -_x[9];
    double helper_31 = -_x[11];
    double helper_32 = -_x[6];
    double helper_33 = -_x[7];
    double helper_34 = std::pow(helper_2 + _x[5], 2) + std::pow(helper_4 + _x[4], 2) + std::pow(helper_9 + _x[3], 2) +
                       1.3333333333333333 * std::pow(helper_23 + 0.5 * _x[1] + 0.5 * _x[4], 2) +
                       1.3333333333333333 * std::pow(helper_30 + 0.5 * _x[0] + 0.5 * _x[3], 2) +
                       1.3333333333333333 * std::pow(helper_31 + 0.5 * _x[2] + 0.5 * _x[5], 2) + 1.5000000000000002 *
                                                                                                 std::pow(helper_26 +
                                                                                                          0.33333333333333326 *
                                                                                                          _x[11] +
                                                                                                          0.33333333333333348 *
                                                                                                          _x[2] +
                                                                                                          0.33333333333333326 *
                                                                                                          _x[5], 2) +
                       1.5000000000000002 * std::pow(
                               helper_32 + 0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] +
                               0.33333333333333326 * _x[9], 2) + 1.5000000000000002 * std::pow(
            helper_33 + 0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] + 0.33333333333333326 * _x[4], 2);
    double helper_35 = std::pow(helper_19, -1.6666666666666665) * helper_34;
    double helper_36 = helper_20 / helper_34;
    double helper_37 = -1.0 * _x[4];
    double helper_38 = helper_31 + _x[2];
    double helper_39 = 0.94280904158206336 * helper_10 * helper_38;
    double helper_40 = 0.94280904158206336 * helper_15;
    double helper_41 = 0.66666666666666663 * helper_1;
    double helper_42 = helper_32 + _x[0];
    double helper_43 = -1.0 * _x[5];
    double helper_44 = 0.94280904158206336 * helper_11;
    double helper_45 = helper_33 + _x[1];
    double helper_46 = 0.94280904158206336 * helper_14 * helper_45;
    double helper_47 = helper_30 + _x[0];
    double helper_48 = -0.99999999999999989 * _x[6];
    double helper_49 = -0.99999999999999989 * _x[7];
    double helper_50 = -0.99999999999999989 * _x[8];

    _g[0] = helper_36 *
            (helper_21 * (helper_0 + 3.0 * _x[0] - 1.0000000000000007 * _x[6] - 0.99999999999999978 * _x[9]) +
             helper_35 *
             (-0.66666666666666663 * helper_17 - 0.66666666666666663 * helper_18 + helper_22 - helper_24 * helper_29 +
              helper_25 - helper_27 * helper_28));
    _g[1] = helper_36 *
            (helper_21 * (helper_37 + 3.0 * _x[1] - 0.99999999999999978 * _x[10] - 1.0000000000000007 * _x[7]) +
             helper_35 *
             (-helper_14 * helper_29 - 0.66666666666666663 * helper_16 - helper_29 * helper_42 - helper_38 * helper_41 +
              helper_39 + helper_40));
    _g[2] = helper_36 *
            (helper_21 * (helper_43 - 0.99999999999999978 * _x[11] + 3.0 * _x[2] - 1.0000000000000007 * _x[8]) +
             helper_35 *
             (-helper_10 * helper_28 - helper_28 * helper_47 - helper_41 * helper_45 - helper_41 * helper_8 +
              helper_44 + helper_46));
    _g[3] = helper_36 *
            (helper_21 * (helper_48 - 1.0 * _x[0] + 3.0 * _x[3] - 1.0 * _x[9]) + helper_35 * (-helper_22 - helper_25));
    _g[4] = helper_36 *
            (helper_21 * (helper_49 - 1.0 * _x[1] - 1.0 * _x[10] + 3.0 * _x[4]) + helper_35 * (-helper_39 - helper_40));
    _g[5] = helper_36 *
            (helper_21 * (helper_50 - 1.0 * _x[11] - 1.0 * _x[2] + 3.0 * _x[5]) + helper_35 * (-helper_44 - helper_46));
    _g[6] = helper_36 * (helper_21 *
                         (-1.0000000000000007 * _x[0] - 0.99999999999999989 * _x[3] + 3.0000000000000004 * _x[6] -
                          0.99999999999999989 * _x[9]) + helper_35 * (-helper_28 * helper_38 - helper_29 * helper_8));
    _g[7] = helper_36 * (helper_21 *
                         (-1.0000000000000007 * _x[1] - 0.99999999999999989 * _x[10] - 0.99999999999999989 * _x[4] +
                          3.0000000000000004 * _x[7]) + helper_35 * (-helper_29 * helper_47 - helper_3 * helper_41));
    _g[8] = helper_36 * (helper_21 *
                         (-0.99999999999999989 * _x[11] - 1.0000000000000007 * _x[2] - 0.99999999999999989 * _x[5] +
                          3.0000000000000004 * _x[8]) + helper_35 * (-helper_14 * helper_28 - helper_24 * helper_41));
    _g[9] = helper_36 * (helper_21 * (helper_0 + helper_48 - 0.99999999999999978 * _x[0] + 3.0 * _x[9]) +
                         helper_35 * (-helper_13 * helper_28 - helper_29 * helper_45));
    _g[10] = helper_36 * (helper_21 * (helper_37 + helper_49 - 0.99999999999999978 * _x[1] + 3.0 * _x[10]) +
                          helper_35 * (-helper_10 * helper_29 - helper_27 * helper_41));
    _g[11] = helper_36 * (helper_21 * (helper_43 + helper_50 + 3.0 * _x[11] - 0.99999999999999978 * _x[2]) +
                          helper_35 * (-helper_28 * helper_42 - helper_41 * helper_5));
  }

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H)
  {
    _H.setZero();
    //trace
    double helper_0 = 3.0 * _x[0];
    double helper_1 = 1.0 * _x[3];
    double helper_2 = 1.0000000000000007 * _x[6];
    double helper_3 = 0.99999999999999978 * _x[9];
    double helper_4 = -helper_1;
    double helper_5 = helper_0 - helper_2 - helper_3 + helper_4;
    double helper_6 = -1.4142135623730951 * _x[0] + 1.4142135623730951 * _x[3];
    double helper_7 = -_x[2];
    double helper_8 = helper_7 + _x[11];
    double helper_9 = -_x[1];
    double helper_10 = helper_9 + _x[7];
    double helper_11 = helper_10 * helper_8;
    double helper_12 = -1.4142135623730951 * _x[2] + 1.4142135623730951 * _x[5];
    double helper_13 = helper_9 + _x[10];
    double helper_14 = -_x[0];
    double helper_15 = helper_14 + _x[6];
    double helper_16 = helper_13 * helper_15;
    double helper_17 = -1.4142135623730951 * _x[1] + 1.4142135623730951 * _x[4];
    double helper_18 = helper_7 + _x[8];
    double helper_19 = helper_14 + _x[9];
    double helper_20 = helper_18 * helper_19;
    double helper_21 = helper_18 * helper_6;
    double helper_22 = helper_17 * helper_8;
    double helper_23 = helper_10 * helper_12;
    double helper_24 = helper_11 * helper_6 + helper_12 * helper_16 - helper_13 * helper_21 - helper_15 * helper_22 +
                       helper_17 * helper_20 - helper_19 * helper_23;
    double helper_25 = std::pow(helper_24, 0.66666666666666663);
    double helper_26 = 1.0 / helper_25;
    double helper_27 = std::pow(helper_14 + _x[3], 2);
    double helper_28 = std::pow(helper_9 + _x[4], 2);
    double helper_29 = std::pow(helper_7 + _x[5], 2);
    double helper_30 = -_x[9];
    double helper_31 = std::pow(helper_30 + 0.5 * _x[0] + 0.5 * _x[3], 2);
    double helper_32 = -_x[10];
    double helper_33 = std::pow(helper_32 + 0.5 * _x[1] + 0.5 * _x[4], 2);
    double helper_34 = -_x[11];
    double helper_35 = std::pow(helper_34 + 0.5 * _x[2] + 0.5 * _x[5], 2);
    double helper_36 = -_x[6];
    double helper_37 = std::pow(
            helper_36 + 0.33333333333333348 * _x[0] + 0.33333333333333326 * _x[3] + 0.33333333333333326 * _x[9], 2);
    double helper_38 = -_x[7];
    double helper_39 = std::pow(
            helper_38 + 0.33333333333333348 * _x[1] + 0.33333333333333326 * _x[10] + 0.33333333333333326 * _x[4], 2);
    double helper_40 = -_x[8];
    double helper_41 = std::pow(
            helper_40 + 0.33333333333333326 * _x[11] + 0.33333333333333348 * _x[2] + 0.33333333333333326 * _x[5], 2);
    double helper_42 =
            helper_27 + helper_28 + helper_29 + 1.3333333333333333 * helper_31 + 1.3333333333333333 * helper_33 +
            1.3333333333333333 * helper_35 + 1.5000000000000002 * helper_37 + 1.5000000000000002 * helper_39 +
            1.5000000000000002 * helper_41;
    double helper_43 = std::pow(helper_24, -1.6666666666666665);
    double helper_44 = helper_32 + _x[1];
    double helper_45 = helper_12 * helper_44;
    double helper_46 = 0.66666666666666663 * helper_45;
    double helper_47 = 0.66666666666666663 * helper_22;
    double helper_48 = helper_40 + _x[2];
    double helper_49 = helper_17 * helper_48;
    double helper_50 = 0.66666666666666663 * helper_49;
    double helper_51 = 0.66666666666666663 * helper_23;
    double helper_52 = 0.94280904158206336 * helper_11;
    double helper_53 = helper_18 * helper_44;
    double helper_54 = 0.94280904158206336 * helper_53;
    double helper_55 = helper_52 + helper_54;
    double helper_56 = -helper_46 - helper_47 - helper_50 - helper_51 + helper_55;
    double helper_57 = helper_43 * helper_56;
    double helper_58 = helper_26 * helper_5 + helper_42 * helper_57;
    double helper_59 = 0.44444444444444431 * helper_25 / std::pow(
            0.66666666666666652 * helper_27 + 0.66666666666666652 * helper_28 + 0.66666666666666652 * helper_29 +
            0.88888888888888873 * helper_31 + 0.88888888888888873 * helper_33 + 0.88888888888888873 * helper_35 +
            helper_37 + helper_39 + helper_41, 2);
    double helper_60 = helper_58 * helper_59;
    double helper_61 = -helper_52 - helper_54;
    double helper_62 = 1.0 / helper_42;
    double helper_63 = std::pow(helper_24, -0.33333333333333337) * helper_62;
    double helper_64 = helper_58 * helper_63;
    double helper_65 = 3.0 * helper_26;
    double helper_66 = 2.3570226039551585 * helper_11;
    double helper_67 = 2.3570226039551585 * helper_53;
    double helper_68 = std::pow(helper_24, -2.6666666666666665) * helper_42;
    double helper_69 = helper_56 * helper_68;
    double helper_70 = helper_25 * helper_62;
    double helper_71 = 3.0 * _x[1];
    double helper_72 = 0.99999999999999978 * _x[10];
    double helper_73 = 1.0 * _x[4];
    double helper_74 = 1.0000000000000007 * _x[7];
    double helper_75 = -helper_71 + helper_72 + helper_73 + helper_74;
    double helper_76 = helper_36 + _x[0];
    double helper_77 = 0.66666666666666663 * helper_12;
    double helper_78 = helper_76 * helper_77;
    double helper_79 = helper_34 + _x[2];
    double helper_80 = helper_6 * helper_79;
    double helper_81 = 0.66666666666666663 * helper_80;
    double helper_82 = 0.66666666666666663 * helper_21;
    double helper_83 = helper_19 * helper_77;
    double helper_84 = helper_15 * helper_79;
    double helper_85 = 0.94280904158206336 * helper_84;
    double helper_86 = 0.94280904158206336 * helper_20;
    double helper_87 = -helper_85 - helper_86;
    double helper_88 = helper_78 + helper_81 + helper_82 + helper_83 + helper_87;
    double helper_89 = helper_85 + helper_86;
    double helper_90 = -helper_78 - helper_81 - helper_82 - helper_83 + helper_89;
    double helper_91 = helper_43 * helper_5;
    double helper_92 = -helper_73;
    double helper_93 = helper_71 - helper_72 - helper_74 + helper_92;
    double helper_94 = 2.3570226039551585 * helper_84;
    double helper_95 = 2.3570226039551585 * helper_20;
    double helper_96 = 1.6666666666666665 * helper_12;
    double helper_97 = -helper_19 * helper_96 - 1.6666666666666665 * helper_21 - helper_76 * helper_96 -
                       1.6666666666666665 * helper_80 + helper_94 + helper_95;
    double helper_98 = 0.99999999999999978 * _x[11];
    double helper_99 = 3.0 * _x[2];
    double helper_100 = 1.0 * _x[5];
    double helper_101 = 1.0000000000000007 * _x[8];
    double helper_102 = helper_100 + helper_101 + helper_98 - helper_99;
    double helper_103 = helper_30 + _x[0];
    double helper_104 = 0.66666666666666663 * helper_17;
    double helper_105 = helper_103 * helper_104;
    double helper_106 = helper_38 + _x[1];
    double helper_107 = 0.66666666666666663 * helper_6;
    double helper_108 = helper_106 * helper_107;
    double helper_109 = helper_107 * helper_13;
    double helper_110 = helper_104 * helper_15;
    double helper_111 = 0.94280904158206336 * helper_16;
    double helper_112 = helper_106 * helper_19;
    double helper_113 = 0.94280904158206336 * helper_112;
    double helper_114 = -helper_111 - helper_113;
    double helper_115 = helper_105 + helper_108 + helper_109 + helper_110 + helper_114;
    double helper_116 = helper_111 + helper_113;
    double helper_117 = -helper_105 - helper_108 - helper_109 - helper_110 + helper_116;
    double helper_118 = -helper_100;
    double helper_119 = -helper_101 + helper_118 - helper_98 + helper_99;
    double helper_120 = 2.3570226039551585 * helper_16;
    double helper_121 = 2.3570226039551585 * helper_112;
    double helper_122 = 1.6666666666666665 * helper_6;
    double helper_123 = 1.6666666666666665 * helper_17;
    double helper_124 =
            -helper_103 * helper_123 - helper_106 * helper_122 + helper_120 + helper_121 - helper_122 * helper_13 -
            helper_123 * helper_15;
    double helper_125 = 1.0 * _x[0];
    double helper_126 = 3.0 * _x[3];
    double helper_127 = 0.99999999999999989 * _x[6];
    double helper_128 = 1.0 * _x[9];
    double helper_129 = helper_125 - helper_126 + helper_127 + helper_128;
    double helper_130 = -1.0 * helper_26;
    double helper_131 = -helper_127;
    double helper_132 = -helper_125 + helper_126 - helper_128 + helper_131;
    double helper_133 = -helper_66 - helper_67;
    double helper_134 = 1.0 * _x[1];
    double helper_135 = 1.0 * _x[10];
    double helper_136 = 3.0 * _x[4];
    double helper_137 = 0.99999999999999989 * _x[7];
    double helper_138 = helper_134 + helper_135 - helper_136 + helper_137;
    double helper_139 = -helper_137;
    double helper_140 = -helper_134 - helper_135 + helper_136 + helper_139;
    double helper_141 = 0.94280904158206336 * _x[11];
    double helper_142 = -helper_141;
    double helper_143 = 0.94280904158206336 * _x[8];
    double helper_144 = helper_42 * helper_43;
    double helper_145 = -helper_94 - helper_95;
    double helper_146 = 1.0 * _x[11];
    double helper_147 = 1.0 * _x[2];
    double helper_148 = 3.0 * _x[5];
    double helper_149 = 0.99999999999999989 * _x[8];
    double helper_150 = helper_146 + helper_147 - helper_148 + helper_149;
    double helper_151 = -helper_149;
    double helper_152 = -helper_146 - helper_147 + helper_148 + helper_151;
    double helper_153 = 0.94280904158206336 * _x[10];
    double helper_154 = 0.94280904158206336 * _x[7];
    double helper_155 = -helper_154;
    double helper_156 = -helper_120 - helper_121;
    double helper_157 = 1.0000000000000007 * _x[0];
    double helper_158 = 0.99999999999999989 * _x[3];
    double helper_159 = 3.0000000000000004 * _x[6];
    double helper_160 = 0.99999999999999989 * _x[9];
    double helper_161 = helper_157 + helper_158 - helper_159 + helper_160;
    double helper_162 = helper_104 * helper_79;
    double helper_163 = helper_13 * helper_77;
    double helper_164 = helper_162 + helper_163;
    double helper_165 = -1.0000000000000007 * helper_26;
    double helper_166 = -helper_162 - helper_163;
    double helper_167 = -helper_157 - helper_158 + helper_159 - helper_160;
    double helper_168 = -helper_123 * helper_79 - helper_13 * helper_96;
    double helper_169 = 1.0000000000000007 * _x[1];
    double helper_170 = 0.99999999999999989 * _x[10];
    double helper_171 = 0.99999999999999989 * _x[4];
    double helper_172 = 3.0000000000000004 * _x[7];
    double helper_173 = helper_169 + helper_170 + helper_171 - helper_172;
    double helper_174 = helper_107 * helper_8;
    double helper_175 = helper_103 * helper_77;
    double helper_176 = helper_174 + helper_175;
    double helper_177 = -helper_174 - helper_175;
    double helper_178 = -helper_169 - helper_170 - helper_171 + helper_172;
    double helper_179 = 0.94280904158206336 * _x[5];
    double helper_180 = -helper_179;
    double helper_181 = -helper_103 * helper_96 - helper_122 * helper_8;
    double helper_182 = 0.99999999999999989 * _x[11];
    double helper_183 = 1.0000000000000007 * _x[2];
    double helper_184 = 0.99999999999999989 * _x[5];
    double helper_185 = 3.0000000000000004 * _x[8];
    double helper_186 = helper_182 + helper_183 + helper_184 - helper_185;
    double helper_187 = helper_107 * helper_44;
    double helper_188 = helper_104 * helper_19;
    double helper_189 = helper_187 + helper_188;
    double helper_190 = -helper_187 - helper_188;
    double helper_191 = -helper_182 - helper_183 - helper_184 + helper_185;
    double helper_192 = -helper_153;
    double helper_193 = 0.94280904158206336 * _x[4];
    double helper_194 = -helper_122 * helper_44 - helper_123 * helper_19;
    double helper_195 = 0.99999999999999978 * _x[0];
    double helper_196 = 3.0 * _x[9];
    double helper_197 = helper_1 + helper_127 + helper_195 - helper_196;
    double helper_198 = helper_104 * helper_18;
    double helper_199 = helper_106 * helper_77;
    double helper_200 = helper_198 + helper_199;
    double helper_201 = -0.99999999999999978 * helper_26;
    double helper_202 = -helper_198 - helper_199;
    double helper_203 = helper_131 - helper_195 + helper_196 + helper_4;
    double helper_204 = -helper_106 * helper_96 - helper_123 * helper_18;
    double helper_205 = 0.99999999999999978 * _x[1];
    double helper_206 = 3.0 * _x[10];
    double helper_207 = helper_137 + helper_205 - helper_206 + helper_73;
    double helper_208 = helper_107 * helper_48;
    double helper_209 = helper_15 * helper_77;
    double helper_210 = helper_208 + helper_209;
    double helper_211 = -helper_208 - helper_209;
    double helper_212 = helper_139 - helper_205 + helper_206 + helper_92;
    double helper_213 = -helper_143;
    double helper_214 = -helper_122 * helper_48 - helper_15 * helper_96;
    double helper_215 = 3.0 * _x[11];
    double helper_216 = 0.99999999999999978 * _x[2];
    double helper_217 = helper_100 + helper_149 - helper_215 + helper_216;
    double helper_218 = helper_10 * helper_107;
    double helper_219 = helper_104 * helper_76;
    double helper_220 = helper_218 + helper_219;
    double helper_221 = -helper_218 - helper_219;
    double helper_222 = helper_118 + helper_151 + helper_215 - helper_216;
    double helper_223 = -helper_193;
    double helper_224 = -helper_10 * helper_122 - helper_123 * helper_76;
    double helper_225 = helper_144 * helper_90 + helper_26 * helper_93;
    double helper_226 = helper_225 * helper_59;
    double helper_227 = helper_225 * helper_63;
    double helper_228 = helper_43 * helper_93;
    double helper_229 = helper_68 * helper_90;
    double helper_230 = helper_43 * helper_90;
    double helper_231 = 0.94280904158206336 * _x[6];
    double helper_232 = 0.94280904158206336 * _x[9];
    double helper_233 = -helper_232;
    double helper_234 = 0.94280904158206336 * _x[3];
    double helper_235 = -helper_234;
    double helper_236 = -helper_231;
    double helper_237 = helper_117 * helper_144 + helper_119 * helper_26;
    double helper_238 = helper_237 * helper_59;
    double helper_239 = helper_237 * helper_63;
    double helper_240 = helper_119 * helper_43;
    double helper_241 = helper_117 * helper_68;
    double helper_242 = helper_117 * helper_43;
    double helper_243 = helper_132 * helper_26 + helper_144 * helper_61;
    double helper_244 = helper_243 * helper_59;
    double helper_245 = helper_243 * helper_63;
    double helper_246 = helper_43 * helper_61;
    double helper_247 = helper_61 * helper_68;
    double helper_248 = helper_132 * helper_43;
    double helper_249 = -0.99999999999999989 * helper_26;
    double helper_250 = 0.94280904158206336 * _x[2];
    double helper_251 = 0.94280904158206336 * _x[1];
    double helper_252 = -helper_251;
    double helper_253 = -helper_250;
    double helper_254 = helper_140 * helper_26 + helper_144 * helper_87;
    double helper_255 = helper_254 * helper_59;
    double helper_256 = helper_254 * helper_63;
    double helper_257 = helper_140 * helper_43;
    double helper_258 = helper_68 * helper_87;
    double helper_259 = helper_43 * helper_87;
    double helper_260 = 0.94280904158206336 * _x[0];
    double helper_261 = -helper_260;
    double helper_262 = helper_114 * helper_144 + helper_152 * helper_26;
    double helper_263 = helper_262 * helper_59;
    double helper_264 = helper_262 * helper_63;
    double helper_265 = helper_114 * helper_43;
    double helper_266 = helper_114 * helper_68;
    double helper_267 = helper_152 * helper_43;
    double helper_268 = helper_144 * helper_166 + helper_167 * helper_26;
    double helper_269 = helper_268 * helper_59;
    double helper_270 = helper_268 * helper_63;
    double helper_271 = 3.0000000000000004 * helper_26;
    double helper_272 = helper_167 * helper_43;
    double helper_273 = helper_166 * helper_68;
    double helper_274 = helper_166 * helper_43;
    double helper_275 = helper_144 * helper_177 + helper_178 * helper_26;
    double helper_276 = helper_275 * helper_59;
    double helper_277 = helper_275 * helper_63;
    double helper_278 = helper_177 * helper_43;
    double helper_279 = helper_177 * helper_68;
    double helper_280 = helper_178 * helper_43;
    double helper_281 = helper_144 * helper_190 + helper_191 * helper_26;
    double helper_282 = helper_281 * helper_59;
    double helper_283 = helper_281 * helper_63;
    double helper_284 = helper_190 * helper_43;
    double helper_285 = helper_190 * helper_68;
    double helper_286 = helper_191 * helper_43;
    double helper_287 = helper_144 * helper_202 + helper_203 * helper_26;
    double helper_288 = helper_287 * helper_59;
    double helper_289 = helper_287 * helper_63;
    double helper_290 = helper_202 * helper_43;
    double helper_291 = helper_202 * helper_68;
    double helper_292 = helper_203 * helper_43;
    double helper_293 = helper_144 * helper_211 + helper_212 * helper_26;
    double helper_294 = helper_293 * helper_59;
    double helper_295 = helper_293 * helper_63;
    double helper_296 = helper_212 * helper_43;
    double helper_297 = helper_211 * helper_68;
    double helper_298 = helper_222 * helper_43;
    double helper_299 = helper_144 * helper_221 + helper_222 * helper_26;

    _H(0, 0) = helper_60 * (-helper_0 + helper_1 + helper_2 + helper_3) +
               helper_64 * (helper_46 + helper_47 + helper_50 + helper_51 + helper_61) + helper_70 *
                                                                                         (2 * helper_5 * helper_57 +
                                                                                          helper_65 + helper_69 *
                                                                                                      (-1.6666666666666665 *
                                                                                                       helper_22 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_23 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_45 -
                                                                                                       1.6666666666666665 *
                                                                                                       helper_49 +
                                                                                                       helper_66 +
                                                                                                       helper_67));
    _H(0, 1) = helper_60 * helper_75 + helper_64 * helper_88 +
               helper_70 * (helper_57 * helper_93 + helper_69 * helper_97 + helper_90 * helper_91);
    _H(0, 2) = helper_102 * helper_60 + helper_115 * helper_64 +
               helper_70 * (helper_117 * helper_91 + helper_119 * helper_57 + helper_124 * helper_69);
    _H(0, 3) = helper_129 * helper_60 + helper_55 * helper_64 +
               helper_70 * (helper_130 + helper_132 * helper_57 + helper_133 * helper_69 + helper_61 * helper_91);
    _H(0, 4) = helper_138 * helper_60 + helper_64 * helper_89 + helper_70 * (helper_140 * helper_57 +
                                                                             helper_144 * (helper_142 + helper_143) +
                                                                             helper_145 * helper_69 +
                                                                             helper_87 * helper_91);
    _H(0, 5) = helper_116 * helper_64 + helper_150 * helper_60 + helper_70 * (helper_114 * helper_91 +
                                                                              helper_144 * (helper_153 + helper_155) +
                                                                              helper_152 * helper_57 +
                                                                              helper_156 * helper_69);
    _H(0, 6) = helper_161 * helper_60 + helper_164 * helper_64 +
               helper_70 * (helper_165 + helper_166 * helper_91 + helper_167 * helper_57 + helper_168 * helper_69);
    _H(0, 7) = helper_173 * helper_60 + helper_176 * helper_64 + helper_70 * (helper_144 * (helper_141 + helper_180) +
                                                                              helper_177 * helper_91 +
                                                                              helper_178 * helper_57 +
                                                                              helper_181 * helper_69);
    _H(0, 8) = helper_186 * helper_60 + helper_189 * helper_64 + helper_70 * (helper_144 * (helper_192 + helper_193) +
                                                                              helper_190 * helper_91 +
                                                                              helper_191 * helper_57 +
                                                                              helper_194 * helper_69);
    _H(0, 9) = helper_197 * helper_60 + helper_200 * helper_64 +
               helper_70 * (helper_201 + helper_202 * helper_91 + helper_203 * helper_57 + helper_204 * helper_69);
    _H(0, 10) = helper_207 * helper_60 + helper_210 * helper_64 + helper_70 * (helper_144 * (helper_179 + helper_213) +
                                                                               helper_211 * helper_91 +
                                                                               helper_212 * helper_57 +
                                                                               helper_214 * helper_69);
    _H(0, 11) = helper_217 * helper_60 + helper_220 * helper_64 + helper_70 * (helper_144 * (helper_154 + helper_223) +
                                                                               helper_221 * helper_91 +
                                                                               helper_222 * helper_57 +
                                                                               helper_224 * helper_69);
    _H(1, 1) = helper_226 * helper_75 + helper_227 * helper_88 +
               helper_70 * (2 * helper_228 * helper_90 + helper_229 * helper_97 + helper_65);
    _H(1, 2) = helper_102 * helper_226 + helper_115 * helper_227 +
               helper_70 * (helper_117 * helper_228 + helper_119 * helper_230 + helper_124 * helper_229);
    _H(1, 3) = helper_129 * helper_226 + helper_227 * helper_55 + helper_70 *
                                                                  (helper_132 * helper_230 + helper_133 * helper_229 +
                                                                   helper_144 * (helper_141 + helper_213) +
                                                                   helper_228 * helper_61);
    _H(1, 4) = helper_138 * helper_226 + helper_227 * helper_89 +
               helper_70 * (helper_130 + helper_140 * helper_230 + helper_145 * helper_229 + helper_228 * helper_87);
    _H(1, 5) = helper_116 * helper_227 + helper_150 * helper_226 + helper_70 * (helper_114 * helper_228 +
                                                                                helper_144 * (helper_231 + helper_233) +
                                                                                helper_152 * helper_230 +
                                                                                helper_156 * helper_229);
    _H(1, 6) = helper_161 * helper_226 + helper_164 * helper_227 + helper_70 * (helper_144 * (helper_142 + helper_179) +
                                                                                helper_166 * helper_228 +
                                                                                helper_167 * helper_230 +
                                                                                helper_168 * helper_229);
    _H(1, 7) = helper_173 * helper_226 + helper_176 * helper_227 +
               helper_70 * (helper_165 + helper_177 * helper_228 + helper_178 * helper_230 + helper_181 * helper_229);
    _H(1, 8) = helper_186 * helper_226 + helper_189 * helper_227 + helper_70 * (helper_144 * (helper_232 + helper_235) +
                                                                                helper_190 * helper_228 +
                                                                                helper_191 * helper_230 +
                                                                                helper_194 * helper_229);
    _H(1, 9) = helper_197 * helper_226 + helper_200 * helper_227 + helper_70 * (helper_144 * (helper_143 + helper_180) +
                                                                                helper_202 * helper_228 +
                                                                                helper_203 * helper_230 +
                                                                                helper_204 * helper_229);
    _H(1, 10) = helper_207 * helper_226 + helper_210 * helper_227 +
                helper_70 * (helper_201 + helper_211 * helper_228 + helper_212 * helper_230 + helper_214 * helper_229);
    _H(1, 11) = helper_217 * helper_226 + helper_220 * helper_227 + helper_70 *
                                                                    (helper_144 * (helper_234 + helper_236) +
                                                                     helper_221 * helper_228 + helper_222 * helper_230 +
                                                                     helper_224 * helper_229);
    _H(2, 2) = helper_102 * helper_238 + helper_115 * helper_239 +
               helper_70 * (2 * helper_117 * helper_240 + helper_124 * helper_241 + helper_65);
    _H(2, 3) = helper_129 * helper_238 + helper_239 * helper_55 + helper_70 *
                                                                  (helper_132 * helper_242 + helper_133 * helper_241 +
                                                                   helper_144 * (helper_154 + helper_192) +
                                                                   helper_240 * helper_61);
    _H(2, 4) = helper_138 * helper_238 + helper_239 * helper_89 + helper_70 * (helper_140 * helper_242 +
                                                                               helper_144 * (helper_232 + helper_236) +
                                                                               helper_145 * helper_241 +
                                                                               helper_240 * helper_87);
    _H(2, 5) = helper_116 * helper_239 + helper_150 * helper_238 +
               helper_70 * (helper_114 * helper_240 + helper_130 + helper_152 * helper_242 + helper_156 * helper_241);
    _H(2, 6) = helper_161 * helper_238 + helper_164 * helper_239 + helper_70 * (helper_144 * (helper_153 + helper_223) +
                                                                                helper_166 * helper_240 +
                                                                                helper_167 * helper_242 +
                                                                                helper_168 * helper_241);
    _H(2, 7) = helper_173 * helper_238 + helper_176 * helper_239 + helper_70 * (helper_144 * (helper_233 + helper_234) +
                                                                                helper_177 * helper_240 +
                                                                                helper_178 * helper_242 +
                                                                                helper_181 * helper_241);
    _H(2, 8) = helper_186 * helper_238 + helper_189 * helper_239 +
               helper_70 * (helper_165 + helper_190 * helper_240 + helper_191 * helper_242 + helper_194 * helper_241);
    _H(2, 9) = helper_197 * helper_238 + helper_200 * helper_239 + helper_70 * (helper_144 * (helper_155 + helper_193) +
                                                                                helper_202 * helper_240 +
                                                                                helper_203 * helper_242 +
                                                                                helper_204 * helper_241);
    _H(2, 10) = helper_207 * helper_238 + helper_210 * helper_239 + helper_70 *
                                                                    (helper_144 * (helper_231 + helper_235) +
                                                                     helper_211 * helper_240 + helper_212 * helper_242 +
                                                                     helper_214 * helper_241);
    _H(2, 11) = helper_217 * helper_238 + helper_220 * helper_239 +
                helper_70 * (helper_201 + helper_221 * helper_240 + helper_222 * helper_242 + helper_224 * helper_241);
    _H(3, 3) = helper_129 * helper_244 + helper_245 * helper_55 +
               helper_70 * (2 * helper_132 * helper_246 + helper_133 * helper_247 + helper_65);
    _H(3, 4) = helper_138 * helper_244 + helper_245 * helper_89 +
               helper_70 * (helper_140 * helper_246 + helper_145 * helper_247 + helper_248 * helper_87);
    _H(3, 5) = helper_116 * helper_245 + helper_150 * helper_244 +
               helper_70 * (helper_114 * helper_248 + helper_152 * helper_246 + helper_156 * helper_247);
    _H(3, 6) = helper_161 * helper_244 + helper_164 * helper_245 +
               helper_70 * (helper_166 * helper_248 + helper_167 * helper_246 + helper_168 * helper_247 + helper_249);
    _H(3, 7) = helper_173 * helper_244 + helper_176 * helper_245 + helper_70 * (helper_144 * (helper_142 + helper_250) +
                                                                                helper_177 * helper_248 +
                                                                                helper_178 * helper_246 +
                                                                                helper_181 * helper_247);
    _H(3, 8) = helper_186 * helper_244 + helper_189 * helper_245 + helper_70 * (helper_144 * (helper_153 + helper_252) +
                                                                                helper_190 * helper_248 +
                                                                                helper_191 * helper_246 +
                                                                                helper_194 * helper_247);
    _H(3, 9) = helper_197 * helper_244 + helper_200 * helper_245 +
               helper_70 * (helper_130 + helper_202 * helper_248 + helper_203 * helper_246 + helper_204 * helper_247);
    _H(3, 10) = helper_207 * helper_244 + helper_210 * helper_245 + helper_70 *
                                                                    (helper_144 * (helper_143 + helper_253) +
                                                                     helper_211 * helper_248 + helper_212 * helper_246 +
                                                                     helper_214 * helper_247);
    _H(3, 11) = helper_217 * helper_244 + helper_220 * helper_245 + helper_70 *
                                                                    (helper_144 * (helper_155 + helper_251) +
                                                                     helper_221 * helper_248 + helper_222 * helper_246 +
                                                                     helper_224 * helper_247);
    _H(4, 4) = helper_138 * helper_255 + helper_256 * helper_89 +
               helper_70 * (helper_145 * helper_258 + 2 * helper_257 * helper_87 + helper_65);
    _H(4, 5) = helper_116 * helper_256 + helper_150 * helper_255 +
               helper_70 * (helper_114 * helper_257 + helper_152 * helper_259 + helper_156 * helper_258);
    _H(4, 6) = helper_161 * helper_255 + helper_164 * helper_256 + helper_70 * (helper_144 * (helper_141 + helper_253) +
                                                                                helper_166 * helper_257 +
                                                                                helper_167 * helper_259 +
                                                                                helper_168 * helper_258);
    _H(4, 7) = helper_173 * helper_255 + helper_176 * helper_256 +
               helper_70 * (helper_177 * helper_257 + helper_178 * helper_259 + helper_181 * helper_258 + helper_249);
    _H(4, 8) = helper_186 * helper_255 + helper_189 * helper_256 + helper_70 * (helper_144 * (helper_233 + helper_260) +
                                                                                helper_190 * helper_257 +
                                                                                helper_191 * helper_259 +
                                                                                helper_194 * helper_258);
    _H(4, 9) = helper_197 * helper_255 + helper_200 * helper_256 + helper_70 * (helper_144 * (helper_213 + helper_250) +
                                                                                helper_202 * helper_257 +
                                                                                helper_203 * helper_259 +
                                                                                helper_204 * helper_258);
    _H(4, 10) = helper_207 * helper_255 + helper_210 * helper_256 +
                helper_70 * (helper_130 + helper_211 * helper_257 + helper_212 * helper_259 + helper_214 * helper_258);
    _H(4, 11) = helper_217 * helper_255 + helper_220 * helper_256 + helper_70 *
                                                                    (helper_144 * (helper_231 + helper_261) +
                                                                     helper_221 * helper_257 + helper_222 * helper_259 +
                                                                     helper_224 * helper_258);
    _H(5, 5) = helper_116 * helper_264 + helper_150 * helper_263 +
               helper_70 * (2 * helper_152 * helper_265 + helper_156 * helper_266 + helper_65);
    _H(5, 6) = helper_161 * helper_263 + helper_164 * helper_264 + helper_70 * (helper_144 * (helper_192 + helper_251) +
                                                                                helper_166 * helper_267 +
                                                                                helper_167 * helper_265 +
                                                                                helper_168 * helper_266);
    _H(5, 7) = helper_173 * helper_263 + helper_176 * helper_264 + helper_70 * (helper_144 * (helper_232 + helper_261) +
                                                                                helper_177 * helper_267 +
                                                                                helper_178 * helper_265 +
                                                                                helper_181 * helper_266);
    _H(5, 8) = helper_186 * helper_263 + helper_189 * helper_264 +
               helper_70 * (helper_190 * helper_267 + helper_191 * helper_265 + helper_194 * helper_266 + helper_249);
    _H(5, 9) = helper_197 * helper_263 + helper_200 * helper_264 + helper_70 * (helper_144 * (helper_154 + helper_252) +
                                                                                helper_202 * helper_267 +
                                                                                helper_203 * helper_265 +
                                                                                helper_204 * helper_266);
    _H(5, 10) = helper_207 * helper_263 + helper_210 * helper_264 + helper_70 *
                                                                    (helper_144 * (helper_236 + helper_260) +
                                                                     helper_211 * helper_267 + helper_212 * helper_265 +
                                                                     helper_214 * helper_266);
    _H(5, 11) = helper_217 * helper_263 + helper_220 * helper_264 +
                helper_70 * (helper_130 + helper_221 * helper_267 + helper_222 * helper_265 + helper_224 * helper_266);
    _H(6, 6) = helper_161 * helper_269 + helper_164 * helper_270 +
               helper_70 * (2 * helper_166 * helper_272 + helper_168 * helper_273 + helper_271);
    _H(6, 7) = helper_173 * helper_269 + helper_176 * helper_270 +
               helper_70 * (helper_177 * helper_272 + helper_178 * helper_274 + helper_181 * helper_273);
    _H(6, 8) = helper_186 * helper_269 + helper_189 * helper_270 +
               helper_70 * (helper_190 * helper_272 + helper_191 * helper_274 + helper_194 * helper_273);
    _H(6, 9) = helper_197 * helper_269 + helper_200 * helper_270 +
               helper_70 * (helper_202 * helper_272 + helper_203 * helper_274 + helper_204 * helper_273 + helper_249);
    _H(6, 10) = helper_207 * helper_269 + helper_210 * helper_270 + helper_70 *
                                                                    (helper_144 * (helper_180 + helper_250) +
                                                                     helper_211 * helper_272 + helper_212 * helper_274 +
                                                                     helper_214 * helper_273);
    _H(6, 11) = helper_217 * helper_269 + helper_220 * helper_270 + helper_70 *
                                                                    (helper_144 * (helper_193 + helper_252) +
                                                                     helper_221 * helper_272 + helper_222 * helper_274 +
                                                                     helper_224 * helper_273);
    _H(7, 7) = helper_173 * helper_276 + helper_176 * helper_277 +
               helper_70 * (2 * helper_178 * helper_278 + helper_181 * helper_279 + helper_271);
    _H(7, 8) = helper_186 * helper_276 + helper_189 * helper_277 +
               helper_70 * (helper_190 * helper_280 + helper_191 * helper_278 + helper_194 * helper_279);
    _H(7, 9) = helper_197 * helper_276 + helper_200 * helper_277 + helper_70 * (helper_144 * (helper_179 + helper_253) +
                                                                                helper_202 * helper_280 +
                                                                                helper_203 * helper_278 +
                                                                                helper_204 * helper_279);
    _H(7, 10) = helper_207 * helper_276 + helper_210 * helper_277 +
                helper_70 * (helper_211 * helper_280 + helper_212 * helper_278 + helper_214 * helper_279 + helper_249);
    _H(7, 11) = helper_217 * helper_276 + helper_220 * helper_277 + helper_70 *
                                                                    (helper_144 * (helper_235 + helper_260) +
                                                                     helper_221 * helper_280 + helper_222 * helper_278 +
                                                                     helper_224 * helper_279);
    _H(8, 8) = helper_186 * helper_282 + helper_189 * helper_283 +
               helper_70 * (2 * helper_191 * helper_284 + helper_194 * helper_285 + helper_271);
    _H(8, 9) = helper_197 * helper_282 + helper_200 * helper_283 + helper_70 * (helper_144 * (helper_223 + helper_251) +
                                                                                helper_202 * helper_286 +
                                                                                helper_203 * helper_284 +
                                                                                helper_204 * helper_285);
    _H(8, 10) = helper_207 * helper_282 + helper_210 * helper_283 + helper_70 *
                                                                    (helper_144 * (helper_234 + helper_261) +
                                                                     helper_211 * helper_286 + helper_212 * helper_284 +
                                                                     helper_214 * helper_285);
    _H(8, 11) = helper_217 * helper_282 + helper_220 * helper_283 +
                helper_70 * (helper_221 * helper_286 + helper_222 * helper_284 + helper_224 * helper_285 + helper_249);
    _H(9, 9) = helper_197 * helper_288 + helper_200 * helper_289 +
               helper_70 * (2 * helper_203 * helper_290 + helper_204 * helper_291 + helper_65);
    _H(9, 10) = helper_207 * helper_288 + helper_210 * helper_289 +
                helper_70 * (helper_211 * helper_292 + helper_212 * helper_290 + helper_214 * helper_291);
    _H(9, 11) = helper_217 * helper_288 + helper_220 * helper_289 +
                helper_70 * (helper_221 * helper_292 + helper_222 * helper_290 + helper_224 * helper_291);
    _H(10, 10) = helper_207 * helper_294 + helper_210 * helper_295 +
                 helper_70 * (2 * helper_211 * helper_296 + helper_214 * helper_297 + helper_65);
    _H(10, 11) = helper_217 * helper_294 + helper_220 * helper_295 +
                 helper_70 * (helper_211 * helper_298 + helper_221 * helper_296 + helper_224 * helper_297);
    _H(11, 11) = helper_217 * helper_299 * helper_59 + helper_220 * helper_299 * helper_63 +
                 helper_70 * (helper_221 * helper_224 * helper_68 + 2 * helper_221 * helper_298 + helper_65);

    for (int i = 1; i < 12; ++i)
      for (int j = 0; j < i; ++j)
        _H(i, j) = _H(j, i);


    // Compute Eigen decomposition H = V D V^T
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(_H);
    const Eigen::MatrixXd &V = solver.eigenvectors();
    Eigen::VectorXd evals = solver.eigenvalues();

    // check eigenvalues against epsilon, for those less that eps (zero)
    const int n = n_unknowns();
    for (int i = 0; i < n; ++i)
    {
      if (evals[i] < 1e-8)
      {
        evals[i] = 1e-8;
      }
    }

    // compute correction matrix M = V * diag(m) * V^T
    _H.noalias() = V * evals.asDiagonal() * V.transpose();
  }
};


//Surface Inverse mean ratio metric energy
class SIMRMEnergy
{
public:
  static int n_unknowns() { return 3; }

  static double eval_f(const double *_x)
  {
    double helper_0 = std::pow(_x[0] - _x[3], 2);
    double helper_1 = std::pow(_x[1] - _x[4], 2);
    double helper_2 = std::pow(_x[2] - _x[5], 2);
    double helper_3 = helper_0 + helper_1 + helper_2;
    double helper_4 = std::sqrt(helper_3);
    double helper_5 = 1.0 / helper_4;
    double helper_6 = -_x[6];
    double helper_7 = std::pow(helper_6 + _x[3], 2);
    double helper_8 = -_x[7];
    double helper_9 = std::pow(helper_8 + _x[4], 2);
    double helper_10 = -_x[8];
    double helper_11 = std::pow(helper_10 + _x[5], 2);
    double helper_12 = std::pow(helper_6 + _x[0], 2);
    double helper_13 = std::pow(helper_8 + _x[1], 2);
    double helper_14 = std::pow(helper_10 + _x[2], 2);
    double helper_15 = helper_12 + helper_13 + helper_14;
    double helper_16 = std::pow(-helper_11 + helper_15 + helper_3 - helper_7 - helper_9, 2) /
                       (4 * helper_0 + 4 * helper_1 + 4 * helper_2);
    double helper_17 =
            1.7320508075688772 * helper_0 + 1.7320508075688772 * helper_1 + 1.7320508075688772 * helper_2;

    return std::log((1.0 / 4.0) * helper_5 * std::pow(helper_15 - helper_16, -0.5) *
                    (2.3094010767585029 * helper_12 + 2.3094010767585029 * helper_13 +
                     2.3094010767585029 * helper_14 -
                     2.3094010767585029 * helper_16 + helper_17 + 0.57735026918962573 * std::pow(helper_4 -
                                                                                                 0.57735026918962573 *
                                                                                                 helper_5 *
                                                                                                 (-1.7320508075688772 *
                                                                                                  helper_11 +
                                                                                                  1.7320508075688772 *
                                                                                                  helper_12 +
                                                                                                  1.7320508075688772 *
                                                                                                  helper_13 +
                                                                                                  1.7320508075688772 *
                                                                                                  helper_14 +
                                                                                                  helper_17 -
                                                                                                  1.7320508075688772 *
                                                                                                  helper_7 -
                                                                                                  1.7320508075688772 *
                                                                                                  helper_9),
                                                                                                 2)));
  }

  static void eval_gradient(const double *_x, Vec3d &_g)
  {
    //log
    double helper_0 = std::pow(_x[0] - _x[3], 2);
    double helper_1 = std::pow(_x[1] - _x[4], 2);
    double helper_2 = std::pow(_x[2] - _x[5], 2);
    double helper_3 = helper_0 + helper_1 + helper_2;
    double helper_4 = std::pow(helper_3, -1.5);
    double helper_5 = 1.0 * _x[0];
    double helper_6 = -helper_5;
    double helper_7 = 1.0 * _x[3];
    double helper_8 = helper_4 * (helper_6 + helper_7);
    double helper_9 = 4 * helper_0 + 4 * helper_1 + 4 * helper_2;
    double helper_10 = 1.0 / helper_9;
    double helper_11 = -_x[6];
    double helper_12 = std::pow(helper_11 + _x[3], 2);
    double helper_13 = -_x[7];
    double helper_14 = std::pow(helper_13 + _x[4], 2);
    double helper_15 = -_x[8];
    double helper_16 = std::pow(helper_15 + _x[5], 2);
    double helper_17 = std::pow(helper_11 + _x[0], 2);
    double helper_18 = std::pow(helper_13 + _x[1], 2);
    double helper_19 = std::pow(helper_15 + _x[2], 2);
    double helper_20 = helper_17 + helper_18 + helper_19;
    double helper_21 = -helper_12 - helper_14 - helper_16 + helper_20 + helper_3;
    double helper_22 = std::pow(helper_21, 2);
    double helper_23 = helper_10 * helper_22;
    double helper_24 = helper_20 - helper_23;
    double helper_25 = std::sqrt(helper_24);
    double helper_26 = 1.0 / helper_25;
    double helper_27 = std::sqrt(helper_3);
    double helper_28 = 1.0 / helper_27;
    double helper_29 =
            1.7320508075688772 * helper_0 + 1.7320508075688772 * helper_1 + 1.7320508075688772 * helper_2;
    double helper_30 =
            -1.7320508075688772 * helper_12 - 1.7320508075688772 * helper_14 - 1.7320508075688772 * helper_16 +
            1.7320508075688772 * helper_17 + 1.7320508075688772 * helper_18 + 1.7320508075688772 * helper_19 +
            helper_29;
    double helper_31 = helper_27 - 0.57735026918962573 * helper_28 * helper_30;
    double helper_32 =
            2.3094010767585029 * helper_17 + 2.3094010767585029 * helper_18 + 2.3094010767585029 * helper_19 -
            2.3094010767585029 * helper_23 + helper_29 + 0.57735026918962573 * std::pow(helper_31, 2);
    double helper_33 = (1.0 / 4.0) * helper_26 * helper_32;
    double helper_34 = 8 * _x[0];
    double helper_35 = helper_22 / std::pow(helper_9, 2);
    double helper_36 = helper_35 * (-helper_34 + 8 * _x[3]);
    double helper_37 = helper_10 * helper_21;
    double helper_38 = helper_37 * (helper_34 - 4 * _x[3] - 4 * _x[6]);
    double helper_39 = (1.0 / 4.0) * helper_28;
    double helper_40 = std::pow(helper_24, -1.5) * helper_32 * helper_39;
    double helper_41 = -3.4641016151377544 * _x[3];
    double helper_42 = 2 * helper_28;
    double helper_43 = 1.1547005383792515 * helper_28;
    double helper_44 = 1.1547005383792515 * helper_30;
    double helper_45 = 0.57735026918962573 * helper_31;
    double helper_46 = helper_26 * helper_39;
    double helper_47 = 4 * helper_25 * helper_27 / helper_32;
    double helper_48 = 1.0 * _x[1];
    double helper_49 = -helper_48;
    double helper_50 = 1.0 * _x[4];
    double helper_51 = helper_49 + helper_50;
    double helper_52 = helper_33 * helper_4;
    double helper_53 = 8 * _x[1];
    double helper_54 = helper_35 * (-helper_53 + 8 * _x[4]);
    double helper_55 = helper_37 * (helper_53 - 4 * _x[4] - 4 * _x[7]);
    double helper_56 = -3.4641016151377544 * _x[4];
    double helper_57 = helper_4 * helper_44;
    double helper_58 = 1.0 * _x[2];
    double helper_59 = -helper_58;
    double helper_60 = 1.0 * _x[5];
    double helper_61 = helper_59 + helper_60;
    double helper_62 = 8 * _x[2];
    double helper_63 = helper_35 * (-helper_62 + 8 * _x[5]);
    double helper_64 = helper_37 * (helper_62 - 4 * _x[5] - 4 * _x[8]);
    double helper_65 = -3.4641016151377544 * _x[5];

    _g[0] = helper_47 *
            (helper_33 * helper_8 + helper_40 * (0.5 * helper_36 + 0.5 * helper_38 + helper_6 + 1.0 * _x[6]) +
             helper_46 *
             (-2.3094010767585029 * helper_36 - 2.3094010767585029 * helper_38 + helper_41 + helper_45 *
                                                                                             (helper_42 *
                                                                                              (helper_5 -
                                                                                               helper_7) -
                                                                                              helper_43 *
                                                                                              (helper_41 +
                                                                                               6.9282032302755088 *
                                                                                               _x[0] -
                                                                                               3.4641016151377544 *
                                                                                               _x[6]) -
                                                                                              helper_44 *
                                                                                              helper_8) +
              8.0829037686547593 * _x[0] - 4.6188021535170058 * _x[6]));
    _g[1] = helper_47 * (helper_40 * (helper_49 + 0.5 * helper_54 + 0.5 * helper_55 + 1.0 * _x[7]) + helper_46 *
                                                                                                     (helper_45 *
                                                                                                      (helper_42 *
                                                                                                       (helper_48 -
                                                                                                        helper_50) -
                                                                                                       helper_43 *
                                                                                                       (helper_56 +
                                                                                                        6.9282032302755088 *
                                                                                                        _x[1] -
                                                                                                        3.4641016151377544 *
                                                                                                        _x[7]) -
                                                                                                       helper_51 *
                                                                                                       helper_57) -
                                                                                                      2.3094010767585029 *
                                                                                                      helper_54 -
                                                                                                      2.3094010767585029 *
                                                                                                      helper_55 +
                                                                                                      helper_56 +
                                                                                                      8.0829037686547593 *
                                                                                                      _x[1] -
                                                                                                      4.6188021535170058 *
                                                                                                      _x[7]) +
                         helper_51 * helper_52);
    _g[2] = helper_47 * (helper_40 * (helper_59 + 0.5 * helper_63 + 0.5 * helper_64 + 1.0 * _x[8]) + helper_46 *
                                                                                                     (helper_45 *
                                                                                                      (helper_42 *
                                                                                                       (helper_58 -
                                                                                                        helper_60) -
                                                                                                       helper_43 *
                                                                                                       (helper_65 +
                                                                                                        6.9282032302755088 *
                                                                                                        _x[2] -
                                                                                                        3.4641016151377544 *
                                                                                                        _x[8]) -
                                                                                                       helper_57 *
                                                                                                       helper_61) -
                                                                                                      2.3094010767585029 *
                                                                                                      helper_63 -
                                                                                                      2.3094010767585029 *
                                                                                                      helper_64 +
                                                                                                      helper_65 +
                                                                                                      8.0829037686547593 *
                                                                                                      _x[2] -
                                                                                                      4.6188021535170058 *
                                                                                                      _x[8]) +
                         helper_52 * helper_61);

  }

  static void eval_hessian(const double *_x, Mat3d &_H)
  {
    _H.setZero();

    //log
    double helper_0 = std::pow(_x[0] - _x[3], 2);
    double helper_1 = std::pow(_x[1] - _x[4], 2);
    double helper_2 = std::pow(_x[2] - _x[5], 2);
    double helper_3 = helper_0 + helper_1 + helper_2;
    double helper_4 = std::sqrt(helper_3);
    double helper_5 = 1.0 / helper_4;
    double helper_6 = 1.0 * _x[0];
    double helper_7 = 1.0 * _x[3];
    double helper_8 = helper_6 - helper_7;
    double helper_9 = helper_5 * helper_8;
    double helper_10 = 4 * helper_0 + 4 * helper_1 + 4 * helper_2;
    double helper_11 = 1.0 / helper_10;
    double helper_12 = -_x[6];
    double helper_13 = std::pow(helper_12 + _x[3], 2);
    double helper_14 = -_x[7];
    double helper_15 = std::pow(helper_14 + _x[4], 2);
    double helper_16 = -_x[8];
    double helper_17 = std::pow(helper_16 + _x[5], 2);
    double helper_18 = std::pow(helper_12 + _x[0], 2);
    double helper_19 = std::pow(helper_14 + _x[1], 2);
    double helper_20 = std::pow(helper_16 + _x[2], 2);
    double helper_21 = helper_18 + helper_19 + helper_20;
    double helper_22 = -helper_13 - helper_15 - helper_17 + helper_21 + helper_3;
    double helper_23 = std::pow(helper_22, 2);
    double helper_24 = helper_11 * helper_23;
    double helper_25 = helper_21 - helper_24;
    double helper_26 = std::sqrt(helper_25);
    double helper_27 = std::pow(helper_3, -1.5);
    double helper_28 = -helper_6;
    double helper_29 = helper_28 + helper_7;
    double helper_30 = helper_27 * helper_29;
    double helper_31 = 1.0 / helper_26;
    double helper_32 =
            1.7320508075688772 * helper_0 + 1.7320508075688772 * helper_1 + 1.7320508075688772 * helper_2;
    double helper_33 =
            -1.7320508075688772 * helper_13 - 1.7320508075688772 * helper_15 - 1.7320508075688772 * helper_17 +
            1.7320508075688772 * helper_18 + 1.7320508075688772 * helper_19 + 1.7320508075688772 * helper_20 +
            helper_32;
    double helper_34 = helper_33 * helper_5;
    double helper_35 = -0.57735026918962573 * helper_34 + helper_4;
    double helper_36 = std::pow(helper_35, 2);
    double helper_37 =
            2.3094010767585029 * helper_18 + 2.3094010767585029 * helper_19 + 2.3094010767585029 * helper_20 -
            2.3094010767585029 * helper_24 + helper_32 + 0.57735026918962573 * helper_36;
    double helper_38 = helper_31 * helper_37;
    double helper_39 = (1.0 / 4.0) * helper_38;
    double helper_40 = std::pow(helper_25, -1.5);
    double helper_41 = 1.0 * _x[6];
    double helper_42 = 8 * _x[0];
    double helper_43 = -helper_42 + 8 * _x[3];
    double helper_44 = std::pow(helper_10, -2);
    double helper_45 = helper_23 * helper_44;
    double helper_46 = helper_43 * helper_45;
    double helper_47 = 0.5 * helper_46;
    double helper_48 = helper_11 * helper_22;
    double helper_49 = helper_42 - 4 * _x[3] - 4 * _x[6];
    double helper_50 = 0.5 * helper_49;
    double helper_51 = helper_48 * helper_50;
    double helper_52 = helper_28 + helper_41 + helper_47 + helper_51;
    double helper_53 = helper_40 * helper_52;
    double helper_54 = (1.0 / 4.0) * helper_5;
    double helper_55 = helper_37 * helper_54;
    double helper_56 = 8.0829037686547593 * _x[0];
    double helper_57 = 3.4641016151377544 * _x[3];
    double helper_58 = -helper_57;
    double helper_59 = 4.6188021535170058 * _x[6];
    double helper_60 = 2.3094010767585029 * helper_46;
    double helper_61 = 2.3094010767585029 * helper_49;
    double helper_62 = helper_48 * helper_61;
    double helper_63 = helper_58 + 6.9282032302755088 * _x[0] - 3.4641016151377544 * _x[6];
    double helper_64 = helper_5 * helper_63;
    double helper_65 = helper_27 * helper_33;
    double helper_66 = 1.1547005383792515 * helper_65;
    double helper_67 = -helper_29 * helper_66 - 1.1547005383792515 * helper_64 + 2 * helper_9;
    double helper_68 = 0.57735026918962573 * helper_35;
    double helper_69 = helper_67 * helper_68;
    double helper_70 = helper_56 + helper_58 - helper_59 - helper_60 - helper_62 + helper_69;
    double helper_71 = helper_31 * helper_54;
    double helper_72 = helper_30 * helper_39 + helper_53 * helper_55 + helper_70 * helper_71;
    double helper_73 = 4 / helper_37;
    double helper_74 = helper_72 * helper_73;
    double helper_75 = helper_26 * helper_74;
    double helper_76 = -helper_41 - helper_47 - helper_51 + helper_6;
    double helper_77 = helper_31 * helper_4;
    double helper_78 = helper_74 * helper_77;
    double helper_79 = -helper_56 + helper_57 + helper_59 + helper_60 + helper_62 - helper_69;
    double helper_80 = helper_26 * helper_4;
    double helper_81 = 0.75000000000000011 * helper_80 /
                       std::pow(0.75 * helper_0 + 0.75 * helper_1 + 0.75 * helper_2 + helper_25 +
                                0.25 * helper_36,
                                2);
    double helper_82 = helper_72 * helper_81;
    double helper_83 = helper_27 * helper_38;
    double helper_84 = -0.25 * helper_83;
    double helper_85 = -3.0 * _x[0];
    double helper_86 = helper_85 + 3.0 * _x[3];
    double helper_87 = std::pow(helper_3, -2.5);
    double helper_88 = helper_29 * helper_87;
    double helper_89 = helper_86 * helper_88;
    double helper_90 = helper_30 * helper_37;
    double helper_91 = (1.0 / 2.0) * helper_70;
    double helper_92 = helper_30 * helper_31;
    double helper_93 = 1.5 * helper_48;
    double helper_94 = 1.5 * helper_46 + helper_49 * helper_93 + helper_85 + 3.0 * _x[6];
    double helper_95 = std::pow(helper_25, -2.5) * helper_55;
    double helper_96 = helper_52 * helper_95;
    double helper_97 = helper_11 * (4 * _x[0] - 2 * _x[3] - 2 * _x[6]);
    double helper_98 = helper_23 / std::pow(helper_10, 3);
    double helper_99 = helper_98 * (-16 * _x[0] + 16 * _x[3]);
    double helper_100 = 0.5 * helper_43;
    double helper_101 = helper_22 * helper_44;
    double helper_102 = helper_101 * helper_43;
    double helper_103 = helper_102 * helper_49;
    double helper_104 = -4.0 * helper_45 + 4.0 * helper_48 - 1.0;
    double helper_105 = helper_40 * helper_55;
    double helper_106 = -0.33333333333333331 * helper_34 + 0.57735026918962573 * helper_4;
    double helper_107 = 2 * helper_8;
    double helper_108 = 1.1547005383792515 * helper_33;
    double helper_109 = -5.9999999999999991 * helper_5 + helper_66;
    double helper_110 = 0.33333333333333331 * helper_65;
    double helper_111 =
            -helper_110 * helper_29 - 0.33333333333333331 * helper_64 + 0.57735026918962573 * helper_9;
    double helper_112 = 2.3094010767585029 * helper_43;
    double helper_113 = 18.475208614068023 * helper_45 - 18.475208614068023 * helper_48 + 8.0829037686547593;
    double helper_114 = helper_73 * helper_80;
    double helper_115 = 1.0 * _x[1];
    double helper_116 = 1.0 * _x[4];
    double helper_117 = helper_115 - helper_116;
    double helper_118 = helper_117 * helper_5;
    double helper_119 = 1.0 * _x[7];
    double helper_120 = 8 * _x[1];
    double helper_121 = -helper_120 + 8 * _x[4];
    double helper_122 = helper_121 * helper_45;
    double helper_123 = 0.5 * helper_122;
    double helper_124 = helper_120 - 4 * _x[4] - 4 * _x[7];
    double helper_125 = helper_124 * helper_48;
    double helper_126 = 0.5 * helper_125;
    double helper_127 = helper_115 - helper_119 - helper_123 - helper_126;
    double helper_128 = 8.0829037686547593 * _x[1];
    double helper_129 = 3.4641016151377544 * _x[4];
    double helper_130 = 4.6188021535170058 * _x[7];
    double helper_131 = 2.3094010767585029 * helper_122;
    double helper_132 = 2.3094010767585029 * helper_125;
    double helper_133 = -helper_129;
    double helper_134 = helper_133 + 6.9282032302755088 * _x[1] - 3.4641016151377544 * _x[7];
    double helper_135 = helper_134 * helper_5;
    double helper_136 = -helper_115;
    double helper_137 = helper_116 + helper_136;
    double helper_138 = 2 * helper_118 - 1.1547005383792515 * helper_135 - helper_137 * helper_66;
    double helper_139 = helper_138 * helper_68;
    double helper_140 = -helper_128 + helper_129 + helper_130 + helper_131 + helper_132 - helper_139;
    double helper_141 = helper_137 * helper_27;
    double helper_142 = -3.0 * _x[1];
    double helper_143 = helper_142 + 3.0 * _x[4];
    double helper_144 = helper_108 * helper_88;
    double helper_145 = 1.1547005383792515 * helper_30;
    double helper_146 = 1.1547005383792515 * helper_63;
    double helper_147 = -helper_134 * helper_145 - helper_141 * helper_146;
    double helper_148 =
            -helper_110 * helper_137 + 0.57735026918962573 * helper_118 - 0.33333333333333331 * helper_135;
    double helper_149 = helper_11 * (4 * _x[1] - 2 * _x[4] - 2 * _x[7]);
    double helper_150 = helper_98 * (-16 * _x[1] + 16 * _x[4]);
    double helper_151 = helper_101 * helper_112;
    double helper_152 = helper_101 * helper_121;
    double helper_153 = -helper_124 * helper_151 - helper_152 * helper_61;
    double helper_154 = 0.5 * helper_102;
    double helper_155 = helper_124 * helper_154 + helper_152 * helper_50;
    double helper_156 = 1.5 * helper_122 + helper_124 * helper_93 + helper_142 + 3.0 * _x[7];
    double helper_157 = helper_39 * helper_88;
    double helper_158 = helper_128 - helper_130 - helper_131 - helper_132 + helper_133 + helper_139;
    double helper_159 = helper_53 * helper_54;
    double helper_160 = helper_119 + helper_123 + helper_126 + helper_136;
    double helper_161 = helper_160 * helper_40;
    double helper_162 = helper_54 * helper_70;
    double helper_163 = (1.0 / 4.0) * helper_92;
    double helper_164 = helper_141 * helper_31;
    double helper_165 = (1.0 / 4.0) * helper_70;
    double helper_166 = (1.0 / 4.0) * helper_90;
    double helper_167 = helper_141 * helper_37;
    double helper_168 = (1.0 / 4.0) * helper_53;
    double helper_169 =
            helper_158 * helper_159 + helper_158 * helper_163 + helper_161 * helper_162 +
            helper_161 * helper_166 +
            helper_164 * helper_165 + helper_167 * helper_168;
    double helper_170 = 1.0 * _x[2];
    double helper_171 = 1.0 * _x[5];
    double helper_172 = helper_170 - helper_171;
    double helper_173 = helper_172 * helper_5;
    double helper_174 = 1.0 * _x[8];
    double helper_175 = 8 * _x[2];
    double helper_176 = -helper_175 + 8 * _x[5];
    double helper_177 = helper_176 * helper_45;
    double helper_178 = 0.5 * helper_177;
    double helper_179 = helper_175 - 4 * _x[5] - 4 * _x[8];
    double helper_180 = helper_179 * helper_48;
    double helper_181 = 0.5 * helper_180;
    double helper_182 = helper_170 - helper_174 - helper_178 - helper_181;
    double helper_183 = 8.0829037686547593 * _x[2];
    double helper_184 = 3.4641016151377544 * _x[5];
    double helper_185 = 4.6188021535170058 * _x[8];
    double helper_186 = 2.3094010767585029 * helper_177;
    double helper_187 = 2.3094010767585029 * helper_180;
    double helper_188 = -helper_184;
    double helper_189 = helper_188 + 6.9282032302755088 * _x[2] - 3.4641016151377544 * _x[8];
    double helper_190 = helper_189 * helper_5;
    double helper_191 = -helper_170;
    double helper_192 = helper_171 + helper_191;
    double helper_193 = 2 * helper_173 - 1.1547005383792515 * helper_190 - helper_192 * helper_66;
    double helper_194 = helper_193 * helper_68;
    double helper_195 = -helper_183 + helper_184 + helper_185 + helper_186 + helper_187 - helper_194;
    double helper_196 = helper_192 * helper_27;
    double helper_197 = -3.0 * _x[2];
    double helper_198 = helper_197 + 3.0 * _x[5];
    double helper_199 = -helper_145 * helper_189 - helper_146 * helper_196;
    double helper_200 =
            -helper_110 * helper_192 + 0.57735026918962573 * helper_173 - 0.33333333333333331 * helper_190;
    double helper_201 = helper_11 * (4 * _x[2] - 2 * _x[5] - 2 * _x[8]);
    double helper_202 = helper_98 * (-16 * _x[2] + 16 * _x[5]);
    double helper_203 = helper_101 * helper_176;
    double helper_204 = -helper_151 * helper_179 - helper_203 * helper_61;
    double helper_205 = helper_154 * helper_179 + helper_203 * helper_50;
    double helper_206 = 1.5 * helper_177 + helper_179 * helper_93 + helper_197 + 3.0 * _x[8];
    double helper_207 = helper_183 - helper_185 - helper_186 - helper_187 + helper_188 + helper_194;
    double helper_208 = helper_174 + helper_178 + helper_181 + helper_191;
    double helper_209 = helper_208 * helper_40;
    double helper_210 = helper_196 * helper_31;
    double helper_211 = helper_196 * helper_37;
    double helper_212 =
            helper_159 * helper_207 + helper_162 * helper_209 + helper_163 * helper_207 +
            helper_165 * helper_210 +
            helper_166 * helper_209 + helper_168 * helper_211;
    double helper_213 = (1.0 / 4.0) * helper_83;
    double helper_214 = helper_105 * helper_160 + helper_137 * helper_213 + helper_158 * helper_71;
    double helper_215 = helper_214 * helper_73;
    double helper_216 = helper_215 * helper_26;
    double helper_217 = helper_215 * helper_77;
    double helper_218 = helper_214 * helper_81;
    double helper_219 = 2 * helper_117;
    double helper_220 = helper_108 * helper_87;
    double helper_221 = helper_137 * helper_220;
    double helper_222 = helper_124 * helper_97;
    double helper_223 = helper_121 * helper_99;
    double helper_224 = helper_160 * helper_95;
    double helper_225 = helper_39 * helper_87;
    double helper_226 = helper_137 * helper_225;
    double helper_227 = (1.0 / 2.0) * helper_158;
    double helper_228 = helper_124 * helper_149;
    double helper_229 = helper_121 * helper_150;
    double helper_230 = helper_124 * helper_152;
    double helper_231 =
            -1.1547005383792515 * helper_134 * helper_196 - 1.1547005383792515 * helper_141 * helper_189;
    double helper_232 = helper_124 * helper_201;
    double helper_233 = helper_121 * helper_202;
    double helper_234 = helper_152 * helper_179;
    double helper_235 = helper_124 * helper_203;
    double helper_236 = -2.3094010767585029 * helper_234 - 2.3094010767585029 * helper_235;
    double helper_237 = 0.5 * helper_234 + 0.5 * helper_235;
    double helper_238 = helper_158 * helper_209 * helper_54 + (1.0 / 4.0) * helper_158 * helper_210 +
                        helper_161 * helper_207 * helper_54 + (1.0 / 4.0) * helper_161 * helper_211 +
                        (1.0 / 4.0) * helper_164 * helper_207 + (1.0 / 4.0) * helper_167 * helper_209;
    double helper_239 = helper_105 * helper_208 + helper_192 * helper_213 + helper_207 * helper_71;
    double helper_240 = helper_239 * helper_73;
    double helper_241 = helper_240 * helper_26;
    double helper_242 = helper_240 * helper_77;
    double helper_243 = helper_239 * helper_81;
    double helper_244 = 2 * helper_172;
    double helper_245 = helper_192 * helper_220;
    double helper_246 = helper_179 * helper_97;
    double helper_247 = helper_176 * helper_99;
    double helper_248 = helper_208 * helper_95;
    double helper_249 = helper_192 * helper_225;
    double helper_250 = helper_149 * helper_179;
    double helper_251 = helper_150 * helper_176;
    double helper_252 = (1.0 / 2.0) * helper_207;
    double helper_253 = helper_179 * helper_201;
    double helper_254 = helper_176 * helper_202;
    double helper_255 = helper_179 * helper_203;

    _H(0, 0) = helper_114 *
               (helper_105 * (helper_100 * helper_99 + 1.0 * helper_103 + helper_104 + helper_50 * helper_97) +
                helper_39 * helper_89 + helper_5 * helper_53 * helper_91 + (1.0 / 2.0) * helper_53 * helper_90 +
                helper_71 * (-4.6188021535170058 * helper_103 + helper_106 *
                                                                (helper_107 * helper_30 -
                                                                 helper_108 * helper_89 +
                                                                 helper_109 -
                                                                 2.3094010767585029 * helper_30 * helper_63) +
                             helper_111 * helper_67 - helper_112 * helper_99 + helper_113 -
                             helper_61 * helper_97) +
                helper_84 + helper_91 * helper_92 + helper_94 * helper_96) + helper_75 * helper_9 +
               helper_76 * helper_78 + helper_79 * helper_82;
    _H(0, 1) = helper_114 * (helper_105 * (helper_100 * helper_150 + helper_149 * helper_50 + helper_155) +
                             helper_143 * helper_157 + helper_156 * helper_96 + helper_169 + helper_71 *
                                                                                             (helper_106 *
                                                                                              (helper_107 *
                                                                                               helper_141 -
                                                                                               helper_143 *
                                                                                               helper_144 +
                                                                                               helper_147) -
                                                                                              helper_112 *
                                                                                              helper_150 +
                                                                                              helper_148 *
                                                                                              helper_67 -
                                                                                              helper_149 *
                                                                                              helper_61 +
                                                                                              helper_153)) +
               helper_118 * helper_75 + helper_127 * helper_78 + helper_140 * helper_82;
    _H(0, 2) = helper_114 * (helper_105 * (helper_100 * helper_202 + helper_201 * helper_50 + helper_205) +
                             helper_157 * helper_198 + helper_206 * helper_96 + helper_212 + helper_71 *
                                                                                             (helper_106 *
                                                                                              (helper_107 *
                                                                                               helper_196 -
                                                                                               helper_144 *
                                                                                               helper_198 +
                                                                                               helper_199) -
                                                                                              helper_112 *
                                                                                              helper_202 +
                                                                                              helper_200 *
                                                                                              helper_67 -
                                                                                              helper_201 *
                                                                                              helper_61 +
                                                                                              helper_204)) +
               helper_173 * helper_75 + helper_182 * helper_78 + helper_195 * helper_82;
    _H(1, 0) = helper_114 * (helper_105 * (helper_155 + 0.5 * helper_222 + 0.5 * helper_223) + helper_169 +
                             helper_224 * helper_94 + helper_226 * helper_86 + helper_71 * (helper_106 *
                                                                                            (helper_147 +
                                                                                             helper_219 *
                                                                                             helper_30 -
                                                                                             helper_221 *
                                                                                             helper_86) +
                                                                                            helper_111 *
                                                                                            helper_138 +
                                                                                            helper_153 -
                                                                                            2.3094010767585029 *
                                                                                            helper_222 -
                                                                                            2.3094010767585029 *
                                                                                            helper_223)) +
               helper_216 * helper_9 + helper_217 * helper_76 + helper_218 * helper_79;
    _H(1, 1) =
            helper_114 * (helper_105 * (helper_104 + 0.5 * helper_228 + 0.5 * helper_229 + 1.0 * helper_230) +
                          helper_143 * helper_226 + helper_156 * helper_224 +
                          (1.0 / 2.0) * helper_161 * helper_167 + helper_161 * helper_227 * helper_5 +
                          helper_164 * helper_227 + helper_71 * (helper_106 * (helper_109 -
                                                                               2.3094010767585029 * helper_134 *
                                                                               helper_141 +
                                                                               helper_141 * helper_219 -
                                                                               helper_143 * helper_221) +
                                                                 helper_113 + helper_138 * helper_148 -
                                                                 2.3094010767585029 * helper_228 -
                                                                 2.3094010767585029 * helper_229 -
                                                                 4.6188021535170058 * helper_230) + helper_84) +
            helper_118 * helper_216 + helper_127 * helper_217 + helper_140 * helper_218;
    _H(1, 2) = helper_114 *
               (helper_105 * (0.5 * helper_232 + 0.5 * helper_233 + helper_237) + helper_198 * helper_226 +
                helper_206 * helper_224 + helper_238 + helper_71 * (helper_106 * (helper_196 * helper_219 -
                                                                                  helper_198 * helper_221 +
                                                                                  helper_231) +
                                                                    helper_138 * helper_200 -
                                                                    2.3094010767585029 * helper_232 -
                                                                    2.3094010767585029 * helper_233 +
                                                                    helper_236)) +
               helper_173 * helper_216 + helper_182 * helper_217 + helper_195 * helper_218;
    _H(2, 0) = helper_114 * (helper_105 * (helper_205 + 0.5 * helper_246 + 0.5 * helper_247) + helper_212 +
                             helper_248 * helper_94 + helper_249 * helper_86 + helper_71 * (helper_106 *
                                                                                            (helper_199 +
                                                                                             helper_244 *
                                                                                             helper_30 -
                                                                                             helper_245 *
                                                                                             helper_86) +
                                                                                            helper_111 *
                                                                                            helper_193 +
                                                                                            helper_204 -
                                                                                            2.3094010767585029 *
                                                                                            helper_246 -
                                                                                            2.3094010767585029 *
                                                                                            helper_247)) +
               helper_241 * helper_9 + helper_242 * helper_76 + helper_243 * helper_79;
    _H(2, 1) = helper_114 *
               (helper_105 * (helper_237 + 0.5 * helper_250 + 0.5 * helper_251) + helper_143 * helper_249 +
                helper_156 * helper_248 + helper_238 + helper_71 * (helper_106 * (helper_141 * helper_244 -
                                                                                  helper_143 * helper_245 +
                                                                                  helper_231) +
                                                                    helper_148 * helper_193 + helper_236 -
                                                                    2.3094010767585029 * helper_250 -
                                                                    2.3094010767585029 * helper_251)) +
               helper_118 * helper_241 + helper_127 * helper_242 + helper_140 * helper_243;
    _H(2, 2) =
            helper_114 * (helper_105 * (helper_104 + 0.5 * helper_253 + 0.5 * helper_254 + 1.0 * helper_255) +
                          helper_198 * helper_249 + helper_206 * helper_248 +
                          (1.0 / 2.0) * helper_209 * helper_211 + helper_209 * helper_252 * helper_5 +
                          helper_210 * helper_252 + helper_71 * (helper_106 * (helper_109 -
                                                                               2.3094010767585029 * helper_189 *
                                                                               helper_196 +
                                                                               helper_196 * helper_244 -
                                                                               helper_198 * helper_245) +
                                                                 helper_113 + helper_193 * helper_200 -
                                                                 2.3094010767585029 * helper_253 -
                                                                 2.3094010767585029 * helper_254 -
                                                                 4.6188021535170058 * helper_255) + helper_84) +
            helper_173 * helper_241 + helper_182 * helper_242 + helper_195 * helper_243;


    // Compute Eigen decomposition H = V D V^T
    Eigen::SelfAdjointEigenSolver<Mat3d> solver(_H);
    const Mat3d &V = solver.eigenvectors();
    Vec3d evals = solver.eigenvalues();

    // check eigenvalues against epsilon, for those less that eps (zero)
    const int n = n_unknowns();
    for (int i = 0; i < n; ++i)
    {
      if (evals[i] < 1e-8)
      {
        evals[i] = 1e-8;
      }
    }

    // compute correction matrix M = V * diag(m) * V^T
    _H.noalias() = V * evals.asDiagonal() * V.transpose();
  }
};
}
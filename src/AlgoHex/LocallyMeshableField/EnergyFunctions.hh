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

  static double eval_f(const Vec3d &_n, const double *_x);

  static void eval_gradient(const Vec3d &_n, const double *_x, double *_g);

  static void eval_hessian(const Vec3d &_n, const double *_x, Eigen::MatrixXd &_H);
};


class ComplexEdgeLength
{
public:
  static int n_unknowns() { return 6; }

  static double eval_f(const double *_x);

  static void eval_gradient(const double *_x, double *_g);

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H);
};


class CurvatureSmooth
{
public:
  static int n_unknowns() { return 9; }

  static double eval_f(const double *_x);

  static void eval_gradient(const double *_x, double *_g);

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H);
};


class RepulsionEnergy
{
public:
  static int n_unknowns() { return 3; }

  static double eval_f(const double *_x, const Vec3d &_target_pt);

  static void eval_gradient(const double *_x, const Vec3d &_target_pt, double *_g);

  static void eval_hessian(Mat3d &_H);
};


//Inverse mean ratio metric energy
class IMRMEnergy
{
public:
  static int n_unknowns() { return 3; }

  static double eval_f(const double *_x);

  static void eval_gradient(const double *_x, Vec3d &_g);

  static void eval_hessian(const double *_x, Mat3d &_H);
};


//Two free vertices
class IMRMEnergy2
{
public:
  static int n_unknowns() { return 6; }

  static double eval_f(const double *_x);

  static void eval_gradient(const double *_x, double *_g);

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H);
};


//Three free vertices
class IMRMEnergy3
{
public:
  static int n_unknowns() { return 9; }

  static double eval_f(const double *_x);

  static void eval_gradient(const double *_x, double *_g);

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H);
};


//Four vertices
class IMRMEnergy4
{
public:
  static int n_unknowns() { return 12; }

  static double eval_f(const double *_x);

  static void eval_gradient(const double *_x, double *_g);

  static void eval_hessian(const double *_x, Eigen::MatrixXd &_H);
};


//Surface Inverse mean ratio metric energy
class SIMRMEnergy
{
public:
  static int n_unknowns() { return 3; }

  static double eval_f(const double *_x);

  static void eval_gradient(const double *_x, Vec3d &_g);

  static void eval_hessian(const double *_x, Mat3d &_H);
};
}
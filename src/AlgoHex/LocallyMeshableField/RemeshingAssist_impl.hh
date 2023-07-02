/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define REMESHINGASSIST_C

#include "RemeshingAssist.hh"

namespace AlgoHex
{
//compare the quality after remeshing with before
template<class MeshT>
bool
RemeshingAssist<MeshT>::is_quality_improved(const std::vector<std::vector<Point>> &_new_cells_points,
                                            const std::vector<std::vector<Point>> &_old_cells_points,
                                            bool _compare_max_energy)
{
  std::vector<double> new_energy, old_energy;
  new_energy.reserve(_new_cells_points.size());
  old_energy.reserve(_old_cells_points.size());
  for (const auto &cvs: _new_cells_points)
    new_energy.push_back(compute_IMRM_energy(cvs));
  for (const auto &cvs: _old_cells_points)
    old_energy.push_back(compute_IMRM_energy(cvs));

  if (_compare_max_energy)
  {
    if (new_energy.empty())
      return false;
    else
      return *std::max_element(new_energy.begin(), new_energy.end())
             < *std::max_element(old_energy.begin(), old_energy.end());
  }
  else
  {//use the sum of energy
    double new_sum = 0., old_sum = 0.;
    int i = 0;
    for (const auto &cvs: _new_cells_points)
    {
      auto volume = (cvs[1] - cvs[0]) % (cvs[2] - cvs[0]) | (cvs[3] - cvs[0]);
      new_sum += new_energy[i] * volume;
      i++;
    }
    i = 0;
    for (const auto &cvs: _old_cells_points)
    {
      auto volume = (cvs[1] - cvs[0]) % (cvs[2] - cvs[0]) | (cvs[3] - cvs[0]);
      old_sum += old_energy[i] * volume;
      i++;
    }

    return new_sum < old_sum;
  }
}


//compute the IMRM energy
template<class MeshT>
double
RemeshingAssist<MeshT>::compute_IMRM_energy(const std::vector<Point> &_points)
{
  auto energy = 0.;

  if (!is_tet_valid(_points))
  {
    return std::numeric_limits<double>::infinity();
  }
  else
  {
    double x[12] = {0.};
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++)
        x[i * 3 + j] = _points[i][j];

    energy = comformalIMRMEnergy(x);
    //handle numerical error
    if (std::isinf(energy) || std::isnan(energy) || energy <= 0)
      energy = std::numeric_limits<double>::infinity();
  }

  return energy;
}


//compute the Gradient
template<class MeshT>
void
RemeshingAssist<MeshT>::compute_IMRM_gradient(const std::vector<Point> &_points, double *_gradient)
{
  double x[12] = {0.};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      x[i * 3 + j] = _points[i][j];

  comformalIMRMGradient(x, _gradient);
}

//compute the hessian
template<class MeshT>
void
RemeshingAssist<MeshT>::compute_IMRM_hessian(const std::vector<Point> &_points, double *_hessian)
{
  double x[12] = {0.};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      x[i * 3 + j] = _points[i][j];

  comformalIMRMHessian(x, _hessian);
}


//check the orientation of the tetrahedron
template<class MeshT>
bool
RemeshingAssist<MeshT>::is_tet_valid(const std::vector<Point> &_points)
{
  assert(_points.size() == 4);

  if (((_points[1] - _points[0]) % (_points[2] - _points[0]) | (_points[3] - _points[0])) <= 0.)
    return false;

  return true;
}


template<class MeshT>
int
RemeshingAssist<MeshT>::boundary_edge_type(double _angle)
{
  int type = 0;
  if ((_angle <= 1.25 * M_PI && _angle >= 0.75 * M_PI))
    type = 0;
  else if (_angle >= 0 && _angle < 0.75 * M_PI)
    type = -1;
  else if (_angle > 1.25 * M_PI)
    type = 1;

  return type;
}

template<class MeshT>
void
RemeshingAssist<MeshT>::sample_points(const Point &_p0, const Point &_p1, const Point &_p2, const int _n, const int _N,
                                      std::vector<Point> &_points)
{
  if (_n == _N)
  {
    _points.push_back((_p0 + _p1 + _p2) / 3);
    return;
  }
  Point p3 = (_p0 + _p1) / 2, p4 = (_p1 + _p2) / 2, p5 = (_p0 + _p2) / 2;
  sample_points(_p0, p3, p5, _n + 1, _N, _points);
  sample_points(_p1, p4, p3, _n + 1, _N, _points);
  sample_points(_p2, p5, p4, _n + 1, _N, _points);
  sample_points(p3, p4, p5, _n + 1, _N, _points);
}

template<class MeshT>
double
RemeshingAssist<MeshT>::comformalIMRMEnergy(const double *_x)
{
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


template<class MeshT>
void
RemeshingAssist<MeshT>::comformalIMRMGradient(const double *_x, double *_g)
{
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

template<class MeshT>
void
RemeshingAssist<MeshT>::comformalIMRMHessian(const double *_x, double *_h)
{
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

  _h[0] = helper_57 * (-helper_0 + helper_1 + helper_2 + helper_3) +
          helper_60 * (-helper_43 - helper_46 + helper_47 + helper_49 + helper_50 + helper_52) +
          helper_64 * (2 * helper_4 * helper_54 + helper_61 + helper_63 * (2.3570226039551585 * helper_10 -
                                                                           1.6666666666666665 * helper_21 -
                                                                           1.6666666666666665 * helper_22 +
                                                                           2.3570226039551585 * helper_45 -
                                                                           1.6666666666666665 * helper_48 -
                                                                           1.6666666666666665 *
                                                                           helper_51));
  _h[1] = helper_57 * helper_69 + helper_60 * helper_81 +
          helper_64 * (helper_54 * helper_84 + helper_63 * helper_86 + helper_82 * helper_83);
  _h[2] = helper_103 * helper_60 + helper_57 * helper_91 +
          helper_64 * (helper_104 * helper_83 + helper_105 * helper_54 + helper_108 * helper_63);
  _h[3] = _h[1];
  _h[4] = helper_111 * helper_69 + helper_112 * helper_81 +
          helper_64 * (2 * helper_109 * helper_84 + helper_113 * helper_86 + helper_61);
  _h[5] = helper_103 * helper_112 + helper_111 * helper_91 +
          helper_64 * (helper_105 * helper_109 + helper_108 * helper_113 + helper_114 * helper_84);
  _h[6] = _h[2];
  _h[7] = _h[5];
  _h[8] = helper_103 * helper_115 * helper_59 + helper_115 * helper_56 * helper_91 +
          helper_64 * (helper_104 * helper_108 * helper_62 + 2 * helper_105 * helper_114 + helper_61);
}

}
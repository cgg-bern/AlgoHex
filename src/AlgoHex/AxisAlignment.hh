/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once


//== INCLUDES =================================================================

#include <AlgoHex/TypeDef.hh>
#include <cassert>

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== ENUM DEFINITION =========================================================



/** \enum AxisAlignment  AxisAlignment.hh

    Brief Description.

    A more elaborate description follows.
*/

enum AxisAlignment
{
  Align_NONE = -1, Align_X = 0, Align_MX = 1, Align_Y = 2, Align_MY = 3, Align_Z = 4, Align_MZ = 5
};

class AxisAlignmentHelpers
{
public:
  static Eigen::Vector3d vector(const AxisAlignment _a)
  {
    switch (_a)
    {
      case Align_NONE: return Eigen::Vector3d(0, 0, 0);
      case Align_X   : return Eigen::Vector3d(1, 0, 0);
      case Align_MX  : return Eigen::Vector3d(-1, 0, 0);
      case Align_Y   : return Eigen::Vector3d(0, 1, 0);
      case Align_MY  : return Eigen::Vector3d(0, -1, 0);
      case Align_Z   : return Eigen::Vector3d(0, 0, 1);
      case Align_MZ  : return Eigen::Vector3d(0, 0, -1);
      default: assert(false);
        return {};
    }
  }

  static AxisAlignment get_dominant_axis(const Eigen::Vector3d &_v)
  {
    int idx(0);
    if (std::abs(_v[1]) > std::abs(_v[idx]))
      idx = 1;
    if (std::abs(_v[2]) > std::abs(_v[idx]))
      idx = 2;

    if (_v[idx] >= 0.0)
      idx = 2 * idx; // positive
    else
      idx = 2 * idx + 1; // negative

    return AxisAlignment(idx);
  }

  static bool is_positive_axis(const AxisAlignment _a)
  {
    return ((int(_a) % 2) == 0);
  }

  static Eigen::Vector3d quaternion_vector(const Quaternion &_q, const AxisAlignment _a)
  {
    // get corresponding unit vector
    int vnr = _a / 2;
    bool negative = _a % 2;

    Vec3d v(0 == vnr, 1 == vnr, 2 == vnr);

    if (negative)
      v = -v;

    return _q * v;
  }


  template<typename Vec3>
  static std::pair<double, AxisAlignment>
  closest_axis(const Quaternion &_q, const Vec3 &_n);
};

extern template ALGOHEX_EXPORT std::pair<double, AxisAlignment> AxisAlignmentHelpers::closest_axis<OVM::Geometry::Vec3d>(const Quaternion &, const OVM::Geometry::Vec3d &);

extern template ALGOHEX_EXPORT std::pair<double, AxisAlignment> AxisAlignmentHelpers::closest_axis<Eigen::Vector3d>(const Quaternion &, const Eigen::Vector3d &);

static AxisAlignment unoriented_axis(const AxisAlignment _a)
{
  switch (_a)
  {
    case Align_MX: return AxisAlignment::Align_X;
    case Align_MY: return AxisAlignment::Align_Y;
    case Align_MZ: return AxisAlignment::Align_Z;
    default: return _a;
  }
}

static AxisAlignment reverse_axis(const AxisAlignment _a)
{
  switch (_a)
  {
    case Align_MX: return AxisAlignment::Align_X;
    case Align_MY: return AxisAlignment::Align_Y;
    case Align_MZ: return AxisAlignment::Align_Z;
    case Align_X : return AxisAlignment::Align_MX;
    case Align_Y : return AxisAlignment::Align_MY;
    case Align_Z : return AxisAlignment::Align_MZ;
    default: return _a;
  }
}

//right-hand rule
//x, -x, y, -y, z, -z: 0, 1, 2, 3, 4, 5
static AxisAlignment the_third_axis(const AxisAlignment _ax0, const AxisAlignment _ax1)
{
  int x0 = _ax0 / 2;
  int x1 = _ax1 / 2;

  if (x0 == x1)
  {
    std::cout << "Warning: input axes are parallel!" << std::endl;
    return _ax0;
  }

  int x2 = 0;
  if (x0 + x1 == 1)
    x2 = 2;
  else if (x0 + x1 == 2)
    x2 = 1;
  else if (x0 + x1 == 3)
    x2 = 0;

  int ax2 = 2 * x2;

  int prod = (_ax1 - _ax0) * (ax2 - _ax1) * (_ax0 - ax2);
  if (prod > 0)
    prod = 1;
  else
    prod = 0;
  prod = (_ax0 % 2 + _ax1 % 2 + prod) % 2;

  if (prod == 1)
    ax2 = ax2 + 1;

  return (AxisAlignment) ax2;
}


//=============================================================================
} // namespace AlgoHex
//=============================================================================

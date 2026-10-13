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

public:
  static Eigen::Vector3d vector(const AxisAlignment _a);

  static AxisAlignment get_dominant_axis(const Eigen::Vector3d &_v);

  static bool is_positive_axis(const AxisAlignment _a);

  static Eigen::Vector3d quaternion_vector(const Quaternion &_q, const AxisAlignment _a);

  template<typename Vec3>
  static std::pair<double, AxisAlignment>
  closest_axis(const Quaternion &_q, const Vec3 &_n);
};

extern template ALGOHEX_EXPORT std::pair<double, AxisAlignment> AxisAlignmentHelpers::closest_axis<OVM::Geometry::Vec3d>(const Quaternion &, const OVM::Geometry::Vec3d &);
extern template ALGOHEX_EXPORT std::pair<double, AxisAlignment> AxisAlignmentHelpers::closest_axis<Eigen::Vector3d>(const Quaternion &, const Eigen::Vector3d &);

AxisAlignment unoriented_axis(const AxisAlignment _a);

AxisAlignment reverse_axis(const AxisAlignment _a);

AxisAlignment the_third_axis(const AxisAlignment _ax0, const AxisAlignment _ax1);

//=============================================================================
} // namespace AlgoHex
//=============================================================================

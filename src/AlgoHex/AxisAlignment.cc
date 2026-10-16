/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "AxisAlignment.hh"

namespace AlgoHex
{

template<typename Vec3>
std::pair<double, AxisAlignment>
AxisAlignmentHelpers::closest_axis(const Quaternion &_q, const Vec3 &_n)
{
  // get local frame
  auto m = _q.toRotationMatrix();

  double dpx = _n.dot(Vec3(m(0, 0), m(1, 0), m(2, 0)));
  double dpy = _n.dot(Vec3(m(0, 1), m(1, 1), m(2, 1)));
  double dpz = _n.dot(Vec3(m(0, 2), m(1, 2), m(2, 2)));

  std::vector<std::pair<double, AxisAlignment>> v;
  v.reserve(6);
  v.emplace_back(dpx, AxisAlignment::Align_X);
  v.emplace_back(-dpx, AxisAlignment::Align_MX);
  v.emplace_back(dpy, AxisAlignment::Align_Y);
  v.emplace_back(-dpy, AxisAlignment::Align_MY);
  v.emplace_back(dpz, AxisAlignment::Align_Z);
  v.emplace_back(-dpz, AxisAlignment::Align_MZ);

  return *std::max_element(v.begin(), v.end(),
                           [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
}

template std::pair<double, AxisAlignment> ALGOHEX_EXPORT AxisAlignmentHelpers::closest_axis<OVM::Geometry::Vec3d>(const Quaternion &, const OVM::Geometry::Vec3d &);
template std::pair<double, AxisAlignment> ALGOHEX_EXPORT AxisAlignmentHelpers::closest_axis<Eigen::Vector3d>(const Quaternion &, const Eigen::Vector3d &);

} // namespace AlgoHex

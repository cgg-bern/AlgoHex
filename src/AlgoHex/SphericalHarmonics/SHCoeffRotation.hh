/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "SHCoeffs.hh"

#include <Eigen/Core>

namespace AlgoHex
{
namespace SHRotation
{

extern const Mat9d RB_x_pi2;
extern const Mat9d RB_x_minus_pi2;
extern const Mat9d RB_y_pi2;
extern const Mat9d RB_y_minus_pi2;

Mat9d RB_z(double radians);

Mat9d RB_x(double radians);

Mat9d RB_y(double radians);

Vec9d rot_z(double radians, Vec9d const &);

Mat9d RB_from_quat(const Quaternion &q);

}
}

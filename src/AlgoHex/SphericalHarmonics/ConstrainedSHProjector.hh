/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "SHProjectorInterface.hh"
#include <AlgoHex/Config/Export.hh>

namespace AlgoHex
{

class ALGOHEX_EXPORT ConstrainedSHProjector
{
public:
  // projector to closest frame that contains `v`
  ConstrainedSHProjector(const Vec3d &v);

  SHProjectionResult project(const SHCoeffs &sh_coeffs);

private:
  Quaternion v_to_z_;
  SHCoeffs r1_, r9_;
  //Vec3d v;
};

} // namespace AlgoHex

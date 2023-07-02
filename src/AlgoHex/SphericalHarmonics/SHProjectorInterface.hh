/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include "SHCoeffs.hh"

namespace AlgoHex
{

struct SHProjectionResult
{
  SHCoeffs shc;
  Quaternion q;
  size_t iterations; // number of iterations used in projection
};

class SHProjectorInterface
{
  virtual SHProjectionResult project(const SHCoeffs &sh_coeffs) = 0;
};


}

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <AlgoHex/TypeDef.hh>
#include "SHCoeffs.hh"


namespace AlgoHex
{

Quaternion ALGOHEX_EXPORT from_projected_shc(const SHCoeffs &shc);

} // namespace AlgoHex


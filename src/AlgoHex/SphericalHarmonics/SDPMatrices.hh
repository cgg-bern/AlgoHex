/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <Eigen/Core>
#include <AlgoHex/Config/Export.hh>

namespace AlgoHex
{

// SdpA and SdpB as computed in
// https://github.com/dpa1mer/arff/blob/2dc37b560840af5052b34ec10b7bf9ba442385ae/src/mbo/OctaMBO.m

extern ALGOHEX_EXPORT Eigen::Map<const Eigen::Matrix<double, 16, 100> > sdpA;
extern ALGOHEX_EXPORT const Eigen::Matrix<double, 16, 1> sdpB;

} // namespace AlgoHex

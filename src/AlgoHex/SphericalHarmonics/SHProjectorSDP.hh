/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include <AlgoHex/Config/Version.hh>

#if ALGOHEX_WITH_MOSEK

#pragma once

#include "SHProjectorInterface.hh"

// forward declaration to avoid slow-to-compile mosek include:
namespace mosek::fusion { class Model; }

namespace AlgoHex {
/// Not threadsafe (contains Mosek model), recommended use
/// is with thread_local storage specifier.
class ALGOHEX_EXPORT SHProjectorSDP : public SHProjectorInterface
{
public:
    SHProjectorSDP();
    ~SHProjectorSDP();
    SHProjectionResult project(const SHCoeffs &sh_coeffs) override;
    SHCoeffs project_only_shc(const SHCoeffs &shc);
private:
    struct impl;
    std::unique_ptr<impl> impl_;
    mosek::fusion::Model *model();
};

} // namespace AlgoHex

#endif // ALGOHEX_WITH_MOSEK


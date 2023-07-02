/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

//== INCLUDES =================================================================

#include <chrono>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================




/** \class StopWatch StopWatch.hh
     
     Brief Description.
     
     A more elaborate description follows.
*/

class StopWatch
{
public:
  // Constructor
  StopWatch() {}

  // Destructor
  ~StopWatch() {}

  // Start time measurement
  void start()
  {
    starttime_ = std::chrono::steady_clock::now();
  }

  // Restart, return time elapsed until now.
  double restart()
  {
    double t = elapsed();
    start();
    return t;
  }

  // Stop time measurement, return time.
  double stop()
  {
    endtime_ = std::chrono::steady_clock::now();
    return elapsed();
  }

  // Get the total time in UnitT (watch has to be stopped).
  template<typename UnitT = std::chrono::milliseconds>
  double elapsed() const
  {
    auto duration = std::chrono::duration_cast<UnitT>(endtime_ - starttime_);
    return (double) duration.count();
  }

private:
  std::chrono::steady_clock::time_point starttime_, endtime_;
};

//=============================================================================
} // namespace AlgoHex
//=============================================================================

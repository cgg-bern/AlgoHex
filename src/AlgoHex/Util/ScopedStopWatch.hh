/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <AlgoHex/Config/Export.hh>
#include <string>
#include <vector>
#include <chrono>
#include <ostream>

namespace AlgoHex
{

class ALGOHEX_EXPORT Watch
{
public:
  using Clock = std::chrono::steady_clock;

  Watch(Watch *parent = nullptr)
          : parent_(parent),
            elapsed_(Clock::duration::zero()) {}

  void reset() { elapsed_ = Clock::duration::zero(); }



  void resume();




  void stop();




  size_t count() const;


  bool running() const { return running_; }

  const Clock::duration &elapsed() const { return elapsed_; }



  double elapsed_ms() const;


private:
  Watch *parent_;
  Clock::time_point start_;
  Clock::duration elapsed_;
  bool running_ = false;
  size_t count_ = 0; // how often have we resumed this stopwatch?
};

class ALGOHEX_EXPORT HierarchicalStopWatch
{
public:
  HierarchicalStopWatch(const std::string &name)
          : name_(name) {}



  HierarchicalStopWatch(const std::string &name,
                        HierarchicalStopWatch &parent);




  void reset();


  const std::string &name() const { return name_; }

  Watch &watch() { return watch_; }

  const Watch &watch() const { return watch_; }

  int level() const { return level_; }



  Watch::Clock::duration unaccounted() const;


  std::vector<HierarchicalStopWatch *>::const_iterator children_begin() const { return children_.cbegin(); }

  std::vector<HierarchicalStopWatch *>::const_iterator children_end() const { return children_.cend(); }

private:
  ALGOHEX_EXPORT friend std::ostream &operator<<(std::ostream &os, const HierarchicalStopWatch &sw);
//    friend std::ostream& print(std::ostream& os, HierarchicalStopWatch &sw);


//    friend std::ostream& print(std::ostream& os, HierarchicalStopWatch &sw);

  void add_child(HierarchicalStopWatch *sw);


  std::string name_;
  Watch watch_;
  int level_ = 0;
  std::vector<HierarchicalStopWatch *> children_;
};

ALGOHEX_EXPORT std::ostream &operator<<(std::ostream &os, const HierarchicalStopWatch &sw);


/// RAII helper to conveniently time the execution of parts of the code.
/// If the given watch is stopped at construction time, it will resume it and stop when it goes out of scope.
/// If it is already running at construction time, do nothing (not upon destruction either).
///
class ALGOHEX_EXPORT ScopedStopWatch
{
public:

public:
  ScopedStopWatch(Watch &watch);




  ScopedStopWatch(HierarchicalStopWatch &hsw);




  ~ScopedStopWatch();


private:


private:
  void maybe_start();


  Watch &watch_;
  bool should_stop_ = true;
};

} //namespace AlgoHex

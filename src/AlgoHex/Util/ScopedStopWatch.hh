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
#include <cassert>
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

  void resume()
  {
    assert(!parent_ || parent_->running());
    assert(!running_);
    running_ = true;
    start_ = Clock::now();
    ++count_;
  }

  void stop()
  {
    assert(running_);
    running_ = false;
    Clock::time_point now = Clock::now();
    Clock::duration d = (now - start_);
    elapsed_ += d;
  }

  size_t count() const
  {
    return count_;
  }

  bool running() const { return running_; }

  const Clock::duration &elapsed() const { return elapsed_; }

  double elapsed_ms() const
  {
    return std::chrono::duration_cast<std::chrono::milliseconds>(elapsed_).count();
  }

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
                        HierarchicalStopWatch &parent)
          : name_(name),
            watch_(parent.watch_),
            level_(parent.level_ + 1)
  {
    parent.add_child(this);
  }

  void reset()
  {
    watch_.reset();
    for (auto ch: children_)
    {
      ch->reset();
    }
  }

  const std::string &name() const { return name_; }

  Watch &watch() { return watch_; }

  const Watch &watch() const { return watch_; }

  int level() const { return level_; }

  Watch::Clock::duration unaccounted() const
  {
    Watch::Clock::duration dur = watch_.elapsed();
    for (auto child: children_)
    {
      dur -= child->watch().elapsed();
    }
    return dur;
  }

  std::vector<HierarchicalStopWatch *>::const_iterator children_begin() const { return children_.cbegin(); }

  std::vector<HierarchicalStopWatch *>::const_iterator children_end() const { return children_.cend(); }

private:
  ALGOHEX_EXPORT friend std::ostream &operator<<(std::ostream &os, const HierarchicalStopWatch &sw);
//    friend std::ostream& print(std::ostream& os, HierarchicalStopWatch &sw);

  void add_child(HierarchicalStopWatch *sw)
  {
    children_.push_back(sw);
  }

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
  ScopedStopWatch(Watch &watch)
          : watch_(watch)
  {
    maybe_start();
  }

  ScopedStopWatch(HierarchicalStopWatch &hsw)
          : watch_(hsw.watch())
  {
    maybe_start();
  }

  ~ScopedStopWatch()
  {
    if (should_stop_)
    {
      watch_.stop();
    }
  }

private:
  void maybe_start()
  {
    if (watch_.running())
    {
      should_stop_ = false;
    }
    else
    {
      watch_.resume();
    }

  }

  Watch &watch_;
  bool should_stop_ = true;
};

} //namespace AlgoHex



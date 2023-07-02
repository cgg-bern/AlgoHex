/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "ScopedStopWatch.hh"

#include <string_view>
#include <iomanip>
#include <memory>
#include <algorithm>

namespace AlgoHex
{


static std::string strrep(const std::string_view s, int n)
{
  if (n <= 0)
    return "";
  std::string indent;
  indent.reserve(s.size() * n);
  for (int i = 0; i < n; ++i)
  {
    indent.append(s);
  }
  return indent;
}

struct WatchResult
{
  WatchResult(const HierarchicalStopWatch &hsw, int _depth = 0)
          : depth(_depth), name(hsw.name()), duration(hsw.watch().elapsed()), unaccounted(hsw.unaccounted()),
            count(hsw.watch().count())
  {
    auto cbegin = hsw.children_begin();
    auto cend = hsw.children_end();
    children.reserve(std::distance(cbegin, cend));
    std::for_each(cbegin, cend,
                  [this](const HierarchicalStopWatch *ch) {
                    children.emplace_back(std::make_unique<WatchResult>(*ch, depth + 1));
                  });
  }

  int depth;
  std::string_view name;
  Watch::Clock::duration duration;
  Watch::Clock::duration unaccounted;
  int count;
  std::vector<std::unique_ptr<WatchResult>> children;
};


/// number of columns used for each indent level,
/// *not* the same as length of the indent string in bytes - unicode!
static const int indent_ncol = 2;

const std::string str_unacc{"(unaccounted)"};

size_t max_print_width(const WatchResult &wr)
{
  size_t my_indent = wr.depth * indent_ncol;
  size_t longest = my_indent + wr.name.size();
  if (!wr.children.empty())
  {
    longest = std::max(longest, my_indent + 2 + str_unacc.size());
  }
  for (auto &ch: wr.children)
  {
    size_t width = max_print_width(*ch);
    if (width > longest)
      longest = width;
  }
  return longest;
}


static std::ostream &print_duration(std::ostream &os, const Watch::Clock::duration &dur)
{
  os << std::chrono::duration_cast<std::chrono::milliseconds>(dur).count() << " ms";
  return os;
}

static std::ostream &print(std::ostream &os, const WatchResult &wr, size_t max_width)
{
  // cf. https://en.wikipedia.org/wiki/Box-drawing_character
  // TODO: properly align all timings
  std::ios oldState(nullptr);
  oldState.copyfmt(os);

  std::string indent = strrep("│ ", wr.depth - 1);
  os << indent;
  int textwidth = max_width;
  if (wr.depth > 0)
  {
    textwidth -= (wr.depth) * 2;
    os << "├ ";
  }

  os << std::setw(textwidth) << std::left
     << wr.name
     << std::setw(10) << std::right;
  print_duration(os, wr.duration);
  os << std::left << ", n = "
     << std::setw(5) << std::right
     << wr.count << "\n";
  os.copyfmt(oldState);

  for (auto &child: wr.children)
  {
    print(os, *child, max_width);
  }
  if (!wr.children.empty())
  {
    if (wr.depth > 0)
      os << indent << "│ ";
    os << "╰ ";
    os << std::setw(textwidth - 2) << std::left
       << str_unacc
       << std::setw(10) << std::right;
    print_duration(os, wr.unaccounted);
    os << "\n";
    os.copyfmt(oldState);
  }
  os << std::flush;
  return os;
}

ALGOHEX_EXPORT std::ostream &operator<<(std::ostream &os, const HierarchicalStopWatch &sw)
{
  WatchResult wr{sw};
  return print(os, wr, max_print_width(wr));
}

} // namespace AlgoHex

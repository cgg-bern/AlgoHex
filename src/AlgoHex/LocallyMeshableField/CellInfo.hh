/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include <ostream>
#include <AlgoHex/TypeDef.hh>

namespace AlgoHex
{
struct CELLINFO
{
public:
  CELLINFO(const CH _ch, const Vec3d &_param) : ch(_ch), param(_param) {}


  bool operator<(const CELLINFO &other) const
  {
    return ch < other.ch;
  }

  friend std::ostream &operator<<(std::ostream &out, const CELLINFO &ci)
  {
    return out << "cell: " << ci.ch << ", param: " << ci.param.transpose() << " distance: " << ci.length << std::endl;
  }

  CH ch;
  Vec3d param;
  //matching from the start cell to the current cell
  int trans = 0;

  double length = 0;
};


struct HALFFACEINFO
{
public:
  explicit HALFFACEINFO(const HFH _hfh, const Vec3d &_param) : hfh(_hfh), param(_param) {}


  bool operator<(const HALFFACEINFO &other) const
  {
    return hfh < other.hfh;
  }

  friend std::ostream &operator<<(std::ostream &out, const HALFFACEINFO &ci)
  {
    return out << "halfface: " << ci.hfh << ", param: " << ci.param.transpose() << " distance: " << ci.length
               << std::endl;
  }

  HFH hfh;
  Vec3d param;
  //matching from the start cell to the current cell
  int trans = 0;

  double length = 0;
};
}
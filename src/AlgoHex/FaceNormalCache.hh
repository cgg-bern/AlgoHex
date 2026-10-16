/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#include "TypeDef.hh"
#include <map>

namespace AlgoHex
{


class ALGOHEX_EXPORT FaceNormalCache
{
public:
  using Normal = TetMesh::PointT;

  FaceNormalCache(const TetMesh &mesh);


  FaceNormalCache(const FaceNormalCache &) = delete;

  Normal operator[](HFH hfh);


private:

  /// Get normal from cache if available, otherwise compute it (and update cache)
  const Normal &face_normal(FH fh);


  const TetMesh &mesh_;
  std::map<FH, Normal> face_normal_;

};

}

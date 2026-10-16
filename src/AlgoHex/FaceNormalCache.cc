/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#include "FaceNormalCache.hh"

namespace AlgoHex
{

  FaceNormalCache::FaceNormalCache(const TetMesh &mesh)
          : mesh_(mesh) {}

  FaceNormalCache::Normal FaceNormalCache::operator[](HFH hfh)
  {
    FH fh = mesh_.face_handle(hfh);
    if (hfh.idx() & 1)
    {
      return -face_normal(fh);
    }
    else
    {
      return face_normal(fh);
    }
  }

  const FaceNormalCache::Normal &FaceNormalCache::face_normal(FH fh)
  {
    auto it = face_normal_.find(fh);
    if (it != face_normal_.end())
    {
      return it->second;
    }
    else
    {
      Normal n = mesh_.normal(mesh_.halfface_handle(fh, 0));
      return face_normal_.emplace(fh, n).first->second;
    }

  }


}

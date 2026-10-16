/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include <AlgoHex/TypeDef.hh>
#include "SingularVertexOptProblem.hh"

namespace AlgoHex
{

template<typename Vec3>
void SingularVertexOptProblemQuaternion::project_on_tangent_plane(const Vec3 &_p0, const Vec3 &_px, const Vec3 &_n,
                                                                  Vec3 &_pp)
{
  _pp = _px - (_px - _p0).dot(_n) * _n;
}

template void ALGOHEX_EXPORT SingularVertexOptProblemQuaternion::project_on_tangent_plane<OVM::Geometry::Vec3d>(const OVM::Geometry::Vec3d &, const OVM::Geometry::Vec3d &, const OVM::Geometry::Vec3d &, OVM::Geometry::Vec3d &);
template void ALGOHEX_EXPORT SingularVertexOptProblemQuaternion::project_on_tangent_plane<Eigen::Vector3d>(const Eigen::Vector3d &, const Eigen::Vector3d &, const Eigen::Vector3d &, Eigen::Vector3d &);

} // namespace AlgoHex

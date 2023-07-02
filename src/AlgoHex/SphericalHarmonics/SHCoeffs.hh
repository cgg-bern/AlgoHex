/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

#include <AlgoHex/TypeDef.hh>
#include <AlgoHex/Config/Export.hh>
#include <ostream>

namespace AlgoHex
{

// It would be easier to derive from Vec9d, however this is currently impossible on MSVC:
// cf https://stackoverflow.com/questions/58630573/why-cant-i-inherit-from-eigenmatrix
//    https://gitlab.com/libeigen/eigen/-/issues/1768
//    https://eigen.tuxfamily.org/dox/TopicCustomizing_InheritingMatrix.html

class ALGOHEX_EXPORT SHCoeffs
{
public:
  using Scalar = double;
  using Vec9d = Eigen::Matrix<Scalar, 9, 1>;

  /// uninitalized!
  SHCoeffs() {}

  template<typename OtherDerived>
  SHCoeffs(const Eigen::MatrixBase<OtherDerived> &other)
          : vec_(other) {}

  SHCoeffs(const std::array<double, 9> &a)
          : vec_((Vec9d() << a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8]).finished())
  {
  }

  SHCoeffs &operator=(const std::array<double, 9> &a)
  {
    *this = SHCoeffs(a);
    return *this;
  }

  template<typename OtherDerived>
  SHCoeffs &operator=(const Eigen::MatrixBase<OtherDerived> &other)
  {
    vec_ = other;
    return *this;
  }

  template<typename OtherDerived>
  SHCoeffs &operator=(const Eigen::MatrixBase<OtherDerived> &&other)
  {
    vec_ = std::move(other);
    return *this;
  }

  template<typename OtherDerived>
  double dot(const Eigen::MatrixBase<OtherDerived> &other) const
  {
    return vec_.dot(other);
  }

  double dot(const SHCoeffs &other) const
  {
    return vec_.dot(other.vec());
  }

  void setZero() { vec_.setZero(); }

  void normalize() { vec_.normalize(); }

  SHCoeffs normalized() { return SHCoeffs(vec_.normalized()); }

  SHCoeffs &operator+=(const SHCoeffs &other)
  {
    vec_ += other.vec();
    return *this;
  }

  SHCoeffs operator+(const SHCoeffs &other) const
  {
    return {vec() + other.vec()};
  }

  SHCoeffs operator-(const SHCoeffs &other) const
  {
    return {vec() - other.vec()};
  }

  double norm() const { return vec_.norm(); }

  double squaredNorm() const { return vec_.squaredNorm(); }

  const Vec9d &vec() const { return vec_; }

  Vec9d &vec() { return vec_; }

  const double &operator[](size_t i) const { return vec()[i]; }

  double &operator[](size_t i) { return vec()[i]; }

  static SHCoeffs axis_aligned_cross()
  {
    return SHCoeffs({0, 0, 0, 0, std::sqrt(7. / 12.), 0, 0, 0, std::sqrt(5. / 12.)});
  }

  const Scalar *data() const
  {
    return vec_.data();
  }

  Scalar *data()
  {
    return vec_.data();
  }

  static SHCoeffs from_quaternion(const Quaternion &q);

  void rotate_x(double radians);

  void rotate_y(double radians);

  void rotate_z(double radians);

  void rotate_xyz(const Vec3d &radians);

  void rotate(const Quaternion &q);

private:
  Vec9d vec_;
};

ALGOHEX_EXPORT std::ostream &operator<<(std::ostream &s, const SHCoeffs &shc);

ALGOHEX_EXPORT inline SHCoeffs operator*(double d, const SHCoeffs &shc)
{
  return {d * shc.vec()};
}

// linearized axis rotations:

ALGOHEX_EXPORT SHCoeffs Ex(const SHCoeffs &a);

ALGOHEX_EXPORT SHCoeffs Ey(const SHCoeffs &a);

ALGOHEX_EXPORT SHCoeffs Ez(const SHCoeffs &a);

} // namespace AlgoHex

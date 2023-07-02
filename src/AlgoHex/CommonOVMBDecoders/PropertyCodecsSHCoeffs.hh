//
// Created by Denis Kalmykov on 16.03.22.
//

#ifndef ALGOHEX_PROPERTYCODECSSHCOEFFS_HH
#define ALGOHEX_PROPERTYCODECSSHCOEFFS_HH
#pragma once

#include "AlgoHex/TypeDef.hh"
#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/PropertyCodecsEigen.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>

#include "AlgoHex/SphericalHarmonics/SHCoeffs.hh"
#include <Eigen/Core>


namespace OpenVolumeMesh::IO::Codecs
{

struct SPHCoeffs
{
  using T = AlgoHex::SHCoeffs;
//        using T = EigenDenseFixedMatrix<double, 9, 1>::T;
  static EigenDenseFixedMatrix<double, 9, 1> mat_codec;

  static void encode(detail::Encoder &enc, const AlgoHex::SHCoeffs &val)
  {
    mat_codec.encode(enc, val.vec());
  }

  static void decode(detail::Decoder &reader, AlgoHex::SHCoeffs &val)
  {
    mat_codec.decode(reader, val.vec());
  }
};

struct Quaternion
{
  using T = AlgoHex::Quaternion;
//        using T = EigenDenseFixedMatrix<double, 9, 1>::T;
  static EigenDenseFixedMatrix<double, 4, 1> mat_codec;

  static void encode(detail::Encoder &enc, const AlgoHex::Quaternion &val)
  {
    mat_codec.encode(enc, val.coeffs());
  }

  static void decode(detail::Decoder &reader, AlgoHex::Quaternion &val)
  {
    mat_codec.decode(reader, val.coeffs());
  }
};

} // namespace OpenVolumeMesh::IO::Codecs

namespace OpenVolumeMesh::IO
{

inline void register_algohex_codecs(PropertyCodecs &_codecs)
{
  using namespace Codecs;
  _codecs.register_codec < SimplePropCodec < SPHCoeffs >> ("SPHCoeffs");
  _codecs.register_codec < SimplePropCodec < Quaternion >> ("Quaternion");
}

} // namespace OpenVolumeMesh::IO

#endif //ALGOHEX_PROPERTYCODECSSHCOEFFS_HH

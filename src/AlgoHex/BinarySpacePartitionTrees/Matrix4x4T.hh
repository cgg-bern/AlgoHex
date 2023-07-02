/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
\*===========================================================================*/





//=============================================================================
//
//  CLASS Matrix4x4T
//
//=============================================================================


#ifndef ACG_MATRIX4X4_HH
#define ACG_MATRIX4X4_HH


//== INCLUDES =================================================================

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <functional>


//== NAMESPACES  ==============================================================


namespace ACG {

    using OpenVolumeMesh::VectorT;

//== MACROS / HELPERS =========================================================

template<typename Scalar> inline bool checkEpsilon(Scalar x) {
  return fabs(x) < (1e-6);
}

template<> inline bool checkEpsilon(float x) {
  return fabs(x) < (1e-4);
}



//== CLASS DEFINITION =========================================================


/** Simple 4x4 Matrix. Inherited by GLMatrix.
*/

template <class Scalar>
class Matrix4x4T
{
public:

  /// constructor: uninitialized values
  Matrix4x4T() {}
 
  /// construct from other matrix type
  template <class OtherScalar>
  inline Matrix4x4T(const Matrix4x4T<OtherScalar>& _rhs) {
    operator=(_rhs);
  }
  
  /** setup matrix using an array of N*N scalar values.
      elements are ordered 'column first' (like OpenGL) */
  inline Matrix4x4T(const Scalar _array[16]) {
    std::copy(_array,_array+16, mat_);
  }

  /// destructor
  ~Matrix4x4T() {}


  /// assignment from other matrix type
  template<typename otherScalar>
  inline Matrix4x4T<Scalar>& operator=(const Matrix4x4T<otherScalar>& _rhs) {
    // produces warning C4244 on msvc (implicit cast double to float)
//    std::copy(_rhs.data(),_rhs.data()+16,mat_);
    for (int i = 0; i < 16; ++i)
      mat_[i] = Scalar(_rhs.data()[i]);

    return *this;
  }



  /// access operator (read and write)
  inline Scalar& operator()(unsigned int row, unsigned int col) {
    return  mat_[(row)+((col)<<2)];
  }

  /// access operator (read only)
  inline const Scalar& operator()(unsigned int row, unsigned int col) const {
    return mat_[(row)+((col)<<2)];
  }


  /// compare two matrices (up to some epsilon)
  inline bool operator== (const Matrix4x4T<Scalar>& _rhs) const {
    int i;
    const Scalar *a = mat_;
    const Scalar *b = _rhs.mat_;
    for(i=0;i< 16;i++,a++,b++)
      if(! checkEpsilon( *a - *b ))
	return false;
    return true;
  }

  /// compare two matrices
  inline bool operator!= (const Matrix4x4T<Scalar>& _rhs) const {
    return !( operator==(_rhs) );
  }


  /// self + _rhs
  inline Matrix4x4T operator+ (Matrix4x4T<Scalar> _rhs) const {
    std::transform(mat_,mat_+16, _rhs.mat_, _rhs.mat_, std::plus<Scalar>());
    return _rhs;
  }

  /// self - _rhs
  inline Matrix4x4T operator- (Matrix4x4T<Scalar> _rhs) const {
    std::transform(mat_,mat_+16, _rhs.mat_, _rhs.mat_, std::minus<Scalar>());
    return _rhs;
  }

  /// self * _rhs
  Matrix4x4T operator*(const Matrix4x4T<Scalar>& inst) const;

  /// self * scalar
  Matrix4x4T operator*(const Scalar& scalar);


  /// self += _rhs
  inline Matrix4x4T& operator+= ( const Matrix4x4T<Scalar>& _rhs) {
    std::transform(mat_,mat_+16,_rhs.mat_, mat_, std::plus<Scalar>());
    return *this;
  }

  /// self -= _rhs
  inline Matrix4x4T& operator-= ( const Matrix4x4T<Scalar>& _rhs) {
    std::transform(mat_,mat_+16, _rhs.mat_, mat_, std::minus<Scalar>());
    return *this;
  }

  /// self *= _rhs
  Matrix4x4T& operator*= (const Matrix4x4T<Scalar>& _rhs);

  /// multiply from left: self = _rhs * self
  Matrix4x4T& leftMult(const Matrix4x4T<Scalar>& _rhs);


  /// matrix by vector multiplication
  template <typename T>
  inline VectorT<T,4> operator*(const VectorT<T,4>& _v) const;

  /// transform point (x',y',z',1) = M * (x,y,z,1)
  template <typename T>
  inline VectorT<T,3> transform_point(const VectorT<T,3>& _v) const;

  /// transform vector (x',y',z',0) = A * (x,y,z,0)
  template <typename T>
  inline VectorT<T,3> transform_vector(const VectorT<T,3>& _v) const;

  /// sets all elements to zero
  inline void clear();

  /// setup an identity matrix
  inline void identity();


  /// check if the matrix is the identity ( up to an epsilon )
  inline bool is_identity() const {
    int i;
    const Scalar *a = mat_;
    Scalar b = 0.0;
    for(i=0;i< 16;i++,a++,b++) {
      if ( ( i == 0) || ( i == 5 ) || ( i == 10 ) || ( i == 15 ) )
         b = 1.0;
      else
         b = 0.0;
      if(! checkEpsilon( *a - b ))
        return false;
    }
    return true;
  }


  /// transpose matrix
  inline void transpose();

  
  /// matrix inversion (returns true on success)
  bool invert();

  Scalar determinant() const {
      return  mat_[12] * mat_[9] * mat_[6] * mat_[3] - mat_[8] * mat_[13] * mat_[6] * mat_[3] -
              mat_[12] * mat_[5] * mat_[10] * mat_[3] + mat_[4] * mat_[13] * mat_[10] * mat_[3] +
              mat_[8] * mat_[5] * mat_[14] * mat_[3] - mat_[4] * mat_[9] * mat_[14] * mat_[3] -
              mat_[12] * mat_[9] * mat_[2] * mat_[7] + mat_[8] * mat_[13] * mat_[2] * mat_[7] +
              mat_[12] * mat_[1] * mat_[10] * mat_[7] - mat_[0] * mat_[13] * mat_[10] * mat_[7] -
              mat_[8] * mat_[1] * mat_[14] * mat_[7] + mat_[0] * mat_[9] * mat_[14] * mat_[7] +
              mat_[12] * mat_[5] * mat_[2] * mat_[11] - mat_[4] * mat_[13] * mat_[2] * mat_[11] -
              mat_[12] * mat_[1] * mat_[6] * mat_[11] + mat_[0] * mat_[13] * mat_[6] * mat_[11] +
              mat_[4] * mat_[1] * mat_[14] * mat_[11] - mat_[0] * mat_[5] * mat_[14] * mat_[11] -
              mat_[8] * mat_[5] * mat_[2] * mat_[15] + mat_[4] * mat_[9] * mat_[2] * mat_[15] +
              mat_[8] * mat_[1] * mat_[6] * mat_[15] - mat_[0] * mat_[9] * mat_[6] * mat_[15] -
              mat_[4] * mat_[1] * mat_[10] * mat_[15] + mat_[0] * mat_[5] * mat_[10] * mat_[15];
  }


  /** access to data array. not very nice, but in case of 4x4 matrices
      this member can be used to pass matrices to OpenGL
      e.g. glLoadMatrixf(m.get_raw_data()); */
  inline const Scalar* get_raw_data() const { return mat_; }
  inline const Scalar* raw() const { return mat_; }
  inline const Scalar* data() const { return mat_; }

protected:

    Scalar mat_[16];
};


/// typedef
typedef Matrix4x4T<float>  Matrix4x4f;
/// typedef
typedef Matrix4x4T<double> Matrix4x4d;




//== IO to/from streams =======================================================


/// output matrix to ostream os
template<typename Scalar>
inline std::ostream& 
operator<<(std::ostream& os, const Matrix4x4T<Scalar>& m)
{
  for(int i=0; i<4; i++)
  {
    for(int j=0; j<4; j++)
      os << m(i,j) << " ";
    os << "\n";
  }
  return os;
}
 
 
/// read the space-separated components of a vector from a stream */
template<typename Scalar>
inline std::istream& 
operator>>(std::istream& is, Matrix4x4T<Scalar>& m) 
{
  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      is >> m(i,j);
  return is;
}

//=============================================================================
} // namespace ACG
//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(ACG_MATRIX4X4_C)
#define ACG_MATRIX4X4_TEMPLATES
#include "Matrix4x4T_impl.hh"
#endif
//=============================================================================
#endif // ACG_MATRIX4X4_HH defined
//=============================================================================


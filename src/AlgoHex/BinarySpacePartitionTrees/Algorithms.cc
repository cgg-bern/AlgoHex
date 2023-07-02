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
//  CLASS Geometry - IMPLEMENTATION
//
//=============================================================================

#define GEO_ALGORITHMS_C

//== INCLUDES =================================================================

#include "Algorithms.hh"


#ifdef max
#  undef max
#endif

#ifdef min
#  undef min
#endif


//----------------------------------------------------------------------------


namespace ACG {
namespace Geometry {


//== IMPLEMENTATION ========================================================== 


//== 3D STUFF ================================================================ 

  
template<typename Scalar>
bool
baryCoord( const VectorT<Scalar,3> & _p,
	   const VectorT<Scalar,3> & _u,
	   const VectorT<Scalar,3> & _v,
	   const VectorT<Scalar,3> & _w,
	   VectorT<Scalar,3> & _result )
{
  VectorT<Scalar,3>  vu = _v - _u,
                          wu = _w - _u,
                          pu = _p - _u;

  
  // find largest absolute coodinate of normal
  Scalar nx = vu[1]*wu[2] - vu[2]*wu[1],
         ny = vu[2]*wu[0] - vu[0]*wu[2],
         nz = vu[0]*wu[1] - vu[1]*wu[0],
         ax = fabs(nx),
         ay = fabs(ny),
         az = fabs(nz);


  unsigned char max_coord;

  if ( ax > ay ) {
    if ( ax > az ) {
      max_coord = 0;
    }
    else {
      max_coord = 2;
    }
  }
  else {
    if ( ay > az ) {
      max_coord = 1;
    }
    else {
      max_coord = 2;
    }
  }


  // solve 2D problem
  switch (max_coord) {

    case 0:
    {      
      if (1.0+ax == 1.0) return false;
      _result[1] = 1.0 + (pu[1]*wu[2]-pu[2]*wu[1])/nx - 1.0;
      _result[2] = 1.0 + (vu[1]*pu[2]-vu[2]*pu[1])/nx - 1.0;
      _result[0] = 1.0 - _result[1] - _result[2];
    }
    break;
      
    case 1:
    {
      if (1.0+ay == 1.0) return false;
      _result[1] = 1.0 + (pu[2]*wu[0]-pu[0]*wu[2])/ny - 1.0;
      _result[2] = 1.0 + (vu[2]*pu[0]-vu[0]*pu[2])/ny - 1.0;
      _result[0] = 1.0 - _result[1] - _result[2];
    }
    break;
      
    case 2:
    {
      if (1.0+az == 1.0) return false;
      _result[1] = 1.0 + (pu[0]*wu[1]-pu[1]*wu[0])/nz - 1.0;
      _result[2] = 1.0 + (vu[0]*pu[1]-vu[1]*pu[0])/nz - 1.0;
      _result[0] = 1.0 - _result[1] - _result[2];
    }
    break;    
  }

  return true;
}


//-----------------------------------------------------------------------------


template <class Vec>
typename Vec::value_type
triangleAreaSquared( const Vec& _v0,
		     const Vec& _v1,
		     const Vec& _v2 )
{
  return 0.25 * ((_v1-_v0)%(_v2-_v0)).sqrnorm();
}


//-----------------------------------------------------------------------------


template<class Vec>
typename Vec::value_type
distPointLineSquared( const Vec& _p,
		      const Vec& _v0,
		      const Vec& _v1,
		      Vec*       _min_v )
{
  Vec d1(_p-_v0), d2(_v1-_v0), min_v(_v0);
  typename Vec::value_type t = (d1|d2)/ d2.sqrnorm();

  if (t >  1.0)       d1 = _p - (min_v = _v1);
  else if (t > 0.0)   d1 = _p - (min_v = _v0 + d2*t);

  if (_min_v) *_min_v = min_v;
  return d1.sqrnorm();
}    

  //-----------------------------------------------------------------------------

template <class Vec>
typename Vec::value_type
distPointTriangleSquared( const Vec& _p,
                          const Vec& _v0,
                          const Vec& _v1,
                          const Vec& _v2,
                          Vec& _nearestPoint,
                          bool _stable)
{
  Vec v0v1 = _v1 - _v0;
  Vec v0v2 = _v2 - _v0;
  Vec n = v0v1 % v0v2; // not normalized !
  typename Vec::value_type d = n.sqrnorm();

  
  // Check if the triangle is degenerated
  if (d < FLT_MIN && d > -FLT_MIN) {
    if (_stable) {
      const double l0 = v0v1.sqrnorm();
      const double l1 = v0v2.sqrnorm();
      const double l2 = (_v2 - _v1).sqrnorm();
      if (l0 > l1 && l0 > l2) {
        return distPointLineSquared(_p, _v0, _v1, &_nearestPoint);
      }
      else if (l1 > l0 && l1 > l2) {
        return distPointLineSquared(_p, _v0, _v2, &_nearestPoint);
      }
      else {
        return distPointLineSquared(_p, _v1, _v2, &_nearestPoint);
      }
    } else {
      std::cerr << "distPointTriangleSquared: Degenerated triangle !\n";
      return -1.0;
    }
  }
  typename Vec::value_type invD = typename Vec::value_type(1.0) / d;

  
  // these are not needed for every point, should still perform
  // better with many points against one triangle
  Vec v1v2 = _v2 - _v1;
  typename Vec::value_type inv_v0v2_2 = typename Vec::value_type(1.0) / v0v2.sqrnorm();
  typename Vec::value_type inv_v0v1_2 = typename Vec::value_type(1.0) / v0v1.sqrnorm();
  typename Vec::value_type inv_v1v2_2 = typename Vec::value_type(1.0) / v1v2.sqrnorm();

  
  Vec v0p = _p - _v0;
  Vec t = v0p % n;
  typename Vec::value_type  s01, s02, s12;
  typename Vec::value_type a = (t | v0v2) * -invD;
  typename Vec::value_type b = (t | v0v1) * invD;

  
  if (a < 0)
  {
    // Calculate the distance to an edge or a corner vertex
    s02 = ( v0v2 | v0p ) * inv_v0v2_2;
    if (s02 < 0.0)
    {
      s01 = ( v0v1 | v0p ) * inv_v0v1_2;
      if (s01 <= 0.0) {
	v0p = _v0;
      } else if (s01 >= 1.0) {
	v0p = _v1;
      } else {
	v0p = _v0 + v0v1 * s01;
      }
    } else if (s02 > 1.0) {
      s12 = ( v1v2 | ( _p - _v1 )) * inv_v1v2_2;
      if (s12 >= 1.0) {
	v0p = _v2;
      } else if (s12 <= 0.0) {
	v0p = _v1;
      } else {
	v0p = _v1 + v1v2 * s12;
      }
    } else {
      v0p = _v0 + v0v2 * s02;
    }
  } else if (b < 0.0) {
    // Calculate the distance to an edge or a corner vertex
    s01 = ( v0v1 | v0p ) * inv_v0v1_2;
    if (s01 < 0.0)
    {
      s02 = ( v0v2 |  v0p ) * inv_v0v2_2;
      if (s02 <= 0.0) {
	v0p = _v0;
      } else if (s02 >= 1.0) {
	v0p = _v2;
      } else {
	v0p = _v0 + v0v2 * s02;
      }
    } else if (s01 > 1.0) {
      s12 = ( v1v2 | ( _p - _v1 )) * inv_v1v2_2;
      if (s12 >= 1.0) {
	v0p = _v2;
      } else if (s12 <= 0.0) {
	v0p = _v1;
      } else {
	v0p = _v1 + v1v2 * s12;
      }
    } else {
      v0p = _v0 + v0v1 * s01;
    }
  } else if (a+b > 1.0) {
    // Calculate the distance to an edge or a corner vertex
    s12 = ( v1v2 | ( _p - _v1 )) * inv_v1v2_2;
    if (s12 >= 1.0) {
      s02 = ( v0v2 | v0p ) * inv_v0v2_2;
      if (s02 <= 0.0) {
	v0p = _v0;
      } else if (s02 >= 1.0) {
	v0p = _v2;
      } else {
	v0p = _v0 + v0v2*s02;
      }
    } else if (s12 <= 0.0) {
      s01 = ( v0v1 |  v0p ) * inv_v0v1_2;
      if (s01 <= 0.0) {
	v0p = _v0;
      } else if (s01 >= 1.0) {
	v0p = _v1;
      } else {
	v0p = _v0 + v0v1 * s01;
      }
    } else {
      v0p = _v1 + v1v2 * s12;
    }
  } else {
    // Calculate the distance to an interior point of the triangle
    _nearestPoint = _p - n*((n|v0p) * invD);
    return (_nearestPoint - _p).sqrnorm();
  }
  
  _nearestPoint = v0p;

  return (_nearestPoint - _p).sqrnorm();
}

//-----------------------------------------------------------------------------


template <class Vec>
typename Vec::value_type
distPointTriangleSquared( const Vec& _p,
                          const Vec& _v0,
                          const Vec& _v1,
                          const Vec& _v2,
                          Vec& _nearestPoint )
{
  return distPointTriangleSquared(_p, _v0, _v1, _v2, _nearestPoint, false);
}


//-----------------------------------------------------------------------------


template <class Vec>
typename Vec::value_type
distPointTriangleSquaredStable( const Vec& _p,
                                const Vec& _v0,
                                const Vec& _v1,
                                const Vec& _v2,
                                Vec& _nearestPoint )
{
  return distPointTriangleSquared(_p, _v0, _v1, _v2, _nearestPoint, true);
}


//-----------------------------------------------------------------------------



//
// Modified code of Dave Eberly (www.magic-software.com)
//

template<typename Scalar>
Scalar
distLineLineSquared( const VectorT<Scalar,3>& _v00,
		     const VectorT<Scalar,3>& _v01,
		     const VectorT<Scalar,3>& _v10,
		     const VectorT<Scalar,3>& _v11,
		     VectorT<Scalar,3>*       _min_v0,
		     VectorT<Scalar,3>*       _min_v1,
		     bool                          _fastApprox )
{
  VectorT<Scalar,3> kDiff = _v00 - _v10;
  VectorT<Scalar,3> d0 = _v01-_v00;
  VectorT<Scalar,3> d1 = _v11-_v10;

  Scalar fA00 = d0.sqrnorm();
  Scalar fA01 = -(d0|d1);
  Scalar fA11 = d1.sqrnorm();
  Scalar fB0  = (kDiff|d0);
  Scalar fC   = kDiff.sqrnorm();
  Scalar fDet = fabs(fA00*fA11-fA01*fA01);
  Scalar fB1, fS, fT, fSqrDist, fTmp;

  
  if ( fDet >= FLT_MIN )
  {
    // line segments are not parallel
    fB1 = -(kDiff|d1);
    fS = fA01*fB1-fA11*fB0;
    fT = fA01*fB0-fA00*fB1;

    
    // conservative approximation only?
    if (_fastApprox)
    {
      if      ( fS < 0.0 )  fS = 0.0;
      else if ( fS > fDet ) fS = fDet;
      if      ( fT < 0.0 )  fT = 0.0;
      else if ( fT > fDet ) fT = fDet;
    }

    
    if ( fS >= 0.0 )
    {
      if ( fS <= fDet )
      {
	if ( fT >= 0.0 )
	{
	  if ( fT <= fDet )  // region 0 (interior)
	  {
	    // minimum at two interior points of 3D lines
	    Scalar fInvDet = 1.0/fDet;
	    fS *= fInvDet;
	    fT *= fInvDet;
	    fSqrDist = fS*(fA00*fS+fA01*fT+2.0*fB0) +
	               fT*(fA01*fS+fA11*fT+2.0*fB1)+fC;
	  }
	  else  // region 3 (side)
	  {
	    fT = 1.0;
	    fTmp = fA01+fB0;
	    if ( fTmp >= 0.0 )
	    {
	      fS = 0.0;
	      fSqrDist = fA11+2.0*fB1+fC;
	    }
	    else if ( -fTmp >= fA00 )
	    {
	      fS = 1.0;
	      fSqrDist = fA00+fA11+fC+2.0*(fB1+fTmp);
	    }
	    else
	    {
	      fS = -fTmp/fA00;
	      fSqrDist = fTmp*fS+fA11+2.0*fB1+fC;
	    }
	  }
	}
	else  // region 7 (side)
	{
	  fT = 0.0;
	  if ( fB0 >= 0.0 )
	  {
	    fS = 0.0;
	    fSqrDist = fC;
	  }
	  else if ( -fB0 >= fA00 )
	  {
	    fS = 1.0;
	    fSqrDist = fA00+2.0*fB0+fC;
	  }
	  else
	  {
	    fS = -fB0/fA00;
	    fSqrDist = fB0*fS+fC;
	  }
	}
      }
      else
      {
	if ( fT >= 0.0 )
	{
	  if ( fT <= fDet )  // region 1 (side)
	  {
	    fS = 1.0;
	    fTmp = fA01+fB1;
	    if ( fTmp >= 0.0 )
	    {
	      fT = 0.0;
	      fSqrDist = fA00+2.0*fB0+fC;
	    }
	    else if ( -fTmp >= fA11 )
	    {
	      fT = 1.0;
	      fSqrDist = fA00+fA11+fC+2.0*(fB0+fTmp);
	    }
	    else
	    {
	      fT = -fTmp/fA11;
	      fSqrDist = fTmp*fT+fA00+2.0*fB0+fC;
	    }
	  }
	  else  // region 2 (corner)
	  {
	    fTmp = fA01+fB0;
	    if ( -fTmp <= fA00 )
	    {
	      fT = 1.0;
	      if ( fTmp >= 0.0 )
	      {
		fS = 0.0;
		fSqrDist = fA11+2.0*fB1+fC;
	      }
	      else
	      {
		fS = -fTmp/fA00;
		fSqrDist = fTmp*fS+fA11+2.0*fB1+fC;
	      }
	    }
	    else
	    {
	      fS = 1.0;
	      fTmp = fA01+fB1;
	      if ( fTmp >= 0.0 )
	      {
		fT = 0.0;
		fSqrDist = fA00+2.0*fB0+fC;
	      }
	      else if ( -fTmp >= fA11 )
	      {
		fT = 1.0;
		fSqrDist = fA00+fA11+fC+2.0*(fB0+fTmp);
	      }
	      else
	      {
		fT = -fTmp/fA11;
		fSqrDist = fTmp*fT+fA00+2.0*fB0+fC;
	      }
	    }
	  }
	}
	else  // region 8 (corner)
	{
	  if ( -fB0 < fA00 )
	  {
	    fT = 0.0;
	    if ( fB0 >= 0.0 )
	    {
	      fS = 0.0;
	      fSqrDist = fC;
	    }
	    else
	    {
	      fS = -fB0/fA00;
	      fSqrDist = fB0*fS+fC;
	    }
	  }
	  else
	  {
	    fS = 1.0;
	    fTmp = fA01+fB1;
	    if ( fTmp >= 0.0 )
	    {
	      fT = 0.0;
	      fSqrDist = fA00+2.0*fB0+fC;
	    }
	    else if ( -fTmp >= fA11 )
	    {
	      fT = 1.0;
	      fSqrDist = fA00+fA11+fC+2.0*(fB0+fTmp);
	    }
	    else
	    {
	      fT = -fTmp/fA11;
	      fSqrDist = fTmp*fT+fA00+2.0*fB0+fC;
	    }
	  }
	}
      }
    }
    else 
    {
      if ( fT >= 0.0 )
      {
	if ( fT <= fDet )  // region 5 (side)
	{
	  fS = 0.0;
	  if ( fB1 >= 0.0 )
	  {
	    fT = 0.0;
	    fSqrDist = fC;
	  }
	  else if ( -fB1 >= fA11 )
	  {
	    fT = 1.0;
	    fSqrDist = fA11+2.0*fB1+fC;
	  }
	  else
	  {
	    fT = -fB1/fA11;
	    fSqrDist = fB1*fT+fC;
	  }
	}
	else  // region 4 (corner)
	{
	  fTmp = fA01+fB0;
	  if ( fTmp < 0.0 )
	  {
	    fT = 1.0;
	    if ( -fTmp >= fA00 )
	    {
	      fS = 1.0;
	      fSqrDist = fA00+fA11+fC+2.0*(fB1+fTmp);
	    }
	    else
	    {
	      fS = -fTmp/fA00;
	      fSqrDist = fTmp*fS+fA11+2.0*fB1+fC;
	    }
	  }
	  else
	  {
	    fS = 0.0;
	    if ( fB1 >= 0.0 )
	    {
	      fT = 0.0;
	      fSqrDist = fC;
	    }
	    else if ( -fB1 >= fA11 )
	    {
	      fT = 1.0;
	      fSqrDist = fA11+2.0*fB1+fC;
	    }
	    else
	    {
	      fT = -fB1/fA11;
	      fSqrDist = fB1*fT+fC;
	    }
	  }
	}
      }
      else   // region 6 (corner)
      {
	if ( fB0 < 0.0 )
	{
	  fT = 0.0;
	  if ( -fB0 >= fA00 )
	  {
	    fS = 1.0;
	    fSqrDist = fA00+2.0*fB0+fC;
	  }
	  else
	  {
	    fS = -fB0/fA00;
	    fSqrDist = fB0*fS+fC;
	  }
	}
	else
	{
	  fS = 0.0;
	  if ( fB1 >= 0.0 )
	  {
	    fT = 0.0;
	    fSqrDist = fC;
	  }
	  else if ( -fB1 >= fA11 )
	  {
	    fT = 1.0;
	    fSqrDist = fA11+2.0*fB1+fC;
	  }
	  else
	  {
	    fT = -fB1/fA11;
	    fSqrDist = fB1*fT+fC;
	  }
	}
      }
    }
  }
  else
  {
    // line segments are parallel
    if ( fA01 > 0.0 )
    {
      // direction vectors form an obtuse angle
      if ( fB0 >= 0.0 )
      {
	fS = 0.0;
	fT = 0.0;
	fSqrDist = fC;
      }
      else if ( -fB0 <= fA00 )
      {
	fS = -fB0/fA00;
	fT = 0.0;
	fSqrDist = fB0*fS+fC;
      }
      else
      {
	fB1 = -(kDiff|d1);
	fS = 1.0;
	fTmp = fA00+fB0;
	if ( -fTmp >= fA01 )
	{
	  fT = 1.0;
	  fSqrDist = fA00+fA11+fC+2.0*(fA01+fB0+fB1);
	}
	else
	{
	  fT = -fTmp/fA01;
	  fSqrDist = fA00+2.0*fB0+fC+fT*(fA11*fT+2.0*(fA01+fB1));
	}
      }
    }
    else
    {
      // direction vectors form an acute angle
      if ( -fB0 >= fA00 )
      {
	fS = 1.0;
	fT = 0.0;
	fSqrDist = fA00+2.0*fB0+fC;
      }
      else if ( fB0 <= 0.0 )
      {
	fS = -fB0/fA00;
	fT = 0.0;
	fSqrDist = fB0*fS+fC;
      }
      else
      {
	fB1 = -(kDiff|d1);
	fS = 0.0;
	if ( fB0 >= -fA01 )
	{
	  fT = 1.0;
	  fSqrDist = fA11+2.0*fB1+fC;
	}
	else
	{
	  fT = -fB0/fA01;
	  fSqrDist = fC+fT*(2.0*fB1+fA11*fT);
	}
      }
    }
  }


  if (_min_v0)  *_min_v0 = _v00 + fS*d0;
  if (_min_v1)  *_min_v1 = _v10 + fT*d1;

  return fabs(fSqrDist);
}

//-----------------------------------------------------------------------------

template < typename VectorT , typename ValueT >
inline
ValueT
distPointPlane(const VectorT& _porigin,
               const VectorT& _pnormal,
               const VectorT& _point)
{
  assert( fabs(_pnormal.norm()) > 0.0000000001) ;
  return( ( (_point - _porigin) | _pnormal ) );
}


//-----------------------------------------------------------------------------

template < typename VectorT >
VectorT projectToEdge(const VectorT& _start , const VectorT& _end , const VectorT& _point )
{
  VectorT direction = ( _end - _start ).normalize();
  assert( fabs(direction.norm()) > 0.0000000001) ;
  const VectorT projected_point = ( ( _point - _start ) | direction ) * direction + _start;
  
  if ( ( ( projected_point - _start ) | direction ) > ( ( _end - _start ) | direction ) )
    return _end;
  
  if ( ( ( projected_point - _start ) | direction ) < 0 )
    return _start;
  
  return projected_point;
}

//-----------------------------------------------------------------------------

template < typename VectorT >
inline
VectorT
projectToPlane(const VectorT& _porigin,
               const VectorT& _pnormal,
               const VectorT& _point)
{
  return (_point - _pnormal * distPointPlane< VectorT , double >( _porigin , _pnormal , _point ) );
}

//-----------------------------------------------------------------------------


template<typename Scalar>
bool
circumCenter( const VectorT<Scalar,3>&  _v0,
	      const VectorT<Scalar,3>&  _v1,
	      const VectorT<Scalar,3>&  _v2,
	      VectorT<Scalar,3>&        _result )
{
  VectorT<Scalar,3>   a(_v1-_v2),
                           b(_v2-_v0),
                           c(_v0-_v1),
                           G((_v0+_v1+_v2)/3.0);
  
  Scalar d0 = -(b|c),
         d1 = -(c|a),
         d2 = -(a|b),
         e0 = d1*d2,
         e1 = d2*d0,
         e2 = d0*d1,
         e  = e0+e1+e2;

  if (e <= std::numeric_limits<Scalar>::min())  return false;

  VectorT<Scalar,3>   H((e0*_v0 + e1*_v1 + e2*_v2)/e);

  _result = (Scalar)0.5 * ((Scalar)3.0*G - H);

  return true;
}



//-----------------------------------------------------------------------------


template<typename Scalar>
Scalar
circumRadiusSquared( const VectorT<Scalar,3>&  _v0,
		     const VectorT<Scalar,3>&  _v1,
		     const VectorT<Scalar,3>&  _v2 )
{
  VectorT<Scalar,3>  v0v1(_v1-_v0),
                          v0v2(_v2-_v0),
                          v1v2(_v2-_v1);

  Scalar denom = 4.0*((v0v1%v0v2).sqrnorm());
  if (denom < std::numeric_limits<Scalar>::min()  && 
      denom > -std::numeric_limits<Scalar>::min())  
    return std::numeric_limits<Scalar>::max();

  return ( v0v1.sqrnorm() *
	   v0v2.sqrnorm() *
	   v1v2.sqrnorm() /
	   denom );
}

//-----------------------------------------------------------------------------

template<typename Scalar>
bool
rotationOfTwoVectors( const VectorT<Scalar,3>&  _v0,
                      const VectorT<Scalar,3>&  _v1,
                      VectorT<Scalar,3>&  _axis,
                      Scalar& _angle,
                      bool _degree ) {

    // Copy axes
    VectorT<Scalar,3> v0 = _v0;
    VectorT<Scalar,3> v1 = _v1;

    // Normalize axes
    v0.normalize();
    v1.normalize();

    // Get rotation axis
    _axis = (v0 % v1).normalize();

    // Is nan?
    if ( std::isnan(_axis[0]) || std::isnan(_axis[1]) || std::isnan(_axis[2])  ) {
        return false;
    }

    // Get rotation angle (in radiant)
    _angle = acos(v0 | v1);

    // Is nan?
    if ( std::isnan(_angle) )
        _angle = 0.0;

    // Convert to degree
    if(_degree) {
        _angle *= 180.0 / M_PI;
    }

    return true;
}

//-----------------------------------------------------------------------------

template<class VectorT>
int isObtuse(const VectorT& _p0,
             const VectorT& _p1,
             const VectorT& _p2 )
{
  const double a0 = ((_p1-_p0)|(_p2-_p0));

  if ( a0<0.0) return 1;
  else
  {
    const double a1 = ((_p0-_p1)|(_p2-_p1));

    if (a1 < 0.0 ) return 2;
    else
    {
      const double a2 = ((_p0-_p2)|(_p1-_p2));

      if (a2 < 0.0 ) return 3;
      else return 0;
    }
  }
}

//-----------------------------------------------------------------------------


template<typename Scalar>
bool
minSphere(const VectorT<Scalar,3>&  _v0,
	  const VectorT<Scalar,3>&  _v1,
	  const VectorT<Scalar,3>&  _v2,
	  VectorT<Scalar,3>&        _center,
	  Scalar&                        _radius)
{
  VectorT<Scalar,3>   a(_v1-_v2),
                           b(_v2-_v0),
                           c(_v0-_v1);

  Scalar d0 = -(b|c),
         d1 = -(c|a),
         d2 = -(a|b);

  
  // obtuse angle ?
  if (d2 < std::numeric_limits<Scalar>::min())
  {
    _center = (_v0+_v1)*0.5;
    _radius = 0.5 * c.norm();
    return true;
  }
  if (d0 < std::numeric_limits<Scalar>::min())
  {
    _center = (_v1+_v2)*0.5;
    _radius = 0.5 * a.norm();
    return true;
  }
  if (d1 < std::numeric_limits<Scalar>::min())
  {
    _center = (_v2+_v0)*0.5;
    _radius = 0.5 * b.norm();
    return true;
  }
  

  // acute angle
  VectorT<Scalar,3>   G((_v0+_v1+_v2)/3.0);
  
  Scalar e0 = d1*d2,
         e1 = d2*d0,
         e2 = d0*d1,
         e  = e0+e1+e2;

  if ( e <= std::numeric_limits<Scalar>::min() )  return false;

  VectorT<Scalar,3>   H((e0*_v0 + e1*_v1 + e2*_v2)/e);

  _center = (Scalar)0.5 * ((Scalar)3.0*G - H);
  _radius = (_center-_v0).norm();
  
  return true;
}


//-----------------------------------------------------------------------------


template<typename Scalar>
Scalar
minRadiusSquared( const VectorT<Scalar,3>&  _v0,
		  const VectorT<Scalar,3>&  _v1,
		  const VectorT<Scalar,3>&  _v2 )
{
  VectorT<Scalar,3>  v0v1(_v1-_v0),
                          v0v2(_v2-_v0),
                          v1v2(_v2-_v1);

  Scalar denom = 4.0*((v0v1%v0v2).sqrnorm());
  if (denom < std::numeric_limits<Scalar>::min() && 
      denom > -std::numeric_limits<Scalar>::min())  
    return std::numeric_limits<Scalar>::max();

  Scalar  l0 = v0v1.sqrnorm(),
          l1 = v0v2.sqrnorm(),
          l2 = v1v2.sqrnorm(),
          l3 = l0*l1*l2/denom;
  
  return  std::max(std::max(l0, l1), std::max(l2, l3));
}


//-----------------------------------------------------------------------------


template<typename Scalar>
bool
circumCenter( const VectorT<Scalar,3>&  _v0,
	      const VectorT<Scalar,3>&  _v1,
	      const VectorT<Scalar,3>&  _v2,
	      const VectorT<Scalar,3>&  _v3,
	      VectorT<Scalar,3>&        _result )
{
  VectorT<Scalar,3>
    v0v1(_v1-_v0),
    v0v2(_v2-_v0),
    v0v3(_v3-_v0);

  Scalar  denom = ((v0v1[0]*v0v2[1]*v0v3[2] +
		    v0v1[1]*v0v2[2]*v0v3[0] +
		    v0v1[2]*v0v2[0]*v0v3[1]) -
		   (v0v1[0]*v0v2[2]*v0v3[1] +
		    v0v1[1]*v0v2[0]*v0v3[2] +
		    v0v1[2]*v0v2[1]*v0v3[0])) * 2.0;

  if (denom < std::numeric_limits<Scalar>::min() && 
      denom > -std::numeric_limits<Scalar>::min())  return false;

  
  _result = _v0 + (( v0v3.sqrnorm()*(v0v1%v0v2) +
		     v0v2.sqrnorm()*(v0v3%v0v1) +
		     v0v1.sqrnorm()*(v0v2%v0v3) ) / denom);
  
  return true;
}


//-----------------------------------------------------------------------------


template <typename Scalar>
VectorT<Scalar,3>
perpendicular( const VectorT<Scalar,3>&  v )
{
  if (fabs(v[0]) < fabs(v[1])) {
    if (fabs(v[0]) < fabs(v[2]))
      return VectorT<Scalar, 3>(Scalar(1.0) - v[0] * v[0], -v[0] * v[1], -v[0] * v[2]).normalize();
  } else {
    if (fabs(v[1]) < fabs(v[2]))
      return VectorT<Scalar, 3>(-v[1] * v[0], Scalar(1.0) - v[1] * v[1], -v[1] * v[2]).normalize();
  }

  return VectorT<Scalar, 3>(-v[2] * v[0], -v[2] * v[1], Scalar(1.0) - v[2] * v[2]).normalize();
}



//== 2D STUFF ================================================================ 



template<typename Scalar>
bool
baryCoord( const VectorT<Scalar,2> & _p,
	   const VectorT<Scalar,2> & _u,
	   const VectorT<Scalar,2> & _v,
	   const VectorT<Scalar,2> & _w,
	   VectorT<Scalar,3> & _result )
{
  Scalar  v0(_v[0]-_u[0]), v1(_v[1]-_u[1]),
          w0(_w[0]-_u[0]), w1(_w[1]-_u[1]),
          p0(_p[0]-_u[0]), p1(_p[1]-_u[1]),
          denom = v0*w1-v1*w0;

  if (1.0+fabs(denom) == 1.0) {
    std::cerr << "degen tri: ("
	      << _u << "), ("
	      << _v << "), ("
	      << _w << ")\n";
    return false;
  }
  
  _result[1] = 1.0 + ((p0*w1-p1*w0)/denom) - 1.0;
  _result[2] = 1.0 + ((v0*p1-v1*p0)/denom) - 1.0;
  _result[0] = 1.0 - _result[1] - _result[2];

  return true;
}


//-----------------------------------------------------------------------------


template<typename Scalar>
bool
lineIntersection( const VectorT<Scalar,2>&  _v0,
		  const VectorT<Scalar,2>&  _v1,
		  const VectorT<Scalar,2>&  _v2,
		  const VectorT<Scalar,2>&  _v3,
		  Scalar& _t1,
		  Scalar& _t2 )
{      
  _t1 = ((_v0[1]-_v2[1])*(_v3[0]-_v2[0])-(_v0[0]-_v2[0])*(_v3[1]-_v2[1])) /
        ((_v1[0]-_v0[0])*(_v3[1]-_v2[1])-(_v1[1]-_v0[1])*(_v3[0]-_v2[0]));
      
  _t2 = ((_v0[1]-_v2[1])*(_v1[0]-_v0[0])-(_v0[0]-_v2[0])*(_v1[1]-_v0[1])) /
    ((_v1[0]-_v0[0])*(_v3[1]-_v2[1])-(_v1[1]-_v0[1])*(_v3[0]-_v2[0]));
      
  return ((_t1>0.0) && (_t1<1.0) && (_t2>0.0) && (_t2<1.0));
}


//-----------------------------------------------------------------------------

template<typename Scalar>
bool
circumCenter( const VectorT<Scalar,2>&  _v0,
	      const VectorT<Scalar,2>&  _v1,
	      const VectorT<Scalar,2>&  _v2,
	      VectorT<Scalar,2>&        _result )
{
  Scalar x0(_v0[0]), y0(_v0[1]), xy0(x0*x0+y0*y0),
         x1(_v1[0]), y1(_v1[1]), xy1(x1*x1+y1*y1),
         x2(_v2[0]), y2(_v2[1]), xy2(x2*x2+y2*y2),
         a(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1),
         b(xy0*y1 - xy0*y2 - xy1*y0 + xy1*y2 + xy2*y0 - xy2*y1),
         c(xy0*x1 - xy0*x2 - xy1*x0 + xy1*x2 + xy2*x0 - xy2*x1);

  if (Scalar(1.0)+a == Scalar(1.0)) {
    std::cerr << "circumcircle: colinear points\n";
    return false;
  }

  _result[0] = 0.5*b/a;
  _result[1] = -0.5*c/a;
  
  return true;
}



//== N-DIM STUFF ============================================================== 


template <typename Scalar, int N>
Scalar
aspectRatio( const VectorT<Scalar, N>& _v0,
	     const VectorT<Scalar, N>& _v1,
	     const VectorT<Scalar, N>& _v2 )
{
  VectorT<Scalar,3> d0 = _v0 - _v1;
  VectorT<Scalar,3> d1 = _v1 - _v2;

  // finds the max squared edge length
  Scalar l2, maxl2 = d0.sqrnorm();
  if ((l2=d1.sqrnorm()) > maxl2)
    maxl2 = l2;
  // keep searching for the max squared edge length
  d1 = _v2 - _v0;
  if ((l2=d1.sqrnorm()) > maxl2)
    maxl2 = l2;

  // squared area of the parallelogram spanned by d0 and d1
  Scalar a2 = (d0 % d1).sqrnorm();

  // the area of the triangle would be
  // sqrt(a2)/2 or length * height / 2
  // aspect ratio = length / height
  //              = length * length / (2*area)
  //              = length * length / sqrt(a2)

  // returns the length of the longest edge
  //         divided by its corresponding height
  return sqrt( (maxl2 * maxl2) / a2 );
}


//-----------------------------------------------------------------------------


template <typename Scalar, int N>
Scalar
roundness( const VectorT<Scalar, N>& _v0,
	   const VectorT<Scalar, N>& _v1,
	   const VectorT<Scalar, N>& _v2 )
{
  const double FOUR_ROOT3 = 6.928203230275509;

  double area = 0.5*((_v1-_v0)%(_v2-_v0)).norm();

  return (Scalar) (FOUR_ROOT3 * area / ((_v0-_v1).sqrnorm() +
					(_v1-_v2).sqrnorm() +
					(_v2-_v0).sqrnorm() ));
}

template<typename Vec>
bool
triangleIntersection( const Vec&  _o,
		      const Vec&  _dir,
		      const Vec&  _v0,
		      const Vec&  _v1,
		      const Vec&  _v2,
		      typename Vec::value_type& _t,
		      typename Vec::value_type& _u,
		      typename Vec::value_type& _v )
{
    //This code effectively replicates the method described by Moeller et al. in "Fast, Minimum Storage Ray-Triangle Intersection".
    Vec edge1, edge2, tvec, pvec, qvec;
    typename Vec::value_type det, inv_det;

    //find vectors for two edges sharing v0
    edge1 = _v1-_v0;
    edge2 = _v2-_v0;

    //begin calculating determinant - also used to calculate u parameter
    pvec = _dir % edge2;

    //if determinant is near zero, the ray lies in plane of triangle
    det = edge1 | pvec;

    static const typename Vec::value_type EPSILON = std::numeric_limits<typename Vec::value_type>::epsilon() * 1e2;
    if (det > -EPSILON && det < EPSILON) {
        return false;
    }
    inv_det = typename Vec::value_type(1.0) / det;

    //calculate distance from vert0 to ray origin
    tvec = _o - _v0;

    //calculate U parameter and test bounds
    _u = (tvec | pvec) * inv_det;
    if (_u < 0.0 || _u > 1.0)
        return false;

    //prepare to test V parameter
    qvec = tvec % edge1;

    //calculate V parameter and test bounds
    _v = (_dir | qvec) * inv_det;
    if (_v < 0.0 || _u + _v > 1.0)
        return false;

    //Intersection found! Calculate t and exit...
    _t = (edge2 | qvec) * inv_det;
    return true;
}

template<typename Vec>
bool
axisAlignedBBIntersection( const Vec&  _o,
                  const Vec&  _dir,
                  const Vec& _bbmin,
                  const Vec& _bbmax,
                  typename Vec::value_type& tmin,
                  typename Vec::value_type& tmax )
{
    /*
    * Ray-box intersection using IEEE numerical properties to ensure that the
    * test is both robust and efficient, as described in:
    *
    *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
    *      "An Efficient and Robust Ray-Box Intersection Algorithm"
    *      Journal of graphics tools, 10(1):49-54, 2005
    *
    */
    typename Vec::value_type tymin, tymax, tzmin, tzmax;
    Vec inv_dir;

    inv_dir[0] = 1/_dir[0];
    inv_dir[1] = 1/_dir[1];
    inv_dir[2] = 1/_dir[2];

    if (inv_dir[0] >= 0) {
        tmin = (_bbmin[0] - _o[0]) * inv_dir[0];
        tmax = (_bbmax[0] - _o[0]) * inv_dir[0];
    }
    else {
        tmin = (_bbmax[0] - _o[0]) * inv_dir[0];
        tmax = (_bbmin[0] - _o[0]) * inv_dir[0];
    }

    if (inv_dir[1] >= 0) {
        tymin = (_bbmin[1] - _o[1]) * inv_dir[1];
        tymax = (_bbmax[1] - _o[1]) * inv_dir[1];
    }
    else {
        tymin = (_bbmax[1] - _o[1]) * inv_dir[1];
        tymax = (_bbmin[1] - _o[1]) * inv_dir[1];
    }

    if ( (tmin > tymax) || (tymin > tmax) )
        return false;
    if (tymin > tmin)
        tmin = tymin;
    if (tymax < tmax)
        tmax = tymax;

    if (inv_dir[2] >= 0) {
        tzmin = (_bbmin[2] - _o[2]) * inv_dir[2];
        tzmax = (_bbmax[2] - _o[2]) * inv_dir[2];
    }
    else {
        tzmin = (_bbmax[2] - _o[2]) * inv_dir[2];
        tzmax = (_bbmin[2] - _o[2]) * inv_dir[2];
    }

    if ( (tmin > tzmax) || (tzmin > tmax) )
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    
    return true;
}

//=============================================================================
} // namespace Geometry
} // namespace ACG
//=============================================================================

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once

//== INCLUDES =================================================================

#include <vector>
#include <algorithm>
#include <gmm/gmm.h>
#include <Eigen/Dense>

#include "AlgoHex/Config/Export.hh"


//== FORWARDDECLARATIONS ======================================================

#define ROUND_TQ(x) (x<0?int((x)-0.5):int((x)+0.5))

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== CLASS DEFINITION =========================================================




/** \class TransitionQuaternion TransitionQuaternion.hh

    Brief Description.

    A more elaborate description follows.
*/

class ALGOHEX_EXPORT TransitionQuaternion
{
public:
  using Quaternion = Eigen::Quaternion<double>;
  using Matrix = Eigen::Matrix3d;

  /// Default constructor
  TransitionQuaternion() { init_transition_quaternions(); }

  /// Destructor
  ~TransitionQuaternion() {}

  // return the index of the inverse transition
  int inverse_transition_idx(const int _i) const { if (_i >= 0) return inverse_transition_[_i]; else return -_i; }

  // return trans_i * trans_j
  int mult_transitions_idx(const int _i, const int _j) const;

  // return the index of the transition closest to _q
  int closest_transition_idx(const Quaternion &_q) const;

  // return the index of the transition t s.t. _i*t is_2d_rotation and t closest to _q
  int closest_transition_idx_generating_2d(const Quaternion &_q, const int _i) const;

  // return the index of the transition t s.t. _iv(j)*t is_2d_rotation and t closest to _q
  bool find_closest_transition_idx_generating_2d(const Quaternion &_q, const std::vector<int> _iv, int &_t) const;

  // get the quaternion of transition _i
  Quaternion transition(const int _i) const
  {
    if (_i >= 0) return transitions_[_i];
    else
      return inverse_transition(-_i);
  }

  // get the quaternion of the inverse transition
  Quaternion inverse_transition(const int _i) const
  {
    if (_i >= 0)
      return transition(inverse_transition_idx(_i));
    else return transition(-_i);
  }

  // notice that _q has to be a unit quaternion!!!
  Quaternion closest_transition(const Quaternion &_q) const { return transition(closest_transition_idx(_q)); }

  // return the index of the transition t s.t. _i*t is_2d_rotation and t closest to _q
  Quaternion closest_transition_generating_2d(const Quaternion &_q, const int _i) const
  {
    return transition(closest_transition_idx_generating_2d(_q, _i));
  }

  // get the matrix (rounded to int) of transition _i
  Matrix transition_matrix_int(const int _i) const;

  // return whether an transition is a simple 2d rotation
  bool is_2d_rotation(const int _i) const { if (_i >= 0) return is_2d_rotation_[_i]; else return is_2d_rotation(-_i); }

  // return the rotation axis
  int rotation_axis_2d(const int _i) const
  {
    if (_i >= 0) return rotation_axis_2d_[_i];
    else
      return rotation_axis_2d(-_i);
  }

  int rotation_angle_2d(const int _i) const
  {
    if (_i >= 0)
      return rotation_angle_deg_2d_[_i];
    else return rotation_angle_2d(-_i);
  }

  // check quaternions
  void check_quaternions() const;

  int axis_after_transition(const int _axis_start, const int _transition) const;

private:

  // initialize vector of transitions
  void init_transition_quaternions();

  void init_mult_table();


  // test method
  void test_round_to_2d_possibilities();

private:

  /// Copy constructor (not used)
  TransitionQuaternion(const TransitionQuaternion &_rhs);

  /// Assignment operator (not used)
  TransitionQuaternion &operator=(const TransitionQuaternion &_rhs);

private:

  // vector of transition quaternions
  std::vector<Quaternion> transitions_;

  // idx of inverse transition
  std::vector<int> inverse_transition_;

  // table for transition multiplication
  gmm::dense_matrix<int> mult_table_;

  // flag whether a rotation is only 2d
  std::vector<bool> is_2d_rotation_;
  // for 2d rotations we have a fixed rotation axis
  // convention [0,1,2,3,4,5] = [x,-x,y,-y,z,-z]
  std::vector<int> rotation_axis_2d_;
  // 2d rotation angle in DEGREE w.r.t. rotation_axis_2d_
  std::vector<int> rotation_angle_deg_2d_;

};


//=============================================================================
} // namespace AlgoHex
//=============================================================================
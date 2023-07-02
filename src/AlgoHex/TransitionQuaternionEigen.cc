/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

//== INCLUDES =================================================================

#include "TransitionQuaternionEigen.hh"

//== NAMESPACES ===============================================================

namespace AlgoHex
{

//== IMPLEMENTATION ==========================================================

void
TransitionQuaternion::
init_transition_quaternions()
{
  transitions_.clear();
  transitions_.reserve(24);

  inverse_transition_.clear();
  inverse_transition_.reserve(24);

  is_2d_rotation_.clear();
  is_2d_rotation_.reserve(24);

  rotation_axis_2d_.clear();
  rotation_axis_2d_.reserve(24);

  rotation_angle_deg_2d_.clear();
  rotation_angle_deg_2d_.reserve(24);

  double h = 0.5;
  double s = 0.5 * sqrt(2);

  // 1: identity
  transitions_.push_back(Quaternion(1, 0, 0, 0)); // 0

  inverse_transition_.push_back(0);

  is_2d_rotation_.push_back(true);

  rotation_axis_2d_.push_back(0);
  rotation_angle_deg_2d_.push_back(0);

  // 2: rotations around {x,y,z} by 180ï¿½
  transitions_.push_back(Quaternion(0, 1, 0, 0)); // 1
  transitions_.push_back(Quaternion(0, 0, 1, 0)); // 2
  transitions_.push_back(Quaternion(0, 0, 0, 1)); // 3

  inverse_transition_.push_back(1);
  inverse_transition_.push_back(2);
  inverse_transition_.push_back(3);

  is_2d_rotation_.push_back(true);
  is_2d_rotation_.push_back(true);
  is_2d_rotation_.push_back(true);

  rotation_axis_2d_.push_back(0);
  rotation_axis_2d_.push_back(2);
  rotation_axis_2d_.push_back(4);

  rotation_angle_deg_2d_.push_back(180);
  rotation_angle_deg_2d_.push_back(180);
  rotation_angle_deg_2d_.push_back(180);


  // 3: rotations around {x,-x,y,-y,z,-z} by 90ï¿½
  transitions_.push_back(Quaternion(s, s, 0, 0)); // 4
  transitions_.push_back(Quaternion(s, -s, 0, 0)); // 5
  transitions_.push_back(Quaternion(s, 0, s, 0)); // 6
  transitions_.push_back(Quaternion(s, 0, -s, 0)); // 7
  transitions_.push_back(Quaternion(s, 0, 0, s)); // 8
  transitions_.push_back(Quaternion(s, 0, 0, -s)); // 9

  inverse_transition_.push_back(5);
  inverse_transition_.push_back(4);
  inverse_transition_.push_back(7);
  inverse_transition_.push_back(6);
  inverse_transition_.push_back(9);
  inverse_transition_.push_back(8);

  is_2d_rotation_.push_back(true);
  is_2d_rotation_.push_back(true);
  is_2d_rotation_.push_back(true);
  is_2d_rotation_.push_back(true);
  is_2d_rotation_.push_back(true);
  is_2d_rotation_.push_back(true);

  rotation_axis_2d_.push_back(0);
  rotation_axis_2d_.push_back(1);
  rotation_axis_2d_.push_back(2);
  rotation_axis_2d_.push_back(3);
  rotation_axis_2d_.push_back(4);
  rotation_axis_2d_.push_back(5);

  rotation_angle_deg_2d_.push_back(90);
  rotation_angle_deg_2d_.push_back(90);
  rotation_angle_deg_2d_.push_back(90);
  rotation_angle_deg_2d_.push_back(90);
  rotation_angle_deg_2d_.push_back(90);
  rotation_angle_deg_2d_.push_back(90);


  // 4: rotations around {(x,y), (y,z), (z,x) (y,-x), (z,-y), (x,-z)} by 180ï¿½
  transitions_.push_back(Quaternion(0, s, s, 0)); // 10
  transitions_.push_back(Quaternion(0, 0, s, s)); // 11
  transitions_.push_back(Quaternion(0, s, 0, s)); // 12
  transitions_.push_back(Quaternion(0, -s, s, 0)); // 13
  transitions_.push_back(Quaternion(0, 0, -s, s)); // 14
  transitions_.push_back(Quaternion(0, s, 0, -s)); // 15

  inverse_transition_.push_back(10);
  inverse_transition_.push_back(11);
  inverse_transition_.push_back(12);
  inverse_transition_.push_back(13);
  inverse_transition_.push_back(14);
  inverse_transition_.push_back(15);

  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);

  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);

  rotation_angle_deg_2d_.push_back(180);
  rotation_angle_deg_2d_.push_back(180);
  rotation_angle_deg_2d_.push_back(180);
  rotation_angle_deg_2d_.push_back(180);
  rotation_angle_deg_2d_.push_back(180);
  rotation_angle_deg_2d_.push_back(180);


  // 5: rotations around
  // {(x,y,z), (x,y,-z), (x,-y, z), (x,-y,-z), (-x,y,z), (-x,y,-z), (-x,-y,z), (-x,-y,-z)}
  // by 120ï¿½
  transitions_.push_back(Quaternion(h, h, h, h)); // 16
  transitions_.push_back(Quaternion(h, h, h, -h)); // 17
  transitions_.push_back(Quaternion(h, h, -h, h)); // 18
  transitions_.push_back(Quaternion(h, h, -h, -h)); // 19
  transitions_.push_back(Quaternion(h, -h, h, h)); // 20
  transitions_.push_back(Quaternion(h, -h, h, -h)); // 21
  transitions_.push_back(Quaternion(h, -h, -h, h)); // 22
  transitions_.push_back(Quaternion(h, -h, -h, -h)); // 23

  inverse_transition_.push_back(23);
  inverse_transition_.push_back(22);
  inverse_transition_.push_back(21);
  inverse_transition_.push_back(20);
  inverse_transition_.push_back(19);
  inverse_transition_.push_back(18);
  inverse_transition_.push_back(17);
  inverse_transition_.push_back(16);

  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);
  is_2d_rotation_.push_back(false);

  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);
  rotation_axis_2d_.push_back(-1);

  rotation_angle_deg_2d_.push_back(120);
  rotation_angle_deg_2d_.push_back(120);
  rotation_angle_deg_2d_.push_back(120);
  rotation_angle_deg_2d_.push_back(120);
  rotation_angle_deg_2d_.push_back(120);
  rotation_angle_deg_2d_.push_back(120);
  rotation_angle_deg_2d_.push_back(120);
  rotation_angle_deg_2d_.push_back(120);

  // init multiplication table
  init_mult_table();

  //  test_round_to_2d_possibilities();
}

//-----------------------------------------------------------------------------

void
TransitionQuaternion::
init_mult_table()
{
  gmm::resize(mult_table_, 24, 24);

  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j)
    {
      Quaternion qi = transition(i);
      Quaternion qj = transition(j);
      Quaternion qij = qi * qj;

      mult_table_(i, j) = closest_transition_idx(qij);
    }
}

//-----------------------------------------------------------------------------

void
TransitionQuaternion::
check_quaternions() const
{
  std::vector<double> angles;
  for (int i = 0; i < (int) transitions_.size(); ++i)
    for (int j = 0; j < (int) transitions_.size(); ++j)
      if (i != j)
      {
        angles.push_back(transition(i).dot(transition(j)));
      }

  std::sort(angles.begin(), angles.end());
  //std::cerr << "pairangles: " << angles << std::endl;

  for (int i = 0; i < (int) transitions_.size(); ++i)
    std::cerr << closest_transition_idx(transition(i)) << " should be " << i << std::endl;

  for (int i = 0; i < (int) transitions_.size(); ++i)
    if (std::abs(std::abs((transition(i) * inverse_transition(i)).dot(transition(0))) - 1.0) > 1e-8)
    {
      auto qi = transition(i) * inverse_transition(i);
      std::cerr << "ERROR: trans*inv_trans != I \n";
      std::cerr << i << ": " << qi.w() << " " << qi.x() << " " << qi.y() << " " << qi.z() << std::endl;
    }
}

//-----------------------------------------------------------------------------

int
TransitionQuaternion::
closest_transition_idx(const Quaternion &_q) const
{
  double max_val = -2.0;
  int max_idx = 0;

  for (int i = 0; i < (int) transitions_.size(); ++i)
  {
    double cur_val = std::abs(_q.dot(transitions_[i]));
    if (cur_val > max_val)
    {
      max_idx = i;
      max_val = cur_val;
    }
  }

  return max_idx;
}

//-----------------------------------------------------------------------------

int
TransitionQuaternion::
closest_transition_idx_generating_2d(const Quaternion &_q, const int _i) const
{
  double max_val = -2.0;
  int max_idx = -1;

  for (int i = 0; i < (int) transitions_.size(); ++i)
  {
    // is valid?
    if (is_2d_rotation(mult_transitions_idx(_i, i)))
    {
      double cur_val = std::abs(_q.dot(transitions_[i]));
      if (cur_val > max_val)
      {
        max_idx = i;
        max_val = cur_val;
      }
    }
  }

  // debug output
  if (max_idx == -1)
    std::cerr << "ERROR: closest_transition_idx_generating_2d did not find valid transition!! "
              << _i << std::endl;

  return max_idx;
}

//-----------------------------------------------------------------------------

bool
TransitionQuaternion::
find_closest_transition_idx_generating_2d(const Quaternion &_q, const std::vector<int> _iv, int &_t) const
{
  double max_val = -2.0;
  int max_idx = -1;

  for (int i = 0; i < (int) transitions_.size(); ++i)
  {
    bool is_valid = true;

    // is valid?
    for (unsigned int j = 0; j < _iv.size(); ++j)
      if (!is_2d_rotation(mult_transitions_idx(_iv[j], i)))
        is_valid = false;

    // is valid?
    if (is_valid)
    {
      double cur_val = std::abs(_q.dot(transitions_[i]));
      if (cur_val > max_val)
      {
        max_idx = i;
        max_val = cur_val;
      }
    }
  }

  // found one?
  if (max_idx != -1)
  {
    _t = max_idx;
    return true;
  }
  else
    return false;
}

//-----------------------------------------------------------------------------

TransitionQuaternion::Matrix
TransitionQuaternion::
transition_matrix_int(int _i) const
{
  Matrix trans = transition(_i).toRotationMatrix();

  for (unsigned int i = 0; i < 3; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      trans(i, j) = ROUND_TQ(trans(i, j));

  return trans;
}

//-----------------------------------------------------------------------------

int
TransitionQuaternion::
mult_transitions_idx(const int _i, const int _j) const
{
  if (_i >= 0 && _j >= 0)
    return mult_table_(_i, _j);
  else if (_i < 0)
    return mult_table_(inverse_transition_idx(-_i), _j);
  else
    return mult_table_(_i, inverse_transition_idx(-_j));
}

//-----------------------------------------------------------------------------

void
TransitionQuaternion::
test_round_to_2d_possibilities()
{
  std::cerr << "test rounding...\n";
  // first check tuples
  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j)
    {
      int c;
      std::vector<int> iv;
      iv.push_back(i);
      iv.push_back(j);
      if (!find_closest_transition_idx_generating_2d(Quaternion(1, 0, 0, 0), iv, c))
        std::cerr << "no possible rounding for tuple " << i << " " << j << std::endl;
    }

  // check triples
  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j)
      for (int k = 0; k < 24; ++k)
      {
        int c;
        std::vector<int> iv;
        iv.push_back(i);
        iv.push_back(j);
        iv.push_back(k);
        if (!find_closest_transition_idx_generating_2d(Quaternion(1, 0, 0, 0), iv, c))
        {
          std::cerr << "no possible rounding for triple "
                    << i << " " << j << " " << k << std::endl;

          for (unsigned int l = 0; l < 24; ++l)
          {
            std::cerr << "transition " << l << " leads to singular edges "
                      << mult_transitions_idx(iv[0], l) << " "
                      << mult_transitions_idx(iv[1], l) << " "
                      << mult_transitions_idx(iv[2], l) << std::endl;
          }
        }
      }

  std::cerr << "test rounding done!\n";
}

//-----------------------------------------------------------------------------

int TransitionQuaternion::axis_after_transition(const int _axis_start, const int _transition) const
{
  if (_transition == -1 || _axis_start > 6 || _axis_start < 0)
  {
    std::cerr << "ERROR: input transition or axis is invalid! Axis start: " << _axis_start << " transition: "
              << _transition << std::endl;
    return -1;
  }

  Matrix m1 = transition_matrix_int(_transition);

  int axes_mapped[3] = {0, 0, 0};
  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
    {
      if (m1(i, j) == 1)
      {
        axes_mapped[j] = 2 * i;
        break;
      }
      if (m1(i, j) == -1)
      {
        axes_mapped[j] = 2 * i + 1;
        break;
      }
    }


  int pos = _axis_start / 2;
  int res = _axis_start % 2;
  int axis = 0;

  if (res == 0)
    axis = axes_mapped[pos];
  else if (res == 1)
    axis = (axes_mapped[pos] % 2 == 0) ? axes_mapped[pos] + 1 : axes_mapped[pos] - 1;

  return axis;
}


//=============================================================================
} // namespace AlgoHex
//=============================================================================


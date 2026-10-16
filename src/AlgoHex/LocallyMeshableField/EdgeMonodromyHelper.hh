/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

#ifdef ALGOHEX_VERBOSE
#define ALGOHEX_DEBUG_ONLY(x) x
#else
#define ALGOHEX_DEBUG_ONLY(x)
#endif

#include <AlgoHex/Geometry.hh>
#include <AlgoHex/TypeDef.hh>
#include <AlgoHex/AxisAlignment.hh>
#include <AlgoHex/TransitionQuaternionEigen.hh>

namespace AlgoHex
{
class EdgeMonodromyHelper
{
public:
  static int halfedge_transition_index(const TetMesh &_mesh, const TransitionQuaternion &_tq, const HFP<int> &_trans_prop,
                                       const HEH _heh, const HFH _hfh_start);


  static int halfedge_transition_index(const TetMesh &_mesh, const TransitionQuaternion &_tq, const HFP<int> &_trans_prop,
                                       const HEH _heh, const HFH _hfh_start, const HFH _hfh_end);


  static int halfedge_rotational_axis_idx(const TetMesh &_mesh, const TransitionQuaternion &_tq,
                                          const CP <Quaternion> &_cell_quaternions,
                                          const HFP<int> &_trans_prop, const EP<int> &_valence, const HEH _heh,
                                          const HFH _hfh_start);


  static int halfedge_transition_index_in_cell(const TetMesh &_mesh, const TransitionQuaternion &_tq,
                                               const CP <Quaternion> &_cell_quaternions,
                                               const HFP<int> &_trans_prop, const EP<int> &_valence, const HEH _heh,
                                               const HFH _hfh_start);



  static int get_dominant_axis_in_cell(const TetMesh &_mesh, const TransitionQuaternion &_tq,
                                       const CP <Quaternion> &_cell_quaternions,
                                       const HFP<int> &_trans_prop, const HEH _heh);


  static int
  halfedge_axis_idx(const TetMesh &_mesh, const CP <Quaternion> &_cell_quaternions, const HEH _heh, const CH _ch,
                    const int _axis);


  static int get_halfedge_aligned_axis_in_cell(const TetMesh &_mesh, const TransitionQuaternion &_tq,
                                               const CP <Quaternion> &_cell_quaternions,
                                               const HFP<int> &_trans_prop, const EP<int> &_valence,
                                               const EP<int> &_feature_edge, const HEH _heh, const HFH _hfh);


  static std::map<CH, int> get_halfedge_aligned_axis_in_cells(const TetMesh &_mesh, const TransitionQuaternion &_tq,
                                                              const CP <Quaternion> &_cell_quaternions,
                                                              const HFP<int> &_trans_prop, const EP<int> &_valence,
                                                              const EP<int> &_feature_edge, const HEH &_heh);


  //TODO: support high valence by checking angle defect
  static std::map<CH, int>
  get_halfedge_aligned_axis_in_cells_wrt_valence(const TetMesh &_mesh, const TransitionQuaternion &_tq,
                                                 const CP <Quaternion> &_cell_quaternions,
                                                 const HFP<int> &_trans_prop, const EP<int> &_valence,
                                                 const EP<int> &_feature_edge, const HEH &_heh);

};
}

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once
#include <AlgoHex/TypeDef.hh>

#include <AlgoHex/TransitionQuaternionEigen.hh>


namespace AlgoHex
{
class QuaternionSmoothing
{
public:
  using VecAxis = std::pair<Vec3d, int>;

  static Quaternion
  ALGOHEX_EXPORT construct_quaternion(const Quaternion &_qt, const std::vector<VecAxis> &_directions);

  static Quaternion project_quaternion(const Quaternion &_qt, const int _axis, const Vec3d &_direction);

  static Quaternion
  ALGOHEX_EXPORT project_quaternion(const Quaternion &_qt, const Vec3d &_direction);

  //align to feature faces, feature edges (and singular edges if added)
  template<class CT>
  static void
  smooth_field_quaternion(const TetMesh &mesh, const HFP<int> &trans_prop, const CT &cells,
                          const HFP<bool> &_is_hf_bdy,
                          const CP<std::vector<Vec3d>> &_ch_alignments, const AlgoHex::TransitionQuaternion &tq,
                          CP<Quaternion> &cell_quaternions);

  //align to feature faces, feature edges (and singular edges if added)
  template<class CT>
  static void
  smooth_field_quaternion(const TetMesh &mesh, const HFP<int> &trans_prop, const CT &cells,
                          const HFP<bool> &_is_hf_bdy,
                          const CP<std::vector<VecAxis>> &_ch_alignments, const AlgoHex::TransitionQuaternion &tq,
                          CP<Quaternion> &cell_quaternions);

  template<class CT>
  static void
  smooth_field_quaternion(const TetMesh &mesh, const HFP<int> &trans_prop, const CT &cells,
                          const std::set<FH> &_bdy_fhs,
                          std::map<CH, std::vector<VecAxis>> &_ch_alignments,
                          const AlgoHex::TransitionQuaternion &tq, CP<Quaternion> &cell_quaternions);

  //
  static void
  ALGOHEX_EXPORT optimize_quaternions(const TetMesh &mesh, const HFP<int> &trans_prop, const AlgoHex::TransitionQuaternion &tq,
                       CP<Quaternion> &cell_quaternions, const EP<int> &valence, const FP<int> &feature_face,
                       const EP<int> &feature_edge, const std::set<CH> &_cells,
                       std::map<CH, int> &_bdy_aligned_axis, const bool _align_sg = false);

  static void
  ALGOHEX_EXPORT optimize_quaternions_in_cells(const TetMesh &mesh, const HFP<int> &trans_prop,
                                            const AlgoHex::TransitionQuaternion &tq,
                                            CP<Quaternion> &cell_quaternions,
                                            const EP<int> &valence, const FP<int> &feature_face,
                                            const EP<int> &feature_edge,
                                            const std::set<CH> &_cells, const bool _align_sg = false);

  static void
  ALGOHEX_EXPORT optimize_quaternions_wrt_prescribed_valence(const TetMesh &mesh, const HFP<int> &trans_prop,
                                                          const AlgoHex::TransitionQuaternion &tq,
                                                          CP<Quaternion> &cell_quaternions,
                                                          const EP<int> &valence, const EP<int> &feature_edge,
                                                          const FP<int> &feature_face, const std::set<CH> &_cells,
                                                          const bool _align_sg = false);

  static void
  ALGOHEX_EXPORT optimize_quaternions_wrt_prescribed_valence(const TetMesh &mesh, const HFP<int> &trans_prop,
                                                          const AlgoHex::TransitionQuaternion &tq,
                                                          CP<Quaternion> &cell_quaternions,
                                                          const EP<int> &valence, const EP<int> &feature_edge,
                                                          const FP<int> &feature_face, const VH _vh_new,
                                                          const bool _align_sg = false);

  //
  static void
  ALGOHEX_EXPORT optimize_quaternions_at_edge(const TetMesh &mesh, const HFP<int> &trans_prop,
                                           const AlgoHex::TransitionQuaternion &tq,
                                           CP<Quaternion> &cell_quaternions, const EP<int> &valence,
                                           const FP<int> &feature_face, const EP<int> &feature_edge, const EH eh);

  static void
  ALGOHEX_EXPORT optimize_quaternions_at_vertex(const TetMesh &mesh, const HFP<int> &trans_prop,
                                             const AlgoHex::TransitionQuaternion &tq,
                                             CP<Quaternion> &cell_quaternions,
                                             const EP<int> &valence, const FP<int> &feature_face,
                                             const EP<int> &feature_edge,
                                             const VH _vh_new, const bool _align_sg = false);

  static void
  ALGOHEX_EXPORT optimize_quaternions_at_edges(const TetMesh &mesh, const HFP<int> &trans_prop,
                                            const AlgoHex::TransitionQuaternion &tq,
                                            CP<Quaternion> &cell_quaternions,
                                            const EP<int> &valence, const FP<int> &feature_face,
                                            const EP<int> &feature_edge, const std::set<EH> &_edges,
                                            const bool _align_sg = false);

  static void
  ALGOHEX_EXPORT optimize_quaternions_with_alignment_constraints(const TetMesh &mesh, const HFP<int> &trans_prop,
                                                              const AlgoHex::TransitionQuaternion &tq,
                                                              CP<Quaternion> &cell_quaternions,
                                                              const std::set<CH> &_cells,
                                                              std::map<CH, std::vector<VecAxis> > &_alignments,
                                                              const int _iter = 20);

  static void
  ALGOHEX_EXPORT optimize_quaternions_with_alignment_constraints(const TetMesh &mesh, const HFP<int> &trans_prop,
                                                              const AlgoHex::TransitionQuaternion &tq,
                                                              CP<Quaternion> &cell_quaternions,
                                                              const VH _vh_new,
                                                              std::map<CH, std::vector<VecAxis> > &_alignments,
                                                              const int _iter = 20);

};

// The container type used at call sites that pass `mesh.cells()` (a pair of cell iterators)
using CellRange = std::decay_t<decltype(std::declval<const TetMesh &>().cells())>;

extern template void ALGOHEX_EXPORT QuaternionSmoothing::smooth_field_quaternion<CellRange>(const TetMesh &, const HFP<int> &, const CellRange &, const HFP<bool> &, const CP<std::vector<VecAxis>> &, const AlgoHex::TransitionQuaternion &, CP<Quaternion> &);
extern template void ALGOHEX_EXPORT QuaternionSmoothing::smooth_field_quaternion<std::set<CH>>(const TetMesh &, const HFP<int> &, const std::set<CH> &, const HFP<bool> &, const CP<std::vector<VecAxis>> &, const AlgoHex::TransitionQuaternion &, CP<Quaternion> &);
extern template void ALGOHEX_EXPORT QuaternionSmoothing::smooth_field_quaternion<std::array<CH, 1>>(const TetMesh &, const HFP<int> &, const std::array<CH, 1> &, const HFP<bool> &, const CP<std::vector<VecAxis>> &, const AlgoHex::TransitionQuaternion &, CP<Quaternion> &);

extern template void ALGOHEX_EXPORT QuaternionSmoothing::smooth_field_quaternion<CellRange>(const TetMesh &, const HFP<int> &, const CellRange &, const std::set<FH> &, std::map<CH, std::vector<VecAxis>> &, const AlgoHex::TransitionQuaternion &, CP<Quaternion> &);
extern template void ALGOHEX_EXPORT QuaternionSmoothing::smooth_field_quaternion<std::set<CH>>(const TetMesh &, const HFP<int> &, const std::set<CH> &, const std::set<FH> &, std::map<CH, std::vector<VecAxis>> &, const AlgoHex::TransitionQuaternion &, CP<Quaternion> &);
extern template void ALGOHEX_EXPORT QuaternionSmoothing::smooth_field_quaternion<std::vector<CH>>(const TetMesh &, const HFP<int> &, const std::vector<CH> &, const std::set<FH> &, std::map<CH, std::vector<VecAxis>> &, const AlgoHex::TransitionQuaternion &, CP<Quaternion> &);


}

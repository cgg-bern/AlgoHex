/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "SplitPaperT.hh"
#include "../DetachAtInvalidSingularNodeT.hh"
#include "../LocalMeshabilityRepairT.hh"
#include "../FrameFieldOptimizationT.hh"

namespace AlgoHex
{
template<class MeshT>
void SplitApproachT<MeshT>::split_approach()
{
  remove_zigzag2();
  split_complex_singular_edges();

  SplitHelperT<MeshT> pp(mesh_);
  pp.split_for_single_alignment();

  FrameFieldOptimizationT<MeshT> ffo(mesh_);
  ffo.align_quaternions_to_feature(100);
}

template<class MeshT>
bool SplitApproachT<MeshT>::check_zigzag()
{
  bool zz_free = true;
  std::cerr << "Zigzag face left: ";
  for (const auto fh: mesh_.faces())
  {
    if (feature_fprop_[fh] > 0)
      continue;
    int n_sge = n_singular_edges_in_face(fh);
    if (n_sge > 1)
    {
      HFH hf0 = mesh_.halfface_handle(fh, 0);
      auto hes = mesh_.halfface(hf0).halfedges();

      std::vector<int> v_trans;
      for (auto hehi: hes)
      {
        auto ehi = mesh_.edge_handle(hehi);
        if (valence_[ehi] != 0)
        {

          int trans = EdgeMonodromyHelperT::halfedge_transition_index_in_cell(mesh_, tq_, cell_quaternions_,
                                                                              trans_prop_, valence_, hehi, hf0);
          v_trans.push_back(trans);
        }
      }

      if (n_sge == 2)
      {
        if (v_trans[0] == v_trans[1])
        {
          std::cerr << " " << fh;
          zz_free = false;
        }
      }
      else
      {
        if (v_trans[0] == v_trans[1] || v_trans[0] == v_trans[2] || v_trans[1] == v_trans[2])
        {
          std::cerr << " " << fh;
          zz_free = false;
        }
      }
    }
  }
  std::cerr << std::endl;

  return zz_free;
}

template<class MeshT>
int SplitApproachT<MeshT>::n_complex_singular_edges()
{
  std::cerr << "Complex singular edges left: ";

  int n = 0;
  //check
  for (const auto eh: mesh_.edges())
  {
    if (valence_[eh] == std::numeric_limits<int>::max())
    {
      std::cerr << " " << eh;
      n++;
    }
  }


  std::cerr << std::endl;

  return n;
}


template<class MeshT>
void SplitApproachT<MeshT>::remove_zigzag()
{
  std::cerr << "###Removing zigzags..." << std::endl;
  std::queue<FH> zz_fhs;

  for (const auto fh: mesh_.faces())
  {
    if (feature_fprop_[fh] == 0 && is_removable_zigzag(fh))
      zz_fhs.push(fh);
  }

  int n = 0;
  while (!zz_fhs.empty() && n < 100000)
  {
    auto fh_cur = zz_fhs.front();
    zz_fhs.pop();

    if (feature_fprop_[fh_cur] > 0)
      continue;

    if (is_removable_zigzag(fh_cur))
    {
      HFH hf0 = mesh_.halfface_handle(fh_cur, 0);
      auto hes = mesh_.halfface(hf0).halfedges();

      std::vector<HEH> sg_hes;
      for (const auto hehi: hes)
        if (valence_[mesh_.edge_handle(hehi)] == 1 || valence_[mesh_.edge_handle(hehi)] == -1)
          sg_hes.push_back(hehi);

      int old_trans = trans_prop_[hf0];
      std::vector<int> e_trans;

      if (!feature_face_edge_[mesh_.edge_handle(sg_hes[0])] && !feature_face_edge_[mesh_.edge_handle(sg_hes[1])])
      {
        for (const auto &hehi: sg_hes)
          e_trans.push_back(halfedge_transition_index(mesh_, tq_, trans_prop_, hehi, hf0));
      }
      else
      {
        //TODO: optimize
        for (int i = 0; i < 24; ++i)
          e_trans.push_back(i);
      }

      bool accept = false;
      for (const auto trans: e_trans)
      {
        trans_prop_[hf0] = tq_.mult_transitions_idx(old_trans, tq_.inverse_transition_idx(trans));
        trans_prop_[mesh_.opposite_halfface_handle(hf0)] = tq_.inverse_transition_idx(trans_prop_[hf0]);

        int n_regular = 0, n_complex = 0;
        HEH sg_heh(-1);
        int sg_val = 0;
        for (const auto &heh: hes)
        {
          int val = sge_.calculate_edge_valence(mesh_.edge_handle(heh));
          if (val == 0)
            n_regular++;
          else if (abs(val) != 1)
            n_complex++;
          else
          {
            sg_heh = heh;
            sg_val = val;
          }
        }
        if (n_complex == 0 && n_regular == 2)
        {
          bool aligned = true;
          for (const auto hehi: hes)
          {
            auto ehi = mesh_.edge_handle(hehi);
            if (feature_edge_[ehi] > 0)
            {
              aligned = fac_.check_field_alignment_at_feature_edge(ehi);
              if (!aligned)
              {
                aligned = false;
                break;
              }
            }
          }

          if (aligned)
          {
            accept = true;
            break;
          }

//                        auto sg_eh = mesh_.edge_handle(sg_heh);
//                        if (feature_edge_[sg_eh]) {
//                            bool aligned = fac_.check_field_alignment_at_feature_edge(sg_eh);
//                            if (aligned) {
//                                accept = true;
//                                break;
//                            }
//                        }
//                        else {
//                            accept = true;
//                            break;
//                        }
        }
      }

      if (!accept)
      {
        trans_prop_[hf0] = old_trans;
        trans_prop_[mesh_.opposite_halfface_handle(hf0)] = tq_.inverse_transition_idx(trans_prop_[hf0]);

        continue;
      }


      n++;

      //update valence
      for (const auto heh: hes)
        sge_.compute_edge_valence(mesh_.edge_handle(heh));

      //update singular vertex
      for (auto fv_it = mesh_.fv_iter(fh_cur); fv_it.valid(); ++fv_it)
        this->update_singular_vertex_property(*fv_it);

      //update quaternions
      std::set<EH> ehs;
      for (auto fe_it = mesh_.fe_iter(fh_cur); fe_it.valid(); ++fe_it)
        ehs.insert(*fe_it);
      QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                         feature_fprop_, feature_edge_, ehs);

      //push neigbour faces
      for (const auto heh: hes)
      {
        for (auto hef_it = mesh_.hef_iter(heh); hef_it.valid(); ++hef_it)
          zz_fhs.push(*hef_it);
      }
    }
  }


  //check
  std::cerr << "Warnings above are harmless... " << std::endl;
  std::cerr << "Removed " << n << "zigzags. " << "Zigzag face left: " << std::endl;
  for (const auto fh: mesh_.faces())
  {
    if (feature_fprop_[fh] == 0 && is_removable_zigzag(fh))
      std::cerr << " " << fh;
  }
  std::cerr << std::endl;
}

template<class MeshT>
void SplitApproachT<MeshT>::remove_zigzag2()
{
  std::cerr << "###Removing zigzags..." << std::endl;
  std::queue<FH> zz_fhs;

  for (const auto fh: mesh_.faces())
  {
    if (n_singular_edges_in_face(fh) >= 2)
      zz_fhs.push(fh);
  }

  int n = 0;
  while (!zz_fhs.empty() && n < 100000)
  {
    auto fh_cur = zz_fhs.front();
    zz_fhs.pop();

    if (feature_fprop_[fh_cur] > 0)
      continue;

    int n_sge = n_singular_edges_in_face(fh_cur);
    if (n_sge >= 2)
    {
      HFH hf0 = mesh_.halfface_handle(fh_cur, 0);
      auto hes = mesh_.halfface(hf0).halfedges();

      std::vector<int> v_trans;
      for (auto hehi: hes)
      {
        auto ehi = mesh_.edge_handle(hehi);
        if (valence_[ehi] != 0)
        {

          int trans = EdgeMonodromyHelperT::halfedge_transition_index_in_cell(mesh_, tq_, cell_quaternions_,
                                                                              trans_prop_, valence_, hehi, hf0);
          v_trans.push_back(trans);

//                        std::cerr << " eh " << ehi << " trans " << trans<<" in cell "<<mesh_.incident_cell(mesh_.opposite_halfface_handle(hf0));
        }
      }

      int adjust_trans = 0;
      if (n_sge == 2)
      {
        if (v_trans[0] != v_trans[1])
          continue;

        adjust_trans = v_trans[0];
      }
      else
      {
        if (v_trans[0] != v_trans[1] && v_trans[2] != v_trans[1] && v_trans[0] != v_trans[2])
          continue;

        if (v_trans[0] == v_trans[1] || v_trans[0] == v_trans[2])
          adjust_trans = v_trans[0];
        else if (v_trans[1] == v_trans[2])
          adjust_trans = v_trans[1];
      }
      int old_trans = trans_prop_[hf0];
      trans_prop_[hf0] = tq_.mult_transitions_idx(old_trans, tq_.inverse_transition_idx(adjust_trans));
      trans_prop_[mesh_.opposite_halfface_handle(hf0)] = tq_.inverse_transition_idx(trans_prop_[hf0]);


      n++;

      //update valence
      for (const auto heh: hes)
        sge_.compute_edge_valence(mesh_.edge_handle(heh));

      //update singular vertex
      for (auto fv_it = mesh_.fv_iter(fh_cur); fv_it.valid(); ++fv_it)
        this->update_singular_vertex_property(*fv_it);

      //update quaternions
      std::set<EH> ehs;
      for (auto fe_it = mesh_.fe_iter(fh_cur); fe_it.valid(); ++fe_it)
        ehs.insert(*fe_it);
      QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                         feature_fprop_, feature_edge_, ehs);

      //push neigbour faces
      for (const auto heh: hes)
      {
        for (auto hef_it = mesh_.hef_iter(heh); hef_it.valid(); ++hef_it)
          zz_fhs.push(*hef_it);
      }
    }
  }


  //check
  std::cerr << "Warnings above are harmless... " << std::endl;
  std::cerr << "Removed " << n << "zigzags." << std::endl;
  check_zigzag();
}


template<class MeshT>
void
SplitApproachT<MeshT>::split_complex_singular_edges()
{
  LocalMeshabilityRepairT<MeshT> flnmv(mesh_);
  DetachAtInvalidSingularNodeT<MeshT> fisn(mesh_);

  std::queue<VH> que;
  for (const auto vh: mesh_.vertices())
  {
    if (sgl_vt_[vh] != 0)
    {
      if (n_incident_complex_singular_edges(vh) == 1)
      {
        que.push(vh);
      }
    }
  }
  int n = 0;

  while (!que.empty() && n < 100000)
  {
    auto vh_cur = que.front();
    que.pop();

    int n_cpe = n_incident_complex_singular_edges(vh_cur);
    if (n_cpe == 1)
    {
      EH cpl_eh(-1);
      for (auto ve_it = mesh_.ve_iter(vh_cur); ve_it.valid(); ++ve_it)
        if (valence_[*ve_it] == std::numeric_limits<int>::max())
          cpl_eh = *ve_it;
      bool suc = eem_.fix_interior_compound_singular_edge(mesh_.halfedge_handle(cpl_eh, 0), true);

      if (suc)
        fisn.fix_interior_invalid_singular_node(vh_cur);
      else
      {
        std::cerr << "Failed to fix complex sg edge " << std::endl;

        continue;
      }

      int n_after = n_incident_complex_singular_edges(vh_cur);

      if (n_after > 0)
        std::cerr << "after fix still complex: " << vh_cur << std::endl;
      if (n_after == 0)
        n++;

      //push candidates
      if (sgl_vt_[vh_cur])
        que.push(vh_cur);
      for (auto vv_it = mesh_.vv_iter(vh_cur); vv_it.valid(); ++vv_it)
        if (sgl_vt_[*vv_it] != 0)
          que.push(*vv_it);
    }
  }

  mesh_.collect_garbage();

  std::cerr << "Removed " << n << " complex singular edges." << std::endl;
}


template<class MeshT>
int
SplitApproachT<MeshT>::n_incident_complex_singular_edges(const VH &_vh) const
{
  int n_sg = 0;
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (valence_[*ve_it] == std::numeric_limits<int>::max())
      n_sg++;

  return n_sg;
}

template<class MeshT>
int
SplitApproachT<MeshT>::n_legal_singular_edges_in_face(const FH &_fh) const
{
  int n_sg = 0;
  for (auto fe_it = mesh_.fe_iter(_fh); fe_it.valid(); ++fe_it)
    if (valence_[*fe_it] == -1 || valence_[*fe_it] == 1)
      n_sg++;

  return n_sg;
}

template<class MeshT>
bool
SplitApproachT<MeshT>::is_removable_zigzag(const FH &_fh) const
{
  if (n_legal_singular_edges_in_face(_fh) == 2)
  {
    for (auto fe_it = mesh_.fe_iter(_fh); fe_it.valid(); ++fe_it)
    {
      if (feature_face_edge_[*fe_it] && !mesh_.is_boundary(*fe_it)) //feature face sector maybe become incorrect
        return false;
    }

    return true;
  }

  return false;
}

template<class MeshT>
int
SplitApproachT<MeshT>::n_singular_edges_in_face(const FH &_fh) const
{
  int n = 0;
  for (auto fe_it = mesh_.fe_iter(_fh); fe_it.valid(); ++fe_it)
  {
    if (valence_[*fe_it] != 0)
      n++;
  }

  return n;
}


}
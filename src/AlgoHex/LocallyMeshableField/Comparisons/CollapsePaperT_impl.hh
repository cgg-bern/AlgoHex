/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#include "CollapsePaperT.hh"

namespace AlgoHex
{
template<class MeshT>
void
CollapseApproachT<MeshT>::collapse_approach()
{
  matching_adjustment();
  collapse_complex_singular_edges();
}

template<class MeshT>
int CollapseApproachT<MeshT>::n_complex_singular_edges() const
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
void
CollapseApproachT<MeshT>::matching_adjustment()
{
  //1. greedily remove complex singular edges
  int n_ce_removed = 0;
  for (const auto eh: mesh_.edges())
  {
    if (valence_[eh] == std::numeric_limits<int>::max())
    {
      auto he0 = mesh_.halfedge_handle(eh, 0);
      for (auto hehf_it = mesh_.hehf_iter(he0); hehf_it.valid(); ++hehf_it)
      {
        auto old_trans = trans_prop_[*hehf_it];

        auto hes = mesh_.halfface(*hehf_it).halfedges();
        std::vector<int> e_trans;

        auto dis = sorted_matchings(*hehf_it);

        bool accept = false;
        for (const auto &di: dis)
        {
          trans_prop_[*hehf_it] = tq_.mult_transitions_idx(old_trans, tq_.inverse_transition_idx(di.second));
          trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = tq_.inverse_transition_idx(trans_prop_[*hehf_it]);

          int n_regular = 0, n_complex = 0;
          for (const auto &heh: hes)
          {
            int val = sge_.calculate_edge_valence(mesh_.edge_handle(heh));
            if (val == 0)
              n_regular++;
            if (val == std::numeric_limits<int>::max())
              n_complex++;
          }
          if (n_complex == 0 && n_regular > 0)
          {
            accept = true;
            break;
          }
        }

        if (!accept)
        {
          trans_prop_[*hehf_it] = old_trans;
          trans_prop_[mesh_.opposite_halfface_handle(*hehf_it)] = tq_.inverse_transition_idx(trans_prop_[*hehf_it]);
        }
        else
        {
          for (const auto &heh: hes)
            sge_.compute_edge_valence(mesh_.edge_handle(heh));

          //update singular vertex
          for (auto hfv_it = mesh_.hfv_iter(*hehf_it); hfv_it.valid(); ++hfv_it)
            this->update_singular_vertex_property(*hfv_it);

          //update quaternions
          std::set<EH> ehs;
          for (auto hfe_it = mesh_.hfe_iter(*hehf_it); hfe_it.valid(); ++hfe_it)
            ehs.insert(*hfe_it);
          QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                             feature_fprop_, feature_edge_, ehs);

          n_ce_removed++;

          break;
        }
      }
    }
  }


  std::cerr << "Removed " << n_ce_removed << "complex edges by matching adjustment. " << "Complex edges left: "
            << std::endl;
  //check
  for (const auto eh: mesh_.edges())
  {
    if (valence_[eh] == std::numeric_limits<int>::max())
    {
      std::cerr << " " << eh;
    }
  }
  std::cerr << std::endl;


  //2. zigzag removal
  remove_zigzag();
}


template<class MeshT>
void CollapseApproachT<MeshT>::remove_zigzag()
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

//                std::cerr<<" adj trans "<<adjust_trans<<std::endl;
      trans_prop_[hf0] = tq_.mult_transitions_idx(old_trans, tq_.inverse_transition_idx(adjust_trans));
      trans_prop_[mesh_.opposite_halfface_handle(hf0)] = tq_.inverse_transition_idx(trans_prop_[hf0]);

      int n_regular = 0, n_complex = 0;
      for (const auto &heh: hes)
      {
        int val = sge_.calculate_edge_valence(mesh_.edge_handle(heh));
        if (val == 0)
          n_regular++;
        if (val == std::numeric_limits<int>::max())
          n_complex++;

      }

      if (n_complex == 0 && n_regular > 0)
      {
        n++;
        //update valence
        for (const auto &heh: hes)
          sge_.compute_edge_valence(mesh_.edge_handle(heh));

        //update singular vertex
        for (auto fv_it = mesh_.fv_iter(fh_cur); fv_it.valid(); ++fv_it)
          this->update_singular_vertex_property(*fv_it);

        //update quaternions
        std::set<EH> ehs;
        for (auto fe_it = mesh_.fe_iter(fh_cur); fe_it.valid(); ++fe_it)
          ehs.insert(*fe_it);
        QuaternionSmoothing::optimize_quaternions_at_edges(mesh_, trans_prop_, tq_, cell_quaternions_, valence_,
                                                           feature_fprop_,
                                                           feature_edge_, ehs);

        //push neigbour faces
        for (const auto &heh: hes)
        {
          for (auto hef_it = mesh_.hef_iter(heh); hef_it.valid(); ++hef_it)
            zz_fhs.push(*hef_it);
        }
      }
      else
      {
        trans_prop_[hf0] = old_trans;
        trans_prop_[mesh_.opposite_halfface_handle(hf0)] = tq_.inverse_transition_idx(trans_prop_[hf0]);
      }
    }
  }


  //check
  std::cerr << "Warnings above are harmless... " << std::endl;
  std::cerr << "Removed " << n << "zigzags." << std::endl;
}


template<class MeshT>
void
CollapseApproachT<MeshT>::collapse_complex_singular_edges()
{

  //find all complex singular edges
  std::queue<HEH> que;
  for (const auto eh: mesh_.edges())
  {
    if (valence_[eh] == std::numeric_limits<int>::max())
    {
      que.push(mesh_.halfedge_handle(eh, 0));
      que.push(mesh_.halfedge_handle(eh, 1));
    }
  }

  int n = 0;
  while (!que.empty() && n < 10000)
  {
    auto he_cur = que.front();
    que.pop();

//            std::cerr<<"collapsing eh "<<mesh_.edge_handle(he_cur)<<" "<<mesh_.halfedge(he_cur).from_vertex()<<" "<<mesh_.halfedge(he_cur).to_vertex()<<std::endl;
    if (mesh_.is_deleted(he_cur))
      continue;

    auto eh_cur = mesh_.edge_handle(he_cur);
    if (valence_[eh_cur] == std::numeric_limits<int>::max())
    {
      //step 1
      std::vector<EH> split_ehs;
      for (auto ef_it = mesh_.ef_iter(eh_cur); ef_it.valid(); ++ef_it)
      {
        std::vector<EH> proper_sg_ehs;
        for (auto fe_it = mesh_.fe_iter(*ef_it); fe_it.valid(); ++fe_it)
          if (valence_[*fe_it] == -1 || valence_[*fe_it] == 1)
            proper_sg_ehs.push_back(*fe_it);

        if (proper_sg_ehs.size() == 2)
        {
          split_ehs.push_back(proper_sg_ehs[0]);
        }
      }
      for (const auto &eh: split_ehs)
      {
        if (es_.is_split_ok(eh) && valence_[eh] != std::numeric_limits<int>::max())
        {
          auto vh_sp = es_.edge_split(eh);
        }
      }


      //step 2
      auto collapse_type = ec_.is_topology_ok(he_cur);

      if (collapse_type == ec_.LinkError)
      {
        std::cerr << "link problem " << std::endl;

        std::set<EH> blocking_ehs;
        auto blk_ehs = ec_.get_blocking_edges(eh_cur);
        blocking_ehs.insert(blk_ehs.begin(), blk_ehs.end());

        //split
        for (const auto &eh: blocking_ehs)
        {
          if (es_.is_split_ok(eh) && valence_[eh] != std::numeric_limits<int>::max())
          {
            auto vh_sp = es_.edge_split(eh);
            std::cerr << "split blocking edge " << eh << std::endl;
          }
        }

        collapse_type = ec_.is_topology_ok(he_cur);
        std::cerr << "after split the collapse type: " << collapse_type << std::endl;
      }

      if (mesh_.is_deleted(he_cur))
        continue;

      if (collapse_type == ec_.CollapseOK)
      {
//                    std::cerr<<"collapse ok, halfedge "<<mesh_.halfedge(he_cur).from_vertex()<<" "<<mesh_.halfedge(he_cur).to_vertex()<<std::endl;
        VH vh = ec_.edge_collapse_allow_inverted(he_cur);

        if (vh.is_valid())
        {
          //update and push
          for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
          {
            if (!mesh_.is_deleted(*ve_it))
            {
              //update valence
//                                auto vhff = mesh_.halfedge(mesh_.halfedge_handle(*ve_it, 0)).from_vertex();
//                                auto vhtt = mesh_.halfedge(mesh_.halfedge_handle(*ve_it, 0)).to_vertex();
              sge_.compute_edge_valence(*ve_it);

              que.push(mesh_.halfedge_handle(*ve_it, 0));
              que.push(mesh_.halfedge_handle(*ve_it, 1));
            }
          }

          //update
          for (auto vv_it = mesh_.vv_iter(vh); vv_it.valid(); ++vv_it)
            this->update_singular_vertex_property(*vv_it);
          this->update_singular_vertex_property(vh);

          //relocate
//                        vo_.vertex_optimize(vh);

          n++;
//                        std::cerr<<"collapsed eh "<<mesh_.edge_handle(he_cur)<<" "<<mesh_.halfedge(he_cur).from_vertex()<<" "<<mesh_.halfedge(he_cur).to_vertex()<<std::endl;

        }
//                    else
//                        std::cerr<<"Failed to collapse edge "<<eh_cur<<std::endl;
      }
      else
      {
        std::cerr << "Error: split didn't solve link error " << eh_cur << std::endl;
      }
    }
  }


  std::cerr << "Collapsed " << n << "complex edges. " << "Complex edges left: " << std::endl;
  //check
  for (const auto eh: mesh_.edges())
  {
    if (valence_[eh] == std::numeric_limits<int>::max())
    {
      std::cerr << " " << eh;
    }
  }
  std::cerr << std::endl;

  mesh_.collect_garbage();
}

template<class MeshT>
void CollapseApproachT<MeshT>::collapse_complex_singular_edge(const HEH &he_cur)
{
  auto eh_cur = mesh_.edge_handle(he_cur);
  if (valence_[eh_cur] == std::numeric_limits<int>::max())
  {
    //step 1
    std::vector<EH> split_ehs;
    for (auto ef_it = mesh_.ef_iter(eh_cur); ef_it.valid(); ++ef_it)
    {
      std::vector<EH> proper_sg_ehs;
      for (auto fe_it = mesh_.fe_iter(*ef_it); fe_it.valid(); ++fe_it)
        if (valence_[*fe_it] == -1 || valence_[*fe_it] == 1)
          proper_sg_ehs.push_back(*fe_it);

      if (proper_sg_ehs.size() == 2)
      {
        split_ehs.push_back(split_ehs[0]);
      }
    }
    for (const auto &eh: split_ehs)
    {
      if (es_.is_split_ok(eh))
      {
        auto vh_sp = es_.edge_split(eh);
      }
    }


    //step 2
    auto collapse_type = ec_.is_collapse_ok(he_cur);

    if (collapse_type == ec_.LinkError)
    {
      std::set<EH> blocking_ehs;
      auto blk_ehs = ec_.get_blocking_edges(eh_cur);
      blocking_ehs.insert(blk_ehs.begin(), blk_ehs.end());

      //split
      for (const auto &eh: blocking_ehs)
      {
        if (es_.is_split_ok(eh))
        {
          auto vh_sp = es_.edge_split(eh);
        }
      }

      collapse_type = ec_.is_collapse_ok(he_cur);
      std::cerr << "after split the collapse type: " << collapse_type << std::endl;
    }


    if (collapse_type == ec_.CollapseOK)
    {
      VH vh = ec_.edge_collapse(he_cur, false);
      if (vh.is_valid())
      {

        //update and push
        for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
        {
          if (!mesh_.is_deleted(*ve_it))
          {
            //update valence
            sge_.compute_edge_valence(*ve_it);
          }
        }

        //update
        for (auto vv_it = mesh_.vv_iter(vh); vv_it.valid(); ++vv_it)
          this->update_singular_vertex_property(*vv_it);
        this->update_singular_vertex_property(vh);

        //relocate
        vo_.vertex_optimize(vh);
      }
    }
  }

  mesh_.collect_garbage();
}

template<class MeshT>
std::vector<std::pair<double, int>>
CollapseApproachT<MeshT>::sorted_matchings(const HFH &_hfh) const
{
  std::vector<std::pair<double, int>> di;
  if (mesh_.is_boundary(mesh_.face_handle(_hfh)))
    return di;

  for (int i = 1; i < 24; ++i)
  {
    auto hfh1 = mesh_.opposite_halfface_handle(_hfh);

    auto ch0 = mesh_.incident_cell(_hfh);
    auto ch1 = mesh_.incident_cell(hfh1);

    Quaternion tranq = tq_.transition(trans_prop_[_hfh]);

    Quaternion q0in1 = cell_quaternions_[ch0] * tranq;
    Mat3d mt0in1 = q0in1.conjugate().toRotationMatrix();
    Mat3d mt1 = cell_quaternions_[ch1].toRotationMatrix();

    Mat3d mt = mt1 - mt0in1;

    di.emplace_back(mt.norm(), i);
  }

  std::sort(di.begin(), di.end());

  return di;
}


template<class MeshT>
int
CollapseApproachT<MeshT>::n_legal_singular_edges_in_face(const FH &_fh) const
{
  int n_sg = 0;
  for (auto fe_it = mesh_.fe_iter(_fh); fe_it.valid(); ++fe_it)
    if (valence_[*fe_it] == -1 || valence_[*fe_it] == 1)
      n_sg++;

  return n_sg;
}

template<class MeshT>
int
CollapseApproachT<MeshT>::n_singular_edges_in_face(const FH &_fh) const
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
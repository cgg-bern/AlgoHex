/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define MESHOPTIMIZATIONT_C

#include "MeshOptimizationT.hh"
#include <AlgoHex/BinarySpacePartitionTrees/Algorithms.hh>
#include "CommonFuncs.hh"
#include "SingularVertexOptProblem.hh"


namespace AlgoHex
{
template<class MeshT>
void MeshOptimizationT<MeshT>::
remesh(const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w, const double _rp_w,
       const int k)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "##Remesh ..." << std::endl;)

  {
    check_cells_volume(mesh_);

    ScopedStopWatch ssw(sw::edge_split);
    auto ns = split_edges(k);
  }

  {
    ScopedStopWatch ssw(sw::edge_collapse);
    auto nc = collapse_edges(k);
  }

  {
    ScopedStopWatch ssw(sw::edge_swap);
    auto nsw = swap_edges(k);
  }


  {
    ScopedStopWatch ssw(sw::singularity_smooth);
    check_cells_volume(mesh_);

    singular_vertices_relocation(_sge_w / scale_, _rgl_w / scale_, _cs_w / scale_, _tps_w / scale_, _rp_w / scale_,
                                 100);

//            fac_.check_field_alignments_at_feature_faces();
//            fac_.check_field_alignment_at_feature_edges();
  }


  iter_++;

  ALGOHEX_DEBUG_ONLY(std::cerr << "\n#Tet mesh data after optimization";
                             std::cerr << "\n#Tet mesh vertex number: " << mesh_.n_vertices();
                             std::cerr << "\n#Tet mesh edge number: " << mesh_.n_edges();
                             std::cerr << "\n#Tet mesh face number: " << mesh_.n_faces();
                             std::cerr << "\n#Tet mesh cell number: " << mesh_.n_cells() << std::endl);
}

template<class MeshT>
std::vector<EH> MeshOptimizationT<MeshT>::
get_singular_edge_pairs(const int _pair_id) const
{
  std::vector<EH> sg_ehs;
  EH eh_s1(-1), eh_s2(-1);
  for (auto ehi: mesh_.edges())
  {
    if (sg_edge_pairs_[ehi] == _pair_id)
    {
      eh_s1 = ehi;
    }
    if (sg_edge_pairs_[ehi] == _pair_id + 1)
    {
      eh_s2 = ehi;
    }

    if (eh_s1.is_valid() && eh_s2.is_valid())
      break;
  }
  if (eh_s1.is_valid())
    sg_ehs = get_edges_on_singular_arc(mesh_, valence_, mesh_.halfedge_handle(eh_s1, 0));

  if (eh_s2.is_valid())
  {
    auto sg_ehs2 = get_edges_on_singular_arc(mesh_, valence_, mesh_.halfedge_handle(eh_s2, 0));
    sg_ehs.insert(sg_ehs.end(), sg_ehs2.begin(), sg_ehs2.end());
  }

  return sg_ehs;
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
remesh(int _pair_id, const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w,
       const double _rp_w, const int k)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "##Remesh ..." << std::endl;)

  {
    std::vector<EH> sg_ehs = get_singular_edge_pairs(_pair_id);
    auto ns = split_edges(sg_ehs, k);
  }

  {
    std::vector<EH> sg_ehs = get_singular_edge_pairs(_pair_id);

    auto nc = collapse_edges(sg_ehs, k);
  }

  {
    std::vector<EH> sg_ehs = get_singular_edge_pairs(_pair_id);

    auto nsw = swap_edges(sg_ehs, k);
  }

  mesh_.collect_garbage();

  {
    std::vector<EH> sg_ehs = get_singular_edge_pairs(_pair_id);
    singular_vertices_relocation(sg_ehs, _sge_w / scale_, _rgl_w / scale_, _cs_w / scale_, _tps_w / scale_,
                                 _rp_w / scale_, 100);
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
remesh(const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w, const double _rp_w,
       const int k, std::vector<MeshT> &meshes)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "##Remesh ..." << std::endl;)

  {
    ScopedStopWatch ssw(sw::edge_split);
    auto ns = split_edges(k);

//            save meshes for visualization
//                meshes.push_back(mesh_);

  }


  {
    ScopedStopWatch ssw(sw::edge_collapse);
    auto nc = collapse_edges(k);
//                save meshes for visualization
//                meshes.push_back(mesh_);
  }


  {
    ScopedStopWatch ssw(sw::edge_swap);
    auto nsw = swap_edges(k);

    //save meshes for visualization
//            meshes.push_back(mesh_);
  }

//meshes.push_back(mesh_);

  {
    check_cells_volume(mesh_);

    ScopedStopWatch ssw(sw::singularity_smooth);
    singular_vertices_relocation(_sge_w / scale_, _rgl_w / scale_, _cs_w / scale_, _tps_w / scale_, _rp_w / scale_,
                                 100);
  }

  iter_++;
  ALGOHEX_DEBUG_ONLY(std::cerr << "\n#Tet mesh data after optimization";
                             std::cerr << "\n#Tet mesh vertex number: " << mesh_.n_vertices();
                             std::cerr << "\n#Tet mesh edge number: " << mesh_.n_edges();
                             std::cerr << "\n#Tet mesh face number: " << mesh_.n_faces();
                             std::cerr << "\n#Tet mesh cell number: " << mesh_.n_cells() << std::endl;)
}


template<class MeshT>
void MeshOptimizationT<MeshT>::remesh_post(const int _remesh_iter)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "\n#Tet mesh data before optimization";
                             std::cerr << "\n#Tet mesh vertex number: " << mesh_.n_vertices();
                             std::cerr << "\n#Tet mesh edge number: " << mesh_.n_edges();
                             std::cerr << "\n#Tet mesh face number: " << mesh_.n_faces();
                             std::cerr << "\n#Tet mesh cell number: " << mesh_.n_cells() << std::endl;)

  auto smoothness_vprop = mesh_.template request_vertex_property<double>("field_smoothness", 0.);
  mesh_.set_persistent(smoothness_vprop, true);

  set_min_target_length(0.3);

  double damping = 1.;

  std::cerr << "###ERROR message on numerical issues in KKT system can be ignored!" << std::endl;

  for (int i = 0; i < _remesh_iter; ++i)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "##Remesh-post: " << i << " ..." << std::endl;)

    //adapt target length of the tetmesh
    if (i < 3)
    {
      damping *= 1.;
      adapt_target_length(damping, 8);
    }
    else
      measure_field_smoothness(smoothness_vprop, 8);

    auto ns = split_edges(-1, true);

    auto nc = collapse_edges(-1, true, true);

    auto nsw = swap_edges(-1, true);

    //            optimize_singular_vertices();

//            optimize_key_vertices();
    optimize_regular_vertices();
  }

  mesh_.collect_garbage();

  check_mesh_quality(mesh_);

  ALGOHEX_DEBUG_ONLY(std::cerr << "\n#Tet mesh data after optimization";
                             std::cerr << "\n#Tet mesh vertex number: " << mesh_.n_vertices();
                             std::cerr << "\n#Tet mesh edge number: " << mesh_.n_edges();
                             std::cerr << "\n#Tet mesh face number: " << mesh_.n_faces();
                             std::cerr << "\n#Tet mesh cell number: " << mesh_.n_cells() << std::endl;)
}

template<class MeshT>
void MeshOptimizationT<MeshT>::remesh_post(const int _remesh_iter, std::vector<MeshT> &_meshes)
{
  double damping = 1.;
  for (int i = 0; i < _remesh_iter; ++i)
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "##Remesh-post: " << i << " ..." << std::endl;)

    auto ns = split_edges(-1, true);
//            _meshes.push_back(mesh_);

    auto nc = collapse_edges(-1, true, true);
//            _meshes.push_back(mesh_);

    auto nsw = swap_edges(-1, true);
//            check_valence_change();
//            _meshes.push_back(mesh_);


//            optimize_singular_vertices();
//            optimize_key_vertices();
//
//            _meshes.push_back(mesh_);


    optimize_regular_vertices();
//            check_valence_change();
//            _meshes.push_back(mesh_);
  }

  mesh_.collect_garbage();

  check_mesh_quality(mesh_);

  ALGOHEX_DEBUG_ONLY(std::cerr << "\n#Tet mesh data after optimization";
                             std::cerr << "\n#Tet mesh vertex number: " << mesh_.n_vertices();
                             std::cerr << "\n#Tet mesh edge number: " << mesh_.n_edges();
                             std::cerr << "\n#Tet mesh face number: " << mesh_.n_faces();
                             std::cerr << "\n#Tet mesh cell number: " << mesh_.n_cells() << std::endl;)
}


template<class MeshT>
void MeshOptimizationT<MeshT>::check_valence_change()
{
  for (const auto ehi: mesh_.edges())
  {
    int val_tmp = sge_.calculate_edge_valence(ehi);
    if (val_tmp != valence_[ehi])
      std::cerr << "Error: valence of edge " << ehi << " changed from " << valence_[ehi] << " to " << val_tmp
                << std::endl;
  }
}


template<class MeshT>
int MeshOptimizationT<MeshT>::
split_edges(const int k, const bool _split_with_opt)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Splitting edges ..." << std::endl;)
  EdgeSplitT<MeshT> es(mesh_);
  double threshold = 5. / 3.;

  //for post remeshing
  if (_split_with_opt)
  {
    es.set_split_with_optimization();
    threshold /= 2.;
  }

  std::set<PairDE, std::greater<PairDE>> kring_ehs;
  if (k == -1)
  {
    for (const auto eh: mesh_.edges())
      kring_ehs.emplace(mesh_.length(eh), eh);
  }
  else if (k > 0)
    kring_ehs = get_k_ring_edges_of_singular_graph<MeshT, std::greater<PairDE>>(mesh_, k);

  int n = 0;
  while (!kring_ehs.empty() && n < 1000000)
  {
    auto de_cur = *kring_ehs.begin();
    kring_ehs.erase(kring_ehs.begin());

    if (es.is_split_ok(de_cur.second))
    {
      if (de_cur.first != mesh_.length(de_cur.second))
        continue;

      auto vh_f = mesh_.edge(de_cur.second).from_vertex();
      auto vh_t = mesh_.edge(de_cur.second).to_vertex();

      if (de_cur.first < threshold * (target_length_[vh_f] + target_length_[vh_t]) / 2.)
        continue;

      auto vh = es.edge_split(de_cur.second, true);
      if (vh.is_valid())
      {
        n++;

        //push
        for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
          kring_ehs.emplace(mesh_.length(*ve_it), *ve_it);

        //push the neigbouring edges in case infinite loop occurs
        for (auto ve_it = mesh_.ve_iter(vh_f); ve_it.valid(); ++ve_it)
        {
          kring_ehs.emplace(mesh_.length(*ve_it), *ve_it);
        }

        for (auto ve_it = mesh_.ve_iter(vh_t); ve_it.valid(); ++ve_it)
        {
          kring_ehs.emplace(mesh_.length(*ve_it), *ve_it);
        }
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "Split " << n << " edges." << std::endl;)
  {
//            ScopedStopWatch ssw(sw::garbage_collect);
    mesh_.collect_garbage();
  }
//        check_mesh_quality(mesh_);

  return n;
}

template<class MeshT>
int MeshOptimizationT<MeshT>::
split_edges(const std::vector<EH> &_sg_ehs, const int k)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Splitting edges ..." << std::endl;)
  EdgeSplitT<MeshT> es(mesh_);
  double threshold = 5. / 3.;

  std::set<PairDE, std::greater<PairDE>> kring_ehs;
  kring_ehs = get_k_ring_edges<MeshT, std::greater<PairDE>>(mesh_, _sg_ehs, k);

  int n = 0;
  while (!kring_ehs.empty() && n < 1000000)
  {
    auto de_cur = *kring_ehs.begin();
    kring_ehs.erase(kring_ehs.begin());

    if (es.is_split_ok(de_cur.second))
    {
      if (de_cur.first != mesh_.length(de_cur.second))
        continue;

      auto vh_f = mesh_.edge(de_cur.second).from_vertex();
      auto vh_t = mesh_.edge(de_cur.second).to_vertex();

      if (de_cur.first < threshold * (target_length_[vh_f] + target_length_[vh_t]) / 2.)
        continue;

      auto vh = es.edge_split(de_cur.second, true);
      if (vh.is_valid())
      {
        n++;

        //push
        for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
          kring_ehs.emplace(mesh_.length(*ve_it), *ve_it);

        //push the neigbouring edges in case infinite loop occurs
        for (auto ve_it = mesh_.ve_iter(vh_f); ve_it.valid(); ++ve_it)
        {
          kring_ehs.emplace(mesh_.length(*ve_it), *ve_it);
        }

        for (auto ve_it = mesh_.ve_iter(vh_t); ve_it.valid(); ++ve_it)
        {
          kring_ehs.emplace(mesh_.length(*ve_it), *ve_it);
        }
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "Split " << n << " edges." << std::endl;)

  return n;
}


template<class MeshT>
int MeshOptimizationT<MeshT>::
collapse_edges(const int k, const bool _post_remesh, const bool _check_energy)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Collapsing edges ..." << std::endl;)

  double threshold = 3. / 5.;

  EdgeCollapseT<MeshT> ec(mesh_);
  if (_post_remesh)
  {
    ec.set_post_remeshing(_post_remesh);
    ec.set_collapse_with_optimization();
//            threshold /= 2.;
  }

  std::set<PairDHE, std::less<PairDHE>> kring_hehs;
  if (k == -1)
  {
    for (const auto eh: mesh_.edges())
    {
      double len = mesh_.length(eh);
      kring_hehs.emplace(len, mesh_.halfedge_handle(eh, 0));
      kring_hehs.emplace(len, mesh_.halfedge_handle(eh, 1));
    }
  }
  else if (k > 0)
  {
    auto kring_ehs = get_k_ring_edges_of_singular_graph<MeshT, std::less<PairDE>>(mesh_, k);
    for (const auto &deh: kring_ehs)
    {
      kring_hehs.emplace(deh.first, mesh_.halfedge_handle(deh.second, 0));
      kring_hehs.emplace(deh.first, mesh_.halfedge_handle(deh.second, 1));
    }
  }

  int n = 0;
  while (!kring_hehs.empty())
  {
    auto dheh_cur = *kring_hehs.begin();
    kring_hehs.erase(kring_hehs.begin());

    auto eh_cur = mesh_.edge_handle(dheh_cur.second);
    if (dheh_cur.first != mesh_.length(eh_cur))
      continue;

    auto vh_f = mesh_.halfedge(dheh_cur.second).from_vertex();
    auto vh_t = mesh_.halfedge(dheh_cur.second).to_vertex();

    if (dheh_cur.first < threshold * (target_length_[vh_f] + target_length_[vh_t]) / 2.)
    {
      auto collapse_type = ec.is_collapse_ok(dheh_cur.second);
      if (collapse_type == ec.CollapseOK)
      {
        VH vh = ec.edge_collapse(dheh_cur.second, _check_energy);
        if (vh.is_valid())
        {
          n++;

          //update edge length and push
          for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
          {
            if (!mesh_.is_deleted(*ve_it))
            {
              //update valence
              //in case valence updated inside edge collapse may be not enough
              if (!mesh_.is_boundary(*ve_it) || valence_[*ve_it] != 0)
              {
//                                    int old_val = valence_[*ve_it];
                sge_.compute_edge_valence(*ve_it);
              }

              double len = mesh_.length(*ve_it);
              kring_hehs.emplace(len, mesh_.halfedge_handle(*ve_it, 0));
              kring_hehs.emplace(len, mesh_.halfedge_handle(*ve_it, 1));
            }
          }


          for (auto vv_it = mesh_.vv_iter(vh); vv_it.valid(); ++vv_it)
            MeshPropertiesT<MeshT>::update_singular_vertex_property(*vv_it);
          MeshPropertiesT<MeshT>::update_singular_vertex_property(vh);
        }
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "Collapsed " << n << " edges." << std::endl;)
  {
//            ScopedStopWatch ssw(sw::garbage_collect);
    mesh_.collect_garbage();
  }

  return n;
}

template<class MeshT>
int MeshOptimizationT<MeshT>::
collapse_edges(const std::vector<EH> &_sg_ehs, const int k, const bool _check_energy)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Collapsing edges ..." << std::endl;)

  double threshold = 3. / 5.;

  EdgeCollapseT<MeshT> ec(mesh_);

  std::set<PairDHE, std::less<PairDHE>> kring_hehs;
  auto kring_ehs = get_k_ring_edges<MeshT, std::less<PairDE>>(mesh_, _sg_ehs, k);
  for (const auto &deh: kring_ehs)
  {
    kring_hehs.emplace(deh.first, mesh_.halfedge_handle(deh.second, 0));
    kring_hehs.emplace(deh.first, mesh_.halfedge_handle(deh.second, 1));
  }


//        std::vector<EH> unclp_ehs;
  int n = 0;
  while (!kring_hehs.empty())
  {
    auto dheh_cur = *kring_hehs.begin();
    kring_hehs.erase(kring_hehs.begin());

    auto eh_cur = mesh_.edge_handle(dheh_cur.second);
    if (dheh_cur.first != mesh_.length(eh_cur))
      continue;

    auto vh_f = mesh_.halfedge(dheh_cur.second).from_vertex();
    auto vh_t = mesh_.halfedge(dheh_cur.second).to_vertex();

    if (dheh_cur.first < threshold * (target_length_[vh_f] + target_length_[vh_t]) / 2.)
    {
      auto collapse_type = ec.is_collapse_ok(dheh_cur.second);
      if (collapse_type == ec.CollapseOK)
      {
        VH vh = ec.edge_collapse(dheh_cur.second, _check_energy);
        if (vh.is_valid())
        {
          n++;

          //update edge length and push
          for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
          {
            if (!mesh_.is_deleted(*ve_it))
            {
              //update valence
              //in case valence updated inside edge collapse may be not enough
              if (!mesh_.is_boundary(*ve_it) || valence_[*ve_it] != 0)
              {
//                                    int old_val = valence_[*ve_it];
                sge_.compute_edge_valence(*ve_it);
              }

              double len = mesh_.length(*ve_it);
              kring_hehs.emplace(len, mesh_.halfedge_handle(*ve_it, 0));
              kring_hehs.emplace(len, mesh_.halfedge_handle(*ve_it, 1));
            }
          }


          for (auto vv_it = mesh_.vv_iter(vh); vv_it.valid(); ++vv_it)
            MeshPropertiesT<MeshT>::update_singular_vertex_property(*vv_it);
          MeshPropertiesT<MeshT>::update_singular_vertex_property(vh);
        }
      }
    }
  }


  ALGOHEX_DEBUG_ONLY(std::cerr << "Collapsed " << n << " edges." << std::endl;)

  return n;
}


template<class MeshT>
int MeshOptimizationT<MeshT>::
swap_edges(const int k, const bool _post_remesh)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Swapping edges..." << std::endl;)


  int n = 0;

  EdgeSwapT<MeshT> es(mesh_);
  es.set_post_remeshing(_post_remesh);

  std::set<PairDE, std::greater<PairDE>> kring_ehs;
  if (k == -1)
  {
    for (const auto eh: mesh_.edges())
      kring_ehs.emplace(mesh_.length(eh), eh);
  }
  else if (k > 0)
    kring_ehs = get_k_ring_edges_of_singular_graph<MeshT, std::greater<PairDE>>(mesh_, k);


  while (!kring_ehs.empty())
  {
    auto de_cur = *kring_ehs.begin();
    kring_ehs.erase(kring_ehs.begin());

    if (es.is_swap_ok(de_cur.second))
    {
      auto new_ehs = es.edge_swap(de_cur.second);

      if (!new_ehs.empty())
      {
        n++;
        //push
        for (const auto eh: new_ehs)
          kring_ehs.emplace(mesh_.length(eh), eh);
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "Swapped " << n << " edges." << std::endl;)

  {
//            ScopedStopWatch ssw(sw::garbage_collect);
    mesh_.collect_garbage();
  }

  return n;
}

template<class MeshT>
int MeshOptimizationT<MeshT>::
swap_edges(const std::vector<EH> &_sg_ehs, const int k)
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Swapping edges..." << std::endl;)

  int n = 0;

  EdgeSwapT<MeshT> es(mesh_);

  std::set<PairDE, std::greater<PairDE>> kring_ehs;
  kring_ehs = get_k_ring_edges<MeshT, std::greater<PairDE>>(mesh_, _sg_ehs, k);


  while (!kring_ehs.empty())
  {
    auto de_cur = *kring_ehs.begin();
    kring_ehs.erase(kring_ehs.begin());

    if (es.is_swap_ok(de_cur.second))
    {
      auto new_ehs = es.edge_swap(de_cur.second);

      if (!new_ehs.empty())
      {
        n++;
        //push
        for (const auto eh: new_ehs)
          kring_ehs.emplace(mesh_.length(eh), eh);
      }
    }
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "Swapped " << n << " edges." << std::endl;)

  return n;
}


template<class MeshT>
void MeshOptimizationT<MeshT>::
optimize_regular_vertices()
{
  ALGOHEX_DEBUG_ONLY(std::cerr << "Optimizing regular vertex position..." << std::endl;)
  VertexOptimizeT<MeshT> vo(mesh_);
  vo.set_energy_type(VertexOptimizeT<MeshT>::WORST_QUALITY);

  int n = 0;
  for (const auto vh: mesh_.vertices())
  {
    //don't move singular vertex position
    if (sgl_vt_[vh] || feature_edge_vertex_[vh])
      continue;

    if (vo.is_vertex_optimize_ok(vh))
    {
      if (vo.vertex_optimize(vh))
        n++;
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << "Optimized " << n << " vertex positions." << std::endl;)
}


template<class MeshT>
std::map<VH, std::map<HFH, int>> MeshOptimizationT<MeshT>::
collect_normal_aligned_axis(const bool _post)
{
  std::map<VH, std::map<HFH, int>> v_aligned_axis;

  //map variable vertices to new indices
  std::set<VH> vhs;

  if (!_post)
  {
    for (const auto &eh: mesh_.edges())
    {
      if (valence_[eh] != 0 && !mesh_.is_boundary(eh))
      {
        auto vh0 = mesh_.edge(eh).from_vertex();
        auto vh1 = mesh_.edge(eh).to_vertex();

        if (feature_face_vertex_[vh0])
          vhs.insert(vh0);

        if (feature_face_vertex_[vh1])
          vhs.insert(vh1);
      }
    }
  }
  else
  {
    for (const auto vhi: mesh_.vertices())
    {
      if (sgl_vt_[vhi] || feature_edge_vertex_[vhi])
      {
        if (feature_face_vertex_[vhi])
          vhs.insert(vhi);
      }
    }
  }

  for (const auto vhi: vhs)
  {
    std::map<HFH, int> aligned_axis;

    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
    {
      for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
      {
        auto fhi = mesh_.face_handle(*chf_it);
        if (feature_fprop_[fhi] > 0)
        {
          auto nm = mesh_.normal(*chf_it);
          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], nm).second;

          aligned_axis.insert(std::make_pair(*chf_it, axis));
        }
      }
    }

    v_aligned_axis.insert(std::make_pair(vhi, aligned_axis));
  }

  return v_aligned_axis;
}

template<class MeshT>
std::map<VH, std::map<HFH, int>> MeshOptimizationT<MeshT>::
collect_normal_aligned_axis(const std::vector<EH> &_sg_ehs)
{
  std::map<VH, std::map<HFH, int>> v_aligned_axis;

  //map variable vertices to new indices
  std::set<VH> vhs;

  for (const auto &eh: _sg_ehs)
  {
    if (valence_[eh] != 0)
    {
      auto vh0 = mesh_.edge(eh).from_vertex();
      auto vh1 = mesh_.edge(eh).to_vertex();

      if (feature_face_vertex_[vh0])
        vhs.insert(vh0);

      if (feature_face_vertex_[vh1])
        vhs.insert(vh1);
    }
  }


  for (const auto &vhi: vhs)
  {
    std::map<HFH, int> aligned_axis;

    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
    {
      for (auto chf_it = mesh_.chf_iter(*vc_it); chf_it.valid(); ++chf_it)
      {
        auto fhi = mesh_.face_handle(*chf_it);
        if (feature_fprop_[fhi] > 0)
        {
          auto nm = mesh_.normal(*chf_it);
          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], nm).second;

          aligned_axis.insert(std::make_pair(*chf_it, axis));
        }
      }
    }

    v_aligned_axis.insert(std::make_pair(vhi, aligned_axis));
  }

  return v_aligned_axis;
}


template<class MeshT>
void MeshOptimizationT<MeshT>::
singular_vertices_relocation(const double _sge_w, const double _rgl_w, const double _cs_w, const double _tps_w,
                             const double _rp_w, const int _max_iter)
{
  //check
  std::vector<int> sgg(mesh_.n_vertices(), 0);
  for (const auto &ehi: mesh_.edges())
  {
    if (valence_[ehi] != 0)
    {
      if (!mesh_.is_boundary(ehi))
      {
        sgg[mesh_.edge(ehi).from_vertex().idx()] = 1;
        sgg[mesh_.edge(ehi).to_vertex().idx()] = 1;
      }
    }
  }
  for (const auto &ehi: mesh_.edges())
  {
    if (valence_[ehi] != 0)
    {
      if (mesh_.is_boundary(ehi))
      {
        sgg[mesh_.edge(ehi).from_vertex().idx()] = 2;
        sgg[mesh_.edge(ehi).to_vertex().idx()] = 2;
      }
    }
  }

  for (const auto &vh: mesh_.vertices())
    if (sgg[vh.idx()] != sgl_vt_[vh])
    {
      std::cout << "Error: singular vertex status not consistent. VH: " << vh << "orig: " << sgl_vt_[vh]
                << "->" << sgg[vh.idx()] << std::endl;
    }

  //check quaternion
  for (const auto chi: mesh_.cells())
  {
    double qtnm = cell_quaternions_[chi].norm();
    if (qtnm < 0.1 || qtnm > 1.1)
      std::cerr << "Error: bad cell quaternion " << cell_quaternions_[chi].coeffs() << std::endl;
  }

  //
  auto v_aligned_axis = collect_normal_aligned_axis();

  //map variable vertices to new indices
  std::map<VH, int> new_vh_idx;
  int num = 0;
  for (auto vhi: mesh_.vertices())
  {
    if (sgl_vt_[vhi] != 0)
      new_vh_idx[vhi] = num++;
  }

  if (new_vh_idx.empty())
    return;

  //initialize variable values
  std::vector<double> x;
  x.resize(num * 3, 0);
  for (auto &vi: new_vh_idx)
  {
    for (int i = 0; i < 3; ++i)
      x[3 * vi.second + i] = mesh_.vertex(vi.first)[i];
  }


  //all original points
  std::vector<double> orig_x;
  orig_x.reserve(mesh_.n_vertices() * 3);
  for (const auto &vh: mesh_.vertices())
  {
    for (int i = 0; i < 3; ++i)
      orig_x.push_back(mesh_.vertex(vh)[i]);
  }

  SingularVertexOptProblemQuaternion svop(orig_x, x);


  //add alignment energy
  add_alignment_energy(svop, new_vh_idx, _sge_w);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added alignment energy" << std::endl;)

  //add shrink energy
  add_edge_shrink_energy(svop, new_vh_idx, _cs_w, _tps_w);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added shrink energy" << std::endl;)

  //add curvature term on singular graph
  add_curvature_smooth_energy(svop, new_vh_idx, _rgl_w);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added curvature_smooth energy" << std::endl;)

  add_repulsion_energy(svop, new_vh_idx, _rp_w);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added repulsion energy" << std::endl;)

//        add_boundary_singular_vertex_damping_energy(svop, new_vh_idx, 0.1);
//        ALGOHEX_DEBUG_ONLY(std::cerr<<"added boundary singular vertex damping energy"<<std::endl;

  add_imrm_energy(svop, new_vh_idx, 0.001);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added imrm energy" << std::endl;)

  add_surface_imrm_energy(svop, new_vh_idx, 0.05);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added surface imrm energy" << std::endl;)

  std::set<VH> bad_vhs;
  add_anti_normal_flip_term(svop, new_vh_idx, bad_vhs, 0.001);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added anti normal flip energy" << std::endl;)


  //hard constraints
  std::vector<VH> movable_vhs2, movable_vhs1, fixed_vhs;
  collect_constrained_vertices(new_vh_idx, bad_vhs, fixed_vhs, movable_vhs1, movable_vhs2);
  ALGOHEX_DEBUG_ONLY(std::cerr << "collected constrained vertices" << std::endl;)

  SMatD A;
  VecD b;
  setup_singular_vertex_constraints(3 * num, movable_vhs2, movable_vhs1, fixed_vhs, new_vh_idx, A, b);
  ALGOHEX_DEBUG_ONLY(std::cerr << "set constraints" << std::endl;)

//        std::vector<COMISO::NConstraintInterface *> constraint_pointers;
//        std::vector<COMISO::LinearConstraint> constraints;
//        setup_singular_vertex_constraints(3 * num, new_vh_idx, constraints, constraint_pointers);

//debug
  ALGOHEX_DEBUG_ONLY(std::cerr << "initial energy: ";
                             double initf = svop.initial_f(x.data());
                             std::cerr << initf << std::endl;)

  COMISO::NPTiming svopt(&svop);


//        COMISO::NPDerivativeChecker dc;
//        dc.check_all(&svop);



  // 2. optimize problem
  if (truncated_newton_)
  {
    COMISO::ConstraintTools::remove_dependent_linear_constraints(A, b);
    //in case LDLT factorization fails
    add_regularizer(svop, new_vh_idx, 1e-6);

    COMISO::TruncatedNewtonPCG tn(1e-4);
    tn.always_update_preconditioner() = true;
//  tn.adaptive_tolerance_modifier() = 0.01;
    tn.allow_warmstart() = true;
    tn.max_pcg_iters() = 300;
    tn.max_iters() = 100;
    tn.tolerance_gdx() = 1e-1;
    tn.solve_projected_normal_equation(&svopt, A, b);
  }
  else
  {
    COMISO::NewtonSolver ns(1e-6, 1e-9, 100);
    ns.solve_infeasible_start(&svopt, A, b);
  }

//        ns.solve(&svopt, A, b);

  //
//        COMISO::IPOPTSolver ipsol;
//
////        ipsol.app().Options()->SetStringValue("derivative_test", "first-order");
//
////        ipsol.app().Options()->SetIntegerValue("print_level", 12);
//
//        ipsol.app().Options()->SetIntegerValue("max_iter", _max_iter);
////      ipsol.app().Options()->SetNumericValue("tolerance", 0.00001);
//        ipsol.solve(&svopt, constraint_pointers);


  //test if ground truth energy decreases in each iter
//        auto coords = sgsp.get_all_coords();
//        coords.erase(coords.begin());
//        coords.erase(coords.begin());
//        _fvals = evaluate_development_of_groud_truth_energy(coords);
//        ALGOHEX_DEBUG_ONLY(std::cerr<<"final energy: ";
//                                   double ff = svop.initial_f(x.data());
//                                   std::cerr<<ff<<std::endl;)


  std::map<VH, Point> orig_ff_sg_vhs, orig_fe_sg_vhs;
//std::cerr<<"mov2 vhs ";
  for (auto vhi: movable_vhs2)
  {
    orig_ff_sg_vhs[vhi] = mesh_.vertex(vhi);
//            std::cerr<<" "<<vhi;
  }

  for (auto vhi: movable_vhs1)
    orig_fe_sg_vhs[vhi] = mesh_.vertex(vhi);
//std::cerr<<std::endl;
//        std::cerr<<"mov1 vhs ";
//        for(auto& vhi : movable_vhs1) {
//            std::cerr<<" "<<vhi;
//        }
//        std::cerr<<std::endl;
//        std::cerr<<"fix vhs ";
//        for(auto& vhi : fixed_vhs) {
//            std::cerr<<" "<<vhi;
//        }
//        std::cerr<<std::endl;

  //set new points
  for (auto &vi: new_vh_idx)
  {
    Point new_pt(0);
    for (int i = 0; i < 3; ++i)
      new_pt[i] = x[3 * vi.second + i];

    mesh_.set_vertex(vi.first, new_pt);
  }

  bool ck_vol = check_cells_volume(mesh_);
//        std::cerr<<"Projecting ";

  //project to original surface
  for (auto vhi: movable_vhs2)
  {
//            std::cerr<<" "<<vhi<<" line "<<mesh_.vertex(vhi)<<" "<<get_normal_on_reference_mesh(mesh_.vertex(vhi)) + mesh_.vertex(vhi)<<std::endl;
    project_vertex_to_original_surface(vhi, orig_ff_sg_vhs[vhi]);
  }


  //project to original feature arcs
  for (auto vhi: movable_vhs1)
  {
    project_vertex_to_original_feature_arc(vhi, orig_fe_sg_vhs[vhi]);
  }

//std::cerr<<std::endl;
  //if has degenerate cell, fall back
  //TODO: only check one-ring cells
  ALGOHEX_DEBUG_ONLY(std::cerr << "After projection";)
  if (ck_vol && !check_cells_volume(mesh_))
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Fall back to original points!" << std::endl;)
    for (auto &vi: new_vh_idx)
    {
      Point new_pt(0);
      for (int i = 0; i < 3; ++i)
        new_pt[i] = orig_x[3 * vi.first.idx() + i];

      mesh_.set_vertex(vi.first, new_pt);
    }
  }

  //align normal to the right axis of quaternions
  for (auto&[vhi, hfa]: v_aligned_axis)
  {
    //smooth quaternions
    std::vector<CH> cells;
    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
      cells.push_back(*vc_it);

    std::set<FH> bdy_fhs;
    for (const auto ch: cells)
      for (auto cf_it = mesh_.cf_iter(ch); cf_it.valid(); ++cf_it)
        if (mesh_.is_boundary(*cf_it))
          bdy_fhs.insert(*cf_it);

    std::map<CH, std::vector<std::pair<Vec3d, int>>> alignments;

    for (auto&[hfi, ax]: hfa)
    {
      CH ch = mesh_.incident_cell(hfi);

      auto nm = mesh_.normal(hfi);

      alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), ax);
    }

    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
    {
      for (auto ce_it = mesh_.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
      {
        if (feature_edge_[*ce_it] > 0)
        {
          auto dir = mesh_.vertex(mesh_.edge(*ce_it).to_vertex()) - mesh_.vertex(mesh_.edge(*ce_it).from_vertex());
          dir.normalize();

          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], dir).second;
          alignments[*vc_it].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
        }
      }
    }

    for (int i = 0; i < 10; ++i)
      QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, cells, bdy_fhs, alignments, tq_,
                                                   cell_quaternions_);
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
singular_vertices_relocation(const std::vector<EH> &_sg_ehs, const double _sge_w, const double _rgl_w,
                             const double _cs_w, const double _tps_w, const double _rp_w, const int _max_iter)
{
  auto v_aligned_axis = collect_normal_aligned_axis(_sg_ehs);

  //map variable vertices to new indices
  std::map<VH, int> new_vh_idx;
  int num = 0;
  for (auto ehi: _sg_ehs)
  {
    if (valence_[ehi] != 0 && valence_[ehi] != std::numeric_limits<int>::max())
    {
      auto vhf = mesh_.edge(ehi).from_vertex();
      auto vht = mesh_.edge(ehi).to_vertex();

      if (new_vh_idx.find(vhf) == new_vh_idx.end())
      {
        new_vh_idx[vhf] = num++;
      }

      if (new_vh_idx.find(vht) == new_vh_idx.end())
      {
        new_vh_idx[vht] = num++;
      }
    }
  }

  if (new_vh_idx.empty())
    return;

  //initialize variable values
  std::vector<double> x;
  x.resize(num * 3, 0);
  for (auto &vi: new_vh_idx)
  {
    for (int i = 0; i < 3; ++i)
      x[3 * vi.second + i] = mesh_.vertex(vi.first)[i];
  }


  //all original points
  std::vector<double> orig_x;
  orig_x.resize(mesh_.n_vertices() * 3, 0.);
  for (const auto &vh: mesh_.vertices())
  {
    const Point &pt = mesh_.vertex(vh);
    for (int i = 0; i < 3; ++i)
      orig_x[3 * vh.idx() + i] = pt[i];
  }

  SingularVertexOptProblemQuaternion svop(orig_x, x);


  //add alignment energy
  add_alignment_energy(svop, _sg_ehs, new_vh_idx, _sge_w);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added alignment energy" << std::endl;)

  //add shrink energy
  add_edge_shrink_energy(svop, _sg_ehs, new_vh_idx, _cs_w, _tps_w);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added shrink energy" << std::endl;)

  //add curvature term on singular graph
  add_curvature_smooth_energy(svop, _sg_ehs, new_vh_idx, _rgl_w);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added curvature_smooth energy" << std::endl;)

//        add_repulsion_energy(svop, _sg_ehs, new_vh_idx, _rp_w);
//        ALGOHEX_DEBUG_ONLY(std::cerr<<"added repulsion energy"<<std::endl;)

  add_imrm_energy2(svop, new_vh_idx, 0.001);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added imrm energy" << std::endl;)

  add_surface_imrm_energy(svop, _sg_ehs, new_vh_idx, 0.05);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added surface imrm energy" << std::endl;)

  std::set<VH> bad_vhs;
  add_anti_normal_flip_term(svop, _sg_ehs, new_vh_idx, bad_vhs, 0.001);
  ALGOHEX_DEBUG_ONLY(std::cerr << "added anti normal flip energy" << std::endl;)


  //hard constraints
  std::vector<VH> movable_vhs2, movable_vhs1, fixed_vhs;
  collect_constrained_vertices(new_vh_idx, bad_vhs, fixed_vhs, movable_vhs1, movable_vhs2);
  ALGOHEX_DEBUG_ONLY(std::cerr << "collected constrained vertices" << std::endl;)

  SMatD A;
  VecD b;
  setup_singular_vertex_constraints(3 * num, movable_vhs2, movable_vhs1, fixed_vhs, new_vh_idx, A, b);
  ALGOHEX_DEBUG_ONLY(std::cerr << "set constraints" << std::endl;)

//        std::vector<COMISO::NConstraintInterface *> constraint_pointers;
//        std::vector<COMISO::LinearConstraint> constraints;
//        setup_singular_vertex_constraints(3 * num, new_vh_idx, constraints, constraint_pointers);

//debug
  ALGOHEX_DEBUG_ONLY(std::cerr << "initial energy: ";
                             double initf = svop.initial_f(x.data());
                             std::cerr << initf << std::endl;)

  COMISO::NPTiming svopt(&svop);


//        COMISO::NPDerivativeChecker dc;
//        dc.check_all(&svop);



  // 2. optimize problem
  if (truncated_newton_)
  {
    COMISO::ConstraintTools::remove_dependent_linear_constraints(A, b);
    //in case LDLT factorization fails
    add_regularizer(svop, new_vh_idx, 1e-6);

    COMISO::TruncatedNewtonPCG tn(1e-4);
    tn.always_update_preconditioner() = true;
//  tn.adaptive_tolerance_modifier() = 0.01;
    tn.allow_warmstart() = true;
    tn.max_pcg_iters() = 300;
    tn.max_iters() = 100;
    tn.tolerance_gdx() = 1e-1;
    tn.solve_projected_normal_equation(&svopt, A, b);
  }
  else
  {
    COMISO::NewtonSolver ns(1e-6, 1e-9, 100);
    ns.solve_infeasible_start(&svopt, A, b);
  }

  std::map<VH, Point> orig_ff_sg_vhs, orig_fe_sg_vhs;

  for (auto vhi: movable_vhs2)
  {
    orig_ff_sg_vhs[vhi] = mesh_.vertex(vhi);
  }

  for (auto vhi: movable_vhs1)
    orig_fe_sg_vhs[vhi] = mesh_.vertex(vhi);

  //set new points
  for (auto &vi: new_vh_idx)
  {
    Point new_pt(0);
    for (int i = 0; i < 3; ++i)
      new_pt[i] = x[3 * vi.second + i];

    mesh_.set_vertex(vi.first, new_pt);
  }

  bool ck_vol = check_cells_volume(mesh_);
//        std::cerr<<"Projecting ";

  //project to original surface
  for (auto vhi: movable_vhs2)
  {
//            std::cerr<<" "<<vhi<<" line "<<mesh_.vertex(vhi)<<" "<<get_normal_on_reference_mesh(mesh_.vertex(vhi)) + mesh_.vertex(vhi)<<std::endl;
    project_vertex_to_original_surface(vhi, orig_ff_sg_vhs[vhi]);
  }


  //project to original feature arcs
  for (auto vhi: movable_vhs1)
  {
    project_vertex_to_original_feature_arc(vhi, orig_fe_sg_vhs[vhi]);
  }

//std::cerr<<std::endl;
  //if has degenerate cell, fall back
  //TODO: only check one-ring cells
  ALGOHEX_DEBUG_ONLY(std::cerr << "After projection";)
  if (ck_vol && !check_cells_volume(mesh_))
  {
    ALGOHEX_DEBUG_ONLY(std::cerr << "Fall back to original points!" << std::endl;)
    for (auto &vi: new_vh_idx)
    {
      Point new_pt(0);
      for (int i = 0; i < 3; ++i)
        new_pt[i] = orig_x[3 * vi.first.idx() + i];

      mesh_.set_vertex(vi.first, new_pt);
    }
  }

  //align normal to the right axis of quaternions
  for (auto&[vhi, hfa]: v_aligned_axis)
  {
    //smooth quaternions
    std::vector<CH> cells;
    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
      cells.push_back(*vc_it);

    std::set<FH> bdy_fhs;
    for (const auto ch: cells)
      for (auto cf_it = mesh_.cf_iter(ch); cf_it.valid(); ++cf_it)
        if (mesh_.is_boundary(*cf_it))
          bdy_fhs.insert(*cf_it);

    std::map<CH, std::vector<std::pair<Vec3d, int>>> alignments;

    for (auto&[hfi, ax]: hfa)
    {
      CH ch = mesh_.incident_cell(hfi);

      auto nm = mesh_.normal(hfi);

      alignments[ch].emplace_back(Vec3d(nm[0], nm[1], nm[2]), ax);
    }

    for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
    {
      for (auto ce_it = mesh_.ce_iter(*vc_it); ce_it.valid(); ++ce_it)
      {
        if (feature_edge_[*ce_it] > 0)
        {
          auto dir = mesh_.vertex(mesh_.edge(*ce_it).to_vertex()) - mesh_.vertex(mesh_.edge(*ce_it).from_vertex());
          dir.normalize();

          int axis = AxisAlignmentHelpers::closest_axis(cell_quaternions_[*vc_it], dir).second;
          alignments[*vc_it].emplace_back(Vec3d(dir[0], dir[1], dir[2]), axis);
        }
      }
    }

    for (int i = 0; i < 10; ++i)
      QuaternionSmoothing::smooth_field_quaternion(mesh_, trans_prop_, cells, bdy_fhs, alignments, tq_,
                                                   cell_quaternions_);
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::
add_alignment_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                     double _alignment_weight) const
{
  for (const auto &eh: mesh_.edges())
  {
    if (valence_[eh] != std::numeric_limits<int>::max() && valence_[eh] != std::numeric_limits<int>::lowest() &&
        valence_[eh] != 0)
    {
      Vec3d dir = best_aligned_rotational_axis(mesh_.halfedge_handle(eh, 0));

      VH vh0 = mesh_.edge(eh).from_vertex();
      VH vh1 = mesh_.edge(eh).to_vertex();

      _problem.add_perpendicular_element(_vh_to_idx[vh0], _vh_to_idx[vh1],
                                         dir, _alignment_weight / std::max(1e-10, mesh_.length(eh)));
    }
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
add_alignment_energy(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                     std::map<VH, int> &_vh_to_idx, double _alignment_weight) const
{
  for (auto eh: _sg_ehs)
  {
    if (valence_[eh] != std::numeric_limits<int>::max() && valence_[eh] != std::numeric_limits<int>::lowest() &&
        valence_[eh] != 0)
    {
      Vec3d dir = best_aligned_rotational_axis(mesh_.halfedge_handle(eh, 0));

      VH vh0 = mesh_.edge(eh).from_vertex();
      VH vh1 = mesh_.edge(eh).to_vertex();

      _problem.add_perpendicular_element(_vh_to_idx[vh0], _vh_to_idx[vh1],
                                         dir, _alignment_weight / std::max(1e-10, mesh_.length(eh)));
    }
  }
}


template<class MeshT>
Vec3d MeshOptimizationT<MeshT>::
best_aligned_rotational_axis(const HEH _heh) const
{
  Vec3d dir_avr(0, 0, 0);

  auto eh = mesh_.edge_handle(_heh);

  if (valence_[eh] % 4 != 0)
  {
    auto hehf_it = mesh_.hehf_iter(_heh);
    CH ch_s(-1);
    bool e_bdy = mesh_.is_boundary(_heh);
    if (!e_bdy)
      ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
    else
      ch_s = mesh_.incident_cell(*hehf_it);

    int axis = EdgeMonodromyHelperT::halfedge_rotational_axis_idx(mesh_, tq_, cell_quaternions_, trans_prop_, valence_,
                                                                  _heh, *hehf_it);
    if (axis == -1)
    {
      std::cout << "Error: singular edge: " << eh << "of valence " << valence_[eh] << " has -1 axis"
                << std::endl;
      return dir_avr;
    }

    dir_avr += AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[ch_s], (AxisAlignment) axis);

    if (e_bdy)
      ++hehf_it;

    for (; hehf_it.valid(); ++hehf_it)
    {
      auto chi = mesh_.incident_cell(*hehf_it);
      if (chi.is_valid())
      {
        if (chi == ch_s)
          break;

        axis = tq_.axis_after_transition(axis, trans_prop_[*hehf_it]);
        dir_avr += AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[chi], (AxisAlignment) axis);
      }
    }
  }
  else if (valence_[eh] % 4 == 0)
  {//TODO: does not consider boundary edge here
    auto hehf_it = mesh_.hehf_iter(_heh);

    CH ch_s = mesh_.incident_cell(mesh_.opposite_halfface_handle(*hehf_it));
    int dm_ax = EdgeMonodromyHelperT::get_dominant_axis_in_cell(mesh_, tq_, cell_quaternions_, trans_prop_, _heh);

    for (; hehf_it.valid(); ++hehf_it)
    {
      CH chi = mesh_.incident_cell(*hehf_it);
      dm_ax = tq_.axis_after_transition(dm_ax, trans_prop_[*hehf_it]);
      dir_avr += AxisAlignmentHelpers::quaternion_vector(cell_quaternions_[chi], (AxisAlignment) dm_ax);
    }
  }


  dir_avr.normalize();

  return dir_avr;
}


template<class MeshT>
void MeshOptimizationT<MeshT>::
add_edge_shrink_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx, double _cs_weight,
                       double _tps_w) const
{
  for (const auto &eh: mesh_.edges())
  {
    if (!feature_face_edge_[eh] && valence_[eh] == std::numeric_limits<int>::max())
    {
      auto vh0 = mesh_.edge(eh).from_vertex();
      auto vh1 = mesh_.edge(eh).to_vertex();

      _problem.add_shrink_edge_element(_vh_to_idx[vh0], _vh_to_idx[vh1], _cs_weight / mesh_.length(eh));
    }
  }


  //shrink zipper node
  for (const auto &vh: mesh_.vertices())
  {
    if (!feature_face_vertex_[vh])
    {
      int n1 = 0, n_1 = 0;
      for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
      {
        if (valence_[*ve_it] == 1)
          n1++;
        if (valence_[*ve_it] == -1)
          n_1++;
      }

      if (n1 == 1 && n_1 == 1)
      {
        for (auto voh_it = mesh_.voh_iter(vh); voh_it.valid(); ++voh_it)
        {
          if (valence_[mesh_.edge_handle(*voh_it)] != 0)
          {
            auto vh0 = mesh_.halfedge(*voh_it).from_vertex();
            auto vh1 = mesh_.halfedge(*voh_it).to_vertex();

            _problem.add_shrink_edge_element(_vh_to_idx[vh0], _vh_to_idx[vh1],
                                             _tps_w / std::max(1e-10, mesh_.length(*voh_it)));
          }
        }
      }
    }
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
add_edge_shrink_energy(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                       std::map<VH, int> &_vh_to_idx, double _cs_weight, double _tps_w) const
{
  //shrink zipper node
  for (auto &[vh, idx]: _vh_to_idx)
  {
    if (!feature_face_vertex_[vh])
    {
      int n1 = 0, n_1 = 0;
      for (auto ve_it = mesh_.ve_iter(vh); ve_it.valid(); ++ve_it)
      {
        if (valence_[*ve_it] == 1)
          n1++;
        if (valence_[*ve_it] == -1)
          n_1++;
      }

      if (n1 == 1 && n_1 == 1)
      {
        for (auto voh_it = mesh_.voh_iter(vh); voh_it.valid(); ++voh_it)
        {
          if (valence_[mesh_.edge_handle(*voh_it)] != 0)
          {
            auto vh1 = mesh_.halfedge(*voh_it).to_vertex();
            if (_vh_to_idx.find(vh1) != _vh_to_idx.end())
              _problem.add_shrink_edge_element(idx, _vh_to_idx[vh1],
                                               _tps_w / std::max(1e-10, mesh_.length(*voh_it)));
          }
        }
      }
    }
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::
add_curvature_smooth_energy(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                            std::map<VH, int> &_vh_to_idx, double _rgl_weight) const
{
  //add curvature term on singular graph
  for (auto &[vh, idx]: _vh_to_idx)
  {
    std::vector<VH> vhs;
    double len = 0.;
    for (auto voh_it = mesh_.voh_iter(vh); voh_it.valid(); ++voh_it)
    {
      if (valence_[mesh_.edge_handle(*voh_it)] != 0)
      {
        auto vht = mesh_.halfedge(*voh_it).to_vertex();
        if (_vh_to_idx.find(vht) != _vh_to_idx.end())
        {
          vhs.push_back(vht);
          len += mesh_.length(*voh_it);
        }
      }
    }

    if (vhs.size() == 2)
    {
      _problem.add_curvature_smooth_element(_vh_to_idx[vhs[0]], idx, _vh_to_idx[vhs[1]],
                                            _rgl_weight / std::max(1e-10, len));
    }
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
add_curvature_smooth_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                            double _rgl_weight) const
{
  //add curvature term on singular graph
  for (const auto &vh: mesh_.vertices())
  {
    if (sgl_vt_[vh] > 0)
    {
      std::vector<VH> vhs;
      double len = 0.;
      for (auto voh_it = mesh_.voh_iter(vh); voh_it.valid(); ++voh_it)
      {
        if (valence_[mesh_.edge_handle(*voh_it)] != 0)
        {
          vhs.push_back(mesh_.halfedge(*voh_it).to_vertex());
          len += mesh_.length(*voh_it);
        }
      }

      if (vhs.size() == 2)
      {
        _problem.add_curvature_smooth_element(_vh_to_idx[vhs[0]], _vh_to_idx[vh], _vh_to_idx[vhs[1]],
                                              _rgl_weight / std::max(1e-10, len));
      }
    }
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::
add_imrm_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx, double _imrm_weight) const
{
  //add IMRM energy
  double w = 10.;
  std::set<CH> chs1, chs2, chs3, chs4;
  for (auto &vi: _vh_to_idx)
  {
    for (auto vc_it = mesh_.vc_iter(vi.first); vc_it.valid(); ++vc_it)
    {
      int n_sgv = n_singular_vertex_in_cell(*vc_it);
      if (n_sgv == 1)
        chs1.insert(*vc_it);
      else if (n_sgv == 2)
        chs2.insert(*vc_it);
      else if (n_sgv == 3)
        chs3.insert(*vc_it);
      else if (n_sgv == 4)
        chs4.insert(*vc_it);
    }
  }
  for (const auto &ch: chs1)
  {
    double weight = _imrm_weight;
    if (imrm_cell_[ch])
      weight *= w;
    auto cvhs = get_cell_vertices(ch);
//            ALGOHEX_DEBUG_ONLY(std::cerr<<"ch1: "<<_vh_to_idx[cvhs[0]]<<" "<<cvhs[1]<<" "<<cvhs[2]<<" "<<cvhs[3]<<std::endl;

    _problem.add_imrm_energy_element_one(_vh_to_idx[cvhs[0]], cvhs[1].idx(), cvhs[2].idx(), cvhs[3].idx(), weight);
  }

  for (const auto &ch: chs2)
  {
    double weight = _imrm_weight;
    if (imrm_cell_[ch])
      weight *= w;

    auto cvhs = get_cell_vertices(ch);
//            std::cerr<<"ch2: "<<_vh_to_idx[cvhs[0]]<<" "<<_vh_to_idx[cvhs[1]]<<" "<<cvhs[2]<<" "<<cvhs[3]<<std::endl;
    _problem.add_imrm_energy_element_two(_vh_to_idx[cvhs[0]], _vh_to_idx[cvhs[1]], cvhs[2].idx(), cvhs[3].idx(),
                                         weight);
  }

  for (const auto &ch: chs3)
  {
    double weight = _imrm_weight;
    if (imrm_cell_[ch])
      weight *= w;

    auto cvhs = get_cell_vertices(ch);
    _problem.add_imrm_energy_element_three(_vh_to_idx[cvhs[0]], _vh_to_idx[cvhs[1]], _vh_to_idx[cvhs[2]], cvhs[3].idx(),
                                           weight);
  }

  for (const auto &ch: chs4)
  {
    double weight = _imrm_weight;
    if (imrm_cell_[ch])
      weight *= w;

    auto cvhs = get_cell_vertices(ch);
    _problem.add_imrm_energy_element_four(_vh_to_idx[cvhs[0]], _vh_to_idx[cvhs[1]], _vh_to_idx[cvhs[2]],
                                          _vh_to_idx[cvhs[3]], weight);
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
add_imrm_energy2(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx, double _imrm_weight) const
{
  //add IMRM energy
  double w = 10.;
  std::set<CH> chs1, chs2, chs3, chs4;
  for (auto &vi: _vh_to_idx)
  {
    for (auto vc_it = mesh_.vc_iter(vi.first); vc_it.valid(); ++vc_it)
    {
      int n_sgv = n_singular_vertex_in_cell2(_vh_to_idx, *vc_it);
      if (n_sgv == 1)
        chs1.insert(*vc_it);
      else if (n_sgv == 2)
        chs2.insert(*vc_it);
      else if (n_sgv == 3)
        chs3.insert(*vc_it);
      else if (n_sgv == 4)
        chs4.insert(*vc_it);
    }
  }
  for (const auto &ch: chs1)
  {
    double weight = _imrm_weight;
    if (imrm_cell_[ch])
      weight *= w;
    auto cvhs = get_cell_vertices2(_vh_to_idx, ch);
//            ALGOHEX_DEBUG_ONLY(std::cerr<<"ch1: "<<_vh_to_idx[cvhs[0]]<<" "<<cvhs[1]<<" "<<cvhs[2]<<" "<<cvhs[3]<<std::endl;

    _problem.add_imrm_energy_element_one(_vh_to_idx[cvhs[0]], cvhs[1].idx(), cvhs[2].idx(), cvhs[3].idx(), weight);
  }

  for (const auto &ch: chs2)
  {
    double weight = _imrm_weight;
    if (imrm_cell_[ch])
      weight *= w;

    auto cvhs = get_cell_vertices2(_vh_to_idx, ch);
//            std::cerr<<"ch2: "<<_vh_to_idx[cvhs[0]]<<" "<<_vh_to_idx[cvhs[1]]<<" "<<cvhs[2]<<" "<<cvhs[3]<<std::endl;
    _problem.add_imrm_energy_element_two(_vh_to_idx[cvhs[0]], _vh_to_idx[cvhs[1]], cvhs[2].idx(), cvhs[3].idx(),
                                         weight);
  }

  for (const auto &ch: chs3)
  {
    double weight = _imrm_weight;
    if (imrm_cell_[ch])
      weight *= w;

    auto cvhs = get_cell_vertices2(_vh_to_idx, ch);
    _problem.add_imrm_energy_element_three(_vh_to_idx[cvhs[0]], _vh_to_idx[cvhs[1]], _vh_to_idx[cvhs[2]], cvhs[3].idx(),
                                           weight);
  }

  for (const auto &ch: chs4)
  {
    double weight = _imrm_weight;
    if (imrm_cell_[ch])
      weight *= w;

    auto cvhs = get_cell_vertices2(_vh_to_idx, ch);
    _problem.add_imrm_energy_element_four(_vh_to_idx[cvhs[0]], _vh_to_idx[cvhs[1]], _vh_to_idx[cvhs[2]],
                                          _vh_to_idx[cvhs[3]], weight);
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
add_surface_imrm_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                        double _sf_imrm_weight) const
{
  for (const auto &vh: mesh_.vertices())
  {
//            if (sgl_vt_[vh] == 1 && mesh_.is_boundary(vh)) {
    if (sgl_vt_[vh] > 0 && feature_face_vertex_[vh])
    {
      if (feature_edge_vertex_[vh] > 0)
        continue;

      Point normal = get_normal_on_reference_mesh(mesh_.vertex(vh));

      //get incident halffaces
      std::set<HFH> ft_hfhs;
      for (auto vf_it = mesh_.vf_iter(vh); vf_it.valid(); ++vf_it)
      {
        if (feature_fprop_[*vf_it] > 0)
        {
          ft_hfhs.insert(mesh_.halfface_handle(*vf_it, 0));
        }
      }

      for (auto hfi: ft_hfhs)
      {
        auto hfvhs = mesh_.get_halfface_vertices(hfi, vh);
        _problem.add_boundary_simrm_element(_vh_to_idx[hfvhs[0]], hfvhs[1].idx(), hfvhs[2].idx(),
                                            Eigen::Vector3d(normal[0], normal[1], normal[2]),
                                            _sf_imrm_weight);
      }
    }
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
add_surface_imrm_energy(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                        std::map<VH, int> &_vh_to_idx, double _sf_imrm_weight) const
{
  for (auto &[vh, idx]: _vh_to_idx)
  {
    if (sgl_vt_[vh] > 0 && feature_face_vertex_[vh])
    {
      if (feature_edge_vertex_[vh] > 0)
        continue;

      Point normal = get_normal_on_reference_mesh(mesh_.vertex(vh));

      //get incident halffaces
      std::set<HFH> ft_hfhs;
      for (auto vf_it = mesh_.vf_iter(vh); vf_it.valid(); ++vf_it)
      {
        if (feature_fprop_[*vf_it] > 0)
        {
          ft_hfhs.insert(mesh_.halfface_handle(*vf_it, 0));
        }
      }

      for (auto hfi: ft_hfhs)
      {
        auto hfvhs = mesh_.get_halfface_vertices(hfi, vh);
        _problem.add_boundary_simrm_element(_vh_to_idx[hfvhs[0]], hfvhs[1].idx(), hfvhs[2].idx(),
                                            Eigen::Vector3d(normal[0], normal[1], normal[2]),
                                            _sf_imrm_weight);
      }
    }
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
add_anti_normal_flip_term(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                          std::set<VH> &_bad_vhs, double _normal_weight) const
{
  std::vector<Point> cell_pts(4);
  ALGOHEX_DEBUG_ONLY(std::cerr << "virtual incident cells are degenerate at vertices: ";)
  for (const auto &vh: mesh_.vertices())
  {
    if (sgl_vt_[vh] != 0
        && !feature_edge_vertex_[vh]//for normal
        && mesh_.is_boundary(vh))
    {
      Point pt = mesh_.vertex(vh);
      auto normal = get_normal_on_reference_mesh(pt);
      auto top_pt = pt + normal * target_length_[vh];

      for (auto vhf_it = mesh_.vhf_iter(vh); vhf_it.valid(); ++vhf_it)
      {
        if (mesh_.is_boundary(*vhf_it))
        {
          auto hfvhs = mesh_.get_halfface_vertices(*vhf_it, vh);

          for (int i = 0; i < 3; ++i)
            cell_pts[i] = mesh_.vertex(hfvhs[i]);
          cell_pts[3] = top_pt;

          if (get_cell_volume(cell_pts) < 1e-16)
          {
            _bad_vhs.insert(vh);
            ALGOHEX_DEBUG_ONLY(std::cerr << " " << vh;)
            break;
          }
          else
            _problem.add_anti_normal_flip_element(_vh_to_idx[hfvhs[0]], hfvhs[1].idx(), hfvhs[2].idx(),
                                                  Eigen::Vector3d(top_pt[0], top_pt[1], top_pt[2]), _normal_weight);
        }
      }
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
}

template<class MeshT>
void MeshOptimizationT<MeshT>::
add_anti_normal_flip_term(SingularVertexOptProblemQuaternion &_problem, const std::vector<EH> &_sg_ehs,
                          std::map<VH, int> &_vh_to_idx, std::set<VH> &_bad_vhs, double _normal_weight) const
{
  std::vector<Point> cell_pts(4);
  ALGOHEX_DEBUG_ONLY(std::cerr << "virtual incident cells are degenerate at vertices: ";)
  for (auto &[vh, idx]: _vh_to_idx)
  {
    if (sgl_vt_[vh] == 1
        && !feature_edge_vertex_[vh]//for normal
        && mesh_.is_boundary(vh))
    {
      Point pt = mesh_.vertex(vh);
      auto normal = get_normal_on_reference_mesh(pt);
      auto top_pt = pt + normal * target_length_[vh];

      for (auto vhf_it = mesh_.vhf_iter(vh); vhf_it.valid(); ++vhf_it)
      {
        if (mesh_.is_boundary(*vhf_it))
        {
          auto hfvhs = mesh_.get_halfface_vertices(*vhf_it, vh);

          for (int i = 0; i < 3; ++i)
            cell_pts[i] = mesh_.vertex(hfvhs[i]);
          cell_pts[3] = top_pt;

          if (get_cell_volume(cell_pts) < 1e-16)
          {
            _bad_vhs.insert(vh);
            ALGOHEX_DEBUG_ONLY(std::cerr << " " << vh;)
            break;
          }
          else
            _problem.add_anti_normal_flip_element(_vh_to_idx[hfvhs[0]], hfvhs[1].idx(), hfvhs[2].idx(),
                                                  Eigen::Vector3d(top_pt[0], top_pt[1], top_pt[2]), _normal_weight);
        }
      }
    }
  }
  ALGOHEX_DEBUG_ONLY(std::cerr << std::endl;)
}


template<class MeshT>
void MeshOptimizationT<MeshT>::
add_regularizer(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx, const double _rgl_weight)
{
  for (auto&[vhi, idx]: _vh_to_idx)
  {
    auto pt = mesh_.vertex(vhi);
    _problem.add_repulsion_term(idx, Eigen::Vector3d(pt[0], pt[1], pt[2]), _rgl_weight);
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::
add_repulsion_energy(SingularVertexOptProblemQuaternion &_problem, std::map<VH, int> &_vh_to_idx,
                     const double _rp_weight)
{
  //add repulsion
  int max_pair_id = *std::max_element(sg_edge_pairs_.begin(), sg_edge_pairs_.end());

  std::vector<std::set<VH>> cd_vhs;
  std::vector<std::set<EH>> cd_ehs;
  std::set<EH> rel_ehs;

  cd_vhs.reserve(max_pair_id);
  cd_ehs.reserve(max_pair_id);

  //need to update edge pair property because of the remeshing
  for (int i = 1; i <= max_pair_id; ++i)
  {
    //find an edge of the id
    EH eh_x(-1);
    for (const auto &ehi: mesh_.edges())
    {
      if (sg_edge_pairs_[ehi] == i && valence_[ehi] >= -2 && valence_[ehi] <= 4)
      {
        eh_x = ehi;
        break;
      }
    }
    if (!eh_x.is_valid())
      continue;

    if (valence_[eh_x] == 0)
    {
      sg_edge_pairs_[eh_x] = 0;
      continue;
    }

    //update
    auto sg_ehs = get_edges_on_singular_arc(mesh_, valence_, mesh_.halfedge_handle(eh_x, 0));
    for (const auto &ehi: sg_ehs)
      sg_edge_pairs_[ehi] = i;

    std::set<VH> vhs_set;
    std::set<EH> ehs_set;
    for (const auto &ehi: sg_ehs)
    {
      auto vh0 = mesh_.edge(ehi).from_vertex();
      auto vh1 = mesh_.edge(ehi).to_vertex();

      vhs_set.insert(vh0);
      vhs_set.insert(vh1);

      ehs_set.insert(ehi);
      rel_ehs.insert(ehi);
    }

    if (!vhs_set.empty())
    {
      cd_vhs.push_back(vhs_set);
      cd_ehs.push_back(ehs_set);
    }
  }

  //cells with higher weight for better tet quality
  for (const auto chi: mesh_.cells())
    imrm_cell_[chi] = false;

  for (auto i = 0u; i < cd_vhs.size(); ++i)
  {
    double max_target_dist = -1, min_target_dist = 100000;
    for (const auto &vhi: cd_vhs[i])
    {
      if (min_target_dist > target_length_[vhi])
        min_target_dist = target_length_[vhi];
      if (max_target_dist < target_length_[vhi])
        max_target_dist = target_length_[vhi];
    }
    ALGOHEX_DEBUG_ONLY(
            std::cerr << "Max target dist: " << max_target_dist << " min dist: " << min_target_dist << " size "
                      << cd_vhs[i].size() << std::endl;)


    for (const auto &vhi: cd_vhs[i])
    {
      //the vertex is not a variable vertex
      if (_vh_to_idx.find(vhi) == _vh_to_idx.end())
        continue;

      if (n_incident_singular_edges(vhi) > 2)
        continue;

      std::set<EH> kr_ehs;
      std::set<VH> kr_vhs{vhi};
      get_k_ring_vertices(mesh_, 2, kr_vhs);
      auto kr_chs = get_k_ring_cells(mesh_, kr_vhs);
      for (const auto chi: kr_chs)
      {
        for (auto ce_it = mesh_.ce_iter(chi); ce_it.valid(); ++ce_it)
        {
          if (valence_[*ce_it] != 0 && rel_ehs.find(*ce_it) != rel_ehs.end())
            kr_ehs.insert(*ce_it);
        }
      }

      auto cls_pt = find_closest_point_on_other_singular_arcs(vhi, cd_ehs[i], kr_ehs);
      if (cls_pt.first < min_target_dist)
      {
        for (auto vc_it = mesh_.vc_iter(vhi); vc_it.valid(); ++vc_it)
          imrm_cell_[*vc_it] = true;

        auto dir = mesh_.vertex(vhi) - cls_pt.second;
        if (dir.norm() == 0)
        {
          std::cerr << "Warning: find the closest point that has the same position. Query vertex " << vhi << std::endl;
          continue;
        }

        dir.normalize();
        auto mid_pt = cls_pt.second + mesh_.vertex(vhi);
        mid_pt /= 2.;
//                  std::cerr<<"vh: "<<vhi<<" move dir: "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<std::endl;
        auto target_pt = mid_pt + min_target_dist * dir;
        //old
        _problem.add_repulsion_term(_vh_to_idx[vhi],
                                    Eigen::Vector3d(target_pt[0], target_pt[1], target_pt[2]),
                                    _rp_weight / std::max(1e-8, cls_pt.first));
      }
    }

    //straigten the arc
    //add curvature term on singular graph
    for (const auto &vhi: cd_vhs[i])
    {
      if (!mesh_.is_boundary(vhi) && _vh_to_idx.find(vhi) != _vh_to_idx.end())
      {
        std::vector<VH> vhs;
        double len = 0.;
        for (auto voh_it = mesh_.voh_iter(vhi); voh_it.valid(); ++voh_it)
        {
          if (valence_[mesh_.edge_handle(*voh_it)] != 0)
          {
            vhs.push_back(mesh_.halfedge(*voh_it).to_vertex());
            len += mesh_.length(*voh_it);
          }
        }

        if (vhs.size() == 2)
        {
          if (_vh_to_idx.find(vhs[0]) != _vh_to_idx.end() && _vh_to_idx.find(vhs[1]) != _vh_to_idx.end())
            _problem.add_curvature_smooth_element(_vh_to_idx[vhs[0]], _vh_to_idx[vhi], _vh_to_idx[vhs[1]],
                                                  10. * _rp_weight / std::max(1e-8, len));
        }
      }
    }
  }
}


template<class MeshT>
void
MeshOptimizationT<MeshT>::collect_constrained_vertices(std::map<VH, int> &_new_vh_idx, const std::set<VH> &_bad_vhs,
                                                       std::vector<VH> &_fixed_vhs, std::vector<VH> &_movable_vhs1,
                                                       std::vector<VH> &_movable_vhs2)
{
  for (auto &vi: _new_vh_idx)
  {
    if (sgl_vt_[vi.first] > 0)
    {
      if (feature_node_[vi.first] > 0)
        _fixed_vhs.push_back(vi.first);
      else if (feature_edge_vertex_[vi.first])
      {
        bool has_nonfe_sge = false;
        for (auto ve_it = mesh_.ve_iter(vi.first); ve_it.valid(); ++ve_it)
          if (valence_[*ve_it] != 0 && feature_edge_[*ve_it] == 0)
          {
            has_nonfe_sge = true;
            break;
          }

        if (has_nonfe_sge)
        {
          if (_bad_vhs.find(vi.first) == _bad_vhs.end())
            _movable_vhs1.push_back(vi.first);
          else
            _fixed_vhs.push_back(vi.first);
        }
        else
          _fixed_vhs.push_back(vi.first);
      }
      else if (feature_face_vertex_[vi.first])
      {
        if (_bad_vhs.find(vi.first) == _bad_vhs.end())
          _movable_vhs2.push_back(vi.first);
        else
          _fixed_vhs.push_back(vi.first);
      }
    }
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::collect_fixed_vertices(std::map<VH, int> &_new_vh_idx, const std::set<VH> &_bad_vhs,
                                                      std::vector<VH> &_fixed_vhs)
{
  for (auto &vi: _new_vh_idx)
  {
    if (sgl_vt_[vi.first] == 2
        || feature_edge_vertex_[vi.first]
        || (sgl_vt_[vi.first] == 1 && mesh_.is_boundary(vi.first) && n_incident_singular_edges(vi.first) > 1)
        || _bad_vhs.find(vi.first) != _bad_vhs.end())
      _fixed_vhs.push_back(vi.first);
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::collect_fixed_vertices(std::map<VH, int> &_new_vh_idx, std::vector<VH> &_fixed_vhs)
{
  for (auto &vi: _new_vh_idx)
  {
    if (mesh_.is_boundary(vi.first))
      _fixed_vhs.push_back(vi.first);
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::collect_fixed_vertices_post(std::map<VH, int> &_new_vh_idx, std::vector<VH> &_fixed_vhs)
{
  for (auto &vi: _new_vh_idx)
  {
    int n_sgvhs = n_incident_singular_edges(vi.first);
    if ((n_sgvhs > 0 && n_sgvhs != 2) || MeshPropertiesT<MeshT>::node_index(vi.first) == 10)
      _fixed_vhs.push_back(vi.first);
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::setup_singular_vertex_constraints(int _n_unknowns, std::map<VH, int> &_new_vh_idx,
                                                                 std::vector<COMISO::LinearConstraint> &_constraints,
                                                                 std::vector<COMISO::NConstraintInterface *> &_constraint_pointers)
{
  COMISO::LinearConstraint::SVectorNC coeffs(_n_unknowns);

  //boundary vertices that can move in the tangent plane
  //boundary vertices
  int n_cons_free = 0;
  for (auto &vi: _new_vh_idx)
    if (sgl_vt_[vi.first] == 1 && mesh_.is_boundary(vi.first) && n_incident_singular_edges(vi.first) == 1)
      n_cons_free++;

  int n_cons_fixed = 0;
  for (auto &vi: _new_vh_idx)
    if (sgl_vt_[vi.first] == 2
        || feature_edge_vertex_[vi.first]
        || (sgl_vt_[vi.first] == 1 && mesh_.is_boundary(vi.first) && n_incident_singular_edges(vi.first) > 1))
      n_cons_fixed++;

  _constraints.reserve(n_cons_free + 3 * n_cons_fixed);
  _constraint_pointers.reserve(_constraints.size());

  //free vertices
  for (auto &vi: _new_vh_idx)
    if (sgl_vt_[vi.first] == 1 && mesh_.is_boundary(vi.first) && n_incident_singular_edges(vi.first) == 1)
    {
      auto pt = mesh_.vertex(vi.first);

      auto n = get_normal_on_reference_mesh(pt);
//                auto nn = vertex_normal_area_weighted(mesh_, vi.first);
//                std::cout<<"vertex: "<<vi.first<<" normal: "<<n<<" vv: "<<nn<<std::endl;

      coeffs.setZero();

      for (int i = 0; i < 3; ++i)
        coeffs.coeffRef(3 * vi.second + i) += n[i];

      _constraints.push_back(
              COMISO::LinearConstraint(coeffs, -n[0] * pt[0] - n[1] * pt[1] - n[2] * pt[2],
                                       COMISO::NConstraintInterface::NC_EQUAL));
    }


  //fixed vertices on boundary singular edges
  for (auto &vi: _new_vh_idx)
    if (sgl_vt_[vi.first] == 2
        || feature_edge_vertex_[vi.first]
        || (sgl_vt_[vi.first] == 1 && mesh_.is_boundary(vi.first) && n_incident_singular_edges(vi.first) > 1))
    {
      auto pt = mesh_.vertex(vi.first);

      for (int i = 0; i < 3; ++i)
      {
        coeffs.setZero();

        coeffs.coeffRef(3 * vi.second + i) += 1.;
        _constraints.push_back(COMISO::LinearConstraint(coeffs, -pt[i], COMISO::NConstraintInterface::NC_EQUAL));
      }
    }


  for (auto &cons: _constraints)
    _constraint_pointers.push_back(&cons);
}


template<class MeshT>
void MeshOptimizationT<MeshT>::setup_singular_vertex_constraints(int _n_unknowns, const std::vector<VH> &_movable_vhs2,
                                                                 const std::vector<VH> &_movable_vhs1,
                                                                 const std::vector<VH> &_fixed_vhs,
                                                                 std::map<VH, int> &_new_vh_idx, SMatD &_A, VecD &_b)
{
  //check
  for (const auto vhi: _movable_vhs2)
    if (_new_vh_idx.find(vhi) == _new_vh_idx.end())
      std::cerr << "Error: the vertex2 " << vhi << " is not in the map" << std::endl;
  for (const auto vhi: _movable_vhs1)
    if (_new_vh_idx.find(vhi) == _new_vh_idx.end())
      std::cerr << "Error: the vertex1 " << vhi << " is not in the map" << std::endl;
  for (const auto vhi: _fixed_vhs)
    if (_new_vh_idx.find(vhi) == _new_vh_idx.end())
      std::cerr << "Error: the vertex0 " << vhi << " is not in the map" << std::endl;

  //vertices on features that can are constrained
  int n_cons_2d = (int) _movable_vhs2.size();
  int n_cons_1d = (int) _movable_vhs1.size();
  int n_cons_fixed = (int) _fixed_vhs.size();


  std::vector<T> v_triplets;
  v_triplets.reserve(3 * n_cons_2d + 6 * n_cons_1d + 3 * n_cons_fixed);

  //cons_1d is not linearly independent
  int n_cons = n_cons_2d + 3 * n_cons_1d + 3 * n_cons_fixed;
  _A.resize(n_cons, _n_unknowns);
  _b.resize(n_cons);

  //row index of A
  int num = 0;
  //constrain vertices that can move on feature surface
  for (const auto vhi: _movable_vhs2)
  {
    auto pt = mesh_.vertex(vhi);
    auto n = get_normal_on_reference_mesh(pt);

    for (int i = 0; i < 3; ++i)
      v_triplets.emplace_back(num, 3 * _new_vh_idx[vhi] + i, n[i]);

    _b[num++] = n[0] * pt[0] + n[1] * pt[1] + n[2] * pt[2];
  }

  //constrain vertices that can move on feature edge
  for (const auto vhi: _movable_vhs1)
  {
    auto pt = mesh_.vertex(vhi);

    //constrained direction
    Point d(0.);
    std::vector<Point> fe_vs;
    for (auto vv_it = mesh_.vv_iter(vhi); vv_it.valid(); ++vv_it)
      if (feature_edge_vertex_[*vv_it])
        fe_vs.push_back(mesh_.vertex(*vv_it));

    if (fe_vs.size() != 2)
      std::cerr << "Error: not a feature edge vertex!" << std::endl;

    d = (fe_vs[0] - pt).normalize() + (pt - fe_vs[1]).normalize();

    //add constraints
    v_triplets.emplace_back(num, 3 * _new_vh_idx[vhi], d[1]);
    v_triplets.emplace_back(num, 3 * _new_vh_idx[vhi] + 1, -d[0]);
    _b[num++] = d[1] * pt[0] - d[0] * pt[1];

    v_triplets.emplace_back(num, 3 * _new_vh_idx[vhi] + 1, d[2]);
    v_triplets.emplace_back(num, 3 * _new_vh_idx[vhi] + 2, -d[1]);
    _b[num++] = d[2] * pt[1] - d[1] * pt[2];

    v_triplets.emplace_back(num, 3 * _new_vh_idx[vhi] + 2, d[0]);
    v_triplets.emplace_back(num, 3 * _new_vh_idx[vhi], -d[2]);
    _b[num++] = d[0] * pt[2] - d[2] * pt[0];
  }

  //constrain fixed vertices
  for (const auto vhi: _fixed_vhs)
  {
    auto pt = mesh_.vertex(vhi);
    for (int i = 0; i < 3; ++i)
    {
      v_triplets.emplace_back(num, 3 * _new_vh_idx[vhi] + i, 1.);
      _b[num++] = pt[i];
    }
  }

  //set up matrix
  _A.setFromTriplets(v_triplets.begin(), v_triplets.end());
}


template<class MeshT>
typename MeshOptimizationT<MeshT>::Point
MeshOptimizationT<MeshT>::get_normal_on_reference_mesh(const Point &_pt) const
{
  auto closest = bsp_->nearest(_pt);
  auto fh = closest.handle;

  auto fv_it = sf_mesh_.fv_iter(fh);
  const Point &pt0 = sf_mesh_.vertex(*fv_it);
  const Point &n0 = vertex_normal_[*fv_it];

  const Point &pt1 = sf_mesh_.vertex(*(++fv_it));
  const Point &n1 = vertex_normal_[*fv_it];

  const Point &pt2 = sf_mesh_.vertex(*(++fv_it));
  const Point &n2 = vertex_normal_[*fv_it];

  if (sharp_face_[fh])
  {
    if (_pt == pt0)
    {
//                std::cerr<<"projection is the same as the face vertex"<<std::endl;
      return n0;
    }
    else if (_pt == pt1)
    {
//                std::cerr<<"projection is the same as the face vertex"<<std::endl;

      return n1;
    }
    else if (_pt == pt2)
    {
//                std::cerr<<"projection is the same as the face vertex"<<std::endl;

      return n2;
    }
    else
      return face_normal_[fh];
  }
  else
  {
    Point p_best;
    double dist = ACG::Geometry::distPointTriangleSquaredStable(_pt, pt0, pt1, pt2, p_best);

    auto bary_cds = compute_bary_coord(p_best, pt0, pt1, pt2);

    auto n = bary_cds[0] * n0 + bary_cds[1] * n1 + bary_cds[2] * n2;
    n.normalize();

    return Point(n[0], n[1], n[2]);
  }
}


template<class MeshT>
std::pair<double, typename MeshT::PointT> MeshOptimizationT<MeshT>::
find_closest_point_on_other_singular_arcs(const VH _vh, const std::set<EH> &_same_arc_ehs,
                                          const std::set<EH> &_all_ehs) const
{
  double min_dist = std::numeric_limits<double>::max();
  Point min_pt(0.);
  const Point &pt = mesh_.vertex(_vh);

  for (const auto ehi: _all_ehs)
  {
    if (_same_arc_ehs.find(ehi) != _same_arc_ehs.end())
      continue;

    const Point &pt0 = mesh_.vertex(mesh_.edge(ehi).from_vertex());
    const Point &pt1 = mesh_.vertex(mesh_.edge(ehi).to_vertex());

    Point pt_m(0.);
    auto dist = ACG::Geometry::distPointLine(pt, pt0, pt1, &pt_m);

    if (dist < min_dist)
    {
      min_pt = pt_m;
      min_dist = dist;
    }
  }

//        std::cerr<<"vh "<<_vh<<" "<<pt<<" "<<min_pt<<std::endl;

  return std::make_pair(min_dist, min_pt);
}


template<class MeshT>
void
MeshOptimizationT<MeshT>::
project_vertex_to_original_surface(const VH _vh, const Point &_orig_pt)
{
  //current point
  auto pt = mesh_.vertex(_vh);
  std::vector<std::vector<Point>> cells_points;
  for (auto ch: mesh_.vertex_cells(_vh))
  {
    auto cell_points = get_cell_points(ch, _vh);
    cells_points.push_back(cell_points);
  }

  //closest to current point
  auto closest_cur = bsp_->nearest(pt);
  auto fh_cur = closest_cur.handle;
  auto normal_cur = face_normal_[fh_cur];

  //closest to original point
  auto normal_orig = get_normal_on_reference_mesh(_orig_pt);
//        auto closest_orig = bsp_->nearest(_orig_pt);
//        auto fh_orig = closest_orig.handle;
//        auto normal_orig = face_normal_[fh_orig];

  Point p_best;
  //if normal is very different from the original normal, it may project to wrong surface
  if (normal_cur.dot(normal_orig) < 0.707)
  {
    p_best = _orig_pt;
    std::cerr << "Warning: prevent wrong projection at vertex " << _vh << std::endl;
  }
  else
  {
    auto fv_it = sf_mesh_.fv_iter(fh_cur);
    const Point &pt0 = sf_mesh_.vertex(*fv_it);
    const Point &pt1 = sf_mesh_.vertex(*(++fv_it));
    const Point &pt2 = sf_mesh_.vertex(*(++fv_it));
    double dist = ACG::Geometry::distPointTriangleSquaredStable(pt, pt0, pt1, pt2, p_best);
  }

  Point dir, query_pt;
  dir = p_best - pt;

  query_pt = p_best;
  int num = 0;
  double coeff = 1.;
  while (num < 10 && (has_invalid_incident_tet(cells_points, query_pt) || has_invalid_incident_triangle(_vh, query_pt)))
  {
    coeff *= 0.8;
    query_pt = pt + coeff * dir;

    num++;
  }

  if (num < 10)
    mesh_.set_vertex(_vh, query_pt);

//        std::cout<<"Projected "<<coeff<<" of to the original "<<std::endl;
}

template<class MeshT>
void
MeshOptimizationT<MeshT>::
project_vertex_to_original_feature_arc(const VH _vh, const Point &_orig_pt)
{
  //current point
  auto pt = mesh_.vertex(_vh);
  std::vector<std::vector<Point>> cells_points;
  for (auto ch: mesh_.vertex_cells(_vh))
  {
    auto cell_points = get_cell_points(ch, _vh);
    cells_points.push_back(cell_points);
  }

  //closest to current point
  //get fe id
  int fe_id = -1;
  for (auto ve: mesh_.vertex_edges(_vh))
    if (feature_edge_[ve] > 0)
    {
      fe_id = feature_edge_[ve];
      break;
    }
  double min_dist = std::numeric_limits<double>::max();
  Point p_best(0.);

  for (auto ehi: sm_fes_[fe_id])
  {
    const Point &pt0 = sf_mesh_.vertex(sf_mesh_.edge(ehi).from_vertex());
    const Point &pt1 = sf_mesh_.vertex(sf_mesh_.edge(ehi).to_vertex());

    Point pt_m(0.);
    auto dist = ACG::Geometry::distPointLine(pt, pt0, pt1, &pt_m);

    if (dist < min_dist)
    {
      p_best = pt_m;
      min_dist = dist;
    }
  }

  Point dir, query_pt;
  dir = p_best - pt;

  query_pt = p_best;
  int num = 0;
  double coeff = 1.;
  while (num < 10 && (has_invalid_incident_tet(cells_points, query_pt) || has_invalid_incident_triangle(_vh, query_pt)))
  {
    coeff *= 0.8;
    query_pt = pt + coeff * dir;

    num++;
  }

  if (num < 10)
    mesh_.set_vertex(_vh, query_pt);
}


template<class MeshT>
bool MeshOptimizationT<MeshT>::
has_invalid_incident_tet(std::vector<std::vector<Point>> &_cells_pts, const Point &_query_pt) const
{
  for (auto &c_pts: _cells_pts)
  {
    c_pts[0] = _query_pt;
    if (get_cell_volume(c_pts) <= 1e-16)
      return true;
  }

  return false;
}


template<class MeshT>
bool MeshOptimizationT<MeshT>::
has_invalid_incident_triangle(const VH _vh, const Point &_query_pt) const
{
  auto normal = get_normal_on_reference_mesh(_query_pt);

  std::vector<std::array<double, 9>> fs_pts;
  fs_pts.reserve(3);

  Point ppx;
  for (auto vhf_it = mesh_.vhf_iter(_vh); vhf_it.valid(); ++vhf_it)
  {
    if (mesh_.is_boundary(*vhf_it))
    {
      auto hfvhs = mesh_.get_halfface_vertices(*vhf_it, _vh);

      std::array<double, 9> fpts;
      for (int i = 0; i < 3; ++i)
        fpts[i] = _query_pt[i];

      for (int i = 1; i < 3; ++i)
      {
        SingularVertexOptProblemQuaternion::project_on_tangent_plane(_query_pt, mesh_.vertex(hfvhs[i]), normal, ppx);
        for (int j = 0; j < 3; ++j)
          fpts[3 * i + j] = ppx[j];
      }

      fs_pts.push_back(fpts);
    }
  }
  for (const auto &f_pts: fs_pts)
  {
    if (!std::isfinite(SIMRMEnergy::eval_f(f_pts.data())))
      return true;
  }

  return false;
}


template<class MeshT>
int
MeshOptimizationT<MeshT>::n_singular_vertex_in_cell(const CH _ch) const
{
  int count = 0;
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
    if (sgl_vt_[*cv_it] > 0)
      count++;

  return count;
}

template<class MeshT>
int
MeshOptimizationT<MeshT>::n_singular_vertex_in_cell2(std::map<VH, int> &_vh_to_idx, const CH _ch) const
{
  int count = 0;
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
    if (sgl_vt_[*cv_it] > 0 && _vh_to_idx.find(*cv_it) != _vh_to_idx.end())
      count++;

  return count;
}

template<class MeshT>
int
MeshOptimizationT<MeshT>::n_incident_singular_edges(const VH _vh) const
{
  int count = 0;
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (valence_[*ve_it] != 0)
      count++;

  return count;
}


template<class MeshT>
bool
MeshOptimizationT<MeshT>::is_movable(const VH _vh) const
{
  int count = 0;
  bool is_bdy = mesh_.is_boundary(_vh);
  for (auto ve_it = mesh_.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (feature_edge_[*ve_it] > 0 || valence_[*ve_it] != 0)
    {
      if (mesh_.is_boundary(*ve_it) != is_bdy)
        return false;
      count++;
    }

  if (count == 2)
    return true;

  return false;
}


template<class MeshT>
std::vector<typename MeshT::PointT>
MeshOptimizationT<MeshT>::get_cell_points(const CH _ch, const VH _vh) const
{
  auto cell_vertices = mesh_.get_cell_vertices(_ch, _vh);
  std::vector<Point> cell_points;
  cell_points.reserve(4);
  for (const auto cv: cell_vertices)
    cell_points.push_back(mesh_.vertex(cv));

  return cell_points;
}


template<class MeshT>
std::vector<VH>
MeshOptimizationT<MeshT>::get_cell_vertices(const CH _ch) const
{
  std::vector<VH> sgl_vhs;
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
    if (sgl_vt_[*cv_it] > 0)
      sgl_vhs.push_back(*cv_it);

  if (sgl_vhs.empty())
    std::cout << "Error: input cell has no singular vertex!" << std::endl;

  if (sgl_vhs.size() == 1)
    return mesh_.get_cell_vertices(_ch, sgl_vhs[0]);
  else if (sgl_vhs.size() == 2)
  {
    auto heh = mesh_.find_halfedge(sgl_vhs[0], sgl_vhs[1]);
    HFH hfh;
    for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
    {
      if (mesh_.incident_cell(*hehf_it) == _ch)
        hfh = *hehf_it;
    }

    return mesh_.get_cell_vertices(hfh, heh);
  }
  else if (sgl_vhs.size() == 3)
  {
    auto hfh = mesh_.find_halfface(sgl_vhs);
    if (mesh_.incident_cell(hfh) != _ch)
      hfh = mesh_.opposite_halfface_handle(hfh);

    return mesh_.get_cell_vertices(hfh);
  }


  return mesh_.get_cell_vertices(_ch);
}

template<class MeshT>
std::vector<VH>
MeshOptimizationT<MeshT>::get_cell_vertices2(std::map<VH, int> &_vh_to_idx, const CH _ch) const
{
  std::vector<VH> sgl_vhs;
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
    if (sgl_vt_[*cv_it] > 0 && _vh_to_idx.find(*cv_it) != _vh_to_idx.end())
      sgl_vhs.push_back(*cv_it);

  if (sgl_vhs.empty())
    std::cout << "Error: input cell has no singular vertex!" << std::endl;

  if (sgl_vhs.size() == 1)
    return mesh_.get_cell_vertices(_ch, sgl_vhs[0]);
  else if (sgl_vhs.size() == 2)
  {
    auto heh = mesh_.find_halfedge(sgl_vhs[0], sgl_vhs[1]);
    HFH hfh;
    for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
    {
      if (mesh_.incident_cell(*hehf_it) == _ch)
        hfh = *hehf_it;
    }

    return mesh_.get_cell_vertices(hfh, heh);
  }
  else if (sgl_vhs.size() == 3)
  {
    auto hfh = mesh_.find_halfface(sgl_vhs);
    if (mesh_.incident_cell(hfh) != _ch)
      hfh = mesh_.opposite_halfface_handle(hfh);

    return mesh_.get_cell_vertices(hfh);
  }


  return mesh_.get_cell_vertices(_ch);
}

template<class MeshT>
std::vector<VH>
MeshOptimizationT<MeshT>::get_cell_vertices_post(const CH _ch) const
{
  std::vector<VH> key_vhs;
  for (auto cv_it = mesh_.cv_iter(_ch); cv_it.valid(); ++cv_it)
    if (sgl_vt_[*cv_it] || feature_edge_vertex_[*cv_it])
      key_vhs.push_back(*cv_it);

  if (key_vhs.empty())
    std::cout << "Error: input cell has no singular vertex!" << std::endl;

  if (key_vhs.size() == 1)
    return mesh_.get_cell_vertices(_ch, key_vhs[0]);
  else if (key_vhs.size() == 2)
  {
    auto heh = mesh_.halfedge(key_vhs[0], key_vhs[1]);
    HFH hfh;
    for (auto hehf_it = mesh_.hehf_iter(heh); hehf_it.valid(); ++hehf_it)
    {
      if (mesh_.incident_cell(*hehf_it) == _ch)
        hfh = *hehf_it;
    }

    return mesh_.get_cell_vertices(hfh, heh);
  }
  else if (key_vhs.size() == 3)
  {
    auto hfh = mesh_.halfface(key_vhs);
    if (mesh_.incident_cell(hfh) != _ch)
      hfh = mesh_.opposite_halfface_handle(hfh);

    return mesh_.get_cell_vertices(hfh);
  }


  return mesh_.get_cell_vertices(_ch);
}


template<class MeshT>
void MeshOptimizationT<MeshT>::save_singular_graph_development(MeshT &_sg_mesh, int _iter)
{
  auto sg_edge_idx = _sg_mesh.template request_edge_property<int>("animation_idx");
  _sg_mesh.set_persistent(sg_edge_idx, true);
  auto edge_color = _sg_mesh.template request_edge_property<OVM::Vec4f>("edge_color");
  _sg_mesh.set_persistent(edge_color, true);
  auto sg_node_idx = _sg_mesh.template request_vertex_property<int>("animation_idx");
  _sg_mesh.set_persistent(sg_node_idx, true);
  auto node_color = _sg_mesh.template request_vertex_property<OVM::Vec4f>("vertex_color");
  _sg_mesh.set_persistent(node_color, true);

  std::map<VH, VH> old_to_new_vh;
  for (const auto &eh: mesh_.edges())
  {
    if (valence_[eh] != 0 || feature_edge_[eh] > 0)
    {
      auto vh0 = mesh_.edge(eh).from_vertex();
      auto vh1 = mesh_.edge(eh).to_vertex();
      if (old_to_new_vh.find(vh0) == old_to_new_vh.end())
      {
        auto newvh0 = _sg_mesh.add_vertex(mesh_.vertex(vh0));
        old_to_new_vh[vh0] = newvh0;
      }

      if (old_to_new_vh.find(vh1) == old_to_new_vh.end())
      {
        auto newvh1 = _sg_mesh.add_vertex(mesh_.vertex(vh1));
        old_to_new_vh[vh1] = newvh1;
      }
    }
  }

  //color edges
  OVM::Vec4f color_vec;
  for (const auto eh: mesh_.edges())
  {
    if (valence_[eh] != 0)
    {
      auto vh0 = mesh_.edge(eh).from_vertex();
      auto vh1 = mesh_.edge(eh).to_vertex();

      auto newvh0 = old_to_new_vh[vh0];
      auto newvh1 = old_to_new_vh[vh1];

      auto neweh = _sg_mesh.add_edge(newvh0, newvh1);

      sg_edge_idx[neweh] = _iter;

      get_edge_color(eh, color_vec);
      edge_color[neweh] = color_vec;

      sg_node_idx[newvh0] = _iter;
      sg_node_idx[newvh1] = _iter;

      get_vertex_color(vh0, color_vec);
      node_color[newvh0] = color_vec;
      get_vertex_color(vh1, color_vec);
      node_color[newvh1] = color_vec;
    }
    else if (feature_edge_[eh] > 0)
    {
      auto vh0 = mesh_.edge(eh).from_vertex();
      auto vh1 = mesh_.edge(eh).to_vertex();

      auto newvh0 = old_to_new_vh[vh0];
      auto newvh1 = old_to_new_vh[vh1];

      auto neweh = _sg_mesh.add_edge(newvh0, newvh1);

      sg_edge_idx[neweh] = _iter;

      edge_color[neweh] = OVM::Vec4f(1, 0, 0, 1);
    }
  }


  //debug: add normal for visualization
//        for(const auto vh : mesh_.vertices()) {
//            if(sgl_vt_[vh] == 1 && mesh_.is_boundary(vh)) {
//                auto normal = get_normal_on_reference_mesh(mesh_.vertex(vh));
//
//                auto newvh0 = _sg_mesh.add_vertex(mesh_.vertex(vh));
//                auto newvh1 = _sg_mesh.add_vertex(mesh_.vertex(vh) + normal);
//                auto neweh = _sg_mesh.add_edge(newvh0, newvh1);
//
//                sg_edge_idx[neweh] = _iter;
//                edge_color_id[neweh] = -1;
//            }
//        }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::save_singular_graph(const std::string &_filename, int _id, const std::set<VH> &_fxb_tps)
{
  if (!_filename.empty())
  {
    auto file_name = _filename + std::to_string(_id) + ".ovm";
    MeshT _sg_mesh;
    auto edge_color = _sg_mesh.template request_edge_property<OVM::Vec4f>("edge_color");
    _sg_mesh.set_persistent(edge_color, true);
    auto node_color = _sg_mesh.template request_vertex_property<OVM::Vec4f>("vertex_color");
    _sg_mesh.set_persistent(node_color, true);

    auto org_vh = _sg_mesh.template request_vertex_property<int>("old_vh_id");
    _sg_mesh.set_persistent(org_vh, true);

    std::map<VH, VH> old_to_new_vh;
    for (const auto &eh: mesh_.edges())
    {
      if (valence_[eh] != 0 || feature_edge_[eh] > 0)
      {
        auto vh0 = mesh_.edge(eh).from_vertex();
        auto vh1 = mesh_.edge(eh).to_vertex();
        if (old_to_new_vh.find(vh0) == old_to_new_vh.end())
        {
          auto newvh0 = _sg_mesh.add_vertex(mesh_.vertex(vh0));
          old_to_new_vh[vh0] = newvh0;
          org_vh[newvh0] = vh0.idx();
        }

        if (old_to_new_vh.find(vh1) == old_to_new_vh.end())
        {
          auto newvh1 = _sg_mesh.add_vertex(mesh_.vertex(vh1));
          old_to_new_vh[vh1] = newvh1;
          org_vh[newvh1] = vh1.idx();
        }
      }
    }

    //color edges
    OVM::Vec4f color_vec;
    for (const auto &eh: mesh_.edges())
    {
      if (valence_[eh] != 0)
      {
        auto vh0 = mesh_.edge(eh).from_vertex();
        auto vh1 = mesh_.edge(eh).to_vertex();

        auto newvh0 = old_to_new_vh[vh0];
        auto newvh1 = old_to_new_vh[vh1];

        auto neweh = _sg_mesh.add_edge(newvh0, newvh1);

        get_edge_color(eh, color_vec);
        edge_color[neweh] = color_vec;

        get_vertex_color(vh0, color_vec);
        node_color[newvh0] = color_vec;
        get_vertex_color(vh1, color_vec);
        node_color[newvh1] = color_vec;
      }
      else if (feature_edge_[eh] > 0)
      {
        auto vh0 = mesh_.edge(eh).from_vertex();
        auto vh1 = mesh_.edge(eh).to_vertex();

        auto newvh0 = old_to_new_vh[vh0];
        auto newvh1 = old_to_new_vh[vh1];

        auto neweh = _sg_mesh.add_edge(newvh0, newvh1);

        edge_color[neweh] = OVM::Vec4f(1, 0, 0, 1);
      }
    }

    for (auto vhi: _fxb_tps)
      node_color[old_to_new_vh[vhi]] = OVM::Vec4f(1, 0, 0, 1);


    OpenVolumeMesh::IO::FileManager fm;
    fm.writeFile(file_name, _sg_mesh);
  }
}

template<class MeshT>
void MeshOptimizationT<MeshT>::get_edge_color(const EH _eh, OVM::Vec4f &_color) const
{
  //colour singular edges
  if (valence_[_eh] == 0)
  {
    _color[0] = 1.;
    _color[1] = 1.;
    _color[2] = 1.;
    _color[3] = 0.5;
  }
  else if (valence_[_eh] == -1)
  {//valance -1, green
    _color[0] = 0.;
    _color[1] = 1.;
    _color[2] = 0.;
    _color[3] = .8;
  }
  else if (valence_[_eh] == 1)
  {//valance 1, blue
    _color[0] = 0.;
    _color[1] = 0.;
    _color[2] = 1.;
    _color[3] = .8;
  }
  else if (valence_[_eh] == std::numeric_limits<int>::lowest())
  {//edge direction is orthogonal to rotational axis
    _color[0] = 1.;
    _color[1] = 0.;
    _color[2] = 0.;
    _color[3] = .8;
  }
  else if (valence_[_eh] == -2)
  {//valence -2
    _color[0] = 1.;
    _color[1] = 1.;
    _color[2] = 0.;
    _color[3] = .8;
  }
  else if (valence_[_eh] == 2)
  {//valence +2
    _color[0] = 0.;
    _color[1] = 1.;
    _color[2] = 1.;
    _color[3] = .8;
  }
  else if (valence_[_eh] == -3)
  {//valence -3
    _color[0] = 1.;
    _color[1] = 0.75;
    _color[2] = 0.;
    _color[3] = .8;
  }
  else if (valence_[_eh] == 3)
  {//valence +3
    _color[0] = 1.;
    _color[1] = .5;
    _color[2] = 0.;
    _color[3] = .8;
  }
  else if (valence_[_eh] == 4)
  {//valence +4
    _color[0] = 1.;
    _color[1] = 0.8;
    _color[2] = 0.8;
    _color[3] = .8;
  }
  else if (valence_[_eh] == std::numeric_limits<int>::max())
  {//3D transition case, black
    _color[0] = 0.;
    _color[1] = 0.;
    _color[2] = 0.;
    _color[3] = .8;
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::get_vertex_color(const VH _vh, OVM::Vec4f &_color) const
{
  int node_type = MeshPropertiesT<MeshT>::node_index(_vh);
  //color singular nodes
  if (node_type == 0)
  {
    _color[0] = 1.;
    _color[1] = 1.;
    _color[2] = 1.;
    _color[3] = 0.5;
  }
  else if (node_type == 1 || node_type == -1)
  {
    _color[0] = 0.25;
    _color[1] = 0.8;
    _color[2] = .0;
    _color[3] = .8;
  }
  else if (node_type == 2)
  {
    _color[0] = 0.;
    _color[1] = 0.5;
    _color[2] = 0.5;
    _color[3] = .8;
  }
  else if (node_type == 3)
  {
    _color[0] = 0.;
    _color[1] = 0.25;
    _color[2] = 0.75;
    _color[3] = .8;
  }
  else if (node_type == 4 || node_type == -7)
  {
    _color[0] = 0.25;
    _color[1] = 0.;
    _color[2] = .8;
    _color[3] = .8;
  }
  else if (node_type == 10)
  {
    _color[0] = 0.5;
    _color[1] = 0.5;
    _color[2] = 0.5;
    _color[3] = .8;
  }
  else if (node_type == -2)
  {
    _color[0] = 1.;
    _color[1] = 1.;
    _color[2] = 0.;
    _color[3] = .8;
  }
  else if (node_type == -3)
  {
    _color[0] = 1.;
    _color[1] = 0.75;
    _color[2] = 0.25;
    _color[3] = .8;
  }
  else if (node_type == -4)
  {
    _color[0] = 1.;
    _color[1] = 0.6;
    _color[2] = 0.4;
    _color[3] = .8;
  }
  else if (node_type == -5)
  {
    _color[0] = 1.;
    _color[1] = 0.5;
    _color[2] = 0.5;
    _color[3] = .8;
  }
  else if (node_type == -6)
  {
    _color[0] = 1.;
    _color[1] = 0.25;
    _color[2] = 0.75;
    _color[3] = .8;
  }
  else if (node_type == -8)
  {
    _color[0] = 1.;
    _color[1] = 0.1;
    _color[2] = 0.9;
    _color[3] = .8;
  }
  else if (node_type == -9)
  {
    _color[0] = 0.;
    _color[1] = 1.;
    _color[2] = 1.;
    _color[3] = .8;
  }
  else if (node_type == -10)
  {
    _color[0] = 0.;
    _color[1] = 1.;
    _color[2] = 0.;
    _color[3] = .8;
  }
  else if (node_type == -11)
  {
    _color[0] = 0.;
    _color[1] = 0.;
    _color[2] = 1.;
    _color[3] = .8;
  }
  else if (node_type == -12)
  {
    _color[0] = 1.;
    _color[1] = .8;
    _color[2] = .6;
    _color[3] = .8;
  }
  else if (node_type == -13)
  {
    _color[0] = 0.8;
    _color[1] = 0.6;
    _color[2] = 1.;
    _color[3] = .8;
  }
  else if (node_type == -14)
  {
    _color[0] = 0.6;
    _color[1] = .8;
    _color[2] = 1.;
    _color[3] = .8;
  }
  else if (node_type == -15)
  {
    _color[0] = 0.4;
    _color[1] = .6;
    _color[2] = .8;
    _color[3] = .8;
  }
  else if (node_type == -16)
  {
    _color[0] = 0.8;
    _color[1] = .6;
    _color[2] = .4;
    _color[3] = .8;
  }
  else if (node_type == -17)
  {
    _color[0] = 0.6;
    _color[1] = .8;
    _color[2] = .4;
    _color[3] = .8;
  }
  else if (node_type == -18)
  {
    _color[0] = 1.;
    _color[1] = .5;
    _color[2] = 0.;
    _color[3] = .8;
  }
  else if (node_type == -19)
  {
    _color[0] = 1.;
    _color[1] = .8;
    _color[2] = .8;
    _color[3] = .8;
  }
  else
  {
    _color[0] = 0.;
    _color[1] = 0.;
    _color[2] = 0.;
    _color[3] = .8;
  }
}


template<class MeshT>
void MeshOptimizationT<MeshT>::adapt_target_length(const double _damping, const int _norm_order)
{
  auto max_tl = *std::max_element(target_length_.begin(), target_length_.end());
  auto min_tl = *std::min_element(target_length_.begin(), target_length_.end());

  ALGOHEX_DEBUG_ONLY(
          std::cerr << "Before max target length: " << max_tl << " min target length: " << min_tl << std::endl;)

  auto smoothness_vprop = mesh_.template request_vertex_property<double>("field_smoothness", 0.);

  measure_field_smoothness(smoothness_vprop, _norm_order);

  for (const auto vhi: mesh_.vertices())
  {
    double factor = std::pow((smoothness_vprop[vhi] + 1.), 2);
    factor *= _damping;

    factor = std::max(1., factor);

    //clamp
    double tl_new = std::max(min_target_len_, 1. / factor * target_length_[vhi]);

    target_length_[vhi] = tl_new;
  }

  auto max_tl_new = *std::max_element(target_length_.begin(), target_length_.end());
  auto min_tl_new = *std::min_element(target_length_.begin(), target_length_.end());

  ALGOHEX_DEBUG_ONLY(
          std::cerr << "Max target length: " << max_tl_new << " min target length: " << min_tl_new << std::endl;)

  MeshPropertiesT<MeshT>::smooth_target_edge_length(3);

  max_tl_new = *std::max_element(target_length_.begin(), target_length_.end());
  min_tl_new = *std::min_element(target_length_.begin(), target_length_.end());
  ALGOHEX_DEBUG_ONLY(
          std::cerr << "After smooth, max target length: " << max_tl_new << " min target length: " << min_tl_new
                    << std::endl;)

  //scale to original max target length
  double scale = (max_tl - min_tl_new) / (max_tl_new - min_tl_new);
  scale = std::max(1., scale);

  for (auto vhi: mesh_.vertices())
  {
    target_length_[vhi] = scale * (target_length_[vhi] - min_tl_new) + min_tl_new;
  }

  max_tl_new = *std::max_element(target_length_.begin(), target_length_.end());
  min_tl_new = *std::min_element(target_length_.begin(), target_length_.end());
  ALGOHEX_DEBUG_ONLY(
          std::cerr << "After scaling, max target length: " << max_tl_new << " min target length: " << min_tl_new
                    << std::endl;)
}


template<class MeshT>
void MeshOptimizationT<MeshT>::measure_field_smoothness(VP<double> &_smth_vprop, const int _norm_order)
{
  auto smth_fprop = mesh_.template request_face_property<double>("field_smoothness_face", 0.);
  mesh_.set_persistent(smth_fprop, true);

  for (const auto fhi: mesh_.faces())
  {
    if (mesh_.is_boundary(fhi))
      continue;

    auto hfh0 = mesh_.halfface_handle(fhi, 0);
    auto hfh1 = mesh_.halfface_handle(fhi, 1);

    auto ch0 = mesh_.incident_cell(hfh0);
    auto ch1 = mesh_.incident_cell(hfh1);

    Quaternion tranq = tq_.transition(trans_prop_[hfh0]);

    Quaternion q0in1 = cell_quaternions_[ch0] * tranq;

    Quaternion rt = cell_quaternions_[ch1] * q0in1.conjugate();
    smth_fprop[fhi] = 2 * std::acos(std::min(std::fabs(rt.w()), 1.));

    if (smth_fprop[fhi] <= M_PI / 180. * 20)
      smth_fprop[fhi] = 0;
  }

  int n10 = 0, n20 = 0, n30 = 0, n45 = 0, nbbb = 0;
  double max_rt_agl = 0.;
  for (const auto fhi: mesh_.faces())
  {
    if (smth_fprop[fhi] <= M_PI / 180. * 10 && smth_fprop[fhi] >= 0.)
      n10++;
    else if (smth_fprop[fhi] <= M_PI / 180. * 20 && smth_fprop[fhi] > M_PI / 180. * 10.)
      n20++;
    else if (smth_fprop[fhi] <= M_PI / 180. * 30 && smth_fprop[fhi] > M_PI / 180. * 20.)
      n30++;
    else if (smth_fprop[fhi] <= M_PI / 180. * 45 && smth_fprop[fhi] > M_PI / 180. * 30.)
      n45++;
    else if (smth_fprop[fhi] > M_PI / 180. * 45.)
    {
      ALGOHEX_DEBUG_ONLY(std::cerr << " face: " << fhi << " angle: " << smth_fprop[fhi] * 180. / M_PI << std::endl;)
      nbbb++;
    }

    if (smth_fprop[fhi] > max_rt_agl)
      max_rt_agl = smth_fprop[fhi];
  }

  ALGOHEX_DEBUG_ONLY(std::cerr << "####Number of faces w.r.t. rotation angle: " << std::endl;
                             std::cerr << "0-10 degree: " << n10 << std::endl;
                             std::cerr << "10-20 degree: " << n20 << std::endl;
                             std::cerr << "20-30 degree: " << n30 << std::endl;
                             std::cerr << "30-45 degree: " << n45 << std::endl;
                             std::cerr << ">45 degree: " << nbbb << std::endl;
                             std::cerr << "Max rotation angle: " << max_rt_agl * 180. / M_PI << " degree."
                                       << std::endl;)


  for (const auto vhi: mesh_.vertices())
  {
    _smth_vprop[vhi] = measure_field_smoothness_at_vertex(vhi, smth_fprop, _norm_order);
  }
}

template<class MeshT>
double MeshOptimizationT<MeshT>::measure_field_smoothness_at_vertex(const VH _vh, const FP<double> &_smth_fprop,
                                                                    int _norm_order) const
{
  double sm = 0.;
  for (auto vf_it = mesh_.vf_iter(_vh); vf_it.valid(); ++vf_it)
  {
    if (mesh_.is_boundary(*vf_it))
      continue;

    sm += std::pow(_smth_fprop[*vf_it], (double) _norm_order);
  }

  return std::pow(sm, 1. / (double) _norm_order);
}


template<class MeshT>
void MeshOptimizationT<MeshT>::set_min_target_length(const double _factor)
{
  double min_tl = *std::min_element(target_length_.begin(), target_length_.end());

  min_target_len_ = _factor * min_tl;
}
}


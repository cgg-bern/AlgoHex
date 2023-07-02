/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#define COMMONFUNCST_C

#include "CommonFuncs.hh"


namespace AlgoHex
{
template<class MeshT>
typename MeshT::PointT face_normal(const MeshT &mesh, const HFH _hfh)
{
  typename MeshT::PointT p0, p1, p2;
  auto hfvs = mesh.get_halfface_vertices(_hfh);
  p0 = mesh.vertex(hfvs[0]);
  p1 = mesh.vertex(hfvs[1]);
  p2 = mesh.vertex(hfvs[2]);

  typename MeshT::PointT vn, v1, v2;
  v1 = p1 - p0;
  v2 = p2 - p1;
  vn = v1 % v2;

  if (vn.norm() < std::numeric_limits<double>::min())
    std::cout << "Warning: zero normal!";

  return vn;
}

template<class MeshT>
typename MeshT::PointT face_normal(const MeshT &mesh, const VH vh0, const VH vh1, const VH vh2)
{
  typename MeshT::PointT p0, p1, p2;
  p0 = mesh.vertex(vh0);
  p1 = mesh.vertex(vh1);
  p2 = mesh.vertex(vh2);

  typename MeshT::PointT vn, v1, v2;
  v1 = p1 - p0;
  v2 = p2 - p1;
  vn = v1 % v2;

  assert(std::isfinite(vn.norm()));

  vn.normalize();

  return vn;
}

template<class MeshT>
typename MeshT::PointT face_normal(const MeshT &mesh, const std::vector<VH> &vhs)
{
  typename MeshT::PointT p0, p1, p2;
  p0 = mesh.vertex(vhs[0]);
  p1 = mesh.vertex(vhs[1]);
  p2 = mesh.vertex(vhs[2]);

  typename MeshT::PointT vn, v1, v2;
  v1 = p1 - p0;
  v2 = p2 - p1;
  vn = v1 % v2;

  assert(std::isfinite(vn.norm()));

  vn.normalize();

  return vn;
}

template<class MeshT>
typename MeshT::PointT vertex_normal_area_weighted(const MeshT &mesh, const VH vh)
{
  typename MeshT::PointT vn(0);
  for (auto vhf_it = mesh.vhf_iter(vh); vhf_it.valid(); ++vhf_it)
    if (mesh.is_boundary(*vhf_it))
    {
      vn += face_normal(mesh, *vhf_it);
    }

  vn.normalize();

  return vn;
}

template<class MeshT>
typename MeshT::PointT vertex_normal_uniform_weighted(const MeshT &mesh, const VH vh)
{
  typename MeshT::PointT vn(0);
  for (auto vhf_it = mesh.vhf_iter(vh); vhf_it.valid(); ++vhf_it)
    if (mesh.is_boundary(*vhf_it))
    {
      vn += face_normal(mesh, *vhf_it).normalize();
    }

  vn.normalize();

  return vn;
}

template<typename Vec3>
std::array<double, 3> compute_bary_coord(const Vec3 &_px,
                                         const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2)
{
  std::array<double, 3> bary_coord;
  auto v0x = _p0 - _px;
  auto v1x = _p1 - _px;
  auto v2x = _p2 - _px;
  auto vn12 = v1x % v2x;
  auto vn20 = v2x % v0x;

  auto v10 = _p1 - _p0;
  auto v20 = _p2 - _p0;
  auto area = (v10 % v20).length();

  bary_coord[0] = vn12.length() / area;
  bary_coord[1] = vn20.length() / area;
  bary_coord[2] = 1.0 - bary_coord[0] - bary_coord[1];

  return bary_coord;
}

template<class MeshT>
double edge_angle(const MeshT &mesh, const HEH heh0, const HEH heh1)
{
  auto pt0 = mesh.vertex(mesh.halfedge(heh0).from_vertex());
  auto pt1 = mesh.vertex(mesh.halfedge(heh0).to_vertex());
  auto pt2 = mesh.vertex(mesh.halfedge(heh1).from_vertex());
  auto pt3 = mesh.vertex(mesh.halfedge(heh1).to_vertex());

  auto v0 = pt1 - pt0;
  v0.normalize();

  auto v1 = pt3 - pt2;
  v1.normalize();

  auto dprd = v0.dot(v1);
  dprd = std::min(1., std::max(-1., dprd));

//        std::cerr<<"edge: "<<mesh_.edge_handle(_heh0)<<" "<<mesh_.edge_handle(_heh1)<<" angle: "<<std::acos(dprd)*180./M_PI<<std::endl;

  return std::acos(dprd);
}


template<class MeshT>
void update_cells_volume(MeshT &mesh, CP<double> &cell_volume, double &total_volume)
{
  total_volume = 0;
  for (const auto ch: mesh.cells())
  {
    cell_volume[ch] = get_cell_volume(mesh, ch);
    total_volume += cell_volume[ch];
  }

  double minvol = *std::min_element(cell_volume.begin(), cell_volume.end());
  std::cout << "Min cell volume: " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << minvol
            << " Max cell volume: " << *std::max_element(cell_volume.begin(), cell_volume.end()) << std::setprecision(6)
            << std::endl;
  if (minvol <= 1e-16)
    std::cout << "Error: degenerate cell or flipped cell. Volume "
              << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << minvol << std::setprecision(6)
              << std::endl;

//        double min_len = DBL_MAX, max_len = DBL_MIN;
//        for(const auto eh : mesh.edges())
//            edge_length_[eh] = mesh.length(eh);

//        std::cout<<"Min edge length: "<<*std::min_element(edge_length_.begin(), edge_length_.end())
//        <<" Max edge length: "<<*std::max_element(edge_length_.begin(), edge_length_.end())<<std::endl;
}

template<class MeshT>
AlgoHex::FH common_face(const MeshT &mesh, const CH _ch0, const CH _ch1)
{
  for (auto chf_it = mesh.chf_iter(_ch0); chf_it.valid(); ++chf_it)
  {
    if (mesh.incident_cell(mesh.opposite_halfface_handle(*chf_it)) == _ch1)
      return mesh.face_handle(*chf_it);
  }

  return FH(-1);
}

template<class MeshT>
void build_surface_mesh(MeshT &mesh, MeshT &sf_mesh)
{
  using Point = typename MeshT::PointT;
  auto feature_fprop = mesh.template request_face_property<int>("AlgoHex::FeatureFaces");
  auto feature_face_vertex = mesh.template request_vertex_property<bool>("AlgoHex::FeatureFaceVertices");
  //get surface mesh of tetmesh
  std::map<VH, VH> vh_to_sfvh;
  for (const auto vhi: mesh.vertices())
    if (feature_face_vertex[vhi])
    {
      auto tri_vh = sf_mesh.add_vertex(mesh.vertex(vhi));
      vh_to_sfvh[vhi] = tri_vh;
    }

  //add trimesh faces
  std::map<HFH, FH> hfh_to_sffh;
  std::vector<bool> visited_fhs(mesh.n_faces(), false);
  //add boundary
  for (const auto hfi: mesh.halffaces())
  {
    auto fhi = mesh.face_handle(hfi);
    if (feature_fprop[fhi] > 0 && mesh.is_boundary(hfi))
    {
      std::vector<VH> tri_vhs;
      tri_vhs.reserve(3);
      for (auto hfv_it = mesh.hfv_iter(hfi); hfv_it.valid(); ++hfv_it)
        tri_vhs.push_back(vh_to_sfvh[*hfv_it]);

      hfh_to_sffh[hfi] = sf_mesh.add_face(tri_vhs);
      visited_fhs[fhi.idx()] = true;
    }
  }

  //add interior feature
  std::stack<HFH> hfh_st;
  for (const auto hfi: mesh.halffaces())
  {
    auto fhi = mesh.face_handle(hfi);
    if (feature_fprop[fhi] > 0 && !visited_fhs[fhi.idx()])
    {
      hfh_st.push(hfi);
    }
  }
  while (!hfh_st.empty())
  {
    auto hf_cur = hfh_st.top();
    hfh_st.pop();

    auto fh_cur = mesh.face_handle(hf_cur);
    if (visited_fhs[fh_cur.idx()])
      continue;
    else
    {
      std::vector<VH> tri_vhs;
      tri_vhs.reserve(3);
      for (auto hfv_it = mesh.hfv_iter(hf_cur); hfv_it.valid(); ++hfv_it)
        tri_vhs.push_back(vh_to_sfvh[*hfv_it]);

      hfh_to_sffh[hf_cur] = sf_mesh.add_face(tri_vhs);
      visited_fhs[fh_cur.idx()] = true;
    }

    //neigbours
    for (auto hfhe_it = mesh.hfhe_iter(hf_cur); hfhe_it.valid(); ++hfhe_it)
    {
      auto heh_opp = mesh.opposite_halfedge_handle(*hfhe_it);
      for (auto hehf_it = mesh.hehf_iter(heh_opp); hehf_it.valid(); ++hehf_it)
      {
        auto fh_next = mesh.face_handle(*hehf_it);
        if (feature_fprop[fh_next] > 0 && !visited_fhs[fh_next.idx()])
          hfh_st.push(*hehf_it);
      }
    }
  }

  auto face_normal = sf_mesh.template request_face_property<Point>("original face normal");
  auto vertex_normal = sf_mesh.template request_vertex_property<Point>("original vertex normal");
  auto sharp_face = sf_mesh.template request_face_property<bool>("original sharp face");
  auto sm_feature_edge = sf_mesh.template request_edge_property<int>("surface mesh: feature edge");

  auto feature_edge = mesh.template request_edge_property<int>("AlgoHex::FeatureEdges");

  for (const auto eh: mesh.edges())
  {
    if (feature_edge[eh] > 0)
    {
      auto vh0 = mesh.edge(eh).from_vertex();
      auto vh1 = mesh.edge(eh).to_vertex();

      for (auto vhf_it = mesh.vhf_iter(vh0); vhf_it.valid(); ++vhf_it)
        if (feature_fprop[mesh.face_handle(*vhf_it)] > 0)
        {
          if (hfh_to_sffh.find(*vhf_it) != hfh_to_sffh.end())
            sharp_face[hfh_to_sffh[*vhf_it]] = true;
        }

      for (auto vhf_it = mesh.vhf_iter(vh1); vhf_it.valid(); ++vhf_it)
        if (feature_fprop[mesh.face_handle(*vhf_it)] > 0)
        {
          if (hfh_to_sffh.find(*vhf_it) != hfh_to_sffh.end())
            sharp_face[hfh_to_sffh[*vhf_it]] = true;
        }

      auto sm_eh = sf_mesh.edge_handle(sf_mesh.find_halfedge(vh_to_sfvh[vh0], vh_to_sfvh[vh1]));
      sm_feature_edge[sm_eh] = feature_edge[eh];
    }
  }

  for (const auto hfi: mesh.halffaces())
    if (feature_fprop[mesh.face_handle(hfi)] > 0)
    {
      if (hfh_to_sffh.find(hfi) != hfh_to_sffh.end())
        face_normal[hfh_to_sffh[hfi]] = mesh.normal(hfi);
    }

  for (const auto vhi: sf_mesh.vertices())
  {
    Point vn(0.);
    for (auto vf_it = sf_mesh.vf_iter(vhi); vf_it.valid(); ++vf_it)
      vn += face_normal[*vf_it];

    vn.normalize();
    vertex_normal[vhi] = vn;
  }
}


template<class MeshT>
std::unique_ptr<AlgoHex::OpenVolumeMeshTriangleBSPT<MeshT>>
build_bsp(const MeshT &sfmesh)
{
  // create Triangle BSP
  auto triangle_bsp = std::make_unique<AlgoHex::OpenVolumeMeshTriangleBSPT<MeshT>>(sfmesh);

  // build Triangle BSP
  triangle_bsp->reserve(sfmesh.n_faces());

  for (const auto fh: sfmesh.faces())
    triangle_bsp->push_back(fh);

  triangle_bsp->build(10, 100); //max vertices per leaf 10, max depth 100

  // return pointer to triangle BinarySpacePartitionTrees
  return triangle_bsp;
}


template<class MeshT>
bool has_interior_singular_edge(const MeshT &mesh, const EP<int> &valence, const VH vh)
{
  for (auto voe_it = mesh.ve_iter(vh); voe_it.valid(); ++voe_it)
    if (valence[*voe_it] != 0 && !mesh.is_boundary(*voe_it))
      return true;

  return false;
}

template<class MeshT>
int n_interior_complex_singular_edge(const MeshT &mesh, const EP<int> &valence, const VH vh)
{
  int n_int_cse = 0;
  for (auto voe_it = mesh.ve_iter(vh); voe_it.valid(); ++voe_it)
    if (std::abs(valence[*voe_it]) > 1 && !mesh.is_boundary(*voe_it))
      n_int_cse++;

  return n_int_cse;
}

template<class MeshT>
int n_invalid_singular_edge(const MeshT &mesh, const EP<int> &valence, const VH vh)
{
  int n_int_cse = 0;
  for (auto voe_it = mesh.ve_iter(vh); voe_it.valid(); ++voe_it)
    if (valence[*voe_it] == std::numeric_limits<int>::max())
      n_int_cse++;

  return n_int_cse;
}

template<class MeshT>
int n_incident_special_edges(const MeshT &mesh, const EP<int> &feature_edge, const EP<int> &valence, const VH vh)
{
  int n = 0;
  for (auto ve_it = mesh.ve_iter(vh); ve_it.valid(); ++ve_it)
  {
    if (feature_edge[*ve_it] > 0 || valence[*ve_it] != 0)
      n++;
  }

  return n;
}

template<class MeshT>
int n_incident_special_edges_on_feature_face(const MeshT &mesh, const EP<bool> &feature_face_edge,
                                             const EP<int> &feature_edge, const EP<int> &valence, const VH vh)
{
  int n = 0;
  for (auto ve_it = mesh.ve_iter(vh); ve_it.valid(); ++ve_it)
  {
    if (feature_face_edge[*ve_it] && (feature_edge[*ve_it] > 0 || valence[*ve_it] != 0))
      n++;
  }

  return n;
}

template<class MeshT>
int n_incident_special_edges_kept_on_feature_face(const MeshT &mesh, const EP<bool> &feature_face_edge,
                                                  const EP<int> &feature_edge, const EP<int> &valence,
                                                  const EP<bool> &keep_on_surface, const VH vh)
{
  int n = 0;
  for (auto ve_it = mesh.ve_iter(vh); ve_it.valid(); ++ve_it)
  {
    if (feature_face_edge[*ve_it] && (feature_edge[*ve_it] > 0 || (keep_on_surface[*ve_it] && valence[*ve_it] != 0)))
      n++;
  }

  return n;
}

template<class MeshT>
int
n_incident_feature_singular_edges(const MeshT &mesh, const EP<int> &feature_edge, const EP<int> &valence, const VH vh)
{
  int n = 0;
  for (auto ve_it = mesh.ve_iter(vh); ve_it.valid(); ++ve_it)
  {
    //boundary feature halfedges
    if (feature_edge[*ve_it] > 0 && valence[*ve_it] != 0)
      n++;
  }

  return n;
}

template<class MeshT>
int n_incident_singular_edges(const MeshT &mesh, const EP<int> &valence, const VH _vh)
{
  int n_sg = 0;
  for (auto ve_it = mesh.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (valence[*ve_it] != 0)
      n_sg++;

  return n_sg;
}

template<class MeshT>
int n_incident_singular_edges_of_valence_one(const MeshT &mesh, const EP<int> &valence, const VH _vh)
{
  int n_sg = 0;
  for (auto ve_it = mesh.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (valence[*ve_it] == -1 || valence[*ve_it] == 1 || valence[*ve_it] == std::numeric_limits<int>::lowest())
      n_sg++;

  return n_sg;
}

template<class MeshT>
int n_incident_feature_edges(const MeshT &mesh, const EP<int> &feature_edge, const VH _vh)
{
  int n_ft = 0;
  for (auto ve_it = mesh.ve_iter(_vh); ve_it.valid(); ++ve_it)
    if (feature_edge[*ve_it] > 0)
      n_ft++;

  return n_ft;
}

template<class MeshT>
void
get_special_halfedges_at_vertex(const MeshT &mesh, const EP<int> &feature_edge, const EP<int> &valence, const VH vh,
                                std::vector<HEH> &ft_hehs, std::vector<HEH> &bdy_hehs, std::vector<HEH> &itr_hehs)
{
  for (auto voh_it = mesh.voh_iter(vh); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh.edge_handle(*voh_it);
    bool is_bdye = mesh.is_boundary(ve);

    //boundary feature halfedges
    if (feature_edge[ve] > 0)
      ft_hehs.push_back(*voh_it);

    //interior singular halfedges
    if (valence[ve] != 0 && !is_bdye)
      itr_hehs.push_back(*voh_it);

    //boudary singular halfedges
    if (is_bdye && valence[ve] != 0)
      bdy_hehs.push_back(*voh_it);
  }
}

template<class MeshT>
std::vector<HEH> get_singular_halfedges(const MeshT &mesh, const EP<int> &valence, const VH _vh)
{
  std::vector<HEH> sg_hehs;
  for (auto voh_it = mesh.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh.edge_handle(*voh_it);
    if (valence[ve] != 0)
      sg_hehs.push_back(*voh_it);
  }

  return sg_hehs;
}

template<class MeshT>
std::vector<HEH>
get_non_feature_singular_halfedges(const MeshT &mesh, const EP<bool> &feature_feace_edge, const EP<int> &valence,
                                   const VH _vh)
{
  std::vector<HEH> sg_hehs;
  for (auto voh_it = mesh.voh_iter(_vh); voh_it.valid(); ++voh_it)
  {
    auto ve = mesh.edge_handle(*voh_it);
    //interior singular halfedges
    if (valence[ve] != 0 && !feature_feace_edge[ve])
      sg_hehs.push_back(*voh_it);
  }

  return sg_hehs;
}

template<class MeshT>
void
count_special_edges_at_vertex(const MeshT &mesh, const EP<bool> &feature_face_edge, const EP<int> &feature_edge,
                              const EP<int> &valence, const VH vh,
                              int &n_nff_val_1, int &n_nff_val1, int &n_nff_val_u, int &n_nff, int &n_bdy_val_1,
                              int &n_bdy_val1, int &n_bdy_val_u, int &n_bdy, int &n_cmpl,
                              int &n_ff_val_1, int &n_ff_val1, int &n_ff_val_u, int &n_ff, int &n_fe)
{
  n_nff_val_1 = 0, n_nff_val1 = 0, n_nff_val_u = 0, n_nff = 0,
  n_bdy_val_1 = 0, n_bdy_val1 = 0, n_bdy_val_u = 0, n_bdy = 0,
  n_cmpl = 0, n_ff_val_1 = 0, n_ff_val1 = 0, n_ff_val_u = 0, n_ff = 0, n_fe = 0;
  for (auto ve_it = mesh.ve_iter(vh); ve_it.valid(); ++ve_it)
  {
    if (feature_edge[*ve_it] > 0)
      n_fe++;

    if (feature_face_edge[*ve_it])
    {
      if (valence[*ve_it] != 0)
      {
        n_ff++;
        if (mesh.is_boundary(*ve_it))
          n_bdy++;
      }

      if (valence[*ve_it] == 1)
      {
        n_ff_val1++;
        if (mesh.is_boundary(*ve_it))
          n_bdy_val1++;
      }

      if (valence[*ve_it] == -1)
      {
        n_ff_val_1++;
        if (mesh.is_boundary(*ve_it))
          n_bdy_val_1++;
      }

      if (valence[*ve_it] == std::numeric_limits<int>::lowest())
      {
        n_ff_val_u++;
        if (mesh.is_boundary(*ve_it))
          n_bdy_val_u++;
      }
    }
    else
    {
      if (valence[*ve_it] == -1)
        n_nff_val_1++;

      if (valence[*ve_it] == 1)
        n_nff_val1++;

      if (valence[*ve_it] != 0)
        n_nff++;

      if (valence[*ve_it] == std::numeric_limits<int>::lowest())
      {
        n_nff_val_u++;
      }
    }

    if (std::abs(valence[*ve_it]) > 1)
      n_cmpl++;
  }
}

template<class MeshT>
std::set<EH>
incident_non_ffe_singular_edges_in_same_region(const MeshT &mesh, const EP<bool> &feature_face_edge,
                                               const EP<int> &valence,
                                               const VH vh, const std::set<CH> &por_chs)
{
  std::set<EH> nffsge;
  for (const auto chi: por_chs)
  {
    for (auto ce_it = mesh.ce_iter(chi); ce_it.valid(); ++ce_it)
    {
      if (valence[*ce_it] != 0 && !feature_face_edge[*ce_it])
      {
        if (mesh.edge(*ce_it).from_vertex() == vh || mesh.edge(*ce_it).to_vertex() == vh)
        {
          nffsge.insert(*ce_it);
        }
      }
    }
  }

  return nffsge;
}

template<class MeshT>
int
n_incident_non_ffe_singular_edges_in_same_region(const MeshT &mesh, const EP<bool> &feature_face_edge,
                                                 const EP<int> &valence,
                                                 const VH vh, const std::set<CH> &por_chs)
{
  std::set<EH> nffsge = incident_non_ffe_singular_edges_in_same_region(mesh, feature_face_edge, valence, vh, por_chs);

  return nffsge.size();
}

template<class MeshT>
std::vector<HEH>
get_boundary_halfedges_in_feature_sector_ccw(const MeshT &mesh, const EP<int> &feature_edge, const EP<int> &valence,
                                             const HEH _heh_s)
{
  assert(mesh.is_boundary(_heh_s));

  std::vector<HEH> hehs;
  hehs.push_back(_heh_s);

  HEH he_it = _heh_s;
  do
  {
    HFH bdy_hfh;
    for (auto hehf_it = mesh.hehf_iter(he_it); hehf_it.valid(); ++hehf_it)
      if (mesh.is_boundary(*hehf_it))
      {
        bdy_hfh = *hehf_it;
        break;
      }

    he_it = mesh.prev_halfedge_in_halfface(he_it, bdy_hfh);
    he_it = mesh.opposite_halfedge_handle(he_it);

    hehs.push_back(he_it);

    auto ehi = mesh.edge_handle(he_it);
    if (feature_edge[ehi] > 0 || valence[ehi] != 0)
      break;
  }
  while (he_it != _heh_s);

  return hehs;
}

template<class MeshT>
std::vector<HEH>
get_boundary_halfedges_in_feature_sector_cw(const MeshT &mesh, const EP<int> &feature_edge, const HEH _heh_s)
{
  assert(mesh.is_boundary(_heh_s));

  std::vector<HEH> hehs;
  hehs.push_back(_heh_s);

  HEH he_it = _heh_s;
  do
  {
    HFH bdy_hfh = *mesh.hehf_iter(he_it);

    he_it = mesh.prev_halfedge_in_halfface(he_it, bdy_hfh);
    he_it = mesh.opposite_halfedge_handle(he_it);

    hehs.push_back(he_it);

    if (feature_edge[mesh.edge_handle(he_it)] > 0)
      break;
  }
  while (he_it != _heh_s);

  return hehs;
}

//view from interior, the sector is ccw
template<class MeshT>
std::vector<HEH>
get_halfedges_in_feature_sector(const MeshT &mesh, const FP<int> &feature_face, const EP<int> &feature_edge,
                                const EP<int> &valence, const HEH _heh_s,
                                std::set<FH> &_visited_fhs, std::set<FH> &_sec_fhs)
{
  _sec_fhs.clear();

  std::vector<HEH> hehs;
  hehs.push_back(_heh_s);

  HEH he_it = _heh_s;
  do
  {
    HFH ft_hfh(-1);
    for (auto hehf_it = mesh.hehf_iter(he_it); hehf_it.valid(); ++hehf_it)
    {
      auto fhi = mesh.face_handle(*hehf_it);
      if (feature_face[fhi] > 0 && _visited_fhs.find(fhi) == _visited_fhs.end())
      {
        if (mesh.incident_cell(*hehf_it).is_valid())
        {
          ft_hfh = *hehf_it;
          _visited_fhs.insert(fhi);
          _sec_fhs.insert(fhi);
          break;
        }
      }
    }

    if (!ft_hfh.is_valid())
    {
      return std::vector<HEH>{};
    }

    he_it = mesh.prev_halfedge_in_halfface(he_it, ft_hfh);
    he_it = mesh.opposite_halfedge_handle(he_it);

    hehs.push_back(he_it);

    auto ehi = mesh.edge_handle(he_it);
    if (feature_edge[ehi] > 0 || valence[ehi] != 0)
      break;
  }
  while (he_it != _heh_s);

  return hehs;
}

template<class MeshT>
void
get_feature_sectors_at_feature_edge_vertex(const MeshT &mesh, const FP<int> &feature_face, const EP<int> &feature_edge,
                                           const VP<bool> &feature_edge_vertex, const EP<int> &valence, const VH vh,
                                           std::vector<std::vector<HEH>> &v_sec_hehs,
                                           std::vector<std::set<FH>> &v_sec_fhs)
{
  if (!feature_edge_vertex[vh])
    return;

  HEH ft_heh_s(-1);
  for (auto voh_it = mesh.voh_iter(vh); voh_it.valid(); ++voh_it)
  {
    if (feature_edge[mesh.edge_handle(*voh_it)] > 0)
    {
      ft_heh_s = *voh_it;
      break;
    }
  }

  if (ft_heh_s.is_valid())
  {
    std::set<FH> visited_fhs;

    HEH he_it = ft_heh_s;
    while (he_it.is_valid())
    {
      std::set<FH> sec_fhs;
      auto sec_hehs = get_halfedges_in_feature_sector(mesh, feature_face, feature_edge,
                                                      valence, he_it, visited_fhs, sec_fhs);


      if (sec_hehs.empty())
        break;

      v_sec_hehs.push_back(sec_hehs);
      v_sec_fhs.push_back(sec_fhs);

      he_it = sec_hehs.back();

      if (he_it == ft_heh_s)
        break;
    }
  }
}

template<class MeshT>
void get_k_ring_vertices(MeshT &mesh, const int k, std::set<VH> &_initial_vertices)
{
  std::set<VH> tmp;
  for (int i = 0; i < k - 1; ++i)
  {
    tmp = _initial_vertices;
    for (const auto vh: tmp)
      for (auto vv_it = mesh.vv_iter(vh); vv_it.valid(); ++vv_it)
        _initial_vertices.insert(*vv_it);
  }
}


template<class MeshT>
std::set<VH> get_k_ring_vertices_of_singular_graph(MeshT &mesh, const int k)
{
  auto sgl_vt = mesh.template request_vertex_property<int>("singular_vertex");
  std::set<VH> kring_vhs;
  for (const auto vh: mesh.vertices())
  {
    if (sgl_vt[vh] > 0)
      kring_vhs.insert(vh);
  }

  get_k_ring_vertices(mesh, k, kring_vhs);

  return kring_vhs;
}


template<class MeshT>
std::set<CH>
get_k_ring_cells(MeshT &mesh, const std::set<VH> &_kring_vhs)
{
  std::set<CH> kring_chs;

  for (const auto vh: _kring_vhs)
  {
    for (auto vc_it = mesh.vc_iter(vh); vc_it.valid(); ++vc_it)
    {
      kring_chs.insert(*vc_it);
    }
  }

  return kring_chs;
}

template<class MeshT>
std::set<CH>
get_k_ring_cells_of_singular_graph(MeshT &mesh, const int k)
{
  auto kring_vhs = get_k_ring_vertices_of_singular_graph(mesh, k);

  return get_k_ring_cells(mesh, kring_vhs);
}


template<class MeshT, typename U>
std::set<std::pair<double, EH>, U>
get_k_ring_edges_of_singular_graph(MeshT &mesh, const int k)
{
  std::set<std::pair<double, EH>, U> edge_set;
  auto kring_chs = get_k_ring_cells_of_singular_graph(mesh, k);

  for (const auto ch: kring_chs)
  {
    for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it)
      edge_set.emplace(mesh.length(*ce_it), *ce_it);
  }

  return edge_set;
}


template<class MeshT>
std::set<EH>
get_k_ring_edges_of_singular_graph(MeshT &mesh, const int k)
{
  std::set<EH> edge_set;

  auto kring_chs = get_k_ring_cells_of_singular_graph(mesh, k);

  for (const auto ch: kring_chs)
  {
    for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it)
      edge_set.insert(*ce_it);
  }

  return edge_set;
}

template<class MeshT, typename U>
std::set<std::pair<double, EH>, U>
get_k_ring_edges(MeshT &mesh, const std::vector<EH> &_sg_ehs, const int k)
{
  std::set<VH> vhs_set;
  for (auto ehi: _sg_ehs)
  {
    vhs_set.insert(mesh.edge(ehi).from_vertex());
    vhs_set.insert(mesh.edge(ehi).to_vertex());
  }
  get_k_ring_vertices(mesh, k, vhs_set);
  auto kring_chs = get_k_ring_cells(mesh, vhs_set);

  std::set<std::pair<double, EH>, U> edge_set;
  for (const auto ch: kring_chs)
  {
    for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it)
      edge_set.emplace(mesh.length(*ce_it), *ce_it);
  }

  return edge_set;
}


//stop at zipper node
template<class MeshT>
std::vector<HEH>
get_halfedges_on_singular_arc_without_zipper_node(const MeshT &mesh, const EP<int> &valence, const HEH _heh_s)
{
  std::vector<HEH> hehs2, hehs1;

  if (valence[mesh.edge_handle(_heh_s)] == 0)
    return hehs2;

  hehs2.push_back(_heh_s);

  HEH heh_it = _heh_s;
  auto eh_s = mesh.edge_handle(_heh_s);
  int val_cur = valence[eh_s];
  bool is_bdy = mesh.is_boundary(eh_s);

  while (heh_it.is_valid())
  {
    VH vh_to = mesh.halfedge(heh_it).to_vertex();
    //if it's an interior arc, stop if touches boundary
    if (!is_bdy && mesh.is_boundary(vh_to))
      break;

    int n = 0;
    HEH heh_next;
    auto eh_cur = mesh.edge_handle(heh_it);
    for (auto voh_it = mesh.voh_iter(vh_to); voh_it.valid(); ++voh_it)
    {
      auto eh = mesh.edge_handle(*voh_it);
      if (valence[eh] != 0 && eh != eh_cur)
      {
        n++;
        if (valence[eh] == val_cur)
          heh_next = *voh_it;
      }
    }

    if (n == 1 && heh_next.is_valid())
    {
      heh_it = heh_next;
      hehs2.push_back(heh_it);
    }
    else
      heh_it = HEH(-1);
  }


  heh_it = mesh.opposite_halfedge_handle(_heh_s);
  while (heh_it.is_valid())
  {
    VH vh_to = mesh.halfedge(heh_it).to_vertex();
    if (!is_bdy && mesh.is_boundary(vh_to))
      break;

    int n = 0;
    HEH heh_next;
    auto eh_cur = mesh.edge_handle(heh_it);
    for (auto voh_it = mesh.voh_iter(vh_to); voh_it.valid(); ++voh_it)
    {
      auto eh = mesh.edge_handle(*voh_it);
      if (valence[eh] != 0 && eh != eh_cur)
      {
        n++;
        if (valence[eh] == val_cur)
          heh_next = *voh_it;
      }
    }

    if (n == 1 && heh_next.is_valid())
    {
      heh_it = heh_next;
      hehs1.push_back(heh_it);
    }
    else
      heh_it = HEH(-1);
  }

  if (!hehs1.empty())
  {
    for (auto i = 0u; i < hehs1.size(); ++i)
      hehs1[i] = mesh.opposite_halfedge_handle(hehs1[i]);

    std::reverse(hehs1.begin(), hehs1.end());
    hehs1.insert(hehs1.end(), hehs2.begin(), hehs2.end());

    return hehs1;
  }

  return hehs2;
}

//singular edges on feature surface
template<class MeshT>
std::vector<HEH>
get_halfedges_of_singular_arc_on_feature_face(const MeshT &mesh, const EP<int> &valence,
                                              const EP<bool> &feature_face_edge, const HEH _heh_s,
                                              const bool _same_valence)
{
  std::vector<HEH> hehs2, hehs1;

  auto eh_s = mesh.edge_handle(_heh_s);
  if (valence[eh_s] == 0 || !feature_face_edge[eh_s])
    return hehs2;

  hehs2.push_back(_heh_s);

  std::set<HEH> visited_hehs;
  visited_hehs.insert(_heh_s);

  HEH heh_it = _heh_s;
  while (heh_it.is_valid())
  {
    VH vh_to = mesh.halfedge(heh_it).to_vertex();

    int n = 0;
    HEH heh_next;
    auto eh_cur = mesh.edge_handle(heh_it);
    for (auto voh_it = mesh.voh_iter(vh_to); voh_it.valid(); ++voh_it)
    {
      auto eh_next = mesh.edge_handle(*voh_it);
      if (feature_face_edge[eh_next] && valence[eh_next] != 0 && eh_next != eh_cur)
      {
        n++;

        if (!_same_valence || (_same_valence && valence[eh_next] == valence[eh_cur]))
        {
          heh_next = *voh_it;
        }
      }
    }

    if (n == 1 && heh_next.is_valid())
    {
      if (visited_hehs.count(heh_next) > 0)
        break;

      heh_it = heh_next;
      hehs2.push_back(heh_it);
      visited_hehs.insert(heh_it);
    }
    else
      heh_it = HEH(-1);
  }


  heh_it = mesh.opposite_halfedge_handle(_heh_s);
  while (heh_it.is_valid())
  {
    VH vh_to = mesh.halfedge(heh_it).to_vertex();

    int n = 0;
    HEH heh_next;
    auto eh_cur = mesh.edge_handle(heh_it);
    for (auto voh_it = mesh.voh_iter(vh_to); voh_it.valid(); ++voh_it)
    {
      auto eh_next = mesh.edge_handle(*voh_it);
      if (feature_face_edge[eh_next] && valence[eh_next] != 0 && eh_next != eh_cur)
      {
        n++;
        if (!_same_valence || (_same_valence && valence[eh_next] == valence[eh_cur]))
        {
          heh_next = *voh_it;
        }
      }
    }

    if (n == 1 && heh_next.is_valid())
    {
      if (visited_hehs.count(heh_next) > 0)
        break;

      heh_it = heh_next;
      hehs1.push_back(heh_it);
      visited_hehs.insert(heh_it);
    }
    else
      heh_it = HEH(-1);
  }

  if (!hehs1.empty())
  {
    for (auto i = 0u; i < hehs1.size(); ++i)
      hehs1[i] = mesh.opposite_halfedge_handle(hehs1[i]);

    std::reverse(hehs1.begin(), hehs1.end());
    hehs1.insert(hehs1.end(), hehs2.begin(), hehs2.end());

    return hehs1;
  }

  return hehs2;
}


template<class MeshT>
std::vector<HEH>
get_halfedges_of_singular_arc_kept_on_feature_face(const MeshT &mesh, const EP<int> &valence,
                                                   const EP<bool> &keep_on_face, const EP<int> &feature_edge,
                                                   const HEH _heh_s)
{
  std::vector<HEH> hehs2, hehs1;

  auto eh_s = mesh.edge_handle(_heh_s);
  if (valence[eh_s] == 0 || !keep_on_face[eh_s])
    return hehs2;

  hehs2.push_back(_heh_s);

  std::set<HEH> visited_hehs;
  visited_hehs.insert(_heh_s);

  HEH heh_it = _heh_s;
  while (heh_it.is_valid())
  {
    VH vh_to = mesh.halfedge(heh_it).to_vertex();

    int n = 0;
    HEH heh_next;
    auto eh_cur = mesh.edge_handle(heh_it);
    for (auto voh_it = mesh.voh_iter(vh_to); voh_it.valid(); ++voh_it)
    {
      auto eh_next = mesh.edge_handle(*voh_it);
      if (keep_on_face[eh_next] && valence[eh_next] != 0 && feature_edge[eh_next] == 0 && eh_next != eh_cur)
      {
        n++;
        heh_next = *voh_it;
      }
    }

    if (n == 1)
    {
      if (visited_hehs.count(heh_next) > 0)
        break;

      heh_it = heh_next;
      hehs2.push_back(heh_it);
      visited_hehs.insert(heh_it);
    }
    else
      heh_it = HEH(-1);
  }


  heh_it = mesh.opposite_halfedge_handle(_heh_s);
  while (heh_it.is_valid())
  {
    VH vh_to = mesh.halfedge(heh_it).to_vertex();

    int n = 0;
    HEH heh_next;
    auto eh_cur = mesh.edge_handle(heh_it);
    for (auto voh_it = mesh.voh_iter(vh_to); voh_it.valid(); ++voh_it)
    {
      auto eh_next = mesh.edge_handle(*voh_it);
      if (keep_on_face[eh_next] && valence[eh_next] != 0 && feature_edge[eh_next] == 0 && eh_next != eh_cur)
      {
        n++;
        heh_next = *voh_it;
      }
    }

    if (n == 1)
    {
      if (visited_hehs.count(heh_next) > 0)
        break;

      heh_it = heh_next;
      hehs1.push_back(heh_it);
      visited_hehs.insert(heh_it);
    }
    else
      heh_it = HEH(-1);
  }

  if (!hehs1.empty())
  {
    for (auto i = 0u; i < hehs1.size(); ++i)
      hehs1[i] = mesh.opposite_halfedge_handle(hehs1[i]);

    std::reverse(hehs1.begin(), hehs1.end());
    hehs1.insert(hehs1.end(), hehs2.begin(), hehs2.end());

    return hehs1;
  }

  return hehs2;
}

template<class MeshT>
std::vector<HEH> get_halfedges_on_singular_arc_directional(const MeshT &mesh, const EP<int> &valence, const HEH heh)
{
  std::vector<HEH> hehs;
  auto sg_vhs = get_vertices_on_singular_arc_directional(mesh, valence, heh);

  if (sg_vhs.empty())
    return hehs;

  for (size_t i = 0; i < sg_vhs.size() - 1; ++i)
    hehs.push_back(mesh.find_halfedge(sg_vhs[i], sg_vhs[i + 1]));

  auto he = mesh.find_halfedge(sg_vhs.back(), sg_vhs[0]);
  if (he.is_valid())
  {
    if (sg_vhs.size() > 2 && valence[mesh.edge_handle(he)] != 0)
    {
      int n_sge_b = n_incident_singular_edges(mesh, valence, sg_vhs.back());
      int n_sge_s = n_incident_singular_edges(mesh, valence, sg_vhs[0]);
      if (!(n_sge_b > 2 && n_sge_s > 2))
        hehs.push_back(he);
    }
  }

  return hehs;
}

//has to be singular edge too
template<class MeshT>
std::vector<HEH>
get_halfedges_on_feature_arc_directional(const MeshT &mesh, const EP<int> &valence, const EP<int> &feature_edge,
                                         const HEH heh)
{
  std::vector<HEH> hehs;

  if (!heh.is_valid())
    return hehs;

  auto eh_s = mesh.edge_handle(heh);

  if (feature_edge[eh_s] == 0)
    return hehs;

  hehs.push_back(heh);

  std::set<HEH> visited_hehs;
  visited_hehs.insert(heh);

  HEH heh_it = heh;

  while (heh_it.is_valid())
  {
    VH vh_to = mesh.halfedge(heh_it).to_vertex();

    int n = 0;
    HEH heh_next;
    auto eh_cur = mesh.edge_handle(heh_it);
    for (auto voh_it = mesh.voh_iter(vh_to); voh_it.valid(); ++voh_it)
    {
      auto eh_next = mesh.edge_handle(*voh_it);
      if (feature_edge[eh_next] > 0 && eh_next != eh_cur)
      {
        if (valence[eh_next] == 0)
          heh_next = *voh_it;
        n++;
      }
    }

    if (n == 1 && heh_next.is_valid())
    {
      if (visited_hehs.count(heh_next) > 0)
        break;
      heh_it = heh_next;
      hehs.push_back(heh_it);
      visited_hehs.insert(heh_it);
    }
    else
      heh_it = HEH(-1);
  }


  return hehs;
}

template<class MeshT>
std::vector<HEH> get_halfedges_on_feature_arc_directional(const MeshT &mesh, const EP<int> &feature_edge, const HEH heh)
{
  std::vector<HEH> hehs;

  if (!heh.is_valid())
    return hehs;

  auto eh_s = mesh.edge_handle(heh);

  if (feature_edge[eh_s] == 0)
    return hehs;

  hehs.push_back(heh);

  std::set<HEH> visited_hehs;
  visited_hehs.insert(heh);

  HEH heh_it = heh;

  while (heh_it.is_valid())
  {
    VH vh_to = mesh.halfedge(heh_it).to_vertex();

    int n = 0;
    HEH heh_next;
    auto eh_cur = mesh.edge_handle(heh_it);
    for (auto voh_it = mesh.voh_iter(vh_to); voh_it.valid(); ++voh_it)
    {
      auto eh_next = mesh.edge_handle(*voh_it);
      if (feature_edge[eh_next] > 0 && eh_next != eh_cur)
      {
        heh_next = *voh_it;
        n++;
      }
    }

    if (n == 1 && heh_next.is_valid())
    {
      if (visited_hehs.count(heh_next) > 0)
        break;
      heh_it = heh_next;
      hehs.push_back(heh_it);
      visited_hehs.insert(heh_it);
    }
    else
      heh_it = HEH(-1);
  }


  return hehs;
}


//keep tracing at zipper node
template<class MeshT>
std::vector<VH>
get_vertices_on_singular_arc(const MeshT &mesh, const EP<int> &valence, const HEH heh)
{
  std::vector<VH> sg_vhs;

  if (!heh.is_valid())
  {
    std::cerr << "Error: input edge is not valid!" << std::endl;
    return std::vector<VH>();
  }
  if (valence[mesh.edge_handle(heh)] == 0)
  {
    std::cerr << "Warning: input edge " << mesh.edge_handle(heh) << " is not singular!" << std::endl;
    return std::vector<VH>();
  }

  VH vh_f = mesh.halfedge(heh).from_vertex();
  VH vh_t = mesh.halfedge(heh).to_vertex();

  std::vector<VH> vhs0, vhs1;
  vhs1.push_back(vh_t);
  vhs0.push_back(vh_f);

  std::set<VH> sgvhs_set;
  sgvhs_set.insert(vh_t);
  sgvhs_set.insert(vh_f);

  std::queue<VH> que_vhs;

  que_vhs.push(vh_f);
  while (!que_vhs.empty())
  {
    VH vh_cur = que_vhs.front();
    que_vhs.pop();


    VH vh_next;
    int n_sg = 0;
    for (auto voh_it = mesh.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
    {
      EH eh_og = mesh.edge_handle(*voh_it);
      VH vh_tmp = mesh.halfedge(*voh_it).to_vertex();
      int val = valence[eh_og];
      if (val != 0 && sgvhs_set.count(vh_tmp) == 0)
      {
        n_sg++;

        if (val != 0)
          vh_next = vh_tmp;
      }
    }

    if (n_sg == 1 && vh_next.is_valid())
    {
      sgvhs_set.insert(vh_next);
      vhs0.push_back(vh_next);

      que_vhs.push(vh_next);
    }
  }

  que_vhs.push(vh_t);

  while (!que_vhs.empty())
  {
    VH vh_cur = que_vhs.front();
    que_vhs.pop();

    VH vh_next;
    int n_sg = 0;
    for (auto voh_it = mesh.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
    {
      EH eh_og = mesh.edge_handle(*voh_it);
      VH vh_tmp = mesh.halfedge(*voh_it).to_vertex();
      int val = valence[eh_og];
      if (val != 0 && sgvhs_set.count(vh_tmp) == 0)
      {
        n_sg++;

        if (val != 0)
          vh_next = vh_tmp;
      }
    }

    if (n_sg == 1 && vh_next.is_valid())
    {
      sgvhs_set.insert(vh_next);
      vhs1.push_back(vh_next);

      que_vhs.push(vh_next);
    }
  }


  sg_vhs.reserve(vhs0.size() + vhs1.size());
  sg_vhs.insert(sg_vhs.end(), vhs0.rbegin(), vhs0.rend());
  sg_vhs.insert(sg_vhs.end(), vhs1.begin(), vhs1.end());

  return sg_vhs;
}

//keep tracing at zipper node
template<class MeshT>
std::vector<VH>
get_vertices_on_singular_arc_with_intertior_surface(const MeshT &mesh, const EP<int> &valence,
                                                    const VP<bool> &feature_face_vertex, const HEH heh)
{
  std::vector<VH> sg_vhs;

  if (!heh.is_valid())
  {
    std::cerr << "Error: input edge is not valid!" << std::endl;
    return std::vector<VH>();
  }
  if (valence[mesh.edge_handle(heh)] == 0)
  {
    std::cerr << "Warning: input edge " << mesh.edge_handle(heh) << " is not singular!" << std::endl;
    return std::vector<VH>();
  }

  VH vh_f = mesh.halfedge(heh).from_vertex();
  VH vh_t = mesh.halfedge(heh).to_vertex();

  std::vector<VH> vhs0, vhs1;
  vhs1.push_back(vh_t);
  vhs0.push_back(vh_f);

  std::set<VH> sgvhs_set;
  sgvhs_set.insert(vh_t);
  sgvhs_set.insert(vh_f);

  std::queue<VH> que_vhs;

  que_vhs.push(vh_f);
  while (!que_vhs.empty())
  {
    VH vh_cur = que_vhs.front();
    que_vhs.pop();

    if (feature_face_vertex[vh_cur])
      break;

    VH vh_next;
    int n_sg = 0;
    for (auto voh_it = mesh.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
    {
      EH eh_og = mesh.edge_handle(*voh_it);
      VH vh_tmp = mesh.halfedge(*voh_it).to_vertex();
      int val = valence[eh_og];
      if (val != 0 && sgvhs_set.count(vh_tmp) == 0)
      {
        n_sg++;

        if (val != 0)
          vh_next = vh_tmp;
      }
    }

    if (n_sg == 1 && vh_next.is_valid())
    {
      sgvhs_set.insert(vh_next);
      vhs0.push_back(vh_next);

      que_vhs.push(vh_next);
    }
  }

  que_vhs.push(vh_t);

  while (!que_vhs.empty())
  {
    VH vh_cur = que_vhs.front();
    que_vhs.pop();

    if (feature_face_vertex[vh_cur])
      break;

    VH vh_next;
    int n_sg = 0;
    for (auto voh_it = mesh.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
    {
      EH eh_og = mesh.edge_handle(*voh_it);
      VH vh_tmp = mesh.halfedge(*voh_it).to_vertex();
      int val = valence[eh_og];
      if (val != 0 && sgvhs_set.count(vh_tmp) == 0)
      {
        n_sg++;

        if (val != 0)
          vh_next = vh_tmp;
      }
    }

    if (n_sg == 1 && vh_next.is_valid())
    {
      sgvhs_set.insert(vh_next);
      vhs1.push_back(vh_next);

      que_vhs.push(vh_next);
    }
  }


  sg_vhs.reserve(vhs0.size() + vhs1.size());
  sg_vhs.insert(sg_vhs.end(), vhs0.rbegin(), vhs0.rend());
  sg_vhs.insert(sg_vhs.end(), vhs1.begin(), vhs1.end());

  return sg_vhs;
}

template<class MeshT>
std::vector<VH> get_vertices_on_singular_arc_directional(const MeshT &mesh, const EP<int> &valence, const HEH heh)
{
  std::vector<VH> sg_vhs;

  if (!heh.is_valid())
  {
    std::cerr << "Error: input edge is not valid!" << std::endl;
    return std::vector<VH>();
  }
  if (valence[mesh.edge_handle(heh)] == 0)
  {
    std::cerr << "Warning: input edge " << mesh.edge_handle(heh) << " is not singular!" << std::endl;
    return std::vector<VH>();
  }

  VH vh_f = mesh.halfedge(heh).from_vertex();
  VH vh_t = mesh.halfedge(heh).to_vertex();

  sg_vhs.push_back(vh_f);
  sg_vhs.push_back(vh_t);

  std::set<VH> sgvhs_set;
  sgvhs_set.insert(vh_t);
  sgvhs_set.insert(vh_f);

  std::queue<VH> que_vhs;
  que_vhs.push(vh_t);

  while (!que_vhs.empty())
  {
    VH vh_cur = que_vhs.front();
    que_vhs.pop();

    VH vh_next;
    int n_sg = 0;
    for (auto voh_it = mesh.voh_iter(vh_cur); voh_it.valid(); ++voh_it)
    {
      EH eh_og = mesh.edge_handle(*voh_it);

      VH vh_tmp = mesh.halfedge(*voh_it).to_vertex();
      int val = valence[eh_og];
      if (val != 0 && sgvhs_set.count(vh_tmp) == 0)
      {
        n_sg++;

        if (val != 0)
          vh_next = vh_tmp;
      }
    }

    if (n_sg == 1 && vh_next.is_valid())
    {
      sgvhs_set.insert(vh_next);
      sg_vhs.push_back(vh_next);

      que_vhs.push(vh_next);
    }
  }


  return sg_vhs;
}

template<class MeshT>
std::vector<HEH>
get_halfedges_on_singular_arc(const MeshT &mesh, const EP<int> &valence, const HEH heh)
{
  std::vector<HEH> hehs;
  auto sg_vhs = get_vertices_on_singular_arc(mesh, valence, heh);

  if (sg_vhs.empty())
    return hehs;

  for (size_t i = 0; i < sg_vhs.size() - 1; ++i)
    hehs.push_back(mesh.find_halfedge(sg_vhs[i], sg_vhs[i + 1]));

  auto he = mesh.find_halfedge(sg_vhs.back(), sg_vhs[0]);
  if (he.is_valid())
  {
    if (sg_vhs.size() > 2 && valence[mesh.edge_handle(he)] != 0)
    {
      int n_sge_b = n_incident_singular_edges(mesh, valence, sg_vhs.back());
      int n_sge_s = n_incident_singular_edges(mesh, valence, sg_vhs[0]);
      if (!(n_sge_b > 2 && n_sge_s > 2))
        hehs.push_back(he);
    }
  }

  return hehs;
}

template<class MeshT>
std::vector<HEH> get_halfedges_on_singular_arc_with_intertior_surface(const MeshT &mesh, const EP<int> &valence,
                                                                      const VP<bool> &feature_face_vertex,
                                                                      const HEH heh)
{
  std::vector<HEH> hehs;
  auto sg_vhs = get_vertices_on_singular_arc_with_intertior_surface(mesh, valence, feature_face_vertex, heh);

  if (sg_vhs.empty())
    return hehs;

  for (size_t i = 0; i < sg_vhs.size() - 1; ++i)
    hehs.push_back(mesh.find_halfedge(sg_vhs[i], sg_vhs[i + 1]));

  auto he = mesh.find_halfedge(sg_vhs.back(), sg_vhs[0]);
  if (he.is_valid())
  {
    if (sg_vhs.size() > 2 && valence[mesh.edge_handle(he)] != 0)
    {
      int n_sge_b = n_incident_singular_edges(mesh, valence, sg_vhs.back());
      int n_sge_s = n_incident_singular_edges(mesh, valence, sg_vhs[0]);
      if (!(n_sge_b > 2 && n_sge_s > 2))
        hehs.push_back(he);
    }
  }

  return hehs;
}


template<class MeshT>
std::vector<EH>
get_edges_on_singular_arc(const MeshT &mesh, const EP<int> &valence, const HEH heh)
{
  std::vector<EH> ehs;
  auto sg_hehs = get_halfedges_on_singular_arc(mesh, valence, heh);

  if (sg_hehs.empty())
    return ehs;

  ehs.reserve(sg_hehs.size());
  for (size_t i = 0; i < sg_hehs.size(); ++i)
    ehs.push_back(mesh.edge_handle(sg_hehs[i]));

  return ehs;
}

template<class MeshT>
bool is_singular_node(const MeshT &mesh, const EP<int> &valence, const VH vh)
{
  int n = 0;
  for (auto ve_it = mesh.ve_iter(vh); ve_it.valid(); ++ve_it)
    if (valence[*ve_it] != 0)
      n++;

  if (n > 2)
    return true;

  return false;
}

template<class MeshT>
std::set<CH> get_onering_cells(const MeshT &mesh, const VH vh)
{
  std::set<CH> onering_chs;
  for (auto vc_it = mesh.vc_iter(vh); vc_it.valid(); ++vc_it)
    onering_chs.insert(*vc_it);

  return onering_chs;
}

template<class MeshT>
std::set<CH>
get_part_of_onering_cells(const MeshT &mesh, const FP<int> &feature_fprop, const VH vh, const std::set<CH> &onering_chs,
                          const CH ch_s, std::set<HFH> &ft_hfhs)
{
  ft_hfhs.clear();

  std::set<CH> por_chs;

  std::queue<CH> ch_que;
  ch_que.push(ch_s);
  por_chs.insert(ch_s);

  while (!ch_que.empty())
  {
    auto ch_cur = ch_que.front();
    ch_que.pop();

    for (auto chf_it = mesh.chf_iter(ch_cur); chf_it.valid(); ++chf_it)
    {
      auto cf = mesh.face_handle(*chf_it);
      if (feature_fprop[cf] > 0)
      {
        ft_hfhs.insert(*chf_it);
        continue;
      }
      else
      {
        CH ch_next = mesh.incident_cell(mesh.opposite_halfface_handle(*chf_it));
        if (ch_next.is_valid() && onering_chs.find(ch_next) != onering_chs.end() &&
            por_chs.find(ch_next) == por_chs.end())
        {
          por_chs.insert(ch_next);
          ch_que.push(ch_next);
        }
      }
    }
  }


  return por_chs;
}

template<class MeshT>
std::map<CH, int>
get_region_id_of_onering_cells(const MeshT &mesh, const FP<int> &feature_fprop, const VH vh)
{
  std::map<CH, int> cell_region_id;
  auto or_chs = get_onering_cells(mesh, vh);
  auto all_chs = or_chs;

  int rg_id = 0;
  while(!all_chs.empty()) {
    auto ch_cur_it = all_chs.begin();

    std::set<HFH> ft_hfhs;
    auto por_chs = get_part_of_onering_cells(mesh, feature_fprop, vh, or_chs, *ch_cur_it, ft_hfhs);

    for(const auto chi : por_chs)
      cell_region_id.insert(std::make_pair(chi, rg_id));

    //erase
    for(const auto chi : por_chs)
      all_chs.erase(chi);

    rg_id++;
  }

  return cell_region_id;
}

}

/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once


#include <assert.h>
#include <AlgoHex/TypeDef.hh>
#include <stack>
#include "AlgoHex/BinarySpacePartitionTrees/TriangleBSPT.hh"


namespace AlgoHex
{
//    template <typename Point>
//    Eigen::Vector3d to_eigen_vec3(const Point& _vec) {return {_vec[0], _vec[1], _vec[2]}; }

template<class MeshT>
typename MeshT::PointT face_normal(const MeshT &_mesh, const HFH _hfh);

template<class MeshT>
typename MeshT::PointT face_normal(const MeshT &_mesh, const VH _vh0, const VH _vh1, const VH _vh2);

template<class MeshT>
typename MeshT::PointT face_normal(const MeshT &_mesh, const std::vector<VH> &_vhs);

template<class MeshT>
typename MeshT::PointT vertex_normal_area_weighted(const MeshT &_mesh, const VH _vh);

template<class MeshT>
typename MeshT::PointT vertex_normal_uniform_weighted(const MeshT &_mesh, const VH _vh);

template<typename PointT>
PointT vertex_normal_uniform_weighted(const PointT &_pt0, const PointT &_pt1, const PointT &_n0, const PointT &_n1);

template<typename Vec3>
std::array<double, 3> compute_bary_coord(const Vec3 &_px, const Vec3 &_p0, const Vec3 &_p1, const Vec3 &_p2);

template<class MeshT>
double edge_angle(const MeshT &mesh, const HEH heh0, const HEH heh1);

template<class MeshT>
void update_cells_volume(MeshT &mesh, CP<double> &cell_volume, double &total_volume);

template<class MeshT>
AlgoHex::FH common_face(const MeshT &mesh, const CH _ch0, const CH _ch1);

template<class MeshT>
void build_surface_mesh(MeshT &mesh, MeshT &sf_mesh);

template<class MeshT>
std::unique_ptr<AlgoHex::OpenVolumeMeshTriangleBSPT<MeshT>> build_bsp(const MeshT &sf_mesh);

template<class MeshT>
bool has_interior_singular_edge(const MeshT &mesh, const EP<int> &valence, const VH vh);

template<class MeshT>
int n_interior_complex_singular_edge(const MeshT &mesh, const EP<int> &valence, const VH vh);

template<class MeshT>
int n_incident_singular_edges_of_valence_one(const MeshT &mesh, const EP<int> &valence, const VH _vh);

template<class MeshT>
int n_incident_special_edges(const VH _vh);

template<class MeshT>
int n_incident_special_edges_on_feature_face(const MeshT &mesh, const EP<bool> &feature_face_edge,
                                             const EP<int> &feature_edge, const EP<int> &valence, const VH vh);

template<class MeshT>
int n_incident_special_edges_kept_on_feature_face(const MeshT &mesh, const EP<bool> &feature_face_edge,
                                                  const EP<int> &feature_edge, const EP<int> &valence,
                                                  const EP<bool> &keep_on_surface, const VH vh);

template<class MeshT>
int n_incident_feature_edges(const MeshT &mesh, const EP<int> &feature_edge, const VH _vh);

template<class MeshT>
int
n_incident_feature_singular_edges(const MeshT &mesh, const EP<int> &feature_edge, const EP<int> &valence, const VH vh);

template<class MeshT>
void
get_special_halfedges_at_vertex(const MeshT &mesh, const EP<int> &feature_edge, const EP<int> &valence, const VH vh,
                                std::vector<HEH> &ft_hehs, std::vector<HEH> &bdy_hehs, std::vector<HEH> &itr_hehs);

template<class MeshT>
std::vector<HEH> get_singular_halfedges(const MeshT &mesh, const EP<int> &valence, const VH _vh);

template<class MeshT>
std::vector<HEH>
get_non_feature_singular_halfedges(const MeshT &mesh, const EP<bool> &feature_feace_edge, const EP<int> &valence,
                                   const VH _vh);

template<class MeshT>
void
count_special_edges_at_vertex(const MeshT &mesh, const EP<bool> &feature_face_edge, const EP<int> &feature_edge,
                              const EP<int> &valence, const VH vh,
                              int &n_nff_val_1, int &n_nff_val1, int &n_nff_val_u, int &n_nff, int &n_bdy_val_1,
                              int &n_bdy_val1, int &n_bdy_val_u, int &n_bdy, int &n_cmpl,
                              int &n_ff_val_1, int &n_ff_val1, int &n_ff_val_u, int &n_ff, int &n_fe);

template<class MeshT>
std::set<EH>
incident_non_ffe_singular_edges_in_same_region(const MeshT &mesh, const EP<bool> &feature_face_edge,
                                               const EP<int> &valence,
                                               const VH vh, const std::set<CH> &por_chs);

template<class MeshT>
int
n_incident_non_ffe_singular_edges_in_same_region(const MeshT &mesh, const EP<bool> &feature_face_edge,
                                                 const EP<int> &valence,
                                                 const VH vh, const std::set<CH> &por_chs);

template<class MeshT>
std::vector<HEH>
get_halfedges_in_feature_sector(const MeshT &mesh, const FP<int> &feature_face, const EP<int> &feature_edge,
                                const EP<int> &valence, const HEH _heh_s,
                                std::set<FH> &_visited_fhs, std::set<FH> &_sec_fhs);

//sectors of feature edges or singular edge on feature surface
template<class MeshT>
void
get_feature_sectors_at_feature_edge_vertex(const MeshT &mesh, const FP<int> &feature_face, const EP<int> &feature_edge,
                                           const VP<bool> &feature_edge_vertex,
                                           const EP<int> &valence, const VH vh,
                                           std::vector<std::vector<HEH>> &v_sec_hehs,
                                           std::vector<std::set<FH>> &v_sec_fhs);

template<class MeshT>
std::vector<HEH>
get_boundary_halfedges_in_feature_sector_ccw(const MeshT &mesh, const EP<int> &feature_edge, const EP<int> &valence,
                                             const HEH _heh_s);

template<class MeshT>
std::vector<HEH>
get_boundary_halfedges_in_feature_sector_cw(const MeshT &mesh, const EP<int> &feature_edge, const HEH _heh_s);

template<class MeshT>
void get_k_ring_vertices(MeshT &mesh, const int k, std::set<VH> &_initial_vertices);

template<class MeshT>
std::set<VH> get_k_ring_vertices_of_singular_graph(MeshT &mesh, const int k);

template<class MeshT>
std::set<CH> get_k_ring_cells(MeshT &mesh, const std::set<VH> &_kring_vhs);

template<class MeshT>
std::set<CH> get_k_ring_cells_of_singular_graph(MeshT &mesh, const int k);

template<class MeshT, typename U>
std::set<std::pair<double, EH>, U> get_k_ring_edges_of_singular_graph(MeshT &mesh, const int k);

template<class MeshT>
std::set<EH> get_k_ring_edges_of_singular_graph(MeshT &mesh, const int k);

template<class MeshT, typename U>
std::set<std::pair<double, EH>, U> get_k_ring_edges(MeshT &mesh, const std::vector<EH> &_sg_ehs, const int k);

//same valence
template<class MeshT>
std::vector<HEH>
get_halfedges_on_singular_arc_without_zipper_node(const MeshT &mesh, const EP<int> &valence, const HEH _heh_s);

//ignore interior singular
template<class MeshT>
std::vector<HEH>
get_halfedges_of_singular_arc_on_feature_face(const MeshT &mesh, const EP<int> &valence,
                                              const EP<bool> &feature_face_edge, const HEH _heh_s,
                                              const bool _same_valence = false);

template<class MeshT>
std::vector<HEH>
get_halfedges_of_singular_arc_kept_on_feature_face(const MeshT &mesh, const EP<int> &valence,
                                                   const EP<bool> &keep_on_face, const EP<int> &feature_edge,
                                                   const HEH _heh_s);

//singular arc that ignore valence
template<class MeshT>
std::vector<VH> get_vertices_on_singular_arc(const MeshT &mesh, const EP<int> &valence, const HEH heh);

template<class MeshT>
std::vector<VH>
get_vertices_on_singular_arc_with_intertior_surface(const MeshT &mesh, const EP<int> &valence,
                                                    const VP<bool> &feature_face_vertex, const HEH heh);

//singular vertices that are on the outgoing singular arc
template<class MeshT>
std::vector<VH> get_vertices_on_singular_arc_directional(const MeshT &mesh, const EP<int> &valence, const HEH heh);

//singular arc that ignore valence
template<class MeshT>
std::vector<HEH> get_halfedges_on_singular_arc(const MeshT &mesh, const EP<int> &valence, const HEH heh);

//singular arc that ignore valence and ends at feature surface
template<class MeshT>
std::vector<HEH> get_halfedges_on_singular_arc_with_intertior_surface(const MeshT &mesh, const EP<int> &valence,
                                                                      const VP<bool> &feature_face_vertex,
                                                                      const HEH heh);

//singular arc that ignore valence
template<class MeshT>
std::vector<EH> get_edges_on_singular_arc(const MeshT &mesh, const EP<int> &valence, const HEH heh);

//halfedges that are on the outgoing singular arc
template<class MeshT>
std::vector<HEH> get_halfedges_on_singular_arc_directional(const MeshT &mesh, const EP<int> &valence, const HEH heh);

//halfedges that are on the outgoing feature arc
template<class MeshT>
std::vector<HEH>
get_halfedges_on_feature_arc_directional(const MeshT &mesh, const EP<int> &valence, const EP<int> &feature_edge,
                                         const HEH heh);

//no distinguish of valence
template<class MeshT>
std::vector<HEH>
get_halfedges_on_feature_arc_directional(const MeshT &mesh, const EP<int> &feature_edge, const HEH heh);

template<class MeshT>
std::set<CH> get_onering_cells(const MeshT &mesh, const VH _vh);

template<class MeshT>
std::set<CH> get_part_of_onering_cells(const MeshT &mesh, const FP<bool> &feature_fprop, const VH vh,
                                       const std::set<CH> &onering_chs, const CH ch_s, std::set<HFH> &ft_hfhs);

template<class MeshT>
std::map<CH, int>
get_region_id_of_onering_cells(const MeshT &mesh, const FP<int> &feature_fprop, const VH vh);
}

//=============================================================================
#if defined(INCLUDE_TEMPLATES) && !defined(COMMONFUNCST_C)
#define COMMONFUNCST_TEMPLATES

#include "CommonFuncs_impl.hh"

#endif
//=============================================================================
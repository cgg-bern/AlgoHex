/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/
#pragma once


#include <AlgoHex/TypeDef.hh>
#include "AlgoHex/BinarySpacePartitionTrees/TriangleBSPT.hh"


namespace AlgoHex
{

TetMesh::PointT ALGOHEX_EXPORT face_normal(const TetMesh &_mesh, const HFH _hfh);

TetMesh::PointT ALGOHEX_EXPORT face_normal(const TetMesh &_mesh, const VH _vh0, const VH _vh1, const VH _vh2);

TetMesh::PointT ALGOHEX_EXPORT face_normal(const TetMesh &_mesh, const std::vector<VH> &_vhs);

TetMesh::PointT ALGOHEX_EXPORT vertex_normal_area_weighted(const TetMesh &_mesh, const VH _vh);

TetMesh::PointT ALGOHEX_EXPORT vertex_normal_uniform_weighted(const TetMesh &_mesh, const VH _vh);

std::array<double, 3> ALGOHEX_EXPORT compute_bary_coord(const Point &_px, const Point &_p0, const Point &_p1, const Point &_p2);

double ALGOHEX_EXPORT edge_angle(const TetMesh &mesh, const HEH heh0, const HEH heh1);

void ALGOHEX_EXPORT update_cells_volume(TetMesh &mesh, CP<double> &cell_volume, double &total_volume);

AlgoHex::FH ALGOHEX_EXPORT common_face(const TetMesh &mesh, const CH _ch0, const CH _ch1);

void ALGOHEX_EXPORT build_surface_mesh(TetMesh &mesh, TetMesh &sf_mesh);

std::unique_ptr<AlgoHex::OpenVolumeMeshTriangleBSPT<TetMesh>> ALGOHEX_EXPORT build_bsp(const TetMesh &sf_mesh);

bool ALGOHEX_EXPORT has_interior_singular_edge(const TetMesh &mesh, const EP<int> &valence, const VH vh);

int ALGOHEX_EXPORT n_interior_complex_singular_edge(const TetMesh &mesh, const EP<int> &valence, const VH vh);

int ALGOHEX_EXPORT n_incident_singular_edges_of_valence_one(const TetMesh &mesh, const EP<int> &valence, const VH _vh);

int ALGOHEX_EXPORT n_incident_singular_edges(const TetMesh &mesh, const EP<int> &valence, const VH _vh);

int ALGOHEX_EXPORT n_invalid_singular_edge(const TetMesh &mesh, const EP<int> &valence, const VH vh);

bool ALGOHEX_EXPORT is_singular_node(const TetMesh &mesh, const EP<int> &valence, const VH vh);

int ALGOHEX_EXPORT n_incident_special_edges(const VH _vh);

int ALGOHEX_EXPORT n_incident_special_edges_on_feature_face(const TetMesh &mesh, const EP<bool> &feature_face_edge,
                                             const EP<int> &feature_edge, const EP<int> &valence, const VH vh);

int ALGOHEX_EXPORT n_incident_special_edges_kept_on_feature_face(const TetMesh &mesh, const EP<bool> &feature_face_edge,
                                                  const EP<int> &feature_edge, const EP<int> &valence,
                                                  const EP<bool> &keep_on_surface, const VH vh);

int ALGOHEX_EXPORT n_incident_feature_edges(const TetMesh &mesh, const EP<int> &feature_edge, const VH _vh);

int
ALGOHEX_EXPORT n_incident_feature_singular_edges(const TetMesh &mesh, const EP<int> &feature_edge, const EP<int> &valence, const VH vh);

void
ALGOHEX_EXPORT get_special_halfedges_at_vertex(const TetMesh &mesh, const EP<int> &feature_edge, const EP<int> &valence, const VH vh,
                                std::vector<HEH> &ft_hehs, std::vector<HEH> &bdy_hehs, std::vector<HEH> &itr_hehs);

std::vector<HEH> ALGOHEX_EXPORT get_singular_halfedges(const TetMesh &mesh, const EP<int> &valence, const VH _vh);

std::vector<HEH>
ALGOHEX_EXPORT get_non_feature_singular_halfedges(const TetMesh &mesh, const EP<bool> &feature_feace_edge, const EP<int> &valence,
                                   const VH _vh);

void
ALGOHEX_EXPORT count_special_edges_at_vertex(const TetMesh &mesh, const EP<bool> &feature_face_edge, const EP<int> &feature_edge,
                              const EP<int> &valence, const VH vh,
                              int &n_nff_val_1, int &n_nff_val1, int &n_nff_val_u, int &n_nff, int &n_bdy_val_1,
                              int &n_bdy_val1, int &n_bdy_val_u, int &n_bdy, int &n_cmpl,
                              int &n_ff_val_1, int &n_ff_val1, int &n_ff_val_u, int &n_ff, int &n_fe);

std::set<EH>
ALGOHEX_EXPORT incident_non_ffe_singular_edges_in_same_region(const TetMesh &mesh, const EP<bool> &feature_face_edge,
                                               const EP<int> &valence,
                                               const VH vh, const std::set<CH> &por_chs);

int
ALGOHEX_EXPORT n_incident_non_ffe_singular_edges_in_same_region(const TetMesh &mesh, const EP<bool> &feature_face_edge,
                                                 const EP<int> &valence,
                                                 const VH vh, const std::set<CH> &por_chs);

std::vector<HEH>
ALGOHEX_EXPORT get_halfedges_in_feature_sector(const TetMesh &mesh, const FP<int> &feature_face, const EP<int> &feature_edge,
                                const EP<int> &valence, const HEH _heh_s,
                                std::set<FH> &_visited_fhs, std::set<FH> &_sec_fhs);

//sectors of feature edges or singular edge on feature surface
void
ALGOHEX_EXPORT get_feature_sectors_at_feature_edge_vertex(const TetMesh &mesh, const FP<int> &feature_face, const EP<int> &feature_edge,
                                           const VP<bool> &feature_edge_vertex,
                                           const EP<int> &valence, const VH vh,
                                           std::vector<std::vector<HEH>> &v_sec_hehs,
                                           std::vector<std::set<FH>> &v_sec_fhs);

std::vector<HEH>
ALGOHEX_EXPORT get_boundary_halfedges_in_feature_sector_ccw(const TetMesh &mesh, const EP<int> &feature_edge, const EP<int> &valence,
                                             const HEH _heh_s);

std::vector<HEH>
ALGOHEX_EXPORT get_boundary_halfedges_in_feature_sector_cw(const TetMesh &mesh, const EP<int> &feature_edge, const HEH _heh_s);

void ALGOHEX_EXPORT get_k_ring_vertices(TetMesh &mesh, const int k, std::set<VH> &_initial_vertices);

std::set<VH> ALGOHEX_EXPORT get_k_ring_vertices_of_singular_graph(TetMesh &mesh, const int k);

std::set<CH> ALGOHEX_EXPORT get_k_ring_cells(TetMesh &mesh, const std::set<VH> &_kring_vhs);

std::set<CH> ALGOHEX_EXPORT get_k_ring_cells_of_singular_graph(TetMesh &mesh, const int k);

template<typename U>
std::set<std::pair<double, EH>, U> get_k_ring_edges_of_singular_graph(TetMesh &mesh, const int k);
extern template ALGOHEX_EXPORT std::set<std::pair<double, EH>, std::greater<std::pair<double, EH>>> get_k_ring_edges_of_singular_graph<std::greater<std::pair<double, EH>>>(TetMesh &, const int);
extern template ALGOHEX_EXPORT std::set<std::pair<double, EH>, std::less<std::pair<double, EH>>> get_k_ring_edges_of_singular_graph<std::less<std::pair<double, EH>>>(TetMesh &, const int);

std::set<EH> ALGOHEX_EXPORT get_k_ring_edges_of_singular_graph(TetMesh &mesh, const int k);

template<typename U>
std::set<std::pair<double, EH>, U> get_k_ring_edges(TetMesh &mesh, const std::vector<EH> &_sg_ehs, const int k);
extern template ALGOHEX_EXPORT std::set<std::pair<double, EH>, std::greater<std::pair<double, EH>>> get_k_ring_edges<std::greater<std::pair<double, EH>>>(TetMesh &, const std::vector<EH> &, const int);
extern template ALGOHEX_EXPORT std::set<std::pair<double, EH>, std::less<std::pair<double, EH>>> get_k_ring_edges<std::less<std::pair<double, EH>>>(TetMesh &, const std::vector<EH> &, const int);

//same valence
std::vector<HEH>
ALGOHEX_EXPORT get_halfedges_on_singular_arc_without_zipper_node(const TetMesh &mesh, const EP<int> &valence, const HEH _heh_s);

//ignore interior singular
std::vector<HEH>
ALGOHEX_EXPORT get_halfedges_of_singular_arc_on_feature_face(const TetMesh &mesh, const EP<int> &valence,
                                              const EP<bool> &feature_face_edge, const HEH _heh_s,
                                              const bool _same_valence = false);

std::vector<HEH>
ALGOHEX_EXPORT get_halfedges_of_singular_arc_kept_on_feature_face(const TetMesh &mesh, const EP<int> &valence,
                                                   const EP<bool> &keep_on_face, const EP<int> &feature_edge,
                                                   const HEH _heh_s);

//singular arc that ignore valence
std::vector<VH> ALGOHEX_EXPORT get_vertices_on_singular_arc(const TetMesh &mesh, const EP<int> &valence, const HEH heh);

std::vector<VH>
ALGOHEX_EXPORT get_vertices_on_singular_arc_with_intertior_surface(const TetMesh &mesh, const EP<int> &valence,
                                                    const VP<bool> &feature_face_vertex, const HEH heh);

//singular vertices that are on the outgoing singular arc
std::vector<VH> ALGOHEX_EXPORT get_vertices_on_singular_arc_directional(const TetMesh &mesh, const EP<int> &valence, const HEH heh);

//singular arc that ignore valence
std::vector<HEH> ALGOHEX_EXPORT get_halfedges_on_singular_arc(const TetMesh &mesh, const EP<int> &valence, const HEH heh);

//singular arc that ignore valence and ends at feature surface
std::vector<HEH> ALGOHEX_EXPORT get_halfedges_on_singular_arc_with_intertior_surface(const TetMesh &mesh, const EP<int> &valence,
                                                                      const VP<bool> &feature_face_vertex,
                                                                      const HEH heh);

//singular arc that ignore valence
std::vector<EH> ALGOHEX_EXPORT get_edges_on_singular_arc(const TetMesh &mesh, const EP<int> &valence, const HEH heh);

//halfedges that are on the outgoing singular arc
std::vector<HEH> ALGOHEX_EXPORT get_halfedges_on_singular_arc_directional(const TetMesh &mesh, const EP<int> &valence, const HEH heh);

//halfedges that are on the outgoing feature arc
std::vector<HEH>
ALGOHEX_EXPORT get_halfedges_on_feature_arc_directional(const TetMesh &mesh, const EP<int> &valence, const EP<int> &feature_edge,
                                         const HEH heh);

//no distinguish of valence
std::vector<HEH>
ALGOHEX_EXPORT get_halfedges_on_feature_arc_directional(const TetMesh &mesh, const EP<int> &feature_edge, const HEH heh);

std::set<CH> ALGOHEX_EXPORT get_onering_cells(const TetMesh &mesh, const VH _vh);

std::set<CH> ALGOHEX_EXPORT get_part_of_onering_cells(const TetMesh &mesh, const FP<int> &feature_fprop, const VH vh,
                                       const std::set<CH> &onering_chs, const CH ch_s, std::set<HFH> &ft_hfhs);

std::map<CH, int>
ALGOHEX_EXPORT get_region_id_of_onering_cells(const TetMesh &mesh, const FP<int> &feature_fprop, const VH vh);
}
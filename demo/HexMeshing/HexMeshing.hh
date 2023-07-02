/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

#pragma once

//== INCLUDES =================================================================

#include <string>
#include <AlgoHex/TypeDef.hh>
#include <AlgoHex/Args.hh>
#include <AlgoHex/SmoothOctahedralFieldGeneratorT.hh>
#include <AlgoHex/SmoothOctahedralFieldGeneratorCellBasedT.hh>
#include <AlgoHex/FieldConstraintGenerator.hh>
#include <AlgoHex/SingularGraphExtractionT.hh>
#include <AlgoHex/Parametrization3DT.hh>
#include <AlgoHex/FrameFieldOptimizer3DT.hh>
#include <AlgoHex/Config/Export.hh>

#include <HexEx/HexExtractor.hh>

#include <AlgoHex/LocallyMeshableField/LocallyMeshableFieldGenerationT.hh>
#include <AlgoHex/LocallyMeshableField/OneringParameterizationChecker.hh>
#include <AlgoHex/LocallyMeshableField/LocalMeshabilityCheckerWithFrames.hh>
#include <AlgoHex/LocallyMeshableField/LocalMeshabilityChecker.hh>

#include <AlgoHex/LocallyMeshableField/Comparisons/SplitPaperT.hh>
#include <AlgoHex/LocallyMeshableField/Comparisons/CollapsePaperT.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{
ALGOHEX_EXPORT
void hexMeshing(Args &args);

void hexMeshing_from_seamless_map(const Args &args);

template<class MeshT>
void load_tetmesh(const Args &args, MeshT &tetmesh);

template<class MeshT>
void get_initial_frame_field(const Args &args, MeshT &tetmesh);

template<class MeshT>
void modify_frame_field(Args &args, MeshT &tetmesh);

template<class MeshT>
void parameterization(const Args &args, MeshT &tetmesh, nlohmann::json &json_data);

template<class MeshT>
void extract_hexmesh(const Args &args, MeshT &tetmesh);

template<class MeshT>
void check_local_meshability(const Args &args, MeshT &tetmesh, nlohmann::json &json_data);

template<typename MeshT>
void write_ovmb_file(const std::string &filename, MeshT &mesh);

template<typename MeshT>
auto read_ovmb_file(const std::string &filename, MeshT &mesh);

template<class MeshT>
void initialize_feature_properties(const Args &args, MeshT &tetmesh);

template<class MeshT>
std::vector<EH> get_feature_edges(MeshT &tetmesh);

template<class MeshT>
void load_cell_quaternions(const Args &args, MeshT &tetmesh);

template<class MeshT>
void save_cell_quaternions(const std::string &_filename, MeshT &_mesh);

//    template <class MeshT>
//    void get_algohex_feature_face_property(MeshT& tetmesh);

template<class MeshT>
void set_locally_hexmeshable_field_coeffs(Args &args, LocallyMeshableFieldGenerationT <MeshT> &lmfg);

double hexmesh_volume(const HexEx::HexahedralMesh &_hexmesh);

void initialize_json_file(const std::string &json_path, nlohmann::json &json_data);

void create_sub_hexmesh(const HexEx::HexahedralMesh &_hexmesh, HexEx::HexahedralMesh &_sub_mesh);

void store_json_file(const Args &_args, const nlohmann::json &_json_data);

template<class MeshT>
void attach_feature_information(MeshT &tetmesh, const std::string &_filename);

//=============================================================================
} // namespace AlgoHex
//=============================================================================

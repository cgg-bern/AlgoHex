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
#include <nlohmann/json.hpp>

#include <Typedefs.hh>

#include <AlgoHex/LocallyMeshableField/LocallyMeshableFieldGeneration.hh>

//== FORWARDDECLARATIONS ======================================================

//== NAMESPACES ===============================================================

namespace AlgoHex
{
ALGOHEX_EXPORT
void hexMeshing(Args &args);

void hexMeshing_from_seamless_map(const Args &args);

void load_tetmesh(const Args &args, TetMesh &tetmesh);

void get_initial_frame_field(const Args &args, TetMesh &tetmesh);

void modify_frame_field(Args &args, TetMesh &tetmesh);

void parameterization(const Args &args, TetMesh &tetmesh, nlohmann::json &json_data);

void extract_hexmesh(const Args &args, TetMesh &tetmesh);

void check_local_meshability(const Args &args, TetMesh &tetmesh, nlohmann::json &json_data);

void write_ovmb_file(const std::string &filename, TetMesh &mesh);

void read_ovmb_file(const std::string &filename, TetMesh &mesh);

void initialize_feature_properties(const Args &args, TetMesh &tetmesh);

std::vector<EH> get_feature_edges(TetMesh &tetmesh);

void load_cell_quaternions(const Args &args, TetMesh &tetmesh);

void save_cell_quaternions(const std::string &_filename, TetMesh &_mesh);

void set_locally_hexmeshable_field_coeffs(Args &args, LocallyMeshableFieldGeneration &lmfg);

double hexmesh_volume(const HexEx::HexahedralMesh &_hexmesh);

void initialize_json_file(const std::string &json_path, nlohmann::json &json_data);

void create_sub_hexmesh(const HexEx::HexahedralMesh &_hexmesh, HexEx::HexahedralMesh &_sub_mesh);

void store_json_file(const Args &_args, const nlohmann::json &_json_data);

void attach_feature_information(TetMesh &tetmesh, const std::string &_filename);

//=============================================================================
} // namespace AlgoHex
//=============================================================================

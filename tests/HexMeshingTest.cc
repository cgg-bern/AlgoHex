/*===========================================================================*\
 *        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
 *         Copyright (C) 2019-2023 David Bommes, University of Bern          *
 *---------------------------------------------------------------------------*
 *  For license information see LICENSE.txt in the AlgoHex root directory.   *
 *  All contributors are stated in CREDITS.txt.                              *
\*===========================================================================*/

// Replicates the HexMeshing demo pipeline (demo/HexMeshing/HexMeshing.cc) on
// demo/HexMeshing/cylinder.ovm using the default command-line arguments of the
// demo (i.e. only "-i <input>" and "-o <output>").
//
// Instead of writing the final hex mesh to disk (the demo's "-o" output), this
// test keeps the extracted mesh in memory and compares the number of
// vertices, edges, faces and polyhedra (cells) against the values the demo
// currently produces:
//
//   Vertices:  11248
//   Edges:     32367
//   Faces:     31056
//   Polyhedra: 9936

#include <gtest/gtest.h>

#include "../demo/HexMeshing/HexMeshing.hh" // pipeline functions from the demo

#include <AlgoHex/TypeDef.hh>
#include <AlgoHex/Geometry.hh>

#include <HexEx/HexExtractor.hh>
#include <HexHex/HexHex.hh>

#include <nlohmann/json.hpp>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace AlgoHex
{
namespace
{

// Element counts of the hex mesh that the demo currently produces for
// demo/HexMeshing/cylinder.ovm with default arguments.
constexpr std::size_t kExpectedVertices = 11248;
constexpr std::size_t kExpectedEdges = 32367;
constexpr std::size_t kExpectedFaces = 31056;
constexpr std::size_t kExpectedCells = 9936;
constexpr double tol_rel = 0.15;

// Same extraction as demo/HexMeshing/HexMeshing.cc::extract_hexmesh, but
// keeps the resulting hex mesh in memory instead of writing it to disk.
std::unique_ptr<HexHex::HexahedralMesh>
extract_hexmesh_in_memory(TetMesh &tetmesh)
{
  // Get the per cell parametrization
  auto c_map_igm =
      tetmesh.request_cell_property<std::map<OpenVolumeMesh::VertexHandle, TetMesh::PointT>>(
          "Parametrization");

  // Convert per cell map to halfface prop used by HexHex
  auto hf_igm = tetmesh.request_halfface_property<OVM::Vec3d>("HexHex::Parametrization");
  for (auto hf_it = tetmesh.hf_iter(); hf_it.is_valid(); ++hf_it)
  {
    OVM::CH ch = tetmesh.incident_cell(*hf_it);
    if (!ch.is_valid())
      continue;
    OVM::VH vh = tetmesh.halfface_opposite_vertex(*hf_it);
    hf_igm[*hf_it] = anyVec2ovm(c_map_igm[ch].at(vh));
  }

  auto config = HexHex::Config();
#ifdef _OPENMP
  config.num_threads = std::min(omp_get_max_threads(), 8);
#else
  config.num_threads = 1;
#endif
  config.assert_igm_validity = true;

  auto hexhex_res = HexHex::extractHexMesh(tetmesh, hf_igm, config);
  if (hexhex_res.success)
  {
    return std::move(hexhex_res.hex_mesh);
  }

  return {};
}

} // namespace

TEST(HexMeshingTest, CylinderHexMeshElementCounts)
{
#if defined(HEXMESHING_TEST_CYLINDER_OVM_PATH)
  const std::string in_file = HEXMESHING_TEST_CYLINDER_OVM_PATH;
#else
  const std::string in_file;
#endif

  if (in_file.empty())
  {
    GTEST_SKIP() << "Path to demo/HexMeshing/cylinder.ovm not configured "
                    "(HEXMESHING_TEST_CYLINDER_OVM_PATH undefined).";
  }

  // Default demo arguments: only the input file is set.
  // outFileName is intentionally left empty so that nothing is written to disk.
  AlgoHex::Args args;
  args.inFileName = in_file;

  TetMesh tetmesh;

  // Pipeline, mirroring demo/HexMeshing/HexMeshing.cc::hexMeshing
  load_tetmesh(args, tetmesh);
  ASSERT_GT(tetmesh.n_cells(), 0u) << "Failed to load input tetmesh: " << in_file;

  get_initial_frame_field(args, tetmesh);
  modify_frame_field(args, tetmesh);

  nlohmann::json json_data;
  json_data["InputTetmeshFile"] = args.inFileName;
  parameterization(args, tetmesh, json_data);

  auto hexmesh = extract_hexmesh_in_memory(tetmesh);
  ASSERT_NE(hexmesh, nullptr);

  EXPECT_NEAR(hexmesh->n_vertices(), kExpectedVertices, tol_rel * kExpectedVertices) << "vertex count changed";
  EXPECT_NEAR(hexmesh->n_edges(), kExpectedEdges, tol_rel * kExpectedEdges) << "edge count changed";
  EXPECT_NEAR(hexmesh->n_faces(), kExpectedFaces, tol_rel * kExpectedFaces) << "face count changed";
  EXPECT_NEAR(hexmesh->n_cells(), kExpectedCells, tol_rel * kExpectedCells) << "hex count changed";
}

} // namespace AlgoHex

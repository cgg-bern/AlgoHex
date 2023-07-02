#include <CLI/CLI.hpp>
#include <AlgoHex/TypeDef.hh>
#include <AlgoHex/MeshGeometry.hh>
#include <AlgoHex/FrameFieldOptimizer3DT.hh>

#include <AlgoHex/LocallyMeshableField/LocalMeshabilityChecker.hh>

#ifdef _WIN32
#  include <windows.h>
#  include <stdlib.h>
#  include <errhandlingapi.h>
#endif


int main(int argc, const char *argv[])
{

#ifdef _WIN32
  SetErrorMode(0);
#endif

  std::string in_filename;

  CLI::App app{"LocalMeshabilityCheck"};
  app.add_option("-i", in_filename, "Input tetrahedral mesh in .ovm/.vtk format (several properties are required!!!).");

  try
  {
    app.parse(argc, argv);
  }
  catch (const CLI::ParseError &e)
  {
    return app.exit(e);
  }

  using Vec3d = OpenVolumeMesh::Geometry::Vec3d;
  using TetMesh = OVM::TetrahedralGeometryKernel<Vec3d, OVM::TetrahedralMeshTopologyKernel>;

  OpenVolumeMesh::IO::FileManager fm;
  TetMesh tetmesh;
//  OpenVolumeMesh::StatusAttrib status(tetmesh);

  if (!in_filename.empty())
  {
    fm.readFile(in_filename, tetmesh);

    AlgoHex::FrameFieldOptimizer3DT <TetMesh> ffopt(tetmesh);
//    ffopt.enable_save_locally_non_meshable(args.locallyNonMeshableFileName);

    ffopt.import_frames_from_vec3d_properties();
    ffopt.import_quaternions_from_frames();
    ffopt.check_valence_consistency();
    ffopt.check_frame_rotation_angles();
//    ffopt.check_local_meshability();

    std::cerr << "------------------------ TEST ALTERNATIVE CHECKER --------------------" << std::endl;
    AlgoHex::LocalMeshabilityChecker lmc(tetmesh);
    lmc.verbose() = true;
    lmc.check_local_meshability(false);
  }

  return 0;
}

#/*===========================================================================*\
#*        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
#*         Copyright (C) 2019-2023 David Bommes, University of Bern          *
#*---------------------------------------------------------------------------*
#*  For license information see LICENSE.txt in the AlgoHex root directory.   *
#*  All contributors are stated in CREDITS.txt.                              *
#\*===========================================================================*/

cmake_minimum_required(VERSION 3.14)

# AutoMOC/AutoUIC
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17")
    cmake_policy(SET CMP0100 NEW)
endif ()

# honor $*_ROOT environment variables:
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()


list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-library/finders)

project(AlgoHex
        VERSION 0.0.1
        LANGUAGES C CXX
        )
set(CMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY ON)
option(BUILD_SHARED_LIBS "Build shared library (*.dll, *.so, *.dylib) instead of static library (*.a, *.lib)" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(ALGOHEX_IS_ROOT TRUE)
    set(ALGOHEX_TARGET_PREFIX "")
else ()
    set(ALGOHEX_IS_ROOT FALSE)
    set(ALGOHEX_TARGET_PREFIX "AlgoHex_")
endif ()


if (ALGOHEX_IS_ROOT)
    if (WIN32)
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build")
    else ()
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
    endif ()

    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    option(STL_DEBUG "Enable STL debug checks. Warning: This must be consistent across the whole build! Not compatible with the libc++ shipped with Apple XCode as of 2019." OFF)
    if (STL_DEBUG)
        # We do not know which STL will be used, defining both is the easist way.
        # https://libcxx.llvm.org/docs/UsingLibcxx.html#libc-configuration-macros
        # https://libcxx.llvm.org/docs/DesignDocs/DebugMode.html#using-debug-mode
        add_definitions(-D_LIBCPP_DEBUG=1)
        # https://gcc.gnu.org/onlinedocs/libstdc%2B%2B/manual/debug_mode_using.html
        add_definitions(-D_GLIBCXX_DEBUG=1)
        add_definitions(-D_GLIBCXX_DEBUG_PEDANTIC=1)
    endif ()
    set(ALGOHEX_CXX_STANDARD 17 CACHE STRING "C++ standard version to use")
    set_property(CACHE ALGOHEX_CXX_STANDARD PROPERTY STRINGS 17 20)
    set(CMAKE_CXX_STANDARD ${ALGOHEX_CXX_STANDARD})
    set(CMAKE_CXX_STANDARD_REQUIRED YES)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif ()

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()


# try to find Mosek before getting the other dependencies, so we don't load
# an incompatible FindMOSEK.cmake from any of them
find_package(Mosek QUIET) # optional dependency


find_package(CHOLMOD CONFIG)
if (NOT CHOLMOD_FOUND)
    if(APPLE)
        message(FATAL_ERROR "Dependency CHOLMOD (from SuiteSparse) not found. You will also require OpenMP. You can install both using Homebrew with 'brew install libomp suitesparse', and setting the environment variable OpenMP_ROOT=/opt/homebrew/opt/libomp/ .")
    else()
        message(FATAL_ERROR "Dependency CHOLMOD (from SuiteSparse) not found. This *may* also be caused by missing OpenMP.")
    endif()
endif()

set(ALGOHEX_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")

set(ACG_COMMON_DO_NOT_COPY_POST_BUILD TRUE)
set(WANT_COMISO_QT OFF CACHE BOOL "Enable Qt support in CoMISo")
include(AlgoHexDependencies)

# add library
add_library(AlgoHex
        src/AlgoHex/Args.hh
        src/AlgoHex/FaceNormalCache.hh
        src/AlgoHex/FaceNormalCache.cc
        src/AlgoHex/FieldConstraints.hh
        src/AlgoHex/FieldConstraints.cc
        src/AlgoHex/FieldConstraintGenerator.hh
        src/AlgoHex/FieldConstraintGenerator.cc
        src/AlgoHex/FrameFieldOptimizer3D.hh
        src/AlgoHex/FrameFieldOptimizer3D.cc
        src/AlgoHex/Geometry.hh
        src/AlgoHex/Geometry.cc
        src/AlgoHex/MeshGeometry.hh
        src/AlgoHex/MeshGeometry.cc
        src/AlgoHex/Parametrization3D.hh
        src/AlgoHex/Parametrization3D.cc
        src/AlgoHex/SingularGraphExtraction.hh
        src/AlgoHex/SingularGraphExtraction.cc
        src/AlgoHex/SmoothOctahedralFieldGenerator.hh
        src/AlgoHex/SmoothOctahedralFieldGenerator.cc
        src/AlgoHex/SmoothOctahedralFieldGeneratorCellBased.hh
        src/AlgoHex/SmoothOctahedralFieldGeneratorCellBased.cc
        src/AlgoHex/Stopwatches.hh
        src/AlgoHex/Stopwatches.cc
        src/AlgoHex/TransitionQuaternionEigen.cc
        src/AlgoHex/TransitionQuaternionEigen.hh
        src/AlgoHex/TypeDef.hh
        src/AlgoHex/AxisAlignment.hh
        src/AlgoHex/AxisAlignment.cc
        src/AlgoHex/AMIPSFrameElement3D.hh
        src/AlgoHex/AMIPSFrameElement3D.cc
        src/AlgoHex/AMIPSFrameElement3DTinyAD.hh
        src/AlgoHex/AMIPSFrameElement3DTinyAD.cc
        src/AlgoHex/BarrierElements.hh
        src/AlgoHex/BarrierElements.cc
        src/AlgoHex/DeformationElements3D.hh
        src/AlgoHex/DeformationElements3D.cc
        src/AlgoHex/DihedralAngleElementTinyAD.hh
        src/AlgoHex/FrameFittingElement3D.hh
        src/AlgoHex/FrameFittingElement3D.cc
        src/AlgoHex/IntegrabilityElements3D.hh
        src/AlgoHex/LinearLeastSquaresElements.hh
        src/AlgoHex/LogDetElement3D.hh
        src/AlgoHex/LogDetElement3D.cc

        src/AlgoHex/Util/StopWatch.hh
        src/AlgoHex/Util/ScopedStopWatch.hh
        src/AlgoHex/Util/ScopedStopWatch.cc

        src/AlgoHex/OpenNL/OpenNL_psm.c
        src/AlgoHex/OpenNL/OpenNL_psm.h

        src/AlgoHex/CommonOVMBDecoders/PropertyCodecsSHCoeffs.hh

        src/AlgoHex/SphericalHarmonics/ConstrainedSHProjector.hh
        src/AlgoHex/SphericalHarmonics/ConstrainedSHProjector.cc
        src/AlgoHex/SphericalHarmonics/OctaToQuat.hh
        src/AlgoHex/SphericalHarmonics/OctaToQuat.cc
        src/AlgoHex/SphericalHarmonics/SDPMatrices.cc
        src/AlgoHex/SphericalHarmonics/SDPMatrices.hh
        src/AlgoHex/SphericalHarmonics/SHCoeffs.hh
        src/AlgoHex/SphericalHarmonics/SHCoeffs.cc
        src/AlgoHex/SphericalHarmonics/SHCoeffRotation.cc
        src/AlgoHex/SphericalHarmonics/SHCoeffRotation.hh
        src/AlgoHex/SphericalHarmonics/SHProjectorInterface.hh
        src/AlgoHex/SphericalHarmonics/SHProjectorRay.hh
        src/AlgoHex/SphericalHarmonics/SHProjectorRay.cc
        src/AlgoHex/SphericalHarmonics/SHProjectorSDP.hh
        src/AlgoHex/SphericalHarmonics/SHProjectorSDP.cc

        src/AlgoHex/LocallyMeshableField/ArcZipping.hh
        src/AlgoHex/LocallyMeshableField/ArcZipping.cc
        src/AlgoHex/LocallyMeshableField/CellInfo.hh
        src/AlgoHex/LocallyMeshableField/Comparisons/CollapsePaper.hh
        src/AlgoHex/LocallyMeshableField/Comparisons/CollapsePaper.cc
        src/AlgoHex/LocallyMeshableField/Comparisons/SplitPaper.hh
        src/AlgoHex/LocallyMeshableField/Comparisons/SplitPaper.cc
        src/AlgoHex/LocallyMeshableField/CommonFuncs.hh
        src/AlgoHex/LocallyMeshableField/CommonFuncs.cc
        src/AlgoHex/LocallyMeshableField/DuplicateOneRingMesh.hh
        src/AlgoHex/LocallyMeshableField/DuplicateOneRingMesh.cc
        src/AlgoHex/LocallyMeshableField/DetachAtInvalidSingularNode.hh
        src/AlgoHex/LocallyMeshableField/DetachAtInvalidSingularNode.cc
        src/AlgoHex/LocallyMeshableField/EdgeMonodromyHelper.hh
        src/AlgoHex/LocallyMeshableField/EdgeMonodromyHelper.cc
        src/AlgoHex/LocallyMeshableField/EnergyFunctions.hh
        src/AlgoHex/LocallyMeshableField/EnergyFunctions.cc
        src/AlgoHex/LocallyMeshableField/EnsureEdgeMeshability.hh
        src/AlgoHex/LocallyMeshableField/EnsureEdgeMeshability.cc
        src/AlgoHex/LocallyMeshableField/ExpDeformationElements3D.hh
        src/AlgoHex/LocallyMeshableField/ExpDeformationElements3D.cc
        src/AlgoHex/LocallyMeshableField/FieldAngleCalculator.hh
        src/AlgoHex/LocallyMeshableField/FieldAngleCalculator.cc
        src/AlgoHex/LocallyMeshableField/FixZipperNode.hh
        src/AlgoHex/LocallyMeshableField/FixZipperNode.cc
        src/AlgoHex/LocallyMeshableField/FixConstrainedZipperNode.hh
        src/AlgoHex/LocallyMeshableField/FixConstrainedZipperNode.cc
        src/AlgoHex/LocallyMeshableField/FrameFieldOptimization.hh
        src/AlgoHex/LocallyMeshableField/FrameFieldOptimization.cc
        src/AlgoHex/LocallyMeshableField/LocallyMeshableFieldGeneration.hh
        src/AlgoHex/LocallyMeshableField/LocallyMeshableFieldGeneration.cc
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityChecker.hh
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityChecker.cc
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityCheckerWithFrames.hh
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityCheckerWithFrames.cc
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityRepair.hh
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityRepair.cc
        src/AlgoHex/LocallyMeshableField/MeshOptimization.hh
        src/AlgoHex/LocallyMeshableField/MeshOptimization.cc
        src/AlgoHex/LocallyMeshableField/MeshProperties.hh
        src/AlgoHex/LocallyMeshableField/MeshProperties.cc
        src/AlgoHex/LocallyMeshableField/OneringParameterizationChecker.hh
        src/AlgoHex/LocallyMeshableField/OneringParameterizationChecker.cc
        src/AlgoHex/LocallyMeshableField/PushSingularVertex.hh
        src/AlgoHex/LocallyMeshableField/PushSingularVertex.cc
        src/AlgoHex/LocallyMeshableField/QuaternionsSmoothing.hh
        src/AlgoHex/LocallyMeshableField/QuaternionsSmoothing.cc
        src/AlgoHex/LocallyMeshableField/RemeshingAssist.hh
        src/AlgoHex/LocallyMeshableField/RemeshingAssist.cc
        src/AlgoHex/LocallyMeshableField/SeparableSingularArcFinder.hh
        src/AlgoHex/LocallyMeshableField/SeparableSingularArcFinder.cc
        src/AlgoHex/LocallyMeshableField/SingularVertexOptProblem.hh
        src/AlgoHex/LocallyMeshableField/SingularVertexOptProblem.cc
        src/AlgoHex/LocallyMeshableField/SplitHelper.hh
        src/AlgoHex/LocallyMeshableField/SplitHelper.cc
        src/AlgoHex/LocallyMeshableField/TetRemeshing.hh
        src/AlgoHex/LocallyMeshableField/TetRemeshing.cc
        src/AlgoHex/LocallyMeshableField/TransversalUnzipping.hh
        src/AlgoHex/LocallyMeshableField/TransversalUnzipping.cc
        src/AlgoHex/LocallyMeshableField/VertexOptimize.hh
        src/AlgoHex/LocallyMeshableField/VertexOptimize.cc
        )

set_source_files_properties(
        src/AlgoHex/SphericalHarmonics/SHCoeffRotation.cc
        src/AlgoHex/SphericalHarmonics/SHCoeffs.cc
        src/AlgoHex/SphericalHarmonics/SHProjectorRay.cc
        PROPERTIES COMPILE_FLAGS
        "$<IF:$<CXX_COMPILER_ID:MSVC>,/O2,-O3>"
)

target_compile_features(AlgoHex PUBLIC cxx_std_17)
add_library(AlgoHex::AlgoHex ALIAS AlgoHex)

target_compile_options(AlgoHex PUBLIC
        "$<$<CXX_COMPILER_ID:MSVC>:/bigobj>"
        # cmake by default adds /W3, that's too much for us:
        "$<$<CXX_COMPILER_ID:MSVC>:/W1>"
        # disable MSVC inline warnings that are enabled with /Ob2, which is in the default cmake release cxxflags:
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4710>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4711>"

        "$<$<CXX_COMPILER_ID:MSVC>:/wd4514>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4571>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5045>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5024>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5025>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5026>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd5027>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4625>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4626>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4820>"
        "$<$<CXX_COMPILER_ID:MSVC>:/wd4464>"
        )
target_compile_definitions(AlgoHex
        PUBLIC
        "INCLUDE_TEMPLATES"
        "$<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_DEPRECATE>" # for GMM
        # SHCoeffs derives from an Eigen Matrix, on MSVC
        # we currently cannot instantiate the template for dllexport.
        # We'll catch the compile time asserts on the other compilers.
        "$<$<CXX_COMPILER_ID:MSVC>:EIGEN_NO_STATIC_ASSERT>"
        PRIVATE
        "GEO_STATIC_LIBS" # for OpenNL
        )

set_target_properties(AlgoHex PROPERTIES
        VERSION ${AlgoHex_VERSION_MAJOR}.${AlgoHex_VERSION_MINOR}
        SOVERSION ${AlgoHex_VERSION_MAJOR}.${AlgoHex_VERSION_MINOR}
        CXX_VISIBILITY_PRESET hidden
        )

if (ALGOHEX_IS_ROOT)
    set_target_properties(AlgoHex PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
            ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
            )
endif ()

if (TARGET MOSEK::FusionCXX)
    set(ALGOHEX_WITH_MOSEK "1") # used by configure_file on Config.hh
    target_link_libraries(AlgoHex PUBLIC MOSEK::FusionCXX)
    message("Compiling AlgoHex with MOSEK support")
else ()
    set(ALGOHEX_WITH_MOSEK "0") # used by configure_file on Config.hh
    message("Compiling AlgoHex without MOSEK support")
endif ()

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(AlgoHex PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_C)
endif ()

target_link_libraries(AlgoHex PUBLIC
        OpenVolumeMesh::OpenVolumeMesh
        CoMISo::CoMISo
        HexEx::HexEx
        HexHex::HexHex
        Eigen3::Eigen
        GMM::GMM
        QGP3D::QGP3D
        tinyad::tinyad
        )


configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/src/AlgoHex/Config/Version.hh.in"
        "${CMAKE_CURRENT_BINARY_DIR}/src/AlgoHex/Config/Version.hh"
)

include(GenerateExportHeader)
generate_export_header(AlgoHex
        BASE_NAME ALGOHEX
        EXPORT_FILE_NAME "src/AlgoHex/Config/Export.hh"
        )


#Set target properties
target_include_directories(AlgoHex
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        #    PRIVATE
        #"${CMAKE_CURRENT_BINARY_DIR}/src"
        )

option(ALGOHEX_BUILD_DEMO "Build AlgoHex demo commandline tool" ${ALGOHEX_IS_ROOT})
if (ALGOHEX_BUILD_DEMO)
    add_subdirectory(demo/HexMeshing)
    add_subdirectory(demo/LocalMeshabilityCheck)
endif ()

option(ALGOHEX_BUILD_TESTS "Build AlgoHex unit tests" ${ALGOHEX_IS_ROOT})
if (ALGOHEX_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif ()


include(CMakePackageConfigHelpers)
include(GNUInstallDirs)


set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/AlgoHex)

configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/AlgoHexConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/AlgoHexConfig.cmake"
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

write_basic_package_version_file(
        AlgoHexConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
)

option(ALGOHEX_NO_INSTALL "Do not try to use cmake install/export functionality" ON)

if(NOT ALGOHEX_NO_INSTALL)

    install(FILES
            "${CMAKE_CURRENT_BINARY_DIR}/AlgoHexConfigVersion.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/AlgoHexConfig.cmake"
            DESTINATION ${INSTALL_CONFIGDIR})

    install(DIRECTORY
            "${CMAKE_CURRENT_BINARY_DIR}/src/AlgoHex/Config"
            DESTINATION include/AlgoHex)

# Install Header Files
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/AlgoHex
            DESTINATION include
            FILES_MATCHING PATTERN "*.hh"
            )

    install(TARGETS AlgoHex
            EXPORT AlgoHexTargets
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            )

    install(EXPORT AlgoHexTargets
            FILE AlgoHexTargets.cmake
            NAMESPACE AlgoHex::
            DESTINATION ${INSTALL_CONFIGDIR}
            )

    export(EXPORT AlgoHexTargets
            NAMESPACE AlgoHex::)
    export(PACKAGE AlgoHex)

endif()

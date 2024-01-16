#/*===========================================================================*\
#*        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
#*         Copyright (C) 2019-2023 David Bommes, University of Bern          *
#*---------------------------------------------------------------------------*
#*  For license information see LICENSE.txt in the AlgoHex root directory.   *
#*  All contributors are stated in CREDITS.txt.                              *
#\*===========================================================================*/

cmake_minimum_required(VERSION 3.9)

# AutoMOC/AutoUIC
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17")
    cmake_policy(SET CMP0100 NEW)
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()


# try to find Mosek before getting the other dependencies, so we don't load
# an incompatible FindMOSEK.cmake from any of them
find_package("MOSEK")

set(ALGOHEX_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")

set(ACG_COMMON_DO_NOT_COPY_POST_BUILD TRUE)
set(WANT_COMISO_QT OFF CACHE BOOL "Enable Qt support in CoMISo")
include(AlgoHexDependencies)

# add library
add_library(AlgoHex
        src/AlgoHex/Args.hh
        src/AlgoHex/FaceNormalCache.hh
        src/AlgoHex/FieldConstraints.hh
        src/AlgoHex/FieldConstraintGenerator.hh
        src/AlgoHex/FrameFieldOptimizer3DT.hh
        src/AlgoHex/FrameFieldOptimizer3DT_impl.hh
        src/AlgoHex/Geometry.hh
        src/AlgoHex/MeshGeometry.hh
        src/AlgoHex/Parametrization3DT.hh
        src/AlgoHex/Parametrization3DT_impl.hh
        src/AlgoHex/SingularGraphExtractionT.hh
        src/AlgoHex/SingularGraphExtractionT_impl.hh
        src/AlgoHex/SmoothOctahedralFieldGeneratorT.hh
        src/AlgoHex/SmoothOctahedralFieldGeneratorT_impl.hh
        src/AlgoHex/SmoothOctahedralFieldGeneratorCellBasedT.hh
        src/AlgoHex/SmoothOctahedralFieldGeneratorCellBasedT_impl.hh
        src/AlgoHex/Stopwatches.hh
        src/AlgoHex/Stopwatches.cc
        src/AlgoHex/TransitionQuaternionEigen.cc
        src/AlgoHex/TransitionQuaternionEigen.hh
        src/AlgoHex/TypeDef.hh
        src/AlgoHex/AxisAlignment.hh
        src/AlgoHex/AMIPSFrameElement3D.hh
        src/AlgoHex/AMIPSFrameElement3DTinyAD.hh
        src/AlgoHex/BarrierElements.hh
        src/AlgoHex/DeformationElements3D.hh
        src/AlgoHex/DihedralAngleElementTinyAD.hh
        src/AlgoHex/FrameFittingElement3D.hh
        src/AlgoHex/IntegrabilityElements3D.hh
        src/AlgoHex/LinearLeastSquaresElements.hh
        src/AlgoHex/LogDetElement3D.hh

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

        src/AlgoHex/LocallyMeshableField/ArcZippingT.hh
        src/AlgoHex/LocallyMeshableField/ArcZippingT_impl.hh
        src/AlgoHex/LocallyMeshableField/CellInfo.hh
        src/AlgoHex/LocallyMeshableField/CommonFuncs.hh
        src/AlgoHex/LocallyMeshableField/CommonFuncs_impl.hh
        src/AlgoHex/LocallyMeshableField/DuplicateOneRingMeshT.hh
        src/AlgoHex/LocallyMeshableField/DuplicateOneRingMeshT_impl.hh
        src/AlgoHex/LocallyMeshableField/DetachAtInvalidSingularNodeT.hh
        src/AlgoHex/LocallyMeshableField/DetachAtInvalidSingularNodeT_impl.hh
        src/AlgoHex/LocallyMeshableField/EdgeMonodromyHelperT.hh
        src/AlgoHex/LocallyMeshableField/EnergyFunctions.hh
        src/AlgoHex/LocallyMeshableField/EnsureEdgeMeshabilityT.hh
        src/AlgoHex/LocallyMeshableField/EnsureEdgeMeshabilityT_impl.hh
        src/AlgoHex/LocallyMeshableField/ExpDeformationElements3D.hh
        src/AlgoHex/LocallyMeshableField/FieldAngleCalculatorT.hh
        src/AlgoHex/LocallyMeshableField/FieldAngleCalculatorT_impl.hh
        src/AlgoHex/LocallyMeshableField/FixZipperNodeT.hh
        src/AlgoHex/LocallyMeshableField/FixZipperNodeT_impl.hh
        src/AlgoHex/LocallyMeshableField/FixConstrainedZipperNodeT.hh
        src/AlgoHex/LocallyMeshableField/FixConstrainedZipperNodeT_impl.hh
        src/AlgoHex/LocallyMeshableField/FrameFieldOptimizationT.hh
        src/AlgoHex/LocallyMeshableField/FrameFieldOptimizationT_impl.hh
        src/AlgoHex/LocallyMeshableField/LocallyMeshableFieldGenerationT.hh
        src/AlgoHex/LocallyMeshableField/LocallyMeshableFieldGenerationT_impl.hh
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityChecker.hh
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityCheckerWithFrames.hh
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityRepairT.hh
        src/AlgoHex/LocallyMeshableField/LocalMeshabilityRepairT_impl.hh
        src/AlgoHex/LocallyMeshableField/MeshOptimizationT.hh
        src/AlgoHex/LocallyMeshableField/MeshOptimizationT_impl.hh
        src/AlgoHex/LocallyMeshableField/MeshProperties.hh
        src/AlgoHex/LocallyMeshableField/OneringParameterizationChecker.hh
        src/AlgoHex/LocallyMeshableField/PushSingularVertexT.hh
        src/AlgoHex/LocallyMeshableField/PushSingularVertexT_impl.hh
        src/AlgoHex/LocallyMeshableField/QuaternionsSmoothing.hh
        src/AlgoHex/LocallyMeshableField/RemeshingAssist.hh
        src/AlgoHex/LocallyMeshableField/RemeshingAssist_impl.hh
        src/AlgoHex/LocallyMeshableField/SeparableSingularArcFinderT.hh
        src/AlgoHex/LocallyMeshableField/SeparableSingularArcFinderT_impl.hh
        src/AlgoHex/LocallyMeshableField/SingularVertexOptProblem.hh
        src/AlgoHex/LocallyMeshableField/SplitHelperT.hh
        src/AlgoHex/LocallyMeshableField/TetRemeshingT.hh
        src/AlgoHex/LocallyMeshableField/TetRemeshingT_impl.hh
        src/AlgoHex/LocallyMeshableField/TransversalUnzippingT.hh
        src/AlgoHex/LocallyMeshableField/TransversalUnzippingT_impl.hh
        src/AlgoHex/LocallyMeshableField/VertexOptimizeT.hh
        src/AlgoHex/LocallyMeshableField/VertexOptimizeT_impl.hh
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

if (TARGET QGP3D::QGP3D)
    target_link_libraries(AlgoHex PUBLIC QGP3D::QGP3D)
endif ()

target_link_libraries(AlgoHex PUBLIC
        OpenVolumeMesh::OpenVolumeMesh
        CoMISo::CoMISo
        HexEx::HexEx
        Eigen3::Eigen
        GMM::GMM
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


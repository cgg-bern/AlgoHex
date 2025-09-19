
include(FetchContent)

set(FETCHCONTENT_QUIET OFF)
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(FETCHCONTENT_UPDATES_DISCONNECTED TRUE)

if(NOT TARGET GMM::GMM)
    FetchContent_Declare(gmm
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        URL http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.4.2.tar.gz
        URL_HASH SHA224=8f4951901a55a1d1987d8199ed0e0b36ff2da26c870a4c3d55694a14
        SOURCE_DIR "${EXTERNAL_DIR}/gmm"
    )
    FetchContent_Populate(gmm)

    message("Downloaded GMM to ${gmm_SOURCE_DIR}")
    set(GMM_INCLUDE_DIR "${gmm_SOURCE_DIR}/include")
    find_package(GMM REQUIRED)
endif()

if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    FetchContent_Declare(openvolumemesh
        GIT_REPOSITORY  https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git
        GIT_TAG v3.3.0
        SOURCE_DIR "${EXTERNAL_DIR}/OpenVolumeMesh"
        )
    FetchContent_MakeAvailable(openvolumemesh)
endif()

if(NOT TARGET HexEx)
    FetchContent_Declare(hexex
        GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/HexEx/libHexEx.git
        GIT_TAG        algohex
        SOURCE_DIR "${EXTERNAL_DIR}/libHexEx"
        )
    FetchContent_MakeAvailable(hexex)
endif()

if(NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(eigen
        #GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        #GIT_TAG 3.4.0
        # temporary fix, 2025-01-16: my (mh) MR is not merged yet, but is required to fix a build error:
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        GIT_TAG 464c1d097891a1462ab28bf8bb763c1683883892 # master 2025-03-10
        SOURCE_DIR "${EXTERNAL_DIR}/eigen"
    )
    FetchContent_Populate(eigen)
    message("Downloaded Eigen3 to ${eigen_SOURCE_DIR}")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE "$<BUILD_INTERFACE:${eigen_SOURCE_DIR}>")
    #target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_HAS_STD_RESULT_OF=0)
endif()


if(NOT TARGET tinyad::tinyad)
    FetchContent_Declare(tinyad
        GIT_REPOSITORY  https://github.com/patr-schm/TinyAD
        GIT_TAG 81fab13c3884b787c0f03bbbbb95b0f794e54434 # main 2023-10-10
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/tinyad"
        )
    #FetchContent_Populate(tinyad)
    #set(TINYAD_DIR "${tinyad_SOURCE_DIR}/include")
    #set(TINYAD_INCLUDE_DIR "${tinyad_SOURCE_DIR}/include")
#[[= # rely on comiso's tinyad finder for now to avoid this issue:
        CMake Error in external/CoMISo/CMakeLists.txt:
        export called with target "CoMISo" which requires target "TinyAD" that is
        not in any export set.
    FetchContent_MakeAvailable(tinyad)

    add_library(tinyad::tinyad ALIAS TinyAD)
#]]
    FetchContent_MakeAvailable(tinyad)
    add_library(tinyad::tinyad ALIAS TinyAD)
endif()

if(NOT TARGET CoMISo)
    FetchContent_Declare(comiso
        GIT_REPOSITORY https://gitlab.vci.rwth-aachen.de:9000/CoMISo/CoMISo.git
        GIT_TAG 95609ed64337d0fd3c528f70c80334a74b81bb51 # cgg2-mh 2024-08-05
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/CoMISo" # case matters
        )
    set(COMISO_NO_INSTALL YES)
    FetchContent_MakeAvailable(comiso)
    if (NOT TARGET CoMISo::CoMISo)
        error(FATAL_ERROR "CoMISo target missing")
    endif()
endif()

if(NOT TARGET CLI11::CLI11)
    FetchContent_Declare(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        v2.5.0
        SOURCE_DIR "${EXTERNAL_DIR}/CLI11"
        )
    FetchContent_MakeAvailable(cli11)
endif()

# important: TS3D before MC3D
if(NOT TARGET TS3D:TS3D)
    FetchContent_Declare(ts3d
        GIT_REPOSITORY https://github.com/cgg-bern/TrulySeamless3D
        GIT_TAG        cgg
        SOURCE_DIR "${EXTERNAL_DIR}/TrulySeamless3D"
        )
    FetchContent_MakeAvailable(ts3d)
endif()

# important: MC3D before QGP3D
if(NOT TARGET MC3D:MC3D)
    FetchContent_Declare(mc3d
        GIT_REPOSITORY https://github.com/cgg-bern/MC3D
        GIT_TAG        gmp-finder
        SOURCE_DIR "${EXTERNAL_DIR}/MC3D"
        )
    FetchContent_MakeAvailable(mc3d)
endif()


if(NOT TARGET QGP3D:QGP3D)
    FetchContent_Declare(qgp3d
        GIT_REPOSITORY https://github.com/HendrikBrueckler/QGP3D
        GIT_TAG        1acfb740547cc2fba38a7466825b1794fb31d7e2 # main on  2025-01-30: "Priority bug fixed"
        SOURCE_DIR "${EXTERNAL_DIR}/QGP3D"
        )
    FetchContent_MakeAvailable(qgp3d)
endif()

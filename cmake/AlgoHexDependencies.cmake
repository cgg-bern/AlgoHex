
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
        EXCLUDE_FROM_ALL
    )
    FetchContent_MakeAvailable(gmm)

    message("Downloaded GMM to ${gmm_SOURCE_DIR}")
    add_library(gmm INTERFACE)
    add_library(GMM::GMM ALIAS gmm)
    target_include_directories(gmm INTERFACE "${gmm_SOURCE_DIR}/include")
endif()

if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    FetchContent_Declare(openvolumemesh
        GIT_REPOSITORY https://github.com/OpenVolumeMesh/OpenVolumeMesh/
        #GIT_REPOSITORY  https://gitlab.vci.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git
        GIT_TAG b70df3adab139d37fbd118d3287577d541b3590a # master on 2026-10-21
        SOURCE_DIR "${EXTERNAL_DIR}/OpenVolumeMesh"
        )
    FetchContent_MakeAvailable(openvolumemesh)
endif()

# before HexHex
if (NOT TARGET Timekeeper::libTimekeeper)
    FetchContent_Declare(timekeeper
        GIT_REPOSITORY  https://github.com/cgg-bern/libTimekeeper
        GIT_TAG f4f486648faac7c740535678b45b387078e26918 # main 2025-12-04 (up-to-date 2026-10-17)
        EXCLUDE_FROM_ALL
        SOURCE_DIR "${EXTERNAL_DIR}/libTimekeeper"
        )
    FetchContent_MakeAvailable(timekeeper)
    message("Downloaded libTimekeeper to ${timekeeper_SOURCE_DIR}")
endif()

if(NOT TARGET HexHex::HexHex)
    FetchContent_Declare(hexhex
        GIT_REPOSITORY https://github.com/cgg-bern/libHexHex.git
        GIT_TAG b1f37fc23938dd69e3c1eeff819902b0a26ca200
        SOURCE_DIR "${EXTERNAL_DIR}/libHexHex"
        )
    FetchContent_MakeAvailable(hexhex)
endif()

if(NOT TARGET HexEx::HexEx)
    FetchContent_Declare(hexex
        GIT_REPOSITORY https://gitlab.vci.rwth-aachen.de:9000/HexEx/libHexEx.git
        GIT_TAG f394a5d638f3bdac92bb74f2502021757d2365d6# algohex-ovm-default-prop branch on 2016-10-16
        SOURCE_DIR "${EXTERNAL_DIR}/libHexEx"
        )
    FetchContent_MakeAvailable(hexex)
endif()


if(NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(eigen
        #URL https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.bz2
        #URL_HASH SHA256=bdca0ec740fb83be21fe038699923f4c589ead9ab904f4058a9c97752e60d50b
        URL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.bz2
        URL_HASH SHA256=e4de6b08f33fd8b8985d2f204381408c660bffa6170ac65b68ae1bd3cd575c0a
        #GIT_REPOSITORY https://gitlab.com/libeigen/eigen
        #GIT_TAG 464c1d097891a1462ab28bf8bb763c1683883892 # master 2025-03-10
        SOURCE_DIR "${EXTERNAL_DIR}/eigen"
        SOURCE_SUBDIR "nonexisting. Do not use Eigen CMake, we use it header-only."
    )
    FetchContent_MakeAvailable(eigen)
    message("Downloaded Eigen3 to ${eigen_SOURCE_DIR}")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories(Eigen3::Eigen INTERFACE "$<BUILD_INTERFACE:${eigen_SOURCE_DIR}>")
    #target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_HAS_STD_RESULT_OF=0)
endif()


if(NOT TARGET tinyad::tinyad)
    FetchContent_Declare(tinyad
        GIT_REPOSITORY https://github.com/cgg-bern/TinyAD
        GIT_TAG dev/cmake-leave-cxx-standard
        #GIT_REPOSITORY  https://github.com/patr-schm/TinyAD
        #GIT_TAG 81fab13c3884b787c0f03bbbbb95b0f794e54434 # main 2023-10-10
        SOURCE_DIR "${EXTERNAL_DIR}/tinyad"
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
endif()

if(NOT TARGET CoMISo::CoMISo)
    FetchContent_Declare(comiso
        GIT_REPOSITORY https://gitlab.vci.rwth-aachen.de:9000/CoMISo/CoMISo.git
        GIT_TAG ccc10fc6dc876d928ac41e56463dc2f0deec5874 # master 2026-10-17
        SOURCE_DIR "${EXTERNAL_DIR}/CoMISo"
        )
    set(COMISO_NO_INSTALL YES)
    set(COMISO_ENABLE_DEFAULT OFF CACHE BOOL "" FORCE)
    set(COMISO_ENABLE_TINYAD ON CACHE BOOL "" FORCE)
    set(COMISO_ENABLE_IPOPT ON CACHE BOOL "" FORCE)
    set(COMISO_ENABLE_GUROBI ON CACHE BOOL "" FORCE)
    set(COMISO_ENABLE_SUITESPARSE_CHOLMOD ON CACHE BOOL "" FORCE)
    #set(COMISO_ENABLE_
    FetchContent_MakeAvailable(comiso)
    #if (NOT TARGET CoMISo::CoMISo)
    #    message(FATAL_ERROR "CoMISo target missing")
    #endif()
endif()

if(NOT TARGET CLI11::CLI11)
    FetchContent_Declare(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        cbd58a3696887b34c70949aef21a71735a0c2ad5 # v2.7.2
        SOURCE_DIR "${EXTERNAL_DIR}/CLI11"
        GIT_SUBMODULES ""
        )
    FetchContent_MakeAvailable(cli11)
endif()

# important: TS3D before MC3D
if(NOT TARGET TS3D::TS3D)
    set(TS3D_SUBMODULES_MANUAL TRUE  CACHE BOOL "" FORCE)
    FetchContent_Declare(ts3d
        GIT_REPOSITORY https://github.com/HendrikBrueckler/TrulySeamless3D
        #GIT_REPOSITORY https://github.com/cgg-bern/TrulySeamless3D
        GIT_TAG        218a4b1ff858eff2ce78e7a0d4b7f8e4f00cc720 # main branch on 2026-10-21
        SOURCE_DIR "${EXTERNAL_DIR}/TrulySeamless3D"
        GIT_SUBMODULES ""
        GIT_SUBMODULES_RECURSE NO
        )
    FetchContent_MakeAvailable(ts3d)
endif()

# important: glog before MC3D
if (NOT TARGET glog::glog)
    set(WITH_GTEST FALSE CACHE BOOL "" FORCE)
    set(WITH_GFLAGS FALSE CACHE BOOL "" FORCE)
    set(BUILD_EXAMPLES FALSE CACHE BOOL "" FORCE)
    set(BUILD_TESTING FALSE CACHE BOOL "" FORCE)
    FetchContent_Declare(glog
        GIT_REPOSITORY https://github.com/google/glog.git
        GIT_TAG        53d58e4531c7c90f71ddab503d915e027432447a # newest commit (archived repo)
        SOURCE_DIR "${EXTERNAL_DIR}/glog"
        GIT_SUBMODULES ""
        GIT_SUBMODULES_RECURSE NO
        )
    FetchContent_MakeAvailable(glog)
endif()

# important: MC3D before QGP3D
if(NOT TARGET MC3D::MC3D)
    set(MC3D_SUBMODULES_MANUAL TRUE  CACHE BOOL "" FORCE)
    FetchContent_Declare(mc3d
        GIT_REPOSITORY https://github.com/HendrikBrueckler/MC3D
        #GIT_REPOSITORY https://github.com/cgg-bern/MC3D
        GIT_TAG        29c321fad52aec2371f035f51ccd2d1e9f0f2c64 # main branch on 2026-10-21
        SOURCE_DIR "${EXTERNAL_DIR}/MC3D"
        GIT_SUBMODULES ""
        GIT_SUBMODULES_RECURSE NO
        )
    FetchContent_MakeAvailable(mc3d)
endif()


if(NOT TARGET QGP3D::QGP3D)
    set(QGP3D_SUBMODULES_MANUAL TRUE  CACHE BOOL "" FORCE)
    FetchContent_Declare(qgp3d
        GIT_REPOSITORY https://github.com/HendrikBrueckler/QGP3D
        #GIT_REPOSITORY https://github.com/cgg-bern/QGP3D
        GIT_TAG        1f15db94be38fbceeea9ad0fe29386916a321666 # main branch on 2026-10-21
        SOURCE_DIR "${EXTERNAL_DIR}/QGP3D"
        GIT_SUBMODULES ""
        GIT_SUBMODULES_RECURSE NO
        )
    FetchContent_MakeAvailable(qgp3d)
endif()


if (NOT TARGET nlohmann_json::nlohmann_json)
    set(JSON_Install TRUE  CACHE BOOL "" FORCE)
    FetchContent_Declare(nlohmann_json
        #GIT_REPOSITORY  https://github.com/nlohmann/json.git
        #GIT_TAG v3.12.0
        URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
        URL_HASH SHA256=42f6e95cad6ec532fd372391373363b62a14af6d771056dbfc86160e6dfff7aa
        SOURCE_DIR "${EXTERNAL_DIR}/nlohmann_json"
        )
    FetchContent_MakeAvailable(nlohmann_json)
    message("Downloaded nlohmann-json to ${nlohmann_json_SOURCE_DIR}")
endif()

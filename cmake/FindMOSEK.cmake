# Once done this will (on success) define the following:
#       MOSEK::MosekC    - only the C interface
#       MOSEK::FusionCXX - C and C++ interface

# The enviroment variable MOSEK_DIR can be used to initialize MOSEK_BASE

# Tested with Mosek 9.2 on Linux

find_path (MOSEK_BASE
    NAMES tools/platform
    PATHS
    $ENV{MOSEK_DIR}
    DOC "Base path of your MOSEK installation")

if (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Windows")
    set(DEFAULT_MOSEK_PLATFORM "win64x86")
elseif (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Darwin")
    set(DEFAULT_MOSEK_PLATFORM "osx64x86")
elseif (${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")
    set(DEFAULT_MOSEK_PLATFORM "linux64x86")
else()
    set(DEFAULT_MOSEK_PLATFORM "unknown")
endif()

set (MOSEK_PLATFORM "${DEFAULT_MOSEK_PLATFORM}" CACHE STRING "mosek platform, e.g. linux64x86 or osx64x86")
mark_as_advanced(MOSEK_PLATFORM)


find_path (MOSEK_INCLUDE_DIR
           NAMES mosek.h
           PATHS "${MOSEK_BASE}/tools/platform/${MOSEK_PLATFORM}/h"
          )

find_path (MOSEK_LIBRARY_DIR
           NAMES libmosek64.dylib
                 libmosek64.so
                 mosek64_9_2.dll
           PATHS "${MOSEK_BASE}/tools/platform/${MOSEK_PLATFORM}/bin")

find_library (MOSEK_LIBRARY
              NAMES mosek64
              PATHS "${MOSEK_LIBRARY_DIR}")

find_library (MOSEK_CXX_LIBRARY
              NAMES fusion64
              PATHS "${MOSEK_LIBRARY_DIR}")


if(MOSEK_LIBRARY AND NOT TARGET MOSEK::MosekC)
    add_library(MOSEK::MosekC SHARED IMPORTED)
    target_include_directories(MOSEK::MosekC INTERFACE ${MOSEK_INCLUDE_DIR})
    set_target_properties(MOSEK::MosekC PROPERTIES IMPORTED_LOCATION ${MOSEK_LIBRARY})
endif()

# TODO: add fusion_cxx target to compile with the current compiler(!), safer than a globally built fusion_cxx
#
if(MOSEK_LIBRARY AND MOSEK_CXX_LIBRARY AND NOT TARGET MOSEK::FusionCXX)
    add_library(MOSEK::FusionCXX SHARED IMPORTED)
    target_include_directories(MOSEK::FusionCXX INTERFACE ${MOSEK_INCLUDE_DIR})
    set_target_properties(MOSEK::FusionCXX PROPERTIES IMPORTED_LOCATION ${MOSEK_CXX_LIBRARY})
    target_link_libraries(MOSEK::FusionCXX INTERFACE MOSEK::MosekC)
endif()



# legacy support:
set(MOSEK_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})
set(MOSEK_LIBRARIES MOSEK::MosekC MOSEK::FusionCXX)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MOSEK  DEFAULT_MSG MOSEK_LIBRARY MOSEK_CXX_LIBRARY MOSEK_INCLUDE_DIR)


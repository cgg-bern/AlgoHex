#/*===========================================================================*\
#*        AlgoHex Library --- https://github.com/cgg-bern/AlgoHex.git        *
#*         Copyright (C) 2019-2023 David Bommes, University of Bern          *
#*---------------------------------------------------------------------------*
#*  For license information see LICENSE.txt in the AlgoHex root directory.   *
#*  All contributors are stated in CREDITS.txt.                              *
#\*===========================================================================*/

cmake_minimum_required(VERSION 3.9)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(HexMeshing
        ./HexMeshing.cc
        ./HexMeshing.hh
        ./main.cc)

target_link_libraries(HexMeshing
        AlgoHex
        CLI11::CLI11)

install(TARGETS HexMeshing
    CONFIGURATIONS RelWithDebInfo Release
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        )

set_target_properties(HexMeshing
    PROPERTIES
        BUILD_RPATH_USE_ORIGIN TRUE
        INSTALL_RPATH_USE_LINK_PATH TRUE
        )

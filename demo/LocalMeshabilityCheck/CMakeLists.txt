#/*===========================================================================*\
#*                                                                           *
#*                            HexMeshing Library                             *
#*     Copyright (C) 2019 by Computer Graphics Group, University of Bern     *
#*                                                                           *
#*---------------------------------------------------------------------------*
#*  This file is part of HexMeshing Library.                                 *
#*                                                                           *
#*  HexMeshing is free software: you can redistribute it and/or modify       *
#*  it under the terms of the GNU Lesser General Public License as           *
#*  published by the Free Software Foundation, either version 3 of           *
#*  the License, or (at your option) any later version with the              *
#*  following exceptions:                                                     *
#*                                                                           *
#*  HexMeshing is distributed in the hope that it will be useful,            *
#*  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
#*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
#*  GNU Lesser General Public License for more details.                      *
#*                                                                           *
#*  You should have received a copy of the GNU LesserGeneral Public          *
#*  License along with HexMeshing.  If not,                              *
#*  see <http://www.gnu.org/licenses/>.                                      *
#*                                                                           *
#\*===========================================================================*/

cmake_minimum_required(VERSION 3.29)
project(AlgoHexLocalMeshabilityCheck)

if (NOT TARGET AlgoHex::AlgoHex)
    find_package (AlgoHex)
endif()

# disable OpenVolumeMesh debug output
#add_definitions(-DNDEBUG)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED YES)
endif()

add_executable(LocalMeshabilityCheck ./main.cc)

target_link_libraries(LocalMeshabilityCheck
        AlgoHex
        CLI11::CLI11)

# Precompiled header for the heavy third-party header stack.
# Clang only: GCC's PCH reuse measured slower than no-PCH
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
target_precompile_headers(LocalMeshabilityCheck
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../../src/AlgoHex/AlgoHexPCH.hpp)
endif ()

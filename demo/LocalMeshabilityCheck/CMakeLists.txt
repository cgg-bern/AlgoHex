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

cmake_minimum_required(VERSION 3.14)

#find_package (AlgoHex)

# disable OpenVolumeMesh debug output
#add_definitions(-DNDEBUG)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_executable(LocalMeshabilityCheck ./main.cc)

target_link_libraries(LocalMeshabilityCheck
        AlgoHex
        CLI11::CLI11)


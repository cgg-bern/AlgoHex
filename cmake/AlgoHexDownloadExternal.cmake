################################################################################
include(DownloadProject)

if(NOT DEFINED ALGOHEX_DOWNLOAD_MISSING_DEPS)
    if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
        set(ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT ON)
    else()
        set(ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT OFF)
    endif()
else()
    set(ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT ${ALGOHEX_DOWNLOAD_MISSING_DEPS})
endif()
option(ALGOHEX_DOWNLOAD_MISSINGS_DEPS "Download missing AlgoHex dependencies to external/ folder" ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})


# Shortcut function
function(algohex_download_project name)
    if (NOT ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS})
        message(WARNING "AlgoHex: not downloading missing dependency '${name}', because ALGOHEX_DOWNLOAD_MISSINGS_DEPS is not set")
        return()
    endif()
    message(STATUS "AlgoHex: downloading missing dependency '${name}'")

    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${ALGOHEX_EXTERNAL}/${name}
        DOWNLOAD_DIR ${ALGOHEX_EXTERNAL}/.cache/${name}
        ${ARGN}
    )
endfunction()

################################################################################
## OpenVolumeMesh
function(algohex_download_openvolumemesh)
    algohex_download_project(OpenVolumeMesh
            GIT_REPOSITORY       https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git
            GIT_TAG              master
        )
    if (${ALGOHEX_USE_LOCAL_OVM})
        add_subdirectory(${ALGOHEX_EXTERNAL}/OpenVolumeMesh)
    endif()
endfunction()


## HexEx
function(algohex_download_hexex)
    algohex_download_project(libHexEx
            GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/HexEx/libHexEx.git
            GIT_TAG        algohex
        )
    if (${ALGOHEX_USE_LOCAL_HEXEX})
        add_subdirectory(${ALGOHEX_EXTERNAL}/libHexEx)
    endif()
endfunction()

## qgp3d
function(algohex_download_qgp3d)
    algohex_download_project(QGP3D
            GIT_REPOSITORY https://github.com/HendrikBrueckler/QGP3D
            GIT_TAG        main
            )
    if (${ALGOHEX_USE_LOCAL_QGP3D})
        add_subdirectory(${ALGOHEX_EXTERNAL}/QGP3D)
    endif()
endfunction()

## Eigen
function(algohex_download_eigen)
    algohex_download_project(eigen
            GIT_REPOSITORY      https://gitlab.com/libeigen/eigen.git
            GIT_TAG             master
    )
    add_library (Eigen3::Eigen INTERFACE IMPORTED)
    target_include_directories (Eigen3::Eigen SYSTEM INTERFACE $<BUILD_INTERFACE:${ALGOHEX_EXTERNAL}/eigen>)
    #    target_include_directories (Eigen3::Eigen SYSTEM INTERFACE "${ALGOHEX_EXTERNAL}/eigen")
    target_compile_definitions(Eigen3::Eigen INTERFACE -DEIGEN_HAS_STD_RESULT_OF=0)
endfunction()


## GMM
function(algohex_download_gmm)
    algohex_download_project(gmm
            URL           http://download-mirror.savannah.gnu.org/releases/getfem/stable/gmm-5.4.tar.gz
            URL_HASH SHA256=7163d5080efbe6893d1950e4b331cd3e9160bb3dcf583d206920fba6af7b1e56
            )
        if (${ALGOHEX_USE_LOCAL_GMM})
            add_library(GMM::GMM INTERFACE IMPORTED)
            target_include_directories(GMM::GMM INTERFACE $<BUILD_INTERFACE:${ALGOHEX_EXTERNAL}/gmm/include>)
    endif()
endfunction()

## TinyAD
function(algohex_download_tinyad)
    algohex_download_project(tinyad
            GIT_REPOSITORY           https://github.com/patr-schm/TinyAD.git
            GIT_TAG                  main
            )
    if (${ALGOHEX_USE_LOCAL_TINYAD})
        add_library(tinyad::tinyad INTERFACE IMPORTED)
        target_include_directories(tinyad::tinyad INTERFACE
            $<BUILD_INTERFACE:${ALGOHEX_EXTERNAL}/tinyad/include>)
    endif()
endfunction()


## CoMISo
function(algohex_download_comiso)
    algohex_download_project(CoMISo
            GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/CoMISo/CoMISo.git
            GIT_TAG        0dde6f43f124fbb32c844be399293209dd5ca855 # cgg2 branch, 2023-11-20
   )
    if (${ALGOHEX_USE_LOCAL_COMISO})
        add_subdirectory(${ALGOHEX_EXTERNAL}/CoMISo)
    endif()
endfunction()


## CLI11
function(algohex_download_cli11)
    algohex_download_project(cli11
            GIT_REPOSITORY    https://github.com/CLIUtils/CLI11.git
            GIT_TAG           v2.0.0
    )
    if (${ALGOHEX_USE_LOCAL_CLI11})
        add_subdirectory(${ALGOHEX_EXTERNAL}/cli11)
    endif()
endfunction()

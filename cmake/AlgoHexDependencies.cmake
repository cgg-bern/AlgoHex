################################################################################
# CMake download helpers
################################################################################

# download external dependencies
include(AlgoHexDownloadExternal)

################################################################################
# Required dependencies
################################################################################

option(ALGOHEX_USE_LOCAL_OVM    "Use OVM from external/"    ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})
option(ALGOHEX_USE_LOCAL_HEXEX  "Use HexEX from external/"  ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})
option(ALGOHEX_USE_LOCAL_EIGEN  "Use EIGEN from external/"  ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})
option(ALGOHEX_USE_LOCAL_GMM    "Use GMM from external/"    ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})
option(ALGOHEX_USE_LOCAL_TINYAD "Use TinyAD from external/" ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})
option(ALGOHEX_USE_LOCAL_COMISO "Use CoMISo from external/" ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})
option(ALGOHEX_USE_LOCAL_CLI11  "Use CLI11 from external/"  ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})
option(ALGOHEX_USE_LOCAL_QGP3D  "Use QGP3D from external/"  ${ALGOHEX_DOWNLOAD_MISSINGS_DEPS_DEFAULT})

find_package(OpenVolumeMesh QUIET)
if(NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    algohex_download_openvolumemesh()
endif()


# HexEx
if(NOT TARGET HexEx)
    algohex_download_hexex()
endif()


# Eigen
find_package(EIGEN3)
if(NOT TARGET Eigen3::Eigen)
	algohex_download_eigen()
endif()


# GMM
find_package(GMM)
if(NOT GMM_FOUND)
	algohex_download_gmm()
endif()

# TINYAD
find_package(TINYAD)
if(NOT TINYAD_FOUND)
    algohex_download_tinyad()
endif()

# CoMiSo
if(NOT TARGET CoMISo)
    algohex_download_comiso()
endif()


# CLI11
if(NOT TARGET CLI11::CLI11)
	algohex_download_cli11()
endif()

# QGP3D
if(NOT TARGET QGP3D)
    algohex_download_qgp3d()
endif()
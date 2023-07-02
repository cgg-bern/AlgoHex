# - Try to find TINYAD
# Once done this will define
#  TINYAD_FOUND          - System has TINYAD
#  TINYAD_INCLUDE_DIRS   - The TINYAD include directories

if (TINYAD_INCLUDE_DIR)
  # in cache already
  set(TINYAD_FOUND TRUE)
  set(TINYAD_INCLUDE_DIRS "${TINYAD_INCLUDE_DIR}" )
else (TINYAD_INCLUDE_DIR)

# Check if the base path is set
if ( NOT CMAKE_WINDOWS_LIBS_DIR )
  # This is the base directory for windows library search used in the finders we shipp.
  set(CMAKE_WINDOWS_LIBS_DIR "c:/libs" CACHE STRING "Default Library search dir on windows." )
endif()


find_path( TINYAD_INCLUDE_DIR 
           NAMES TinyAD/Scalar.hh
           PATHS $ENV{TINYAD_DIR}
                 $ENV{TINYAD_DIR}/include
                 /usr/include
                 /usr/include
                 /usr/local/include
          )

set(TINYAD_INCLUDE_DIRS "${TINYAD_INCLUDE_DIR}" )


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TINYAD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TINYAD  DEFAULT_MSG
                                  TINYAD_INCLUDE_DIR)

mark_as_advanced(TINYAD_INCLUDE_DIR)

endif(TINYAD_INCLUDE_DIR)


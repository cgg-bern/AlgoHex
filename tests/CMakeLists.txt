set(TARGET_NAME "${ALGOHEX_TARGET_PREFIX}unittests")

enable_testing()
find_package(GTest QUIET)

if (NOT TARGET GTest::gtest)
    message(STATUS "AlgoHex tests enabled, but GTest not found - downloading using FetchContent.")
    include(FetchContent)

    set(FETCHCONTENT_QUIET OFF)
    set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../external")

    FetchContent_Declare(googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        v1.18.x
        EXCLUDE_FROM_ALL # do not install
        SOURCE_DIR "${EXTERNAL_DIR}/googletest"
    )
    if(WIN32)
        # avoid linking errors, cf https://stackoverflow.com/questions/12540970/how-to-make-gtest-build-mdd-instead-of-mtd-by-default-using-cmake
        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    endif()
    FetchContent_MakeAvailable(googletest)
    # gtest targets have their output directories explicitly set to something different
    # than what we use. Fix it so the tests run on windows (without rpath).
    foreach (_target gtest gtest_main)
        set_target_properties(${_target} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
        )
    endforeach()
endif()

if (TARGET GTest::gtest)
    include(GoogleTest)

    set(TEST_SOURCES
            SphericalHarmonicsTest.cc
            SHProjectorTest.cc
            ConstrainedSHProjectorTest.cc
            )

    if (ALGOHEX_BUILD_EXPENSIVE_TESTS)
        # HexMeshingTest replicates the demo pipeline, so compile the
        # actual demo source (HexMeshing.cc, not main.cc) into the test.
        list(APPEND TEST_SOURCES
                HexMeshingTest.cc
                ${CMAKE_SOURCE_DIR}/demo/HexMeshing/HexMeshing.cc
                ${CMAKE_SOURCE_DIR}/demo/HexMeshing/HexMeshing.hh
                )
    endif ()

    add_executable(${TARGET_NAME} ${TEST_SOURCES})
    target_link_libraries(${TARGET_NAME}
            AlgoHex::AlgoHex
            GTest::gtest
            GTest::gtest_main
            )

    if (ALGOHEX_BUILD_EXPENSIVE_TESTS)
        # Input tetmesh for the HexMeshing pipeline test
        target_compile_definitions(${TARGET_NAME} PRIVATE
                HEXMESHING_TEST_CYLINDER_OVM_PATH="${CMAKE_SOURCE_DIR}/demo/HexMeshing/cylinder.ovm"
                )
    endif ()

    gtest_add_tests(TARGET ${TARGET_NAME}
            EXTRA_ARGS "--gtest_color=yes;--gtest_output=xml"
            )

    set(TEST_DIR "${CMAKE_BINARY_DIR}/Build/tests")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/TestFiles/"
            DESTINATION "${TEST_DIR}"
            )

    set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${TEST_DIR}")
else ()
    message(WARNING "Unable to find/compile gtest, building without tests.")
endif ()

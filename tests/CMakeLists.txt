set(TARGET_NAME "${ALGOHEX_TARGET_PREFIX}unittests")

enable_testing()
find_package(GTest)

if (GTEST_FOUND)
    include(GoogleTest)

    add_executable(${TARGET_NAME}
            SphericalHarmonicsTest.cc
            SHProjectorTest.cc
            ConstrainedSHProjectorTest.cc
            )
    target_link_libraries(${TARGET_NAME}
            AlgoHex::AlgoHex
            GTest::GTest
            GTest::Main
            )

    gtest_add_tests(TARGET ${TARGET_NAME}
            EXTRA_ARGS "--gtest_color=yes;--gtest_output=xml"
            )

    set(TEST_DIR "${CMAKE_BINARY_DIR}/Build/tests")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/TestFiles/"
            DESTINATION "${TEST_DIR}"
            )

    set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${TEST_DIR}")
else ()
    message(WARNING "GoogleTest not found, not building AlgoHex tests.")
endif ()

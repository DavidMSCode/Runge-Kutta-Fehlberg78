cmake_minimum_required(VERSION 3.20)
# Set the project name and language
project(RKF78 LANGUAGES CXX)

# Set the C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the source files for the RKF78 library
add_library(RKF78 INTERFACE)


# Include the header files for the RKF78 library
target_include_directories(RKF78 INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_executable(RKF78_example
    src/test.cpp
)
target_link_libraries(RKF78_example PRIVATE RKF78)
target_include_directories(RKF78_example PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Tests
enable_testing()

add_executable(rkf78_tests
    tests/test_rkf78.cpp
)
target_link_libraries(rkf78_tests PRIVATE RKF78)
target_include_directories(rkf78_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
add_test(NAME rkf78_tests COMMAND rkf78_tests)

# Documentation (Doxygen)
find_package(Doxygen QUIET)
if (DOXYGEN_FOUND)
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM)
else()
    message(STATUS "Doxygen not found. 'docs' target will be unavailable.")
endif()
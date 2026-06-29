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
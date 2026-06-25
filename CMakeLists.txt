cmake_minimum_required(VERSION 3.20)
# Set the project name and language
project(RKF78 LANGUAGES CXX)

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the source files for the RKF78 library
add_library(RKF78 STATIC
    src/RKF78.cpp
)

# Include the header files for the RKF78 library
target_include_directories(RKF78 PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
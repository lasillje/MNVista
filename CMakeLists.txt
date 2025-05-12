cmake_minimum_required(VERSION 3.14)
project(MNVista LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set(CMAKE_BUILD_TYPE Release)

# Find the HTSlib include directory by searching for the header file.
# Adjust "htslib/hts.h" if your header is located elsewhere.
find_path(HTSLIB_INCLUDE_DIR
  NAMES htslib/hts.h
  PATHS ENV(HTSLIB_ROOT) ${HTSLIB_INCLUDE_DIR} /usr/include /usr/local/include
)

# Find the HTSlib library (typically named hts)
find_library(HTSLIB_LIBRARY
  NAMES hts
  PATHS ENV(HTSLIB_ROOT) /usr/lib /usr/local/lib
)

# Ensure both HTSLIB_INCLUDE_DIR and HTSLIB_LIBRARY were found
if (NOT HTSLIB_INCLUDE_DIR)
  message(FATAL_ERROR "Could not find the HTSlib include directory!")
endif()

if (NOT HTSLIB_LIBRARY)
  message(FATAL_ERROR "Could not find the HTSlib library!")
endif()

message(STATUS "Found HTSlib include directory: ${HTSLIB_INCLUDE_DIR}")
message(STATUS "Found HTSlib library: ${HTSLIB_LIBRARY}")

# Add HTSlib include directories
include_directories(${HTSLIB_INCLUDE_DIR})

#Add argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

include_directories(src)
file(GLOB SOURCES "src/*.cpp")

add_executable(MNVista ${SOURCES})

set_target_properties(MNVista PROPERTIES OUTPUT_NAME "MNVista" RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

set_target_properties(MNVista PROPERTIES 
    BUILD_RPATH "/apps/software/OpenSSL/1.1/lib;/apps/software/HTSlib/1.19.1-GCCcore-11.3.0/lib"
    INSTALL_RPATH "/apps/software/OpenSSL/1.1/lib;/apps/software/HTSlib/1.19.1-GCCcore-11.3.0/lib"
)

#Link HTSlib, argparse
find_package(OpenSSL REQUIRED)
target_link_libraries(MNVista ${HTSLIB_LIBRARY} OpenSSL::Crypto argparse)
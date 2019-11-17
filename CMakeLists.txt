# BSD 3-Clause License; see https://github.com/jpivarski/awkward-1.0/blob/master/LICENSE

cmake_minimum_required(VERSION 2.8.12.2)
project(fastjet_binding)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#file(READ "VERSION_INFO" VERSION_INFO)
#string(STRIP ${VERSION_INFO} VERSION_INFO)
#add_definitions(-DVERSION_INFO="${VERSION_INFO}")

# Use RPath on Mac OS.
set(CMAKE_MACOSX_RPATH 1)

# Path to additonal modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Resolve warning
#cmake_policy(SET CMP0074 NEW)

# Require fastjet
find_package(Boost REQUIRED)
find_package(FASTJET REQUIRED)

message(STATUS "Fastjet libs: ${FASTJET_LIBS}, dirs: ${FASTJET_LIBS_DIR}")

include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
link_directories(${FASTJET_LIBS_DIR} ${BOOST_LIBS_DIR})
add_definitions(${FASTJET_DEFINITIONS})

add_subdirectory(pybind11)
pybind11_add_module(fastjet_binding src/fastjet_binding.cpp)
target_link_libraries(fastjet_binding PRIVATE ${FASTJET_LIBS})

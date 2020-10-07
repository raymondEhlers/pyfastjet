# BSD 3-Clause License; see https://github.com/jpivarski/awkward-1.0/blob/master/LICENSE

cmake_minimum_required(VERSION 3.15)
project(pyfastjet)

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
#find_package(Boost REQUIRED)
find_package(FASTJET REQUIRED)

# Setup for awkward1
execute_process(COMMAND python -m awkward1.config --incdir OUTPUT_VARIABLE AWKWARD_INCLUDE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND python -m awkward1.config --libdir OUTPUT_VARIABLE AWKWARD_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
include_directories(BEFORE "${AWKWARD_INCLUDE}")
# Statically link to awkward. This way, we don't need to worry about updating the LD_LIBRARY_PATH.
find_library(CPU-KERNELS awkward-cpu-kernels-static REQUIRED HINTS ${AWKWARD_LIBRARIES})
find_library(LIBAWKWARD awkward-static REQUIRED HINTS ${AWKWARD_LIBRARIES})

message(STATUS "Fastjet libs: ${FASTJET_LIBS}, dirs: ${FASTJET_LIBS_DIR}")

include_directories(SYSTEM ${FASTJET_INCLUDE_DIR})
link_directories(${FASTJET_LIBS_DIR} ${BOOST_LIBS_DIR})
add_definitions(${FASTJET_DEFINITIONS})

add_subdirectory(pybind11)
pybind11_add_module(_src pyfastjet/src/binding.cpp)
# The fastjet cmake module apparently will drop the gmp argument if it's static. So we need to include
# it explicitly for the library to load successfully. There must be a better way to do this, but it works,
# so fine for now.
target_link_libraries(_src PRIVATE ${FASTJET_LIBS} gmp ${CPU-KERNELS} ${LIBAWKWARD})
# For adding this as a submodule to the package
#set_target_properties(pyfastjet PROPERTIES OUTPUT_NAME pyfastjet/binding)

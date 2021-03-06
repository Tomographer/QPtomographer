cmake_minimum_required(VERSION 3.1)

project(QPtomographer)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type (Release, RelWithDebInfo, Debug)")
endif()
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g" CACHE STRING "C++ flags for relwithdebinfo")

# DEBUG: cause eigen_assert() failures to abort(): immediate crash & core dump (which may
# be more useful for debugging the cause)
#add_definitions(-DTOMOGRAPHERPY_DEBUG_EIGEN_ASSERT_CAUSES_ABORT)


message(WARNING "
  ============================================================================
  This CMake script is only intended for development and debugging. To install
  the `QPtomographer` python package, you should use the setup script, e.g.:
  
    > python setup.py install
  
  Check out the README.md file for more information.
  ============================================================================
")


find_package(SCS REQUIRED)

find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

find_package(Threads REQUIRED)

find_package(TomographerPy 5.4 REQUIRED)


include_directories(${SCS_INCLUDE_DIRS})

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories(${MatIO_INCLUDE_DIR})

include_directories(${TOMOGRAPHER_INCLUDE_DIR})


# === PYTHON PACKAGE ===

# copy the skeleton module & all python sources into custom `QPtomographer`
# directory:
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/QPtomographer)
file(COPY ${CMAKE_SOURCE_DIR}/QPtomographer/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/QPtomographer
  FILES_MATCHING PATTERN *.py)

find_package(pybind11 2.2 REQUIRED)

# === Python interface -- pydnormbistates
pybind11_add_module(bistates MODULE cxx/pydnormbistates.cxx NO_EXTRAS)

set_target_properties(bistates PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/QPtomographer")

target_include_directories(bistates PUBLIC  ${SCS_INCLUDE_DIRS})
target_link_libraries(bistates PUBLIC ${SCS_DIR_LIBRARIES})
target_link_libraries(bistates PUBLIC TomographerPy::TomographerPy Threads::Threads)


# === Python interface -- pyfigofmerit
pybind11_add_module(figofmerit MODULE cxx/pyfigofmerit.cxx NO_EXTRAS)

set_target_properties(figofmerit PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/QPtomographer")

target_include_directories(figofmerit PUBLIC  ${SCS_INCLUDE_DIRS})
target_link_libraries(figofmerit PUBLIC ${SCS_DIR_LIBRARIES})
target_link_libraries(figofmerit PUBLIC TomographerPy::TomographerPy)


# === Python interface -- pydnormchannelspace
pybind11_add_module(channelspace MODULE cxx/pydnormchannelspace.cxx NO_EXTRAS)

set_target_properties(channelspace PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/QPtomographer")

target_include_directories(channelspace PUBLIC  ${SCS_INCLUDE_DIRS})
target_link_libraries(channelspace PUBLIC ${SCS_DIR_LIBRARIES})
target_link_libraries(channelspace PUBLIC TomographerPy::TomographerPy Threads::Threads)




# === Enable Testing.

enable_testing()

add_subdirectory(test)

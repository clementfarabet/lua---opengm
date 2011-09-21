cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
project(OpenGM)

find_package(HDF5)
if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    add_definitions(${HDF5_CPPFLAGS})
    message(STATUS "HDF5 found.")
else()
  message(STATUS "HDF5 *NOT* found.")
endif()

find_package(Boost 1.41)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
    message(STATUS "Boost found.")
else()
    message(STATUS "Boost *NOT* found.")
endif()

# disable deprecation warnings
if(WIN32)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
endif()

include_directories(include)
add_subdirectory(src)

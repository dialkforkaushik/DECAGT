cmake_minimum_required(VERSION 2.0)

project(decagt)

IF(NOT DEFINED MULTICORE OR MULTICORE)
    ADD_DEFINITIONS(-DMULTICORE)
ENDIF(MULTICORE)
UNSET(MULTICORE CACHE)

SET(PYTHON true)
ADD_DEFINITIONS(-DPYTHON)

SET(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fopenmp")

# Set source directory
SET(SOURCE_DIR "src/core")

# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})
file(GLOB SOURCES
    ${CMAKE_SOURCE_DIR}/src/core/*.h
    ${CMAKE_SOURCE_DIR}/src/core/*.cc
)

include_directories(${CMAKE_SOURCE_DIR}/lib)
link_directories(decagt ${CMAKE_SOURCE_DIR}/lib)

add_subdirectory(lib/pybind11)
pybind11_add_module(decagt ${SOURCES} "${SOURCE_DIR}/bindings.cc")

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(prefix='', plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE
  OUTPUT_STRIP_TRAILING_WHITESPACE)


if(NOT (${CMAKE_CURRENT_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR}))
add_custom_command(TARGET decagt
           POST_BUILD
           COMMAND cp -r ${CMAKE_CURRENT_BINARY_DIR}/*.so ${CMAKE_SOURCE_DIR})
endif()

install(TARGETS decagt
  COMPONENT python3
  LIBRARY DESTINATION "${PYTHON_SITE}"
  ARCHIVE DESTINATION "${PYTHON_SITE}")



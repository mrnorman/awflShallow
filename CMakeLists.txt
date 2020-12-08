cmake_minimum_required(VERSION 3.0)
project(driver)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

set(DRIVER_SRC driver.cpp)

set(YAKL_HOME ${CMAKE_CURRENT_SOURCE_DIR}/YAKL)
set(YAKL_BIN  ${CMAKE_CURRENT_BINARY_DIR}/yakl)

add_subdirectory(${YAKL_HOME} ${YAKL_BIN})

# Main driver
add_executable(driver ${DRIVER_SRC})
include_directories(${YAKL_HOME})
include_directories(${YAKL_BIN})
target_link_libraries(driver yakl ${NCFLAGS} -lyaml-cpp)

set_source_files_properties(${DRIVER_SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
if ("${ARCH}" STREQUAL "CUDA")
  set_source_files_properties(${DRIVER_SRC} PROPERTIES LANGUAGE CUDA)
  include_directories(${YAKL_HOME}/cub)
endif()


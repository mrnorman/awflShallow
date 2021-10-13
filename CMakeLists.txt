cmake_minimum_required(VERSION 3.0)
project(driver)

enable_language(Fortran)
enable_language(CXX)
enable_language(C)
if ("${YAKL_ARCH}" STREQUAL "CUDA")
  enable_language(CUDA)
endif()

set(DRIVER_SRC driver.cpp)

set(YAKL_HOME ${CMAKE_CURRENT_SOURCE_DIR}/YAKL)
set(YAKL_BIN  ${CMAKE_CURRENT_BINARY_DIR}/yakl)
include_directories(${YAKL_BIN})

add_subdirectory(${YAKL_HOME} ${YAKL_BIN})

# Main driver
add_executable(driver ${DRIVER_SRC})

include(${YAKL_HOME}/yakl_utils.cmake)
yakl_process_target(driver)

target_link_libraries(driver ${NCFLAGS} -lyaml-cpp)



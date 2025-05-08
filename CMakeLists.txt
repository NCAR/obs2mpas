################################################################################
# obs2mpas
################################################################################

cmake_minimum_required( VERSION 3.12 )
project( obs2mpas VERSION 1.0.0 LANGUAGES C CXX Fortran )

## Ecbuild integration
find_package( ecbuild 3.3.2 REQUIRED )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()
list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Configuration options
include( GNUInstallDirs )

## Dependencies
find_package( atlas 0.20.2 REQUIRED )
find_package(NetCDF REQUIRED COMPONENTS Fortran C)

add_subdirectory( src )

## Package Config
ecbuild_install_project( NAME ${PROJECT_NAME} )

ecbuild_print_summary()

cmake_minimum_required(VERSION 2.6)

project(log-demons)

# Set runtime, library, and archive output directories
SET ( CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/bin  CACHE PATH "Single Directory for all executables"       )
SET ( CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/lib  CACHE PATH "Single Directory for all dynamic libraries" )
SET ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/lib  CACHE PATH "Single Directory for all static libraries"  )
MESSAGE ( STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" )
MESSAGE ( STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" )
MESSAGE ( STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}" )

option(LOG_DEMONS_BUILD_PLUGIN "Create LogDemons medInria plugin" OFF)
#TODO option(LOG_DEMONS_ENABLE_FFTW "Add FFTW support")

include_directories(${PROJECT_SOURCE_DIR}/dependencies/TCLAP)
set(TCLAP_INCLUDEDIR ${PROJECT_SOURCE_DIR}/dependencies/TCLAP)
set(TCLAP_DIR ${PROJECT_SOURCE_DIR}/dependencies/TCLAP)

add_subdirectory(src)

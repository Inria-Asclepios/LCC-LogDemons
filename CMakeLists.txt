cmake_minimum_required(VERSION 2.6)

project(LogDemons)

# Set runtime, library, and archive output directories
SET ( CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/bin  CACHE PATH "Single Directory for all executables"       )
SET ( CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/lib  CACHE PATH "Single Directory for all dynamic libraries" )
SET ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/lib  CACHE PATH "Single Directory for all static libraries"  )

option(LOG_DEMONS_BUILD_PLUGIN "Create LogDemons medInria plugin" OFF)
#TODO option(LOG_DEMONS_ENABLE_FFTW "Add FFTW support")

include_directories(${PROJECT_SOURCE_DIR}/dependencies/TCLAP)
set(TCLAP_INCLUDEDIR ${PROJECT_SOURCE_DIR}/dependencies/TCLAP)
set(TCLAP_DIR ${PROJECT_SOURCE_DIR}/dependencies/TCLAP)

add_subdirectory(src)






##################################################################
# Setup use file
##################################################################

if(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in)
configure_file( ## Common use file
  ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Use.cmake.in
  ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Use.cmake
  @ONLY IMMEDIATE)
endif()


###############################################################################
# Setup Configuration files
###############################################################################

if(EXISTS ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in)
configure_file( ## Build tree configure file
  ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
  ${${PROJECT_NAME}_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  @ONLY IMMEDIATE)
endif()
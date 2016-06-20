cmake_minimum_required(VERSION 3.1.0)

set(TARGET_NAME LCCLogDemonsPlugin)

## #################################################################
## Setup version numbering
## #################################################################

set(${TARGET_NAME}_VERSION_MAJOR 3 CACHE STRING "${TARGET_NAME} major version number." FORCE)
set(${TARGET_NAME}_VERSION_MINOR 0 CACHE STRING "${TARGET_NAME} minor version number." FORCE)
set(${TARGET_NAME}_VERSION_BUILD 0 CACHE STRING "${TARGET_NAME} build version number." FORCE)
set(${TARGET_NAME}_VERSION
  "${${TARGET_NAME}_VERSION_MAJOR}.${${TARGET_NAME}_VERSION_MINOR}.${${TARGET_NAME}_VERSION_BUILD}")

mark_as_advanced(${TARGET_NAME}_VERSION_BUILD)

string(TOUPPER ${TARGET_NAME} TARGET_NAME_UP)
add_definitions(-D${TARGET_NAME_UP}_VERSION="${${TARGET_NAME}_VERSION}")

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

## #################################################################
## Resolve dependencies
## #################################################################

if (UNIX AND NOT APPLE)
# Don't export unnecessary symbols, otherwise possible crash on exit
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -fvisibility-inlines-hidden")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--exclude-libs,ALL")
endif()

find_package(dtk REQUIRED)
include_directories(${dtk_INCLUDE_DIRS})

find_package(medInria REQUIRED)

include(list_source_files)
include(list_header_directories_to_include)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

SET(ITKIO_LIBRARIES
  ITKIOImageBase
  ITKIOBMP
  ITKIOBioRad
  ITKIOHDF5
  ITKIOGDCM
  ITKIOGE
  ITKIOGIPL
  ITKIOJPEG
  ITKIOLSM
  ITKIOMeta
  ITKIONIFTI
  ITKIONRRD
  ITKIOPNG
  ITKIOStimulate
  ${ITKIOPhilipsREC_LIBRARIES}
  ITKIOTIFF
  ITKIOVTK
  ITKIOTransformBase
  ITKIOTransformInsightLegacy
  ITKIOTransformHDF5
  ITKIOTransformMatlab
  ITKIOMRC
)

find_package(RPI REQUIRED)
include(${RPI_USE_FILE})

## #################################################################
## Input
## #################################################################

list_source_files(${TARGET_NAME}
  .
  )

list_header_directories_to_include(${TARGET_NAME}
  ${${TARGET_NAME}_HEADERS}
  )

## #################################################################
## Build rules
## #################################################################

add_library(${TARGET_NAME} SHARED
  ${${TARGET_NAME}_CFILES}
  ${${TARGET_NAME}_QRC}
  )

## #############################################################################
## include directories.
## #############################################################################

target_include_directories(${TARGET_NAME}
  PRIVATE ${${TARGET_NAME}_INCLUDE_DIRS}
  )
# PRIVATE because plugins are not supposed to be used somewhere else

## #############################################################################
## Link
## #############################################################################

target_link_libraries(${TARGET_NAME}
  Qt5::Core
  dtkLog
  dtkCoreSupport
  medCoreLegacy
  medRegistration
  ${ITKIO_LIBRARIES}
  libLCClogDemons
)

if (PLUGINS_LEGACY_PATH)
  set_property(TARGET ${TARGET_NAME} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${PLUGINS_LEGACY_PATH})
  set_property(TARGET ${TARGET_NAME} PROPERTY LIBRARY_OUTPUT_DIRECTORY ${PLUGINS_LEGACY_PATH})
endif()

## #################################################################
## Install rules
## #################################################################

install(TARGETS ${TARGET_NAME}
  RUNTIME DESTINATION plugins_legacy
  LIBRARY DESTINATION plugins_legacy
  ARCHIVE DESTINATION plugins_legacy)

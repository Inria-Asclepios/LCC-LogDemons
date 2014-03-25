# Check the CMake version
CMAKE_MINIMUM_REQUIRED ( VERSION 2.6 )
IF ( COMMAND CMAKE_POLICY )
    CMAKE_POLICY( SET CMP0003 NEW )
ENDIF ( COMMAND CMAKE_POLICY )


# Project name
PROJECT ( LCC_LOG_DEMONS )

# Set runtime, library, and archive output directories
SET ( CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/bin  CACHE PATH "Single Directory for all executables"       )
SET ( CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/lib  CACHE PATH "Single Directory for all dynamic libraries" )
SET ( CMAKE_ARCHIVE_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/lib  CACHE PATH "Single Directory for all static libraries"  )
MESSAGE ( STATUS "CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" )
MESSAGE ( STATUS "CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" )
MESSAGE ( STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}" )


# Specify the include directories
INCLUDE_DIRECTORIES (
    ${PROJECT_SOURCE_DIR}/TCLAP
    ${PROJECT_SOURCE_DIR}/RegistrationProgrammingInterface
    ${PROJECT_SOURCE_DIR}/LogDemons
    ${PROJECT_SOURCE_DIR}/LCClogDemons
    ${PROJECT_SOURCE_DIR}/SVFLogJacobian
)


# Specify the directory where libraries are found
LINK_DIRECTORIES ( ${PROJECT_BINARY_DIR}/lib )

set(BUILD_SHARED_LIBS FALSE)


# Add sub-directories
ADD_SUBDIRECTORY( LCClogDemons )
#ADD_SUBDIRECTORY( SpatioTemporalLocalLogDemons            )
#ADD_SUBDIRECTORY( RegistrationProgrammingInterface               )
ADD_SUBDIRECTORY( SVFLogJacobian )

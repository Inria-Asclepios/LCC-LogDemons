# Project name
PROJECT( SSDLog_Demons )

# Define the minimum CMake version needed
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )



#======================================
#
# CMake for standard project
#
#======================================

# Check if ITK is found and include useful files
FIND_PACKAGE( ITK )
IF( NOT ITK_FOUND )
    MESSAGE( "Project ${PROJECT_NAME} requires ITK and ITK was not found. ${PROJECT_NAME} will not be built." )
    RETURN()
ENDIF()
INCLUDE( ${ITK_USE_FILE} )


# Set used libraries
SET( LIBRARIES
${ITK_LIBRARIES}
)

# Create old log demons executable
#ADD_EXECUTABLE(        oldExeLogDemons LogDomainDemonsRegistration.cxx )
#TARGET_LINK_LIBRARIES( oldExeLogDemons ${LIBRARIES} )
#SET_TARGET_PROPERTIES( oldExeLogDemons PROPERTIES OUTPUT_NAME "LogDemons" )

# Install rules
#INSTALL( TARGETS oldExeLogDemons RUNTIME DESTINATION bin )



#======================================
#
# CMake for RPI project
#
#======================================

# Check if TCLAP is found and include directory
FIND_PATH(TCLAP_INCLUDEDIR ${SOURCE_BASE_DIR}/TCLAP/tclap/CmdLine.h
  ${SOURCE_BASE_DIR}/TCLAP/tclap/include
  ${SOURCE_BASE_DIR}/TCLAP/
  )
IF(TCLAP_INCLUDEDIR)
SET(TCLAP_FOUND TRUE)
ELSE()
SET(TCLAP_FOUND FALSE)
ENDIF()


IF( NOT TCLAP_FOUND )
    MESSAGE( "Project ${PROJECT_NAME} requires TCLAP and TCLAP was not found. ${PROJECT_NAME} will not be built." )
    RETURN()
    ENDIF()
INCLUDE_DIRECTORIES( ${TCLAP_INCLUDE_DIR} )

# Set used libraries
SET( LIBRARIES
     ${ITK_LIBRARIES}
)

# Create log demons library
ADD_LIBRARY(           libLogDemons rpiLogDemons.cxx )
TARGET_LINK_LIBRARIES( libLogDemons ${LIBRARIES} )
SET_TARGET_PROPERTIES( libLogDemons PROPERTIES OUTPUT_NAME "RPILogDemons" )

# Create log demons executable
ADD_EXECUTABLE(        exeLogDemons rpiLogDemonsExecutable.cxx )
TARGET_LINK_LIBRARIES( exeLogDemons libLogDemons )
SET_TARGET_PROPERTIES( exeLogDemons PROPERTIES OUTPUT_NAME "rpiLogDemons" )

# Install rules
INSTALL( TARGETS libLogDemons exeLogDemons
         RUNTIME DESTINATION bin
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib/static)

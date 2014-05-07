# Project name
PROJECT ( LCC_LOG_DEMONS )


FIND_PACKAGE( ITK )
IF( NOT ITK_FOUND )
    MESSAGE( "Project ${PROJECT_NAME} requires ITK and ITK was not found. ${PROJECT_NAME} will not be built." )
    RETURN()
ENDIF()
INCLUDE( ${ITK_USE_FILE} )

# Select useful libraries
#SET(Libraries
#ITKCommon
#ITKStatistics
#ITKNumerics
#ITKAlgorithms
#)


# Create executable

ADD_LIBRARY(           libLCClogDemons rpiLCClogDemons.cxx )
TARGET_LINK_LIBRARIES ( libLCClogDemons ${ITK_LIBRARIES} )
SET_TARGET_PROPERTIES ( libLCClogDemons PROPERTIES OUTPUT_NAME "RPILCClogDemons" )



ADD_EXECUTABLE        ( exeLCClogDemons rpiLCClogDemonsExecutable.cxx )
TARGET_LINK_LIBRARIES ( exeLCClogDemons libLCClogDemons )
SET_TARGET_PROPERTIES ( exeLCClogDemons PROPERTIES OUTPUT_NAME "rpiLCClogDemonsExe" )


# Install rules
INSTALL( TARGETS libLCClogDemons exeLCClogDemons
         RUNTIME DESTINATION bin
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib/static)

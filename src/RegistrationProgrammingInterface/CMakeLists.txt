# Project name
PROJECT( RPI )

# Define the minimum CMake version needed
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )


# Check if ITK is found and include useful files
FIND_PACKAGE( ITK )
IF( NOT ITK_FOUND )
    MESSAGE( "Project ${PROJECT_NAME} requires ITK and ITK was not found. ${PROJECT_NAME} will not be built." )
    RETURN()
ENDIF()
INCLUDE( ${ITK_USE_FILE} )


# Check if TCLAP is found and include directory
IF( NOT TCLAP_FOUND )
    MESSAGE( "Project ${PROJECT_NAME} requires TCLAP and TCLAP was not found. ${PROJECT_NAME} will not be built." )
    RETURN()
ENDIF()
INCLUDE_DIRECTORIES( ${TCLAP_INCLUDE_DIR} )


# Check if TinyXML is found and include directory
IF( NOT TinyXML_FOUND )
    MESSAGE( "Project ${PROJECT_NAME} requires TinyXML and TinyXML was not found. ${PROJECT_NAME} will not be built." )
    RETURN()
ENDIF()
INCLUDE_DIRECTORIES( ${TinyXML_INCLUDE_DIR} )


# Add MIPS include directory
INCLUDE_DIRECTORIES( ${MIPS_IMAGEPKG_INCLUDE_DIRS} )


# Set used libraries
SET(LIBRARIES
    ITKIO
    ITKAlgorithms
    ITKCommon
    ITKBasicFilters
    ITKStatistics
    ITKNumerics
    mipsItkImageIO
)


# Create ConvertLinearToDF executable
ADD_EXECUTABLE(        exeConvertLinearToDF rpiConvertLinearToDF.cxx )
TARGET_LINK_LIBRARIES( exeConvertLinearToDF ${LIBRARIES} )
SET_TARGET_PROPERTIES( exeConvertLinearToDF PROPERTIES OUTPUT_NAME "rpiConvertLinearToDF" )


# Create ConvertLinearToSVF executable
ADD_EXECUTABLE(        exeConvertLinearToSVF rpiConvertLinearToSVF.cxx )
TARGET_LINK_LIBRARIES( exeConvertLinearToSVF ${LIBRARIES} )
SET_TARGET_PROPERTIES( exeConvertLinearToSVF PROPERTIES OUTPUT_NAME "rpiConvertLinearToSVF" )


# Create rpiFuseTransformations executable
ADD_EXECUTABLE(        exeFuseTransformations rpiFuseTransformations.cxx )
TARGET_LINK_LIBRARIES( exeFuseTransformations ${LIBRARIES} TinyXML )
SET_TARGET_PROPERTIES( exeFuseTransformations PROPERTIES OUTPUT_NAME "rpiFuseTransformations" )


# Create rpiResampleImage executable
ADD_EXECUTABLE(        exeResampleImage rpiResampleImage.cxx )
TARGET_LINK_LIBRARIES( exeResampleImage ${LIBRARIES} )
SET_TARGET_PROPERTIES( exeResampleImage PROPERTIES OUTPUT_NAME "rpiResampleImage" )

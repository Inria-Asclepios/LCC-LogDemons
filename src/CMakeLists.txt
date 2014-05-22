# Specify the include directories
INCLUDE_DIRECTORIES (LogDemons
                     LCClogDemons
                     SVFLogJacobian
)


# Specify the directory where libraries are found
LINK_DIRECTORIES ( ${PROJECT_BINARY_DIR}/lib )

set(BUILD_SHARED_LIBS FALSE)

if (LOG_DEMONS_BUILD_PLUGIN)
  add_subdirectory(medInriaPlugin)
endif()

# If medInria plugin is enabled, this will retrieve the superproject packages
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(RPI REQUIRED)
include(${RPI_USE_FILE})

# Add sub-directories
#ADD_SUBDIRECTORY( LogDemons )
ADD_SUBDIRECTORY( LCClogDemons )
ADD_SUBDIRECTORY( SVFLogJacobian )

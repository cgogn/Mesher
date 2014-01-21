cmake_minimum_required(VERSION 2.8)
project(CGoGN)

SET ( CGoGN_ROOT_DIR ${CMAKE_SOURCE_DIR} )
IF(WIN32)
	SET ( CMAKE_PREFIX_PATH ${CGoGN_ROOT_DIR}/windows_dependencies CACHE STRING "path to dependencies" )
ENDIF(WIN32)

find_package(Boost COMPONENTS regex thread system REQUIRED)
find_package(ZLIB REQUIRED)
find_package(LibXml2 REQUIRED)

# define includes of external libs 
SET (CGoGN_EXT_INCLUDES
	${ZLIB_INCLUDE_DIRS}
	${LIBXML2_INCLUDE_DIR}
	${Boost_INCLUDE_DIRS}
)

# define libs for external libs
SET (CGoGN_EXT_LIBS
	${ZLIB_LIBRARIES}
	${LIBXML2_LIBRARIES}
	${Boost_SYSTEM_LIBRARY}
	${Boost_REGEX_LIBRARY}
	${Boost_THREAD_LIBRARY}
)

IF(WIN32)
	add_subdirectory(Release)
	set(CMAKE_CONFIGURATION_TYPES Release Debug)
	set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING "Only Release or Debug" FORCE)
ELSE(WIN32)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Woverride-init -fPIC")
	add_subdirectory(Release)
	add_subdirectory(Debug)
ENDIF(WIN32)


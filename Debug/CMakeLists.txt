cmake_minimum_required(VERSION 2.8)

project(CGoGNDebug)

SET(CMAKE_BUILD_TYPE Debug)

IF(WIN32)
	SET(LIBRARY_OUTPUT_PATH ${CGoGN_ROOT_DIR}/lib)
ELSE(WIN32)
	SET(LIBRARY_OUTPUT_PATH ${CGoGN_ROOT_DIR}/lib/Debug)
ENDIF(WIN32)


include_directories(
	${CMAKE_CURRENT_BINARY_DIR}
	${CGoGN_ROOT_DIR}/include
	${CGoGN_EXT_INCLUDES}
)


file(
	GLOB_RECURSE
	files_topology
	${CGoGN_ROOT_DIR}/src/Topology/*.cpp
	${CGoGN_ROOT_DIR}/include/Topology/*.hpp
	${CGoGN_ROOT_DIR}/include/Topology/*.h
)

file(
	GLOB_RECURSE
	files_container
	${CGoGN_ROOT_DIR}/src/Container/*.cpp
	${CGoGN_ROOT_DIR}/include/Container/*.hpp
	${CGoGN_ROOT_DIR}/include/Container/*.h
)

file(
	GLOB_RECURSE
	files_algo
	${CGoGN_ROOT_DIR}/src/Algo/*.cpp
	${CGoGN_ROOT_DIR}/src/Algo/*.c
	${CGoGN_ROOT_DIR}/include/Algo/*.hpp
	${CGoGN_ROOT_DIR}/include/Algo/*.h
)



file(
	GLOB		# WARNING NO MORE RECURSE TO AVOID TAKING QT FILES
	files_utils
	${CGoGN_ROOT_DIR}/src/Utils/*.cpp
	${CGoGN_ROOT_DIR}/include/Utils/*.hpp
	${CGoGN_ROOT_DIR}/include/Utils/*.h
)


add_library( cgognD ${files_topology} ${files_container} ${files_algo} ${files_utils} )



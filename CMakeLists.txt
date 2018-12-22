cmake_minimum_required(VERSION 3.9)

set(APPLICATION coopti)
project(${APPLICATION})

# -std=c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -fPIC")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(lib/Eigen)
find_package(OpenGL REQUIRED) # define OPENGL_LIBRARIES

include_directories(
	${CMAKE_SOURCE_DIR}/src
	lib/Eigen
)

add_library(BRDFReader
	src/BRDFReader/brdfreader.cpp
	src/BRDFReader/brdfreader.h
	)
add_library(PCA
	src/PCA/PCASolver.cpp
	src/PCA/PCASolver.h
	)
add_library(Optimisation
	src/Optimisation/OptimisationSolver.cpp
	src/Optimisation/OptimisationSolver.h
	)

add_executable(${APPLICATION}
	src/main.cpp
)

set(EXT_LIBS ${OPENGL_gl_LIBRARY} )
target_link_libraries(${APPLICATION}
	${EXT_LIBS}
	BRDFReader
	PCA
	Optimisation
	)

add_subdirectory(src/BRDFMap)

cmake_minimum_required(VERSION 3.5)

set(APPLICATION coopti)
project(${APPLICATION})

# -std=c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -fPIC -fopenmp")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(lib/stxxl)

include_directories(
	${CMAKE_SOURCE_DIR}/src
	lib/Eigen
	${STXXL_INCLUDE_DIRS}
)

add_subdirectory(src/BRDFReader)
add_subdirectory(src/Parametrisation)
add_subdirectory(src/Optimisation)

add_executable(${APPLICATION}
	src/main.cpp
)
target_link_libraries(${APPLICATION}
	BRDFReader
	Parametrisation
	Optimisation
	${STXXL_LIBRARIES}
)

add_subdirectory(tests)

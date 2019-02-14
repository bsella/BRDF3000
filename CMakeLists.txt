cmake_minimum_required(VERSION 3.5)

set(APPLICATION coopti)
project(${APPLICATION})

# -std=c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -fPIC -fopenmp -O3")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

#add_subdirectory(lib/stxxl)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/lib/Eigen
	${CMAKE_CURRENT_SOURCE_DIR}/lib/bitmap
	${CMAKE_CURRENT_SOURCE_DIR}/lib/boost/libs/multiprecision/include/boost/multiprecision
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/boost/libs/math/include/
	${CMAKE_CURRENT_SOURCE_DIR}/lib/boost/libs/numeric/odeint/include/
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
	stdc++fs
	quadmath
	gmp
)

add_subdirectory(tests)

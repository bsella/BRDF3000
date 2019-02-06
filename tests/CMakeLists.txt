cmake_minimum_required(VERSION 3.5)

set(APPLICATION coopti_tests)
project(${APPLICATION})

# -std=c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/tests)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/../src
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/../lib/Eigen	
)

add_executable(${APPLICATION}
	main.cpp
	BaseTest.h
	BaseTest.cpp
	OptimisationTest.h
	OptimisationTest.cpp
	ParametrisationTest.h
	ParametrisationTest.cpp
	BRDFReaderTest.h
	BRDFReaderTest.cpp
)
target_link_libraries(${APPLICATION}
	BRDFReader
	Parametrisation
	Optimisation
	stdc++fs
)

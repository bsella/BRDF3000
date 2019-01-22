cmake_minimum_required(VERSION 3.5)

project(Parametrisation)

# -std=c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -fPIC")

include_directories(
	${CMAKE_SOURCE_DIR}/src
	lib/Eigen
)

file(GLOB SOURCES
	*.cpp
)
file(GLOB HEADERS
	*.h
)

add_library(Parametrisation
	${SOURCES}
	${HEADERS}
)

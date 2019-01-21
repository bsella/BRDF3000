cmake_minimum_required(VERSION 3.5)
project(BRDFMap)

find_package(OpenGL REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)

# -std=c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")

include_directories(
	${Qt5Widgets_INCLUDE_DIRS}
	${Qt5OpenGL_INCLUDE_DIRS}
	)

add_definitions(${Qt5Widgets_DEFINITIONS})
add_definitions(${Qt5OpenGL_DEFINITIONS})

add_library(BRDFMap
	BRDFMapWidget.h
	BRDFMapWidget.cpp
	)

qt5_use_modules(BRDFMap Widgets OpenGL)

set(EXT_LIBS ${QT_LIBRARIES} ${OPENGL_LIBRARIES})
target_link_libraries(BRDFMap ${EXT_LIBS})

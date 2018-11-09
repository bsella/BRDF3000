cmake_minimum_required(VERSION 3.9)

set(APPLICATION coopti)
project(${APPLICATION})

# -std=c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(lib/Eigen)
find_package(OpenGL REQUIRED) # define OPENGL_LIBRARIES


# Add Sources
FILE(GLOB_RECURSE
    folder_source
    ${CMAKE_SOURCE_DIR}/src/*.cpp
)

FILE(GLOB_RECURSE
    folder_header
    ${CMAKE_SOURCE_DIR}/src/*.h

)

include_directories(
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/glm
    lib/Eigen
)

SOURCE_GROUP("Source Files" FILES ${folder_source})
SOURCE_GROUP("Header Files" FILES ${folder_header})
#SOURCE_GROUP("Shader Files" FILES ${folder_shader})

add_executable(${APPLICATION}
               ${folder_source}
               ${folder_header}
               )

set(EXT_LIBS ${OPENGL_gl_LIBRARY} )
target_link_libraries(${APPLICATION} ${EXT_LIBS} )


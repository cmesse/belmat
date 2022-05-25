cmake_minimum_required(VERSION 3.20)
project(belmat)

set(CMAKE_CXX_STANDARD 11)

# make sure that the library and the demo are
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}  )


include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_subdirectory( src )

add_executable(demo main.cpp)
add_dependencies(demo libbelmat.a )
target_link_libraries(demo libbelmat.a )
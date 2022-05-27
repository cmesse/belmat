cmake_minimum_required(VERSION 3.20)
project(belmat)

set(CMAKE_CXX_STANDARD 11)

option( BUILD_SHARED "Build Shared Library instead of static one" OFF )
option( USE_HDF5 "Link to HDF5 libraries" ON )
option( USE_DEBUG "Use Debug Flags" OFF )

# temporary settings
set( HDF5_INCLUDE_PATH "/opt/tpls/sw_gcc/include" )
set( HDF5_LIBRARY_PATH "/opt/tpls/sw_gcc/lib64/" )

# make sure that the library and the demo are
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

if( USE_DEBUG )
    add_definitions( -DDEBUG )
    set( CMAKE_CXX_FLAGS "-Og -g" )
else()
    set( CMAKE_CXX_FLAGS "-O2" )
endif()

if( USE_HDF5 )
    add_definitions( -DHDF5 )
    include_directories( ${HDF5_INCLUDE_PATH} )
    set( BELMAT_LIBS "-Wl,-rpath,{HDF5_LIBRARY_PATH} -L${HDF5_LIBRARY_PATH} -lhdf5_hl -lhdf5 -lz -ldl" )
endif()



include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src  )
add_subdirectory( src )

add_executable( demo main.cpp )
add_dependencies( demo libbelmat.a )
target_include_directories(demo PUBLIC ${BELMAT_INCLUDES} )
target_link_libraries( demo PUBLIC libbelmat.a ${BELMAT_LIBS} )
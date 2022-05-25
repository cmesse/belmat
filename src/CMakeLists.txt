set( LIBNAME belmat )

set( SOURCES
        database.cpp
        )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )

# create library
add_library( libbelmat.a STATIC ${SOURCES} )
set_target_properties( libbelmat.a PROPERTIES OUTPUT_NAME belmat )
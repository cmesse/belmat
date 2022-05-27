set( LIBNAME belmat )

set( SOURCES
        database.cpp
        Hdf5File.cpp
        )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )

# create library

if( BUILD_SHARED )
    if( APPLE )
        add_library( libbelmat.dylibSHARED ${SOURCES} )
        set_target_properties( libbelmat.dylib PROPERTIES OUTPUT_NAME belmatt )
    else()
        add_library( libbelmat.so SHARED ${SOURCES} )
        set_target_properties( libbelmat.so PROPERTIES OUTPUT_NAME belmat )
    endif()
else()
    add_library( libbelmat.a STATIC ${SOURCES} )
    set_target_properties( libbelmat.a PROPERTIES OUTPUT_NAME belmat )
endif()
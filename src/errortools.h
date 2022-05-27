//
// Created by Christian Messe on 12/6/21.
//

#ifndef BELMAT_HDF5TOOLS_ERROR_HPP
#define BELMAT_HDF5TOOLS_ERROR_HPP

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <memory>

namespace errortools {
//------------------------------------------------------------------------------

    /**
     * converts a variable input into a string
     *
     * /tparam Args
     * /param format
     * /param args
     * /return
     */
    template < typename ... Args >
    std::string sprint( const char * format, const Args ... args )
    {
        // Determine size of string.
        auto size = std::snprintf( nullptr, 0, format, args ... );

        // create unique pointer with length tSize. Add Extra space for '\0'.
        std::unique_ptr< char[] > buffer( new char[ size + 1 ] );

        // write formatted string into buffer
        std::snprintf( buffer.get(), size + 1, format, args ... );

        // return formatted string
        return std::string( buffer.get(), buffer.get() + size );
    }

//------------------------------------------------------------------------------

    /**
     * throw an error if condition is not fulfilled
     * /param check  condition that is to be checked
     * /param args information that is to be printed
     */
    template < typename ... Args > void
    error( const bool condition, const   Args ...      args ) {
        // throw error if condition is not fulfilled
        if( ! condition )
        {
            throw std::runtime_error( sprint( args ... ).c_str() );
        }
    }

//------------------------------------------------------------------------------
}
#endif //BELMAT_HDF5TOOLS_ERROR_HPP

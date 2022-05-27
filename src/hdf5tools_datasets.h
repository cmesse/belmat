//
// Created by Christian Messe on 12/6/21.
//

#ifndef STRUMPACK_HDF5TOOLS_DATASETS_HPP
#define STRUMPACK_HDF5TOOLS_DATASETS_HPP

#include "hdf5tools_datatypes.h"

namespace hdf5tools
{
//------------------------------------------------------------------------------

    /**
     * test if a dataset exists
     */
    inline bool
    dataset_exists( hid_t file_id, const std::string & label )
    {
#ifdef HDF5
        hid_t data_set = 0;
        return H5Lexists( file_id, label.c_str(), data_set );
#else
        return false;
#endif
    }

//------------------------------------------------------------------------------

    /**
     * test if a group exists
     */
    inline bool
    group_exists( hid_t file_id, const std::string & label )
    {
#ifdef HDF5
        std::string l = "/" + label;
        hid_t group = 0;
        return H5Lexists( file_id, l.c_str(), group );
#else
        return false;
#endif
    }

//------------------------------------------------------------------------------
}
#endif // HDF5TOOLS_DATASETS_HPP

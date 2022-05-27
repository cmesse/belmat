//
// Created by Christian Messe on 11/23/21.
//





#ifndef STRUMPACK_HDF5TOOLS_DATATYPES_HPP
#define STRUMPACK_HDF5TOOLS_DATATYPES_HPP

#ifdef HDF5
#include <hdf5.h>
#else
typedef long int hid_t;
typedef int herr_t;
typedef long long unsigned int hsize_t;
#endif

#include "errortools.h"

namespace hdf5tools
{
//------------------------------------------------------------------------------

    /**
     *
     * /brief                 returns an HDF5 enum defining the
     *                        data type that is to be communicated.
     *
     * see also https://support.hdfgroup.org/HDF5/doc/H5.user/Datatypes.html
     */
    template<typename T> hid_t
    hdf5_type() {
        return T::hdf5_type();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef HDF5
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline hid_t
    hdf5_type< bool > () {
        return H5T_NATIVE_HBOOL;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    template<> inline hid_t
    hdf5_type< int >(){
        return H5T_NATIVE_INT ;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline hid_t
    hdf5_type< long int >() {
        return H5T_NATIVE_LONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline hid_t
    hdf5_type < unsigned int> () {
        return H5T_NATIVE_UINT;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline hid_t
    hdf5_type < long unsigned int >() {
        return H5T_NATIVE_ULONG;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline hid_t
    hdf5_type< float >() {
        return H5T_NATIVE_FLOAT;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline hid_t
    hdf5_type< double > () {
        return H5T_NATIVE_DOUBLE;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    template<> inline hid_t
    hdf5_type< long double >() {
        return H5T_NATIVE_LDOUBLE;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#endif // HDF5
//------------------------------------------------------------------------------
}

#endif // STRUMPACK_HDF5TOOLS_DATATYPES_HPP

//
// Created by Christian Messe on 12/6/21.
//

#ifndef STRUMPACK_HDTFTOOLS_WRITE_HPP
#define STRUMPACK_HDTFTOOLS_WRITE_HPP

#include "hdf5tools_datatypes.hpp"
#include "hdf5tools_datasets.hpp"
#include "errortools.hpp"

namespace hdf5tools
{
//----------------------------------------------------------------------

    /**
      * writes a boolean value to a file
      *
      * \param loc_id  handler to group or file
      * \param label   name of dataset
      * \param value   the data
      * \param status  error handler
      */
    inline void
    write_bool(
            hid_t               loc_id,
            const std::string & label,
            const bool          value,
            herr_t            & status ) {
#ifdef HDF5
        // test if dataset exists
        errortools::error( ! dataset_exists( loc_id, label ),
               "dataset %s already exists",
               label.c_str() );

        // select data type for matrix to save
        hid_t data_type = H5Tcopy( H5T_NATIVE_HBOOL );

        // matrix dimensions
        hsize_t dims[ 1 ] = { 1 };

        // create data space
        hid_t data_space = H5Screate_simple( 1, dims, nullptr );

        // create new dataset
        hid_t data_set = H5Dcreate(
                loc_id,
                label.c_str(),
                data_type,
                data_space,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT);

        // value to cast bool to
        hbool_t hvalue = ( hbool_t ) value;

        // write data into dataset
        status = H5Dwrite(
                data_set,
                data_type,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                & hvalue );

        // close resources
        H5Sclose( data_space );
        H5Tclose( data_type );
        H5Dclose( data_set );


        // check for error
        errortools::error( status == 0,
           "error %d occured while trying to write bool %s",
           ( int ) status, label.c_str() );
#endif
    }

//----------------------------------------------------------------------

      /**
      * writes a string to a file
      *
      * \param loc_id  handler to group or file
      * \param label   name of dataset
      * \param data   the data
      * \param status  error handler
      */
    inline void
    write_string(
            hid_t                 loc_id,
            const std::string   & label,
            const std::string   & value,
            herr_t              & status ) {
#ifdef HDF5
          errortools::error( ! dataset_exists( loc_id, label ),
               "dataset %s already exists",
               label.c_str() );

        // select data type for string
        hid_t data_type = H5Tcopy( H5T_C_S1 );

        // set size of output type
        status  = H5Tset_size( data_type, value.length() );

        // matrix dimensions
        hsize_t dims[ 1 ] = { 1 };

        // create data space
        hid_t data_space = H5Screate_simple( 1, dims, nullptr );

        // create new dataset
        hid_t data_set = H5Dcreate(
                loc_id,
                label.c_str(),
                data_type,
                data_space,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // write data into dataset
        status = H5Dwrite(
                data_set,
                data_type,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                value.c_str() );

        // close open hids
        H5Sclose( data_space );
        H5Tclose( data_type );
        H5Dclose( data_set );

        // check for error
          errortools::error( status == 0,
           "error %d occured while trying to write string %s",
           ( int ) status, label.c_str() );
#endif
    }

//----------------------------------------------------------------------

    /**
    * saves a scalar value to a file
    *
    * \param loc_id  handler to hdf5 file
    * \param label   label of matrix to save
    * \param value   value that is to be stored
    * \param status  error handler
    */
    template < typename T >
    void
    write_scalar(
            hid_t                 loc_id,
            const std::string   & label,
            const T               value,
            herr_t              & status ) {
#ifdef  HDF5
        // test if dataset exists
        errortools::error( ! dataset_exists( loc_id, label ),
               "dataset %s already exists",
               label.c_str() );

        // select datatype for data to save
        hid_t data_type = H5Tcopy( hdf5tools::hdf5_type<T>() );

        // set data type to little endian
        status = H5Tset_order( data_type, H5T_ORDER_LE );

        // matrix dimensions
        hsize_t dims[ 1 ] = { 1 };

        // create data space
        hid_t data_space = H5Screate_simple( 1, dims, nullptr );

        // create new dataset
        hid_t data_set = H5Dcreate(
                loc_id,
                label.c_str(),
                data_type,
                data_space,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // write data into dataset
        status = H5Dwrite(
                data_set,
                data_type,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                & value );

        // close open hids
        H5Sclose( data_space );
        H5Tclose( data_type );
        H5Dclose( data_set );

        // check for error
        errortools::error( status == 0,
           "error %d occured while trying to write scalar %s",
           ( int ) status, label.c_str() );
#endif
    }

//----------------------------------------------------------------------

    /**
    * saves a raw array to a file
    *
    * \param loc_id  handler to hdf5 file
    * \param label   label of matrix to save
    * \param value   value that is to be stored
    * \param size    length of array
    * \param status  error handler
    */
    template < typename T >
    void
    write_array(
            hid_t               loc_id,
            const std::string & label,
            const T           * data,
            const hsize_t       size,
            herr_t            & status ) {
#ifdef  HDF5
        // test if dataset exists
        errortools::error( ! dataset_exists( loc_id, label ),
               "dataset %s already exists",
               label.c_str() );

        // select datatype for data to save
        hid_t data_type = H5Tcopy( hdf5tools::hdf5_type<T>() );

        // set data type to little endian
        status = H5Tset_order( data_type, H5T_ORDER_LE );

        // matrix dimensions
        hsize_t dims[ 1 ];
        dims[ 0 ] = size;

        // create data space
        hid_t  data_space
                = H5Screate_simple( 1, dims, nullptr );

        // create new dataset
        hid_t data_set = H5Dcreate(
                loc_id,
                label.c_str(),
                data_type,
                data_space,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT );

        // test if vector is not empty
        if( size > 0 )
        {
            // write data into dataset
            status = H5Dwrite(
                    data_set,
                    data_type,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &data[ 0 ] );
        }

        // close open hids
        H5Sclose( data_space );
        H5Tclose( data_type );
        H5Dclose( data_set );

        // check for error
        errortools::error( status == 0,
           "error %d occured while trying to write array %s",
           ( int ) status, label.c_str() );
#endif
    }

//----------------------------------------------------------------------
}

#endif //STRUMPACK_HDTFTOOLS_WRITE_HPP

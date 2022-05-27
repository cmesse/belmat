//
// Created by Christian Messe on 12/6/21.
//

#ifndef BELMAT_HDFT5OOLS_READ_HPP
#define BELMAT_HDFT5OOLS_READ_HPP

#include "hdf5tools_datatypes.h"
#include "hdf5tools_datasets.h"
#include "errortools.h"


namespace hdf5tools
{
//------------------------------------------------------------------------------

     /**
      * reads a boolean value from a file
      *
      * \param loc_id  handler to group or file
      * \param label   name of dataset
      * \param value   the data
      * \param status  error handler
      */
     inline void
     read_bool(
             hid_t               loc_id,
             const std::string & label,
             bool              & value,
             herr_t            & status )
    {
#ifdef HDF5
        // test if dataset exists
         errortools::error( dataset_exists( loc_id, label ),
               "dataset %s does not exist",
               label.c_str() );

        // open the data set
        hid_t data_set = H5Dopen1( loc_id, label.c_str() );

        // get the data type of the set
        hid_t data_type = H5Dget_type( data_set );

        // get the handler of the data space
        hid_t data_space = H5Dget_space( data_set );

        // value to cast bool to
        hbool_t hvalue ;

        // read data from file
        status = H5Dread(
                data_set,
                data_type,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &hvalue );

         // cast output value
        value = ( bool ) hvalue;

        // close resources
        H5Sclose( data_space );
        H5Tclose( data_type );
        H5Dclose( data_set );


        // check for error
        errortools::error( status == 0,
           "error %d occured while trying to read bool %s",
           ( int ) status, label.c_str() );
#endif
    }

//------------------------------------------------------------------------------

    /**
      * reads a bool value from a file
      *
      * \param loc_id  handler to group or file
      * \param label   name of dataset
      * \param data    the data
      * \param status  error handler
      */
    inline void
    read_string(
            hid_t               loc_id,
            const std::string & label,
            std::string       & data,
            herr_t            & status ) {
#ifdef HDF5
        // test if dataset exists
        errortools::error( dataset_exists( loc_id, label ),
               "dataset %s does not exist",
               label.c_str() );

        // open the data set
        hid_t data_set = H5Dopen1( loc_id, label.c_str() );

        // get the data type of the set
        hid_t data_type = H5Dget_type( data_set );

        // get the handler of the data space
        hid_t data_space = H5Dget_space( data_set );

        // get length of string
        hsize_t length = H5Dget_storage_size( data_set );

        // allocate buffer for string
        char * buffer = ( char * ) malloc( length * sizeof( char ) );

        // load string into buffer
        status = H5Dread(
                data_set,
                data_type,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                buffer );

        // create string from buffer
        data.assign( buffer, length );

        // delete buffer
        free( buffer );

        // close resources
        H5Sclose( data_space );
        H5Tclose( data_type );
        H5Dclose( data_set );

        // check for error
        errortools::error( status == 0,
           "error %d occured while trying to read string %s",
           ( int ) status, label.c_str() );
#endif
    }

//------------------------------------------------------------------------------

      /**
       * reads a scalar from a file
       *
       * \param loc_id  handler to group or file
       * \param label   name of dataset
       * \param value   the data
       * \param status  error handler
       */
     template < typename T > void
     read_scalar(
             hid_t               loc_id,
             const std::string & label,
             T                 & value,
             herr_t            & status ) {
#ifdef HDF5
        // test if dataset exists
        errortools::error( dataset_exists( loc_id, label ),
               "dataset %s does not exist",
               label.c_str() );

        // open the data set
        hid_t data_set = H5Dopen1( loc_id, label.c_str() );

        // get the data type of the set
        hid_t data_type = H5Dget_type( data_set );

        // get the handler of the data space
        hid_t data_space = H5Dget_space( data_set );

        // read data from file
        status = H5Dread(
                data_set,
                data_type,
                H5S_ALL,
                H5S_ALL,
                H5P_DEFAULT,
                &value );

        // close resources
        H5Sclose( data_space );
        H5Tclose( data_type );
        H5Dclose( data_set );

        // check for error
        errortools::error( status == 0,
           "error %d occured while trying to read scalar %s",
           ( int ) status, label.c_str() );
#endif
    }

//------------------------------------------------------------------------------

      /**
      * reads a raw array from a file
      *
      * \param loc_id  handler to group or file
      * \param label   name of dataset
      * \param data    the data
      * \param size    the length of the array
      * \param status  error handler
      */
    template < typename T >
    void
    read_array(
            hid_t               loc_id,
            const std::string & label,
            T                 * data,
            const hsize_t       size,
            herr_t            & status ) {
#ifdef HDF5
        // test if dataset exists
        errortools::error( dataset_exists( loc_id, label ),
               "dataset %s does not exist",
               label.c_str() );


        // check for error
        errortools::error( status == 0,
           "error %d occured while trying to read array %s",
           ( int ) status, label.c_str() );

        // open the data set
        hid_t data_set = H5Dopen1( loc_id, label.c_str() );

        // get the data type of the set
        hid_t data_type = H5Dget_type( data_set );

        // get the handler of the data space
        hid_t data_space = H5Dget_space( data_set );

        // matrix dimensions
        hsize_t dims[ 1 ];

        // ask hdf for dimensions
        status  = H5Sget_simple_extent_dims( data_space, dims, nullptr );

        // check length
        errortools::error( dims[ 0 ] == size, "length of array %s to not match: is %lu, expect %lu",
               ( long unsigned int ) dims[ 0 ], ( long unsigned int ) size );

        // test if array is empty
        if( dims[ 0 ] > 0 )
        {
           // read data from file
            status = H5Dread(
                    data_set,
                    data_type,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    &data[ 0 ] );
        }
        else
        {
            // do nothing
            status = 0;
        }

        // close resources
        H5Sclose( data_space );
        H5Tclose( data_type );
        H5Dclose( data_set );

        // check for error
        errortools::error( status == 0,
           "error %d occured while trying to read array %s",
           ( int ) status, label.c_str() );
#endif
    }
//------------------------------------------------------------------------------
}
#endif //BELMAT_HDFT5OOLS_READ_HPP

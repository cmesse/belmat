//
// Created by Christian Messe on 12/6/21.
//

#ifndef MATREAD_HD5FILE_HPP
#define MATREAD_HD5FILE_HPP

#include <string>

#include "hdf5tools_datatypes.hpp"

//------------------------------------------------------------------------------


//! the file mode is passed to the constructor in order to tell
//! how the file should be handled
enum class Hdf5FileMode
{
    NEW,
    OPEN_RDONLY,
    OPEN_RDWR
};

//------------------------------------------------------------------------------

class Hdf5File
{
    //! path to data file
    std::string mPath ;

    //! HDF5-pointer to the file
    hid_t mFile ;

    //! HDF5-pointer to active group
    hid_t mActiveGroup = -1;

    //! error status of file
    herr_t mStatus = 0 ;

    //! flag telling if the file is still open, needed by destructor
    bool mFileIsOpen = false ;

//------------------------------------------------------------------------------
public:
//------------------------------------------------------------------------------

    /**
     * constructor
     */
    Hdf5File(
            const std::string & path,
            const enum Hdf5FileMode mode=Hdf5FileMode::NEW );

//------------------------------------------------------------------------------

    /**
     * destructor
     */
    ~Hdf5File();

//------------------------------------------------------------------------------

    /**
     * close the file
     */
     void
     close();

//------------------------------------------------------------------------------
// group functions
//------------------------------------------------------------------------------

    /**
     * create a data group
     *
     * /param label  name of new group
     * /param parent parent object, default: file
     */
    hid_t
    create_group( const std::string & label, const hid_t parent=-1 );

//------------------------------------------------------------------------------

    /**
     * select a data group
     *
     * /param label  name of group to select
     * /param parent parent object, default: file
     */
    hid_t
    select_group( const std::string & label, const hid_t parent=-1 );

//------------------------------------------------------------------------------

    /**
     * close the active group
     */
    void
    close_active_group();

//------------------------------------------------------------------------------
// special type read functions
//------------------------------------------------------------------------------

    /**
     * read a boolean from the active group
     *
     * \param label   name of boolean
     * \param value   value of boolean
     */
     void
     read( const std::string & label, bool & value );


//------------------------------------------------------------------------------

    /**
     * read a string from the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    read( const std::string & label, std::string & value );

//------------------------------------------------------------------------------
// scalar read functions
//------------------------------------------------------------------------------

    /**
     * read an int from the active group
     *
     * \param label   name of boolean
     * \param value   content of int
     */
    void
    read( const std::string & label, int & value );


//------------------------------------------------------------------------------

    /**
     * read a long int from the active group
     *
     * \param label   name of boolean
     * \param value   content of int
     */
    void
    read( const std::string & label, long int &  value );

//------------------------------------------------------------------------------

    /**
     * read an unsigned int from the active group
     *
     * \param label   name of boolean
     * \param value   content of int
     */
    void
    read( const std::string & label, unsigned int &  value );

//------------------------------------------------------------------------------

    /**
     * read a long  unsigned int from the active group
     *
     * \param label   name of boolean
     * \param value   content of int
     */
    void
    read( const std::string & label, long unsigned int &  value );

//------------------------------------------------------------------------------

    /**
     * read a float from the active group
     *
     * \param label   name of boolean
     * \param value   content of int
     */
    void
    read( const std::string & label, float &  value );

//------------------------------------------------------------------------------

    /**
     * read a double from the active group
     *
     * \param label   name of boolean
     * \param value   content of int
     */
    void
    read( const std::string & label, double &  value );

//------------------------------------------------------------------------------

    /**
     * read a long double from the active group
     *
     * \param label   name of boolean
     * \param value   content of int
     */
    void
    read( const std::string & label, long double &  value );

//------------------------------------------------------------------------------
// array read functions
//------------------------------------------------------------------------------

    /**
     * read an int array from the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array with data
     * \param size    size of array
     */
    void
    read( const std::string & label, int * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * read a long int array from the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array with data
     * \param size    size of array
     */
    void
    read( const std::string & label, long int * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * read an unsigned int array from the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array with data
     * \param size    size of array
     */
    void
    read( const std::string & label, unsigned int * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * read a long unsigned int array from the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array with data
     * \param size    size of array
     */
    void
    read( const std::string & label, long unsigned int * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * read a float array from the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array with data
     * \param size    size of array
     */
    void
    read( const std::string & label, float * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * read a double array from the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array with data
     * \param size    size of array
     */
    void
    read( const std::string & label, double * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * read a double array from the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array with data
     * \param size    size of array
     */
    void
    read( const std::string & label, long double * data, const hsize_t size );

//------------------------------------------------------------------------------
// special type write functions
//------------------------------------------------------------------------------

    /**
     * write a boolean into the active group
     *
     * \param label  name of boolean
     * \param value  value of boolean
     */
    void
    write( const std::string & label, const bool value );

//------------------------------------------------------------------------------

    /**
     * write a string into the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    write( const std::string & label, const std::string & value );

//------------------------------------------------------------------------------
// scalar write functions
//------------------------------------------------------------------------------

    /**
     * write an int into the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    write( const std::string & label, const int value );

//------------------------------------------------------------------------------

    /**
     * write a long int into the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    write( const std::string & label, const long int value );

//------------------------------------------------------------------------------

    /**
     * write an unsigned int into the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    write( const std::string & label, const unsigned int value );

//------------------------------------------------------------------------------

    /**
     * write a long unsigned int into the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    write( const std::string & label, const long unsigned int value );

//------------------------------------------------------------------------------

    /**
     * write a float into the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    write( const std::string & label, const float value );

//------------------------------------------------------------------------------

    /**
     * write a double into the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    write( const std::string & label, const double value );

//------------------------------------------------------------------------------

    /**
     * write a long double into the active group
     *
     * \param label   name of boolean
     * \param value   content of string
     */
    void
    write( const std::string & label, const long double value );

//------------------------------------------------------------------------------
// array write functions
//------------------------------------------------------------------------------

    /**
     * write an int array into the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array for data
     * \param size    size of array
     */
    void
    write( const std::string & label, const int * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * write a long int array into the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array for data
     * \param size    size of array
     */
    void
    write( const std::string & label, const long int * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * write an unsigned int array into the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array for data
     * \param size    size of array
     */
    void
    write( const std::string & label, const unsigned int * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * write an long unsigned int array into the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array for data
     * \param size    size of array
     */
    void
    write( const std::string & label, const long unsigned int * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * write a float array into the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array for data
     * \param size    size of array
     */
    void
    write( const std::string & label, const float * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * write a float array into the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array for data
     * \param size    size of array
     */
    void
    write( const std::string & label, const double * data, const hsize_t size );

//------------------------------------------------------------------------------

    /**
     * write a float array into the active group
     *
     * \param label   name of boolean
     * \param data    preallocated array for data
     * \param size    size of array
     */
    void
    write( const std::string & label, const long double * data, const hsize_t size );

//------------------------------------------------------------------------------
private:
//------------------------------------------------------------------------------

    /**
     * check if a file exists
     */
     bool
     file_exists( const std::string & path );

//------------------------------------------------------------------------------
};

#endif //MATREAD_HD5FILE_HPP

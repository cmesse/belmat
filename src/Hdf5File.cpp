//
// Created by Christian Messe on 12/6/21.
//

#include <fstream>
#include "Hdf5File.h"
#include "errortools.h"
#include "hdf5tools_datasets.h"
#include "hdf5tools_read.h"
#include "hdf5tools_write.h"

//------------------------------------------------------------------------------

Hdf5File::Hdf5File(
        const std::string & path,
        const enum Hdf5FileMode mode ) :
                mPath( path )
{
#ifdef HDF5

    // first we check if a path was given
    errortools::error( ! mPath.empty(), "There was no file path given." );

    // check the file mode
    switch( mode )
    {
        case( Hdf5FileMode::NEW ) :
        {
            // call HDF5 API
            mFile = H5Fcreate(
                    mPath.c_str(),
                    H5F_ACC_TRUNC, // If file exists, erasing all existing data.
                    H5P_DEFAULT,   // File creation property list identifier
                    H5P_DEFAULT);  // Access property list identifier.

            break;
        }
        case( Hdf5FileMode::OPEN_RDONLY ) :
        {
            // make sure that file exists
            errortools::error( this->file_exists( mPath ),
                        "File %s does not exist.",
                        mPath.c_str() );

            // call HDF5 API
            mFile = H5Fopen(
                    mPath.c_str(),
                    H5F_ACC_RDONLY,   // File creation property list identifier
                    H5P_DEFAULT);  // Access property list identifier.

            break ;
        }
        case( Hdf5FileMode::OPEN_RDWR ) :
        {
            // make sure that file exists
            errortools::error( this->file_exists( mPath ),
                   "File %s does not exist.",
                   mPath.c_str() );

            // call HDF5 API
            mFile = H5Fopen(
                    mPath.c_str(),
                    H5F_ACC_RDWR,   // File creation property list identifier
                    H5P_DEFAULT);  // Access property list identifier.
                    break ;
        }
        default:
        {
            errortools::error( false, "Illegal file mode");
        }
    }

    // check status of file
    errortools::error(mFile > 0,
          "Something went wrong while trying to create file %s\nIs it in use?",
          mPath.c_str());

    mFileIsOpen = true ;

    // copy file into active group
    mActiveGroup = mFile ;
#endif
}

//------------------------------------------------------------------------------

Hdf5File::~Hdf5File()
{
    // close the file, if it is still open
    this->close();
}

//------------------------------------------------------------------------------

bool
Hdf5File::file_exists( const std::string & path )
{
    // create a stream object to check if the file exist
    std::ifstream file( path );

    // check status
    if( file )
    {
        // close file
        file.close() ;

        // return that it exists
        return true ;
    }
    else
    {
        // file does not exist
        return false ;
    }
}

//------------------------------------------------------------------------------

void
Hdf5File::close()
{
#ifdef HDF5
    // check if file is open
    if( mFileIsOpen ) {
        // check if an a group is open
        if ( mActiveGroup != mFile )
        {
            // close this group
            H5Gclose(mActiveGroup );
        }

        // close the file
        H5Fclose(mFile);

        // unset flag
        mFileIsOpen = false;
    }
#endif
}

//------------------------------------------------------------------------------

hid_t
Hdf5File::create_group( const std::string & label, const hid_t parent )
{
#ifdef HDF5
    // location of new group
    hid_t target = ( parent == -1 ) ? ( mFile ) : ( parent ) ;

    // check if group exists
    errortools::error( ! hdf5tools::group_exists( target, label ),
                      "group %s already exists in file %s",
                      label.c_str(),
                      mPath.c_str() );

    // add backslash to label
    std::string l = "/" + label;

    // check if any group is open
    if( mActiveGroup != mFile )
    {
        // close active group
        H5Gclose( mActiveGroup );
    }

    // create group and close it
    mActiveGroup =
            H5Gcreate2(
                    target,
                    label.c_str(),
                    H5P_DEFAULT,
                    H5P_DEFAULT,
                    H5P_DEFAULT );

    // return the group
    return mActiveGroup ;
#else
    return -1;
#endif
}

//------------------------------------------------------------------------------

hid_t
Hdf5File::select_group( const std::string & label, const hid_t parent )
{
#ifdef HDF5
    // location of new group
    hid_t target = ( parent == -1 ) ? ( mFile ) : ( parent ) ;

    // check if group exists
    errortools::error( hdf5tools::group_exists( target, label ),
                      "group %s does not exists in file %s",
                      label.c_str(),
                      mPath.c_str() );

    // check if any group is open
    if( mActiveGroup != mFile )
    {
        // close active group
        H5Gclose( mActiveGroup );
    }

    // load the grou-
    mActiveGroup = H5Gopen2(
            target,
            label.c_str(),
            H5P_DEFAULT );

    return mActiveGroup;
#else
    return -1;
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::close_active_group() {
#ifdef HDF5
    if( mActiveGroup != mFile )
    {
        // close active group
        H5Gclose( mActiveGroup );

        mActiveGroup = mFile;
    }
#endif
}

//------------------------------------------------------------------------------
// special type read functions
//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, bool & value ) {
#ifdef HDF5
   hdf5tools::read_bool( mActiveGroup, label, value, mStatus );
#endif
}


//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, std::string & value ) {
#ifdef HDF5
    hdf5tools::read_string( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------
// scalar read functions
//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, int & value ) {
#ifdef HDF5
    hdf5tools::read_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, long int & value ) {
#ifdef HDF5
    hdf5tools::read_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, unsigned int & value ) {
#ifdef HDF5
    hdf5tools::read_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, long unsigned int & value ) {
#ifdef HDF5
    hdf5tools::read_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, float & value ) {
#ifdef HDF5
    hdf5tools::read_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, double & value ) {
#ifdef HDF5
    hdf5tools::read_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, long double & value ) {
#ifdef HDF5
    hdf5tools::read_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------
// array read functions
//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, int * data, const hsize_t size ) {
#ifdef HDF5
    hdf5tools::read_array( mActiveGroup, label, data, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, long int * data, const hsize_t size ) {
#ifdef HDF5
    hdf5tools::read_array( mActiveGroup, label, data, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, unsigned int * data, const hsize_t size ) {
#ifdef HDF5
    hdf5tools::read_array( mActiveGroup, label, data, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, long unsigned int * data, const hsize_t size ) {
#ifdef HDF5
    hdf5tools::read_array( mActiveGroup, label, data, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, float * data, const hsize_t size ) {
#ifdef HDF5
    hdf5tools::read_array( mActiveGroup, label, data, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, double * data, const hsize_t size ) {
#ifdef HDF5
    hdf5tools::read_array( mActiveGroup, label, data, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::read( const std::string & label, long double  * data, const hsize_t size ) {
#ifdef HDF5
    hdf5tools::read_array( mActiveGroup, label, data, size, mStatus );
#endif
}

//------------------------------------------------------------------------------
// special type write functions
//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const bool value ) {
#ifdef HDF5
    hdf5tools::write_bool( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const std::string & value ) {
#ifdef HDF5
    hdf5tools::write_string( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------
// scalar write functions
//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const int value ) {
#ifdef HDF5
    hdf5tools::write_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const long int value ) {
#ifdef HDF5
    hdf5tools::write_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const unsigned int value ) {
#ifdef HDF5
    hdf5tools::write_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const long unsigned int value ) {
#ifdef HDF5
    hdf5tools::write_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const float value ) {
#ifdef HDF5
    hdf5tools::write_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const double value ) {
#ifdef HDF5
    hdf5tools::write_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const long double value ) {
#ifdef HDF5
    hdf5tools::write_scalar( mActiveGroup, label, value, mStatus );
#endif
}

//------------------------------------------------------------------------------
// array write functions
//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const int * value, const hsize_t size ){
#ifdef HDF5
    hdf5tools::write_array( mActiveGroup, label, value, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const long int * value, const hsize_t size ){
#ifdef HDF5
    hdf5tools::write_array( mActiveGroup, label, value, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const unsigned int * value, const hsize_t size ){
#ifdef HDF5
    hdf5tools::write_array( mActiveGroup, label, value, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const long unsigned int * value, const hsize_t size ){
#ifdef HDF5
    hdf5tools::write_array( mActiveGroup, label, value, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const float * value, const hsize_t size ){
#ifdef HDF5
    hdf5tools::write_array( mActiveGroup, label, value, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const double * value, const hsize_t size ){
#ifdef HDF5
    hdf5tools::write_array( mActiveGroup, label, value, size, mStatus );
#endif
}

//------------------------------------------------------------------------------

void
Hdf5File::write( const std::string & label, const long double * value, const hsize_t size ){
#ifdef HDF5
    hdf5tools::write_array( mActiveGroup, label, value, size, mStatus );
#endif
}

//
// Created by Christian Messe on 31.07.22.
//
#include <limits>
#include "material.h"
#include "Hdf5File.h"

namespace belmat
{
//------------------------------------------------------------------------------

    material::material( const std::string & filename, const std::string & label, const uint rrr) :
        mLabel( label ),
        mRRR( rrr ),
        mCpLabel( label + "Cp" ),
        mKLabel( rrr == 0 ? label+"K" : label+"K_RRR"+ std::to_string( rrr ) ),
        mRhoLabel( rrr == 0 ? label+"Rho" : label+"Rho_RRR"+ std::to_string( rrr ) )
    {

        // open the file
        Hdf5File hdf5file( filename, Hdf5FileMode::OPEN_RDONLY );

        // read density
        if( hdf5file.dataset_exists( label + "Density" ) )
        {
            hdf5file.read( label + "Density", mDensity );
        }
        else
        {
            mDensity = std::numeric_limits<double>::quiet_NaN();
        }

        // test if datasets exists
        if( hdf5file.group_exists( mCpLabel ) )
        {
            mCp = new database( hdf5file, mCpLabel );
        }
        if( hdf5file.group_exists( mKLabel ) )
        {
            mK = new database( hdf5file, mKLabel );
        }
        if( hdf5file.group_exists( mRhoLabel ) )
        {
            mRho = new database( hdf5file, mRhoLabel );
        }

        hdf5file.close();
    }

//------------------------------------------------------------------------------

    material::~material()
    {
        if( this->has_cp() )
        {
            delete mCp ;
        }
        if( this->has_k() )
        {
            delete mK ;
        }
        if ( this->has_rho() )
        {
            delete mRho ;
        }
    }

//------------------------------------------------------------------------------
}
//
// Created by christian on 4/28/22.
//
#include <cmath>
#include <numeric>
#include "errortools.h"

#include "mpitools_comm.h"
#include "Hdf5File.h"
#include "database.h"

namespace belmat
{
//----------------------------------------------------------------------------

    database::database( const std::string & filename,
                        const std::string & tablename )
    {
#ifdef HDF5
        if( mpitools::comm_rank() == 0 )
        {
            this->load_from_file(filename, tablename);
        }

        // send data from master to others
        this->distribute_data() ;

        // link the shape function
        this->link_interpolation_function();
#endif
    }

//----------------------------------------------------------------------------

    // the database constructor (passing database)
    database::database( Hdf5File & file ,
        const std::string & tablename )
    {
#ifdef HDF5
        if( mpitools::comm_rank() == 0 )
        {
            this->load_from_file(file, tablename);
        }

        // send data from master to others
        this->distribute_data() ;

        // link the shape function
        this->link_interpolation_function();
#endif
    }

//----------------------------------------------------------------------------

    database::~database()
    {
        this->free();
    }

//----------------------------------------------------------------------------

    void
    database::free()
    {
        if( mNumPoints != nullptr )
        {
            ::free( mNumPoints );
            mNumPoints = nullptr ;
        }
        if( mXmax != nullptr )
        {
            ::free( mXmax );
            mXmax = nullptr ;
        }
        if( mXmin != nullptr )
        {
            ::free( mXmin );
            mXmin = nullptr ;
        }
        if( mX != nullptr )
        {
            ::free( mX );
            mX = nullptr ;
        }
        if ( mStep != nullptr )
        {
            ::free( mStep );
            mStep = nullptr ;
        }
        if ( mInvStep != nullptr )
        {
            ::free( mInvStep );
            mInvStep = nullptr ;
        }
        if ( mValues != nullptr )
        {
            ::free( mValues );
            mValues = nullptr ;
        }
        if( mWork != nullptr )
        {
            ::free( mWork );
            mWork = nullptr ;
        }
        if( mCoeffs != nullptr )
        {
            ::free( mCoeffs );
            mCoeffs = nullptr ;
        }

        mInterpolationOrder = 0 ;
        mNumberOfDimensions = 0 ;
        mMemorySize = 0 ;
        mNumberOfPointsPerCell = 0 ;

    }

//----------------------------------------------------------------------------

    void
    database::distribute_data()
    {
        // wait for other procs
        mpitools::comm_barrier() ;

        // prepare info tag
        std::vector< uint > iinfo( 8, 0 );
        std::vector< double > dinfo ;

        if( mpitools::comm_rank() == 0 )
        {
            std::size_t count = 0 ;

            iinfo[ count++ ] = mInterpolationOrder ;
            iinfo[ count++ ] = mNumberOfDimensions ;
            iinfo[ count++ ] = mMemorySize ;

            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                iinfo[ count++ ] = mNumPoints[ k ];
            }

            mpitools::distribute( iinfo );

            dinfo.resize( 3 * mNumberOfDimensions );

            // reset counter
            count = 0 ;

            // send min and max data
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                dinfo[ count++ ] = mXmin[ k ];
            }
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                dinfo[ count++ ] = mXmax[ k ];
            }
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                dinfo[ count++ ] = mStep[ k ];
            }
            mpitools::distribute( dinfo );

            // catch special case
            if( mOrder == 0 )
            {
                dinfo.resize( 4 * mMemorySize );
                std::copy(mCoeffs, mCoeffs + 4 * mMemorySize, dinfo.begin() );
            }
            else
            {
                // prepare data
                dinfo.resize(mMemorySize);
                std::copy(mValues, mValues + mMemorySize, dinfo.begin());
            }
            // send data
            mpitools::distribute( dinfo );

        }
        else
        {
            std::size_t count = 0 ;

            mpitools::distribute( iinfo );

            // unpack info
            mInterpolationOrder    = iinfo[ count++ ];
            mNumberOfDimensions    = iinfo[ count++ ];
            mMemorySize            = iinfo[ count++ ] ;

            mNumPoints = ( uint * ) malloc( mNumberOfDimensions * sizeof ( uint  ) ) ;
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mNumPoints[ k ] = iinfo[ count++ ];
            }

            // order, but as double
            mOrder = double( mInterpolationOrder );

            // compute the inverse of the order
            mInvOrder = 1.0 / mOrder ;

            // receive double values
            mpitools::distribute( dinfo );

            // allocate data
            mXmin = ( double * ) malloc( mNumberOfDimensions * sizeof ( double  ) ) ;
            mXmax = ( double * ) malloc( mNumberOfDimensions * sizeof ( double  ) ) ;
            mStep = ( double * ) malloc( mNumberOfDimensions * sizeof ( double  ) ) ;
            mInvStep = ( double * ) malloc( mNumberOfDimensions * sizeof ( double  ) ) ;

            // allocate work array
            mX = ( double * ) malloc( mNumberOfDimensions * sizeof ( double ) ) ;

            // populate data
            count = 0 ;
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mXmin[ k ] = dinfo[ count++ ];
            }
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mXmax[ k ] = dinfo[ count++ ];
            }
            for( uint k=0; k<mNumberOfDimensions; ++k )
            {
                mStep[ k ] = dinfo[ count++ ];
                mInvStep[ k ] = 1.0 / mStep[ k ];
            }

            mpitools::distribute( dinfo );

            if( mOrder == 0 )
            {
                // allocate memory for coefficients
                mCoeffs = ( double * ) malloc(4 * mMemorySize * sizeof( double ) );

                // copy data from vector
                std::copy(dinfo.begin(), dinfo.end(), mCoeffs );
            }
            else
            {
                // allocate value memory
                mValues = ( double * ) malloc(mMemorySize * sizeof( double ) );

                // copy data from vector
                std::copy(dinfo.begin(), dinfo.end(), mValues );
            }
        }

        // wait for other procs
        mpitools::comm_barrier() ;
    }

//----------------------------------------------------------------------------

    void
    database::link_interpolation_function()
    {
        switch( mNumberOfDimensions )
        {
            case( 1 ) :
            {
                switch ( mInterpolationOrder )
                {
                    case ( 0 ) :
                    {
                        mInterpolate = & database::interpolate_1dspline ;
                        mNumberOfPointsPerCell = 0 ;
                        break ;
                    }
                    case ( 1 ) :
                    {
                        mInterpolate = &database::interpolate_line2;
                        mNumberOfPointsPerCell = 2;
                        break;
                    }
                    case ( 2 ) :
                    {
                        mInterpolate = &database::interpolate_line3;
                        mNumberOfPointsPerCell = 3;
                        break;
                    }
                    case ( 3 ) :
                    {
                        mInterpolate = &database::interpolate_line4;
                        mNumberOfPointsPerCell = 4;
                        break;
                    }
                    default :
                    {
                        errortools::error(false, "Invalid interpolation order: %u",
                                          (unsigned int) mInterpolationOrder);
                    }
                }
                break ;
            }
            case( 2 ) :
            {
                switch ( mInterpolationOrder )
                {
                    case( 1 ) :
                    {
                        mInterpolate = & database::interpolate_quad4 ;
                        mNumberOfPointsPerCell = 4 ;
                        break ;
                    }
                    case( 2 ) :
                    {
                        mInterpolate = & database::interpolate_quad9 ;
                        mNumberOfPointsPerCell = 9 ;
                        break ;
                    }
                    case( 3 ) :
                    {
                        mInterpolate = & database::interpolate_quad16 ;
                        mNumberOfPointsPerCell = 16 ;
                        break ;
                    }
                    default :
                    {
                        errortools::error( false, "Invalid interpolation order: %u",
                                           ( unsigned int ) mInterpolationOrder );
                    }
                }
                break ;
            }
            case( 3 ) :
            {
                switch ( mInterpolationOrder )
                {
                    case( 1 ) :
                    {
                        mInterpolate = &  database::interpolate_hex8 ;
                        mNumberOfPointsPerCell = 8 ;
                        break ;
                    }
                    case( 2 ) :
                    {
                        mInterpolate = & database::interpolate_hex27 ;
                        mNumberOfPointsPerCell = 27 ;
                        break ;
                    }
                    case( 3 ) :
                    {
                        mInterpolate = & database::interpolate_hex64 ;
                        mNumberOfPointsPerCell = 64 ;
                        break ;
                    }
                    default :
                    {
                        errortools::error( false, "Invalid interpolation order: %u",
                                           ( unsigned int ) mInterpolationOrder );
                    }
                }
                break ;
            }
            default :
            {
                errortools::error( false, "Invalid dimension: %u",
                                   ( unsigned int ) mNumberOfDimensions );
            }
        }

        if( mNumberOfPointsPerCell > 0 )
        {
            // allocate work vector
            mWork = (double *)
                    malloc(mNumberOfPointsPerCell * sizeof(double));
        }
    }

//----------------------------------------------------------------------------

    void
    database::load_from_file(
            const std::string & filename,
            const std::string & tablename )
    {
        // open the database
        Hdf5File file( filename, Hdf5FileMode::OPEN_RDONLY );

        this->load_from_file( file, tablename );
        file.close();
    }

//----------------------------------------------------------------------------

    void
    database::load_from_file(
            Hdf5File & file,
            const std::string & tablename )
    {
        // reset data
        this->free();

        // select the group
        file.select_group( tablename );

        // get the number of dimensions
        file.read( "dimension", mNumberOfDimensions );

        // get the interpolation order
        file.read( "order", mInterpolationOrder );

        // order, but as double
        mOrder = double( mInterpolationOrder );

        // compute the inverse of the order
        mInvOrder = 1.0 / mOrder ;

        // read the number of samples
        mNumPoints = ( uint * ) malloc( mNumberOfDimensions * sizeof ( uint  ) ) ;
        file.read( "numpoints", mNumPoints, mNumberOfDimensions );

        // read the stepsize
        mStep =  ( double * ) malloc( mNumberOfDimensions * sizeof ( double ) ) ;
        file.read( "step", mStep, mNumberOfDimensions );

        // inverse stepsize
        mInvStep = ( double * ) malloc( mNumberOfDimensions * sizeof ( double ) ) ;
        for( uint k = 0; k<mNumberOfDimensions; ++k )
        {
            mInvStep[ k ] = 1.0 / mStep[ k ];
        }

        // read the offset
        mXmin = ( double * ) malloc( mNumberOfDimensions * sizeof ( double ) ) ;
        file.read( "offset", mXmin, mNumberOfDimensions );

        // compute the maximum values
        mXmax = ( double * ) malloc( mNumberOfDimensions * sizeof ( double ) ) ;
        for( uint d=0; d<mNumberOfDimensions; ++d )
        {
            mXmax[ d ] = mXmin[ d ] + mStep[ d ] *  double( mNumPoints[ d ] - 1 ) * mInvOrder ;
        }

        // allocate the work vector for the coordinates
        mX = ( double * ) malloc( mNumberOfDimensions * sizeof ( double ) ) ;

        for( uint d=0; d<mNumberOfDimensions; ++d )
        {
            mMemorySize *= mNumPoints[ d ];
        }

        // compute the memory size
        mMemorySize = 1 ;
        for( uint d=0; d<mNumberOfDimensions; ++d )
        {
            mMemorySize *= mNumPoints[ d ];
        }

        if( mOrder == 0 )
        {
            // read the precomputed coefficients
            std::size_t tSize = 4 * mMemorySize ;
            mCoeffs = ( double * ) malloc( tSize * sizeof ( double  ) ) ;
            file.read( "coeffs", mCoeffs, tSize );
        }
        else
        {
            // read the data
            mValues = ( double * ) malloc(mMemorySize * sizeof( double ) );
            file.read("values", mValues, mMemorySize );
        }
        file.close_active_group();
    }

//----------------------------------------------------------------------------

    double
    database::compute( const double T )
    {
        mX[ 0 ] = T ;
        return ( this->*mInterpolate )( mX );
    }

//----------------------------------------------------------------------------

    double
    database::compute( const double T, const double B )
    {
        mX[ 0 ] = T ;
        mX[ 1 ] = B ;
        return ( this->*mInterpolate )( mX );
    }

//----------------------------------------------------------------------------

    double
    database::compute( const double T, const double B, const double theta )
    {
#ifdef HDF5
        mX[ 0 ] = theta ;
        mX[ 1 ] = T ;
        mX[ 2 ] = B ;
        return ( this->*mInterpolate )( mX );
#else
        double Jc0 = 3.5e10/46.5 ;
        double alpha=1.2 ;
        double beta=1 ;
        double gamma=8.02 ;
        double Bc0=240e-3;

        double eps= std::sqrt(std::pow(std::sin(theta),2)/gamma
                    +std::pow(std::cos(theta),2) );

        return Jc0/std::pow(1+eps*std::pow(B/Bc0,beta),alpha);
#endif
    }

//----------------------------------------------------------------------------

    double
    database::interpolate_line2( const double * x )
    {
        // identify the indices of the first point of the cell
        const uint i = find_cell_linear( x, 0 );

        // compute parameter coordinates
        const double xi  = this->xtoxi_linear( x, i, 0 );

        // compute the interpolation
        return 0.5 * (  mValues[ i ] * ( 1.0 - xi )
                      + mValues[ i + 1 ] * ( 1.0 + xi ) );
    }

//----------------------------------------------------------------------------

    double
    database::interpolate_line3( const double * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = find_cell_higher_order( x, 0 );

        // step 2: collect the node data
        mWork[ 0 ] = mValues[ i ];
        mWork[ 1 ] = mValues[ i+2 ];
        mWork[ 2 ] = mValues[ i+1 ];

        // step 3: compute parameter coordinate
        const double xi  = this->xtoxi_higher_order( x, i, 0 );
        double xi2 = xi*xi ;

        // step 4: compute the interpolation
        mWork[ 0 ] *=  xi2 - xi ;
        mWork[ 1 ] *= xi2 + xi  ;
        mWork[ 2 ] *= 1.0 - xi2;

        return 0.5 * ( mWork[ 0 ] +  mWork[ 1 ] ) +  mWork[ 2 ] ;
    }

//----------------------------------------------------------------------------

    double
    database::interpolate_line4( const double * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = find_cell_higher_order( x, 0 );

        // step : compute parameter coordinate
        const double xi  = this->xtoxi_higher_order( x, i, 0 );

        // compute helpers
        mWork[ 0 ] = 1.0 - xi ;
        mWork[ 1 ] = 1.0 + xi ;
        mWork[ 2 ] = 3.0 * xi - 1.0 ;
        mWork[ 3 ] = 3.0 * xi + 1.0 ;

        // compute interpolation
        return  (   mWork[ 2 ] * mWork[ 3 ] *
                  ( mWork[ 0 ] * mValues[ i ] + mWork[ 1 ] * mValues[ i+3 ] )
                + 9.0 * mWork[ 0 ] * mWork[ 1 ] * ( mWork[ 3 ] * mValues[ i+2 ]
                 - mWork[ 2 ] * mValues[ i+1 ] ) ) *  0.06250 ;
    }

//----------------------------------------------------------------------------

    double
    database::interpolate_quad4( const double * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = find_cell_linear( x, 0 );
        const uint j = find_cell_linear( x, 1 );

        // step 2: collect the node data
        mWork[ 0 ] = mValues[ this->index( i, j ) ];
        mWork[ 1 ] = mValues[ this->index( i+1, j ) ];
        mWork[ 2 ] = mValues[ this->index( i+1, j+1 ) ];
        mWork[ 3 ] = mValues[ this->index( i, j+1 ) ];

        // step 3: compute parameter coordinates
        const double xi  = this->xtoxi_linear( x, i, 0 );
        const double eta = this->xtoxi_linear( x, j, 1 );

        // step 4: compute the interpolation
        mWork[ 0 ] *= ( ( 1.0 - xi ) * ( 1.0 - eta ) ) ;
        mWork[ 1 ] *= ( ( 1.0 + xi ) * ( 1.0 - eta ) ) ;
        mWork[ 2 ] *= ( ( 1.0 + xi ) * ( 1.0 + eta ) ) ;
        mWork[ 3 ] *= ( ( 1.0 - xi ) * ( 1.0 + eta ) ) ;

        // compute and return the value
        return 0.25 * std::accumulate( mWork, mWork + 4, 0.0 ) ;
    }

//----------------------------------------------------------------------------

    double
    database::interpolate_quad9( const double * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = find_cell_higher_order( x, 0 );
        const uint j = find_cell_higher_order( x, 1 );

        // step 2: collect the node data
		mWork[ 0 ] = mValues[ this->index( i  , j   ) ];
		mWork[ 1 ] = mValues[ this->index( i+2, j   ) ];
		mWork[ 2 ] = mValues[ this->index( i+2, j+2 ) ];
		mWork[ 3 ] = mValues[ this->index( i  , j+2 ) ];
		mWork[ 4 ] = mValues[ this->index( i+1, j   ) ];
		mWork[ 5 ] = mValues[ this->index( i+2, j+1 ) ];
		mWork[ 6 ] = mValues[ this->index( i+1, j+2 ) ];
		mWork[ 7 ] = mValues[ this->index( i  , j+1 ) ];
		mWork[ 8 ] = mValues[ this->index( i+1, j+1 ) ];

        // step 3: compute parameter coordinates
        const double xi  = this->xtoxi_quadratic( x, i, 0 );
        const double eta = this->xtoxi_quadratic( x, j, 1 );

        // step 4: compute the interpolation
        const double    c = xi * eta * 0.25;
        const double  xi2 = xi*xi;
        const double eta2 = eta*eta;
        mWork[ 0 ] *= ( c * ( eta - 1.0 ) * (xi - 1.0) );
        mWork[ 1 ] *= ( c * ( eta - 1.0 ) * (xi + 1.0) );
        mWork[ 2 ] *= ( c * ( eta + 1.0 ) * (xi + 1.0) );
        mWork[ 3 ] *= ( c * ( eta + 1.0 ) * (xi - 1.0) );
        mWork[ 4 ] *= ( eta * ( 1.0 - xi2 ) * ( eta - 1.0 ) ) * 0.5;
        mWork[ 5 ] *= ( xi * ( 1.0 - eta2)*( xi + 1.0 ) )*0.5;
        mWork[ 6 ] *= ( eta * (1.0 - xi2)*( eta + 1.0 ) )*0.5;
        mWork[ 7 ] *= ( xi*( 1.0 - eta2 )*( xi - 1.0 ) )*0.5;
        mWork[ 8 ] *= ( eta2 - 1.0 )*( xi2 - 1.0 );

        return std::accumulate( mWork, mWork + 9, 0.0 );
    }

//----------------------------------------------------------------------------

    double
    database::interpolate_quad16( const double * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = find_cell_higher_order( x, 0 );
        const uint j = find_cell_higher_order( x, 1 );

        // step 2: collect the node data
        mWork[  0 ] = mValues[ this->index( i  , j   ) ];
		mWork[  1 ] = mValues[ this->index( i+3, j   ) ];
		mWork[  2 ] = mValues[ this->index( i+3, j+3 ) ];
		mWork[  3 ] = mValues[ this->index( i  , j+3 ) ];
		mWork[  4 ] = mValues[ this->index( i+1, j   ) ];
		mWork[  5 ] = mValues[ this->index( i+2, j   ) ];
		mWork[  6 ] = mValues[ this->index( i+3, j+1 ) ];
		mWork[  7 ] = mValues[ this->index( i+3, j+2 ) ];
		mWork[  8 ] = mValues[ this->index( i+2, j+3 ) ];
		mWork[  9 ] = mValues[ this->index( i+1, j+3 ) ];
		mWork[ 10 ] = mValues[ this->index( i  , j+2 ) ];
		mWork[ 11 ] = mValues[ this->index( i  , j+1 ) ];
		mWork[ 12 ] = mValues[ this->index( i+1, j+1 ) ];
		mWork[ 13 ] = mValues[ this->index( i+2, j+1 ) ];
		mWork[ 14 ] = mValues[ this->index( i+2, j+2 ) ];
		mWork[ 15 ] = mValues[ this->index( i+1, j+2 ) ];

        // step 3: compute parameter coordinates
        const double xi  = this->xtoxi_higher_order( x, i, 0 );
        const double eta = this->xtoxi_higher_order( x, j, 1 );

        // step 4: compute the interpolation
        const double a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
        const double a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
        const double a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
        const double a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

        const double b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
        const double b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
        const double b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
        const double b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

        mWork[ 0 ] *= a0*b0;
        mWork[ 1 ] *= a3*b0;
        mWork[ 2 ] *= a3*b3;
        mWork[ 3 ] *= a0*b3;
        mWork[ 4 ] *= a1*b0;
        mWork[ 5 ] *= a2*b0;
        mWork[ 6 ] *= a3*b1;
        mWork[ 7 ] *= a3*b2;
        mWork[ 8 ] *= a2*b3;
        mWork[ 9 ] *= a1*b3;
        mWork[ 10 ] *= a0*b2;
        mWork[ 11 ] *= a0*b1;
        mWork[ 12 ] *= a1*b1;
        mWork[ 13 ] *= a2*b1;
        mWork[ 14 ] *= a2*b2;
        mWork[ 15 ] *= a1*b2;

        return std::accumulate( mWork, mWork + 16, 0.0 );
    }


//----------------------------------------------------------------------------

    double
    database::interpolate_hex8( const double * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = find_cell_linear( x, 0 );
        const uint j = find_cell_linear( x, 1 );
        const uint k = find_cell_linear( x, 2 );

        // step 2: collect the node data
        mWork[ 0 ] = mValues[ this->index( i, j, k ) ];
        mWork[ 1 ] = mValues[ this->index( i+1, j, k ) ];
        mWork[ 2 ] = mValues[ this->index( i+1, j+1, k ) ];
        mWork[ 3 ] = mValues[ this->index( i, j+1, k ) ];
        mWork[ 4 ] = mValues[ this->index( i, j, k+1 ) ];
        mWork[ 5 ] = mValues[ this->index( i+1, j, k+1 ) ];
        mWork[ 6 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
        mWork[ 7 ] = mValues[ this->index( i, j+1, k+1 ) ];

        // step 3: compute parameter coordinates
        const double  xi  = this->xtoxi_linear( x, i, 0 );
        const double  eta = this->xtoxi_linear( x, j, 1 );
        const double zeta = this->xtoxi_linear( x, k, 2 );

        // step 4: compute the interpolation
        mWork[ 0 ] *=  - ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[ 1 ] *=    ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[ 2 ] *=  - ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[ 3 ] *=    ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[ 4 ] *=    ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        mWork[ 5 ] *=  - ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[ 6 ] *=    ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[ 7 ] *=  - ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );

        // compute and return the value
        return 0.125 * std::accumulate( mWork, mWork + 8, 0.0 );
    }

//----------------------------------------------------------------------------

    double
    database::interpolate_hex27( const double * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = find_cell_higher_order( x, 0 );
        const uint j = find_cell_higher_order( x, 1 );
        const uint k = find_cell_higher_order( x, 2 );

        // step 2: collect the node data
		mWork[  0 ] = mValues[ this->index( i  , j  , k   ) ];
		mWork[  1 ] = mValues[ this->index( i+2, j  , k   ) ];
		mWork[  2 ] = mValues[ this->index( i+2, j+2, k+2 ) ];
		mWork[  3 ] = mValues[ this->index( i  , j+2, k+2 ) ];
		mWork[  4 ] = mValues[ this->index( i  , j  , k   ) ];
		mWork[  5 ] = mValues[ this->index( i+2, j  , k   ) ];
		mWork[  6 ] = mValues[ this->index( i+2, j+2, k+2 ) ];
		mWork[  7 ] = mValues[ this->index( i  , j+2, k+2 ) ];
		mWork[  8 ] = mValues[ this->index( i+1, j  , k   ) ];
		mWork[  9 ] = mValues[ this->index( i+2, j+1, k+1 ) ];
		mWork[ 10 ] = mValues[ this->index( i+1, j+2, k+2 ) ];
		mWork[ 11 ] = mValues[ this->index( i  , j+1, k+1 ) ];
		mWork[ 12 ] = mValues[ this->index( i  , j  , k   ) ];
		mWork[ 13 ] = mValues[ this->index( i+2, j  , k   ) ];
		mWork[ 14 ] = mValues[ this->index( i+2, j+2, k+2 ) ];
		mWork[ 15 ] = mValues[ this->index( i  , j+2, k+2 ) ];
		mWork[ 16 ] = mValues[ this->index( i+1, j  , k   ) ];
		mWork[ 17 ] = mValues[ this->index( i+2, j+1, k+1 ) ];
		mWork[ 18 ] = mValues[ this->index( i+1, j+2, k+2 ) ];
		mWork[ 19 ] = mValues[ this->index( i  , j+1, k+1 ) ];
		mWork[ 20 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
		mWork[ 21 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
		mWork[ 22 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
		mWork[ 23 ] = mValues[ this->index( i  , j+1, k+1 ) ];
		mWork[ 24 ] = mValues[ this->index( i+2, j+1, k+1 ) ];
		mWork[ 25 ] = mValues[ this->index( i+1, j  , k   ) ];
		mWork[ 26 ] = mValues[ this->index( i+1, j+2, k+2 ) ];

        // step 3: compute parameter coordinates
        const double   xi = this->xtoxi_quadratic( x, i, 0 );
        const double  eta = this->xtoxi_quadratic( x, j, 1 );
        const double zeta = this->xtoxi_quadratic( x, k, 2 );

        // step 4: compute the interpolation

        const double   xi2 = xi*xi;
        const double  eta2 = eta*eta;
        const double zeta2 = zeta*zeta;

        const double a = -0.25 * eta * zeta;
        const double b = -0.25 * xi * zeta;
        const double c = -0.25 * xi * eta;
        const double d = 0.125 * xi * eta * zeta;

        mWork[  0 ] *=  d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[  1 ] *=  d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[  2 ] *=  d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[  3 ] *=  d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[  4 ] *=  d * ( eta - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        mWork[  5 ] *=  d * ( eta - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[  6 ] *=  d * ( eta + 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[  7 ] *=  d * ( eta + 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        mWork[  8 ] *=  a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta - 1.0 );
        mWork[  9 ] *=  b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta - 1.0 );
        mWork[ 10 ] *=  a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta - 1.0 );
        mWork[ 11 ] *=  b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta - 1.0 );
        mWork[ 12 ] *=  c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi - 1.0 );
        mWork[ 13 ] *=  c * ( zeta2 - 1.0 ) * ( eta - 1.0 ) * ( xi + 1.0 );
        mWork[ 14 ] *=  c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi + 1.0 );
        mWork[ 15 ] *=  c * ( zeta2 - 1.0 ) * ( eta + 1.0 ) * ( xi - 1.0 );
        mWork[ 16 ] *=  a * ( xi2 - 1.0 ) * ( eta - 1.0 ) * ( zeta + 1.0 );
        mWork[ 17 ] *=  b * ( eta2 - 1.0 ) * ( xi + 1.0 ) * ( zeta + 1.0 );
        mWork[ 18 ] *=  a * ( xi2 - 1.0 ) * ( eta + 1.0 ) * ( zeta + 1.0 );
        mWork[ 19 ] *=  b * ( eta2 - 1.0 ) * ( xi - 1.0 ) * ( zeta + 1.0 );
        mWork[ 20 ] *=  -( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta2 - 1.0 );
        mWork[ 21 ] *=  ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta - 1.0 ) ) * 0.5;
        mWork[ 22 ] *=  ( zeta * ( eta2 - 1.0 ) * ( xi2 - 1.0 ) * ( zeta + 1.0 ) ) * 0.5;
        mWork[ 23 ] *=  ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi - 1.0 ) ) * 0.5;
        mWork[ 24 ] *=  ( xi * ( eta2 - 1.0 ) * ( zeta2 - 1.0 ) * ( xi + 1.0 ) ) * 0.5;
        mWork[ 25 ] *=  ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta - 1.0 ) ) * 0.5;
        mWork[ 26 ] *=  ( eta * ( xi2 - 1.0 ) * ( zeta2 - 1.0 ) * ( eta + 1.0 ) ) * 0.5;

        return std::accumulate( mWork, mWork + 27, 0.0 );
    }
//----------------------------------------------------------------------------

    double
    database::interpolate_hex64( const double * x )
    {
        // step 1: identify the indices of the first point of the cell
        const uint i = find_cell_higher_order( x, 0 );
        const uint j = find_cell_higher_order( x, 1 );
        const uint k = find_cell_higher_order( x, 2 );

        // step 2: collect the node data
        mWork[  0 ] = mValues[ this->index( i  , j  , k   ) ];
        mWork[  1 ] = mValues[ this->index( i+3, j  , k   ) ];
        mWork[  2 ] = mValues[ this->index( i+3, j+3, k   ) ];
        mWork[  3 ] = mValues[ this->index( i  , j+3, k   ) ];
        mWork[  4 ] = mValues[ this->index( i  , j  , k+3 ) ];
        mWork[  5 ] = mValues[ this->index( i+3, j  , k+3 ) ];
        mWork[  6 ] = mValues[ this->index( i+3, j+3, k+3 ) ];
        mWork[  7 ] = mValues[ this->index( i  , j+3, k+3 ) ];
        mWork[  8 ] = mValues[ this->index( i+1, j  , k   ) ];
        mWork[  9 ] = mValues[ this->index( i+2, j  , k   ) ];
        mWork[ 10 ] = mValues[ this->index( i  , j+1, k   ) ];
        mWork[ 11 ] = mValues[ this->index( i  , j+2, k   ) ];
        mWork[ 12 ] = mValues[ this->index( i  , j  , k+1 ) ];
        mWork[ 13 ] = mValues[ this->index( i  , j  , k+2 ) ];
        mWork[ 14 ] = mValues[ this->index( i+3, j+1, k   ) ];
        mWork[ 15 ] = mValues[ this->index( i+3, j+2, k   ) ];
        mWork[ 16 ] = mValues[ this->index( i+3, j  , k+1 ) ];
        mWork[ 17 ] = mValues[ this->index( i+3, j  , k+2 ) ];
        mWork[ 18 ] = mValues[ this->index( i+2, j+3, k   ) ];
        mWork[ 19 ] = mValues[ this->index( i+1, j+3, k   ) ];
        mWork[ 20 ] = mValues[ this->index( i+3, j+3, k+1 ) ];
        mWork[ 21 ] = mValues[ this->index( i+3, j+3, k+2 ) ];
        mWork[ 22 ] = mValues[ this->index( i  , j+3, k+1 ) ];
        mWork[ 23 ] = mValues[ this->index( i  , j+3, k+2 ) ];
        mWork[ 24 ] = mValues[ this->index( i+1, j  , k+3 ) ];
        mWork[ 25 ] = mValues[ this->index( i+2, j  , k+3 ) ];
        mWork[ 26 ] = mValues[ this->index( i  , j+1, k+3 ) ];
        mWork[ 27 ] = mValues[ this->index( i  , j+2, k+3 ) ];
        mWork[ 28 ] = mValues[ this->index( i+3, j+1, k+3 ) ];
        mWork[ 29 ] = mValues[ this->index( i+3, j+2, k+3 ) ];
        mWork[ 30 ] = mValues[ this->index( i+2, j+3, k+3 ) ];
        mWork[ 31 ] = mValues[ this->index( i+1, j+3, k+3 ) ];
        mWork[ 32 ] = mValues[ this->index( i+1, j+1, k   ) ];
        mWork[ 33 ] = mValues[ this->index( i+1, j+2, k   ) ];
        mWork[ 34 ] = mValues[ this->index( i+2, j+2, k   ) ];
        mWork[ 35 ] = mValues[ this->index( i+2, j+1, k   ) ];
        mWork[ 36 ] = mValues[ this->index( i+1, j  , k+1 ) ];
        mWork[ 37 ] = mValues[ this->index( i+2, j  , k+1 ) ];
        mWork[ 38 ] = mValues[ this->index( i+2, j  , k+2 ) ];
        mWork[ 39 ] = mValues[ this->index( i+1, j  , k+2 ) ];
        mWork[ 40 ] = mValues[ this->index( i  , j+1, k+1 ) ];
        mWork[ 41 ] = mValues[ this->index( i  , j+1, k+2 ) ];
        mWork[ 42 ] = mValues[ this->index( i  , j+2, k+2 ) ];
        mWork[ 43 ] = mValues[ this->index( i  , j+2, k+1 ) ];
        mWork[ 44 ] = mValues[ this->index( i+3, j+1, k+1 ) ];
        mWork[ 45 ] = mValues[ this->index( i+3, j+2, k+1 ) ];
        mWork[ 46 ] = mValues[ this->index( i+3, j+2, k+2 ) ];
        mWork[ 47 ] = mValues[ this->index( i+3, j+1, k+2 ) ];
        mWork[ 48 ] = mValues[ this->index( i+2, j+3, k+1 ) ];
        mWork[ 49 ] = mValues[ this->index( i+1, j+3, k+1 ) ];
        mWork[ 50 ] = mValues[ this->index( i+1, j+3, k+2 ) ];
        mWork[ 51 ] = mValues[ this->index( i+2, j+3, k+2 ) ];
        mWork[ 52 ] = mValues[ this->index( i+1, j+1, k+3 ) ];
        mWork[ 53 ] = mValues[ this->index( i+2, j+1, k+3 ) ];
        mWork[ 54 ] = mValues[ this->index( i+2, j+2, k+3 ) ];
        mWork[ 55 ] = mValues[ this->index( i+1, j+2, k+3 ) ];
        mWork[ 56 ] = mValues[ this->index( i+1, j+1, k+1 ) ];
        mWork[ 57 ] = mValues[ this->index( i+2, j+1, k+1 ) ];
        mWork[ 58 ] = mValues[ this->index( i+2, j+2, k+1 ) ];
        mWork[ 59 ] = mValues[ this->index( i+1, j+2, k+1 ) ];
        mWork[ 60 ] = mValues[ this->index( i+1, j+1, k+2 ) ];
        mWork[ 61 ] = mValues[ this->index( i+2, j+1, k+2 ) ];
        mWork[ 62 ] = mValues[ this->index( i+2, j+2, k+2 ) ];
        mWork[ 63 ] = mValues[ this->index( i+1, j+2, k+2 ) ];

        // step 3: compute parameter coordinates
        const double  xi  = this->xtoxi_quadratic( x, i, 0 );
        const double  eta = this->xtoxi_quadratic( x, j, 1 );
        const double zeta = this->xtoxi_quadratic( x, k, 2 );

        const double a0 =  ( xi*( 1.0 + 9.0 * xi * ( 1.0 - xi ) ) - 1.0 )*0.0625;
        const double a1 =  ( 9.0 - xi * ( 27.0 + xi*( 9.0 - 27.0*xi ) ) )*0.0625;
        const double a2 =  ( 9.0 + xi * ( 27.0 - xi*( 9.0 + 27.0*xi ) ) )*0.0625;
        const double a3 = ( -xi*( 1.0 - 9.0 * xi * ( 1.0 + xi ) ) - 1.0 )*0.0625;

        const double b0 =  ( eta*( 1.0 + 9.0 * eta * ( 1.0 - eta ) ) - 1.0 )*0.0625;
        const double b1 =  ( 9.0 - eta * ( 27.0 + eta*( 9.0 - 27.0*eta ) ) )*0.0625;
        const double b2 =  ( 9.0 + eta * ( 27.0 - eta*( 9.0 + 27.0*eta ) ) )*0.0625;
        const double b3 = ( -eta*( 1.0 - 9.0 * eta * ( 1.0 + eta ) ) - 1.0 )*0.0625;

        const double c0 =  ( zeta*( 1.0 + 9.0 * zeta * ( 1.0 - zeta ) ) - 1.0 )*0.0625;
        const double c1 =  ( 9.0 - zeta * ( 27.0 + zeta*( 9.0 - 27.0*zeta ) ) )*0.0625;
        const double c2 =  ( 9.0 + zeta * ( 27.0 - zeta*( 9.0 + 27.0*zeta ) ) )*0.0625;
        const double c3 = ( -zeta*( 1.0 - 9.0 * zeta * ( 1.0 + zeta ) ) - 1.0 )*0.0625;

        mWork[  0 ] *= a0 * b0 * c0;
        mWork[  1 ] *= a3 * b0 * c0;
        mWork[  2 ] *= a3 * b3 * c0;
        mWork[  3 ] *= a0 * b3 * c0;
        mWork[  4 ] *= a0 * b0 * c3;
        mWork[  5 ] *= a3 * b0 * c3;
        mWork[  6 ] *= a3 * b3 * c3;
        mWork[  7 ] *= a0 * b3 * c3;
        mWork[  8 ] *= a1 * b0 * c0;
        mWork[  9 ] *= a2 * b0 * c0;
        mWork[ 10 ] *= a0 * b1 * c0;
        mWork[ 11 ] *= a0 * b2 * c0;
        mWork[ 12 ] *= a0 * b0 * c1;
        mWork[ 13 ] *= a0 * b0 * c2;
        mWork[ 14 ] *= a3 * b1 * c0;
        mWork[ 15 ] *= a3 * b2 * c0;
        mWork[ 16 ] *= a3 * b0 * c1;
        mWork[ 17 ] *= a3 * b0 * c2;
        mWork[ 18 ] *= a2 * b3 * c0;
        mWork[ 19 ] *= a1 * b3 * c0;
        mWork[ 20 ] *= a3 * b3 * c1;
        mWork[ 21 ] *= a3 * b3 * c2;
        mWork[ 22 ] *= a0 * b3 * c1;
        mWork[ 23 ] *= a0 * b3 * c2;
        mWork[ 24 ] *= a1 * b0 * c3;
        mWork[ 25 ] *= a2 * b0 * c3;
        mWork[ 26 ] *= a0 * b1 * c3;
        mWork[ 27 ] *= a0 * b2 * c3;
        mWork[ 28 ] *= a3 * b1 * c3;
        mWork[ 29 ] *= a3 * b2 * c3;
        mWork[ 30 ] *= a2 * b3 * c3;
        mWork[ 31 ] *= a1 * b3 * c3;
        mWork[ 32 ] *= a1 * b1 * c0;
        mWork[ 33 ] *= a1 * b2 * c0;
        mWork[ 34 ] *= a2 * b2 * c0;
        mWork[ 35 ] *= a2 * b1 * c0;
        mWork[ 36 ] *= a1 * b0 * c1;
        mWork[ 37 ] *= a2 * b0 * c1;
        mWork[ 38 ] *= a2 * b0 * c2;
        mWork[ 39 ] *= a1 * b0 * c2;
        mWork[ 40 ] *= a0 * b1 * c1;
        mWork[ 41 ] *= a0 * b1 * c2;
        mWork[ 42 ] *= a0 * b2 * c2;
        mWork[ 43 ] *= a0 * b2 * c1;
        mWork[ 44 ] *= a3 * b1 * c1;
        mWork[ 45 ] *= a3 * b2 * c1;
        mWork[ 46 ] *= a3 * b2 * c2;
        mWork[ 47 ] *= a3 * b1 * c2;
        mWork[ 48 ] *= a2 * b3 * c1;
        mWork[ 49 ] *= a1 * b3 * c1;
        mWork[ 50 ] *= a1 * b3 * c2;
        mWork[ 51 ] *= a2 * b3 * c2;
        mWork[ 52 ] *= a1 * b1 * c3;
        mWork[ 53 ] *= a2 * b1 * c3;
        mWork[ 54 ] *= a2 * b2 * c3;
        mWork[ 55 ] *= a1 * b2 * c3;
        mWork[ 56 ] *= a1 * b1 * c1;
        mWork[ 57 ] *= a2 * b1 * c1;
        mWork[ 58 ] *= a2 * b2 * c1;
        mWork[ 59 ] *= a1 * b2 * c1;
        mWork[ 60 ] *= a1 * b1 * c2;
        mWork[ 61 ] *= a2 * b1 * c2;
        mWork[ 62 ] *= a2 * b2 * c2;
        mWork[ 63 ] *= a1 * b2 * c2;

        return std::accumulate( mWork, mWork + 64, 0.0 );
    }
//----------------------------------------------------------------------------
}

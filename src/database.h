//
// Created by christian on 4/28/22.
//

#ifndef BELMAT_DATABASE_H
#define BELMAT_DATABASE_H

#include <cmath>
#include <string>
#include "Hdf5File.h"

namespace belmat
{
//------------------------------------------------------------------------

    typedef unsigned int        uint;

//------------------------------------------------------------------------

    class database
    {
        uint mInterpolationOrder = 0 ;
        uint mNumberOfDimensions = 0 ;
        uint mNumberOfPointsPerCell = 0 ;

        uint mMemorySize = 0 ;

        uint * mNumPoints = nullptr ;

        // inverse of order
        double mOrder = 0 ;
        double mInvOrder = 0 ;

        // miniumum value in dataset
        double * mXmin = nullptr ;

        // maximum value in dataset
        double * mXmax = nullptr ;

        // stepwidth
        double * mStep = nullptr ;

        // inverse stepwidth
        double * mInvStep = nullptr ;

        // the actual data of the set
        double * mValues = nullptr ;

        // coefficients, only for 1D cubic
        double * mCoeffs = nullptr ;

        // work vector for the interpolation
        double * mWork = nullptr ;

        // work vector for the coordinates
        double  * mX = nullptr ;

        // pointer to the function that collects the data
        double
        ( database::*mInterpolate )( const double * x ) ;

//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        // the database constructor (standalone via filename)
        database( const std::string & filename,
                  const std::string & tablename );

//----------------------------------------------------------------------------

        // the database constructor (passing database)
        database( Hdf5File & file ,
                  const std::string & tablename );

//----------------------------------------------------------------------------

        // the destructor
        ~database() ;

//------------------------------------------------------------------------

        void
        load_from_file( const std::string & filename,
                        const std::string & tablename );

//------------------------------------------------------------------------

        void
        load_from_file( Hdf5File          & file,
                        const std::string & tablename );

//----------------------------------------------------------------------------

        /**
         * a dummy function that returns the critical current density
         *
         * @param T        temperature in K
         * @param B        magnetic field strength in T
         * @return
         */
        double
        compute( const double T ) ;

//----------------------------------------------------------------------------

        /**
         * a dummy function that returns the critical current density
         *
         * @param T        temperature in K
         * @param B        magnetic field strength in T
         * @return
         */
        double
        compute( const double T, const double B ) ;

//----------------------------------------------------------------------------

        /**
         * a dummy function that returns the critical current density
         *
         * @param T        temperature in K
         * @param B        magnetic field strength in T
         * @param theta    angle in rad
         * @return
         */
        double
        compute( const double T, const double B, const double theta ) ;

//----------------------------------------------------------------------------

        /**
         * returns the number of dimensions: 1, 2 or 3
         */
        inline uint
        numparams() const
        {
            return mNumberOfDimensions ;
        }

//----------------------------------------------------------------------------

        /**
         * direct access to interpolation function (for material)
         */
        inline double
        interpolate( const double *x )
        {
            return ( this->*mInterpolate )( x );
        }

//------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------

        /**
         * clear the memory
         */
         void
         free();

//----------------------------------------------------------------------------

        void
        distribute_data();

//----------------------------------------------------------------------------

        void
        link_interpolation_function();

//----------------------------------------------------------------------------

        double
        interpolate_line2( const double * x );

//----------------------------------------------------------------------------

        double
        interpolate_line3( const double * x );

//----------------------------------------------------------------------------

        double
        interpolate_line4( const double * x );

//----------------------------------------------------------------------------

        double
        interpolate_quad4( const double * x );

//----------------------------------------------------------------------------

        double
        interpolate_quad9( const double * x );

//----------------------------------------------------------------------------

        double
        interpolate_quad16( const double * x );

//----------------------------------------------------------------------------

        double
        interpolate_hex8( const double * x );

//----------------------------------------------------------------------------

        double
        interpolate_hex27( const double * x );

//----------------------------------------------------------------------------

        double
        interpolate_hex64( const double * x );

//----------------------------------------------------------------------------

        /**
         * the index function for the 2d case
         */
        inline uint
        index( const uint i, const uint j ) const
        {
            return mNumPoints[ 0 ] * j + i ;
        }

//----------------------------------------------------------------------------

        /**
         * the index function for the 3d case
         */
        inline uint
        index( const uint i, const uint j, const uint k ) const
        {
            return mNumPoints[ 0 ] *  ( k * mNumPoints[ 1 ] +  j ) + i ;
        }

//----------------------------------------------------------------------------

        // only for 1d cubic spline
        inline uint
        find_cell_1dspline( const double x )
        {
            long int pivot = std::floor( ( x - mXmin[ 0 ] ) * mInvStep[ 0 ] ) ;

            if( pivot < 0 )
            {
                return 0 ;
            }
            else if ( pivot >= mNumPoints[ 0 ] )
            {
                return 4 * ( mNumPoints[ 0 ] - 1 );
            }
            else
            {
                return 4 * pivot ;
            }
        }

//------------------------------------------------------------------------

        inline uint
        find_cell_linear( const double * x, const uint dimension ) const
        {
            if( x[ dimension ] < mXmin[ dimension ] )
            {
                return 0 ;
            }
            else if ( mXmax[ dimension ] < x[ dimension ] )
            {
                return mNumPoints[ dimension ] - 2 ;
            }
            else
            {
                return uint( (x[ dimension ] - mXmin[ dimension ] ) * mInvStep[ dimension ] );
            }
        }

//------------------------------------------------------------------------

        inline uint
        find_cell_higher_order( const double * x, const unsigned int dimension ) const
        {
            if( x[ dimension ] < mXmin[ dimension ] )
            {
                return 0 ;
            }
            else if ( mXmax[ dimension ] < x[ dimension ] )
            {
                return mNumPoints[ dimension ] - mInterpolationOrder - 1 ;
            }
            else
            {
                return mInterpolationOrder * uint( ( x[ dimension ] - mXmin[ dimension ] ) * mInvStep[ dimension ] );
            }
        }

//------------------------------------------------------------------------

        inline double
        xtoxi_linear( const double * x, const uint index, const unsigned int dimension ) const
        {
            return ( x[ dimension ] + x[ dimension ]
                - mXmin[ dimension ] - mXmin[ dimension ] ) * mInvStep[ dimension ]
                         - index - index - 1 ;
        }

//------------------------------------------------------------------------

        inline double
        xtoxi_quadratic( const double * x, const uint index, const unsigned int dimension ) const
        {
            return ( x[ dimension ] + x[ dimension ]
                     - mXmin[ dimension ] - mXmin[ dimension ] ) * mInvStep[ dimension ]
                   - index - 1 ;
        }

//------------------------------------------------------------------------

        inline double
        xtoxi_higher_order( const double * x, const uint index, const unsigned int dimension ) const
        {
            return 2 * ( x[ dimension ] - mXmin[ dimension ]
                    - mInvOrder * index * mStep[ dimension ] )
                       * mInvStep[ dimension ] - 1 ;
        }


//------------------------------------------------------------------------


        inline double
        interpolate_1dspline( const double * x )
        {
            double X = *x ;

            // get index
            std::size_t k = this->find_cell_1dspline( X );

            return ( (   mCoeffs[ k ]     * X
                       + mCoeffs[ k+1 ] ) * X
                       + mCoeffs[ k+2 ] ) * X
                       + mCoeffs[ k+3 ]  ;
        }
//------------------------------------------------------------------------
    };
//------------------------------------------------------------------------
}



#endif //BELMAT_DATABASE_H

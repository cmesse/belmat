//
// Created by christian on 4/28/22.
//

#ifndef BELMAT_DATABASE_H
#define BELMAT_DATABASE_H

#include <cmath>
#include <string>

namespace belmat
{
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

        // the actual data of the set
        double * mValues = nullptr ;

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

        // the database constructor
        database();

//----------------------------------------------------------------------------

        // the destructor
        ~database() ;

//------------------------------------------------------------------------

        void
        load_from_file( const std::string & filename,
                        const std::string & materialname );

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
        jc( const double T, const double B, const double theta ) ;

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
         link_interpolation_function();

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
                return uint( (x[ dimension ] - mXmin[ dimension ] ) / mStep[ dimension ] );
            }
        }

//------------------------------------------------------------------------

        inline std::size_t
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
                return mInterpolationOrder * std::size_t ( (x[ dimension ] - mXmin[ dimension ] ) / mStep[ dimension ] );
            }
        }

//------------------------------------------------------------------------

        inline double
        xtoxi_linear( const double * x, const uint index, const unsigned int dimension ) const
        {
            return ( x[ dimension ] + x[ dimension ]
                - mXmin[ dimension ] - mXmin[ dimension ] ) / mStep[ dimension ]
                         - index - index - 1 ;
        }

//------------------------------------------------------------------------

        inline double
        xtoxi_quadratic( const double * x, const uint index, const unsigned int dimension ) const
        {
            return ( x[ dimension ] + x[ dimension ]
                     - mXmin[ dimension ] - mXmin[ dimension ] ) / mStep[ dimension ]
                   - index - 1 ;
        }

//------------------------------------------------------------------------

        inline double
        xtoxi_higher_order( const double * x, const uint index, const unsigned int dimension ) const
        {
            return 2 * ( x[ dimension ] - mXmin[ dimension ]
                    - mInvOrder * index * mStep[ dimension ] )
                    / mStep[ dimension ] - 1 ;
        }

//------------------------------------------------------------------------
    };
//------------------------------------------------------------------------
}



#endif //BELMAT_DATABASE_H

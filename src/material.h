//
// Created by Christian Messe on 31.07.22.
//

#ifndef BELMAT_MATERIAL_HPP
#define BELMAT_MATERIAL_HPP

#include <cassert>
#include <string>
#include "database.h"

namespace belmat
{
    class material
    {
        const std::string mLabel ;

        const uint        mRRR ;
        const std::string mCpLabel ;
        const std::string mKLabel ;
        const std::string mRhoLabel ;

        double       mDensity ;

        database * mCp      = nullptr ;
        database * mK       = nullptr ;
        database * mRho     = nullptr ;

        // data containers
        double mKX[ 3 ];
        double mRhoX[ 3 ];

//------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------


        material( const std::string & filename, const std::string & label, const uint rrr=0 );

        ~material();

//------------------------------------------------------------------------

        /**
         * return the label
         */
         const std::string &
         label() const ;

//------------------------------------------------------------------------

        /**
         * return the density for 298.15 K
         */
        double
        density() const ;

//------------------------------------------------------------------------

        /**
         * return the specific heat capacity
         */
        double
        cp( const double T ) ;

//------------------------------------------------------------------------

        /**
         * return the thermal conductivity
         * @param T     in K
         * @param B     in T
         * @param Theta in rad
         */
        double
        k( const double T, const double B=0.0, const double Theta=0.0 ) ;

//------------------------------------------------------------------------

        double
        rho( const double T, const double B=0.0, const double Theta=0.0 ) ;

//------------------------------------------------------------------------

        /**
         * tell if specific heat exisis
         */
        bool
        has_cp() const ;

//------------------------------------------------------------------------

        /**
         * tell if specific heat exisis
         */
        bool
        has_k() const ;

//------------------------------------------------------------------------

        /**
         * tell if specific heat exisis
         */
        bool
        has_rho() const ;

//------------------------------------------------------------------------
    };

//------------------------------------------------------------------------

    inline const std::string &
    material::label() const
    {
        return mLabel ;
    }

//------------------------------------------------------------------------
    inline bool
    material::has_cp() const
    {
        return mCp != nullptr ;
    }

//------------------------------------------------------------------------

    inline bool
    material::has_k() const
    {
        return mK != nullptr ;
    }

//------------------------------------------------------------------------

    inline bool
    material::has_rho() const
    {
        return mRho != nullptr ;
    }

//------------------------------------------------------------------------

    inline double
    material::density() const
    {
        return mDensity;
    }

//------------------------------------------------------------------------

    inline double
    material::cp( const double T )
    {
        assert( this->has_cp() );
        return mCp->interpolate( &T );
    }

//------------------------------------------------------------------------

    inline double
    material::k( const double T, const double B, const double Theta )
    {
        assert( this->has_k() );
        mKX[ 0 ] = T ;
        mKX[ 1 ] = B ;
        mKX[ 2 ] = Theta ;
        return mK->interpolate( mKX );
    }

//------------------------------------------------------------------------

    inline double
    material::rho( const double T, const double B, const double Theta )
    {
        assert( this->has_rho() );
        mRhoX[ 0 ] = T ;
        mRhoX[ 1 ] = B ;
        mRhoX[ 2 ] = Theta ;
        return mRho->interpolate( mRhoX );
    }

//------------------------------------------------------------------------
}
#endif //BELMAT_MATERIAL_HPP

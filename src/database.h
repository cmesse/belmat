//
// Created by christian on 4/28/22.
//

#ifndef BELMAT_DATABASE_H
#define BELMAT_DATABASE_H

namespace belmat
{
//------------------------------------------------------------------------
    class database
    {
//----------------------------------------------------------------------------
    public:
//----------------------------------------------------------------------------

        // the database constructor
        database();

//----------------------------------------------------------------------------

        // the destructor
        ~database() = default ;

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
        jc( const double T, const double B, const double theta ) const ;

//------------------------------------------------------------------------
    };
//------------------------------------------------------------------------
}



#endif //BELMAT_DATABASE_H
